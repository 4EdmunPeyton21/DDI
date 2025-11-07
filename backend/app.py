from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from sqlalchemy import create_engine, text
from transformers import T5ForConditionalGeneration, T5Tokenizer
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import json
import os
import torch

# --- Configuration ---
DB_URL = "postgresql://dduser:devpass@localhost:5432/ddidb"
# -- Model 1: Generative (T5) --
T5_MODEL_PATH = "models/biot5_finetuned"
# -- Model 2: Predictive (SMILES) --
SMILES_MODEL_PATH = "models/smiles_MLP_model.joblib"
SMILES_SCALER_PATH = "models/smiles_MLP_scaler.joblib"
FINGERPRINT_SIZE = 2048

# --- Global Resources (Loaded at Startup) ---
engine = None
tokenizer_t5 = None
model_t5 = None
model_smiles = None
scaler_smiles = None
device = "cuda" if torch.cuda.is_available() else "cpu"

# --- Initialize FastAPI app ---
app = FastAPI()

# --- Helper Function for Fingerprinting ---
def get_fingerprint(smiles_string):
    """Converts a SMILES string to a Morgan Fingerprint."""
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None: return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=FINGERPRINT_SIZE)
        return list(fp)
    except:
        return None

# --- Startup Event ---
@app.on_event("startup")
def load_resources():
    global engine, tokenizer_t5, model_t5, model_smiles, scaler_smiles
    
    print("--- SERVER STARTING UP ---")
    
    # 1. Connect to Database
    try:
        engine = create_engine(DB_URL)
        with engine.connect() as conn:
            conn.execute(text("SELECT 1"))
        print("Database connection successful.")
    except Exception as e:
        print(f"WARNING: Database connection failed: {e}")

    # 2. Load Generative T5 Model (The "Historian")
    try:
        print(f"Loading T5 model from '{T5_MODEL_PATH}'...")
        tokenizer_t5 = T5Tokenizer.from_pretrained(T5_MODEL_PATH)
        model_t5 = T5ForConditionalGeneration.from_pretrained(T5_MODEL_PATH).to(device)
        print("T5 model loaded successfully.")
    except Exception as e:
        print(f"FATAL ERROR: Could not load T5 model: {e}")

    # 3. Load Predictive SMILES Model (The "Analyst")
    try:
        print(f"Loading SMILES model from '{SMILES_MODEL_PATH}'...")
        model_smiles = joblib.load(SMILES_MODEL_PATH)
        print("SMILES model loaded successfully.")
        
        print(f"Loading SMILES scaler from '{SMILES_SCALER_PATH}'...")
        scaler_smiles = joblib.load(SMILES_SCALER_PATH)
        print("SMILES scaler loaded successfully.")
    except Exception as e:
        print(f"FATAL ERROR: Could not load SMILES model or scaler: {e}")

    print("--- SERVER READY ---")

# --- API Models ---
class CheckIn(BaseModel):
    drugs: list[str]

class InteractionResponse(BaseModel):
    pair: list[str]
    description: str
    source: str
    probability: float | None = None

# --- API Endpoint ---
@app.post("/check", response_model=list[InteractionResponse])
def check(payload: CheckIn):
    if len(payload.drugs) < 2:
        raise HTTPException(status_code=400, detail="Please provide at least two drug names.")

    results = []
    
    with engine.connect() as conn:
        for i in range(len(payload.drugs)):
            for j in range(i + 1, len(payload.drugs)):
                name_a = payload.drugs[i]
                name_b = payload.drugs[j]
                
                print(f"Checking database for: {name_a} + {name_b}")
                
                # Query for both drugs at once
                db_rows = conn.execute(
                    text("SELECT drugbank_id, name, smiles FROM drugs WHERE name = :name_a OR name = :name_b"),
                    [{"name_a": name_a, "name_b": name_b}]
                ).fetchall()
                
                # --- THIS IS THE FIX ---
                # Use the correct index (row[1] is the name) as the dictionary key
                id_map = {row[1]: {'id': row[0], 'smiles': row[2]} for row in db_rows}
                
                info_a = id_map.get(name_a)
                info_b = id_map.get(name_b)

                if not info_a or not info_b:
                    print("  -> One or both drugs not found in database. Skipping pair.")
                    continue
                
                id_a = info_a['id']
                id_b = info_b['id']

                interaction = conn.execute(
                    text("SELECT * FROM interactions WHERE (drug_a = :a AND drug_b = :b) OR (drug_a = :b AND drug_b = :a)"),
                    [{"a": id_a, "b": id_b}]
                ).fetchone()
                
                # --- 1. Database Check (Highest Confidence) ---
                if interaction:
                    print("  -> Found in database.")
                    evidence_list = interaction.evidence if isinstance(interaction.evidence, list) else []
                    description = evidence_list[0].get('description', 'Interaction found.') if evidence_list else 'Interaction found.'
                    
                    results.append(InteractionResponse(
                        pair=[name_a, name_b],
                        description=description,
                        source="DrugBank (Curated Database)",
                        probability=1.0
                    ))
                    continue 

                # --- 2. Predictive Model Check (The "Analyst") ---
                print("  -> Not in database. Checking predictive model...")
                if not model_smiles or not scaler_smiles:
                    print("  -> Predictive model not loaded. Skipping.")
                    continue

                fp_a = get_fingerprint(info_a['smiles'])
                fp_b = get_fingerprint(info_b['smiles'])
                
                if not fp_a or not fp_b:
                    print("  -> Could not generate fingerprints (missing SMILES?). Skipping.")
                    continue
                
                features = np.array(fp_a + fp_b, dtype=np.float32).reshape(1, -1)
                features_scaled = scaler_smiles.transform(features)
                
                proba = model_smiles.predict_proba(features_scaled)[0][1] # Get probability for class 1
                
                if proba < 0.5: # 50% threshold
                    print(f"  -> Model predicts NO interaction (Prob: {proba*100:.2f}%).")
                    continue
                
                print(f"  -> Model predicts YES interaction (Prob: {proba*100:.2f}%).")

                # --- 3. Generative Model Check (The "Historian") ---
                print("  -> Asking T5 model to describe the interaction...")
                if not model_t5 or not tokenizer_t5:
                    print("  -> T5 model not loaded. Skipping description.")
                    results.append(InteractionResponse(
                        pair=[name_a, name_b],
                        description="Interaction predicted based on chemical structure, but description is unavailable.",
                        source="AI Prediction (SMILES Model)",
                        probability=proba
                    ))
                    continue

                input_text = f"describe interaction: {name_a} + {name_b}"
                inputs = tokenizer_t5(input_text, return_tensors="pt", max_length=128, truncation=True).to(device)
                
                with torch.no_grad():
                    outputs = model_t5.generate(
                        **inputs, max_length=256, num_beams=4, early_stopping=True
                    )
                generated_text = tokenizer_t5.decode(outputs[0], skip_special_tokens=True)
                
                results.append(InteractionResponse(
                    pair=[name_a, name_b],
                    description=generated_text,
                    source="AI Prediction (Hybrid Model)",
                    probability=proba
                ))

    return results

@app.get("/")
def read_root():
    return {"status": "DDI Detector Hybrid API is running"}