from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from sqlalchemy import create_engine, text
import os
import requests 

# --- Configuration ---
# Vercel automatically sets POSTGRES_URL
DB_URL = os.getenv("POSTGRES_URL") 
if DB_URL and DB_URL.startswith("postgres://"):
    DB_URL = DB_URL.replace("postgres://", "postgresql://", 1)

# Hugging Face API (Replaces local model)
HF_API_URL = "https://api-inference.huggingface.co/models/google/flan-t5-base"
HF_API_KEY = os.getenv("HF_API_KEY")

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- Database Connection ---
engine = None
try:
    if DB_URL:
        # SSL is required for Vercel/Neon
        engine = create_engine(DB_URL, connect_args={"sslmode": "require"})
    else:
        print("Warning: No Database URL found.")
except Exception as e:
    print(f"DB Error: {e}")

# --- Models ---
class DrugRequest(BaseModel):
    drugs: list[str]

class InteractionResponse(BaseModel):
    pair: list[str]
    description: str
    source: str
    probability: float = 1.0

# --- Endpoints ---
@app.get("/")
def read_root():
    return {"status": "Vercel API is Live"}

@app.get("/autocomplete")
def autocomplete(query: str = ""):
    if not query or len(query.strip()) < 2: return []
    if not engine: return []

    try:
        with engine.connect() as conn:
            # Fast, case-insensitive search
            sql = text("""
                SELECT DISTINCT name FROM drugs 
                WHERE name ILIKE :p 
                ORDER BY name ASC LIMIT 10
            """)
            result = conn.execute(sql, {"p": f"{query}%"}).fetchall()
            return [row[0] for row in result]
    except Exception as e:
        print(f"Autocomplete Error: {e}")
        return []

@app.post("/check", response_model=list[InteractionResponse])
def check(payload: DrugRequest):
    if len(payload.drugs) < 2:
        raise HTTPException(status_code=400, detail="Provide at least two drugs.")

    results = []
    
    # 1. DATABASE CHECK
    if engine:
        try:
            with engine.connect() as conn:
                for i in range(len(payload.drugs)):
                    for j in range(i + 1, len(payload.drugs)):
                        d1, d2 = payload.drugs[i], payload.drugs[j]
                        
                        # Look for interaction in either direction (A->B or B->A)
                        query = text("""
                            SELECT evidence FROM interactions 
                            WHERE (drug_a_name ILIKE :a AND drug_b_name ILIKE :b)
                            OR (drug_a_name ILIKE :b AND drug_b_name ILIKE :a)
                            LIMIT 1
                        """)
                        row = conn.execute(query, {"a": d1, "b": d2}).fetchone()
                        
                        if row:
                            # Parse JSON evidence if it exists
                            desc = "Interaction found."
                            if row[0]:
                                # Check if it's a dict or list (depends on ingestion)
                                if isinstance(row[0], dict):
                                    desc = row[0].get('description', desc)
                                elif isinstance(row[0], str): # Sometimes stored as stringified JSON
                                     desc = row[0]
                            
                            results.append(InteractionResponse(
                                pair=[d1, d2],
                                description=desc,
                                source="Database (Vercel Postgres)",
                                probability=1.0
                            ))
                            continue
                        
                        # 2. AI PREDICTION (API Call)
                        # Only call AI if NO database record found
                        if HF_API_KEY:
                            try:
                                response = requests.post(
                                    HF_API_URL,
                                    headers={"Authorization": f"Bearer {HF_API_KEY}"},
                                    json={"inputs": f"describe interaction: {d1} and {d2}"}
                                )
                                if response.status_code == 200:
                                    data = response.json()
                                    if isinstance(data, list) and len(data) > 0:
                                        generated = data[0].get('generated_text', "Potential interaction detected.")
                                        results.append(InteractionResponse(
                                            pair=[d1, d2],
                                            description=generated,
                                            source="AI Prediction (Hugging Face API)",
                                            probability=0.85 # Placeholder confidence
                                        ))
                            except:
                                pass 
        except Exception as e:
            print(f"Query Error: {e}")

    return results