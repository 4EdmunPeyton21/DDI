from fastapi import FastAPI
from pydantic import BaseModel
from sqlalchemy import create_engine, text
import json
import os

# --- Configuration ---
DB_URL = "postgresql://dduser:devpass@localhost:5432/ddidb"
NAME_INDEX_PATH = "data/processed/name_to_drugbank.json"

# --- Initialize Database Engine and FastAPI app ---
engine = create_engine(DB_URL)
app = FastAPI()

# --- Load the name index ---
print("--- SERVER STARTING UP ---")
if not os.path.exists(NAME_INDEX_PATH):
    print(f"FATAL ERROR: Name index file not found at {NAME_INDEX_PATH}")
    name_to_id_map = {}
else:
    with open(NAME_INDEX_PATH, "r", encoding="utf-8") as f:
        name_to_id_map = json.load(f)
    print(f"Successfully loaded name index with {len(name_to_id_map)} entries.")
print("--- SERVER READY ---")


class CheckIn(BaseModel):
    drugs: list[str]

@app.post("/check")
def check(payload: CheckIn):
    """
    Checks for interactions between a list of provided drug names.
    """
    print(f"Received request for drugs: {payload.drugs}")
    
    # Use lowercase for case-insensitive lookup
    drug_ids = [name_to_id_map.get(name.lower()) for name in payload.drugs]
    print(f"Mapped to DrugBank IDs: {drug_ids}")
    
    results = []
    
    # Filter out any names that weren't found
    valid_ids = [id for id in drug_ids if id is not None]
    
    if len(valid_ids) < 2:
        print("Result: Could not find at least two valid drugs. Sending error response.")
        return {"message": "Could not find at least two valid drugs to check."}

    # Create a map from original names to their found IDs to handle the payload.drugs indexing correctly
    original_names_map = {}
    for i, drug_id in enumerate(drug_ids):
        if drug_id:
            original_names_map[drug_id] = payload.drugs[i]

    with engine.connect() as conn:
        for i in range(len(valid_ids)):
            for j in range(i + 1, len(valid_ids)):
                id_a = valid_ids[i]
                id_b = valid_ids[j]

                drug_a = min(id_a, id_b)
                drug_b = max(id_a, id_b)
                
                print(f"Checking interaction between {drug_a} and {drug_b}...")

                interaction = conn.execute(
                    text("SELECT * FROM interactions WHERE drug_a = :da AND drug_b = :db"),
                    {"da": drug_a, "db": drug_b}
                ).fetchone()

                if interaction:
                    print(f"  -> Interaction FOUND!")
                    # Directly access 'interaction.evidence' which is already a Python list
                    evidence_list = interaction.evidence if isinstance(interaction.evidence, list) else []
                    description = evidence_list[0].get('description', 'No description available.') if evidence_list else 'No description available.'
                    
                    # Use the map to get the original names for the pair
                    name_a = original_names_map.get(id_a)
                    name_b = original_names_map.get(id_b)

                    results.append({
                        "pair": [name_a, name_b],
                        "description": description
                    })
                else:
                    print(f"  -> No interaction found in DB.")

    if not results:
        print("\nResult: No interactions were found for any pairs. Sending 'not found' response.")
        return {"message": "No interactions found for the given drugs."}

    print("\nResult: Found interactions. Sending successful response.")
    return {"interactions": results}