import pandas as pd
from sqlalchemy import create_engine, text
import os

# --- Configuration ---
DB_URL = "postgresql://dduser:devpass@localhost:5432/ddidb"
OUTPUT_DIR = "data/processed"
DRUGS_FILE = os.path.join(OUTPUT_DIR, "drugs_with_smiles.csv")
INTERACTIONS_FILE = os.path.join(OUTPUT_DIR, "positive_interactions.csv")

def main():
    engine = create_engine(DB_URL)
    
    # 1. Fetch drugs with SMILES
    print("Fetching drugs and SMILES...")
    query_drugs = text("SELECT drugbank_id, smiles FROM drugs WHERE smiles IS NOT NULL")
    with engine.connect() as conn:
        df_drugs = pd.read_sql(query_drugs, conn)
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df_drugs.to_csv(DRUGS_FILE, index=False)
    print(f"Saved {len(df_drugs)} drugs to '{DRUGS_FILE}'.")

    # 2. Fetch positive interaction pairs
    print("Fetching positive interaction pairs...")
    query_interactions = text("SELECT drug_a, drug_b FROM interactions")
    with engine.connect() as conn:
        df_positive = pd.read_sql(query_interactions, conn)
    
    df_positive.to_csv(INTERACTIONS_FILE, index=False)
    print(f"Saved {len(df_positive)} interactions to '{INTERACTIONS_FILE}'.")
    
    print("\nExport for Colab complete.")

if __name__ == "__main__":
    main()