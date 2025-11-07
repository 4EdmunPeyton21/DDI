import pandas as pd
import numpy as np
from sqlalchemy import create_engine, text
from rdkit import Chem
from rdkit.Chem import AllChem
import random
import os

# --- Configuration ---
DB_URL = "postgresql://dduser:devpass@localhost:5432/ddidb"
OUTPUT_DIR = "data/processed"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "structure_dataset.parquet")
FINGERPRINT_SIZE = 2048 # Standard size for Morgan fingerprints

def get_fingerprint(smiles_string):
    """Converts a SMILES string to a Morgan Fingerprint."""
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None: return None
        # Get a 2048-bit fingerprint
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=FINGERPRINT_SIZE)
        return list(fp.ToBitString()) # Convert to a list of '0's and '1's
    except:
        return None

def main():
    engine = create_engine(DB_URL)
    
    # 1. Fetch all drugs with SMILES and create fingerprints
    print("Fetching drugs and generating fingerprints...")
    query_drugs = text("SELECT drugbank_id, smiles FROM drugs WHERE smiles IS NOT NULL")
    with engine.connect() as conn:
        df_drugs = pd.read_sql(query_drugs, conn)

    df_drugs['fingerprint'] = df_drugs['smiles'].apply(get_fingerprint)
    # Drop drugs where fingerprint generation failed
    df_drugs.dropna(subset=['fingerprint'], inplace=True)
    
    # Create a fast lookup dictionary for fingerprints
    fp_dict = pd.Series(df_drugs.fingerprint.values, index=df_drugs.drugbank_id).to_dict()
    all_drug_ids = set(df_drugs['drugbank_id'])
    print(f"Successfully generated fingerprints for {len(fp_dict)} drugs.")

    # 2. Fetch positive interaction pairs
    print("Fetching positive interaction pairs...")
    query_interactions = text("SELECT drug_a, drug_b FROM interactions")
    with engine.connect() as conn:
        df_positive = pd.read_sql(query_interactions, conn)
    
    # Filter out pairs if one of the drugs doesn't have a fingerprint
    df_positive = df_positive[
        df_positive['drug_a'].isin(all_drug_ids) & df_positive['drug_b'].isin(all_drug_ids)
    ]
    df_positive['label'] = 1
    num_positives = len(df_positive)
    print(f"Found {num_positives} valid positive interaction pairs.")

    # 3. Generate negative interaction pairs
    print("Generating negative interaction pairs...")
    positive_pairs = set(zip(df_positive['drug_a'], df_positive['drug_b']))
    drug_id_list = list(all_drug_ids)
    negative_pairs = set()

    # Create a balanced dataset by generating an equal number of negative pairs
    while len(negative_pairs) < num_positives:
        d1 = random.choice(drug_id_list)
        d2 = random.choice(drug_id_list)
        if d1 == d2: continue
        
        # Ensure consistent order (A,B) vs (B,A)
        pair = tuple(sorted((d1, d2)))
        
        # Add to set if it's not a known positive interaction
        if pair not in positive_pairs and pair not in negative_pairs:
            negative_pairs.add(pair)
            
    df_negative = pd.DataFrame(list(negative_pairs), columns=['drug_a', 'drug_b'])
    df_negative['label'] = 0
    print(f"Generated {len(df_negative)} negative interaction pairs.")

    # 4. Combine datasets
    print("Combining positive and negative datasets...")
    final_dataset_df = pd.concat([df_positive, df_negative], ignore_index=True)
    
    # 5. Map fingerprints to pairs
    print("Mapping fingerprints to interaction pairs...")
    # This step is slow but necessary. We'll create two new columns.
    final_dataset_df['fp_a'] = final_dataset_df['drug_a'].map(fp_dict)
    final_dataset_df['fp_b'] = final_dataset_df['drug_b'].map(fp_dict)
    
    # Concatenate the two fingerprints to create the final feature vector
    # This creates a 4096-bit vector (2048 + 2048)
    final_dataset_df['features'] = final_dataset_df['fp_a'] + final_dataset_df['fp_b']
    
    # 6. Save final dataset
    final_df_to_save = final_dataset_df[['features', 'label']]
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    # Using Parquet is much more efficient for this type of data than CSV
    final_df_to_save.to_parquet(OUTPUT_FILE, index=False)
    
    print(f"\nSuccessfully created and saved structural dataset to '{OUTPUT_FILE}'.")
    print(f"Total examples: {len(final_df_to_save)} (50% positive, 50% negative)")

if __name__ == "__main__":
    main()