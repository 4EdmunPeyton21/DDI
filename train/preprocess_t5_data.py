import pandas as pd
from sqlalchemy import create_engine, text
from sklearn.model_selection import train_test_split
from datasets import Dataset
import json
import os

# --- Configuration ---
DB_URL = "postgresql://dduser:devpass@localhost:5432/ddidb"
OUTPUT_DIR = "data/processed"

def fetch_data_from_db(engine):
    """
    Connects to the database and fetches interaction pairs with their names and evidence.
    """
    print("Connecting to database and fetching interaction data...")
    query = text("""
        SELECT 
            d1.name AS drug_a_name, 
            d2.name AS drug_b_name, 
            i.evidence 
        FROM 
            interactions AS i
        JOIN 
            drugs AS d1 ON i.drug_a = d1.drugbank_id
        JOIN 
            drugs AS d2 ON i.drug_b = d2.drugbank_id
        WHERE 
            d1.name IS NOT NULL AND d2.name IS NOT NULL;
    """)
    with engine.connect() as conn:
        df = pd.read_sql(query, conn)
    print(f"Fetched {len(df)} interaction records from the database.")
    return df

def main():
    """
    Main function to fetch, process, and save the dataset for T5 training.
    """
    engine = create_engine(DB_URL)
    df = fetch_data_from_db(engine)

    # 1. Format the data into "input" and "target" columns
    print("Formatting data into 'input' and 'target' pairs...")
    df['input'] = df['drug_a_name'] + ' + ' + df['drug_b_name']
    
    # Safely extract the description from the JSONB evidence column
    df['target'] = df['evidence'].apply(
        lambda evidence: evidence[0]['description'] if (evidence and isinstance(evidence, list) and len(evidence) > 0) else None
    )
    
    # Drop rows where the target description could not be extracted
    df.dropna(subset=['target'], inplace=True)
    
    # 2. Select the final columns
    final_df = df[['input', 'target']]

    # 3. Split the data into training (80%), validation (10%), and test (10%) sets
    print("Splitting data into train, validation, and test sets...")
    train_df, temp_df = train_test_split(final_df, test_size=0.2, random_state=42)
    val_df, test_df = train_test_split(temp_df, test_size=0.5, random_state=42)

    print(f"Train set size: {len(train_df)}")
    print(f"Validation set size: {len(val_df)}")
    print(f"Test set size: {len(test_df)}")

    # 4. Convert pandas DataFrames to Hugging Face Dataset objects
    train_dataset = Dataset.from_pandas(train_df, preserve_index=False)
    val_dataset = Dataset.from_pandas(val_df, preserve_index=False)
    test_dataset = Dataset.from_pandas(test_df, preserve_index=False)

    # 5. Save the datasets as JSON Lines files
    print("Saving datasets to JSONL files...")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    train_dataset.to_json(os.path.join(OUTPUT_DIR, "train_dataset.jsonl"))
    val_dataset.to_json(os.path.join(OUTPUT_DIR, "val_dataset.jsonl"))
    test_dataset.to_json(os.path.join(OUTPUT_DIR, "test_dataset.jsonl"))

    print("\nPreprocessing complete!")
    print(f"Datasets saved in '{OUTPUT_DIR}'")

if __name__ == "__main__":
    main()