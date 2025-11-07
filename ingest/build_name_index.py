import json
from sqlalchemy import create_engine, text

DB_URL = "postgresql://dduser:devpass@localhost:5432/ddidb"
OUTPUT_PATH = "data/processed/name_to_drugbank.json"

def main():
    engine = create_engine(DB_URL)
    
    print("Connecting to database to fetch all drug names and synonyms...")
    with engine.connect() as conn:
        # Fetch the primary name and the list of synonyms
        result = conn.execute(text("SELECT name, synonyms, drugbank_id FROM drugs WHERE name IS NOT NULL"))
        rows = result.fetchall()

    print(f"Found {len(rows)} drugs. Building comprehensive index...")
    
    name_map = {}
    for name, synonyms, drugbank_id in rows:
        # Add the primary name (and its lowercase version)
        if name:
            name_map[name] = drugbank_id
            name_map[name.lower()] = drugbank_id
        
        # Add all synonyms (and their lowercase versions)
        if synonyms:
            for synonym in synonyms:
                if synonym:
                    name_map[synonym] = drugbank_id
                    name_map[synonym.lower()] = drugbank_id
    
    with open(OUTPUT_PATH, "w", encoding="utf-8") as f:
        json.dump(name_map, f, indent=2)
        
    print(f"Successfully created comprehensive name index with {len(name_map)} entries at: {OUTPUT_PATH}")

if __name__ == "__main__":
    main()