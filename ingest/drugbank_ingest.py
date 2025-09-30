import xmltodict
import json
import sys
from sqlalchemy import create_engine, text

# Database connection URL
DB_URL = "postgresql://dduser:devpass@localhost:5432/ddidb"

def parse_drugbank(xml_path, max_items=None):
    """A generator to parse the DrugBank XML file incrementally."""
    with open(xml_path, "r", encoding="utf-8") as f:
        doc = xmltodict.parse(f.read())
    
    for i, drug in enumerate(doc['drugbank']['drug']):
        if max_items and i >= max_items:
            break
        yield drug

def extract_simple(drug):
    """Extracts a simplified record from a single DrugBank drug entry."""
    drugbank_id = None
    # Handle the drugbank-id field, which can be a list or a single dictionary
    id_field = drug.get('drugbank-id')
    if isinstance(id_field, list):
        for d in id_field:
            if isinstance(d, dict) and d.get('@primary') == 'true':
                drugbank_id = d.get('#text')
                break
    elif isinstance(id_field, dict):
        drugbank_id = id_field.get('#text')

    if not drugbank_id and id_field:
        if isinstance(id_field, list) and id_field[0]:
            drugbank_id = id_field[0].get('#text')

    name = drug.get('name')
    
    # Robustly handle synonyms, which can be a list of dicts, a single dict, or a single string
    synonyms = []
    synonyms_field = drug.get('synonyms', {}).get('synonym')
    if synonyms_field:
        if isinstance(synonyms_field, list):
            for s in synonyms_field:
                if isinstance(s, dict):
                    synonyms.append(s.get('#text'))
                elif isinstance(s, str):
                    synonyms.append(s)
        elif isinstance(synonyms_field, dict):
            synonyms.append(synonyms_field.get('#text'))
        elif isinstance(synonyms_field, str):
            synonyms.append(synonyms_field)

    # Extract SMILES string from calculated properties
    smiles = None
    if drug.get('calculated-properties'):
        for prop in drug.get('calculated-properties', {}).get('property', []):
            if isinstance(prop, dict) and prop.get('kind') == 'SMILES':
                smiles = prop.get('value')
                break
    
    raw = drug
    properties = {}

    return dict(
        drugbank_id=drugbank_id,
        name=name,
        synonyms=synonyms,
        smiles=smiles,
        properties=properties,
        raw=raw
    )

def upsert_drug(engine, rec):
    """Inserts or updates a drug record in the database."""
    with engine.begin() as conn:
        stmt = text("""
            INSERT INTO drugs (drugbank_id, name, synonyms, smiles, properties, raw)
            VALUES (:drugbank_id, :name, :synonyms, :smiles, :properties, :raw)
            ON CONFLICT (drugbank_id) DO UPDATE SET 
                name = EXCLUDED.name, 
                synonyms = EXCLUDED.synonyms, 
                smiles = EXCLUDED.smiles, 
                properties = EXCLUDED.properties, 
                raw = EXCLUDED.raw
        """)
        
        # Prepare the record for insertion, converting dicts to JSON strings
        params = {
            "drugbank_id": rec['drugbank_id'],
            "name": rec['name'],
            "synonyms": rec['synonyms'],
            "smiles": rec['smiles'],
            "properties": json.dumps(rec['properties']),
            "raw": json.dumps(rec['raw'])
        }
        
        conn.execute(stmt, params)

def main(xml_path, max_items=None):
    """Main function to run the ingestion process."""
    engine = create_engine(DB_URL)
    
    print(f"Starting ingestion of '{xml_path}'...")
    count = 0
    for drug in parse_drugbank(xml_path, max_items=int(max_items) if max_items else None):
        rec = extract_simple(drug)
        if rec and rec['drugbank_id']:
            upsert_drug(engine, rec)
            count += 1
            if count % 50 == 0:
                print(f"  ...processed {count} drugs")
    
    print(f"Ingestion complete. Total drugs processed: {count}.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python ingest/drugbank_ingest.py <path_to_drugbank.xml> [max_items]")
        sys.exit(1)
        
    xml_file = sys.argv[1]
    max_items_arg = sys.argv[2] if len(sys.argv) > 2 else None
    main(xml_file, max_items_arg)