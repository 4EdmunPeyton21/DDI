import sys
import json
from lxml import etree
from sqlalchemy import create_engine, text

# --- Configuration ---
DB_URL = "postgresql://dduser:devpass@localhost:5432/ddidb"
BATCH_SIZE = 1000  # Process 1000 drugs per database transaction

# --- Helper function to convert XML element to a simpler dict ---
def elem_to_dict(elem):
    """Converts an lxml element to a dictionary."""
    d = {}
    for child in elem:
        # Extract tag name without the namespace
        tag = etree.QName(child).localname
        if child.text and child.text.strip():
            value = child.text.strip()
            # Handle attributes if they exist
            if child.attrib:
                value = {'#text': value, **child.attrib}
            
            if tag in d:
                if not isinstance(d[tag], list):
                    d[tag] = [d[tag]]
                d[tag].append(value)
            else:
                d[tag] = value
    return d

def parse_drugbank_optimized(xml_path):
    """
    Parses the DrugBank XML iteratively, yielding one drug at a time
    without loading the entire file into memory.
    """
    context = etree.iterparse(xml_path, events=('end',), tag='{http://www.drugbank.ca}drug')
    for event, elem in context:
        yield elem_to_dict(elem)
        # Clear the element and its predecessors to free up memory
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]

def extract_optimized(drug_dict):
    """Extracts a simplified record from a dictionary created from an XML element."""
    drugbank_id = None
    id_field = drug_dict.get('drugbank-id')

    if isinstance(id_field, list):
        # It's a list, so we look for the primary ID
        for item in id_field:
            if isinstance(item, dict) and item.get('@primary') == 'true':
                drugbank_id = item.get('#text')
                break
        # Fallback if no primary is found
        if not drugbank_id:
            drugbank_id = id_field[0].get('#text') if isinstance(id_field[0], dict) else id_field[0]
    elif isinstance(id_field, dict):
        # It's a single dictionary
        drugbank_id = id_field.get('#text')
    elif isinstance(id_field, str):
        # It's a plain string
        drugbank_id = id_field

    name = drug_dict.get('name')
    
    synonyms_field = drug_dict.get('synonym', [])
    if isinstance(synonyms_field, list):
        synonyms = [s if isinstance(s, str) else s.get('#text') for s in synonyms_field]
    else: # Handle single synonym case
        synonyms = [synonyms_field if isinstance(synonyms_field, str) else synonyms_field.get('#text')]

    smiles = None
    # This logic can be expanded to parse the nested calculated-properties
    # For now, we'll leave it as a placeholder to keep the script simpler.

    return dict(
        drugbank_id=drugbank_id,
        name=name,
        synonyms=synonyms,
        smiles=smiles,
        properties={},
        raw={} # Raw JSON is omitted for performance
    )
    
def upsert_drug_batch(engine, batch):
    """Inserts or updates a batch of drug records in a single transaction."""
    if not batch:
        return

    with engine.begin() as conn:
        stmt = text("""
            INSERT INTO drugs (drugbank_id, name, synonyms, smiles, properties, raw)
            VALUES (:drugbank_id, :name, :synonyms, :smiles, :properties, :raw)
            ON CONFLICT (drugbank_id) DO UPDATE SET 
                name = EXCLUDED.name, 
                synonyms = EXCLUDED.synonyms
        """)
        # Note: We're only updating a subset of fields for simplicity
        
        conn.execute(stmt, batch)

def main(xml_path):
    """Main function to run the optimized ingestion process."""
    engine = create_engine(DB_URL)
    batch = []
    count = 0
    
    print(f"Starting optimized ingestion of '{xml_path}'...")
    
    for drug_dict in parse_drugbank_optimized(xml_path):
        rec = extract_optimized(drug_dict)
        if rec and rec['drugbank_id']:
            params = {
                "drugbank_id": rec['drugbank_id'],
                "name": rec['name'],
                "synonyms": rec['synonyms'],
                "smiles": rec['smiles'],
                "properties": json.dumps(rec['properties']),
                "raw": json.dumps(rec['raw'])
            }
            batch.append(params)
            count += 1
            
            if len(batch) >= BATCH_SIZE:
                upsert_drug_batch(engine, batch)
                print(f"  ...processed {count} drugs")
                batch = []

    # Insert any remaining records in the last batch
    if batch:
        upsert_drug_batch(engine, batch)
    
    print(f"Ingestion complete. Total drugs processed: {count}.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python ingest/ingest_optimized.py <path_to_drugbank.xml>")
        sys.exit(1)
        
    xml_file = sys.argv[1]
    main(xml_file)