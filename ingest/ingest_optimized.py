import sys
import json
from lxml import etree
from sqlalchemy import create_engine, text

# --- Configuration ---
DB_URL = "postgresql://dduser:devpass@localhost:5432/ddidb"
BATCH_SIZE = 1000

def parse_drugbank_optimized(xml_path):
    """
    Parses the DrugBank XML iteratively, yielding one drug lxml element at a time.
    This is memory-efficient for very large files.
    """
    # The namespace is required by lxml to find tags in an XML with a default namespace.
    ns = {'db': 'http://www.drugbank.ca'}
    context = etree.iterparse(xml_path, events=('end',), tag='{http://www.drugbank.ca}drug')
    for event, elem in context:
        yield elem, ns # Yield both the element and the namespace map
        # Clear the element from memory after it's processed.
        elem.clear()
        # Also clear its preceding siblings to free up more memory.
        while elem.getprevious() is not None:
            del elem.getparent()[0]

def extract_definitive(drug_elem, ns):
    """
    Extracts a complete and robust record from a single lxml drug element,
    with paranoid checks for missing data to prevent silent failures.
    """
    drugbank_id = None
    name = None
    synonyms = []
    smiles = None

    # --- Safely extract DrugBank ID ---
    # First, try to find the ID marked as primary.
    primary_id_elem = drug_elem.find('db:drugbank-id[@primary="true"]', namespaces=ns)
    if primary_id_elem is not None and primary_id_elem.text:
        drugbank_id = primary_id_elem.text
    else:
        # If no primary is found, fall back to the very first ID listed.
        first_id_elem = drug_elem.find('db:drugbank-id', namespaces=ns)
        if first_id_elem is not None and first_id_elem.text:
            drugbank_id = first_id_elem.text

    # --- Safely extract the drug name ---
    name_elem = drug_elem.find('db:name', namespaces=ns)
    if name_elem is not None and name_elem.text:
        name = name_elem.text
    
    # --- Safely extract all synonyms ---
    synonyms_elems = drug_elem.findall('db:synonyms/db:synonym', namespaces=ns)
    if synonyms_elems:
        # Ensure we only add synonyms that actually have text content.
        synonyms = [s.text for s in synonyms_elems if s.text]
    
    # --- Safely and correctly extract the SMILES string ---
    properties_elem = drug_elem.find('db:calculated-properties', namespaces=ns)
    if properties_elem is not None:
        # Loop through all properties to find the one for SMILES.
        for prop in properties_elem.findall('db:property', namespaces=ns):
            kind_elem = prop.find('db:kind', namespaces=ns)
            # Check if the 'kind' element exists and its text is 'SMILES'.
            if kind_elem is not None and kind_elem.text == 'SMILES':
                value_elem = prop.find('db:value', namespaces=ns)
                # If the 'value' element also exists and has text, we've found our SMILES string.
                if value_elem is not None and value_elem.text:
                    smiles = value_elem.text
                break # We found it, so we can stop looping through properties.

    return dict(
        drugbank_id=drugbank_id, name=name, synonyms=synonyms, smiles=smiles,
        properties={}, raw={} # Placeholders for now.
    )

def upsert_drug_batch(engine, batch):
    """Inserts or updates a batch of drug records in a single database transaction."""
    if not batch: return

    with engine.begin() as conn:
        stmt = text("""
            INSERT INTO drugs (drugbank_id, name, synonyms, smiles, properties, raw)
            VALUES (:drugbank_id, :name, :synonyms, :smiles, :properties, :raw)
            ON CONFLICT (drugbank_id) DO UPDATE SET 
                name = EXCLUDED.name, 
                synonyms = EXCLUDED.synonyms,
                smiles = EXCLUDED.smiles; -- Ensure smiles is also updated.
        """)
        conn.execute(stmt, batch)

def main(xml_path):
    """Main function to run the definitive, corrected ingestion process."""
    engine = create_engine(DB_URL)
    batch = []
    count = 0
    
    print(f"Starting definitive corrected ingestion of '{xml_path}'...")
    
    for drug_elem, ns in parse_drugbank_optimized(xml_path):
        rec = extract_definitive(drug_elem, ns)
        # Only process records that have a valid drugbank_id.
        if rec and rec['drugbank_id']:
            params = {
                "drugbank_id": rec['drugbank_id'], "name": rec['name'],
                "synonyms": rec['synonyms'], "smiles": rec['smiles'],
                "properties": json.dumps(rec['properties']), "raw": json.dumps(rec['raw'])
            }
            batch.append(params)
            count += 1
            
            if len(batch) >= BATCH_SIZE:
                upsert_drug_batch(engine, batch)
                print(f"  ...processed {count} drugs")
                batch = []

    # Don't forget to save the last, partially filled batch.
    if batch:
        upsert_drug_batch(engine, batch)
    
    print(f"Ingestion complete. Total drugs processed: {count}.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python ingest/ingest_optimized.py <path_to_drugbank.xml>")
        sys.exit(1)
        
    xml_file = sys.argv[1]
    main(xml_file)