import sys
import json
from lxml import etree
from sqlalchemy import create_engine, text

# --- Configuration ---
DB_URL = "postgresql://dduser:devpass@localhost:5432/ddidb"
BATCH_SIZE = 1000  # Process 1000 interactions per database transaction

# --- Helper function to convert XML element to a simpler dict ---
def elem_to_dict(elem):
    d = {}
    for child in elem:
        tag = etree.QName(child).localname
        if child.text and child.text.strip():
            value = child.text.strip()
            d[tag] = value
    return d

def parse_drugbank_optimized(xml_path):
    """Parses the DrugBank XML iteratively, yielding one drug element at a time."""
    context = etree.iterparse(xml_path, events=('end',), tag='{http://www.drugbank.ca}drug')
    for event, elem in context:
        yield elem
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]

def extract_interactions(drug_elem):
    """A generator that yields interaction records from a single drug element."""
    # Find the primary DrugBank ID for the source drug
    primary_id_elem = drug_elem.find('{http://www.drugbank.ca}drugbank-id[@primary="true"]')
    if primary_id_elem is None:
        return
    source_drug_id = primary_id_elem.text

    # Find the interactions section
    interactions_elem = drug_elem.find('{http://www.drugbank.ca}drug-interactions')
    if interactions_elem is None:
        return

    # Yield each interaction found
    for interaction in interactions_elem.findall('{http://www.drugbank.ca}drug-interaction'):
        interaction_dict = elem_to_dict(interaction)
        yield {
            'source_drug_id': source_drug_id,
            'target_drug_id': interaction_dict.get('drugbank-id'),
            'description': interaction_dict.get('description')
        }

def upsert_interactions_batch(engine, batch):
    """Inserts or updates a batch of interaction records."""
    if not batch:
        return

    with engine.begin() as conn:
        stmt = text("""
            INSERT INTO interactions (drug_a, drug_b, evidence)
            VALUES (:drug_a, :drug_b, :evidence)
            ON CONFLICT (drug_a, drug_b) DO UPDATE SET
                evidence = interactions.evidence || CAST(:evidence AS JSONB);
        """)
        # Note: The conflict update now explicitly casts the new evidence to JSONB
        # to ensure the || operator works correctly.
        
        conn.execute(stmt, batch)
        
def main(xml_path):
    """Main function to run the interactions ingestion process."""
    engine = create_engine(DB_URL)
    batch = []
    count = 0
    
    print(f"Starting interactions ingestion from '{xml_path}'...")
    
    for drug_elem in parse_drugbank_optimized(xml_path):
        for interaction in extract_interactions(drug_elem):
            # Ensure consistent order to avoid duplicate pairs (A,B) and (B,A)
            drug_a = min(interaction['source_drug_id'], interaction['target_drug_id'])
            drug_b = max(interaction['source_drug_id'], interaction['target_drug_id'])

            params = {
                'drug_a': drug_a,
                'drug_b': drug_b,
                'evidence': json.dumps([{'source': 'DrugBank', 'description': interaction['description']}])
            }
            batch.append(params)
            count += 1
            
            if len(batch) >= BATCH_SIZE:
                upsert_interactions_batch(engine, batch)
                print(f"  ...processed {count} interactions")
                batch = []

    if batch:
        upsert_interactions_batch(engine, batch)
    
    print(f"Ingestion complete. Total interactions found: {count}.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python ingest/interactions_ingest.py <path_to_drugbank.xml>")
        sys.exit(1)
        
    xml_file = sys.argv[1]
    main(xml_file)