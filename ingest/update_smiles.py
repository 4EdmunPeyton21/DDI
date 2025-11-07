import sys
from lxml import etree
from sqlalchemy import create_engine, text

DB_URL = "postgresql://dduser:devpass@localhost:5432/ddidb"

def main(xml_path):
    """
    Iterates through the XML and runs a targeted UPDATE for each drug with a SMILES string.
    """
    engine = create_engine(DB_URL)
    ns = {'db': 'http://www.drugbank.ca'}
    context = etree.iterparse(xml_path, events=('end',), tag='{http://www.drugbank.ca}drug')
    count = 0
    updated_count = 0

    print("Starting targeted SMILES update...")

    with engine.connect() as conn:
        for event, drug_elem in context:
            count += 1
            if count % 5000 == 0:
                print(f"  ...scanned {count} drugs, found and updated {updated_count} SMILES strings.")

            # Safely get the primary DrugBank ID
            primary_id_elem = drug_elem.find('db:drugbank-id[@primary="true"]', namespaces=ns)
            if primary_id_elem is None or not primary_id_elem.text:
                drug_elem.clear() # Clear memory and skip
                continue
            drugbank_id = primary_id_elem.text
            
            # Safely get the SMILES string
            smiles = None
            properties_elem = drug_elem.find('db:calculated-properties', namespaces=ns)
            if properties_elem is not None:
                for prop in properties_elem.findall('db:property', namespaces=ns):
                    kind_elem = prop.find('db:kind', namespaces=ns)
                    if kind_elem is not None and kind_elem.text == 'SMILES':
                        value_elem = prop.find('db:value', namespaces=ns)
                        if value_elem is not None and value_elem.text:
                            smiles = value_elem.text
                        break
            
            # If a SMILES string was found, execute an immediate UPDATE
            if smiles:
                stmt = text("UPDATE drugs SET smiles = :smiles WHERE drugbank_id = :dbid")
                conn.execute(stmt, {"smiles": smiles, "dbid": drugbank_id})
                conn.commit() # Commit the change immediately
                updated_count += 1

            # Clear memory
            drug_elem.clear()
            while drug_elem.getprevious() is not None:
                del drug_elem.getparent()[0]

    print(f"\nUpdate complete. Scanned {count} drugs and updated {updated_count} SMILES strings.")

if __name__ == "__main__":
    xml_file = sys.argv[1]
    main(xml_file)