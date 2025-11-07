import sys
from lxml import etree

# The DrugBank ID for Acetylsalicylic acid
TARGET_DRUG_ID = 'DB00945'

def main(xml_path):
    """
    Finds a single drug in the XML and prints its extracted data and raw XML block.
    """
    ns = {'db': 'http://www.drugbank.ca'}
    context = etree.iterparse(xml_path, events=('end',), tag='{http://www.drugbank.ca}drug')

    print(f"Searching for drug {TARGET_DRUG_ID} in '{xml_path}'...")

    for event, drug_elem in context:
        # Find the primary ID for the current drug
        primary_id_elem = drug_elem.find('db:drugbank-id[@primary="true"]', namespaces=ns)
        if primary_id_elem is None:
            continue
        
        current_drug_id = primary_id_elem.text
        
        # If we found our target drug, process it
        if current_drug_id == TARGET_DRUG_ID:
            print(f"\n--- FOUND {TARGET_DRUG_ID} ---")
            
            # 1. Run the extraction logic
            name_elem = drug_elem.find('db:name', namespaces=ns)
            name = name_elem.text if name_elem is not None else "Name not found"
            
            smiles = "SMILES not found"
            properties_elem = drug_elem.find('db:calculated-properties', namespaces=ns)
            if properties_elem is not None:
                for prop in properties_elem.findall('db:property', namespaces=ns):
                    kind_elem = prop.find('db:kind', namespaces=ns)
                    kind = kind_elem.text if kind_elem is not None else ""
                    if kind == 'SMILES':
                        value_elem = prop.find('db:value', namespaces=ns)
                        smiles = value_elem.text if value_elem is not None else "SMILES value is empty"
                        break
            
            print("\n[1] EXTRACTED DATA:")
            print(f"  - DrugBank ID: {current_drug_id}")
            print(f"  - Name: {name}")
            print(f"  - SMILES: {smiles}")
            
            # 2. Print the raw XML block for calculated properties
            print("\n[2] RAW XML FOR 'calculated-properties':")
            if properties_elem is not None:
                # Pretty-print the XML block to make it readable
                print(etree.tostring(properties_elem, pretty_print=True).decode('utf-8'))
            else:
                print("  -> 'calculated-properties' section not found in XML.")
            
            break # Stop after finding the drug

    print("\nDebug script finished.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python debug_smiles.py <path_to_drugbank.xml>")
        sys.exit(1)
        
    xml_file = sys.argv[1]
    main(xml_file)