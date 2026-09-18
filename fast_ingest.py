import csv
import os
import sys
import psycopg2
from lxml import etree
from tqdm import tqdm

# ---------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------
# Keep your existing DB_URL
DB_URL = "postgresql://neondb_owner:npg_AY5PZWVCQ8jJ@ep-green-snow-ad00zclq-pooler.c-2.us-east-1.aws.neon.tech/neondb?sslmode=require"

XML_FILE = "data/raw/drugbank/full database.xml"
DRUGS_CSV = "drugs_temp.csv"
INTER_CSV = "interactions_temp.csv"
NS = "{http://www.drugbank.ca}"

def parse_xml_to_csv(xml_path):
    print(f"🔨 Parsing XML and writing CSVs locally...")
    
    with open(DRUGS_CSV, 'w', newline='', encoding='utf-8') as f_drug, \
         open(INTER_CSV, 'w', newline='', encoding='utf-8') as f_inter:
        
        # Use QUOTE_MINIMAL to let Python handle the CSV structure
        drug_writer = csv.writer(f_drug, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
        inter_writer = csv.writer(f_inter, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
        
        context = etree.iterparse(xml_path, events=('end',), tag=f'{NS}drug')
        
        count = 0
        for event, elem in tqdm(context, desc="Parsing XML", unit=" drugs"):
            if elem.get('type') not in ['biotech', 'small molecule']:
                elem.clear()
                continue
            
            try:
                db_id = elem.findtext(f'{NS}drugbank-id[@primary="true"]')
                name = elem.findtext(f'{NS}name')
                
                # --- FIX: ROBUST ARRAY ESCAPING ---
                syns = []
                for s in elem.findall(f'{NS}synonyms/{NS}synonym'):
                    if s.text:
                        # 1. Escape backslashes first (Python needs \\ to mean literal \)
                        val = s.text.replace('\\', '\\\\')
                        # 2. Escape double quotes with backslash (Postgres Array Syntax: "Name \"A\"")
                        val = val.replace('"', '\\"')
                        syns.append(f'"{val}"')
                
                # Wrap in braces for Postgres Array: {"Syn1","Syn2"}
                synonyms_str = "{" + ",".join(syns) + "}"
                
                # SMILES
                smiles = elem.findtext(f'.//{NS}calculated-properties/{NS}property[{NS}kind="SMILES"]/{NS}value')
                
                # Write Drug Row
                drug_writer.writerow([db_id, name, synonyms_str, smiles, "{}", "{}", "{}"])
                
                # Interactions
                interactions = elem.findall(f'{NS}drug-interactions/{NS}drug-interaction')
                for inter in interactions:
                    target_id = inter.findtext(f'{NS}drugbank-id')
                    target_name = inter.findtext(f'{NS}name')
                    description = inter.findtext(f'{NS}description') or ""
                    
                    # Clean description for JSON
                    clean_desc = description.replace('\\', '\\\\').replace('"', '\\"')
                    evidence_json = f'{{"description": "{clean_desc}"}}'
                    
                    inter_writer.writerow([db_id, target_id, name, target_name, evidence_json, "Unknown"])
                
                count += 1
                elem.clear()
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
                    
            except Exception as e:
                # Print error but keep going to avoid crashing 99% of progress
                # print(f"Skipping row due to error: {e}")
                pass 

    print(f"✅ CSV Generation Complete. Processed {count} drugs.")

def upload_csv_to_db():
    print("🚀 Connecting to Cloud Database for Bulk Upload...")
    
    conn = None
    try:
        conn = psycopg2.connect(DB_URL)
        cur = conn.cursor()
        
        print("🧹 Wiping old data to prevent duplicates...")
        cur.execute("TRUNCATE TABLE interactions, drugs RESTART IDENTITY CASCADE;")
        conn.commit()
        
        print("⚡ Dropping indexes to speed up write...")
        cur.execute("DROP INDEX IF EXISTS idx_drugs_name;")
        cur.execute("DROP INDEX IF EXISTS idx_inter_a;")
        cur.execute("DROP INDEX IF EXISTS idx_inter_b;")
        conn.commit()
        
        print("📤 Uploading Drugs (Stream Copy)...")
        with open(DRUGS_CSV, 'r', encoding='utf-8') as f:
            cur.copy_expert("COPY drugs (drugbank_id, name, synonyms, smiles, atc_codes, properties, raw) FROM STDIN WITH CSV DELIMITER E'\t' QUOTE '\"' ESCAPE '\"'", f)
        conn.commit()
        
        print("📤 Uploading Interactions (Stream Copy)...")
        with open(INTER_CSV, 'r', encoding='utf-8') as f:
            cur.copy_expert("COPY interactions (drug_a, drug_b, drug_a_name, drug_b_name, evidence, severity) FROM STDIN WITH CSV DELIMITER E'\t' QUOTE '\"' ESCAPE '\"'", f)
        conn.commit()
        
        print("🔧 Re-creating indexes...")
        cur.execute("CREATE INDEX idx_drugs_name ON drugs(name);")
        cur.execute("CREATE INDEX idx_inter_a ON interactions(drug_a);")
        cur.execute("CREATE INDEX idx_inter_b ON interactions(drug_b);")
        conn.commit()
        
        print("✅ SUCCESS! Database is hydrated.")
        
    except Exception as e:
        print(f"❌ Upload Error: {e}")
    finally:
        if conn: conn.close()
        # Clean up temp files
        if os.path.exists(DRUGS_CSV): os.remove(DRUGS_CSV)
        if os.path.exists(INTER_CSV): os.remove(INTER_CSV)

if __name__ == "__main__":
    parse_xml_to_csv(XML_FILE)
    upload_csv_to_db()