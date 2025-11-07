import json

NAME_INDEX_PATH = "data/processed/name_to_drugbank.json"
DRUGS_TO_CHECK = ["Warfarin", "Aspirin"]

print(f"Loading name index from: {NAME_INDEX_PATH}")
with open(NAME_INDEX_PATH, "r", encoding="utf-8") as f:
    name_map = json.load(f)

print("\n--- Checking for drugs ---")
for drug in DRUGS_TO_CHECK:
    if drug in name_map:
        print(f"✅ Found '{drug}': ID is {name_map[drug]}")
    else:
        print(f"❌ NOT FOUND: '{drug}' does not exist as an exact key in the index.")