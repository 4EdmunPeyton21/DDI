import requests
import time
import concurrent.futures
from sqlalchemy import create_engine, text

# --- Configuration ---
DB_URL = "postgresql://dduser:devpass@localhost:5432/ddidb"
RXNAV_API_URL = "https://rxnav.nlm.nih.gov/REST/rxcui.json?name={}"
MAX_WORKERS = 10  # Number of parallel requests to make. Adjust based on your connection.

def get_rxcui(name):
    """Fetches the RxCUI for a given drug name from the RxNav API."""
    try:
        r = requests.get(RXNAV_API_URL.format(name), timeout=10)
        r.raise_for_status()
        data = r.json()
        return data.get('idGroup', {}).get('rxnormId', [None])[0]
    except requests.exceptions.RequestException:
        # We can ignore failed requests for this bulk process
        return None

def fetch_worker(drug):
    """A single worker's task: take a drug, fetch its RxCUI."""
    drugbank_id, name = drug
    if not name:
        return None
    rxcui = get_rxcui(name)
    if rxcui:
        return {'dbid': drugbank_id, 'rxcui': rxcui}
    return None

def update_rxnorm_id_batch(engine, results):
    """Updates the database with a batch of results."""
    if not results:
        return
    with engine.begin() as conn:
        stmt = text("UPDATE drugs SET rxnorm_id = :rxcui WHERE drugbank_id = :dbid")
        conn.execute(stmt, results)

def main():
    """Main function to run the optimized RxNorm lookup process."""
    engine = create_engine(DB_URL)
    
    with engine.connect() as conn:
        result = conn.execute(text("SELECT drugbank_id, name FROM drugs WHERE rxnorm_id IS NULL"))
        drugs_to_process = result.fetchall()

    print(f"Found {len(drugs_to_process)} drugs to process. Starting parallel API requests...")
    
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # Use executor.map to run fetch_worker on all drugs concurrently
        future_to_drug = {executor.submit(fetch_worker, drug): drug for drug in drugs_to_process}
        
        for i, future in enumerate(concurrent.futures.as_completed(future_to_drug)):
            res = future.result()
            if res:
                results.append(res)
            
            if (i + 1) % 100 == 0:
                print(f"  ...requests completed: {i + 1} / {len(drugs_to_process)}")

    print(f"\nAPI requests complete. Found {len(results)} RxNorm IDs.")
    print("Updating database in a single batch...")
    
    update_rxnorm_id_batch(engine, results)
    
    print("RxNorm lookup complete.")

if __name__ == "__main__":
    main()