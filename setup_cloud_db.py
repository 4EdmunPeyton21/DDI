import sqlalchemy
from sqlalchemy import create_engine, text

# ------------------------------------------------------------------
# PASTE YOUR NEON / VERCEL DATABASE URL HERE
# Example: "postgresql://neondb_owner:..."
# ------------------------------------------------------------------
DB_URL = "postgresql://neondb_owner:npg_AY5PZWVCQ8jJ@ep-green-snow-ad00zclq-pooler.c-2.us-east-1.aws.neon.tech/neondb?sslmode=require" 

# Define the Tables (The Schema)
SCHEMA_SQL = """
-- 1. Create Drugs Table
CREATE TABLE IF NOT EXISTS drugs (
    drugbank_id VARCHAR(50) PRIMARY KEY,
    name VARCHAR(255),
    synonyms TEXT[],
    smiles TEXT,
    atc_codes TEXT[],
    properties JSONB DEFAULT '{}',
    raw JSONB DEFAULT '{}'
);

-- 2. Create Interactions Table
CREATE TABLE IF NOT EXISTS interactions (
    id SERIAL PRIMARY KEY,
    drug_a VARCHAR(50) REFERENCES drugs(drugbank_id),
    drug_b VARCHAR(50) REFERENCES drugs(drugbank_id),
    drug_a_name VARCHAR(255),
    drug_b_name VARCHAR(255),
    evidence JSONB,
    severity VARCHAR(50)
);

-- 3. Create Indexes for Speed
CREATE INDEX IF NOT EXISTS idx_drugs_name ON drugs(name);
CREATE INDEX IF NOT EXISTS idx_inter_a ON interactions(drug_a);
CREATE INDEX IF NOT EXISTS idx_inter_b ON interactions(drug_b);
"""

def init_db():
    print(f"🔌 Connecting to: {DB_URL.split('@')[1] if '@' in DB_URL else 'Database'}...")
    
    # Connect with SSL (Required for Neon/Vercel)
    try:
        # Check if we need to fix the postgres:// prefix for SQLAlchemy
        clean_url = DB_URL.replace("postgres://", "postgresql://")
        
        engine = create_engine(clean_url, connect_args={"sslmode": "require"})
        
        with engine.connect() as conn:
            # Run the SQL commands
            conn.execute(text(SCHEMA_SQL))
            conn.commit()
            
        print("✅ SUCCESS! Tables 'drugs' and 'interactions' created.")
        print("🚀 You can now run the ingestion script.")
        
    except Exception as e:
        print(f"❌ ERROR: {e}")

if __name__ == "__main__":
    init_db()