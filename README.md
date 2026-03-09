

-----

Hybrid AI-Powered Drug Interaction Detector
Deterministic Database Lookup + Probabilistic ML Prediction + LLM Explanations
📌 Overview
Most drug interaction checkers are limited to static databases—if a pair isn't in the records, the system returns "No Interaction," which can be clinically dangerous.

This project implements a Waterfall Architecture that combines a verified database of 1.4 million interactions with a Random Forest Classifier to predict risks for novel drugs, and a BioT5 LLM to generate human-readable clinical explanations.

🚀 Key Features
Massive Dataset: Indexes 73,000+ drugs and 1.4M interactions from DrugBank.

Hybrid Waterfall Logic: 1. Level 1 (Deterministic): Ultra-fast PostgreSQL lookup for verified data.
2. Level 2 (Predictive): ML model predicts interaction probability for undocumented pairs.
3. Level 3 (Generative): BioT5 generates biological mechanistic explanations.

High Performance: Sub-20ms query latency via optimized indexing and in-memory caching.

Production Ready: Fully containerized with Docker and served via FastAPI.

🛠️ Technical Stack
Backend: FastAPI (Asynchronous Python)

Database: PostgreSQL (Cloud Neon / Dockerized)

Data Engineering: lxml for iterative XML parsing (handling 5GB+ datasets)

Machine Learning: Scikit-learn (Random Forest), RDKit (Cheminformatics)

Generative AI: BioT5-base via Hugging Face Inference API

Deployment: Docker, Vercel/Render

🧠 The Machine Learning Pipeline
For novel drug detection, we convert chemical SMILES strings into Morgan Fingerprints (2048-bit vectors).

Featurization: Drugs are converted into bit-vectors capturing molecular substructures.

Inference: A Random Forest model (92.5% accuracy) identifies potential clashing structures.

Explanation: If a risk is flagged, the system prompts BioT5 to describe the mechanism (e.g., "CYP3A4 inhibition").

📂 Project Structure
Bash
├── data_pipeline/         # ETL scripts for DrugBank XML parsing
├── ml_models/             # Random Forest training & Morgan Fingerprint logic
├── api/                   # FastAPI endpoints and business logic
├── db/                    # PostgreSQL schema and migration scripts
├── docker-compose.yml     # Container orchestration
└── requirements.txt       # Project dependencies
⚙️ Installation & Setup
1. Clone the Repository
Bash
git clone https://github.com/your-username/drug-interaction-detector.git
cd drug-interaction-detector
2. Environment Variables
Create a .env file:

Code snippet
DATABASE_URL=your_postgres_url
HF_TOKEN=your_hugging_face_token
3. Run via Docker
Bash
docker-compose up --build
The API will be available at http://localhost:8000

📈 Performance Metrics
Database Query Latency: <20ms

ML Prediction Accuracy: 92.5%

Recall (Sensitivity): 93.2% (Optimized for clinical safety)

Data Processing: 5GB+ XML ingested with <500MB RAM usage.

🤝 Contributing
Contributions are welcome! If you'd like to improve the predictive model or add more data sources, please open an issue or submit a pull request.

-----

## Setup and Installation

### Prerequisites

  * Python 3.10+
  * Git
  * Docker Desktop

### 1\. Clone the Repository

```bash
git clone https://github.com/4EdmunPeyton21/DDI.git
cd DDI

# Create and activate a Python virtual environment
python -m venv .venv
.venv\Scripts\Activate.ps1

# Install required dependencies
pip install fastapi uvicorn sqlalchemy psycopg2-binary pydantic lxml requests
2. Prepare Data
Download the DrugBank XML file (full database.xml) and place it inside the data/raw/drugbank/ directory.

3. Start & Set Up Database
Ensure Docker Desktop is running.

PowerShell

# Run the PostgreSQL container
docker run --name ddidb -e POSTGRES_PASSWORD=devpass -e POSTGRES_USER=dduser -e POSTGRES_DB=ddidb -p 5432:5432 -d postgres:15

# Create the database tables from the schema file
type db\schema.sql | docker exec -i ddidb psql -U dduser -d ddidb
4. Run Data Ingestion Pipeline
Run these scripts in order to populate the database. This is a one-time setup process.

PowerShell

# 1. Ingest drugs
python ingest/ingest_optimized.py "data/raw/drugbank/full database.xml"

# 2. Ingest interactions
python ingest/interactions_ingest.py "data/raw/drugbank/full database.xml"

# 3. Get RxNorm IDs
python ingest/rxnorm_lookup.py

# 4. Build the name index for the API
python ingest/build_name_index.py
Usage
1. Run the API Server
Start the FastAPI server from the project's root directory.

PowerShell

uvicorn backend.app:app --reload --port 8000
The server will be available at http://127.0.0.1:8000.

2. Test the Endpoint
In a new terminal, send a POST request to the /check endpoint.

Example using PowerShell:

PowerShell

Invoke-RestMethod -Uri "http://127.0.0.1:8000/check" -Method Post -ContentType "application/json" -Body '{"drugs":["Warfarin", "Aspirin"]}'
Expected Response:

JSON

{
  "interactions": [
    {
      "pair": [
        "Warfarin",
        "Aspirin"
      ],
      "description": "Acetylsalicylic acid may increase the anticoagulant activities of Warfarin."
    }
  ]
}






