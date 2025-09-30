
DDI Detector: A Drug-Drug Interaction Checker
DDI Detector is a backend system and data pipeline for identifying potential drug-drug interactions (DDIs). It processes raw data from biomedical sources like DrugBank, stores it in a structured PostgreSQL database, and exposes a real-time REST API to check for interactions between a given list of medications.

Overview
The project's core is a robust ETL (Extract, Transform, Load) pipeline that efficiently handles multi-gigabyte XML data files. The processed data is served via a high-performance FastAPI application, designed for quick and responsive interaction checks. This serves as a proof-of-concept for a clinical decision support tool.

Features
High-Performance Data Ingestion: Optimized Python scripts using iterative parsing and batch inserts to process large datasets in hours, not days.

Containerized Database: A PostgreSQL database running in a Docker container for a consistent and portable development environment.

Data Standardization: Enriches the dataset with universal drug identifiers from the public RxNorm API.

Real-Time API: A fast and responsive REST API built with FastAPI, capable of handling on-demand interaction checks.

Optimized Lookups: Utilizes a pre-computed, in-memory name index that includes millions of drug names and synonyms for near-instant, case-insensitive lookups.

Technology Stack
Backend: Python, FastAPI

Database: PostgreSQL, SQLAlchemy

Data Processing: lxml, pandas

Environment: Docker, Python venv

API Client: requests

Project Structure
ddidetector/
├── .venv/                  # Python virtual environment (ignored by Git)
├── backend/
│   └── app.py              # FastAPI application and the /check endpoint
├── data/
│   ├── raw/                # Raw, untouched data sources (ignored by Git)
│   │   └── drugbank/
│   │       └── full database.xml
│   └── processed/
│       └── name_to_drugbank.json # The pre-computed name index for the API
├── db/
│   └── schema.sql          # SQL script to create the database tables
├── ingest/
│   ├── ingest_optimized.py     # Script to ingest drugs into the 'drugs' table
│   ├── interactions_ingest.py  # Script to ingest interactions into the 'interactions' table
│   ├── rxnorm_lookup.py        # Script to enrich data with RxNorm IDs
│   └── build_name_index.py     # Script to create the name_to_drugbank.json index
├── .gitignore              # Specifies files and folders for Git to ignore
└── README.md               # Project documentation (this file)
Setup and Installation
Prerequisites
Python 3.10+

Git

Docker Desktop

1. Clone & Set Up Environment
PowerShell

# Clone the repository
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






