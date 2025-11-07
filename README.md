

-----

# DDI Detector (Drug-Drug Interaction Detector)

DDI Detector is a proof-of-concept application designed to identify and report potential interactions between a list of drugs. It uses a comprehensive database built from publicly available biomedical data sources and exposes a simple REST API for checking interactions.

## Features

  * **Comprehensive Database:** Ingests and processes data from sources like DrugBank to build a robust PostgreSQL database of drugs and their known interactions.
  * **Data Ingestion Pipeline:** Includes optimized Python scripts for parsing large XML files and performing bulk data lookups.
  * **Name Normalization:** Uses the RxNorm API to map various drug names to standardized identifiers.
  * **REST API:** A simple FastAPI backend with a `/check` endpoint to query for interactions in real-time.
  * **Containerized Database:** Uses Docker to run a PostgreSQL database for easy setup and consistency.

## Technology Stack

  * **Backend:** Python, FastAPI
  * **Database:** PostgreSQL
  * **Data Libraries:** lxml, SQLAlchemy, pandas
  * **Environment:** Docker, Python `venv`

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






