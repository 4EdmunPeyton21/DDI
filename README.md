Of course. Here is a complete `README.md` file content for your project. You can copy and paste this directly into your `README.md` file.

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
```

### 2\. Set Up the Environment

Create and activate a Python virtual environment.

```powershell
# For PowerShell
python -m venv .venv
.venv\Scripts\Activate.ps1
```

Install the required dependencies.

```powershell
pip install -r requirements.txt
```

*(Note: You will need to create a `requirements.txt` file by running `pip freeze > requirements.txt`)*

### 3\. Prepare Data

Download the necessary raw data files (e.g., the DrugBank XML) and place them in the appropriate subdirectories within `data/raw/`.

### 4\. Start the Database

Ensure Docker Desktop is running, then start the PostgreSQL container.

```bash
docker run --name ddidb -e POSTGRES_PASSWORD=devpass -e POSTGRES_USER=dduser -e POSTGRES_DB=ddidb -p 5432:5432 -d postgres:15
```

Create the database schema.

```powershell
# For PowerShell
type db\schema.sql | docker exec -i ddidb psql -U dduser -d ddidb
```

### 5\. Run Ingestion Scripts

Run the scripts in the following order to populate the database.

```powershell
# 1. Ingest drugs from DrugBank
python ingest/ingest_optimized.py "data/raw/drugbank/your_drugbank_file.xml"

# 2. Ingest interactions from DrugBank
python ingest/interactions_ingest.py "data/raw/drugbank/your_drugbank_file.xml"

# 3. Normalize names with RxNorm API
python ingest/rxnorm_lookup.py

# 4. Build the name index for the API
python ingest/build_name_index.py
```

-----

## Usage

### 1\. Run the API Server

Start the FastAPI server from the project's root directory.

```powershell
uvicorn backend.app:app --reload --port 8000
```

The server will be available at `http://127.0.0.1:8000`.

### 2\. Test the Endpoint

In a new terminal, you can send a POST request to the `/check` endpoint to find interactions.

**Example using PowerShell:**

```powershell
Invoke-RestMethod -Uri "http://127.0.0.1:8000/check" -Method Post -ContentType "application/json" -Body '{"drugs":["Warfarin", "Aspirin"]}'
```

**Expected Response:**

```json
{
  "interactions": [
    {
      "pair": [
        "Warfarin",
        "Acetylsalicylic acid"
      ],
      "description": "Acetylsalicylic acid may increase the anticoagulant activities of Warfarin."
    }
  ]
}
```