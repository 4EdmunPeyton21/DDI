# Hybrid AI-Powered Drug Interaction Detector

**Deterministic Database Lookup · Predictive Machine Learning · LLM-Based Explanations**

---

## Overview

Traditional drug interaction systems rely heavily on static databases. If a drug pair is not recorded, they often return *“No Interaction”*, which can be misleading and potentially unsafe in clinical scenarios.

This project addresses that limitation by implementing a **hybrid, multi-stage architecture** that combines:

* A large-scale verified interaction database
* A machine learning model for predicting unknown interactions
* A language model for generating clear, human-readable explanations

The system is designed to provide both **high reliability** (through verified data) and **intelligent inference** (for unseen drug combinations).

---

## Key Features

* **Extensive Dataset**
  Covers over **73,000 drugs** and **1.4 million known interactions** sourced from DrugBank.

* **Three-Level Hybrid Architecture**

  1. **Deterministic Layer**
     Fast PostgreSQL lookup for confirmed interactions
  2. **Predictive Layer**
     Machine learning model estimates risk for unknown drug pairs
  3. **Generative Layer**
     BioT5 generates biological and clinical explanations

* **High Performance**
  Optimized queries with sub-20ms latency using indexing and caching.

* **Production-Ready Design**
  Containerized using Docker and exposed via FastAPI.

---

## Technical Stack

| Component        | Technology Used                                 |
| ---------------- | ----------------------------------------------- |
| Backend          | FastAPI (asynchronous Python)                   |
| Database         | PostgreSQL (Neon / Docker)                      |
| Data Processing  | lxml (efficient XML parsing for large datasets) |
| Machine Learning | Scikit-learn, RDKit                             |
| Generative AI    | BioT5 (Hugging Face Inference API)              |
| Deployment       | Docker, Vercel / Render                         |

---

## Machine Learning Pipeline

To detect interactions between previously unseen drugs, the system processes chemical structures using SMILES notation.

### Workflow

1. **Featurization**
   SMILES strings are converted into **Morgan fingerprints** (2048-bit vectors), capturing molecular substructures.

2. **Prediction**
   A trained **Random Forest classifier** (92.5% accuracy) identifies potential interaction risks.

3. **Explanation Generation**
   If an interaction is predicted, the system uses BioT5 to generate mechanistic explanations such as enzyme inhibition or metabolic interference.

---

## Project Structure

```
DDIDETECTOR/
├── api/                           # API routes and configurations
├── backend/
│   └── app.py                     # FastAPI entry point
├── data/
│   ├── external-mappings/         # External ID mappings (e.g., RxNorm)
│   ├── processed/                 # Cleaned datasets
│   └── raw/                       # Raw DrugBank XML files
├── db/
│   └── schema.sql                 # Database schema
├── docs/
│   └── screenshots/               # Documentation assets
├── ingest/                        # Data ingestion and ETL pipeline
│   ├── build_name_index.py
│   ├── drugbank_ingest.py
│   ├── ingest_optimized.py
│   ├── rxnorm_lookup.py
│   └── update_smiles.py
├── models/                        # Trained models and artifacts
│   ├── biot5_finetuned/
│   ├── smiles_MLP_model.joblib
│   └── smiles_MLP_scaler.joblib
├── train/                         # Model training scripts
│   ├── 1_prepare_data.py
│   ├── 2_train_model.py
│   ├── prepare_structure_dataset.py
│   ├── preprocess_t5_data.py
│   ├── run_evaluation.py
│   ├── test_model.py
│   └── test_smiles_model.py
├── ui/                            # Frontend components
├── requirements.txt
├── setup_cloud_db.py
├── fast_ingest.py
├── check_index.py
├── debug_smiles.py
├── export_for_colab.py
├── vercel.json
└── .gitignore
```

---

## Installation and Setup

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/drug-interaction-detector.git
cd drug-interaction-detector
```

### 2. Configure Environment Variables

Create a `.env` file in the root directory:

```env
DATABASE_URL=your_postgres_url
HF_TOKEN=your_hugging_face_token
```

### 3. Run with Docker

```bash
docker-compose up --build
```

The API will be available at:

```
http://localhost:8000
```

---

## Performance Metrics

* **Database Query Latency:** < 20 ms
* **Model Accuracy:** 92.5%
* **Recall (Sensitivity):** 93.2%
* **Data Processing Capability:**
  Handles 5GB+ XML datasets with under 500MB RAM usage

---

## Contributing

Contributions are welcome. You can:

* Improve model performance
* Add new datasets or sources
* Enhance system architecture or APIs

Please open an issue or submit a pull request with your proposed changes.

---
