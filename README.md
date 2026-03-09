
---

# **Hybrid AI-Powered Drug Interaction Detector**

### **Deterministic Database Lookup + Probabilistic ML Prediction + LLM Explanations**

## **📌 Overview**

Most drug interaction checkers are limited to static databases—if a pair isn't in the records, the system returns "No Interaction," which can be clinically dangerous.

This project implements a **Waterfall Architecture** that combines a verified database of **1.4 million interactions** with a **Random Forest Classifier** to predict risks for novel drugs, and a **BioT5 LLM** to generate human-readable clinical explanations.

---

## **🚀 Key Features**

* **Massive Dataset:** Indexes **73,000+ drugs** and **1.4M interactions** from DrugBank.
* **Hybrid Waterfall Logic:** 1. **Level 1 (Deterministic):** Ultra-fast PostgreSQL lookup for verified data.
2. **Level 2 (Predictive):** ML model predicts interaction probability for undocumented pairs.
3. **Level 3 (Generative):** BioT5 generates biological mechanistic explanations.
* **High Performance:** Sub-20ms query latency via optimized indexing and in-memory caching.
* **Production Ready:** Fully containerized with Docker and served via FastAPI.

---

## **🛠️ Technical Stack**

* **Backend:** FastAPI (Asynchronous Python)
* **Database:** PostgreSQL (Cloud Neon / Dockerized)
* **Data Engineering:** `lxml` for iterative XML parsing (handling 5GB+ datasets)
* **Machine Learning:** Scikit-learn (Random Forest), RDKit (Cheminformatics)
* **Generative AI:** BioT5-base via Hugging Face Inference API
* **Deployment:** Docker, Vercel/Render

---

## **🧠 The Machine Learning Pipeline**

For novel drug detection, we convert chemical **SMILES** strings into **Morgan Fingerprints** (2048-bit vectors).

1. **Featurization:** Drugs are converted into bit-vectors capturing molecular substructures.
2. **Inference:** A Random Forest model (92.5% accuracy) identifies potential clashing structures.
3. **Explanation:** If a risk is flagged, the system prompts **BioT5** to describe the mechanism (e.g., "CYP3A4 inhibition").

---

## **📂 Project Structure**

DDIDETECTOR/
├── api/                           # API related configurations/routes
├── backend/
│   └── app.py                     # Main FastAPI application entry point
├── data/                          # Datasets (ignored in version control usually)
│   ├── external-mappings/         # External ID mappings (e.g., RxNorm)
│   ├── processed/                 # Cleaned and processed data ready for training
│   └── raw/                       # Raw DrugBank XML files
├── db/
│   └── schema.sql                 # PostgreSQL database schema definitions
├── docs/                          # Project documentation and assets
│   └── Screenshot...png           
├── ingest/                        # ETL and Data Ingestion Pipeline
│   ├── build_name_index.py        # In-memory index builder for fast lookups
│   ├── drugbank_ingest.py         # Standard XML parsing script
│   ├── ingest_optimized.py        # Optimized iterparse script for memory efficiency
│   ├── rxnorm_lookup.py           # Drug mapping logic
│   └── update_smiles.py           # SMILES extraction and formatting
├── models/                        # Serialized ML Models
│   ├── biot5_finetuned/           # Fine-tuned BioT5 model weights/configs
│   ├── smiles_MLP_model.joblib    # Trained Multi-Layer Perceptron model
│   └── smiles_MLP_scaler.joblib   # Data scaler for Morgan Fingerprints
├── train/                         # Machine Learning Training Scripts
│   ├── 1_prepare_data.py          # Data splitting and balancing
│   ├── 2_train_model.py           # MLP Classifier training logic
│   ├── prepare_structure_dataset.py # SMILES to Morgan Fingerprint logic
│   ├── preprocess_t5_data.py      # LLM prompt/completion formatting
│   ├── run_evaluation.py          # Metrics calculation (Accuracy, ROC, ROUGE)
│   ├── test_model.py              # Inference testing scripts
│   └── test_smiles_model.py       # Specific SMILES prediction tests
├── ui/                            # Frontend interface components
├── .gitignore                     # Git ignore rules
├── check_index.py                 # Utility to verify DB indexing
├── debug_smiles.py                # Utility to troubleshoot structural hashes
├── export_for_colab.py            # Script to package data/models for GPU training
├── fast_ingest.py                 # High-speed data loading utility
├── requirements.txt               # Python dependencies
├── setup_cloud_db.py              # Cloud PostgreSQL (Neon) initialization
└── vercel.json                    # Deployment configuration for Vercel
```

---

## **⚙️ Installation & Setup**

### **1. Clone the Repository**

```bash
git clone https://github.com/your-username/drug-interaction-detector.git
cd drug-interaction-detector

```

### **2. Environment Variables**

Create a `.env` file:

```env
DATABASE_URL=your_postgres_url
HF_TOKEN=your_hugging_face_token

```

### **3. Run via Docker**

```bash
docker-compose up --build

```

*The API will be available at `http://localhost:8000*`

---

## **📈 Performance Metrics**

* **Database Query Latency:** <20ms
* **ML Prediction Accuracy:** 92.5%
* **Recall (Sensitivity):** 93.2% (Optimized for clinical safety)
* **Data Processing:** 5GB+ XML ingested with <500MB RAM usage.

---

## **🤝 Contributing**

Contributions are welcome! If you'd like to improve the predictive model or add more data sources, please open an issue or submit a pull request.

---
