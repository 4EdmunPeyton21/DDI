import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# --- Configuration ---
MODEL_PATH = "models/smiles_MLP_model.joblib"
SCALER_PATH = "models/smiles_MLP_scaler.joblib"
FINGERPRINT_SIZE = 2048

# --- Example SMILES strings ---
WARFARIN_SMILES = "CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3O)C(=O)OC2"
ASPIRIN_SMILES = "CC(=O)OC1=CC=CC=C1C(O)=O"
ACETAMINOPHEN_SMILES = "CC(=O)NC1=CC=C(O)C=C1"

def get_fingerprint(smiles_string):
    """Converts a SMILES string to a Morgan Fingerprint."""
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None: return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=FINGERPRINT_SIZE)
        return list(fp)
    except:
        return None

def main():
    print(f"Loading model from {MODEL_PATH}...")
    model = joblib.load(MODEL_PATH)

    print(f"Loading scaler from {SCALER_PATH}...")
    scaler = joblib.load(SCALER_PATH)

    print("Models loaded successfully.")

    # --- Test 1: Warfarin + Aspirin (Should be 1 - Interacts) ---
    fp_a = get_fingerprint(WARFARIN_SMILES)
    fp_b = get_fingerprint(ASPIRIN_SMILES)

    # Combine features and convert to 2D array
    features = np.array(fp_a + fp_b, dtype=np.float32).reshape(1, -1)
    # Scale the features
    features_scaled = scaler.transform(features)

    prediction = model.predict(features_scaled)[0]
    proba = model.predict_proba(features_scaled)[0][1] # Probability of class 1

    print(f"\nTest 1: Warfarin + Aspirin")
    print(f"  -> Prediction: {prediction} ('1' means interacts)")
    print(f"  -> Interaction Probability: {proba * 100:.2f}%")

    # --- Test 2: Lisinopril + Acetaminophen (Should be 0 - No Interaction) ---
    fp_a = get_fingerprint("Lisinopril_SMILES_placeholder") # You'd need the real SMILES
    fp_b = get_fingerprint(ACETAMINOPHEN_SMILES)

    # Note: This test will fail without the real Lisinopril SMILES.
    # It's here to show the process.

if __name__ == "__main__":
    main()