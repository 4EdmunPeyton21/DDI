import os
from datasets import load_dataset
from transformers import T5Tokenizer

# --- Configuration ---
MODEL_NAME = "QizhiPei/biot5-base" # Use the specialized biomedical model
DATA_DIR = "data/processed"
TOKENIZED_DATA_DIR = "data/processed/tokenized_data"

def main():
    """
    Loads raw datasets and saves tokenized versions to disk.
    """
    # 1. Load the tokenizer for the specified model
    print(f"Loading tokenizer for '{MODEL_NAME}'...")
    tokenizer = T5Tokenizer.from_pretrained(MODEL_NAME)

    # 2. Load your full datasets from the JSONL files
    print("Loading datasets from JSONL files...")
    train_dataset_full = load_dataset('json', data_files=os.path.join(DATA_DIR, 'train_dataset.jsonl'), split='train')
    val_dataset_full = load_dataset('json', data_files=os.path.join(DATA_DIR, 'val_dataset.jsonl'), split='train')

    # --- Create a medium-sized subset for efficient, high-quality training ---
    print("Creating a medium subset for efficient, high-quality training...")
    train_dataset = train_dataset_full.select(range(100000)) # Use 100,000 examples
    val_dataset = val_dataset_full.select(range(10000))   # Use 10,000 for validation
    # -----------------------------------------------------------------------

    # 3. Create a function to tokenize the data
    def tokenize_function(examples):
        inputs = ["describe interaction: " + doc for doc in examples['input']]
        model_inputs = tokenizer(inputs, max_length=128, truncation=True, padding="max_length")
        labels = tokenizer(text_target=examples['target'], max_length=256, truncation=True, padding="max_length")
        model_inputs["labels"] = labels["input_ids"]
        return model_inputs

    # 4. Apply the tokenization to the datasets
    print("Tokenizing datasets...")
    tokenized_train_dataset = train_dataset.map(tokenize_function, batched=True, remove_columns=['input', 'target'])
    tokenized_val_dataset = val_dataset.map(tokenize_function, batched=True, remove_columns=['input', 'target'])

    # 5. Save the tokenized datasets to disk
    print(f"Saving tokenized datasets to '{TOKENIZED_DATA_DIR}'...")
    os.makedirs(TOKENIZED_DATA_DIR, exist_ok=True)
    tokenized_train_dataset.save_to_disk(os.path.join(TOKENIZED_DATA_DIR, 'train'))
    tokenized_val_dataset.save_to_disk(os.path.join(TOKENIZED_DATA_DIR, 'validation'))

    print("Data preparation complete.")

if __name__ == "__main__":
    main()