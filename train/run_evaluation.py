import os
import torch
from datasets import load_dataset
from transformers import T5ForConditionalGeneration, T5Tokenizer
from evaluate import load # Use the new evaluate library
import nltk # NLTK is needed for BLEU calculation

# --- Configuration ---
MODEL_PATH = "models/biot5_finetuned" # Path to your local fine-tuned model
DATA_DIR = "data/processed"
TEST_DATA_FILE = os.path.join(DATA_DIR, 'test_dataset.jsonl')
BATCH_SIZE = 16 # Adjust based on your GPU memory
SUBSET_SIZE = 1000 # Set to None to use the full test set, or a number for a subset

def main():
    """
    Loads the fine-tuned model, generates predictions on the test set (or a subset),
    and calculates ROUGE scores.
    """
    # --- 1. Load Model and Tokenizer ---
    print(f"Loading model and tokenizer from '{MODEL_PATH}'...")
    try:
        tokenizer = T5Tokenizer.from_pretrained(MODEL_PATH)
        model = T5ForConditionalGeneration.from_pretrained(MODEL_PATH)
    except Exception as e:
        print(f"Error loading model: {e}. Ensure model files are in {MODEL_PATH}")
        return

    # Determine device (CPU or GPU)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model.to(device)
    print(f"Model loaded successfully onto {device}.")

    # --- 2. Load Test Dataset ---
    print(f"Loading test dataset from '{TEST_DATA_FILE}'...")
    try:
        test_dataset_full = load_dataset('json', data_files=TEST_DATA_FILE, split='train')
        print(f"Loaded {len(test_dataset_full)} total examples for evaluation.")

        # --- SELECT SUBSET (IF CONFIGURED) ---
        if SUBSET_SIZE is not None:
            if SUBSET_SIZE < len(test_dataset_full):
                test_dataset = test_dataset_full.select(range(SUBSET_SIZE))
                print(f"Using a subset of {len(test_dataset)} examples for quick evaluation.")
            else:
                test_dataset = test_dataset_full # Subset size is larger than dataset
                print(f"Subset size ({SUBSET_SIZE}) >= dataset size. Using full test set.")
        else:
            test_dataset = test_dataset_full # Use the full test set
            print("Using the full test set for evaluation.")
        # ------------------------------------

    except Exception as e:
        print(f"Error loading test dataset: {e}")
        return

    # --- 3. Prepare Data Processing Function ---
    def prepare_data(examples):
        inputs = ["describe interaction: " + doc for doc in examples['input']]
        model_inputs = tokenizer(inputs, max_length=128, truncation=True, padding="longest")
        return {
            "input_ids": model_inputs.input_ids,
            "attention_mask": model_inputs.attention_mask,
            "labels": examples["target"] # Keep original target text for comparison
        }

    print("Preprocessing test dataset...")
    processed_test_dataset = test_dataset.map(prepare_data, batched=True, batch_size=BATCH_SIZE)

    # --- 4. Generate Predictions ---
    print("Generating predictions on the test set (this may take a while)...")
    all_predictions = []
    all_references = []

    # Calculate total batches needed based on the potentially smaller test_dataset
    total_batches = (len(processed_test_dataset) + BATCH_SIZE - 1) // BATCH_SIZE

    for i in range(0, len(processed_test_dataset), BATCH_SIZE):
        batch = processed_test_dataset[i : i + BATCH_SIZE]
        input_ids = torch.tensor(batch['input_ids']).to(device)
        attention_mask = torch.tensor(batch['attention_mask']).to(device)

        try:
            with torch.no_grad():
                outputs = model.generate(
                    input_ids=input_ids,
                    attention_mask=attention_mask,
                    max_length=256,
                    num_beams=4,
                    early_stopping=True
                )
            
            predictions = tokenizer.batch_decode(outputs, skip_special_tokens=True)
            all_predictions.extend(predictions)
            all_references.extend(batch['labels'])

            current_batch_num = i // BATCH_SIZE + 1
            if current_batch_num % 10 == 0 or current_batch_num == total_batches: # Print every 10 batches or on the last batch
                 print(f"  ...processed batch {current_batch_num} / {total_batches}")

        except Exception as e:
             print(f"Error during generation for batch starting at index {i}: {e}")
             all_predictions.extend(["<error>"] * len(batch['input_ids']))
             all_references.extend(batch['labels'])


    print("Prediction generation complete.")

    # --- 5. Calculate Metrics ---
    print("Calculating ROUGE scores...")
    try:
        rouge_metric = load("rouge")
        rouge_results = rouge_metric.compute(predictions=all_predictions, references=all_references)

        print("\n--- Evaluation Results ---")
        print("ROUGE Scores:")
        for key, value in rouge_results.items():
            print(f"  {key}: {value * 100:.2f}") # Print scores as percentages
    except Exception as e:
        print(f"Could not calculate ROUGE scores: {e}")

    # --- Optional: Calculate BLEU Score ---
    # try:
    #     print("\nCalculating BLEU score...")
    #     nltk.download('punkt', quiet=True) # Required for BLEU
    #     bleu_metric = load("bleu")
    #     # BLEU expects references to be a list of lists
    #     bleu_results = bleu_metric.compute(predictions=all_predictions, references=[[ref] for ref in all_references])
    #     print("\nBLEU Score:")
    #     print(f"  BLEU: {bleu_results['bleu'] * 100:.2f}")
    # except Exception as e:
    #     print(f"\nCould not calculate BLEU score: {e}")

    print("\nEvaluation complete.")

if __name__ == "__main__":
    main()