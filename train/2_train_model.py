import os
from datasets import load_from_disk
from transformers import T5ForConditionalGeneration, T5Tokenizer, Trainer, TrainingArguments

# --- Configuration ---
MODEL_NAME = "QizhiPei/biot5-base" # Use the specialized biomedical model
TOKENIZED_DATA_DIR = "data/processed/tokenized_data"
MODEL_OUTPUT_DIR = "models/biot5_finetuned" # Output directory for the new model

def main():
    """
    Loads tokenized data and fine-tunes the BioT5 model.
    """
    # 1. Load the tokenizer and a pre-trained BioT5 model
    print(f"Loading tokenizer and model for '{MODEL_NAME}'...")
    tokenizer = T5Tokenizer.from_pretrained(MODEL_NAME)
    model = T5ForConditionalGeneration.from_pretrained(MODEL_NAME)

    # 2. Load the pre-processed datasets from disk
    print("Loading tokenized datasets from disk...")
    tokenized_train_dataset = load_from_disk(os.path.join(TOKENIZED_DATA_DIR, 'train'))
    tokenized_val_dataset = load_from_disk(os.path.join(TOKENIZED_DATA_DIR, 'validation'))

    # 3. Define the training arguments with speed optimizations
    training_args = TrainingArguments(
        output_dir=MODEL_OUTPUT_DIR,
        
        # --- SPEED OPTIMIZATIONS FOR POWERFUL GPUs ---
        fp16=True,                          # Enable mixed-precision training
        per_device_train_batch_size=16,     # Increase batch size
        per_device_eval_batch_size=16,
        gradient_accumulation_steps=4,      # Simulate an even larger effective batch size
        
        # --- Standard Arguments ---
        num_train_epochs=3,
        warmup_steps=200,
        weight_decay=0.01,
        logging_dir='./logs',
        logging_steps=100,
        eval_strategy="steps",
        eval_steps=200,                     # Evaluate/save periodically
        save_steps=200,
        load_best_model_at_end=True,
    )

    # 4. Create the Trainer instance
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=tokenized_train_dataset,
        eval_dataset=tokenized_val_dataset,
    )

    # 5. Start the fine-tuning process (will resume from checkpoint if available)
    print("\n--- Starting Fine-Tuning ---")
    # Use trainer.train() for the very first run
    # Use trainer.train(resume_from_checkpoint=True) for subsequent runs to resume
    trainer.train() # Or trainer.train(resume_from_checkpoint=True)
    print("--- Fine-Tuning Complete ---")

    # 6. Save the final model and tokenizer
    print(f"Saving the final model to '{MODEL_OUTPUT_DIR}'...")
    trainer.save_model(MODEL_OUTPUT_DIR)
    tokenizer.save_pretrained(MODEL_OUTPUT_DIR)
    print("Model saved successfully.")

if __name__ == "__main__":
    main()