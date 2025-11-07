from transformers import T5ForConditionalGeneration, T5Tokenizer

# --- Configuration ---
# This path points to where you saved your downloaded model files
MODEL_PATH = "models/biot5_finetuned"

# --- The drug pair you want to test ---
drug_a = "Glimepiride"
drug_b = "acetylsalicylic acid"
# The input format must match what the model was trained on
input_text = f"describe interaction: {drug_a} + {drug_b}"

def main():
    """
    Loads the fine-tuned model from a local directory and generates
    an interaction description.
    """
    print(f"Loading fine-tuned model and tokenizer from '{MODEL_PATH}'...")
    try:
        # Load the tokenizer and model from your local directory
        tokenizer = T5Tokenizer.from_pretrained(MODEL_PATH)
        model = T5ForConditionalGeneration.from_pretrained(MODEL_PATH)
        print("Model loaded successfully.")
    except Exception as e:
        print(f"Error loading model from {MODEL_PATH}: {e}")
        print("Please ensure the model files are correctly placed in the directory.")
        return

    print(f"\nInput prompt: '{input_text}'")

    # 1. Tokenize the input text
    # Convert the text into numerical IDs the model understands
    inputs = tokenizer(input_text, return_tensors="pt", max_length=128, truncation=True)

    # 2. Generate the output sequence of IDs
    # `max_length` controls how long the generated description can be
    print("Generating interaction description...")
    outputs = model.generate(**inputs, max_length=256, num_beams=4, early_stopping=True)

    # 3. Decode the output IDs back into text
    # Convert the numerical IDs back into a human-readable sentence
    generated_text = tokenizer.decode(outputs[0], skip_special_tokens=True)

    print("\n--- Generated Interaction ---")
    print(generated_text)
    print("---------------------------\n")

if __name__ == "__main__":
    main()