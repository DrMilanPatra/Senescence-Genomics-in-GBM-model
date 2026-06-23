import scrublet as scr
import scipy.io
import pandas as pd
import os

input_dir = "/home/user/scrublet_input"
output_file = "/home/user/scrublet_combined_results.csv"
combined = []

for sample in os.listdir(input_dir):
    sample_dir = os.path.join(input_dir, sample)
    if not os.path.isdir(sample_dir):
        continue

    print(f"Processing {sample}...")

    # Load matrix and metadata
    matrix = scipy.io.mmread(os.path.join(sample_dir, "matrix.mtx")).T.tocsc()
    barcodes = pd.read_csv(os.path.join(sample_dir, "barcodes.tsv"), header=None)[0]

    # Run Scrublet
    scrub = scr.Scrublet(matrix, expected_doublet_rate=0.06)
    scores, preds = scrub.scrub_doublets()

    # Create DataFrame
    df = pd.DataFrame({
        'barcode': barcodes,
        'GSMID': sample,
        'doublet_score': scores,
        'predicted_doublet': preds
    })

    combined.append(df)

# Combine all results
result = pd.concat(combined)
result.to_csv(output_file, index=False)
print(f"\n Saved all results to {output_file}")
