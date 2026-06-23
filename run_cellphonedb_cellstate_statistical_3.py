from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import os

# ------------------------------
# Input files
# ------------------------------
cpdb_file_path = "/home/user/CellphoneDB-master/v5.0.0/cellphonedb.zip"

patients = [
    "rGBM-01", "ndGBM-01", "ndGBM-11", "ndGBM-02", "rGBM-02",
    "rGBM-03", "rGBM-04", "ndGBM-03", "rGBM-05", "ndGBM-10",
    "LGG-04", "ndGBM-04", "ndGBM-05", "ndGBM-06", "ndGBM-07", "ndGBM-08"
]

# Optional
active_tf_path = None
microenvs_file_path = None

# Paths
input_path = "/home/user/CellPhoneDB_cellstate_outputs/"
base_out_path = "/home/user/CellPhoneDB_cellstate_outputs/results"
os.makedirs(base_out_path, exist_ok=True)

# ------------------------------
# Loop over patients ONLY
# ------------------------------
for p in patients:
    
    print(f"Processing patient: {p}")
    
    meta_file_path = os.path.join(input_path, f"meta_{p}_cellstate.txt")
    counts_file_path = os.path.join(input_path, f"counts_{p}_cellstate.txt")

    # Skip if files missing
    if not os.path.exists(meta_file_path) or not os.path.exists(counts_file_path):
        print(f"Skipping {p}: input files not found")
        continue

    # Output folder
    out_path = os.path.join(base_out_path, p)
    os.makedirs(out_path, exist_ok=True)

    # Run CellPhoneDB
    cpdb_results = cpdb_statistical_analysis_method.call(
        cpdb_file_path = cpdb_file_path,
        meta_file_path = meta_file_path,
        counts_file_path = counts_file_path,
        counts_data = 'hgnc_symbol',
        active_tfs_file_path = active_tf_path,
        microenvs_file_path = microenvs_file_path,
        score_interactions = True,
        iterations = 1000,
        threshold = 0.1,
        threads = 28,
        debug_seed = 42,
        result_precision = 3,
        pvalue = 0.05,
        subsampling = False,
        separator = '|',
        output_path = out_path,
        output_suffix = "cellstate"
    )

    print(f"Finished patient: {p}, results in {out_path}")
