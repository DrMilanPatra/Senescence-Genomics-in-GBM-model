from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import os

cpdb_file_path = "/home/user/CellphoneDB-master/v5.0.0/cellphonedb.zip"
input_path = "/home/user/CellPhoneDB_patient_outputs_final"
base_out_path = os.path.join(input_path, "results")
os.makedirs(base_out_path, exist_ok=True)

# ------------------------------
# AUTO-DETECT PATIENTS
# ------------------------------
patients = sorted(list(set([
    f.split("_")[1]
    for f in os.listdir(input_path)
    if f.startswith("meta_")
])))

sen_cols = [f"celltype_sen_{i}" for i in range(1, 8)]

print("Patients detected:", len(patients))

# ------------------------------
# RUN CELL PHONE DB
# ------------------------------
for p in patients:
    for sen_col in sen_cols:

        meta_file_path = os.path.join(input_path, f"meta_{p}_{sen_col}.txt")
        counts_file_path = os.path.join(input_path, f"counts_{p}_{sen_col}.txt")

        # ------------------------------
        # SAFETY CHECKS
        # ------------------------------
        if not os.path.exists(meta_file_path):
            print(f"SKIP {p} {sen_col}: missing meta")
            continue

        if not os.path.exists(counts_file_path):
            print(f"SKIP {p} {sen_col}: missing counts")
            continue

        # check file size (VERY IMPORTANT)
        if os.path.getsize(meta_file_path) < 100:
            print(f"SKIP {p} {sen_col}: empty meta")
            continue

        out_path = os.path.join(base_out_path, p, sen_col)
        os.makedirs(out_path, exist_ok=True)

        print(f"Running CPDB: {p} | {sen_col}")

        cpdb_statistical_analysis_method.call(
            cpdb_file_path=cpdb_file_path,
            meta_file_path=meta_file_path,
            counts_file_path=counts_file_path,
            counts_data='hgnc_symbol',
            score_interactions=True,
            iterations=1000,
            threshold=0.1,
            threads=28,
            debug_seed=42,
            pvalue=0.05,
            subsampling=False,
            output_path=out_path,
            output_suffix="statistical"
        )

        print(f"DONE: {p}, {sen_col} -> {out_path}")
