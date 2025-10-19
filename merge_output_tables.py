import pandas as pd
import numpy as np
import os

# Common settings
OUTPUT_DIR = "OUTPUT"

# Table families and suffixes
table_families = ["table_lotss_closed_dipole", "table_lotss_pfss","table_vlass_closed_dipole", "table_vlass_pfss"]
suffixes = ["0.1", "1.0", "10.0"]

# Columns that differ between tables
diff_cols = [
    "Flux_r_S(muJy)",
    "Flux_reconnect(muJy)",
    "Flux_sb(muJy)",
    "Flux_r_S_NOABS",
    "Flux_reconnect_NOABS",
    "Flux_sb_NOABS",
    "M_A_nominal"
]

for family in table_families:
    print(f"🔄 Processing {family}...")

    # Read all CSVs for this family
    dfs = []
    for suffix in suffixes:
        file_path = os.path.join(OUTPUT_DIR, f"{family}_M_star_dot_{suffix}.csv")
        dfs.append(pd.read_csv(file_path))

    # Use the first as base
    df_base = dfs[0]
    common_cols = [c for c in df_base.columns if c not in diff_cols]
    df_common = df_base[common_cols].copy()

    # Add differing columns for each suffix
    for col in diff_cols:
        for suffix, df in zip(suffixes, dfs):
            df_common[f"{col}_{suffix}"] = df[col]
    '''
    # Format numeric values in scientific notation with two decimals
    numeric_cols = df_common.select_dtypes(include=[np.number]).columns
    df_common[numeric_cols] = df_common[numeric_cols].applymap(
        lambda x: f"{x:.2e}" if pd.notnull(x) else ""
    )
    '''
    # Save the merged result
    output_file = os.path.join(OUTPUT_DIR, f"{family}_merged.csv")
    df_common.to_csv(output_file, index=False, na_rep='')

    print(f"✅ Saved merged table: {output_file}")



print("🎉 All tables processed successfully!")

