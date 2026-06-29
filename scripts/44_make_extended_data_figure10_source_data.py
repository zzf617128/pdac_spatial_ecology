from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
TABLE_DIR = PROJECT / "results" / "tables"
SOURCE_DIR = PROJECT / "results" / "source_data"
SOURCE_DIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    gradients = pd.read_csv(TABLE_DIR / "gse274557_full_caf_core_gradients.csv")
    context = pd.read_csv(TABLE_DIR / "gse274557_full_caf_core_context_summary.csv")

    sample_table = (
        gradients[
            [
                "dataset_id",
                "geo_accession",
                "sample_id",
                "title",
                "tissue",
                "treatment",
                "patient_id",
                "n_spots",
                "n_caf_core_spots",
            ]
        ]
        .drop_duplicates()
        .sort_values(["tissue", "sample_id"])
    )
    composition = (
        sample_table.groupby("tissue", as_index=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            total_spots=("n_spots", "sum"),
            median_spots_per_sample=("n_spots", "median"),
            total_caf_core_spots=("n_caf_core_spots", "sum"),
        )
        .sort_values("tissue")
    )
    composition.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_10A.csv", index=False)

    context = context.sort_values(["tissue", "target_program"])
    context.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_10B.csv", index=False)
    sample_table.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_10_Sample_Metadata.csv", index=False)

    print("Wrote Extended Data Figure 10 source data tables.")


if __name__ == "__main__":
    main()
