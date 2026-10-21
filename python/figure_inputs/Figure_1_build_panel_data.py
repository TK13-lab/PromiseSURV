#!/usr/bin/env python3
"""Build the single tidy input table consumed by Figure_1_cohort_construction.R."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def read_csv(input_root: Path, relative_path: str) -> pd.DataFrame:
    return pd.read_csv(input_root / relative_path)


def read_parquet(input_root: Path, relative_path: str) -> pd.DataFrame:
    return pd.read_parquet(input_root / relative_path)


def sex_from_one_hot(
    frame: pd.DataFrame,
    source_label: str,
    female_column: str,
    male_column: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, row in frame.iterrows():
        if row.get(female_column, 0) == 1:
            sex = "Female"
        elif row.get(male_column, 0) == 1:
            sex = "Male"
        else:
            sex = "Missing"
        rows.append(
            {
                "source": source_label,
                "cancer_label": row["cancer_label"],
                "patient_id": row["patient_id"],
                "sex": sex,
            }
        )
    return pd.DataFrame(rows)


def build_panel_data(input_root: Path, cancer_order: tuple[str, ...]) -> pd.DataFrame:
    task1 = read_csv(input_root, "outputs/task1/table1_discovery_validation_cohorts.csv")
    locked_missingness = read_csv(
        input_root,
        "outputs/task2/tables/supp_table_s4_missingness_rate_locked_mvp.csv",
    )
    cptac = read_parquet(
        input_root,
        "outputs/task3/model_inputs/cptac_mvp_complete_case_input.parquet",
    )
    tcga_clinical = read_csv(
        input_root,
        "outputs/task2/processed_modalities/clinical/clinical_features_tcga.csv",
    )
    tcga_predictions = read_csv(
        input_root,
        "outputs/task4/external_tcga/task4_tcga_common_clinical_predictions.csv",
    )

    task1_mvp = task1.loc[task1["cancer_label"].isin(cancer_order)].copy()
    cptac_summary = (
        cptac.groupby("cancer_label", as_index=False)
        .agg(n=("OS", "size"), events=("OS", "sum"), median_followup_days=("OS_time", "median"))
        .assign(source="CPTAC")
    )
    tcga_summary = (
        tcga_predictions.groupby("cancer_label", as_index=False)
        .agg(n=("OS", "size"), events=("OS", "sum"), median_followup_days=("OS_time", "median"))
        .assign(source="TCGA")
    )
    cohort_summary = pd.concat([cptac_summary, tcga_summary], ignore_index=True)
    cohort_summary["non_events"] = cohort_summary["n"] - cohort_summary["events"]
    cohort_summary["event_rate"] = cohort_summary["events"] / cohort_summary["n"]
    cohort_summary["cancer_label"] = pd.Categorical(
        cohort_summary["cancer_label"], categories=cancer_order, ordered=True
    )
    cohort_summary = cohort_summary.sort_values(["source", "cancer_label"]).reset_index(drop=True)

    locked_clinical = locked_missingness.query("modality == 'clinical'").copy()
    count_rows: list[dict[str, object]] = []
    for cancer in cancer_order:
        task1_row = task1_mvp.loc[task1_mvp["cancer_label"] == cancer].iloc[0]
        cptac_clinical_row = locked_clinical.loc[
            locked_clinical["cancer_label"] == cancer
        ].iloc[0]
        count_rows.extend(
            [
                {
                    "cancer_label": cancer,
                    "data_layer": "TCGA CDR",
                    "patient_count": int(task1_row["tcga_cdr_patients"]),
                },
                {
                    "cancer_label": cancer,
                    "data_layer": "TCGA deploy",
                    "patient_count": int(
                        tcga_summary.loc[tcga_summary["cancer_label"] == cancer, "n"].iloc[0]
                    ),
                },
                {
                    "cancer_label": cancer,
                    "data_layer": "CPTAC clinical",
                    "patient_count": int(cptac_clinical_row["available_count"]),
                },
                {
                    "cancer_label": cancer,
                    "data_layer": "CPTAC complete",
                    "patient_count": int(
                        cptac_summary.loc[cptac_summary["cancer_label"] == cancer, "n"].iloc[0]
                    ),
                },
            ]
        )
    count_matrix = pd.DataFrame(count_rows)

    tcga_clinical = tcga_clinical.loc[
        tcga_clinical["patient_id"].isin(set(tcga_predictions["patient_id"]))
        & tcga_clinical["cancer_label"].isin(cancer_order)
    ].copy()
    sex = pd.concat(
        [
            sex_from_one_hot(cptac, "CPTAC", "gender_clean_FEMALE", "gender_clean_MALE"),
            sex_from_one_hot(tcga_clinical, "TCGA", "gender_FEMALE", "gender_MALE"),
        ],
        ignore_index=True,
    )
    sex["sex"] = pd.Categorical(
        sex["sex"], categories=("Female", "Male", "Missing"), ordered=True
    )
    sex_composition = (
        sex.groupby(["source", "cancer_label", "sex"], observed=False)
        .size()
        .reset_index(name="n")
    )
    sex_composition["total"] = sex_composition.groupby(
        ["source", "cancer_label"], observed=False
    )["n"].transform("sum")
    sex_composition["percent"] = np.where(
        sex_composition["total"] > 0,
        100 * sex_composition["n"] / sex_composition["total"],
        0,
    )

    cptac_age_columns = ["cancer_label", "patient_id", "age_at_diagnosis_z"]
    cptac_age = cptac.loc[
        cptac["age_at_diagnosis_z_missing"].fillna(0) == 0,
        cptac_age_columns,
    ].assign(source="CPTAC")
    tcga_age = tcga_clinical.loc[
        tcga_clinical["age_at_diagnosis_z_missing"].fillna(0) == 0,
        ["source", "cancer_label", "patient_id", "age_at_diagnosis_z"],
    ]
    age_distribution = pd.concat([cptac_age, tcga_age], ignore_index=True)

    followup_distribution = pd.concat(
        [
            cptac[["cancer_label", "patient_id", "OS_time"]].assign(source="CPTAC"),
            tcga_predictions[["source", "cancer_label", "patient_id", "OS_time"]],
        ],
        ignore_index=True,
    )
    followup_distribution["followup_years"] = followup_distribution["OS_time"] / 365.25

    tcga_modality_rows: list[dict[str, object]] = []
    for cancer in cancer_order:
        patient_count = int((tcga_predictions["cancer_label"] == cancer).sum())
        for modality, available_count in (
            ("clinical", patient_count),
            ("rna", patient_count),
            ("proteomics", 0),
            ("phosphoproteomics", 0),
        ):
            tcga_modality_rows.append(
                {
                    "source": "TCGA_DEPLOYMENT",
                    "cancer_label": cancer,
                    "cohort_id": f"TCGA_{cancer}",
                    "modality": modality,
                    "patient_count": patient_count,
                    "available_count": available_count,
                    "missing_count": patient_count - available_count,
                    "missing_rate": (patient_count - available_count) / patient_count,
                }
            )
    missingness = pd.concat(
        [locked_missingness, pd.DataFrame(tcga_modality_rows)], ignore_index=True
    )
    missingness = missingness.loc[
        missingness["modality"].isin(
            ("clinical", "rna", "proteomics", "phosphoproteomics")
        )
    ].copy()
    missingness["cohort_label"] = np.where(
        missingness["source"].str.contains("TCGA", case=False, na=False),
        "TCGA " + missingness["cancer_label"].astype(str),
        "CPTAC " + missingness["cancer_label"].astype(str),
    )

    panel_frames = []
    for panel, frame in (
        ("A", count_matrix),
        ("B", cohort_summary),
        ("C", followup_distribution),
        ("D", age_distribution),
        ("E", sex_composition),
        ("F", missingness),
    ):
        panel_frame = frame.copy()
        panel_frame.insert(0, "figure_panel", panel)
        panel_frames.append(panel_frame)
    return pd.concat(panel_frames, ignore_index=True, sort=False)


def main() -> None:
    package_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, default=package_root)
    parser.add_argument(
        "--output",
        type=Path,
        default=package_root
        / "outputs"
        / "figures"
        / "figure2"
        / "figure2_cohort_construction_panel_data.csv",
    )
    args = parser.parse_args()

    # Locked manuscript values are intentionally local and visible here.
    cancer_order = ("LUAD", "KIRC", "HNSC")
    panel_data = build_panel_data(args.input_root.resolve(), cancer_order)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    panel_data.to_csv(args.output, index=False)
    print(f"Wrote {args.output.resolve()}")


if __name__ == "__main__":
    main()
