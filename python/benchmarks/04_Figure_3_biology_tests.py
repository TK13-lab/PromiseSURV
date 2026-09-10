#!/usr/bin/env python3
"""Recompute two-sided correlation tests underlying Figure 3."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    values = p_values.to_numpy(dtype=float)
    adjusted = np.full(len(values), np.nan)
    finite_indices = np.flatnonzero(np.isfinite(values))
    if len(finite_indices) == 0:
        return pd.Series(adjusted, index=p_values.index)
    finite_values = values[finite_indices]
    order = np.argsort(finite_values)
    ranked = finite_values[order]
    raw_adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    monotone_adjusted = np.minimum.accumulate(raw_adjusted[::-1])[::-1]
    restored = np.empty(len(ranked))
    restored[order] = np.minimum(monotone_adjusted, 1.0)
    adjusted[finite_indices] = restored
    return pd.Series(adjusted, index=p_values.index)


def safe_spearman(x: pd.Series, y: pd.Series) -> tuple[int, float, float]:
    paired = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(paired) < 3 or paired["x"].nunique() < 2 or paired["y"].nunique() < 2:
        return len(paired), float("nan"), float("nan")
    result = spearmanr(paired["x"], paired["y"], alternative="two-sided")
    return len(paired), float(result.statistic), float(result.pvalue)


def teacher_student_tests(
    scores: pd.DataFrame,
    cancer_order: tuple[str, ...],
) -> pd.DataFrame:
    rows = []
    for cancer in cancer_order:
        cohort = scores.loc[scores["cancer_label"] == cancer]
        sample_count, rho, p_value = safe_spearman(
            cohort["teacher_risk_score"], cohort["rna_clinical_student_risk_score"]
        )
        rows.append(
            {
                "cancer_label": cancer,
                "n": sample_count,
                "spearman_rho": rho,
                "p_value_two_sided": p_value,
            }
        )
    output = pd.DataFrame(rows)
    output["fdr_bh"] = benjamini_hochberg(output["p_value_two_sided"])
    return output


def program_tests(
    frame: pd.DataFrame,
    cancer_order: tuple[str, ...],
    program_columns: dict[str, str | tuple[str, ...]],
    layer: str,
) -> pd.DataFrame:
    rows = []
    for cancer in cancer_order:
        cohort = frame.loc[frame["cancer_label"] == cancer]
        for program_label, configured_columns in program_columns.items():
            score_columns = (
                (configured_columns,)
                if isinstance(configured_columns, str)
                else configured_columns
            )
            if not all(column in cohort.columns for column in score_columns):
                continue
            program_score = cohort[list(score_columns)].mean(axis=1, skipna=True)
            sample_count, rho, p_value = safe_spearman(
                cohort["risk_score"], program_score
            )
            rows.append(
                {
                    "layer": layer,
                    "cancer_label": cancer,
                    "program": program_label,
                    "score_column": ";".join(score_columns),
                    "n": sample_count,
                    "spearman_rho": rho,
                    "p_value_two_sided": p_value,
                }
            )
    output = pd.DataFrame(rows)
    output["fdr_bh_within_layer_cancer"] = output.groupby(
        ["layer", "cancer_label"], group_keys=False
    )["p_value_two_sided"].apply(benjamini_hochberg)
    return output


def main() -> None:
    package_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, default=package_root)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=package_root / "outputs" / "benchmarks" / "figure3_biology",
    )
    args = parser.parse_args()

    # Exact Figure 3 cancers/programs are local and visible in this file.
    cancer_order = ("HNSC", "KIRC", "LUAD")
    rna_program_columns = {
        "PI3K/AKT": "rna_pi3k_akt_mtor_score",
        "Angiogenesis": "rna_angiogenesis_score",
        "EMT/Stromal": ("rna_emt_score", "rna_stromal_ecm_score"),
        "Proliferation": "rna_proliferation_score",
        "Immune": "rna_immune_cytotoxic_score",
        "DNA damage": "rna_dna_damage_score",
    }
    rppa_program_columns = {
        "PI3K/AKT": "rppa_pi3k_akt_mtor",
        "EMT/Stromal": "rppa_emt_stromal",
        "Proliferation": "rppa_proliferation",
        "Immune": "rppa_immune_stromal",
        "DNA damage": "rppa_dna_damage",
    }
    fdr_method = "Benjamini-Hochberg within each layer and cancer"
    all_tests_two_sided = True

    input_root = args.input_root.resolve()
    teacher_student_file = (
        input_root / "outputs/task6/cptac/task6_cptac_teacher_student_risk_scores.csv"
    )
    rna_scores_file = input_root / "outputs/task6/tcga/task6_tcga_rna_pathway_scores.csv"
    rppa_scores_file = input_root / "outputs/task6b/rppa/task6b_rppa_risk_input.csv"
    teacher_student = pd.read_csv(teacher_student_file)
    rna_scores = pd.read_csv(rna_scores_file)
    rppa_scores = pd.read_csv(rppa_scores_file)

    transfer_tests = teacher_student_tests(teacher_student, cancer_order)
    rna_tests = program_tests(
        rna_scores, cancer_order, rna_program_columns, "RNA"
    )
    rppa_tests = program_tests(
        rppa_scores, cancer_order, rppa_program_columns, "RPPA"
    )
    program_results = pd.concat([rna_tests, rppa_tests], ignore_index=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    transfer_tests.to_csv(
        args.output_dir / "figure3_teacher_student_spearman_tests.csv", index=False
    )
    program_results.to_csv(
        args.output_dir / "figure3_program_spearman_tests.csv", index=False
    )
    settings = {
        "cancers": cancer_order,
        "rna_program_columns": rna_program_columns,
        "rppa_program_columns": rppa_program_columns,
        "test": "Spearman rank correlation",
        "all_tests_two_sided": all_tests_two_sided,
        "multiple_testing": fdr_method,
    }
    (args.output_dir / "figure3_biology_test_settings.json").write_text(
        json.dumps(settings, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote Figure 3 biology tests to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
