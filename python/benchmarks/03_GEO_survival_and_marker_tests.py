#!/usr/bin/env python3
"""Recompute the GEO survival tests and marker-shift tests used in Figures 4 and 6."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
from scipy.stats import ttest_ind


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


def survival_statistics(
    frame: pd.DataFrame,
    time_column: str,
    event_column: str,
    score_column: str,
) -> dict[str, float | int]:
    analysis = frame[[time_column, event_column, score_column]].dropna().copy()
    analysis.columns = ["time", "event", "score"]
    standard_deviation = float(analysis["score"].std(ddof=0))
    analysis["score_z"] = (
        analysis["score"] - float(analysis["score"].mean())
    ) / standard_deviation
    cox = CoxPHFitter(penalizer=0.0)
    cox.fit(analysis[["time", "event", "score_z"]], duration_col="time", event_col="event")
    cox_row = cox.summary.loc["score_z"]

    median_score = float(analysis["score"].median())
    high = analysis["score"] > median_score
    logrank = logrank_test(
        analysis.loc[high, "time"],
        analysis.loc[~high, "time"],
        event_observed_A=analysis.loc[high, "event"],
        event_observed_B=analysis.loc[~high, "event"],
    )
    group_frame = analysis.assign(high_score=high.astype(int))
    group_cox = CoxPHFitter(penalizer=0.0)
    group_cox.fit(
        group_frame[["time", "event", "high_score"]],
        duration_col="time",
        event_col="event",
    )
    return {
        "n": len(analysis),
        "events": int(analysis["event"].sum()),
        "c_index": float(concordance_index(analysis["time"], -analysis["score"], analysis["event"])),
        "cox_hr_per_sd": float(cox_row["exp(coef)"]),
        "cox_ci_lower": float(cox_row["exp(coef) lower 95%"]),
        "cox_ci_upper": float(cox_row["exp(coef) upper 95%"]),
        "cox_p_two_sided": float(cox_row["p"]),
        "median_high_vs_low_hr": float(np.exp(group_cox.params_["high_score"])),
        "logrank_p_two_sided": float(logrank.p_value),
    }


def standardized_marker_tests(
    frame: pd.DataFrame,
    dataset_label: str,
    marker_genes: tuple[str, ...],
    positive_markers: tuple[str, ...],
) -> pd.DataFrame:
    analysis = frame.loc[frame["risk_group"].isin(("low", "high"))].copy()
    rows = []
    for gene in marker_genes:
        if gene not in analysis.columns:
            continue
        values = pd.to_numeric(analysis[gene], errors="coerce")
        standardized = (values - values.mean()) / values.std(ddof=1)
        high = standardized.loc[analysis["risk_group"] == "high"].dropna()
        low = standardized.loc[analysis["risk_group"] == "low"].dropna()
        raw_shift = float(high.mean() - low.mean()) if len(high) > 0 and len(low) > 0 else float("nan")
        if len(high) >= 2 and len(low) >= 2:
            test_statistic, test_p_value = ttest_ind(
                high, low, equal_var=False, alternative="two-sided"
            )
        else:
            test_statistic, test_p_value = float("nan"), float("nan")
        rows.append(
            {
                "dataset": dataset_label,
                "gene_symbol": gene,
                "high_n": len(high),
                "low_n": len(low),
                "high_minus_low_z_shift": raw_shift,
                "direction_consistent_shift": raw_shift if gene in positive_markers else -raw_shift,
                "welch_t_statistic": float(test_statistic),
                "p_value_two_sided": float(test_p_value),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    package_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, default=package_root)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=package_root / "outputs" / "benchmarks" / "geo",
    )
    args = parser.parse_args()

    # Locked manuscript values are explicit within this script.
    geo_accessions = ("GSE31210", "GSE68465", "GSE41613", "GSE65858", "GSE29609")
    kirc_analysis_sets = ("primary_all", "primary_m0_sensitivity")
    marker_genes = ("NAA38", "RAN", "MRPL42", "CDC45", "MAF", "FMO1", "ZBTB4", "STAT5B")
    positive_markers = ("NAA38", "RAN", "MRPL42", "CDC45")
    negative_markers = ("MAF", "FMO1", "ZBTB4", "STAT5B")
    confidence_alpha = 0.05
    all_tests_two_sided = True

    input_root = args.input_root.resolve()
    geo_file = input_root / "outputs/task4/geo_external_validation/plan2_geo_marker_score_predictions.csv"
    kirc_file = input_root / "outputs/task4/geo_kirc_dfs_progression/plan2_kirc_gse22541_predictions.csv"
    geo = pd.read_csv(geo_file)
    kirc = pd.read_csv(kirc_file)

    survival_rows = []
    marker_frames = []
    for accession in geo_accessions:
        cohort = geo.loc[geo["accession"] == accession].copy()
        cancer_label = str(cohort["cancer_label"].dropna().iloc[0])
        survival_rows.append(
            {
                "dataset": accession,
                "cancer_label": cancer_label,
                "endpoint": "OS",
                **survival_statistics(cohort, "OS_time", "OS", "promise_marker_score"),
            }
        )
        marker_frames.append(
            standardized_marker_tests(
                cohort,
                f"{cancer_label} {accession}",
                marker_genes,
                positive_markers,
            )
        )

    for analysis_set in kirc_analysis_sets:
        cohort = kirc.loc[kirc["analysis_set"] == analysis_set].copy()
        dataset_label = (
            "KIRC GSE22541 primary"
            if analysis_set == "primary_all"
            else "KIRC GSE22541 M0-only"
        )
        survival_rows.append(
            {
                "dataset": analysis_set,
                "cancer_label": "KIRC",
                "endpoint": "DFS/progression",
                **survival_statistics(
                    cohort,
                    "dfs_time_months",
                    "dfs_event",
                    "promise_marker_score",
                ),
            }
        )
        marker_frames.append(
            standardized_marker_tests(
                cohort,
                dataset_label,
                marker_genes,
                positive_markers,
            )
        )

    survival_results = pd.DataFrame(survival_rows)
    marker_results = pd.concat(marker_frames, ignore_index=True)
    marker_results["fdr_bh"] = marker_results.groupby("dataset", group_keys=False)[
        "p_value_two_sided"
    ].apply(benjamini_hochberg)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    survival_results.to_csv(args.output_dir / "geo_survival_statistics.csv", index=False)
    marker_results.to_csv(args.output_dir / "geo_marker_shift_tests.csv", index=False)
    settings = {
        "geo_accessions": geo_accessions,
        "kirc_analysis_sets": kirc_analysis_sets,
        "positive_markers": positive_markers,
        "negative_markers": negative_markers,
        "confidence_alpha": confidence_alpha,
        "all_tests_two_sided": all_tests_two_sided,
        "marker_test": "Welch two-sample t test within each dataset; BH FDR within dataset",
    }
    (args.output_dir / "geo_benchmark_settings.json").write_text(
        json.dumps(settings, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote GEO benchmark outputs to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
