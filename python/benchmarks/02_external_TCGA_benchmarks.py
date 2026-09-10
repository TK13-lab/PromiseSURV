#!/usr/bin/env python3
"""Recompute external TCGA survival metrics and paired model comparisons."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
from sklearn.metrics import roc_auc_score


def c_index(frame: pd.DataFrame) -> float:
    if frame["OS"].nunique() < 2:
        return float("nan")
    return float(concordance_index(frame["OS_time"], -frame["risk_score"], frame["OS"]))


def rank_to_unit_interval(values: pd.Series) -> np.ndarray:
    if len(values) == 1:
        return np.array([0.5])
    ranks = values.rank(method="average").to_numpy()
    return (ranks - 1.0) / (len(values) - 1.0)


def horizon_metrics(frame: pd.DataFrame, horizon_days: int) -> dict[str, float | int]:
    usable = frame.loc[
        (frame["OS_time"] > horizon_days)
        | ((frame["OS_time"] <= horizon_days) & (frame["OS"] == 1))
    ].copy()
    outcome = ((usable["OS_time"] <= horizon_days) & (usable["OS"] == 1)).astype(int)
    auc = float("nan") if outcome.nunique() < 2 else float(roc_auc_score(outcome, usable["risk_score"]))
    rank_risk = rank_to_unit_interval(usable["risk_score"])
    brier = float(np.mean((outcome.to_numpy() - rank_risk) ** 2))
    return {
        f"auc_{horizon_days}d": auc,
        f"rank_brier_{horizon_days}d": brier,
        f"usable_{horizon_days}d": len(usable),
        f"events_{horizon_days}d": int(outcome.sum()),
    }


def cox_and_logrank(frame: pd.DataFrame) -> dict[str, float | int]:
    analysis = frame[["OS_time", "OS", "risk_score"]].dropna().copy()
    risk_standard_deviation = float(analysis["risk_score"].std(ddof=0))
    analysis["risk_z"] = (
        analysis["risk_score"] - float(analysis["risk_score"].mean())
    ) / risk_standard_deviation
    cox = CoxPHFitter(penalizer=0.0)
    cox.fit(analysis[["OS_time", "OS", "risk_z"]], duration_col="OS_time", event_col="OS")
    cox_row = cox.summary.loc["risk_z"]

    median_risk = float(analysis["risk_score"].median())
    high = analysis["risk_score"] > median_risk
    logrank = logrank_test(
        analysis.loc[high, "OS_time"],
        analysis.loc[~high, "OS_time"],
        event_observed_A=analysis.loc[high, "OS"],
        event_observed_B=analysis.loc[~high, "OS"],
    )
    group_frame = analysis.assign(high_risk=high.astype(int))
    group_cox = CoxPHFitter(penalizer=0.0)
    group_cox.fit(
        group_frame[["OS_time", "OS", "high_risk"]],
        duration_col="OS_time",
        event_col="OS",
    )
    return {
        "n": len(analysis),
        "events": int(analysis["OS"].sum()),
        "cox_hr_per_sd": float(cox_row["exp(coef)"]),
        "cox_ci_lower": float(cox_row["exp(coef) lower 95%"]),
        "cox_ci_upper": float(cox_row["exp(coef) upper 95%"]),
        "cox_p_two_sided": float(cox_row["p"]),
        "median_high_vs_low_hr": float(np.exp(group_cox.params_["high_risk"])),
        "logrank_p_two_sided": float(logrank.p_value),
    }


def summarize_predictions(
    predictions: pd.DataFrame,
    cancer_order: tuple[str, ...],
    horizons_days: tuple[int, ...],
) -> pd.DataFrame:
    rows = []
    for model_name, model_frame in predictions.groupby("model", sort=True):
        for cancer in (*cancer_order, "Pooled"):
            subset = (
                model_frame.copy()
                if cancer == "Pooled"
                else model_frame.loc[model_frame["cancer_label"] == cancer].copy()
            )
            row: dict[str, object] = {
                "model": model_name,
                "cancer_label": cancer,
                "c_index": c_index(subset),
                **cox_and_logrank(subset),
            }
            for horizon in horizons_days:
                row.update(horizon_metrics(subset, horizon))
            rows.append(row)
    return pd.DataFrame(rows)


def holm_adjust(p_values: list[float]) -> list[float]:
    count = len(p_values)
    order = np.argsort(p_values)
    adjusted = np.empty(count, dtype=float)
    running_maximum = 0.0
    for rank, index in enumerate(order):
        candidate = min(1.0, (count - rank) * p_values[index])
        running_maximum = max(running_maximum, candidate)
        adjusted[index] = running_maximum
    return adjusted.tolist()


def paired_bootstrap_difference(
    predictions: pd.DataFrame,
    reference_model: str,
    comparator_model: str,
    cancer_order: tuple[str, ...],
    iterations: int,
    random_seed: int,
    confidence_alpha: float,
) -> pd.DataFrame:
    rows = []
    for cancer_index, cancer in enumerate(cancer_order):
        reference = predictions.loc[
            (predictions["model"] == reference_model)
            & (predictions["cancer_label"] == cancer),
            ["patient_id", "OS_time", "OS", "risk_score"],
        ].rename(columns={"risk_score": "reference_risk"})
        comparator = predictions.loc[
            (predictions["model"] == comparator_model)
            & (predictions["cancer_label"] == cancer),
            ["patient_id", "risk_score"],
        ].rename(columns={"risk_score": "comparator_risk"})
        paired = reference.merge(comparator, on="patient_id", validate="one_to_one")
        random_generator = np.random.default_rng(random_seed + cancer_index)
        differences = []
        for _ in range(iterations):
            sampled = paired.iloc[random_generator.integers(0, len(paired), size=len(paired))]
            if sampled["OS"].nunique() < 2:
                continue
            try:
                reference_c = concordance_index(
                    sampled["OS_time"], -sampled["reference_risk"], sampled["OS"]
                )
                comparator_c = concordance_index(
                    sampled["OS_time"], -sampled["comparator_risk"], sampled["OS"]
                )
                differences.append(reference_c - comparator_c)
            except ZeroDivisionError:
                continue
        lower, upper = np.quantile(differences, [confidence_alpha / 2, 1 - confidence_alpha / 2])
        two_sided_p = min(
            1.0,
            2 * min(np.mean(np.asarray(differences) <= 0), np.mean(np.asarray(differences) >= 0)),
        )
        rows.append(
            {
                "cancer_label": cancer,
                "reference_model": reference_model,
                "comparator_model": comparator_model,
                "paired_patients": len(paired),
                "bootstrap_iterations_requested": iterations,
                "bootstrap_iterations_used": len(differences),
                "mean_c_index_difference": float(np.mean(differences)),
                "ci_lower": float(lower),
                "ci_upper": float(upper),
                "p_value_two_sided": float(two_sided_p),
            }
        )
    output = pd.DataFrame(rows)
    output["p_value_holm"] = holm_adjust(output["p_value_two_sided"].tolist())
    return output


def main() -> None:
    package_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, default=package_root)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=package_root / "outputs" / "benchmarks" / "external_tcga",
    )
    args = parser.parse_args()

    # Locked values are intentionally repeated here instead of imported globally.
    cancer_order = ("LUAD", "KIRC", "HNSC")
    horizons_days = (365, 1095, 1825)
    paired_bootstrap_iterations = 2000
    random_seed = 20260623
    confidence_alpha = 0.05
    reference_model = "promise_rna_clinical_student"
    comparator_model = "common_clinical_cox"

    input_root = args.input_root.resolve()
    student_file = (
        input_root
        / "outputs/task4/external_tcga/task4_tcga_rna_clinical_student_predictions.csv"
    )
    clinical_file = (
        input_root
        / "outputs/task4/external_tcga/task4_tcga_common_clinical_predictions.csv"
    )
    student = pd.read_csv(student_file)
    student["model"] = reference_model
    clinical = pd.read_csv(clinical_file)
    clinical["model"] = comparator_model
    predictions = pd.concat([student, clinical], ignore_index=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    metrics = summarize_predictions(predictions, cancer_order, horizons_days)
    paired_tests = paired_bootstrap_difference(
        predictions,
        reference_model,
        comparator_model,
        cancer_order,
        paired_bootstrap_iterations,
        random_seed,
        confidence_alpha,
    )
    metrics.to_csv(args.output_dir / "external_tcga_survival_metrics.csv", index=False)
    paired_tests.to_csv(
        args.output_dir / "external_tcga_paired_bootstrap_tests.csv", index=False
    )
    settings = {
        "cancers": cancer_order,
        "horizons_days": horizons_days,
        "paired_bootstrap_iterations": paired_bootstrap_iterations,
        "random_seed": random_seed,
        "confidence_alpha": confidence_alpha,
        "reference_model": reference_model,
        "comparator_model": comparator_model,
        "all_p_values_are_two_sided": True,
    }
    (args.output_dir / "external_tcga_benchmark_settings.json").write_text(
        json.dumps(settings, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote external TCGA benchmark outputs to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
