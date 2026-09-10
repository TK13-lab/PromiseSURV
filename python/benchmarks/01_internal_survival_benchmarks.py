#!/usr/bin/env python3
"""Run the locked CPTAC model benchmark, repeated splits, and paired tests."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy.stats import wilcoxon
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler


def feature_groups(frame: pd.DataFrame) -> dict[str, list[str]]:
    id_columns = {"cancer_label", "cohort_id", "patient_id", "OS", "OS_time", "source"}
    numeric = [
        column
        for column in frame.columns
        if column not in id_columns and pd.api.types.is_numeric_dtype(frame[column])
    ]
    clinical = [
        column
        for column in numeric
        if not column.startswith(("rna__", "proteomics__", "phosphoproteomics__"))
    ]
    rna = [column for column in numeric if column.startswith("rna__")]
    protein = [column for column in numeric if column.startswith("proteomics__")]
    phosphoprotein = [
        column for column in numeric if column.startswith("phosphoproteomics__")
    ]
    return {
        "clinical_cox": clinical,
        "rna_cox": rna,
        "rna_clinical_cox": clinical + rna,
        "protein_cox": protein,
        "phospho_cox": phosphoprotein,
        "protein_phospho_cox": protein + phosphoprotein,
        "full_multiomics_cox": clinical + rna + protein + phosphoprotein,
    }


def split_discovery_data(
    frame: pd.DataFrame,
    test_fraction: float,
    random_seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    strata = frame["cancer_label"].astype(str) + "_event" + frame["OS"].astype(int).astype(str)
    stratify = strata if strata.value_counts().min() >= 2 else frame["cancer_label"]
    train_index, test_index = train_test_split(
        frame.index,
        test_size=test_fraction,
        random_state=random_seed,
        stratify=stratify,
    )
    return frame.loc[train_index].copy(), frame.loc[test_index].copy()


def keep_variable_features(
    train: pd.DataFrame,
    features: list[str],
    minimum_variance: float,
) -> list[str]:
    kept = []
    for feature in features:
        values = pd.to_numeric(train[feature], errors="coerce")
        if values.nunique(dropna=True) >= 2 and float(values.var(skipna=True)) > minimum_variance:
            kept.append(feature)
    if not kept:
        raise RuntimeError("Every requested feature was removed by variance filtering.")
    return kept


def c_index(frame: pd.DataFrame, risk_score: np.ndarray) -> float:
    return float(concordance_index(frame["OS_time"], -risk_score, frame["OS"]))


def fit_cox(
    train: pd.DataFrame,
    test: pd.DataFrame,
    features: list[str],
    penalizer: float,
    l1_ratio: float,
    minimum_variance: float,
) -> tuple[dict[str, float | int], np.ndarray, np.ndarray, list[str]]:
    kept_features = keep_variable_features(train, features, minimum_variance)
    scaler = StandardScaler()
    train_x = pd.DataFrame(
        scaler.fit_transform(train[kept_features]), columns=kept_features, index=train.index
    )
    test_x = pd.DataFrame(
        scaler.transform(test[kept_features]), columns=kept_features, index=test.index
    )
    train_fit = pd.concat([train[["OS_time", "OS"]], train_x], axis=1)
    test_fit = pd.concat([test[["OS_time", "OS"]], test_x], axis=1)
    model = CoxPHFitter(penalizer=penalizer, l1_ratio=l1_ratio)
    model.fit(train_fit, duration_col="OS_time", event_col="OS", show_progress=False)
    train_risk = model.predict_log_partial_hazard(train_fit[kept_features]).to_numpy()
    test_risk = model.predict_log_partial_hazard(test_fit[kept_features]).to_numpy()
    metrics: dict[str, float | int] = {
        "train_c_index": c_index(train, train_risk),
        "test_c_index": c_index(test, test_risk),
        "feature_count": len(kept_features),
        "train_n": len(train),
        "test_n": len(test),
        "train_events": int(train["OS"].sum()),
        "test_events": int(test["OS"].sum()),
    }
    return metrics, train_risk, test_risk, kept_features


def fit_distilled_student(
    train: pd.DataFrame,
    test: pd.DataFrame,
    student_features: list[str],
    train_teacher_risk: np.ndarray,
    test_teacher_risk: np.ndarray,
    ridge_alpha: float,
    minimum_variance: float,
) -> tuple[dict[str, float | int], np.ndarray, np.ndarray]:
    kept_features = keep_variable_features(train, student_features, minimum_variance)
    scaler = StandardScaler()
    train_x = scaler.fit_transform(train[kept_features])
    test_x = scaler.transform(test[kept_features])
    student = Ridge(alpha=ridge_alpha)
    student.fit(train_x, train_teacher_risk)
    train_risk = student.predict(train_x)
    test_risk = student.predict(test_x)
    metrics: dict[str, float | int] = {
        "train_c_index": c_index(train, train_risk),
        "test_c_index": c_index(test, test_risk),
        "feature_count": len(kept_features),
        "train_n": len(train),
        "test_n": len(test),
        "train_events": int(train["OS"].sum()),
        "test_events": int(test["OS"].sum()),
        "teacher_test_c_index": c_index(test, test_teacher_risk),
        "test_teacher_student_mse": float(np.mean((test_risk - test_teacher_risk) ** 2)),
    }
    return metrics, train_risk, test_risk


def bootstrap_c_index(
    predictions: pd.DataFrame,
    iterations: int,
    random_seed: int,
    alpha: float,
) -> pd.DataFrame:
    random_generator = np.random.default_rng(random_seed)
    rows = []
    for model_name, model_frame in predictions.groupby("model", sort=True):
        values = model_frame[["OS_time", "OS", "risk_score"]].to_numpy()
        bootstrap_scores = []
        for _ in range(iterations):
            sampled = values[random_generator.integers(0, len(values), size=len(values))]
            if len(np.unique(sampled[:, 1])) < 2:
                continue
            try:
                bootstrap_scores.append(
                    concordance_index(sampled[:, 0], -sampled[:, 2], sampled[:, 1])
                )
            except ZeroDivisionError:
                continue
        lower, upper = np.quantile(bootstrap_scores, [alpha / 2, 1 - alpha / 2])
        rows.append(
            {
                "model": model_name,
                "bootstrap_iterations_requested": iterations,
                "bootstrap_iterations_used": len(bootstrap_scores),
                "c_index_mean": float(np.mean(bootstrap_scores)),
                "ci_lower": float(lower),
                "ci_upper": float(upper),
                "ci_level": 1 - alpha,
            }
        )
    return pd.DataFrame(rows)


def bootstrap_repeated_split_means(
    repeated_metrics: pd.DataFrame,
    iterations: int,
    random_seed: int,
    alpha: float,
) -> pd.DataFrame:
    random_generator = np.random.default_rng(random_seed)
    rows = []
    for model_name, model_frame in repeated_metrics.groupby("model", sort=True):
        split_scores = model_frame["test_c_index"].dropna().to_numpy()
        bootstrap_means = np.asarray(
            [
                random_generator.choice(split_scores, size=len(split_scores), replace=True).mean()
                for _ in range(iterations)
            ]
        )
        lower, upper = np.quantile(bootstrap_means, [alpha / 2, 1 - alpha / 2])
        rows.append(
            {
                "model": model_name,
                "repeated_splits": len(split_scores),
                "bootstrap_iterations": iterations,
                "mean_c_index": float(split_scores.mean()),
                "bootstrap_mean_c_index": float(bootstrap_means.mean()),
                "ci_lower": float(lower),
                "ci_upper": float(upper),
                "ci_level": 1 - alpha,
            }
        )
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


def paired_model_tests(
    repeated_metrics: pd.DataFrame,
    reference_model: str,
    comparator_models: tuple[str, ...],
    bootstrap_iterations: int,
    random_seed: int,
    alpha: float,
) -> pd.DataFrame:
    pivot = repeated_metrics.pivot(index="repeat", columns="model", values="test_c_index")
    rows = []
    for comparator_index, comparator in enumerate(comparator_models):
        paired = pivot[[reference_model, comparator]].dropna()
        differences = (paired[reference_model] - paired[comparator]).to_numpy()
        random_generator = np.random.default_rng(random_seed + comparator_index)
        bootstrap_mean_differences = np.asarray(
            [
                random_generator.choice(differences, size=len(differences), replace=True).mean()
                for _ in range(bootstrap_iterations)
            ]
        )
        difference_lower, difference_upper = np.quantile(
            bootstrap_mean_differences, [alpha / 2, 1 - alpha / 2]
        )
        result = wilcoxon(
            paired[reference_model],
            paired[comparator],
            alternative="two-sided",
            zero_method="wilcox",
            method="auto",
        )
        rows.append(
            {
                "reference_model": reference_model,
                "comparator_model": comparator,
                "paired_splits": len(paired),
                "mean_c_index_difference": float(np.mean(differences)),
                "median_c_index_difference": float(np.median(differences)),
                "bootstrap_iterations": bootstrap_iterations,
                "difference_ci_lower": float(difference_lower),
                "difference_ci_upper": float(difference_upper),
                "wilcoxon_statistic": float(result.statistic),
                "p_value_two_sided": float(result.pvalue),
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
        default=package_root / "outputs" / "benchmarks" / "internal",
    )
    parser.add_argument(
        "--quick-check",
        action="store_true",
        help="Run one split and 20 bootstrap samples to verify an installation; full locked values remain the default.",
    )
    args = parser.parse_args()

    # Locked values are local and explicit; no shared/global configuration is read.
    random_seed = 20260623
    test_fraction = 0.30
    repeated_splits = 10
    repeated_mean_bootstrap_iterations = 5000
    patient_bootstrap_iterations = 2000
    confidence_alpha = 0.05
    cox_penalizer = 0.50
    cox_l1_ratio = 0.00
    ridge_alpha = 10.0
    minimum_variance = 1e-8
    student_model_name = "rna_clinical_student_ridge_distilled_from_full_teacher"
    comparator_models = ("clinical_cox", "rna_clinical_cox", "full_multiomics_cox")
    if args.quick_check:
        repeated_splits = 1
        repeated_mean_bootstrap_iterations = 20
        patient_bootstrap_iterations = 20

    input_file = (
        args.input_root.resolve()
        / "outputs"
        / "task3"
        / "model_inputs"
        / "cptac_mvp_complete_case_input.parquet"
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_parquet(input_file)
    groups = feature_groups(data)

    metric_rows: list[dict[str, object]] = []
    fixed_split_predictions: list[pd.DataFrame] = []
    for repeat_index in range(repeated_splits):
        split_seed = random_seed + repeat_index
        train, test = split_discovery_data(data, test_fraction, split_seed)
        split_risk: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for model_name, model_features in groups.items():
            print(
                f"Split {repeat_index + 1}/{repeated_splits}: fitting {model_name}",
                flush=True,
            )
            metrics, train_risk, test_risk, _ = fit_cox(
                train,
                test,
                model_features,
                cox_penalizer,
                cox_l1_ratio,
                minimum_variance,
            )
            metric_rows.append(
                {"repeat": repeat_index + 1, "seed": split_seed, "model": model_name, **metrics}
            )
            split_risk[model_name] = (train_risk, test_risk)

        teacher_train_risk, teacher_test_risk = split_risk["full_multiomics_cox"]
        print(
            f"Split {repeat_index + 1}/{repeated_splits}: fitting {student_model_name}",
            flush=True,
        )
        student_metrics, student_train_risk, student_test_risk = fit_distilled_student(
            train,
            test,
            groups["rna_clinical_cox"],
            teacher_train_risk,
            teacher_test_risk,
            ridge_alpha,
            minimum_variance,
        )
        metric_rows.append(
            {
                "repeat": repeat_index + 1,
                "seed": split_seed,
                "model": student_model_name,
                **student_metrics,
            }
        )
        split_risk[student_model_name] = (student_train_risk, student_test_risk)

        if repeat_index == 0:
            for model_name, (_, test_risk) in split_risk.items():
                fixed_split_predictions.append(
                    test[["cancer_label", "patient_id", "OS_time", "OS"]]
                    .copy()
                    .assign(model=model_name, risk_score=test_risk)
                )

    repeated_metrics = pd.DataFrame(metric_rows)
    repeated_summary = (
        repeated_metrics.groupby("model", as_index=False)
        .agg(
            repeats=("test_c_index", "count"),
            train_c_index_mean=("train_c_index", "mean"),
            test_c_index_mean=("test_c_index", "mean"),
            test_c_index_sd=("test_c_index", "std"),
            test_c_index_median=("test_c_index", "median"),
            test_c_index_q025=("test_c_index", lambda values: values.quantile(0.025)),
            test_c_index_q975=("test_c_index", lambda values: values.quantile(0.975)),
            feature_count_median=("feature_count", "median"),
        )
        .sort_values("test_c_index_mean", ascending=False)
    )
    fixed_predictions = pd.concat(fixed_split_predictions, ignore_index=True)
    patient_bootstrap = bootstrap_c_index(
        fixed_predictions,
        patient_bootstrap_iterations,
        random_seed,
        confidence_alpha,
    )
    repeated_mean_bootstrap = bootstrap_repeated_split_means(
        repeated_metrics,
        repeated_mean_bootstrap_iterations,
        random_seed,
        confidence_alpha,
    )
    paired_tests = paired_model_tests(
        repeated_metrics,
        student_model_name,
        comparator_models,
        repeated_mean_bootstrap_iterations,
        random_seed,
        confidence_alpha,
    )

    repeated_metrics.to_csv(args.output_dir / "internal_repeated_split_metrics.csv", index=False)
    repeated_summary.to_csv(args.output_dir / "internal_repeated_split_summary.csv", index=False)
    fixed_predictions.to_csv(args.output_dir / "internal_fixed_split_predictions.csv", index=False)
    patient_bootstrap.to_csv(
        args.output_dir / "internal_fixed_split_patient_bootstrap_cindex.csv", index=False
    )
    repeated_mean_bootstrap.to_csv(
        args.output_dir / "internal_repeated_split_bootstrap_mean_cindex.csv", index=False
    )
    paired_tests.to_csv(args.output_dir / "internal_paired_wilcoxon_tests.csv", index=False)
    settings = {
        "random_seed": random_seed,
        "test_fraction": test_fraction,
        "repeated_splits": repeated_splits,
        "repeated_mean_bootstrap_iterations": repeated_mean_bootstrap_iterations,
        "patient_bootstrap_iterations": patient_bootstrap_iterations,
        "confidence_alpha": confidence_alpha,
        "cox_penalizer": cox_penalizer,
        "cox_l1_ratio": cox_l1_ratio,
        "ridge_alpha": ridge_alpha,
        "minimum_variance": minimum_variance,
        "quick_check": args.quick_check,
    }
    (args.output_dir / "internal_benchmark_settings.json").write_text(
        json.dumps(settings, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote internal benchmark outputs to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
