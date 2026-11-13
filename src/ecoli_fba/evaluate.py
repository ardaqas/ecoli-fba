"""Comparing predicted essentiality with experimental labels.

The positive class is always "essential":

    TP  predicted essential,      observed essential
    FP  predicted essential,      observed non-essential
    FN  predicted non-essential,  observed essential
    TN  predicted non-essential,  observed non-essential
"""

import math

import pandas as pd

OUTCOMES = ("TP", "FP", "FN", "TN")


def classify_outcomes(predicted: pd.Series, observed: pd.Series) -> pd.Series:
    """TP / FP / FN / TN for each gene, from two aligned boolean series."""
    if not predicted.index.equals(observed.index):
        raise ValueError("predicted and observed must have the same index")
    predicted = predicted.astype(bool)
    observed = observed.astype(bool)
    outcome = pd.Series("TN", index=predicted.index, name="outcome")
    outcome[predicted & observed] = "TP"
    outcome[predicted & ~observed] = "FP"
    outcome[~predicted & observed] = "FN"
    return outcome


def confusion_counts(outcome: pd.Series) -> dict[str, int]:
    """Number of genes in each confusion-matrix cell."""
    counts = outcome.value_counts()
    return {name: int(counts.get(name, 0)) for name in OUTCOMES}


def _ratio(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else math.nan


def metrics(counts: dict[str, int]) -> dict[str, float]:
    """Summary metrics from confusion-matrix counts.

    Essential genes are a minority, so accuracy alone flatters a model that
    calls everything non-essential. Precision, recall and the Matthews
    correlation coefficient (MCC) are reported next to it. A metric whose
    denominator is zero is NaN.
    """
    tp, fp, fn, tn = (counts[k] for k in OUTCOMES)
    n = tp + fp + fn + tn
    recall = _ratio(tp, tp + fn)
    specificity = _ratio(tn, tn + fp)
    precision = _ratio(tp, tp + fp)
    mcc_denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return {
        "n_genes": n,
        "n_observed_essential": tp + fn,
        "n_predicted_essential": tp + fp,
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
        "accuracy": _ratio(tp + tn, n),
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
        "balanced_accuracy": (recall + specificity) / 2,
        "f1": _ratio(2 * tp, 2 * tp + fp + fn),
        "mcc": _ratio(tp * tn - fp * fn, mcc_denominator),
    }


def compare(
    predictions: pd.DataFrame,
    labels: pd.DataFrame,
    predicted_column: str = "predicted_essential",
    observed_column: str = "observed_essential",
) -> pd.DataFrame:
    """Join predictions and labels on ``gene_id`` and add the outcome.

    Only genes present in both tables with a known label are kept; genes whose
    label is missing (NaN) cannot be scored and are dropped.
    """
    for name, table in [("predictions", predictions), ("labels", labels)]:
        if table["gene_id"].duplicated().any():
            raise ValueError(f"duplicate gene_id in {name}")
    merged = predictions.merge(
        labels[["gene_id", observed_column]].dropna(subset=[observed_column]),
        on="gene_id",
        how="inner",
    )
    merged[observed_column] = merged[observed_column].astype(bool)
    merged["outcome"] = classify_outcomes(merged[predicted_column], merged[observed_column])
    return merged
