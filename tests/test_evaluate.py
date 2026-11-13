import math

import pandas as pd
import pytest
from sklearn.metrics import matthews_corrcoef

from ecoli_fba import evaluate


def test_outcomes_use_essential_as_positive_class():
    predicted = pd.Series([True, True, False, False])
    observed = pd.Series([True, False, True, False])

    outcome = evaluate.classify_outcomes(predicted, observed)

    assert outcome.tolist() == ["TP", "FP", "FN", "TN"]


def test_outcomes_reject_misaligned_series():
    with pytest.raises(ValueError, match="same index"):
        evaluate.classify_outcomes(pd.Series([True], index=["a"]), pd.Series([True], index=["b"]))


def test_metrics_against_hand_calculation():
    # 10 genes: TP=3, FP=1, FN=2, TN=4
    counts = {"TP": 3, "FP": 1, "FN": 2, "TN": 4}

    m = evaluate.metrics(counts)

    assert m["n_genes"] == 10
    assert m["n_observed_essential"] == 5
    assert m["n_predicted_essential"] == 4
    assert m["accuracy"] == pytest.approx(7 / 10)
    assert m["precision"] == pytest.approx(3 / 4)
    assert m["recall"] == pytest.approx(3 / 5)
    assert m["specificity"] == pytest.approx(4 / 5)
    assert m["balanced_accuracy"] == pytest.approx((3 / 5 + 4 / 5) / 2)
    assert m["f1"] == pytest.approx(6 / 9)
    # (3*4 - 1*2) / sqrt(4 * 5 * 5 * 6)
    assert m["mcc"] == pytest.approx(10 / math.sqrt(600))


def test_mcc_agrees_with_scikit_learn():
    predicted = pd.Series([True] * 3 + [True] + [False] * 2 + [False] * 4)
    observed = pd.Series([True] * 3 + [False] + [True] * 2 + [False] * 4)

    counts = evaluate.confusion_counts(evaluate.classify_outcomes(predicted, observed))

    assert evaluate.metrics(counts)["mcc"] == pytest.approx(matthews_corrcoef(observed, predicted))


def test_undefined_metrics_are_nan_not_errors():
    m = evaluate.metrics({"TP": 0, "FP": 0, "FN": 2, "TN": 5})

    assert math.isnan(m["precision"])
    assert math.isnan(m["mcc"])
    assert m["recall"] == 0.0


def test_compare_keeps_only_genes_with_a_known_label():
    predictions = pd.DataFrame(
        {"gene_id": ["fake001", "fake002", "fake003", "fake004"],
         "predicted_essential": [True, True, False, False]}
    )  # fmt: skip
    labels = pd.DataFrame(
        {"gene_id": ["fake002", "fake001", "fake003", "fake999"],
         "observed_essential": [False, True, None, True]}
    )  # fmt: skip

    merged = evaluate.compare(predictions, labels)

    # fake003 has no label, fake004 is not in the label table, fake999 is not in the model
    assert merged["gene_id"].tolist() == ["fake001", "fake002"]
    assert merged["outcome"].tolist() == ["TP", "FP"]


def test_compare_rejects_duplicate_genes():
    predictions = pd.DataFrame({"gene_id": ["fake001", "fake001"], "predicted_essential": [True, False]})
    labels = pd.DataFrame({"gene_id": ["fake001"], "observed_essential": [True]})

    with pytest.raises(ValueError, match="duplicate gene_id in predictions"):
        evaluate.compare(predictions, labels)
