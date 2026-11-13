import pandas as pd
import pytest

from ecoli_fba import knockout, media
from ecoli_fba.model import MODEL_PATH, load_model, metabolic_genes

ESSENTIAL_ON_A = {"gT", "g3a", "g3b"}
ALL_TOY_GENES = {"gT", "g1", "g1alt", "g2a", "g2b", "g3a", "g3b"}


def essential_set(result: pd.DataFrame) -> set[str]:
    return set(result.loc[knockout.classify_essential(result["growth_ratio"]), "gene_id"])


def test_pseudo_gene_is_not_screened(toy_model):
    assert "s0001" in toy_model.genes
    assert set(metabolic_genes(toy_model)) == ALL_TOY_GENES


def test_toy_model_essential_genes_match_hand_derivation(toy_model):
    result = knockout.run_knockouts(toy_model, processes=1)

    assert set(result["gene_id"]) == ALL_TOY_GENES
    assert essential_set(result) == ESSENTIAL_ON_A


def test_isozymes_and_parallel_reactions_keep_full_growth(toy_model):
    result = knockout.run_knockouts(toy_model, processes=1).set_index("gene_id")

    for gene in ["g1", "g1alt", "g2a", "g2b"]:
        assert result.loc[gene, "growth_ratio"] == pytest.approx(1.0)
    for gene in ESSENTIAL_ON_A:
        assert result.loc[gene, "growth"] == pytest.approx(0.0)


def test_essentiality_depends_on_medium(toy_model):
    with toy_model:
        media.apply_medium(toy_model, {"EX_A_e": 10.0, "EX_C_e": 10.0})
        result = knockout.run_knockouts(toy_model, processes=1)

    # C enters downstream of the A transporter, so gT is rescued.
    assert essential_set(result) == {"g3a", "g3b"}


def test_rows_follow_requested_gene_order(toy_model):
    genes = ["g3b", "g1", "gT"]
    result = knockout.run_knockouts(toy_model, genes=genes, processes=1)

    assert result["gene_id"].tolist() == genes
    assert list(result.columns) == ["gene_id", "growth", "growth_ratio", "status"]


def test_knockouts_do_not_modify_the_model(toy_model):
    before = toy_model.slim_optimize()
    knockout.run_knockouts(toy_model, processes=1)

    assert toy_model.slim_optimize() == pytest.approx(before)
    assert all(g.functional for g in toy_model.genes)


def test_screen_refuses_medium_without_growth(toy_model):
    with toy_model:
        media.apply_medium(toy_model, {})
        with pytest.raises(ValueError, match="wild type does not grow"):
            knockout.run_knockouts(toy_model, processes=1)


def test_classify_uses_strict_threshold():
    ratios = pd.Series([0.0, 0.049, 0.05, 0.5, 1.0])

    assert knockout.classify_essential(ratios, threshold=0.05).tolist() == [
        True, True, False, False, False,
    ]  # fmt: skip


def test_threshold_sensitivity_counts():
    ratios = pd.Series([0.0, 0.0, 0.3, 0.8, 1.0])

    table = knockout.threshold_sensitivity(ratios, (t for t in [0.01, 0.5, 0.9]))

    assert table["threshold"].tolist() == [0.01, 0.5, 0.9]
    assert table["n_essential"].tolist() == [2, 3, 4]


@pytest.mark.slow
@pytest.mark.skipif(not MODEL_PATH.exists(), reason="data/iML1515.json not downloaded")
def test_full_model_screen():
    model = load_model(download=False)

    # Same wild-type growth rate as notebook 01.
    assert model.slim_optimize() == pytest.approx(0.877, abs=1e-3)

    result = knockout.run_knockouts(model, processes=1)
    assert len(result) == 1515
    assert "s0001" not in set(result["gene_id"])
    assert result["growth_ratio"].between(0, 1 + 1e-6).all()
