import pytest

from ecoli_fba import media


def test_lb_list_is_the_published_65_compounds():
    assert len(media.LB_COMPOUNDS) == 65
    assert len(set(media.LB_COMPOUNDS)) == 65


def test_glucose_minimal_medium_is_a_copy(toy_model):
    medium = media.glucose_minimal_medium(toy_model)
    assert medium == {"EX_A_e": 10.0}

    medium["EX_C_e"] = 5.0
    assert toy_model.medium == {"EX_A_e": 10.0}


def test_lb_like_medium_adds_only_compounds_the_model_can_exchange(toy_model, monkeypatch):
    monkeypatch.setattr(media, "LB_COMPOUNDS", ("A", "C", "not_in_model"))
    base = {"EX_A_e": 10.0}

    medium, missing = media.lb_like_medium(toy_model, base, uptake=3.0)

    assert medium == {"EX_A_e": 10.0, "EX_C_e": 3.0}  # base bound kept, new one added
    assert missing == ["not_in_model"]
    assert base == {"EX_A_e": 10.0}  # input not modified


def test_apply_medium_closes_everything_else(toy_model):
    with toy_model:
        media.apply_medium(toy_model, {"EX_C_e": 4.0})
        assert toy_model.reactions.EX_C_e.lower_bound == -4.0
        assert toy_model.reactions.EX_A_e.lower_bound == 0.0

    # reverted by the context manager
    assert toy_model.reactions.EX_A_e.lower_bound == -10.0


def test_apply_medium_rejects_unknown_exchange(toy_model):
    with pytest.raises(KeyError, match="EX_nope_e"):
        media.apply_medium(toy_model, {"EX_nope_e": 1.0})
