import pytest
from cobra import Metabolite, Reaction

from ecoli_fba import errors


@pytest.fixture
def two_precursor_model(toy_model):
    """Toy model whose biomass needs D_c and a second ingredient E_c (made by gE)."""
    toy_model.add_metabolites([Metabolite("E_c", compartment="c")])
    make_e = Reaction("R_E", lower_bound=0.0, upper_bound=1000.0)
    toy_model.add_reactions([make_e])
    make_e.add_metabolites({toy_model.metabolites.A_c: -1, toy_model.metabolites.E_c: 1})
    make_e.gene_reaction_rule = "gE"
    toy_model.reactions.BIOMASS.add_metabolites({toy_model.metabolites.E_c: -1})
    return toy_model


def test_biomass_reaction_is_the_objective(toy_model):
    assert errors.biomass_reaction(toy_model).id == "BIOMASS"


def test_blocked_precursors_match_hand_derivation(two_precursor_model):
    table = errors.blocked_biomass_precursors(
        two_precursor_model, ["gE", "g3a", "gT", "g2a"]
    ).set_index("gene_id")

    assert table.loc["gE", "blocked_precursors"] == "E_c"  # only its own product
    assert table.loc["g3a", "blocked_precursors"] == "D_c"
    assert table.loc["gT", "blocked_precursors"] == "D_c; E_c"  # no A uptake: nothing can be made
    assert table.loc["gT", "n_blocked"] == 2
    assert table.loc["g2a", "blocked_precursors"] == ""  # isozyme g2b takes over
    assert table.loc["g2a", "n_blocked"] == 0


def test_model_is_left_unchanged(two_precursor_model):
    n_reactions = len(two_precursor_model.reactions)
    growth = two_precursor_model.slim_optimize()

    errors.blocked_biomass_precursors(two_precursor_model, ["gE", "gT"])

    assert len(two_precursor_model.reactions) == n_reactions
    assert errors.biomass_reaction(two_precursor_model).id == "BIOMASS"
    assert two_precursor_model.slim_optimize() == pytest.approx(growth)
    assert all(g.functional for g in two_precursor_model.genes)


def test_single_nutrient_rescue_matches_hand_derivation(toy_model):
    table = errors.single_nutrient_rescue(toy_model, ["gT", "g3a"]).set_index("gene_id")

    # C enters downstream of the A transporter, so it rescues gT ...
    assert table.loc["gT", "rescuing_nutrients"] == "EX_C_e"
    assert table.loc["gT", "n_rescuing"] == 1
    # ... but not g3a, which acts after C.
    assert table.loc["g3a", "rescuing_nutrients"] == ""
    assert table.loc["g3a", "n_rescuing"] == 0


def test_rescue_screen_leaves_medium_unchanged(toy_model):
    before = dict(toy_model.medium)

    errors.single_nutrient_rescue(toy_model, ["gT"])

    assert toy_model.medium == before
    assert all(g.functional for g in toy_model.genes)
