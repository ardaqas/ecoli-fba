"""Shared fixtures.

The toy model is small enough to work out every knockout by hand:

    A_e --T_A (gT)--> A_c --R1 (g1) / R1_ALT (g1alt)--> B_c
    B_c --R2 (g2a or g2b)--> C_c --R3 (g3a and g3b)--> D_c --BIOMASS-->
    C_e --T_C (no gene)--> C_c          (only if C is in the medium)
    R3 also makes a by-product W_c, removed by SPONT (pseudo-gene s0001)

On the default medium (A only):
    essential:      gT, g3a, g3b
    non-essential:  g1, g1alt (parallel reactions), g2a, g2b (isozymes)
With C added to the medium, gT is no longer essential.
"""

import cobra
import pytest
from cobra import Metabolite, Model, Reaction


def _reaction(model, rxn_id, stoichiometry, gpr="", bounds=(0.0, 1000.0)):
    rxn = Reaction(rxn_id)
    rxn.bounds = bounds
    model.add_reactions([rxn])
    rxn.add_metabolites(
        {model.metabolites.get_by_id(m): c for m, c in stoichiometry.items()}
    )
    rxn.gene_reaction_rule = gpr
    return rxn


@pytest.fixture
def toy_model() -> cobra.Model:
    model = Model("toy")
    model.add_metabolites(
        [
            Metabolite("A_e", compartment="e"),
            Metabolite("C_e", compartment="e"),
            Metabolite("A_c", compartment="c"),
            Metabolite("B_c", compartment="c"),
            Metabolite("C_c", compartment="c"),
            Metabolite("D_c", compartment="c"),
            Metabolite("W_c", compartment="c"),
        ]
    )
    _reaction(model, "EX_A_e", {"A_e": -1}, bounds=(-10.0, 1000.0))
    _reaction(model, "EX_C_e", {"C_e": -1}, bounds=(0.0, 1000.0))
    _reaction(model, "T_A", {"A_e": -1, "A_c": 1}, gpr="gT")
    _reaction(model, "T_C", {"C_e": -1, "C_c": 1})
    _reaction(model, "R1", {"A_c": -1, "B_c": 1}, gpr="g1")
    _reaction(model, "R1_ALT", {"A_c": -1, "B_c": 1}, gpr="g1alt")
    _reaction(model, "R2", {"B_c": -1, "C_c": 1}, gpr="g2a or g2b")
    _reaction(model, "R3", {"C_c": -1, "D_c": 1, "W_c": 1}, gpr="g3a and g3b")
    _reaction(model, "SPONT", {"W_c": -1}, gpr="s0001")
    _reaction(model, "BIOMASS", {"D_c": -1})
    model.objective = "BIOMASS"
    return model
