import pandas as pd

from ecoli_fba import annotation


def test_gpr_structure_of_toy_genes(toy_model):
    toy_model.reactions.R2.subsystem = "Pathway X"
    toy_model.reactions.R3.subsystem = "Pathway Y"
    table = annotation.annotate_genes(toy_model).set_index("gene_id")

    # isozymes: the reaction survives losing either one
    assert table.loc["g2a", ["n_reactions", "n_reactions_lost"]].tolist() == [1, 0]
    assert table.loc["g2a", "has_backup"]
    assert not table.loc["g2a", "in_complex"]
    assert table.loc["g2a", "subsystems"] == "Pathway X"

    # complex subunits: losing either one kills the reaction
    assert table.loc["g3a", "n_reactions_lost"] == 1
    assert not table.loc["g3a", "has_backup"]
    assert table.loc["g3a", "in_complex"]

    # single gene with a parallel reaction: its own reaction is lost,
    # the backup is a different reaction, not an isozyme
    assert table.loc["g1", "n_reactions_lost"] == 1
    assert not table.loc["g1", "has_backup"]
    assert not table.loc["g1", "in_complex"]


def test_complex_inside_isozyme_rule(toy_model):
    toy_model.reactions.R2.gene_reaction_rule = "g2a or (g2b and g3a)"
    table = annotation.annotate_genes(toy_model).set_index("gene_id")

    assert not table.loc["g2a", "in_complex"]
    assert table.loc["g2b", "in_complex"]
    assert table.loc["g2b", "has_backup"]  # g2a still catalyses R2
    assert table.loc["g3a", ["n_reactions", "n_reactions_lost"]].tolist() == [2, 1]


def test_count_by_subsystem_counts_a_gene_once_per_subsystem():
    table = pd.DataFrame(
        {"gene_id": ["x", "y", "z"], "subsystems": ["P1; P2", "P1", ""]}
    )

    counts = annotation.count_by_subsystem(table)

    assert counts.to_dict() == {"P1": 2, "P2": 1}
    assert counts.index[0] == "P1"
