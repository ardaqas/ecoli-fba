"""Per-gene annotation taken from the model: name, pathways and GPR structure."""

import ast
from collections.abc import Iterable

import cobra
import pandas as pd

from ecoli_fba.model import metabolic_genes


def _genes_in_complexes(rxn: cobra.Reaction) -> set[str]:
    """Genes that sit under an ``and`` in the reaction's gene-reaction rule."""
    found: set[str] = set()
    if rxn.gpr.body is None:
        return found
    for node in ast.walk(rxn.gpr.body):
        if isinstance(node, ast.BoolOp) and isinstance(node.op, ast.And):
            found |= {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}
    return found


def annotate_genes(model: cobra.Model, genes: Iterable[str] | None = None) -> pd.DataFrame:
    """One row per gene.

    Columns:
        gene_name          common name, e.g. ``thrA``
        n_reactions        reactions whose rule mentions the gene
        n_reactions_lost   reactions that stop working when only this gene is deleted
        has_backup         at least one of its reactions survives the deletion
                           (an isozyme or alternative complex takes over)
        in_complex         the gene is a subunit of a protein complex (``and``)
        subsystems         model pathways of its reactions, ``"; "``-separated
    """
    gene_ids = list(genes) if genes is not None else metabolic_genes(model)
    rows = []
    for gene_id in gene_ids:
        gene = model.genes.get_by_id(gene_id)
        reactions = sorted(gene.reactions, key=lambda r: r.id)
        lost = [r for r in reactions if not r.gpr.eval(knockouts={gene_id})]
        rows.append(
            {
                "gene_id": gene_id,
                "gene_name": gene.name,
                "n_reactions": len(reactions),
                "n_reactions_lost": len(lost),
                "has_backup": len(lost) < len(reactions),
                "in_complex": any(gene_id in _genes_in_complexes(r) for r in reactions),
                "subsystems": "; ".join(sorted({r.subsystem for r in reactions if r.subsystem})),
            }
        )
    return pd.DataFrame(rows)


def count_by_subsystem(annotation: pd.DataFrame) -> pd.Series:
    """Number of genes per subsystem, most frequent first.

    A gene with reactions in several subsystems is counted once in each, so
    the counts can add up to more than the number of genes.
    """
    counts = (
        annotation["subsystems"]
        .str.split("; ")
        .explode()
        .loc[lambda s: s != ""]
        .value_counts()
    )
    return counts.rename("n_genes").rename_axis("subsystem")
