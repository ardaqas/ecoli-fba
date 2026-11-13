"""Evidence for the error analysis: why does the model call a gene essential?

A knockout fails to grow in FBA when at least one ingredient of the biomass
reaction can no longer be made. Finding out *which* ingredient points to the
pathway the model believes has no alternative.
"""

from collections.abc import Iterable

import cobra
import pandas as pd
from cobra.util.solver import linear_reaction_coefficients


def biomass_reaction(model: cobra.Model) -> cobra.Reaction:
    """The reaction the model currently maximises."""
    (reaction,) = linear_reaction_coefficients(model)
    return reaction


def blocked_biomass_precursors(
    model: cobra.Model, genes: Iterable[str], tolerance: float = 1e-6
) -> pd.DataFrame:
    """For each gene, the biomass ingredients its knockout can no longer make.

    Every metabolite consumed by the biomass reaction gets a temporary drain
    ("demand") reaction. With the gene deleted, each drain is maximised on its
    own; if nothing can flow, that ingredient cannot be produced. Uses the
    medium currently set on the model and leaves the model unchanged.

    Returns one row per gene with ``n_blocked`` and ``blocked_precursors``
    (``"; "``-separated metabolite IDs, empty if the knockout can still make
    everything).
    """
    precursors = [m for m, coefficient in biomass_reaction(model).metabolites.items() if coefficient < 0]
    rows = []
    with model:
        demands = {m.id: model.add_boundary(m, type="demand") for m in precursors}
        for gene_id in genes:
            blocked = []
            with model:
                model.genes.get_by_id(gene_id).knock_out()
                for metabolite_id, demand in demands.items():
                    model.objective = demand
                    if model.slim_optimize(error_value=0.0) < tolerance:
                        blocked.append(metabolite_id)
            rows.append(
                {
                    "gene_id": gene_id,
                    "n_blocked": len(blocked),
                    "blocked_precursors": "; ".join(sorted(blocked)),
                }
            )
    return pd.DataFrame(rows, columns=["gene_id", "n_blocked", "blocked_precursors"])
