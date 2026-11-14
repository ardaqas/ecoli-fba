"""Evidence for the error analysis: why does the model call a gene essential?

A knockout fails to grow in FBA when at least one ingredient of the biomass
reaction can no longer be made. Finding out *which* ingredient points to the
pathway the model believes has no alternative.
"""

from collections.abc import Iterable

import cobra
import pandas as pd
from cobra.util.solver import linear_reaction_coefficients

from ecoli_fba.knockout import DEFAULT_THRESHOLD


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


def single_nutrient_rescue(
    model: cobra.Model,
    genes: Iterable[str],
    uptake: float = 10.0,
    threshold: float = DEFAULT_THRESHOLD,
) -> pd.DataFrame:
    """For each gene, the single extra nutrients that let its knockout grow.

    Every exchange reaction that is closed in the current medium is opened in
    turn (one at a time, up to ``uptake``) and the knockout is simulated
    again. A nutrient "rescues" the knockout if growth reaches ``threshold``
    times the wild-type growth on the unsupplemented medium.

    A rescue shows what the knockout is missing. Whether the cells in the
    experiment really had access to that compound is a separate question.

    Returns one row per gene with ``n_rescuing`` and ``rescuing_nutrients``
    (``"; "``-separated exchange reaction IDs).
    """
    wild_type = model.slim_optimize(error_value=0.0)
    if wild_type <= 0:
        raise ValueError("wild type does not grow on this medium")
    closed = [r for r in model.exchanges if r.lower_bound >= 0]
    rows = []
    for gene_id in genes:
        rescuing = []
        with model:
            model.genes.get_by_id(gene_id).knock_out()
            for exchange in closed:
                exchange.lower_bound = -uptake
                if model.slim_optimize(error_value=0.0) >= threshold * wild_type:
                    rescuing.append(exchange.id)
                exchange.lower_bound = 0.0
        rows.append(
            {
                "gene_id": gene_id,
                "n_rescuing": len(rescuing),
                "rescuing_nutrients": "; ".join(rescuing),
            }
        )
    return pd.DataFrame(rows, columns=["gene_id", "n_rescuing", "rescuing_nutrients"])


# Evidence groups. They describe what the model and the data show for a
# mispredicted gene; they are not a verdict on the biological cause.
FP_DUPLICATION = "deletion strain keeps an intact gene copy (Yamamoto 2009)"
FP_RESCUED = "one extra nutrient restores growth in the model"
FP_UNRESCUED = "no single nutrient restores growth in the model"
FP_GROUPS = (FP_DUPLICATION, FP_RESCUED, FP_UNRESCUED)

FN_BACKUP_GENE = "another gene can catalyse all of its reactions"
FN_REROUTED_COST = "reactions lost, model reroutes with slower growth"
FN_REROUTED_FREE = "reactions lost, model growth unaffected"
FN_GROUPS = (FN_BACKUP_GENE, FN_REROUTED_COST, FN_REROUTED_FREE)

UNAFFECTED_RATIO = 0.99


def group_false_positives(table: pd.DataFrame) -> pd.Series:
    """Evidence group of each false positive (model: essential, mutant viable).

    Needs the columns ``duplication_2009`` and ``n_rescuing``. The duplication
    flag takes precedence because it questions the label itself.
    """
    group = pd.Series(FP_UNRESCUED, index=table.index, name="evidence_group")
    group[table["n_rescuing"] > 0] = FP_RESCUED
    group[table["duplication_2009"].astype(bool)] = FP_DUPLICATION
    return group


def group_false_negatives(table: pd.DataFrame) -> pd.Series:
    """Evidence group of each false negative (model: dispensable, no viable mutant).

    Needs the columns ``n_reactions_lost`` and ``growth_ratio``.
    """
    group = pd.Series(FN_REROUTED_FREE, index=table.index, name="evidence_group")
    group[table["growth_ratio"] < UNAFFECTED_RATIO] = FN_REROUTED_COST
    group[table["n_reactions_lost"] == 0] = FN_BACKUP_GENE
    return group
