"""Single-gene knockout screen with flux balance analysis."""

from collections.abc import Iterable

import cobra
import pandas as pd
from cobra.flux_analysis import single_gene_deletion

from ecoli_fba.model import metabolic_genes

# A knockout is called essential when it grows slower than this fraction of
# the wild type. See notebook 02 for how sensitive the calls are to it.
DEFAULT_THRESHOLD = 0.05


def run_knockouts(
    model: cobra.Model,
    genes: Iterable[str] | None = None,
    processes: int = 1,
) -> pd.DataFrame:
    """Delete each gene in turn and maximise growth.

    Uses the medium and objective currently set on the model. Returns one row
    per gene: ``gene_id``, ``growth`` (1/h), ``growth_ratio`` (relative to the
    wild type) and the solver ``status``. An infeasible knockout cannot make
    biomass at all and is reported as zero growth. The full iML1515 screen
    takes about 15 s on one core.
    """
    gene_ids = list(genes) if genes is not None else metabolic_genes(model)
    wild_type = model.slim_optimize(error_value=0.0)
    if wild_type <= 0:
        raise ValueError("wild type does not grow on this medium")

    raw = single_gene_deletion(model, gene_list=gene_ids, processes=processes)
    result = pd.DataFrame(
        {
            "gene_id": [next(iter(ids)) for ids in raw["ids"]],
            "growth": raw["growth"].fillna(0.0).clip(lower=0.0).to_numpy(),
            "status": raw["status"].to_numpy(),
        }
    )
    result["growth_ratio"] = result["growth"] / wild_type
    result = result.set_index("gene_id").loc[gene_ids].reset_index()
    return result[["gene_id", "growth", "growth_ratio", "status"]]


def classify_essential(
    growth_ratio: pd.Series, threshold: float = DEFAULT_THRESHOLD
) -> pd.Series:
    """True where the knockout grows slower than ``threshold`` x wild type."""
    return (growth_ratio < threshold).rename("predicted_essential")


def threshold_sensitivity(
    growth_ratio: pd.Series, thresholds: Iterable[float]
) -> pd.DataFrame:
    """Number of genes called essential at each threshold."""
    thresholds = list(thresholds)
    return pd.DataFrame(
        {
            "threshold": thresholds,
            "n_essential": [int((growth_ratio < t).sum()) for t in thresholds],
        }
    )
