"""In-silico growth media.

A medium is a dict ``{exchange_reaction_id: maximum_uptake}`` in mmol/gDW/h,
the same format as ``cobra.Model.medium``.

Two media are used in this project:

* glucose minimal: the medium iML1515 ships with (glucose, oxygen, ammonium
  and inorganic salts).
* LB-like: the glucose minimal medium plus the organic nutrients of a
  published in-silico LB definition. LB is made from tryptone and yeast
  extract and is not chemically defined, so this is an approximation.
"""

import cobra

# In-silico LB (lysogeny broth) composition, BiGG metabolite IDs.
# Source: CarveMe media database, medium "LB" (Machado et al. 2018,
# Nucleic Acids Res. 46:7542), carveme/data/input/media_db.tsv at commit
# 9dbe654 (Apache-2.0). Copied unchanged, 65 compounds.
LB_COMPOUNDS = (
    "adn", "ala__L", "amp", "arg__L", "aso3", "asp__L", "ca2", "cbl1", "cd2",
    "cl", "cmp", "cobalt2", "cro4", "cu2", "cys__L", "dad_2", "dcyt", "fe2",
    "fe3", "fol", "glc__D", "glu__L", "gly", "gmp", "gsn", "h2o", "h2s", "h",
    "hg2", "his__L", "hxan", "ile__L", "ins", "k", "leu__L", "lipoate",
    "lys__L", "met__L", "mg2", "mn2", "mobd", "na1", "nac", "nh4", "ni2",
    "o2", "phe__L", "pheme", "pi", "pnto__R", "pro__L", "pydx", "ribflv",
    "ser__L", "so4", "thm", "thr__L", "thymd", "trp__L", "tyr__L", "ump",
    "ura", "uri", "val__L", "zn2",
)  # fmt: skip

# Uptake limit for nutrients that LB adds to the minimal medium. Real uptake
# rates in LB are not known; the screen is repeated with other values to
# check that the essentiality calls do not depend on this choice.
DEFAULT_LB_UPTAKE = 10.0


def exchange_id(compound: str) -> str:
    """BiGG exchange reaction ID of an extracellular compound."""
    return f"EX_{compound}_e"


def glucose_minimal_medium(model: cobra.Model) -> dict[str, float]:
    """The medium currently set on the model.

    Call this on a freshly loaded iML1515 to get its default aerobic glucose
    minimal medium.
    """
    return dict(model.medium)


def lb_like_medium(
    model: cobra.Model,
    base_medium: dict[str, float],
    uptake: float = DEFAULT_LB_UPTAKE,
) -> tuple[dict[str, float], list[str]]:
    """Add the LB compounds to a base medium.

    Compounds already in ``base_medium`` keep their bound; new ones get
    ``uptake``. LB compounds without an exchange reaction in the model cannot
    be supplied and are returned as the second element.
    """
    medium = dict(base_medium)
    missing = []
    for compound in LB_COMPOUNDS:
        rxn_id = exchange_id(compound)
        if rxn_id not in model.reactions:
            missing.append(compound)
        elif rxn_id not in medium:
            medium[rxn_id] = uptake
    return medium, missing


def apply_medium(model: cobra.Model, medium: dict[str, float]) -> None:
    """Set the medium: listed exchanges open for uptake, all others closed.

    Use inside ``with model:`` to have the change reverted afterwards.
    """
    unknown = [rxn_id for rxn_id in medium if rxn_id not in model.reactions]
    if unknown:
        raise KeyError(f"exchange reactions not in model: {unknown}")
    model.medium = medium
