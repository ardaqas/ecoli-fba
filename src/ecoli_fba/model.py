"""Loading the iML1515 genome-scale metabolic model."""

from pathlib import Path
from urllib.request import urlretrieve

import cobra
from cobra.io import load_json_model

from ecoli_fba import DATA_DIR

MODEL_URL = "https://bigg.ucsd.edu/static/models/iML1515.json"
MODEL_PATH = DATA_DIR / "iML1515.json"

# s0001 is not a real gene: iML1515 uses it as the "gene" of spontaneous
# reactions. It has no knockout strain, so it is excluded from every screen.
PSEUDO_GENES = frozenset({"s0001"})


def load_model(path: Path = MODEL_PATH, download: bool = True) -> cobra.Model:
    """Load iML1515 from disk, downloading it from BiGG Models first if needed."""
    path = Path(path)
    if not path.exists():
        if not download:
            raise FileNotFoundError(f"{path} not found and download=False")
        path.parent.mkdir(parents=True, exist_ok=True)
        urlretrieve(MODEL_URL, path)
    return load_json_model(path)


def metabolic_genes(model: cobra.Model) -> list[str]:
    """IDs of the real genes in the model (pseudo-genes removed)."""
    return [g.id for g in model.genes if g.id not in PSEUDO_GENES]
