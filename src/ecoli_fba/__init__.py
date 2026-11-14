"""Gene essentiality prediction in E. coli iML1515 with flux balance analysis."""

from pathlib import Path

# Repository root, so code works the same from notebooks/, tests/ and scripts/,
# whatever the current working directory. This relies on the editable install
# (`pip install -e .`), which keeps the package inside the repository.
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
RESULTS_DIR = REPO_ROOT / "results"
