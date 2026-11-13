# E. coli gene essentiality: Flux Balance Analysis vs. experiment

> Work in progress.

## Question
How well does a genome-scale metabolic model (iML1515) predict which *E. coli* K-12 genes are essential,
using only network structure and Flux Balance Analysis (FBA)? Where does it fail, and why?

## Approach
1. Load the iML1515 genome-scale metabolic model with COBRApy and explore its structure.
2. Simulate growth under different conditions (carbon source, oxygen).
3. Run single-gene knockout simulations and classify genes as essential / non-essential.
4. Validate predictions against the Keio collection (Baba et al., 2006).
5. Analyse false positives / false negatives and categorise the causes
   (regulation, missing pathways, annotation errors, medium mismatch, label noise).

## Repository structure
```
notebooks/   one notebook per stage
src/         reusable analysis code (the `ecoli_fba` package)
tests/       pytest tests for the package
data/        external data (Keio essentiality list); see data/README.md
results/     figures and tables
```

## Setup
With conda:
```bash
conda env create -f environment.yml
conda activate ecoli-fba
jupyter lab
```
Or with a plain virtual environment (Python 3.10 or newer):
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
jupyter lab
```
Both install the code in `src/ecoli_fba/` as an editable package, so notebooks and tests can `import ecoli_fba`.

Run the tests with `pytest` (fast unit tests) or `pytest -m slow` (full-model tests, need `data/iML1515.json`).
The first run of `notebooks/01_model_exploration.ipynb` downloads `iML1515.json` from BiGG Models (https://bigg.ucsd.edu/static/models/iML1515.json) into `data/` and needs internet access.

## Results
_To be filled in as the analysis progresses._

## Limitations
_To be filled in. FBA assumes steady state, ignores gene regulation and enzyme kinetics, and depends on the chosen objective function and medium._

## References
- Monk et al. (2017) *Nature Biotechnology*: iML1515 model.
- Baba et al. (2006) *Molecular Systems Biology*: Keio collection.
- Ebrahim et al. (2013) *BMC Systems Biology*: COBRApy.
