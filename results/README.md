# Results

Every file here is written by one of the notebooks; rerunning notebooks 02 to 04 in order reproduces them (byte for byte with the software versions listed in the main README). Genes are identified by their MG1655 b-number (`gene_id`) and common name (`gene_name`). Growth rates are in h⁻¹; `growth_ratio` is the knockout growth rate divided by the wild-type growth rate on the same medium.

## Tables

| File | Notebook | Content |
|---|---|---|
| `knockout_glucose_minimal.csv` | 02 | One row per gene (1,515): predicted growth, `growth_ratio`, solver status and `predicted_essential` (ratio below 0.05) on glucose minimal medium |
| `knockout_lb_like.csv` | 02 | The same on the LB-like medium |
| `threshold_sensitivity.csv` | 02 | Number of predicted essential genes for essentiality thresholds from 0.000001 to 0.9 |
| `lb_uptake_sensitivity.csv` | 02 | Wild-type growth and essential genes on LB-like medium for uptake limits of 1, 10 and 100 mmol/gDW/h |
| `predicted_essential_by_subsystem.csv` | 02 | Predicted essential genes per model pathway, essential on both media or only on glucose minimal medium |
| `metrics.csv` | 03 | Confusion matrix and metrics (accuracy, precision, recall, specificity, balanced accuracy, F1, MCC) for comparisons A, B and C |
| `comparison_glucose.csv` | 03 | Per-gene prediction, label and outcome (TP/FP/FN/TN) for comparison C (1,508 genes) |
| `comparison_lb.csv` | 03 | The same for comparison B |
| `mops_cutoff_sensitivity.csv` | 03 | Comparison C recomputed with other no-growth cut-offs (OD600 0.04 to 0.15) and with the 24 h reading |
| `errors_glucose.csv` | 04 | One row per false positive and false negative of comparison C, with all evidence columns (see below) |
| `errors_lb.csv` | 04 | The same for comparison B |
| `error_groups.csv` | 04 | Number of mispredicted genes per evidence group, comparisons B and C |
| `scenarios.csv` | 04 | What-if scenarios: confusion counts, precision, recall, MCC, and the genes whose errors were resolved or newly created |

Evidence columns of `errors_*.csv`:

| Column | Source | Meaning |
|---|---|---|
| `n_reactions`, `n_reactions_lost`, `lost_reactions` | model | reactions of the gene, and those that stop working when only this gene is deleted |
| `has_backup`, `backup_genes` | model | other genes in the rules of the reactions that survive the deletion (isoenzymes, alternative complexes) |
| `in_complex` | model | the gene is a subunit of a complex (`and` in a gene rule) |
| `subsystems` | model | pathways of the gene's reactions |
| `blocked_precursors` | model, false positives only | biomass components the knockout can no longer make |
| `rescuing_nutrients` | model, false positives only | single nutrients that restore growth of the knockout when added to the medium |
| `od_lb_22h`, `od_mops_24h`, `od_mops_48h` | Baba et al. 2006 | OD600 of the Keio mutant in LB and in glucose MOPS medium (empty if no mutant was obtained) |
| `pec` | Baba et al. 2006 | essentiality call of the PEC database as listed in the Keio table (E essential, N non-essential, u unknown) |
| `duplication_2009` | Yamamoto et al. 2009 | both isolates of the mutant still carry an intact copy of the gene |
| `evidence_group` | derived | the evidence group used in the error analysis |

## Figures

| File | Notebook | Shows |
|---|---|---|
| `knockout_growth_ratio.png` | 02 | Distribution of knockout growth ratios on both media |
| `predicted_essential_by_subsystem.png` | 02 | Predicted essential genes per pathway |
| `keio_growth_lb_vs_mops.png` | 03 | OD600 of Keio mutants in LB and in glucose MOPS medium |
| `validation_metrics.png` | 03 | Precision, recall and MCC of comparisons A, B and C |
| `error_evidence_groups.png` | 04 | Evidence groups of false positives and false negatives |

The figures are explained in the [main README](../README.md#results).
