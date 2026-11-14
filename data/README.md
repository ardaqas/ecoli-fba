# Data

Nothing in this folder except this file is committed to git. Everything is
downloaded by code, and each file is checked against a SHA-256 so that every
run uses exactly the same inputs.

## iML1515 model

| | |
|---|---|
| File | `iML1515.json` |
| Source | BiGG Models, https://bigg.ucsd.edu/static/models/iML1515.json |
| Reference | Monk et al. (2017) *Nat Biotechnol* 35:904–908 |
| How | Downloaded on first use by notebook 01 or `ecoli_fba.model.load_model()` |

## Keio collection tables

Download with:

```bash
python scripts/fetch_keio.py
```

The files land in `raw/keio/`. The 2006 archive is served slowly by Europe PMC
(about six minutes); the script prints its progress.

| File | Content | Used for |
|---|---|---|
| `msb4100050-s5.xls` | Baba et al. 2006, Supplementary Table 3: every targeted ORF, number of correct deletion colonies, OD600 after 22 h in LB and after 24 h / 48 h in glucose MOPS medium | essentiality labels |
| `msb4100050-s8.xls` | Baba et al. 2006, Supplementary Table 6: the 303 essential gene candidates | cross-check of the list derived from Table 3 |
| `msb200992-s5.xls` | Yamamoto et al. 2009, Supplementary Table V: PCR re-examination of both isolates of each mutant | flagging mutants that still carry an intact copy of the gene |

Sources (retrieved 2026-11-14 through the Europe PMC REST API):

- Baba T. et al. (2006) Construction of *Escherichia coli* K-12 in-frame, single-gene knockout mutants: the Keio collection. *Mol Syst Biol* 2:2006.0008. PMC1681482.
  `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC1681482/supplementaryFiles`
- Yamamoto N. et al. (2009) Update on the Keio collection of *Escherichia coli* single-gene deletion mutants. *Mol Syst Biol* 5:335. PMC2824493.
  `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2824493/supplementaryFiles`

SHA-256:

```
a1a9c30b6cd703ff835f2809b420afe38b1c139865c9af87566a410669a30d2e  msb4100050-s5.xls
ded72d1f890d47b0544fa52edcb491b187a010fa6b938437fc04418e08e31638  msb4100050-s8.xls
8b7225b7357d4a3a119acf9dc0c59a2e628c56a86b0a9d1a2cf419a6296ff95e  msb200992-s5.xls
```

The raw tables stay out of the repository because the 2006 article carries a
publisher copyright notice and no open licence statement.

## How the labels are defined

The experimental conditions decide what "essential" means, so two labels are
derived (code: `src/ecoli_fba/keio.py`).

| Label | Definition | Experimental condition |
|---|---|---|
| `lb_essential` | The ORF was targeted and fewer than two colonies carried the correct deletion. This is the paper's own rule and reproduces its 303 candidates exactly. | Mutants selected on LB agar with kanamycin, 37 °C, aerobic |
| `glucose_essential` | `lb_essential`, or the mutant exists but its OD600 after 48 h in glucose MOPS medium is below 0.08 | 0.4 % glucose MOPS medium, 96-well plate, 37 °C, no shaking |

The 0.08 cut-off is the one the paper uses in its Table III (groups VI and VII);
it equals the mean minus two standard deviations of the 48 h readings of all
mutants, which the test suite recomputes from the data.

Genes that were never targeted have no label and are left out of every
comparison. Seven of them are in the model (`araA`, `araB`, `araD`, `rhaA`,
`rhaB`, `rhaD`, `lacZ`): they are already deleted or disrupted in BW25113, the
parent strain of the collection.
