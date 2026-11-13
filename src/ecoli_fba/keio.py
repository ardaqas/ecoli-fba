"""Experimental essentiality labels from the Keio collection.

Sources (see data/README.md for download details):

* Baba et al. 2006, Mol Syst Biol 2:2006.0008, Supplementary Table 3.
  One row per targeted ORF: how many of the picked colonies carried the
  correct deletion, and the OD600 of the mutant after 22 h in LB and after
  24 h and 48 h in 0.4 % glucose MOPS minimal medium.
* Yamamoto et al. 2009, Mol Syst Biol 5:335, Supplementary Table V.
  Re-examination of the collection: for each mutant, whether its two
  isolates were verified or still carry a second, intact copy of the gene.

Two labels are derived:

``lb_essential``
    No deletion mutant could be obtained on LB agar. This is the paper's
    definition of an essential gene candidate (303 ORFs).
``glucose_essential``
    Essential on LB, or the mutant exists but does not grow in glucose MOPS
    minimal medium (OD600 at 48 h below 0.08, the cut-off of groups VI and VII
    in Table III of the paper: mean - 2 s.d. of all mutants).
"""

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ecoli_fba import DATA_DIR

RAW_DIR = DATA_DIR / "raw" / "keio"


@dataclass(frozen=True)
class RawFile:
    source: str
    sha256: str
    description: str


BABA_TABLE_S3 = "msb4100050-s5.xls"
BABA_TABLE_S6 = "msb4100050-s8.xls"
YAMAMOTO_TABLE_S5 = "msb200992-s5.xls"

RAW_FILES = {
    BABA_TABLE_S3: RawFile(
        source="baba2006",
        sha256="a1a9c30b6cd703ff835f2809b420afe38b1c139865c9af87566a410669a30d2e",
        description="Baba et al. 2006, Supplementary Table 3 (targets, mutants, growth)",
    ),
    BABA_TABLE_S6: RawFile(
        source="baba2006",
        sha256="ded72d1f890d47b0544fa52edcb491b187a010fa6b938437fc04418e08e31638",
        description="Baba et al. 2006, Supplementary Table 6 (essential gene candidates; cross-check only)",
    ),
    YAMAMOTO_TABLE_S5: RawFile(
        source="yamamoto2009",
        sha256="8b7225b7357d4a3a119acf9dc0c59a2e628c56a86b0a9d1a2cf419a6296ff95e",
        description="Yamamoto et al. 2009, Supplementary Table V (mutant re-examination)",
    ),
}

SOURCE_ARCHIVES = {
    "baba2006": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC1681482/supplementaryFiles",
    "yamamoto2009": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2824493/supplementaryFiles",
}

# "Mutants were scored as correct if two or more colonies had the expected
# structure" (Baba et al. 2006).
MIN_CORRECT_COLONIES = 2

# OD600 below which a mutant counts as not growing in glucose MOPS medium.
MOPS_NO_GROWTH_OD = 0.08
MOPS_TIMEPOINT = "od_mops_48h"

_B_NUMBER = r"b\d{4}"

# Position of the columns used from Supplementary Table 3 (the sheet has a
# four-row merged header, so columns are addressed by position).
_BABA_COLUMNS = {
    0: "eck_id",
    1: "gene_name",
    2: "jw_id",
    8: "gene_id",
    14: "n_correct_colonies",
    15: "n_tested_colonies",
    17: "od_lb_22h",
    18: "od_mops_24h",
    19: "od_mops_48h",
    21: "pec",
    28: "comment",
}
_BABA_HEADER_ROWS = 4
_BABA_NUMERIC = ["n_correct_colonies", "n_tested_colonies", "od_lb_22h", "od_mops_24h", "od_mops_48h"]


def read_baba_table_s3(path: Path = RAW_DIR / BABA_TABLE_S3) -> pd.DataFrame:
    """Read and parse Supplementary Table 3 of Baba et al. 2006."""
    return parse_baba_table_s3(pd.read_excel(path, header=None))


def parse_baba_table_s3(raw: pd.DataFrame) -> pd.DataFrame:
    """Tidy the raw sheet of Supplementary Table 3, one row per ORF (JW id).

    Text placeholders in numeric columns ("N.A." for ORFs that were not
    targeted or gave no mutant, "NT" for not tested) become NaN.
    ``targeted`` is False for ORFs the authors did not try to delete.
    """
    header = raw.iloc[:_BABA_HEADER_ROWS].astype(str)
    if "b number" not in header[8].values or "MOPS_48hr" not in header[19].values:
        raise ValueError("sheet does not look like Baba et al. 2006 Supplementary Table 3")

    table = raw.iloc[_BABA_HEADER_ROWS:, list(_BABA_COLUMNS)].copy()
    table.columns = list(_BABA_COLUMNS.values())
    for column in _BABA_NUMERIC:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    for column in ["eck_id", "gene_name", "jw_id", "gene_id"]:
        table[column] = table[column].astype(str).str.strip()
    table["targeted"] = table["n_tested_colonies"].notna()
    table["mutant_obtained"] = table["n_correct_colonies"] >= MIN_CORRECT_COLONIES
    return table.reset_index(drop=True)


def read_baba_table_s6(path: Path = RAW_DIR / BABA_TABLE_S6) -> pd.DataFrame:
    """Supplementary Table 6 of Baba et al. 2006: the essential gene candidates.

    Only used to cross-check the list derived from Supplementary Table 3.
    """
    raw = pd.read_excel(path, header=None)
    table = raw.iloc[_BABA_HEADER_ROWS:, [0, 1, 2, 6]].copy()
    table.columns = ["eck_id", "gene_name", "jw_id", "gene_id"]
    return table.astype(str).apply(lambda column: column.str.strip()).reset_index(drop=True)


def read_yamamoto_table_s5(path: Path = RAW_DIR / YAMAMOTO_TABLE_S5) -> pd.DataFrame:
    """Read and parse Supplementary Table V of Yamamoto et al. 2009."""
    return parse_yamamoto_table_s5(pd.read_excel(path, header=None))


def parse_yamamoto_table_s5(raw: pd.DataFrame) -> pd.DataFrame:
    """Tidy the raw sheet of Supplementary Table V, one row per mutant.

    ``isolate_1`` / ``isolate_2`` are the status of the two independent
    isolates of each mutant: "verified", "duplication" (an intact copy of the
    gene is still present), "mixture" or "use other isolate".
    """
    starts = raw.index[raw[0].astype(str).str.fullmatch(r"ECK\d{4}")]
    if len(starts) == 0 or "Keio collection" not in str(raw.iat[0, 0]):
        raise ValueError("sheet does not look like Yamamoto et al. 2009 Supplementary Table V")

    table = raw.loc[starts[0] :, [0, 1, 2, 3, 7, 8]].copy()
    table.columns = ["eck_id", "gene_id", "jw_id", "gene_name", "isolate_1", "isolate_2"]
    for column in table.columns:
        table[column] = table[column].astype("string").str.strip()
    return table.reset_index(drop=True)


def duplication_genes(yamamoto: pd.DataFrame) -> set[str]:
    """Genes whose two isolates both still carry an intact copy of the gene.

    Yamamoto et al. 2009 list these 25 genes as candidates for essential
    genes: the "deletion mutant" exists only because the gene was duplicated.
    """
    both = (yamamoto["isolate_1"] == "duplication") & (yamamoto["isolate_2"] == "duplication")
    return set(yamamoto.loc[both, "gene_id"])


def build_labels(baba: pd.DataFrame, yamamoto: pd.DataFrame | None = None) -> pd.DataFrame:
    """One row per gene (b-number) with the two essentiality labels.

    A few genes were targeted as two separate ORFs. Such a gene counts as
    deleted if either ORF gave a mutant, and its OD is the highest one
    measured, because one growing deletion mutant is enough to show the gene
    is dispensable.

    Labels are missing (``pd.NA``) for genes that were never targeted, and
    ``glucose_essential`` is also missing if a mutant exists but its growth in
    MOPS medium was not measured.
    """
    genes = baba[baba["gene_id"].str.fullmatch(_B_NUMBER)]
    od_columns = ["od_lb_22h", "od_mops_24h", "od_mops_48h"]
    labels = genes.groupby("gene_id", sort=True).agg(
        gene_name=("gene_name", "first"),
        jw_ids=("jw_id", "; ".join),
        targeted=("targeted", "any"),
        mutant_obtained=("mutant_obtained", "any"),
        **{column: (column, "max") for column in od_columns},
        pec=("pec", "first"),
    )

    lb_essential = (~labels["mutant_obtained"]).astype("boolean").where(labels["targeted"])
    no_growth = (labels[MOPS_TIMEPOINT] < MOPS_NO_GROWTH_OD).astype("boolean")
    no_growth = no_growth.where(labels[MOPS_TIMEPOINT].notna())
    labels["lb_essential"] = lb_essential
    labels["glucose_essential"] = (lb_essential | no_growth).where(labels["targeted"])

    if yamamoto is not None:
        labels["duplication_2009"] = labels.index.isin(duplication_genes(yamamoto))
    return labels.reset_index()
