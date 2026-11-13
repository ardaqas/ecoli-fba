"""Tests for the Keio label parser.

The parser tests build a small table with the same layout as the real
supplementary sheets but with obviously fake identifiers, in memory, at test
time. No synthetic label data is stored in the repository. The tests at the
bottom run on the real files and are skipped until
``python scripts/fetch_keio.py`` has been run.
"""

import numpy as np
import pandas as pd
import pytest

from ecoli_fba import keio

needs_raw_data = pytest.mark.skipif(
    not all((keio.RAW_DIR / name).is_file() for name in keio.RAW_FILES),
    reason="Keio tables not downloaded (run scripts/fetch_keio.py)",
)


def fake_baba_sheet(rows: list[dict]) -> pd.DataFrame:
    """Raw sheet shaped like Supplementary Table 3: 4 header rows, 29 columns."""
    sheet = pd.DataFrame(np.nan, index=range(4 + len(rows)), columns=range(29), dtype=object)
    sheet.iat[0, 0] = "Supplementary Table 3. Information on Keio collection deletion mutants"
    sheet.iat[3, 2], sheet.iat[3, 8] = "JW id", "b number"
    sheet.iat[3, 17], sheet.iat[3, 18], sheet.iat[3, 19] = "LB_22hr", "MOPS_24hr", "MOPS_48hr"
    for i, row in enumerate(rows, start=4):
        sheet.iat[i, 0] = row.get("eck", f"FAKE{i:04d}")
        sheet.iat[i, 1] = row["name"]
        sheet.iat[i, 2] = row["jw"]
        sheet.iat[i, 8] = row["b"]
        sheet.iat[i, 14] = row["verified"]
        sheet.iat[i, 15] = row["tested"]
        sheet.iat[i, 17] = row.get("lb", "N.A.")
        sheet.iat[i, 18] = row.get("m24", "N.A.")
        sheet.iat[i, 19] = row.get("m48", "N.A.")
        sheet.iat[i, 21] = row.get("pec", "u")
        sheet.iat[i, 28] = row.get("comment", np.nan)
    return sheet


FAKE_ROWS = [
    # grows everywhere
    dict(name="fakA", jw="JW9001", b="b9001", verified=8, tested=8, lb=0.80, m24=0.30, m48=0.35),
    # mutant exists but does not grow on minimal medium (auxotroph)
    dict(name="fakB", jw="JW9002", b="b9002", verified=6, tested=8, lb=0.75, m24=0.01, m48=0.03),
    # no mutant: essential candidate
    dict(name="fakC", jw="JW9003", b="b9003", verified=0, tested=8),
    # a single correct colony is below the paper's "two or more" rule
    dict(name="fakD", jw="JW9004", b="b9004", verified=1, tested=8),
    # never targeted: no label at all
    dict(name="fakE", jw="JW9005", b="b9005", verified="N.A.", tested="N.A.", comment="not targeted (IS ORF)"),
    # one gene targeted as two ORFs, only one of which gave a mutant
    dict(name="fakF", jw="JW9006", b="b9006", verified=0, tested=8),
    dict(name="fakF", jw="JW9007", b="b9006", verified=4, tested=8, lb=0.60, m24=0.25, m48=0.20),
    # mutant exists, 48 h reading missing
    dict(name="fakG", jw="JW9008", b="b9008", verified=8, tested=8, lb=0.70, m24=0.30, m48="NT"),
    # ORF without an MG1655 b-number
    dict(name="fakH", jw="JW9009", b="none", verified=8, tested=8, lb=0.70, m24=0.30, m48=0.30),
    # exactly at the cut-off counts as growth (the rule is strictly "below")
    dict(name="fakI", jw="JW9010", b="b9010", verified=8, tested=8, lb=0.70, m24=0.05, m48=0.08),
]


@pytest.fixture
def fake_baba() -> pd.DataFrame:
    return keio.parse_baba_table_s3(fake_baba_sheet(FAKE_ROWS))


def test_parse_keeps_one_row_per_orf_and_converts_placeholders(fake_baba):
    assert len(fake_baba) == len(FAKE_ROWS)
    by_jw = fake_baba.set_index("jw_id")

    assert by_jw.loc["JW9001", "od_mops_48h"] == pytest.approx(0.35)
    assert np.isnan(by_jw.loc["JW9003", "od_lb_22h"])  # "N.A."
    assert np.isnan(by_jw.loc["JW9008", "od_mops_48h"])  # "NT"
    assert not by_jw.loc["JW9005", "targeted"]
    assert by_jw.loc["JW9003", "targeted"]


def test_mutant_needs_two_correct_colonies(fake_baba):
    by_jw = fake_baba.set_index("jw_id")["mutant_obtained"]

    assert by_jw["JW9001"]
    assert not by_jw["JW9003"]
    assert not by_jw["JW9004"]


def test_parse_rejects_a_different_sheet():
    with pytest.raises(ValueError, match="Supplementary Table 3"):
        keio.parse_baba_table_s3(pd.DataFrame(np.nan, index=range(6), columns=range(29), dtype=object))


def test_labels_per_gene(fake_baba):
    labels = keio.build_labels(fake_baba).set_index("gene_id")

    assert "none" not in labels.index
    assert labels.index.is_unique

    expected = {
        # gene:    (lb_essential, glucose_essential)
        "b9001": (False, False),
        "b9002": (False, True),
        "b9003": (True, True),
        "b9004": (True, True),
        "b9005": (pd.NA, pd.NA),
        "b9006": (False, False),
        "b9008": (False, pd.NA),
        "b9010": (False, False),
    }
    for gene, (lb, glucose) in expected.items():
        for column, value in [("lb_essential", lb), ("glucose_essential", glucose)]:
            got = labels.loc[gene, column]
            assert (pd.isna(got) and pd.isna(value)) or got == value, (gene, column, got)


def test_split_orf_gene_uses_the_growing_mutant(fake_baba):
    row = keio.build_labels(fake_baba).set_index("gene_id").loc["b9006"]

    assert row["jw_ids"] == "JW9006; JW9007"
    assert row["mutant_obtained"]
    assert row["od_mops_48h"] == pytest.approx(0.20)


def test_duplication_flag_requires_both_isolates(fake_baba):
    sheet = pd.DataFrame(
        [
            ["Supplementary Table V  Information on Keio collection deletion mutants"] + [np.nan] * 9,
            ["gene id"] + [np.nan] * 9,
            ["ECK9001", "b9001", "JW9001", "fakA", "1, 2", "A", 1, "verified", "verified", np.nan],
            ["ECK9002", "b9002", "JW9002", "fakB", "1, 2", "A", 2, "duplication", "duplication", np.nan],
            ["ECK9010", "b9010", "JW9010", "fakI", "1, 2", "A", 3, "duplication", "verified", np.nan],
        ]
    )
    yamamoto = keio.parse_yamamoto_table_s5(sheet)

    assert len(yamamoto) == 3
    assert keio.duplication_genes(yamamoto) == {"b9002"}

    labels = keio.build_labels(fake_baba, yamamoto).set_index("gene_id")
    assert labels["duplication_2009"].to_dict()["b9002"]
    assert labels["duplication_2009"].sum() == 1


# --- contract checks on the real files ---------------------------------------


@needs_raw_data
def test_real_table_reproduces_the_numbers_in_the_paper():
    baba = keio.read_baba_table_s3()

    # Baba et al. 2006: 4288 ORFs targeted, 3985 mutants, 303 essential candidates.
    assert baba["targeted"].sum() == 4288
    assert baba["mutant_obtained"].sum() == 3985
    assert (baba["targeted"] & ~baba["mutant_obtained"]).sum() == 303
    assert baba["jw_id"].is_unique


@needs_raw_data
def test_derived_essential_list_equals_supplementary_table_6():
    baba = keio.read_baba_table_s3()
    derived = set(baba.loc[baba["targeted"] & ~baba["mutant_obtained"], "jw_id"])

    assert derived == set(keio.read_baba_table_s6()["jw_id"])


@needs_raw_data
def test_mops_cutoff_matches_table_iii_of_the_paper():
    # Table III groups mutants by mean +/- 2 s.d.; its lower MOPS bound is 0.08.
    mutants = keio.read_baba_table_s3().query("mutant_obtained")
    lower = mutants["od_mops_48h"].mean() - 2 * mutants["od_mops_48h"].std()

    assert round(lower, 2) == keio.MOPS_NO_GROWTH_OD


@needs_raw_data
def test_real_labels_are_well_formed():
    labels = keio.build_labels(keio.read_baba_table_s3(), keio.read_yamamoto_table_s5())

    assert labels["gene_id"].is_unique
    assert labels["gene_id"].str.fullmatch(r"b\d{4}").all()
    # Yamamoto et al. 2009: partial duplications in both isolates of 25 mutants.
    assert labels["duplication_2009"].sum() == 25
    # every LB-essential gene is also essential on minimal medium
    essential = labels[labels["lb_essential"].fillna(False)]
    assert essential["glucose_essential"].all()
