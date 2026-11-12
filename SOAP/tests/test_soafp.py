# -*- coding: utf-8 -*-
"""Tests for the soafp module.  Run from SOAP/:  python -m pytest tests"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import soafp  # noqa: E402
from soafp import yields  # noqa: E402

rdkit = pytest.importorskip("rdkit")


# --- yield tables ---------------------------------------------------------

def test_tables_cover_c5_to_c30_for_every_column():
    for case in yields.CASES:
        t = yields.table(case)
        assert set(t) == set(yields.COLUMNS)
        for col in yields.COLUMNS:
            assert sorted(t[col]) == list(range(5, 31))


def test_avg_is_mean_of_cases_and_matches_original_avg_sheet():
    v, flags = yields.mass_yield("n-Alkane", 5)
    assert v["avg"] == pytest.approx((v["noVWL"] + v["VWL"]) / 2)
    assert v["avg"] == pytest.approx(0.000990985, rel=1e-5)   # avg_py, n-Alkane, C5
    assert flags == []


def test_out_of_range_carbon_numbers():
    v, flags = yields.mass_yield("Aromatic", 4)
    assert v == {"noVWL": 0.0, "VWL": 0.0, "avg": 0.0} and "below" in flags[0]
    v35, flags = yields.mass_yield("Aromatic", 35)
    v30, _ = yields.mass_yield("Aromatic", 30)
    assert v35 == v30 and "C30" in flags[0]


def test_oxygen_column_caps_at_five():
    assert yields.oxygen_column(1) == "O1"
    assert yields.oxygen_column(7) == "O5"


# --- structure route ------------------------------------------------------

@pytest.mark.parametrize("smiles, column", [
    ("CCCCCCCCCC", "n-Alkane"),
    ("CC(C)CC(C)(C)C", "b-Alkane"),
    ("CC1CCCCC1", "Cycloalkane"),
    ("[H][C@@]12CCCC[C@]1([H])CCCC2", "Cycloalkane"),   # cis-decalin, stereo SMILES
    ("CC[C@@H]1CCC[C@@H]1C", "Cycloalkane"),
    ("c1ccccc1", "Aromatic"),                          # aromatic (lower-case) SMILES
    ("CC1=CC=CC=C1", "Aromatic"),                      # Kekule SMILES
    ("C1CC2=CC=CC=C2C1", "Aromatic"),                  # indane: one aromatic ring
    ("C1=CC2=CC=CC=C2C=C1", "PAH"),
    ("C1=CC=C(C=C1)C2=CC=CC=C2", "PAH"),               # biphenyl
    ("CCCCCCCCC=C", "Alkenes"),
    ("C1CCC=CC1", "Cycloalkene"),
    ("CC1=CCC2CC1C2(C)C", "Cycloalkene"),              # alpha-pinene (terpene rule)
    ("CC(=C)C1CCC(C)=CC1", "Cycloalkene"),             # limonene
    ("CC(C)=CCC=C(C)C=C", "Cycloalkene"),              # ocimene: acyclic terpene, terpene rule
    ("CCCCCCC(C)=O", "O1"),
    ("OCC(O)CO", "O3"),
])
def test_structure_classes(smiles, column):
    assert soafp.soafp_from_structure(smiles)["column"] == column


def test_structure_value_and_bounds():
    r = soafp.soafp_from_structure("CC1=CC=CC=C1")      # toluene, C7
    assert r["n_C"] == 7
    assert r["soafp_low"] <= r["soafp_avg"] <= r["soafp_high"]
    assert r["soafp_avg"] == pytest.approx(0.0895767, rel=1e-4)


def test_heteroatoms():
    assert soafp.soafp_from_structure("Clc1ccccc1")["soafp_avg"] == 0.0
    d5 = soafp.soafp_from_structure(
        "C[Si]1(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O1")
    assert d5["column"] == "O5" and any("Si" in f for f in d5["flags"])


def test_formula_mismatch_is_flagged():
    r = soafp.soafp_from_structure("CCCCCCCCC", formula="C9H18")
    assert any("does not match" in f for f in r["flags"])


def test_carbon_number_override_for_lumped_species():
    r = soafp.soafp_from_structure("CCCCCCCCC1=CC=CC=C1", n_C=18)
    assert r["column"] == "Aromatic" and r["n_C"] == 18
    assert any("given value used" in f for f in r["flags"])
    r = soafp.soafp_from_structure("CCCCCCC(C)=O", n_O=3)
    assert r["column"] == "O3"


def test_bad_smiles_raises():
    with pytest.raises(ValueError):
        soafp.soafp_from_structure("---")


# --- formula route --------------------------------------------------------

def test_formula_ranges():
    r = soafp.soafp_range_from_formula("C10H22")
    assert r["column"] == "n-Alkane; b-Alkane"
    assert r["soafp_avg"] is None
    assert r["soafp_low"] <= r["soafp_avg_low"] <= r["soafp_avg_high"] <= r["soafp_high"]
    assert soafp.soafp_range_from_formula("C9H12")["column"] == "Aromatic"
    assert soafp.soafp_range_from_formula("C10H8")["column"] == "PAH"


def test_structure_value_lies_in_formula_range():
    for smiles in ("CCCCCCCCCC", "CC(C)CCCCCCC", "CC1CCCCCCCC1", "CCCCCCCCC=C",
                   "CC1=CC=CC=C1", "CC1=CCC2CC1C2(C)C", "CCCCCCCC(C)=O"):
        s = soafp.soafp_from_structure(smiles)
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        f = soafp.soafp_range_from_formula(
            rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smiles)))
        assert f["soafp_avg_low"] - 1e-12 <= s["soafp_avg"] <= f["soafp_avg_high"] + 1e-12


def test_terpene_formulas():
    assert soafp.soafp_range_from_formula("C10H16")["column"] == "Cycloalkene"
    # C15H24 is a sesquiterpene or an alkylbenzene (e.g. tripropylbenzene)
    assert soafp.soafp_range_from_formula("C15H24")["column"] == "Cycloalkene; Aromatic"


def test_formula_heteroatom_and_small():
    assert soafp.soafp_range_from_formula("C3H6O")["soafp_high"] == 0.0     # < C5
    r = soafp.soafp_range_from_formula("C8H16O")                            # DBE 1
    assert r["column"] == "O1" and r["n_candidates"] == 1
    # halocarbon: from 0 to the H-substituted analogue (benzene)
    r = soafp.soafp_range_from_formula("C6H5Cl")
    assert r["column"] == "none; Aromatic" and r["soafp_low"] == 0.0
    assert r["soafp_high"] == soafp.soafp_range_from_formula("C6H6")["soafp_high"]
    # N / S / Si: from 0 to an oxygenate with the heteroatoms as O-groups
    assert soafp.soafp_range_from_formula("C10H21NO")["column"] == "none; O2"
    assert soafp.soafp_range_from_formula("C10H30O5Si5")["column"] == "none; O5"


def test_oxygenated_aromatic_formula():
    r = soafp.soafp_range_from_formula("C7H8O")       # cresols, benzyl alcohol, anisole
    assert r["column"] == "O1; Aromatic"
    o1 = soafp.mass_yield("O1", 7)[0]["avg"]
    ar = soafp.mass_yield("Aromatic", 7)[0]["avg"]
    assert r["soafp_central"] == pytest.approx((o1 + ar) / 2)
    assert soafp.soafp_range_from_formula("C7H14O")["column"] == "O1"   # DBE 1: aliphatic


def test_central_and_uncertainty():
    r = soafp.soafp_range_from_formula("C10H22")
    a = soafp.mass_yield("n-Alkane", 10)[0]["avg"]
    b = soafp.mass_yield("b-Alkane", 10)[0]["avg"]
    assert r["soafp_central"] == pytest.approx((a + b) / 2)
    assert r["soafp_low"] <= r["soafp_central"] <= r["soafp_high"]
    assert r["soafp_minus"] == pytest.approx(r["soafp_central"] - r["soafp_low"])
    assert r["soafp_plus"] == pytest.approx(r["soafp_high"] - r["soafp_central"])
    s = soafp.soafp_from_structure("CC1=CC=CC=C1")
    assert s["soafp_central"] == s["soafp_avg"] and s["n_candidates"] == 1


def test_weights():
    w = soafp.soafp_range_from_formula("C10H22", weights={"b-Alkane": 3})
    b = soafp.mass_yield("b-Alkane", 10)[0]["avg"]
    a = soafp.mass_yield("n-Alkane", 10)[0]["avg"]
    assert w["soafp_central"] == pytest.approx((a + 3 * b) / 4)
    assert any("weighted" in f for f in w["flags"])


def test_soafp_from_formulas():
    pytest.importorskip("pandas")
    out = soafp.soafp_from_formulas(["C10H16", "C7H8O", "not a formula", None])
    assert list(out["formula"])[:2] == ["C10H16", "C7H8O"]
    assert out.loc[0, "soafp_central"] > 0 and out.loc[1, "n_candidates"] == 2
    assert out.loc[2, "soafp_method"] == "none" and out["soafp_central"].isna()[2]
    assert "not calculated" in out.loc[3, "soafp_flags"]


# --- table helper ---------------------------------------------------------

def test_soafp_table_routes():
    pd = pytest.importorskip("pandas")
    df = pd.DataFrame({"smiles": ["CC1=CC=CC=C1", "---", None, "CCCCCCCCCCCC"],
                       "formula": ["C7H8", "C10H16", None, "C12H26"],
                       "cls": [None, None, None, "Cycloalkane"],
                       "nC": [None, None, None, 20]})
    out = soafp.soafp_table(df, smiles="smiles", formula="formula", column="cls", n_C="nC")
    assert list(out["soafp_method"]) == ["structure", "formula", "none", "class"]
    assert out.loc[3, "n_C"] == 20
    assert "cannot parse" in out.loc[1, "soafp_flags"]


# --- calc_soafp script -----------------------------------------------------

def test_calc_soafp_input_detection():
    import calc_soafp as cs
    assert cs.resolve("C10H16")[0] == "formula"
    assert cs.resolve("CCCCCC")[0] == "smiles"          # no H: not read as a formula
    assert cs.resolve("CO")[0] == "smiles"              # methanol
    assert cs.resolve("CS2", kind="formula")[0] == "formula"
    assert cs.resolve("C10H17+")[0] == "unknown"        # ion formula
    assert cs.resolve("toluene")[0] == "unknown"        # names need online=True
    r = cs.soafp_of("O=Cc1ccccc1")                      # benzaldehyde: same as its formula
    assert r["column"] == "O1; Aromatic"
    assert r["soafp_central"] == pytest.approx(
        soafp.soafp_range_from_formula("C7H6O")["soafp_central"])
    out = cs.soafp_of_list(["C10H16", "nonsense!"])
    assert out["soafp_central"].notna().tolist() == [True, False]
