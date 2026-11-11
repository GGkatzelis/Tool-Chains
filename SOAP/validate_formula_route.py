# -*- coding: utf-8 -*-
"""
How good is the formula-only SOAFP? Test it on species whose structure is known.

For every species in SOAP/profiles with a readable SMILES, the formula is
taken from that SMILES (so placeholder formulas in the tables do not
matter) and two SOAFP values are compared:

    reference   structure route (one SOM class, avg of no-VWL and VWL)
    formula     formula route (DBE candidates): central, low, high

Reported per profile and for all profiles together:
    coverage      share of species whose reference lies in [low, high]
    within x2     share with formula central within a factor 2 of the reference
    median ratio  median of central / reference (species with reference > 0)
    weighted      sum(wt % x SOAFP) / sum(wt %) for reference and formula central

Writes results/formula_validation.csv (one row per species).
Usage:  python validate_formula_route.py
"""
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import soafp  # noqa: E402

from rdkit import Chem  # noqa: E402
from rdkit.Chem import rdMolDescriptors  # noqa: E402

#  profile:  (SMILES column, name column, weight column)
PROFILES = {
    "gasoline": ("Speciate SMILES", "Compound", "wt %"),
    "diesel":   ("Speciate SMILES", "Compound", "wt %"),
    "vcp":      ("Speciate SMILES", "Compound", "wt %.6"),   # total VCP column
    "cooking":  ("SMILES", "Assigned Species", "Mass Fraction Determined by PMF (wt %)"),
}


def species(profile):
    smiles_col, name_col, wt_col = PROFILES[profile]
    df = pd.read_csv(HERE / "profiles" / f"{profile}.csv")
    rows = []
    for _, rec in df.iterrows():
        mol = Chem.MolFromSmiles(str(rec[smiles_col])) if isinstance(rec[smiles_col], str) else None
        if mol is None:
            continue
        formula = rdMolDescriptors.CalcMolFormula(mol).rstrip("+-")
        s = soafp.soafp_from_structure(rec[smiles_col])
        f = soafp.soafp_range_from_formula(formula)
        rows.append({"profile": profile, "name": rec[name_col], "smiles": rec[smiles_col],
                     "formula": formula, "wt": pd.to_numeric(rec[wt_col], errors="coerce"),
                     "class": s["column"], "reference": s["soafp_central"],
                     "formula_column": f["column"], "central": f["soafp_central"],
                     "low": f["soafp_low"], "high": f["soafp_high"],
                     "n_candidates": f["n_candidates"]})
    return pd.DataFrame(rows)


def summary(d, label):
    eps = 1e-12
    covered = (d["reference"] >= d["low"] - eps) & (d["reference"] <= d["high"] + eps)
    pos = d[d["reference"] > 0]
    ratio = pos["central"] / pos["reference"]
    w = d["wt"].fillna(0)
    wsum = w.sum() or float("nan")
    return {"set": label, "species": len(d),
            "coverage": covered.mean(),
            "within_x2": ((ratio >= 0.5) & (ratio <= 2)).mean() if len(pos) else None,
            "median_ratio": ratio.median() if len(pos) else None,
            "weighted_reference": (w * d["reference"]).sum() / wsum,
            "weighted_formula": (w * d["central"]).sum() / wsum}


if __name__ == "__main__":
    data = pd.concat([species(p) for p in PROFILES], ignore_index=True)
    (HERE / "results").mkdir(exist_ok=True)
    data.to_csv(HERE / "results" / "formula_validation.csv", index=False)
    table = [summary(data[data["profile"] == p], p) for p in PROFILES]
    table.append(summary(data, "all"))
    for cls, d in data.groupby("class", dropna=False):
        table.append(summary(d, f"class {cls}"))
    pd.set_option("display.width", 140)
    print(pd.DataFrame(table).round(4).to_string(index=False))
