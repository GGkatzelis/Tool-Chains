# -*- coding: utf-8 -*-
"""
SOAFP for the four speciation profiles in SOAP/profiles (CSV), written to
SOAP/results/<profile>_soafp.csv.

Each species gets
  * the structure-based value (SMILES via RDKit; formula range if no SMILES),
  * the formula-only range in the formula_* columns, for comparison with what
    a mass-spectrometer formula alone would give.

Lumped species whose SMILES/formula are placeholders (e.g. "C20 bicycloalkanes"
entered with an n-alkane SMILES) are corrected through NAME_HINTS: a match on
the compound name sets the class by hand and is flagged. This is a property of
these input tables, not of the method.

Usage:  python run_profiles.py [gasoline diesel vcp cooking]
"""
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import soafp  # noqa: E402

PROFILES = {
    #  name       SMILES column      formula column     name column         n_C / n_O columns  use them?
    "gasoline": ("Speciate SMILES", "Species Formula", "Compound",         None, None,        None),
    "diesel":   ("Speciate SMILES", "Species Formula", "Compound",         "n_C", "n_O",      "override"),
    "vcp":      ("Speciate SMILES", "Species Formula", "Compound",         "No. C", "No. O",  "check"),
    "cooking":  ("SMILES",          "Formula",         "Assigned Species", None, None,        None),
}
# The diesel and VCP tables give carbon / oxygen numbers per row.
#  diesel  "override": most mismatches are lumped species whose SMILES is only
#          a representative ("Aromatic" C15...C23 all with a C14 SMILES, "UCM"
#          n_C = 26 with dodecane), so the table's numbers are used.
#  vcp     "check": mismatches are mostly errors in the numbers for named
#          compounds (3-methylpentane given C4, benzaldehyde 0 O), so the
#          structure is used and the mismatch is only flagged.

# compound-name pattern -> SOM column (hydrocarbon lumps only)
NAME_HINTS = [
    (r"branched alkane", "b-Alkane"),
    (r"cycloalkane", "Cycloalkane"),       # also bi-/tricycloalkanes
    (r"cycloalkene", "Cycloalkene"),
]


# flags that point at a problem in the input table rather than an assumption
ISSUE_PATTERN = ("does not match|cannot parse|not calculated|class set by user|"
                 "given value used|structure kept|differs, not used")


def name_hint(name):
    for pattern, column in NAME_HINTS:
        if re.search(pattern, str(name), flags=re.IGNORECASE):
            return column
    return None


def run(profile):
    smiles_col, formula_col, name_col, nC_col, nO_col, use = PROFILES[profile]
    df = pd.read_csv(HERE / "profiles" / f"{profile}.csv")
    df["soafp_class_from_name"] = df[name_col].map(name_hint)
    override = use == "override"
    out = soafp.soafp_table(df, smiles=smiles_col, formula=formula_col,
                            column="soafp_class_from_name",
                            n_C=nC_col if override else None,
                            n_O=nO_col if override else None)
    if use == "check":
        for given, ours in ((nC_col, "n_C"), (nO_col, "n_O")):
            g = pd.to_numeric(df[given], errors="coerce")
            bad = g.notna() & out[ours].notna() & (g != out[ours])
            out.loc[bad, "soafp_flags"] = [
                (f + "; " if f else "") + f"table '{given}' = {int(v)} differs, not used"
                for f, v in zip(out.loc[bad, "soafp_flags"], g[bad])]
    by_formula = soafp.soafp_table(df[[formula_col]], formula=formula_col, route="formula")
    for col in ("soafp_column", "soafp_central", "soafp_low", "soafp_high",
                "n_candidates", "soafp_flags"):
        out["formula_" + col.replace("soafp_", "", 1)] = by_formula[col].values
    (HERE / "results").mkdir(exist_ok=True)
    path = HERE / "results" / f"{profile}_soafp.csv"
    out.to_csv(path, index=False)
    issues = out["soafp_flags"].str.contains(ISSUE_PATTERN, regex=True, na=False)
    out.loc[issues, [name_col, formula_col, smiles_col, "n_C", "n_O", "soafp_column",
                     "soafp_flags"]].to_csv(HERE / "results" / f"{profile}_input_issues.csv",
                                            index=False)
    return out, path


if __name__ == "__main__":
    for name in sys.argv[1:] or PROFILES:
        table, path = run(name)
        n_zero = int((table["soafp_avg_high"].fillna(0) == 0).sum())
        n_flag = int((table["soafp_flags"] != "").sum())
        print(f"{name:9s} {len(table):4d} species, {n_zero:3d} with SOAFP 0, "
              f"{n_flag:3d} flagged -> {path.name}")
