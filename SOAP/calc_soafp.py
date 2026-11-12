# -*- coding: utf-8 -*-
"""
calc_soafp.py: SOAFP (SOM SOA mass yield, g SOA / g reacted) with uncertainty
for any list of species.

WHAT YOU PROVIDE (per species, one of)
  formula   neutral molecular formula, e.g. C10H16, C7H8O, C6H5Cl, CH3CH2OH.
            Elements C H O N S P Si F Cl Br I; no charges, no isotopes, no
            "+H" ion formulas (give the neutral molecule, e.g. C10H16, not
            C10H17+). -> central value and range over the classes the DBE allows.
  SMILES    e.g. CC1=CC=CC=C1. Needs RDKit (pip install rdkit).
            -> one SOM class, range = vapor-wall-loss range only.
  name      e.g. limonene. Only with --online: the name is looked up on
            PubChem (only the name is sent) and the SMILES found is used.
            Check the "smiles" column: a trade or mixture name may resolve to
            the wrong compound.

Auto-detection (--kind auto, default): a string that is a valid formula AND
contains H is a formula; else it is tried as SMILES; else as a name. Force
the type with --kind formula|smiles|name when this guesses wrong (e.g. "CO"
is read as the SMILES of methanol, "CS2" needs --kind formula).

LIMITS (every row also lists its assumptions in soafp_flags)
  * SOM covers C5-C30: fewer than 5 C -> 0; more than 30 C -> C30 value.
  * Oxygenates use alkane-based O-columns (O1..O5); aromatic oxygenates
    range from that column up to the Aromatic column.
  * N, S, P, Si, halogens have no SOM fit: their range starts at 0.
  * Lumped mixtures (UCM, "C18 aromatics", plant extracts) are not single
    species: give a representative formula with the right carbon number.
  * Values are per mass reacted at the SOM run conditions, not per mass emitted.

USAGE
  python calc_soafp.py C10H16 C7H8O "CC1=CC=CC=C1" limonene
  python calc_soafp.py --file peaks.csv --column Formula            # CSV or .xlsx
  python calc_soafp.py --file list.xlsx --column Name --kind name --online --out out.csv

  from calc_soafp import soafp_of, soafp_of_list
  soafp_of("C10H16")              # one species -> dict
  soafp_of_list(["C10H16", "toluene"])   # many -> pandas DataFrame
"""
import argparse
import json
import sys
import urllib.parse
import urllib.request
from functools import lru_cache
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import soafp  # noqa: E402
from soafp.classify import ce  # noqa: E402   (ClassEstimationDBE)

KINDS = ("auto", "formula", "smiles", "name")
SHOW = ["input", "interpreted_as", "smiles", "formula", "soafp_central", "soafp_low",
        "soafp_high", "soafp_minus", "soafp_plus", "n_candidates", "column",
        "n_C", "n_O", "dbe", "flags"]


def _is_smiles(text):
    try:
        from rdkit import Chem, RDLogger
    except ImportError:
        return False
    RDLogger.DisableLog("rdApp.*")
    return Chem.MolFromSmiles(text) is not None


@lru_cache(maxsize=None)
def _pubchem(name):
    """SMILES of a compound name from PubChem, or None."""
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           + urllib.parse.quote(name) + "/property/SMILES/JSON")
    req = urllib.request.Request(url, headers={"User-Agent": "soafp-tool-chain"})
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            return json.load(r)["PropertyTable"]["Properties"][0]["SMILES"]
    except Exception:
        return None


def resolve(text, kind="auto", online=False):
    """(kind found, SMILES or None, formula or None, notes) for one input string."""
    text = str(text).strip()
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}")
    if text[-1:] in "+-" and ce.is_formula(text[:-1]):
        return "unknown", None, None, ["ion formula: give the neutral molecule "
                                       "(e.g. C10H16 for the C10H17+ ion)"]
    if kind == "formula" or (kind == "auto" and ce.is_formula(text) and "H" in ce.counts(text)):
        return "formula", None, text, []
    if kind in ("auto", "smiles") and _is_smiles(text):
        return "smiles", text, None, []
    if kind in ("auto", "name"):
        if online:
            smi = _pubchem(text)
            if smi:
                return "name", smi, None, [f"name found on PubChem: {smi}"]
        return "unknown", None, None, ["not a formula or SMILES; name not found on PubChem"
                                       if online else "not a formula or SMILES "
                                       "(for a compound name use --online)"]
    return "unknown", None, None, [f"not a valid {kind}"]


def soafp_of(text, kind="auto", online=False, weights=None):
    """SOAFP of one species given as formula, SMILES or name -> dict (see SHOW)."""
    found, smi, formula, notes = resolve(text, kind, online)
    out = dict.fromkeys(SHOW)
    out.update(input=text, interpreted_as=found, smiles=smi, formula=formula)
    try:
        if found == "formula":
            r = soafp.soafp_range_from_formula(formula, weights)
        elif smi:
            r = soafp.soafp_from_structure(smi)
            out["formula"] = soafp.core._formula(soafp.classify.column_from_smiles(smi)[1])
        else:
            out["flags"] = "; ".join(notes + ["SOAFP not calculated"])
            return out
    except (ValueError, ImportError) as err:
        out["flags"] = "; ".join(notes + [str(err), "SOAFP not calculated"])
        return out
    out.update({k: r.get(k) for k in SHOW if k in r})
    out["flags"] = "; ".join(notes + r["flags"])
    return out


def soafp_of_list(items, kind="auto", online=False, weights=None):
    """SOAFP for many species -> pandas DataFrame, one row per input."""
    import pandas as pd
    out = pd.DataFrame([soafp_of(x, kind, online, weights) for x in items], columns=SHOW)
    out["n_candidates"] = out["n_candidates"].astype("Int64")
    return out


def _read_column(path, column):
    import pandas as pd
    path = Path(path)
    df = pd.read_excel(path) if path.suffix.lower() in (".xlsx", ".xls") else pd.read_csv(path)
    if column not in df.columns:
        raise SystemExit(f"column {column!r} not in {path.name}; columns: {list(df.columns)}")
    return df, [x for x in df[column]]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog="Run with no arguments to see this help.")
    p.add_argument("species", nargs="*", help="formulas, SMILES or names")
    p.add_argument("--file", help="CSV or Excel file with one species per row")
    p.add_argument("--column", help="column of --file holding the species")
    p.add_argument("--kind", default="auto", choices=KINDS)
    p.add_argument("--online", action="store_true", help="look up unknown names on PubChem")
    p.add_argument("--out", help="write the result table to this CSV")
    a = p.parse_args(argv)
    if not a.species and not a.file:
        p.print_help()
        print(__doc__)
        return
    import pandas as pd
    if a.file:
        if not a.column:
            raise SystemExit("--file needs --column")
        df, items = _read_column(a.file, a.column)
        rows = [soafp_of(x, a.kind, a.online) if pd.notna(x) and str(x).strip()
                else {**dict.fromkeys(SHOW), "flags": "empty cell"} for x in items]
        res = pd.DataFrame(rows, columns=SHOW, index=df.index).drop(columns="input")
        res.columns = [c if c.startswith("soafp_") else "soafp_" + c for c in res.columns]
        res = pd.concat([df, res], axis=1)
    else:
        res = soafp_of_list(a.species, a.kind, a.online)
    if a.out:
        res.to_csv(a.out, index=False)
        print(f"{len(res)} rows -> {a.out}")
    else:
        cols = ["interpreted_as", "formula", "soafp_central", "soafp_low",
                "soafp_high", "n_candidates", "column"]
        cols = [a.column or "input"] + [c if c in res else "soafp_" + c for c in cols]
        print(res[cols].round(4).to_string(index=False))
        flags = "flags" if "flags" in res else "soafp_flags"
        for _, row in res[res["soafp_central"].isna()].iterrows():
            print(f"  no value for {row[cols[0]]!r}: {row[flags]}")
        print("Assumptions for every row are in the flags column (use --out to save it).")


if __name__ == "__main__":
    main()
