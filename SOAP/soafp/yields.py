# -*- coding: utf-8 -*-
"""
SOM mass-yield tables (plain CSV in soafp/data, no spreadsheet needed).

Two cases are tabulated, both as SOA mass yield (g SOA / g precursor reacted):

    noVWL   SOM parameters fit without vapor-wall-loss correction (lower bound)
    VWL     SOM parameters fit with vapor-wall-loss correction (upper bound)

"avg" is their arithmetic mean, the value used by the original SOAFP scripts.
"""
import csv
from functools import lru_cache
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data"
CASES = ("noVWL", "VWL")
ALL_CASES = CASES + ("avg",)

HYDROCARBON_COLUMNS = ("n-Alkane", "b-Alkane", "Cycloalkane", "Aromatic", "PAH",
                       "Alkenes", "Cycloalkene")
OXYGEN_COLUMNS = ("O1", "O2", "O3", "O4", "O5")
COLUMNS = HYDROCARBON_COLUMNS + OXYGEN_COLUMNS

MIN_C, MAX_C = 5, 30      # carbon range of the tables


@lru_cache(maxsize=None)
def table(case):
    """{column: {n_C: yield}} for case 'noVWL' or 'VWL'."""
    if case not in CASES:
        raise ValueError(f"case must be one of {CASES}, not {case!r}")
    path = DATA_DIR / f"som_mass_yield_{case}.csv"
    with open(path, encoding="utf-8") as f:
        rows = csv.DictReader(line for line in f if not line.startswith("#"))
        out = {col: {} for col in COLUMNS}
        for row in rows:
            n_C = int(row["n_C"])
            for col in COLUMNS:
                out[col][n_C] = float(row[col])
    return out


def oxygen_column(n_O):
    """'O1' ... 'O5' for n_O >= 1 (five or more O share 'O5')."""
    if n_O < 1:
        raise ValueError("oxygen columns need n_O >= 1")
    return f"O{min(int(n_O), 5)}"


def mass_yield(column, n_C):
    """
    SOM mass yield of one table column at carbon number n_C, for every case.

    Returns ({'noVWL': y, 'VWL': y, 'avg': y}, flags). Outside the table:
    n_C < 5 gives 0 (as in the original method); n_C > 30 uses the C30 value
    (the yields plateau above ~C23). Both are reported in flags.
    """
    if column not in COLUMNS:
        raise ValueError(f"unknown column {column!r}; expected one of {COLUMNS}")
    flags = []
    n_C = int(n_C)
    if n_C < MIN_C:
        flags.append(f"n_C={n_C} < {MIN_C}: below the SOM table, SOAFP set to 0")
        return {case: 0.0 for case in ALL_CASES}, flags
    if n_C > MAX_C:
        flags.append(f"n_C={n_C} > {MAX_C}: C{MAX_C} value used")
        n_C = MAX_C
    values = {case: table(case)[column][n_C] for case in CASES}
    values["avg"] = (values["noVWL"] + values["VWL"]) / 2
    return values, flags
