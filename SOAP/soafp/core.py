# -*- coding: utf-8 -*-
"""
SOAFP (SOA formation potential = SOM SOA mass yield) for single species.

Every function returns a dict with the same keys, so results from the
structure, class-override and formula routes stack into one table:

    method          'structure' | 'class' | 'formula'
    column          SOM column(s) used, '; '-joined when the formula is ambiguous
    n_C, n_O, dbe   composition used for the lookup
    soafp_avg       mean of the no-VWL and VWL yields (single-column routes only)
    soafp_avg_low / soafp_avg_high   min / max of soafp_avg over the candidates
    soafp_low / soafp_high           min / max over candidates and both VWL cases
    flags           list of every assumption or problem met on the way
"""
from . import classify
from .yields import CASES, OXYGEN_COLUMNS, mass_yield, oxygen_column

import class_estimation as ce   # made importable by classify


def _result(method, candidates, c, dbe, flags):
    """Look up every (label, column) candidate and build the common result dict."""
    n_C, n_O = c.get("C", 0), c.get("O", 0)
    out = {"method": method, "column": None, "n_C": n_C, "n_O": n_O, "dbe": dbe,
           "soafp_noVWL": 0.0, "soafp_VWL": 0.0, "soafp_avg": 0.0,
           "soafp_avg_low": 0.0, "soafp_avg_high": 0.0,
           "soafp_low": 0.0, "soafp_high": 0.0, "candidates": {}, "flags": list(flags)}
    if not candidates:
        return out
    for label, column in candidates:
        values, more = mass_yield(column, n_C)
        out["candidates"][label] = values
        out["flags"] += [f for f in more if f not in out["flags"]]
    vals = list(out["candidates"].values())
    out["column"] = "; ".join(column for _, column in candidates)
    out["soafp_avg_low"] = min(v["avg"] for v in vals)
    out["soafp_avg_high"] = max(v["avg"] for v in vals)
    out["soafp_low"] = min(v[case] for v in vals for case in CASES)
    out["soafp_high"] = max(v[case] for v in vals for case in CASES)
    if len(vals) == 1:
        out["soafp_noVWL"], out["soafp_VWL"], out["soafp_avg"] = (
            vals[0]["noVWL"], vals[0]["VWL"], vals[0]["avg"])
    else:
        out["soafp_noVWL"] = out["soafp_VWL"] = out["soafp_avg"] = None
    return out


def soafp_from_structure(smiles, formula=None, n_C=None, n_O=None):
    """
    SOAFP from a SMILES string: exactly one SOM column. If a formula is given
    too, it is only checked against the structure (a mismatch is flagged and
    the structure wins).

    n_C / n_O override the carbon / oxygen numbers of the structure, for
    lumped species whose SMILES is only a representative of the class
    (e.g. "C18 aromatics" entered with a C14 SMILES). The class still comes
    from the structure; overrides are flagged.
    """
    column, c, flags = classify.column_from_smiles(smiles)
    if formula is not None and ce.is_formula(str(formula).strip()):
        given = ce.counts(str(formula).strip())
        if given != c:
            flags.append(f"formula {formula} does not match SMILES ({_formula(c)}); "
                         "structure used")
    c, column = _override(c, column, n_C, n_O, flags)
    candidates = [(column, column)] if column else []
    return _result("structure", candidates, c, ce.dbe(c), flags)


def _override(c, column, n_C, n_O, flags):
    """Apply given carbon / oxygen numbers (flagged) to a composition and its column."""
    c = dict(c)
    for element, value in (("C", n_C), ("O", n_O)):
        if value is None or value != value:      # None or NaN
            continue
        value = int(value)
        if element == "O" and value == 0 and column in OXYGEN_COLUMNS:
            flags.append(f"n_O=0 given, but the structure has {c['O']} O: structure kept")
            continue
        if value != c.get(element, 0):
            flags.append(f"n_{element}={value} given, structure has {c.get(element, 0)}: "
                         "given value used")
            c[element] = value
    if column in OXYGEN_COLUMNS and c.get("O", 0) >= 1:
        column = oxygen_column(c["O"])
    return c, column


def soafp_from_class(column, n_C, n_O=0):
    """SOAFP for a class assigned by hand (e.g. lumped species such as 'C20 bicycloalkanes')."""
    c = {"C": int(n_C), "O": int(n_O)}
    return _result("class", [(column, column)], c, None,
                   [f"class set by user: {column}"])


def soafp_range_from_formula(formula):
    """
    SOAFP range from a neutral molecular formula alone: every SOM column the
    formula allows, with soafp_low/high spanning candidates and VWL cases.
    """
    candidates, c, dbe, flags = classify.columns_from_formula(str(formula).strip())
    return _result("formula", candidates, c, dbe, flags)


def _formula(c):
    order = ["C", "H"] + sorted(e for e in c if e not in ("C", "H"))
    return "".join(f"{e}{c[e] if c[e] != 1 else ''}" for e in order if c.get(e))
