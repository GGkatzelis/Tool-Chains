# -*- coding: utf-8 -*-
"""SOAFP for a whole species table (needs pandas)."""
from .core import soafp_from_class, soafp_from_structure, soafp_range_from_formula
from . import classify

OUTPUT_COLUMNS = ["soafp_method", "soafp_column", "n_C", "n_O", "dbe",
                  "soafp_noVWL", "soafp_VWL", "soafp_avg",
                  "soafp_avg_low", "soafp_avg_high", "soafp_low", "soafp_high",
                  "soafp_flags"]


def _present(value):
    return value is not None and value == value and str(value).strip() not in ("", "n/a")


def soafp_one(smiles=None, formula=None, column=None, route="structure",
              n_C=None, n_O=None):
    """
    One species, choosing the route:

    * `column` given (a SOM column name): that class, with n_C / n_O from the
      SMILES or, failing that, the formula;
    * route="structure" (default): the SMILES if it parses, else the formula range;
    * route="formula": the formula range even when a SMILES is available.

    n_C / n_O, when given, override the carbon / oxygen numbers of the
    structure or the class route (lumped species); the formula route ignores them.
    """
    smiles = str(smiles).strip() if _present(smiles) else None
    formula = str(formula).strip() if _present(formula) else None
    n_C = _number(n_C)
    n_O = _number(n_O)
    note = []
    if column is not None and _present(column):
        c = None
        if smiles:
            try:
                _, c, _ = classify.column_from_smiles(smiles)
            except ValueError:
                pass
        if c is None and formula and classify.ce.is_formula(formula):
            c = classify.ce.counts(formula)
        if c is None and n_C is None:
            return _empty("class", ["no SMILES, formula or n_C to get n_C from"])
        c = c or {}
        return soafp_from_class(column, n_C if n_C is not None else c.get("C", 0),
                                n_O if n_O is not None else c.get("O", 0))
    if smiles and route == "structure":
        try:
            return soafp_from_structure(smiles, formula, n_C, n_O)
        except ValueError as err:
            note.append(f"{err}; formula used instead")
    if formula and classify.ce.is_formula(formula):
        r = soafp_range_from_formula(formula)
        r["flags"] = note + r["flags"]
        return r
    return _empty("none", note + ["no usable SMILES or formula: SOAFP not calculated"])


def _number(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None                 # None, NaN-like strings, 'n/a'
    return None if value != value else int(value)


def _empty(method, flags):
    r = dict.fromkeys(OUTPUT_COLUMNS)
    r.update({"method": method, "column": None, "flags": flags})
    return r


def soafp_table(df, smiles=None, formula=None, column=None, route="structure",
                n_C=None, n_O=None):
    """
    A copy of DataFrame `df` with the SOAFP output columns added.
    `smiles`, `formula`, `column`, `n_C`, `n_O` are the names of the columns of
    `df` to use (any may be None); see soafp_one for what each does.
    """
    rows = []
    for _, rec in df.iterrows():
        r = soafp_one(rec[smiles] if smiles else None,
                      rec[formula] if formula else None,
                      rec[column] if column else None, route,
                      rec[n_C] if n_C else None, rec[n_O] if n_O else None)
        rows.append({"soafp_method": r["method"], "soafp_column": r["column"],
                     **{k: r.get(k) for k in OUTPUT_COLUMNS[2:-1]},
                     "soafp_flags": "; ".join(r["flags"])})
    out = df.copy()
    for col in OUTPUT_COLUMNS:
        out[col] = [row[col] for row in rows]
    return out
