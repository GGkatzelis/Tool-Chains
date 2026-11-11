# -*- coding: utf-8 -*-
"""
Compound class from a neutral molecular formula alone: element counts and the
double bond equivalent (DBE), no library lookup and no assumption about how
the formula was obtained.

Our classes, every rule reading the formula:

    nitrogen    any N                                   (exclusive)
    terpene     C/H only, C10H16, C15H24 or C20H32      (exclusive, except
                aromatic for C15H24 / C20H32)
    aromatic    >= 6 C, DBE >= 4; a C6 with DBE = 4 is aromatic only
    alkene      C/H only, >= 2 C, 1 <= DBE <= 3
    ketone      >= 3 C, 1 O, DBE >= 1
    aldehyde    >= 1 C, 1 O, DBE >= 1
    ester       >= 2 C, >= 2 O, DBE >= 1
    acid        >= 1 C, >= 2 O, DBE >= 1
    alcohol     1-2 O, DBE = 0
    other       nothing above fits

The rules rule classes out; they do not pick one. Several can remain and all
of them are returned - that ambiguity is information for the caller, not
something to resolve by guessing.

reference_classes() translates a formula into the hydrocarbon classes of the
SOA-yield reference (n-Alkane, b-Alkane, Cycloalkane, Aromatic, PAH, Alkenes,
Cycloalkene), again listing every class the formula allows. Turning several
candidates into a lower and upper yield is left to the caller; see README.md.

@author: exp.voc
"""
import re

ELEMENTS = ("H", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Br", "I")
HALOGENS = ("F", "Cl", "Br", "I")

CLASSES = ["ketone", "aldehyde", "ester", "alcohol", "acid", "aromatic",
           "alkene", "terpene", "nitrogen"]
MIN_C = {"ketone": 3, "aldehyde": 1, "ester": 2, "acid": 1, "alcohol": 1,
         "aromatic": 6, "alkene": 2}
MIN_O = {"ketone": 1, "aldehyde": 1, "ester": 2, "acid": 2, "alcohol": 1}
MAX_O = {"ketone": 1, "aldehyde": 1, "alcohol": 2}
MIN_DBE = {"ketone": 1, "aldehyde": 1, "ester": 1, "acid": 1, "aromatic": 4}

REFERENCE = ["n-Alkane", "b-Alkane", "Cycloalkane", "Aromatic", "PAH",
             "Alkenes", "Cycloalkene"]
PAH_MIN_DBE, PAH_MIN_C = 7, 10

_ELEMENT = re.compile(r"([A-Z][a-z]?)(\d*)")
_FORMULA = re.compile(r"(?:[A-Z][a-z]?\d*)+")


# --- formulas -----------------------------------------------------------

def is_formula(formula):
    """True for a plain formula of known elements, e.g. 'C7H8' or 'C2H6O'."""
    if not isinstance(formula, str):
        return False
    f = formula.strip()
    return bool(_FORMULA.fullmatch(f)) and all(e in ELEMENTS for e, _ in _ELEMENT.findall(f))


def counts(formula):
    """'C2H6O' -> {'C': 2, 'H': 6, 'O': 1}; repeated elements are summed."""
    out = {}
    for element, n in _ELEMENT.findall(formula.strip()):
        out[element] = out.get(element, 0) + int(n or 1)
    return {e: n for e, n in out.items() if n}


def dbe(c):
    """
    Double bond equivalent of a composition, None without carbon:
    C - (H + halogens)/2 + N/2 + 1; O, S, Si and P do not enter.
    """
    if c.get("C", 0) == 0:
        return None
    halogens = sum(c.get(e, 0) for e in HALOGENS)
    return c["C"] - (c.get("H", 0) + halogens) / 2 + c.get("N", 0) / 2 + 1


def is_hydrocarbon(c):
    return c.get("C", 0) > 0 and set(c) <= {"C", "H"}


def is_terpene_formula(c):
    """C10H16, C15H24 or C20H32: the (C5H8)n formulas of mono-, sesqui- and diterpenes."""
    C = c.get("C", 0)
    return is_hydrocarbon(c) and C in (10, 15, 20) and c.get("H", 0) * 5 == C * 8


# --- our classes ----------------------------------------------------------

def possible_classes(c):
    """The classes a composition c does not rule out."""
    C, O, d = c.get("C") or None, c.get("O", 0), dbe(c)
    hetero = any(n > 0 for e, n in c.items() if e not in ("C", "H"))
    if C is None:
        return ["other"]                # inorganic: no class applies
    if c.get("N", 0) > 0:
        return ["nitrogen"]
    if is_terpene_formula(c):
        # C15H24 (DBE 4) and C20H32 (DBE 5) are also alkylbenzene formulas,
        # e.g. tripropylbenzene; C10H16 (DBE 3) cannot be aromatic.
        if C >= MIN_C["aromatic"] and d >= MIN_DBE["aromatic"]:
            return ["terpene", "aromatic"]
        return ["terpene"]
    excluded = {"nitrogen", "terpene"}
    if d is not None and abs(d) > 0.5:
        excluded.add("alcohol")
    if d is not None and not (0.5 <= d <= 3.5):
        excluded.add("alkene")
    if hetero:
        excluded.add("alkene")
    if C == 6 and d is not None and abs(d - 4) <= 0.05:
        excluded |= set(CLASSES) - {"aromatic"}
    for cls, n in MIN_C.items():
        if C is not None and C < n:
            excluded.add(cls)
    for cls, n in MIN_O.items():
        if O < n:
            excluded.add(cls)
    for cls, n in MAX_O.items():
        if O > n:
            excluded.add(cls)
    for cls, n in MIN_DBE.items():
        if d is not None and d < n - 0.01:
            excluded.add(cls)
    return [cls for cls in CLASSES if cls not in excluded] or ["other"]


# --- SOA-yield reference classes ----------------------------------------

def reference_classes(c):
    """
    The SOA-yield reference classes a composition c allows. Only hydrocarbons
    map onto the reference; anything with a heteroatom returns an empty list.
    Terpene formulas (which have their own yields) return only the
    non-terpene isomers the DBE allows: none for C10H16, Aromatic for the
    alkylbenzenes that share the C15H24 / C20H32 formulas.

        DBE 0                   n-Alkane, b-Alkane
        DBE 1                   Alkenes, Cycloalkane
        DBE 2                   Alkenes, Cycloalkene, Cycloalkane (bicyclic)
        DBE 3                   Alkenes, Cycloalkene
        DBE >= 7, >= 10 C       PAH
        DBE >= 4, >= 6 C        Aromatic (all other cases)
    """
    if not is_hydrocarbon(c):
        return []
    C, d = c["C"], dbe(c)
    if is_terpene_formula(c):
        return ["Aromatic"] if "aromatic" in possible_classes(c) else []
    if abs(d - round(d)) > 0.01 or d < 0:
        return []                       # not a closed-shell molecule
    d = int(round(d))
    if d == 0:
        return ["n-Alkane", "b-Alkane"]
    if d == 1:
        return ["Alkenes", "Cycloalkane"]
    if d == 2:
        return ["Alkenes", "Cycloalkene", "Cycloalkane"]
    if d == 3:
        return ["Alkenes", "Cycloalkene"]
    if C >= PAH_MIN_C and d >= PAH_MIN_DBE:
        return ["PAH"]
    if C >= 6:
        return ["Aromatic"]
    return []


# --- one call per formula ------------------------------------------------

def classify(formula):
    """
    The DBE of a neutral formula and every class it allows, ours and the
    reference's. Raises ValueError for anything that is not a formula.
    """
    if not is_formula(formula):
        raise ValueError(f"{formula!r} is not a neutral formula of known elements")
    c = counts(formula)
    return {"formula": formula, "dbe": dbe(c),
            "possible_classes": possible_classes(c),
            "reference_classes": reference_classes(c)}


def classify_all(formulas):
    """classify() for many formulas; wrap in pandas.DataFrame(...) for a table."""
    return [classify(f) for f in formulas]
