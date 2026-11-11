# -*- coding: utf-8 -*-
"""
Mapping a species onto a column of the SOM yield table.

Two routes:

* structure (SMILES, parsed with RDKit) -> exactly one column;
* formula only (element counts + DBE, via ClassEstimationDBE) -> every
  column the formula allows.

All assumptions are module constants so they can be changed in one place
and are listed in SOAP/README.md.
"""
import sys
from pathlib import Path

# ClassEstimationDBE lives next to SOAP in the Tool-Chains repository.
_DBE_DIR = Path(__file__).resolve().parents[2] / "ClassEstimationDBE"
if str(_DBE_DIR) not in sys.path:
    sys.path.insert(0, str(_DBE_DIR))
import class_estimation as ce  # noqa: E402

from .yields import oxygen_column  # noqa: E402

# --- assumptions -------------------------------------------------------
TERPENE_COLUMN = "Cycloalkene"   # terpenes use the alpha-pinene surrogate
# Species with heteroatoms other than O (N, S, P, Si, halogens):
#   "oxygen_column": with >= 1 O use the O-number column like any oxygenate
#                    (the original scripts); without O they get 0;
#   "zero":          always 0.
HETEROATOM_POLICY = "oxygen_column"
OTHER_HETEROATOMS = ("N", "S", "P", "Si", "F", "Cl", "Br", "I")


# --- shared helpers -----------------------------------------------------

def _heteroatom_column(c, flags):
    """Column for a species that is not a pure hydrocarbon, or None (-> SOAFP 0)."""
    others = [e for e in OTHER_HETEROATOMS if c.get(e, 0)]
    n_O = c.get("O", 0)
    if others:
        tag = "+".join(others)
        if HETEROATOM_POLICY == "zero" or n_O == 0:
            flags.append(f"contains {tag}: no SOM parameterization, SOAFP set to 0")
            return None
        flags.append(f"contains {tag}: treated as an oxygenate with {n_O} O (assumption)")
    if n_O >= 1:
        flags.append("oxygenate: O-number column (alkane surrogate), independent of structure")
        return oxygen_column(n_O)
    return None


# --- structure route ----------------------------------------------------

def _rdkit():
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import rdMolDescriptors
    except ImportError as err:          # pragma: no cover
        raise ImportError("the structure route needs RDKit: pip install rdkit") from err
    RDLogger.DisableLog("rdApp.*")
    return Chem, rdMolDescriptors


def hydrocarbon_class(mol):
    """
    SOM hydrocarbon column of an RDKit molecule made of C and H only.

        >= 2 aromatic rings                     PAH
        1 aromatic ring                         Aromatic  (indane, indene, tetralin, styrene)
        terpene formula (C10H16, C15H24, C20H32) Cycloalkene (alpha-pinene surrogate)
        C=C or C#C and a ring                   Cycloalkene
        C=C or C#C, no ring                     Alkenes
        saturated with a ring                   Cycloalkane
        saturated, a C bonded to >= 3 C         b-Alkane
        otherwise                               n-Alkane
    """
    Chem, _ = _rdkit()
    ring_info = mol.GetRingInfo()
    aromatic_rings = sum(
        all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
        for ring in ring_info.AtomRings())
    flags = []
    if aromatic_rings >= 2:
        return "PAH", flags
    if aromatic_rings == 1:
        return "Aromatic", flags
    if ce.is_terpene_formula(_mol_counts(mol)):
        flags.append(f"terpene: {TERPENE_COLUMN} (alpha-pinene) surrogate (assumption)")
        return TERPENE_COLUMN, flags
    bonds = [b.GetBondType() for b in mol.GetBonds() if not b.GetIsAromatic()]
    has_triple = Chem.BondType.TRIPLE in bonds
    unsaturated = has_triple or Chem.BondType.DOUBLE in bonds
    if has_triple:
        flags.append("alkyne: treated as an alkene (assumption)")
    if unsaturated:
        return ("Cycloalkene" if ring_info.NumRings() else "Alkenes"), flags
    if ring_info.NumRings():
        return "Cycloalkane", flags
    branched = any(sum(n.GetSymbol() == "C" for n in a.GetNeighbors()) >= 3
                   for a in mol.GetAtoms() if a.GetSymbol() == "C")
    return ("b-Alkane" if branched else "n-Alkane"), flags


def _mol_counts(mol):
    _, rdMolDescriptors = _rdkit()
    return ce.counts(rdMolDescriptors.CalcMolFormula(mol).rstrip("+-"))


def column_from_smiles(smiles):
    """
    (column or None, counts, flags) for one SMILES string. None means no SOM
    parameterization applies (SOAFP 0). Raises ValueError for an unreadable SMILES.
    """
    Chem, _ = _rdkit()
    mol = Chem.MolFromSmiles(str(smiles).strip()) if smiles else None
    if mol is None:
        raise ValueError(f"cannot parse SMILES {smiles!r}")
    c = _mol_counts(mol)
    flags = []
    if c.get("C", 0) == 0:
        flags.append("no carbon: SOAFP set to 0")
        return None, c, flags
    if ce.is_hydrocarbon(c):
        column, more = hydrocarbon_class(mol)
        return column, c, flags + more
    return _heteroatom_column(c, flags), c, flags


# --- formula route ------------------------------------------------------

def columns_from_formula(formula):
    """
    ([(label, column), ...], counts, dbe, flags) for a neutral formula: every
    SOM column the formula allows. An empty list means SOAFP 0.
    """
    if not ce.is_formula(formula):
        raise ValueError(f"{formula!r} is not a neutral formula of known elements")
    c = ce.counts(formula)
    d = ce.dbe(c)
    flags = []
    if c.get("C", 0) == 0:
        flags.append("no carbon: SOAFP set to 0")
        return [], c, d, flags
    if not ce.is_hydrocarbon(c):
        column = _heteroatom_column(c, flags)
        return ([(column, column)] if column else []), c, d, flags
    candidates = [(cls, cls) for cls in ce.reference_classes(c)]
    if ce.is_terpene_formula(c):
        flags.append(f"terpene formula: {TERPENE_COLUMN} (alpha-pinene) surrogate (assumption)")
        candidates.insert(0, ("terpene", TERPENE_COLUMN))
    if not candidates:
        flags.append(f"DBE={d}: no hydrocarbon class fits, SOAFP set to 0")
    elif len(candidates) > 1:
        flags.append("formula is ambiguous: range over "
                     + ", ".join(label for label, _ in candidates))
    return candidates, c, d, flags
