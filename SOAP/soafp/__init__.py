# -*- coding: utf-8 -*-
"""
soafp: SOA formation potential (SOM SOA mass yield) of individual species.

    import sys; sys.path.insert(0, r"<Tool-Chains>/SOAP")
    import soafp

    soafp.soafp_from_structure("CC1=CC=CC=C1")        # toluene -> one value
    soafp.soafp_range_from_formula("C10H18")          # formula -> range
    soafp.soafp_table(df, smiles="SMILES", formula="Formula")   # pandas table

Yield tables are plain CSV in soafp/data; no spreadsheet is needed.
The structure route needs RDKit; the formula route needs only the standard
library (and ClassEstimationDBE from this repository).
See SOAP/README.md for the method and every assumption.
"""
from .core import soafp_from_class, soafp_from_structure, soafp_range_from_formula
from .table import soafp_one, soafp_table
from .yields import CASES, COLUMNS, mass_yield

__all__ = ["soafp_from_structure", "soafp_range_from_formula", "soafp_from_class",
           "soafp_one", "soafp_table", "mass_yield", "COLUMNS", "CASES"]
