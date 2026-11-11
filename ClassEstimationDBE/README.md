# ClassEstimationDBE

Estimates a compound's class from its **neutral molecular formula** alone,
using element counts and the double bond equivalent (DBE). It makes no
library lookup and no assumption about how the formula was obtained. It
returns **every class the formula allows**, not a best guess. Deciding what to
do when several are possible is up to the calling code.

The module is a single file with no dependencies: `class_estimation.py`.

## Usage

```python
import sys
sys.path.insert(0, r"C:\Users\exp.voc\PythoRepos\Tool-Chains\ClassEstimationDBE")
import class_estimation as ce

ce.classify("C4H8O")
# {'formula': 'C4H8O', 'dbe': 1.0,
#  'possible_classes': ['ketone', 'aldehyde'], 'reference_classes': []}

import pandas as pd
table = pd.DataFrame(ce.classify_all(["C7H8", "C10H8", "C9H18", "C4H8O"]))
```

**Input:** the formula of the neutral molecule, written with element symbols
and counts (`C2H6O`). Repeated elements are summed, so `CH3CH2OH` reads the
same as `C2H6O`. Supported elements: C, H, N, O, F, Si, P, S, Cl, Br, I.
Anything else raises a `ValueError`: a charge sign, an unknown element, an
isotope label or a bare mass.

**Output:**

| field | meaning |
|---|---|
| `dbe` | C − (H + halogens)/2 + N/2 + 1; O, S, Si and P do not count. `None` without carbon |
| `possible_classes` | every one of **our** classes that the formula doesn't rule out; `["other"]` if none |
| `reference_classes` | every **SOA-reference** class the formula allows; empty if none fits |

The building blocks are also available on their own: `counts()` turns a
formula into an element-count dict such as `{"C": 7, "H": 8}`, and `dbe()`,
`possible_classes()` and `reference_classes()` take such a dict.

## Our classes

Each rule rules classes *out*. Whatever is left is returned.

| class | rule |
|---|---|
| nitrogen | any N (if present, nothing else is returned) |
| terpene | C/H only and exactly C10H16, C15H24 or C20H32 (only `aromatic` can also be returned: C15H24 and C20H32 are alkylbenzene formulas too) |
| aromatic | ≥ 6 C, DBE ≥ 4; a C6 with DBE = 4 is aromatic only |
| alkene | C/H only, ≥ 2 C, 1 ≤ DBE ≤ 3 |
| ketone | ≥ 3 C, 1 O, DBE ≥ 1 |
| aldehyde | ≥ 1 C, 1 O, DBE ≥ 1 |
| ester | ≥ 2 C, ≥ 2 O, DBE ≥ 1 |
| acid | ≥ 1 C, ≥ 2 O, DBE ≥ 1 |
| alcohol | 1–2 O, DBE = 0 |
| other | nothing above fits; this **includes saturated hydrocarbons** (DBE 0) and anything without carbon |

Si, S and P are left out of the DBE. So organosilicon compounds such as
siloxanes get a low or negative DBE (C6H18OSi2 → −2) and come out as
`other`, which is intended.

## Mapping to the SOA-yield reference classes

The reference uses seven hydrocarbon classes: n-Alkane, b-Alkane,
Cycloalkane, Aromatic, PAH, Alkenes and Cycloalkene. A formula only gives the
DBE (rings + double bonds), not where they are, so several reference classes
often share one formula. `reference_classes()` returns all of them:

| formula (C/H only) | DBE | our class | reference classes returned | resolved? |
|---|---|---|---|---|
| CnH2n+2 | 0 | other | n-Alkane, b-Alkane | no: these are isomers |
| CnH2n | 1 | alkene | Alkenes, Cycloalkane | no |
| CnH2n−2 | 2 | alkene | Alkenes (dienes), Cycloalkene, Cycloalkane (bicyclic, e.g. decalin) | no |
| CnH2n−4 | 3 | alkene | Alkenes (trienes), Cycloalkene | no |
| ≥ 10 C | ≥ 7 | aromatic | PAH | **yes** |
| ≥ 6 C, not PAH | ≥ 4 | aromatic | Aromatic | **yes** |
| C10H16 | 3 | terpene | none: terpenes have their own SOA yields | — |
| C15H24, C20H32 | 4–5 | terpene, aromatic | Aromatic (the alkylbenzene isomers, e.g. tripropylbenzene); the terpene isomers need their own yields | no |
| any heteroatom (O, N, S, …) | any | ketone … nitrogen | none: the reference covers hydrocarbons only | — |

Notes on the cut-offs:

- **Aromatic vs PAH.** The PAH cut is DBE ≥ 7 with ≥ 10 C, set by naphthalene,
  C10H8. Aromatics with one extra ring or double bond (styrene, indane, indene,
  DBE 5–6) stay as Aromatic. Both thresholds are module constants:
  `PAH_MIN_DBE` and `PAH_MIN_C`.
- **n- vs b-alkane can never be separated by formula.** Every isomer has the
  same formula.
- **Oxygenated aromatics** such as phenols, cresols and benzaldehyde fall
  under our `aromatic` class but get **no** reference class, because they
  contain O.
- **Isoprene (C5H8)** isn't covered by the terpene rule. It is a DBE 2
  hydrocarbon, so it comes back as Alkenes, Cycloalkene, Cycloalkane. If your
  reference gives isoprene its own yield, handle it before using this mapping.
- A tricyclic alkane at DBE 3, a hydrocarbon with DBE ≥ 4 and fewer than 6 C,
  and any formula with a non-integer or negative DBE are not mapped to any
  reference class.

## Building SOA estimates on top (not implemented here)

What a caller needs in order to turn the candidate lists into bounds:

1. **A lookup table.** One row per reference class. Columns: the SOA yield,
   plus whatever it depends on in the reference (for example carbon number).
   If the yields depend on carbon number, key the lookup on class and carbon
   number; the carbon number comes from elsewhere, not from this module.
2. **Lower and upper bound per compound.** Look up the yield of every entry
   in `reference_classes`:
   - lower = the minimum over the candidates;
   - upper = the maximum over the candidates;
   - when there is only one candidate, the two bounds are equal.

   Never average the candidates: the formula gives no information about
   which one applies.
3. **Empty `reference_classes`.** Terpenes and compounds with heteroatoms get
   no reference class. The caller has to decide what happens to them, either
   with its own yields or by excluding them. Report that choice rather than
   dropping these compounds silently.
