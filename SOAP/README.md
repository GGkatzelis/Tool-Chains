# SOAFP: SOA formation potential from structure or formula

Assigns each organic species a **secondary organic aerosol formation potential
(SOAFP)**. Here SOAFP is the **SOA mass yield** (g SOA per g precursor reacted)
from Statistical Oxidation Model (SOM) box-model runs, looked up by compound
class and carbon number (or oxygen number).

There are two routes:

| You know | Route | Result |
|---|---|---|
| the structure (SMILES) | `soafp_from_structure` | one SOM class → one value, plus its vapor-wall-loss (VWL) range |
| only the molecular formula (e.g. PTR-ToF-MS, CIMS) | `soafp_range_from_formula` | DBE → every class the formula allows → a min–max range |

The module needs **no spreadsheet**: the yield tables are plain CSV in
`soafp/data/`. The structure route needs [RDKit](https://www.rdkit.org/)
(`pip install rdkit`). The formula route needs only the standard library and
`ClassEstimationDBE` from this repository.

---

## 1. Quick start

```python
import sys
sys.path.insert(0, r"<path to Tool-Chains>/SOAP")
import soafp

r = soafp.soafp_from_structure("CC1=CC=CC=C1")       # toluene
r["column"], r["soafp_avg"], r["soafp_low"], r["soafp_high"]
# ('Aromatic', 0.0896, 0.0670, 0.1122)

r = soafp.soafp_range_from_formula("C10H18")
r["column"], r["soafp_avg_low"], r["soafp_avg_high"]
# ('Alkenes; Cycloalkene; Cycloalkane', 0.0867, 0.2005)

import pandas as pd
out = soafp.soafp_table(df, smiles="SMILES", formula="Formula")   # whole table
```

To run the four source profiles (gasoline, diesel, VCP, cooking):

```
python run_profiles.py            # -> results/<profile>_soafp.csv, results/<profile>_input_issues.csv
python -m pytest tests            # 31 tests
```

## 2. Files

| Path | Content |
|---|---|
| `soafp/yields.py` | reads the yield tables, `mass_yield(column, n_C)` |
| `soafp/classify.py` | structure → class (RDKit); formula → candidate classes (DBE). **All assumptions are constants here.** |
| `soafp/core.py` | `soafp_from_structure`, `soafp_range_from_formula`, `soafp_from_class` |
| `soafp/table.py` | `soafp_one`, `soafp_table` (pandas) |
| `soafp/data/som_mass_yield_noVWL.csv`, `..._VWL.csv` | SOM mass yields, C5–C30 × 12 columns, with provenance header |
| `profiles/*.csv` | the four speciation profiles (exported from `legacy/Tables for SOAFP calculation.xlsx`) |
| `run_profiles.py` | SOAFP for the four profiles, with the formula-only range alongside |
| `tests/test_soafp.py` | tests |
| `legacy/` | the original scripts and workbooks by Yi Wu (they still run, see §8) |

## 3. The yield tables

Two SOM box-model runs, both from the workbook `Copy of SOM_boxmodel(for Shantanu).xlsx`:

| Case | SOM parameters fit to chamber data… | Use |
|---|---|---|
| `noVWL` | without vapor-wall-loss correction | lower bound |
| `VWL` | with vapor-wall-loss correction | upper bound |
| `avg` | mean of the two (= sheets `avg` / `avg_py`, checked to be identical) | central value, as in the original scripts |

Columns and surrogate precursors, as in the workbook headers:

| Column | Surrogate |
|---|---|
| `n-Alkane` | n-dodecane |
| `b-Alkane`, `Cycloalkane` | not stated in the workbook. Probably the C12 alkane fits of Cappa et al. (2013): 2-methylundecane, hexylcyclohexane. **To confirm.** |
| `Aromatic` | benzene (n_C ≤ 6), toluene (7), m-xylene (≥ 8). The column therefore jumps between C6, C7 and C8. |
| `PAH` | naphthalene |
| `Alkenes` | isoprene |
| `Cycloalkene` | α-pinene |
| `O1` … `O5` | alkane surrogate with 1…5 (≥ 5) O atoms. Labelled m-undecane in `no VWL` and dodecane in `VWL`. Constant above C20, so extrapolated. |

**Provenance.** The workbook is labelled "SOM II (Updated 11/23/2016)" and
contains a block "Table S7b. Zhao (n-alkanes, SOA mass yield)". The SOM
framework is Cappa & Wilson (2012) and Cappa et al. (2013). The VWL-corrected
fits follow Zhang et al. (2014) and Cappa et al. (2016, Table S1). The same
family of yields is used by McDonald et al. (2018, Table S12) for VCPs. See
`LITERATURE_REVIEW.md`.

**Not documented in the workbook:**
- the OH exposure, OA loading, NOx regime and temperature of the runs. A margin
  note reads "correct for OH exposure at the time considered";
- which wall-loss rate the `VWL` case uses. Cappa et al. (2016) give 1×10⁻⁴
  and 2.5×10⁻⁴ s⁻¹.

## 4. Method

### 4.1 Structure route (SMILES → one class)

The SMILES is parsed with RDKit, so aromatic (`c1ccccc1`), Kekulé and stereo
SMILES all work. For C/H-only species the rules are applied in this order:

| Rule | Column |
|---|---|
| ≥ 2 aromatic rings | `PAH` (naphthalenes, biphenyl) |
| 1 aromatic ring | `Aromatic`. This includes indane, indene, tetralin and styrene. |
| terpene formula C10H16 / C15H24 / C20H32, not aromatic | `Cycloalkene` (α-pinene). This includes limonene and acyclic terpenes. |
| C=C or C≡C, with a ring | `Cycloalkene` |
| C=C or C≡C, no ring | `Alkenes`. Alkynes are flagged. |
| saturated, with a ring | `Cycloalkane` |
| saturated, a carbon bonded to ≥ 3 carbons | `b-Alkane` |
| otherwise | `n-Alkane` |

Species with oxygen use the `O<n_O>` column at their carbon number,
independent of the carbon skeleton. Species with N, S, P, Si or halogens follow
`HETEROATOM_POLICY` (see §5).

A formula given together with the SMILES is only checked against it: a
mismatch is flagged and the structure is used. Given `n_C` / `n_O` values
override the structure (for lumped species, see §7).

### 4.2 Formula route (formula → candidate classes → range)

`ClassEstimationDBE` computes the DBE and lists every class the formula
allows:

| C/H formula | DBE | Candidates |
|---|---|---|
| CₙH₂ₙ₊₂ | 0 | n-Alkane, b-Alkane |
| CₙH₂ₙ | 1 | Alkenes, Cycloalkane |
| CₙH₂ₙ₋₂ | 2 | Alkenes, Cycloalkene, Cycloalkane (bicyclic) |
| CₙH₂ₙ₋₄ | 3 | Alkenes, Cycloalkene |
| ≥ 6 C | ≥ 4 | Aromatic; PAH if DBE ≥ 7 and ≥ 10 C |
| C10H16 | 3 | Cycloalkene (terpene) |
| C15H24, C20H32 | 4–5 | Cycloalkene (terpene), Aromatic (alkylbenzenes) |
| with O | any | `O<n_O>`: a single value, because the O-number columns ignore structure |

The range is reported two ways:
- `soafp_avg_low` – `soafp_avg_high`: over the candidates, using the `avg` table;
- `soafp_low` – `soafp_high`: over the candidates **and** both VWL cases.

The candidates are never averaged, because the formula carries no information
on which one applies.

## 5. Assumptions (all in `soafp/classify.py` and `soafp/yields.py`)

| # | Assumption | Default | Why / alternative |
|---|---|---|---|
| A1 | SOAFP = SOM mass yield at the (undocumented) box-model conditions | — | not reactivity- or exposure-weighted; see LITERATURE_REVIEW.md |
| A2 | central value = mean of no-VWL and VWL | `avg` | the two cases are reported as the bounds |
| A3 | n_C < 5 → 0 | | the SOM tables start at C5 |
| A4 | n_C > 30 → C30 value | | yields plateau above about C23. The legacy scripts gave 0. |
| A5 | oxygenates: `O<min(n_O,5)>` column, regardless of structure | | aromatic or alkene oxygenates (phenols, aldehydes…) get the alkane-based value |
| A6 | terpenes → `Cycloalkene` (α-pinene) | `TERPENE_COLUMN` | applies to sesquiterpenes too, scaled only by carbon number |
| A7 | species with N, S, P, Si or halogens: with O they are treated as oxygenates (e.g. siloxanes → `O5`); without O they get 0 | `HETEROATOM_POLICY = "oxygen_column"` | `"zero"` gives them all 0 |
| A8 | one aromatic ring = `Aromatic`; only ≥ 2 aromatic rings = `PAH` | | the legacy SMILES rules put indane, indene and tetralin in PAH, which is about ×2 higher at C9–C10 |
| A9 | alkynes → as alkenes | | flagged; mostly < C5 anyway |
| A10 | formula route: DBE ≥ 4 with ≥ 6 C is aromatic | `ce.PAH_MIN_DBE`, `PAH_MIN_C` | non-aromatic polyenes or polycyclics at DBE ≥ 4 are not considered |

Every assumption that applies to a species is written to its `soafp_flags`.

## 6. Output columns

| Column | Meaning |
|---|---|
| `soafp_method` | `structure`, `class` (set by hand), `formula`, or `none` |
| `soafp_column` | SOM column(s) used; `; `-separated when there are several candidates |
| `n_C`, `n_O`, `dbe` | composition used for the lookup |
| `soafp_noVWL`, `soafp_VWL`, `soafp_avg` | single-class values (empty for ambiguous formulas) |
| `soafp_avg_low`, `soafp_avg_high` | range of `soafp_avg` over the candidates |
| `soafp_low`, `soafp_high` | range over the candidates and both VWL cases |
| `soafp_flags` | every assumption or problem, `; `-separated |
| `formula_*` (run_profiles only) | the same, computed from the formula alone |

## 7. The four profiles (`run_profiles.py`)

**Given carbon and oxygen numbers.** The diesel and VCP tables carry per-row
carbon and oxygen numbers (`n_C`/`n_O`, `No. C`/`No. O`). They are handled
differently:

| Profile | Policy | Why |
|---|---|---|
| diesel | **override** the structure | the mismatches are lumped species whose SMILES is only a representative: "Aromatic" C15…C23 all carry the SMILES of a C14 alkylbenzene, and "UCM" has `n_C` = 26 with the SMILES of dodecane |
| VCP | **check only**, structure used | the mismatches are mostly wrong numbers for named compounds: 3-methylpentane is given C4, benzaldehyde 0 O, dipropylene glycol methyl ether acetate C9 |

The legacy VCP script used `No. C`/`No. O` and the legacy diesel script used
the formula, so both policies differ from legacy for some rows. Compound names
matching "branched alkane", "cycloalkane" or "cycloalkene" set the class by
hand (`NAME_HINTS`), as the legacy scripts did. Every such row is listed in
`results/<profile>_input_issues.csv` for review.

**Weighted mean SOAFP, Σ(wt % × SOAFP)/Σ wt %:**

| Profile | Species changed | Rows to review | Legacy | This module (avg) | VWL range | Formula only, avg | Formula only, incl. VWL |
|---|---|---|---|---|---|---|---|
| gasoline | 5/137 | 2 | 0.0284 | 0.0283 | 0.0220–0.0345 | 0.0280–0.0309 | 0.0190–0.0414 |
| diesel | 20/217 | 98 | 0.2270 | **0.3592** | 0.3163–0.4020 | 0.2089–0.2604 | 0.1640–0.3462 |
| VCP | 35/512 | 31 | 0.0379 | **0.0332** | 0.0222–0.0442 | 0.0310–0.0386 | 0.0196–0.0550 |
| cooking | 1/46 | 11 | 0.0029 | 0.0029 | 0.0014–0.0044 | 0.0029 | 0.0014–0.0044 |

**Two rows dominate the differences:**
- **diesel "UCM"** (11.9 wt %): the legacy script parsed the placeholder
  formula C12H26 and got 0.069. The table's `n_C` = 26 gives 1.05, which
  explains almost all of the diesel increase. *Please confirm that 26 is meant.*
- **VCP "Witch hazel"** (0.54 wt %): `No. C` = 20 and `No. O` = 14 in the
  table, but the SMILES and formula are decane. Legacy took C20 (0.66); this
  module takes decane (0.033). Neither is a good representation of a plant
  extract.

The other changes across the profiles:
- indane, indene and tetralin → Aromatic (A8);
- terpenes → Cycloalkene (limonene, pinene, pine oil and turpentine were
  Aromatic or PAH before);
- the diesel "Aromatic" and "Polycyclic Aromatic" lumps now use their given
  carbon numbers;
- some rows where the SMILES, the formula or the given numbers disagree.

The formula-only range does not reach the diesel value, because the lumps'
placeholder formulas carry the wrong carbon numbers.

**Input rows to check** (all are in `results/*_input_issues.csv`):
- gasoline "napthalene": formula C10H8, but the SMILES is methylnaphthalene,
  which is used;
- cooking undecanal/undecenal: the formulas and SMILES are swapped or wrong;
- cooking: eight species have no SMILES (`---`), so the formula range is used;
- VCP: Witch hazel (above); Texanol `No. O` is 4 for a molecule with 3 O;
  "C9 cycloalkanes" and "Cis-bicyclo[4.3.0]nonane" carry the SMILES of other
  compounds; siloxane lumps have inconsistent `No. C`.

## 8. Legacy scripts

`legacy/` holds the original per-source scripts (`SOAFP_calculation_*.py`) and
both workbooks unchanged. The only edit is to the file names in the scripts,
which now match the uploaded workbooks. Run them from that folder; they write
`WRF Speciation with SOAFP_*.xlsx` / `SOAFP_cooking.xlsx` next to themselves
(ignored by git). Compared with this module, they:
- classify from SMILES with string matching, so aromatic or stereo SMILES fail
  and indane etc. come out as PAH;
- distinguish cycloalkenes only through the `Speciate CRACMM` column;
- in the cooking script, use no structure at all (only the O number and the
  terpene flag);
- give 0 above C30.

## 9. Open items

1. Document the box-model conditions (OH exposure, OA loading, NOx, T) and the
   VWL rate. Confirm the b-Alkane and Cycloalkane surrogates.
2. Confirm the diesel UCM carbon number and fix the input rows listed in §7.
   Better still, give lumped species an explicit class and carbon number
   instead of a placeholder SMILES.
3. Possible method upgrades are ranked in `LITERATURE_REVIEW.md`.

## References

- Cappa, C. D. and Wilson, K. R.: Multi-generation gas-phase oxidation, equilibrium partitioning, and the formation and evolution of secondary organic aerosol, Atmos. Chem. Phys., 12, 9505–9528, 2012.
- Cappa, C. D., et al.: Application of the Statistical Oxidation Model (SOM) to SOA formation from photooxidation of C12 alkanes, Atmos. Chem. Phys., 13, 1591–1606, doi:10.5194/acp-13-1591-2013, 2013.
- Cappa, C. D., et al.: Simulating secondary organic aerosol in a regional air quality model using the statistical oxidation model – Part 2: Assessing the influence of vapor wall losses, Atmos. Chem. Phys., 16, 3041–3059, doi:10.5194/acp-16-3041-2016, 2016 (Table S1: SOM parameters).
- Zhang, X., et al.: Influence of vapor wall loss in laboratory chambers on yields of secondary organic aerosol, PNAS, 111, 5802–5807, doi:10.1073/pnas.1404727111, 2014.
- McDonald, B. C., et al.: Volatile chemical products emerging as largest petrochemical source of urban organic emissions, Science, 359, 760–764, doi:10.1126/science.aaq0524, 2018 (Table S12: SOA yields).

Original scripts: Yi Wu. Module: Tool-Chains.
