# SOAFP: SOA formation potential from structure or formula

Assigns each organic species a **secondary organic aerosol formation potential
(SOAFP)**. Here SOAFP is the **SOA mass yield** (g SOA per g precursor reacted)
from Statistical Oxidation Model (SOM) box-model runs, looked up by compound
class and carbon number (or oxygen number). The SOM parameters are fits to
chamber experiments, so the tool adds no chemistry of its own: it only decides
which SOM class a species belongs to, and reports how uncertain that choice is.

| You know | Function | Result |
|---|---|---|
| a list of molecular formulas (e.g. PTR-ToF-MS, CIMS peak list) | `soafp_from_formulas` | per formula: central value, low–high range, number of candidate classes, assumptions |
| one formula | `soafp_range_from_formula` | the same, as a dict |
| the structure (SMILES) | `soafp_from_structure` | one SOM class → one value, with its vapor-wall-loss (VWL) range |

Like a log C* estimate from a formula, every value comes with an uncertainty
and a list of the assumptions behind it.

The module needs **no spreadsheet**: the yield tables are plain CSV in
`soafp/data/`. The formula route needs only the standard library (pandas for
`soafp_from_formulas`) and `ClassEstimationDBE` from this repository. The
structure route needs [RDKit](https://www.rdkit.org/) (`pip install rdkit`).

---

## 1. Quick start

### From the command line: `calc_soafp.py`

```
python calc_soafp.py C10H16 C7H8O "CC1=CC=CC=C1"           # formulas and SMILES, mixed
python calc_soafp.py limonene "benzyl alcohol" --online     # compound names (PubChem)
python calc_soafp.py --file peaks.csv --column Formula --out peaks_soafp.csv   # CSV or .xlsx
python calc_soafp.py                                        # help: inputs and limits
```

With `--file`, every column of the input is kept and the `soafp_*` columns
are appended, so peak lists, inventories etc. pass straight through. In
Python: `from calc_soafp import soafp_of, soafp_of_list`.

**What each entry can be**

| Input | Example | Notes |
|---|---|---|
| neutral formula | `C10H16`, `C7H8O`, `C6H5Cl`, `CH3CH2OH` | elements C H O N S P Si F Cl Br I. **Not** ion formulas: give C10H16, not C10H17+ (refused with a message). Isotopes and charges are not supported. |
| SMILES | `CC1=CC=CC=C1`, `O=Cc1ccccc1` | needs RDKit; gives the class directly, so the range is narrower |
| compound name | `limonene` | only with `--online`: looked up on PubChem (only the name is sent). Check the returned `smiles`: trade and mixture names can resolve to the wrong compound. |

Auto-detection: a valid formula that contains H is read as a formula;
anything else is tried as SMILES, then as a name. `CCCCCC` is therefore
hexane (SMILES), and `CO` is methanol. Force the type with
`--kind formula|smiles|name` (e.g. `CS2 --kind formula`).

**Limits of the input.** Lumped mixtures (UCM, "C18 aromatics", plant
extracts, siloxane blends) are not single species: give a representative
formula with the right carbon number. Species below C5 get 0, above C30 the
C30 value. Entries that are not understood get no value and the reason in
`flags`; they are never silently set to 0.

### From Python: the `soafp` package

```python
import sys
sys.path.insert(0, r"<path to Tool-Chains>/SOAP")
import soafp

soafp.soafp_from_formulas(["C10H16", "C10H18", "C7H8O", "C8H16O2"])
#  formula  soafp_column                       soafp_central  soafp_low  soafp_high  n_candidates
#  C10H16   Cycloalkene                               0.1323     0.0930      0.1716             1
#  C10H18   Alkenes; Cycloalkene; Cycloalkane         0.1399     0.0512      0.3304             3
#  C7H8O    O1; Aromatic                              0.0453     0.0001      0.1122             2
#  C8H16O2  O2                                        0.0058     0.0012      0.0104             1
#  (plus soafp_minus, soafp_plus, dbe, n_C, n_O, soafp_flags, ...)

r = soafp.soafp_from_structure("CC1=CC=CC=C1")       # toluene
r["column"], r["soafp_central"], r["soafp_low"], r["soafp_high"]
# ('Aromatic', 0.0896, 0.0670, 0.1122)

out = soafp.soafp_table(df, smiles="SMILES", formula="Formula")   # whole pandas table
```

```
python run_profiles.py            # the four source profiles -> results/<profile>_soafp.csv
python validate_formula_route.py  # formula route vs known structures -> results/formula_validation.csv
python -m pytest tests            # 36 tests
```

## 2. Files

| Path | Content |
|---|---|
| `soafp/yields.py` | reads the yield tables, `mass_yield(column, n_C)` |
| `soafp/classify.py` | structure → class (RDKit); formula → candidate classes (DBE). **All assumptions are constants here.** |
| `soafp/core.py` | `soafp_from_structure`, `soafp_range_from_formula`, `soafp_from_class` |
| `soafp/table.py` | `soafp_from_formulas`, `soafp_one`, `soafp_table` (pandas) |
| `calc_soafp.py` | command-line / one-call front end: formulas, SMILES or names, single or from a CSV/Excel column |
| `soafp/data/som_mass_yield_noVWL.csv`, `..._VWL.csv` | SOM mass yields, C5–C30 × 12 columns, with provenance header |
| `profiles/*.csv` | the four speciation profiles (exported from `legacy/Tables for SOAFP calculation.xlsx`) |
| `run_profiles.py` | SOAFP for the four profiles, with the formula-only value alongside |
| `validate_formula_route.py` | tests the formula route on the ~900 profile species of known structure (§7) |
| `tests/test_soafp.py` | tests |
| `LITERATURE_REVIEW.md` | provenance of the yield tables and background on the method |
| `legacy/` | the original scripts and workbooks by Yi Wu (they still run, see §10) |

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

### 4.1 Formula route (formula → DBE → candidate classes → value ± range)

`ClassEstimationDBE` computes the double-bond equivalent,
DBE = C − H/2 + N/2 + 1 (halogens count as H), and the SOM classes the formula
allows:

| Formula | Condition | Candidate SOM columns |
|---|---|---|
| CₙH₂ₙ₊₂ | DBE 0 | n-Alkane, b-Alkane |
| CₙH₂ₙ | DBE 1 | Alkenes, Cycloalkane |
| CₙH₂ₙ₋₂ | DBE 2 | Alkenes, Cycloalkene, Cycloalkane (bicyclic) |
| CₙH₂ₙ₋₄ | DBE 3 | Alkenes, Cycloalkene |
| CₓHᵧ | DBE ≥ 4, ≥ 6 C | Aromatic; PAH if DBE ≥ 7 and ≥ 10 C |
| C10H16 | terpene | Cycloalkene (α-pinene) |
| C15H24, C20H32 | terpene or alkylbenzene | Cycloalkene, Aromatic |
| CₓHᵧOz | DBE < 4 or < 6 C | `O<n_O>` (aliphatic oxygenate) |
| CₓHᵧOz | DBE ≥ 4, ≥ 6 C | `O<n_O>` and Aromatic (aromatic oxygenate: phenols, benzaldehyde, benzyl alcohol…) |
| with F, Cl, Br, I | | 0, and the candidates of the formula with each halogen replaced by H |
| with N, S, P, Si | | 0, and `O<n>` with n = n_O + n_N + n_S + n_P (each heteroatom one O-bearing group; at least 1); plus Aromatic if DBE ≥ 4 and ≥ 6 C |
| < 5 C | | 0 (below the SOM table) |
| > 30 C | | the C30 value |
| unreadable formula | | no value, flagged (never a silent 0) |

From the candidates:

- **`soafp_central`** = mean over the candidates of the `avg` yield (mean of
  no-VWL and VWL). Equal weights by default: the formula alone says nothing
  about which isomer it is. If a project knows more (e.g. from an emission
  inventory), pass `weights={"b-Alkane": 3, "n-Alkane": 1}`; labels not
  given keep weight 1. The arithmetic mean is used so that weight-fraction
  sums over many species stay unbiased.
- **`soafp_low`, `soafp_high`** = min and max over the candidates **and**
  both VWL cases. The range therefore combines the two uncertainties the SOM
  table carries: which class (structure), and the chamber wall loss.
- **`soafp_minus`, `soafp_plus`** = central − low and high − central, for
  reporting as central (−minus / +plus).
- **`n_candidates`** = how many classes the formula allows. With 1 the range is
  the wall-loss range alone.

### 4.2 Structure route (SMILES → one class)

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
independent of the carbon skeleton. If they contain a benzene ring (six
aromatic carbons; furans and pyridines do not count), the Aromatic column is
added as a second candidate, as on the formula route (A6). Species with N, S, P
or Si and at least one O are treated as oxygenates; halocarbons and heteroatom
species without O get 0 (`HETEROATOM_POLICY`). With one candidate,
`soafp_central` equals `soafp_avg` and the range is the wall-loss range.

A formula given together with the SMILES is only checked against it: a
mismatch is flagged and the structure is used. Given `n_C` / `n_O` values
override the structure (for lumped species, see §9).

## 5. Assumptions (all in `soafp/classify.py` and `soafp/yields.py`)

| # | Assumption | Setting | Why / alternative |
|---|---|---|---|
| A1 | SOAFP = SOM mass yield at the (undocumented) box-model conditions | — | per mass reacted; not reactivity- or exposure-weighted |
| A2 | central value = mean of no-VWL and VWL; the two runs are the bounds | `avg` | |
| A3 | n_C < 5 → 0 | `MIN_C` | the SOM tables start at C5 |
| A4 | n_C > 30 → C30 value | `MAX_C` | yields plateau above about C23. The legacy scripts gave 0. |
| A5 | oxygenates use the `O<min(n_O,5)>` column (alkane surrogate) | | the only oxygenate fits in the table |
| A6 | aromatic oxygenates: Aromatic column as a second candidate. Formula route: CHO (and N/S) with DBE ≥ 4 and ≥ 6 C; structure route: an oxygenate with a benzene ring | `OXY_AROMATIC_MIN_C`, `OXY_AROMATIC_MIN_DBE` | the O-columns are far below measured yields of aromatic oxygenates; the Aromatic column of the same C bounds them from above |
| A7 | terpenes → `Cycloalkene` (α-pinene) | `TERPENE_COLUMN` | applies to sesquiterpenes too, scaled only by carbon number |
| A8 | formula route, halogens: 0 to the H-substituted analogue | `HALOGENS` | no SOM fit for halocarbons |
| A9 | formula route, N/S/P/Si: 0 to an oxygenate with each N, S, P as one O-group (and Aromatic if DBE ≥ 4) | `GROUP_ATOMS` | no SOM fit; amines, nitro, thio… groups taken as comparable to one oxygenated group |
| A10 | structure route, N/S/P/Si with O: as an oxygenate (siloxanes → `O5`); without O, and halocarbons: 0 | `HETEROATOM_POLICY = "oxygen_column"` | `"zero"` gives them all 0 |
| A11 | one aromatic ring = `Aromatic`; only ≥ 2 aromatic rings = `PAH` | | the legacy SMILES rules put indane, indene and tetralin in PAH, which is about ×2 higher at C9–C10 |
| A12 | alkynes → as alkenes | | flagged; mostly < C5 anyway |
| A13 | formula route: candidates equally likely | `weights=None` | pass `weights` when the isomer mix is known |

Every assumption that applies to a species is written to its `soafp_flags`.

## 6. Output columns

| Column | Meaning |
|---|---|
| `soafp_method` | `formula`, `structure`, `class` (set by hand), or `none` (no usable input) |
| `soafp_column` | SOM column(s) of the candidates, `; `-separated; `none` = a candidate with SOAFP 0 |
| `n_C`, `n_O`, `dbe` | composition used for the lookup |
| `soafp_central` | best estimate (§4.1) |
| `soafp_low`, `soafp_high` | range over the candidates and both VWL cases |
| `soafp_minus`, `soafp_plus` | central − low, high − central |
| `n_candidates` | number of candidate classes |
| `soafp_noVWL`, `soafp_VWL`, `soafp_avg` | single-class values (empty when there are several candidates) |
| `soafp_avg_low`, `soafp_avg_high` | range of `avg` over the candidates (wall loss averaged out) |
| `soafp_flags` | every assumption or problem, `; `-separated |
| `formula_*` (run_profiles only) | central, low, high, candidates from the formula alone |

## 7. How good is the formula route? (`validate_formula_route.py`)

The 903 profile species with a readable SMILES were run twice: through the
structure route (reference) and through the formula route, with the formula
taken from the SMILES.

| Set | Species | Reference in [low, high] | Central within ×2 of reference | Median central / reference | Weighted reference | Weighted formula central |
|---|---|---|---|---|---|---|
| gasoline | 137 | 100 % | 83 % | 1.34 | 0.0284 | 0.0296 |
| diesel | 217 | 100 % | 99.5 % | 1.01 | 0.2139 | 0.2350 |
| VCP | 512 | 100 % | 94 % | 1.00 | 0.0336 | 0.0344 |
| cooking | 37 | 100 % | 96 % | 1.00 | 0.0013 | 0.0023 |
| all | 903 | 100 % | 94 % | 1.00 | | |

Weighted = Σ(wt % × SOAFP)/Σ wt % over these species (the diesel lumps are
evaluated at their SMILES here, not at their given carbon numbers).

Reading the table:
- The range always contains the structure-based value. This is by
  construction: the structure's class is always one of the formula's
  candidates. The check confirms that the DBE rules and the structure rules
  agree.
- For 94 % of species the central value is within ×2. Of the 49 misses, 42
  are branched alkanes (the n-Alkane candidate is higher, median ratio 1.34),
  6 are N- or S-containing species (the formula route adds a 0 candidate that
  the structure route does not), and one is decatrienal. Cycloalkanes are
  biased high too (median 1.6×, from the Alkenes / Cycloalkene candidates at
  DBE 1–2) but stay within ×2.
- On profile totals the formula route is 2–10 % high for gasoline, diesel and
  VCP. Cooking is higher (0.0023 vs 0.0013), almost entirely from one
  species: decatrienal (C10H14O, DBE 4), whose formula also allows aromatic
  isomers.

## 8. Limitations

- **Per mass reacted, at fixed conditions.** The yields are for the SOM run
  conditions (OH exposure, OA loading, NOx, temperature; see §3), not per mass
  emitted. To get SOA from emissions, multiply by the fraction reacted.
- **The range is not a confidence interval.** It is min–max over the classes
  the formula allows and over the two wall-loss fits. It covers the true value
  only if the species really belongs to one of the candidate classes, and only
  as far as SOM reproduces the chamber data for that class.
- **Oxygenates rest on an alkane surrogate.** The O-columns are fits for
  oxygenated alkanes and ignore structure. Aromatic oxygenates get the Aromatic
  column as an upper candidate (A6); furans, oxygenated terpenoids
  (C10H16O, C10H18O…) and unsaturated aliphatic oxygenates do not, and are
  likely underestimated.
- **No SOM fits exist for N-, S-, Si- or halogen-containing species.** Their
  ranges (A8–A10) are placeholders that start at 0.
- **Each class has one surrogate**, scaled only by carbon number: all PAHs are
  naphthalene-like, all terpenes α-pinene-like, all alkenes isoprene-like.
- **Equal candidate weights** (A13): a formula shared by a common and a rare
  isomer is still split evenly unless weights are supplied.
- **Lumped or placeholder species** (UCM, "C18 aromatics") carry the carbon
  number of their placeholder formula unless the true one is supplied.

## 9. The four profiles (`run_profiles.py`)

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

**Weighted mean SOAFP, Σ(wt % × SOAFP)/Σ wt %** (VCP: total column `wt %.6`;
cooking: PMF mass fraction; species without SMILES use the formula route):

| Profile | Legacy | This module, central | VWL range | Formula only, central | Formula only, range |
|---|---|---|---|---|---|
| gasoline | 0.0284 | 0.0283 | 0.0220–0.0347 | 0.0295 | 0.0190–0.0416 |
| diesel | 0.2270 | **0.3595** | 0.3163–0.4029 | 0.2350 | 0.1640–0.3471 |
| VCP | 0.0379 | 0.0340 | 0.0222–0.0469 | 0.0344 | 0.0169–0.0606 |
| cooking | 0.0029 | 0.0035 | 0.0014–0.0059 | 0.0044 | 0.0014–0.0082 |

**Two rows dominate the differences:**
- **diesel "UCM"** (11.9 wt %): the legacy script parsed the placeholder
  formula C12H26 and got 0.069. The table's `n_C` = 26 gives 1.05, which
  explains almost all of the diesel increase. *Please confirm that 26 is meant.*
  The formula-only column cannot reach this value, because the placeholder
  formulas carry the wrong carbon numbers.
- **VCP "Witch hazel"** (0.54 wt %): `No. C` = 20 and `No. O` = 14 in the
  table, but the SMILES and formula are decane. Legacy took C20 (0.66); this
  module takes decane (0.033). Neither is a good representation of a plant
  extract.

Other changes against legacy: indane, indene and tetralin → Aromatic (A11);
terpenes → Cycloalkene (limonene, pinene, pine oil and turpentine were
Aromatic or PAH before); the diesel "Aromatic" and "Polycyclic Aromatic" lumps
use their given carbon numbers; aromatic oxygenates (benzaldehyde, benzyl
alcohol, phenols) range up to the Aromatic column; cooking species without
SMILES now get the formula central value instead of 0.

**Input rows to check** (all are in `results/*_input_issues.csv`):
- gasoline "napthalene": formula C10H8, but the SMILES is methylnaphthalene,
  which is used;
- cooking undecanal/undecenal: the formulas and SMILES are swapped or wrong;
- cooking: eight species have no SMILES (`---`), so the formula route is used;
- VCP: Witch hazel (above); Texanol `No. O` is 4 for a molecule with 3 O;
  "C9 cycloalkanes" and "Cis-bicyclo[4.3.0]nonane" carry the SMILES of other
  compounds; siloxane lumps have inconsistent `No. C`.

## 10. Legacy scripts

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

## 11. Open items

1. Document the box-model conditions (OH exposure, OA loading, NOx, T) and the
   VWL rate. Confirm the b-Alkane and Cycloalkane surrogates.
2. Confirm the diesel UCM carbon number and fix the input rows listed in §9.
   Better still, give lumped species an explicit class and carbon number
   instead of a placeholder SMILES.

## References

- Cappa, C. D. and Wilson, K. R.: Multi-generation gas-phase oxidation, equilibrium partitioning, and the formation and evolution of secondary organic aerosol, Atmos. Chem. Phys., 12, 9505–9528, 2012.
- Cappa, C. D., et al.: Application of the Statistical Oxidation Model (SOM) to SOA formation from photooxidation of C12 alkanes, Atmos. Chem. Phys., 13, 1591–1606, doi:10.5194/acp-13-1591-2013, 2013.
- Cappa, C. D., et al.: Simulating secondary organic aerosol in a regional air quality model using the statistical oxidation model – Part 2: Assessing the influence of vapor wall losses, Atmos. Chem. Phys., 16, 3041–3059, doi:10.5194/acp-16-3041-2016, 2016 (Table S1: SOM parameters).
- Zhang, X., et al.: Influence of vapor wall loss in laboratory chambers on yields of secondary organic aerosol, PNAS, 111, 5802–5807, doi:10.1073/pnas.1404727111, 2014.
- McDonald, B. C., et al.: Volatile chemical products emerging as largest petrochemical source of urban organic emissions, Science, 359, 760–764, doi:10.1126/science.aaq0524, 2018 (Table S12: SOA yields).

Original scripts: Yi Wu. Module: Tool-Chains.
