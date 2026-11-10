# SOAFP calculation (SOA formation potential)

Scripts and lookup tables used to assign a **secondary organic aerosol formation potential (SOAFP)** to each species in an emission speciation profile. Each species gets a value taken from a Statistical Oxidation Model (SOM) box-model table, according to its **carbon number** and its **compound class** (or, for oxygenated compounds, its **oxygen number**).

Author: yi.wu

---

## 1. Files

| File | Purpose |
|---|---|
| `Tables_for_SOAFP_calculation.xlsx` | Input speciation profiles. One sheet per source: `gasoline`, `diesel`, `VCPs`, `cooking`. |
| `Copy_of_SOM_boxmodel_for_Shantanu_.xlsx` | SOM box-model results. **Only the sheet `avg_py` is read by the scripts.** |
| `SOAFP_calculation_gasoline.py` | SOAFP for the `gasoline` sheet. |
| `SOAFP_calculation_diesel.py` | SOAFP for the `diesel` sheet. |
| `SOAFP_calculation_VCP.py` | SOAFP for the `VCPs` sheet (volatile chemical products). |
| `SOAFP_calculation_cooking.py` | SOAFP for the `cooking` sheet. |

Requirements: Python 3, `numpy`, `pandas`, `openpyxl`.

---

## 2. Method in brief

For every species the scripts do the following:

1. **Count atoms.** Parse the molecular formula (e.g. `C11H10`) to get the number of carbons (`n_C`) and oxygens (`n_O`).
2. **Classify the species.** Assign a `family` (`CH`, `CHO` or `Other`) from the formula. For `CH` species, assign a `subfamily` from the SMILES string.
3. **Look up SOAFP** in the SOM table (`avg_py`) at the row where `No_C == n_C`:
   - **Oxygenated species** use the column `O = 1` … `O = 5` (5 means 5 or more O atoms), based on `n_O`.
   - **Hydrocarbon species** use the column for their subfamily. This overrides the O-based lookup.
4. **Cut-off.** Species with `n_C < 5` get `SOAFP = 0`.
5. **Write** the table, with the new columns `family`, `subfamily`, `n_C`, `n_O`, `SOAFP`, to an Excel file.

### Subfamily rules (SMILES-based, `CH` species only)

| Subfamily | Rule on the `Speciate SMILES` string | SOM column |
|---|---|---|
| n-Alkane | no `=`, no `(`, no ring labels `C1`/`C2` | `n-Alkane` |
| b-Alkane (branched) | no `=`, no ring labels, but contains `(` | `b-Alkane` |
| Cycloalkane | no `=`, contains `C1`; or compound name contains "cyclohexane" | `Cycloalkane` |
| Aromatic | contains `=` and `C1`, but not `C2` | `Aromatic` |
| PAH | contains `=`, `C1` and `C2` | `PAH` |
| Alkene | contains `=`, no `C1`/`C2` | `Alkenes` |
| Cycloalkene | contains `=` and `C1` but not `C2`, and `Speciate CRACMM` is `OLT` or `OLI`; or compound name contains "cycloalkene" | `Cycloalkene` |

Rules are applied in the order above, so later rules overwrite earlier ones. Aromatic and Cycloalkene overlap in their SMILES pattern, so the CRACMM class is what separates them.

### The SOM table (`avg_py`)

- Rows: `No_C` = 5 to 30.
- Columns: `n-Alkane`, `b-Alkane`, `Cycloalkane`, `Aromatic`, `PAH`, `Alkenes`, `Cycloalkene`, `O = 1`, `O = 2`, `O = 3`, `O = 4`, `O=5`.
- It is the average of the "No VWL" and "VWL" (vapor wall loss) SOM runs. The header rows of the `avg` sheet document which surrogate represents each column, e.g. benzene / toluene / m-xylene for aromatics, naphthalene for PAH, isoprene for alkenes, α-pinene for cycloalkenes, m-undecane for the oxygenated columns.
- `avg_py` is the clean, header-free copy of the `avg` sheet made for reading into Python. If you change the SOM numbers, update `avg_py` (the sheet `exp` is an identical copy).

---

## 3. Differences between the scripts

| | gasoline | diesel | VCP | cooking |
|---|---|---|---|---|
| Input sheet | `gasoline` | `diesel` | `VCPs` | `cooking` |
| Formula column | `Species Formula` | `Species Formula` | `Species Formula` | `Formula` |
| `n_C` / `n_O` | parsed from formula | parsed from formula | taken from `No. C` / `No. O`; parsed only where they are `n/a` | parsed from formula |
| Subfamily from SMILES | yes | yes (+ name-based rules for branched alkanes and cycloalkanes) | yes | **no**: only O-count lookup |
| Special case | none | none | none | `Componud Class == "Terpene"` uses the `Cycloalkene` column |
| Output file | `WRF Speciation with SOAFP_gasoline.xlsx` | `WRF Speciation with SOAFP_diesel.xlsx` | see note 4 below | `SOAFP_cooking.xlsx` |

---

## 4. How to run

1. Put the two Excel files in a folder.
2. In each script, change the **hard-coded paths** to your own. These are the `pd.read_excel(...)` lines (profile file and SOM file) and the `to_excel(...)` line at the bottom. File names in the scripts use spaces, e.g. `Tables for SOAFP calculation.xlsx` and `Copy of SOM_boxmodel(for Shantanu).xlsx`, so either rename the files or edit the scripts to match.
3. Run the script, e.g. `python SOAFP_calculation_gasoline.py`.
4. Open the output file and check the `family`, `subfamily` and `SOAFP` columns.

Scripts are written in Spyder cell style (`#%%`), so they can also be run cell by cell.

---

## 5. Known issues and things to check

These came up when reading the scripts against the Excel files. Please review them before relying on the numbers.

1. **Column name `O=5`.** In `avg_py` the last column is named `O=5` (no spaces), but all scripts look up `'O = 5'`. Any species with 5 or more O atoms and `n_C > 4` will raise a `KeyError`. Fix by renaming the column in the Excel file or the string in the scripts.
2. **VCP script is not directly runnable as written.**
   - It reads the columns `n_C` / `n_O`, but the `VCPs` sheet has `No. C` / `No. O`.
   - After the `'n/a'` entries in `No. C` are filled in, the mask for filling `No. O` is evaluated again and no longer matches, so `No. O` stays `'n/a'`.
   - The output path points to the **gasoline** file name and would overwrite that result.
3. **Diesel script: `'ch'` vs `'CH'`.** The conditions for branched alkane (SMILES part), cycloalkane (SMILES part) use `family == 'ch'` (lower case), which never matches. For diesel, branched alkanes and cycloalkanes are therefore only identified through the `Compound` name ("Branched Alkane", "cyclohexane", "cycloalkane").
4. **Cooking script.** The `family` column is derived with a regex applied to the `SMILES` column, whereas the other scripts use the formula column. It is not used for the lookup. Species that have no oxygen and are not terpenes get `SOAFP = 0`. The column header `Componud Class` is misspelled in the Excel file, and the script matches that spelling.
5. **Species with 0 oxygen and no hydrocarbon subfamily** (e.g. halogenated or N-containing species, family `Other`) get `SOAFP = 0`, because no lookup rule applies to them.
6. **Carbon range.** The SOM table covers C5 to C30. Species with `n_C > 30` would fail the lookup (`IndexError`). Current inputs go up to C26 (diesel) and C20 (VCPs).
7. **Units / meaning of SOAFP.** *(Add here the definition and units of the SOM values, e.g. yield per mass or per mole of VOC reacted, and the OH exposure / conditions used in the box model.)*

---

## 6. Output columns added

| Column | Meaning |
|---|---|
| `family` | `CH`, `CHO` or `Other` |
| `subfamily` | `n-Alkane`, `b-Alkane`, `Cycloalkane`, `Aromatic`, `PAH`, `Alkene`, `Cycloalkene` or `Other` |
| `n_C`, `n_O` (`n_H` for gasoline/diesel) | Atom counts from the formula |
| `SOAFP` | Value from the SOM table (0 if `n_C < 5` or no rule applies) |
