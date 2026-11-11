# SOA formation potential (SOAFP) from speciated emission profiles: a literature review of yield parameterizations, provenance and uncertainties

*Prepared for the SOAP/SOAFP tool chain (gasoline, diesel, VCP and cooking
profiles); literature checked up to November 2026. This is an AI-assisted
literature search. Numbers were read from the primary sources where possible.
Anything marked **[not verified]** comes from memory or a secondary citation
and must be checked before publication. Section 9 lists the open checks.*

---

## 0. What the current table actually contains

From the workbook `legacy/Copy of SOM_boxmodel(for Shantanu).xlsx` (sheets
`no VWL`, `VWL`, `avg`, `avg_py`):

- **Header and surrogates.** The header reads "SOM II (Updated 11/23/2016)".
  The surrogates are benzene (≤ C6), toluene (C7), m-xylene (≥ C8),
  naphthalene, isoprene and α-pinene. The oxygenated columns (O = 1–5) are
  labelled "m-undecane" (no-VWL sheet) or "dodecane" (VWL sheet). The alkane
  columns use n-dodecane.
- **The oxygenated columns are capped:** constant from C20 upward.
- **Unphysical entries.**
  - The Aromatic column has values at C5, although no C5 aromatic exists.
  - The no-VWL sheet has two unlabelled scratch columns; one goes negative
    (−0.007 at C14).
  - The O = 5 value at C5 is identical in both sheets (0.00371).
- **Mass yields above 1.** Above about C22 the alkane yields exceed 1 (up to
  1.06–1.09). That is consistent with mass yields (oxygen is added), but it
  means those species are treated as fully converted to SOA even though they
  are already semivolatile or low-volatility.
- **An embedded Zhao block.** A block "Table S7b. Zhao (n-alkanes, SOA mass
  yield)" gives log C\* for C5–C30, from 9.3 down to −0.6 (about −0.4 per
  carbon). The SOM n-alkane yields next to it are identical to the n-Alkane
  column.
- **The box-model conditions are not recorded.** The workbook gives no OH
  concentration or time, OA loading, NOx branch or k_wall. The only note is
  "correct for OH exposure at the time considered".
- **Reference values from the avg sheet** (g SOA per g reacted):

  | Species | Yield |
  |---|---|
  | n-alkane C12 / C15 / C18 / C21 | 0.069 / 0.144 / 0.34 / 0.86 |
  | toluene (C7) | 0.090 |
  | m-xylene (C8) | 0.049 |
  | α-pinene (C10) | 0.132 |
  | isoprene (C5) | 0.026 |
  | O = 3 column, C6 / C8 | 0.0034 / 0.0083 |

- **VWL / no-VWL ratio:**

  | Species | Ratio |
  |---|---|
  | C12 n-alkane | ≈ 1.2 |
  | toluene | 1.7 |
  | α-pinene | 1.85 |
  | isoprene | 4 |
  | O = 1 at C10 | ≈ 11 |

---

## 1. Summary of recommendations (ranked)

| # | Change | Expected impact | Effort | Key reference |
|---|---|---|---|---|
| 1 | **Replace the SOM-only O = 1–5 columns for acyclic OVOCs.** Use a bracket: lower bound ≈ 0 (Humes et al. 2022, OFR); central = Sasidharan et al. (2023), y = 0.002 + 0.032·(N_C − N_O) (VWL-corrected, 12 h at 1.5×10⁶ cm⁻³ OH, C_OA = 10 µg m⁻³, high NOx). The current columns are 10–25× below the yields McDonald et al. (2018) actually used. | Large for VCP SOAFP: OVOCs are about 25–40% of VCP SOA in recent work | Low | Sasidharan 2023; Humes 2022; McDonald 2018 |
| 2 | **Document and harmonize the reference condition:** SOAFP as a mass yield per g reacted at C_OA = 10 µg m⁻³, 298 K, high NOx, VWL-corrected (the convention of Gentner 2012, Lu 2020, Sasidharan 2023, CRACMM). Ask the Jathar/Cappa group for the actual "SOM II 2016" run settings. | Makes the numbers interpretable and comparable | Low | Cappa 2016; Akherati 2019 |
| 3 | **Report a per-emitted-mass SOAFP next to the per-reacted yield:** weight by the fraction reacted F = 1 − exp(−k_OH·[OH]·Δt) at a stated exposure (e.g. 12 h × 1.5×10⁶ cm⁻³). The profiles already carry k_OH. | Large: benzene and D5 drop by about 10×, C12+ alkanes and terpenes barely change | Low | Sasidharan 2023; Grosjean 1992; Derwent 2010 |
| 4 | **Treat VWL as a range, not an average:** noVWL as the lower bound, VWL as the central/upper value. In Cappa (2016), whole-model VOC-SOA rises 2–5× (low k_wall) and 5–10× (high k_wall), and the high-VWL case matched observations best. | Moderate to large; makes the range honest | Low | Cappa 2016; Zhang 2014; Krechmer 2016 |
| 5 | **Siloxanes (D4–D6): use ≈ 0.01** (range 0.005–0.05) at ambient OH, not the O = 5 column (≈ 0.17 for D5) | Removes a roughly 10× overestimate for personal-care VCPs | Low | Charan 2022; Janechek 2019 |
| 6 | **Monoterpenes, sesquiterpenes, fragrance terpenoids: compound-specific values.** Limonene is ≈ 0.30 by mass (CRACMM) versus 0.13 now. Sesquiterpenes need NOx-dependent values. Flag linalool and other oxygenated terpenoids rather than running them through the O = 1 column. | Moderate for VCP and cooking (herbs, spices) | Low | Pye 2023 (CRACMM); Klein 2016 |
| 7 | **Oxygenated aromatics, furans, benzyl alcohol: explicit values.** Phenoxyethanols ≈ 0.15 (Humes), benzyl alcohol ≈ 0.13 (Sasidharan), phenols/guaiacols per Akherati 2020, furan ≈ 0.05 (CRACMM). | Moderate for VCP, cooking and burning-like sources | Low–moderate | Akherati 2020; Pye 2023 |
| 8 | **Formula-only species: add estimates to the DBE bracket.** Add a volatility-based estimate (Li et al. 2016 formula volatility, ±0.5–1 log unit for isomers) and Sasidharan's (N_C − N_O) estimate. | Tighter, better-defined ranges for PTR/CIMS species | Moderate | Li 2016; Isaacman-VanWertz & Aumont 2021 |
| 9 | **Cross-check alkane, IVOC and aromatic yields** against Lu et al. (2020), Gentner et al. (2012) and CRACMM ROC values, and report the spread as structural uncertainty. Do not count SVOC/LVOC species (C\* ≲ 10² µg m⁻³) as SOA precursors unless POA evaporation is modeled. | Moderate for diesel | Moderate | Lu 2020; Gentner 2012 |
| 10 | **Sensitivity factors** for C_OA (1–50 µg m⁻³), temperature (ΔH_vap from Epstein et al. 2010) and NOx (low- vs high-NOx SOM sets), applied as multiplicative ranges | Moderate | Moderate | Gentner 2012; Epstein 2010; Cappa 2016 |

---

## 2. Provenance of the SOM carbon-number table

### 2.1 The SOM framework

Each species is defined by its carbon and oxygen numbers (N_C, N_O). The
tunable parameters are:
- p_1O–p_4O: the probability of adding 1–4 oxygen atoms per reaction;
- ΔLVP: the decrease in log vapor pressure per added oxygen;
- m_frag: the fragmentation exponent, with P_frag = (O:C)^m_frag.

Partitioning is absorptive (Pankow 1994), with dynamic gas–particle transfer
(Cappa & Wilson 2012; Cappa et al. 2016, Sect. 2.2).

### 2.2 How carbon-number tables arise

A surrogate's fitted parameters set the oxidation trajectory. Running the same
parameters with the precursor placed at a different N_C gives a different
yield. Sasidharan (2023, SI S1) states it directly: "the SOA mass yields for a
precursor, with a different carbon and oxygen number than the assigned
surrogate, will be different from those for the surrogate."

### 2.3 The surrogate sets match the Jathar/Cappa group's conventions

- **Cappa et al. (2016), Table S1: VWL-corrected SOM parameters.**
  - Precursors: long alkanes (ALK5, n-dodecane; data Cappa 2013 / Loza 2014),
    benzene, toluene (ARO1), m-xylene (ARO2) (Ng 2007; Zhang 2014), isoprene
    and α-pinene (TRP1/SESQ) (Chhabra 2011).
  - Separate low- and high-NOx fits, for k_wall = 1×10⁻⁴ s⁻¹ (SOM-low) and
    2.5×10⁻⁴ s⁻¹ (SOM-high).
  - The no-VWL (SOM-no, k_wall = 0) parameters are in Jathar et al. (2015,
    Part 1).
  - The fits use α_particle = 1 (conservative).
  - Footnote: "the ALK5 class is separated into long alkane species grouped
    according to carbon number". Sesquiterpenes reuse the α-pinene parameters
    and get larger yields only through their carbon number.
- **Akherati et al. (2019), Table 1** has the same structure, with VWL
  "No"/"Yes" and low/high NOx for each precursor. Example: n-dodecane, high
  NOx, no VWL has ΔLVP = 1.39 and m_frag = 0.098.
- **Sasidharan (2023, thesis/ES&T), SI S1** lists the group's surrogate
  assignments:

  | Species | Surrogate |
  |---|---|
  | linear alkanes | n-dodecane |
  | branched alkanes | methylundecane |
  | cyclic alkanes | hexylcyclohexane |
  | alkenes C5–C9 | isoprene |
  | alkenes C10+ | α-pinene |
  | single-ring aromatics C8+ | m-xylene |
  | multi-ring aromatics | naphthalene |

  These are exactly the columns of our table. The b-Alkane and Cycloalkane
  parameters come from Cappa et al. (2013) (Caltech C12 alkane experiments,
  Loza 2014).

**McDonald et al. (2018, Science)** states that "the aerosol yields used in
this study (table S12) are mostly estimated from the Statistical Oxidation
Model (SOM) … along with a one-dimensional volatility basis set for OVOCs".
Co-authors include Jathar, Akherati and Cappa. The "for Shantanu" file name
and the "SOM II 2016" date point to our table being the SOM input to that
Table S12.

### 2.4 McDonald's OVOC yields were not SOM-only

Seltzer et al. (2021, ACP, Sect. 2.5) describe McDonald's "other non-cyclic
oxygenates" as the **arithmetic average of SOM and a 1-D VBS approach**.
Humes et al. (2022) quote McDonald's estimated yields:

| Compound | McDonald 2018 | Our O column (avg) |
|---|---|---|
| methyl carbitol | 0.066 | — |
| carbitol (C6H14O3) | 0.087 | 0.0034 (O = 3, C6) |
| butyl carbitol (C8H18O3) | 0.16 | 0.0083 (O = 3, C8) |
| phenoxyethanol | 0.023 | — |
| 1-phenoxy-2-propanol | 0.036 | — |

So our O columns are probably only the SOM half of what McDonald used, and
they are 10–25× lower.

### 2.5 Zhao "Table S7b"

Zhao et al. (2015 ES&T, diesel IVOCs; 2016 ES&T, gasoline IVOCs) built a
79-species IVOC mechanism at C_OA = 10 µg m⁻³ under high NOx. Lu et al.
(2020) note that "C18−22 n-alkanes are assumed to have the same SOA yields as
C17 n-alkanes". The log C\* column in our sheet is consistent with Zhao's
volatility-binned n-alkane treatment. **[not verified]**: the ACS supplements
were not accessible, so whether "Table S7b" is from Zhao 2014, 2015 or 2016,
and what it contains, is unconfirmed.

### 2.6 Jathar et al. 2014 PNAS

**[full text not read]** As described by Akherati et al. (2019), unspeciated
IVOCs were modeled as an equivalent linear alkane (C13 for gasoline, C15 for
diesel and biomass burning) using n-dodecane SOM parameters.

### 2.7 Box-model conditions

The exact settings behind our table are **not documented in any accessible
source**. Group-standard conditions in comparable work:

| Source | [OH] | Time | C_OA | Other |
|---|---|---|---|---|
| Sasidharan (2023) | 1.5×10⁶ cm⁻³ | 12 h (some runs to 48 h) | 0.1, 1 and 10 µg m⁻³ | high-NOx parameters, VWL on |
| Akherati (2019), Fig. S5 VWL enhancement | — | — | 9 µg m⁻³ | — |
| Lu (2020), Fig. 2 | 3×10⁶ cm⁻³ | — | 5 µg m⁻³ (also 1–50) | 298 K |

Sasidharan reports toluene 0.096 at 12 h, 10 µg m⁻³, VWL on, which is close to
our VWL toluene value (0.112). That suggests similar conditions (about
10 µg m⁻³, half a day to a day of aging), but **this is an inference**. Our
m-xylene (0.035–0.06) and toluene (0.07–0.11) values are closer to high-NOx
chamber yields (Ng et al. 2007), so the table is probably high-NOx or a
low/high average **[inference]**. The yields are mass yields: some exceed 1,
and the header says so.

**Action:** ask S. Jathar / C. Cappa for [OH], duration, C_OA, NOx branch and
k_wall (1.0 or 2.5×10⁻⁴ s⁻¹).

---

## 3. Alternative parameterizations and magnitude comparison

**Gentner et al. 2012 (PNAS, SI Table S8).**
- High-NOx yields at C_OA = 10 µg m⁻³ by carbon number for 8 classes.
  n-Alkanes come from the Jordan et al. model. Branched alkanes assume 30%
  alkoxy decomposition. The other classes come from a Monte Carlo over
  surrogate scenarios.
- Lowering C_OA from 10 to 5 µg m⁻³ reduces C16–22 alkane yields by about 30%
  and C6–8 aromatic yields by about 25%.
- Bulk yields: diesel 0.15 ± 0.05, gasoline 0.023 ± 0.007.

**Lu et al. 2020 (ACP, Table 2).** A six-group VBS fit to Zhao 2015/2016
(high NOx, 298 K). Yields at 10 µg m⁻³:

| Group | Molecular weight | Yield | k_OH |
|---|---|---|---|
| IVOCP6-ALK | 184 | 0.15 | 1.55×10⁻¹¹ cm³ s⁻¹ |
| IVOCP5-ALK | 219 | 0.35 | — |
| IVOCP4-ALK | 255 | 0.43 | — |
| IVOCP3-ALK | 297 | 0.43 | — |
| IVOCP6-ARO | 162 | 0.25 | — |
| IVOCP5-ARO | 197 | 0.36 | 7.56×10⁻¹¹ cm³ s⁻¹ |

Lu et al. (2018) give the matching IVOC emission fractions: 4.6% of gasoline,
28.5% of aircraft and 55.5% of diesel NMOG.

**CRACMM v1.0 (Pye et al. 2023, ACP).** Yields assume population-weighted RO₂
fates, a 10 s bimolecular RO₂ lifetime and C_OA = 10 µg m⁻³. A SAR uses
yield-vs-log C\* fits from Seltzer et al. (2021), mostly from high-NOx data.
Single-ring aromatics are scaled for VWL following Zhang (2014). Values:

| Species | Yield |
|---|---|
| ROCP3ALK / ROCP4ALK | SAR ≈ 0.83 / 0.55 by mass |
| ROCP5ALK, ROCP6ALK, HC10 (mechanism prompt yields) | 0.28, 0.18, 0.06 by mass |
| ROCP5ARO, ROCP6ARO, NAPH | 37%, 21%, 21% molar |
| VROCIOXY (primary oxygenated IVOCs) | 14.9% molar; k_OH-limited, ≈ 40 h lifetime |
| FURAN | ≈ 5% by mass |
| API | ≈ 21% by mass typical, up to ≈ 37% at high NOx |
| LIM | ≈ 30% by mass, up to 50% |
| Sesquiterpenes, RO₂ + HO₂ | 50–91% molar prompt (C_OA 1–10 µg m⁻³) |
| Sesquiterpenes, RO₂ + NO | ≈ 12% molar |
| HC3, HC5, SLOWROC | 2.8×10⁻⁵, 1.3×10⁻³, 0.10% molar |

Because our profiles carry a CRACMM column, these are a direct cross-check.

**Seltzer et al. 2021 (ACP, VCPs).** A hybrid set:

| Species group | Yield source |
|---|---|
| linear alkanes | polynomial fit to Presto (2010) at 10 µg m⁻³ |
| cyclic alkanes | the linear alkane with 3 more carbons (Tkacik 2012) |
| branched alkanes, alkenes | SOM as reported in McDonald (2018) |
| aromatics | Ng (2007), high NOx |
| siloxanes | Janechek (2019) two-product fit |
| glycol ethers | Li & Cocker (2018) |
| other non-cyclic OVOCs | average of SOM and a 1-D VBS |
| halocarbons, species with < 5 C | 0 |
| anything else | 0 if C\* ≥ 3×10⁶ µg m⁻³, otherwise the n-dodecane yield |

The resulting sector effective yield is 5.3%, ranging from about 0% (soaps)
to 14.8% (printing inks).

**Akherati et al. 2020 (ES&T).** SOM parameters for oxygenated aromatics
(phenols, methoxyphenols) and heterocyclics (furans). Oxygenated aromatics
explain just under 60% of biomass-burning SOA and heterocyclics about 20%
**[abstract only; per-compound values not read]**.

**Other sources.**
- Pye & Pouliot (2012, ES&T) added long-alkane and PAH SOA (with oligomers)
  to CMAQ; alkanes contribute more SOA than PAHs **[yield values not read]**.
- Reviews: Charan et al. (2019, Chem Rev) on chamber simulation and wall-loss
  correction; Srivastava et al. (2022, npj Clim Atmos Sci) on anthropogenic
  yields and their condition dependence **[tables not read]**.

**Magnitude comparison** (mass yield at about 10 µg m⁻³, high NOx):

| Species | Ours (avg) | Gentner 2012 | Lu 2020 / CRACMM | Sasidharan 2023 |
|---|---|---|---|---|
| n-C13 | 0.088 | 0.019 | ≈ 0.15 (P6-ALK) | — |
| n-C15 / C16 | 0.14 / 0.19 | 0.055 / 0.089 | ≈ 0.35 (P5-ALK) | — |
| n-C18 | 0.34 | 0.23 | ≈ 0.43 (P4-ALK) | — |
| n-C21 | 0.86 | 0.77 | ≈ 0.43 (P3-ALK) | — |
| branched C12 | 0.069 | 0.0017 | — | — |
| toluene | 0.090 | — | — | 0.096 |
| C10 aromatic | 0.10 | 0.12 ± 0.08 | — | — |
| α-pinene | 0.13 | — | ≈ 0.21 (CRACMM) | — |
| carbitol C6H14O3 | 0.0034 | — | — | 0.098 (formula) |

What the comparison shows:
- Our alkane yields at C ≤ 16 lie between Gentner (lower) and Zhao/Lu
  (higher).
- Our branched-alkane yields are far above Gentner's at C ≤ 16.
- At C ≥ 20 we exceed Lu, because Zhao capped yields at C17.
- The structural spread between parameterizations is a factor of 2–5 for
  IVOC alkanes, and more than 10 for OVOCs.

---

## 4. Missing or poorly treated groups

### 4.1 Oxygenated VCPs (glycol ethers, esters, alcohols, ketones)

**Sasidharan (2023), chamber experiments with SOM-TOMAS fits** for propylene
glycol, diethylene glycol, carbitol, triethanolamine, benzyl alcohol, dimethyl
glutarate, butyl carbitol, DPGMEA, Texanol and methyl palmitate:
- Yields at 12 h and 10 µg m⁻³ are comparable to those of alkanes and
  monoterpenes.
- They fit y = 0.002 + 0.032·(N_C − N_O), with R² = 0.57 (0.88 without
  Texanol, whose yield is 0.02). Benzyl alcohol gives 0.13.
- Suggested bounds:

  | Group | Lower bound | Upper bound |
  |---|---|---|
  | glycols | diethylene glycol | propylene glycol |
  | glycol ethers | butyl carbitol | carbitol |
  | esters / acetates | Texanol | DPGMEA |

**Humes (2022), oxidation flow reactor.**
- No significant SOA from carbitol, methyl/butyl/hexyl carbitol, butyl
  butyrate or butyl acetate.
- The aromatic glycol ethers phenoxyethanol and phenoxypropanol gave ≈ 0.15
  (C_OA ≈ 10–80 µg m⁻³).

**Other chamber data.** Li & Cocker (2018) found that the position of –OH
controls glycol-ether SOA. Li et al. found 10–55% for some glycol ethers and
acetates at high loading (as summarized by Sasidharan).

**Recommendation:** for acyclic OVOCs, use the range 0 to the Sasidharan
formula. Give oxygenated aromatics explicit values.

### 4.2 Siloxanes

| Source | Yield | Conditions |
|---|---|---|
| Wu & Johnston (2017), as cited in Janechek (2019) | 0.08–0.16 | 1.2–12 µg m⁻³, OH ≈ 10⁸ cm⁻³ |
| Janechek et al. (2019), OFR | 0.22–0.50 | 68–220 µg m⁻³, OH ≈ 10¹⁰ cm⁻³ |
| Charan et al. (2022) | "largely < 5% and usually ∼1%" | OH ≲ 10⁷ cm⁻³, or exposure ≲ 2×10¹¹ molec cm⁻³ s; NOx-independent |
| Charan et al. (2022) | > 90% | exposure > 10¹² molec cm⁻³ s only |

Our O = 5 column gives ≈ 0.17 for D5, which is too high for ambient
conditions.

### 4.3 Terpenoid fragrances and sesquiterpenes

- Limonene and sesquiterpenes are underestimated by the α-pinene column (see
  the CRACMM values in Section 3).
- Klein et al. (2016): terpenes from cooking with herbs and pepper had an
  effective SOA yield of about 40% (herbs 45%, pepper 29%).
- Linalool: **no verified yield found**. Don't use the O = 1 column (≈ 0.01);
  flag it instead.

### 4.4 Cooking

- Takhar et al. (2021), heated canola oil at OH exposure 6.4×10¹⁰ cm⁻³ s:
  aldehydes make up about 84% of particle-phase products (alkadienals ≈ 46%,
  alkenals ≈ 31%, heptanal ≈ 7%). The molar yields are in the SI
  **[not read]**.
- Liu et al. (2018): alkenals explain 5–34% of observed cooking SOA; the rest
  is attributed to unidentified S/IVOCs.
- Liu et al. (2017): SOA potential scales with the content of
  monounsaturated and ω-6 fats.
- Chacon-Madrid et al. (2010): n-aldehydes fragment more than n-alkanes and
  give lower SOA **[values not read]**.
- **What this means for us:** setting < C5 species to 0 is fine (acetaldehyde,
  propanal, acetic acid). Cooking SOA is dominated by C7–C10 alkenals and
  dienals and by unmeasured S/IVOCs, so the completeness of the profile
  matters more than the yield table.

### 4.5 Furans and oxygenated aromatics

- CRACMM: FURAN ≈ 5% by mass (Bruns 2016); FURANONE 8% by mass.
- A 2026 EGUsphere preprint **[not peer reviewed]** reports, under NOx:
  furfural 9.8%, 2-methylfuran 1.4%, 3-methylfuran 1.1%.
- Phenols, guaiacols and syringol are outliers in Sasidharan's fit; syringol
  gives 0.15.

### 4.6 Halogenated and N-containing species, alkynes

- Seltzer (2021) and CRACMM set halocarbons to 0.
- Triethanolamine is parameterized by Sasidharan.
- No verified alkyne yields were found. Alkyne emissions are small and mostly
  below C5.

### 4.7 Isoprene and small (< C5) species

- The SOM isoprene yields (C5, 0.01–0.04) do not include IEPOX aqueous
  uptake. CRACMM treats IEPOX separately (Henry's law constant
  3.0×10⁷ M atm⁻¹; organosulfate rate 8.83×10⁻³ M⁻² s⁻¹).
- Glyoxal and methylglyoxal (C2–C3) form aqueous SOA that a gas-phase yield
  table cannot represent. State this as an omission.

---

## 5. Formula-only (mass-spectrometry) approaches

### 5.1 Volatility from formula (Li et al. 2016, ACP)

log10 C₀ = (n_C⁰ − n_C)·b_C − n_O·b_O − 2·(n_C·n_O/(n_C + n_O))·b_CO − n_N·b_N − n_S·b_S

Coefficients (Table 1):

| Class | n_C⁰ | b_C | b_O | b_CO | b_N | b_S |
|---|---|---|---|---|---|---|
| CH | 23.80 | 0.4861 | — | — | — | — |
| CHO | 22.66 | 0.4481 | 1.656 | −0.7790 | — | — |
| CHN | 24.59 | 0.4066 | — | — | 0.9619 | — |
| CHON | 24.13 | 0.3667 | 0.7732 | −0.07790 | 1.114 | — |
| CHOS | 24.06 | 0.3637 | 1.327 | −0.3988 | — | 0.7579 |
| CHONS | 28.50 | 0.3848 | 1.011 | 0.2921 | 1.053 | 1.316 |

Check: C12H26 gives log C₀ ≈ 5.7. Earlier formula-based approaches: Donahue
et al. (2011, 2D-VBS), Daumit et al. (2013).

### 5.2 Isomer uncertainty (Isaacman-VanWertz & Aumont 2021)

- Isomers of one formula differ in vapor pressure by half to one order of
  magnitude more than estimation methods differ from each other, and in k_OH
  by about a factor of 2.
- For isomer mixtures, formula-based estimates are unbiased, with roughly
  normal errors.
- In practice: use σ ≈ 0.5–1 log unit on log C\*, and take k_OH as roughly
  constant per element set.

### 5.3 Yield from formula

**No published, validated general formula-to-yield SAR was found.** The
closest options are:
- Sasidharan (2023), with N_C − N_O, for acyclic OVOCs;
- the CRACMM SAR (yield vs log C\*, from the Seltzer 2021 fits)
  **[coefficients not read]**;
- Akherati (2019): the VWL enhancement decreases with carbon number.

### 5.4 Proposed bounding recipe for formula-only species

1. Keep the DBE-derived candidate classes (the current approach).
2. Compute log C\* with Li et al. (2016), ± 1σ (≈ 0.5–1 log unit).
3. For each candidate, take the class yield at the carbon number of
   equivalent volatility, not at the raw N_C. This treats oxygen content
   consistently.
4. For formulas with O, add the Sasidharan value and 0 (the Humes lower bound)
   to the candidates.
5. Report the minimum, central value and maximum over the candidates and over
   C\* ± σ.
6. Report formulas with C\* ≤ 10² µg m⁻³ separately: they are SVOCs, and
   their SOA depends on partitioning.

---

## 6. Sensitivities and how to express them

| Factor | Evidence | Suggested treatment |
|---|---|---|
| Vapor wall loss | Zhang et al. (2014): chamber toluene yields underestimated by up to ≈ 4× **[factor from a secondary citation]**. Cappa (2016): whole-model VOC-SOA ×2–5 (low k_wall) and ×5–10 (high k_wall); high VWL fit observations best. La et al. (2016, GECKO-A): up to 2×. Krechmer (2016): equivalent wall mass 0.3–10 mg m⁻³. Akherati (2019): the enhancement falls with carbon number. | lower = noVWL, central = VWL; do not average |
| NOx | Cappa (2016) and Akherati (2019) give low- and high-NOx SOM sets. Aromatic yields are higher at low NOx (Ng 2007). Zhao et al. (2017): raising NMOG:NOx from 4 to 10 ppbC/ppb gave ≈ 8× higher gasoline-exhaust SOA yield. D5 is NOx-independent (Charan 2022). | state "high NOx"; low-NOx SOM runs, if available, as the upper range for aromatics and terpenes |
| OA loading | Gentner: 10 → 5 µg m⁻³ lowers yields by ≈ 25–30%. Lu (2020): 1–50 µg m⁻³. Sasidharan: 0.1, 1 and 10 µg m⁻³. | report at 10 µg m⁻³; ×0.7 for 5 µg m⁻³; stronger reduction at ≤ 1 µg m⁻³ (remote) |
| OH exposure | D5: ≈ 1% yield at ambient exposure but > 90% at 10¹² molec cm⁻³ s (Charan 2022). Aromatic yields rise with OH exposure per molecule (Chen et al. 2026). CRACMM yields decline at very long ages. | fix the exposure (e.g. 6.5×10¹⁰ molec cm⁻³ s ≈ 12 h at 1.5×10⁶ cm⁻³) and use the same exposure for the fraction reacted |
| Temperature | ΔH_vap = −11·log10 C\*₃₀₀ + 129 kJ mol⁻¹ (Epstein 2010). Lu (2020) uses 19–52 kJ mol⁻¹. With ΔH = 50 kJ mol⁻¹, a 10 K drop lowers C\* by about 2×. | state 298 K; optionally scale C\* via ΔH_vap |
| RH / aqueous phase | IEPOX and glyoxal pathways are outside SOM. Limonene SOA reportedly increases with RH **[not verified]**. | qualitative omission |

### 6.1 Weighting by reactivity

The classic SOAFP definitions include the fraction of the VOC that reacts:
- Grosjean (1992): fractional aerosol coefficient applied to emitted VOC.
- Derwent et al. (2010): SOAP from MCM trajectory runs for 113 VOCs; 15
  compounds give 97% of UK SOAP. The values are relative to toluene
  **[toluene = 100 not verified]**.
- Wu & Xie (2018): SOAFP for China, led by toluene, n-dodecane and
  m/p-xylene.
- Sasidharan and McDonald report SOA from the VOC reacted at a fixed age.

Our values are yields per mass reacted. A per-emitted-mass SOAFP is

SOAFP_emit = Y × [1 − exp(−k_OH · OH_exp)].

At an exposure of 6.5×10¹⁰ molec cm⁻³ s, the fraction reacted is (to be
checked against our k_OH column):

| Species | k_OH (cm³ s⁻¹) | Fraction reacted |
|---|---|---|
| benzene | 1.2×10⁻¹² | ≈ 0.075 |
| D5 | 1.55×10⁻¹² | ≈ 0.10 |
| toluene | 5.6×10⁻¹² | ≈ 0.31 |
| n-dodecane | 1.32×10⁻¹¹ | ≈ 0.57 |

Report both: the per-reacted yield for comparing mechanisms, and the
per-emitted value for ranking sources.

---

## 7. Developments 2020–2026

- **Updated SOM fits:**
  - Akherati (2019): S/IVOCs, VWL, NOx;
  - Akherati (2020): oxygenated aromatics and furans;
  - Sasidharan (2023): SOM-TOMAS for 10 OVOCs, with the formula relation and
    VBS-SOM parameters for CTMs;
  - SOM-TOMAS has also been used to reconcile chamber and OFR α-pinene data.
- **Mechanism-integrated SOA.** CRACMM v1.0 (Pye 2023) ties SOA yields to the
  RO₂ fate and lumps species by volatility into ROC classes. Since our
  profiles already map to CRACMM, it is the most direct upgrade path.
- **GECKO-A.** La et al. (2016) simulated 41 high-NOx alkane and alkene
  chamber experiments (C7–C17), with up to 2× VWL bias. No comprehensive
  GECKO-A yield table by carbon number was found.
- **Machine learning.**
  - Lumiaro et al. (2021): kernel ridge regression predicts p_sat and
    partitioning coefficients to within 0.3 log units of COSMOtherm.
  - tgBoost (Galeazzo & Shiraiwa): phase state and viscosity.
  - Chen et al. (2026): ML analysis of aromatic chamber data; OH exposure per
    molecule drives the yield.
  - **No validated ML model was found that predicts SOA mass yield directly
    from precursor structure.** This remains a gap.

---

## 8. Candidate parameterization upgrades

| Compound group | Current treatment | Proposed treatment | Reference | Data availability |
|---|---|---|---|---|
| Acyclic OVOCs (glycol ethers, esters, alcohols, ketones) | SOM O = 1–5 columns (e.g. carbitol 0.0034) | range 0 to 0.002 + 0.032·(N_C − N_O); compound values where available | Sasidharan 2023; Humes 2022 | formula in thesis/paper; compound values in ES&T SI |
| Oxygenated aromatics (phenoxyethanol, benzyl alcohol, phenols, guaiacols) | O column | phenoxy-alcohols ≈ 0.15; benzyl alcohol 0.13; phenols/guaiacols per Akherati 2020 | Humes 2022; Sasidharan 2023; Akherati 2020 | good for a few compounds |
| Siloxanes (D4–D6) | O = 5 column (≈ 0.17) | ≈ 0.01 (0.005–0.05) at ambient exposure | Charan 2022; Janechek 2019 | good for D5; D4/D6 by analogy |
| Monoterpenes / limonene | α-pinene SOM by N_C (0.13) | API ≈ 0.21, LIM ≈ 0.30 by mass (high NOx up to 0.37 / 0.50) | Pye 2023 | good |
| Sesquiterpenes | α-pinene SOM at C15 (≈ 0.52) | NOx-dependent: RO₂+HO₂ 50–91% molar, RO₂+NO ≈ 12% molar | Pye 2023; Cappa 2016 | moderate |
| Oxygenated terpenoids (linalool etc.) | O = 1 column (≈ 0.01) | flag; terpene value as the upper bound, O column as the lower | — | poor |
| Furans | < C5 → 0; methylfurans by class | furan 0.05; methylfurans 0.01–0.02; furfural ≈ 0.10 | Pye 2023; 2026 preprint | moderate |
| IVOC alkanes (diesel) | SOM n/b/cyclo by N_C | cross-check with the Lu 2020 VBS; report the spread | Lu 2020; Gentner 2012 | good |
| IVOC aromatics / PAH | m-xylene / naphthalene SOM | Lu 2020 P5/P6-ARO (0.36 / 0.25); CRACMM ROCP5/6ARO | Lu 2020; Pye 2023 | good |
| SVOC/LVOC (C\* ≲ 10² µg m⁻³) | yields ≥ 1 | exclude from VOC-SOA, or treat as a POA-evaporation scenario | Akherati 2019 | conceptual |
| Cooking aldehydes / alkenals / dienals | O = 1 column | Takhar 2021 SI molar yields; flag the S/IVOC gap | Takhar 2021; Liu 2018 | moderate (SI) |
| Halocarbons, N species | 0 or O column | 0 for halocarbons; amines case by case | Seltzer 2021; Sasidharan 2023 | poor |
| Isoprene / glyoxal / methylglyoxal | SOM isoprene; < C5 → 0 | keep; document the aqueous-SOA omission | Pye 2023 | — |
| Formula-only species | DBE class bracket | add Li 2016 volatility (±0.5–1 log), the Sasidharan formula and 0 to the bracket | Li 2016; Isaacman-VanWertz & Aumont 2021 | good |

---

## 9. Not verified (check before publication)

- the contents of McDonald 2018 Table S12 (the Science supplement was not
  accessible);
- the contents of Zhao "Table S7b";
- the Jathar 2014 PNAS full text;
- the Pye & Pouliot 2012 yields;
- the Srivastava 2022 tables;
- the Takhar 2021 SI molar yields;
- linalool yields;
- the Derwent toluene = 100 normalization;
- the "factor 4" attributed to Zhang 2014 (secondary citation);
- the Sasidharan formula, taken from the 2023 CSU thesis (check it against
  the ES&T version);
- DOIs marked "from memory" below.

## 10. References

- Akherati, A. et al. (2019) ACP 19, 4561–4594. https://doi.org/10.5194/acp-19-4561-2019
- Akherati, A. et al. (2020) ES&T 54, 8568–8579. https://doi.org/10.1021/acs.est.0c01345
- Cappa, C. D. & Wilson, K. R. (2012) ACP 12, 9505–9528. https://doi.org/10.5194/acp-12-9505-2012
- Cappa, C. D. et al. (2013) ACP 13, 1591–1606. https://doi.org/10.5194/acp-13-1591-2013
- Cappa, C. D. et al. (2016) ACP 16, 3041–3059. https://doi.org/10.5194/acp-16-3041-2016
- Chacon-Madrid, H. J. et al. (2010) PCCP 12, 13975. https://doi.org/10.1039/c0cp00200c
- Charan, S. M., Huang, Y., Seinfeld, J. H. (2019) Chem. Rev. 119, 11912–11944. https://doi.org/10.1021/acs.chemrev.9b00358
- Charan, S. M. et al. (2022) ACP 22, 917–928. https://doi.org/10.5194/acp-22-917-2022
- Chen, T. et al. (2026) Science Bulletin. https://doi.org/10.1016/j.scib.2026.01.016
- Daumit, K. E. et al. (2013) Faraday Discuss. 165, 181. https://doi.org/10.1039/c3fd00045a
- Derwent, R. G. et al. (2010) Sci. Total Environ. 408, 3374–3381. https://doi.org/10.1016/j.scitotenv.2010.04.013
- Donahue, N. M. et al. (2011) ACP 11, 3303–3318. https://doi.org/10.5194/acp-11-3303-2011
- Epstein, S. A., Riipinen, I., Donahue, N. M. (2010) ES&T 44, 743–748. https://doi.org/10.1021/es902497z
- Gentner, D. R. et al. (2012) PNAS 109, 18318–18323. https://doi.org/10.1073/pnas.1212272109
- Grosjean, D. (1992) Atmos. Environ. 26A, 953. https://doi.org/10.1016/0960-1686(92)90027-I
- Humes, M. B. et al. (2022) ES&T 56, 4806–4815. https://doi.org/10.1021/acs.est.1c07354
- Isaacman-VanWertz, G. & Aumont, B. (2021) ACP 21, 6541–6563. https://doi.org/10.5194/acp-21-6541-2021
- Janechek, N. J. et al. (2019) ACP 19, 1649–1664. https://doi.org/10.5194/acp-19-1649-2019
- Jathar, S. H. et al. (2014) PNAS 111, 10473–10478. https://doi.org/10.1073/pnas.1323740111
- Jathar, S. H. et al. (2015) ACP, SOM in UCD/CIT, Part 1. https://doi.org/10.5194/acp-15-9253-2015 [from memory]
- Klein, F. et al. (2016) Sci. Rep. 6, 36623. https://doi.org/10.1038/srep36623
- Krechmer, J. E. et al. (2016) ES&T 50, 5757–5765. https://doi.org/10.1021/acs.est.6b00606
- La, Y. S. et al. (2016) ACP 16, 1417–1431. https://doi.org/10.5194/acp-16-1417-2016
- Li, L. & Cocker, D. R. (2018) Atmos. Environ. https://doi.org/10.1016/j.atmosenv.2017.12.025
- Li, Y., Pöschl, U., Shiraiwa, M. (2016) ACP 16, 3327–3344. https://doi.org/10.5194/acp-16-3327-2016
- Liu, T. et al. (2017) ACP 17, 7333. https://doi.org/10.5194/acp-17-7333-2017
- Liu, T. et al. (2018) ACP 18, 11363. https://doi.org/10.5194/acp-18-11363-2018
- Loza, C. L. et al. (2014) ACP 14, 1423–1439. https://doi.org/10.5194/acp-14-1423-2014
- Lu, Q., Zhao, Y., Robinson, A. L. (2018) ACP 18, 17637–17654. https://doi.org/10.5194/acp-18-17637-2018
- Lu, Q. et al. (2020) ACP 20, 4313–4332. https://doi.org/10.5194/acp-20-4313-2020
- Lumiaro, E. et al. (2021) ACP 21, 13227. https://doi.org/10.5194/acp-21-13227-2021
- McDonald, B. C. et al. (2018) Science 359, 760–764. https://doi.org/10.1126/science.aaq0524
- Ng, N. L. et al. (2007) ACP 7, 3909–3922. https://doi.org/10.5194/acp-7-3909-2007
- Pankow, J. F. (1994) Atmos. Environ. 28, 185. https://doi.org/10.1016/1352-2310(94)90093-0
- Presto, A. A. et al. (2010) ES&T 44, 2029. https://doi.org/10.1021/es903712r
- Pye, H. O. T. & Pouliot, G. A. (2012) ES&T 46, 6041–6047. https://doi.org/10.1021/es300409w
- Pye, H. O. T. et al. (2023) ACP 23, 5043–5099 (CRACMM). https://doi.org/10.5194/acp-23-5043-2023
- Sasidharan, S. et al. (2023) ES&T 57, 11891–11902. https://doi.org/10.1021/acs.est.3c00683 ; thesis: https://mountainscholar.org/items/2729120b-69f1-4176-bf63-44ec5ed1afb3
- Seltzer, K. M. et al. (2021) ACP 21, 5079–5100. https://doi.org/10.5194/acp-21-5079-2021
- Srivastava, D. et al. (2022) npj Clim. Atmos. Sci. 5, 22. https://doi.org/10.1038/s41612-022-00238-6
- Takhar, M., Li, Y., Chan, A. W. H. (2021) ACP 21, 5137. https://doi.org/10.5194/acp-21-5137-2021
- Wu, W. & Xie, S. (2018) ES&T 52, 8146–8156. https://doi.org/10.1021/acs.est.8b01269
- Zhang, X. et al. (2014) PNAS 111, 5802. https://doi.org/10.1073/pnas.1404727111
- Zhao, Y. et al. (2014) ES&T 48, 13743. https://doi.org/10.1021/es5035188
- Zhao, Y. et al. (2015) ES&T 49, 11516. https://doi.org/10.1021/acs.est.5b02841
- Zhao, Y. et al. (2016) ES&T 50, 4554–4563. https://doi.org/10.1021/acs.est.5b06247
- Zhao, Y. et al. (2017) PNAS. https://doi.org/10.1073/pnas.1620911114
- Furan preprint (2026), EGUsphere-2026-4881, not peer reviewed. https://egusphere.copernicus.org/preprints/2026/egusphere-2026-4881/
