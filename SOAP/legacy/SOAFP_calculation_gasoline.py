# -*- coding: utf-8 -*-
"""
SOAFP calculation - gasoline
Created on Tue Apr 30 15:13:50 2024
@author: yi.wu

Assigns a SOA formation potential (SOAFP) to every species of the "gasoline"
speciation profile, using the SOM box-model table (sheet 'avg_py').

Logic
-----
1. n_C / n_O are parsed from the column "Species Formula".
2. family / subfamily are derived from the formula and the SMILES string.
3. SOAFP is read from the SOM table at the row No_C == n_C:
     - column 'O = 1' ... 'O = 5' (5 = five or more O) for oxygenated species,
     - column of the subfamily (n-Alkane, b-Alkane, ...) (ignored for n_C < 5).
   The subfamily lookup overrides the O-based lookup.
4. Species with n_C < 5 get SOAFP = 0.

Input files are expected in the same folder as this script (change the
settings below if they are somewhere else).
"""
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# %% ---------------------------------------------------------------- settings
try:
    DATA_DIR = Path(__file__).resolve().parent
except NameError:                       # e.g. running cell by cell in Spyder
    DATA_DIR = Path.cwd()

PROFILE_FILE = DATA_DIR / 'Tables for SOAFP calculation.xlsx'
SOM_FILE = DATA_DIR / 'Copy of SOM_boxmodel(for Shantanu).xlsx'
OUTPUT_FILE = DATA_DIR / 'WRF Speciation with SOAFP_gasoline.xlsx'
PROFILE_SHEET = 'gasoline'
SOM_SHEET = 'avg_py'
MIN_CARBON = 5                          # species with fewer C get SOAFP = 0

# %% ----------------------------------------------------------------- helpers
def parse_formula(formula):
    """Return {element: count} for a formula string such as 'C11H10' or 'C2H3Cl'."""
    counts = {}
    for element, number in re.findall(r'([A-Z][a-z]?)(\d*)', str(formula)):
        counts[element] = counts.get(element, 0) + (int(number) if number else 1)
    return counts


def contains(series, pattern, case=True):
    """Boolean str.contains that treats missing values as False."""
    return series.astype(str).str.contains(pattern, case=case, regex=True, na=False)


def lookup_soafp(SOM, n_C, column):
    """SOM value for a given carbon number and column; 0 (with a warning) if not tabulated."""
    values = SOM.loc[SOM['No_C'] == n_C, column].values
    if len(values) == 0 or pd.isna(values[0]):
        warnings.warn('No SOM value for n_C = {}, column "{}": SOAFP set to 0'.format(n_C, column))
        return 0.0
    return values[0]


# subfamily -> column name in the SOM table
SUBFAMILY_COLUMN = {
    'n-Alkane': 'n-Alkane',
    'b-Alkane': 'b-Alkane',
    'Cycloalkane': 'Cycloalkane',
    'Aromatic': 'Aromatic',
    'PAH': 'PAH',
    'Alkene': 'Alkenes',
    'Cycloalkene': 'Cycloalkene',
}

# %% ------------------------------------------------------------- read input
masslist = pd.read_excel(PROFILE_FILE, sheet_name=PROFILE_SHEET)
nVOC = len(masslist)

SOM = pd.read_excel(SOM_FILE, sheet_name=SOM_SHEET)
SOM = SOM.rename(columns={'O=5': 'O = 5'})      # header is 'O=5' in the Excel file

# %% ---------------------------------------------------------- atom numbers
n_C = np.zeros(nVOC)
n_H = np.zeros(nVOC)
n_O = np.zeros(nVOC)
for n in range(nVOC):
    if pd.notna(masslist.loc[n, 'Species Formula']):
        atoms = parse_formula(masslist.loc[n, 'Species Formula'])
        n_C[n] = atoms.get('C', 0)
        n_H[n] = atoms.get('H', 0)
        n_O[n] = atoms.get('O', 0)
masslist['n_C'] = n_C
masslist['n_H'] = n_H
masslist['n_O'] = n_O

# %% ---------------------------------------------------- family and subfamily
masslist['family'] = 'Other'
masslist.loc[contains(masslist['Species Formula'], r'^C[0-9]*H[0-9]*$'), 'family'] = 'CH'
masslist.loc[contains(masslist['Species Formula'], r'^C[0-9]*H[0-9]*O[0-9]*$'), 'family'] = 'CHO'

smiles = masslist['Speciate SMILES']
name = masslist['Compound']
is_CH = masslist['family'] == 'CH'
has_double = contains(smiles, '=')
has_ring1 = contains(smiles, 'C1')
has_ring2 = contains(smiles, 'C2')
has_branch = contains(smiles, r'\(')

masslist['subfamily'] = 'Other'
# n-alkane: no double bond, no branch, no ring
masslist.loc[is_CH & ~has_double & ~has_branch & ~has_ring1 & ~has_ring2, 'subfamily'] = 'n-Alkane'
# branched alkane
masslist.loc[(is_CH & ~has_double & ~has_ring1 & ~has_ring2 & has_branch)
             | contains(name, 'Branched Alkane', case=False), 'subfamily'] = 'b-Alkane'
# cycloalkane
masslist.loc[(is_CH & ~has_double & has_ring1)
             | contains(name, 'cyclohexane|cycloalkane', case=False), 'subfamily'] = 'Cycloalkane'
# aromatic (one ring) and PAH (two rings)
masslist.loc[is_CH & has_double & has_ring1 & ~has_ring2, 'subfamily'] = 'Aromatic'
masslist.loc[is_CH & has_double & has_ring1 & has_ring2, 'subfamily'] = 'PAH'
# alkene
masslist.loc[is_CH & has_double & ~has_ring1 & ~has_ring2, 'subfamily'] = 'Alkene'
# cycloalkene: one ring + double bond and CRACMM class OLT/OLI (or named cycloalkene)
is_olefin_cracmm = masslist['Speciate CRACMM'].isin(['OLT', 'OLI'])
masslist.loc[(is_CH & has_double & has_ring1 & ~has_ring2 & is_olefin_cracmm)
             | (is_CH & contains(name, 'cycloalkene', case=False)), 'subfamily'] = 'Cycloalkene'

# %% ------------------------------------------------ SOA formation potential
SOAFP = np.zeros(nVOC)
for i in range(nVOC):
    n_C = masslist.at[i, 'n_C']
    n_O = masslist.at[i, 'n_O']
    if n_C < MIN_CARBON:
        continue
    if n_O >= 1:
        SOAFP[i] = lookup_soafp(SOM, n_C, 'O = {}'.format(int(min(n_O, 5))))
    column = SUBFAMILY_COLUMN.get(masslist.at[i, 'subfamily'])
    if column is not None:                      # subfamily overrides the O-based value
        SOAFP[i] = lookup_soafp(SOM, n_C, column)
masslist['SOAFP'] = SOAFP

print(masslist['subfamily'].value_counts().to_string())
print('{} of {} species have SOAFP = 0'.format(int((masslist['SOAFP'] == 0).sum()), nVOC))

# %% ------------------------------------------------------------------- save
masslist.to_excel(OUTPUT_FILE, index=False)
print('Saved', OUTPUT_FILE)
