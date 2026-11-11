# Tool-Chains

Reusable analysis building blocks for atmospheric-chemistry projects. Each
folder is self-contained, with its own README listing its method and
assumptions.

| Folder | What it does | Needs |
|---|---|---|
| [ClassEstimationDBE](ClassEstimationDBE/) | Compound classes allowed by a neutral molecular formula (element counts + DBE) | standard library |
| [SOAP](SOAP/) | SOA formation potential (SOM SOA mass yield) for a list of formulas (central value ± range, with assumptions) or from structure (SMILES) | standard library; pandas for tables; RDKit for SMILES |

The SOAP formula route builds on ClassEstimationDBE, so keep the two folders
side by side.
