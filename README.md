# Tool-Chains

Reusable analysis building blocks for atmospheric-chemistry projects. Each
folder is self-contained, with its own README listing its method and
assumptions.

| Folder | What it does | Needs |
|---|---|---|
| [ClassEstimationDBE](ClassEstimationDBE/) | Compound classes allowed by a neutral molecular formula (element counts + DBE) | standard library |
| [SOAP](SOAP/) | SOA formation potential (SOM SOA mass yield) from structure (SMILES) or, as a range, from formula alone | standard library; RDKit for SMILES; pandas for tables |

The SOAP formula route builds on ClassEstimationDBE, so keep the two folders
side by side.
