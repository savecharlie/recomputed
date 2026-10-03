#!/usr/bin/env python3
"""Every number NOTES.md quotes that no other script here computes.

`swear check NOTES.md` listed ten numbers as UNSOURCED and was right about all
ten: I had read them with grep and a curl, not through a recorded run.  So this
records them, by name, from the files they came from.  Names make them pinnable;
a pin resolves by NAME and so does not care about the written precision or how
crowded the ledger is, which is what the value scan cannot do.
"""
import csv
import json
from pathlib import Path

D = Path(__file__).parent / "data"
R = json.loads((Path(__file__).parent / "reresolve.json").read_text())

bxa = next(x for x in R["rows"] if x["lookup"] == "baloxavir acid")
for i, c in enumerate(bxa["live_cids"]):
    print(f"bxa_cid_{i} = {c}")
print(f"bxa_n_cids = {len(bxa['live_cids'])}")

row = next(r for r in csv.DictReader((D / "ligand_resolution.csv").open())
           if "baloxavir" in r["ligand_name"].lower())
print(f"bxa_author_cid = {row['pubchem_cids']}")

c0125 = next(r for r in csv.DictReader((D / "reproduction_local.csv").open())
             if r["claim_id"] == "C0125")
print(f"c0125_reported_score = {c0125['reported_score']}")
print(f"c0125_rerun_score = {c0125['rerun_score']}")
print(f"c0125_abs_delta = {c0125['abs_delta']}")
print(f"c0125_smiles_source_cid = {c0125['smiles_source'].split(':')[-1]}")

auth = json.loads((D / "RESULT_CHECKS.author.json").read_text())
print(f"author_resolvable = {auth['ligand_resolution']['resolvable']}")
print(f"author_unresolved = {auth['ligand_resolution']['unresolved']}")
print(f"author_n_strings = {auth['ligand_resolution']['n_unique_reported_ligands']}")
print(f"resolved_today = {R['resolved_today']}")
print(f"verdict_flips_vs_paper = {R['flips']}")
print(f"query_failures = {R['query_failures']}")
