#!/usr/bin/env python3
"""ANCHORS: re-derive the paper's published numbers from its own CSVs, with my
own code, before computing anything new.

arXiv:2609.37542v1, Giap Duc Ha, "Which reported inputs govern molecular docking
reproducibility?"  The author shipped data/, verify_results.py and
RESULT_CHECKS.json.  His verifier passes 13/13 on this machine.  That verifies
HIS arithmetic.  These anchors verify that I can reach the same numbers by a
route that shares no code with his -- which is the only thing that licenses me
to compute a NEW number from the same tables.

Every anchor compares against the value printed in his RESULT_CHECKS.json or
TABLE_SOURCE_MAP.md, which I treat as the published claim.
"""
import csv
import json
import statistics as st
from pathlib import Path

D = Path(__file__).parent / "data"
AUTH = json.loads((D / "RESULT_CHECKS.author.json").read_text())

OK = FAIL = 0


def chk(name, got, want, tol=5e-4):
    global OK, FAIL
    ok = (got == want) if isinstance(want, int) and isinstance(got, int) \
        else abs(got - want) <= tol
    print(("  PASS  " if ok else "  FAIL  ") + f"{name:<52} got = {got!r}   want = {want!r}")
    if ok:
        OK += 1
    else:
        FAIL += 1


def rows(fn):
    with (D / fn).open() as fh:
        return list(csv.DictReader(fh))


def med(xs):
    return st.median(xs)


print("ANCHOR 1  reporting-field counts over 50 audited papers")
# LONG format: one row per (paper, field).  My first pass read it as wide and
# reported 650 papers.  The file said otherwise.
aud = rows("reporting_audit.csv")
chk("audit rows", len(aud), 650)
chk("n papers", len({r["openalex_id"] for r in aud}), AUTH["reporting"]["n_papers"])
# PRESENCE IS span_verified == "yes", NOT reported == True.  Those differ:
# `reported` is what the regex/claim pass found, `span_verified` is what a
# human could point at in the text.  Using `reported` gave me 49 for
# pdb_receptor against the paper's 47 and 46 against 43 for receptor
# preparation -- i.e. a handful of papers LOOK like they report a field and
# the span could not be verified.  Table 1 counts the verified ones.
counts = {}
claimed = {}
for r in aud:
    counts.setdefault(r["field"], 0)
    claimed.setdefault(r["field"], 0)
    if r["span_verified"] == "yes":
        counts[r["field"]] += 1
    if r["reported"].strip().lower() == "true":
        claimed[r["field"]] += 1
for f, want in sorted(AUTH["reporting"]["counts"].items()):
    chk(f"reported: {f}", counts.get(f, 0), want)
extra = sorted(set(counts) - set(AUTH["reporting"]["counts"]))
print(f"          fields in the audit but NOT in the 12-field denominator: {extra}")
gap = {f: claimed[f] - counts[f] for f in counts if claimed[f] != counts[f]}
print(f"          claimed-minus-span-verified, by field: {gap}")
for f in sorted(gap):
    if f.startswith("__"):
        continue
    print(f"          claimed_{f} = {claimed[f]}   span_verified_{f} = {counts[f]}")
print(f"          total claims not span-verifiable = "
      f"{sum(v for k, v in gap.items() if not k.startswith('__'))}")

print("\nANCHOR 2  search-parameter perturbations (casf_scaleup.csv)")
sc = rows("casf_scaleup.csv")
for p, key in (("box_size", "box_size"), ("exhaustiveness", "exhaustiveness")):
    d = [float(r["abs_delta"]) for r in sc if r["parameter"] == p]
    chk(f"{p} n", len(d), AUTH["search_scaleup"][key]["n"])
    chk(f"{p} median", med(d), AUTH["search_scaleup"][key]["median"])
    chk(f"{p} max", max(d), AUTH["search_scaleup"][key]["max"], tol=5e-4)

print("\nANCHOR 3  retained decomposition (casf_decomposition.csv)")
dec = rows("casf_decomposition.csv")
for p in ("box_size", "box_center", "exhaustiveness"):
    d = [float(r["abs_delta"]) for r in dec if r["parameter"] == p]
    w = AUTH["retained_decomposition"][p]
    chk(f"{p} n", len(d), w["n"])
    chk(f"{p} median", med(d), w["median"])
    chk(f"{p} max", max(d), w["max"])

print("\nANCHOR 4  protonation: the retained inclusion rule excludes non-ionizable")
pr = [r for r in dec if r["parameter"] == "protonation"]
# the retained rule is tag == "ionizable"; the other 7 rows are tagged
# "no-change".  My first guess ("non_ionizable" not in tag) kept all 44.
ion = [float(r["abs_delta"]) for r in pr if r["tag"] == "ionizable"]
print(f"          protonation rows total = {len(pr)}, after the tag rule = {len(ion)}")
chk("protonation n (ionizable)", len(ion), AUTH["retained_decomposition"]["protonation"]["n"])
chk("protonation median", med(ion), AUTH["retained_decomposition"]["protonation"]["median"])

print("\nANCHOR 5  ligand identity, TWO numbers for one effect")
idp = AUTH["identity_pairs"]
ident = rows("identity_outcomes.csv")
base = {r["complex"]: float(r["baseline_score"]) for r in dec
        if r["parameter"] == "identity" and r["baseline_score"]}
pairs = [abs(float(r["rerun_score"]) - base[r["target_pdb"]]) for r in ident
         if r["target_pdb"] in base and r["rerun_score"]]
chk("identity contrasts n", len(pairs), idp["n"])
chk("identity contrast median", med(pairs), idp["median"])
spreads = [float(r["abs_delta"]) for r in dec if r["parameter"] == "identity"]
chk("identity spreads n", len(spreads), idp["spreads_n"])
chk("identity spread median", med(spreads), idp["median_spread"])
# All 124 rows are smiles_source=rdkit_variant with box_size_assumed and
# center_assumed EMPTY, 124 distinct ligand names over 28 complexes -- so these
# are controlled variant swaps, not re-executions.  I guessed otherwise from
# the shared column names and the data refused it.
print(f"          124 rows: smiles sources = "
      f"{sorted({r['smiles_source'].split(':')[0] for r in ident if r['target_pdb'] in base and r['rerun_score']})}, "
      f"assumed boxes = "
      f"{sorted({r['box_size_assumed'] for r in ident if r['target_pdb'] in base and r['rerun_score']})}")

print("\nANCHOR 6  Vinardo cross-docking: 85 contrasts, median 0.955")
vin = rows("vinardo_crossdock.csv")
ok = [r for r in vin if r["status"] == "ok"]
self_ = {(r["target"], r["ligand_name"]): float(r["score"])
         for r in ok if r["kind"] == "self"}
contrasts = [abs(float(r["score"]) - self_[(r["target"], r["ligand_name"])])
             for r in ok if r["kind"] == "cross"
             and (r["target"], r["ligand_name"]) in self_]
chk("vinardo contrasts n", len(contrasts), AUTH["vinardo_crossdock"]["n"])
chk("vinardo median", med(contrasts), AUTH["vinardo_crossdock"]["median"])

print("\nANCHOR 7  self-consistency: 15 ligand ranges, median 0.040")
sel = rows("self_consistency.csv")
rcol = next((c for c in sel[0] if "range" in c.lower()), None)
print(f"          self_consistency.csv range column: {rcol!r}  "
      f"(columns: {list(sel[0])})")
if rcol:
    r_ = [float(x[rcol]) for x in sel if x[rcol] not in ("", None)]
    chk("self-consistency n", len(r_), AUTH["self_consistency"]["n"])
    chk("self-consistency median", med(r_), AUTH["self_consistency"]["median_range"])
    chk("self-consistency max", max(r_), AUTH["self_consistency"]["max_range"])

print("\nANCHOR 8  re-execution within tolerance: local 29/37, cloud 27/37")
for which, fn in (("local", "reproduction_local.csv"), ("cloud", "reproduction_cloud.csv")):
    rr = rows(fn)
    a = AUTH[which]
    chk(f"{which} manifest rows", len(rr), a["manifest_rows"])
    chk(f"{which} obtained scores",
        sum(1 for r in rr if r["status"] == "docked"), a["obtained_scores"])
    # QC: a docked pose scoring worse than -2 kcal/mol is not a usable pose.
    # My first rule (a parseable within_2 column) gave 40 locally and 0 in the
    # cloud file, which has no such column -- the ruler, not the data.
    qc = [r for r in rr if r["status"] == "docked" and r["rerun_score"]
          and float(r["rerun_score"]) <= -2]
    chk(f"{which} QC n", len(qc), a["qc_n"])
    chk(f"{which} distinct papers", len({r["paper_id"] for r in qc}), a["paper_n"])
    chk(f"{which} within 2.0",
        sum(1 for r in qc if float(r["abs_delta"]) <= 2), a["within2_n"])
    chk(f"{which} median |error|",
        med([abs(float(r["abs_delta"])) for r in qc]), a["median_abs_error"], tol=5e-4)

print("\nANCHOR 9  ligand resolvability: 90 of 116 resolved")
lr = rows("ligand_resolution.csv")
lrs = AUTH["ligand_resolution"]
chk("unique reported ligand strings", len(lr), lrs["n_unique_reported_ligands"])
ccol = next((c for c in lr[0] if "class" in c.lower()), None)
if ccol:
    from collections import Counter
    cc = Counter(r[ccol] for r in lr)
    for k, want in sorted(lrs["classification_counts"].items()):
        chk(f"class: {k}", cc.get(k, 0), want)
    resolved = cc.get("explicit_identifier", 0) + cc.get("pubchem_unique", 0)
    chk("resolved", resolved, lrs["resolvable"])
    chk("unresolved/ambiguous", len(lr) - resolved, lrs["unresolved"])

print("\nANCHOR 10  cross-environment: local vs cloud on the 35 shared claims")
def qcrows(fn):
    return [r for r in rows(fn) if r["status"] == "docked" and r["rerun_score"]
            and float(r["rerun_score"]) <= -2]
L = {r["claim_id"]: r for r in qcrows("reproduction_local.csv")}
V = {r["claim_id"]: r for r in qcrows("reproduction_cloud.csv")}
common = sorted(L.keys() & V.keys())
d = [abs(float(L[c]["rerun_score"]) - float(V[c]["rerun_score"])) for c in common]
ce = AUTH["cross_environment"]
chk("shared claims n", len(d), ce["n"])
chk("median |local-cloud|", med(d), ce["median_abs_delta"])
chk("max |local-cloud|", max(d), ce["max_abs_delta"])
chk("tolerance verdict flips",
    sum((float(L[c]["abs_delta"]) <= 2) != (float(V[c]["abs_delta"]) <= 2)
        for c in common), ce["verdict_flips"])

print("\nANCHOR 11  what the package ASSERTS rather than recomputes")
# verify_results.py reads ligand_resolution_summary.json straight off disk for
# the 90/116 claim.  So the resolver output has no independent check anywhere in
# the archive.  Confirm that the per-string table at least AGREES with the
# summary it is supposed to support.
import collections
lrr = rows("ligand_resolution.csv")
cls = collections.Counter(r[ccol] for r in lrr) if (ccol := next(
    (c for c in lrr[0] if "class" in c.lower()), None)) else {}
summ = json.loads((D / "ligand_resolution_summary.json").read_text())
agree = all(cls.get(k, 0) == v for k, v in summ["classification_counts"].items())
chk("per-string table agrees with its own summary", int(agree), 1)
print("          ...but neither is an independent resolution: no query was re-run.")

print(f"\n{'='*72}\nANCHORS  pass = {OK}   fail = {FAIL}")
