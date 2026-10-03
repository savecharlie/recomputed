#!/usr/bin/env python3
"""Re-resolve the paper's 116 reported ligand strings against LIVE PubChem.

This is the one quantitative claim in the author's archive with no
recomputation behind it anywhere.  `verify_results.py` does:

    resolver = json.loads((DATA/'ligand_resolution_summary.json').read_text())
    result['ligand_resolution'] = resolver

-- it reads the answer off disk.  The source map says why, and says it plainly:
"This is the retained deterministic resolver output; online chemistry queries
were not repeated."  So 90 resolved / 26 unresolved of 116 is an assertion
supported by a per-string table, and the table and the summary are consistent
with each other, and neither is a resolution.

So I resolve them.  PubChem PUG-REST, name -> CIDs, throttled under the 5 req/s
courtesy limit, with the author's own classification rules applied to whatever
comes back today:

    explicit_identifier        a non-empty reported_identifier -> resolvable
    pubchem_unique             exactly one CID                 -> resolvable
    pubchem_ambiguous          more than one CID               -> not
    pubchem_not_found          zero CIDs                       -> not
    descriptive_or_series_label  a label like "11a"            -> not

AND THE POINT THAT IS NOT A GOTCHA: a resolver against a live database is not
deterministic across time.  PubChem gains synonyms.  So a disagreement between
today and the author's run is not necessarily his error -- it may be the
database, which is itself a reproducibility fact about name-based ligand
identity, and a sharper one than the headline rate.
"""
import csv
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

D = Path(__file__).parent / "data"
OUT = Path(__file__).parent / "reresolve.json"
BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/JSON"
UA = "iris-maker/1.0 (reproduction check; savecharlie@gmail.com)"
SLEEP = 0.25          # under PubChem's 5 req/s courtesy limit


def cids(name: str, tries: int = 3):
    """CIDs for a name today.  Returns (list_or_None, note)."""
    url = BASE.format(urllib.parse.quote(name, safe=""))
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=30) as fh:
                js = json.loads(fh.read().decode())
            return js.get("IdentifierList", {}).get("CID", []), "ok"
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return [], "ok"          # a real "not found"
            if e.code in (429, 503):
                time.sleep(2 + 3 * a)
                continue
            return None, f"http {e.code}"
        except Exception as e:            # noqa: BLE001
            if a == tries - 1:
                return None, type(e).__name__
            time.sleep(1 + 2 * a)
    return None, "retries exhausted"


def classify(row, live):
    if (row["reported_identifier"] or "").strip():
        return "explicit_identifier", True
    if row["classification"] == "descriptive_or_series_label":
        return "descriptive_or_series_label", False
    if live is None:
        return "query_failed", None
    if len(live) == 1:
        return "pubchem_unique", True
    if len(live) > 1:
        return "pubchem_ambiguous", False
    return "pubchem_not_found", False


rowsx = list(csv.DictReader((D / "ligand_resolution.csv").open()))
summ = json.loads((D / "ligand_resolution_summary.json").read_text())
print(f"strings to resolve: {len(rowsx)}")
print(f"author's summary:   {summ['resolvable']} resolvable / "
      f"{summ['unresolved']} not, of {summ['n_unique_reported_ligands']}")
print(f"author's classes:   {summ['classification_counts']}\n")

res = []
t0 = time.time()
for i, r in enumerate(rowsx, 1):
    lookup = (r["lookup_name"] or r["ligand_name"]).strip()
    need_query = not (r["reported_identifier"] or "").strip() and \
        r["classification"] != "descriptive_or_series_label"
    live, note = (cids(lookup) if need_query else ([], "not queried"))
    if need_query:
        time.sleep(SLEEP)
    cls, ok = classify(r, live if need_query else [])
    was_ok = r["resolvable"].strip().lower() == "true"
    res.append({"ligand_name": r["ligand_name"], "lookup": lookup,
                "author_class": r["classification"], "author_resolvable": was_ok,
                "live_cids": (live if live is None else live[:8]),
                "n_live_cids": (None if live is None else len(live)),
                "live_class": cls, "live_resolvable": ok, "note": note})
    if i % 20 == 0:
        print(f"  ...{i}/{len(rowsx)}  ({time.time()-t0:.0f}s)", flush=True)

fails = [x for x in res if x["note"] not in ("ok", "not queried")]
print(f"\nqueries that failed outright: {len(fails)}")
for f in fails[:8]:
    print(f"   {f['lookup']!r}: {f['note']}")

usable = [x for x in res if x["live_resolvable"] is not None]
print(f"\nstrings with a usable verdict today: {len(usable)} of {len(res)}")
live_ok = sum(1 for x in usable if x["live_resolvable"])
print(f"RESOLVED TODAY        {live_ok} / {len(usable)}  "
      f"({100*live_ok/len(usable):.1f}%)")
auth_ok = sum(1 for x in usable if x["author_resolvable"])
print(f"RESOLVED BY THE PAPER {auth_ok} / {len(usable)}  "
      f"({100*auth_ok/len(usable):.1f}%)  (same subset)")

print(f"\nclass counts today: {dict(Counter(x['live_class'] for x in usable))}")
print(f"class counts paper: {dict(Counter(x['author_class'] for x in usable))}")

flips = [x for x in usable if x["live_resolvable"] != x["author_resolvable"]]
print(f"\nverdict changes: {len(flips)}")
for x in flips:
    print(f"   {x['lookup']!r:<46} paper={x['author_class']:<28}"
          f" today={x['live_class']} (n={x['n_live_cids']})")

json.dump({"n": len(res), "usable": len(usable), "resolved_today": live_ok,
           "resolved_paper_same_subset": auth_ok, "flips": len(flips),
           "query_failures": len(fails), "rows": res},
          OUT.open("w"), indent=1)
print(f"\nwrote {OUT.name}")
