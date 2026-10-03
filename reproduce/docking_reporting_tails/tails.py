#!/usr/bin/env python3
"""The paper ranks reporting fields by MEDIAN effect.  Reproducibility is
decided in the TAIL.  So: rank them again by tail.

arXiv:2609.37542v1 concludes that search-related fields (box centre, box size,
exhaustiveness, random seed) are the ones to deprioritise in reporting, because
perturbing them moves the docking score by a median of no more than
0.08 kcal/mol, while receptor structure moves it by 1.02.

The same paper judges a re-execution successful when it lands within
2.0 kcal/mol of the published score, and measures its own floors: 15 ligands
re-docked with different random seeds spread by a median of 0.04 (max 0.38),
and the same 35 claims run on local and cloud hardware differ by a median of
0.18 (max 0.80).

A median says what happens on a typical day.  A reporting checklist is not for
typical days -- it exists for the case where somebody cannot reproduce you.  So
the quantity that matters is how often a factor pushes a result past the bar
the paper itself uses, and that is a tail probability, not a median.

Everything here is computed from the author's redistributed CSVs.  anchors.py
reproduces all 65 of his published numbers from those same files with code that
shares nothing with his, which is what licenses these.
"""
import csv
import json
import random
import statistics as st
from pathlib import Path

D = Path(__file__).parent / "data"
random.seed(20261003)


def rows(fn):
    with (D / fn).open() as fh:
        return list(csv.DictReader(fh))


def boot_ci(xs, fn=st.median, n=4000, lo=2.5, hi=97.5):
    k = len(xs)
    reps = sorted(fn([xs[random.randrange(k)] for _ in range(k)]) for _ in range(n))
    return reps[int(lo / 100 * n)], reps[int(hi / 100 * n) - 1]


# ------------------------------------------------- assemble every factor
dec = rows("casf_decomposition.csv")
sc = rows("casf_scaleup.csv")
vin = [r for r in rows("vinardo_crossdock.csv") if r["status"] == "ok"]
selfc = rows("self_consistency.csv")

ident_base = {r["complex"]: float(r["baseline_score"]) for r in dec
              if r["parameter"] == "identity" and r["baseline_score"]}
ident_pairs = [abs(float(r["rerun_score"]) - ident_base[r["target_pdb"]])
               for r in rows("identity_outcomes.csv")
               if r["target_pdb"] in ident_base and r["rerun_score"]]

vself = {(r["target"], r["ligand_name"]): float(r["score"])
         for r in vin if r["kind"] == "self"}
receptor = [abs(float(r["score"]) - vself[(r["target"], r["ligand_name"])])
            for r in vin if r["kind"] == "cross"
            and (r["target"], r["ligand_name"]) in vself]

FACTORS = {
    "random seed":        [float(r["range"]) for r in selfc],
    "exhaustiveness":     [float(r["abs_delta"]) for r in sc
                           if r["parameter"] == "exhaustiveness"],
    "box size":           [float(r["abs_delta"]) for r in sc
                           if r["parameter"] == "box_size"],
    "box centre":         [float(r["abs_delta"]) for r in dec
                           if r["parameter"] == "box_center"],
    "ligand protonation": [float(r["abs_delta"]) for r in dec
                           if r["parameter"] == "protonation"
                           and r["tag"] == "ionizable"],
    "ligand identity":    ident_pairs,
    "ligand identity (within-complex spread)":
                          [float(r["abs_delta"]) for r in dec
                           if r["parameter"] == "identity"],
    "receptor structure (Vinardo)": receptor,
}

# the paper's own reference levels, in kcal/mol
TOL = 2.0          # its re-execution success criterion
SEED_MAX = 0.38    # the largest spread it measured from random seed alone
HW_MAX = 0.803     # the largest local-vs-cloud difference it measured

print("TABLE A  the same factors, ranked by MEDIAN -- the paper's ordering")
print(f"  {'factor':<40} {'n':>4} {'median':>8}  {'95% CI':>18}")
for name, xs in sorted(FACTORS.items(), key=lambda kv: st.median(kv[1])):
    lo, hi = boot_ci(xs)
    print(f"  {name:<40} {len(xs):>4} {st.median(xs):>8.3f}  [{lo:.3f}, {hi:.3f}]")

print("\nTABLE B  the same factors, ranked by TAIL -- how often they cross a bar")
print(f"  {'factor':<40} {'n':>4} {'p90':>7} {'p95':>7} {'max':>7}"
      f" {'>0.38':>7} {'>0.80':>7} {'>2.0':>7}")
tail = {}
for name, xs in FACTORS.items():
    s = sorted(xs)
    def q(p):
        return s[min(len(s) - 1, int(p * len(s)))]
    frac = lambda t: sum(1 for x in xs if x > t) / len(xs)
    tail[name] = (frac(SEED_MAX), frac(HW_MAX), frac(TOL))
    print(f"  {name:<40} {len(xs):>4} {q(.90):>7.3f} {q(.95):>7.3f} "
          f"{max(xs):>7.3f} {frac(SEED_MAX):>6.1%} {frac(HW_MAX):>6.1%} "
          f"{frac(TOL):>6.1%}")

print("\nTABLE C  the ordering changes.  Rank by median vs rank by P(>0.80):")
by_med = [n for n, _ in sorted(FACTORS.items(), key=lambda kv: st.median(kv[1]))]
by_tail = [n for n in sorted(FACTORS, key=lambda n: tail[n][1])]
print(f"  {'by median':<42} {'by P(exceeds hardware noise)':<42}")
for a, b in zip(by_med, by_tail):
    mark = "   " if a == b else "  *"
    print(f" {mark}{a:<41} {b:<42}")

print("\nKEY NUMBERS")
bc = FACTORS["box centre"]
pro = FACTORS["ligand protonation"]
print(f"  box_centre        median = {st.median(bc):.3f}   "
      f"P(>0.80) = {tail['box centre'][1]:.1%}   "
      f"P(>2.0) = {tail['box centre'][2]:.1%}   max = {max(bc):.3f}")
print(f"  protonation       median = {st.median(pro):.3f}   "
      f"P(>0.80) = {tail['ligand protonation'][1]:.1%}   "
      f"P(>2.0) = {tail['ligand protonation'][2]:.1%}   max = {max(pro):.3f}")
print(f"  ratio of medians  protonation/box_centre = {st.median(pro)/st.median(bc):.2f}")
print(f"  ratio of maxima   box_centre/protonation = {max(bc)/max(pro):.2f}")

# is the box-centre tail distinguishable from the protonation tail?
from scipy import stats as sps
u = sps.mannwhitneyu(bc, pro, alternative="two-sided")
print(f"\n  Mann-Whitney box_centre vs protonation: U = {u.statistic:.1f}, "
      f"p = {u.pvalue:.4g}   (medians differ: protonation is LARGER)")
n_above = sum(1 for x in bc if x > HW_MAX)
print(f"  box_centre values above the hardware-noise maximum: {n_above} of {len(bc)}")
print(f"  box_centre values above the 2.0 tolerance:          "
      f"{sum(1 for x in bc if x > TOL)} of {len(bc)}")
print(f"  protonation values above the 2.0 tolerance:         "
      f"{sum(1 for x in pro if x > TOL)} of {len(pro)}")

print("\nREPORTING RATE vs TAIL RISK  (Table 1 span-verified counts over 50 papers)")
aud = rows("reporting_audit.csv")
papers = {r["openalex_id"] for r in aud}
ver = {}
for r in aud:
    ver.setdefault(r["field"], 0)
    if r["span_verified"] == "yes":
        ver[r["field"]] += 1
MAP = {"random seed": "random_seed", "exhaustiveness": "exhaustiveness",
       "box size": "grid_size", "box centre": "grid_center",
       "ligand protonation": "ligand_preparation",
       "receptor structure (Vinardo)": "pdb_receptor"}
print(f"  {'factor':<32} {'field':<22} {'reported':>9} {'P(>0.80)':>9} {'P(>2.0)':>8}")
for name, field in MAP.items():
    print(f"  {name:<32} {field:<22} {ver.get(field,0):>4}/{len(papers)} "
          f"{tail[name][1]:>8.1%} {tail[name][2]:>7.1%}")

json.dump({k: {"n": len(v), "median": st.median(v), "max": max(v),
               "p_gt_0.38": tail[k][0], "p_gt_0.80": tail[k][1],
               "p_gt_2.0": tail[k][2]} for k, v in FACTORS.items()},
          open(Path(__file__).parent / "tails.json", "w"), indent=1)
print("\nwrote tails.json")


# ---------------------------------------------------------------- honesty
print("\n" + "=" * 72)
print("TABLE D  the tail probabilities WITH intervals, and the clustering")
print("A tail estimated from 37 rows is not a tail.  And these rows are not")
print("independent: each factor is perturbed within a small set of complexes,")
print("so the effective sample size is the number of COMPLEXES, not of rows.")

# The right cluster variable differs by factor, and I had to look rather than
# assume.  In casf_decomposition and casf_scaleup, `complex` is a DISTINCT PDB
# per row -- 105 box-centre rows are 105 different complexes, each perturbed
# once -- so those are independent and my first worry was unfounded.  The
# Vinardo receptor table is the opposite: 85 contrasts over THREE protein
# families (THR, TRYP, CAII), 17 ligands, 18 PDB entries.  The author says in
# the source map that "shared ligands/receptors limit pooled inference"; the
# number is three.  And the random-seed floor is 15 ranges over TWO PDBs.
CLUSTER = {
    "box centre": [r["complex"] for r in dec if r["parameter"] == "box_center"],
    "box size": [r["complex"] for r in sc if r["parameter"] == "box_size"],
    "exhaustiveness": [r["complex"] for r in sc if r["parameter"] == "exhaustiveness"],
    "ligand protonation": [r["complex"] for r in dec
                           if r["parameter"] == "protonation" and r["tag"] == "ionizable"],
    "ligand identity": [r["target_pdb"] for r in rows("identity_outcomes.csv")
                        if r["target_pdb"] in ident_base and r["rerun_score"]],
    "ligand identity (within-complex spread)":
        [r["complex"] for r in dec if r["parameter"] == "identity"],
    "receptor structure (Vinardo)": [r["target"] for r in vin if r["kind"] == "cross"
                                     and (r["target"], r["ligand_name"]) in vself],
    # ^ `target` is the protein FAMILY label, which is the honest cluster here
    "random seed": [r["pdb"] for r in selfc],
}

print(f"\n  {'factor':<40} {'rows':>5} {'clusters':>9} "
      f"{'P(>2.0)':>9} {'95% CI (Clopper-Pearson)':>26}")
for name, xs in FACTORS.items():
    k = sum(1 for x in xs if x > TOL)
    n = len(xs)
    lo, hi = sps.beta.ppf(0.025, k, n - k + 1) if k else 0.0, \
             sps.beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    nclust = len(set(CLUSTER[name]))
    print(f"  {name:<40} {n:>5} {nclust:>9} {k/n:>8.1%} "
          f"  [{lo:.1%}, {hi:.1%}]")

print("\n  So: 'protonation never exceeds 2.0' rests on 0 of 37, whose upper")
print("  95% bound is the number printed above -- it is not a demonstration")
print("  that protonation is safe, only that 37 rows did not show it.")

# cluster-aware version of the box-centre tail: one value per complex
import collections
bycx = collections.defaultdict(list)
for r in dec:
    if r["parameter"] == "box_center":
        bycx[r["complex"]].append(float(r["abs_delta"]))
worst = {c: max(v) for c, v in bycx.items()}
print(f"\n  box centre: each row is already a different complex, so the collapse is")
print(f"    complexes = {len(worst)}, median of worst = {st.median(list(worst.values())):.3f}")
print(f"    complexes whose worst box-centre shift exceeds 0.80 = "
      f"{sum(1 for v in worst.values() if v > HW_MAX)} of {len(worst)}")
print(f"    complexes whose worst box-centre shift exceeds 2.00 = "
      f"{sum(1 for v in worst.values() if v > TOL)} of {len(worst)}")
prob = collections.defaultdict(list)
for r in dec:
    if r["parameter"] == "protonation" and r["tag"] == "ionizable":
        prob[r["complex"]].append(float(r["abs_delta"]))
pworst = {c: max(v) for c, v in prob.items()}
print(f"  protonation, same collapse:")
print(f"    complexes = {len(pworst)}, median of worst = "
      f"{st.median(list(pworst.values())):.3f}")
print(f"    complexes whose worst protonation shift exceeds 0.80 = "
      f"{sum(1 for v in pworst.values() if v > HW_MAX)} of {len(pworst)}")

print("\n  WHERE THE EFFECTIVE SAMPLE SIZE REALLY BITES")
print(f"    receptor structure: {len(receptor)} contrasts over "
      f"{len({r['target'] for r in vin if r['kind']=='cross'})} protein families "
      f"({sorted({r['target'] for r in vin})}), "
      f"{len({r['ligand_name'] for r in vin})} ligands, "
      f"{len({r['target_pdb'] for r in vin})} PDB entries.")
print(f"    random seed:        {len(FACTORS['random seed'])} ranges over "
      f"{len({r['pdb'] for r in selfc})} PDBs "
      f"({sorted({r['pdb'] for r in selfc})}).")
print("    So of the two floors I used as reference levels, the hardware one")
print("    (35 shared claims, max 0.80) is the better grounded and the seed one")
print("    (2 proteins) should not carry weight.  The ordering below uses the")
print("    hardware floor for that reason.")

print("\n  THE CLAIM I AM WILLING TO MAKE, and nothing wider:")
print("  Ranked by median, box centre sits BELOW ligand protonation and the")
print("  paper deprioritises it.  Ranked by the chance of a shift larger than")
print("  the paper's own hardware-noise ceiling, box centre sits ABOVE it.")
print("  Both orderings are computed from the author's own redistributed data.")
print("  Which one a reporting checklist should follow is a question about what")
print("  a checklist is for, and the paper does not ask it.")
