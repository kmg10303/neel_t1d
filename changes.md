# Changes — code for every item in todo.md

Each section names the file, where to put the code, and why. Ordered to match
todo.md: Phase 1 (offline, no API) → Phase 2 (pipeline hardening) → Phase 3
(needs a renewed token).

---

## Phase 1 — stats engine (works right now, no API needed)

### 1. MR-Egger regression

**File:** `mr_t1d.py`. Add near `plot_mr`, before `mr_analysis`.

Orients every SNP to a positive exposure effect (the standard MR-Egger
convention — flip sign of both beta_exp and beta_out where beta_exp < 0),
then regresses outcome on exposure with an intercept. A significant
intercept is evidence of directional pleiotropy — i.e. the IVW no-pleiotropy
assumption is violated and the Egger slope should be preferred over IVW.

```python
def mr_egger(merged):
    """Weighted MR-Egger regression. A significant intercept indicates
    directional pleiotropy -- the IVW estimate may be biased."""
    d = merged.copy()
    flip = d["beta_exp"] < 0
    d.loc[flip, "beta_exp"] *= -1
    d.loc[flip, "beta_out"] *= -1

    X = sm.add_constant(d["beta_exp"])
    weights = 1 / (d["se_out"] ** 2)
    egger = sm.WLS(d["beta_out"], X, weights=weights).fit()

    intercept, slope = egger.params["const"], egger.params["beta_exp"]
    intercept_se, slope_se = egger.bse["const"], egger.bse["beta_exp"]
    intercept_p, slope_p = egger.pvalues["const"], egger.pvalues["beta_exp"]

    print("\n----- MR-Egger -----")
    print(f"  Slope (causal estimate): {slope:.4f}  SE {slope_se:.4f}  p={slope_p:.4f}")
    print(f"  Intercept (pleiotropy):  {intercept:.4f}  SE {intercept_se:.4f}  p={intercept_p:.4f}")
    if intercept_p < 0.05:
        print("  WARNING: intercept significantly different from 0 -- evidence of directional")
        print("  pleiotropy. Prefer the Egger slope over IVW.")
    else:
        print("  Intercept not significant -- no strong evidence against IVW's no-pleiotropy")
        print("  assumption from this test alone.")
    return slope, slope_se, slope_p, intercept, intercept_se, intercept_p
```

**Caveat to say out loud in the session:** with only 3 instruments, MR-Egger
has 1 degree of freedom (n − 2). It will run and print a number, but it is
underpowered to detect anything short of gross pleiotropy. Don't oversell it
with n=3 — say so in the write-up.

---

### 2. Weighted median estimator

**File:** `mr_t1d.py`. Add after `mr_egger`.

Consistent even if up to 50% of the instrument *weight* comes from invalid
(pleiotropic) SNPs — a weaker assumption than IVW's "every instrument is
valid" or Egger's InSIDE assumption. SE is bootstrapped (no closed form),
following Bowden et al. 2016.

```python
def weighted_median(merged, n_boot=1000, seed=0):
    """Weighted median of per-SNP Wald ratios, with a bootstrap SE."""
    def _wm(d):
        d = d.sort_values("wald")
        cw = d["weight"].cumsum() - 0.5 * d["weight"]
        cw /= d["weight"].sum()
        below = d[cw.values <= 0.5]
        above = d[cw.values > 0.5]
        if below.empty:
            return d["wald"].iloc[0]
        if above.empty:
            return d["wald"].iloc[-1]
        i0, i1 = below.index[-1], above.index[0]
        w0, w1 = cw.loc[i0], cw.loc[i1]
        v0, v1 = d.loc[i0, "wald"], d.loc[i1, "wald"]
        return v0 + (0.5 - w0) / (w1 - w0) * (v1 - v0)

    d = merged.copy()
    d["wald"] = d["beta_out"] / d["beta_exp"]
    d["wald_se"] = d["se_out"] / d["beta_exp"].abs()
    d["weight"] = 1 / (d["wald_se"] ** 2)

    median = _wm(d)

    rng = np.random.default_rng(seed)
    boots = [_wm(d.sample(len(d), replace=True, weights=d["weight"],
                           random_state=rng.integers(1_000_000_000)))
             for _ in range(n_boot)]
    se = float(np.std(boots, ddof=1))
    z = median / se
    p = 2 * (1 - stats.norm.cdf(abs(z)))

    print("\n----- Weighted median -----")
    print(f"  Estimate: {median:.4f}  bootstrap SE {se:.4f}  p={p:.4f}")
    return median, se, p
```

**Caveat:** with n=3, "median" is really just the middle Wald ratio — this
estimator wants more instruments to be meaningful. Still worth running and
reporting for completeness / consistency with the sensitivity-trio standard.

---

### 3. Cochran's Q heterogeneity test

**File:** `mr_t1d.py`. Add after `weighted_median`.

Tests whether the per-SNP Wald ratios agree with each other, given the
overall IVW estimate. Large Q relative to (n_snps − 1) df ⇒ heterogeneity ⇒
at least one instrument is likely invalid. This is the real-data analogue of
the Sargan test already implemented (on simulated data) in
`archived_outputs/causality.py`.

```python
def cochrans_q(merged, ivw_beta):
    """Cochran's Q heterogeneity test across per-SNP Wald ratios."""
    d = merged.copy()
    d["wald"] = d["beta_out"] / d["beta_exp"]
    d["wald_var"] = (d["se_out"] / d["beta_exp"]) ** 2
    q = float(((d["wald"] - ivw_beta) ** 2 / d["wald_var"]).sum())
    dfree = len(d) - 1
    p = float(1 - stats.chi2.cdf(q, dfree)) if dfree > 0 else float("nan")

    print("\n----- Cochran's Q (heterogeneity) -----")
    if dfree <= 0:
        print("  Not testable with fewer than 2 instruments.")
    else:
        print(f"  Q = {q:.3f} on {dfree} df, p = {p:.3f}")
        if p < 0.05:
            print("  WARNING: significant heterogeneity -- IVW's equal-effect assumption is")
            print("  violated. At least one instrument is likely pleiotropic or invalid.")
        else:
            print("  No significant heterogeneity detected.")
    return q, dfree, p
```

---

### 4. Wire Phase 1 estimators into `mr_analysis()`

**File:** `mr_t1d.py`. In `mr_analysis`, right after the existing IVW block:

```python
        print("\n----- MR Results (IVW via WLS) -----")
        print(wls.summary())
```

add:

```python
        mr_egger(merged)
        weighted_median(merged)
        cochrans_q(merged, wls.params.iloc[0])
```

Now every run prints the full IVW / Egger / weighted-median / Cochran's Q
sensitivity block — the standard trio the write-up currently lacks.

---

### 5. Restore the falsification block into `causal.py`

`causal.py` is a stripped-down copy of `archived_outputs/causality.py` — the
refutation/falsification section (negative control outcome, Sargan
overidentification, leave-one-out) got dropped. `PROJECT_STATUS.md` already
describes these as done; they currently aren't, in the file you have open.
This is the block that PROJECT_STATUS.md calls "the concrete causal AI
component the write-up has so far only gestured at" — worth having back.

**File:** `causal.py`. Add these functions after `run()` (or before it —
placement doesn't matter, just keep them top-level):

```python
def _iv_wald(df, snps, outcome="T1D", exposure="IL6"):
    """Two-stage least squares by hand, independent of DoWhy's refuter
    internals -- so these tests don't depend on what DoWhy's IV estimator
    is doing under the hood."""
    import statsmodels.api as sm
    Z = sm.add_constant(df[snps])
    first = sm.OLS(df[exposure], Z).fit()
    second = sm.OLS(df[outcome], sm.add_constant(first.fittedvalues)).fit()
W

def _negative_control_outcome(observed, full):
    """A negative control outcome driven by the same confounders/U as T1D,
    but NOT by IL-6. A valid instrument set must return ~0 here. If it
    doesn't, the SNPs reach the outcome through something other than IL-6 --
    the exclusion restriction is broken."""
    rng = np.random.default_rng(SEED + 1)
    nc = (0.01 * (full["BMI"] - 25) + 0.03 * full["SMOKING"]
          + 0.40 * full["U"] + rng.normal(0, 1.0, len(full)))
    d = observed.copy()
    d["NEG_CONTROL"] = nc
    snps = list(INSTRUMENTS["rsid"])
    est, _, _ = _iv_wald(d, snps, outcome="NEG_CONTROL")
    print("\n[negative control outcome]")
    print(f"  IV estimate on an outcome IL-6 does not cause: {est:+.4f}")
    print(f"  {'PASS' if abs(est) < 0.5 else 'FAIL'} (expect ~0; "
          f"compare to the real estimate of {IVW_BETA:.2f})")


def _overidentification_test(observed):
    """Sargan test. With 3 instruments and 1 exposure there are 2
    overidentifying restrictions -- if the SNPs disagree about the causal
    effect, at least one is invalid. Real-data analogue: cochrans_q() in
    mr_t1d.py."""
    from scipy import stats as st
    import statsmodels.api as sm
    snps = list(INSTRUMENTS["rsid"])
    _, resid, _ = _iv_wald(observed, snps)
    aux = sm.OLS(resid, sm.add_constant(observed[snps])).fit()
    n = len(observed)
    sargan = n * aux.rsquared
    dfree = len(snps) - 1
    p = 1 - st.chi2.cdf(sargan, dfree)
    print("\n[Sargan overidentification test]")
    print(f"  statistic = {sargan:.3f} on {dfree} df, p = {p:.3f}")
    print(f"  {'PASS' if p > 0.05 else 'FAIL'} -- a small p says the instruments disagree, "
          f"implicating pleiotropy or a broken exclusion restriction.")


def _leave_one_out(observed):
    """Is the estimate driven by one influential SNP?"""
    snps = list(INSTRUMENTS["rsid"])
    print("\n[leave-one-out]")
    full_est, _, _ = _iv_wald(observed, snps)
    print(f"  all {len(snps)} instruments : {full_est:+.4f}")
    for drop in snps:
        keep = [s for s in snps if s != drop]
        est, _, _ = _iv_wald(observed, keep)
        print(f"  dropping {drop:<12}: {est:+.4f}  (shift {est - full_est:+.4f})")
```

Then, at the end of `run()` in `causal.py` (after the existing IV-vs-IVW
comparison block, before the `# Real estimate` stub — which you can now
delete), add:

```python
    print("\n----- Falsification tests -----")
    _negative_control_outcome(observed, df)
    _overidentification_test(observed)
    _leave_one_out(observed)
```

`df` and `observed` are already defined earlier in `run()` — no new
variables needed.

---

## Phase 2 — harden `mr_t1d.py` so a token expiry can't cost work again

### 1. Fix the harmonization flip bug

**File:** `mr_t1d.py`, around line 123–131. Current code builds `swap`
against the pre-filter index, then applies it via `.loc` after `merged` has
already been subset with `merged[same | swap]` — this only works because no
rows get dropped by the filter. The first time a SNP fails allele matching,
`.loc[swap, ...]` will raise `IndexingError: Unalignable boolean Series`.

Replace:

```python
        same = (merged["effect_allele_exp"] == merged["effect_allele_out"]) & \
               (merged["other_allele_exp"] == merged["other_allele_out"])
        
        swap = (merged["effect_allele_exp"] == merged["other_allele_out"]) & \
               (merged["other_allele_exp"] == merged["effect_allele_out"])
        
        merged = merged[same | swap].copy()
        merged.loc[swap, "beta_out"] *= -1
```

with:

```python
        same = (merged["effect_allele_exp"] == merged["effect_allele_out"]) & \
               (merged["other_allele_exp"] == merged["other_allele_out"])

        swap = (merged["effect_allele_exp"] == merged["other_allele_out"]) & \
               (merged["other_allele_exp"] == merged["effect_allele_out"])

        merged = merged[same | swap].copy()

        # Recompute against the POST-filter index -- the mask above was built
        # before rows were dropped and is not guaranteed to align afterward.
        swap = (merged["effect_allele_exp"] == merged["other_allele_out"]) & \
               (merged["other_allele_exp"] == merged["effect_allele_out"])
        merged.loc[swap, "beta_out"] *= -1
```

---

### 2. Palindromic SNP handling

**File:** `mr_t1d.py`. Insert right after the merge (before the `same`/`swap`
harmonization block above). A/T and C/G SNPs are strand-ambiguous — `same`
and `swap` can both evaluate true for them, so they need to be resolved
(via EAF) or excluded outright. Excluding is the safe default without a
reliable EAF field from this endpoint:

```python
        def _is_palindromic(a1, a2):
            return {a1, a2} in ({"A", "T"}, {"C", "G"})

        merged["palindromic"] = merged.apply(
            lambda r: _is_palindromic(r["effect_allele_exp"], r["other_allele_exp"]), axis=1
        )
        n_palindromic = int(merged["palindromic"].sum())
        if n_palindromic:
            print(f"WARNING: excluding {n_palindromic} palindromic SNP(s) (A/T or C/G) -- "
                  f"strand cannot be resolved without a trustworthy EAF field.")
            merged = merged[~merged["palindromic"]].copy()
```

---

### 3. Remove the dead default argument

**File:** `mr_t1d.py`, line 76. `snps=HLA_SNPs` is never used — `snps` is
unconditionally overwritten inside the function (`snps =
exposure["rsid"].dropna().unique().tolist()`, line 100).

Replace:

```python
def mr_analysis(snps=HLA_SNPs):
```

with:

```python
def mr_analysis():
```

(and drop the now-unused `snps` reference in the two commented-out
docstring lines if you want it fully clean — not required).

---

### 4. Stop swallowing exceptions

**File:** `mr_t1d.py`, end of `mr_analysis`. Currently:

```python
    except Exception as e:
        print(f"An error occurred: {e}")
```

This hides the real failure behind a one-liner — exactly what made the
tophits-empty bug hard to diagnose last session. Replace with:

```python
    except Exception:
        import traceback
        traceback.print_exc()
        raise
```

You still see a clean traceback, and the failure isn't silently absorbed —
callers (or a future test) can tell the run failed instead of reading `None`
as success.

---

### 5. Batch instead of truncate

**File:** `mr_t1d.py`. `get_instruments()` currently does `return
snp_list[:63]` — a workaround for OpenGWAS's `N(id) * N(variant) <= 64`
limit that silently discards instruments past #63. Add a batching helper and
use it everywhere `gwas.associations` is called with a variant list:

```python
def _batch_associations(variants, ids, batch_size=63):
    """OpenGWAS caps N(id) * N(variant) <= 64 per call. Batch instead of
    truncating the instrument list."""
    results = []
    for i in range(0, len(variants), batch_size):
        chunk = variants[i:i + batch_size]
        res = gwas.associations(variant=chunk, id=ids)
        if res and not isinstance(res, dict):
            results.extend(res)
    return results
```

Then in `get_instruments()`, change:

```python
    print(f"Tophits returned {len(snp_list)} clumped instruments: {snp_list}")
    return snp_list[:63]
```

to:

```python
    print(f"Tophits returned {len(snp_list)} clumped instruments: {snp_list}")
    return snp_list
```

And in `mr_analysis`, change:

```python
        exposure_raw = gwas.associations(variant=instruments, id=[GWASid])
```

to:

```python
        exposure_raw = _batch_associations(instruments, [GWASid])
```

and:

```python
        outcome_raw = gwas.associations(variant=snps, id=[OUTCOME_ID])
```

to:

```python
        outcome_raw = _batch_associations(snps, [OUTCOME_ID])
```

---

### 6. Cache harmonized data to CSV

**File:** `mr_t1d.py`. Add near the top:

```python
import os
CACHE_DIR = "outputs"

def _cache_path(exposure_id, outcome_id):
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, f"harmonized_{exposure_id}_{outcome_id}.csv")
```

In `mr_analysis`, right after `num_snps = len(merged)` (i.e. as soon as
harmonization has succeeded and produced a non-empty frame), add:

```python
        cache_file = _cache_path(GWASid, OUTCOME_ID)
        merged.to_csv(cache_file, index=False)
        print(f"Cached harmonized data -> {cache_file}")
```

This is the direct fix for "no cached data" — per PROJECT_STATUS.md this has
already cost the project a full session twice (an expired token blocked
*all* work, including analysis that needed no network, because nothing had
ever been written to disk).

**Optional but recommended (5 more minutes):** fall back to the cache when
the API call fails, so an expired token during the session degrades to
"working on last session's data" instead of "blocked":

```python
        try:
            exposure_raw = _batch_associations(instruments, [GWASid])
            if not exposure_raw or isinstance(exposure_raw, dict):
                raise RuntimeError(f"API returned no usable exposure data: {exposure_raw}")
        except Exception as e:
            if os.path.exists(_cache_path(GWASid, OUTCOME_ID)):
                print(f"Live API call failed ({e}); falling back to cached harmonized data.")
                merged = pd.read_csv(_cache_path(GWASid, OUTCOME_ID))
                # skip straight to the F-stat / IVW block using this merged df
            else:
                raise
```

(Wire-up note: this changes control flow more than the other fixes — if
you're short on time in the session, do the plain write-cache version above
first and only add the read-fallback if time remains.)

---

## Phase 3 — once the token is renewed

**Your action first:** renew at https://api.opengwas.io, paste the new JWT
into `.ieugwaspy.json` (`jwt` field). Confirmed dead as of this writing —
`gwas.gwasinfo(["prot-a-1538"])` returns `{"message": "Invalid token..."}`;
decoding the current JWT shows it expired 2026-08-15.

### 1. Verify what `prot-a-1538` actually is

The write-up currently calls it "IL-6" without this having been checked.
Run once, by hand, and paste the output into the write-up / PROJECT_STATUS.md:

```python
import ieugwaspy as gwas
info = gwas.gwasinfo(["prot-a-1538"])
row = info[0]
print(row["trait"], "| n =", row["sample_size"], "| author:", row["author"], "| pmid:", row.get("pmid"))
```

### 2. Re-run tophits with LD clumping

`USE_TOPHITS = True` is already set in `mr_t1d.py`, so just re-run:

```
source venv/bin/activate
python mr_t1d.py
```

Archived outputs show 4 different rsIDs (`rs1800795`, `rs1800797`,
`rs7802307`, `rs1554606`) returning identical betas to 6 significant
figures — they are proxies in near-perfect LD, not independent instruments.
That means the current headline **3.483 / p=0.038** understates its own
uncertainty (IVW assumes independent instruments). This re-run, now going
through `tophits(clump=1, kb=10000, r2=0.001)`, is the fix — report whatever
instrument set survives clumping, not the old hardcoded list.

Because of the Phase 1 wiring, this run will now also print MR-Egger,
weighted median, and Cochran's Q on the *real* (LD-clumped) data — the
first non-simulated sensitivity analysis in the project.

### 3. Reconcile DoWhy IV vs hand-rolled 2SLS

PROJECT_STATUS.md flags a ~3% gap: DoWhy's `iv.instrumental_variable`
returns 3.254, the hand-rolled 2SLS in `causality.py`/`causal.py` returns
3.349. With `_iv_wald` now in `causal.py` (Phase 1, item 5), add a direct
comparison at the end of `run()`:

```python
    hand_2sls, _, _ = _iv_wald(observed, list(INSTRUMENTS["rsid"]))
    print(f"\n  Hand-rolled 2SLS       : {hand_2sls:.4f}")
    print(f"  DoWhy IV estimate      : {iv_effect:.4f}")
    print(f"  Difference             : {hand_2sls - iv_effect:+.4f} "
          f"({100 * (hand_2sls - iv_effect) / iv_effect:+.1f}%)")
```

Likely cause to check first: DoWhy's `iv.instrumental_variable` with 3
instruments defaults to a specific IV estimator (check whether it's doing
2SLS or a different combination rule / weighting under the hood — inspect
`iv.estimator` or the DoWhy docs for `iv.instrumental_variable`'s default
method). If it's using a different combination of the 3 instruments than
plain 2SLS, that alone would explain a few-percent gap without either
number being "wrong."
