# Project status — Causal AI × T1D MR

Last updated: 2026-08-01

## Where the project actually stands

**Working and trustworthy**

- OpenGWAS query → harmonization → IVW/WLS pipeline (`mr_t1d.py`). Real, runs, produces sane output.
- Instrument extraction with `tophits` + LD clumping (`clump=1, kb=10000, r2=0.001`). This was the big unblock of the last session.
- Per-instrument F-statistics with a weak-instrument flag.
- Scatter plot with SNP error bars (`plot_mr`).
- **New (this session):** a working DoWhy causal model (`causality.py`) — DAG, identification, IV estimation, refutation, and IV-specific falsification tests that actually pass.

**Headline numbers so far**

| Exposure | Instruments | Causal estimate | p | Read |
|---|---|---|---|---|
| prot-a-1538 (IL-6) | 3 | **3.483** (SE 0.702) | 0.038 | Nominally positive |
| IL-1β | 5 | −27.7 | 0.243 | Weak instruments (F=1.87); keep as contrast only |
| White blood cell count | 60 | 0.005 | 0.889 | Null — matches literature. **This is your best result.** |

## The single most important thing to understand right now

**DoWhy needs individual-level data. Two-sample MR summary statistics are not individual-level data.**

The sketch in the old `causality.py` had `treatment="beta_exp"`, `outcome="beta_out"`. Those columns
have one row per *SNP*, not per *person*. Running DoWhy on that would return a number, and the number
would be meaningless — there is no person to confound, no distribution to adjust over, no BMI or
smoking column that could exist.

So `causality.py` now builds a **synthetic cohort calibrated to the real SNP–exposure betas**, with a
known true effect set to the IVW estimate. That makes it a *method validation*, and it is labelled as
such in every run. It is not new evidence about T1D. To do DoWhy on real data you need
individual-level genotype + phenotype (UK Biobank, FinnGen individual-level) — a different data
access problem, not a coding problem.

## Current results from `causality.py`

```
Naive regression (adjusts BMI/AGE/SMOKING)  : 4.158   <- biased by unmeasured U
DoWhy IV estimate                           : 3.254   <- recovers the truth
IVW/WLS from mr_t1d.py                      : 3.483   (95% CI 2.11–4.86)
```

The IV estimate lands inside the IVW confidence interval — DoWhy's IV estimator and our summary-level
IVW regression target the same estimand. The naive regression is off by +0.68, which is the concrete
argument for why this project uses MR at all. **That contrast is the most publishable thing in the
repo.**

IV-specific falsification, all passing:

- Negative control outcome: **+0.014** (expect ~0, vs. real estimate 3.48)
- Sargan overidentification: χ²=0.001, 2 df, **p=0.999**
- Leave-one-out: estimate moves by at most 0.013 across all three SNPs

## ⚠ The IL-6 headline result is probably not as strong as it looks

The 3.483 / p=0.038 run used the hardcoded `IL6_INSTRUMENTS` list, which was picked from the
literature and **never LD-clumped**. Evidence from your own archived outputs that these instruments
are not independent:

```
output6.md   rs1800795   beta -0.015203   se 0.001863
output7.md   rs1800797   beta -0.015203   se 0.001863   <- identical
output9.md   rs7802307   beta -0.015203   se 0.001863   <- identical
output10.md  rs1554606   beta -0.015253   se 0.001862   <- ~identical
```

Four different rsIDs returning the same beta to 6 significant figures means the API is serving
proxies in near-perfect LD — they are one signal, not four. IVW assumes *independent* instruments;
correlated instruments understate the standard error, so p=0.038 is optimistic and the true evidence
is weaker than reported. With effectively 1–2 independent instruments, the IL-6 analysis is close to
a single-SNP Wald ratio.

**Action:** re-run IL-6 through the same `tophits(clump=1)` path used for the WBC analysis, and
report whatever survives clumping. Do not put 3.483 / p=0.038 in a write-up until this is done.

## Blockers

1. **OpenGWAS token is expired.** Issued 2026-06-14, expired **2026-06-28**. Today is 2026-08-01.
   Every live API call currently returns `Invalid token`. Renew at https://api.opengwas.io and paste
   into `.ieugwaspy.json`. Nobody but the account holder can do this.
2. **No cached data.** Because nothing was ever written to disk, an expired token blocks *all* work,
   including analysis that needs no network. Harmonized data should be cached to CSV on every
   successful run. This has now bitten the project at least twice.

## Known bugs / debt in `mr_t1d.py`

- **Harmonization flip bug** (line ~131): `merged.loc[swap, "beta_out"] *= -1` uses a boolean mask
  built against the *pre-filter* index, after `merged` has already been subset. It works only because
  no rows have been dropped yet; the first time a SNP fails allele matching this raises
  `IndexingError: Unalignable boolean Series`. Recompute `swap` after the filter.
- **Palindromic SNPs are not handled.** A/T and C/G SNPs are silently dropped rather than
  strand-resolved. Standard practice is to resolve by allele frequency or exclude explicitly.
- `snps=HLA_SNPs` default arg on `mr_analysis` is dead — the parameter is overwritten on line 100.
- `return snp_list[:63]` — the 63 is a workaround for the API's `N(id) * N(variant) <= 64` limit.
  Should batch instead of truncate; right now instruments are silently discarded.
- Bare `except Exception` swallows real errors into a one-line print.
- `DoWhy IV` (3.254) vs. hand-rolled 2SLS (3.349) differ by ~3%. Worth one hour to reconcile before
  either number goes in a write-up.

## Still missing for a defensible MR paper

- **MR-Egger** (pleiotropy intercept) and **weighted median** — the standard sensitivity trio with
  IVW. Currently absent. The Sargan test in `causality.py` is the same idea but on simulated data.
- **Cochran's Q** heterogeneity on the real summary data.
- **Steiger filtering** (is the SNP really exposure-first, not outcome-first?).
- Sample overlap check between exposure and outcome GWAS — matters for bias direction.
- Confirmation of what `prot-a-1538` actually measures. Run `gwas.gwasinfo(["prot-a-1538"])` once the
  token works and record the trait string. The write-up should not say "IL-6" until this is verified.

## Suggested next three sessions

1. Renew token → verify `prot-a-1538` trait → cache all harmonized data to CSV → fix the
   harmonization bug.
2. Implement MR-Egger + weighted median + Cochran's Q on the real data. This is the highest-value
   remaining work and it is well-defined.
3. Decide the honest framing of "causal AI" (see below) and write the methods section around what
   was actually built.

## On the "causal AI" framing

As of this session it is no longer a hand-wave — there is a real DAG, a real identification step, and
real falsification tests. But be precise in the write-up about what was done:

- **Honest:** "We formalized the assumptions as a DAG, verified identification, and validated the
  IVW estimator against a DoWhy IV estimator on a calibrated simulation, with negative-control and
  overidentification falsification."
- **Not honest:** "We used causal AI to discover that IL-6 causes T1D." No causal *discovery*
  algorithm has been run, and the DoWhy numbers come from simulated individuals.

If you want genuine causal discovery, `causal-learn` (PC / GES) shipped as a DoWhy dependency and is
already installed — but it needs individual-level data too. That is the same access problem as above.
