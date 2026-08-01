"""
DoWhy causal model for the IL-6 -> T1D question.

IMPORTANT — read this before interpreting any number this script prints.

DoWhy operates on INDIVIDUAL-LEVEL data: one row per person, with columns for
genotype, exposure, outcome and confounders. What the MR pipeline in mr_t1d.py
has is two-sample SUMMARY statistics: one row per SNP, with beta_exp / beta_out.
Those are not the same thing. Handing DoWhy a dataframe whose rows are SNPs and
calling beta_exp the "treatment" would produce a number, but the number would be
meaningless -- there is no person, no confounder, no distribution to adjust over.

So this script builds a synthetic cohort that is CALIBRATED to the real summary
statistics (real SNP-exposure betas and standard errors from prot-a-1538) and
that has a KNOWN true causal effect, set to the IVW estimate from mr_t1d.py.
That makes this a method-validation study, not new evidence about T1D:

  - It verifies the DAG we are claiming is identifiable, and by what strategy.
  - It verifies DoWhy's IV estimator recovers the same estimand as our IVW/WLS
    regression, so the two analyses can be reported as one coherent story.
  - It shows how badly a naive (confounded) regression misses, which is the
    argument for doing MR at all.
  - It runs refutation tests, which is the concrete "causal AI" component the
    write-up has so far only gestured at.

To turn this into a real-data analysis you need individual-level genotype +
phenotype data (UK Biobank, FinnGen individual level, etc). Summary statistics
alone cannot do it. See PROJECT_STATUS.md.

Run:  ./venv/bin/python causality.py
"""

import numpy as np
import pandas as pd
from dowhy import CausalModel

# ---------------------------------------------------------------------------
# Calibration from the real MR run (archived_outputs/output10.md)
# exposure prot-a-1538, outcome ebi-a-GCST90014023 (T1D)
# ---------------------------------------------------------------------------

# Real SNP -> exposure effects. These fix how strong our instruments are.
INSTRUMENTS = pd.DataFrame({
    "rsid":     ["rs1554606", "rs1524107", "rs2069852"],
    "beta_exp": [-0.015253,   -0.025149,   -0.026594],
    "se_exp":   [0.001862,     0.004419,    0.005321],
    # PLACEHOLDER effect-allele frequencies. OpenGWAS returns these in the
    # 'eaf' field -- swap them in once the API token is renewed. They only
    # affect instrument variance, not the direction of the result.
    "eaf":      [0.40,         0.30,        0.25],
})

# The IVW/WLS causal estimate we are trying to reproduce with DoWhy.
IVW_BETA = 3.4829
IVW_SE = 0.702

N_PEOPLE = 50_000
SEED = 20260801


def simulate_cohort(true_effect=IVW_BETA, n=N_PEOPLE, seed=SEED):
    """Individual-level cohort consistent with the real SNP-exposure betas.

    Structure (this is the DAG we assert in causal_graph()):

        G1,G2,G3 --> IL6 --> T1D
        BMI, AGE, SMOKING --> IL6, T1D      (measured confounders)
        U --> IL6, T1D                      (unmeasured confounding)

    The instruments touch T1D ONLY through IL6. That is the exclusion
    restriction, assumed here by construction and untestable in the real data.
    """
    rng = np.random.default_rng(seed)

    # --- genotypes: additive allele count, 0/1/2 under Hardy-Weinberg ---
    genotypes = {}
    for _, snp in INSTRUMENTS.iterrows():
        genotypes[snp["rsid"]] = rng.binomial(2, snp["eaf"], size=n)

    # --- confounders ---
    age = rng.normal(45, 14, n)
    bmi = rng.normal(26.5, 4.5, n) + 0.03 * (age - 45)
    smoking = rng.binomial(1, 0.22, n)

    # U is real but never recorded: shared inflammatory / autoimmune liability,
    # chronic infection burden, population structure. This is the whole reason
    # a naive IL6 ~ T1D regression cannot be read causally.
    u = rng.normal(0, 1, n)

    # --- exposure: IL-6, in the normalised units the pQTL GWAS reports ---
    il6 = np.zeros(n)
    for _, snp in INSTRUMENTS.iterrows():
        il6 += snp["beta_exp"] * genotypes[snp["rsid"]]
    il6 += 0.030 * (bmi - 26.5)      # adiposity raises circulating IL-6
    il6 += 0.004 * (age - 45)
    il6 += 0.060 * smoking
    il6 += 0.250 * u
    il6 += rng.normal(0, 0.30, n)

    # --- outcome ---
    # Simulated on a LINEAR scale so the IV estimand is directly comparable to
    # IVW_BETA. The real outcome GWAS is binary (log-odds); see the scale note
    # printed at the end of the run.
    t1d = (
        true_effect * il6
        + 0.010 * (bmi - 26.5)
        - 0.002 * (age - 45)
        + 0.030 * smoking
        + 0.400 * u                  # U pushes the naive estimate off
        + rng.normal(0, 1.0, n)
    )

    df = pd.DataFrame(genotypes)
    df["BMI"] = bmi
    df["AGE"] = age
    df["SMOKING"] = smoking
    df["U"] = u
    df["IL6"] = il6
    df["T1D"] = t1d
    return df


def causal_graph():
    """DAG in DOT. U is declared but will be withheld from the data."""
    snps = list(INSTRUMENTS["rsid"])
    edges = []
    for s in snps:
        edges.append(f"{s} -> IL6;")
    edges.append("IL6 -> T1D;")
    for c in ["BMI", "AGE", "SMOKING", "U"]:
        edges.append(f"{c} -> IL6;")
        edges.append(f"{c} -> T1D;")
    edges.append("AGE -> BMI;")
    return "digraph { " + " ".join(edges) + " }"


def naive_estimate(df):
    """What you would get regressing T1D on IL6 adjusting only for what you
    can measure. Biased by U -- this is the comparison that motivates MR."""
    import statsmodels.api as sm
    X = sm.add_constant(df[["IL6", "BMI", "AGE", "SMOKING"]])
    return sm.OLS(df["T1D"], X).fit().params["IL6"]


def run():
    print("=" * 72)
    print("DoWhy causal model: IL-6 -> Type 1 Diabetes")
    print("SYNTHETIC cohort calibrated to prot-a-1538 / ebi-a-GCST90014023")
    print("=" * 72)

    df = simulate_cohort()
    # U is unmeasured by construction: drop it before anything sees the data.
    observed = df.drop(columns=["U"])
    print(f"\nCohort: {len(observed):,} individuals, "
          f"{len(INSTRUMENTS)} instruments, 3 measured confounders + 1 unmeasured (U)")
    print(f"True causal effect baked into the simulation: {IVW_BETA:.4f}")

    # ---- 1. model ----
    model = CausalModel(
        data=observed,
        treatment="IL6",
        outcome="T1D",
        graph=causal_graph(),
        instruments=list(INSTRUMENTS["rsid"]),
    )

    # ---- 2. identify ----
    print("\n" + "-" * 72)
    print("STEP 1 — IDENTIFICATION (what is estimable given the DAG?)")
    print("-" * 72)
    estimand = model.identify_effect(proceed_when_unidentifiable=False)
    print(estimand)

    # ---- 3. estimate ----
    print("\n" + "-" * 72)
    print("STEP 2 — ESTIMATION")
    print("-" * 72)

    iv = model.estimate_effect(estimand, method_name="iv.instrumental_variable")
    iv_effect = float(iv.value)

    naive = naive_estimate(observed)

    print(f"\n  Naive regression (adjusts BMI/AGE/SMOKING only) : {naive:.4f}")
    print(f"  DoWhy IV estimate                              : {iv_effect:.4f}")
    print(f"  IVW/WLS from mr_t1d.py (summary-level MR)       : {IVW_BETA:.4f} "
          f"(95% CI {IVW_BETA - 1.96 * IVW_SE:.3f} to {IVW_BETA + 1.96 * IVW_SE:.3f})")
    print(f"\n  Naive bias vs truth : {naive - IVW_BETA:+.4f}")
    print(f"  DoWhy IV bias vs truth : {iv_effect - IVW_BETA:+.4f}")
    if abs(iv_effect - IVW_BETA) < 1.96 * IVW_SE:
        print("  -> IV estimate lands inside the IVW confidence interval. The DoWhy")
        print("     IV estimator and our IVW regression target the same estimand.")
    else:
        print("  -> IV estimate is OUTSIDE the IVW interval; check calibration.")

    # ---- 4. refute ----
    print("\n" + "-" * 72)
    print("STEP 3 — REFUTATION (does the estimate behave the way a causal")
    print("         estimate should when we attack it?)")
    print("-" * 72)

    # --- DoWhy's built-in refuters ---
    for name, why, kwargs in [
        ("data_subset_refuter",
         "Re-estimate on random 80% subsets. Should be stable.",
         {"subset_fraction": 0.8, "num_simulations": 20}),
        ("random_common_cause",
         "Add an irrelevant common cause. A valid estimate should NOT move.",
         {}),
        ("placebo_treatment_refuter",
         "Replace IL6 with a permutation of itself. See caveat below.",
         {"placebo_type": "permute", "num_simulations": 20}),
    ]:
        print(f"\n[{name}]\n  {why}")
        try:
            res = model.refute_estimate(estimand, iv, method_name=name, **kwargs)
            print("  " + str(res).replace("\n", "\n  "))
        except Exception as e:
            print(f"  DID NOT RUN -- {type(e).__name__}: {e}")
            print("  (DoWhy's generic refuters are written against backdoor/adjustment")
            print("   estimators; several do not compose with iv.instrumental_variable.")
            print("   Do not report this as a passed test.)")

    print("""
  CAVEAT ON placebo_treatment_refuter: this test is NOT valid for an IV
  estimator and its output above should be ignored. Permuting IL6 destroys the
  first stage (SNP -> IL6), so the Wald ratio's denominator goes to ~0 and the
  ratio explodes instead of collapsing to 0. A non-zero "new effect" here is an
  artefact of dividing by noise, not evidence against our estimate. The IV
  analogues that DO work are below.""")

    # --- IV-appropriate falsification tests (these are the ones that matter) ---
    print("\n" + "-" * 72)
    print("STEP 3b — IV-SPECIFIC FALSIFICATION")
    print("-" * 72)
    _negative_control_outcome(observed, df)
    _overidentification_test(observed)
    _leave_one_out(observed)


def _iv_wald(df, snps, outcome="T1D", exposure="IL6"):
    """Two-stage least squares by hand, so the tests below don't depend on
    DoWhy's refuter internals."""
    import statsmodels.api as sm
    Z = sm.add_constant(df[snps])
    first = sm.OLS(df[exposure], Z).fit()
    second = sm.OLS(df[outcome], sm.add_constant(first.fittedvalues)).fit()
    return second.params.iloc[1], second.resid, first


def _negative_control_outcome(observed, full):
    """A negative control outcome: something driven by the SAME confounders and
    the SAME unmeasured U as T1D, but NOT by IL-6. A valid instrument set must
    return ~0 here. If it doesn't, the SNPs are reaching the outcome by some
    path other than IL-6 -- i.e. the exclusion restriction is broken."""
    rng = np.random.default_rng(SEED + 1)
    nc = (0.01 * (full["BMI"] - 26.5) + 0.03 * full["SMOKING"]
          + 0.40 * full["U"] + rng.normal(0, 1.0, len(full)))
    d = observed.copy()
    d["NEG_CONTROL"] = nc
    snps = list(INSTRUMENTS["rsid"])
    est, _, _ = _iv_wald(d, snps, outcome="NEG_CONTROL")
    print(f"\n[negative control outcome]")
    print(f"  IV estimate on an outcome IL-6 does not cause: {est:+.4f}")
    print(f"  {'PASS' if abs(est) < 0.5 else 'FAIL'} (expect ~0; "
          f"compare to the real estimate of {IVW_BETA:.2f})")


def _overidentification_test(observed):
    """Sargan test. With 3 instruments and 1 exposure we have 2 overidentifying
    restrictions, so instrument validity becomes partially testable. This is the
    same idea as Cochran's Q / the MR-Egger intercept in the MR literature:
    if the SNPs disagree with each other about the causal effect, at least one
    of them is invalid."""
    from scipy import stats as st
    import statsmodels.api as sm
    snps = list(INSTRUMENTS["rsid"])
    _, resid, _ = _iv_wald(observed, snps)
    aux = sm.OLS(resid, sm.add_constant(observed[snps])).fit()
    n = len(observed)
    sargan = n * aux.rsquared
    dfree = len(snps) - 1
    p = 1 - st.chi2.cdf(sargan, dfree)
    print(f"\n[Sargan overidentification test]")
    print(f"  statistic = {sargan:.3f} on {dfree} df, p = {p:.3f}")
    print(f"  {'PASS' if p > 0.05 else 'FAIL'} — a small p would say the instruments")
    print(f"  disagree, implicating pleiotropy or a broken exclusion restriction.")


def _leave_one_out(observed):
    """Is the estimate driven by one influential SNP? Standard MR practice."""
    snps = list(INSTRUMENTS["rsid"])
    print(f"\n[leave-one-out]")
    full_est, _, _ = _iv_wald(observed, snps)
    print(f"  all {len(snps)} instruments : {full_est:+.4f}")
    for drop in snps:
        keep = [s for s in snps if s != drop]
        est, _, _ = _iv_wald(observed, keep)
        print(f"  dropping {drop:<12}: {est:+.4f}  (shift {est - full_est:+.4f})")

    # ---- 5. caveats, printed every run on purpose ----
    print("\n" + "=" * 72)
    print("HOW TO READ THIS")
    print("=" * 72)
    print("""
1. SYNTHETIC DATA. The cohort is simulated. This run validates the method and
   the DAG; it is not independent evidence that IL-6 causes T1D. The only real
   inputs are the three SNP-exposure betas and the IVW estimate.

2. OUTCOME SCALE. T1D is simulated on a linear scale so the IV estimand matches
   IVW_BETA exactly. The real T1D GWAS is a case/control study reporting
   log-odds. Linear IV on a binary outcome estimates a risk difference, not a
   log-OR, so on real individual-level data these two numbers would NOT be
   expected to match without a scale correction.

3. THE EXCLUSION RESTRICTION IS ASSUMED, NOT TESTED. Here it is true by
   construction. In real data, a SNP affecting T1D through anything other than
   IL-6 breaks the whole analysis, and no refutation test in DoWhy can detect
   it. MR-Egger and weighted median are the tools for that, and they are still
   not implemented -- see PROJECT_STATUS.md.

4. EAFs ARE PLACEHOLDERS. Replace INSTRUMENTS['eaf'] with the real 'eaf' field
   from OpenGWAS once the API token is renewed.
""")


if __name__ == "__main__":
    run()
