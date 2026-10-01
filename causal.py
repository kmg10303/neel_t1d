"""
DoWhy causal model for the IL-6 -> T1D question.

DoWhy needs individual-level data (one row per person). mr_t1d.py's IVW/WLS
pipeline operates on two-sample SUMMARY statistics (one row per SNP). Those are
not the same thing, so this script builds a synthetic cohort CALIBRATED to the
real SNP-exposure betas from prot-a-1538, with a known true causal effect set
to the IVW estimate from mr_t1d.py. This is a method-validation study, not new
evidence about T1D -- see PROJECT_STATUS.md for what that does and doesn't mean.

MR-Egger / weighted median / Cochran's Q are NOT run here: those operate on
summary-stats (beta_exp, beta_out, se_out per SNP), which this individual-level
simulation doesn't have. They live in mr_t1d.py, where the real summary data is.

Run:  ./venv/bin/python causal.py
"""

import numpy as np
import pandas as pd
from dowhy import CausalModel
import statsmodels.api as sm

INSTRUMENTS = pd.DataFrame({
    "rsid":     ["rs1554606", "rs1524107", "rs2069852"],
    "beta_exp": [-0.015253,   -0.025149,   -0.026594],
    "se_exp":   [0.001862,     0.004419,    0.005321],

    # These are example values
    "eaf":      [0.40,         0.30,        0.25],
})

IVW_BETA = 3.4829
IVW_SE = 0.702

N_PEOPLE = 50_000
SEED = 20260801

def simulate(n=N_PEOPLE, seed=SEED):
    rng = np.random.default_rng(seed)

    # Confounding variables
    genotypes = {}
    for _, snp in INSTRUMENTS.iterrows():
        genotypes[snp["rsid"]] = rng.binomial(2, snp["eaf"], size=n)

    # confounders -> These are all values that might need tweaking in order to create a better representation of the world.
    age = rng.normal(45, 14, n)
    bmi = rng.normal(25, 5, n) + 0.03 * (age - 45)
    smoking = rng.binomial(1, 0.22, n)

    u = rng.normal(0, 1, n)

    # Exposure
    il6 = np.zeros(n)

    for _, snp in INSTRUMENTS.iterrows():
        il6 += snp["beta_exp"] * genotypes[snp["rsid"]]

    # Add in our confounders, weighting each.

    il6 += (bmi - 25) * 0.030
    il6 += (age - 45) * 0.004
    il6 += (smoking) * 0.060
    il6 += u * 0.250
    il6 += rng.normal(0, 0.30, n)

    # Outcome

    t1d = (
        IVW_BETA * il6
        + 0.010 * (bmi - 25)
        - .002 * (age - 45)
        + 0.030 * smoking
        + 0.4 * u
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


def _iv_wald(df, snps, outcome="T1D", exposure="IL6"):
    """Two-stage least squares by hand, independent of DoWhy's refuter
    internals -- so these tests don't depend on what DoWhy's IV estimator is
    doing under the hood."""
    Z = sm.add_constant(df[snps])
    first = sm.OLS(df[exposure], Z).fit()
    second = sm.OLS(df[outcome], sm.add_constant(first.fittedvalues)).fit()
    return second.params.iloc[1], second.resid, first


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


# Might have to fiddle around with the DAG here. Adding/removing confounding variables
def causal_graph():
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
    X = sm.add_constant(df[["IL6", "BMI", "AGE", "SMOKING"]])
    return sm.OLS(df["T1D"], X).fit().params["IL6"]

def run():

    df = simulate()

    observed = df.drop(columns=["U"])

    # 1. Model
    model = CausalModel(
        data=observed,
        treatment="IL6",
        outcome="T1D",
        graph=causal_graph(),
        instruments=list(INSTRUMENTS["rsid"]),
    )

    # 2. Identify effects
    print("STEP 2 - Identification (what is estimable given the DAG?)")
    estimand = model.identify_effect(proceed_when_unidentifiable=False)
    print(estimand)

    iv = model.estimate_effect(estimand, method_name="iv.instrumental_variable")
    iv_effect = float(iv.value)

    # 3. Naive Estimate

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

    # 4. Falsification tests
    print("\n----- Falsification tests -----")
    _negative_control_outcome(observed, df)
    _overidentification_test(observed)
    _leave_one_out(observed)

    # 5. Reconcile DoWhy's IV estimator against hand-rolled 2SLS
    hand_2sls, _, _ = _iv_wald(observed, list(INSTRUMENTS["rsid"]))
    print(f"\n  Hand-rolled 2SLS       : {hand_2sls:.4f}")
    print(f"  DoWhy IV estimate      : {iv_effect:.4f}")
    print(f"  Difference             : {hand_2sls - iv_effect:+.4f} "
          f"({100 * (hand_2sls - iv_effect) / iv_effect:+.1f}%)")


if __name__ == "__main__":
    run()
