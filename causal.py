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

def iv_wald(df, snps, outcome="T1D", exposure="IL6"):
    Z = sm.add_constant(df[snps])
    
    first = sm.OLS(df[exposure], Z).fit()
    second = sm.OLD(df[outcome], sm.add_constant(first.fittedvalues)).fit()
    return second.params[1], second.pvalues[1], first

# We already have an outcome.
# What if we took our exposure and other data, ignored the outcome, and tried to create our own outcome.
# So now we can compare our two outcomes to see how alike they are. 
def negative_control_outcome(____, _____): 
    # we need some randomness
    # rng = np.random... 
    # nc = weight * (full["BMI"] - 25) + weight 
    # Note, we do want to include U in this. 
    # Now we can compare the negative control with the ovserved.
    # d = observed.copy() [NOTE: You will have to pass observed into the function as parameter]. 
    # Sequence of prints
    # print("negative control outcome test")
    # PASS = abs(est - obs) < 0.5 else FAIL (we want something really close to 0). 

    pass

def leave_one_out(____):
    # First grab the snps
    # print ("leave one out test")
    # full_est_, _, _ = iv_wald(observed, snps) [NOTE: observed is passed into the function as a parameter. SNPS is grabbed at the top line of the function.]
    # for loop
    # keep snps except one
    # 


    pass


def mr_egger(merged):
    d = merged.copy()

    X = sm.add_constant(d["beta_exp"])

    egger = sm.WLS(d["beta_out"], X, weights=weights).fit()

    intercept, slope = egger.params["const"], egger.params["beta_exp"]
    intercept_p, slope_p = egger.pvalues["const"], egger.pvalues["beta_exp"]

    print("---- MR-Egger ----")
    print(f"Slope (causal estimate): {slope:.4f}")
    print(f" Intercept (pleiotropy): {intercept:.4f}")

    if intercept_p < .1:
        print("The intercept is significant, suggesting pleiotropy. Prefer to use the Egger slope")
    else:
        print("The intercept is not significant, suggesting no pleiotropy")
    return slope, slope_p, intercept, intercept_p

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


    # MR Egger regression

    mr_egger(observed)

    # Wald Regression
    iv_wald(observed, list(INSTRUMENTS["rsid"]))




if __name__ == "__main__":
    run()



