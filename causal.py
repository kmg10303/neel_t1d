import numpy as np  
import pandas as pd 
from dowhy import CausalModel

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
    import statsmodels.api as sm
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


    # Real estimate


    





if __name__ == "__main__":
    run()



