

#HLA-DRB1
#HLA-DQA1
#HLA-DQB1

#rs2187668

import pandas as pd
import ieugwaspy as gwas
import statsmodels.api as sm
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

# Ahola-Olli 2017 IL-6 pQTL GWAS (ieu-b- batch, 41 cytokines, N=8293)
# Verify ID with: gwas.gwasinfo(["ieu-b-30"])
# https://opengwas.io/datasets/
#ebi-a-GCST004446
#ieu-b-30
# prot-a-1538
# ukb-e-30000_CSA - White blood count
GWASid = "prot-a-1538"

# Outcome ID
OUTCOME_ID = "ebi-a-GCST90014023"

# Literature-curated IL-6 pQTL instruments (Ahola-Olli 2017 + IL6R literature)
IL6_PRIMARY = "rs4537545"
IL6_INSTRUMENTS = [
    "rs1524107",
    "rs1554606",
    "rs2069852",
]

# Set True to fetch instruments from OpenGWAS tophits with LD clumping (requires valid JWT).
# Set False to use the hardcoded IL6_INSTRUMENTS list above.
USE_TOPHITS = True

def get_instruments():
    if not USE_TOPHITS:
        return IL6_INSTRUMENTS
    print("Fetching tophits + LD clumping from OpenGWAS (requires valid JWT)...")
    hits = gwas.tophits(id=[GWASid], pval=5e-5, clump=1, kb=10000, r2=0.001, pop="EUR")
    # API returns a dict with "message" key on auth errors instead of a list
    if not hits or isinstance(hits, dict):
        msg = hits.get("message", "empty response") if isinstance(hits, dict) else "empty response"
        print(f"WARNING: tophits failed ({msg}). Falling back to hardcoded IL6_INSTRUMENTS.")
        return IL6_INSTRUMENTS
    snp_list = [h["rsid"] for h in hits if "rsid" in h]
    if not snp_list:
        print("WARNING: tophits returned hits but no rsids. Falling back to hardcoded IL6_INSTRUMENTS.")
        return IL6_INSTRUMENTS
    print(f"Tophits returned {len(snp_list)} clumped instruments: {snp_list}")
    return snp_list[:63]

# Retained for backward-compatibility with gwas_id_check()
HLA_SNPs = ["rs12722495-A", "rs61839660-C", "rs12722496-G"]


# x = []; y = []
# plt.plot(x, y, label="")
# plt.xlabel("")
# plt.ylabel("")

def plot_mr(x, y, merged):
    plt.errorbar(merged["beta_exp"], merged["beta_out"],
             xerr=merged["se_exp"], yerr=merged["se_out"],
             fmt='o', color='black', label='SNPs')
    plt.plot(x, y, label="IVW slope")
    plt.xlabel("input")
    plt.ylabel("outcome")
    plt.legend()
    plt.savefig("scatterplot.png")

def mr_analysis(snps=HLA_SNPs):
    try:
        print("Loading data")
        # Standardizing column mapping for ieugwaspy
        col_map = {
            "ea": "effect_allele",
            "nea": "other_allele",
            "beta": "beta",
            "se": "se",
            "rsid": "rsid"
        }


        instruments = get_instruments()
        print(f"Fetching IL-6 exposure data for {len(instruments)} instruments...")
        exposure_raw = gwas.associations(variant=instruments, id=[GWASid])

        if not exposure_raw or isinstance(exposure_raw, dict):
            print("No exposure data found. Check GWASid or API token.")
            print(type(exposure_raw), exposure_raw)
            return
        

        exposure = pd.DataFrame(exposure_raw).rename(columns=col_map)
        snps = exposure["rsid"].dropna().unique().tolist()
        print(f"Found {len(snps)} genetic instruments.")

        # 2. Get Outcome Data
        outcome_raw = gwas.associations(variant=snps, id=[OUTCOME_ID])
        if not outcome_raw:
            print("No outcome data found.")
            return
            
        outcome = pd.DataFrame(outcome_raw).rename(columns=col_map)

        # 3. Merge and Harmonize
        merged = pd.merge(
            exposure[["rsid", "beta", "se", "effect_allele", "other_allele"]].rename(
                columns={"beta": "beta_exp", "se": "se_exp"}
            ),
            outcome[["rsid", "beta", "se", "effect_allele", "other_allele"]].rename(
                columns={"beta": "beta_out", "se": "se_out"}
            ),
            on="rsid",
            suffixes=("_exp", "_out")
        )

        # Harmonization Logic
        same = (merged["effect_allele_exp"] == merged["effect_allele_out"]) & \
               (merged["other_allele_exp"] == merged["other_allele_out"])
        
        swap = (merged["effect_allele_exp"] == merged["other_allele_out"]) & \
               (merged["other_allele_exp"] == merged["effect_allele_out"])
        
        merged = merged[same | swap].copy()
        merged.loc[swap, "beta_out"] *= -1

        # summary_df = merged.rename(columns={"beta_exp": "IL6", "beta_out": "T1D"})
        
        if merged.empty:
            print("No harmonized SNPs left.")
            return
        
        num_snps = len(merged)

        # F-statistic per instrument
        merged["F_stat"] = (merged["beta_exp"] / merged["se_exp"]) ** 2
        merged["weak"] = merged["F_stat"] < 10
        print("\n----- Instrument Strength (F-statistics) -----")
        print(merged[["rsid", "beta_exp", "se_exp", "F_stat", "weak"]].to_string(index=False))
        weak_count = int(merged["weak"].sum())
        if weak_count:
            print(f"WARNING: {weak_count} instrument(s) have F < 10 (weak instrument bias risk).")

        # If only one SNP, perform Wald ratio
        if num_snps < 2:
            row = merged.iloc[0]
            b_exp = row["beta_exp"]
            b_out = row["beta_out"]
            se_out = row["se_out"]
            
            wald_beta = b_out / b_exp
            wald_se = se_out / abs(b_exp)  # Simplified SE calculation
            z_score = wald_beta / wald_se
            p_val = 2 * (1 - stats.norm.cdf(abs(z_score)))

            print(f"SNP: {row['rsid']}")
            print(f"Causal Estimate (Beta): {wald_beta:.4f}")
            print(f"Standard Error: {wald_se:.4f}")
            print(f"P-value: {p_val:.4f}")
            return

        # 4. Regression (IVW)
        X = merged["beta_exp"]
        y = merged["beta_out"]
        # Inverse Variance Weighting
        weights = 1 / (merged["se_out"]**2)
        wls = sm.WLS(y, X, weights=weights).fit()

        print("\n----- MR Results (IVW via WLS) -----")
        print(wls.summary())

        # print(wls.params)
        x_line = np.linspace(merged["beta_exp"].min(), merged["beta_exp"].max())
        plot_mr(x_line, wls.params.iloc[0] * x_line, merged)
        
    except Exception as e:
        print(f"An error occurred: {e}")

def gwas_id_check(ids, snps):
    for id in ids:
        info = gwas.gwasinfo([id])
        
        if info:
            print("Found info in id", id)
        else:
            print("None found")

        row = info[0]

        # Found info in id ebi-a-GCST90002009
        # {'id': 'ebi-a-GCST90002009', 'trait': 'HLA DR on CD14- CD16-', 'build': 'HG19/GRCh37', 'group_name': 'public', 'category': 'NA', 'subcategory': 'NA', 'population': 'European', 'sex': 'NA', 'author': 'Orr<U+00F9> V', 'nsnp': 15034296, 'sample_size': 3629, 'year': 2020, 'ontology': 'NA', 'unit': 'NA', 'consortium': 'NA', 'pmid': 32929287, 'mr': 1, 'priority': 0, 'note': 'NA'}
        # Found info in id ebi-a-GCST90002010
        # {'id': 'ebi-a-GCST90002010', 'trait': 'HLA DR on monocyte', 'build': 'HG19/GRCh37', 'group_name': 'public', 'category': 'NA', 'subcategory': 'NA', 'population': 'European', 'sex': 'NA', 'author': 'Orr<U+00F9> V', 'nsnp': 15034296, 'sample_size': 3629, 'year': 2020, 'ontology': 'NA', 'unit': 'NA', 'consortium': 'NA', 'pmid': 32929287, 'mr': 1, 'priority': 0, 'note': 'NA'}
        # Found info in id ebi-a-GCST004448
        # {'id': 'ebi-a-GCST004448', 'trait': 'Interleukin-1-beta levels', 'build': 'HG19/GRCh37', 'group_name': 'public', 'category': 'NA', 'subcategory': 'NA', 'population': 'European', 'sex': 'NA', 'author': 'Ahola-Olli AV', 'nsnp': 9983642, 'sample_size': 3309, 'year': 2016, 'ontology': 'NA', 'unit': 'NA', 'consortium': 'NA', 'pmid': 27989323, 'mr': 1, 'priority': 0, 'note': 'NA'}
        # Done
        # Trait, Sample Size, # SNPs, Population, Build
        print(row.get('author'))
        print(row.get('trait'))
        print(row.get('sample_size'))
        print(row.get('nsnp'))
        print(row.get('population'))
        print(row.get('build'))


    for id in ids:
        # Check tophits for each. 
        hits = gwas.associations(variant=snps, id=[id])
        print(f"Top hits for {id}: {len(hits)}")

        # Top hits for ebi-a-GCST90002009: 4
        # Top hits for ebi-a-GCST90002010: 4
        # Top hits for ebi-a-GCST004448: 5 --> Ahola-Olli AV --> https://opengwas.io/datasets/
        # Done

if __name__ == "__main__":
    mr_analysis()
    print("Done")



