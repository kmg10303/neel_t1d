

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

GWASid = "ebi-a-GCST004448"
OUTCOME_ID = "ebi-a-GCST90014023"

HLA_SNPs = ["rs2187668", "rs9273368", "rs9273363", "rs9272346", "rs2647044"]


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

def mr_analysis(variants = [], snps=HLA_SNPs):
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


        # variants = ["rs9272346"]
        exposure_raw = gwas.tophits(GWASid, pval=5e-8, clump=1) # Clump=1 is usually better for MR
        print(type(exposure_raw))
        print(exposure_raw)
        if not exposure_raw:
            print("No significant SNPs found via tophits. Trying specific variant...")
            exposure_raw = gwas.associations(variant=variants, id=[GWASid])
        
        if not exposure_raw:
            print("No exposure data found.")
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
        
        if merged.empty:
            print("No harmonized SNPs left.")
            return
        
        num_snps = len(merged)

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
        x_line = np.linspace(merged["beta_exp"].min(), merged["beta_out"].max())
        plot_mr(x_line, wls.params[0] * x_line, merged)
        
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
    gwas_id_to_check = ["rs2187668", "rs9273368", "rs9273363", "rs9272346", "rs2647044"]
    mr_analysis(gwas_id_to_check, HLA_SNPs)
    print("Done")
