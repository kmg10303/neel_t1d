

#HLA-DRB1
#HLA-DQA1
#HLA-DQB1

#rs2187668

import os
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

CACHE_DIR = "outputs"


def _cache_path(exposure_id, outcome_id):
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, f"harmonized_{exposure_id}_{outcome_id}.csv")


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


def get_instruments():
    if not USE_TOPHITS:
        return IL6_INSTRUMENTS
    print("Fetching tophits + LD clumping from OpenGWAS (requires valid JWT)...")
    try:
        hits = gwas.tophits(id=[GWASid], pval=5e-5, clump=1, kb=10000, r2=0.001, pop="EUR")
    except Exception as e:
        # When the API is fully unreachable (not just an auth-error dict), the
        # underlying ieugwaspy call raises instead of returning a dict/empty list.
        print(f"WARNING: tophits raised {type(e).__name__}: {e}. Falling back to hardcoded IL6_INSTRUMENTS.")
        return IL6_INSTRUMENTS
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
    return snp_list

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


def weighted_median(merged, n_boot=1000, seed=0):
    """Weighted median of per-SNP Wald ratios, with a bootstrap SE."""
    def _wm(d):
        # reset_index is required: bootstrap resamples (replace=True) create
        # duplicate index labels, and the .loc[i0]/.loc[i1] lookups below
        # return a Series instead of a scalar when the index isn't unique.
        d = d.sort_values("wald").reset_index(drop=True)
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


def mr_analysis():
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

        cache_file = _cache_path(GWASid, OUTCOME_ID)
        merged = None

        try:
            instruments = get_instruments()
            print(f"Fetching IL-6 exposure data for {len(instruments)} instruments...")
            exposure_raw = _batch_associations(instruments, [GWASid])

            if not exposure_raw or isinstance(exposure_raw, dict):
                raise RuntimeError(f"No exposure data found. Check GWASid or API token. "
                                    f"({type(exposure_raw)}, {exposure_raw})")

            exposure = pd.DataFrame(exposure_raw).rename(columns=col_map)
            snps = exposure["rsid"].dropna().unique().tolist()
            print(f"Found {len(snps)} genetic instruments.")

            # 2. Get Outcome Data
            outcome_raw = _batch_associations(snps, [OUTCOME_ID])
            if not outcome_raw:
                raise RuntimeError("No outcome data found.")

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

            # Palindromic SNPs (A/T, C/G) are strand-ambiguous: "same" and "swap"
            # can both evaluate true for them without a trustworthy EAF to resolve
            # strand, so exclude them outright rather than risk a silent sign flip.
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

            # Harmonization Logic
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

            if merged.empty:
                raise RuntimeError("No harmonized SNPs left.")

            merged.to_csv(cache_file, index=False)
            print(f"Cached harmonized data -> {cache_file}")

        except Exception as e:
            if os.path.exists(cache_file):
                print(f"Live API call failed ({e}); falling back to cached harmonized data ({cache_file}).")
                merged = pd.read_csv(cache_file)
            else:
                raise

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

        mr_egger(merged)
        weighted_median(merged)
        cochrans_q(merged, wls.params.iloc[0])

        # print(wls.params)
        x_line = np.linspace(merged["beta_exp"].min(), merged["beta_exp"].max())
        plot_mr(x_line, wls.params.iloc[0] * x_line, merged)

    except Exception:
        import traceback
        traceback.print_exc()
        raise

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
