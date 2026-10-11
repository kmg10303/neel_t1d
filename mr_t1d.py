

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

# https://opengwas.io/datasets/
# Exposure datasets (all confirmed with gwas.gwasinfo on 2026-10-09):
#   prot-a-1538         Interleukin-6, Sun BB 2018 (PMID 29875488), N=3,301
#                       NOTE: this is Sun 2018, NOT Ahola-Olli. It has zero hits at 5e-8.
#   ebi-a-GCST90012005  Interleukin-6 levels, Folkersen L 2020 (PMID 33067605), N=21,758
#   ebi-a-GCST004446    Interleukin-6 levels, Ahola-Olli AV 2016, N=8,189 (zero hits at 5e-8)
# ukb-e-30000_CSA - White blood count
SUN_IL6 = "prot-a-1538"
FOLKERSEN_IL6 = "ebi-a-GCST90012005"

# Outcome ID: Type 1 diabetes, Chiou J 2021 (PMID 34012112), 18,942 cases / 501,638 controls
OUTCOME_ID = "ebi-a-GCST90014023"

# Each entry is one MR analysis. "mode" is "genome" (all clumped hits) or "cis"
# (only hits inside CIS_REGIONS).
RUNS = [
    {"label": "sun_5e-5",       "exposure": SUN_IL6,       "pval": 5e-5, "mode": "genome"},  # original weak-instrument run, kept as a contrast
    {"label": "sun_5e-8",       "exposure": SUN_IL6,       "pval": 5e-8, "mode": "genome"},
    {"label": "folkersen_5e-8", "exposure": FOLKERSEN_IL6, "pval": 5e-8, "mode": "genome"},
    {"label": "folkersen_cis",  "exposure": FOLKERSEN_IL6, "pval": 5e-8, "mode": "cis"},
]

# Gene coordinates (GRCh37) plus a flanking window: (chromosome, start, end)
CIS_REGIONS = {
    "IL6":  ("7", 22766766 - 500_000, 22771621 + 500_000),
    "IL6R": ("1", 154377669 - 500_000, 154441926 + 500_000),
}

# The MHC is the strongest T1D locus in the genome and has long-range LD, so any
# instrument there is almost certainly pleiotropic. Extended MHC, GRCh37.
EXCLUDE_MHC = True
MHC_REGION = ("6", 25_000_000, 34_000_000)

# Drop instruments that explain more variance in T1D than in IL-6 (wrong causal direction).
STEIGER_FILTER = True
# Population lifetime prevalence of T1D, used to put the outcome r2 on the liability scale.
T1D_PREVALENCE = 0.005

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


def _cache_path(label, exposure_id, outcome_id):
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, f"harmonized_{label}_{exposure_id}_{outcome_id}.csv")


def _batch_associations(variants, ids, batch_size=63):
    """OpenGWAS caps N(id) * N(variant) <= 64 per call. Batch instead of
    truncating the instrument list.

    proxies=0: by default the API silently substitutes an LD proxy (a different
    SNP, with different betas) but still labels it with the rsid we asked for."""
    results = []
    for i in range(0, len(variants), batch_size):
        chunk = variants[i:i + batch_size]
        res = gwas.associations(variant=chunk, id=ids, proxies=0)
        if res and not isinstance(res, dict):
            results.extend(res)
    return results


def _in_region(df, region):
    chrom, start, end = region
    pos = pd.to_numeric(df["position"], errors="coerce")
    return (df["chr"].astype(str) == chrom) & (pos >= start) & (pos <= end)


def get_instruments(exposure_id, pval):
    if not USE_TOPHITS:
        return IL6_INSTRUMENTS
    print(f"Fetching tophits (p < {pval:g}) + LD clumping from OpenGWAS (requires valid JWT)...")
    try:
        hits = gwas.tophits(id=[exposure_id], pval=pval, clump=1, kb=10000, r2=0.001, pop="EUR")
    except Exception as e:
        # When the API is fully unreachable (not just an auth-error dict), the
        # underlying ieugwaspy call raises instead of returning a dict/empty list.
        print(f"WARNING: tophits raised {type(e).__name__}: {e}. Falling back to hardcoded IL6_INSTRUMENTS.")
        return IL6_INSTRUMENTS
    # API returns a dict with "message" key on auth errors instead of a list
    if isinstance(hits, dict):
        print(f"WARNING: tophits failed ({hits.get('message', 'no message')}). Falling back to hardcoded IL6_INSTRUMENTS.")
        return IL6_INSTRUMENTS
    # An empty list is a real answer: nothing in this GWAS passes the threshold.
    snp_list = [h["rsid"] for h in hits if "rsid" in h]
    print(f"Tophits returned {len(snp_list)} clumped instruments: {snp_list}")
    return snp_list


def describe_datasets(exposure_id, outcome_id):
    """Print what each GWAS ID actually is, and return the two info rows."""
    info = {row["id"]: row for row in gwas.gwasinfo([exposure_id, outcome_id])}
    print("\n----- Datasets -----")
    for role, gid in (("Exposure", exposure_id), ("Outcome", outcome_id)):
        row = info[gid]
        size = f"N={row.get('sample_size')}"
        if row.get("ncase"):
            size += f" ({row['ncase']} cases / {row['ncontrol']} controls)"
        print(f"  {role}: {gid} -- {row.get('trait')}, {row.get('author')} {row.get('year')}, "
              f"PMID {row.get('pmid')}, {size}, {row.get('population')}")
    return info[exposure_id], info[outcome_id]

# Retained for backward-compatibility with gwas_id_check()
HLA_SNPs = ["rs12722495-A", "rs61839660-C", "rs12722496-G"]


# x = []; y = []
# plt.plot(x, y, label="")
# plt.xlabel("")
# plt.ylabel("")

def plot_mr(x, y, merged, label):
    plt.figure()
    plt.errorbar(merged["beta_exp"], merged["beta_out"],
             xerr=merged["se_exp"], yerr=merged["se_out"],
             fmt='o', color='black', label='SNPs')
    plt.plot(x, y, label="IVW slope")
    plt.xlabel("SNP effect on IL-6 (SD)")
    plt.ylabel("SNP effect on T1D (log odds)")
    plt.title(label)
    plt.legend()
    plt.savefig(f"scatterplot_{label}.png")
    plt.close()


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


def steiger(merged, n_exp, out_info):
    """Steiger directionality test. An instrument for IL-6 should explain more
    variance in IL-6 than in T1D; if it's the other way round, the SNP probably
    acts on T1D first. Outcome r2 uses the liability-scale formula from
    TwoSampleMR's get_r_from_lor (logit model)."""
    d = merged
    # Exposure (continuous): r2 from the F statistic
    f = (d["beta_exp"] / d["se_exp"]) ** 2
    d["r2_exp"] = f / (f + n_exp - 2)

    # Outcome (case/control): turn the log odds ratio into variance explained
    ncase, ncontrol = out_info["ncase"], out_info["ncontrol"]
    lor, af = d["beta_out"], d["eaf_out"]
    case_frac = ncase / (ncase + ncontrol)
    # The GWAS has more cases than the population does, so reweight the allele frequency
    af_case = af * np.exp(lor) / (1 - af + af * np.exp(lor))
    af_control = (af - case_frac * af_case) / (1 - case_frac)
    pop_af = T1D_PREVALENCE * af_case + (1 - T1D_PREVALENCE) * af_control
    vg = lor ** 2 * pop_af * (1 - pop_af)
    d["r2_out"] = vg / (vg + np.pi ** 2 / 3)

    z = (np.arctanh(np.sqrt(d["r2_exp"])) - np.arctanh(np.sqrt(d["r2_out"]))) / \
        np.sqrt(1 / (n_exp - 3) + 1 / (ncase + ncontrol - 3))
    d["steiger_p"] = 2 * (1 - stats.norm.cdf(np.abs(z)))
    d["steiger_ok"] = d["r2_exp"] > d["r2_out"]

    print("\n----- Steiger directionality -----")
    print(d[["rsid", "r2_exp", "r2_out", "steiger_ok", "steiger_p"]].to_string(index=False))
    n_bad = int((~d["steiger_ok"]).sum())
    print(f"  {n_bad} of {len(d)} instrument(s) explain more variance in T1D than in IL-6.")
    return d


def overlap_check(exp_info, out_info, mean_f):
    """Sample overlap biases two-sample MR toward the confounded observational
    estimate by roughly (overlap fraction) x (1 / mean F). The API can't tell us
    which cohorts are shared, so report the worst case and what to check by hand."""
    print("\n----- Sample overlap -----")
    print(f"  Exposure: {exp_info.get('author')} {exp_info.get('year')}, N={exp_info.get('sample_size')}, "
          f"{exp_info.get('population')}  (PMID {exp_info.get('pmid')})")
    print(f"  Outcome:  {out_info.get('author')} {out_info.get('year')}, N={out_info.get('sample_size')}, "
          f"{out_info.get('population')}  (PMID {out_info.get('pmid')})")
    print(f"  Mean F = {mean_f:.1f}. Worst case (every exposure participant is also in the outcome GWAS):")
    print(f"  relative bias toward the confounded estimate ~ 1/F = {100 / mean_f:.1f}%.")
    print("  TODO by hand: compare the cohort lists in the two papers to get the real overlap fraction.")


def _report(label, method, beta, se, p, nsnp):
    """Print the estimate as an odds ratio and return a row for the summary table."""
    lo, hi = np.exp(beta - 1.96 * se), np.exp(beta + 1.96 * se)
    print(f"  {method}: OR per SD higher IL-6 = {np.exp(beta):.3f} "
          f"(95% CI {lo:.3f} to {hi:.3f}), p = {p:.4g}, {nsnp} SNP(s)")
    return {"run": label, "method": method, "nsnp": nsnp, "beta": beta, "se": se,
            "or_per_sd": np.exp(beta), "ci_low": lo, "ci_high": hi, "p": p}


def mr_analysis(run):
    label, exposure_id, pval, mode = run["label"], run["exposure"], run["pval"], run["mode"]
    rows = []
    try:
        print(f"\n\n========== RUN: {label} (exposure {exposure_id}, p < {pval:g}, {mode}) ==========")
        # Standardizing column mapping for ieugwaspy
        col_map = {
            "ea": "effect_allele",
            "nea": "other_allele",
            "beta": "beta",
            "se": "se",
            "rsid": "rsid"
        }

        cache_file = _cache_path(label, exposure_id, OUTCOME_ID)
        merged = None

        try:
            exp_info, out_info = describe_datasets(exposure_id, OUTCOME_ID)

            instruments = get_instruments(exposure_id, pval)
            if not instruments:
                print(f"No instruments pass p < {pval:g} in {exposure_id}. Nothing to analyse for this run.")
                return rows
            print(f"Fetching IL-6 exposure data for {len(instruments)} instruments...")
            exposure_raw = _batch_associations(instruments, [exposure_id])

            if not exposure_raw or isinstance(exposure_raw, dict):
                raise RuntimeError(f"No exposure data found. Check GWASid or API token. "
                                    f"({type(exposure_raw)}, {exposure_raw})")

            exposure = pd.DataFrame(exposure_raw).rename(columns=col_map)

            if EXCLUDE_MHC:
                in_mhc = _in_region(exposure, MHC_REGION)
                if in_mhc.any():
                    print(f"Excluding {int(in_mhc.sum())} MHC instrument(s): {exposure.loc[in_mhc, 'rsid'].tolist()}")
                    exposure = exposure[~in_mhc].copy()

            if mode == "cis":
                in_cis = pd.Series(False, index=exposure.index)
                for region in CIS_REGIONS.values():
                    in_cis |= _in_region(exposure, region)
                print(f"Cis-only: keeping {int(in_cis.sum())} of {len(exposure)} instrument(s) "
                      f"inside {list(CIS_REGIONS)}.")
                exposure = exposure[in_cis].copy()

            snps = exposure["rsid"].dropna().unique().tolist()
            print(f"Found {len(snps)} genetic instruments.")
            if not snps:
                print("No instruments left after filtering. Nothing to analyse for this run.")
                return rows

            # 2. Get Outcome Data
            outcome_raw = _batch_associations(snps, [OUTCOME_ID])
            if not outcome_raw:
                raise RuntimeError("No outcome data found.")

            outcome = pd.DataFrame(outcome_raw).rename(columns=col_map)

            # 3. Merge and Harmonize
            merged = pd.merge(
                exposure[["rsid", "chr", "position", "beta", "se", "effect_allele", "other_allele", "eaf", "n"]].rename(
                    columns={"beta": "beta_exp", "se": "se_exp", "eaf": "eaf_exp", "n": "n_exp"}
                ),
                outcome[["rsid", "beta", "se", "effect_allele", "other_allele", "eaf"]].rename(
                    columns={"beta": "beta_out", "se": "se_out", "eaf": "eaf_out"}
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
            merged.loc[swap, "eaf_out"] = 1 - merged.loc[swap, "eaf_out"]

            if merged.empty:
                raise RuntimeError("No harmonized SNPs left.")

            # F-statistic per instrument
            merged["F_stat"] = (merged["beta_exp"] / merged["se_exp"]) ** 2
            merged["weak"] = merged["F_stat"] < 10

            n_exp = merged["n_exp"].fillna(exp_info.get("sample_size")).astype(float)
            merged = steiger(merged, n_exp, out_info)
            overlap_check(exp_info, out_info, merged["F_stat"].mean())

            merged.to_csv(cache_file, index=False)
            print(f"Cached harmonized data -> {cache_file}")

        except Exception as e:
            if os.path.exists(cache_file):
                print(f"Live API call failed ({e}); falling back to cached harmonized data ({cache_file}).")
                merged = pd.read_csv(cache_file)
            else:
                raise

        print("\n----- Instrument Strength (F-statistics) -----")
        print(merged[["rsid", "chr", "position", "beta_exp", "se_exp", "F_stat", "weak"]].to_string(index=False))
        weak_count = int(merged["weak"].sum())
        if weak_count:
            print(f"WARNING: {weak_count} instrument(s) have F < 10 (weak instrument bias risk).")

        if STEIGER_FILTER and not merged["steiger_ok"].all():
            dropped = merged.loc[~merged["steiger_ok"], "rsid"].tolist()
            print(f"Steiger filter: dropping {len(dropped)} instrument(s): {dropped}")
            merged = merged[merged["steiger_ok"]].copy()
            if merged.empty:
                print("No instruments left after Steiger filtering.")
                return rows

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

            print(f"\n----- MR Results (Wald ratio, single SNP {row['rsid']}) -----")
            print(f"Causal Estimate (Beta): {wald_beta:.4f}")
            print(f"Standard Error: {wald_se:.4f}")
            rows.append(_report(label, "Wald ratio", wald_beta, wald_se, p_val, 1))
            print("  Egger / weighted median / Cochran's Q need 3+ instruments -- not run.")
            return rows

        # 4. Regression (IVW)
        X = merged["beta_exp"]
        y = merged["beta_out"]
        # Inverse Variance Weighting
        weights = 1 / (merged["se_out"]**2)
        wls = sm.WLS(y, X, weights=weights).fit()

        print("\n----- MR Results (IVW via WLS) -----")
        print(wls.summary())
        rows.append(_report(label, "IVW", wls.params.iloc[0], wls.bse.iloc[0], wls.pvalues.iloc[0], num_snps))

        if num_snps >= 3:
            slope, slope_se, slope_p, _, _, _ = mr_egger(merged)
            rows.append(_report(label, "MR-Egger", slope, slope_se, slope_p, num_snps))
            median, median_se, median_p = weighted_median(merged)
            rows.append(_report(label, "Weighted median", median, median_se, median_p, num_snps))
        else:
            print("  Egger / weighted median need 3+ instruments -- not run.")
        cochrans_q(merged, wls.params.iloc[0])

        # print(wls.params)
        x_line = np.linspace(merged["beta_exp"].min(), merged["beta_exp"].max())
        plot_mr(x_line, wls.params.iloc[0] * x_line, merged, label)
        return rows

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
    summary = []
    for run in RUNS:
        summary.extend(mr_analysis(run))

    print("\n\n========== SUMMARY (odds ratio for T1D per SD higher IL-6) ==========")
    summary = pd.DataFrame(summary)
    print(summary.round(4).to_string(index=False))
    summary.to_csv(os.path.join(CACHE_DIR, "mr_summary.csv"), index=False)
    print("Done")
