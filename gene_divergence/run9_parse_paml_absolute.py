import os
import re
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, mannwhitneyu, binomtest
import matplotlib.pyplot as plt
import seaborn as sns

# --- Configuración ---
INSIDE_DIR = "bp/paml_inputs"
OUTSIDE_DIR = "out/paml_inputs"
OUTPUT_FILE_INSIDE = "results_omega_inside_acc.csv"
OUTPUT_FILE_OUTSIDE = "results_omega_outside_acc.csv"
PLOTS_PDF = "omega_distributions_acc.pdf"

DS_MIN_THRESHOLD = 0.001
DN_MIN_THRESHOLD = 0
DS_MAX_THRESHOLD = 3

GENE_LIST_INSIDE_AGGREGATE = "gene_ids_aggregate_inside.txt"
GENE_LIST_OUTSIDE_AGGREGATE = "gene_ids_aggregate_outside.txt"
GENE_LIST_INSIDE_STRICT = "gene_ids_strict_inside.txt"
GENE_LIST_OUTSIDE_STRICT = "gene_ids_strict_outside.txt"

# ---------------------

def safe_float(x):
    try:
        return float(x)
    except:
        return 0

def extract_paml_metrics(mlc_file_path):
    metrics = {}
    with open(mlc_file_path, "r") as f:
        text = f.read()

    m = re.search(r"dN\s*&\s*dS\s*for\s*each\s*branch(.+?)tree length", text, re.S)
    if not m:
        print(f"No se encontró la tabla en {mlc_file_path}")
        return None

    table = m.group(1)
    pattern = re.compile(
        r"(\d+\.\.\d+)\s+"      # branch
        r"([-.]?\d*\.?\d+)\s+"  # t
        r"(\d*\.?\d+)\s+"       # N
        r"(\d*\.?\d+)\s+"       # S
        r"([-.]?\d*\.?\d+)\s+"  # omega
        r"([-.]?\d*\.?\d+)\s+"  # dN
        r"([-.]?\d*\.?\d+)",    # dS
        re.M
    )
    matches = pattern.findall(table)
    if not matches:
        print(f"No se pudieron leer líneas válidas en la tabla de {mlc_file_path}")
        return None

    total_N_sites = total_S_sites = total_dN_changes = total_dS_changes = 0
    branches = []
    for branch, t, N, S, omega, dN, dS in matches:
        N = safe_float(N)
        S = safe_float(S)
        dN = safe_float(dN)
        dS = safe_float(dS)
        omega = safe_float(omega)

        branches.append({
            "branch": branch,
            "N_sites": N,
            "S_sites": S,
            "dN": dN,
            "dS": dS,
            "omega": omega
        })

        total_N_sites += N
        total_S_sites += S
        total_dN_changes += dN * N
        total_dS_changes += dS * S

    dN_agg = total_dN_changes / total_N_sites if total_N_sites > 0 else 0
    dS_agg = total_dS_changes / total_S_sites if total_S_sites > 0 else 0
    omega_agg = dN_agg / dS_agg if dS_agg > 0 else 0

    return {
        "branches": branches,
        "N_sites": total_N_sites,
        "S_sites": total_S_sites,
        "N_changes": total_dN_changes,
        "S_changes": total_dS_changes,
        "dN_agg": dN_agg,
        "dS_agg": dS_agg,
        "omega_agg": omega_agg
    }

def process_directory(input_dir, output_file):
    results = []
    for filename in os.listdir(input_dir):
        if filename.endswith(".mlc"):
            gene_id = filename.replace(".mlc", "")
            file_path = os.path.join(input_dir, filename)
            metrics = extract_paml_metrics(file_path)
            if metrics:
                results.append({
                    "gene_id": gene_id,
                    "N_sites": metrics["N_sites"],
                    "S_sites": metrics["S_sites"],
                    "N_changes": metrics["N_changes"],
                    "S_changes": metrics["S_changes"],
                    "dN": metrics["dN_agg"],
                    "dS": metrics["dS_agg"],
                    "omega": metrics["omega_agg"]
                })
    if results:
        df = pd.DataFrame(results)
        df.to_csv(output_file, index=False)
        print(f"✅ Successfully wrote {len(results)} raw results to {output_file}")
        return df
    else:
        print(f"No valid results found in {input_dir}")
        return pd.DataFrame()

def get_aggregate_table(df_inside, df_outside):
    inside_totals = {
        "CDS_number": len(df_inside),
        "N_sites": df_inside["N_sites"].sum(),
        "S_sites": df_inside["S_sites"].sum(),
        "N_changes": df_inside["N_changes"].sum(),
        "S_changes": df_inside["S_changes"].sum()
    }
    outside_totals = {
        "CDS_number": len(df_outside),
        "N_sites": df_outside["N_sites"].sum(),
        "S_sites": df_outside["S_sites"].sum(),
        "N_changes": df_outside["N_changes"].sum(),
        "S_changes": df_outside["S_changes"].sum()
    }

    inside_totals["dN_agg"] = inside_totals["N_changes"] / inside_totals["N_sites"]
    inside_totals["dS_agg"] = inside_totals["S_changes"] / inside_totals["S_sites"]
    inside_totals["omega_agg"] = inside_totals["dN_agg"] / inside_totals["dS_agg"]

    outside_totals["dN_agg"] = outside_totals["N_changes"] / outside_totals["N_sites"]
    outside_totals["dS_agg"] = outside_totals["S_changes"] / outside_totals["S_sites"]
    outside_totals["omega_agg"] = outside_totals["dN_agg"] / outside_totals["dS_agg"]

    table = pd.DataFrame([
        ["Inside Inversion", inside_totals["CDS_number"], inside_totals["N_sites"], inside_totals["S_sites"], inside_totals["N_changes"], inside_totals["S_changes"], inside_totals["dN_agg"], inside_totals["dS_agg"], inside_totals["omega_agg"]],
        ["Outside Inversion", outside_totals["CDS_number"], outside_totals["N_sites"], outside_totals["S_sites"], outside_totals["N_changes"], outside_totals["S_changes"], outside_totals["dN_agg"], outside_totals["dS_agg"], outside_totals["omega_agg"]]
    ], columns=["Group", "CDS_number", "N_sites", "S_sites", "N_changes", "S_changes", "dN_agg", "dS_agg", "omega_agg"])
    return table, inside_totals, outside_totals

def chi2_test_aggregates(inside_totals, outside_totals):
    table = np.array([
        [int(round(inside_totals["N_changes"])), int(round(inside_totals["S_changes"]))],
        [int(round(outside_totals["N_changes"])), int(round(outside_totals["S_changes"]))]
    ])
    chi2, p, dof, expected = chi2_contingency(table)
    print("\n--- Chi-Squared Test on Aggregated Changes ---")
    print("Contingency Table (N_changes vs S_changes):")
    print(table)
    print(f"Chi2 = {chi2:.4f}, p = {p:.5f}")
    print(f"Interpretation: {'Significant difference' if p<0.05 else 'No significant difference'}")

def mann_whitney_aggregates(df_inside, df_outside):
    df_inside = df_inside.copy()
    df_outside = df_outside.copy()
    df_inside["omega_agg"] = df_inside["dN"] / df_inside["dS"].replace(0, np.nan)
    df_outside["omega_agg"] = df_outside["dN"] / df_outside["dS"].replace(0, np.nan)

    metrics = ["dN", "dS", "omega_agg"]
    print("\n--- Mann-Whitney U Test on Aggregated Metrics ---")
    for metric in metrics:
        inside_vals = df_inside[metric].dropna()
        outside_vals = df_outside[metric].dropna()
        stat, p_two = mannwhitneyu(inside_vals, outside_vals, alternative="two-sided", method="auto")
        stat_less, p_less = mannwhitneyu(inside_vals, outside_vals, alternative="less", method="auto")
        stat_greater, p_greater = mannwhitneyu(inside_vals, outside_vals, alternative="greater", method="auto")
        def q(x):
            return np.quantile(x, [0.1, 0.25, 0.5, 0.75, 0.9])

        q_in = q(inside_vals)
        q_out = q(outside_vals)
        
        print(f"\nMetric: {metric}")
        print(f"  Inside  quantiles (10%,25%,50%,75%,90%): "
              f"{', '.join(f'{v:.4f}' for v in q_in)}")
        print(f"  Outside quantiles (10%,25%,50%,75%,90%): "
              f"{', '.join(f'{v:.4f}' for v in q_out)}")
        print(f"  Two-sided p: {p_two:.5f}")
        print(f"  One-sided Inside<Outside p: {p_less:.5f}")
        print(f"  One-sided Inside>Outside p: {p_greater:.5f}")

def bootstrap_aggregate_metrics(df_inside, df_outside, n_boot=10000, seed=123):
    np.random.seed(seed)

    def compute_agg(df):
        N_sites = df["N_sites"].sum()
        S_sites = df["S_sites"].sum()
        N_changes = df["N_changes"].sum()
        S_changes = df["S_changes"].sum()

        dN = N_changes / N_sites if N_sites > 0 else np.nan
        dS = S_changes / S_sites if S_sites > 0 else np.nan
        omega = dN / dS if dS > 0 else np.nan
        return dN, dS, omega

    # valores observados
    obs_dN_in, obs_dS_in, obs_omega_in = compute_agg(df_inside)
    obs_dN_out, obs_dS_out, obs_omega_out = compute_agg(df_outside)

    obs_diff = {
        "dN": obs_dN_in - obs_dN_out,
        "dS": obs_dS_in - obs_dS_out,
        "omega": obs_omega_in - obs_omega_out
    }

    diffs = {"dN": [], "dS": [], "omega": []}

    for _ in range(n_boot):
        samp_in = df_inside.sample(len(df_inside), replace=True)
        samp_out = df_outside.sample(len(df_outside), replace=True)

        dN_in, dS_in, om_in = compute_agg(samp_in)
        dN_out, dS_out, om_out = compute_agg(samp_out)

        diffs["dN"].append(dN_in - dN_out)
        diffs["dS"].append(dS_in - dS_out)
        diffs["omega"].append(om_in - om_out)

    ci = {
        k: (np.nanpercentile(diffs[k], 5),
            np.nanpercentile(diffs[k], 95))
        for k in diffs
    }

    return obs_diff, ci

from scipy.stats import binomtest

def enrichment_test_rates(inside_totals, outside_totals):
    """
    Enrichment-style tests for dN and dS
    Outside is treated as background rate
    """

    print("\n--- Enrichment Tests for Substitution Rates (Outside as background) ---")

    tests = {
        "dN": {
            "k_in": inside_totals["N_changes"],
            "n_sites_in": inside_totals["N_sites"],
            "rate_out": outside_totals["N_changes"] / outside_totals["N_sites"]
        },
        "dS": {
            "k_in": inside_totals["S_changes"],
            "n_sites_in": inside_totals["S_sites"],
            "rate_out": outside_totals["S_changes"] / outside_totals["S_sites"]
        }
    }

    for name, t in tests.items():
        res = binomtest(
            k=int(round(t["k_in"])),
            n=int(round(t["n_sites_in"])),
            p=t["rate_out"],
            alternative="two-sided"
        )

        obs_rate = t["k_in"] / t["n_sites_in"]

        print(f"\n{name} enrichment:")
        print(f"  Outside background rate = {t['rate_out']:.6f}")
        print(f"  Inside observed rate    = {obs_rate:.6f}")
        print(f"  Binomial p-value        = {res.pvalue:.5e}")

        if res.pvalue < 0.05:
            print("  Interpretation: Significant enrichment/depletion")
        else:
            print("  Interpretation: No significant enrichment")

def subsample_outside_aggregate_test(
    df_inside,
    df_outside,
    n_iter=10000,
    seed=123
):
    """
    Empirical enrichment test:
    - Inside is fixed
    - Outside is subsampled to same number of genes
    - Aggregated dN, dS, omega are computed
    """

    np.random.seed(seed)

    def compute_agg(df):
        N_sites = df["N_sites"].sum()
        S_sites = df["S_sites"].sum()
        N_changes = df["N_changes"].sum()
        S_changes = df["S_changes"].sum()

        dN = N_changes / N_sites if N_sites > 0 else np.nan
        dS = S_changes / S_sites if S_sites > 0 else np.nan
        omega = dN / dS if dS > 0 else np.nan
        return dN, dS, omega

    # Observed inside aggregates
    dN_in, dS_in, om_in = compute_agg(df_inside)

    n_inside = len(df_inside)

    null_dist = {
        "dN": [],
        "dS": [],
        "omega": []
    }

    for _ in range(n_iter):
        samp_out = df_outside.sample(n_inside, replace=False)
        dN_o, dS_o, om_o = compute_agg(samp_out)

        null_dist["dN"].append(dN_o)
        null_dist["dS"].append(dS_o)
        null_dist["omega"].append(om_o)

    # Convert to arrays
    for k in null_dist:
        null_dist[k] = np.array(null_dist[k])

    # Empirical CIs
    ci = {
        k: (
            np.nanpercentile(null_dist[k], 2.5),
            np.nanpercentile(null_dist[k], 97.5)
        )
        for k in null_dist
    }

    # Empirical p-values (two-sided)
    pvals = {
        k: np.mean(
            np.abs(null_dist[k] - np.nanmean(null_dist[k]))
            >= np.abs(val - np.nanmean(null_dist[k]))
        )
        for k, val in zip(
            ["dN", "dS", "omega"],
            [dN_in, dS_in, om_in]
        )
    }

    obs = {"dN": dN_in, "dS": dS_in, "omega": om_in}

    return obs, ci, pvals, null_dist


if __name__ == "__main__":
    df_inside = process_directory(INSIDE_DIR, OUTPUT_FILE_INSIDE)
    df_outside = process_directory(OUTSIDE_DIR, OUTPUT_FILE_OUTSIDE)

    if len(df_inside)==0 or len(df_outside)==0:
        print("No se encontraron datos para analizar.")
        exit()

    # Filtrar genes según dS y dN máximos
    df_inside_filt = df_inside[(df_inside["dS"]<=DS_MAX_THRESHOLD) & (df_inside["dN"]<=DS_MAX_THRESHOLD)
           & (df_inside["dS"]>=DS_MIN_THRESHOLD) & (df_inside["dN"]>=DN_MIN_THRESHOLD) ]
    df_outside_filt = df_outside[(df_outside["dS"]<=DS_MAX_THRESHOLD) & (df_outside["dN"]<=DS_MAX_THRESHOLD) & (df_outside["dS"]>=DS_MIN_THRESHOLD) & (df_outside["dN"]>=DN_MIN_THRESHOLD)]
    
    # Exportar listas de genes
    df_inside_filt["gene_id"].to_csv(GENE_LIST_INSIDE_AGGREGATE, index=False, header=False)
    df_outside_filt["gene_id"].to_csv(GENE_LIST_OUTSIDE_AGGREGATE, index=False, header=False)

    # Tabla de agregados
    table, inside_totals, outside_totals = get_aggregate_table(df_inside_filt, df_outside_filt)
    print("\n--- Aggregated Table ---")
    print(table.to_string(index=False))

    # Chi-cuadrado sobre agregados
    chi2_test_aggregates(inside_totals, outside_totals)

    # Enrichment explícito (outside = background)
    enrichment_test_rates(inside_totals, outside_totals)
    
    # Mann-Whitney sobre agregados por gen
    mann_whitney_aggregates(df_inside_filt, df_outside_filt)

    print("\n--- Subsampling Enrichment Test (Outside as empirical null) ---")

    obs_sub, ci_sub, pvals_sub, null_dist = subsample_outside_aggregate_test(
    df_inside_filt,
    df_outside_filt,
    n_iter=10000
    )

    for metric in ["dN", "dS", "omega"]:
        lo, hi = ci_sub[metric]
        print(
        f"{metric}_agg Inside = {obs_sub[metric]:.6f} "
        f"[Outside subsample 95% CI: {lo:.6f}, {hi:.6f}] "
        f"empirical p = {pvals_sub[metric]:.5f}"
        )   


    #Bootstrap test por omega, dS and dN por gen

    print("\n--- Bootstrap on Aggregated Metrics ---")
    obs_diff, ci = bootstrap_aggregate_metrics(
    df_inside_filt,
    df_outside_filt,
    n_boot=10000
    )

print(f"Observed dN_agg difference (Inside - Outside): {obs_diff['dN']:.6f}")
print(f"Observed dS_agg difference (Inside - Outside): {obs_diff['dS']:.6f}")
print(f"Observed omega_agg difference (Inside - Outside): {obs_diff['omega']:.6f}")
print("\n--- Bootstrap 90% Confidence Intervals (Aggregated Metrics) ---")

for metric in ["dN", "dS", "omega"]:
    lo, hi = ci[metric]
    print(
        f"{metric}_agg difference (Inside - Outside): "
        f"{obs_diff[metric]:.6f} "
        f"[95% CI: {lo:.6f}, {hi:.6f}]"
    )


