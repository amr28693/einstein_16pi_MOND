"""
MOND as Bose-Einstein Occupation of Entropy Degrees of Freedom
Rodriguez 2026 

The MOND interpolation function:
    g_obs = g_bar / (1 - exp(-θ))     where θ = sqrt(g_bar / a₀)

can be rewritten as:
    g_obs / g_bar = 1 + n_BE(θ)       where n_BE = 1 / (exp(θ) - 1)

n_BE is the Bose-Einstein occupation number.

USAGE: python3 mond_bose_einstein.py
       python3 mond_bose_einstein.py path/to/rar_galaxies.txt
"""
import numpy as np
import sys, os

def parse_rar_file(filepath):
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line[0].isalpha() or line.startswith('=') or line.startswith('---'):
                continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    vals = [float(parts[i]) for i in range(4)]
                    if -14 < vals[0] < -6 and -14 < vals[2] < -6:
                        data.append(vals)
                except (ValueError, IndexError):
                    continue
    return np.array(data)

# ============================================================
# LOAD DATA
# ============================================================
script_dir = os.path.dirname(os.path.abspath(__file__))
filepath = sys.argv[1] if len(sys.argv) > 1 else os.path.join(script_dir, 'rar_galaxies.txt')

if not os.path.exists(filepath):
    print(f"ERROR: File not found: {filepath}")
    sys.exit(1)

print(f"Reading: {filepath}")
data = parse_rar_file(filepath)
N = len(data)
print(f"Parsed {N} data points\n")

gbar_log = data[:, 0]
gbar_err = data[:, 1]
gobs_log = data[:, 2]
gobs_err = data[:, 3]

gbar = 10**gbar_log
gobs = 10**gobs_log
a0 = 1.2e-10
log_a0 = np.log10(a0)

# ============================================================
# THE KEY COMPUTATION
# ============================================================

# Entropy parameter
theta = np.sqrt(gbar / a0)

# Observed mass discrepancy
D_obs = gobs / gbar

# Bose-Einstein prediction
n_BE = 1.0 / (np.exp(theta) - 1.0)  # occupation number
D_pred = 1.0 + n_BE                  # predicted mass discrepancy

# Observed "occupation number" (dark matter as excitation)
n_obs = D_obs - 1.0  # the "dark matter" contribution

print(f"{'='*65}")
print(f"  MOND AS BOSE-EINSTEIN DISTRIBUTION")
print(f"{'='*65}")
print(f"")
print(f"  MOND:  g_obs/g_bar = 1/(1 - exp(-theta))")
print(f"       = 1 + 1/(exp(theta) - 1)")
print(f"       = 1 + n_BE(theta)")
print(f"")
print(f"  where theta = sqrt(g_bar / a0)")
print(f"  and n_BE = Bose-Einstein occupation number")
print(f"")

# ============================================================
# COMPARISON: PREDICTED vs OBSERVED
# ============================================================

# Residuals
residual = np.log10(D_obs) - np.log10(D_pred)

print(f"{'='*65}")
print(f"  PREDICTED vs OBSERVED MASS DISCREPANCY")
print(f"{'='*65}")
print(f"")
print(f"  Residual = log10(D_obs) - log10(D_pred)")
print(f"  Mean residual:   {np.mean(residual):+.4f} dex")
print(f"  Median residual: {np.median(residual):+.4f} dex")
print(f"  Std residual:    {np.std(residual):.4f} dex")
print(f"  RMS residual:    {np.sqrt(np.mean(residual**2)):.4f} dex")
print(f"")
print(f"  For comparison:")
print(f"  RAR intrinsic scatter (McGaugh+2016): 0.13 dex")
print(f"  Our residual std:                     {np.std(residual):.2f} dex")

if np.std(residual) < 0.15:
    print(f"  --> CONSISTENT: residuals within known RAR scatter")
else:
    print(f"  --> TENSION: residuals exceed known RAR scatter")

# ============================================================
# OCCUPATION NUMBER ANALYSIS
# ============================================================

print(f"\n{'='*65}")
print(f"  OCCUPATION NUMBER n = D - 1 ('DARK MATTER' AS EXCITATION)")
print(f"{'='*65}")

# Bin by theta
theta_bins = np.array([0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.693, 1.0, 1.5, 2.0, 3.0, 5.0])
print(f"\n{'theta':>8} | {'n_BE(pred)':>10} | {'n_obs(med)':>10} | {'N':>5} | {'note':>20}")
print(f"{'-'*8}-+-{'-'*10}-+-{'-'*10}-+-{'-'*5}-+-{'-'*20}")

for i in range(len(theta_bins) - 1):
    lo, hi = theta_bins[i], theta_bins[i+1]
    mask = (theta >= lo) & (theta < hi)
    n = mask.sum()
    if n < 2:
        continue
    mid = (lo + hi) / 2
    n_be_mid = 1.0 / (np.exp(mid) - 1.0)
    n_obs_med = np.median(n_obs[mask])
    
    note = ""
    if abs(mid - np.log(2)) < 0.1:
        note = "<-- theta = ln(2)"
    elif abs(n_be_mid - 1.0) < 0.15:
        note = "<-- n_BE = 1"
    elif n_be_mid < 0.1:
        note = "baryons dominate"
    elif n_be_mid > 5:
        note = "DM dominates"
    
    print(f"{mid:>8.3f} | {n_be_mid:>9.3f}  | {n_obs_med:>9.3f}  | {n:>5} | {note:>20}")

# ============================================================
# THE CROSSOVER: n_BE = 1
# ============================================================

print(f"\n{'='*65}")
print(f"  THE CROSSOVER: n_BE = 1")
print(f"{'='*65}")
print(f"")
print(f"  n_BE(theta) = 1  when  exp(theta) - 1 = 1")
print(f"                    when  exp(theta) = 2")
print(f"                    when  theta = ln(2) = {np.log(2):.6f}")
print(f"")
print(f"  theta = sqrt(g_bar/a0) = ln(2)")
print(f"  => g_bar/a0 = (ln 2)^2 = {np.log(2)**2:.6f}")
print(f"")
print(f"  At this point:")
print(f"    n_BE = 1  (exactly one quantum of occupation)")
print(f"    D = 1 + n_BE = 2  (mass discrepancy = 2)")
print(f"    g_obs = 2 * g_bar  (PG velocity = observed velocity)")
print(f"")
print(f"  Above: n_BE < 1, entropy sector mostly empty, baryons enough")
print(f"  Below: n_BE > 1, entropy sector occupied, 'dark matter' appears")
print(f"  AT:    n_BE = 1, exactly one bit of entropy excitation")


# ============================================================
# ADDITIONAL RESULTS
# ============================================================

results = {
    'N': N,
    'residual_mean': float(np.mean(residual)),
    'residual_std': float(np.std(residual)),
    'residual_rms': float(np.sqrt(np.mean(residual**2))),
    'crossover_theta': float(np.log(2)),
    'crossover_g_over_a0': float(np.log(2)**2),
    'n_BE_at_crossover': 1.0
}

import json
outfile = os.path.join(script_dir, 'bose_einstein_results.json')
with open(outfile, 'w') as f:
    json.dump(results, f, indent=2)
print(f"Results saved: {outfile}")

# Also save the full data for plotting
plotfile = os.path.join(script_dir, 'bose_einstein_plotdata.csv')
with open(plotfile, 'w') as f:
    f.write("theta,n_BE_pred,n_obs,D_pred,D_obs,gbar_log,gobs_log\n")
    for i in range(N):
        f.write(f"{theta[i]:.6f},{n_BE[i]:.6f},{n_obs[i]:.6f},"
                f"{D_pred[i]:.6f},{D_obs[i]:.6f},"
                f"{gbar_log[i]:.4f},{gobs_log[i]:.4f}\n")
print(f"Plot data saved: {plotfile}")
print(f"\nRan Successfully.")
