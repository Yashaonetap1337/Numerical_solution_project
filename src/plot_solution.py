"""
Visualization for CJ detonation verification (Test 2).

Plots:
  Figure 1 — x_front(t): detonation front position vs time.
             Reference line with slope D_CJ is overlaid.

  Figure 2 — Final profiles P(x), rho(x), W(x) at t = t_final.
             Vertical dashed line marks the detected front position.
             CJ reference values are marked with horizontal lines.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

BASE = os.path.dirname(os.path.abspath(__file__))

# ── CJ reference parameters (must match config.txt) ──────────────────────────
GAMMA   = 3.0
RHO0    = 1.84        # g/cm³  — unreacted explosive
D_CJ    = 0.88        # cm/µs
P_CJ    = D_CJ**2 * RHO0 / (GAMMA + 1)   # ≈ 0.3564 Mbar
U_CJ    = D_CJ / (GAMMA + 1)              # ≈ 0.22  cm/µs
RHO_CJ  = RHO0 * (GAMMA + 1) / GAMMA     # ≈ 2.453 g/cm³

# ── Load data ─────────────────────────────────────────────────────────────────
front_csv = os.path.join(BASE, "results", "front_position.csv")
sol_csv   = os.path.join(BASE, "results", "numerical_solution.csv")

if not os.path.exists(front_csv):
    sys.exit(f"Not found: {front_csv}  (run the simulation first)")
if not os.path.exists(sol_csv):
    sys.exit(f"Not found: {sol_csv}")

front = pd.read_csv(front_csv)
sol   = pd.read_csv(sol_csv)

# quasi-1D: keep unique x, sorted
sol1d = sol.drop_duplicates(subset="x").sort_values("x").reset_index(drop=True)
has_omega = "omega" in sol1d.columns

# ── Figure 1: x_front(t) only (no pressure) ────────────────────────────────────
fig1, ax1 = plt.subplots(figsize=(8, 5))
fig1.suptitle("Detonation front tracking", fontsize=13)

ax1.scatter(front["t"], front["x_front"], s=18, color="tab:blue",
            label="$x_{front}(t)$ simulated", zorder=3)

# Reference line D_CJ through first data point
t0, x0 = front["t"].iloc[0], front["x_front"].iloc[0]
t_ref = np.array([front["t"].min(), front["t"].max()])
ax1.plot(t_ref, x0 + D_CJ * (t_ref - t0), "--", color="tab:red", linewidth=1.5,
         label=f"$D_{{CJ}}$ = {D_CJ} cm/µs  (reference)")

# Estimate measured D from linear fit (at least 3 points)
if len(front) >= 3:
    coeffs = np.polyfit(front["t"], front["x_front"], 1)
    D_meas = coeffs[0]
    ax1.set_title(f"$D_{{measured}}$ = {D_meas:.4f} cm/µs  "
                  f"(error {abs(D_meas - D_CJ)/D_CJ*100:.1f}%  vs D_CJ={D_CJ})",
                  fontsize=10)

ax1.set_xlabel("t, µs")
ax1.set_ylabel("$x_{front}$, cm")
ax1.tick_params(axis="y")
ax1.legend(loc="upper left")
ax1.grid(True, linestyle="--", alpha=0.4)

fig1.tight_layout()
out1 = os.path.join(BASE, "results", "front_tracking.png")
fig1.savefig(out1, dpi=150)
print(f"Saved: {out1}")

# ── Figure 2: final profiles P(x), rho(x), W(x) ─────────────────────────────
n_rows = 3 if has_omega else 2
fig2, axes = plt.subplots(n_rows, 1, figsize=(10, 3.5 * n_rows), sharex=True)
fig2.suptitle("Final-time profiles", fontsize=13)

# Detect front in final profile: cell with max P
idx_front = sol1d["p"].idxmax()
x_front_final = sol1d.loc[idx_front, "x"]
P_max_final   = sol1d.loc[idx_front, "p"]

# — Pressure —
ax = axes[0]
ax.plot(sol1d["x"], sol1d["p"], color="tab:red", linewidth=1.5, label="P(x)")
ax.axvline(x_front_final, color="gray", linestyle="--", linewidth=0.8,
           label=f"front x = {x_front_final:.4f} cm")
ax.annotate(f"$P_{{max}}$ = {P_max_final:.4f}", xy=(x_front_final, P_max_final),
            xytext=(x_front_final + 0.02, P_max_final),
            arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
            fontsize=8)
ax.set_ylabel("P, Mbar")
ax.legend(fontsize=8)
ax.grid(True, linestyle="--", alpha=0.4)

# — Density —
ax = axes[1]
ax.plot(sol1d["x"], sol1d["rho"], color="tab:blue", linewidth=1.5, label="ρ(x)")
ax.axhline(RHO_CJ, color="tab:blue", linestyle=":", linewidth=1,
           label=f"$ρ_{{CJ}}$ = {RHO_CJ:.4f} g/cm³")
ax.axhline(RHO0,   color="tab:cyan", linestyle=":", linewidth=1,
           label=f"$ρ_0$ = {RHO0} g/cm³")
ax.axvline(x_front_final, color="gray", linestyle="--", linewidth=0.8)
ax.set_ylabel("ρ, g/cm³")
ax.legend(fontsize=8)
ax.grid(True, linestyle="--", alpha=0.4)

# — Omega (W) — only if saved —
if has_omega:
    ax = axes[2]
    ax.plot(sol1d["x"], sol1d["omega"], color="tab:green", linewidth=1.5, label="W(x)")
    ax.axvline(x_front_final, color="gray", linestyle="--", linewidth=0.8)
    ax.set_ylabel("W (unreacted fraction)")
    ax.set_ylim(-0.05, 1.1)
    ax.axhline(0, color="k", linewidth=0.5)
    ax.axhline(1, color="k", linewidth=0.5, linestyle=":")
    ax.legend(fontsize=8)
    ax.grid(True, linestyle="--", alpha=0.4)
else:
    axes[-1].text(0.5, 0.5, "omega not saved — recompile with Mader method",
                  ha="center", va="center", transform=axes[-1].transAxes,
                  fontsize=9, color="gray")

axes[-1].set_xlabel("x, cm")
fig2.tight_layout()
out2 = os.path.join(BASE, "results", "solution_profiles.png")
fig2.savefig(out2, dpi=150)
print(f"Saved: {out2}")

# ── Console summary ───────────────────────────────────────────────────────────
print("\n── Verification summary ──────────────────────────────────")
if len(front) >= 3:
    print(f"  D_measured = {D_meas:.4f} cm/µs   D_CJ = {D_CJ:.4f}   "
          f"error = {abs(D_meas-D_CJ)/D_CJ*100:.2f}%")
print(f"  Front position (final) = {x_front_final:.4f} cm")
if has_omega:
    w_behind = sol1d.loc[sol1d["x"] < x_front_final, "omega"]
    if not w_behind.empty:
        print(f"  W behind front: min={w_behind.min():.4f}  mean={w_behind.mean():.4f}  "
              f"(target → 0)")
print("─────────────────────────────────────────────────────────\n")

plt.show()