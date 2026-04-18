# =============================================================================
# Dynamic Modelling of an Alkaline Water Electrolysis (AWE) Cell
# =============================================================================

import numpy as np                          # numerical arrays and math
import matplotlib.pyplot as plt             # plotting
from scipy.integrate import solve_ivp       # ODE solver for dynamic model
import os                                   # to create folders

# --- Plot style: clean, publication-ready ---
plt.rcParams.update({
    "figure.facecolor":  "white",
    "axes.facecolor":    "#F9F9F9",
    "axes.grid":         True,
    "grid.color":        "white",
    "grid.linewidth":    1.0,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "font.size":         11,
    "axes.titlesize":    13,
    "axes.titleweight":  "bold",
    "lines.linewidth":   2.0,
    "savefig.dpi":       200,
    "savefig.bbox":      "tight",
})

# --- Colours we'll reuse across all plots ---
BLUE   = "#1A5EA8"
RED    = "#C0392B"
GREEN  = "#1E8449"
ORANGE = "#E67E22"
PURPLE = "#7D3C98"
GRAY   = "#555555"

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

R_GAS = 8.314    # universal gas constant  [J mol⁻¹ K⁻¹]
F     = 96485    # Faraday constant        [C mol⁻¹]
V_HHV = 1.481    # thermoneutral voltage   [V]

# =============================================================================
# EQUATION 1 — Reversible voltage  E_rev(T)
# =============================================================================
# Minimum voltage needed to split water at a given temperature.
# Derived from Gibbs free energy: H2O → H2 + ½O2
# Empirical fit valid 20–90°C  (Ulleberg, 2003)

def E_rev(T_C):
    T = T_C + 273.15                         # convert Celsius → Kelvin
    return (1.5184
            - 1.5421e-3 * T
            + 9.523e-5  * T * np.log(T)
            + 9.84e-8   * T**2)

# =============================================================================
# EQUATION 2 — Activation overpotential  eta_act(I, T)
# =============================================================================
# Extra voltage to overcome the energy barrier at the electrode surface.
# Butler-Volmer kinetics (Tafel approximation):
#   eta_act = (R·T) / (alpha·F) × ln(I / i0)

def eta_act(I, T_C):
    T      = T_C + 273.15
    i0     = 1.0e-3 * np.exp(0.06 * (T_C - 25.0))  # exchange current [A m⁻²]
    alpha  = 0.40                                     # charge transfer coeff [-]
    I_safe = np.maximum(I, i0 * 1.001)               # avoid log(0) at I≈0
    return (R_GAS * T) / (alpha * F) * np.log(I_safe / i0)

# =============================================================================
# EQUATION 3 — Ohmic overpotential  eta_ohm(I, T)
# =============================================================================
# Voltage lost to resistance — membrane + KOH electrolyte.
# R_ohm decreases with temperature: hot KOH conducts better.

def eta_ohm(I, T_C):
    R_ohm = 0.18 * np.exp(-0.01 * (T_C - 25.0))    # area resistance [Ω m²]
    return R_ohm * I / 1000.0

# =============================================================================
# TOTAL CELL VOLTAGE  V(I, T)
# =============================================================================
# Master equation:  V = E_rev + eta_act + eta_ohm

def V_cell(I, T_C):
    return E_rev(T_C) + eta_act(I, T_C) + eta_ohm(I, T_C)

# =============================================================================
# FIGURE 1 — Polarisation Curve
# =============================================================================

def plot_polarisation_curve():

    I_range = np.linspace(10, 600, 500)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # --- Left panel: V-I at 3 temperatures ---
    ax           = axes[0]
    temperatures = [40,   60,     80  ]
    colors       = [BLUE, ORANGE, RED ]

    for T_C, col in zip(temperatures, colors):
        V = V_cell(I_range, T_C)
        ax.plot(I_range, V, color=col, label=f"T = {T_C} °C")

    ax.axhline(V_HHV, color=GRAY, linestyle="--", linewidth=1.2,
               label=f"Thermoneutral  V_HHV = {V_HHV} V")
    ax.set_xlabel("Current density   i   (A m⁻²)")
    ax.set_ylabel("Cell voltage   V   (V)")
    ax.set_title("Polarisation Curves — Temperature Effect")
    ax.legend()
    ax.set_xlim(0, 620)
    ax.set_ylim(1.1, 2.8)

    # --- Right panel: stacked loss breakdown at 60°C ---
    ax2  = axes[1]
    T_C  = 60
    Erev = np.full_like(I_range, E_rev(T_C))
    Eact = eta_act(I_range, T_C)
    Eohm = eta_ohm(I_range, T_C)

    ax2.stackplot(I_range, Erev, Eact, Eohm,
                  labels=["Reversible  E_rev",
                           "Activation loss  η_act",
                           "Ohmic loss  η_ohm"],
                  colors=["#AED6F1", "#F1948A", "#A9DFBF"],
                  alpha=0.85)
    ax2.set_xlabel("Current density   i   (A m⁻²)")
    ax2.set_ylabel("Voltage contribution  (V)")
    ax2.set_title("Voltage Loss Breakdown   (T = 60 °C)")
    ax2.legend(loc="upper left")
    ax2.set_xlim(0, 620)
    ax2.set_ylim(0, 2.8)

    fig.tight_layout()
    fig.savefig("plots/fig1_polarisation_curve.png")
    print("  Saved → plots/fig1_polarisation_curve.png")
    plt.show()

# =============================================================================
# FIGURE 2 — System Efficiency
# =============================================================================

def plot_efficiency():

    I_range      = np.linspace(10, 600, 500)
    fig, ax      = plt.subplots(figsize=(8, 5))
    temperatures = [40,   60,     80  ]
    colors       = [BLUE, ORANGE, RED ]

    for T_C, col in zip(temperatures, colors):
        eta_V   = (V_HHV / V_cell(I_range, T_C)) * 100   # voltage eff [%]
        eta_sys = eta_V * 0.98                             # system eff  [%]
        ax.plot(I_range, eta_sys, color=col, linestyle="-",  label=f"T = {T_C} °C")
        ax.plot(I_range, eta_V,   color=col, linestyle="--", linewidth=1.0, alpha=0.4)

    ax.axvspan(100, 400, alpha=0.07, color=GREEN, label="Typical operating window")
    ax.axhline(70, color=GRAY, linestyle=":", linewidth=1.0)
    ax.text(570, 70.6, "70% threshold", ha="right", fontsize=9, color=GRAY)
    ax.text(420, 87,
            "Solid = system eff.\nDashed = voltage eff.",
            fontsize=8, color=GRAY,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#CCCCCC"))
    ax.set_xlabel("Current density   i   (A m⁻²)")
    ax.set_ylabel("Efficiency   η   (%)")
    ax.set_title("AWE System Efficiency vs Current Density")
    ax.legend(loc="upper right")
    ax.set_xlim(0, 620)
    ax.set_ylim(55, 95)

    fig.tight_layout()
    fig.savefig("plots/fig2_efficiency.png")
    print("  Saved → plots/fig2_efficiency.png")
    plt.show()

# =============================================================================
# DYNAMIC MODEL — First-order ODE
# =============================================================================
# When current changes, voltage responds gradually — not instantly.
# We model this as:   tau * dV/dt = V_ss(I,T) - V(t)

def dV_dt(t, V, I_before, I_after, T_C, tau, t_step):
    I_now = I_before if t < t_step else I_after
    V_ss  = V_cell(I_now, T_C)
    return [(V_ss - V[0]) / tau]


def simulate_step(I_before, I_after, T_C, tau, t_step=5.0, t_end=40.0):

    t_eval   = np.linspace(0, t_end, 1000)
    V0       = [V_cell(I_before, T_C)]

    solution = solve_ivp(
        dV_dt,
        t_span=(0, t_end),
        y0=V0,
        t_eval=t_eval,
        args=(I_before, I_after, T_C, tau, t_step),
        method="RK45",
        rtol=1e-8,
        atol=1e-10,
    )

    V_ss_line = np.where(t_eval < t_step,
                         V_cell(I_before, T_C),
                         V_cell(I_after,  T_C))

    return {
        "t":           solution.t,
        "V":           solution.y[0],
        "V_ss":        V_ss_line,
        "V_before":    V_cell(I_before, T_C),
        "V_after":     V_cell(I_after,  T_C),
        "delta_V":     V_cell(I_after, T_C) - V_cell(I_before, T_C),
        "settle_time": t_step + 5 * tau,
        "tau":         tau,
        "t_step":      t_step,
        "I_before":    I_before,
        "I_after":     I_after,
    }

# =============================================================================
# FIGURE 3 — Dynamic Voltage Response
# =============================================================================

def plot_dynamic_response():

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # --- Left panel: single step, fully annotated ---
    ax  = axes[0]
    res = simulate_step(I_before=150, I_after=400, T_C=60, tau=5.0)

    ax.plot(res["t"], res["V"],   color=BLUE, linewidth=2.5, label="Cell voltage  V(t)")
    ax.plot(res["t"], res["V_ss"], color=RED,  linewidth=1.4, linestyle="--",
            label="Steady-state  V_ss")
    ax.axvline(res["t_step"], color=GRAY, linestyle=":", linewidth=1.2)

    ax.annotate(
        f"Step input\ni: {res['I_before']} → {res['I_after']} A m⁻²",
        xy=(res["t_step"], res["V_before"]),
        xytext=(res["t_step"] + 1.5, res["V_before"] + 0.05),
        arrowprops=dict(arrowstyle="->", color=GRAY, lw=1.0),
        fontsize=9, color=GRAY
    )

    ax.axvline(res["settle_time"], color=GREEN, linestyle="--", linewidth=1.0, alpha=0.8)
    ax.text(res["settle_time"] + 0.4, res["V_after"] - 0.07,
            f"≈99% settled\nt = {res['settle_time']:.0f}s  (= 5τ)",
            fontsize=9, color=GREEN)

    ax.annotate("",
        xy=(2.0, res["V_after"]), xytext=(2.0, res["V_before"]),
        arrowprops=dict(arrowstyle="<->", color=PURPLE, lw=1.4))
    ax.text(2.5, (res["V_before"] + res["V_after"]) / 2,
            f"ΔV = {res['delta_V']:.3f} V",
            fontsize=9, color=PURPLE, va="center")

    ax.set_xlabel("Time  (s)")
    ax.set_ylabel("Cell voltage  V  (V)")
    ax.set_title("Dynamic Response to Step Load   (τ = 5 s)")
    ax.legend(loc="lower right")
    ax.set_ylim(res["V_before"] - 0.15, res["V_after"] + 0.15)

    # --- Right panel: parametric tau study ---
    ax2  = axes[1]
    taus = [2,    5,      10,     20  ]
    cols = [BLUE, GREEN,  ORANGE, RED ]

    for tau, col in zip(taus, cols):
        r = simulate_step(150, 400, T_C=60, tau=tau)
        ax2.plot(r["t"], r["V"], color=col, label=f"τ = {tau} s")

    ax2.plot(r["t"], r["V_ss"], color=GRAY, linestyle="--",
             linewidth=1.2, label="V_ss  reference")
    ax2.axvline(r["t_step"], color=GRAY, linestyle=":", linewidth=1.0)
    ax2.set_xlabel("Time  (s)")
    ax2.set_ylabel("Cell voltage  V  (V)")
    ax2.set_title("Parametric Study — Effect of Time Constant τ")
    ax2.legend()
    ax2.set_ylim(res["V_before"] - 0.15, res["V_after"] + 0.15)

    fig.tight_layout()
    fig.savefig("plots/fig3_dynamic_response.png")
    print("  Saved → plots/fig3_dynamic_response.png")
    plt.show()

# =============================================================================
# FIGURE 4 — 2D Operating Map
# =============================================================================

def plot_operating_map():

    I_arr = np.linspace(20,  600, 300)
    T_arr = np.linspace(25,  90,  300)
    I_grid, T_grid = np.meshgrid(I_arr, T_arr)

    eff_grid = (V_HHV / V_cell(I_grid, T_grid)) * 100 * 0.98
    V_grid   = V_cell(I_grid, T_grid)
    P_grid   = V_cell(I_grid, T_grid) * I_grid / 1000.0

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # --- Left panel: efficiency colour map + voltage isolines ---
    ax = axes[0]
    cf = ax.contourf(I_grid, T_grid, eff_grid,
                     levels=np.linspace(60, 88, 28), cmap="RdYlGn")
    fig.colorbar(cf, ax=ax, label="System efficiency  η  (%)")

    cs = ax.contour(I_grid, T_grid, V_grid,
                    levels=[1.6, 1.8, 2.0, 2.2, 2.4],
                    colors="white", linewidths=0.8, alpha=0.6)
    ax.clabel(cs, fmt="%.1f V", fontsize=8, colors="white")

    I_opt_list = []
    for T_C in T_arr:
        I_test   = np.linspace(20, 600, 2000)
        eff_test = (V_HHV / V_cell(I_test, T_C)) * 100 * 0.98
        I_opt_list.append(I_test[np.argmax(eff_test)])

    ax.plot(I_opt_list, T_arr, color="white", linestyle="--",
            linewidth=1.8, label="Peak efficiency locus")
    ax.legend(loc="upper right", fontsize=9)
    ax.set_xlabel("Current density   i   (A m⁻²)")
    ax.set_ylabel("Temperature   T   (°C)")
    ax.set_title("Efficiency Map  (η contours + V isolines)")

    # --- Right panel: power density + efficiency contours ---
    ax2 = axes[1]
    cf2 = ax2.contourf(I_grid, T_grid, P_grid, levels=30, cmap="Blues")
    fig.colorbar(cf2, ax=ax2, label="Power density  (kW m⁻²)")

    ax2.contour(I_grid, T_grid, eff_grid,
                levels=[70, 75, 80],
                colors=[RED, ORANGE, GREEN],
                linewidths=1.8)

    from matplotlib.lines import Line2D
    legend_lines = [
        Line2D([0], [0], color=GREEN,  linewidth=1.8, label="η = 80%"),
        Line2D([0], [0], color=ORANGE, linewidth=1.8, label="η = 75%"),
        Line2D([0], [0], color=RED,    linewidth=1.8, label="η = 70%"),
    ]
    ax2.legend(handles=legend_lines, loc="upper left", fontsize=9)
    ax2.set_xlabel("Current density   i   (A m⁻²)")
    ax2.set_ylabel("Temperature   T   (°C)")
    ax2.set_title("Power Density Map + Efficiency Contours")

    fig.tight_layout()
    fig.savefig("plots/fig4_operating_map.png")
    print("  Saved → plots/fig4_operating_map.png")
    plt.show()

# =============================================================================
# SUMMARY — print key numbers to terminal
# =============================================================================

def print_summary(T_C=60.0, I_op=200.0, tau=5.0):

    V_op    = V_cell(I_op, T_C)
    Erev    = E_rev(T_C)
    Eact    = eta_act(I_op, T_C)
    Eohm    = eta_ohm(I_op, T_C)
    eta_V   = (V_HHV / V_op) * 100
    eta_sys = eta_V * 0.98
    P_den   = V_op * I_op / 1000.0
    settle  = 5 * tau

    print("\n" + "="*58)
    print("  AWE CELL — MODEL SUMMARY")
    print("="*58)
    print(f"  Operating temperature      : {T_C} °C")
    print(f"  Operating current density  : {I_op} A m⁻²")
    print(f"  Dynamic time constant      : {tau} s")
    print("-"*58)
    print(f"  Reversible voltage E_rev   : {Erev:.4f}  V")
    print(f"  Activation loss  eta_act   : {Eact:.4f}  V  ({Eact/V_op*100:.1f}% of total)")
    print(f"  Ohmic loss       eta_ohm   : {Eohm:.4f}  V  ({Eohm/V_op*100:.1f}% of total)")
    print(f"  Total cell voltage         : {V_op:.4f}  V")
    print("-"*58)
    print(f"  Voltage efficiency         : {eta_V:.2f} %")
    print(f"  System efficiency (HHV)    : {eta_sys:.2f} %")
    print(f"  Power density              : {P_den:.3f}  kW m⁻²")
    print(f"  Settling time (5τ)         : {settle:.1f}  s")
    print("="*58 + "\n")

# =============================================================================
# MAIN — single RUN block at the very bottom, runs everything once
# =============================================================================

if __name__ == "__main__":

    os.makedirs("plots", exist_ok=True)

    print("\n  Alkaline Water Electrolysis — Dynamic Simulation")
    print("  Running model...\n")

    print_summary(T_C=60.0, I_op=200.0, tau=5.0)

    print("  [1/4] Polarisation curves...")
    plot_polarisation_curve()

    print("  [2/4] Efficiency curves...")
    plot_efficiency()

    print("  [3/4] Dynamic step response...")
    plot_dynamic_response()

    print("  [4/4] 2D operating map...")
    plot_operating_map()

    print("\n  Done. All figures saved to ./plots/\n")