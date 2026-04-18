# Dynamic Modelling of an Alkaline Water Electrolysis (AWE) Cell

A physics-based electrochemical simulation of an AWE cell capturing static
polarisation behaviour, voltage loss mechanisms, and transient dynamic response
under variable load conditions.

Built to demonstrate core competencies for PhD research in green hydrogen
systems: dynamic physics-based modelling · numerical ODE integration ·
parametric analysis · scientific visualisation.

---

## What this project does

When current is pushed through an AWE cell to split water into hydrogen and
oxygen, the cell voltage is not constant — it depends on the current density,
the temperature, and how fast the load changes. This project models all of that.

**Static model** — decomposes cell voltage into three physical terms:
V(I, T)  =  E_rev(T)  +  η_act(I, T)  +  η_ohm(I, T)

| Term | Physical meaning |
|------|-----------------|
| `E_rev(T)` | Reversible (Nernst) voltage — thermodynamic minimum |
| `η_act(I, T)` | Activation overpotential — Butler-Volmer kinetics |
| `η_ohm(I, T)` | Ohmic overpotential — membrane + electrolyte resistance |

**Dynamic model** — first-order ODE capturing transient response to load steps:
τ · dV/dt  =  V_ss(I, T)  −  V(t)

Solved numerically using `scipy` RK45 integrator.

---

# Check the figures in pdf file as i have shared

---

## How to run

**Install dependencies:**
```bash
pip install numpy matplotlib scipy
```

**Run the simulation:**
```bash
python awe_simulation.py
```

All 4 figures are saved automatically to the `plots/` folder.

**Requirements:** Python 3.9 or newer. No non-standard dependencies.

---

## Project structure
awe-dynamic-model/
├── awe_simulation.py          # full model — physics, ODE solver, all plots
├── README.md
└── plots/
├── fig1_polarisation_curve.png
├── fig2_efficiency.png
├── fig3_dynamic_response.png
└── fig4_operating_map.png

---

## Key parameters

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Gas constant | R | 8.314 | J mol⁻¹ K⁻¹ |
| Faraday constant | F | 96 485 | C mol⁻¹ |
| Thermoneutral voltage | V_HHV | 1.481 | V |
| Charge transfer coefficient | α | 0.40 | — |
| Exchange current density | i₀ (25 °C) | 1×10⁻³ | A m⁻² |
| Area-specific resistance | R_ohm (25 °C) | 0.18 | Ω m² |
| Faradaic efficiency | η_faradaic | 0.98 | — |
| Time constant (typical) | τ | 1–30 | s |

---

## Sample output

Running `python awe_simulation.py` prints this summary before the plots open:
===========================================================
AWE CELL — MODEL SUMMARY
Operating temperature      : 60.0 °C
Operating current density  : 200.0 A m⁻²
Dynamic time constant      : 5.0 s
Reversible voltage E_rev   : 1.1999  V
Activation loss  eta_act   : 0.7253  V   (37.2% of total)
Ohmic loss       eta_ohm   : 0.0254  V   (1.3% of total)
Total cell voltage         : 1.9505  V
Voltage efficiency         : 75.93 %
System efficiency (HHV)    : 74.41 %
Power density              : 0.390  kW m⁻²
Settling time (5τ)         : 25.0   s

---

## References

1. Ulleberg, Ø. (2003). Modeling of advanced alkaline electrolyzers.
   *Int. J. Hydrogen Energy*, 28(1), 21–33.
   https://doi.org/10.1016/S0360-3199(02)00033-2

2. Hammoudi, M. et al. (2012). New multi-physics approach for modelling
   and design of alkaline electrolyzers.
   *Int. J. Hydrogen Energy*, 37(19), 13895–13913.
   https://doi.org/10.1016/j.ijhydene.2012.07.015

3. Ursua, A., Gandia, L. M., & Sanchis, P. (2011). Hydrogen production
   from water electrolysis: current status and future trends.
   *Proceedings of the IEEE*, 100(2), 410–426.
   https://doi.org/10.1109/JPROC.2011.2156750

4. Abdin, Z., Webb, C. J., & Gray, E. M. (2015). Modelling and simulation
   of a PEM electrolyser cell.
   *Int. J. Hydrogen Energy*, 40(39), 13243–13257.
   https://doi.org/10.1016/j.ijhydene.2015.07.129

5. Zeng, K., & Zhang, D. (2010). Recent progress in alkaline water
   electrolysis for hydrogen production.
   *Progress in Energy and Combustion Science*, 36(3), 307–326.
   https://doi.org/10.1016/j.pecs.2009.11.002

---

## Author

**Rohit Dogra**
Department of Mechatronics and Cyber Physical Systems
Mechanical Engineering
Deggendorf Institute of Technology, Germany
