"""
E. coli (iML1515) + P. putida (iJN1462) co-culture on M9 minimal medium
========================================================================
Objective  : maximize total consortium biomass after 100 COMETS cycles
Optimizer  : Optuna (Windows-compatible, Bayesian TPE sampler)
Simulator  : COMETS via COMETSSimulator

Parameters explored
-------------------
  glc_conc             D-glucose initial concentration  [1 – 20 mmol/gridbox]
  nh4_conc             Ammonium  initial concentration  [1 – 50 mmol/gridbox]
  biomass_ecoli_init   E. coli   initial biomass        [0.001 – 0.2 g DW]
  biomass_pputida_init P. putida initial biomass        [0.001 – 0.2 g DW]

How it works
------------
Optuna minimizes the objective f(params) = 1 / total_final_biomass.
Minimizing 1/x is equivalent to maximizing x, so the optimizer steers
toward parameter combinations that produce the largest consortium biomass.
"""

import warnings

from cobra.io import load_model

from FLYCOP.COMETSSimulator import COMETSSimulator
from FLYCOP.OptunaOptimizer import OptunaOptimizer


# ---------------------------------------------------------------------------
# 1. M9 minimal medium  (BiGG exchange-metabolite IDs, mmol/gridbox)
#    Concentrations of glucose and ammonium are tuned by the optimizer.
# ---------------------------------------------------------------------------
M9_BASE: dict = {
    "glc__D_e":   10.0,   # D-glucose        — carbon & energy source
    "nh4_e":      10.0,   # ammonium         — nitrogen source
    "pi_e":        1.0,   # inorganic phosphate
    "so4_e":       2.0,   # sulfate
    "mg2_e":       2.0,   # magnesium
    "k_e":         3.0,   # potassium
    "na1_e":      10.0,   # sodium
    "cl_e":        0.5,   # chloride
    "o2_e":       20.0,   # oxygen (aerobic conditions)
    "h2o_e":     100.0,   # water
    "h_e":         0.0,   # proton (pH-buffered)
}


# ---------------------------------------------------------------------------
# 2. Load metabolic models once — shared across all optimization trials
# ---------------------------------------------------------------------------
print("Loading iML1515  (E. coli K-12 MG1655) …")
ecoli = load_model("iML1515")

print("Loading iJN1462  (P. putida KT2440)    …")
pputida = load_model("iJN1462")

print(f"  iML1515  — {len(ecoli.reactions)} reactions, {len(ecoli.metabolites)} metabolites")
print(f"  iJN1462  — {len(pputida.reactions)} reactions, {len(pputida.metabolites)} metabolites")
print()


# ---------------------------------------------------------------------------
# 3. Parameter search space  (canonical FLYCOP list-of-dicts)
# ---------------------------------------------------------------------------
PARAM_SPACE = [
    {
        "type": "float",
        "name": "glc_conc",
        "min":  1.0,
        "max": 20.0,
    },
    {
        "type": "float",
        "name": "nh4_conc",
        "min":  1.0,
        "max": 50.0,
    },
    {
        "type": "float",
        "name": "biomass_ecoli_init",
        "min": 0.001,
        "max": 0.2,
    },
    {
        "type": "float",
        "name": "biomass_pputida_init",
        "min": 0.001,
        "max": 0.2,
    },
]


# ---------------------------------------------------------------------------
# 4. Objective function
#    Called by Optuna with a plain dict {param_name: value} per trial.
#    Returns 1 / total_final_biomass so Optuna can *minimize* it.
# ---------------------------------------------------------------------------
def evaluate_consortium(params: dict) -> float:
    """
    Instantiate a fresh COMETSSimulator, configure M9 medium with the
    trial's parameter values, run the simulation, and return the
    reciprocal of the total final biomass (Optuna minimizes → biomass is
    maximized).

    A penalty of 1e6 is returned for trials where the simulation fails or
    produces zero/negative biomass.
    """
    sim = COMETSSimulator()

    # --- Build M9 medium for this trial ---
    media = dict(M9_BASE)
    media["glc__D_e"] = params["glc_conc"]
    media["nh4_e"]    = params["nh4_conc"]

    # --- Initial biomass for each organism ---
    initial_biomasses = {
        ecoli.id:   params["biomass_ecoli_init"],
        pputida.id: params["biomass_pputida_init"],
    }

    # --- Load consortium and medium ---
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")   # suppress BiGG namespace warnings
        sim.load_consortia([ecoli, pputida], initial_biomases=initial_biomasses)
    sim.load_media(media)

    # --- Run COMETS (100 cycles by default) ---
    try:
        sim.simulate()
    except Exception as exc:
        print(f"  [warning] simulation failed ({exc})")
        return 1e6

    # --- Extract total biomass at the final cycle ---
    # cometspy.total_biomass DataFrame: columns = ['cycle', <model.id>, ...]
    bm_df = sim.results.total_biomass
    model_cols = [c for c in bm_df.columns if c != "cycle"]
    total_final_biomass = float(bm_df.iloc[-1][model_cols].sum())

    if total_final_biomass <= 0:
        return 1e6

    # Return reciprocal: minimizing 1/B ↔ maximizing B
    return 1.0 / total_final_biomass


# ---------------------------------------------------------------------------
# 5. Configure and run the Optuna optimizer
# ---------------------------------------------------------------------------
optimizer = OptunaOptimizer(
    configuration_space=PARAM_SPACE,
    objective_function=evaluate_consortium,
    parameters={
        "n_trials":  50,          # increase for a more thorough search
        "direction": "minimize",  # minimize 1/biomass = maximize biomass
    },
)

print("=" * 65)
print("  FLYCOP optimization: E. coli + P. putida  |  M9  |  Biomass")
print(f"  Models   : {ecoli.id}  +  {pputida.id}")
print(f"  Simulator: COMETS  (100 cycles per trial)")
print(f"  Optimizer: Optuna TPE sampler  (50 trials)")
print("=" * 65)
print()

best_params = optimizer.run_optimization()


# ---------------------------------------------------------------------------
# 6. Report results
# ---------------------------------------------------------------------------
best_val = optimizer.study.best_value

print()
print("=" * 65)
print("  Optimization complete")
print("=" * 65)
print(f"  Best 1/biomass    : {best_val:.6f}")
print(f"  Max total biomass : {1.0 / best_val:.4f} g DW / gridbox")
print()
print("  Best parameter values:")
for name, value in best_params.items():
    print(f"    {name:<28s} = {value:.5f}")

# --- Optional: visualize Optuna convergence history ---
try:
    import optuna.visualization as vis
    import matplotlib
    matplotlib.use("Agg")   # non-interactive backend for scripts

    fig = vis.plot_optimization_history(optimizer.study)
    fig.write_html("optim_history_ecoli_pputida_m9.html")
    print("\n  Convergence plot saved to: optim_history_ecoli_pputida_m9.html")

    fig2 = vis.plot_param_importances(optimizer.study)
    fig2.write_html("param_importances_ecoli_pputida_m9.html")
    print("  Parameter importance plot saved to: param_importances_ecoli_pputida_m9.html")
except Exception:
    pass   # visualization is optional
