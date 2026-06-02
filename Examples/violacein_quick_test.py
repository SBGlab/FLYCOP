"""
violacein_quick_test.py
=======================
Rapid test of the violacein_flycop2 example:
  - FVA computed only for glucose and 3 flux percentages [0, 50, 100]
    (vs. 8 carbon sources × 11 percentages in the full run)
  - carbon_source fixed to 'glc__D_e'
  - n_trials = 10

Run from the Examples/ directory:
    python violacein_quick_test.py
"""

import copy
import os
import warnings

# Required so that COMETS' JVM can locate the Gurobi Java interface
os.environ.setdefault("GUROBI_COMETS_HOME", r"C:\gurobi901\win64")

import cobra
import cometspy as c
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")

from FLYCOP.OptunaOptimizer import OptunaOptimizer

# ---------------------------------------------------------------------------
# Reuse helper and model-building code from violacein_flycop2
# ---------------------------------------------------------------------------

CARBON_SOURCES = ["glc__D_e"]          # quick test: glucose only
PERCENTAGES    = [0, 50, 100]          # quick test: 3 flux levels

_CARBON_ATOMS = {
    "glc__D_e": 6, "gal_e": 6, "fru_e": 6, "sucr_e": 12,
    "glyc_e": 3, "xyl__D_e": 5, "mal__D_e": 4, "lcts_e": 12,
}

def carbon_normalization(carbon_source: str) -> int:
    return _CARBON_ATOMS[carbon_source]


def make_df_and_graph(strains, metabolites, comets, max_cycles):
    df_media   = copy.deepcopy(comets.media.loc[comets.media["cycle"] < max_cycles])
    df_biomass = copy.deepcopy(
        comets.total_biomass.loc[comets.total_biomass["cycle"] < max_cycles]
    )
    columns = ["cycle"] + strains
    df_biomass.columns = columns

    for d in metabolites:
        columns.append(d)
        met  = df_media.loc[df_media["metabolite"] == d]
        temp = np.zeros(max_cycles)
        df_biomass[d] = temp
        for j in range(1, max_cycles):
            if (met["cycle"] == j).any():
                df_biomass.loc[j - 1, d] = float(met[met["cycle"] == j]["conc_mmol"])

    df_biomass.columns = columns

    # --- save data ---
    np.savetxt(
        "biomass_vs_violacein_e_template.txt",
        df_biomass.values,
        fmt="%s", delimiter="\t",
        header="\t".join(columns),
    )

    # --- save figure ---
    plt.ioff()
    fig, ax = plt.subplots()
    ax.set_xlabel("time (h)")
    ax.set_ylabel("biomass (g/L)")
    colours = ["k", "m"]
    for i, s in enumerate(strains):
        ax.plot(df_biomass["cycle"] * 0.1, df_biomass[s],
                label=s, color=colours[i % len(colours)])
    ax2 = ax.twinx()
    ax2.set_ylabel("violacein_e (mM)")
    ax2.plot(df_biomass["cycle"] * 0.1, df_biomass["violacein_e"],
             label="violacein_e", color="purple", linestyle="--")
    handles, labels = [], []
    for src in (ax, ax2):
        for h, l in zip(*src.get_legend_handles_labels()):
            if l not in labels:
                handles.append(h); labels.append(l)
    plt.legend(handles, labels)
    plt.savefig("biomass_vs_violacein_e_template_plot.pdf")
    plt.close(fig)
    return df_biomass


# ---------------------------------------------------------------------------
# Build models
# ---------------------------------------------------------------------------

def build_violacein_producer():
    model = cobra.io.load_model("iWFL_1372")

    def _add(rxn_id, name, bounds, mets):
        rxn = cobra.Reaction(rxn_id); rxn.name = name; rxn.bounds = bounds
        rxn.add_metabolites(mets); model.add_reactions([rxn])

    i3i3ppa = cobra.Metabolite("2i3i3ppa_c", formula="C11H9N2O2",
                               name="2-imine-3-(indol-3-yl)propanoate",
                               compartment="c", charge=-1)
    _add("vioA", "Flavin-dependent L-tryptophan oxidase", (-1000, 1000), {
        model.metabolites.trp__L_c: -1, model.metabolites.o2_c: -1,
        model.metabolites.h2o2_c: 1, i3i3ppa: 1})

    i3pyridm = cobra.Metabolite("i3pyridm_c", formula="C22H18N4O4",
                                name="indole-3-pyruvate imine dimer",
                                compartment="c", charge=0)
    _add("vioB", "2-imino-3-(indol-3yl)propanoate dimerase", (-1000, 1000), {
        model.metabolites.get_by_id("2i3i3ppa_c"): -2, i3pyridm: 1})

    ptdeovio = cobra.Metabolite("ptdeovio_c", formula="C21H15N3O2",
                                name="protodeoxyviolaceinic acid",
                                compartment="c", charge=0)
    _add("vioE", "Prodeoxyviolacein synthase", (0, 1000), {
        model.metabolites.i3pyridm_c: -1, ptdeovio: 1})

    ptvio = cobra.Metabolite("ptvio_c", formula="C21H15N3O3",
                             name="protoviolaceinic acid",
                             compartment="c", charge=0)
    _add("vioD", "Protodeoxyviolaceinate monooxygenase", (0, 1000), {
        model.metabolites.ptdeovio_c: -1, model.metabolites.o2_c: -1,
        model.metabolites.nadph_c: -1, model.metabolites.h_c: -1,
        model.metabolites.nadp_c: 1, model.metabolites.h2o_c: 1, ptvio: 1})

    violacein_c = cobra.Metabolite("violacein_c", formula="C20H13N3O3",
                                   name="Violacein", compartment="c", charge=0)
    _add("vioC", "Violacein synthase", (0, 1000), {
        model.metabolites.ptvio_c: -1, model.metabolites.o2_c: -1,
        model.metabolites.nadph_c: -1, model.metabolites.h_c: -1,
        model.metabolites.nadp_c: 1, model.metabolites.h2o_c: 1, violacein_c: 1})

    violacein_e = cobra.Metabolite("violacein_e", formula="C20H13N3O3",
                                   name="violacein", compartment="e", charge=0)
    _add("VIOtr", "Violacein diffusion", (0, 1000), {
        model.metabolites.violacein_c: -1, violacein_e: 1})
    model.add_boundary(model.metabolites.get_by_id("violacein_e"), type="exchange")

    model.reactions.SERD_L.bounds  = (0, 0)
    model.reactions.CYSDS.bounds   = (0, 0)
    model.reactions.TRPAS2.bounds  = (0, 0)
    model.reactions.EX_trp__L_e.bounds = (-5, 5)

    # Open all carbon sources (the medium controls what is actually available)
    _ALL_CS = ["glc__D_e","gal_e","fru_e","sucr_e","glyc_e","xyl__D_e","mal__D_e","lcts_e"]
    for met in _ALL_CS:
        model.reactions.get_by_id(f"EX_{met}").lower_bound = -1.72
    return model


def build_trp_transformer():
    model = cobra.io.load_model("iEC1372_W3110")
    model.reactions.EX_glc__D_e.bounds  = (0, 1000)
    model.reactions.EX_trp__L_e.bounds  = (0, 1000)
    _ALL_CS = ["glc__D_e","gal_e","fru_e","sucr_e","glyc_e","xyl__D_e","mal__D_e","lcts_e"]
    for met in _ALL_CS:
        model.reactions.get_by_id(f"EX_{met}").lower_bound = -1.72
    return model


# ---------------------------------------------------------------------------
# Step 1 — load models
# ---------------------------------------------------------------------------
print("\n[1/5] Loading metabolic models …")
_model_producer    = build_violacein_producer()
_model_transformer = build_trp_transformer()
print(f"  Producer    : {len(_model_producer.reactions)} reactions")
print(f"  Transformer : {len(_model_transformer.reactions)} reactions")

# clean copies (before FVA modifies bounds)
_model_producer_clean    = _model_producer.copy()
_model_transformer_clean = _model_transformer.copy()

# ---------------------------------------------------------------------------
# Step 2 — FVA (glucose only, 3 percentages)
# ---------------------------------------------------------------------------
print("\n[2/5] Pre-computing FVA flux bounds (glucose, 3 levels) …")

trp_flux_transformer_dict:    dict = {}
violacein_flux_producer_dict: dict = {}

_ALL_CS = ["glc__D_e","gal_e","fru_e","sucr_e","glyc_e","xyl__D_e","mal__D_e","lcts_e"]

for cs in CARBON_SOURCES:
    for pct in PERCENTAGES:
        for met in _ALL_CS:
            lb = -1.72 if met == cs else 0.0
            _model_transformer.reactions.get_by_id(f"EX_{met}").lower_bound = lb
            _model_producer.reactions.get_by_id(f"EX_{met}").lower_bound    = lb

        fva_t = cobra.flux_analysis.flux_variability_analysis(
            _model_transformer, reaction_list=["TRPAS2"],
            fraction_of_optimum=pct / 100)
        fva_p = cobra.flux_analysis.flux_variability_analysis(
            _model_producer, reaction_list=["vioC"],
            fraction_of_optimum=pct / 100)

        trp_flux = fva_t.loc["TRPAS2", "minimum"]
        vio_flux = fva_p.loc["vioC",   "maximum"]
        trp_flux_transformer_dict[(cs, pct)]    = trp_flux
        violacein_flux_producer_dict[(cs, pct)] = vio_flux
        print(f"  {cs}  {pct:3d}%  TRP={trp_flux:.4f}  VIO={vio_flux:.4f}")

# ---------------------------------------------------------------------------
# Step 3 — COMETS simulation parameters
# ---------------------------------------------------------------------------
print("\n[3/5] Configuring COMETS parameters …")

_comets_params = c.params()
_comets_params.all_params["maxCycles"]           = 240
_comets_params.all_params["timeStep"]            = 0.1
_comets_params.all_params["spaceWidth"]          = 0.05
_comets_params.all_params["allowCellOverlap"]    = True
_comets_params.all_params["deathRate"]           = 0.0
_comets_params.all_params["maxSpaceBiomass"]     = 1000
_comets_params.all_params["defaultVmax"]         = 20
_comets_params.all_params["showCycleTime"]       = False
_comets_params.all_params["useLogNameTimeStamp"] = False
_comets_params.all_params["FluxLogRate"]         = 1
_comets_params.all_params["MediaLogRate"]        = 1
_comets_params.all_params["exchangestyle"]       = "Standard FBA"
_comets_params.all_params["writeTotalBiomassLog"]= True
_comets_params.all_params["writeMediaLog"]       = True


def _build_layout():
    m_t = c.model(copy.deepcopy(_model_transformer_clean))
    m_p = c.model(copy.deepcopy(_model_producer_clean))
    m_t.id = "transformer_strain"
    m_p.id = "producer_strain"
    layout = c.layout([m_t, m_p])
    layout.add_typical_trace_metabolites()
    return layout, m_t, m_p


# ---------------------------------------------------------------------------
# Step 4 — Parameter space (quick test: glucose only)
# ---------------------------------------------------------------------------
PARAM_SPACE = [
    {"type": "int",         "name": "ratio",               "min": -10,  "max": 10},
    {"type": "categorical", "name": "carbon_source",       "values": CARBON_SOURCES},
    {"type": "categorical", "name": "trp_ratio_flux",      "values": PERCENTAGES},
    {"type": "categorical", "name": "violacein_ratio_flux","values": PERCENTAGES},
]

_INITIAL_BIOMASS = 0.05


# ---------------------------------------------------------------------------
# Step 5 — Objective function
# ---------------------------------------------------------------------------
trial_log = []   # collects per-trial results for final display

def evaluate_violacein(params: dict) -> float:
    cs    = params["carbon_source"]
    t_pct = params["trp_ratio_flux"]
    v_pct = params["violacein_ratio_flux"]
    ratio = params["ratio"]

    layout, m_t, m_p = _build_layout()

    trp_ub = trp_flux_transformer_dict[(cs, t_pct)]
    vio_ub = violacein_flux_producer_dict[(cs, v_pct)]
    # Force transformer to degrade tryptophan (TRPAS2 backward, upper_bound = FVA min)
    m_t.change_bounds("TRPAS2", -1000, trp_ub)
    # Force producer to produce EXACTLY vio_ub mmol/gDW/h of violacein
    # (setting lb=ub forces the FBA to always route through this pathway;
    #  without this, a biomass-maximising FBA will never use vioC)
    m_p.change_bounds("vioC", vio_ub, vio_ub)
    layout.update_models()

    if ratio >= 0:
        bm_t = _INITIAL_BIOMASS * max(ratio, 0.01)
        bm_p = _INITIAL_BIOMASS
    else:
        bm_t = _INITIAL_BIOMASS
        bm_p = _INITIAL_BIOMASS * max(-ratio, 0.01)

    layout.initial_pop = [[0.0, 0.0, bm_t, bm_p]]
    conc = 27.7 * (6 / carbon_normalization(cs))
    layout.set_specific_metabolite(cs, conc)

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sim = c.comets(layout, _comets_params)
            sim.run()
    except Exception as exc:
        print(f"  [warn] COMETS failed: {exc}")
        trial_log.append({**params, "violacein_e": 0.0, "fitness": 1e6})
        return 1e6

    try:
        results = make_df_and_graph(
            [m_t.id, m_p.id], ["violacein_e"], sim, 170
        )
        # First trial: show which metabolites were tracked
        if len(trial_log) == 0:
            tracked = sim.media["metabolite"].unique().tolist()
            vio_tracked = "violacein_e" in tracked
            print(f"  [diag] violacein_e tracked in media: {vio_tracked}  "
                  f"(total metabolites: {len(tracked)})")
        violacein = results.loc[165, "violacein_e"]
        bm1       = results.iloc[165, 1]
        bm2       = results.iloc[165, 2]
    except Exception as exc:
        print(f"  [warn] result extraction failed: {exc}")
        trial_log.append({**params, "violacein_e": 0.0, "fitness": 1e6})
        return 1e6

    fitness = 1e6 if violacein <= 0 else 1.0 / violacein
    trial_log.append({**params, "violacein_e": round(violacein, 6),
                      "bm_transformer": round(bm1, 5), "bm_producer": round(bm2, 5),
                      "fitness": round(fitness, 6)})
    print(f"  ratio={ratio:3d}  trp={t_pct:3d}%  vio={v_pct:3d}%  "
          f"ET={bm1:.4f}  EV={bm2:.4f}  violacein_e={violacein:.6f} mM")
    return fitness


# ---------------------------------------------------------------------------
# Run optimisation — 10 trials
# ---------------------------------------------------------------------------
print("\n[4/5] Running Optuna optimisation  (n_trials=10) …\n")

optimizer = OptunaOptimizer(
    configuration_space=PARAM_SPACE,
    objective_function=evaluate_violacein,
    parameters={"n_trials": 10, "direction": "minimize"},
)
best_params = optimizer.run_optimization()

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
print("\n[5/5] Results\n" + "=" * 60)
best_val = optimizer.study.best_value

print(f"  Best 1/violacein : {best_val:.6f}")
if best_val < 1e6:
    print(f"  Max violacein_e  : {1.0 / best_val:.6f} mM  (cycle 165)")
else:
    print("  No violacein produced in any trial.")

print("\n  Best parameter values:")
for k, v in best_params.items():
    print(f"    {k:<24s} = {v}")

print("\n  All trials:")
print(f"  {'ratio':>6}  {'trp%':>5}  {'vio%':>5}  "
      f"{'ET':>8}  {'EV':>8}  {'violacein_e':>12}  {'fitness':>12}")
print("  " + "-" * 70)
for t in trial_log:
    print(f"  {t.get('ratio','-'):>6}  {t.get('trp_ratio_flux','-'):>5}  "
          f"{t.get('violacein_ratio_flux','-'):>5}  "
          f"{t.get('bm_transformer', 0.0):>8.4f}  "
          f"{t.get('bm_producer', 0.0):>8.4f}  "
          f"{t.get('violacein_e', 0.0):>12.6f}  "
          f"{t.get('fitness', 1e6):>12.6f}")

# Optional: Optuna visualisation
try:
    import optuna.visualization as vis
    vis.plot_optimization_history(optimizer.study).write_html(
        "violacein_quick_optim_history.html")
    print("\n  Plot saved: violacein_quick_optim_history.html")
except Exception:
    pass
