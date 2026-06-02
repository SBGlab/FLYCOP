"""
Violacein biosynthesis optimization with FLYCOP (v2 / Optuna backend)
======================================================================
Two-strain *E. coli* co-culture on a defined minimal medium:

  • Transformer strain  (iEC1372_W3110) — converts carbon to L-tryptophan
  • Producer   strain  (iWFL_1372 + vio pathway) — converts trp to violacein

Objective : maximize extracellular violacein at cycle 165.
            Optuna *minimizes*  f(params) = 1 / violacein_e   (or 1e6 if 0).

Parameters optimised
--------------------
  ratio                int  [-10, 10]      relative initial biomass (+ → more transformer)
  carbon_source        categorical         one of 8 carbon sources
  trp_ratio_flux       categorical [%]     TRPAS2 upper-bound fraction of FVA maximum
  violacein_ratio_flux categorical [%]     vioC   upper-bound fraction of FVA maximum

Differences from the original violacein.py
-------------------------------------------
  • Imports: ConfigSpace / smac removed; OptunaOptimizer used instead.
  • Parameter space: canonical list-of-dicts (compatible with all FLYCOP backends).
  • violacein_scenario → evaluate_violacein(params: dict)
    – No longer receives a ``self`` FLYCOP instance; accesses the shared
      cometspy layout and FVA dicts through a closure.
    – Builds a fresh deep-copy of the layout per trial to avoid state leaks.
  • FLYCOP / SMAC3 orchestration replaced by a direct OptunaOptimizer call.
  • Clean model copies saved before FVA so per-trial layouts start from the
    correct exchange-bound state.
"""

import copy
import warnings

import cobra
import cometspy as c
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")   # non-interactive backend for scripts

from FLYCOP.OptunaOptimizer import OptunaOptimizer

###############################################################################
#  Utility — visualisation / data export
###############################################################################

def make_df_and_graph(strains, metabolites, comets, max_cycles):
    """Create a biomass-vs-metabolite DataFrame and save a PDF plot.

    Parameters
    ----------
    strains    : list[str]  model IDs
    metabolites: list[str]  BiGG metabolite IDs to track
    comets     : cometspy.comets  finished simulation object
    max_cycles : int
    Returns
    -------
    pd.DataFrame  columns = [cycle, strain1…, met1…]
    """
    file_name = "_".join(metabolites)
    df_media   = copy.deepcopy(comets.media.loc[comets.media["cycle"] < max_cycles])
    df_biomass = copy.deepcopy(comets.total_biomass.loc[comets.total_biomass["cycle"] < max_cycles])

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
    np.savetxt(
        f"biomass_vs_{file_name}_template.txt",
        df_biomass.values,
        fmt="%s",
        delimiter="\t",
        header="\t".join(columns),
    )

    # --- figure ---
    plt.ioff()
    fig, ax = plt.subplots()
    ax.set_xlabel("time (h)")
    ax.set_ylabel("biomass (g/L)")
    colours = ["k", "m", "b", "g", "r"]
    for i, strain in enumerate(strains):
        ax.plot(df_biomass["cycle"] * 0.1, df_biomass[strain], label=strain,
                color=colours[i % len(colours)])
    ax2 = ax.twinx()
    ax2.set_ylabel("metabolite conc (mM)")
    for m in metabolites:
        ax2.plot(df_biomass["cycle"] * 0.1, df_biomass[m], label=m)
    # unified legend
    handles, labels = [], []
    for src in (ax, ax2):
        for h, l in zip(*src.get_legend_handles_labels()):
            if l not in labels:
                handles.append(h)
                labels.append(l)
    plt.legend(handles, labels)
    plt.savefig(f"biomass_vs_{file_name}_template_plot.pdf")
    plt.close(fig)
    return df_biomass


###############################################################################
#  Model construction — Producer (iWFL_1372 + vio pathway)
###############################################################################

def build_violacein_producer():
    """Load iWFL_1372 and add the five violacein biosynthesis reactions."""
    model = cobra.io.load_model("iWFL_1372")

    # ── vioA : Flavin-dependent L-tryptophan oxidase ──────────────────────
    rxn = cobra.Reaction("vioA")
    rxn.name = "Flavin-dependent L-tryptophan oxidase"
    rxn.bounds = (-1000, 1000)
    i3i3ppa = cobra.Metabolite("2i3i3ppa_c", formula="C11H9N2O2",
                               name="2-imine-3-(indol-3-yl)propanoate",
                               compartment="c", charge=-1)
    rxn.add_metabolites({
        model.metabolites.trp__L_c: -1,
        model.metabolites.o2_c:     -1,
        model.metabolites.h2o2_c:    1,
        i3i3ppa:                     1,
    })
    model.add_reactions([rxn])

    # ── vioB : 2-imino-3-(indol-3yl)propanoate dimerase ───────────────────
    rxn = cobra.Reaction("vioB")
    rxn.name = "2-imino-3-(indol-3yl)propanoate dimerase"
    rxn.bounds = (-1000, 1000)
    i3pyridm = cobra.Metabolite("i3pyridm_c", formula="C22H18N4O4",
                                name="indole-3-pyruvate imine dimer",
                                compartment="c", charge=0)
    rxn.add_metabolites({
        model.metabolites.get_by_id("2i3i3ppa_c"): -2,
        i3pyridm:                                    1,
    })
    model.add_reactions([rxn])

    # ── vioE : Prodeoxyviolacein synthase ─────────────────────────────────
    rxn = cobra.Reaction("vioE")
    rxn.name = "Prodeoxyviolacein synthase"
    rxn.bounds = (0, 1000)
    ptdeovio = cobra.Metabolite("ptdeovio_c", formula="C21H15N3O2",
                                name="protodeoxyviolaceinic acid",
                                compartment="c", charge=0)
    rxn.add_metabolites({
        model.metabolites.i3pyridm_c: -1,
        ptdeovio:                      1,
    })
    model.add_reactions([rxn])

    # ── vioD : Protodeoxyviolaceinate monooxygenase ───────────────────────
    rxn = cobra.Reaction("vioD")
    rxn.name = "Protodeoxyviolaceinate monooxygenase"
    rxn.bounds = (0, 1000)
    ptvio = cobra.Metabolite("ptvio_c", formula="C21H15N3O3",
                             name="protoviolaceinic acid",
                             compartment="c", charge=0)
    rxn.add_metabolites({
        model.metabolites.ptdeovio_c: -1,
        model.metabolites.o2_c:       -1,
        model.metabolites.nadph_c:    -1,
        model.metabolites.h_c:        -1,
        model.metabolites.nadp_c:      1,
        model.metabolites.h2o_c:       1,
        ptvio:                         1,
    })
    model.add_reactions([rxn])

    # ── vioC : Violacein synthase ─────────────────────────────────────────
    rxn = cobra.Reaction("vioC")
    rxn.name = "Violacein synthase"
    rxn.bounds = (0, 1000)
    violacein_c = cobra.Metabolite("violacein_c", formula="C20H13N3O3",
                                   name="Violacein", compartment="c", charge=0)
    rxn.add_metabolites({
        model.metabolites.ptvio_c:  -1,
        model.metabolites.o2_c:     -1,
        model.metabolites.nadph_c:  -1,
        model.metabolites.h_c:      -1,
        model.metabolites.nadp_c:    1,
        model.metabolites.h2o_c:     1,
        violacein_c:                 1,
    })
    model.add_reactions([rxn])

    # ── VIOtr : Violacein diffusion to extracellular space ────────────────
    rxn = cobra.Reaction("VIOtr")
    rxn.name = "Violacein diffusion"
    rxn.bounds = (0, 1000)
    violacein_e = cobra.Metabolite("violacein_e", formula="C20H13N3O3",
                                   name="violacein", compartment="e", charge=0)
    rxn.add_metabolites({
        model.metabolites.violacein_c: -1,
        violacein_e:                    1,
    })
    model.add_reactions([rxn])
    model.add_boundary(model.metabolites.get_by_id("violacein_e"), type="exchange")

    # ── Knock-outs and exchange adjustments ───────────────────────────────
    model.reactions.SERD_L.bounds  = (0, 0)
    model.reactions.CYSDS.bounds   = (0, 0)
    model.reactions.TRPAS2.bounds  = (0, 0)
    model.reactions.EX_trp__L_e.bounds = (-5, 5)

    for met in CARBON_SOURCES:
        model.reactions.get_by_id(f"EX_{met}").lower_bound = -1.72

    cobra.io.write_sbml_model(model, "violacein_producer.xml")
    return model


###############################################################################
#  Model construction — Transformer (iEC1372_W3110)
###############################################################################

def build_trp_transformer():
    """Load iEC1372_W3110 and configure tryptophan export."""
    model = cobra.io.load_model("iEC1372_W3110")
    model.reactions.EX_glc__D_e.bounds  = (0,  1000)
    model.reactions.EX_trp__L_e.bounds  = (0,  1000)
    for met in CARBON_SOURCES:
        model.reactions.get_by_id(f"EX_{met}").lower_bound = -1.72
    cobra.io.write_sbml_model(model, "trp_transformer.xml")
    return model


def test_violacein_model(model_producer, violacein_flux, carbon_source):
    """Sanity-check: fix vioC lower bound and call optimize()."""
    for met in CARBON_SOURCES:
        lb = -1.72 if met == carbon_source else 0
        model_producer.reactions.get_by_id(f"EX_{met}").lower_bound = lb
    model_producer.reactions.vioC.lower_bound = violacein_flux
    sol = model_producer.optimize()
    print(f"[sanity] carbon={carbon_source}  violacein_flux={violacein_flux:.4f}")
    print(f"         growth={sol.objective_value:.4f}  vioC={sol.fluxes['vioC']:.4f}")


###############################################################################
#  Constants
###############################################################################

CARBON_SOURCES = [
    "glc__D_e", "gal_e", "fru_e", "sucr_e",
    "glyc_e", "xyl__D_e", "mal__D_e", "lcts_e",
]
PERCENTAGES = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

_CARBON_ATOMS = {
    "glc__D_e": 6, "gal_e": 6, "fru_e": 6, "sucr_e": 12,
    "glyc_e": 3, "xyl__D_e": 5, "mal__D_e": 4, "lcts_e": 12,
}


def carbon_normalization(carbon_source: str) -> int:
    """Return number of carbon atoms in one molecule of *carbon_source*."""
    return _CARBON_ATOMS[carbon_source]


###############################################################################
#  Step 1 — Build cobra models
###############################################################################

print("Building violacein producer model (iWFL_1372 + vio pathway) …")
_model_producer = build_violacein_producer()

print("Building tryptophan transformer model (iEC1372_W3110) …")
_model_transformer = build_trp_transformer()

###############################################################################
#  Step 2 — Save clean copies BEFORE FVA modifies exchange bounds
###############################################################################

_model_producer_clean    = _model_producer.copy()
_model_transformer_clean = _model_transformer.copy()

###############################################################################
#  Step 3 — Pre-compute FVA flux bounds for every (carbon_source, %) pair
###############################################################################

print("Pre-computing FVA flux bounds … (this may take a few minutes)")

trp_flux_transformer_dict: dict   = {}
violacein_flux_producer_dict: dict = {}

for carbon_source in CARBON_SOURCES:
    for pct in PERCENTAGES:
        # Set all carbon sources: open only the current one
        for met in CARBON_SOURCES:
            lb = -1.72 if met == carbon_source else 0.0
            _model_transformer.reactions.get_by_id(f"EX_{met}").lower_bound = lb
            _model_producer.reactions.get_by_id(f"EX_{met}").lower_bound    = lb

        fva_t = cobra.flux_analysis.flux_variability_analysis(
            _model_transformer,
            reaction_list=["TRPAS2"],
            fraction_of_optimum=pct / 100,
        )
        fva_p = cobra.flux_analysis.flux_variability_analysis(
            _model_producer,
            reaction_list=["vioC"],
            fraction_of_optimum=pct / 100,
        )

        trp_flux  = fva_t.loc["TRPAS2", "minimum"]
        vio_flux  = fva_p.loc["vioC",   "maximum"]
        trp_flux_transformer_dict[(carbon_source, pct)]    = trp_flux
        violacein_flux_producer_dict[(carbon_source, pct)] = vio_flux

        print(f"  {carbon_source:12s}  {pct:3d}%  "
              f"TRP={trp_flux:.4f}  VIO={vio_flux:.4f}")

# Sanity check with glucose at 80 %
test_violacein_model(
    _model_producer,
    violacein_flux_producer_dict[("glc__D_e", 80)],
    "glc__D_e",
)

###############################################################################
#  Step 4 — COMETS simulation parameters (shared, read-only per trial)
###############################################################################

_comets_params = c.params()
_comets_params.all_params["maxCycles"]          = 240
_comets_params.all_params["timeStep"]           = 0.1
_comets_params.all_params["spaceWidth"]         = 0.05
_comets_params.all_params["allowCellOverlap"]   = True
_comets_params.all_params["deathRate"]          = 0.0
_comets_params.all_params["maxSpaceBiomass"]    = 1000
_comets_params.all_params["defaultVmax"]        = 20
_comets_params.all_params["showCycleTime"]      = True
_comets_params.all_params["useLogNameTimeStamp"]= False
_comets_params.all_params["FluxLogRate"]        = 1
_comets_params.all_params["MediaLogRate"]       = 1
_comets_params.all_params["exchangestyle"]      = "Standard FBA"
_comets_params.all_params["writeTotalBiomassLog"] = True
_comets_params.all_params["writeMediaLog"]      = True


def _build_layout():
    """Return a fresh (deep-copied) cometspy layout ready for one trial.

    Using the *clean* cobra model copies (saved before FVA) guarantees that
    each trial starts from correct, uncontaminated exchange-bound state.
    """
    m_t = c.model(copy.deepcopy(_model_transformer_clean))
    m_p = c.model(copy.deepcopy(_model_producer_clean))
    m_t.id = "transformer_strain"
    m_p.id = "producer_strain"
    layout = c.layout([m_t, m_p])
    layout.add_typical_trace_metabolites()
    return layout, m_t, m_p


###############################################################################
#  Step 5 — Parameter search space  (canonical FLYCOP list-of-dicts)
###############################################################################

PARAM_SPACE = [
    # Relative initial biomass: positive → more transformer; negative → more producer
    {
        "type": "int",
        "name": "ratio",
        "min": -10,
        "max":  10,
    },
    # Carbon source available in the medium
    {
        "type":   "categorical",
        "name":   "carbon_source",
        "values": CARBON_SOURCES,
    },
    # TRPAS2 flux upper-bound expressed as % of FVA maximum in the transformer
    {
        "type":   "categorical",
        "name":   "trp_ratio_flux",
        "values": PERCENTAGES,
    },
    # vioC flux upper-bound expressed as % of FVA maximum in the producer
    {
        "type":   "categorical",
        "name":   "violacein_ratio_flux",
        "values": PERCENTAGES,
    },
]


###############################################################################
#  Step 6 — Objective function
###############################################################################

_INITIAL_BIOMASS = 0.05   # g DW / gridbox baseline


def evaluate_violacein(params: dict) -> float:
    """Run one COMETS trial and return 1/violacein_e (Optuna minimizes).

    Parameters
    ----------
    params : dict
        Keys: ratio, carbon_source, trp_ratio_flux, violacein_ratio_flux
    Returns
    -------
    float
        1 / violacein_e at cycle 165, or 1e6 if production is zero / simulation fails.
    """
    carbon_source        = params["carbon_source"]
    trp_pct              = params["trp_ratio_flux"]
    vio_pct              = params["violacein_ratio_flux"]
    ratio                = params["ratio"]

    # ── Fresh layout for this trial ────────────────────────────────────────
    layout, m_t, m_p = _build_layout()

    # ── Apply FVA-derived flux bounds ──────────────────────────────────────
    trp_ub = trp_flux_transformer_dict[(carbon_source, trp_pct)]
    vio_ub = violacein_flux_producer_dict[(carbon_source, vio_pct)]
    m_t.change_bounds("TRPAS2", -1000, trp_ub)
    # Force lb=ub so FBA must route exactly vio_ub flux through violacein synthesis.
    # With lb=0 the biomass-maximising FBA never activates this pathway.
    m_p.change_bounds("vioC", vio_ub, vio_ub)
    layout.update_models()

    # ── Initial biomass (ratio controls transformer : producer balance) ────
    if ratio >= 0:
        bm_transformer = _INITIAL_BIOMASS * max(ratio, 0.01)
        bm_producer    = _INITIAL_BIOMASS
    else:
        bm_transformer = _INITIAL_BIOMASS
        bm_producer    = _INITIAL_BIOMASS * max(-ratio, 0.01)

    layout.initial_pop = [[0.0, 0.0, bm_transformer, bm_producer]]

    # ── Medium: one carbon source, normalised to glucose carbon content ────
    conc = 27.7 * (6 / carbon_normalization(carbon_source))
    layout.set_specific_metabolite(carbon_source, conc)

    # ── Run COMETS ─────────────────────────────────────────────────────────
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sim = c.comets(layout, _comets_params)
            sim.run()
    except Exception as exc:
        print(f"  [warning] COMETS failed: {exc}")
        return 1e6

    # ── Extract violacein at cycle 165 ─────────────────────────────────────
    try:
        results = make_df_and_graph(
            [m_t.id, m_p.id], ["violacein_e"], sim, 170
        )
        violacein = results.loc[165, "violacein_e"]
        bm1       = results.iloc[165, 1]
        bm2       = results.iloc[165, 2]
    except Exception as exc:
        print(f"  [warning] result extraction failed: {exc}")
        return 1e6

    config_str = ", ".join(f"{k}={v}" for k, v in params.items())
    print(f"  {config_str}  "
          f"ET={bm1:.4f}  EV={bm2:.4f}  violacein={violacein:.6f}")

    if violacein <= 0:
        return 1e6

    return 1.0 / violacein


###############################################################################
#  Step 7 — Run Optuna optimisation
###############################################################################

optimizer = OptunaOptimizer(
    configuration_space=PARAM_SPACE,
    objective_function=evaluate_violacein,
    parameters={
        "n_trials":  1500,        # same trial budget as the original SMAC3 run
        "direction": "minimize",  # minimize 1/violacein = maximize violacein
    },
)

print()
print("=" * 70)
print("  FLYCOP — Violacein optimisation")
print("  Transformer: iEC1372_W3110   Producer: iWFL_1372 + vio pathway")
print("  Simulator : COMETS (240 cycles)   Optimizer: Optuna TPE (1500 trials)")
print("=" * 70)
print()

best_params = optimizer.run_optimization()

###############################################################################
#  Step 8 — Report results
###############################################################################

best_val = optimizer.study.best_value

print()
print("=" * 70)
print("  Optimisation complete")
print("=" * 70)
print(f"  Best 1/violacein  : {best_val:.6f}")
print(f"  Max violacein_e   : {1.0 / best_val:.6f} mM  (cycle 165)")
print()
print("  Best parameter values:")
for name, value in best_params.items():
    print(f"    {name:<24s} = {value}")

# --- Optional: save Optuna visualisations --------------------------------
try:
    import optuna.visualization as vis

    vis.plot_optimization_history(optimizer.study).write_html(
        "violacein_optim_history.html"
    )
    if len(PARAM_SPACE) > 1:
        vis.plot_param_importances(optimizer.study).write_html(
            "violacein_param_importances.html"
        )
    print("\n  Plots saved: violacein_optim_history.html, "
          "violacein_param_importances.html")
except Exception:
    pass
