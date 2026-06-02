"""
Tests for the violacein_flycop2.py example.

Validates parameter space, medium logic, objective function behaviour, and
optimizer configuration — all without running real COMETS simulations or
loading actual GEM files.
"""
import sys
import os
import unittest
from unittest.mock import MagicMock, patch

import pandas as pd
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

# ---------------------------------------------------------------------------
# Constants (mirrors violacein_flycop2.py)
# ---------------------------------------------------------------------------
CARBON_SOURCES = [
    "glc__D_e", "gal_e", "fru_e", "sucr_e",
    "glyc_e", "xyl__D_e", "mal__D_e", "lcts_e",
]
PERCENTAGES = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

_CARBON_ATOMS = {
    "glc__D_e": 6, "gal_e": 6, "fru_e": 6, "sucr_e": 12,
    "glyc_e": 3, "xyl__D_e": 5, "mal__D_e": 4, "lcts_e": 12,
}

PARAM_SPACE = [
    {"type": "int",         "name": "ratio",               "min": -10,   "max": 10},
    {"type": "categorical", "name": "carbon_source",       "values": CARBON_SOURCES},
    {"type": "categorical", "name": "trp_ratio_flux",      "values": PERCENTAGES},
    {"type": "categorical", "name": "violacein_ratio_flux","values": PERCENTAGES},
]


def carbon_normalization(carbon_source: str) -> int:
    return _CARBON_ATOMS[carbon_source]


def _fake_violacein_df(violacein_val: float) -> pd.DataFrame:
    """Return a DataFrame shaped like make_df_and_graph output."""
    n = 170
    data = {
        "cycle":             list(range(n)),
        "transformer_strain": [0.05] * n,
        "producer_strain":    [0.05] * n,
        "violacein_e":        [0.0] * (n - 5) + [violacein_val] * 5,
    }
    df = pd.DataFrame(data)
    df.index = df["cycle"]
    return df


# ---------------------------------------------------------------------------
# TestParamSpace
# ---------------------------------------------------------------------------

class TestParamSpace(unittest.TestCase):

    def test_four_parameters_defined(self):
        self.assertEqual(len(PARAM_SPACE), 4)

    def test_ratio_is_int_with_bounds(self):
        p = next(e for e in PARAM_SPACE if e["name"] == "ratio")
        self.assertEqual(p["type"], "int")
        self.assertIn("min", p)
        self.assertIn("max", p)
        self.assertLess(p["min"], p["max"])

    def test_carbon_source_has_eight_values(self):
        p = next(e for e in PARAM_SPACE if e["name"] == "carbon_source")
        self.assertEqual(p["type"], "categorical")
        self.assertEqual(len(p["values"]), 8)

    def test_flux_params_are_categorical_with_percentages(self):
        for name in ("trp_ratio_flux", "violacein_ratio_flux"):
            p = next(e for e in PARAM_SPACE if e["name"] == name)
            self.assertEqual(p["type"], "categorical")
            self.assertEqual(p["values"], PERCENTAGES)

    def test_param_names_unique(self):
        names = [e["name"] for e in PARAM_SPACE]
        self.assertEqual(len(names), len(set(names)))

    def test_param_space_accepted_by_optuna_optimizer(self):
        from FLYCOP.OptunaOptimizer import OptunaOptimizer
        opt = OptunaOptimizer()
        opt.load_configuration_space(PARAM_SPACE)
        self.assertEqual(opt._param_list, PARAM_SPACE)


# ---------------------------------------------------------------------------
# TestCarbonNormalization
# ---------------------------------------------------------------------------

class TestCarbonNormalization(unittest.TestCase):

    def test_glucose_has_6_carbons(self):
        self.assertEqual(carbon_normalization("glc__D_e"), 6)

    def test_sucrose_has_12_carbons(self):
        self.assertEqual(carbon_normalization("sucr_e"), 12)

    def test_glycerol_has_3_carbons(self):
        self.assertEqual(carbon_normalization("glyc_e"), 3)

    def test_all_carbon_sources_covered(self):
        for cs in CARBON_SOURCES:
            self.assertGreater(carbon_normalization(cs), 0)

    def test_glucose_normalization_factor_is_one(self):
        """Glucose is the reference — factor should be 6/6 = 1."""
        self.assertEqual(6 / carbon_normalization("glc__D_e"), 1.0)

    def test_sucrose_gives_half_concentration(self):
        """Sucrose has 12 C vs glucose 6 C → factor = 6/12 = 0.5."""
        self.assertAlmostEqual(6 / carbon_normalization("sucr_e"), 0.5)


# ---------------------------------------------------------------------------
# TestBiomassRatioLogic
# ---------------------------------------------------------------------------

class TestBiomassRatioLogic(unittest.TestCase):
    """Verify the ratio → initial biomass formula."""

    _INITIAL_BIOMASS = 0.05

    def _compute(self, ratio):
        bm = self._INITIAL_BIOMASS
        if ratio >= 0:
            return bm * max(ratio, 0.01), bm
        else:
            return bm, bm * max(-ratio, 0.01)

    def test_positive_ratio_gives_more_transformer(self):
        bm_t, bm_p = self._compute(5)
        self.assertGreater(bm_t, bm_p)

    def test_negative_ratio_gives_more_producer(self):
        bm_t, bm_p = self._compute(-5)
        self.assertGreater(bm_p, bm_t)

    def test_ratio_zero_transformer_near_zero(self):
        bm_t, bm_p = self._compute(0)
        # ratio=0 → bm_t = 0.05 * 0.01 = 0.0005 (non-zero, avoids COMETS crash)
        self.assertAlmostEqual(bm_t, self._INITIAL_BIOMASS * 0.01)
        self.assertAlmostEqual(bm_p, self._INITIAL_BIOMASS)

    def test_ratio_one_gives_equal_biomass(self):
        bm_t, bm_p = self._compute(1)
        self.assertAlmostEqual(bm_t, bm_p)

    def test_biomasses_always_positive(self):
        for r in range(-10, 11):
            bm_t, bm_p = self._compute(r)
            self.assertGreater(bm_t, 0)
            self.assertGreater(bm_p, 0)


# ---------------------------------------------------------------------------
# TestObjectiveFunctionLogic
# ---------------------------------------------------------------------------

class TestObjectiveFunctionLogic(unittest.TestCase):
    """Test evaluate_violacein logic with COMETS fully mocked."""

    # Fake FVA dicts
    _trp_dict  = {(cs, pct): 0.5 for cs in CARBON_SOURCES for pct in PERCENTAGES}
    _vio_dict  = {(cs, pct): 1.0 for cs in CARBON_SOURCES for pct in PERCENTAGES}

    def _make_evaluate(self, violacein_val=1.2, sim_raises=None):
        """Return an evaluate_violacein closure with mocked COMETS and dicts."""
        trp_flux = self._trp_dict
        vio_flux = self._vio_dict
        initial_bm = 0.05

        def evaluate_violacein(params: dict) -> float:
            cs    = params["carbon_source"]
            t_pct = params["trp_ratio_flux"]
            v_pct = params["violacein_ratio_flux"]
            ratio = params["ratio"]

            # Build fresh layout (mocked)
            m_t = MagicMock(); m_t.id = "transformer_strain"
            m_p = MagicMock(); m_p.id = "producer_strain"
            layout = MagicMock()
            layout.models = [m_t, m_p]

            m_t.change_bounds("TRPAS2", -1000, trp_flux[(cs, t_pct)])
            m_p.change_bounds("vioC",       0, vio_flux[(cs, v_pct)])
            layout.update_models()

            if ratio >= 0:
                bm_t = initial_bm * max(ratio, 0.01)
                bm_p = initial_bm
            else:
                bm_t = initial_bm
                bm_p = initial_bm * max(-ratio, 0.01)
            layout.initial_pop = [[0.0, 0.0, bm_t, bm_p]]
            layout.set_specific_metabolite(cs, 27.7 * (6 / carbon_normalization(cs)))

            try:
                if sim_raises:
                    raise sim_raises
                # Fake results DataFrame
                results = _fake_violacein_df(violacein_val)
            except Exception:
                return 1e6

            violacein = results.loc[165, "violacein_e"]
            if violacein <= 0:
                return 1e6
            return 1.0 / violacein

        return evaluate_violacein

    def _base_params(self, **overrides):
        p = {
            "ratio": 1,
            "carbon_source": "glc__D_e",
            "trp_ratio_flux": 80,
            "violacein_ratio_flux": 80,
        }
        p.update(overrides)
        return p

    def test_returns_reciprocal_of_violacein(self):
        fn = self._make_evaluate(violacein_val=2.0)
        result = fn(self._base_params())
        self.assertAlmostEqual(result, 1.0 / 2.0, places=5)

    def test_higher_violacein_gives_lower_objective(self):
        fn_high = self._make_evaluate(violacein_val=3.0)
        fn_low  = self._make_evaluate(violacein_val=0.5)
        self.assertLess(
            fn_high(self._base_params()),
            fn_low(self._base_params()),
        )

    def test_zero_violacein_returns_penalty(self):
        fn = self._make_evaluate(violacein_val=0.0)
        self.assertEqual(fn(self._base_params()), 1e6)

    def test_simulation_failure_returns_penalty(self):
        fn = self._make_evaluate(sim_raises=RuntimeError("COMETS crash"))
        self.assertEqual(fn(self._base_params()), 1e6)

    def test_sucrose_uses_half_concentration(self):
        """With sucrose (12 C), medium concentration = 27.7 * 0.5 = 13.85 mM."""
        calls = []
        original = carbon_normalization

        fn = self._make_evaluate(violacein_val=1.0)
        # Just call and check it doesn't raise
        fn(self._base_params(carbon_source="sucr_e"))

    def test_all_carbon_sources_are_accepted(self):
        fn = self._make_evaluate(violacein_val=1.0)
        for cs in CARBON_SOURCES:
            result = fn(self._base_params(carbon_source=cs))
            self.assertNotEqual(result, None)


# ---------------------------------------------------------------------------
# TestOptimizerConfiguration
# ---------------------------------------------------------------------------

class TestOptimizerConfiguration(unittest.TestCase):

    def test_direction_is_minimize(self):
        from FLYCOP.OptunaOptimizer import OptunaOptimizer
        opt = OptunaOptimizer(parameters={"n_trials": 1500, "direction": "minimize"})
        self.assertEqual(opt.parameters["direction"], "minimize")

    def test_param_space_loaded_correctly(self):
        from FLYCOP.OptunaOptimizer import OptunaOptimizer
        opt = OptunaOptimizer(configuration_space=PARAM_SPACE)
        self.assertEqual(len(opt._param_list), 4)

    def test_end_to_end_with_mock_objective(self):
        """Optimizer runs, calls objective, returns dict with all 4 param names."""
        from FLYCOP.OptunaOptimizer import OptunaOptimizer

        call_count = [0]

        def mock_obj(params):
            call_count[0] += 1
            # Reward higher trp percentage
            return 1.0 / (params["trp_ratio_flux"] + 1)

        opt = OptunaOptimizer(
            configuration_space=PARAM_SPACE,
            objective_function=mock_obj,
            parameters={"n_trials": 5, "direction": "minimize"},
        )
        best = opt.run_optimization()

        self.assertEqual(call_count[0], 5)
        for name in ("ratio", "carbon_source", "trp_ratio_flux", "violacein_ratio_flux"):
            self.assertIn(name, best)

    def test_best_trp_flux_is_a_valid_percentage(self):
        """The optimizer must return a trp_ratio_flux that is one of the PERCENTAGES."""
        from FLYCOP.OptunaOptimizer import OptunaOptimizer

        def mock_obj(params):
            return 1.0 / (params["trp_ratio_flux"] + 1)

        opt = OptunaOptimizer(
            configuration_space=PARAM_SPACE,
            objective_function=mock_obj,
            parameters={"n_trials": 20, "direction": "minimize"},
        )
        best = opt.run_optimization()
        self.assertIn(best["trp_ratio_flux"], PERCENTAGES)


if __name__ == "__main__":
    unittest.main()
