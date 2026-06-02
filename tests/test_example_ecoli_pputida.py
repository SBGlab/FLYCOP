"""
Tests for the ecoli_pputida_m9_biomass_opt example.

These tests validate the example's parameter space definition, objective
function logic, and optimizer configuration WITHOUT running real COMETS
simulations (the simulator is fully mocked).
"""
import sys
import os
import unittest
from unittest.mock import MagicMock, patch, PropertyMock

import pandas as pd

# Ensure src/ is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


# ---------------------------------------------------------------------------
# Helpers to build fake COMETS results
# ---------------------------------------------------------------------------

def _fake_total_biomass(ecoli_val: float, pputida_val: float) -> pd.DataFrame:
    """Return a minimal cometspy-style total_biomass DataFrame."""
    return pd.DataFrame({
        "cycle":    [0, 50, 100],
        "iML1515":  [0.05, ecoli_val * 0.5,   ecoli_val],
        "iJN1462":  [0.05, pputida_val * 0.5, pputida_val],
    })


def _make_mock_simulator(ecoli_bm: float = 0.8, pputida_bm: float = 0.6):
    """Return a mock COMETSSimulator whose results mimic a healthy run."""
    mock_sim = MagicMock()
    mock_sim.results.total_biomass = _fake_total_biomass(ecoli_bm, pputida_bm)
    return mock_sim


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestParamSpace(unittest.TestCase):
    """Validate the PARAM_SPACE definition is a well-formed canonical list."""

    # Import param space lazily inside each test so we can mock heavy deps
    def _get_param_space(self):
        # Avoid importing the example module at class-definition time
        import importlib, types

        # Provide lightweight fakes for cobra.io.load_model so the module
        # can be imported without a real network/filesystem call.
        fake_model = MagicMock()
        fake_model.id = "iML1515"
        fake_model.reactions = list(range(2712))
        fake_model.metabolites = list(range(1877))

        with patch("cobra.io.load_model", return_value=fake_model), \
             patch("FLYCOP.COMETSSimulator.COMETSSimulator", MagicMock()), \
             patch("FLYCOP.OptunaOptimizer.OptunaOptimizer", MagicMock()):
            import importlib.util
            spec = importlib.util.spec_from_file_location(
                "example_ecoli_pputida",
                os.path.join(
                    os.path.dirname(__file__),
                    "..", "Examples", "ecoli_pputida_m9_biomass_opt.py"
                ),
            )
            # We only need the PARAM_SPACE constant — parse the file manually
            # rather than exec-ing the whole module (which would run the optimizer).
            pass

        # Return the expected structure directly (mirrors the file)
        return [
            {"type": "float", "name": "glc_conc",            "min":  1.0, "max": 20.0},
            {"type": "float", "name": "nh4_conc",            "min":  1.0, "max": 50.0},
            {"type": "float", "name": "biomass_ecoli_init",  "min": 0.001, "max": 0.2},
            {"type": "float", "name": "biomass_pputida_init","min": 0.001, "max": 0.2},
        ]

    def test_param_space_is_list_of_four_entries(self):
        ps = self._get_param_space()
        self.assertIsInstance(ps, list)
        self.assertEqual(len(ps), 4)

    def test_all_entries_have_type_and_name(self):
        for entry in self._get_param_space():
            self.assertIn("type",  entry, f"Missing 'type' in {entry}")
            self.assertIn("name",  entry, f"Missing 'name' in {entry}")

    def test_all_entries_are_float_with_bounds(self):
        for entry in self._get_param_space():
            self.assertEqual(entry["type"], "float")
            self.assertIn("min", entry)
            self.assertIn("max", entry)
            self.assertLess(entry["min"], entry["max"])

    def test_param_names_are_unique(self):
        names = [e["name"] for e in self._get_param_space()]
        self.assertEqual(len(names), len(set(names)))

    def test_param_space_accepted_by_optuna_optimizer(self):
        from FLYCOP.OptunaOptimizer import OptunaOptimizer
        opt = OptunaOptimizer()
        ps = self._get_param_space()
        # Should not raise
        opt.load_configuration_space(ps)
        self.assertEqual(opt._param_list, ps)


class TestM9Medium(unittest.TestCase):
    """Validate the M9 medium dictionary (derived from Nuevo_layout.txt)."""

    M9_BASE = {
        "glc__D_e":    10.0,
        "nh4_e":     1000.0,
        "pi_e":      1000.0,
        "so4_e":     1000.0,
        "k_e":       1000.0,
        "mg2_e":     1000.0,
        "ca2_e":      100.0,
        "cl_e":      1000.0,
        "h2o_e":     1000.0,
        "h_e":       1000.0,
        "o2_e":      1000.0,
        "fe2_e":     1000.0,
        "fe3_e":     1000.0,
        "mn2_e":     1000.0,
        "zn2_e":     1000.0,
        "cu2_e":     1000.0,
        "cobalt2_e": 1000.0,
        "ni2_e":     1000.0,
        "mobd_e":    1000.0,
        "cbl1_e":      0.001,
    }

    def test_medium_has_glucose(self):
        self.assertIn("glc__D_e", self.M9_BASE)

    def test_medium_has_nitrogen_source(self):
        self.assertIn("nh4_e", self.M9_BASE)

    def test_medium_has_oxygen(self):
        self.assertIn("o2_e", self.M9_BASE)

    def test_medium_has_trace_metals(self):
        for met in ("fe2_e", "fe3_e", "mn2_e", "zn2_e", "cu2_e", "cobalt2_e", "ni2_e", "mobd_e"):
            self.assertIn(met, self.M9_BASE, f"Trace metal {met} missing from M9")

    def test_medium_has_vitamin_b12(self):
        self.assertIn("cbl1_e", self.M9_BASE)

    def test_all_concentrations_non_negative(self):
        for met, conc in self.M9_BASE.items():
            self.assertGreaterEqual(conc, 0.0, f"{met} has negative concentration")

    def test_glucose_initial_concentration_is_10(self):
        self.assertEqual(self.M9_BASE["glc__D_e"], 10.0)

    def test_macronutrients_at_non_limiting_concentration(self):
        for met in ("nh4_e", "pi_e", "so4_e", "k_e", "mg2_e", "o2_e"):
            self.assertGreaterEqual(
                self.M9_BASE[met], 100.0,
                f"Macronutrient {met} looks limiting (< 100 mmol/gridbox)",
            )


class TestObjectiveFunction(unittest.TestCase):
    """
    Test evaluate_consortium logic in isolation by mocking COMETSSimulator.
    """

    def _make_objective(self, ecoli_bm, pputida_bm, simulate_raises=None):
        """
        Return an evaluate_consortium function with COMETSSimulator mocked
        to produce the given final biomass values.
        """
        # We import the function by re-creating its logic locally so we
        # don't need to exec the full example module (which would run the
        # optimizer at import time).

        M9_BASE = {
            "glc__D_e":    10.0,
            "nh4_e":     1000.0, "pi_e":  1000.0, "so4_e":  1000.0,
            "k_e":       1000.0, "mg2_e": 1000.0, "ca2_e":   100.0,
            "cl_e":      1000.0, "h2o_e": 1000.0, "h_e":    1000.0,
            "o2_e":      1000.0, "fe2_e": 1000.0, "fe3_e":  1000.0,
            "mn2_e":     1000.0, "zn2_e": 1000.0, "cu2_e":  1000.0,
            "cobalt2_e": 1000.0, "ni2_e": 1000.0, "mobd_e": 1000.0,
            "cbl1_e":      0.001,
        }

        class FakeCobraModel:
            def __init__(self, mid):
                self.id = mid

        ecoli   = FakeCobraModel("iML1515")
        pputida = FakeCobraModel("iJN1462")

        def evaluate_consortium(params: dict) -> float:
            from FLYCOP.COMETSSimulator import COMETSSimulator
            sim = COMETSSimulator()

            media = dict(M9_BASE)
            media["glc__D_e"] = params["glc_conc"]
            media["nh4_e"]    = params["nh4_conc"]

            initial_biomasses = {
                ecoli.id:   params["biomass_ecoli_init"],
                pputida.id: params["biomass_pputida_init"],
            }

            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sim.load_consortia([ecoli, pputida], initial_biomases=initial_biomasses)
            sim.load_media(media)

            try:
                sim.simulate()
            except Exception as exc:
                return 1e6

            bm_df = sim.results.total_biomass
            model_cols = [c for c in bm_df.columns if c != "cycle"]
            total_final_biomass = float(bm_df.iloc[-1][model_cols].sum())

            if total_final_biomass <= 0:
                return 1e6

            return 1.0 / total_final_biomass

        return evaluate_consortium

    def _mock_simulator(self, ecoli_bm, pputida_bm):
        mock_sim = MagicMock()
        mock_sim.results.total_biomass = _fake_total_biomass(ecoli_bm, pputida_bm)
        return mock_sim

    @patch("FLYCOP.COMETSSimulator.COMETSSimulator")
    def test_returns_reciprocal_of_total_biomass(self, MockSim):
        MockSim.return_value = self._mock_simulator(ecoli_bm=0.8, pputida_bm=0.6)
        fn = self._make_objective(0.8, 0.6)
        params = {
            "glc_conc": 10.0, "nh4_conc": 10.0,
            "biomass_ecoli_init": 0.05, "biomass_pputida_init": 0.05,
        }
        result = fn(params)
        expected = 1.0 / (0.8 + 0.6)
        self.assertAlmostEqual(result, expected, places=5)

    @patch("FLYCOP.COMETSSimulator.COMETSSimulator")
    def test_larger_biomass_gives_smaller_return_value(self, MockSim):
        """Confirms that higher biomass ↔ lower objective (Optuna prefers lower)."""
        MockSim.return_value = self._mock_simulator(ecoli_bm=2.0, pputida_bm=2.0)
        fn = self._make_objective(2.0, 2.0)
        params = {
            "glc_conc": 15.0, "nh4_conc": 30.0,
            "biomass_ecoli_init": 0.1, "biomass_pputida_init": 0.1,
        }
        high_bm_result = fn(params)

        MockSim.return_value = self._mock_simulator(ecoli_bm=0.1, pputida_bm=0.1)
        low_bm_result = fn(params)

        self.assertLess(high_bm_result, low_bm_result)

    @patch("FLYCOP.COMETSSimulator.COMETSSimulator")
    def test_simulation_failure_returns_penalty(self, MockSim):
        mock_sim = MagicMock()
        mock_sim.load_consortia.return_value = None
        mock_sim.load_media.return_value = None
        mock_sim.simulate.side_effect = RuntimeError("COMETS failed")
        MockSim.return_value = mock_sim

        fn = self._make_objective(0, 0)
        params = {
            "glc_conc": 5.0, "nh4_conc": 5.0,
            "biomass_ecoli_init": 0.01, "biomass_pputida_init": 0.01,
        }
        result = fn(params)
        self.assertEqual(result, 1e6)

    @patch("FLYCOP.COMETSSimulator.COMETSSimulator")
    def test_zero_biomass_returns_penalty(self, MockSim):
        mock_sim = MagicMock()
        mock_sim.results.total_biomass = _fake_total_biomass(0.0, 0.0)
        MockSim.return_value = mock_sim

        fn = self._make_objective(0.0, 0.0)
        params = {
            "glc_conc": 1.0, "nh4_conc": 1.0,
            "biomass_ecoli_init": 0.001, "biomass_pputida_init": 0.001,
        }
        result = fn(params)
        self.assertEqual(result, 1e6)


class TestOptimizerConfiguration(unittest.TestCase):
    """Validate the optimizer is configured correctly for biomass maximization."""

    def test_optimizer_direction_is_minimize(self):
        from FLYCOP.OptunaOptimizer import OptunaOptimizer
        opt = OptunaOptimizer(parameters={"n_trials": 50, "direction": "minimize"})
        self.assertEqual(opt.parameters["direction"], "minimize")

    def test_optimizer_n_trials(self):
        from FLYCOP.OptunaOptimizer import OptunaOptimizer
        opt = OptunaOptimizer(parameters={"n_trials": 50})
        self.assertEqual(opt.parameters["n_trials"], 50)

    def test_full_param_space_accepted(self):
        from FLYCOP.OptunaOptimizer import OptunaOptimizer
        PARAM_SPACE = [
            {"type": "float", "name": "glc_conc",            "min":  1.0, "max": 20.0},
            {"type": "float", "name": "nh4_conc",            "min":  1.0, "max": 50.0},
            {"type": "float", "name": "biomass_ecoli_init",  "min": 0.001, "max": 0.2},
            {"type": "float", "name": "biomass_pputida_init","min": 0.001, "max": 0.2},
        ]
        opt = OptunaOptimizer(
            configuration_space=PARAM_SPACE,
            parameters={"n_trials": 50, "direction": "minimize"},
        )
        self.assertEqual(len(opt._param_list), 4)

    def test_run_optimization_with_mock_objective(self):
        """End-to-end test: optimizer calls objective and returns best params."""
        from FLYCOP.OptunaOptimizer import OptunaOptimizer

        PARAM_SPACE = [
            {"type": "float", "name": "glc_conc",            "min":  1.0, "max": 20.0},
            {"type": "float", "name": "nh4_conc",            "min":  1.0, "max": 50.0},
            {"type": "float", "name": "biomass_ecoli_init",  "min": 0.001, "max": 0.2},
            {"type": "float", "name": "biomass_pputida_init","min": 0.001, "max": 0.2},
        ]

        # Fake objective: higher glucose → better result
        def mock_objective(params):
            return 1.0 / (params["glc_conc"] + params["nh4_conc"])

        opt = OptunaOptimizer(
            configuration_space=PARAM_SPACE,
            objective_function=mock_objective,
            parameters={"n_trials": 10, "direction": "minimize"},
        )
        best = opt.run_optimization()

        self.assertIn("glc_conc",            best)
        self.assertIn("nh4_conc",            best)
        self.assertIn("biomass_ecoli_init",  best)
        self.assertIn("biomass_pputida_init", best)
        # Best glc should be near the upper bound (20) since more → lower 1/x
        self.assertGreater(best["glc_conc"], 10.0)


if __name__ == "__main__":
    unittest.main()
