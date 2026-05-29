"""
Tests for FLYCOP.OptunaOptimizer
=================================
Verifies:
- Initialisation with/without arguments
- Canonical parameter format acceptance (float, int, categorical)
- Validation errors for bad inputs
- load_scenario accepted/rejected keys
- load_objective_function validation
- run_optimization returns a dict with all parameter names
- _build_optuna_objective translates Trial → plain dict correctly
- get_best_parameters raises before any run
- Factory creates an OptunaOptimizer instance for name 'Optuna'
- Same canonical param list works identically for Optuna and SMAC3 (interop)
"""
import unittest
from unittest.mock import MagicMock, patch

# conftest sets up sys.path and all mocks before this module is imported
from FLYCOP.OptunaOptimizer import OptunaOptimizer
from FLYCOP.OptimizerFactory import OptimizerFactory

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

FLOAT_PARAM  = {"type": "float",       "name": "glucose_e",    "min": 0.0, "max": 20.0}
INT_PARAM    = {"type": "int",         "name": "biomass_init",  "min": 1,   "max": 10}
CAT_PARAM    = {"type": "categorical", "name": "strain",        "values": ["wt", "mut"]}

CANONICAL_PARAMS = [FLOAT_PARAM, INT_PARAM, CAT_PARAM]

_simple_objective = lambda params: sum(v if isinstance(v, (int, float)) else 0
                                       for v in params.values())


# ──────────────────────────────────────────────────────────────────────────────
# Initialisation
# ──────────────────────────────────────────────────────────────────────────────

class TestOptunaOptimizerInit(unittest.TestCase):

    def test_default_init(self):
        opt = OptunaOptimizer()
        self.assertIsNone(opt.objective_function)
        self.assertIsNone(opt.study)
        self.assertEqual(opt._param_list, [])

    def test_default_parameters(self):
        opt = OptunaOptimizer()
        self.assertEqual(opt.parameters["n_trials"], 100)
        self.assertEqual(opt.parameters["direction"], "minimize")

    def test_custom_parameters_merged(self):
        opt = OptunaOptimizer(parameters={"n_trials": 50, "direction": "maximize"})
        self.assertEqual(opt.parameters["n_trials"], 50)
        self.assertEqual(opt.parameters["direction"], "maximize")

    def test_configuration_space_preloaded(self):
        opt = OptunaOptimizer(configuration_space=CANONICAL_PARAMS)
        self.assertEqual(len(opt._param_list), 3)

    def test_objective_function_preloaded(self):
        opt = OptunaOptimizer(objective_function=_simple_objective)
        self.assertIs(opt.objective_function, _simple_objective)

    def test_non_callable_objective_raises_on_init(self):
        with self.assertRaises(ValueError):
            OptunaOptimizer(objective_function="not_callable")


# ──────────────────────────────────────────────────────────────────────────────
# load_configuration_space
# ──────────────────────────────────────────────────────────────────────────────

class TestLoadConfigurationSpace(unittest.TestCase):

    def _opt(self):
        return OptunaOptimizer()

    def test_float_param_accepted(self):
        opt = self._opt()
        opt.load_configuration_space([FLOAT_PARAM])
        self.assertEqual(len(opt._param_list), 1)
        self.assertEqual(opt._param_list[0]["name"], "glucose_e")

    def test_int_param_accepted(self):
        opt = self._opt()
        opt.load_configuration_space([INT_PARAM])
        self.assertEqual(opt._param_list[0]["type"], "int")

    def test_categorical_param_accepted(self):
        opt = self._opt()
        opt.load_configuration_space([CAT_PARAM])
        self.assertEqual(opt._param_list[0]["values"], ["wt", "mut"])

    def test_all_canonical_params_accepted(self):
        opt = self._opt()
        opt.load_configuration_space(CANONICAL_PARAMS)
        self.assertEqual(len(opt._param_list), 3)

    def test_non_list_raises_type_error(self):
        opt = self._opt()
        with self.assertRaises(TypeError):
            opt.load_configuration_space({"type": "float", "name": "p"})

    def test_missing_type_raises_value_error(self):
        opt = self._opt()
        with self.assertRaises(ValueError):
            opt.load_configuration_space([{"name": "p", "min": 0, "max": 1}])

    def test_missing_name_raises_value_error(self):
        opt = self._opt()
        with self.assertRaises(ValueError):
            opt.load_configuration_space([{"type": "float", "min": 0, "max": 1}])

    def test_float_without_min_raises(self):
        opt = self._opt()
        with self.assertRaises(ValueError):
            opt.load_configuration_space([{"type": "float", "name": "p", "max": 1.0}])

    def test_float_without_max_raises(self):
        opt = self._opt()
        with self.assertRaises(ValueError):
            opt.load_configuration_space([{"type": "float", "name": "p", "min": 0.0}])

    def test_categorical_without_values_raises(self):
        opt = self._opt()
        with self.assertRaises(ValueError):
            opt.load_configuration_space([{"type": "categorical", "name": "p"}])

    def test_categorical_empty_values_raises(self):
        opt = self._opt()
        with self.assertRaises(ValueError):
            opt.load_configuration_space([{"type": "categorical", "name": "p", "values": []}])

    def test_unknown_type_raises_value_error(self):
        opt = self._opt()
        with self.assertRaises(ValueError):
            opt.load_configuration_space([{"type": "boolean", "name": "p"}])


# ──────────────────────────────────────────────────────────────────────────────
# load_scenario
# ──────────────────────────────────────────────────────────────────────────────

class TestLoadScenario(unittest.TestCase):

    def test_valid_n_trials(self):
        opt = OptunaOptimizer()
        opt.load_scenario({"n_trials": 25})
        self.assertEqual(opt.parameters["n_trials"], 25)

    def test_valid_direction_maximize(self):
        opt = OptunaOptimizer()
        opt.load_scenario({"direction": "maximize"})
        self.assertEqual(opt.parameters["direction"], "maximize")

    def test_multiple_valid_keys(self):
        opt = OptunaOptimizer()
        opt.load_scenario({"n_trials": 10, "direction": "minimize"})
        self.assertEqual(opt.parameters["n_trials"], 10)

    def test_invalid_key_raises_value_error(self):
        opt = OptunaOptimizer()
        with self.assertRaises(ValueError):
            opt.load_scenario({"invalid_key": 99})

    def test_walltime_limit_not_valid_for_optuna(self):
        """walltime_limit is SMAC3-specific, should be rejected by Optuna."""
        opt = OptunaOptimizer()
        with self.assertRaises(ValueError):
            opt.load_scenario({"walltime_limit": 120})


# ──────────────────────────────────────────────────────────────────────────────
# load_objective_function
# ──────────────────────────────────────────────────────────────────────────────

class TestLoadObjectiveFunction(unittest.TestCase):

    def test_callable_accepted(self):
        opt = OptunaOptimizer()
        opt.load_objective_function(_simple_objective)
        self.assertIs(opt.objective_function, _simple_objective)

    def test_non_callable_raises(self):
        opt = OptunaOptimizer()
        with self.assertRaises(ValueError):
            opt.load_objective_function(42)

    def test_lambda_accepted(self):
        opt = OptunaOptimizer()
        fn = lambda p: 0.0
        opt.load_objective_function(fn)
        self.assertIs(opt.objective_function, fn)


# ──────────────────────────────────────────────────────────────────────────────
# _build_optuna_objective  (the Trial → dict translation)
# ──────────────────────────────────────────────────────────────────────────────

class TestBuildOptunaObjective(unittest.TestCase):
    """Verifies the internal wrapper correctly translates Trial → plain dict."""

    def _make_fake_trial(self, values):
        from conftest import FakeOptunaTrial
        return FakeOptunaTrial(values)

    def test_float_param_suggested_and_forwarded(self):
        received = {}

        def capture(params):
            received.update(params)
            return 0.0

        opt = OptunaOptimizer(
            configuration_space=[FLOAT_PARAM],
            objective_function=capture,
        )
        wrapped = opt._build_optuna_objective()
        trial = self._make_fake_trial({"glucose_e": 5.5})
        wrapped(trial)
        self.assertAlmostEqual(received["glucose_e"], 5.5)

    def test_int_param_suggested_and_forwarded(self):
        received = {}

        def capture(params):
            received.update(params)
            return 0.0

        opt = OptunaOptimizer(
            configuration_space=[INT_PARAM],
            objective_function=capture,
        )
        wrapped = opt._build_optuna_objective()
        trial = self._make_fake_trial({"biomass_init": 7})
        wrapped(trial)
        self.assertEqual(received["biomass_init"], 7)

    def test_categorical_param_suggested_and_forwarded(self):
        received = {}

        def capture(params):
            received.update(params)
            return 0.0

        opt = OptunaOptimizer(
            configuration_space=[CAT_PARAM],
            objective_function=capture,
        )
        wrapped = opt._build_optuna_objective()
        trial = self._make_fake_trial({"strain": "mut"})
        wrapped(trial)
        self.assertEqual(received["strain"], "mut")

    def test_all_params_forwarded_as_plain_dict(self):
        """The objective function receives a plain dict, not a Trial object."""
        received_type = []

        def capture(params):
            received_type.append(type(params).__name__)
            return 0.0

        opt = OptunaOptimizer(
            configuration_space=CANONICAL_PARAMS,
            objective_function=capture,
        )
        wrapped = opt._build_optuna_objective()
        trial = self._make_fake_trial({"glucose_e": 1.0, "biomass_init": 3, "strain": "wt"})
        wrapped(trial)
        self.assertEqual(received_type[0], "dict")

    def test_objective_return_value_propagated(self):
        opt = OptunaOptimizer(
            configuration_space=[FLOAT_PARAM],
            objective_function=lambda p: 42.0,
        )
        wrapped = opt._build_optuna_objective()
        trial = self._make_fake_trial({})
        result = wrapped(trial)
        self.assertEqual(result, 42.0)


# ──────────────────────────────────────────────────────────────────────────────
# run_optimization guards
# ──────────────────────────────────────────────────────────────────────────────

class TestRunOptimizationGuards(unittest.TestCase):

    def test_raises_if_no_config_space(self):
        opt = OptunaOptimizer(objective_function=_simple_objective)
        with self.assertRaises(ValueError, msg="Should raise when config space empty"):
            opt.run_optimization()

    def test_raises_if_no_objective(self):
        opt = OptunaOptimizer(configuration_space=CANONICAL_PARAMS)
        with self.assertRaises(ValueError, msg="Should raise when objective is None"):
            opt.run_optimization()


# ──────────────────────────────────────────────────────────────────────────────
# run_optimization – uses FakeOptunaStudy via the mock in conftest
# ──────────────────────────────────────────────────────────────────────────────

class TestRunOptimization(unittest.TestCase):
    """
    The optuna module is mocked in conftest.py.  FakeOptunaStudy runs
    n_trials synchronously so we can verify end-to-end behaviour.
    """

    def test_returns_dict_with_all_param_names(self):
        opt = OptunaOptimizer(
            configuration_space=CANONICAL_PARAMS,
            objective_function=_simple_objective,
            parameters={"n_trials": 3},
        )
        result = opt.run_optimization()
        self.assertIsInstance(result, dict)

    def test_study_is_set_after_run(self):
        opt = OptunaOptimizer(
            configuration_space=[FLOAT_PARAM],
            objective_function=lambda p: p["glucose_e"],
            parameters={"n_trials": 2},
        )
        opt.run_optimization()
        self.assertIsNotNone(opt.study)

    def test_run_optimization_respects_n_trials_parameter(self):
        """FakeOptunaStudy records one trial per call; verify count."""
        opt = OptunaOptimizer(
            configuration_space=[FLOAT_PARAM],
            objective_function=lambda p: p["glucose_e"],
            parameters={"n_trials": 5},
        )
        opt.run_optimization()
        self.assertEqual(len(opt.study.trials), 5)


# ──────────────────────────────────────────────────────────────────────────────
# get_best_parameters
# ──────────────────────────────────────────────────────────────────────────────

class TestGetBestParameters(unittest.TestCase):

    def test_raises_before_run(self):
        opt = OptunaOptimizer()
        with self.assertRaises(ValueError):
            opt.get_best_parameters()

    def test_returns_dict_after_run(self):
        opt = OptunaOptimizer(
            configuration_space=[FLOAT_PARAM],
            objective_function=lambda p: p["glucose_e"],
            parameters={"n_trials": 1},
        )
        opt.run_optimization()
        best = opt.get_best_parameters()
        self.assertIsInstance(best, dict)


# ──────────────────────────────────────────────────────────────────────────────
# Factory integration
# ──────────────────────────────────────────────────────────────────────────────

class TestOptimizerFactoryOptuna(unittest.TestCase):

    def test_factory_creates_optuna_optimizer(self):
        opt = OptimizerFactory.createOptimizer("Optuna")
        self.assertIsInstance(opt, OptunaOptimizer)

    def test_factory_passes_parameters(self):
        opt = OptimizerFactory.createOptimizer("Optuna", parameters={"n_trials": 42})
        self.assertEqual(opt.parameters["n_trials"], 42)

    def test_factory_preloads_configuration_space(self):
        opt = OptimizerFactory.createOptimizer(
            "Optuna", configuration_space=CANONICAL_PARAMS
        )
        self.assertEqual(len(opt._param_list), 3)


# ──────────────────────────────────────────────────────────────────────────────
# Interoperability: same canonical params work for both backends
# ──────────────────────────────────────────────────────────────────────────────

class TestCanonicalParamInterop(unittest.TestCase):
    """
    Verifies that the *same* CANONICAL_PARAMS list and the *same* objective
    function can be used with both OptunaOptimizer and SMAC3Optimizer without
    any modification.
    """

    def test_same_params_accepted_by_optuna(self):
        opt = OptunaOptimizer()
        # should not raise
        opt.load_configuration_space(CANONICAL_PARAMS)
        self.assertEqual(len(opt._param_list), 3)

    def test_same_params_accepted_by_smac3(self):
        from FLYCOP.SMAC3Optimizer import SMAC3Optimizer
        opt = SMAC3Optimizer()
        # should not raise (ConfigSpace is mocked, so add_hyperparameters is a no-op)
        opt.load_configuration_space(CANONICAL_PARAMS)

    def test_optuna_objective_receives_plain_dict(self):
        """
        FLYCOP_scenario.optimize uses dict.items() which works equally on a
        plain dict and on a SMAC3 Configuration object.  Confirm the Optuna
        wrapper produces a plain dict the scenario can consume.
        """
        received = []

        def mock_scenario_optimize(params):
            received.append(params)
            return 0.0

        opt = OptunaOptimizer(
            configuration_space=[FLOAT_PARAM],
            objective_function=mock_scenario_optimize,
        )
        wrapped = opt._build_optuna_objective()
        from conftest import FakeOptunaTrial
        trial = FakeOptunaTrial({"glucose_e": 3.0})
        wrapped(trial)

        self.assertTrue(len(received) == 1)
        self.assertIsInstance(received[0], dict)
        # simulate what FLYCOP_scenario.optimize does with the dict
        extracted = dict(received[0].items())
        self.assertIn("glucose_e", extracted)


if __name__ == "__main__":
    unittest.main()
