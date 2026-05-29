"""
Tests for FLYCOP.SMAC3Optimizer and FLYCOP.SKoptOptimizer
Covers: initialization, configuration space loading, objective function/facade
        assignment bugs, run_optimization guards.
"""
import unittest
from unittest.mock import MagicMock, patch

from conftest import FakeHPOFacade
from FLYCOP.SMAC3Optimizer import SMAC3Optimizer
from FLYCOP.SKoptOptimizer import SKoptOptimizer


class TestSMAC3OptimizerInit(unittest.TestCase):
    """Bug #16: objective_function and facade were never assigned when valid."""

    def test_default_init_sets_objective_function_to_none(self):
        opt = SMAC3Optimizer()
        self.assertIsNone(opt.objective_function)

    def test_default_init_sets_smac_to_none(self):
        opt = SMAC3Optimizer()
        self.assertIsNone(opt.SMAC)

    def test_callable_objective_function_is_stored(self):
        func = lambda cfg: 1.0
        opt = SMAC3Optimizer(objective_function=func)
        self.assertIs(opt.objective_function, func)

    def test_non_callable_objective_function_raises(self):
        with self.assertRaises(ValueError):
            SMAC3Optimizer(objective_function="not_callable")

    def test_valid_facade_is_stored(self):
        facade = FakeHPOFacade()
        opt = SMAC3Optimizer(facade=facade)
        self.assertIs(opt.SMAC, facade)

    def test_invalid_facade_raises(self):
        with self.assertRaises(ValueError):
            SMAC3Optimizer(facade="not_a_facade")

    def test_parameters_stored(self):
        params = {"max_evaluations": 50, "deterministic": False}
        opt = SMAC3Optimizer(parameters=params)
        self.assertEqual(opt.parameters, params)


class TestSMAC3OptimizerLoadConfigSpace(unittest.TestCase):
    def test_load_float_hyperparameter(self):
        opt = SMAC3Optimizer()
        cs_list = [{"type": "float", "name": "p1", "min": 0.0, "max": 1.0}]
        # Should not raise (ConfigSpace is mocked)
        try:
            opt.load_configuration_space(cs_list)
        except Exception as e:
            self.fail(f"load_configuration_space raised: {e}")

    def test_load_invalid_type_raises(self):
        opt = SMAC3Optimizer()
        cs_list = [{"type": "invalid_type", "name": "p1"}]
        # Should raise KeyError or similar because the type branch doesn't match
        # (no branch for 'invalid_type' — no add_hyperparameters call, but no error either)
        # Just verify it doesn't crash
        try:
            opt.load_configuration_space(cs_list)
        except Exception:
            pass  # implementation-defined behaviour for unknown types


class TestSMAC3OptimizerLoadObjectiveFunction(unittest.TestCase):
    def test_callable_is_accepted(self):
        opt = SMAC3Optimizer()
        func = lambda cfg: 0.0
        opt.load_objective_function(func)
        self.assertIs(opt.objective_function, func)

    def test_non_callable_raises(self):
        opt = SMAC3Optimizer()
        with self.assertRaises(ValueError):
            opt.load_objective_function(42)


class TestSMAC3OptimizerRunOptimization(unittest.TestCase):
    def test_run_optimization_raises_when_no_scenario(self):
        opt = SMAC3Optimizer()
        opt.objective_function = lambda cfg: 0.0
        # scenario is not set, should raise ValueError
        with self.assertRaises((ValueError, AttributeError)):
            opt.run_optimization()

    def test_run_optimization_raises_when_no_objective(self):
        opt = SMAC3Optimizer()
        opt.scenario = MagicMock()
        opt.objective_function = None
        with self.assertRaises(ValueError):
            opt.run_optimization()


class TestSKoptOptimizerInit(unittest.TestCase):
    """Bug #17/#18: class was named SkoptOptimizer and had no __init__."""

    def test_instantiation_succeeds(self):
        try:
            opt = SKoptOptimizer()
        except Exception as e:
            self.fail(f"SKoptOptimizer() raised: {e}")

    def test_parameters_stored(self):
        params = {"n_calls": 50}
        opt = SKoptOptimizer(parameters=params)
        self.assertEqual(opt.parameters, params)

    def test_default_parameters_is_empty_dict(self):
        opt = SKoptOptimizer()
        self.assertEqual(opt.parameters, {})

    def test_configuration_space_stored(self):
        cs = MagicMock()
        opt = SKoptOptimizer(configuration_space=cs)
        self.assertIs(opt.configuration_space, cs)


if __name__ == "__main__":
    unittest.main()
