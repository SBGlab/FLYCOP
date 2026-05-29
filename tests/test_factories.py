"""
Tests for FLYCOP.SimulatorFactory and FLYCOP.OptimizerFactory
Covers: valid simulator/optimizer creation and error handling for unknown types.
"""
import unittest
from unittest.mock import MagicMock, patch

from FLYCOP.SimulatorFactory import SimulatorFactory
from FLYCOP.OptimizerFactory import OptimizerFactory


class TestSimulatorFactory(unittest.TestCase):
    @patch("FLYCOP.SimulatorFactory.COMETSSimulator", create=True)
    def test_create_comets_simulator(self, _):
        with patch.dict("sys.modules", {"FLYCOP.COMETSSimulator": MagicMock()}):
            # createSimulator returns the instance from COMETSSimulator()
            sim = SimulatorFactory.createSimulator("COMETS")
            self.assertIsNotNone(sim)

    def test_unknown_simulator_raises_value_error(self):
        with self.assertRaises(ValueError):
            SimulatorFactory.createSimulator("UNKNOWN_SIM")

    def test_dFBA_simulator_creation(self):
        sim = SimulatorFactory.createSimulator("dFBA")
        self.assertIsNotNone(sim)

    def test_all_known_types_do_not_raise(self):
        """Each known simulator type should return something without raising."""
        for sim_type in ("dFBA",):
            with self.subTest(sim_type=sim_type):
                sim = SimulatorFactory.createSimulator(sim_type)
                self.assertIsNotNone(sim)


class TestOptimizerFactory(unittest.TestCase):
    def test_create_smac3_optimizer(self):
        opt = OptimizerFactory.createOptimizer("SMAC3")
        self.assertIsNotNone(opt)

    def test_create_skopt_optimizer(self):
        opt = OptimizerFactory.createOptimizer("SKopt")
        self.assertIsNotNone(opt)

    def test_create_optuna_optimizer(self):
        from FLYCOP.OptunaOptimizer import OptunaOptimizer
        opt = OptimizerFactory.createOptimizer("Optuna")
        self.assertIsInstance(opt, OptunaOptimizer)

    def test_unknown_optimizer_raises_value_error(self):
        with self.assertRaises(ValueError):
            OptimizerFactory.createOptimizer("UNKNOWN_OPT")


if __name__ == "__main__":
    unittest.main()
