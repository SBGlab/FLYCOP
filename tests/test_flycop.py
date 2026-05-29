"""
Tests for FLYCOP.FLYCOP (FLYCOP class) and FLYCOP.FLYCOP_scenario
Covers: proper initialization of attributes, create_simulator/optimizator,
        FLYCOP_scenario.optimize using dict instead of locals(), set_scenario.
"""
import unittest
from unittest.mock import MagicMock, patch, call

from FLYCOP.FLYCOP import FLYCOP, FLYCOP_scenario


class TestFLYCOPInit(unittest.TestCase):
    """FLYCOP.__init__ must correctly set all attributes (bugs #1, #2)."""

    @patch("FLYCOP.FLYCOP.SimulatorFactory")
    @patch("FLYCOP.FLYCOP.OptimizerFactory")
    def test_simulator_is_not_none_after_init(self, mock_opt_factory, mock_sim_factory):
        fake_sim = MagicMock()
        mock_sim_factory.return_value.createSimulator.return_value = fake_sim
        fake_opt = MagicMock()
        mock_opt_factory.createOptimizer.return_value = fake_opt

        flycop = FLYCOP(simulator_type="COMETS", optimizator_type="SMAC3")

        self.assertIsNotNone(flycop.simulator)
        self.assertIs(flycop.simulator, fake_sim)

    @patch("FLYCOP.FLYCOP.SimulatorFactory")
    @patch("FLYCOP.FLYCOP.OptimizerFactory")
    def test_optimizator_is_not_none_after_init(self, mock_opt_factory, mock_sim_factory):
        mock_sim_factory.return_value.createSimulator.return_value = MagicMock()
        fake_opt = MagicMock()
        mock_opt_factory.createOptimizer.return_value = fake_opt

        flycop = FLYCOP(simulator_type="COMETS", optimizator_type="SMAC3")

        self.assertIsNotNone(flycop.optimizator)
        self.assertIs(flycop.optimizator, fake_opt)

    @patch("FLYCOP.FLYCOP.SimulatorFactory")
    @patch("FLYCOP.FLYCOP.OptimizerFactory")
    def test_objective_function_is_initialized(self, mock_opt_factory, mock_sim_factory):
        mock_sim_factory.return_value.createSimulator.return_value = MagicMock()
        mock_opt_factory.createOptimizer.return_value = MagicMock()

        flycop = FLYCOP()

        self.assertIsNotNone(flycop.objective_function)

    @patch("FLYCOP.FLYCOP.SimulatorFactory")
    @patch("FLYCOP.FLYCOP.OptimizerFactory")
    def test_scenario_created_with_simulator(self, mock_opt_factory, mock_sim_factory):
        fake_sim = MagicMock()
        mock_sim_factory.return_value.createSimulator.return_value = fake_sim
        mock_opt_factory.createOptimizer.return_value = MagicMock()

        flycop = FLYCOP()

        self.assertIsInstance(flycop.scenario, FLYCOP_scenario)
        self.assertIs(flycop.scenario.simulator, fake_sim)


class TestFLYCOPMethods(unittest.TestCase):
    @patch("FLYCOP.FLYCOP.SimulatorFactory")
    @patch("FLYCOP.FLYCOP.OptimizerFactory")
    def setUp(self, mock_opt_factory, mock_sim_factory):
        self.fake_sim = MagicMock()
        mock_sim_factory.return_value.createSimulator.return_value = self.fake_sim
        self.fake_opt = MagicMock()
        mock_opt_factory.createOptimizer.return_value = self.fake_opt
        self.flycop = FLYCOP()

    def test_set_simulator_replaces_simulator(self):
        new_sim = MagicMock()
        self.flycop.set_simulator(new_sim)
        self.assertIs(self.flycop.simulator, new_sim)

    def test_load_parameters_stores_dict(self):
        params = {"p1": 1.0, "p2": 2.0}
        self.flycop.load_parameters(params)
        self.assertEqual(self.flycop.parameters, params)

    def test_simulate_delegates_to_simulator(self):
        params = {"x": 1}
        self.flycop.simulate(params)
        self.fake_sim.simulate.assert_called_once_with(params)


class TestFLYCOPScenarioInit(unittest.TestCase):
    def test_consortia_initialized_to_none(self):
        """Bug #5: FLYCOP_scenario must initialize self.consortia."""
        sim = MagicMock()
        scenario = FLYCOP_scenario(sim)
        self.assertIsNone(scenario.consortia)

    def test_simulator_stored(self):
        sim = MagicMock()
        scenario = FLYCOP_scenario(sim)
        self.assertIs(scenario.simulator, sim)

    def test_fitness_defaults_to_none(self):
        sim = MagicMock()
        scenario = FLYCOP_scenario(sim)
        self.assertIsNone(scenario.fitness)


class TestFLYCOPScenarioSetScenario(unittest.TestCase):
    def test_set_scenario_assigns_directly(self):
        """Bug #4: set_scenario was using __get__ descriptor incorrectly."""
        sim = MagicMock()
        scenario = FLYCOP_scenario(sim)
        new_fn = MagicMock()
        scenario.set_scenario(new_fn)
        self.assertIs(scenario.scenario, new_fn)


class TestFLYCOPScenarioOptimize(unittest.TestCase):
    """Bug #3: optimize() was using locals() to extract config values (no-op)."""

    def _make_scenario(self):
        sim = MagicMock()
        fitness = MagicMock()
        fitness.execute_fitness_function.return_value = 0.5
        s = FLYCOP_scenario(sim, fitness=fitness)
        s.consortia = []
        return s

    def test_optimize_calls_simulate(self):
        scenario = self._make_scenario()
        config = {"biomass_ecoli": 0.1, "glc_e": 10.0}
        scenario.optimize(config)
        scenario.simulator.simulate.assert_called_once()

    def test_optimize_passes_biomass_variables_to_load_consortia(self):
        scenario = self._make_scenario()
        config = {"biomass_ecoli": 0.1, "glc_e": 10.0}
        scenario.optimize(config)
        # biomass_ecoli starts with 'biomass' → biomass_variables = [0.1]
        scenario.simulator.load_consortia.assert_called_once_with([], [0.1])

    def test_optimize_passes_media_variables_to_load_media(self):
        scenario = self._make_scenario()
        config = {"glc_e": 10.0, "ac_e": 5.0}
        scenario.optimize(config)
        scenario.simulator.load_media.assert_called_once_with({"glc_e": 10.0, "ac_e": 5.0})

    def test_optimize_no_biomass_calls_load_consortia_without_biomass(self):
        scenario = self._make_scenario()
        config = {"glc_e": 10.0}
        scenario.optimize(config)
        scenario.simulator.load_consortia.assert_called_once_with([])

    def test_optimize_returns_fitness_value(self):
        scenario = self._make_scenario()
        config = {}
        result = scenario.optimize(config)
        self.assertEqual(result, 0.5)


if __name__ == "__main__":
    unittest.main()
