"""
Tests for FLYCOP.Fitness
Covers: initialization, add/remove functions, default management,
        metabolite_maximization, biomass_maximization, execute helpers.
"""
import unittest
import pandas as pd

from FLYCOP.Fitness import Fitness


class TestFitnessInit(unittest.TestCase):
    def setUp(self):
        self.fitness = Fitness()

    def test_default_functions_registered(self):
        self.assertIn("metabolite_maximization", self.fitness.fitness_functions)
        self.assertIn("biomass_maximization", self.fitness.fitness_functions)

    def test_default_fitness_is_metabolite_maximization(self):
        self.assertEqual(self.fitness.default_fitness, "metabolite_maximization")


class TestFitnessManagement(unittest.TestCase):
    def setUp(self):
        self.fitness = Fitness()

    def test_add_fitness_function_stores_it(self):
        func = lambda x: x
        self.fitness.add_fitness_function("custom", func)
        self.assertIn("custom", self.fitness.fitness_functions)
        self.assertIs(self.fitness.fitness_functions["custom"], func)

    def test_add_fitness_function_updates_default(self):
        self.fitness.add_fitness_function("custom", lambda x: x)
        self.assertEqual(self.fitness.default_fitness, "custom")

    def test_remove_existing_fitness_function(self):
        self.fitness.remove_fitness_function("metabolite_maximization")
        self.assertNotIn("metabolite_maximization", self.fitness.fitness_functions)

    def test_remove_nonexistent_does_not_raise(self):
        # Bug was: calling .remove() on dict → should now use del without error
        try:
            self.fitness.remove_fitness_function("does_not_exist")
        except Exception as e:
            self.fail(f"remove_fitness_function raised unexpectedly: {e}")

    def test_set_default_fitness_existing(self):
        self.fitness.set_default_fitness("biomass_maximization")
        self.assertEqual(self.fitness.default_fitness, "biomass_maximization")

    def test_set_default_fitness_nonexistent_does_not_change(self):
        self.fitness.set_default_fitness("nonexistent")
        self.assertEqual(self.fitness.default_fitness, "metabolite_maximization")

    def test_list_fitness_functions_does_not_raise(self):
        try:
            self.fitness.list_fitness_functions()
        except Exception as e:
            self.fail(f"list_fitness_functions raised unexpectedly: {e}")


class TestFitnessFunctions(unittest.TestCase):
    def setUp(self):
        self.fitness = Fitness()

    def test_metabolite_maximization_returns_inverse_of_second_to_last(self):
        df = pd.DataFrame({"met_A": [1.0, 2.0, 4.0, 8.0]})
        result = self.fitness.metabolite_maximization(df, "met_A")
        # iloc[-2] is index 2, value 4.0  →  1/4.0 = 0.25
        self.assertAlmostEqual(result, 1.0 / 4.0)

    def test_metabolite_maximization_with_minimum_two_rows(self):
        df = pd.DataFrame({"met_A": [5.0, 10.0]})
        result = self.fitness.metabolite_maximization(df, "met_A")
        # iloc[-2] is index 0, value 5.0  →  1/5.0
        self.assertAlmostEqual(result, 1.0 / 5.0)

    def test_biomass_maximization_returns_inverse_of_sum(self):
        # Bug was: NameError 'columnas_biomass' — should use 'columns_biomass'
        df = pd.DataFrame({"biomass_A": [1.0, 2.0, 3.0], "biomass_B": [1.0, 2.0, 1.0]})
        result = self.fitness.biomass_maximization(df)
        # iloc[-2] is index 1  →  biomass_A=2.0, biomass_B=2.0  →  1/4.0
        self.assertAlmostEqual(result, 1.0 / 4.0)

    def test_biomass_maximization_single_strain(self):
        df = pd.DataFrame({"biomass_X": [1.0, 5.0, 10.0]})
        result = self.fitness.biomass_maximization(df)
        # iloc[-2] value = 5.0  →  1/5.0
        self.assertAlmostEqual(result, 1.0 / 5.0)

    def test_execute_default_fitness(self):
        df = pd.DataFrame({"met_A": [1.0, 3.0, 9.0]})
        result = self.fitness.execute_default_fitness(df, "met_A")
        # default is metabolite_maximization, iloc[-2]=3.0  →  1/3.0
        self.assertAlmostEqual(result, 1.0 / 3.0)

    def test_execute_fitness_function_by_name(self):
        df = pd.DataFrame({"met_B": [2.0, 8.0, 32.0]})
        result = self.fitness.execute_fitness_function("metabolite_maximization", df, "met_B")
        # iloc[-2] = 8.0  →  1/8.0
        self.assertAlmostEqual(result, 1.0 / 8.0)


if __name__ == "__main__":
    unittest.main()
