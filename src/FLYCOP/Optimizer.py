# Abstract base class for consortium parameter optimizers.
# Concrete backends: SMAC3Optimizer, OptunaOptimizer, SKoptOptimizer.
# All share the same canonical parameter list format:
#   [{'type': 'float'|'int'|'categorical', 'name': ..., 'min'/'max' or 'values': ...}, ...]
from abc import ABC, abstractmethod


class ParameterOptimizer(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def load_configuration_space(self, configuration_space):
        """Load parameters from the canonical list-of-dicts format."""
        pass

    @abstractmethod
    def load_scenario(self, scenario):
        """Configure optimizer-specific settings (n_trials, direction, etc.)."""
        pass

    @abstractmethod
    def load_objective_function(self, objective_function):
        """Store the callable objective function to minimise/maximise."""
        pass

    @abstractmethod
    def run_optimization(self):
        """Execute the search and return the best parameters found."""
        pass

    @abstractmethod
    def get_best_parameters(self):
        """Return the best parameters found as a plain dict {name: value}."""
        pass

    @abstractmethod
    def visualize_results(self):
        """Visualize optimization results (implementation-specific)."""
        pass