from skopt import gp_minimize
from .Optimizer import ParameterOptimizer


class SKoptOptimizer(ParameterOptimizer):
    """scikit-optimize backend (stub). Implements the ParameterOptimizer ABC."""

    def __init__(self, parameters=None, configuration_space=None):
        super().__init__()
        self.parameters = parameters or {"n_trials": 100}
        self._param_list = []
        self.objective_function = None
        self._result = None
        if configuration_space is not None:
            self.load_configuration_space(configuration_space)

    def load_configuration_space(self, configuration_space):
        """Store the canonical param list and derive skopt search space."""
        if not isinstance(configuration_space, list):
            raise TypeError("configuration_space must be a list of parameter dicts")
        self._param_list = list(configuration_space)

    def load_scenario(self, scenario):
        """Update optimizer settings (n_trials, etc.)."""
        valid_keys = {"n_trials", "n_initial_points", "acq_func"}
        for key in scenario:
            if key not in valid_keys:
                raise ValueError(f"'{key}' is not a valid scenario key. Valid: {valid_keys}")
        self.parameters.update(scenario)

    def load_objective_function(self, objective_function):
        if not callable(objective_function):
            raise ValueError("objective_function must be callable")
        self.objective_function = objective_function

    def run_optimization(self):
        """Run optimisation using gp_minimize.

        Translates the canonical parameter list to skopt space tuples.
        The user objective is called with a plain list of values (skopt
        convention); override or wrap as needed for full FLYCOP integration.
        """
        if not self._param_list:
            raise ValueError("Configuration space not loaded.")
        if self.objective_function is None:
            raise ValueError("Objective function not loaded.")

        space = []
        for p in self._param_list:
            if p["type"] in ("float", "int"):
                space.append((p["min"], p["max"]))
            elif p["type"] == "categorical":
                space.append(p["values"])

        n_calls = self.parameters.get("n_trials", 100)
        self._result = gp_minimize(self.objective_function, space, n_calls=n_calls)
        return self.get_best_parameters()

    def get_best_parameters(self):
        if self._result is None:
            raise ValueError("No optimization has been run yet.")
        names = [p["name"] for p in self._param_list]
        return dict(zip(names, self._result.x))

    def visualize_results(self):
        if self._result is None:
            raise ValueError("No optimization has been run yet.")
        # skopt plotting not implemented here
        pass