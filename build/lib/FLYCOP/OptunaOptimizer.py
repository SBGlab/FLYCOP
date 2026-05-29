"""
OptunaOptimizer
===============
Optuna-based implementation of ParameterOptimizer.

Uses the same **canonical parameter list** format as SMAC3Optimizer so
callers do not need to change anything when switching between backends:

    [
        {'type': 'float',       'name': 'glucose_e',    'min': 0.0, 'max': 20.0},
        {'type': 'int',         'name': 'biomass_ecoli','min': 1,   'max': 10},
        {'type': 'categorical', 'name': 'strain',       'values': ['wt', 'mut']},
    ]

The user objective function is called with a plain ``dict`` ``{name: value}``.
Because ``dict.items()`` and ``ConfigSpace.Configuration.items()`` behave the
same way, ``FLYCOP_scenario.optimize`` works unchanged with both backends.

Scenario keys accepted by ``load_scenario``:
    - ``n_trials``   (int,  default 100)
    - ``direction``  (str,  'minimize' | 'maximize', default 'minimize')
    - ``sampler``    (optuna.samplers.BaseSampler | None)
    - ``pruner``     (optuna.pruners.BasePruner  | None)
"""
import optuna

from .Optimizer import ParameterOptimizer

# Silence Optuna's verbose INFO logging by default.
optuna.logging.set_verbosity(optuna.logging.WARNING)

_VALID_SCENARIO_KEYS = {"n_trials", "direction", "sampler", "pruner"}
_VALID_PARAM_TYPES = {"float", "int", "categorical"}


class OptunaOptimizer(ParameterOptimizer):
    """Optuna-backed hyperparameter optimizer following the FLYCOP contract."""

    def __init__(
        self,
        configuration_space=None,
        objective_function=None,
        parameters=None,
    ):
        super().__init__()
        self._param_list = []
        self.study = None
        self.objective_function = None
        self.parameters = {"n_trials": 100, "direction": "minimize"}
        if parameters:
            self.parameters.update(parameters)

        if configuration_space is not None:
            self.load_configuration_space(configuration_space)

        if objective_function is not None:
            self.load_objective_function(objective_function)

    # ------------------------------------------------------------------
    # ParameterOptimizer interface
    # ------------------------------------------------------------------

    def load_configuration_space(self, configuration_space):
        """Load parameters from the canonical list-of-dicts format.

        Each entry must have at minimum ``type`` and ``name``.
        Float/int entries must also have ``min`` and ``max``.
        Categorical entries must have ``values``.
        """
        if not isinstance(configuration_space, list):
            raise TypeError(
                "configuration_space must be a list of parameter dicts"
            )
        for entry in configuration_space:
            if "type" not in entry or "name" not in entry:
                raise ValueError(
                    f"Each parameter dict must contain 'type' and 'name'. Got: {entry}"
                )
            if entry["type"] not in _VALID_PARAM_TYPES:
                raise ValueError(
                    f"Unknown parameter type '{entry['type']}'. "
                    f"Valid types: {_VALID_PARAM_TYPES}"
                )
            if entry["type"] in ("float", "int"):
                if "min" not in entry or "max" not in entry:
                    raise ValueError(
                        f"Float/int parameter '{entry['name']}' must have 'min' and 'max'"
                    )
            if entry["type"] == "categorical":
                if "values" not in entry or not entry["values"]:
                    raise ValueError(
                        f"Categorical parameter '{entry['name']}' must have a non-empty 'values' list"
                    )
        self._param_list = list(configuration_space)

    def load_scenario(self, scenario):
        """Update optimizer settings from a dict.

        Valid keys: ``n_trials``, ``direction``, ``sampler``, ``pruner``.
        """
        for key in scenario:
            if key not in _VALID_SCENARIO_KEYS:
                raise ValueError(
                    f"'{key}' is not a valid scenario key. "
                    f"Valid keys: {_VALID_SCENARIO_KEYS}"
                )
        self.parameters.update(scenario)

    def load_objective_function(self, objective_function):
        """Store the user objective. Must be callable."""
        if not callable(objective_function):
            raise ValueError("objective_function must be callable")
        self.objective_function = objective_function

    def run_optimization(self):
        """Create an Optuna study and run the optimization.

        Returns
        -------
        dict
            Best parameter values found (``{name: value}``).
        """
        if not self._param_list:
            raise ValueError("Configuration space not loaded. Call load_configuration_space first.")
        if self.objective_function is None:
            raise ValueError("Objective function not loaded. Call load_objective_function first.")

        direction = self.parameters.get("direction", "minimize")
        n_trials = self.parameters.get("n_trials", 100)
        sampler = self.parameters.get("sampler", None)
        pruner = self.parameters.get("pruner", None)

        self.study = optuna.create_study(
            direction=direction,
            sampler=sampler,
            pruner=pruner,
        )
        self.study.optimize(self._build_optuna_objective(), n_trials=n_trials)
        return self.get_best_parameters()

    def get_best_parameters(self):
        """Return the best parameters found as a plain dict.

        Raises ``ValueError`` if no optimization has been run yet.
        """
        if self.study is None:
            raise ValueError("No optimization has been run yet. Call run_optimization first.")
        return self.study.best_params

    def visualize_results(self):
        """Plot optimization history and parameter importances.

        Requires ``plotly`` to be installed (``pip install plotly``).
        """
        if self.study is None:
            raise ValueError("No optimization has been run yet. Call run_optimization first.")
        try:
            from optuna.visualization import (
                plot_optimization_history,
                plot_param_importances,
            )
            plot_optimization_history(self.study).show()
            if len(self._param_list) > 1:
                plot_param_importances(self.study).show()
        except ImportError:
            print(
                "Install plotly for Optuna visualizations: pip install plotly"
            )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_optuna_objective(self):
        """Return a closure that Optuna calls on each trial.

        The closure translates the ``Trial`` API into a plain dict and
        forwards it to the user's objective function, keeping the FLYCOP
        interface backend-agnostic.
        """
        param_list = self._param_list
        user_objective = self.objective_function

        def _optuna_objective(trial):
            params = {}
            for p in param_list:
                if p["type"] == "float":
                    params[p["name"]] = trial.suggest_float(
                        p["name"], p["min"], p["max"]
                    )
                elif p["type"] == "int":
                    params[p["name"]] = trial.suggest_int(
                        p["name"], p["min"], p["max"]
                    )
                elif p["type"] == "categorical":
                    params[p["name"]] = trial.suggest_categorical(
                        p["name"], p["values"]
                    )
            return user_objective(params)

        return _optuna_objective
