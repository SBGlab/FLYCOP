"""
pytest conftest.py - Sets up sys.path and mocks all heavy external dependencies
so tests can run without requiring cobra, cometspy, smac, etc. to be installed.
"""
import sys
import os
from unittest.mock import MagicMock

# Add src/ to the Python path so packages are importable as FLYCOP.xxx
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


# ---------------------------------------------------------------------------
# Real stub classes needed for isinstance() checks inside the source code
# ---------------------------------------------------------------------------

class FakeCobraModel:
    """Minimal stub that satisfies isinstance(x, cobra.Model) checks."""
    def __init__(self, model_id="model", name="model"):
        self.id = model_id
        self.name = name
        self.metabolites = []
        self.reactions = []


class FakeHPOFacade:
    """Minimal stub that satisfies isinstance(x, HPOFacade) checks."""
    pass


# ---------------------------------------------------------------------------
# Build mocks for each external namespace
# ---------------------------------------------------------------------------

# cobra
cobra_mock = MagicMock()
cobra_mock.Model = FakeCobraModel
cobra_io_mock = MagicMock()
cobra_io_mock.load_model = MagicMock(return_value=FakeCobraModel())
cobra_flux_mock = MagicMock()
cobra_mock.flux_analysis = cobra_flux_mock
cobra_mock.io = cobra_io_mock

sys.modules.setdefault("cobra", cobra_mock)
sys.modules.setdefault("cobra.io", cobra_io_mock)
sys.modules.setdefault("cobra.flux_analysis", cobra_flux_mock)

# smac
smac_mock = MagicMock()
smac_mock.HyperparameterOptimizationFacade = FakeHPOFacade
smac_mock.RunHistory = MagicMock()
smac_mock.Scenario = MagicMock()
sys.modules.setdefault("smac", smac_mock)

# ConfigSpace (used by FLYCOP.py and SMAC3Optimizer.py at import time)
cs_mock = MagicMock()
cs_mock.Configuration = MagicMock()
cs_mock.ConfigurationSpace = MagicMock()
cs_mock.Float = MagicMock()
cs_mock.Integer = MagicMock()
cs_mock.Categorical = MagicMock()
sys.modules.setdefault("ConfigSpace", cs_mock)

# cometspy
sys.modules.setdefault("cometspy", MagicMock())

# skopt
sys.modules.setdefault("skopt", MagicMock())

# networkx
sys.modules.setdefault("networkx", MagicMock())
sys.modules.setdefault("networkx.drawing", MagicMock())
sys.modules.setdefault("networkx.drawing.nx_agraph", MagicMock())

# matplotlib
sys.modules.setdefault("matplotlib", MagicMock())
sys.modules.setdefault("matplotlib.pyplot", MagicMock())

# metconsin (SurfinFBASimulator)
sys.modules.setdefault("metconsin", MagicMock())

# ---------------------------------------------------------------------------
# optuna – real-enough stubs for unit testing OptunaOptimizer without
# requiring the actual optuna package to be installed.
# ---------------------------------------------------------------------------

class FakeOptunaTrial:
    """Stub for optuna.Trial that records suggestions and returns fixed values."""
    def __init__(self, param_values):
        # param_values: dict {name: value} to return for each suggest call
        self._values = param_values

    def suggest_float(self, name, low, high):
        return float(self._values.get(name, (low + high) / 2.0))

    def suggest_int(self, name, low, high):
        return int(self._values.get(name, (low + high) // 2))

    def suggest_categorical(self, name, choices):
        return self._values.get(name, choices[0])


class FakeOptunaStudy:
    """Stub for optuna.Study that runs n_trials calls synchronously."""
    def __init__(self, direction="minimize", param_values=None):
        self.direction = direction
        self._param_values = param_values or {}
        self.trials = []
        self.best_value = None
        self.best_params = {}

    def optimize(self, objective, n_trials=1, **kwargs):
        for _ in range(n_trials):
            trial = FakeOptunaTrial(self._param_values)
            value = objective(trial)
            self.trials.append({"value": value, "params": self._param_values.copy()})
        # pick best
        if self.direction == "minimize":
            best = min(self.trials, key=lambda t: t["value"])
        else:
            best = max(self.trials, key=lambda t: t["value"])
        self.best_value = best["value"]
        self.best_params = best["params"]


def _fake_create_study(direction="minimize", sampler=None, pruner=None, **kwargs):
    return FakeOptunaStudy(direction=direction)


optuna_mock = MagicMock()
optuna_mock.create_study = _fake_create_study
optuna_mock.logging = MagicMock()
optuna_mock.logging.WARNING = 30
optuna_mock.logging.set_verbosity = MagicMock()
sys.modules.setdefault("optuna", optuna_mock)
sys.modules.setdefault("optuna.visualization", MagicMock())
sys.modules.setdefault("optuna.samplers", MagicMock())
sys.modules.setdefault("optuna.pruners", MagicMock())

# Expose the stubs so individual test files can reuse them
FakeCobraModel = FakeCobraModel
FakeHPOFacade = FakeHPOFacade
FakeOptunaTrial = FakeOptunaTrial
FakeOptunaStudy = FakeOptunaStudy
