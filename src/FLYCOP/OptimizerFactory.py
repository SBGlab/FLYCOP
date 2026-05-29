

class OptimizerFactory:
    def __init__(self):
        pass

    @staticmethod
    def createOptimizer(optimizer_name, parameters=None, configuration_space=None):
        """Create an optimizer by name.

        Supported names: 'SMAC3', 'Optuna', 'SKopt'.

        All optimizers accept the same *canonical* parameter list format via
        ``load_configuration_space``:
            [{'type': 'float'|'int'|'categorical', 'name': ...,
              'min'/'max' (float/int) or 'values' (categorical)}, ...]

        Parameters
        ----------
        optimizer_name : str
            Backend identifier.
        parameters : dict | None
            Keyword arguments forwarded to the optimizer constructor.
        configuration_space : list | None
            Canonical param list to pre-load into the optimizer.
        """
        if optimizer_name == 'SMAC3':
            from .SMAC3Optimizer import SMAC3Optimizer
            return SMAC3Optimizer(parameters=parameters,
                                  configuration_space=configuration_space)
        elif optimizer_name == 'Optuna':
            from .OptunaOptimizer import OptunaOptimizer
            return OptunaOptimizer(parameters=parameters,
                                   configuration_space=configuration_space)
        elif optimizer_name == 'SKopt':
            from .SKoptOptimizer import SKoptOptimizer
            return SKoptOptimizer(parameters=parameters,
                                  configuration_space=configuration_space)
        else:
            raise ValueError(
                f"Unknown optimizer '{optimizer_name}'. "
                "Supported: 'SMAC3', 'Optuna', 'SKopt'."
            )
