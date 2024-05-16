

class OptimizerFactory:
    def __init__(self):
        pass
    @staticmethod
    def createOptimizer(optimizer_name, parameters=None,configuration_space=None):
        if optimizer_name == 'SMAC3':
            from SMAC3Optimizer import SMAC3Optimizer
            return SMAC3Optimizer(parameters=parameters)
        elif optimizer_name == 'SKopt':
            from SKoptOptimizer import SKoptOptimizer
            return SKoptOptimizer(parameters=parameters,configuration_space=configuration_space)
        else:
            raise ValueError('Invalid optimizer name')
