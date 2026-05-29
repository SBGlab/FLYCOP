# Clase para definir la función objetivo basada en funciones con código de cobra

class  ObjectiveFunction:
    def __init__(self, model, objective_function):
        self.model = model
        self.objective_function = objective_function

    def evaluate(self, parameters):
            # Implementación específica para evaluar la función objetivo
        pass
    def set_objective_function(self, objective_function):
        # Implementación específica para establecer la función objetivo
        # La función objetivo tiene que ser una función de python con código de cobra
        self.objective_function = objective_function
        pass