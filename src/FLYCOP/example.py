# Python code to demonstrate the working of FLYCOP2


import os
import sys
import time
import json
import numpy as np

from . import FLYCOP2 as flycop
from . import ParameterOptimizer
from . import SMAC3Optimizer

def main():
    # Create a FLYCOP2 object
    flycop2 = FLYCOP2()

    # Load the consortia
    flycop2.consortia = flycop.Consortia()
    flycop2.consortia.load_from_file('consortia.json')

    # Load the parameters to optimize
    flycop2.parameters = flycop.Parameters()
    flycop2.parameters.load_from_file('parameters.json')

    # Load the objective function
    flycop2.objective_function = flycop.ObjectiveFunction()
    flycop2.objective_function.load_from_file('objective_function.json')

    # Load the parameter optimizer
    flycop2.optimizator = flycop.SMAC3Optimizer(flycop2.parameters)

    # Run the optimization
    flycop2.run_optimization()

    # Save the results
    flycop2.save_results()

    # Visualize the results
    flycop2.visualize_results()