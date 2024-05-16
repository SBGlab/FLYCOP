# Class Factory of Simulators like COMETSSimualtor and dFBASimulator
class SimulatorFactory:
    def __init__(self):
        pass
    @staticmethod
    def createSimulator(simulator_type):
        if simulator_type == "dFBA":
            from dFBASimulator import dFBASimulator
            return dFBASimulator()
        elif simulator_type == "COMETS":
            from COMETSSimulator import COMETSSimulator
            return COMETSSimulator()
        elif simulator_type == "surfin_FBA":
            from SurfinFBASimulator import SurfinFBASimulator
            return SurfinFBASimulator()
        elif simulator_type == "MMODES":
            from MMODESSimulator import MMODESSimulator
            return MMODESSimulator()
        else:
            raise ValueError("Simulator type not supported")