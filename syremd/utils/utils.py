import os
import subprocess



def get_intertrator():
    """
        generate a integrator
        
    """
    import simtk.openmm as mm
    from simtk.unit import picoseconds
    integrator = mm.VerletIntegrator(0.002*picoseconds)
    return integrator

def get_AndersenThermostat(temp):
    """
        This is just a quick way to generate a system, its positions.
        It is used for testing
        
    """
    from openmmtools import testsystems
    import simtk.openmm as mm
    themostat = mm.AndersenThermostat(temp*kelvin, 1.0/picosecond)
    return themostat




def get_intertrator_langevin():
    """
        generate a integrator
        
    """
    import simtk.openmm as mm
    from simtk.unit import picoseconds
    integrator = mm.VerletIntegrator(0.002*picoseconds)
    return integrator

def get_context(system,integrator,positions):
    """
        generate a context
        
    """
    context = mm.Context(system, integrator)
    context.setPositions(positions)
    return context

def temp_range(tmin,tmax,num):
    """
        generate a range of temperature for REMD
        
    """
    import numpy as np
    trange = [tmin * np.exp(np.log(1.0*tmax/tmin)*(1.0*i/(num-1))) for i in range(0,num)]
    return trange


def run_linux_process(command):
    p=subprocess.Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    p.wait()
    output, err=p.communicate()
    return output, err

def run_linux_process2(command):
    results = os.popen(command).read()
    return results



