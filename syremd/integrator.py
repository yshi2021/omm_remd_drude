import numpy as np
import os
import sys
import copy
from simtk.unit import *
from simtk import unit
from simtk.openmm import app
import simtk.openmm as mm

#from openmmtools import testsystems


class my_integrator(object):
    """
        the class now only have integrator, not other things
    """
    def __init__(self,drude=False,**kwargs):
        """
        This is simplify version of for reducing unnesseary easily debug
        """
        if "integrator" in kwargs.keys():
            self.integrator_keyword = kwargs['integrator'].lower()
        else:
            self.integrator_keyword = None
        if 'temp' in kwargs.keys():
            self.temp = kwargs['temp']
        else:
            self.temp = 300
        self.integrator = None
        if drude == True:
            self.integrator = mm.DrudeLangevinIntegrator(float(self.temp)*kelvin, 5/picosecond, 1*kelvin, 20/picosecond, 0.001*picoseconds)
            self.integrator.setMaxDrudeDistance(0.025)
        elif self.integrator_keyword == "langevin":
            self.integrator = mm.LangevinIntegrator(float(self.temp)*kelvin, 1/picosecond, 0.002*picoseconds)
        elif self.integrator_keyword == "verlet":
            self.integrator = mm.VerletIntegrator(0.002*picoseconds)
        else:
            raise "This is the deubg mod, only langevin, drudeLangevin, andersen, verlet is supported"






class my_integrator2(object):
    def __init__(self,drude=False,**kwargs):
        self._parsearg(self,drude,**kwargs)
        self._initialize(self)
    def _parsearg(self):
        pass
    def _initialize(self):
        pass
    
           
            
            
        
if __name__ == "__main__":
    # just test
    i = my_integrator(temp=300,integrator="langevin")
    print(i)
    print(i.integrator)
    print("OK")
