from mpi4py import MPI
import logging
from os.path import abspath
#from MPIFileYHandler.MPIFileYHandler import MPIFileHandler


import signal
import numpy as np
import os
import sys
import copy
from simtk.unit import *
from simtk import unit
from simtk.openmm import app
import simtk.openmm as mm

#from openmmtools import testsystems


class mySimulation(mm.app.simulation.Simulation):
    """
        This class is for a single simulation, it can provide units of replicas in REMD.

    """

    def __init__(self,topology,system,integrator,dcdfrep=5000,repNamePrefix="run",simulationSeg=1,outputDir=None,initReporter=True,drude=True,dcdstepsize=20000, **kwargs):
#       super(mySimulation,self).__init__(self,topology,system,integrator,platform=None,platformProperties=None,state=None):
        if "platform" in kwargs.keys():
            self.platform = kwargs.pop("platform")
        if self.platform.getName() == "CPU":
            self.platformProperties = dict()
        elif "platformProperties" in kwargs.keys():
            self.platformProperties = kwargs.pop("platformProperties")
        else:
            self.platformProperties = dict(CudaPrecision='mixed')
#       platformProperties = kwargs.pop("platformProperties")
#       state = kwargs.pop("state")
        self.drude = drude
        super().__init__(topology,system,integrator,platform=self.platform,platformProperties=self.platformProperties,state=None)
        self.integrator  = integrator
        self.alltemp = kwargs.pop("temps")
        self.context.setPositions(kwargs.pop("positions"))
        self.kineticEnergy = None
        self.potentialEnergy = self._calcE()
        self.langevin = self._isLangevin()
        self.temp = self._getTemp()
        self.temps = self.alltemp
        self.scaleVel = True
        # state of the simulation
        self.dcdfrep  = dcdfrep
        self.restartFrep = dcdfrep*100
        self.repNamePrefix = repNamePrefix
        self.simulationSeg = simulationSeg
        self.outputDir = os.getcwd() if outputDir == None else outputDir
        self.myreporters = []
        # for the restart file, this is for the restart and checkpoint
        self.restart = False
        self.restart_path = None  # this specifies the output path
        if initReporter == True:
            self._reporter_init()
        self._isLangevin()
        self.dcdstepsize = dcdstepsize

    def system_info(self,getPos=False,getVel=False,getEne=True):
        state = self.context.getState(getPositions=getPos,
                                      getVelocities=getVel,
                                      getEnergy=getEne
                                      )
        self.kineticEnergy = state.getKineticEnergy()
        self.potentialEnergy = state.getPotentialEnergy()
        if getPos == True and getVel != True:
            return state.getPositions()
        if getPos == True and getVel == True:
            return(state.getPositions(),state.getVelocities())
        return self.potentialEnergy


#   def reduced_energy(self):
#       """
#       """
#       beta = 1.0 / (unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA)
#       assert self.potentialEnergy != None
#       reduced_potential = self.potentialEnergy * unit.AVOGADRO_CONSTANT_NA
#       # not consider pressure for now
#       betas =  beta * np.array(self.temps)
#       return betas * reduced_potential

    def reduced_energy(self):
        """
        calculate the reduced potential for each temperature
        """
#       beta = 1.0 / (unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA)
        betas = [1.0 / (unit.BOLTZMANN_CONSTANT_kB * temperature) for temperature in self.alltemp]
#       reduced_potential = potential_energy / unit.AVOGADRO_CONSTANT_NA
        assert self.potentialEnergy != None
        reduced_potential = self.potentialEnergy / unit.AVOGADRO_CONSTANT_NA
        # not consider pressure for now
#        betas =  beta * np.array(temps)
        return [beta * reduced_potential for beta in betas]




    def run(self,steps):
        """
        run simulation
        """
        self.reduced_energy()
        self.step(steps)
        state = self.context.getState(getPositions=True,getVelocities=True)
#       self.myreporter[dcd].writeModel(state.getPositions(), periodicBoxVectors=state.getPeriodicBoxVectors())
#       state = self.context.getState( getPositions=True, getVelocities=True )
        if self.restart:
            if self.restart_path != None:
                assert self.restart_path.endswith(".rst")
                with open(self.restart_path, 'w') as f:
                    f.write(mm.XmlSerializer.serialize(state))
            with open(os.path.join(self.outputDir,"{jobname}_{replicaID}_{seg}.rst".format(jobname=self.repNamePrefix,replicaID=self.replicaID,seg=self.simulationSeg)), 'w') as f:
                f.write(mm.XmlSerializer.serialize(state))
        self.potentialEnergy =  self._calcE()

    def change_temperature(self,targetTemperature):
        '''
        To change temperature, coordinates are not modified, one way is to rescale the temperature, this is for Andersen etc method, (seems to have a bug, find time to solve it)
        For the langevin method, only using  the setTemperature is ok
        '''
        if not self.langevin:
            systemInContext   = copy.deepcopy(self.context.getSystem())
            for forcePosition, tmp in enumerate(systemInContext.getForces()):
                if  'Thermostat' in tmp.__class__.__name__:
                    break
            systemInContext.removeForce(forcePosition)
            newThermostat =  mm.AndersenThermostat(targetTemperature*kelvin, 1.0/picosecond)
            systemInContext.addForce(newThermostat)
            newIntegrator = mm.VerletIntegrator(0.002*picoseconds)
            state = self.context.getState(getPositions=True, getVelocities=True)
            PositionsInContext = state.getPositions()
            VelocitiesInContext =state.getVelocities()
            if self.scaleVel:
                scale = np.sqrt(targetTemperature/self.temp._value)
                VelocitiesInContext = VelocitiesInContext * scale
            self.context = mm.Context(systemInContext, newIntegrator)
            self.context.setPositions(PositionsInContext)
            self.context.setVelocities(VelocitiesInContext)
            self.temp = self._getTemp()
            self.integrator = newIntegrator
        else:
            # langevin
#           state = self.context.getState(getPositions=True, getVelocities=True)
#           PositionsInContext = state.getPositions()
#           VelocitiesInContext = state.getVelocities()
#           if self.scaleVel:
#               scale = np.sqrt(targetTemperature/self.temp._value)
#               VelocitiesInContext = VelocitiesInContext * scale
            
            self.integrator.setTemperature(targetTemperature)
#           self.context.setVelocities(VelocitiesInContext)
            self.temp = self._getTemp()
            

    def _getTemp(self):
        print(self.integrator)
        if self.drude:
            te =  self.integrator.getTemperature()
#           te =  self.integrator.getTemperature
            return te
        elif not self.langevin:
            for tmp in self.context.getSystem().getForces():
                if  'Thermostat' in tmp.__class__.__name__:
                    break
            return tmp.getDefaultTemperature()
        else:
            return self.integrator.getTemperature()

    def _calcE(self):
        return self.context.getState(getEnergy=True).getPotentialEnergy()

    def _reporter_init(self):
        dcd_file = os.path.join(self.outputDir,"{jobname}_{seg}.dcd".format(jobname=self.repNamePrefix,seg=self.simulationSeg))
        out_file = os.path.join(self.outputDir,"{jobname}_{seg}.out".format(jobname=self.repNamePrefix,seg=self.simulationSeg))
#       dcd = app.DCDFile(open(dcd_file,"wb"), self.topology, self.integrator.getStepSize(), 0, self.dcdstepsize)
#       dcd = app.DCDFile(open(dcd_file,"wb"), self.topology, self.integrator.getStepSize(), 0, 20)
##
        self.reporters.append(app.DCDReporter(dcd_file, self.dcdfrep))
        out = app.StateDataReporter(out_file, 1000, step=True, time=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True)
        self.reporters.append(out)
    
    def _isLangevin(self):
        if self.drude:
            self.langevin = True
        if "Langevin" in str(self.integrator.__class__):
            self.langevin = True

    def read_rst(self,rst):
        with open(rst, 'r') as f:
            self.context.setState(mm.XmlSerializer.deserialize(f.read()))

    def read_chk(self,chk):
        with open(chk, 'r') as f:
            self.context.setState(mm.XmlSerializer.deserialize(f.read()))

    def write_rst(self,rst_id):
        state = self.context.getState( getPositions=True, getVelocities=True )
        with open("{jobname}".format(jobname=self.repNamePrefix) + '.rst', 'w') as f:
            f.write(mm.XmlSerializer.serialize(state))

    def write_pdb(self,pdbname):
        crd = self.context.getState(getPositions=True).getPositions()
        mm.app.PDBFile(self.topology,crd,open(pdbname,'w'))


