import os
from simtk import openmm
from simtk import unit
from simtk.openmm import app

from simtk.unit import *
from simtk.openmm.app.charmmpsffile import *
from simtk.openmm.app.charmmpsffile import CharmmPsfFile
from simtk.openmm.app.charmmpsffile import CharmmPsfFile

import numpy as np

#from .. drudeutility import *





class TestSystem(object):
    def __init__(self, **kwargs):
        self._system = openmm.System()
        self._positions = unit.Quantity(np.zeros([0, 3], np.float), unit.nanometers)
        self._topology = app.Topology()
        self._mdtraj_topology = None
        return

    @property
    def system(self):
        return self._system
    @system.setter
    def system(self, value):
        self._system = value
    @system.deleter
    def system(self):
        del self._system
    @property
    def positions(self):
        """The simtk.unit.Quantity object containing the particle positions, with units compatible with simtk.unit.nanometers."""
        return self._positions

    @positions.setter
    def positions(self, value):
        self._positions = value

    @positions.deleter
    def positions(self):
        del self._positions

    @property
    def topology(self):
        """The simtk.openmm.app.Topology object corresponding to the test system."""
        return self._topology

    @topology.setter
    def topology(self, value):
        self._topology = value
        self._mdtraj_topology = None

    @topology.deleter
    def topology(self):
        del self._topology


    @property
    def mdtraj_topology(self):
        """The mdtraj.Topology object corresponding to the test system (read-only)."""
        import mdtraj as md
        if self._mdtraj_topology is None:
            self._mdtraj_topology = md.Topology.from_openmm(self._topology)
        return self._mdtraj_topology



class villin(TestSystem):
    def __init__(self,psf,crd,params):
        self._psf = CharmmPsfFile(psf)
        self.size = 54
        size = self.size
        boxsize = (float(size) + float(size) * 0.01) / 10 # PBC size in nanometer
        self._psf.setBox(boxsize*nanometer,boxsize*nanometer,boxsize*nanometer)
        self._crd = app.CharmmCrdFile(crd)
        self.params = app.CharmmParameterSet(*params)
        self._system = self._psf.createSystem(self.params,  nonbondedMethod=ff.PME, nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0001, constraints=ff.HBonds)
        self.platform = mm.Platform.getPlatformByName('CPU')
        self._positions = self._crd.positions
        self._topology = self._psf.topology
        self._mdtraj_topology = None
        self._nbfix()

    def _nbfix(self):
        nbforce = [self._system.getForce(i) for i in range(self._system.getNumForces()) if isinstance(self._system.getForce(i), mm.NonbondedForce)][0]
        nbfix = [self._system.getForce(i) for i in range(self._system.getNumForces()) if isinstance(self._system.getForce(i), mm.CustomNonbondedForce)][0]

        nbforce.setNonbondedMethod(mm.NonbondedForce.PME)
        nbforce.setEwaldErrorTolerance(0.0001)
        nbforce.setCutoffDistance(1.2*nanometer)
        nbforce.setUseSwitchingFunction(True)
        nbforce.setSwitchingDistance(1.0*nanometer)

        nbfix.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        nbfix.setCutoffDistance(1.2*nanometer)
        nbfix.setUseSwitchingFunction(True)
        nbfix.setSwitchingDistance(1.0*nanometer)





class helix(TestSystem):
    def __init__(self,psf,crd,params):
        self._psf = CharmmPsfFile(psf)
        self.size = 45
        size = self.size
        boxsize = (float(size) + float(size) * 0.01) / 10 # PBC size in nanometer
        self._psf.setBox(boxsize*nanometer,boxsize*nanometer,boxsize*nanometer)
        self._crd = app.CharmmCrdFile(crd)
        self.params = app.CharmmParameterSet(*params)
        self._system = self._psf.createSystem(self.params,  nonbondedMethod=ff.PME, nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0001, constraints=ff.HBonds)
        self.platform = mm.Platform.getPlatformByName('CPU')
        self._positions = self._crd.positions
        self._topology = self._psf.topology
        self._mdtraj_topology = None
        self._nbfix()

    def _nbfix(self):
        nbforce = [self._system.getForce(i) for i in range(self._system.getNumForces()) if isinstance(self._system.getForce(i), mm.NonbondedForce)][0]
        nbfix = [self._system.getForce(i) for i in range(self._system.getNumForces()) if isinstance(self._system.getForce(i), mm.CustomNonbondedForce)][0]

        nbforce.setNonbondedMethod(mm.NonbondedForce.PME)
        nbforce.setEwaldErrorTolerance(0.0001)
        nbforce.setCutoffDistance(1.2*nanometer)
        nbforce.setUseSwitchingFunction(True)
        nbforce.setSwitchingDistance(1.0*nanometer)

        nbfix.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        nbfix.setCutoffDistance(1.2*nanometer)
        nbfix.setUseSwitchingFunction(True)
        nbfix.setSwitchingDistance(1.0*nanometer)



