import os
import os.path
import numpy as np
import numpy.random
import itertools
import inspect

import scipy
import scipy.special
import scipy.integrate

from simtk import openmm
from simtk import unit
from simtk.openmm import app

from simtk import openmm as mm
from simtk.openmm import app
from simtk.unit import *
from sys import stdout, exit, stderr
from subprocess import call
#from . import drudeutility.CharmmDrudePsfFile as CharmmPsfFile

#import openmmtools as mmtools

#from openmmtools.mcmc import  BaseIntegratorMove
#from openmmtools.utils import *
#from openmmtools.integrators import PrettyPrintableIntegrator
from simtk.openmm.app import forcefield as ff

print(os.getcwd())
#from . drudeutility import *
#from .drudeutility import *
#from . import drudeutility

pi = np.pi



def test_system_AlanineDipeptideImplicit(temp):
    """
        This is just a quick way to generate a system, its positions.
        It is used for testing
        
    """
    from openmmtools import testsystems
    testsystem = testsystems.AlanineDipeptideImplicit()
    themostat = mm.AndersenThermostat(temp*kelvin, 1.0/picosecond)
    testsystem.system.addForce(themostat)
    return testsystem.system,testsystem.positions,testsystem.topology

class TestSystem(object):
    """
    part copy from openmm tools
    """
    def __init__(self, **kwargs):
        self._system = openmm.System()
        self._positions = unit.Quantity(np.zeros([0, 3], np.float), unit.nanometers)
        self._topology = app.Topology()
        self._mdtraj_topology = None
        self._mdtraj_pdb = None # if start from pdb
        self._psf = None
        self._crd = None
        self._param = None
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
        if self._mdtraj_topology is None:
            self._mdtraj_topology = md.loadpdb("self._mdtraj_pdb")
        return self._mdtraj_topology
    def _initialize(self):
        """
        parse the parameters
        """
        pass

    def buildfromPsf(self):
        self._initialize(self)
        psf = app.CharmmPsfFile(self._psf)
        crd = app.CharmmCrdFile(self._crd)
        pdb =       app.PDBFile(self._pdb)
        assert self._pre_size != None
        # if _pre_size is nanometer in units
        boxsize = self._pre_size
        psf.setBox(boxsize*nanometer,boxsize*nanometer,boxsize*nanometer)
        params = app.CharmmParameterSet(*self.params)
        self._system = self.psf.createSystem(self.params,  nonbondedMethod=ff.PME, nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0001, constraints=ff.HBonds)
        self._positions = crd.positions
    
    def buildfromPsfDrude(self):
        self._initialize(self)
        psf = app.CharmmPsfFile(self._psf)
        crd = app.CharmmCrdFile(self._crd)
        pdb =       app.PDBFile(self._pdb)
        assert self._pre_size != None
        # if _pre_size is nanometer in units
        boxsize = self._pre_size
        psf.setBox(boxsize*nanometer,boxsize*nanometer,boxsize*nanometer)
        params = app.CharmmParameterSet(*self.params)
        self._system = self.psf.createSystem(self.params,  nonbondedMethod=ff.PME, nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0001, constraints=ff.HBonds)
        self._positions = crd.positions
    def _nbfix():
        """
        for modify Drude
        """
        

# copy from openmmtools
    def serialize(self):
        """Return the System and positions in serialized XML form.

        Returns
        -------

        system_xml : str
            Serialized XML form of System object.

        state_xml : str
            Serialized XML form of State object containing particle positions.

        """

        from simtk.openmm import XmlSerializer

        # Serialize System.
        system_xml = XmlSerializer.serialize(self._system)

        # Serialize positions via State.
        if self._system.getNumParticles() == 0:
            # Cannot serialize the State of a system with no particles.
            state_xml = None
        else:
            platform = openmm.Platform.getPlatformByName('Reference')
            integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
            context = openmm.Context(self._system, integrator, platform)
            context.setPositions(self._positions)
            state = context.getState(getPositions=True)
            del context, integrator
            state_xml = XmlSerializer.serialize(state)

        return (system_xml, state_xml)




class drudeSystem(TestSystem):
    def __init__(self,psf,crd,params):
        """
        This overwrite the TestSystem init function
        """
        self._psf = app.CharmmPsfFile(psf)
        self.size = 54
        size = self.size
        boxsize = (float(size) + float(size) * 0.01) / 10 # PBC size in nanometer
        self._psf.setBox(boxsize*nanometer,boxsize*nanometer,boxsize*nanometer)
        self._crd = app.CharmmCrdFile(crd)
        self.params = app.CharmmParameterSet(*params)
        self._system = self._psf.createSystem(self.params,  nonbondedMethod=ff.PME, nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0001, constraints=ff.HBonds)
        self.platform = mm.Platform.getPlatformByName('CPU')
        self._positions = self._crd.positions
        self._topology = app.Topology()
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



class charmm36System(TestSystem):
    def __init__(self,psf,crd,params):
        """
        This overwrite the TestSystem init function
        """
        self._psf = app.CharmmPsfFile(psf)
        self.size = 53
        size = self.size
        boxsize = (float(size) + float(size) * 0.01) / 10 # PBC size in nanometer
        self._psf.setBox(boxsize*nanometer,boxsize*nanometer,boxsize*nanometer)
        self._crd = app.CharmmCrdFile(crd)
        self.params = app.CharmmParameterSet(*params)
        self._system = self._psf.createSystem(self.params,  nonbondedMethod=ff.PME, nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0001, constraints=ff.HBonds)
        self.platform = mm.Platform.getPlatformByName('CPU')
        self._positions = self._crd.positions
        self._topology = app.Topology()
        self._mdtraj_topology = None
    def read_rst(self,rst):
        pass
    
    def calculate_size_from_crd(self):
        pass



class aaqaa(charmm36System):
    """
    params = ["/home/shiyi/yank_CHARMM_druce_delicate/aaqaa/toppar/par_all36m_prot.prm","/home/shiyi/yank_CHARMM_druce_delicate/aaqaa/toppar/toppar_water_ions.str"]
    """
    def _get_aaqaa_sys(self):
        from . data.get_data import get_drude_params,get_36m_params,get_aaqaa,get_villin
        psf,crd,pdb,rst = get_aaqaa(drude=False)
        params_tmp = get_36m_params()
        d = os.path.dirname(params_tmp[0])
        params = [os.path.join(d,x) for x in ["par_all36m_prot.prm","toppar_water_ions.str"]]
        return psf[0], crd[0], pdb[0], rst[0], params   # params become one list
    def __init__(self):
        psf,crd,pdb,rst,params = self._get_aaqaa_sys()
        super().__init__(psf,crd,params)
        self._topology = self._psf.topology




def get_aaqaa_drude_sys(temps):
    pass

def get_villin_sys(temps):
    pass




# I have failed for now, use this file as a scripts to directly run the test, check it out later for further reason
if __name__ == "__main__":
    # get all data files
    filepath = os.path.abspath(__file__)
    datapath = os.path.dirname(filepath,"data")
    allfile = []
    for dirpath,_,filenames in os.walk(datapath):
        for f in filenames:
            allfile.append(os.path.join(dirpath,f))
    import re
    
    # get all data files end

#   villin = villin(psf,crd,params)
#   villin = villinDrude(psf,crd,params)
#   psf = "/home/shiyi/my/lastfile/inputfile/step2_drude_2019f2.psf"
#   psf = "/home/shiyi/my/test60/step2_drude_2019f2.psf"
#   #psf2 = CharmmDrudePsfFile("vill_wild_durde.psf")
#   crd = "/home/shiyi/my/test60/step2_drude_2019f2.crd"
#   params = ["/home/shiyi/my/lastfile/inputfile/toppar_drude_master_protein_2019f.str"]
#   params = ["/home/shiyi/my/lastfile/inputfile/toppar_drude_master_protein_2019f.str"]

