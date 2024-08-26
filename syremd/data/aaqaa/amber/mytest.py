#!/home/shiyi/anaconda3/bin/python
#import math
#from simtk import unit
#from openmmtools import testsystems, states, mcmc
#testsystem = testsystems.AlanineDipeptideImplicit()
#import os
#import tempfile
#
#
#n_replicas = 40  # Number of temperature replicas.
#T_min = 298.0 * unit.kelvin  # Minimum temperature.
#T_max = 600.0 * unit.kelvin  # Maximum temperature.
#temperatures = [T_min + (T_max - T_min) * (math.exp(float(i) / float(n_replicas-1)) - 1.0) / (math.e - 1.0)
#                for i in range(n_replicas)]
#thermodynamic_states = [states.ThermodynamicState(system=testsystem.system, temperature=T)
#                        for T in temperatures]
#
#
#move = mcmc.GHMCMove(timestep=2.0*unit.femtoseconds, n_steps=50)
#simulation = ReplicaExchangeSampler(mcmc_moves=move, number_of_iterations=2)
#
#te simulation with its storage file (in a temporary directory) and run.
#
#storage_path = tempfile.NamedTemporaryFile(delete=False).name + '.nc'
#reporter = multistate.MultiStateReporter(storage_path, checkpoint_interval=1)
#simulation.create(thermodynamic_states=thermodynamic_states,
#                  sampler_states=states.SamplerState(testsystem.positions),
#                  storage=reporter)
#simulation.run()  # This runs for a maximum of 2 iterations.
#simulation.iteration
#
#simulation.run(n_iterations=1)
#simulation.iteration


import numpy as np
from simtk import unit
from openmmtools import testsystems, states, mcmc
import tempfile

from simtk import openmm
from simtk import unit
from simtk.openmm import app
from simtk.openmm.app import (forcefield as ff, Topology, element) 

from yank.multistate import ParallelTemperingSampler, ParallelTemperingAnalyzer, MultiStateReporter, MultiStateSampler


#testsystem = testsystems.AlanineDipeptideImplicit()
import os

prmtop_filename = "/home/shiyi/yank_CHARMM_druce_delicate/amber_aapaa/aaqaa.prmtop"
crd_filename = "/home/shiyi/yank_CHARMM_druce_delicate/amber_aapaa/aaqaa.inpcrd"

prmtop = app.AmberPrmtopFile(prmtop_filename)
prmtop.topology.setPeriodicBoxVectors(vectors=np.diag([4.5,4.5,4.5]))
#system = prmtop.createSystem(implicitSolvent=None, constraints=constraints, nonbondedCutoff=None, hydrogenMass=hydrogenMass)
# system = psf.createSystem(self.params,  nonbondedMethod=ff.PME, nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer, ewaldErrorTolerance = 0.0001, constraints=ff.HBonds)
system = prmtop.createSystem(implicitSolvent=None, constraints=ff.HBonds, nonbondedCutoff=None)
system = prmtop.createSystem(implicitSolvent=None, constraints=ff.HBonds,nonbondedMethod=ff.PME, nonbondedCutoff=None)
distance_unit = unit.nanometer
#system.setDefaultPeriodicBoxVectors(*(np.identity(3) * distance_unit * 4.58018620))

testsystem = system

topology = prmtop.topology

# Read positions.
inpcrd = app.AmberInpcrdFile(crd_filename)
positions = inpcrd.getPositions(asNumpy=True)

#system, positions = system, positions




n_replicas = 40  # Number of temperature replicas.
T_min = 280.0 * unit.kelvin  # Minimum temperature.
T_max = 400.0 * unit.kelvin  # Maximum temperature.

def temp_range(tmin,tmax,num):
        import numpy as np
        trange = [tmin * np.exp(np.log(1.0*tmax/tmin)*(1.0*i/(num-1))) for i in range(0,num)]
        return trange



temperatures = [float(x) for x in "300.00 302.61 305.24 307.89 310.57 313.27 316.00 318.74 321.52 324.31 327.13 329.98 332.85 335.74 338.66 341.60 344.57 347.57 350.59 353.64 356.72 359.82 362.95 366.10 369.29 372.50 375.74 379.01 382.30 385.63 388.98 392.36 395.77 399.21 402.69 406.19 409.72 413.28 416.88 420.50 424.16 427.85 431.57 435.32 439.10 442.92 446.77 450.66 454.58 458.53 462.52 466.54 470.60 474.69 478.82 482.98 487.18 491.42 495.69 500.00".split()] * unit.kelvin
temperatures = [float(x) for x in "300.00 302.61 305.24 307.89 310.57 313.27 316.00 318.74 321.52 324.31 327.13 329.98 332.85 335.74 338.66 341.60 344.57 347.57 350.59 353.64 356.72 359.82 362.95 366.10 369.29 372.50 375.74 379.01 382.30 385.63 388.98 392.36 395.77 399.21 402.69 406.19 409.72 413.28 416.88 420.50".split()] * unit.kelvin
temperatures = temp_range(280,400,40) * unit.kelvin
                                                          


reference_state = states.ThermodynamicState(system=testsystem, temperature=T_min)


move = mcmc.GHMCMove(timestep=2.0*unit.femtoseconds, n_steps=2000)
simulation = ParallelTemperingSampler(replica_mixing_scheme='swap-neighbors',mcmc_moves=move, number_of_iterations=2)


#storage_path = tempfile.NamedTemporaryFile(delete=False).name + '.nc'
storage_path = "run1_amber_4.nc"
reporter = MultiStateReporter(storage_path, checkpoint_interval=1)
#simulation.create(reference_state,
#                  states.SamplerState(positions),
#                  reporter, min_temperature=T_min,
#                  max_temperature=T_max, n_temperatures=n_replicas)
#
simulation.create(reference_state,states.SamplerState(positions),reporter,n_temperatures=40,temperatures=temperatures)


simulation.minimize()
simulation.equilibrate(10)
simulation.run()



