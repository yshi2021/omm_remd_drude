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

prmtop_filename = "/home/shiyi/yank_CHARMM_druce_delicate/amber/sams-explicit/complex.prmtop"
crd_filename = "/home/shiyi/yank_CHARMM_druce_delicate/amber/sams-explicit/complex.inpcrd"

prmtop = app.AmberPrmtopFile(prmtop_filename)
#system = prmtop.createSystem(implicitSolvent=None, constraints=constraints, nonbondedCutoff=None, hydrogenMass=hydrogenMass)
system = prmtop.createSystem(implicitSolvent=None, constraints=ff.HBonds, nonbondedCutoff=None)
distance_unit = unit.nanometer
system.setDefaultPeriodicBoxVectors(*(np.identity(3) * distance_unit * 4.55))

testsystem = system

topology = prmtop.topology

# Read positions.
inpcrd = app.AmberInpcrdFile(crd_filename)
positions = inpcrd.getPositions(asNumpy=True)

#system, positions = system, positions




n_replicas = 40  # Number of temperature replicas.
T_min = 298.0 * unit.kelvin  # Minimum temperature.
T_max = 400.0 * unit.kelvin  # Maximum temperature.
reference_state = states.ThermodynamicState(system=testsystem, temperature=T_min)


move = mcmc.GHMCMove(timestep=2.0*unit.femtoseconds, n_steps=2000)
simulation = ParallelTemperingSampler(replica_mixing_scheme='swap-neighbors',mcmc_moves=move, number_of_iterations=20)


#storage_path = tempfile.NamedTemporaryFile(delete=False).name + '.nc'
storage_path = "run1_amber_xs.nc"
reporter = MultiStateReporter(storage_path, checkpoint_interval=1)
simulation.create(reference_state,
                  states.SamplerState(positions),
                  reporter, min_temperature=T_min,
                  max_temperature=T_max, n_temperatures=n_replicas)
#simulation.minimize()
#simulation.equilibrate(10)
#simulation.run()



