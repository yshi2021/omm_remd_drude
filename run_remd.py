#!/home/shiyi/anaconda3/bin/python
from syremd.mysimulation import mySimulation
from syremd.utils.utils import get_intertrator,temp_range
from syremd.mylogger import myformatter
from syremd.data import get_testsystem,get_data
import simtk.openmm as mm
from mpi4py import MPI
import logging
import os

from simtk.unit import *
from simtk import unit


from collections import OrderedDict

class SlicableOrderedDict(OrderedDict):
    def __getitem__(self, k):
        if not isinstance(k, slice):
            return OrderedDict.__getitem__(self, k)
        x = SlicableOrderedDict()
        for idx, key in enumerate(self.keys()):
            if k.start <= idx < k.stop:
                x[key] = self[key]
        return x




#####
platform = mm.Platform.getPlatformByName('CUDA')


######
# setting parameters
num_replica = 45
tmin = 298
tmax = 450
temps = temp_range(tmin,tmax,num_replica)



######
# input files

psf = "aaqaa/step3_charmm2omm.psf"
crd = "aaqaa/step3_charmm2omm.crd"
rst = "aaqaa/replica_20_4.rst"
params = ["aaqaa/toppar_drude_master_protein_2019g.str"]

prop = dict(CudaPrecision='mixed')


# set up and initialize the system
from syremd.myREMD import myREMD

task_groups_infos = SlicableOrderedDict() 
for x,t in enumerate(temps):
    task_groups_infos["replica_"+str(x)] = t

task_size = len(temps)
iteration = 1
maxitIteration = 10

temps_unit = [x * kelvin for x in temps]
remd_sim = myREMD(task_groups_infos,task_size,iteration,maxitIteration,temps_unit,crd=crd,psf=psf,size=4.9,params=params,runtype="drude",platform="CUDA",platformProperties=prop)


# setting up MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

remd_sim.node_job_detect()
remd_sim.cycle_per_steps = 2000

# set up the systemd




# read chk or rst file
## if restart, here should specify which replica is which temperature. such as storing the temperature for each replica in a file called sort_replica_temp.dat, 
## for l in  open("sort_replica_temp.dat","r"):
##     a,b = l.split()
##     a = int(a.strip("replica_"))
##     b = float(b) * kelvin
##     temp_units_new[a] = b
## then using sim_replica.change_temperature(temp_units_new[k]), see below

for k,sim_replica in remd_sim.node_task.items():
    # k represent replica id (the replica coordinates do not change, while the temperature can change)
    k = int(k.strip("replica_"))

    name = rst
#   name = chk
    sim_replica.read_rst("{}".format(name))
#   sim_replica.loadCheckpoint(chk)
#   sim_replica.change_temperature(temp_units_new[k])


remd_sim.comm.barrier()



#os.chdir("output")

# setting the dcd segments
for i in range(0,21):
    # each segment is about 10 ns
    iteration = 1 + i * 2500
    maxitIteration = (i+1) * 2500
    # run simulation
    remd_sim.iteration = iteration
    remd_sim.maxitIteration = maxitIteration
    remd_sim.executor()
    # save restart files
    for k,xxx in remd_sim.node_task.items():
        xxx.saveState("{}_{}.rst".format(k,str(i)))
    remd_sim.comm.barrier()





