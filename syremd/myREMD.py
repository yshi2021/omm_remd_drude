import os
from mpi4py import MPI
import logging
from os.path import abspath
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

class bcolors:
    """
    output color for debug
    """
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
#from MPIFileYHandler.MPIFileYHandler import MPIFileHandler


import inspect
import signal
import numpy as np
#import pandas as pd
import os
import sys
import copy
from simtk.unit import *
from simtk import unit
from simtk.openmm import app
import simtk.openmm as mm




from . mylogger import myformatter
from .mytestsystem import test_system_AlanineDipeptideImplicit, drudeSystem, charmm36System
from .integrator import my_integrator
from .mysimulation import mySimulation


def mylogger_hand(fi,line):
#   with open(fi,'a') as f:
    with open(fi,'a') as f:
        f.write(line)


def mylogger_hand(fi,line):
#   with open(fi,'a') as f:
    with open(fi,'a') as f:
        f.write(line)





def get_mpicomm(logger, disable_mpi=False):
    """
        copy from mpiplus module
        generate a MPI commincator which is mainly used for handling error
        
    """
    from mpi4py import MPI
    mpicomm = MPI.COMM_WORLD
    # Override sys.excepthook to abort MPI on exception.
    def mpi_excepthook(type, value, traceback):
        sys.__excepthook__(type, value, traceback)
        node_name = '{}/{}'.format(MPI.COMM_WORLD.rank+1, MPI.COMM_WORLD.size)
        # logging.exception() automatically print the sys.exc_info(), but here
        # we may want to save the exception traceback of another MPI node so
        # we pass the traceback manually.
        logger.critical('MPI node {} raised an exception and called Abort()! The '
                        'exception traceback follows'.format(node_name), exc_info=value)
        # Flush everything.
        sys.stdout.flush()
        sys.stderr.flush()
        for logger_handler in logger.handlers:
            logger_handler.flush()
        # Abort MPI execution.
        if MPI.COMM_WORLD.size > 1:
            MPI.COMM_WORLD.Abort(1)
    # Use our exception handler.
    sys.excepthook = mpi_excepthook
    def handle_signal(signal, frame):
        if mpicomm.size > 1:
            mpicomm.Abort(1)
    for sig in [signal.SIGINT, signal.SIGTERM, signal.SIGABRT]:
        signal.signal(sig, handle_signal)

    # Cache and return the MPI communicator.

    # Report initialization
    logger.debug("MPI initialized on node {}/{}".format(mpicomm.rank+1, mpicomm.size))

    return mpicomm

class myREMD(object):
    def __init__(self,task_groups_infos,task_size,iteration,maxitIteration,temps,**kwargs):
        '''
            REMD simulation
        '''
        self.handle_kargs(kwargs)
        if len(task_groups_infos) == 0:
            type_sim = kwargs.pop("type_sim")

        if "template" in kwargs.keys():
            self.template = kwargs["template"]
        self.logger = myformatter()
        self.comm = get_mpicomm(self.logger)
        self.rank_size = self.comm.Get_size()
        self.rank = self.comm.Get_rank()
        self.task_size = task_size
        self.temps = temps
        print("temps is {}".format(self.temps))
        self.task_groups = self._initialize(task_groups_infos)
        self.iteration = iteration
        self.maxitIteration = maxitIteration
        self.cycle_per_steps = 2000
        self.results_tmp = None
        self.all_system_energy = None
        self.modify_state = None
        self.exchage_pattern = self._splitList(list(range(self.task_size))) # this corresponding to temperature, not replica id
        self.replica_temp_pair = {}  # this monitor the replica id and its corresponding temperature
        self.recvbuf = None
        self.templog = "temp_traj.log"
        self.to_update = {}


    def __repr__(self):
        """this is facilitate the state monitor of the object , it should be enriched later after """
        fmt = "in node {}\
               the tasks are {}\n".format(self.rank,self.task_groups)
        return(fmt)

    def executor(self):
        """
            worke flow
        """
        for itertime in range(self.iteration ,self.maxitIteration+1):
            # run each replica at first, send results to root rank
            self.check_node_task()
            self.run_one_cycle()
            self.check_rank1_info()
            # attempt exchange,and send back resutls to each rank
            self.root_to_other_node()
            # record exchange process to a file
            for task_k,task_v in self.node_task.items():
                old_temp = task_v._getTemp()._value
                if self.to_update[old_temp] != old_temp:
                    self.logger.info("In iteration {}, The replica {} in rank {} before switch is {}\n".format(self.iteration, task_k,self.rank,task_v._getTemp()))
                    task_v.change_temperature(self.to_update[old_temp])
                    self.logger.info("In iteration {}, The replica {} in rank {} switch temperature from {} to {}\n".format(self.iteration,task_k,self.rank,old_temp,self.to_update[old_temp]))
                    # check a lillte
                    self.logger.info("In iteration {}, Check: The replica {} in rank {} temperature {} to {},indeed is {}\n".format(self.iteration,task_k,self.rank,old_temp,self.to_update[old_temp],task_v._getTemp()))
                    mylogger_hand("hand_log{}.log".format(self.rank),"In iteration {}, The replica {} in rank {} before switch is {}\n".format(self.iteration, task_k,self.rank,task_v._getTemp()))
                    mylogger_hand("hand_log{}.log".format(self.rank),"In iteration {}, The replica {} in rank {} switch temperature from {} to {}\n".format(self.iteration,task_k,self.rank,old_temp,self.to_update[old_temp]))
                    mylogger_hand("hand_log{}.log".format(self.rank),"In iteration {}, Check: The replica {} in rank {} temperature {} to {},indeed is {}\n".format(self.iteration,task_k,self.rank,old_temp,self.to_update[old_temp],task_v._getTemp()))
                else:
                    self.logger.info("The replica {} in rank {} before switch is {}\n".format(task_k,self.rank,task_v._getTemp()))
                    self.logger.info("The replica {} in rank {} do not change temperature {}\n".format(task_k,self.rank,self.to_update[old_temp]))
                    self.logger.info("In iteration {}, Check: The replica {} in rank {} temperature {} to {},indeed is {}\n".format(self.iteration,task_k,self.rank,old_temp,self.to_update[old_temp],task_v._getTemp()))
                    mylogger_hand("hand_log{}.log".format(self.rank),"The replica {} in rank {} before switch is {}\n".format(task_k,self.rank,task_v._getTemp()))
                    mylogger_hand("hand_log{}.log".format(self.rank),"The replica {} in rank {} do not change temperature {}\n".format(task_k,self.rank,self.to_update[old_temp]))
                    mylogger_hand("hand_log{}.log".format(self.rank),"In iteration {}, Check: The replica {} in rank {} temperature {} to {},indeed is {}\n".format(self.iteration,task_k,self.rank,old_temp,self.to_update[old_temp],task_v._getTemp()))
            self.iteration += 1
            self.logger.info("The rank {} is in iteration of {}\n".format(self.rank,self.iteration))
            mylogger_hand("hand_log{}.log".format(self.rank),"The rank {} is in iteration of {}\n".format(self.rank,self.iteration))


    def run_one_cycle(self):
        # store some results in this list
        results = []
        # each node doing different jobs
        for replica_id, node_task_job in self.node_task.items():
            print("the rank {} is runing:".format(self.rank))
            node_task_job.run(self.cycle_per_steps)
            temp = node_task_job.temp._value
            reduceEnergy = [x for x in node_task_job.reduced_energy()]
            results.append([replica_id,temp,reduceEnergy])

        self.results_tmp = results
#       self.logger.debug("the results:{}".format(self.results_tmp))
#       mylogger_hand("hand_log{}.log".format(self.rank),"the results:{}".format(self.results_tmp))

        allresults = self.comm.gather(results,root=0)
#       self.logger.debug("all the results:{}".format(allresults))
#       mylogger_hand("hand_log{}.log".format(self.rank),"all the results:{}".format(allresults))

        if self.rank == 0 and self.rank_size > 1:
            # flatten the list
            allresults = [y for x in allresults for y in x]

#           self.logger.debug("all the results2:{}".format(allresults))
#           mylogger_hand("hand_log{}.log".format(self.rank),"all the results2:{}\n".format(allresults))
            # this make sure the sequence is sorted, the temrature is sorted from low to high
            # sort by second element, which is temperature (also because it is a number, while replica_id is a string)
            index_sort = np.argsort([x[1] for x in allresults])
#           self.logger.debug("index sort:{}".format(index_sort))
#           mylogger_hand("hand_log{}.log".format(self.rank),"index sort:{}\n".format(index_sort))
            # construct a new array for storage of every reduced energy for every state
            tmp_array = np.zeros([self.task_size,self.task_size])
            # assign the array
            replica_temp_tmp = []
            # replica_id do not change, which is corresponding to the array_index
            # result_index correspinding to the temperature
            for array_index,result_index in enumerate(index_sort):
#               self.logger.debug("allresults[result_index][2]:{}".format(allresults[result_index][2]))
#               mylogger_hand("hand_log{}.log".format(self.rank),"allresults[result_index][2]:{}\n".format(allresults[result_index][2]))
#               tmp_array[array_index] = allresults[result_index][2]._value
                # tmp_array store the data by sorted temperature, one raw for one temp
                # it will resutls a matrix or ndarray, the row represents temp(each simulation), each column represents estimated in diff temp
                tmp_array[array_index] = [x for x in allresults[result_index][2]]
                replica_temp_tmp.append([allresults[result_index][0],allresults[result_index][1]])
#           self.logger.debug("the tmp_array is:{}".format(tmp_array))
#           mylogger_hand("hand_log{}.log".format(self.rank),"the tmp_array is:{}\n".format(tmp_array))
#           self.logger.debug("the exchage pattern is:{}".format(self.exchage_pattern))
#           mylogger_hand("hand_log{}.log".format(self.rank),"the exchage pattern is:{}\n".format(self.exchage_pattern))
            self.all_system_energy = (self.iteration,tmp_array)
            self.replica_temp_pair = replica_temp_tmp
#           self.logger.debug("the replica_temp_pair is:{}".format(self.replica_temp_pair))
#           mylogger_hand("hand_log{}.log".format(self.rank),"the replica_temp_pair is:{}\n".format(self.replica_temp_pair))
        else:
            allresults = results
            index_sort = np.argsort([x[1] for x in allresults])
            tmp_array = np.zeros([self.task_size,self.task_size])
            replica_temp_tmp = []
            for array_index,result_index in enumerate(index_sort):
#               self.logger.debug("allresults[result_index][2]:{}".format(allresults[result_index][2]))
#               mylogger_hand("hand_log{}.log".format(self.rank),"allresults[result_index][2]:{}\n".format(allresults[result_index][2]))
#               tmp_array[array_index] = allresults[result_index][2]._value
                # tmp_array store the data by sorted temperature, one raw for one temp
                # it will resutls a matrix or ndarray, the row represents temp(each simulation), each column represents estimated in diff temp
                tmp_array[array_index] = [x for x in allresults[result_index][2]]
                replica_temp_tmp.append([allresults[result_index][0],allresults[result_index][1]])
            self.logger.debug("the tmp_array is:{}".format(tmp_array))
#           mylogger_hand("hand_log{}.log".format(self.rank),"the tmp_array is:{}\n".format(tmp_array))
            self.logger.debug("the exchage pattern is:{}".format(self.exchage_pattern))
#           mylogger_hand("hand_log{}.log".format(self.rank),"the exchage pattern is:{}\n".format(self.exchage_pattern))
            self.all_system_energy = (self.iteration,tmp_array)
            self.replica_temp_pair = replica_temp_tmp
#           self.logger.debug("the replica_temp_pair is:{}".format(self.replica_temp_pair))
#           mylogger_hand("hand_log{}.log".format(self.rank),"the replica_temp_pair is:{}\n".format(self.replica_temp_pair))

            
        # may be add a barrier here to make sure self.replica_temp_pair is updated
        self.comm.barrier()

    def root_to_other_node(self):
        """ This is used for the swap attempted """
        #check that the energy is equal to the current iteration cycle
        self.to_update = {}
        if self.rank == 0:
            assert self.all_system_energy[0] == self.iteration
            swap_results = self.swap()
#           self.logger.debug("the swap_results:{}".format(swap_results))
#           mylogger_hand("hand_log{}.log".format(self.rank),"the swap_results:{}\n".format(swap_results))
        else:
            swap_results = None

#       self.recvbuf = self.comm.scatter(swap_results, root=0)
        self.recvbuf = self.comm.bcast(swap_results, root=0)
        self.check_node_recvbuf()

#       self.logger.debug("the swap_results back:{}".format(self.recvbuf))
#       mylogger_hand("hand_log{}.log".format(self.rank),"the swap_results back:{}".format(self.recvbuf))
        for a,b,c in self.recvbuf:
            # c is to change or not
            # self.temps has unit in openmm
            # in self.to_update, the key represents old temperature, while the value represents the temperature in next cycle
            if c:
                self.to_update[self.temps[a]._value] =  self.temps[b]._value
                self.to_update[self.temps[b]._value] =  self.temps[a]._value
            else:
                self.to_update[self.temps[a]._value] =  self.temps[a]._value 
                self.to_update[self.temps[b]._value] =  self.temps[b]._value
        # if two terminal replica is not involed
        if not self.temps[0]._value in self.to_update: self.to_update[self.temps[0]._value] = self.temps[0]._value
        if not self.temps[-1]._value in self.to_update: self.to_update[self.temps[-1]._value] = self.temps[-1]._value
#       self.logger.debug("In iteration {} the to_update is:{}".format(self.iteration,self.to_update))
#       mylogger_hand("hand_log{}.log".format(self.rank),"In iteration {} the to_update is:{}".format(self.iteration,self.to_update))


    def swap(self):
        """ This is used for the swap attempted """
        # neighbour exchange for now, but it is easy for extenxiable
        success = []
        if np.random.choice([0,1]):
            for p in self.exchage_pattern[0]:
                # x is a list
                i,j = p
#               self.logger.debug("the swap i:{}".format(i))
#               self.logger.debug("the swap j:{}".format(j))
#               mylogger_hand("hand_log{}.log".format(self.rank),"the swap i:{}".format(i))
#               mylogger_hand("hand_log{}.log".format(self.rank),"the swap j:{}".format(j))
                success.append((i,j,self._swap_attempts(i,j)))
        else:
            for p in self.exchage_pattern[1]:
                i,j = p
#               self.logger.debug("the swap i:{}".format(i))
#               self.logger.debug("the swap j:{}".format(j))
#               mylogger_hand("hand_log{}.log".format(self.rank),"the swap i:{}".format(i))
#               mylogger_hand("hand_log{}.log".format(self.rank),"the swap j:{}".format(j))
                success.append((i,j,self._swap_attempts(i,j)))
        return success

    def _swap_attempts(self,i,j):
        energy_ij = self.all_system_energy[1][i,j] #replica i in temperature i estimated in temperture j
        energy_ji = self.all_system_energy[1][j,i]
        energy_ii = self.all_system_energy[1][i,i]
        energy_jj = self.all_system_energy[1][j,j]
        log_p_accept = - (energy_ij + energy_ji) + energy_ii + energy_jj
        if log_p_accept >= 0.0 or np.random.rand() < np.exp(log_p_accept):
            return True
        else:
            return False



    def accept(self):
        decision = []
        for i,j in proposal:
            decision.append(swap_attempts(i,j))
        return decision

    def _initialize(self, task_groups_infos):
        """
        This is the initialzation of REMD objects, it should be privided by single MD simulation projects, _initialize_temp_job should be handle the what kind of the projects
        """
        if isinstance(task_groups_infos, SlicableOrderedDict):
            raise "This is not a dict, the dict should be to store the objects that corresponding to the node"
        if self.rank_size == 1:
            print("This is just one node, all replica will load to one node, be carefull of your gpu memory")
        elif self.comm.Get_size() > 1:
            print("This is more than one node, all replica will load to one node,but also be carefull of your gpu memory")
        else:
            raise("Something is wrong about mpi, Please check")
            
        # to make sure the number of task size is equal to the number of privided temps
        assert self.task_size == len(self.temps)

        # attribute tasks into nodes corresponding to the env
        self._initialize_temp_job(task_groups_infos)


    def get_gpu_mem_info(self):
        import pycuda
        import pycuda.autoinit
        import pycuda.driver as cuda
        (free,total)=cuda.mem_get_info()
        print("Global memory occupancy:%f%% free"%(free*100/total))
        return free*100/total

#       for devicenum in range(cuda.Device.count()):
#           device=cuda.Device(devicenum)
#           attrs=device.get_attributes()

#           #Beyond this point is just pretty printing
#           print("\n===Attributes for device %d"%devicenum)
#           for (key,value) in attrs.iteritems():
#               print("%s:%s"%(str(key),str(value)))

    def _initialize_temp_job(self,task_groups_infos):
        if self.rank_size > 1:
            each_distri = np.ceil(self.task_size / self.rank_size)
            each_left = self.task_size % self.rank_size
            if each_distri == int(each_distri):
                print("the replicas is equlity distributed into all the ranks")
                rank_job_temp = task_groups_infos[self.rank*each_distri:self.rank*each_distri+each_distri]
            else:
                print("the replicas is not equlity distributed into all the ranks")
                print("it would be not better to do so!")
                if self.rank < each_left:
                    rank_job_temp = task_groups_infos[self.rank*each_distri:self.rank*each_distri+each_distri]
                else:
                    rank_job_temp = task_groups_infos[self.rank*each_distri:self.rank*each_distri+each_distri-1]
#               rank_job_temp = task_groups_infos[self.rank:self.rank+each_distri] if self.rank < each_left else task_groups_infos[self.rank:self.rank+each_distri-1]
            self.task_groups = {"rank":self.rank, "job":rank_job_temp}
            self.node_task = {}
            for replica_id, temp_now in self.task_groups["job"].items():
                temps = self.temps
                s_mod = self._getsystem()
                i = self._getIntegrator(temp_now)
                t = s_mod._psf.topology
                s = s_mod._system
                p = self._get_position()
                platform = self._getPlatform()
                self.node_task[replica_id] = mySimulation(t,s,i,dcdfrep=50000,repNamePrefix="run_{}".format(replica_id),temps=temps,platform=platform,positions=p,initReporter=True)

        else:
            rank_job_temp = task_groups_infos
            self.task_groups = {"rank":self.rank, "job":rank_job_temp}
            self.node_task = {}
            for replica_id, temp_now in self.task_groups["job"].items():
                temps = self.temps
                s_mod = self._getsystem()
                i = self._getIntegrator(temp_now)
                t = s_mod._psf.topology
                s = s_mod._system
                p = self._get_position()
                platform = self._getPlatform()
                self.node_task[replica_id] = mySimulation(t,s,i,temps=temps,platform=platform,positions=p,initReporter=True,platformProperties=self.platformProperties)
            

    def node_job_detect(self):
        print("{} job is: {}".format(self.rank,self.node_task.keys()))
        return "{} job is: {}".format(self.rank,self.node_task.keys())
                
            


    def _get_position(self):
        """
            This function were added for if the number of system is too many, use cache
        """
        if not self.restart:
            try:
                crd = self.crd
                positions = app.CharmmCrdFile(self.crd).positions
            except:
                pdb = self.pdb
                positions = app.pdbfile.PDBFile(pdb).positions
            return positions
            
        else:
            restart = self.restart_list
            return restart

    def _get_restart(self,**kwargs):
        pass

    def _getIntegrator(self,temp=300):
        if self.runtype == "drude":
            integrator = my_integrator(drude=True,temp=temp)
        elif self.runtype == "langevin":
            integrator = my_integrator(drude=False,temp=temp,integrator='langevin')
        return integrator.integrator


    def _getsystem(self):
        crd = self.crd
        psf = self.psf
        size = self.size
        params = self.params
        try:
             self.runtype
        except NameError: 
            raise "you have to input what force field type of systems is. Now avabilable is charmm36m and drude"
        if self.runtype == "drude":
            system = drudeSystem(psf,crd,params)
        elif self.runtype == "langevin":
            system = charmm36System(psf,crd,params)
        return system

    def _getPlatform(self):
        """
        the platform is GPU for default
        """
        try:
            self.platform
        except:
            self.platform = "CPU"
        platform = mm.Platform.getPlatformByName(self.platform)
        return platform
        


#   def _restart_temp_job():
#       pass

#   def _initialize_rest_job():
#       pass

#   def _initialize_reus_job():
#       pass

    def _generate_drude(self,system,temps):
        import copy 
        systems = {}
        for t in temps:
            systemtmp = copy.deepcopy(system)
#           t = t * kelvin
#           systemtmp.setDefaultTemperature(t)
            intergrator.setTemperature(t)
            systems[t] = systemtmp
        return systems
            

    def _generate_36m(self):
        import copy 
        systems = {}
        for t in temps:
            systemtmp = copy.deepcopy(system)
            for forcePosition, tmp in enumerate(systemInContext.getForces()):
                if  'Thermostat' in tmp.__class__.__name__:
                    break
            systemtmp.removeForce(forcePosition)
            newThermostat =  mm.AndersenThermostat(t*kelvin, 1.0/picosecond)
            systemtmp.addForce(newThermostat)
#           t = t * kelvin
#           systemtmp.setDefaultTemperature(t)
            systems[t] = systemtmp
        return systems

    def _generate_36m_langevin(self):
        import copy 
        systems = {}
        for t in temps:
            systemtmp = copy.deepcopy(system)
            for forcePosition, tmp in enumerate(systemInContext.getForces()):
                if  'Thermostat' in tmp.__class__.__name__:
                    break
            systemtmp.removeForce(forcePosition)
            newThermostat =  mm.AndersenThermostat(t*kelvin, 1.0/picosecond)
            systemtmp.addForce(newThermostat)
#           t = t * kelvin
#           systemtmp.setDefaultTemperature(t)
            systems[t] = systemtmp
        return systems

    def _task_divide(self):
        node_job_ids = [[self.task_groups["replica_"+str(x)] for x in j] for j in  node_job_index]
        return node_job_ids
    def _splitList(self,lst):
        # split the list by 2 elements
        f = lambda a:map(lambda b:a[b:b+2],range(0,len(a),2))
        tmp1 = list(f(lst))
        tmp2 = list(f(lst[1:]))
        # remove the list that has only one element
        for i,l in enumerate(tmp1):
            if len(l) == 1:
                tmp1.remove(l)
        for i,l in enumerate(tmp2):
            if len(l) == 1:
                tmp2.remove(l)
        # return the two choices
        return tmp1, tmp2

    def handle_args(self,*args):
        pass

    def handle_kargs(self,kwargs):
#       print(kwargs)
        print("run\n\n\n")
        for k,v in kwargs.items():
            if k == "crd":
                self.crd = v
            if k == "psf":
                self.psf = v
            if k == "pdb":
                self.pdb = v
            if k == "params":
                self.params = v
            if k == "size":
                self.size = v
            if k == "runtype":
                self.runtype = v

            if k == "restart":
                self.restart = v
            else:
                self.restart = None
            if k == "platform":
                self.platform = v
            if k == "platformProperties":
                self.platformProperties = v
            else:
                self.platformProperties = dict(CudaPrecision='single')

    def check_node_task(self):
        print("In rank {}, the node_task id is {}\n".format(self.rank, self.node_task.keys()))

#   def check_node_energy(self):
#       print("In rank {}, the node_task id is {}\n".format(self.rank, self.node_task.keys()))

    def check_rank1_info(self):
        if self.rank == 0:
            print("the system energy that is combine by rank 1 is:\n{}".format(self.all_system_energy))
            print("the system pair is now is: ".format(self.replica_temp_pair))
#print(bcolors.WARNING + "Warning: No active frommets remain. Continue?" + bcolors.ENDC)

    def check_node_recvbuf(self):
        self.logger.info("In rank: {} data is {}".format(self.rank,self.recvbuf))



    def check_run_type(self):
        if self.runtype == "drude": self.langevin = True
        if self.runtype == "langevin": self.langevin = True
        for k,v in self.node_task:
            assert v.langevin == self.runtype

    def check_PropertyValue(self):
        for replicaid, replica_task in self.node_task():
            assert self.platform.getPropertyValue(replica_task.context,'Precision') == self.platformProperties["CudaPrecision"]
            print("the presision in replicaid {} is {}".format(replicaid, self.platform.getPropertyValue(replica_task.context,'Precision')))


    def check_state(self):
        pass

