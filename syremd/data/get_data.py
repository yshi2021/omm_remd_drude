import os
import re

def get_data(*args,**kwargs):
#   villin = kwargs.get('villlin')
#   if kwargs.get('drude'):
#       villin =
#   
#   
    pass


def get_villin(drude=True):
    directory, _ = os.path.split(os.path.abspath(__file__))
    if drude:
        d = os.path.join(directory,"villin","drude")
    files = os.listdir(d)
    pdb = [os.path.join(d,x) for x in files if x.endswith(".pdb")]
    crd = [os.path.join(d,x) for x in files if x.endswith(".crd")]
    psf = [os.path.join(d,x) for x in files if x.endswith(".psf")]
    rst = [os.path.join(d,x) for x in files if x.endswith(".rst")]
    return psf[0], crd[0], pdb[0],rst[0]
    


def get_drude_params(version="2019f"):
    directory, _ = os.path.split(os.path.abspath(__file__))
    d = os.path.join(directory,"drude_toppar")
    if version == "2019f":
        files = os.listdir(d)
    files = [os.path.join(d,x) for x in files]
    return files



def get_36m_params():
    directory, _ = os.path.split(os.path.abspath(__file__))
    d = os.path.join(directory,"toppar")
    files = os.listdir(d)
    files = [os.path.join(d,x) for x in files]
    return files



def get_aaqaa(drude=True):
    directory, _ = os.path.split(os.path.abspath(__file__))
    if drude:
        d = os.path.join(directory,"aaqaa","aaqaadrude2019")
    else:
        d = os.path.join(directory,"aaqaa","aaqaa36m")
    files = os.listdir(d)
    pdb = [os.path.join(d,x) for x in files if x.endswith(".pdb")]
    crd = [os.path.join(d,x) for x in files if x.endswith(".crd")]
    psf = [os.path.join(d,x) for x in files if x.endswith(".psf")]
    rst = [os.path.join(d,x) for x in files if x.endswith(".rst")]
    return psf, crd, pdb, rst


def get_amber_aaqaa():
    directory, _ = os.path.split(os.path.abspath(__file__))
    d = os.path.join(directory,"aaqaa","amber")
    files = os.listdir(d)
    prmtop = [os.path.join(d,x) for x in files if x.endswith(".prmtop")]
    inpcrd = [os.path.join(d,x) for x in files if x.endswith(".inpcrd")]
    return prmtop,inpcrd


