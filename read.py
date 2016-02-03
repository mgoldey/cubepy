import numpy as np
from utility import *
from constants import bohr_to_angstrom,dict_of_atomic_numbers
def read_cube(fl):
    stream=open(fl,'r')
    stream.readline()
    stream.readline()
    nat,xoff,yoff,zoff=stream.readline().split()
    cell=np.array([stream.readline().split() for i in range(3)],dtype=float).reshape((3,4))
    npts=cell[:,0]
    cell=cell[:,1:]*npts*bohr_to_angstrom
    mol=[]
    nat=int(nat)
    for i in range(nat):
        line=stream.readline().split()
        xyz=np.array(line,dtype=float)
        mol.append([xyz[0],xyz[-3],xyz[-2],xyz[-1]])
    mol=np.array(mol)
    mol=molecule(mol[:,0],bohr_to_angstrom*mol[:,1:],cell)
    data=np.array(stream.read().split(),dtype=float).reshape(np.array(npts,dtype=int))
    return mol,data

def read_xyz(fl):
    stream=open(fl,'r')
    nat=int(stream.readline())
    try:
        cell=np.array(stream.readline(),dtype=float).reshape((3,3))
    except:
        cell=np.ones(3)
    mol=np.array(stream.read().split()).reshape((nat,4))
    atomids=mol[:,0]
    coords=np.array(mol[:,1:],dtype=float)
    if (type(atomids[0][0])==str):
        atomids=np.array([dict_of_atomic_numbers[i] for i in atomids],dtype=int)
    mol=molecule(atomids,coords,cell)
    return mol

def read(fl):
    if (".cub" in fl):
        return read_cube(fl)
    if (".xyz" in fl):
        return read_xyz(fl)
