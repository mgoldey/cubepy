import numpy as np
from cubeutils.utility import molecule
from cubeutils.cons import bohr_to_angstrom,dict_of_atomic_numbers
from cubeutils.cons import angstrom_to_bohr

def read_cube(fl): 
    """ Reads cube files"""
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

def read_pp(fl):
    """ Reads tmp.pp files from quantum espresso and converts them to the same ordering as cube files"""
    stream=open(fl,'r')
    stream.readline()

    data=np.array(stream.readline().split(),dtype=int)
    npts=data[:3]
    ntyp=data[-1]
    nat=data[-2]
    types=np.arange(ntyp)+1
    data=np.array(stream.readline().split(),dtype=float)
    if (int(data[0])==0):
        A=data[1]*bohr_to_angstrom
        cell=np.array([stream.readline().split() for i in range(3)],dtype=float)
        cell*=A
    elif (int(data[0])==1):
        A=data[1]*bohr_to_angstrom
        cell=np.identity(3)*A
    stream.readline()


    atoms=[]
    for iat in range(ntyp):
        atoms.append(stream.readline().split()[1])
    atom_dict=dict(zip(types,atoms))

    atomids=[]
    coords=[]
    for iat in range(nat):
        data=np.array(stream.readline().split(),dtype=float)
        atomids.append(dict_of_atomic_numbers[atom_dict[data[-1]]])
        coords.append((A*data[1:4]).tolist())
    coords=np.array(coords)

    mol=molecule(atomids,coords,cell)
    data=np.array(stream.read().split(),dtype=float).reshape(np.array(npts,dtype=int).T).T
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
    if ("pp" in fl):  # DANGEROUS DUE TO LACK OF FILE ENDING UNIFORMITY
        return read_pp(fl)

def write_cube(mol,data,file_name):
    """ Writes cubefile with mol and data as passed to it """
    outf=open(file_name+".cub",'w')
    outf.writelines("CUBE FILE WRITTEN USING MBG'S CUBEUTILS CODE VERSION 0.1\n")
    outf.writelines("YOUR MILEAGE MAY VARY\n")
    outf.writelines(str(mol.natoms)+"\t0.0\t0.0\t0.0\n") # natoms x y z (offset)
    npts=np.array(data.shape,dtype=float)
    cell=npts**-1*mol.cell*angstrom_to_bohr
    coords=mol.coords*angstrom_to_bohr
    for i in range(3):
        ic=cell[i]
        outf.writelines(str(int(data.shape[i]))+"\t"+str(ic[0])+"\t"+str(ic[1])+"\t"+str(ic[2])+"\n")
    for i,iat in enumerate(coords):
        outf.writelines(str(int(mol.atomids[i]))+"\t"+str(float(mol.atomids[i]))+"\t"+str(iat[0])+"\t"+str(iat[1])+"\t"+str(iat[2])+"\n")
    ictr=0
    for i in data:
        for j in i:
            for k in j:
                outf.write(str(k))
                if (ictr==6):
                    outf.write("\n")
                    ictr=0
                else:
                    outf.write("\t")
                    ictr+=1

    outf.close()

def write(mol,file_name):
    """ Writes filename using mol intrinsic function """
    mol.write(file_name)

