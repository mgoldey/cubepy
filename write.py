import numpy as np
from constants import angstrom_to_bohr
def write_cube(mol,data,file_name):
    """ Writes cubefile with mol and data as passed to it """
    outf=open(file_name+".cub",'w')
    outf.writelines("CUBE FILE WRITTEN USING MBG CUBEPY CODE VERSION 0.01\n")
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

