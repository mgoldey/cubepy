import numpy as np
from copy import deepcopy
import constants
import view

def get_grid(mol,data):
    l1,l2,l3=mol.cell.diagonal()
    n1,n2,n3=mol.cell.diagonal()/data.shape
    return np.mgrid[0:l1:n1,0:l2:n2,0:l3:n3]


class molecule:
    """ Helper molecule class for holding data and organizing dependent functions """
    def __init__(self,atomids,coords,cell=np.identity(3)):
        self.natoms=len(atomids)
        self.atomids=atomids
        self.coords=coords
        self.cell=cell
        self.__center_of_mass()
    def __str__(self):
        ret_str=str(self.natoms)+"\n"
        ret_str+="\t".join([str(ic) for ic in self.cell.flatten()])+"\n"
        for i in xrange(self.natoms):
            ret_str+=constants.dict_of_atomic_abbr[self.atomids[i]]+"\t"+"\t".join([str(self.coords[i][j]) for j in range(3)])+"\n"
        return ret_str
    def __center_of_mass(self):
        """ This computes the centroid and center of mass using standard atomic masses """
        self.com=np.zeros(3)
        self.centroid=np.zeros(3)
        if len(self.coords)==0:
            return
        self.centroid=np.sum(self.coords)/len(self.coords)
        wts=np.array([constants.dict_of_atomic_masses[self.atomids[i]]  for i in range(self.natoms)])
        self.com=wts.dot(self.coords)/np.sum(wts)
    def write(self,fl):
        """ Write to xyz """
        os=open(fl,'w')
        odata=self.__str__()
        os.writelines(odata)
        os.close()
    def view(self):
        view.view_mol(self)
    def fixpbc(self):
        self.coords=np.imag(np.log(np.exp(2*np.pi*np.complex(1j)*self.coords.dot(np.linalg.inv(self.cell))))).dot(self.cell)/(2*np.pi)
        self.__center_of_mass()
    def move(self,vec):
        self.coords=self.coords+vec
    def copy(self,offset):
        imol=deepcopy(self)
        icoords=imol.coords
        rx=(offset[1]-offset[0])
        ry=(offset[3]-offset[2])
        rz=(offset[5]-offset[4])
        imol.natoms*=rx*ry*rz
        imol.cell[0]*=rx
        imol.cell[1]*=ry
        imol.cell[2]*=rz
        imol.atomids=np.array(imol.atomids.tolist()*rx*ry*rz)
        imol.coords=np.array(imol.coords.tolist()*rx*ry*rz).reshape(rx,ry,rz,self.natoms,3)

        for irx,rx in enumerate(range(offset[0],offset[1])):
            for iry,ry in enumerate(range(offset[2],offset[3])):
                for irz,rz in enumerate(range(offset[4],offset[5])):
                    imol.coords[irx,iry,irz]+=rx*self.cell[0]+ry*self.cell[1]+rz*self.cell[2]
        imol.coords=imol.coords.reshape(imol.natoms,3)
        return imol

    def center(self):
        self.fixpbc()
        self.move(-self.com)
        self.__center_of_mass()
