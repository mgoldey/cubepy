import numpy as np
from copy import deepcopy
import constants
import view

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
        view.view(self)
    def fixpbc(self):
        self.coords=np.imag(np.log(np.exp(2*np.pi*np.complex(1j)*self.coords.dot(np.linalg.inv(self.cell))))).dot(self.cell)/(2*np.pi)
        self.__center_of_mass()
    def move(self,vec):
        self.coords=self.coords+vec
    def center(self):
        self.fixpbc()
        self.move(-self.com)
        self.__center_of_mass()
