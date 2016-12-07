from __future__ import print_function

import numpy as np
from copy import deepcopy
import cubeutils.cons as cons

#import view

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
            ret_str+=cons.dict_of_atomic_abbr[self.atomids[i]]+"\t"+"\t".join([str(self.coords[i][j]) for j in range(3)])+"\n"
        return ret_str
    def __center_of_mass(self):
        """ This computes the centroid and center of mass using standard atomic masses """
        self.com=np.zeros(3)
        self.centroid=np.zeros(3)
        if len(self.coords)==0:
            return
        self.centroid=np.sum(self.coords)/len(self.coords)
        wts=np.array([cons.dict_of_atomic_masses[self.atomids[i]]  for i in range(self.natoms)])
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

def ipr(mol,data):
    val=(data**2).sum()/data.sum()**2
    dv=(mol.cell.diagonal()/data.shape).prod()
    return val/dv

def spread(mol,data):
    data2=deepcopy(data)
    denom=np.sum(data2)
    data2/=denom
    l1,l2,l3=mol.cell.diagonal()
    a,b,c=mol.cell.diagonal()/data2.shape
    nr1,nr2,nr3=data2.shape
    x,y,z = np.ogrid[0:l1:a,0:l2:b,0:l3:c]
    x=x.flatten()[:nr1]
    y=y.flatten()[:nr2]
    z=z.flatten()[:nr3]
    
    fac=2.0*np.pi/(l1)
    valx=np.dot(np.sum(np.sum(data2,axis=1),axis=1),(np.exp(x*fac*np.complex(1j))))
    valx=np.imag(np.log(valx))/fac

    fac=2.0*np.pi/(l2)
    valy=np.dot(np.sum(np.sum(data2,axis=2),axis=0),(np.exp(y*fac*np.complex(1j))))
    valy=np.imag(np.log(valy))/fac

    fac=2.0*np.pi/(l3)
    valz=np.dot(np.sum(np.sum(data2,axis=0),axis=0),(np.exp(z*fac*np.complex(1j))))
    valz=np.imag(np.log(valz))/fac

    if (valx<0.0):
        valx+=l1
    if (valy<0.0):
        valy+=l2
    if (valz<0.0):
        valz+=l3

    x1=x-valx
    y1=y-valy
    z1=z-valz


    r = lambda x,y,z: np.sqrt(x**2+y**2+z**2)

    dim=[l1,l2,l3]
    def box(x,y,z):
        x[np.where((x)>dim[0]/2.)]-=dim[0]
        y[np.where((y)>dim[1]/2.)]-=dim[1]
        z[np.where((z)>dim[2]/2.)]-=dim[2]

        x[np.where((x)<-dim[0]/2.)]+=dim[0]
        y[np.where((y)<-dim[1]/2.)]+=dim[1]
        z[np.where((z)<-dim[2]/2.)]+=dim[2]
        return

    box(x1,y1,z1)

    x2=x1**2
    y2=y1**2
    z2=z1**2
    val2x=np.sum(np.tensordot(data2,x2,axes=([0],[0])))
    val2y=np.sum(np.tensordot(data2,y2,axes=([1],[0])))
    val2z=np.sum(np.tensordot(data2,z2,axes=([2],[0])))

    print("x is centered at {0:6.3f} with spread {1:6.3f}".format(valx,np.sqrt(val2x)))
    print("y is centered at {0:6.3f} with spread {1:6.3f}".format(valy,np.sqrt(val2y)))
    print("z is centered at {0:6.3f} with spread {1:6.3f}".format(valz,np.sqrt(val2z)))
    print("gridded data is centered at ({0:6.3f},{1:6.3f},{2:6.3f}) with spread {3:6.3f}".format(valx,valy,valz,np.sqrt(val2x+val2y+val2z)))
