import numpy as np
import scipy.special
import scipy.misc
from mayavi import mlab
from constants import *
from mayavi.tools.mlab_scene_model import MlabSceneModel
from traits.api import HasTraits, Range, Instance, on_trait_change,Int,Float,Button
import itertools
from traitsui.api import View, Item, HGroup,VGroup
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.core.ui.mayavi_scene import MayaviScene

def plotmol(scene,mol,offset=[0,1,0,1,0,1]):
    imol=mol.copy(offset)
    atoms=np.array(imol.atomids,dtype=int)
    fac=1.0
    sizes=covalent_radii[atoms]*fac
    x,y,z=imol.coords.T
    plts=[]
    for iat,i in enumerate(atoms):
        iplt=scene.mlab.points3d(x[iat],y[iat],z[iat],scale_factor=sizes[iat],resolution=40,color=atomic_colors[i],scale_mode='none')
        plts.append(iplt)
    return plts


def plotbonds(scene,mol,offset=[0,1,0,1,0,1]):
    imol=mol.copy(offset)
    thresh=1.2
    atoms=np.array(imol.atomids,dtype=int)
    fac=1.0
    sizes=covalent_radii[atoms]
    x,y,z=imol.coords.T
    plts=[]
    for i,j in itertools.combinations(range(imol.natoms),2):
        dist=np.linalg.norm(imol.coords[i]-imol.coords[j])
        dmin=(sizes[i]+sizes[j])*thresh
        if (dist<dmin):
            iplt=scene.mlab.plot3d(x[[i,j]],y[[i,j]],z[[i,j]] ,tube_radius=fac*0.08, color=tuple(0.5*(np.array(atomic_colors[atoms[i]])+np.array(atomic_colors[atoms[j]]))))
            plts.append(iplt)
    return plts


def plotdata(scene,data,dim=np.ones(3),contours=[],offset=[0,1,0,1,0,1],extent=[]):
    npts=np.array(data.shape,dtype=int)
    nx,ny,nz=np.complex(1j)*npts
    ix1,ix2=[(offset[0])*dim[0],(offset[1])*dim[0]]
    iy1,iy2=[(offset[2])*dim[1],(offset[3])*dim[1]]
    iz1,iz2=[(offset[4])*dim[2],(offset[5])*dim[2]]
    rx=(offset[1]-offset[0])
    ry=(offset[3]-offset[2])
    rz=(offset[5]-offset[4])
    nx*=rx
    ny*=ry
    nz*=rz
    x,y,z = np.mgrid[ix1:ix2:nx,iy1:iy2:ny,iz1:iz2:nz]
    ox,oy,oz=npts
    d2=np.zeros(x.shape)
    d2[:ox,:oy,:oz]=data
    
    for irx in range(rx):
        for iry in range(ry):
            for irz in range(rz):
                d2[(irx)*ox:(irx+1)*ox,(iry)*oy:(iry+1)*oy,(irz)*oz:(irz+1)*oz]=data

    if (len(contours)==0):
        vmax=data.max()*.5
        vmin=data.min()*.5
        contours=[vmin,vmax]
    if (len(contours)==2):
        if (len(extent)==0):
            plot1=scene.mlab.contour3d(x,y,z,d2.real,contours=[contours[0]],transparent=True,color=(1,0,0),opacity=0.3)
            plot2=scene.mlab.contour3d(x,y,z,d2.real,contours=[contours[1]],transparent=True,color=(0,0,1),opacity=0.3)
        else:
            plot1=scene.mlab.contour3d(x,y,z,d2.real,contours=[contours[0]],transparent=True,color=(1,0,0),opacity=0.3,extent=extent)
            plot2=scene.mlab.contour3d(x,y,z,d2.real,contours=[contours[1]],transparent=True,color=(0,0,1),opacity=0.3,extent=extent)
        return [plot1,plot2]
    
    plots=[]
    for ictr in contours:
        if (len(extent)==0):
            iplt=scene.mlab.contour3d(x,y,z,d2.real,contours=[ictr],transparent=True,opacity=0.3)
        else:
            iplt=scene.mlab.contour3d(x,y,z,d2.real,contours=[ictr],transparent=True,opacity=0.3,extent=extent)
        plots.append(iplt)
    return plots

def plotscalar(scene,data,dim=np.ones(3),offset=[0,1,0,1,0,1],lims=[]):
    if (lims==[]):
        lims=[0.9*data.min(),0.9*data.max()]
    #print offset, lims
    npts=np.array(data.shape,dtype=int)
    nx,ny,nz=np.complex(1j)*npts
    ix1,ix2=[(offset[0])*dim[0],(offset[1])*dim[0]]
    iy1,iy2=[(offset[2])*dim[1],(offset[3])*dim[1]]
    iz1,iz2=[(offset[4])*dim[2],(offset[5])*dim[2]]
    rx=(offset[1]-offset[0])
    ry=(offset[3]-offset[2])
    rz=(offset[5]-offset[4])
    nx*=rx
    ny*=ry
    nz*=rz
    x,y,z = np.mgrid[ix1:ix2:nx,iy1:iy2:ny,iz1:iz2:nz]
    ox,oy,oz=npts
    d2=np.zeros(x.shape)
    for irx in range(rx):
        for iry in range(ry):
            for irz in range(rz):
                d2[(irx)*ox:(irx+1)*ox,(iry)*oy:(iry+1)*oy,(irz)*oz:(irz+1)*oz]=data
    source = mlab.pipeline.scalar_field(x,y,z,d2)
    vol = scene.mlab.pipeline.volume(source,vmin=lims[0],vmax=lims[1])
    return [vol]


def view_data(data):
    mi=float(data.min())
    ma=float(data.max())
    class Visualization(HasTraits):
        upiso=Float(0.9*ma)
        downiso=Float(0.9*mi)

        upx,upy,upz=[Int(1) for i in range(3)]
        downx,downy,downz=[Int(0) for i in range(3)]

        scene = Instance(MlabSceneModel, ())
        plots=[]
        def __init__(self):
            HasTraits.__init__(self)
            for iplt in plotdata(self.scene,data,np.ones(3),[mi*0.9,ma*0.9]):
                self.plots.append(iplt)

        @on_trait_change('upiso,downiso',post_init=True)
        def update_plot(self):
            plt1,plt2=self.plots
            plt1.contour.contours=[self.downiso]
            plt2.contour.contours=[self.upiso]        
        
        @on_trait_change('upx,upy,upz,downx,downy,downz',post_init=True)
        def redo_plots(self):
            for iplt in self.plots:
                iplt.stop()
            nx1=self.downx
            nx2=self.upx
            ny1=self.downy
            ny2=self.upy
            nz1=self.downz
            nz2=self.upz
            offset=[nx1,nx2,ny1,ny2,nz1,nz2]
            self.plots=[]
            for iplt in plotdata(self.scene,data,np.ones[3],[self.downiso,self.upiso],offset):
                self.plots.append(iplt)
        # the layout of the dialog created
        view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),height=800, width=800, show_label=False),
               HGroup('upx','upy','upz','downx','downy','downz'),
               HGroup('downiso','upiso,fixpbc'), resizable=True,)

    visualization = Visualization()
    visualization.configure_traits()


def view_mol(mol):
    class Visualization(HasTraits):
        upx,upy,upz=[Int(1) for i in range(3)]
        downx,downy,downz=[Int(0) for i in range(3)]
        fixpbc=Button(label="Fix pbc")
        scene = Instance(MlabSceneModel, ())
        mplts=[]
        bplts=[]
        def __init__(self):
            HasTraits.__init__(self)
            self.mplts=plotmol(self.scene,mol)
            self.bplts=plotbonds(self.scene,mol)

        @on_trait_change('fixpbc',post_init=True)
        def fix_pbc(self):
            mol.fixpbc()
            plotmol(self.scene,mol)
            plotbonds(self.scene,mol)
        
        @on_trait_change('upx,upy,upz,downx,downy,downz',post_init=True)
        def redo_plots(self):
            for iplt in self.mplts:
                iplt.stop()
            for iplt in self.bplts:
                iplt.stop()
            nx1=self.downx
            nx2=self.upx
            ny1=self.downy
            ny2=self.upy
            nz1=self.downz
            nz2=self.upz
            offset=[nx1,nx2,ny1,ny2,nz1,nz2]
            
            self.mplts=plotmol(self.scene,mol,offset)
            self.bplts=plotbonds(self.scene,mol,offset)

        # the layout of the dialog created
        view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),height=800, width=800, show_label=False),
               HGroup('upx','upy','upz','downx','downy','downz'),
               HGroup('fixpbc'), resizable=True,)

    visualization = Visualization()
    visualization.configure_traits()

def view_vol(mol,data):
    mi=float(data.min())
    ma=float(data.max())
    class Visualization(HasTraits):
        upiso=Float(ma*0.9)
        downiso=Float(mi*0.9)
        fixpbc=Button(label="Fix pbc")

        upx,upy,upz=[Int(1) for i in range(3)]
        downx,downy,downz=[Int(0) for i in range(3)]

        scene = Instance(MlabSceneModel, ())
        plots=[]
        def __init__(self):
            HasTraits.__init__(self)
            self.mplts=plotmol(self.scene,mol)
            self.bplts=plotbonds(self.scene,mol)
            for iplt in plotscalar(self.scene,data,mol.cell.diagonal(),[0,1,0,1,0,1],[self.downiso,self.upiso]):
                self.plots.append(iplt)

        @on_trait_change('fixpbc',post_init=True)
        def fix_pbc(self):
            mol.fixpbc()
            plotmol(self.scene,mol)
            plotbonds(self.scene,mol)

        @on_trait_change('upiso,downiso,upx,upy,upz,downx,downy,downz',post_init=True)
        def redo_plots(self):
            for iplt in self.plots:
                iplt.stop()
            for iplt in self.mplts:
                iplt.stop()
            for iplt in self.bplts:
                iplt.stop()
            nx1=self.downx
            nx2=self.upx
            ny1=self.downy
            ny2=self.upy
            nz1=self.downz
            nz2=self.upz
            offset=[nx1,nx2,ny1,ny2,nz1,nz2]

            self.mplts=plotmol(self.scene,mol,offset)
            self.bplts=plotbonds(self.scene,mol,offset)
            
            self.plots=[]
            for iplt in plotscalar(self.scene,data,mol.cell.diagonal(),offset,[self.downiso,self.upiso]):
                self.plots.append(iplt)
        # the layout of the dialog created
        view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),height=800, width=800, show_label=False),
               HGroup('upx','upy','upz','downx','downy','downz'),
               HGroup('downiso','upiso'), resizable=True,)

    visualization = Visualization()
    visualization.configure_traits()


    scene = Instance(MlabSceneModel, ())
    plotmol(scene,mol)
    plotbonds(scene,mol)
    plotscalar(scene,data,mol.cell.diagonal())
    scene.mlab.colorbar()
    scene.mlab.show()

def iview(mol,data):
    mi=float(data.min())
    ma=float(data.max())
    class Visualization(HasTraits):
        upiso=Float(0.9*ma)
        downiso=Float(0.9*mi)

        upx,upy,upz=[Int(1) for i in range(3)]
        downx,downy,downz=[Int(0) for i in range(3)]
        fixpbc=Button(label="Fix pbc")

        scene = Instance(MlabSceneModel, ())
        plots=[]
        mplts=[]
        bplts=[]

        def __init__(self):
            HasTraits.__init__(self)
            self.mplts=plotmol(self.scene,mol)
            self.bplts=plotbonds(self.scene,mol)
            for iplt in plotdata(self.scene,data,mol.cell.diagonal(),[mi*0.9,ma*0.9]):
                self.plots.append(iplt)

        @on_trait_change('fixpbc',post_init=True)
        def fix_pbc(self):
            mol.fixpbc()
            plotmol(self.scene,mol)
            plotbonds(self.scene,mol)

        @on_trait_change('upiso,downiso',post_init=True)
        def update_plot(self):
            plt1,plt2=self.plots
            plt1.contour.contours=[self.downiso]
            plt2.contour.contours=[self.upiso]        
        
        @on_trait_change('upx,upy,upz,downx,downy,downz',post_init=True)
        def redo_plots(self):
            for iplt in self.plots:
                iplt.stop()
            for iplt in self.mplts:
                iplt.stop()
            for iplt in self.bplts:
                iplt.stop()
            nx1=self.downx
            nx2=self.upx
            ny1=self.downy
            ny2=self.upy
            nz1=self.downz
            nz2=self.upz
            offset=[nx1,nx2,ny1,ny2,nz1,nz2]
            
            self.mplts=plotmol(self.scene,mol,offset)
            self.bplts=plotbonds(self.scene,mol,offset)

            self.plots=[]
            for iplt in plotdata(self.scene,data,mol.cell.diagonal(),[self.downiso,self.upiso],offset):
                self.plots.append(iplt)
        # the layout of the dialog created
        view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),height=800, width=800, show_label=False),
               HGroup('upx','upy','upz','downx','downy','downz'),
               HGroup('downiso','upiso','fixpbc'), resizable=True,)

    visualization = Visualization()
    visualization.configure_traits()

