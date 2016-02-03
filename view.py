import numpy as np
import scipy.special
import scipy.misc
from mayavi import mlab
from constants import *

def plotbonds(mol):
    import itertools
    thresh=1.2
    atoms=np.array(mol.atomids,dtype=int)
    fac=1.0
    sizes=covalent_radii[atoms]
    x,y,z=mol.coords.T
    for i,j in itertools.combinations(range(mol.natoms),2):
        dist=np.linalg.norm(mol.coords[i]-mol.coords[j])
        dmin=(sizes[i]+sizes[j])*thresh
        if (dist<dmin):
            mlab.plot3d(x[[i,j]],y[[i,j]],z[[i,j]] ,tube_radius=fac*0.08, color=tuple(0.5*(np.array(atomic_colors[atoms[i]])+np.array(atomic_colors[atoms[j]]))))


def plotdata(data,dim=np.ones(3)):
    npts=np.array(data.shape,dtype=int)
    x,y,z = np.mgrid[0:dim[0]:np.complex(1j)*npts[0],0:dim[1]:np.complex(1j)*npts[1],0:dim[2]:np.complex(1j)*npts[2]]
    vmax=data.max()*.5
    vmin=data.min()*.5
    mlab.contour3d(x,y,z,data.real,contours=[vmax],transparent=True,color=(.5,0,0))
    mlab.contour3d(x,y,z,data.real,contours=[vmin],transparent=True,color=(0,0,.5))

def plotscalar(data):
    source = mlab.pipeline.scalar_field(data)
    min = data.min()
    max = data.max()
    vol = mlab.pipeline.volume(source)

def plotmol(mol):
    atoms=np.array(mol.atomids,dtype=int)
    fac=1.0
    sizes=covalent_radii[atoms]*fac
    x,y,z=mol.coords.T
    for iat,i in enumerate(atoms):
        mlab.points3d(x[iat],y[iat],z[iat],scale_factor=sizes[iat],resolution=40,color=atomic_colors[i],scale_mode='none')

def viewdata(data):
    mlab.figure()
    plotdata(data)
    mlab.colorbar()
    mlab.show()

def view(mol,data):
    mlab.figure()
    plotmol(mol)
    plotbonds(mol)
    plotdata(data,mol.cell.diagonal())
    mlab.colorbar()
    mlab.show()

def viewmol(mol):
    mlab.figure()
    plotmol(mol)
    plotbonds(mol)
    mlab.show()

def viewvol(mol,data):
    mlab.figure()
    plotmol(mol,np.array(data.shape))
    plotbonds(mol,np.array(data.shape))
    plotscalar(data)
    mlab.colorbar()
    mlab.show()

def iview(mol,data):
    import itertools
    from traits.api import HasTraits, Range, Instance, \
                        on_trait_change
    from traitsui.api import View, Item, HGroup
    from tvtk.pyface.scene_editor import SceneEditor
    from mayavi.tools.mlab_scene_model import \
                        MlabSceneModel
    from mayavi.core.ui.mayavi_scene import MayaviScene
    npts=np.array(data.shape,dtype=int)
    dim=mol.cell.diagonal()
    x,y,z = np.mgrid[0:dim[0]:np.complex(1j)*npts[0],0:dim[1]:np.complex(1j)*npts[1],0:dim[2]:np.complex(1j)*npts[2]]
    mi=float(data.min())
    ma=float(data.max())
    class Visualization(HasTraits):
        upiso=Range(0.0,ma,0.9*ma)
        downiso=Range(mi,0.0,0.9*mi)
        scene = Instance(MlabSceneModel, ())

        def __init__(self):
            # Do not forget to call the parent's __init__
            HasTraits.__init__(self)
            atoms=np.array(mol.atomids,dtype=int)
            fac=1.0
            sizes=covalent_radii[atoms]*fac
            mx,my,mz=mol.coords.T
            for iat,i in enumerate(atoms):
                self.scene.mlab.points3d(mx[iat],my[iat],mz[iat],scale_factor=sizes[iat],resolution=40,color=atomic_colors[i],scale_mode='none')
            for i,j in itertools.combinations(range(mol.natoms),2):
                dist=np.linalg.norm(mol.coords[i]-mol.coords[j])
                dmin=(sizes[i]+sizes[j])*1.2
                if (dist<dmin):
                    self.scene.mlab.plot3d(mx[[i,j]],my[[i,j]],mz[[i,j]] ,tube_radius=fac*0.08, color=tuple(0.5*(np.array(atomic_colors[atoms[i]])+np.array(atomic_colors[atoms[j]]))))

            vmax=ma*.9
            vmin=mi*.9
            self.plot=self.scene.mlab.contour3d(x,y,z,data.real,contours=[vmax,vmin],transparent=True,vmin=mi,vmax=ma)
            #self.plot2=self.scene.mlab.contour3d(x,y,z,data.real,contours=[vmin],transparent=True,color=(0,0,.5))

        @on_trait_change('upiso')
        def update_plot1(self):
            self.plot.contour.contours=[self.upiso,self.downiso]
            self.scene.mlab.colorbar(orientation='vertical')
        
        @on_trait_change('downiso')
        def update_plot2(self):
            self.plot.contour.contours=[self.upiso,self.downiso]
            self.scene.mlab.colorbar(orientation='vertical')

        # the layout of the dialog created
        view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                        height=500, width=500, show_label=False),
                    HGroup('_','downiso',),
                    HGroup('_','upiso',),
                    #HGroup('_', 'upiso', 'downiso',),
                    )

    visualization = Visualization()
    visualization.configure_traits()

