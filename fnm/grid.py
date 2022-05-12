import numpy as np
import addpaths
from dicts import dotdict

class grid():

  def __init__(self,*args,**kwargs):
    self.opt = dotdict({'nx' : 100,
                        'ny' : 1,
                        'nz' : 100,
                        'dx' : 1e-3,
                        'dy' : 1e-3,
                        'dz' : 1e-3,
                        'offset_x' : 0,
                        'offset_y' : 0,
                        'offset_z' : 0})
    self.opt.update(**kwargs)
  def values(self):
    nx = self.opt.nx
    ny = self.opt.ny
    nz = self.opt.nz
    dx = self.opt.dx
    dy = self.opt.dy
    dz = self.opt.dz
    offset_x = self.opt.offset_x
    offset_y = self.opt.offset_y
    offset_z = self.opt.offset_z

    # TODO: Support any plane
    xs1 = (np.r_[0:nx] - (nx-1.0)/2.0 + offset_x) * dx
    zs1 = (np.r_[0:nz] - (nz-1.0)/2.0 + offset_z) * dz

    xs,zs = np.meshgrid(xs1,zs1,indexing='ij')

    pos = np.c_[xs.flatten(), np.zeros(nx*nz), zs.flatten()].astype(np.float32)
    return pos
