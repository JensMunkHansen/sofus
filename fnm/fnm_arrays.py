# Note: Not tested
import addpaths
import numpy as np
import swig_fnm as fnm

from dicts import dotdict
from euler import euler2rot

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import art3d
import matplotlib.colors as colors
plt.ion()

def show_rect(ax,corners):
    rect = art3d.Poly3DCollection([np.roll(corners,-2,axis=0)])
    rect.set_color(colors.rgb2hex([0.60,0.05,0.65]))
    rect.set_edgecolor('k')
    ax.add_collection3d(rect)

class rect (dotdict):
    def __init__(self,*args, **kwargs):
        opt = dotdict({'hw' : 0,
                       'hh' : 0,
                       'center'    : [0,0,0],
                       'euler'     : [0,0,0],
                       'conv'      : 'yxz',
                       'intrinsic' : True})
        opt.update(*args,**kwargs)
        super(rect,self).update(**opt)

        # Only 2 conventions are supported by Sofus
        conventionsSupported3 = ['yxy']
        conventionsSupported2 = ['yxz', 'yxy']
        if not(conventionsSupported2.count(opt.conv) > 0):
            raise Exception('Convention not supported')
        if (opt.euler[2] != 0.0):
            if not(conventionsSupported3.count(opt.conv) > 0):
                raise Exception('If convention yxz is used, z must be zero')

        self.center = np.array(self.center).flatten()
        self.euler  = np.array(self.euler).flatten()

    def corners(self):
        """
        Corners of the rectangle ordered clock-wise in the
        orientation dictated by the 'Euler' angles, the convention and
        the intrinsic property.
        """
        rotm = euler2rot(self.euler[0],self.euler[1],self.euler[2],
                         conv = self.conv,
                         intrinsic = self.intrinsic)
        # Ordered clock-wise
        a = np.array([[1,1],[-1,1],[-1,-1],[1,-1]])

        corners = np.c_[a * np.r_[self.hw,self.hh],[0,0,0,0]]
        corners = np.r_[[np.dot(rotm,corners[i]) + self.center for i in range(4)]]
        return corners

    def element(self):
        return np.r_[np.r_[self.hw], np.r_[self.hh], self.center, self.euler]

def field_array(*args,**kwargs):

    f2data = args[0]

    (nSubElements, nValues) = f2data.shape

    assert(nValues == 26)

    nElements    = int(f2data[nSubElements-1,0] + 1)
    nSubElements = int(max(f2data[:,1]) + 1)
    
    f2data = f2data.reshape((nElements,nSubElements,26))

    widths  = f2data[:,:,2]
    heights = f2data[:,:,3]
    tan_xz  = f2data[:,:,5]
    tan_yz  = f2data[:,:,6]

    apodization = f2data[:,:,4]
    center      = f2data[:,:,7:10]
    corners     = f2data[:,:,10:22].reshape(nElements,nSubElements,4,3)

    delays      = f2data[:,:,22]
    if 0:
        fig = plt.figure()
        ax = fig.add_subplot(121, projection='3d')

    rects = []
    for iElement in range(nElements):
        for jElement in range(nSubElements):
            if 0:
                show_rect(ax,corners[iElement,jElement])
  
    # Work for focused linear array
    if 0:
        ax = fig.add_subplot(122, projection='3d')

    # TODO: Figure out how (tan_yz,tan_xz) -> euler
    for iElement in range(nElements):
        for jElement in range(nSubElements):
            r = rect(hh=heights[iElement,jElement]/2.0,
                     hw=widths[iElement,jElement]/2.0,
                     center = center[iElement,jElement],
                     euler = [-tan_xz[iElement,jElement],tan_yz[iElement,jElement],0])
            rects.append(r)
            if 0:
                show_rect(ax,r.corners())


    elements = np.r_[[rects[i].element() for i in range(len(rects))]]
    elements = elements.reshape((nElements, nSubElements,8))

    a = fnm.ApertureFloat()
    a.subelements = elements.astype(np.float32)

    return a

#    widths = f2data[2

def linear_array(*args,**kwargs):
    """
    Key-value arguments: nElements, nSubH, pitch, kerf, height, focus
    """
    opt = dotdict({'nElements' : 192,
                   'nSubH'     : 1, # Elevation
                   'nSubW'     : 1,
                   'pitch'     : 0.2e-3,
                   'kerf'      : 0.2e-4,
                   'height'    : 1.0e-2,
                   'focus'     : None})

    opt.update(**kwargs)

    half_width  = (opt.pitch - opt.kerf) / 2.0
    half_height = opt.height / 2.0
    rects = []

    focus = 1.0
    R     = 1.0

    bFocused = opt.focus != None and opt.nSubH > 1
    
    if bFocused:
        focus = opt.focus
        elSector = 2 * np.arctan2(half_height, focus)
        R = np.sqrt(focus**2 + (opt.height / 2.0)**2)
    else:
        elSector = 0

    dEl =  elSector / max((opt.nSubH-1),1)
    elAngles = (np.r_[0:opt.nSubH] - (opt.nSubH - 1.0)/2) * dEl

    if bFocused:
        chordLength = 2 * focus * np.sin(dEl/2)
    else:
        chordLength = opt.height/opt.nSubH

    for iElement in range(opt.nElements):
        for iSubH in range(opt.nSubH):
            if bFocused:
                center = np.r_[(iElement - (opt.nElements-1.0)/2)*opt.pitch,
                               focus * np.tan(elAngles[iSubH]),
                               -(R * np.cos(elAngles[iSubH]) - focus)]
            else:
                center = np.r_[(iElement - (opt.nElements-1.0)/2)*opt.pitch,
                               (iSubH - (opt.nSubH - 1.0)/2.0)*chordLength,
                               0]
            for iSubW in range(opt.nSubW):
                center1 = center + 2 * half_width/opt.nSubW * np.r_[1,0,0] * (iSubW - (opt.nSubW-1)/2.0)
                r = rect(hw=half_width/opt.nSubW,hh=chordLength/2.0,
                         center=center1,
                         euler=[0,-elAngles[iSubH],0],conv='yxz',intrinsic=True)
                rects.append(r)

    elements = np.r_[[rects[i].element() for i in range(len(rects))]]
    elements = elements.reshape((opt.nElements, opt.nSubH * opt.nSubW, 8))

    a = fnm.ApertureFloat()
    a.subelements = elements.astype(np.float32)
    return a

# TODO: Set element on the inside instead, so no gaps
def convex_array(*args,**kwargs):
    opt = dotdict({'radius'    : 6.1e-2,
                   'nElements' : 192,
                   'nSubH'     : 1, # Elevation
                   'nSubW'     : 1,
                   'pitch'     : None,
                   'kerf'      : None,
                   'height'    : 1.0e-2,
                   'focus'     : 0.02,
                   'sector'    : 60.0/180 * np.pi})

    opt.updateset(**kwargs)

    if (opt.sector != None):
        dAz = opt.sector / max((opt.nElements - 1),1)
        if (opt.pitch == None):
            # We compute a pitch based on sector size
            opt.pitch = 2.0*opt.radius*np.sin(dAz/2.0)
            print('pitch is %f' % opt.pitch)
    else:
        # We compute sector based on pitch
        opt.sector = \
            (opt.nElements - 1) * 2.0 * np.arcsin(0.5 * opt.pitch / opt.radius)
        dAz = opt.sector / max((opt.nElements - 1),1)

    azAngles = (np.r_[0:opt.nElements] - (opt.nElements - 1.0)/2) * dAz

    rects = []
    el = 0.0

    elSector = 2 * np.arctan2(opt.height/2.0, opt.focus)
    dEl =  elSector / max((opt.nSubH-1),1)
    elAngles = (np.r_[0:opt.nSubH] - (opt.nSubH - 1.0)/2) * dEl

    chordLength = \
        {True : 2 * opt.focus * np.sin(dEl/2), False : opt.height}[opt.nSubH > 1]

    R = np.sqrt(opt.focus**2 + (opt.height / 2.0)**2)

    if opt.kerf == None:
        opt.kerf = 0.0

    half_width = (opt.pitch - opt.kerf) / 2.0

    for iAz in range(opt.nElements):
        for iEl in range(opt.nSubH):
            # Center position on element
            center = np.r_[0,0,opt.radius]
            
            # Candidate (no rotation is done around elevation focus)
            center = center + np.r_[0, opt.focus * np.tan(elAngles[iEl]), -(R * np.cos(elAngles[iEl]) - opt.focus) ]

            rotm = \
                euler2rot(azAngles[iAz],0,0,conv='yxz',intrinsic=True)

            # Rotate about origin
            center = np.dot(rotm,center)

            # Translate
            center = center - [0,0,opt.radius]

            # Translate sub-elements
            rots = euler2rot(azAngles[iAz],elAngles[iEl],0,conv='yxz',intrinsic=True)

            for iSubW in range(opt.nSubW):
                # Shift center coordinate
                center1 = center + 2 * half_width/opt.nSubW * rots[:,0] * (iSubW - (opt.nSubW-1)/2.0)
                
                r = rect(hw=half_width/opt.nSubW,hh=chordLength/2.0,
                         center=center1,
                         euler=[azAngles[iAz],elAngles[iEl],0],conv='yxz',intrinsic=True)
                rects.append(r)

    elements = np.r_[[rects[i].element() for i in range(len(rects))]]
    elements = elements.reshape((opt.nElements, opt.nSubH * opt.nSubW, 8))

    a = fnm.ApertureFloat()
    a.subelements = elements.astype(np.float32)
    return a
