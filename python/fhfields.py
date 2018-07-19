import os, sys
import abc
import numpy as np

import matplotlib.pyplot as plt

from dicts import dotdict

filedir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(filedir, '../fnm'))

plt.ion()

# In Python3, we can do
#
# class FraunhoferApertureIF(dotdict, metaclass=abc.ABCMeta):
#   @abc.abstractmethod
#   def extent(self):
#     raise NotImplementedError("Implement extent returning a 3x2 array")

class FraunhoferApertureIF(dotdict):
  __metaclass__ = abc.ABCMeta
  @abc.abstractmethod
  def ExtentGet(self):
    raise NotImplementedError("Implement extent returning a 3x2 array")

  @abc.abstractmethod
  def F0Get(self):
    raise NotImplementedError

  @abc.abstractmethod
  def FocusGet(self):
    raise NotImplementedError

  @abc.abstractproperty
  def c(self):
    raise NotImplementedError

  @property
  def extent(self):
    return self.ExtentGet()

  @property
  def f0(self):
    return self.F0Get()

  @property
  def focus(self):
    return self.FocusGet()

  @property
  def bandwidth(self):
    raise NotImplementedError

def trondwin(N,beta):
  x=np.linspace(-beta,beta,N)
  w=np.exp(-x**4)
  return w

def tukeywin(window_length, alpha=0.5):
  '''
  The Tukey window, also known as the tapered cosine window, can be
  regarded as a cosine lobe of width \alpha * N / 2 that is convolved
  with a rectangle window of width (1 - \alpha / 2). At \alpha = 1 it
  becomes rectangular, and at \alpha = 0 it becomes a Hann window.

  We use the same reference as MATLAB to provide the same results in
  case users compare a MATLAB output to this function output

  Reference
  ---------
  http://www.mathworks.com/access/helpdesk/help/toolbox/signal/tukeywin.html

  '''
  # Special cases
  if alpha <= 0:
    return np.ones(window_length) #rectangular window
  elif alpha >= 1:
    return np.hanning(window_length)

  # Normal case
  x = np.linspace(0, 1, window_length)
  w = np.ones(x.shape)

  # first condition 0 <= x < alpha/2
  first_condition = x<alpha/2
  w[first_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[first_condition] - alpha/2) ))

  # second condition already taken care of

  # third condition 1 - alpha / 2 <= x <= 1
  third_condition = x>=(1 - alpha/2)
  w[third_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[third_condition] - 1 + alpha/2)))
  return w

class Concrete(FraunhoferApertureIF):
  def __init__(self, **kwargs):
    opt = dotdict({'pitch0' : 0.2e-3,
                   'N0'     : 64,
                   'pitch1' : 0.2e-3,
                   'N1'     : 64,
                   'f0'     : 7e6,
                   'focus'  : np.r_[0.0,0.0,3.0e-2]})
    opt.update(**kwargs)
    self._f0 = opt.f0
    self._focus = opt.focus
    super(Concrete, self).update(**opt)

  def ExtentGet(self):
    retval = np.zeros((3,2))
    retval[0,0] = - self.N0*self.pitch0 / 2.0
    retval[0,1] = self.N0*self.pitch0 / 2.0
    retval[1,0] = - self.N1*self.pitch1 / 2.0
    retval[1,1] = self.N1*self.pitch1 / 2.0
    return retval

  def F0Get(self):
    return self._f0

  def FocusGet(self):
    return self._focus

  @property
  def c(self):
    return 1540.0

  @property
  def bandwidth(self):
    return 0.0

class FraunhoferField:
  @staticmethod
  def FieldAtFocusNorway(xdc, **kwargs):
    opt = dotdict({'Apod'         : np.r_[0.0,0.0], # Rectangular windows
                   'oversampling' : 1,
                   'Unused'       : None,
                   'P'            : None})

    oversampling = opt.oversampling

    c = xdc.c

    # Focus distance (azimuth, elevation) - assumed equal
    Rx = xdc.focus[2]
    Ry = xdc.focus[2]

    # Diameter (azimuth, elevation)
    [Dx,Dy,_] = xdc.extent[:,1] - xdc.extent[:,0]

    # Resolution
    resolution = np.r_[6e-5, 6e-5, 6e-5]

    # Range (x,y,z)
    rng=np.array([0.006, 0.006, 0.0006]) #in x, y, z direction

    # (Rxy, Dxy, rng, bwm f0)
    Rxy = np.r_[Rx, Ry]
    R = np.mean(Rxy)
    [betax, betay] = opt.Apod

    circsymm = 1

    bw = xdc.bandwidth
    f0 = xdc.f0

    if opt.P == None:
      df=c/2.0/rng[2]
      Nf=np.round(bw/df) + 1
      pData = dotdict({})
      pData.f1 = np.linspace(f0 - bw/2, f0 + bw/2, Nf)
      pData.P  = trondwin(len(pData.f1),1.5).T # Trond W. variant
      res = pData.P
      opt.f = pData.f1
      opt.P = pData.P

    e = None
    f = opt.f
    P = opt.P

    # TODO(JMH): Enable this
    ROCcorrection = 0
    sphericalCorrection = 0

    testPlot = False

    if (len(f) == 1):
      fmax = f[0]
    else:
      fmax = 1.1*f[-1]

    # Always odd number FFT (TODO: Specify freely)
    Nx=1+2*np.round(rng[0]*fmax*Dx/R/c)
    Ny=1+2*np.round(rng[1]*fmax*Dy/R/c)
    fx=np.linspace(-fmax*Dx/R,fmax*Dx/R,Nx)
    fy=np.linspace(-fmax*Dy/R,fmax*Dy/R,Ny)
    fxi=np.linspace(-fmax*Dx/R,fmax*Dx/R,1+oversampling*(Nx-1))[None,:]
    fyi=np.linspace(-fmax*Dy/R,fmax*Dy/R,1+oversampling*(Ny-1))[None,:]
    fxm=np.dot(fxi.T,np.ones((1,fyi.shape[1])))
    fym=np.dot(np.ones((fxi.shape[1],1)),fyi)

    currfig=plt.gcf()
    hlp2=np.ones((oversampling,oversampling))/(oversampling)**2#2D lpfilter
    P1=np.zeros((len(fx),len(fy),len(f)),dtype=np.complex128)
    P1i0=np.zeros((len(fxi.T),len(fyi.T)),dtype=np.complex128)

    for k in range(len(f)):
      fxmax=0.5*Dx/R*f[k]
      fymax=0.5*Dy/R*f[k]
      Ix=np.where(np.abs(fxi)<fxmax)[1]
      Iy=np.where(np.abs(fyi)<fymax)[1]
      wx=tukeywin(len(Ix),betax) # wx[0]=0.5*wx[0]wx[-1]=0.5*wx[-1]
      wy=tukeywin(len(Iy),betay) # wy[0]=0.5*wy[0]wy[-1]=0.5*wy[-1]
      w=np.dot(wx[:,None],wy[None,:])
      if circsymm:
        w=tukey2(len(Ix),betax)
      P1k=P1i0
      P1k[Ix[0]:Ix[-1]+1,Iy[0]:Iy[-1]+1] = w*P[k]
      if e != None:
        Iblank=np.where((np.abs(fxm/fxmax)**e+np.abs(fym/fymax)**e)>1)[0]
        P1k[Iblank]=0
      if oversampling>1:
        P1k=convolve2d(P1k,hlp2,'same')
      if 0:
        figure(7)
        #subplot(2,2,1)
        #imagesc(np.abs(P1k))axis('image')colormap(gray)gridpause

      if ROCcorrection:
        Ha=apertureCorrection(R,Rxy[0],Rxy[1],f[k],fxi,fyi)
        P1k=P1k*Ha
      if sphericalCorrection:
        Hs=focalcorrection(R,f[k],fxi,fyi)
        P1u=P1k
        P1k=conv2(P1k,Hs,'same')
        if testPlot:
          figure(7)
          #subplot(2,2,1)imagesc(np.abs(P1u))axis('image')colormap(gray)grid
          #subplot(2,2,3)imagesc(np.abs(P1k))axis('image')
          #subplot(2,2,2)imagesc(np.abs(Hs))axis('image')
          #subplot(2,2,4)imagesc(np.abs(ifft(fftshift(P1k))))axis('image')
      P1k=P1k[0::oversampling,0::oversampling]#decimation to correct sampling rate
      P1[:,:,k]=P1k
    #k=1:len(f),
    if testPlot:
      if currfig!=None:
        figure(currfig)
    return (P1,fx,fy,f)
  def FieldAtFocus(xdc, **kwargs):
    opt = dotdict({'nx' : 100,
                   'ny' : 100,
                   'dx' : 6e-5,
                   'dy' : 6e-5})
    # Maximum frequency
    fmax = None
    if xdc.bandwidth == 0.0:
      fmax = xdc.f0




# Local variables: #
# tab-width: 2 #
# python-indent: 2 #
# indent-tabs-mode: nil #
# End: #
