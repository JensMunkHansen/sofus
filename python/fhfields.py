import os, sys
import abc
import numpy as np

from scipy.signal import convolve2d

import matplotlib.pyplot as plt

from scipy.fftpack import (fft, ifft, ifftn, fft2, ifft2)

filedir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(filedir, '../fnm'))

from dicts import dotdict

plt.ion()

talk = True

# In Python3, we can do
#
# class FraunhoferApertureIF(dotdict, metaclass=abc.ABCMeta):
#   @abc.abstractmethod
#   def extent(self):
#     raise NotImplementedError("Implement extent returning a 3x2 array")

class FraunhoferApertureIF(object):
  """
  Abstract base class for aperture
  """
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

#  @property
#  def bandwidth(self):
#    raise NotImplementedError

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

def tukey2(N,R):
  return tukey2ht(N,R)
#  return tukey2jmh(N,R)

def tukey2ht(N,R):
  """
  Circular Tukey window - implementation from Hans Torp
  """
  Nl=512;
  ic=np.round(N/2);
  x = np.dot(np.ones((N,1)), np.linspace(-1,1,N)[None,:])
  y = np.dot(np.linspace(-1,1,N)[:,None], np.ones((1,N)))
  r = np.sqrt(x**2 + y**2)
  Ir = np.where(r<1)
  wN = tukeywin(Nl,R)
  w = 0*x
  w[Ir] = wN[np.round(0.49*Nl*(1+r[Ir])).astype(np.int64)]
  return w

def tukey2jmh(N,alpha):
  """
  Circular Tukey windows - analytic (does not work)
  """
  if alpha <= 0:
    return np.ones((N,N))
  x = np.dot(np.ones((N,1)), np.linspace(-1,1,N)[None,:])
  y = np.dot(np.linspace(-1,1,N)[:,None], np.ones((1,N)))
  r = np.sqrt(x**2 + y**2)
  w = 0*x

  # first condition 0 <= r < alpha/2
  first_condition = np.where(np.logical_and((r<1), (0.5*(1-r) < alpha/2)))
  w[first_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * ( 0.5*(1.0 - r[first_condition]) - alpha/2) ))

  # second condition already taken care of

  # third condition 1 - alpha / 2 <= r <= 1
  third_condition = np.where(np.logical_and((r<1),(0.5*(1-r) >= (1 - alpha/2))), (0.5*(1-r) >= alpha/2)  )
  w[third_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * ( 0.5*(1.0 - r[third_condition]) - 1 + alpha/2)))

  return w

def apertureCorrection(R,Rx,Ry,f,fxv,fyv):
  c = 1540.0
  fx=np.dot(fxv.reshape((fxv.shape[1],1)),np.ones(fyv.shape))
  fy=np.dot(np.ones((fxv.shape[1],1)),fyv)
  fx2=(fx/f)**2;
  fy2=(fy/f)**2;
  delta=R*(1-np.sqrt(1-fx2-fy2)) - Rx*(1-np.sqrt(1-(R/Rx)**2*fx2)) - Ry*(1-np.sqrt(1-(R/Ry)**2*fy2))
  Ha=np.exp(-0.5*1j*2*np.pi*f/c*delta)
  return Ha

def focalcorrection(R,f0,fx,fy):
  # correction of phase due to curvature of focal "plane" for single frequency f0
  c = 1540.0
  dfx = fx[0,1] - fx[0,0]
  dfy = fy[0,1] - fy[0,0]
  xmax = 1.0/dfx * c/2.0
  ymax = 1.0/dfy * c/2.0
  # This causes asymmetry
  xax=np.linspace(-xmax,xmax,fx.shape[1]) - 2*xmax/fx.shape[1] #gjengs konstant trukket fra
  yax=np.linspace(-ymax,ymax,fy.shape[1])- 2*ymax/fy.shape[1]

  [x,y] = np.meshgrid(xax,yax,indexing='ij')
  #x=xax'*ones(1,length(yax));
  #y=ones(length(xax),1)*yax;

  tau=(R-np.sqrt(R**2 - x**2 - y**2))/c
  h = np.exp(-1j*2*np.pi*f0*tau)
  tind=0.8
  w=np.dot(tukeywin(len(xax),tind), tukeywin(len(yax),tind).T)
  # subplot(2,2,4);imagesc(w);
  h=h*w # skip weighting
  h=np.fft.fftshift(h)
  Hs=np.fft.fft2(h)
  Hs=np.fft.fftshift(Hs)
  return Hs

class FraunhoferAperture(FraunhoferApertureIF):
  def __init__(self, **kwargs):
    opt = dotdict({'pitch0' : 0.2e-3,
                   'N0'     : 64,
                   'pitch1' : 0.2e-3,
                   'N1'     : 64,
                   'f0'     : 7e6,
                   'focus'  : np.r_[0.0,0.0,3.0e-2]})
    opt.update(**kwargs)

    self._f0        = opt.f0
    self._focus     = opt.focus
    self._bandwidth = 1.0
    self.N0 = opt.N0
    self.pitch0 = opt.pitch0
    self.N1 = opt.N1
    self.pitch1 = opt.pitch1

  def ExtentGet(self):
    retval = np.zeros((3,2))
    retval[0,0] = - self.N0*self.pitch0 / 2.0
    retval[0,1] =   self.N0*self.pitch0 / 2.0
    retval[1,0] = - self.N1*self.pitch1 / 2.0
    retval[1,1] =   self.N1*self.pitch1 / 2.0
    return retval

  def F0Get(self):
    return self._f0

  def FocusGet(self):
    return self._focus

  @property
  def c(self):
    return 1540.0

#   def getbandwidth(self):
#     return self._bandwidth
#   def setbandwidth(self, value):
#     self._bandwidth = value
#   def delbandwidth(self):
#     del self._bandwidth
#   bandwidth = property(getbandwidth, setbandwidth, delbandwidth, "I'm the 'x' property.")
  @property
  def bandwidth(self):
    return self._bandwidth

  @bandwidth.setter
  def bandwidth(self, bandwidth):
    self._bandwidth = bandwidth

class FraunhoferField:
  @staticmethod
  def FieldAtFocusNorway(xdc, **kwargs):
    """
    Field computed in the focal plane using the Fraunhofer approximation.

    This implementation is taken from Hans Torp, Trondheim.

    K-space is not computed properly to be used with FFT's.
    """
    opt = dotdict({'Apod'         : np.r_[0.0,0.0], # Rectangular windows
                   'oversampling' : 1,
                   'P'            : None,
                   'rng'          : np.array([0.006, 0.006, 0.0006])})

    oversampling = opt.oversampling

    c = xdc.c

    # Focus distance (azimuth, elevation) - assumed equal for now
    Rx = xdc.focus[2]
    Ry = xdc.focus[2]

    # Diameter (azimuth, elevation)
    [Dx,Dy,_] = xdc.extent[:,1] - xdc.extent[:,0]

    # Range (x,y,z)
    rng = opt.rng

    Rxy = np.r_[Rx, Ry]
    R = np.mean(Rxy)
    [betax, betay] = opt.Apod

    circsymm = 1

    bw = xdc.f0 * xdc.bandwidth
    f0 = xdc.f0

    if opt.P == None:
      df=c/2.0/rng[2]
      Nf=np.round(bw/df) + 1
      pData = dotdict({})
      pData.f1 = np.linspace(f0 - bw/2, f0 + bw/2, Nf)
      # Trond W. variant - exp(x**4), -1.5 < x < 1.5
      pData.P  = trondwin(len(pData.f1),1.5).T
      res = pData.P
      opt.f = pData.f1
      opt.P = pData.P

    e = None
    f = opt.f
    P = opt.P

    # TODO(JMH): Enable this
    ROCcorrection = 1
    sphericalCorrection = 1

    if (len(f) == 1):
      fmax = f[0]
    else:
      fmax = 1.1*f[-1]

    # Always odd number FFT (TODO: Specify freely instead of using Dx, Dy)
    Nx=1+2*np.round(rng[0]*fmax*Dx/R/c)
    Ny=1+2*np.round(rng[1]*fmax*Dy/R/c)
    fx=np.linspace(-fmax*Dx/R,fmax*Dx/R,Nx)
    fy=np.linspace(-fmax*Dy/R,fmax*Dy/R,Ny)
    fxi=np.linspace(-fmax*Dx/R,fmax*Dx/R,1+oversampling*(Nx-1))[None,:]
    fyi=np.linspace(-fmax*Dy/R,fmax*Dy/R,1+oversampling*(Ny-1))[None,:]
    fxm=np.dot(fxi.T,np.ones((1,fyi.shape[1])))
    fym=np.dot(np.ones((fxi.shape[1],1)),fyi)

    # 2D low-pass filter
    hlp2=np.ones((oversampling,oversampling))/(oversampling)**2
    P1=np.zeros((len(fx),len(fy),len(f)),dtype=np.complex128)
    P1i0=np.zeros((len(fxi.T),len(fyi.T)),dtype=np.complex128)

    # Loop over frequencies
    for k in range(len(f)):
      fxmax=0.5*Dx/R*f[k]
      fymax=0.5*Dy/R*f[k]
      Ix=np.where(np.abs(fxi)<fxmax)[1]
      Iy=np.where(np.abs(fyi)<fymax)[1]
      wx=tukeywin(len(Ix),betax)
      wy=tukeywin(len(Iy),betay)
      w=np.dot(wx[:,None],wy[None,:])
      if circsymm:
        w=tukey2(len(Ix),betax)
      P1k=P1i0
      P1k[Ix[0]:Ix[-1]+1,Iy[0]:Iy[-1]+1] = w*P[k]
      if oversampling>1:
        # Low-pass filter
        P1k=convolve2d(P1k,hlp2,'same')
      if ROCcorrection:
        Ha=apertureCorrection(R,Rxy[0],Rxy[1],f[k],fxi,fyi)
        P1k=P1k*Ha
      if sphericalCorrection:
        Hs=focalcorrection(R,f[k],fxi,fyi)
        P1u=P1k
        P1k=convolve2d(P1k,Hs,'same')
      P1k=P1k[0::oversampling,0::oversampling] # Decimation to correct sampling rate
      P1[:,:,k] = P1k
    return (P1,fx,fy,f)

  @staticmethod
  def FourierToSpaceTimeNorway(P12,fx1,fy1,f1,**kwargs):

    opt = dotdict({'resolution' : np.r_[6e-5, 6e-5, 6e-5],
                   'c'          : 1540.0})
    opt.update(**kwargs)

    resolution = opt.resolution
    c          = opt.c

    testPlot = False

    dfx=(fx1[1]-fx1[0]);
    dfy=(fy1[1]-fy1[0]);
    df=(f1[1]-f1[0]);

    Nfx=int(2*np.round(c/dfx/resolution[0]/2)+1)
    Nfy=int(2*np.round(c/dfy/resolution[1]/2)+1)
    Nfs=int(2*np.round(c/2/df/resolution[2]/2)+1)

    xax=np.linspace(-1/dfx*c/2,1/dfx*c/2,Nfx);
    yax=np.linspace(-1/dfy*c/2,1/dfy*c/2,Nfy);
    zax=np.linspace(-1/df*c/4,1/df*c/4,Nfs);

    [Nx,Ny,Nf]=P12.shape
    P12zp=np.zeros((Nfx,Nfy,Nfs),dtype=np.complex128)
    ix1=int(np.round(1+(Nfx-1)/2-(Nx-1)/2))
    iy1=int(np.round(1+(Nfy-1)/2-(Ny-1)/2))
    if1=int(np.round(1+(Nfs-1)/2+1+f1[0]/df) - 1)
    P12zp[ix1:ix1+Nx,iy1:iy1+Ny,if1:if1+Nf]=P12;
    P12zp=np.fft.fftshift(P12zp);
    p12=ifftn(P12zp);
    p12=np.fft.fftshift(p12);

    return (p12,xax,yax,zax)


  def FieldAtFocus(xdc, **kwargs):
    opt = dotdict({'nx' : 100,
                   'ny' : 100,
                   'dx' : 6e-5,
                   'dy' : 6e-5})
    # Maximum frequency
    fmax = None
    if xdc.bandwidth == 0.0:
      fmax = xdc.f0

class OurClass:

    def __init__(self, a):
        self.OurAtt = a

    @property
    def OurAtt(self):
        return self.__OurAtt

    @OurAtt.setter
    def OurAtt(self, val):
        if val < 0:
            self.__OurAtt = 0
        elif val > 1000:
            self.__OurAtt = 1000
        else:
            self.__OurAtt = val


x = OurClass(10)
print(x.OurAtt)


# Local variables: #
# tab-width: 2 #
# python-indent: 2 #
# indent-tabs-mode: nil #
# End: #
