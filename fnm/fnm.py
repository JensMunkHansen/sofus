import numpy as np
import sys

from dicts import dotdict
from euler import euler2rot

from timeit import default_timer as timer

# Test integrating x**2 from [-1,1]

def create_time_string(seconds):
  m, s = divmod(seconds, 60)
  h, m = divmod(m, 60)
  ms = np.round(s * 1000) % 1000
  timeString = "%d:%02d:%02d,%03d" % (h, m, s, (1000*ms))
  return timeString;


class rect(dotdict):
    def __init__(self,*args, **kwargs):
        """
        Initialize a rectangular element

        hRect = rect(args,*,hw=1.0, hh=1.0, center=[0,0,0],euler=[0,0,0],conv='yxz',intrinsic=True,nAbcissa=[1,1])

        A rect object is returned, which can be used for computing a continuous wave field (ARFI pressure)

        hw        - half width of the element
        hh        - half height of the element
        center    - center of element
        euler     - orientation of element
        conv      - convention used for the orientation, order of rotations
        intrinsic - if true, the coordinate system is rotating
        nAbcissa  - number of Gaussian abcissas used for integration [nWidth,nHeight]
        """
        opt = dotdict({'hw' : 1,
                       'hh' : 1,
                       'center'    : [0,0,0],
                       'euler'     : [0,0,0],
                       'conv'      : 'yxz',
                       'intrinsic' : True,
                       'nAbcissa'  : 1})
        opt.update(*args,**kwargs)
        opt.nAbcissa = np.array(opt.nAbcissa).flatten()
        
        super(rect,self).update(**opt)

        if len(self.nAbcissa)==1:
            self.nAbcissa = np.r_[self.nAbcissa[0], self.nAbcissa[0]]

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

    def Hz_ref(self,s,l,z,k):
      """
      Field above a corner of an element with dimensions s and l (just a function)
      
      H_{s,l}(z;k) =
       - i/(2\pi k) ( s \int_0^l (exp(-ik\sqrt{\sigma^2+z^2+s^2}) - exp(-ikz)) / (\sigma^2+s^2) d\sigma +
                      l \int_0^s (exp(-ik\sqrt{\sigma^2+z^2+l^2}) - exp(-ikz)) / (\sigma^2+l^2) d\sigma)

      Reference implementation
      """

      ndiv = self.nAbcissa[0]
      xs, ws = np.polynomial.legendre.leggauss(ndiv)

      ndiv = self.nAbcissa[1]
      ys, ws1 = np.polynomial.legendre.leggauss(ndiv)

      ls = l/2.0*xs + l/2.0
      ss = s/2.0*ys + s/2.0

      result = -1j / (2*np.pi*k) * (
          l/2.0 * s * np.sum(((np.exp(-k*np.sqrt(ls**2+z**2+s**2)*1j)-np.ones(ndiv)*np.exp(-k*z*1j))/(ls**2+s**2))*ws) +
          s/2.0 * l * np.sum(((np.exp(-k*np.sqrt(ss**2+z**2+l**2)*1j)-np.ones(ndiv)*np.exp(-k*z*1j))/(ss**2+l**2))*ws1))
      return result

    def H_single(self,u1,u2,v,z,k,us,ws):
      ndiv = len(us)
      ss = (u2-u1)/2.0 * us + (u2+u1)/2.0
      ss2 = ss**2
      return (u2-u1)/2.0 * np.sum(((np.exp(-k*np.sqrt(z**2+ss2+v**2)*1j)-np.ones(ndiv)*np.exp(-k*z*1j))/(ss2+v**2))*ws)

    def H_ref(self,point,k):
      """ 
      Field at an arbitrary position 'point' of the continuous
      wave with wave-number k. Superposition is used as first
      demonstrated by Lockwood and Wilette (1973)

      Reference implementation
      """
      vertices = self.corners()

      # Changed sign of roll
      edges = vertices - np.roll(vertices,-1,axis=0)
      edges = edges / np.r_[3*[np.sqrt(np.sum(edges**2,axis=1))]].transpose()

      w = np.cross(edges[0],edges[1])
      ref = self.center
      ref2point = point - ref
      z = np.dot(w,ref2point)
      
      x = np.dot(edges[0],ref2point)
      y = np.dot(edges[1],ref2point)

      s = np.zeros(4)
      l = np.zeros(4)
      a = self.hw
      b = self.hh

      s[0] = np.abs(x)+a
      s[1] = a - np.abs(x)
      s[2] = np.abs(x)+a
      s[3] = a - np.abs(x)

      l[0] = np.abs(y)+b
      l[1] = np.abs(y)+b
      l[2] = b - np.abs(y)
      l[3] = b - np.abs(y)

      result = 0j

      # Four contributions are computed and superposed to give the right value
      for i in range(4):
          result1 = self.Hz_ref(np.abs(s[i]),np.abs(l[i]),z,k)*np.sign(s[i])*np.sign(l[i])
          result = result + result1

      return result

    def H9(self,xs,ys,zs,k):
      """
      CW field of element on a grid specified using nx,dx,nz and the wave-number k. 
      
      Reference implementation
      """

      nx = xs.shape[0]
      nz = zs.shape[1]
      results = np.zeros((nx,nz),dtype=np.complex64)

      ref = self.center
      vertices = self.corners()

      edges = vertices - np.roll(vertices,-1,axis=0)
      edges = edges / np.r_[3*[np.sqrt(np.sum(edges**2,axis=1))]].transpose()
      w = np.cross(edges[0],edges[1])

      s = np.zeros(4)
      l = np.zeros(4)
      a = self.hw
      b = self.hh

      us, ws = np.polynomial.legendre.leggauss(ndiv)
      
      for ix in range(nx):
        for iz in range(nz):
          point = [xs[ix,0],0,zs[0,iz]]
          ref2point = point - ref
          x = np.dot(edges[0],ref2point)
          y = np.dot(edges[1],ref2point)
          z = np.dot(w,ref2point)

          s[0] = s[2] = np.abs(x)+a
          s[1] = s[3] = a - np.abs(x) # reverse and combine with s[0]

          l[0] = l[1] = np.abs(y)+b
          l[2] = l[3] = b - np.abs(y)
          ndiv = self.nAbcissa[0]

          result1 = 0
          if (np.abs(x) > a):
            _l = l[0]
            # 2 integrals
            result1 = result1 + _l * self.H_single(np.abs(x),np.abs(x)+a,_l,z,k,us,ws)
            result1 = result1 + _l * self.H_single(np.abs(x)-a,np.abs(x),_l,z,k,us,ws)
            _l = np.abs(l[2]) # Check if sign is needed
            # 2 integrals
            result1 = result1 + _l * self.H_single(np.abs(x),np.abs(x)+a,_l,z,k,us,ws)
            result1 = result1 + _l * self.H_single(np.abs(x)-a,np.abs(x),_l,z,k,us,ws)
          else:
            _s = s[0]
            _l = l[0]
            result1 = result1 + _l * self.H_single(0,_s,_l,z,k,us,ws)
            _s = s[1]
            _l = l[1]
            result1 = result1 + _l * self.H_single(0,_s,_l,z,k,us,ws)
            _s = s[2]
            _l = l[2] # Gives sign
            result1 = result1 + _l * self.H_single(0,_s,_l,z,k,us,ws)
            _s = s[3]
            _l = l[3] # Gives sign
            result1 = result1 + _l * self.H_single(0,_s,_l,z,k,us,ws)
          if (np.abs(y) > b):
            _s = s[0]
            # 2 integrals, 0 and 2
            result1 = result1 + _s * self.H_single(np.abs(y)-b,np.abs(y),_s,z,k,us,ws)
            result1 = result1 + _s * self.H_single(np.abs(y),np.abs(y)+b,_s,z,k,us,ws)
            _s = np.abs(s[1]) # Is sign needed?
            # 2 integrals, 1 and 3
            result1 = result1 + _s * self.H_single(np.abs(y)-b,np.abs(y),_s,z,k,us,ws)
            result1 = result1 + _s * self.H_single(np.abs(y),np.abs(y)+b,_s,z,k,us,ws)
          else:
            _s = s[0]
            _l = l[0]
            result1 = result1 + _s*self.H_single(0,_l,_s,z,k,us,ws)
            _s = s[1] # Gives sign
            _l = l[1]
            result1 = result1 + _s*self.H_single(0,_l,_s,z,k,us,ws)
            _s = s[2]
            _l = l[2]
            result1 = result1 + _s*self.H_single(0,_l,_s,z,k,us,ws)
            _s = s[3] # Gives sign
            _l = l[3]
            result1 = result1 + _s*self.H_single(0,_l,_s,z,k,us,ws)
          results[ix,iz] = result1

      return results / (2*np.pi * k * 1j)

    def H4(self,xs,ys,zs,k):
      return self.H_accurate(xs,ys,zs,k)

    def H_accurate(self,xs,ys,zs,k):
      """
      CW field of element on a grid specified using nx,dx,nz and the wave-number k. 
      
      Reference implementation
      """

      nx = xs.shape[0]
      nz = zs.shape[1]
      results = np.zeros((nx,nz),dtype=np.complex64)

      ref = self.center
      vertices = self.corners()

      edges = vertices - np.roll(vertices,-1,axis=0)
      edges = edges / np.r_[3*[np.sqrt(np.sum(edges**2,axis=1))]].transpose()
      w = np.cross(edges[0],edges[1])

      s = np.zeros(4)
      l = np.zeros(4)
      a = self.hw
      b = self.hh

      ndiv0 = self.nAbcissa[0]
      us0, ws0 = np.polynomial.legendre.leggauss(ndiv0)

      ndiv1 = self.nAbcissa[1]
      us1, ws1 = np.polynomial.legendre.leggauss(ndiv1)
      
      for ix in range(nx):
        for iz in range(nz):
          point = [xs[ix,0],ys[ix,iz],zs[0,iz]]
          ref2point = point - ref
          x = np.dot(edges[0],ref2point)
          y = np.dot(edges[1],ref2point)
          z = np.dot(w,ref2point)

          s[0] = s[2] = np.abs(x)+a
          s[1] = s[3] = a - np.abs(x) # reverse and combine with s[0]

          l[0] = l[1] = np.abs(y)+b
          l[2] = l[3] = b - np.abs(y)

          result1 = 0
          if (np.abs(x) > a):
            _l = l[0]
            # 2 integrals
            result1 = result1 + _l * self.H_single(np.abs(x),np.abs(x)+a,_l,z,k,us1,ws1)
            result1 = result1 + _l * self.H_single(np.abs(x)-a,np.abs(x),_l,z,k,us1,ws1)
            _l = l[2]#np.abs(l[2]) # Sign is not needed (verify)
            # 2 integrals
            result1 = result1 + _l * self.H_single(np.abs(x),np.abs(x)+a,_l,z,k,us1,ws1)
            result1 = result1 + _l * self.H_single(np.abs(x)-a,np.abs(x),_l,z,k,us1,ws1)
          else:
            _s = s[0]
            _l = l[0]
            result1 = result1 + _l * self.H_single(0,_s,_l,z,k,us1,ws1)
            _s = s[1]
            _l = l[1]
            result1 = result1 + _l * self.H_single(0,_s,_l,z,k,us1,ws1)
            _s = s[2]
            _l = l[2] # Gives sign
            result1 = result1 + _l * self.H_single(0,_s,_l,z,k,us1,ws1)
            _s = s[3]
            _l = l[3] # Gives sign
            result1 = result1 + _l * self.H_single(0,_s,_l,z,k,us1,ws1)
          if (np.abs(y) > b):
            _s = s[0]
            # 2 integrals, 0 and 2
            result1 = result1 + _s * self.H_single(np.abs(y)-b,np.abs(y),_s,z,k,us0,ws0)
            result1 = result1 + _s * self.H_single(np.abs(y),np.abs(y)+b,_s,z,k,us0,ws0)
            _s = s[1]#np.abs(s[1]) # Sign is not needed (verify)
            # print('s[1], abs(s[1]): %f %f' % (s[1],np.abs(s[1])))
            # 2 integrals, 1 and 3
            result1 = result1 + _s * self.H_single(np.abs(y)-b,np.abs(y),_s,z,k,us0,ws0)
            result1 = result1 + _s * self.H_single(np.abs(y),np.abs(y)+b,_s,z,k,us0,ws0)
          else:
            _s = s[0]
            _l = l[0]
            result1 = result1 + _s*self.H_single(0,_l,_s,z,k,us0,ws0)
            _s = s[1] # Gives sign
            _l = l[1]
            result1 = result1 + _s*self.H_single(0,_l,_s,z,k,us0,ws0)
            _s = s[2]
            _l = l[2]
            result1 = result1 + _s*self.H_single(0,_l,_s,z,k,us0,ws0)
            _s = s[3] # Gives sign
            _l = l[3]
            result1 = result1 + _s*self.H_single(0,_l,_s,z,k,us0,ws0)
          results[ix,iz] = result1

      return results / (2*np.pi * k * 1j)

    def HN(self,xs,ys,zs,k):
      """
      Fast version of H_accurate. For performance reasons, it is preferred that n1 is larger than n0
      """
      vertices = self.corners()
      edges = vertices - np.roll(vertices,-1,axis=0)
      edges = edges / np.r_[3*[np.sqrt(np.sum(edges**2,axis=1))]].transpose()
      
      w = np.cross(edges[0],edges[1])
      
      ref = self.center
      
      # Do not say anything about how many different x- or z-coordinates we have
      n0 = xs.shape[0]
      n1 = zs.shape[1]
      
      results = np.zeros((n0,n1),dtype=np.complex64)
      
      a = self.hw
      b = self.hh

      # Projection
      xs1 = np.abs((xs - ref[0]*np.ones(xs.shape))*edges[0,0] +
                   (ys - ref[1]*np.ones(ys.shape))*edges[0,1] +
                   (zs - ref[2]*np.ones(zs.shape))*edges[0,2])
      ys1 = np.abs((xs - ref[0]*np.ones(xs.shape))*edges[1,0] +
                   (ys - ref[1]*np.ones(zs.shape))*edges[1,1] +
                   (zs - ref[2]*np.ones(zs.shape))*edges[1,2])

      zs1 = (xs - ref[0]*np.ones(xs.shape))*w[0] + \
            (ys - ref[1]*np.ones(xs.shape))*w[1] + \
            (zs - ref[2]*np.ones(xs.shape))*w[2] 

      del ref,edges,vertices,w
      
      ndiv0   = self.nAbcissa[0]
      ndiv1   = self.nAbcissa[1]

      us0, ws0 = np.polynomial.legendre.leggauss(ndiv0)
      us0 = us0.reshape((ndiv0,1))
      
      us1, ws1 = np.polynomial.legendre.leggauss(ndiv1)
      us1 = us1.reshape((ndiv1,1))

      ws0 = np.r_[n1*[ws0]].transpose()
      ws1 = np.r_[n1*[ws1]].transpose()

      # Needed when used as not the integrand
      l0 = np.abs(ys1)+b
      l2 = b - np.abs(ys1)
      s0 = np.abs(xs1)+a
      s1 = a - np.abs(xs1)
      
      # s-integral
      mask  = np.abs(xs1) > a
      slow0  = np.zeros(xs1.shape)
      slow1  = np.zeros(xs1.shape)
      shigh0 = s0.copy()
      shigh1 = s1.copy()

      tmp = np.abs(xs1)-a
      slow0[mask]  = tmp[mask]
      tmp = np.abs(xs1)
      slow1[mask]  = tmp[mask]
      shigh0[mask] = tmp[mask]
      tmp = s0.copy()
      shigh1[mask] = tmp[mask]

      # l-integral
      mask = np.abs(ys1) > b
      llow0 = np.zeros(ys1.shape)
      llow1 = np.zeros(ys1.shape)
      lhigh0 = l0.copy()
      lhigh1 = l2.copy()

      tmp = np.abs(ys1)-b
      llow0[mask]  = tmp[mask]
      tmp = np.abs(ys1)
      llow1[mask]  = tmp[mask]
      lhigh0[mask] = tmp[mask]
      tmp = l0.copy()
      lhigh1[mask] = tmp[mask]

      del tmp, mask
      
      for i0 in range(n0):
        z2 = (zs1[i0]**2).reshape((1,n1))
        expz = np.exp(-k*zs1[i0]*1j)
        expz0 = np.r_[ndiv0 * [expz]]
        expz1 = np.r_[ndiv1 * [expz]]

        sscale0 = (shigh0[i0] - slow0[i0])/2.0
        ss02 = (np.dot(us1, sscale0.reshape((1,n1))) + np.r_[ndiv1 * [(shigh0[i0] + slow0[i0])/2.0]])**2

        sscale1 = (shigh1[i0] - slow1[i0])/2.0
        ss12 = (np.dot(us1, sscale1.reshape((1,n1))) + np.r_[ndiv1 * [(shigh1[i0] + slow1[i0])/2.0]])**2
        
        l = l0[i0]
        l22 = l**2
        
        num0 = np.exp(-k*np.sqrt(ss02 + z2 + l22)*1j) - expz1
        num1 = np.exp(-k*np.sqrt(ss12 + z2 + l22)*1j) - expz1

        if False:#i0==0:
          print('sscale0.shape'+str(sscale0.shape)) # (n1,)
          print('l.shape'+str(l.shape))             # (n1,)
          print('num0.shape'+str(num0.shape))       # (ndiv1,n1)
          print('ss02.shape'+str(ss02.shape))       # (ndiv1,n1)
          print('l22.shape'+str(l22.shape))         # (n1)
          print('ws1.shape'+str(ws1.shape))         # (ndiv1,n1), n1 = 250
          
        tmp = np.sum((l * ( sscale0 * (num0 / (ss02 + l22)) +
                            sscale1 * (num1 / (ss12 + l22))) * ws1), axis=0)

        l = l2[i0]
        l22 = l**2
        
        num0 = np.exp(-k*np.sqrt(ss02 + z2 + l22)*1j) - expz1
        num1 = np.exp(-k*np.sqrt(ss12 + z2 + l22)*1j) - expz1

        tmp = tmp + np.sum((l * ( sscale0 * (num0 / (ss02 + l22)) +
                                  sscale1 * (num1 / (ss12 + l22))) * ws1), axis=0)

        # l-integral
        
        lscale0 = (lhigh0[i0] - llow0[i0])/2.0
        ls02 = (np.dot(us0, lscale0.reshape((1,n1))) + np.r_[ndiv0 * [(lhigh0[i0] + llow0[i0])/2.0]])**2

        lscale1 = (lhigh1[i0] - llow1[i0])/2.0
        ls12 = (np.dot(us0, lscale1.reshape((1,n1))) + np.r_[ndiv0 * [(lhigh1[i0] + llow1[i0])/2.0]])**2
        
        s = s0[i0]
        s22 = s**2
        
        num0 = np.exp(-k*np.sqrt(ls02 + z2 + s22)*1j) - expz0
        num1 = np.exp(-k*np.sqrt(ls12 + z2 + s22)*1j) - expz0

        tmp = tmp + np.sum((s * ( lscale0 * (num0 / (ls02 + s22)) +
                                  lscale1 * (num1 / (ls12 + s22))) * ws0), axis=0)
        
        s = s1[i0]
        s22 = s**2

        num0 = np.exp(-k*np.sqrt(ls02 + z2 + s22)*1j) - expz0
        num1 = np.exp(-k*np.sqrt(ls12 + z2 + s22)*1j) - expz0

        results[i0] = (-1j / (2*np.pi*k)) * \
                      (tmp + np.sum((s * ( lscale0 * (num0 / (ls02 + s22)) +
                                           lscale1 * (num1 / (ls12 + s22))) * ws0), axis=0))
        
      return results
    
      
    def H(self,xs,ys,zs,k):
      """
      Compute the impulse response at positions specified using 3
      matrices of x-,y-, and z-coordinates for a continuous wave
      with wave-number k
      
      This is a fast implementation, but less accurate than H_accurate or HN
      """
      vertices = self.corners()
      
      # Changed sign of roll
      edges = vertices - np.roll(vertices,-1,axis=0)
      edges = edges / np.r_[3*[np.sqrt(np.sum(edges**2,axis=1))]].transpose()
      
      w = np.cross(edges[0],edges[1])
      
      ref = self.center
      
      # Doesn't say anything about how many different x- or z-coordinates we have
      nx = xs.shape[0]
      nz = zs.shape[1]
      
      results = np.zeros((nx,nz),dtype=np.complex64)
      
      s = np.zeros((4,nx))
      l = np.zeros((4,nz))
      a = self.hw
      b = self.hh
      
      xs1 = np.abs((xs - ref[0]*np.ones(xs.shape))*edges[0,0] +
                   (ys - ref[1]*np.ones(ys.shape))*edges[0,1] +
                   (zs - ref[2]*np.ones(zs.shape))*edges[0,2])
      ys1 = np.abs((xs - ref[0]*np.ones(xs.shape))*edges[1,0] +
                   (ys - ref[1]*np.ones(zs.shape))*edges[1,1] +
                   (zs - ref[2]*np.ones(zs.shape))*edges[1,2])
      
      s0_s2 = xs1 + a
      s1_s3 = a - xs1
      
      l0_l1 = ys1 + b
      l2_l3 = b - ys1
      
      ndiv   = self.nAbcissa[0]
      ndiv_1 = self.nAbcissa[1]
      
      xsi, ws = np.polynomial.legendre.leggauss(ndiv)
      xsi = xsi.reshape((ndiv,1))
      
      xsi_1, ws_1 = np.polynomial.legendre.leggauss(ndiv_1)
      xsi_1 = xsi_1.reshape((ndiv_1,1))
      
      zs1 = (xs - ref[0]*np.ones(xs.shape))*w[0] + (zs - ref[2]*np.ones(zs.shape))*w[2] 
      ws = np.r_[nz*[ws]].transpose()
      ws_1 = np.r_[nz*[ws_1]].transpose()
      
      for ix in range(nx):
        z = zs1[ix].reshape((1,nz))
        z2 = z**2
        con1   = np.dot(np.ones((ndiv,1)),np.exp(-k*z*1j))
        con1_1 = np.dot(np.ones((ndiv_1,1)),np.exp(-k*z*1j))
      
        # l0, s0
        l0_l1_ix = l0_l1[ix]
        s0_s2_ix = s0_s2[ix]
      
        ls01_2 = np.dot(1+xsi_1,(l0_l1_ix/2.0).reshape((1,nz)))
        ls01_2 = ls01_2**2
        ss02_2 = np.dot(1+xsi,(s0_s2_ix/2.0).reshape((1,nz)))
        ss02_2 = ss02_2**2
      
        s0_s2_2 = s0_s2_ix**2
        l0_l1_2 = l0_l1_ix**2
        
        bum  = np.exp(-k*np.sqrt(ls01_2 + z2 + s0_s2_2)*1j) - con1_1 
        bum1 = np.exp(-k*np.sqrt(ss02_2 + z2 + l0_l1_2)*1j) - con1
        result2 = -1j / (2*np.pi*k) * (np.sum((
            l0_l1_ix/2.0 * s0_s2_ix * (bum / (ls01_2 + s0_s2_2))) * ws_1,axis=0) +
            np.sum((s0_s2_ix/2.0 * l0_l1_ix * (bum1/  (ss02_2 + l0_l1_2))) * ws,axis=0))
        results[ix] = results[ix] + result2
        
        # l0, s1
        s1_s3_ix = s1_s3[ix]
      
        ss13_2 = np.dot(1+xsi,(s1_s3_ix/2.0).reshape((1,nz)))
        ss13_2 = ss13_2**2
      
        s1_s3_2 = s1_s3_ix**2
        
        bum  = np.exp(-k*np.sqrt(ls01_2 + z2 + s1_s3_2)*1j) - con1_1 
        bum1 = np.exp(-k*np.sqrt(ss13_2 + z2 + l0_l1_2)*1j) - con1
        result2 = -1j / (2*np.pi*k) * (np.sum((
            l0_l1_ix/2 * s1_s3_ix * (bum / (ls01_2 + s1_s3_2)))*ws_1,axis=0) +
            np.sum(( s1_s3_ix/2 * l0_l1_ix * (bum1/  (ss13_2 + l0_l1_2)))* ws,axis=0))
        results[ix] = results[ix] + result2
      
        # s0, l1
        l2_l3_ix = l2_l3[ix]
        
        l23_2 = np.dot(1+xsi_1,(l2_l3_ix/2.0).reshape((1,nz)))
        l23_2 = l23_2**2 
      
        l2_l3_2 = l2_l3_ix**2
      
        bum  = np.exp(-k*np.sqrt(l23_2 + z2 + s0_s2_2)*1j) - con1_1 
        bum1 = np.exp(-k*np.sqrt(ss02_2 + z2 + l2_l3_2)*1j) - con1
        result2 = -1j / (2*np.pi*k) * (np.sum((
            l2_l3_ix/2.0 * s0_s2_ix * (bum / (l23_2 + s0_s2_2)))*ws_1,axis=0) +
            np.sum(( s0_s2_ix/2 * l2_l3_ix * (bum1/  (ss02_2 + l2_l3_2))) * ws,axis=0))
        results[ix] = results[ix] + result2
        
        # s1, l1
        bum  = np.exp(-k*np.sqrt(l23_2 + z2 + s1_s3_2)*1j) - con1_1 
        bum1 = np.exp(-k*np.sqrt(ss13_2 + z2 + l2_l3_2)*1j) - con1
      
        result2 = -1j / (2*np.pi*k) * (np.sum((
            l2_l3_ix/2.0 * s1_s3_ix * (bum / (l23_2 + s1_s3_2)))*ws_1,axis=0)+
            np.sum((s1_s3_ix/2.0 * l2_l3_ix * (bum1/  (ss13_2 + l2_l3_2)))*ws,axis=0))
      
        results[ix] = results[ix] + result2
      
      return results
    

class linear_array(dotdict):
    def __init__(self,*args, **kwargs):

      opt = dotdict({'nElements' : 192,
                     'nSubH'     : 1,
                     'nAbcissa'  : [1,1],
                     'pitch'     : 0.2e-3,
                     'kerf'      : 0.2e-4,
                     'height'    : 1.0e-2,
                     'c'         : 1540.0,
                     'focus'     : None})

      opt.update(**kwargs)

      half_width  = (opt.pitch - opt.kerf) / 2.0
      half_height = opt.height / 2.0
      self.rects = []
      self.nElements = opt.nElements
      self.nSubH     = opt.nSubH

      assert(self.nSubH % 2 == 1)
      
      focus = 1.0
      R     = 1.0

      if (opt.focus != None):
          focus = opt.focus
          elSector = 2 * np.arctan2(half_height, focus)
          R = np.sqrt(focus**2 + (opt.height / 2.0)**2)
      else:
          elSector = 0
      
      dEl =  elSector / max((opt.nSubH-1),1)
      elAngles = (np.r_[0:opt.nSubH] - (opt.nSubH - 1.0)/2) * dEl
      
      if (opt.focus != None):
          chordLength = 2 * focus * np.sin(dEl/2)
      else:
          chordLength = opt.height/opt.nSubH
      
      for iElement in range(opt.nElements):
          for iSubH in range(opt.nSubH):
              if opt.focus != None:
                  center = np.r_[(iElement - (opt.nElements-1.0)/2)*opt.pitch,
                                 focus * np.tan(elAngles[iSubH]),
                                 -(R * np.cos(elAngles[iSubH]) - focus)]
              else:
                  center = np.r_[(iElement - (opt.nElements-1.0)/2)*opt.pitch,
                                 (iSubH - (opt.nSubH - 1.0)/2.0)*chordLength,
                                 0]
              for iSubW in range(1):
                  center1 = center + 2 * half_width/1.0 * np.r_[1,0,0] * (iSubW - (1-1)/2.0)
                  r = rect(hw=half_width/1.0,hh=chordLength/2.0,
                           center=center1,
                           euler=[0,elAngles[iSubH],0],
                           conv='yxz',
                           intrinsic=True,
                           nAbcissa=opt.nAbcissa)
                  self.rects.append(r)

      self.phases = np.ones(opt.nElements,dtype=np.complex128)
      self.c = opt.c

    def set_focus(self,focus,f0):
      """
      Set electronic focus
      """
      focus = np.array(focus).flatten()

      k = 2*np.pi / (self.c / f0)
      for i in range(self.nElements):
        self.phases[i] = np.exp(-1j * np.angle(self.rects[i].HN(np.r_[focus[0]].reshape((1,1)),
                                                                np.r_[focus[1]].reshape((1,1)),
                                                                np.r_[focus[2]].reshape((1,1)),k)[0,0]))
            
    def cw_pressure(self,xs,ys,zs,k):
      """
      The CW pressure is the absolute result of this function multiplied by
      
      -i \omega * \rho * v_n 

      where \omega is the angular frequency, \rho is the density and v_n is normal component of the velocity.
      """
      start = timer()
      nSubElements = len(self.rects)
      #      print('Element: %d' % 0)  
      result = self.phases[0] * self.rects[0].HN(xs,ys,zs,k)
      for i in range(1,nSubElements):
        # print('Element: %d' % (i))
        result = result + self.phases[i / self.nSubH] * self.rects[i].HN(xs,ys,zs,k)
      end = timer()
      timeString = create_time_string(end-start)
      print(timeString)
      return result
        
