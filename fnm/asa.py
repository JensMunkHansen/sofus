from dicts import dotdict
import numpy as np

# TODO: Implement propagation of convolution of fund pressure

def asa(y0,dx,dy,zs,**kwargs):
    """
    TODO: Add time-like propagators
    """
    from scipy.fftpack import (fft2,ifft2)
    opt = dotdict({'c'          : 1540.0,
                   'f0'         : 1e6,
                   'beta'       : 0.5, # dB/(cm*MHz)
                   'Nx'         : None,
                   'Ny'         : None,
                   'N'          : 128,
                   'propagator' : 'P', # 'p', 'P', 'v' or 'V'
                   'restrict'   : True,
               })

    opt.update(**kwargs)

    dtype = y0.dtype
    
    dBperNeper = 20 * np.log10(np.e)
    eps = np.finfo(np.float64).resolution

    if (opt.Nx == None):
        opt.Nx = opt.N
    if (opt.Ny == None):
        opt.Ny = opt.N

    Nx = opt.Nx
    Ny = opt.Ny

    lambd = opt.c / opt.f0
    
    alpha = opt.beta / dBperNeper * 100 * opt.f0 / 1e6
    
    [nx,ny] = y0.shape

    assert((Nx >= nx) and (Ny >= ny))

    if ((dx >= lambd/2) or (dy >= lambd/2)):
        print('Input is sampled coarser than lambda/2')
            
    Y0 = fft2(y0,(Nx,Ny))
    
    nz = len(zs)
    zs = zs - zs[0]
    dzs = (np.roll(zs,-1) - zs)[:-1]

    if np.any(dzs >= lambd/2):
        print('Radial steps are larger than lambda/2')
    
    output = np.zeros((nz,nx,ny),dtype=dtype)
    output[0,:,:] = y0

    k = 2*np.pi/lambd

    # Note asymmetric, also for N odd
    kx = np.r_[-Nx//2:Nx//2:1]*lambd/(Nx*dx);
    ky = np.r_[-Ny//2:Ny//2:1]*lambd/(Ny*dy);
        
    [kxspace,kyspace] = np.meshgrid(kx,ky,indexing='ij')
    kxsq_ysq = np.fft.fftshift(kxspace**2 + kyspace**2).astype(dtype)
    kzspace = k * np.sqrt(1.0 - kxsq_ysq)
    
    for iz in range(nz-1):
      # Distance (relative to origing)
      z     = zs[iz+1]

      # Delta increment in z (should it be dx and dy?)
      delta = dzs[iz]

      # Basic spectral propagator
      if opt.propagator[0]=='P':
        if zs[iz+1] > zs[0]:
            # Forward
            H = np.conj(np.exp(1j * z * kzspace)); # exp (i*(z-z0)*K(kx,ky))
        else:
            # Backward
            H = np.exp(-1j*z * kzspace) * kxsq_ysq <= 1;
      elif opt.propagator[0]=='V':
        H = k*np.conj(np.exp(1j*kzspace*z)) / (1j*kzspace)
        H[np.isnan(H)] = eps

      # Attenuation
      if opt.beta > 0:
        evans_mode = np.sqrt(kxsq_ysq)<1.0
        # Needs to be complex
        H = H * np.exp(- alpha * z / np.cos(np.arcsin(np.sqrt(kxsq_ysq))) * evans_mode) * evans_mode

      # Restricted angular content
      if opt.restrict:
        # For a square field area
        # D = (N-1)*delta
        # thres = np.sqrt(0.5*D**2/(0.5*D**2+z**2))
        Dx = (Nx-1)*dx
        Dy = (Ny-1)*dy
        D2xy = 0.25*Dx**2 + 0.25*Dy**2
        thres = np.sqrt(D2xy/(D2xy+z**2))
        H = H * (np.sqrt(kxsq_ysq) <= thres)
      newpress = ifft2(Y0*H,(Nx,Ny))
      output[iz+1,:,:] = newpress[0:nx,0:ny]
    return output
