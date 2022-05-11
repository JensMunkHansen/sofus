# -*- coding: utf-8 -*-
"""Python reference implementation for computation of arc length of sphere
intersecting a rectange

Example
-------
TODO: How to run or use this

    $ python example_numpy.py

Notes
-----
    No notes for this file

"""
import os, sys
import numpy as np
from copy import deepcopy

filedir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(filedir, '../../fnm'))

from dicts import dotdict
from euler import euler2rot

debug = False

class rect(dotdict):
  def __init__(self,*args, **kwargs):
    opt = dotdict({'hw' : 1.0,
                   'hh' : 1.0,
                   'center'    : [0,0,0],
                   'euler'     : [0,0,0],
                   'conv'      : 'yxy',
                   'ndiv'      : 18})

    opt.update(*args,**kwargs)
    super(rect,self).update(**opt)

    self.center = np.array(self.center).flatten()
    self.euler  = np.array(self.euler).flatten()

    # This introduces errors. Instead use basisvectors
    rotm = euler2rot(self.euler[0], self.euler[1], self.euler[2],
                     conv = 'yxy',
                     intrinsic = True)

    # Vertices ordered clock-wise
    a = np.array([[1,1],[-1,1],[-1,-1],[1,-1]])
    vertices = np.c_[a * np.r_[self.hw,self.hh],[0,0,0,0]]
    self.vertices = np.r_[[np.dot(rotm,vertices[i]) + self.center for i in range(4)]]

    edges = self.vertices - np.roll(self.vertices,-1,axis=0)
    edges = edges / np.r_[3*[np.sqrt(np.sum(edges**2,axis=1))]].transpose()
    self.normal = np.cross(edges[0],edges[1])
    self.uvector = edges[0]
    self.vvector = edges[1]

  def element(self):
    """Return variables defining the element as used in SOFUS

    """
    return np.r_[np.r_[self.hw], np.r_[self.hh], self.center, self.euler]

  def local(self, point):
    """Projection of point onto local rectangle coordinates. Projection,
    distances to vertices and distance to the center of the element
    are returned

    """
    point = np.array(point).flatten()
    r2p = point - self.center

    dist = np.sqrt(np.sum(r2p**2))
    u = np.dot(self.uvector,r2p)
    v = np.dot(self.vvector,r2p)
    w = np.abs(np.dot(self.normal,r2p))

    vdists = np.sqrt(np.sum((self.vertices - np.r_[4*[point]])**2,axis=1))

    return (np.r_[u,v,w], vdists, dist)

  def direct(self, proj):
    """
    Direct term using Gauss-Legendre. TODO: Verify this
    """
    xs,ws = np.polynomial.legendre.leggauss(self.ndiv)

    # Experiment introducing adjacent.
    # TODO: Split up in 8 integrals

    # U-integral (TODO: Should we scale by sm, i.e. self.hw
    us = self.hw * xs + np.fabs(proj[0])
    us2 = us**2
    adj = (self.hh - np.fabs(proj[1]))
    denom = (us2 + adj**2)
    result = adj * np.sum(ws * (1.0 / denom))

    adj = (self.hh + np.fabs(proj[1]))
    denom = (us2 + adj**2)
    result = result + adj*np.sum(ws * (1.0 / denom))
    result = result * self.hw

    # V-integral
    vs = self.hh * xs + np.fabs(proj[1])
    vs2 = vs**2
    adj = (self.hw - np.fabs(proj[0]))
    denom = (vs2 + adj**2)
    result = result + self.hh * adj * np.sum(ws * (1.0 / denom))
    adj = (self.hw + np.fabs(proj[0]))
    denom = (vs2 + adj**2)
    result = result + self.hh * adj * np.sum(ws * (1.0 / denom))

    return result

  def fstart(self, pulse, proj, vdists):
    c = pulse.c
    fs = pulse.fs

    z = proj[2]

    fSampleStop  = fs * (np.max(vdists) / c)
    fSampleStart = fs * (np.min(vdists) / c)

    if (np.fabs(proj[0]) < self.hw):
      if (np.fabs(proj[1]) < self.hh):
        fSampleStart = fs * np.fabs(z) / c
      else:
        fSampleStart = fs * np.sqrt(proj[2]**2 + (np.fabs(proj[1]) - self.hh)**2) / c
    elif (np.fabs(proj[1]) < self.hh):
      fSampleStart = fs * np.sqrt(proj[2]**2 + (np.fabs(proj[0]) - self.hw)**2) / c
    return fSampleStart

  def calc_direct_response(self, pulse, point, delay=0.0):
    """
    Compute direct response. Consider rewriting this as a sum 4 integrals
    """
    proj, vdists, dist = self.local(point)

    fSampleSignalStart = 0.0
    w = pulse.w
    c = pulse.c
    fs = pulse.fs

    invfs = 1.0 / fs

    rho = 1000.0
    factor = - 1.0 * rho * c / (2*np.pi)

    iSampleSignalStart = int(np.floor(fSampleSignalStart))

    z = proj[2]

    fSampleStop  = fs * (np.max(vdists) / c)
    fSampleStart = self.fstart(pulse, proj, vdists)

    iSampleStart = int(fSampleStart + delay*fs) + 1
    iSampleStop  = int(fSampleStop + delay*fs + w*fs) + 1

    t1 = z / c

    spatial = self.direct(proj)

    signal = np.zeros(iSampleStop)
    for iSample in range(iSampleStart, iSampleStop):
      iSampleSignal = iSample - iSampleSignalStart
      t = invfs * iSample - delay# Figure this out
      direct = - pulse.eval(t-t1) * spatial
      signal[iSampleSignal] = factor*direct

    return np.array(signal)

  def calc_edge_response(self, pulse, iEdge, point, delay=0.0):
    """
    Compute edge response
    """
    proj, vdists, dist = self.local(point)

    # TODO: Compute this globally
    fSampleSignalStart = 0.0

    w = pulse.w
    c = pulse.c
    fs = pulse.fs

    rho = 1000.0
    factor = - 1.0 * rho * c / (2*np.pi)

    invfs = 1.0 / fs

    iSampleSignalStart = np.floor(fSampleSignalStart)

    z = proj[2]

    t1 = z / c

    vertices = []

    # Figure out where to integrate and what is the adjacent edge
    if (iEdge == 0):
      vertices.append(0)
      vertices.append(3)
      hl = self.hh # integration
      pl = proj[1]
      nl = self.ndiv
      hs = self.hw
      ps = proj[0]
      if (ps > 0.0):
        adjacent = hs - np.fabs(ps)
      else:
        adjacent = hs + np.fabs(ps)
    elif (iEdge == 1):
      vertices.append(0)
      vertices.append(1)
      hl = self.hw # integration
      pl = proj[0]
      nl = self.ndiv
      hs = self.hh
      ps = proj[1]
      if (ps > 0.0):
        adjacent = hs - np.fabs(ps)
      else:
        adjacent = hs + np.fabs(ps)
    elif (iEdge == 2):
      vertices.append(1) # 0
      vertices.append(2) # 3
      hl = self.hh # integration
      pl = proj[1]
      nl = self.ndiv
      hs = self.hw
      ps = proj[0]
      if (ps < 0.0):
        adjacent = hs - np.fabs(ps)
      else:
        adjacent = hs + np.fabs(ps)
    elif (iEdge == 3):
      vertices.append(3)
      vertices.append(2)
      hl = self.hw # integration
      pl = proj[0]
      nl = self.ndiv
      hs = self.hh
      ps = proj[1]
      if (ps < 0.0):
        adjacent = hs - np.fabs(ps)
      else:
        adjacent = hs + np.fabs(ps)

    t12min = np.min(vdists[vertices]) / c
    fSampleStart = fs * t12min
    t12max = np.max(vdists[vertices]) / c
    fSampleStop  = fs * t12max

    if (np.fabs(pl) < hl):
      t12min = np.sqrt(z**2 + adjacent**2) / c
      fSampleStart = fs * t12min

    iSampleStart = int(fSampleStart + delay*fs) + 1
    iSampleStop  = int(fSampleStop + delay*fs + w*fs) + 1

    fSigma = dotdict()
    iSigma = dotdict()

    print("ps: %f, pl: %f, adjacent: %f" % (ps,pl,adjacent))

    xs, ws = np.polynomial.legendre.leggauss(self.ndiv)

    # l's
    ls = hl * xs

    iSigma.low.max = nl
    iSigma.low.min = 0
    iSigma.up.max = nl
    iSigma.up.min = 0

    # TODO: Figure out if it is enough to control these
    iSigma.low.low = 0
    iSigma.low.high = 0
    iSigma.up.low = 0
    iSigma.up.high = 0

    if (np.fabs(pl) < hl):
      indices = np.where(ls > pl)[0]
      if (len(indices) > 0):
        iSigma.low.low  = indices[0]
        iSigma.low.high = indices[0]
        iSigma.up.low  = indices[0]
        iSigma.up.high = indices[0]
      else:
        iSigma.low.low  = nl
        iSigma.low.high = nl
    elif (pl > hl):
      iSigma.low.low  = nl
      iSigma.low.high = nl
    else:
      iSigma.up.low  = 0
      iSigma.up.high = 0

    debug_low = []
    debug_high = []
    debug_count = []

    #gen = (iSample for iSample in range(iSampleStart, iSampleStop))
    gen = (iSample for iSample in range(0, iSampleStop))

    signal = np.zeros(iSampleStop)
    E = np.zeros(2)
    for iSample in gen:
      if (iSample < iSampleStart):
        continue
      iSampleSignal = iSample - iSampleSignalStart;
      t = invfs * iSample - delay;
      # Compute two sigma segments
      deltaSigma = np.sqrt(max(0.0, (c*t)**2 - z**2 - adjacent**2))

      # Lower sigma range
      fSigma.low.max = np.clip(-hl, hl, pl)
      fSigma.low.min = np.clip(pl - deltaSigma, -hl, hl)

      # Upper sigma range
      fSigma.up.max = np.clip(pl + deltaSigma, -hl, hl)
      fSigma.up.min = np.clip(pl, -hl, hl)

      if (t > t12min + w) and (t < t12max + w):
        # Shorten ranges
        deltaSigma = np.sqrt(max(0.0, (c*(t-w))**2 - z**2 - adjacent**2))
        fSigma.up.min  = np.clip(pl + deltaSigma, -hl, hl)
        fSigma.low.max = np.clip(pl - deltaSigma, -hl, hl);

      # Remove sigma from upper
      for isigma in range(iSigma.up.low, iSigma.up.max):
        if (isigma == iSigma.up.high):
          break
        sigma = ls[isigma]
        if (sigma > fSigma.up.min):
          break
        iSigma.up.low = isigma + 1
        weight = - hl * ws[isigma-1] # We remove
        pulse.adjust(pl-sigma, weight, adjacent, z, E) # WAS sigma + pl

      # Add sigma to upper (iSigma.up.high was not previously included)
      for isigma in range(iSigma.up.high, iSigma.up.max):
        sigma = ls[isigma]
        if (sigma > fSigma.up.max):
          break
        iSigma.up.high = isigma + 1
        weight = hl * ws[isigma-1]
        pulse.adjust(pl - sigma, weight, adjacent, z, E)

      # Remove sigma from lower (excluding iSigma.low.min)
      for isigma in range(iSigma.low.high, iSigma.low.min, -1):
        if (sigma == iSigma.low.low):
          break
        sigma = ls[isigma-1]
        if (sigma < fSigma.low.max):
          break
        iSigma.low.high = isigma - 1
        weight = -hl * ws[isigma-1]
        pulse.adjust(pl - sigma, weight, adjacent, z, E)

      # Add sigma to lower
      for isigma in range(iSigma.low.low, iSigma.low.min, -1):
        sigma = ls[isigma-1]
        if (sigma < fSigma.low.min):
          break
        iSigma.low.low = isigma - 1
        weight = hl * ws[isigma-1]
        pulse.adjust(pl - sigma, weight, adjacent, z, E)

      edge = 0.0
      for iTerm in range(2):
        edge = edge + E[iTerm] * pulse.tbf(iTerm, t)
      signal[iSample] = factor * adjacent * edge

      if (debug):
        debug_high.append("fSigmaMin[Upper]: %1.5f, fSigmaMax[Upper]: %1.5f" % (fSigma.up.min, fSigma.up.max))
        debug_high.append("iSigmaMin[Upper]: %d, iSigmaMax[Upper]: %d" % (iSigma.up.low, iSigma.up.high))
        debug_low.append("fSigmaMin[Lower]: %1.5f, fSigmaMax[Lower]: %1.5f" % (fSigma.low.min, fSigma.low.max))
        debug_low.append("iSigmaMin[Lower]: %d, iSigmaMax[Lower]: %d" % (iSigma.low.low, iSigma.low.high))
        debug_count.append([iSigma.low.high - iSigma.low.low, iSigma.up.high - iSigma.up.low])
    if (debug):
      debug_count = np.array(debug_count).T
      print('\n'.join(debug_low))
      print('\n'.join(debug_high))
      print(debug_count)
    return signal

  def calc_fnm_value(self, realdist, u, v, z, vdists):
    """Compute response using fast nearfield method.

    """
    proj = np.array([u,v,z])
    spatial = self.direct(proj)

    # Edge terms
    edge = 0
    # Edge 0:
    if ((realdist < vdists[0]) or realdist < vdists[3]):
      if (u > 0):
        adjacent = hs - np.fabs(ps)
        hl = self.hh # integration
        pl = proj[1]
        hs = self.hw

        # TODO: Figure out how to compute strip
        deltaSigma = np.sqrt(max(0.0, realdist**2 - z**2 - adjacent**2))
    return spatial

  def calc_sir_fnm(self, point, pulse, delay=0):
    """
    Compute response using FNM. TODO: FIX ME
    """
    c  = pulse.c
    fs = pulse.fs
    ds = c / fs

    proj, vdists, dist = self.local(point)
    z = np.fabs(proj[2])
    u = proj[0]
    v = proj[1]

    edges = self.vertices - np.roll(self.vertices,-1,axis=0)
    edges = edges / np.r_[3*[np.sqrt(np.sum(edges**2,axis=1))]].transpose()
    p2edges = self.vertices - np.r_[4*[point]]
    edists = np.sum(np.cross(p2edges,edges)**2,axis=1)
    edists = np.sqrt(edists)

    # Names with capitals are float samples
    Dists        = edists / ds + delay * fs
    iDists       = np.ceil(Dists).astype(np.int)
    Vdists       = vdists / ds + delay * fs
    iCornerDists = np.ceil(Vdists).astype(np.int)

    # Check boundary
    iSampleInside = min(iDists)

    iSampleEnd    = np.ceil(max(Vdists)).astype(np.int) + 1

    Dist2plane       = z / ds + delay*fs
    iDist2plane      = np.ceil(Dist2plane).astype(np.int)
    iDist2planeStore = iDist2plane

    iDist2plane   = iSampleEnd

    outside = False
    if (abs(u) > self.hw or abs(v) > self.hh):
      outside = True
    if abs(u) > self.hw and abs(v) > self.hh:
      # Outside and in the corner
      iSampleInside = min(iCornerDists)
    elif (u > self.hw):
      iSampleInside = iDists[3]
    elif (u < -self.hw):
      iSampleInside = iDists[1]
    else:
      # Including edges
      if v > self.hh:
        iSampleInside = iDists[0]
      elif v < -self.hh:
        iSampleInside = iDists[2]
      else:
        # Inside - include the first integral
        iDist2plane   = iDist2planeStore

    results = []

    # Will be left out if outside
    for iSample in range(iDist2plane,iSampleInside):
      realdist = ds*iSample - delay*c
      bla = self.calc_fnm_value(realdist, u, v, z, vdists)
      results.append(bla)

    for iSample in range(iSampleInside,iSampleEnd):
      realdist = ds*iSample - delay*c
      bla = self.calc_fnm_value(realdist, u, v, z, vdists)
      results.append(bla)

    return results


  def calc_sir(self, point, pulse, delay=0):
    """
    Returns arc lengths divided with distance and radii for all samples
    """
    c  = pulse.c
    fs = pulse.fs
    ds = c / fs

    proj, vdists, dist = self.local(point)
    dist2plane = np.fabs(proj[2])
    u = proj[0]
    v = proj[1]

    # Do we really need distances to edges
    edges = self.vertices - np.roll(self.vertices,-1,axis=0)
    edges = edges / np.r_[3*[np.sqrt(np.sum(edges**2,axis=1))]].transpose()
    p2edges = self.vertices - np.r_[4*[point]]
    dists = np.sum(np.cross(p2edges,edges)**2,axis=1)
    dists = np.sqrt(dists)

    resp       = []
    radii      = []
    realdists  = []

    Dists        = dists / ds + delay * fs
    iDists       = np.ceil(Dists).astype(np.int)
    Vdists       = vdists / ds + delay * fs
    iCornerDists = np.ceil(Vdists).astype(np.int)

    arcLengths = []

    # Check boundary
    iSampleInside = min(iDists)

    iSampleEnd    = np.ceil(max(Vdists)).astype(np.int) + 1

    Dist2plane       = dist2plane / ds + delay*fs
    iDist2plane      = np.ceil(Dist2plane).astype(np.int)
    iDist2planeStore = iDist2plane

    arrival_samples      = 1e6*np.ones(9)
    # Floating point
    arrival_samples[0:4] = Vdists

    iDist2plane   = iSampleEnd

    nDiscontinuities = 4

    outside = False
    if (abs(u) > self.hw or abs(v) > self.hh):
      outside = True
    if abs(u) > self.hw and abs(v) > self.hh:
      # Outside and in the corner
      iSampleInside = min(iCornerDists)
    elif (u > self.hw):
      iSampleInside = iDists[3]
      arrival_samples[4:6] = [Dists[1], Dists[3]]
      nDiscontinuities = 6
    elif (u < -self.hw):
      iSampleInside = iDists[1]
      arrival_samples[4:6] = [Dists[1], Dists[3]]
      nDiscontinuities = 6
    else:
      # Including edges
      if v > self.hh:
        iSampleInside = iDists[0]
        arrival_samples[4:6] = [Dists[0], Dists[2]]
        nDiscontinuities = 6
      elif v < -self.hh:
        iSampleInside = iDists[2]
        arrival_samples[4:6] = [Dists[0], Dists[2]]
        nDiscontinuities = 6
      else:
        # Inside - include the first
        iSampleInside = iSampleInside
        iDist2plane   = iDist2planeStore
        arrival_samples[1:5] = Vdists
        arrival_samples[5:9] = Dists
        arrival_samples[0] = Dist2plane
        nDiscontinuities = 9

    # Will be left out if outside
    for iSample in range(iDist2plane,iSampleInside):
      realdist = ds*iSample - delay*c
      if outside:
        radius = 0
      else:
        radius = np.sqrt(realdist**2 - dist2plane**2)
      arcLengths.append(2*np.pi)
      resp.append(radius*2*np.pi / realdist)
      radii.append(radius)
    for iSample in range(iSampleInside,iSampleEnd):
      realdist = ds*iSample - delay*c
      radius = np.sqrt(max(realdist**2 - dist2plane**2,0))

      arc, angles = self.calc_arc_length(radius,realdist,u,v,vdists)

      arcLengths.append(arc)
      radii.append(radius)
      realdists.append(realdist)

    radii      = np.array(radii)
    arcs       = np.array(arcLengths)
    dists      = np.array(realdists)

    # We compute radii * arcs / dists afterwards
    return (arcs, radii, dists, np.sort(arrival_samples), nDiscontinuities)


  def calc_sir_diff(self,point,pulse,delay=0):
    """
    Compute derivative of spatial impulse response
    """
    c  = pulse.c
    fs = pulse.fs
    ds = c/fs

    arcs, radii, dists, arrival_samples, nDiscontinuities = self.calc_sir(point, pulse, delay)

    resp = arcs * radii / dists

    # To match Field II
    resp  = resp / (2*np.pi)

    # Wrong for the first and last sample (simple)
    radii = np.r_[radii[0],radii,radii[-1]]

    # Radii in plane doesn't make sense before touch
    dr = 0.5*(np.roll(radii,-1) - np.roll(radii,1))

    # We don't use the first and last result
    ints = np.cumsum(resp*dr[1:-1])

    diffs = ints
    diffs = deepcopy(diffs)
    for index in range(len(diffs)-1):
      diffs[index] = (diffs[index+1] - diffs[index]) * fs

    diffs = np.r_[np.zeros(2),diffs]
    return diffs

  def calc_arc_length(self, radius, realdist, u, v, vdists):
    """Compute sum of arc lengths. Result range is [0,2pi]. The arc
    length as well as the four contributions from each edge are
    returned.

    """

    angles = np.array([0,0,
                       np.pi/2,np.pi/2,
                       np.pi,np.pi,
                       3*np.pi/2,3*np.pi/2])

    if (radius < np.finfo(np.float32).eps):
      return (0, angles)

    args = np.r_[(self.hw - u)/radius, (self.hh - v)/radius,
                 (self.hw + u)/radius, (self.hh + v)/radius]

    args = np.maximum(args,-1)
    args = np.minimum(args,1)

    if (radius + u > self.hw):
      if (realdist < vdists[3]):
        if realdist < vdists[0] and abs(v) > self.hh:
          angles[0] = 0.0
        else:
          angles[0] = - np.arccos(args[0])
      else:
        angles[0] = 0.0
      if (realdist < vdists[0]):
        if (realdist < vdists[3]) and abs(v) > self.hh:
          angles[1] = 0.0
        else:
          angles[1] = + np.arccos(args[0])
      else:
        angles[1] = np.pi/2
    if (radius + v > self.hh):
      if (realdist < vdists[0]):
        if (realdist < vdists[1] and abs(u) > self.hw):
          angles[2] = np.pi/2
        else:
          angles[2] = np.pi/2 - np.arccos(args[1])
      else:
        angles[2] = np.pi/2
      if (realdist < vdists[1]):
        if (realdist < vdists[0] and abs(u) > self.hw):
          angles[3] = np.pi/2
        else:
          angles[3] = np.pi/2 + np.arccos(args[1])
      else:
          angles[3] = np.pi
    if (radius - u > self.hw):
      if (realdist < vdists[1]):
        if (realdist < vdists[2]) and (abs(v) > self.hh):
          angles[4] = np.pi
        else:
          angles[4] = np.pi - np.arccos(args[2])
      else:
        angles[4] = np.pi
      if (realdist < vdists[2]):
        # Python here,
        if realdist < vdists[1] and (abs(v) > self.hh):
          angles[5] = np.pi
        else:
          angles[5] = np.pi + np.arccos(args[2])
      else:
        angles[5] = 3*np.pi/2.0

    if (radius - v > self.hh):
      if (realdist < vdists[2]):
        if (realdist < vdists[3] and (abs(u) > self.hw)):
          angles[6] = 3*np.pi/2
        else:
          angles[6] = 3*np.pi/2 - np.arccos(args[3])
      else:
        angles[6] = 3*np.pi/2.0
      if (realdist < vdists[3]):
        if (realdist < vdists[2]) and abs(u) > self.hw:
          angles[7] = 3*np.pi/2
        else:
          angles[7] = 3*np.pi/2 + np.arccos(args[3])
      else:
        angles[7] = 2*np.pi

    arc = 2*np.pi - sum(angles[1::2] - angles[0:-1:2])

    return arc, angles

# Local variables: #
# indent-tab-mode: nil #
# tab-width: 2 #
# python-indent: 2 #
# indent-tabs-mode: nil #
# End: #
