# -*- coding: utf-8; tab-width: 2; python-indent: 2; indent-tabs-mode: nil -*-

import numpy as np

from scipy.optimize import fsolve

def grating_lobe_angle(steering_angle, pitch, f0, c=1540, **kwargs):
  """
  Return the angle at which grating lobes occur.

  Parameters
  ----------
  steering_angle:
      Steering angle for the beam. [rad]
  pitch:
      Transducer pitch                       [m]
  f0:
      Center frequency                       [m]
  c:
      Speed of sound                         [m/s]


  kwargs:

  ============= =============================================================
   Name           Description
  ============= =============================================================
  return_denom    Return denominator and angle used to find position.
  ============= =============================================================

  Returns
  -------
  phi:
      An array with angles for the grating lobes. Grating lobes
      are considered only when they appear at angles < 90 deg.
      A non-steered lambda-pitch array has grating lobes at 90
      deg. These will not be reported.      [rad]

  den:
      The denominator of the one-way array response

  theta:
      The angles at which den has been evaluated.

  `den` and `theta` are returned only if `return_denom` is set to `True`

  """

# The radiation pattern of an array is:
#
#
#          sin(pi/lambda * pitch * NumElem * sin(theta + steering_angle))
# P = ------------------------------------------------------------------
#          sin(pi/lambda*pitch*sin(theta+steering_angle))
#
#
#  A grating lobe appears when the denominator is equal to 0 at an
#  angle that is different than the steering angle.
#  This algorithm finds only the first grating lobes closest to the
#  transducer.
#


# The algorithm works as follows:
# 1. Create a fine-spaced vector for the angle theta and
#    split the vector into two subvectors - left and right.
#    The left sub-vector is for the negative angles
#
# 2. Get the derivative of the response from the steering angle - outward
#
# 3. Find the points where the derivative crosses the zero from a
#    a negative to a positive direction
#

  opts = {
      'return_denom': False,
      }

  opts.update(kwargs)
  return_denom = opts['return_denom']
  wavelen = c / f0

  theta = np.r_[-910:911] * np.pi / 1800   # Very, very fine spacing
  theta_d = np.r_[-910:911] / 10.


  steering_angle_d = np.round(steering_angle * 180 / np.pi * 10) / 10
  steering_angle = steering_angle_d * np.pi / 180
  den = np.abs(np.sin(np.pi * (pitch / wavelen) * np.sin(theta - steering_angle) ))
  steer_idx = np.where(theta_d == steering_angle_d)[0]

  if len(steer_idx) == 0:
      raise RuntimeError('Could not find position of steering angle')
  else:
      steer_idx = steer_idx[0]


  left = den[1: steer_idx]
  right = den[steer_idx + 1:-1]
  left_theta = theta[: steer_idx]
  right_theta = theta[steer_idx + 1:]

  phi_left = np.array([])
  phi_right= np.array([])

  idx_grating = np.array([])
  if len(left) > 0:
      left = left[-1::-1]
      left_theta = left_theta[-1::-1]
      deriv = np.diff(left)
      left_theta = left_theta[1:]
      sel_grating = np.logical_and(deriv[:-1] <= 0, deriv[1:] > 0)
      idx_grating = np.where(sel_grating)[0]
      if len(idx_grating) > 0:
          phi_left = left_theta[idx_grating]
          phi_left = phi_left[ left[idx_grating] < 0.3]
  idx_left = idx_grating

  if len(right) > 0:
       deriv = np.diff(right)
       right_theta = right_theta[1:]
       sel_grating = np.logical_and(deriv[:-1] <= 0, deriv[1:] > 0)
       idx_grating = np.where(sel_grating)[0]
       if len(idx_grating) > 0:
           phi_right = right_theta[idx_grating]
           phi_right = phi_right[right[idx_grating] < 0.3]

  phi = np.r_[phi_left, phi_right]


  return (phi, den, theta) if return_denom else phi

def beam_pattern(theta, pitch, nelem, f0, kerf=0,
                 steering=0, c=1540., monostatic=False, weight=None):
    '''
    Continuous-wave approximation of beam pattern.
    '''

    u = np.sin(theta - steering)

    wavelen = c / f0

    K = 2 * np.pi / wavelen if monostatic else np.pi / wavelen


    nom = np.sin(K * nelem * pitch * u)
    den = np.sin(K * pitch * u)

    arrayPattern = np.ones_like(nom)
    sel = (den != 0)
    #sel = (np.abs(den) > 10*np.finfo(np.float32).eps)

    arrayPattern[sel] = nom[sel] / den[sel] / nelem

    elemFactor = 1 if kerf == 0 else elem_factor(theta, pitch-kerf, f0, c)

    beamPattern = arrayPattern * elemFactor

    return beamPattern




def elem_factor(theta, width, f0, c=1540, weight=None):
  '''
  Element factor as function of angle.

  Paramaters:
  -----------
  theta: array_like or scalar
      Angle at which to calculate the element factor
  width: scalar
      Width of element   [m]
  f0: scalar
      Center frequency   [Hz]
  c: scalar
      Speed of sound
  weight: callable (function)
      Weighting function. Single argument - steering angle.

  Returns
  -------
  h: Element sensitivity at steering angle (theta)

  Selfridge, Kino, Khuri-Yakub, A theory for the radiation pattern
  of a narrow-strip acoustic transducer, Appl. Phys. Lett 37(1),
  pp 35-36, 1980

  By default we use $\cos^2$ instead of the $\cos$ used in the paper
  '''

  if (weight == None):
    weight = lambda x: np.cos(x)**2
  wavelen = c / f0

  u = np.sin(theta)
  h = np.sinc(width / wavelen * u) * weight(theta)
  return h

def accept_angle(level, width, f0, c=1540, **kwargs):
  '''
  Acceptance angle (half) for transducer elements. Continuous-wave approx

  This function solves the equation `elem_factor() == level`

  Parameters
  ----------
  level: scalar
      Level at which to calculate the acceptance level. linear (0 .. 1) or dB
  width: scalar
      Width of the transducer element
  f0: scalar
      Center frequency of the transducer elements
  c: scalar
      Speef of sound

  kwargs: named parameters

  ============= =============================================================
   Name            Description
  ============= =============================================================
  level_scale    Scale of `level`. One of 'linear', 'log'. Default 'linear'.
  baffle_type    String: 'soft', 'hard', 'default'
  weight         Callable weight(angle). Used when `baffle_type` == `default`
  ============= =============================================================

  Returns
  -------
  The steering angle (positive) where the sensitivity equals `level`


  Notes
  -----
  `weight` is defined as `lambda x: np.cos(x)**2`, and is used when
           `baffle_type` is `default`.

  '''

  opts = {
      'level_scale': 'linear',
      'baffle_type': 'default',
      'weight': lambda x: np.cos(x)**2,
      }
  opts.update(kwargs)

  level_scale = opts['level_scale'].lower()[:3]
  baffle_type = opts['baffle_type'].lower().strip()

  level = 10 ** (level / 20) if level_scale == 'log' else level

  if baffle_type == 'hard':
    weight = lambda x: 1
  elif baffle_type == 'soft':
    weight = lambda x: np.cos(x)
  else:
    weight = opts['weight']

  theta = fsolve(lambda x: elem_factor(x, width, f0, c, weight) - level, np.pi / 4)
  return theta


# Local variables: #
# tab-width: 2 #
# python-indent: 2 #
# indent-tabs-mode: nil #
# End: #
