import numpy as np

# This feature can be used to write out formulas using SymPy
test = False
if (test):
    from sympy import *
else:
    from numpy import (sin,cos)

def euler2rot(a,b,c,conv='zxz',intrinsic=True):
    """
    The angles a,b,c are used for rotations based on the given
    convention

    Classical Euler angles:
    z-x-z, x-y-x, y-z-y, z-y-z, x-z-x, y-x-y

    Tait-Bryan angles (aka Cardan angles, nautical angles, heading,
    elevation, and bank, or yaw, pitch, and roll):
    x-y-z, y-z-x, z-x-y, x-z-y, z-y-x, y-x-z

    All rotation may be intrinsic or extrinsic.

    For intrinsic rotations, the rotations are about the axes of the
    rotating coordinate system XYZ, which changes its orientation
    after each elemental rotation.

    For extrinsic rotations, the rotations are about the original
    axes of the coordinate system (motionless)
    """
    c1 = cos(a)
    s1 = sin(a)
    c2 = cos(b)
    s2 = sin(b)
    c3 = cos(c)
    s3 = sin(c)

    # For testing
    if (test):
        [c1,s1,c2,s2,c3,s3] = symbols(['c1','s1','c2','s2','c3','s3'])

    # We may need different ordering
    _c1 = c1
    _s1 = s1
    _c2 = c2
    _s2 = s2
    _c3 = c3
    _s3 = s3
    if (['zxz','yxy','xzx'].count(conv)) > 0:
        if (conv=='xzx'):
            _c1 = c1
            _s1 = s1
            _c3 = c2
            _s3 = s2
            _c2 = c3
            _s2 = s3
        elif (conv=='yxy'):
            _c2 = c1
            _s2 = s1
            _c1 = c2
            _s1 = s2
            _c3 = c3
            _s3 = s3
        elif (conv=='zxz'):
            _c3 = c1
            _s3 = s1
            _c1 = c2
            _s1 = s2
            _c2 = c3
            _s2 = s3
        Rx = np.array(
            [[1,    0,   0],
             [0,  _c1,-_s1],
             [0,  _s1, _c1]])
        Ry = np.array(
            [[_c2,  0, _s2],
             [0,    1,   0],
             [-_s2, 0, _c2]])
        Rz = np.array(
            [[_c3, -_s3, 0],
             [_s3,  _c3, 0],
             [0,      0, 1]])
        if (conv=='xzx'):
            R3 = np.array(
                [[1,    0,   0],
                 [0,  _c2,-_s2],
                 [0,  _s2, _c2]])
        elif (conv=='yxy'):
            R3 = np.array(
                [[_c3,  0, _s3],
                 [0,    1,   0],
                 [-_s3, 0, _c3]])
        elif (conv=='zxz'):
            R3 = np.array(
                [[_c2, -_s2, 0],
                 [_s2,  _c2, 0],
                 [0,      0, 1]])
        if (intrinsic):
            R = {'xzx' : np.dot(np.dot(Rx,Rz),R3),
                 'yxy' : np.dot(np.dot(Ry,Rx),R3),
                 'zxz' : np.dot(np.dot(Rz,Rx),R3)}[conv]
        else:
            R = {'xzx' : np.dot(np.dot(R3,Rz),Rx),
                 'yxy' : np.dot(np.dot(R3,Rx),Ry),
                 'zxz' : np.dot(np.dot(R3,Rx),Rz)}[conv]
        return R
    elif (['xzy','xyz','yxz','yzx','zyx','zxy'].count(conv)) > 0:
        # Tait-Bryan angles
        if conv=='xzy':
            _c2 = c3
            _s2 = s3
            _c3 = c2
            _s3 = s2
        elif conv=='yxz':
            _c1 = c2
            _s1 = s2
            _c2 = c1
            _s2 = s1
        elif conv=='zyx':
            _c1 = c3
            _s1 = s3
            _c3 = c1
            _s3 = s1
        elif conv=='yzx':
            _c1 = c2
            _s1 = s2
            _c2 = c3
            _s2 = s3
            _c3 = c1
            _s3 = s1
        elif conv=='zxy':
            _c1 = c3
            _s1 = s3
            _c2 = c1
            _s2 = s1
            _c3 = c2
            _s3 = s2

        Rx = np.array(
            [[1,  0, 0],
             [0,  _c1, -_s1],
             [0,  _s1, _c1]])
        Rz = np.array(
            [[_c3, -_s3, 0],
             [_s3,  _c3, 0],
             [0,    0, 1]])
        Ry = np.array(
            [[_c2,  0, _s2],
             [0,   1, 0],
             [-_s2,  0, _c2]])
        if (intrinsic):
            # Note the Ambiguity xyz (intrinsic) = zyx (extrinsic)
            R = {'xzy' : np.dot(np.dot(Rx,Rz),Ry),
                 'xyz' : np.dot(np.dot(Rx,Ry),Rz),
                 'yxz' : np.dot(np.dot(Ry,Rx),Rz),
                 'yzx' : np.dot(np.dot(Ry,Rz),Rx),
                 'zyx' : np.dot(np.dot(Rz,Ry),Rx),
                 'zxy' : np.dot(np.dot(Rz,Rx),Ry)}[conv]
        else:
            R = {'xzy' : np.dot(np.dot(Rz,Ry),Rx),
                 'xyz' : np.dot(np.dot(Rz,Ry),Rx),
                 'yxz' : np.dot(np.dot(Rz,Rx),Ry),
                 'yzx' : np.dot(np.dot(Rx,Rz),Ry),
                 'zyx' : np.dot(np.dot(Rx,Ry),Rz),
                 'zxy' : np.dot(np.dot(Ry,Rx),Rz)}[conv]
        return R
    else:
        raise Exception('Convention %s not supported' % conv)
