import numpy as np


def angle_between(v1, v2, normal=None):
    """
    Returns the angle in radians between vectors 'v1' and 'v2'::

        >>> angle_between((1, 0, 0), (0, 1, 0))
        1.5707963267948966
        >>> angle_between((1, 0, 0), (1, 0, 0))
        0.0
        >>> angle_between((1, 0, 0), (-1, 0, 0))
        3.141592653589793

    https://stackoverflow.com/questions/2827393/
    angles-between-two-n-dimensional-vectors-in-python/
    13849249#13849249

    If normal is given, the angle polarity is determined using the
    cross product of the two vectors.

    """

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    if normal is not None:
        # Get normal vector and cross product to determine sign.
        cross = np.cross(v1_u, v2_u)
        if np.dot(normal, cross) < 0:
            angle = -angle
    return angle


def unit_vector(vector):
    """
    Returns the unit vector of the vector.

    https://stackoverflow.com/questions/2827393/
    angles-between-two-n-dimensional-vectors-in-python/
    13849249#13849249

    """
    return vector / np.linalg.norm(vector)


def get_dihedral(pt1, pt2, pt3, pt4):
    """
    Calculate the dihedral between four points.

    Uses Praxeolitic formula --> 1 sqrt, 1 cross product
    Output in range (-pi to pi).
    From: https://stackoverflow.com/questions/20305272/
    dihedral-torsion-angle-from-four-points-in-cartesian-
    coordinates-in-python
    (new_dihedral(p))

    """

    p0 = np.asarray(pt1)
    p1 = np.asarray(pt2)
    p2 = np.asarray(pt3)
    p3 = np.asarray(pt4)

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))
