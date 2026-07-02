#  -*- coding: utf-8 -*-
# Copyright (C) 2026  CEA, EDF
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#


def MEDCouplingHasPyVistaSupport():
    try:
        import pyvista
    except:
        return False
    return True


def DADfitPlaneL2(self):
    """
    Adjust plan : a*x + b*y + c*z + d = 0

    by minimizing sum of squared ortho distances contained in self (with a,b,c normalized)

    :return: tuple (a, b, c, d)
    """
    import numpy as np

    if self.getNumberOfComponents() != 3:
        raise ValueError("Nb of components must be 3")
    if self.getNumberOfTuples() < 3:
        raise ValueError("At least 3 tuples explected")

    points = self.toNumPyArray()

    P = np.asarray(points, dtype=float)

    centroid = P.mean(axis=0)

    _, _, vh = np.linalg.svd(P - centroid, full_matrices=False)

    normal = vh[-1]

    a, b, c = normal

    s = np.sqrt(a * a + b * b + c * c)

    normal = normal / s

    d = -normal @ centroid

    a, b, c = normal
    return float(a), float(b), float(c), float(d)


def DADmaxDistanceToPlane(self, a, b, c, d):
    """
    Return max distance of this to plan
        a*x + b*y + c*z + d = 0
    """
    import numpy as np

    points = self.toNumPyArray()
    P = np.asarray(points, dtype=float)

    if self.getNumberOfComponents() != 3:
        raise ValueError("Nb of components must be 3")

    normal_norm = np.sqrt(a * a + b * b + c * c)
    if normal_norm == 0:
        raise ValueError("Normal vector is nul.")

    distances = np.abs(P @ np.array([a, b, c]) + d) / normal_norm
    return float(np.max(distances))
