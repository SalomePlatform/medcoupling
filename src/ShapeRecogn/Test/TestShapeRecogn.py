#! /usr/bin/env python3
#  -*- coding: utf-8 -*-
# Copyright (C) 2024  CEA, EDF
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

import ShapeRecogn as sr
import MEDLoader as ml
from numpy.testing import assert_allclose
import unittest


def getResourceFile(filename, levelUp=2):
    ROOT_DIR_0 = "MEDCOUPLING_ROOT_DIR"
    ROOT_DIR_1 = "MEDCOUPLING_RESOURCE_DIR"
    import os
    from pathlib import Path

    if ROOT_DIR_0 in os.environ:
        resourceFile = (
            Path(os.environ[ROOT_DIR_0]).absolute()
            / "share"
            / "resources"
            / "med"
            / filename
        )
        if resourceFile.exists():
            return f"{resourceFile}"
    if ROOT_DIR_1 in os.environ:
        resourceFile = Path(os.environ[ROOT_DIR_1].split(":")[-1]).absolute() / filename
        if resourceFile.exists():
            return f"{resourceFile}"
    p = Path.cwd()
    for i in range(levelUp):
        p = p.parent
    resourceFile = p / "resources" / filename
    if not resourceFile.exists():
        raise RuntimeError(
            f"getResourceFile: could not open resource test file {filename}"
        )
    return f"{resourceFile}"


class MEDCouplingIterativeStatisticsTest(unittest.TestCase):
    def testPlane(self):
        """
        Direct translation of PlaneTest::testArea
        """
        self.assertEqual(
            sr.AllManagedPrimitivesStr(),
            ("Plane", "Sphere", "Cylinder", "Cone", "Torus", "Unknown"),
        )
        fname = "ShapeRecognPlane.med"
        mesh = ml.ReadUMeshFromFile(getResourceFile(fname, 3), 0)
        srMesh = sr.ShapeRecognMeshBuilder(mesh)
        rem = srMesh.recognize()
        areas = srMesh.getAreas()
        self.assertEqual(areas.getNumberOfAreas(), 1)
        self.assertEqual(areas.getNumberOfNodes(0), 36)
        nodeIds = areas.getNodeIds(0)
        f = rem.getAreaPrimitiveType()
        self.assertEqual(areas.getPrimitiveType(0), 0)
        self.assertEqual(areas.getPrimitiveTypeName(0), "Plane")
        self.assertEqual(sr.ConvertStringToPrimitive("Plane"), 0)
        self.assertEqual(sr.ConvertPrimitiveToString(0), "Plane")
        f_normal = rem.getNodeNormal()
        affinePoint = srMesh.buildAreaAffinePoint().getArray()[nodeIds[0]]
        normal = f_normal.getArray()[nodeIds[0]]
        normalRef = sr.DataArrayDouble([0.781525, 0.310606, -0.541056], 1, 3)
        proportion0 = normal[0, 0] / normalRef[0, 0]
        proportion1 = normal[0, 1] / normalRef[0, 1]
        proportion2 = normal[0, 2] / normalRef[0, 2]
        proportion3 = (
            sr.DataArrayDouble.Dot(normal, affinePoint)[0]
            / sr.DataArrayDouble.Dot(normalRef, affinePoint)[0]
        )
        assert_allclose(
            [proportion0, proportion1, proportion2, proportion3],
            [1.0, 1.0, 1.0, 1.0],
            rtol=1e-5,
        )
        angle = sr.DataArrayDouble.CrossProduct(normal, normalRef).magnitude()[0]
        assert_allclose([angle], [0.0], atol=1e-6)


if __name__ == "__main__":
    unittest.main()
