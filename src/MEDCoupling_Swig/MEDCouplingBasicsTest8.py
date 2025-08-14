#  -*- coding: utf-8 -*-
# Copyright (C) 2007-2025  CEA, EDF
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

import medcoupling as mc
import unittest


class MEDCouplingBasicsTest8(unittest.TestCase):
    def test_MEDCoupling_DistanceToPoint_0(self):
        """
        [EDF33340] : check distance returned by MEDCouplingUMesh.distanceToPoint
        """
        coo = mc.DataArrayDouble(
            [
                [0.0, -3.4, 2.30],
                [0.0, -3.7, 1.61],
                [0.0, -4.7, 1.64],
                [0.0, -4.5, 2.44],
                [0.0, -2.8, 0.82],
                [0.0, -2.5, 2.14],
            ]
        )

        m1 = mc.MEDCouplingUMesh.New("Test_distanceToPoint", 2)
        m1.setCoords(coo)
        m1.allocateCells(2)
        m1.insertNextCell(mc.NORM_QUAD4, [0, 1, 2, 3])
        m1.insertNextCell(mc.NORM_QUAD4, [4, 1, 0, 5])
        m1.finishInsertingCells()

        # point is on the plane of QUAD4 (plane is X = 0) inside cell #0
        point = mc.DataArrayDouble([0.0, -3.9, 1.8], 1, 3)

        distance, cell = m1.distanceToPoint(point)
        self.assertEqual(cell, 0)
        self.assertAlmostEqual(distance, 0.0, 12)
        distance0, _ = m1[0].distanceToPoint(point)
        self.assertAlmostEqual(distance0, 0.0, 12)
        distance1, _ = m1[1].distanceToPoint(point)
        self.assertAlmostEqual(distance1, 0.2591719724193925, 10)
        #
        point2 = mc.DataArrayDouble([1.0, -3.9, 1.8], 1, 3)
        distance3, cell = m1.distanceToPoint(point2)
        self.assertAlmostEqual(distance3, 1.0, 12)
        self.assertEqual(cell, 0)
        distance0, _ = m1[0].distanceToPoint(point2)
        self.assertAlmostEqual(distance0, 1.0, 12)
        distance4, _ = m1[1].distanceToPoint(point2)
        self.assertAlmostEqual(distance4, 1.0330392593158102, 12)
        pass


if __name__ == "__main__":
    unittest.main()
