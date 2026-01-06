#  -*- coding: utf-8 -*-
# Copyright (C) 2007-2026  CEA, EDF
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

    def test_UMesh_buildDescendingCrude_explodeCrude_0(self):
        """
        [EDF33593] : fine control on buildDescendingConnectivity
        """
        NX = 4
        arr = mc.DataArrayDouble(NX)
        arr.iota()
        m = mc.MEDCouplingCMesh()
        m.setCoords(arr, arr)
        m2D = m.buildUnstructured()
        m1D, descIndx, revNodalIndex, revDesc = m2D.buildDescendingConnectivityCrude()
        # revDesc
        self.assertTrue(len(revDesc), (NX - 1) * (NX - 1) + 4)
        revDesc.rearrange(4)
        self.assertTrue(revDesc[:, 0].isIota((NX - 1) * (NX - 1)))
        for i in range(1, 4):
            self.assertTrue(revDesc[:, 0].isEqual(revDesc[:, i]))
        # revNodalIndex : for each node gives # of seg2 sharing it
        self.assertEqual(len(revNodalIndex), NX * NX + 1)
        self.assertTrue(
            revNodalIndex.isEqual(
                mc.DataArrayInt(
                    [0, 2, 6, 10, 12, 16, 24, 32, 36, 40, 48, 56, 60, 62, 66, 70, 72]
                )
            )
        )
        #
        self.assertTrue((descIndx // 4).isIota((NX - 1) * (NX - 1) + 1))
        m1D.checkConsistency()
        m1D = mc.MEDCoupling1SGTUMesh(m1D)
        self.assertEqual(
            m1D.getCoords().getHiddenCppPointer(), m2D.getCoords().getHiddenCppPointer()
        )
        self.assertEqual(m1D.getNumberOfCells(), 4 * (NX - 1) * (NX - 1))
        self.assertEqual(m1D.getCellModelEnum(), mc.NORM_SEG2)
        #
        m = mc.MEDCouplingCMesh()
        m.setCoords(arr, arr, arr)
        m3D = m.buildUnstructured()
        m1D, descIndx, revNodalIndex, revDesc = m3D.explode3DMeshTo1DCrude()
        m1D.checkConsistency()
        m1D = mc.MEDCoupling1SGTUMesh(m1D)
        self.assertEqual(
            m1D.getCoords().getHiddenCppPointer(), m3D.getCoords().getHiddenCppPointer()
        )
        self.assertEqual(m1D.getNumberOfCells(), 12 * (NX - 1) * (NX - 1) * (NX - 1))
        self.assertEqual(m1D.getCellModelEnum(), mc.NORM_SEG2)
        pass

    def test_UMesh_findCommonCells_ct89_0(self):
        """
        [EDF33593] : add compType 8 and 9 to findCommonCells
        """
        coo = mc.DataArrayDouble([(0, 0), (1, 0), (2, 0), (1, 0)])
        m = mc.MEDCouplingUMesh("mesh", 1)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(mc.NORM_SEG3, [0, 2, 1])
        m.insertNextCell(mc.NORM_SEG3, [0, 2, 3])
        m.insertNextCell(mc.NORM_SEG3, [0, 2, 1])
        m.checkConsistency()
        self.assertTrue(m.findCommonCells(8)[0].isEqual(mc.DataArrayInt([0, 1, 2])))
        self.assertTrue(m.findCommonCells(9)[0].isEqual(mc.DataArrayInt([0, 1])))
        self.assertTrue(m.findCommonCells(2)[0].isEqual(mc.DataArrayInt([0, 2])))

    def test_1GTUMesh_computeEulerCharacteristic(self):
        """
        [EDF31086] : MEDCoupling1DGTUMesh.computeEulerCharacteristic
        """
        coo = mc.DataArrayDouble(
            [
                1.9677875479706612,
                -1.4344454267558717,
                148.5520867356951,
                1.7155516292042539,
                -1.7263125116735436,
                147.87223133920384,
                0.9481350791305305,
                -1.06561677318218,
                147.87618431809045,
                0.9189192758350372,
                -0.9279708487895826,
                148.0268677360358,
                1.200491438203907,
                -0.7771242691723659,
                148.55166214174324,
                1.0757582767339655,
                -1.1537872957619573,
                149.20365350652625,
                1.7304879222512763,
                -1.71284091932755,
                149.20154119526893,
                0.63060069877032,
                -1.268946912758893,
                148.02192605406867,
                0.46386721206230075,
                -1.6481588568186534,
                148.5386824372086,
                0.802582743014336,
                -1.479984836084472,
                149.20790682169851,
                0.7337474432330666,
                -1.3091423701933405,
                147.88130443833788,
                1.084107912422434,
                -2.184169361649895,
                147.8848784283266,
                0.7375968039933086,
                -2.323604033133538,
                148.53613406856113,
                1.0893277565584536,
                -2.180837664804089,
                149.20973021713132,
            ],
            14,
            3,
        )
        conn_poly = mc.DataArrayInt(
            [
                0,
                1,
                2,
                3,
                4,
                -1,
                4,
                5,
                6,
                0,
                -1,
                3,
                7,
                8,
                9,
                5,
                4,
                -1,
                10,
                11,
                12,
                8,
                7,
                -1,
                7,
                3,
                2,
                10,
                -1,
                2,
                1,
                11,
                10,
                -1,
                0,
                6,
                13,
                12,
                11,
                1,
                -1,
                12,
                13,
                9,
                8,
                -1,
                13,
                6,
                5,
                9,
            ]
        )
        mesh_3D = mc.MEDCouplingUMesh("poly", 3)
        mesh_3D.setCoords(coo)
        mesh_3D.allocateCells()
        mesh_3D.insertNextCell(mc.NORM_POLYHED, conn_poly)
        mesh_3D.finishInsertingCells()
        ec = mc.MEDCoupling1DGTUMesh(mesh_3D).computeEulerCharacteristic()
        self.assertTrue(ec.isEqual(mc.DataArrayInt32([2])))


if __name__ == "__main__":
    unittest.main()
