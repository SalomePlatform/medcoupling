#  -*- coding: iso-8859-1 -*-
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
# Author : Aymeric SONOLET (CEA)

import unittest

from medcoupling import (
    DataArrayDouble,
    DataArrayInt,
    MEDCouplingCMesh,
    MEDCouplingUMesh,
)
from MEDLoader import MEDFileUMesh


class CrackAlongTest(unittest.TestCase):
    def testBuildInnerBoundaryAlongM1Group1(self):
        fname = "Pyfile44.med"
        m = MEDCouplingCMesh.New()
        m.setCoordsAt(0, DataArrayDouble.New([0.0, 1.1, 2.3, 3.6, 5.0, 6.5]))
        m.setCoordsAt(1, DataArrayDouble.New([0.0, 1.1, 2.3, 3.6, 5.0]))
        m = m.buildUnstructured()
        m.setName("AnthonyDuplicate")
        m.getCoords().setInfoOnComponents(["X [km]", "Z [mm]"])
        m2 = m.buildDescendingConnectivity()[0][
            [8, 11, 14, 20, 21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 37]
        ]
        m2.setName(m.getName())
        grp = DataArrayInt.New([4, 6, 8])
        grp.setName("Grp")
        grp2 = DataArrayInt.New([9, 16])
        grp2.setName("Grp2")
        mm = MEDFileUMesh.New()
        mm.setMeshAtLevel(0, m)
        mm.setMeshAtLevel(-1, m2)
        mm.setGroupsAtLevel(-1, [grp, grp2])
        grpNode = DataArrayInt.New([4, 21, 23])
        grpNode.setName("GrpNode")
        mm.setGroupsAtLevel(1, [grpNode])
        ref0 = [4, 15, 14, 20, 21, 4, 16, 15, 21, 22, 4, 17, 16, 22, 23]
        ref1 = [4, 9, 8, 14, 15, 4, 10, 9, 15, 16, 4, 11, 10, 16, 17]
        ref2 = [4, 30, 14, 20, 21, 4, 31, 30, 21, 22, 4, 32, 31, 22, 23]
        #
        self.assertEqual(30, mm.getNumberOfNodes())
        self.assertEqual(
            ref0, mm.getMeshAtLevel(0)[[12, 13, 14]].getNodalConnectivity().getValues()
        )
        self.assertEqual(
            ref1, mm.getMeshAtLevel(0)[[7, 8, 9]].getNodalConnectivity().getValues()
        )
        #
        c2o2nN = mm.crackAlong("Grp")
        self.assertEqual(
            {12: {15: 30}, 13: {15: 30, 16: 31}, 14: {16: 31, 17: 32}}, c2o2nN
        )
        self.assertEqual(33, mm.getNumberOfNodes())
        self.assertEqual([4, 6, 8, 17, 18, 19], mm.getGroupArr(-1, "Grp").getValues())
        self.assertEqual([9, 16], mm.getGroupArr(-1, "Grp2").getValues())
        self.assertEqual([4, 21, 23], mm.getGroupArr(1, "GrpNode").getValues())
        self.assertEqual(
            ref2, mm.getMeshAtLevel(0)[[12, 13, 14]].getNodalConnectivity().getValues()
        )  # cells 7,8,9 and 12,13,14 are lying on "Grp" but only 12,13,14 are renumbered
        self.assertEqual(
            ref1, mm.getMeshAtLevel(0)[[7, 8, 9]].getNodalConnectivity().getValues()
        )  #
        # fmt: off
        refValues = DataArrayDouble.New(
            [1.21, 1.32, 1.43, 1.54, 1.65, 1.32, 1.44, 1.56, 1.68, 1.80,
             1.43, 1.56, 1.69, 1.82, 1.95, 1.54, 1.68, 1.82, 1.96, 2.10]
        )
        # fmt: on
        valsToTest = mm.getMeshAtLevel(0).getMeasureField(True).getArray()
        delta = valsToTest - refValues
        delta.abs()
        self.assertTrue(delta.getMaxValue()[0] < 1e-12)

    def testBuildInnerBoundaryAlongM1Group2(self):
        fname = "Pyfile45.med"
        m = MEDCouplingCMesh.New()
        m.setCoordsAt(0, DataArrayDouble.New([0.0, 1.1, 2.3, 3.6, 5.0, 6.5]))
        m.setCoordsAt(1, DataArrayDouble.New([0.0, 1.1, 2.3, 3.6, 5.0]))
        m = m.buildUnstructured()
        m.setName("AnthonyDuplicate")
        m.getCoords().setInfoOnComponents(["X [km]", "Z [mm]"])
        m2 = m.buildDescendingConnectivity()[0][
            [8, 11, 14, 20, 21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 37]
        ]
        m2.setName(m.getName())
        grp = DataArrayInt.New([4, 6])
        grp.setName("Grp")
        grp2 = DataArrayInt.New([9, 16])
        grp2.setName("Grp2")
        mm = MEDFileUMesh.New()
        mm.setMeshAtLevel(0, m)
        mm.setMeshAtLevel(-1, m2)
        mm.setGroupsAtLevel(-1, [grp, grp2])
        grpNode = DataArrayInt.New([4, 21, 23])
        grpNode.setName("GrpNode")
        mm.setGroupsAtLevel(1, [grpNode])
        ref0 = [4, 15, 14, 20, 21, 4, 16, 15, 21, 22, 4, 17, 16, 22, 23]
        ref1 = [4, 9, 8, 14, 15, 4, 10, 9, 15, 16]
        ref2 = [4, 30, 14, 20, 21, 4, 16, 30, 21, 22, 4, 17, 16, 22, 23]
        #
        self.assertEqual(30, mm.getNumberOfNodes())
        self.assertEqual(
            ref0, mm.getMeshAtLevel(0)[[12, 13, 14]].getNodalConnectivity().getValues()
        )
        self.assertEqual(
            ref1, mm.getMeshAtLevel(0)[[7, 8]].getNodalConnectivity().getValues()
        )
        #
        c2o2nN = mm.crackAlong("Grp")
        self.assertEqual({13: {15: 30}, 12: {15: 30}}, c2o2nN)
        self.assertEqual(31, mm.getNumberOfNodes())
        self.assertEqual([4, 6, 17, 18], mm.getGroupArr(-1, "Grp").getValues())
        self.assertEqual([9, 16], mm.getGroupArr(-1, "Grp2").getValues())
        self.assertEqual([4, 21, 23], mm.getGroupArr(1, "GrpNode").getValues())
        self.assertEqual(
            ref2, mm.getMeshAtLevel(0)[[12, 13, 14]].getNodalConnectivity().getValues()
        )  # cells 7,8,9 and 12,13,14 are lying on "Grp" but only 12 and 13 are renumbered
        self.assertEqual(
            ref1, mm.getMeshAtLevel(0)[[7, 8]].getNodalConnectivity().getValues()
        )
        # fmt: off
        refValues = DataArrayDouble.New(
            [1.21, 1.32, 1.43, 1.54, 1.65, 1.32, 1.44, 1.56, 1.68, 1.80,
             1.43, 1.56, 1.69, 1.82, 1.95, 1.54, 1.68, 1.82, 1.96, 2.10]
        )
        # fmt: on
        valsToTest = mm.getMeshAtLevel(0).getMeasureField(True).getArray()
        delta = valsToTest - refValues
        delta.abs()
        self.assertTrue(delta.getMaxValue()[0] < 1e-12)

    def testBuildInnerBoundaryAlongM1Group3(self):
        """Test buildInnerBoundaryAlongM1Group() with *non-connex* cracks"""
        fname = "Pyfile73.med"
        m = MEDCouplingCMesh.New()
        m.setCoordsAt(0, DataArrayDouble([0.0, 1.1, 2.3, 3.6, 5.0]))
        m.setCoordsAt(1, DataArrayDouble([0.0, 1.0, 2.0]))
        m = m.buildUnstructured()
        m.setName("simple")
        m2 = m.buildDescendingConnectivity()[0]
        m2.setName(m.getName())

        # A crack in two non connected parts of the mesh:
        grpSeg = DataArrayInt([3, 19])
        grpSeg.setName("Grp")

        mm = MEDFileUMesh.New()
        mm.setMeshAtLevel(0, m)
        mm.setMeshAtLevel(-1, m2)
        mm.setGroupsAtLevel(-1, [grpSeg])
        c2o2nN = mm.crackAlong("Grp")
        self.assertEqual({1: {1: 15}, 7: {13: 16}}, c2o2nN)
        # self.assertEqual([1, 13], nodes.getValues())
        # self.assertEqual([0, 6], cellsMod.getValues())
        # self.assertEqual([1, 7], cellsNotMod.getValues())
        self.assertEqual(17, mm.getNumberOfNodes())
        self.assertEqual({3, 19, 22, 23}, set(mm.getGroupArr(-1, "Grp").getValues()))

        refValues = DataArrayDouble([1.1, 1.2, 1.3, 1.4, 1.1, 1.2, 1.3, 1.4])
        valsToTest = mm.getMeshAtLevel(0).getMeasureField(True).getArray()
        delta = valsToTest - refValues
        delta.abs()
        self.assertTrue(delta.getMaxValue()[0] < 1e-10)

    def testBuildInnerBoundaryAlongM1Group4(self):
        """Test case where cells touch the M1 group on some nodes only and not
        on full edges (triangle mesh for ex)
        """
        # fmt: off
        coo = DataArrayDouble(
            [0.0, 0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 0.0, 1.0, 1.0, 1.0,
             2.0, 1.0, 3.0, 1.0, 0.0, 2.0, 1.0, 2.0, 2.0, 2.0, 3.0, 2.0],
            12,
            2,
        )
        conn = [3, 0, 4, 1,
                3, 1, 4, 5,
                3, 5, 9, 10,
                3, 5, 10, 6,
                3, 2, 6, 7,
                3, 2, 7, 3,
                3, 4, 8, 9,
                3, 4, 9, 5,
                3, 1, 5, 6,
                3, 1, 6, 2,
                3, 6, 10, 11,
                3, 6, 11, 7]
        # fmt: on
        # Only TRI3:
        connI = DataArrayInt()
        connI.alloc(13, 1)
        connI.iota()
        connI *= 4
        m2 = MEDCouplingUMesh("2D", 2)
        m2.setCoords(coo)
        m2.setConnectivity(DataArrayInt(conn), connI)
        m2.checkConsistency()
        m1, _, _, _, _ = m2.buildDescendingConnectivity()
        grpIds = DataArrayInt([9, 11])
        grpIds.setName("group")
        grpIds2 = DataArrayInt([0, 1])
        grpIds2.setName("group2")
        mfu = MEDFileUMesh()
        mfu.setMeshAtLevel(0, m2)
        mfu.setMeshAtLevel(-1, m1)
        mfu.setGroupsAtLevel(-1, [grpIds, grpIds2])
        nNod = m2.getNumberOfNodes()

        ref0 = [3, 5, 10, 6, 3, 2, 7, 3, 3, 6, 10, 11]
        ref1 = [3, 2, 6, 7, 3, 1, 5, 6, 3, 1, 6, 2, 3, 6, 11, 7]
        ref2 = [3, 2, 12, 7, 3, 1, 5, 12, 3, 1, 12, 2, 3, 6, 11, 13]

        self.assertEqual(ref0, m2[[3, 5, 10]].getNodalConnectivity().getValues())
        self.assertEqual(ref1, m2[[4, 8, 9, 11]].getNodalConnectivity().getValues())

        c2o2nN = mfu.crackAlong("group")
        self.assertEqual({4: {6: 12}, 8: {6: 12}, 9: {6: 12}, 11: {7: 13}}, c2o2nN)

        self.assertEqual(ref0, m2[[3, 5, 10]].getNodalConnectivity().getValues())
        self.assertEqual(ref2, m2[[4, 8, 9, 11]].getNodalConnectivity().getValues())

        m2_bis = mfu.getMeshAtLevel(0)
        m2_bis.checkConsistency()
        m1_bis = mfu.getMeshAtLevel(-1)
        m1_bis.checkConsistency()
        self.assertEqual(nNod + 2, mfu.getNumberOfNodes())
        self.assertEqual(nNod + 2, m2_bis.getNumberOfNodes())
        self.assertEqual(nNod + 2, m1_bis.getNumberOfNodes())
        self.assertEqual([9, 11, 23, 24], mfu.getGroupArr(-1, "group").getValues())
        self.assertEqual([0, 1], mfu.getGroupArr(-1, "group2").getValues())

        m_bis0 = mfu.getMeshAtLevel(-1)
        m_desc, _, _, _, _ = m_bis0.buildDescendingConnectivity()
        m_bis0.checkDeepEquivalOnSameNodesWith(mfu.getMeshAtLevel(-1), 2, 9.9999999)

    def testBuildInnerBoundary5(self):
        """Full 3D test with tetras only. In this case a tri from the group is not duplicated because it is made only
        of non duplicated nodes. The tri in question is hence not part of the final new "dup" group."""
        # fmt: off
        coo = DataArrayDouble(
            [200.0, 200.0, 0.0, 200.0, 200.0, 200.0, 200.0, 0.0, 200.0, 200.0, 0.0, 0.0, 0.0, 200.0, 0.0, 0.0, 200.0, 200.0, 0.0, 0.0, 0.0, 0.0, 0.0, 200.0, 400.0, 200.0, 0.0, 400.0, 200.0, 200.0, 400.0, 0.0, 0.0, 400.0, 0.0, 200.0, 0.0, 100.00000000000016, 200.0, 63.15203310314546, 200.0, 200.0, 134.45205700643342, 200.0, 200.0, 200.0, 100.00000000000016, 200.0, 63.15203310314546, 0.0, 200.0, 134.45205700643342, 0.0, 200.0, 0.0, 100.00000000000016, 0.0, 63.15203310314546, 200.0, 0.0, 134.45205700643342, 200.0, 0.0, 200.0, 100.00000000000016, 0.0, 63.15203310314546, 0.0, 0.0, 134.45205700643342, 0.0, 0.0, 200.0, 200.0, 100.02130053568538, 0.0, 200.0, 100.00938163175135, 200.0, 0.0, 100.02130053568538, 0.0, 0.0, 100.00938163175135, 299.3058739933347, 200.0, 200.0, 400.0, 98.68100542924483, 200.0, 302.8923433403344, 0.0, 200.0, 302.8923433403344, 200.0, 0.0, 400.0, 100.00000000000016, 0.0, 302.8923433403344, 0.0, 0.0, 400.0, 200.0, 98.55126825835082, 400.0, 0.0, 100.02162286181577, 99.31624553977466, 99.99999998882231, 200.0, 99.31624576683302, 100.00000010178034, 0.0, 99.31624560596512, 200.0, 100.0050761312483, 99.31624560612883, 0.0, 100.00507613125338, 200.0, 99.99999995813045, 100.00950673487786, 0.0, 99.99999989928207, 100.0041870621175, 301.29063354383015, 100.0000000093269, 0.0, 301.29063360689975, 0.0, 100.00957769061164, 140.52853868782435, 99.99999963972768, 100.00509135751312, 297.87779091770784, 97.16750463405486, 97.18018457127863],
            46,
            3,
        )
        c0 = [14, 45, 31, 21, 42, 14, 37, 38, 20, 44, 14, 39, 36, 41, 44, 14, 5, 25, 12, 13, 14, 38, 36, 44, 41, 14, 21, 20, 24, 44, 14, 38, 25, 41, 19, 14, 37, 38, 44, 41, 14, 16, 27, 39, 41, 14, 21, 45, 26, 40, 14, 39, 37, 44, 41, 14, 14, 15, 24, 44, 14, 25, 38, 41, 13, 14, 27, 18, 6, 22, 14, 38, 36, 41, 13, 14, 44, 14, 15, 36, 14, 44, 23, 39, 26, 14, 21, 26, 23, 44, 14, 38, 44, 14, 24, 14, 39, 37, 41, 22, 14, 21, 33, 45, 42, 14, 27, 22, 39, 41, 14, 23, 26, 21, 3, 14, 27, 18, 22, 41, 14, 39, 36, 44, 17, 14, 21, 26, 44, 40, 14, 39, 37, 22, 23, 14, 37, 38, 41, 19, 14, 25, 12, 13, 41, 14, 30, 26, 43, 45, 14, 38, 36, 13, 14, 14, 12, 36, 13, 41, 14, 20, 44, 21, 37, 14, 16, 36, 12, 41, 14, 39, 36, 17, 16, 14, 44, 20, 24, 38, 14, 27, 16, 12, 41, 14, 26, 15, 17, 44, 14, 19, 18, 41, 37, 14, 40, 45, 26, 15, 14, 37, 38, 19, 20, 14, 17, 15, 26, 2, 14, 39, 36, 16, 41, 14, 24, 21, 44, 40, 14, 16, 7, 27, 12, 14, 22, 18, 37, 41, 14, 21, 31, 45, 24, 14, 44, 40, 15, 24, 14, 24, 45, 15, 28, 14, 44, 40, 26, 15, 14, 24, 20, 21, 0, 14, 38, 36, 14, 44, 14, 39, 37, 23, 44, 14, 45, 31, 42, 32, 14, 25, 18, 19, 4, 14, 36, 44, 17, 15, 14, 25, 19, 18, 41, 14, 24, 15, 14, 1, 14, 45, 24, 34, 28, 14, 35, 45, 30, 43, 14, 17, 44, 39, 26, 14, 44, 23, 21, 37, 14, 30, 45, 29, 15, 14, 45, 35, 33, 43, 14, 30, 15, 26, 45, 14, 31, 21, 0, 24, 14, 33, 35, 32, 10, 14, 29, 45, 34, 28, 14, 32, 45, 34, 29, 14, 45, 31, 32, 34, 14, 33, 26, 45, 43, 14, 45, 31, 34, 24, 14, 33, 26, 21, 45, 14, 11, 30, 35, 29, 14, 33, 35, 45, 32, 14, 33, 45, 42, 32, 14, 32, 8, 34, 31, 14, 21, 26, 33, 3, 14, 35, 45, 32, 29, 14, 29, 34, 9, 28, 14, 15, 45, 24, 40, 14, 29, 45, 28, 15, 14, 21, 24, 45, 40, 14, 24, 15, 1, 28, 14, 35, 45, 29, 30, 14, 26, 15, 30, 2]
        cI0 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300, 305, 310, 315, 320, 325, 330, 335, 340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 390, 395, 400, 405, 410, 415, 420, 425, 430]
        # fmt: on
        m3 = MEDCouplingUMesh("3D", 3)
        m3.setCoords(coo)
        m3.setConnectivity(DataArrayInt(c0), DataArrayInt(cI0))
        m3.checkConsistency()
        m2, _, _, _, _ = m3.buildDescendingConnectivity()
        grpIds = DataArrayInt([36, 74])
        grpIds.setName("group")
        mfu = MEDFileUMesh()
        mfu.setMeshAtLevel(0, m3)
        mfu.setMeshAtLevel(-1, m2)
        grpIds3D = DataArrayInt([0, 1])
        grpIds3D.setName("group_3d")
        mfu.setGroupsAtLevel(0, [grpIds3D])  # just to check preservation of 3D group
        mfu.setGroupsAtLevel(-1, [grpIds])
        nNod = m3.getNumberOfNodes()

        c2o2nN = mfu.crackAlong("group", grpMustBeFullyDup=False)
        self.assertEqual({77: {3: 46}}, c2o2nN)

        m3_bis = mfu.getMeshAtLevel(0)
        m3_bis.checkConsistency()
        m2_bis = mfu.getMeshAtLevel(-1)
        m2_bis.checkConsistency()
        self.assertEqual(nNod + 1, mfu.getNumberOfNodes())
        self.assertEqual(nNod + 1, m3_bis.getNumberOfNodes())
        self.assertEqual(nNod + 1, m2_bis.getNumberOfNodes())
        self.assertEqual(
            m3_bis.getCoords()[3].getValues(), m3_bis.getCoords()[nNod:].getValues()
        )
        self.assertEqual([36, 74, 213, 214], mfu.getGroupArr(-1, "group").getValues())
        self.assertEqual([0, 1], mfu.getGroupArr(0, "group_3d").getValues())
        m_bis0 = mfu.getMeshAtLevel(-1)
        m_desc, _, _, _, _ = m_bis0.buildDescendingConnectivity()
        m_bis0.checkDeepEquivalOnSameNodesWith(mfu.getMeshAtLevel(-1), 2, 9.9999999)

    def testBuildInnerBoundary6(self):
        """3D test where the crack has a funny shape with a singular point (i.e. two faces of the M1 group are only connected by one point, not a full segment)
        The singular point was wrongly duplicated.
        """
        # fmt: off
        coo = DataArrayDouble(
            [(-1.38778e-17, 0.226, 0), (-1.38778e-17, -1.38778e-17, 0), (0.226, 0.226, 0), (0.226, -1.38778e-17, 0), (0.452, 0.226, 0), (0.452, -1.38778e-17, 0), (-1.38778e-17, 0.452, 0), (0.226, 0.452, 0), (0.452, 0.452, 0), (-1.38778e-17, 0.226, 0.25), (0.226, 0.226, 0.25), (0.226, -1.38778e-17, 0.25), (-1.38778e-17, -1.38778e-17, 0.25), (-1.38778e-17, 0.226, 0.779375), (0.226, 0.226, 0.779375), (0.226, -1.38778e-17, 0.779375), (-1.38778e-17, -1.38778e-17, 0.779375), (-1.38778e-17, 0.226, 1.30875), (0.226, 0.226, 1.30875), (0.226, -1.38778e-17, 1.30875), (-1.38778e-17, -1.38778e-17, 1.30875), (0.452, 0.226, 0.25), (0.452, -1.38778e-17, 0.25), (0.452, 0.226, 0.779375), (0.452, -1.38778e-17, 0.779375), (0.452, 0.226, 1.30875), (0.452, -1.38778e-17, 1.30875), (-1.38778e-17, 0.452, 0.25), (0.226, 0.452, 0.25), (-1.38778e-17, 0.452, 0.779375), (0.226, 0.452, 0.779375), (-1.38778e-17, 0.452, 1.30875), (0.226, 0.452, 1.30875), (0.452, 0.452, 0.25), (0.452, 0.452, 0.779375), (0.452, 0.452, 1.30875), (0.146, 0.226, 0.779375), (0.146, -1.38778e-17, 0.779375), (0.146, 0.226, 1.30875), (0.146, -1.38778e-17, 1.30875), (0.146, 0.452, 0.779375), (0.146, 0.452, 1.30875)]
        )
        c0 = [18, 0, 2, 3, 1, 9, 10, 11, 12, 18, 9, 10, 11, 12, 13, 36, 37, 16, 18, 13, 36, 37, 16, 17, 38, 39, 20, 18, 2, 4, 5, 3, 10, 21, 22, 11, 18, 10, 21, 22, 11, 14, 23, 24, 15, 18, 14, 23, 24, 15, 18, 25, 26, 19, 18, 6, 7, 2, 0, 27, 28, 10, 9, 18, 27, 28, 10, 9, 29, 40, 36, 13, 18, 29, 40, 36, 13, 31, 41, 38, 17, 18, 7, 8, 4, 2, 28, 33, 21, 10, 18, 28, 33, 21, 10, 30, 34, 23, 14, 18, 30, 34, 23, 14, 32, 35, 25, 18]
        # fmt: on
        cI0 = [0, 9, 18, 27, 36, 45, 54, 63, 72, 81, 90, 99, 108]
        m3 = MEDCouplingUMesh("3D", 3)
        m3.setCoords(coo)
        m3.setConnectivity(DataArrayInt(c0), DataArrayInt(cI0))
        m3.checkConsistency()
        m2, _, _, _, _ = m3.buildDescendingConnectivity()
        grpIds = DataArrayInt([7, 12, 22, 27])
        grpIds.setName("group")
        mfu = MEDFileUMesh()
        mfu.setMeshAtLevel(0, m3)
        mfu.setMeshAtLevel(-1, m2)
        mfu.setGroupsAtLevel(-1, [grpIds])
        nNod = m3.getNumberOfNodes()
        c2o2nN = mfu.crackAlong("group")
        self.assertEqual(
            {
                7: {13: 42, 36: 48},
                8: {13: 42, 36: 48, 17: 44, 38: 49},
                10: {23: 46, 14: 43},
                11: {23: 46, 25: 47, 14: 43, 18: 45},
            },
            c2o2nN,
        )
        m3_bis = mfu.getMeshAtLevel(0)
        m3_bis.checkConsistency()
        m2_bis = mfu.getMeshAtLevel(-1)
        m2_bis.checkConsistency()
        self.assertEqual(nNod + 8, mfu.getNumberOfNodes())
        self.assertEqual(nNod + 8, m3_bis.getNumberOfNodes())
        self.assertEqual(nNod + 8, m2_bis.getNumberOfNodes())
        # self.assertEqual([13, 14, 17, 18, 23, 25, 36, 38], nodesDup.getValues())
        # self.assertEqual(
        #     m3_bis.getCoords()[nodesDup].getValues(),
        #     m3_bis.getCoords()[nNod:].getValues(),
        # )
        # self.assertEqual(set([1, 2, 4, 5]), set(cells1.getValues()))
        # self.assertEqual(set([7, 8, 10, 11]), set(cells2.getValues()))
        self.assertEqual(
            [7, 12, 22, 27, 56, 57, 58, 59], mfu.getGroupArr(-1, "group").getValues()
        )
        m_desc, _, _, _, _ = m3_bis.buildDescendingConnectivity()
        m_desc.checkDeepEquivalOnSameNodesWith(m2_bis, 2, 9.9999)

    def testBuildInnerBoundary7(self):
        """3D test where the crack has another funny shape with another singular point (i.e. two faces of the M1 group are only connected by one point, not a full segment)
        Once the crack is inserted, the cells on either side of the crack do not necessarily form a connex spread zone. This was not properly handled either.
        """
        # fmt: off
        coo = DataArrayDouble(
            [(5, 17, 0), (0, 17, 0), (0, 12, 0), (5, 12, 0), (15, 17, 0), (15, 12, 0), (20, 12, 0), (20, 17, 0), (20, 2, 0), (15, 2, 0), (15, -3, 0), (20, -3, 0), (5, -3, 0), (5, 2, 0), (0, -3, 0), (0, 2, 0), (5, 17, 10), (5, 17, 20), (5, 17, 30), (5, 17, 40), (0, 17, 10), (0, 17, 20), (0, 17, 30), (0, 17, 40), (0, 12, 10), (0, 12, 20), (0, 12, 30), (0, 12, 40), (5, 12, 10), (5, 12, 20), (5, 12, 30), (5, 12, 40), (15, 17, 10), (15, 17, 20), (15, 17, 30), (15, 17, 40), (15, 12, 10), (15, 12, 20), (15, 12, 30), (15, 12, 40), (20, 12, 10), (20, 12, 20), (20, 12, 30), (20, 12, 40), (20, 17, 10), (20, 17, 20), (20, 17, 30), (20, 17, 40), (20, 2, 10), (20, 2, 20), (20, 2, 30), (20, 2, 40), (15, 2, 10), (15, 2, 20), (15, 2, 30), (15, 2, 40), (15, -3, 10), (15, -3, 20), (15, -3, 30), (15, -3, 40), (20, -3, 10), (20, -3, 20), (20, -3, 30), (20, -3, 40), (5, -3, 10), (5, -3, 20), (5, -3, 30), (5, -3, 40), (5, 2, 10), (5, 2, 20), (5, 2, 30), (5, 2, 40), (0, -3, 10), (0, -3, 20), (0, -3, 30), (0, -3, 40), (0, 2, 10), (0, 2, 20), (0, 2, 30), (0, 2, 40), (20, 8, 0), (0, 8, 0), (20, 8, 10), (20, 8, 20), (20, 8, 30), (20, 8, 40), (15, 8, 30), (15, 8, 40), (5, 8, 30), (5, 8, 40), (0, 8, 10), (0, 8, 20), (0, 8, 30), (0, 8, 40)]
        )
        c = DataArrayInt(
            [31, 0, 3, 2, 1, -1, 16, 20, 24, 28, -1, 0, 16, 28, 3, -1, 3, 28, 24, 2, -1, 2, 24, 20, 1, -1, 1, 20, 16, 0, 31, 16, 28, 24, 20, -1, 17, 21, 25, 29, -1, 16, 17, 29, 28, -1, 28, 29, 25, 24, -1, 24, 25, 21, 20, -1, 20, 21, 17, 16, 31, 17, 29, 25, 21, -1, 18, 22, 26, 30, -1, 17, 18, 30, 29, -1, 29, 30, 26, 25, -1, 25, 26, 22, 21, -1, 21, 22, 18, 17, 31, 18, 30, 26, 22, -1, 19, 23, 27, 31, -1, 18, 19, 31, 30, -1, 30, 31, 27, 26, -1, 26, 27, 23, 22, -1, 22, 23, 19, 18, 31, 4, 5, 3, 0, -1, 32, 16, 28, 36, -1, 4, 32, 36, 5, -1, 5, 36, 28, 3, -1, 3, 28, 16, 0, -1, 0, 16, 32, 4, 31, 32, 36, 28, 16, -1, 33, 17, 29, 37, -1, 32, 33, 37, 36, -1, 36, 37, 29, 28, -1, 28, 29, 17, 16, -1, 16, 17, 33, 32, 31, 33, 37, 29, 17, -1, 34, 18, 30, 38, -1, 33, 34, 38, 37, -1, 37, 38, 30, 29, -1, 29, 30, 18, 17, -1, 17, 18, 34, 33, 31, 34, 38, 30, 18, -1, 35, 19, 31, 39, -1, 34, 35, 39, 38, -1, 38, 39, 31, 30, -1, 30, 31, 19, 18, -1, 18, 19, 35, 34, 31, 6, 5, 4, 7, -1, 40, 44, 32, 36, -1, 6, 40, 36, 5, -1, 5, 36, 32, 4, -1, 4, 32, 44, 7, -1, 7, 44, 40, 6, 31, 40, 36, 32, 44, -1, 41, 45, 33, 37, -1, 40, 41, 37, 36, -1, 36, 37, 33, 32, -1, 32, 33, 45, 44, -1, 44, 45, 41, 40, 31, 41, 37, 33, 45, -1, 42, 46, 34, 38, -1, 41, 42, 38, 37, -1, 37, 38, 34, 33, -1, 33, 34, 46, 45, -1, 45, 46, 42, 41, 31, 42, 38, 34, 46, -1, 43, 47, 35, 39, -1, 42, 43, 39, 38, -1, 38, 39, 35, 34, -1, 34, 35, 47, 46, -1, 46, 47, 43, 42, 31, 80, 9, 5, 6, -1, 82, 40, 36, 52, -1, 80, 82, 52, 9, -1, 9, 52, 36, 5, -1, 5, 36, 40, 6, -1, 6, 40, 82, 80, 31, 82, 52, 36, 40, -1, 83, 41, 37, 53, -1, 82, 83, 53, 52, -1, 52, 53, 37, 36, -1, 36, 37, 41, 40, -1, 40, 41, 83, 82, 31, 83, 53, 37, 41, -1, 84, 42, 38, 86, -1, 83, 84, 86, 53, -1, 53, 86, 38, 37, -1, 37, 38, 42, 41, -1, 41, 42, 84, 83, 31, 84, 86, 38, 42, -1, 85, 43, 39, 87, -1, 84, 85, 87, 86, -1, 86, 87, 39, 38, -1, 38, 39, 43, 42, -1, 42, 43, 85, 84, 31, 10, 9, 8, 11, -1, 56, 60, 48, 52, -1, 10, 56, 52, 9, -1, 9, 52, 48, 8, -1, 8, 48, 60, 11, -1, 11, 60, 56, 10, 31, 56, 52, 48, 60, -1, 57, 61, 49, 53, -1, 56, 57, 53, 52, -1, 52, 53, 49, 48, -1, 48, 49, 61, 60, -1, 60, 61, 57, 56, 31, 57, 53, 49, 61, -1, 58, 62, 50, 54, -1, 57, 58, 54, 53, -1, 53, 54, 50, 49, -1, 49, 50, 62, 61, -1, 61, 62, 58, 57, 31, 58, 54, 50, 62, -1, 59, 63, 51, 55, -1, 58, 59, 55, 54, -1, 54, 55, 51, 50, -1, 50, 51, 63, 62, -1, 62, 63, 59, 58, 31, 12, 13, 9, 10, -1, 64, 56, 52, 68, -1, 12, 64, 68, 13, -1, 13, 68, 52, 9, -1, 9, 52, 56, 10, -1, 10, 56, 64, 12, 31, 64, 68, 52, 56, -1, 65, 57, 53, 69, -1, 64, 65, 69, 68, -1, 68, 69, 53, 52, -1, 52, 53, 57, 56, -1, 56, 57, 65, 64, 31, 65, 69, 53, 57, -1, 66, 58, 54, 70, -1, 65, 66, 70, 69, -1, 69, 70, 54, 53, -1, 53, 54, 58, 57, -1, 57, 58, 66, 65, 31, 66, 70, 54, 58, -1, 67, 59, 55, 71, -1, 66, 67, 71, 70, -1, 70, 71, 55, 54, -1, 54, 55, 59, 58, -1, 58, 59, 67, 66, 31, 14, 15, 13, 12, -1, 72, 64, 68, 76, -1, 14, 72, 76, 15, -1, 15, 76, 68, 13, -1, 13, 68, 64, 12, -1, 12, 64, 72, 14, 31, 72, 76, 68, 64, -1, 73, 65, 69, 77, -1, 72, 73, 77, 76, -1, 76, 77, 69, 68, -1, 68, 69, 65, 64, -1, 64, 65, 73, 72, 31, 73, 77, 69, 65, -1, 74, 66, 70, 78, -1, 73, 74, 78, 77, -1, 77, 78, 70, 69, -1, 69, 70, 66, 65, -1, 65, 66, 74, 73, 31, 74, 78, 70, 66, -1, 75, 67, 71, 79, -1, 74, 75, 79, 78, -1, 78, 79, 71, 70, -1, 70, 71, 67, 66, -1, 66, 67, 75, 74, 31, 2, 3, 13, 81, -1, 24, 90, 68, 28, -1, 2, 24, 28, 3, -1, 3, 28, 68, 13, -1, 13, 68, 90, 81, -1, 81, 90, 24, 2, 31, 24, 28, 68, 90, -1, 25, 91, 69, 29, -1, 24, 25, 29, 28, -1, 28, 29, 69, 68, -1, 68, 69, 91, 90, -1, 90, 91, 25, 24, 31, 25, 29, 69, 91, -1, 26, 92, 88, 30, -1, 25, 26, 30, 29, -1, 29, 30, 88, 69, -1, 69, 88, 92, 91, -1, 91, 92, 26, 25, 31, 26, 30, 88, 92, -1, 27, 93, 89, 31, -1, 26, 27, 31, 30, -1, 30, 31, 89, 88, -1, 88, 89, 93, 92, -1, 92, 93, 27, 26, 31, 13, 3, 5, 9, -1, 68, 52, 36, 28, -1, 13, 68, 28, 3, -1, 3, 28, 36, 5, -1, 5, 36, 52, 9, -1, 9, 52, 68, 13, 31, 68, 28, 36, 52, -1, 69, 53, 37, 29, -1, 68, 69, 29, 28, -1, 28, 29, 37, 36, -1, 36, 37, 53, 52, -1, 52, 53, 69, 68, 31, 69, 29, 37, 53, -1, 88, 86, 38, 30, -1, 69, 88, 30, 29, -1, 29, 30, 38, 37, -1, 37, 38, 86, 53, -1, 53, 86, 88, 69, 31, 88, 30, 38, 86, -1, 89, 87, 39, 31, -1, 88, 89, 31, 30, -1, 30, 31, 39, 38, -1, 38, 39, 87, 86, -1, 86, 87, 89, 88]
        )
        cI = DataArrayInt(
            [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570, 600, 630, 660, 690, 720, 750, 780, 810, 840, 870, 900, 930, 960, 990, 1020, 1050, 1080]
        )
        # fmt: on
        m3 = MEDCouplingUMesh("box", 3)
        m3.setCoords(coo)
        m3.setConnectivity(c, cI)
        m3.checkConsistency()
        m2, _, _, _, _ = m3.buildDescendingConnectivity()
        grpIds = DataArrayInt([2, 7, 12, 17, 95, 99, 103, 107, 129, 133, 137, 141])
        grpIds.setName("group")
        mfu = MEDFileUMesh()
        mfu.setMeshAtLevel(0, m3)
        mfu.setMeshAtLevel(-1, m2)
        mfu.setGroupsAtLevel(-1, [grpIds])
        nNod = m3.getNumberOfNodes()

        c2o2nN = mfu.crackAlong("group")

        self.assertEqual(
            {
                4: {0: 94, 3: 95, 16: 99, 28: 103},
                5: {16: 99, 17: 100, 28: 103, 29: 104},
                6: {17: 100, 18: 101, 29: 104, 30: 105},
                7: {18: 101, 19: 102, 30: 105, 31: 106},
                24: {12: 96, 13: 97, 64: 107, 68: 111},
                25: {64: 107, 65: 108, 68: 111, 69: 113},
                26: {65: 108, 66: 109, 69: 113, 70: 115},
                27: {66: 109, 67: 110, 70: 115, 71: 116},
                28: {13: 98, 68: 112},
                29: {68: 112, 69: 114},
                30: {69: 114},
                32: {3: 95, 28: 103},
                33: {28: 103, 29: 104},
                34: {29: 104, 30: 105, 88: 117},
                35: {30: 105, 31: 106, 88: 117, 89: 118},
            },
            c2o2nN,
        )
        m3_bis = mfu.getMeshAtLevel(0)
        m3_bis.checkConsistency()
        m2_bis = mfu.getMeshAtLevel(-1)
        m2_bis.checkConsistency()
        self.assertEqual(nNod + 25, mfu.getNumberOfNodes())
        self.assertEqual(nNod + 25, m3_bis.getNumberOfNodes())
        self.assertEqual(nNod + 25, m2_bis.getNumberOfNodes())
        self.assertEqual(
            [2  , 7  , 12 , 17 , 95 , 99 , 103, 107, 129, 133, 137, 141,
             151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162],
            mfu.getGroupArr(-1, "group").getValues(),
        )  # fmt: skip
        m_desc, _, _, _, _ = m3_bis.buildDescendingConnectivity()
        m_desc.checkDeepEquivalOnSameNodesWith(m2_bis, 2, 9.9999)

    # def testBuildInnerBoundary8(self):
    #     """3D test where the crack leaves 'naked' cells. If we call a 'close-to-crack cell' a cell which shares a face with the M1 group,
    #     a 'naked cell' is a cell that has some node duplicated, but which do not share any face with a 'close-to-crack cell'. In this case
    #     it is tricky to decide whether this cell should be renumbered or not ...
    #     Warning: on the mesh below some points have already been doubled by a previous cut.
    #     """
    #     # fmt: off
    #     coo = DataArrayDouble(
    #         [(0, 15, 0), (0, 5, 0), (3, 5, 0), (5, 5, 0), (5, 15, 0), (5, 20, 0), (0, 20, 0), (15, 20, 0), (15, 15, 0), (20, 15, 0), (20, 20, 0), (20, 5, 0), (15, 5, 0), (15, 0, 0), (20, 0, 0), (5, -1.60551e-25, 0), (5, 3, 0), (3, 0, 0), (3, 3, 0), (0, 0, 0), (0, 3, 0), (0, 15, 10), (0, 15, 20), (0, 15, 30), (0, 15, 40), (0, 5, 10), (0, 5, 20), (0, 5, 30), (0, 5, 40), (3, 5, 10), (3, 5, 20), (3, 5, 30), (3, 5, 40), (5, 5, 10), (5, 5, 20), (5, 5, 30), (5, 5, 40), (5, 15, 10), (5, 15, 20), (5, 15, 30), (5, 15, 40), (5, 20, 10), (5, 20, 20), (5, 20, 30), (5, 20, 40), (0, 20, 10), (0, 20, 20), (0, 20, 30), (0, 20, 40), (15, 20, 10), (15, 20, 20), (15, 20, 30), (15, 20, 40), (15, 15, 10), (15, 15, 20), (15, 15, 30), (15, 15, 40), (20, 15, 10), (20, 15, 20), (20, 15, 30), (20, 15, 40), (20, 20, 10), (20, 20, 20), (20, 20, 30), (20, 20, 40), (20, 5, 10), (20, 5, 20), (20, 5, 30), (20, 5, 40), (15, 5, 10), (15, 5, 20), (15, 5, 30), (15, 5, 40), (15, 0, 10), (15, 0, 20), (15, 0, 30), (15, 0, 40), (20, 0, 10), (20, 0, 20), (20, 0, 30), (20, 0, 40), (5, -1.60551e-25, 10), (5, -1.60551e-25, 20), (5, -1.60551e-25, 30), (5, -1.60551e-25, 40), (5, 3, 10), (5, 3, 20), (5, 3, 30), (5, 3, 40), (3, 0, 10), (3, 0, 20), (3, 0, 30), (3, 0, 40), (3, 3, 10), (3, 3, 20), (3, 3, 30), (3, 3, 40), (0, 0, 10), (0, 0, 20), (0, 0, 30), (0, 0, 40), (0, 3, 10), (0, 3, 20), (0, 3, 30), (0, 3, 40), (0, 9, 0), (3, 9, 0), (20, 9, 0), (0, 9, 10), (0, 9, 20), (0, 9, 30), (0, 9, 40), (3, 9, 10), (3, 9, 20), (3, 9, 30), (3, 9, 40), (5, 9, 30), (5, 9, 40), (20, 9, 10), (20, 9, 20), (20, 9, 30), (20, 9, 40), (15, 9, 30), (15, 9, 40), (0, 15, 0), (20, 15, 0), (0, 15, 10), (0, 15, 20), (0, 15, 30), (0, 15, 40), (5, 15, 30), (5, 15, 40), (15, 15, 30), (15, 15, 40), (20, 15, 10), (20, 15, 20), (20, 15, 30), (20, 15, 40)]
    #     )
    #     c = DataArrayInt(
    #         [31, 5, 4, 124, 6, -1, 41, 45, 126, 37, -1, 5, 41, 37, 4, -1, 4, 37, 126, 124, -1, 124, 126, 45, 6, -1, 6, 45, 41, 5, 31, 41, 37, 126, 45, -1, 42, 46, 127, 38, -1, 41, 42, 38, 37, -1, 37, 38, 127, 126, -1, 126, 127, 46, 45, -1, 45, 46, 42, 41, 31, 42, 38, 127, 46, -1, 43, 47, 128, 130, -1, 42, 43, 130, 38, -1, 38, 130, 128, 127, -1, 127, 128, 47, 46, -1, 46, 47, 43, 42, 31, 43, 130, 128, 47, -1, 44, 48, 129, 131, -1, 43, 44, 131, 130, -1, 130, 131, 129, 128, -1, 128, 129, 48, 47, -1, 47, 48, 44, 43, 31, 7, 8, 4, 5, -1, 49, 41, 37, 53, -1, 7, 49, 53, 8, -1, 8, 53, 37, 4, -1, 4, 37, 41, 5, -1, 5, 41, 49, 7, 31, 49, 53, 37, 41, -1, 50, 42, 38, 54, -1, 49, 50, 54, 53, -1, 53, 54, 38, 37, -1, 37, 38, 42, 41, -1, 41, 42, 50, 49, 31, 50, 54, 38, 42, -1, 51, 43, 130, 132, -1, 50, 51, 132, 54, -1, 54, 132, 130, 38, -1, 38, 130, 43, 42, -1, 42, 43, 51, 50, 31, 51, 132, 130, 43, -1, 52, 44, 131, 133, -1, 51, 52, 133, 132, -1, 132, 133, 131, 130, -1, 130, 131, 44, 43, -1, 43, 44, 52, 51, 31, 125, 8, 7, 10, -1, 134, 61, 49, 53, -1, 125, 134, 53, 8, -1, 8, 53, 49, 7, -1, 7, 49, 61, 10, -1, 10, 61, 134, 125, 31, 134, 53, 49, 61, -1, 135, 62, 50, 54, -1, 134, 135, 54, 53, -1, 53, 54, 50, 49, -1, 49, 50, 62, 61, -1, 61, 62, 135, 134, 31, 135, 54, 50, 62, -1, 136, 63, 51, 132, -1, 135, 136, 132, 54, -1, 54, 132, 51, 50, -1, 50, 51, 63, 62, -1, 62, 63, 136, 135, 31, 136, 132, 51, 63, -1, 137, 64, 52, 133, -1, 136, 137, 133, 132, -1, 132, 133, 52, 51, -1, 51, 52, 64, 63, -1, 63, 64, 137, 136, 31, 107, 12, 8, 9, -1, 118, 57, 53, 69, -1, 107, 118, 69, 12, -1, 12, 69, 53, 8, -1, 8, 53, 57, 9, -1, 9, 57, 118, 107, 31, 118, 69, 53, 57, -1, 119, 58, 54, 70, -1, 118, 119, 70, 69, -1, 69, 70, 54, 53, -1, 53, 54, 58, 57, -1, 57, 58, 119, 118, 31, 119, 70, 54, 58, -1, 120, 59, 55, 122, -1, 119, 120, 122, 70, -1, 70, 122, 55, 54, -1, 54, 55, 59, 58, -1, 58, 59, 120, 119, 31, 120, 122, 55, 59, -1, 121, 60, 56, 123, -1, 120, 121, 123, 122, -1, 122, 123, 56, 55, -1, 55, 56, 60, 59, -1, 59, 60, 121, 120, 31, 13, 12, 11, 14, -1, 73, 77, 65, 69, -1, 13, 73, 69, 12, -1, 12, 69, 65, 11, -1, 11, 65, 77, 14, -1, 14, 77, 73, 13, 31, 73, 69, 65, 77, -1, 74, 78, 66, 70, -1, 73, 74, 70, 69, -1, 69, 70, 66, 65, -1, 65, 66, 78, 77, -1, 77, 78, 74, 73, 31, 74, 70, 66, 78, -1, 75, 79, 67, 71, -1, 74, 75, 71, 70, -1, 70, 71, 67, 66, -1, 66, 67, 79, 78, -1, 78, 79, 75, 74, 31, 75, 71, 67, 79, -1, 76, 80, 68, 72, -1, 75, 76, 72, 71, -1, 71, 72, 68, 67, -1, 67, 68, 80, 79, -1, 79, 80, 76, 75, 31, 17, 18, 16, 15, -1, 89, 81, 85, 93, -1, 17, 89, 93, 18, -1, 18, 93, 85, 16, -1, 16, 85, 81, 15, -1, 15, 81, 89, 17, 31, 89, 93, 85, 81, -1, 90, 82, 86, 94, -1, 89, 90, 94, 93, -1, 93, 94, 86, 85, -1, 85, 86, 82, 81, -1, 81, 82, 90, 89, 31, 90, 94, 86, 82, -1, 91, 83, 87, 95, -1, 90, 91, 95, 94, -1, 94, 95, 87, 86, -1, 86, 87, 83, 82, -1, 82, 83, 91, 90, 31, 91, 95, 87, 83, -1, 92, 84, 88, 96, -1, 91, 92, 96, 95, -1, 95, 96, 88, 87, -1, 87, 88, 84, 83, -1, 83, 84, 92, 91, 31, 19, 20, 18, 17, -1, 97, 89, 93, 101, -1, 19, 97, 101, 20, -1, 20, 101, 93, 18, -1, 18, 93, 89, 17, -1, 17, 89, 97, 19, 31, 97, 101, 93, 89, -1, 98, 90, 94, 102, -1, 97, 98, 102, 101, -1, 101, 102, 94, 93, -1, 93, 94, 90, 89, -1, 89, 90, 98, 97, 31, 98, 102, 94, 90, -1, 99, 91, 95, 103, -1, 98, 99, 103, 102, -1, 102, 103, 95, 94, -1, 94, 95, 91, 90, -1, 90, 91, 99, 98, 31, 99, 103, 95, 91, -1, 100, 92, 96, 104, -1, 99, 100, 104, 103, -1, 103, 104, 96, 95, -1, 95, 96, 92, 91, -1, 91, 92, 100, 99, 31, 1, 2, 18, 20, -1, 25, 101, 93, 29, -1, 1, 25, 29, 2, -1, 2, 29, 93, 18, -1, 18, 93, 101, 20, -1, 20, 101, 25, 1, 31, 25, 29, 93, 101, -1, 26, 102, 94, 30, -1, 25, 26, 30, 29, -1, 29, 30, 94, 93, -1, 93, 94, 102, 101, -1, 101, 102, 26, 25, 31, 26, 30, 94, 102, -1, 27, 103, 95, 31, -1, 26, 27, 31, 30, -1, 30, 31, 95, 94, -1, 94, 95, 103, 102, -1, 102, 103, 27, 26, 31, 27, 31, 95, 103, -1, 28, 104, 96, 32, -1, 27, 28, 32, 31, -1, 31, 32, 96, 95, -1, 95, 96, 104, 103, -1, 103, 104, 28, 27, 31, 3, 4, 8, 12, -1, 33, 69, 53, 37, -1, 3, 33, 37, 4, -1, 4, 37, 53, 8, -1, 8, 53, 69, 12, -1, 12, 69, 33, 3, 31, 33, 37, 53, 69, -1, 34, 70, 54, 38, -1, 33, 34, 38, 37, -1, 37, 38, 54, 53, -1, 53, 54, 70, 69, -1, 69, 70, 34, 33, 31, 34, 38, 54, 70, -1, 116, 122, 55, 39, -1, 34, 116, 39, 38, -1, 38, 39, 55, 54, -1, 54, 55, 122, 70, -1, 70, 122, 116, 34, 31, 116, 39, 55, 122, -1, 117, 123, 56, 40, -1, 116, 117, 40, 39, -1, 39, 40, 56, 55, -1, 55, 56, 123, 122, -1, 122, 123, 117, 116, 31, 16, 18, 2, 3, -1, 85, 33, 29, 93, -1, 16, 85, 93, 18, -1, 18, 93, 29, 2, -1, 2, 29, 33, 3, -1, 3, 33, 85, 16, 31, 85, 93, 29, 33, -1, 86, 34, 30, 94, -1, 85, 86, 94, 93, -1, 93, 94, 30, 29, -1, 29, 30, 34, 33, -1, 33, 34, 86, 85, 31, 86, 94, 30, 34, -1, 87, 35, 31, 95, -1, 86, 87, 95, 94, -1, 94, 95, 31, 30, -1, 30, 31, 35, 34, -1, 34, 35, 87, 86, 31, 87, 95, 31, 35, -1, 88, 36, 32, 96, -1, 87, 88, 96, 95, -1, 95, 96, 32, 31, -1, 31, 32, 36, 35, -1, 35, 36, 88, 87, 31, 4, 3, 106, 105, 0, -1, 37, 21, 108, 112, 33, -1, 3, 4, 37, 33, -1, 106, 3, 33, 112, -1, 105, 106, 112, 108, -1, 0, 105, 108, 21, -1, 4, 0, 21, 37, 31, 37, 33, 112, 108, 21, -1, 38, 22, 109, 113, 34, -1, 33, 37, 38, 34, -1, 112, 33, 34, 113, -1, 108, 112, 113, 109, -1, 21, 108, 109, 22, -1, 37, 21, 22, 38, 31, 38, 34, 113, 109, 22, -1, 39, 23, 110, 114, 116, -1, 34, 38, 39, 116, -1, 113, 34, 116, 114, -1, 109, 113, 114, 110, -1, 22, 109, 110, 23, -1, 38, 22, 23, 39, 31, 39, 116, 114, 110, 23, -1, 40, 24, 111, 115, 117, -1, 116, 39, 40, 117, -1, 114, 116, 117, 115, -1, 110, 114, 115, 111, -1, 23, 110, 111, 24, -1, 39, 23, 24, 40, 31, 16, 3, 12, 13, 15, -1, 85, 81, 73, 69, 33, -1, 3, 16, 85, 33, -1, 12, 3, 33, 69, -1, 13, 12, 69, 73, -1, 15, 13, 73, 81, -1, 16, 15, 81, 85, 31, 85, 33, 69, 73, 81, -1, 86, 82, 74, 70, 34, -1, 33, 85, 86, 34, -1, 69, 33, 34, 70, -1, 73, 69, 70, 74, -1, 81, 73, 74, 82, -1, 85, 81, 82, 86, 31, 86, 34, 70, 74, 82, -1, 87, 83, 75, 71, 35, -1, 34, 86, 87, 35, -1, 70, 34, 35, 71, -1, 74, 70, 71, 75, -1, 82, 74, 75, 83, -1, 86, 82, 83, 87, 31, 87, 35, 71, 75, 83, -1, 88, 84, 76, 72, 36, -1, 35, 87, 88, 36, -1, 71, 35, 36, 72, -1, 75, 71, 72, 76, -1, 83, 75, 76, 84, -1, 87, 83, 84, 88]
    #     )
    #     cI = DataArrayInt(
    #         [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570, 600, 630, 660, 690, 720, 750, 780, 810, 840, 870, 900, 930, 960, 990, 1020, 1050, 1080, 1110, 1140, 1170, 1200, 1237, 1274, 1311, 1348, 1385, 1422, 1459, 1496]
    #     )
    #     # fmt: on
    #     m3 = MEDCouplingUMesh("box", 3)
    #     m3.setCoords(coo)
    #     m3.setConnectivity(c, cI)
    #     m3.mergeNodes(1e-12)
    #     m3.checkConsistency()
    #     m2, _, _, _, _ = m3.buildDescendingConnectivity()
    #     # grpIds = DataArrayInt(
    #     #     [2, 7, 12, 17, 101, 106, 111, 116, 160, 164, 170, 173, 176, 179]
    #     # )
    #     grpIds = DataArrayInt(
    #         [16, 149, 12, 146, 8, 143, 137, 108, 134, 103, 131, 98, 128, 93]
    #     )
    #     grpIds.setName("group")
    #     mfu = MEDFileUMesh()
    #     mfu.setMeshAtLevel(0, m3)
    #     mfu.setMeshAtLevel(-1, m2)
    #     mfu.setGroupsAtLevel(-1, [grpIds])
    #     nNod = m3.getNumberOfNodes()

    #     mfu.write("case8dupnode_in.med", 2)

    #     c2o2nN = mfu.crackAlong("group", False)

    #     mfu.write("case8dupnode_out.med", 2)

    #     mfu.OpenCrack(c2o2nN)

    #     mfu.write("case8dupnode_cracked.med", 2)

    #     self.assertEqual(
    #         {
    #             44: {85: 166, 15: 163, 81: 159, 16: 158},
    #             6: {130: 162, 43: 161, 38: 160, 42: 157},
    #             41: {38: 160, 37: 154},
    #             47: {83: 146, 84: 143, 87: 148, 35: 139, 88: 145, 36: 138},
    #             46: {86: 150, 82: 149, 83: 146, 87: 148, 35: 139},
    #             7: {131: 165, 44: 164, 130: 162, 43: 161},
    #             36: {33: 155, 3: 140},
    #             43: {116: 147, 117: 144, 39: 142, 40: 141},
    #             33: {38: 160, 37: 154},
    #             4: {4: 156, 37: 154, 41: 153, 5: 152},
    #             42: {38: 160, 116: 147, 39: 142},
    #             45: {85: 166, 81: 159, 86: 150, 82: 149},
    #             32: {4: 156, 37: 154},
    #             37: {33: 155, 34: 151},
    #             38: {34: 151},
    #             34: {38: 160},
    #             5: {38: 160, 42: 157, 37: 154, 41: 153},
    #             40: {4: 156, 37: 154},
    #         },
    #         c2o2nN,
    #     )
    #     m3_bis = mfu.getMeshAtLevel(0)
    #     m3_bis.checkConsistency()
    #     m2_bis = mfu.getMeshAtLevel(-1)
    #     m2_bis.checkConsistency()
    #     self.assertEqual(nNod + 23, mfu.getNumberOfNodes())
    #     self.assertEqual(nNod + 23, m3_bis.getNumberOfNodes())
    #     self.assertEqual(nNod + 23, m2_bis.getNumberOfNodes())
    #     self.assertEqual(
    #         [5, 15, 16, 35, 36, 39, 40, 41, 42, 43, 44, 81, 82, 83, 84, 85, 86, 87, 88, 116, 117, 130, 131],
    #         nodesDup.getValues(),
    #     )  # fmt: skip
    #     self.assertEqual(
    #         m3_bis.getCoords()[nodesDup].getValues(),
    #         m3_bis.getCoords()[nNod:].getValues(),
    #     )
    #     self.assertEqual(
    #         set([0, 1, 2, 3, 20, 21, 22, 23, 34, 35, 36, 37, 38, 39]),
    #         set(cells1.getValues()),
    #     )
    #     self.assertEqual(
    #         set([4, 5, 6, 7, 42, 43, 44, 45, 46, 47]), set(cells2.getValues())
    #     )
    #     self.assertEqual(
    #         [2, 7, 12, 17, 101, 106, 111, 116, 160, 164, 170, 173, 176, 179],
    #         mfu.getGroupArr(-1, "group").getValues(),
    #     )
    #     self.assertEqual(
    #         [212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225],
    #         mfu.getGroupArr(-1, "group_dup").getValues(),
    #     )  # fmt: skip
    #     m_desc, _, _, _, _ = m3_bis.buildDescendingConnectivity()
    #     m_desc.checkDeepEquivalOnSameNodesWith(m2_bis, 2, 9.9999)

    def testBuildInnerBoundary9(self):
        """3D test where the crack is performed so that two non-connex parts are found facing one single connex part on the other side
        of the crack.
        """
        # fmt: off
        coo = DataArrayDouble(
            [(0, 4.6, 0), (3, 4.6, 0), (5, 4.6, 0), (15, 4.6, 0), (15, 0, 0), (5, -1.60551e-25, 0), (5, 3, 0), (3, 0, 0), (3, 3.8, 0), (0, 0, 0), (0, 3.8, 0), (0, 4.6, 10), (0, 4.6, 20), (3, 4.6, 10), (3, 4.6, 20), (5, 4.6, 10), (5, 4.6, 20), (15, 4.6, 10), (15, 4.6, 20), (15, 0, 10), (15, 0, 20), (5, -1.60551e-25, 10), (5, -1.60551e-25, 20), (5, 3, 10), (5, 3, 20), (3, 0, 10), (3, 0, 20), (3, 3.8, 10), (3, 3.8, 20), (0, 0, 10), (0, 0, 20), (0, 3.8, 10), (0, 3.8, 20), (3, 3, 0), (0, 3, 0), (3, 3, 10), (3, 3, 20), (0, 3, 10), (0, 3, 20), ]
        )
        c = DataArrayInt(
            [31, 7, 33, 6, 5, -1, 25, 21, 23, 35, -1, 7, 25, 35, 33, -1, 33, 35, 23, 6, -1, 6, 23, 21, 5, -1, 5, 21, 25, 7, 31, 25, 35, 23, 21, -1, 26, 22, 24, 36, -1, 25, 26, 36, 35, -1, 35, 36, 24, 23, -1, 23, 24, 22, 21, -1, 21, 22, 26, 25, 31, 9, 34, 33, 7, -1, 29, 25, 35, 37, -1, 9, 29, 37, 34, -1, 34, 37, 35, 33, -1, 33, 35, 25, 7, -1, 7, 25, 29, 9, 31, 29, 37, 35, 25, -1, 30, 26, 36, 38, -1, 29, 30, 38, 37, -1, 37, 38, 36, 35, -1, 35, 36, 26, 25, -1, 25, 26, 30, 29, 31, 0, 1, 8, 10, -1, 11, 31, 27, 13, -1, 0, 11, 13, 1, -1, 1, 13, 27, 8, -1, 8, 27, 31, 10, -1, 10, 31, 11, 0, 31, 11, 13, 27, 31, -1, 12, 32, 28, 14, -1, 11, 12, 14, 13, -1, 13, 14, 28, 27, -1, 27, 28, 32, 31, -1, 31, 32, 12, 11, 31, 6, 8, 1, 2, -1, 23, 15, 13, 27, -1, 6, 23, 27, 8, -1, 8, 27, 13, 1, -1, 1, 13, 15, 2, -1, 2, 15, 23, 6, 31, 23, 27, 13, 15, -1, 24, 16, 14, 28, -1, 23, 24, 28, 27, -1, 27, 28, 14, 13, -1, 13, 14, 16, 15, -1, 15, 16, 24, 23, 31, 6, 2, 3, 4, 5, -1, 23, 21, 19, 17, 15, -1, 2, 6, 23, 15, -1, 3, 2, 15, 17, -1, 4, 3, 17, 19, -1, 5, 4, 19, 21, -1, 6, 5, 21, 23, 31, 23, 15, 17, 19, 21, -1, 24, 22, 20, 18, 16, -1, 15, 23, 24, 16, -1, 17, 15, 16, 18, -1, 19, 17, 18, 20, -1, 21, 19, 20, 22, -1, 23, 21, 22, 24]
        )
        # fmt: on
        cI = DataArrayInt([0, 30, 60, 90, 120, 150, 180, 210, 240, 277, 314])
        m3 = MEDCouplingUMesh("box", 3)
        m3.setCoords(coo)
        m3.setConnectivity(c, cI)
        m3.checkConsistency()
        m2, _, _, _, _ = m3.buildDescendingConnectivity()
        grpIds = DataArrayInt([4, 9, 35, 39])
        grpIds.setName("group")
        mfu = MEDFileUMesh()
        mfu.setMeshAtLevel(0, m3)
        mfu.setMeshAtLevel(-1, m2)
        mfu.setGroupsAtLevel(-1, [grpIds])
        m2, _, _, _, _ = m3.buildDescendingConnectivity()
        grpIds = DataArrayInt([4, 9, 35, 39])
        grpIds.setName("group")
        mfu = MEDFileUMesh()
        mfu.setMeshAtLevel(0, m3)
        mfu.setMeshAtLevel(-1, m2)
        mfu.setGroupsAtLevel(-1, [grpIds])
        nNod = m3.getNumberOfNodes()

        c2o2nN = mfu.crackAlong("group")

        self.assertEqual(
            {
                6: {6: 41, 23: 47},
                7: {23: 47, 24: 49},
                8: {6: 42, 23: 48, 21: 45, 5: 40, 2: 39, 15: 43},
                9: {23: 48, 21: 45, 24: 50, 22: 46, 15: 43, 16: 44},
            },
            c2o2nN,
        )

        m3_bis = mfu.getMeshAtLevel(0)
        m3_bis.checkConsistency()
        m2_bis = mfu.getMeshAtLevel(-1)
        m2_bis.checkConsistency()
        self.assertEqual(nNod + 12, mfu.getNumberOfNodes())
        self.assertEqual(nNod + 12, m3_bis.getNumberOfNodes())
        self.assertEqual(nNod + 12, m2_bis.getNumberOfNodes())
        # self.assertEqual([2, 5, 6, 15, 16, 21, 22, 23, 24], nodesDup.getValues())
        # self.assertEqual(
        #     m3_bis.getCoords()[nodesDup].getValues(),
        #     m3_bis.getCoords()[nNod:].getValues(),
        # )
        # self.assertEqual(set([0, 1, 6, 7]), set(cells1.getValues()))
        # self.assertEqual(set([8, 9]), set(cells2.getValues()))
        self.assertEqual(
            [4, 9, 35, 39, 49, 50, 51, 52], mfu.getGroupArr(-1, "group").getValues()
        )
        m_desc, _, _, _, _ = m3_bis.buildDescendingConnectivity()
        m_desc.checkDeepEquivalOnSameNodesWith(m2_bis, 2, 9.9999)

    def testBuildInnerBoundary10(self):
        """2D tests where some cells are touching the M1 group with just a node, and are **not** neighbor
        with any cells touching the M1 group by a face.
        """
        # fmt: off
        coo = DataArrayDouble(
            [0.0, 0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 1.0, 3.0, 1.0, 0.0, 2.0, 1.0, 2.0, 2.0, 2.0, 3.0, 2.0, 2.5, 0.0, 3.0, 0.5, ],
            14,
            2,
        )
        c = DataArrayInt(
            [3, 12, 2, 6, 3, 3, 12, 6, 3, 13, 3, 6, 3, 7, 13, 6, 4, 1, 0, 4, 5, 4, 2, 1, 5, 6, 4, 5, 4, 8, 9, 4, 6, 5, 9, 10, 4, 7, 6, 10, 11]
        )
        # fmt: on
        cI = DataArrayInt([0, 4, 8, 12, 16, 21, 26, 31, 36, 41])
        m2 = MEDCouplingUMesh("mesh", 2)
        m2.setCoords(coo)
        m2.setConnectivity(c, cI)
        m2.checkConsistency()
        m1, _, _, _, _ = m2.buildDescendingConnectivity()
        grpIds = DataArrayInt([8, 14])
        grpIds.setName("group")
        mfu = MEDFileUMesh()
        mfu.setMeshAtLevel(0, m2)
        mfu.setMeshAtLevel(-1, m1)
        mfu.setGroupsAtLevel(-1, [grpIds])
        nNod = m2.getNumberOfNodes()

        c2o2nN = mfu.crackAlong("group")

        self.assertEqual({7: {6: 14}, 8: {6: 14, 7: 15}}, c2o2nN)

        m2_bis = mfu.getMeshAtLevel(0)
        m2_bis.checkConsistency()
        m1_bis = mfu.getMeshAtLevel(-1)
        m1_bis.checkConsistency()
        self.assertEqual(nNod + 2, mfu.getNumberOfNodes())
        self.assertEqual(nNod + 2, m2_bis.getNumberOfNodes())
        self.assertEqual(nNod + 2, m1_bis.getNumberOfNodes())
        # self.assertEqual([6, 7], nodesDup.getValues())
        # self.assertEqual([2.0, 1.0, 3.0, 1.0], m2_bis.getCoords()[nNod:].getValues())
        # self.assertEqual(set([0, 1, 2, 3, 5]), set(cells1.getValues()))
        # self.assertEqual(set([7, 8]), set(cells2.getValues()))
        self.assertEqual([8, 14, 22, 23], mfu.getGroupArr(-1, "group").getValues())


if __name__ == "__main__":
    unittest.main()
