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

    def test_bug26801_sortToHaveConsecutivePairs(self):
        """
        EDF26801 : fix bug in sortToHaveConsecutivePairs + sortToHaveConsecutivePairs
        """
        import itertools

        data = [(11, 10), (10, 14), (15, 13), (12, 11), (13, 12)]
        ref = mc.DataArrayInt([15, 13, 12, 11, 10, 14])

        def check(permut):
            a = mc.DataArrayInt(data)[list(permut)]
            a.sortToHaveConsecutivePairs()
            self.assertTrue(a.fromLinkedListOfPairToList().isEqual(ref))

        for permutation in itertools.permutations(range(5)):
            check(permutation)

        #
        permutation = mc.DataArrayInt([4, 2, 1, 3, 5, 0])
        data2 = mc.DataArrayInt(
            [(11, 10), (10, 14), (12, 11), (21, 20), (20, 24), (22, 21)]
        )
        data2 = data2[permutation]
        # 2 chains -> raises
        self.assertRaises(mc.InterpKernelException, data2.sortToHaveConsecutivePairs)
        ret = data2.splitPairsIntoChains()
        self.assertEqual(len(ret), 2)
        self.assertTrue(mc.DataArrayInt([(22, 21), (21, 20), (20, 24)]).isEqual(ret[0]))
        self.assertTrue(mc.DataArrayInt([(12, 11), (11, 10), (10, 14)]).isEqual(ret[1]))
        # 2 chains with cycle
        data3 = mc.DataArrayInt(
            [(11, 10), (10, 12), (12, 11), (21, 20), (20, 22), (22, 21)]
        )
        data3 = data3[permutation]
        ret = data3.splitPairsIntoChains()
        self.assertEqual(len(ret), 2)
        self.assertTrue(mc.DataArrayInt([(21, 20), (20, 22), (22, 21)]).isEqual(ret[0]))
        self.assertTrue(mc.DataArrayInt([(10, 12), (12, 11), (11, 10)]).isEqual(ret[1]))

    def test_bug35017_distanceToPoint(self):
        """
        [EDF35017] : Fix distanceToPoint
        """
        # fmt: off
        from math import sqrt, pi
        def rotate( mesh, ptsMesh, rotationInRad ):
            retMesh = mesh.deepCopy()
            retMesh.rotate( [0,0,0], [0,0,1], rotationInRad )
            retPtsMesh = ptsMesh.deepCopy()
            retPtsMesh.rotate( [0,0,0], [0,0,1], rotationInRad )
            return retMesh, retPtsMesh

        coo = mc.DataArrayDouble( [(0,0),(2,0),(2,1),(0,1)] )
        n2o = mc.DataArrayInt( [0,2,3,1] )
        o2n = n2o.invertArrayN2O2O2N(4)
        self.assertTrue( o2n.isEqual( mc.DataArrayInt( [0,3,1,2] ) ) )
        m = mc.MEDCouplingUMesh("",2)
        m.setCoords( coo[n2o] )
        m.allocateCells()
        m.insertNextCell( mc.NORM_QUAD4, o2n.getValues() )
        m.changeSpaceDimension(3,0.)
        # First testing pts outside polygon
        ptsOut = mc.DataArrayDouble( [(-1,-1), (0,-1), (1.8,-1), (3,-1), (3,0.7), (3,2), (0.8, 2), (-1,2), (-1, 0.9) ] )
        ptsOut = ptsOut.changeNbOfComponents(3,0.)
        pts = mc.MEDCouplingUMesh.Build0DMeshFromCoords( ptsOut )
        rotations = [ 0., pi/6, pi/3, pi/2, 3*pi/4, pi, 5*pi/4, 11 * pi /6 ]
        ref_ptsOut = mc.DataArrayDouble([ sqrt(2.0), 1.0, 1.0, sqrt(2.0), 1.0, sqrt(2.0), 1.0, sqrt(2.0), 1.0 ])
        for curRot in rotations:
            mCpy, ptsCpy = rotate( m, pts, curRot )
            self.assertTrue( mCpy.distanceToPoints( ptsCpy.getCoords() )[0].isEqual( ref_ptsOut, 1e-12 ) )
        # Now manage Z
        for zLev in [2.0, -3.0]:
            pts.getCoords()[:,2] = zLev
            for curRot in rotations:
                mCpy, ptsCpy = rotate( m, pts, curRot )
                mCpy.distanceToPoints( ptsCpy.getCoords() )
                curRef = ( ref_ptsOut**2 + zLev**2 )**0.5
                self.assertTrue( mCpy.distanceToPoints( ptsCpy.getCoords() )[0].isEqual( curRef, 1e-12 ) )
        # Now testing points inside polygon
        ptsIn = mc.DataArrayDouble( [(0.1,0.1), (1.9,0.1), (1.0, 0.5), (0.1,0.9), (0,0),(2,0),(2,1),(0,1) ] )
        ptsIn = ptsIn.changeNbOfComponents(3,0.)
        pts = mc.MEDCouplingUMesh.Build0DMeshFromCoords( ptsIn )
        for curRot in rotations:
            mCpy, ptsCpy = rotate( m, pts, curRot )
            self.assertTrue( mCpy.distanceToPoints( ptsCpy.getCoords() )[0].isUniform( 0, 1e-12 ) )
        # Now manage Z
        for zLev in [2.0, -3.0]:
            pts.getCoords()[:,2] = zLev
            for curRot in rotations:
                mCpy, ptsCpy = rotate( m, pts, curRot )
                mCpy.distanceToPoints( ptsCpy.getCoords() )
                self.assertTrue( mCpy.distanceToPoints( ptsCpy.getCoords() )[0].isUniform( abs( zLev ), 1e-12 ) )
        # And to finish : a test in 3D
        coords = mc.DataArrayDouble( [(-9.529956, -29.38396651, 8.2), (-9.529956, -30.37755672, 8.2), (-9.329956, -30.24157895, 8.2), (-9.329956, -29.26605263, 8.2), (-9.329956, -30.125, 9.224), (-9.329956, -29.18833333, 9.224)] )
        m3D = mc.MEDCouplingUMesh("",2)
        m3D.setCoords( coords )
        m3D.allocateCells()
        m3D.insertNextCell( mc.NORM_QUAD4, [0,1,2,3] )
        m3D.insertNextCell( mc.NORM_QUAD4, [4,1,0,5] )
        pt = mc.DataArrayDouble( [-9.35495563, -30.0, 8.2] , 1, 3 )
        a,b = m3D.distanceToPoints(pt)
        self.assertTrue( a.isEqual( mc.DataArrayDouble([0.0]), 1e-12 ) )
        self.assertTrue( b.isEqual( mc.DataArrayInt([0]) ) )
        # fmt: on


if __name__ == "__main__":
    unittest.main()
