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


import sys
from medcoupling import *
import unittest
from math import pi, e, sqrt, cos, sin
from datetime import datetime
from MEDCouplingDataForTest import MEDCouplingDataForTest
import rlcompleter, readline  # this line has to be here,to ensure a usability of MEDCoupling/MEDLoader. B4 removing it please notify to anthony.geay@cea.fr


class MEDCouplingBasicsTest7(unittest.TestCase):
    def testDAIGetIdsEqual1(self):
        tab1 = [5, -2, -4, -2, 3, 2, -2]
        da = DataArrayInt64.New()
        da.setValues(tab1, 7, 1)
        da2 = da.findIdsEqual(-2)
        self.assertEqual(3, da2.getNumberOfTuples())
        self.assertEqual(1, da2.getNumberOfComponents())
        expected1 = [1, 3, 6]
        self.assertEqual(expected1, da2.getValues())
        pass

    def testDAIGetIdsEqualList1(self):
        tab1 = [5, -2, -4, -2, 3, 2, -2]
        da = DataArrayInt64.New()
        da.setValues(tab1, 7, 1)
        da2 = da.findIdsEqualList([3, -2, 0])
        self.assertEqual(4, da2.getNumberOfTuples())
        self.assertEqual(1, da2.getNumberOfComponents())
        expected1 = [1, 3, 4, 6]
        self.assertEqual(expected1, da2.getValues())
        pass

    def testDAIsUniform1(self):
        tab1 = [1, 1, 1, 1, 1]
        da = DataArrayInt64.New()
        da.setValues(tab1, 5, 1)
        self.assertTrue(da.isUniform(1))
        da.setIJ(2, 0, 2)
        self.assertTrue(not da.isUniform(1))
        da.setIJ(2, 0, 1)
        self.assertTrue(da.isUniform(1))
        da2 = da.convertToDblArr()
        self.assertTrue(da2.isUniform(1.0, 1.0e-12))
        da2.setIJ(1, 0, 1.0 + 1.0e-13)
        self.assertTrue(da2.isUniform(1.0, 1.0e-12))
        da2.setIJ(1, 0, 1.0 + 1.0e-11)
        self.assertTrue(not da2.isUniform(1.0, 1.0e-12))
        pass

    def testDAIBuildComplement1(self):
        a = DataArrayInt64.New()
        tab = [3, 1, 7, 8]
        a.setValues(tab, 4, 1)
        b = a.buildComplement(12)
        self.assertEqual(8, b.getNumberOfTuples())
        self.assertEqual(1, b.getNumberOfComponents())
        expected1 = [0, 2, 4, 5, 6, 9, 10, 11]
        for i in range(8):
            self.assertEqual(expected1[i], b.getIJ(0, i))
            pass
        pass

    def testDAIBuildUnion1(self):
        a = DataArrayInt64.New()
        tab1 = [3, 1, 7, 8]
        a.setValues(tab1, 4, 1)
        c = DataArrayInt64.New()
        tab2 = [5, 3, 0, 18, 8]
        c.setValues(tab2, 5, 1)
        b = a.buildUnion(c)
        self.assertEqual(7, b.getNumberOfTuples())
        self.assertEqual(1, b.getNumberOfComponents())
        expected1 = [0, 1, 3, 5, 7, 8, 18]
        for i in range(7):
            self.assertEqual(expected1[i], b.getIJ(0, i))
            pass
        b = DataArrayInt64.BuildUnion([a, c])
        self.assertEqual(7, b.getNumberOfTuples())
        self.assertEqual(1, b.getNumberOfComponents())
        expected1 = [0, 1, 3, 5, 7, 8, 18]
        for i in range(7):
            self.assertEqual(expected1[i], b.getIJ(0, i))
            pass
        pass

    def testDAIBuildIntersection1(self):
        a = DataArrayInt64.New()
        tab1 = [3, 1, 7, 8]
        a.setValues(tab1, 4, 1)
        c = DataArrayInt64.New()
        tab2 = [5, 3, 0, 18, 8]
        c.setValues(tab2, 5, 1)
        b = a.buildIntersection(c)
        self.assertEqual(2, b.getNumberOfTuples())
        self.assertEqual(1, b.getNumberOfComponents())
        expected1 = [3, 8]
        for i in range(2):
            self.assertEqual(expected1[i], b.getIJ(0, i))
            pass
        b = DataArrayInt64.BuildIntersection([a, c])
        self.assertEqual(2, b.getNumberOfTuples())
        self.assertEqual(1, b.getNumberOfComponents())
        expected1 = [3, 8]
        for i in range(2):
            self.assertEqual(expected1[i], b.getIJ(0, i))
            pass
        pass

    def testDAIDeltaShiftIndex1(self):
        a = DataArrayInt64.New()
        tab = [1, 3, 6, 7, 7, 9, 15]
        a.setValues(tab, 7, 1)
        b = a.deltaShiftIndex()
        self.assertEqual(6, b.getNumberOfTuples())
        self.assertEqual(1, b.getNumberOfComponents())
        expected1 = [2, 3, 1, 0, 2, 6]
        for i in range(6):
            self.assertEqual(expected1[i], b.getIJ(0, i))
            pass
        pass

    def testDAIBuildSubstraction1(self):
        a = DataArrayInt64.New()
        aa = [2, 3, 6, 8, 9]
        a.setValues(aa, 5, 1)
        b = DataArrayInt64.New()
        bb = [1, 3, 5, 9, 11]
        b.setValues(bb, 5, 1)
        self.assertEqual([2, 6, 8], a.buildSubstraction(b).getValues())
        pass

    def testDAIBuildPermutationArr1(self):
        a = DataArrayInt64.New()
        a.setValues([4, 5, 6, 7, 8], 5, 1)
        b = DataArrayInt64.New()
        b.setValues([5, 4, 8, 6, 7], 5, 1)
        c = a.buildPermutationArr(b)
        self.assertEqual([1, 0, 4, 2, 3], c.getValues())
        self.assertTrue(a.isEqualWithoutConsideringStrAndOrder(b))
        b.setIJ(0, 0, 9)
        self.assertTrue(not a.isEqualWithoutConsideringStrAndOrder(b))
        self.assertRaises(InterpKernelException, a.buildPermutationArr, b)
        a.setIJ(3, 0, 4)
        b.setIJ(0, 0, 5)
        b.setIJ(4, 0, 4)  # a==[4,5,6,4,8] and b==[5,4,8,6,4]
        self.assertTrue(a.isEqualWithoutConsideringStrAndOrder(b))
        c = a.buildPermutationArr(b)
        self.assertEqual([1, 3, 4, 2, 3], c.getValues())
        d = b.convertToDblArr()
        expect3 = [4, 4, 5, 6, 8]
        b.sort()
        self.assertEqual(expect3, b.getValues())
        d.sort()
        self.assertEqual(5, d.getNumberOfTuples())
        self.assertEqual(1, d.getNumberOfComponents())
        for i in range(5):
            self.assertAlmostEqual(float(expect3[i]), d.getIJ(i, 0), 14)
            pass
        pass

    def testDAIAggregateMulti1(self):
        a = DataArrayInt64.New()
        a.setValues(list(range(4)), 2, 2)
        a.setName("aa")
        b = DataArrayInt64.New()
        b.setValues(list(range(6)), 3, 2)
        c = DataArrayInt64.Aggregate([a, b])
        self.assertEqual(list(range(4)) + list(range(6)), c.getValues())
        self.assertEqual("aa", c.getName())
        self.assertEqual(5, c.getNumberOfTuples())
        self.assertEqual(2, c.getNumberOfComponents())
        pass

    def testDAICheckAndPreparePermutation1(self):
        vals1 = [9, 10, 0, 6, 4, 11, 3, 7]
        expect1 = [5, 6, 0, 3, 2, 7, 1, 4]
        vals2 = [9, 10, 0, 6, 10, 11, 3, 7]
        da = DataArrayInt64.New()
        da.setValues(vals1, 8, 1)
        da2 = da.checkAndPreparePermutation()
        self.assertEqual(8, da2.getNumberOfTuples())
        self.assertEqual(1, da2.getNumberOfComponents())
        for i in range(8):
            self.assertEqual(expect1[i], da2.getIJ(i, 0))
            pass
        #
        da = DataArrayInt64.New()
        da.alloc(8, 1)
        da.iota(0)
        da2 = da.checkAndPreparePermutation()
        self.assertEqual(1, da2.getNumberOfComponents())
        self.assertTrue(da2.isIota(8))
        #
        da = DataArrayInt64.New()
        da.alloc(8, 1)
        da.setValues(vals2, 8, 1)
        self.assertRaises(InterpKernelException, da.checkAndPreparePermutation)
        pass

    def testDAIChangeSurjectiveFormat1(self):
        vals1 = [0, 3, 2, 3, 2, 2, 1, 2]
        expected1 = [0, 1, 2, 6, 8]
        expected2 = [0, 6, 2, 4, 5, 7, 1, 3]
        da = DataArrayInt64.New()
        da.setValues(vals1, 8, 1)
        #
        da2, da2I = da.changeSurjectiveFormat(4)
        self.assertEqual(5, da2I.getNumberOfTuples())
        self.assertEqual(8, da2.getNumberOfTuples())
        self.assertEqual(expected1, da2I.getValues())
        self.assertEqual(expected2, da2.getValues())
        #
        self.assertRaises(InterpKernelException, da.changeSurjectiveFormat, 3)
        #
        pass

    def testDAIGetIdsNotEqual1(self):
        d = DataArrayInt64.New()
        vals1 = [2, 3, 5, 6, 8, 5, 5, 6, 1, -5]
        d.setValues(vals1, 10, 1)
        d2 = d.findIdsNotEqual(5)
        self.assertEqual(7, d2.getNumberOfTuples())
        self.assertEqual(1, d2.getNumberOfComponents())
        expected1 = [0, 1, 3, 4, 7, 8, 9]
        for i in range(7):
            self.assertEqual(expected1[i], d2.getIJ(0, i))
            pass
        d.rearrange(2)
        self.assertRaises(InterpKernelException, d.findIdsNotEqual, 5)
        vals2 = [-4, 5, 6]
        vals3 = vals2
        d.rearrange(1)
        d3 = d.findIdsNotEqualList(vals3)
        self.assertEqual(5, d3.getNumberOfTuples())
        self.assertEqual(1, d3.getNumberOfComponents())
        expected2 = [0, 1, 4, 8, 9]
        for i in range(5):
            self.assertEqual(expected2[i], d3.getIJ(0, i))
            pass
        pass

    def testDAIComputeOffsets1(self):
        d = DataArrayInt64.New()
        vals1 = [3, 5, 1, 2, 0, 8]
        expected1 = [0, 3, 8, 9, 11, 11]
        d.setValues(vals1, 6, 1)
        d.computeOffsets()
        self.assertEqual(6, d.getNumberOfTuples())
        self.assertEqual(1, d.getNumberOfComponents())
        for i in range(6):
            self.assertEqual(expected1[i], d.getIJ(0, i))
            pass
        pass

    def testDAITransformWithIndArr1(self):
        if not MEDCouplingUse64BitIDs():
            return
        tab1 = [17, 18, 22, 19]
        tab2 = [0, 1, 1, 3, 3, 0, 1, 3, 2, 2, 3, 0]
        expected = [17, 18, 18, 19, 19, 17, 18, 19, 22, 22, 19, 17]
        d = DataArrayInt64.New()
        d.setValues(tab1, 4, 1)
        d1 = DataArrayInt64.New()
        d1.setValues(tab2, 12, 1)
        d2 = d1[:]
        #
        d1.transformWithIndArr(d)
        self.assertEqual(12, d1.getNumberOfTuples())
        self.assertEqual(1, d1.getNumberOfComponents())
        for i in range(12):
            self.assertEqual(expected[i], d1.getIJ(i, 0))
            pass
        #
        d1 = d2
        d1.transformWithIndArr(tab1)
        self.assertEqual(12, d1.getNumberOfTuples())
        self.assertEqual(1, d1.getNumberOfComponents())
        for i in range(12):
            self.assertEqual(expected[i], d1.getIJ(i, 0))
            pass
        pass

    def testDAIBuildPermArrPerLevel1(self):
        arr = [2, 0, 1, 1, 0, 1, 2, 0, 1, 1, 0, 0]
        expected1 = [10, 0, 5, 6, 1, 7, 11, 2, 8, 9, 3, 4]
        da = DataArrayInt64.New()
        da.setValues(arr, 12, 1)
        da2 = da.buildPermArrPerLevel()
        self.assertEqual(12, da2.getNumberOfTuples())
        self.assertEqual(1, da2.getNumberOfComponents())
        for i in range(12):
            self.assertEqual(expected1[i], da2.getIJ(i, 0))
            pass
        pass

    def testDAIOperations1(self):
        arr1 = [-1, -2, 4, 7, 3, 2, 6, 6, 4, 3, 0, 1]
        da = DataArrayInt64.New()
        da.setValues(arr1, 4, 3)
        da1 = DataArrayInt64.New()
        da1.alloc(12, 1)
        da1.iota(2)
        self.assertRaises(
            InterpKernelException, DataArrayInt64.Add, da, da1
        )  # not same number of tuples/Components
        da1.rearrange(3)
        da2 = DataArrayInt64.Add(da, da1)
        self.assertEqual(4, da2.getNumberOfTuples())
        self.assertEqual(3, da2.getNumberOfComponents())
        expected1 = [1, 1, 8, 12, 9, 9, 14, 15, 14, 14, 12, 14]
        for i in range(12):
            self.assertEqual(expected1[i], da2.getIJ(0, i))
            pass
        da1.substractEqual(da)
        expected2 = [3, 5, 0, -2, 3, 5, 2, 3, 6, 8, 12, 12]
        for i in range(12):
            self.assertEqual(expected2[i], da1.getIJ(0, i))
            pass
        da1.rearrange(1)
        da1.iota(2)
        da1.rearrange(3)
        da1.addEqual(da)
        for i in range(12):
            self.assertEqual(expected1[i], da1.getIJ(0, i))
            pass
        da1.rearrange(1)
        da1.iota(2)
        da1.rearrange(3)
        da2 = DataArrayInt64.Multiply(da, da1)
        self.assertEqual(4, da2.getNumberOfTuples())
        self.assertEqual(3, da2.getNumberOfComponents())
        expected3 = [-2, -6, 16, 35, 18, 14, 48, 54, 40, 33, 0, 13]
        for i in range(12):
            self.assertEqual(expected3[i], da2.getIJ(0, i))
            pass
        da.divideEqual(da1)
        self.assertEqual(4, da.getNumberOfTuples())
        self.assertEqual(3, da.getNumberOfComponents())
        expected4 = [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        for i in range(12):
            self.assertEqual(expected4[i], da.getIJ(0, i))
            pass
        da.setValues(arr1, 4, 3)
        da1.multiplyEqual(da)
        self.assertEqual(4, da1.getNumberOfTuples())
        self.assertEqual(3, da1.getNumberOfComponents())
        for i in range(12):
            self.assertEqual(expected3[i], da1.getIJ(0, i))
            pass
        da1.rearrange(1)
        da1.iota(2)
        da1.rearrange(3)
        da2 = DataArrayInt64.Divide(da, da1)
        self.assertEqual(4, da2.getNumberOfTuples())
        self.assertEqual(3, da2.getNumberOfComponents())
        for i in range(12):
            self.assertEqual(expected4[i], da2.getIJ(0, i))
            pass
        da1.applyInv(321)
        self.assertEqual(4, da1.getNumberOfTuples())
        self.assertEqual(3, da1.getNumberOfComponents())
        expected5 = [160, 107, 80, 64, 53, 45, 40, 35, 32, 29, 26, 24]
        for i in range(12):
            self.assertEqual(expected5[i], da1.getIJ(0, i))
            pass
        da1.applyDivideBy(2)
        self.assertEqual(4, da1.getNumberOfTuples())
        self.assertEqual(3, da1.getNumberOfComponents())
        expected6 = [80, 53, 40, 32, 26, 22, 20, 17, 16, 14, 13, 12]
        for i in range(12):
            self.assertEqual(expected6[i], da1.getIJ(0, i))
            pass
        expected7 = [3, 4, 5, 4, 5, 1, 6, 3, 2, 0, 6, 5]
        da1.applyModulus(7)
        for i in range(12):
            self.assertEqual(expected7[i], da1.getIJ(0, i))
            pass
        da1.applyLin(1, 1)
        expected8 = [3, 3, 3, 3, 3, 1, 3, 3, 0, 0, 3, 3]
        da1.applyRModulus(3)
        for i in range(12):
            self.assertEqual(expected8[i], da1.getIJ(0, i))
            pass
        pass

    def testDAITransformWithIndArrR1(self):
        tab1 = [2, 4, 5, 3, 6, 7]
        tab2 = [-1, -1, 0, 1, 2, 3, 4, 5, -1, -1, -1, -1]
        expected = [0, 3, 1, 2, 4, 5]
        d = DataArrayInt64.New()
        d.setValues(tab1, 6, 1)
        d1 = DataArrayInt64.New()
        d1.setValues(tab2, 12, 1)
        d2 = d1[:]
        #
        d3 = d.transformWithIndArrR(d1)
        self.assertEqual(6, d3.getNumberOfTuples())
        self.assertEqual(1, d3.getNumberOfComponents())
        for i in range(6):
            self.assertEqual(expected[i], d3.getIJ(i, 0))
            pass
        #
        d1 = d2
        d3 = d.transformWithIndArrR(tab2)
        self.assertEqual(6, d3.getNumberOfTuples())
        self.assertEqual(1, d3.getNumberOfComponents())
        for i in range(6):
            self.assertEqual(expected[i], d3.getIJ(i, 0))
            pass
        pass

    def testDAISplitByValueRange1(self):
        val1 = [6, 5, 0, 3, 2, 7, 8, 1, 4]
        val2 = [0, 4, 9]
        d = DataArrayInt64.New()
        d.setValues(val1, 9, 1)
        e, f, g = d.splitByValueRange(val2)
        self.assertEqual(9, e.getNumberOfTuples())
        self.assertEqual(1, e.getNumberOfComponents())
        self.assertEqual(9, f.getNumberOfTuples())
        self.assertEqual(1, f.getNumberOfComponents())
        self.assertEqual(2, g.getNumberOfTuples())
        self.assertEqual(1, g.getNumberOfComponents())
        #
        expected1 = [1, 1, 0, 0, 0, 1, 1, 0, 1]
        expected2 = [2, 1, 0, 3, 2, 3, 4, 1, 0]
        for i in range(9):
            self.assertEqual(expected1[i], e.getIJ(i, 0))
            self.assertEqual(expected2[i], f.getIJ(i, 0))
            pass
        self.assertEqual(0, g.getIJ(0, 0))
        self.assertEqual(1, g.getIJ(1, 0))
        #
        d.setIJ(6, 0, 9)
        self.assertRaises(InterpKernelException, d.splitByValueRange, val2)
        # non regression test in python wrapping
        rg = DataArrayInt64(
            [
                0,
                10,
                29,
                56,
                75,
                102,
                121,
                148,
                167,
                194,
                213,
                240,
                259,
                286,
                305,
                332,
                351,
                378,
                397,
                424,
                443,
                470,
                489,
                516,
            ]
        )
        a, b, c = DataArrayInt64([75]).splitByValueRange(rg)
        self.assertTrue(a.isEqual(DataArrayInt64([4])))
        self.assertTrue(b.isEqual(DataArrayInt64([0])))
        self.assertTrue(c.isEqual(DataArrayInt64([4])))
        pass

    def testDAIBuildExplicitArrByRanges1(self):
        d = DataArrayInt64.New()
        vals1 = [0, 2, 3]
        d.setValues(vals1, 3, 1)
        e = DataArrayInt64.New()
        vals2 = [0, 3, 6, 10, 14, 20]
        e.setValues(vals2, 6, 1)
        #
        f = d.buildExplicitArrByRanges(e)
        self.assertEqual(11, f.getNumberOfTuples())
        self.assertEqual(1, f.getNumberOfComponents())
        expected1 = [0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 13]
        for i in range(11):
            self.assertEqual(expected1[i], f.getIJ(i, 0))
            pass
        pass

    def testDAIComputeOffsets2(self):
        d = DataArrayInt64.New()
        vals1 = [3, 5, 1, 2, 0, 8]
        expected1 = [0, 3, 8, 9, 11, 11, 19]
        d.setValues(vals1, 6, 1)
        d.computeOffsetsFull()
        self.assertEqual(7, d.getNumberOfTuples())
        self.assertEqual(1, d.getNumberOfComponents())
        for i in range(7):
            self.assertEqual(expected1[i], d.getIJ(0, i))
            pass
        pass

    def testDAIBuildOld2NewArrayFromSurjectiveFormat2(self):
        arr = [0, 3, 5, 7, 9]
        arrI = [0, 2, 5]
        a = DataArrayInt.New()
        a.setValues(arr, 5, 1)
        b = DataArrayInt.New()
        b.setValues(arrI, 3, 1)
        ret, newNbTuple = DataArrayInt64.ConvertIndexArrayToO2N(10, a, b)
        expected = [0, 1, 2, 0, 3, 4, 5, 4, 6, 4]
        self.assertEqual(10, ret.getNbOfElems())
        self.assertEqual(7, newNbTuple)
        self.assertEqual(1, ret.getNumberOfComponents())
        self.assertEqual(expected, ret.getValues())
        self.assertRaises(
            InterpKernelException, DataArrayInt64.ConvertIndexArrayToO2N, 9, a, b
        )
        pass

    def testDAIBuildUnique1(self):
        d = DataArrayInt64([1, 2, 2, 3, 3, 3, 3, 4, 5, 5, 7, 7, 7, 19])
        e = d.buildUnique()
        self.assertTrue(e.isEqual(DataArrayInt64([1, 2, 3, 4, 5, 7, 19])))
        pass

    def testDAIPartitionByDifferentValues1(self):
        d = DataArrayInt64([1, 0, 1, 2, 0, 2, 2, -3, 2])
        expected = [[-3, [7]], [0, [1, 4]], [1, [0, 2]], [2, [3, 5, 6, 8]]]
        for i, elt in enumerate(zip(*d.partitionByDifferentValues())):
            self.assertEqual(expected[i][0], elt[1])
            self.assertEqual(expected[i][1], elt[0].getValues())
            pass
        pass

    def testDAICheckMonotonic1(self):
        data1 = [-1, 0, 2, 2, 4, 5]
        data2 = [6, 2, 0, -8, -9, -56]
        data3 = [-1, 0, 3, 2, 4, 6]
        data4 = [7, 5, 2, 3, 0, -6]
        d = DataArrayInt64.New(data1)
        self.assertTrue(d.isMonotonic(True))
        self.assertTrue(not d.isMonotonic(False))
        d.checkMonotonic(True)
        self.assertRaises(InterpKernelException, d.checkMonotonic, False)
        d = DataArrayInt64.New(data2)
        self.assertTrue(d.isMonotonic(False))
        self.assertTrue(not d.isMonotonic(True))
        d.checkMonotonic(False)
        self.assertRaises(InterpKernelException, d.checkMonotonic, True)
        d = DataArrayInt64.New(data3)
        self.assertTrue(not d.isMonotonic(False))
        self.assertTrue(not d.isMonotonic(True))
        self.assertRaises(InterpKernelException, d.checkMonotonic, True)
        self.assertRaises(InterpKernelException, d.checkMonotonic, False)
        d = DataArrayInt64.New(data4)
        self.assertTrue(not d.isMonotonic(False))
        self.assertTrue(not d.isMonotonic(True))
        self.assertRaises(InterpKernelException, d.checkMonotonic, True)
        self.assertRaises(InterpKernelException, d.checkMonotonic, False)
        d = DataArrayInt64.New(0, 1)
        self.assertTrue(d.isMonotonic(True))
        self.assertTrue(d.isMonotonic(False))
        d.checkMonotonic(True)
        d.checkMonotonic(False)
        d = DataArrayInt64.New(data4, 3, 2)  # throw because nbComp!=1
        self.assertRaises(InterpKernelException, d.isMonotonic, True)
        self.assertRaises(InterpKernelException, d.isMonotonic, False)
        self.assertRaises(InterpKernelException, d.checkMonotonic, True)
        self.assertRaises(InterpKernelException, d.checkMonotonic, False)
        pass

    def testDAIBuildSubstractionOptimized1(self):
        da1 = DataArrayInt64.New([1, 3, 5, 6, 7, 9, 13])
        da2 = DataArrayInt64.New([3, 5, 9])
        da3 = DataArrayInt64.New([1, 3, 5])
        da4 = DataArrayInt64.New([1, 3, 5, 6, 7, 9, 13])
        #
        a = da1.buildSubstractionOptimized(da2)
        self.assertTrue(a.isEqual(DataArrayInt64([1, 6, 7, 13])))
        #
        a = da1.buildSubstractionOptimized(da3)
        self.assertTrue(a.isEqual(DataArrayInt64([6, 7, 9, 13])))
        #
        a = da1.buildSubstractionOptimized(da4)
        self.assertTrue(a.isEqual(DataArrayInt64([])))
        pass

    def testDAIIsStrictlyMonotonic1(self):
        da1 = DataArrayInt64.New([1, 3, 5, 6, 7, 9, 13])
        self.assertTrue(da1.isStrictlyMonotonic(True))
        da1.checkStrictlyMonotonic(True)
        self.assertTrue(da1.isMonotonic(True))
        da1.checkMonotonic(True)
        self.assertTrue(not da1.isStrictlyMonotonic(False))
        self.assertRaises(InterpKernelException, da1.checkStrictlyMonotonic, False)
        self.assertTrue(not da1.isMonotonic(False))
        self.assertRaises(InterpKernelException, da1.checkMonotonic, False)
        #
        da1 = DataArrayInt64.New([1, 3, 5, 6, 6, 9, 13])
        self.assertTrue(not da1.isStrictlyMonotonic(True))
        self.assertRaises(InterpKernelException, da1.checkStrictlyMonotonic, True)
        self.assertTrue(da1.isMonotonic(True))
        da1.checkMonotonic(True)
        self.assertTrue(not da1.isStrictlyMonotonic(False))
        self.assertRaises(InterpKernelException, da1.checkStrictlyMonotonic, False)
        self.assertTrue(not da1.isMonotonic(False))
        self.assertRaises(InterpKernelException, da1.checkMonotonic, False)
        #
        da1 = DataArrayInt64.New([1, 3, 5, 6, 5, 9, 13])
        self.assertTrue(not da1.isStrictlyMonotonic(True))
        self.assertRaises(InterpKernelException, da1.checkStrictlyMonotonic, True)
        self.assertTrue(not da1.isMonotonic(True))
        self.assertRaises(InterpKernelException, da1.checkMonotonic, True)
        self.assertTrue(not da1.isStrictlyMonotonic(False))
        self.assertRaises(InterpKernelException, da1.checkStrictlyMonotonic, False)
        self.assertTrue(not da1.isMonotonic(False))
        self.assertRaises(InterpKernelException, da1.checkMonotonic, False)
        #
        da1 = DataArrayInt64.New([13, 9, 7, 6, 5, 3, 1])
        self.assertTrue(not da1.isStrictlyMonotonic(True))
        self.assertRaises(InterpKernelException, da1.checkStrictlyMonotonic, True)
        self.assertTrue(not da1.isMonotonic(True))
        self.assertRaises(InterpKernelException, da1.checkMonotonic, True)
        self.assertTrue(da1.isStrictlyMonotonic(False))
        da1.checkStrictlyMonotonic(False)
        self.assertTrue(da1.isMonotonic(False))
        da1.checkMonotonic(False)
        #
        da1 = DataArrayInt64.New([13, 9, 6, 6, 5, 3, 1])
        self.assertTrue(not da1.isStrictlyMonotonic(True))
        self.assertRaises(InterpKernelException, da1.checkStrictlyMonotonic, True)
        self.assertTrue(not da1.isMonotonic(True))
        self.assertRaises(InterpKernelException, da1.checkMonotonic, True)
        self.assertTrue(not da1.isStrictlyMonotonic(False))
        self.assertRaises(InterpKernelException, da1.checkStrictlyMonotonic, False)
        self.assertTrue(da1.isMonotonic(False))
        da1.checkMonotonic(False)
        #
        da1 = DataArrayInt64.New([13, 9, 5, 6, 5, 3, 1])
        self.assertTrue(not da1.isStrictlyMonotonic(True))
        self.assertRaises(InterpKernelException, da1.checkStrictlyMonotonic, True)
        self.assertTrue(not da1.isMonotonic(True))
        self.assertRaises(InterpKernelException, da1.checkMonotonic, True)
        self.assertTrue(not da1.isStrictlyMonotonic(False))
        self.assertRaises(InterpKernelException, da1.checkStrictlyMonotonic, False)
        self.assertTrue(not da1.isMonotonic(False))
        self.assertRaises(InterpKernelException, da1.checkMonotonic, False)
        #
        da1 = DataArrayInt64.New([])
        self.assertTrue(da1.isStrictlyMonotonic(True))
        da1.checkStrictlyMonotonic(True)
        self.assertTrue(da1.isMonotonic(True))
        da1.checkMonotonic(True)
        self.assertTrue(da1.isStrictlyMonotonic(False))
        da1.checkStrictlyMonotonic(False)
        self.assertTrue(da1.isMonotonic(False))
        da1.checkMonotonic(False)
        #
        da1 = DataArrayInt64.New([13])
        self.assertTrue(da1.isStrictlyMonotonic(True))
        da1.checkStrictlyMonotonic(True)
        self.assertTrue(da1.isMonotonic(True))
        da1.checkMonotonic(True)
        self.assertTrue(da1.isStrictlyMonotonic(False))
        da1.checkStrictlyMonotonic(False)
        self.assertTrue(da1.isMonotonic(False))
        da1.checkMonotonic(False)
        pass

    def testDAIIndicesOfSubPart(self):
        a = DataArrayInt64([9, 10, 0, 6, 4, 11, 3, 8])
        b = DataArrayInt64([6, 0, 11, 8])
        c = a.indicesOfSubPart(b)
        self.assertTrue(c.isEqual(DataArrayInt([3, 2, 5, 7])))
        #
        d = DataArrayInt64([9, 10, 0, 6, 4, 11, 0, 8])
        self.assertRaises(
            InterpKernelException, d.indicesOfSubPart, b
        )  # 0 appears twice in the d array
        f = DataArrayInt64([6, 0, 11, 8, 12])
        self.assertRaises(
            InterpKernelException, a.indicesOfSubPart, f
        )  # 12 in f does not exist in a
        pass

    def testDAIsortToHaveConsecutivePairs(self):
        dref = DataArrayInt64(
            [(6, 216), (216, 218), (218, 220), (220, 222), (222, 224), (224, 226)]
        )
        dtest = DataArrayInt64(
            [(6, 216), (218, 216), (224, 226), (222, 220), (218, 220), (222, 224)]
        )
        dtest.sortToHaveConsecutivePairs()
        self.assertTrue(dtest.isEqual(dref))

    def testDAIFromLinkedListOfPairToList1(self):
        d = DataArrayInt64([(5, 7), (7, 3), (3, 12), (12, 17)])
        zeRes = DataArrayInt64([5, 7, 3, 12, 17])
        self.assertTrue(d.fromLinkedListOfPairToList().isEqual(zeRes))
        d.rearrange(1)
        self.assertRaises(InterpKernelException, d.fromLinkedListOfPairToList)
        d.rearrange(2)
        self.assertTrue(d.fromLinkedListOfPairToList().isEqual(zeRes))
        d2 = DataArrayInt64([(5, 7)])
        self.assertTrue(d2.fromLinkedListOfPairToList().isEqual(DataArrayInt64([5, 7])))
        d3 = DataArrayInt64([(5, 7), (7, 3), (4, 12), (12, 17)])
        self.assertRaises(
            InterpKernelException, d3.fromLinkedListOfPairToList
        )  # not a linked list of pair
        d4 = DataArrayInt64([(5, 7), (7, 3), (12, 3), (12, 17)])
        self.assertRaises(
            InterpKernelException, d4.fromLinkedListOfPairToList
        )  # not a linked list of pair, but can be repaired !
        d4.sortEachPairToMakeALinkedList()
        self.assertTrue(d4.fromLinkedListOfPairToList().isEqual(zeRes))
        pass

    def testDAIfindIdsExt1(self):
        d = DataArrayInt64([4, 6, -2, 3, 7, 0, 10])
        self.assertTrue(
            d.findIdsGreaterOrEqualTo(3).isEqual(DataArrayInt([0, 1, 3, 4, 6]))
        )
        self.assertTrue(d.findIdsGreaterThan(3).isEqual(DataArrayInt([0, 1, 4, 6])))
        self.assertTrue(d.findIdsLowerThan(3).isEqual(DataArrayInt([2, 5])))
        self.assertTrue(d.findIdsLowerOrEqualTo(3).isEqual(DataArrayInt([2, 3, 5])))
        pass

    def testDAICheckUniformAndGuess1(self):
        d = DataArrayInt64([3, 3], 1, 2)
        self.assertRaises(
            InterpKernelException, d.checkUniformAndGuess
        )  # non single compo
        d = DataArrayInt64([])
        self.assertRaises(InterpKernelException, d.checkUniformAndGuess)  # empty
        d = DataArrayInt64()
        self.assertRaises(
            InterpKernelException, d.checkUniformAndGuess
        )  # non allocated
        d = DataArrayInt64([3, 3, 3])
        self.assertEqual(3, d.checkUniformAndGuess())
        d = DataArrayInt64([7])
        self.assertEqual(7, d.checkUniformAndGuess())
        d = DataArrayInt64([3, 4, 3])
        self.assertRaises(InterpKernelException, d.checkUniformAndGuess)  # non uniform
        pass

    def testDAIFindIdForEach1(self):
        a1 = DataArrayInt64([17, 27, 2, 10, -4, 3, 12, 27, 16])
        b1 = DataArrayInt64([3, 16, -4, 27, 17])
        ret = a1.findIdForEach(b1)
        self.assertTrue(ret.isEqual(DataArrayInt([5, 8, 4, 7, 0])))
        self.assertTrue(a1[ret].isEqual(b1))
        b2 = DataArrayInt64([3, 16, 22, 27, 17])
        self.assertRaises(InterpKernelException, a1.findIdForEach, b2)  # 22 not in a1 !
        a1.rearrange(3)
        self.assertRaises(
            InterpKernelException, a1.findIdForEach, b1
        )  # a1 is not single component
        pass

    def testGlobalHelpers(self):
        arr0 = vtk2med_cell_types()
        self.assertEqual(len(arr0), 43)
        arr1 = med2vtk_cell_types()
        self.assertEqual(len(arr1), 34)
        arr2 = AllGeometricTypes()
        self.assertEqual(len(arr2), 25)
        for elt in arr2:
            MEDCouplingUMesh.GetReprOfGeometricType(elt)
            self.assertNotEqual(MEDCouplingUMesh.GetDimensionOfGeometricType(elt), -1)
        pass

    def testVoronoi2D_3(self):
        """
        Non regression test for EDF20418 : After 8.5.0 MEDCouplingUMesh.Interset2DMeshes method (called by voronoize) is sensible to cell orientation of 2 input meshes. This test check correct behavior in
        case of non standard orientation input cell.
        """
        coo = DataArrayDouble(
            [
                0.018036113896685007,
                0.030867224641316506,
                0.019000000000000003,
                0.030833333333333407,
                0.018518056948342503,
                0.030850278987324904,
                0.018773068345659904,
                0.031180320157635305,
                0.018546136691319805,
                0.031527306981937307,
                0.018291125294002404,
                0.031197265811626906,
            ],
            6,
            2,
        )
        m = MEDCouplingUMesh("mesh", 2)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_TRI3, [0, 1, 4])
        f = MEDCouplingFieldDouble(ON_GAUSS_PT)
        f.setMesh(m)
        f.setArray(
            DataArrayDouble(
                [
                    12613576.708019681,
                    18945164.734307285,
                    22385248.637775388,
                    17074219.938821714,
                    19361929.467164982,
                    19258841.562907547,
                ]
            )
        )
        f.setGaussLocalizationOnType(
            NORM_TRI3,
            [0, 0, 1, 0, 0, 1],
            [
                0.0915762,
                0.0915762,
                0.816848,
                0.0915762,
                0.0915762,
                0.816848,
                0.445948,
                0.108103,
                0.445948,
                0.445948,
                0.108103,
                0.445948,
            ],
            [0.0549759, 0.0549759, 0.0549759, 0.111691, 0.111691, 0.111691],
        )
        f.setName("field")
        f_voro = f.voronoize(1e-13)
        ref_area = DataArrayDouble(
            [
                4.6679303278867127,
                4.2514546761810138,
                4.2514546761809337,
                6.6206415950989804,
                6.2643538685231039,
                6.6206415950989884,
            ]
        )
        area = f_voro.buildMeasureField(True).getArray() * 1e8
        self.assertTrue(ref_area.isEqual(area, 1e-6))
        ref_bary = DataArrayDouble(
            [
                (0.018231625096313969, 0.030950287685553721),
                (0.018826045778781105, 0.030916927013719033),
                (0.018533397739746087, 0.031364396601025746),
                (0.018541498169815956, 0.030944333493252929),
                (0.018660195622447071, 0.031132366117047686),
                (0.018400646702087166, 0.031159700554391174),
            ]
        )
        bary = f_voro.getMesh().computeCellCenterOfMass()
        self.assertTrue(ref_bary.isEqual(bary, 1e-8))
        self.assertTrue(f_voro.getArray().isEqual(f.getArray(), 1e-8))
        pass

    def testDAIOccurenceRankInThis(self):
        arr = DataArrayInt([5, 3, 2, 1, 4, 5, 2, 1, 0, 11, 5, 4])
        self.assertTrue(
            arr.occurenceRankInThis().isEqual(
                DataArrayInt([0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 2, 1])
            )
        )

    def testDAIFindPermutationFromFirstToSecondDuplicate(self):
        arr0 = DataArrayInt([5, 3, 2, 1, 4, 5, 2, 1, 0, 11, 5, 4])
        arr1 = DataArrayInt([0, 1, 1, 2, 2, 3, 4, 4, 5, 5, 5, 11])
        self.assertTrue(
            DataArrayInt.FindPermutationFromFirstToSecondDuplicate(arr0, arr1).isEqual(
                DataArrayInt([8, 5, 3, 1, 6, 9, 4, 2, 0, 11, 10, 7])
            )
        )
        self.assertTrue(
            DataArrayInt.FindPermutationFromFirstToSecondDuplicate(arr1, arr0).isEqual(
                DataArrayInt([8, 3, 7, 2, 6, 1, 4, 11, 0, 5, 10, 9])
            )
        )

    def testDAIIndexOfSameConsecutiveValueGroups(self):
        arr = DataArrayInt([0, 1, 1, 2, 2, 3, 4, 4, 5, 5, 5, 11])
        self.assertTrue(
            arr.indexOfSameConsecutiveValueGroups().isEqual(
                DataArrayInt([0, 1, 3, 5, 6, 8, 11, 12])
            )
        )

    def testSkyLineGroupPacks(self):
        arr = DataArrayInt(
            [
                1,
                4,
                5,
                0,
                2,
                4,
                5,
                6,
                1,
                3,
                5,
                6,
                7,
                2,
                6,
                7,
                0,
                1,
                5,
                8,
                9,
                0,
                1,
                2,
                4,
                6,
                8,
                9,
                10,
                1,
                2,
                3,
                5,
                7,
                9,
                10,
                11,
                2,
                3,
                6,
                10,
                11,
                4,
                5,
                9,
                12,
                13,
                4,
                5,
                6,
                8,
                10,
                12,
                13,
                14,
                5,
                6,
                7,
                9,
                11,
                13,
                14,
                15,
                6,
                7,
                10,
                14,
                15,
                8,
                9,
                13,
                8,
                9,
                10,
                12,
                14,
                9,
                10,
                11,
                13,
                15,
                10,
                11,
                14,
            ]
        )
        arrI = DataArrayInt(
            [0, 3, 8, 13, 16, 21, 29, 37, 42, 47, 55, 63, 68, 71, 76, 81, 84]
        )
        sk = MEDCouplingSkyLineArray(arrI, arr)
        part = DataArrayInt([0, 3, 4, 7, 16])
        sk2 = sk.groupPacks(part)
        self.assertTrue(sk2.getValuesArray().isEqual(arr))
        self.assertTrue(sk2.getIndexArray().isEqual(DataArrayInt([0, 13, 16, 37, 84])))

    def testSkyLineUniqueNotSortedByPack(self):
        arrI = DataArrayInt([0, 3, 9, 15, 18, 24, 36, 48, 54])
        arr = DataArrayInt(
            [
                1,
                4,
                5,
                0,
                4,
                5,
                2,
                5,
                6,
                3,
                6,
                7,
                1,
                5,
                6,
                2,
                6,
                7,
                0,
                1,
                5,
                5,
                8,
                9,
                0,
                1,
                4,
                6,
                9,
                10,
                1,
                2,
                4,
                6,
                8,
                9,
                2,
                3,
                5,
                7,
                9,
                10,
                1,
                2,
                5,
                7,
                10,
                11,
                2,
                3,
                6,
                6,
                10,
                11,
            ]
        )
        sk = MEDCouplingSkyLineArray(arrI, arr)
        sk2 = sk.uniqueNotSortedByPack()
        self.assertTrue(
            sk2.getIndexArray().isEqual(DataArrayInt([0, 3, 8, 13, 16, 21, 29, 37, 42]))
        )
        self.assertTrue(
            sk2.getValuesArray().isEqual(
                DataArrayInt(
                    [
                        1,
                        4,
                        5,
                        0,
                        2,
                        4,
                        5,
                        6,
                        1,
                        3,
                        5,
                        6,
                        7,
                        2,
                        6,
                        7,
                        0,
                        1,
                        5,
                        8,
                        9,
                        0,
                        1,
                        2,
                        4,
                        6,
                        8,
                        9,
                        10,
                        1,
                        2,
                        3,
                        5,
                        7,
                        9,
                        10,
                        11,
                        2,
                        3,
                        6,
                        10,
                        11,
                    ]
                )
            )
        )

    def testSkyLineAggregatePacks1(self):
        arr = DataArrayDouble(3)
        arr.iota()
        m = MEDCouplingCMesh()
        m.setCoords(arr, arr)
        m = m.buildUnstructured()
        a, b = m.computeEnlargedNeighborsOfNodes()
        sk = MEDCouplingSkyLineArray(b, a)
        sk1 = sk.deepCopy()
        sk1.getValuesArray()[:] *= 2
        sk2 = sk.deepCopy()
        sk2.getValuesArray()[:] *= 3
        skOut = MEDCouplingSkyLineArray.AggregatePacks([sk, sk1, sk2])
        self.assertTrue(
            skOut.getIndexArray().isEqual(
                DataArrayInt([0, 9, 24, 33, 48, 72, 87, 96, 111, 120])
            )
        )
        self.assertTrue(
            skOut.getValuesArray().isEqual(
                DataArrayInt(
                    [
                        1,
                        3,
                        4,
                        2,
                        6,
                        8,
                        3,
                        9,
                        12,
                        0,
                        2,
                        3,
                        4,
                        5,
                        0,
                        4,
                        6,
                        8,
                        10,
                        0,
                        6,
                        9,
                        12,
                        15,
                        1,
                        4,
                        5,
                        2,
                        8,
                        10,
                        3,
                        12,
                        15,
                        0,
                        1,
                        4,
                        6,
                        7,
                        0,
                        2,
                        8,
                        12,
                        14,
                        0,
                        3,
                        12,
                        18,
                        21,
                        0,
                        1,
                        2,
                        3,
                        5,
                        6,
                        7,
                        8,
                        0,
                        2,
                        4,
                        6,
                        10,
                        12,
                        14,
                        16,
                        0,
                        3,
                        6,
                        9,
                        15,
                        18,
                        21,
                        24,
                        1,
                        2,
                        4,
                        7,
                        8,
                        2,
                        4,
                        8,
                        14,
                        16,
                        3,
                        6,
                        12,
                        21,
                        24,
                        3,
                        4,
                        7,
                        6,
                        8,
                        14,
                        9,
                        12,
                        21,
                        3,
                        4,
                        5,
                        6,
                        8,
                        6,
                        8,
                        10,
                        12,
                        16,
                        9,
                        12,
                        15,
                        18,
                        24,
                        4,
                        5,
                        7,
                        8,
                        10,
                        14,
                        12,
                        15,
                        21,
                    ]
                )
            )
        )

    def testDACopySorted1(self):
        d = DataArrayInt32([5, 1, 100, 20])
        self.assertTrue(d.copySorted().isEqual(DataArrayInt32([1, 5, 20, 100])))
        d = DataArrayInt64([5, 1, 100, 20])
        self.assertTrue(d.copySorted().isEqual(DataArrayInt64([1, 5, 20, 100])))
        d = DataArrayInt([5, 1, 100, 20])
        self.assertTrue(d.copySorted().isEqual(DataArrayInt([1, 5, 20, 100])))
        d = DataArrayDouble([5, 1, 100, 20])
        self.assertTrue(d.copySorted().isEqual(DataArrayDouble([1, 5, 20, 100]), 1e-10))

    def testFieldAreStrictlyCompatible(self):
        arr = DataArrayDouble(10)
        arr.iota()
        m = MEDCouplingCMesh()
        m.setCoords(arr, arr)
        m = m.buildUnstructured()
        f = MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        f2 = MEDCouplingFieldDouble(ON_CELLS)
        f2.setMesh(m)
        self.assertTrue(f.areStrictlyCompatible(f2))
        self.assertTrue(f.areStrictlyCompatibleForMulDiv(f2))
        f2.setMesh(f2.getMesh().deepCopy())
        self.assertTrue(not f.areStrictlyCompatible(f2))
        self.assertTrue(not f.areStrictlyCompatibleForMulDiv(f2))
        f3 = MEDCouplingFieldDouble(ON_NODES)
        f3.setMesh(m)
        self.assertTrue(not f.areStrictlyCompatible(f3))
        self.assertTrue(not f.areStrictlyCompatibleForMulDiv(f3))

    def testBugZipConnectivityTraducer(self):
        """
        Non regression test : here cell #1 and cell #2 are nearly the same but not the same. zipConnectivityTraducer called by areCellsIncludedIn
        failed to capture that.
        """
        coo = DataArrayDouble([0, 1, 2, 3], 4, 1)
        m = MEDCouplingUMesh("", 1)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_SEG2, [0, 1])
        m.insertNextCell(NORM_SEG2, [1, 2])
        m.insertNextCell(NORM_SEG2, [2, 1])
        #
        a, b = m.areCellsIncludedIn(m, 0)
        self.assertTrue(a)
        self.assertTrue(b.isIota(3))
        #
        self.assertTrue(m.deepCopy().zipConnectivityTraducer(0).isIota(3))
        self.assertTrue(m.deepCopy().zipConnectivityTraducer(1).isIota(3))
        self.assertTrue(
            m.deepCopy().zipConnectivityTraducer(2).isEqual(DataArrayInt([0, 1, 1]))
        )

    def testBugAreCellsIncludedIn1(self):
        """
        Non regression test: a.areCellsIncludedIn(b) was buggy when some cells in b were duplicated into a following specified policy.
        """
        coo = DataArrayDouble([0, 1, 2, 3, 4], 5, 1)
        m = MEDCouplingUMesh("", 1)
        m.setCoords(coo)
        m.allocateCells()
        # m contains several duplicated cells - some of those duplicated cells will be in m1
        for i in range(3):
            m.insertNextCell(NORM_SEG2, [0, 1])
        for i in range(4):
            m.insertNextCell(NORM_SEG2, [1, 2])
        for i in range(2):
            m.insertNextCell(NORM_SEG2, [2, 3])
        for i in range(2):
            m.insertNextCell(NORM_SEG2, [3, 4])
        #
        bexp = DataArrayInt([0, 1, 2, 3, 4, 5, 6, 9, 10])
        m1 = m[bexp]
        #
        a, b = m.areCellsIncludedIn(m1, 0)
        self.assertTrue(a)
        self.assertTrue(b.isEqual(DataArrayInt([2, 2, 2, 6, 6, 6, 6, 10, 10])))
        #
        bexp2 = DataArrayInt([0, 1, 2, 3, 4, 0, 6, 9, 10])
        m2 = m[bexp2]
        a, b = m.areCellsIncludedIn(m2, 0)
        self.assertTrue(a)
        self.assertTrue(b.isEqual(DataArrayInt([2, 2, 2, 6, 6, 2, 6, 10, 10])))

    def testSkyLineArrayThreshold(self):
        x = DataArrayInt(
            [0, 1, 2, 11, 12, 13, 3, 4, 5, 6, 14, 15, 16, 17, 9, 10, 18, 19]
        )
        xi = DataArrayInt([0, 6, 14, 18])
        sk = MEDCouplingSkyLineArray(xi, x)
        lsk, rsk = sk.thresholdPerPack(11)
        self.assertTrue(
            lsk.getValuesArray().isEqual(DataArrayInt([0, 1, 2, 3, 4, 5, 6, 9, 10]))
        )
        self.assertTrue(lsk.getIndexArray().isEqual(DataArrayInt([0, 3, 7, 9])))
        self.assertTrue(
            rsk.getValuesArray().isEqual(
                DataArrayInt([11, 12, 13, 14, 15, 16, 17, 18, 19])
            )
        )
        self.assertTrue(rsk.getIndexArray().isEqual(DataArrayInt([0, 3, 7, 9])))

    def testPenta18GaussNE(self):
        conn = [1, 0, 2, 4, 3, 5, 6, 7, 8, 9, 13, 14, 11, 10, 15, 12, 17, 16]
        coo = DataArrayDouble(
            [
                (27.237499999999997, 9.8, 0.0),
                (26.974999999999994, 9.8, 0.0),
                (27.111517409545634, 9.532083869948877, 0.0),
                (27.237499999999997, 9.8, 0.5000000000000001),
                (26.974999999999994, 9.8, 0.5000000000000002),
                (27.111517409545634, 9.532083869948877, 0.5),
                (27.106249999999996, 9.8, 0.0),
                (27.17450870477282, 9.666041934974439, 0.0),
                (27.04325870477281, 9.666041934974439, 0.0),
                (27.106249999999996, 9.8, 0.5000000000000001),
                (27.237499999999997, 9.8, 0.25),
                (26.974999999999994, 9.8, 0.2500000000000001),
                (27.106249999999996, 9.8, 0.2500000000000001),
                (27.174508704772816, 9.666041934974439, 0.5),
                (27.043258704772814, 9.666041934974439, 0.5000000000000001),
                (27.111517409545634, 9.532083869948877, 0.25),
                (27.043258704772818, 9.666041934974436, 0.25000000000000006),
                (27.174508704772816, 9.666041934974436, 0.25),
            ]
        )
        m = MEDCouplingUMesh("mesh", 3)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_PENTA18, conn)
        f = MEDCouplingFieldDouble(ON_GAUSS_NE)
        f.setMesh(m)
        f.setArray(DataArrayDouble(18 * [0.0]))
        self.assertTrue(f.getLocalizationOfDiscr().isEqual(coo[conn], 1e-10))

    def testDADEigenValuesPb(self):
        """EDF22126 : eigen values with Identity matrix returned nan. Now it returns correct eigen values 1.0"""
        valuesExp = DataArrayDouble(
            [(1.0, 1.0, 1.0), (2.0, -1.0, 0.0), (2.0, 0.0, 1.0), (3.0, 0.0, 0.0)]
        )
        d = DataArrayDouble(4, 6)
        for i, (
            v0,
            v1,
            v2,
            v3,
            v4,
            v5,
        ) in enumerate(
            [
                (1, 1, 1, 0, 0, 0),
                (1, 0, 0, 1, 0, 1),
                (1, 1, 1, 0, 1, 0),
                (1, 1, 1, 1, 1, 1),
            ]
        ):
            d[i] = [v0, v1, v2, v3, v4, v5]
        self.assertTrue(d.eigenValues().isEqual(valuesExp, 1e-12))
        pass

    def testBugOnReprOf1SGTUMesh(self):
        """Non reg bug on repr of empty MEDCoupling1SGTUMesh instance"""
        m = MEDCoupling1SGTUMesh()
        m.simpleRepr()
        str(m)
        m.advancedRepr()
        repr(m)
        m = MEDCoupling1DGTUMesh()
        m.simpleRepr()
        str(m)
        m.advancedRepr()
        repr(m)

    def testCheckGeomConsistency0(self):
        """Test of MEDCouplingUMesh.checkGeomConsistency"""
        m = MEDCouplingUMesh("", 2)
        m.setCoords(DataArrayDouble([(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)]))
        m.allocateCells()
        m.insertNextCell(NORM_TRI6, [0, 1, 2, 3, 4, 5])
        m.insertNextCell(NORM_TRI6, [0, 1, 3, 3, 4, 5])
        m.checkConsistency()
        self.assertRaises(
            InterpKernelException, m.checkGeomConsistency
        )  # cell1 is incorrect because node 3 is repeated twice
        m.getNodalConnectivity()[10] = 2  # replace 3 by 2 for cell#1 to fix the problem
        m.checkConsistency()
        m.checkGeomConsistency()  # now m is OK

    def testInt32Int64Arr0(self):
        n = 30
        arr = DataArrayInt32(n)
        arr.iota()
        arr.rearrange(3)
        comps = ["a", "bb", "ccc"]
        name = "aaa"
        arr.setName(name)
        arr.setInfoOnComponents(comps)
        self.assertEqual(arr.accumulate(), [135, 145, 155])
        arr2 = arr.convertToInt64Arr()  # test is here
        self.assertEqual(arr2.accumulate(), [135, 145, 155])
        self.assertTrue(isinstance(arr2, DataArrayInt64))
        self.assertEqual(arr2.getName(), name)
        self.assertEqual(arr2.getInfoOnComponents(), comps)
        arr3 = arr2.convertToInt32Arr()  # test is here
        self.assertEqual(arr3.accumulate(), [135, 145, 155])
        self.assertTrue(isinstance(arr3, DataArrayInt32))
        self.assertEqual(arr3.getName(), name)
        self.assertEqual(arr3.getInfoOnComponents(), comps)
        self.assertTrue(arr3.isEqual(arr))

    def testComputeMeshCenterOfMass0(self):
        # 2D
        arr = DataArrayDouble(5)
        arr.iota()
        m = MEDCouplingCMesh()
        m.setCoords(arr, arr)
        m = m.buildUnstructured()
        self.assertTrue(
            m.computeMeshCenterOfMass().isEqual(DataArrayDouble([2, 2], 1, 2), 1e-12)
        )
        # 3D
        m = MEDCouplingCMesh()
        m.setCoords(arr, arr, arr)
        m = m.buildUnstructured()
        self.assertTrue(
            m.computeMeshCenterOfMass().isEqual(DataArrayDouble([2, 2, 2], 1, 3), 1e-12)
        )

    def testBugPenta15_0(self):
        """
        Non regression test from Roberto Da Via pointing error in connectivity of 5th sub face of penta15 cell.
        """
        coo = DataArrayDouble(
            [
                (0, 1, 1),
                (0, 0, 1),
                (1, 0, 1),
                (0, 1, 0),
                (0, 0, 0),
                (1, 0, 0),
                (0, 0.5, 1),
                (0.5, 0, 1),
                (0.5, 0.5, 1),
                (0, 0.5, 0),
                (0.5, 0, 0),
                (0.5, 0.5, 0),
                (0, 1, 0.5),
                (0, 0, 0.5),
                (1, 0, 0.5),
            ]
        )

        m = MEDCouplingUMesh("penta15", 3)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_PENTA15, list(range(15)))
        bm = m.buildBoundaryMesh(True)
        # bm.writeVTK("boundary.vtu")
        conn_expected = [
            [6, 0, 1, 2, 6, 7, 8],
            [6, 3, 5, 4, 11, 10, 9],
            [8, 0, 3, 4, 1, 12, 9, 13, 6],
            [8, 1, 4, 5, 2, 13, 10, 14, 7],
            [8, 2, 5, 3, 0, 14, 11, 12, 8],  # old = [8, 2, 4, 5, 0, 14, 11, 12, 8]
        ]
        self.assertTrue(
            bm.getNodalConnectivity().isEqual(DataArrayInt(sum(conn_expected, [])))
        )

    def testBugWithPolyhedInterpWithMoreThan255Nodes(self):
        """
        [EDF25207] : Check interpolation containing polyhedron with more than 255 nodes is OK at bbox computation stage
        """
        n = 8
        arr = DataArrayDouble(n)
        arr.iota()
        m = MEDCouplingCMesh()
        m.setCoords(arr, arr, arr)
        m = m.buildUnstructured()
        skin = m.computeSkin()
        skin.zipCoords()
        # check that skin contains more than 2**8-1 node to reveal bug
        self.assertTrue(skin.getNumberOfNodes() > 255)
        # Build a single polyhedron cell from skin
        skin1 = MEDCoupling1SGTUMesh(skin)
        conn = skin1.getNodalConnectivity()
        conn.rearrange(
            MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(skin1.getCellModelEnum())
        )
        connPolyhed = conn.changeNbOfComponents(
            MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(skin1.getCellModelEnum())
            + 1,
            -1,
        )
        connPolyhed.rearrange(1)
        connPolyhed.popBackSilent()
        meshSinglePolyhed = MEDCouplingUMesh("", 3)
        meshSinglePolyhed.allocateCells()
        meshSinglePolyhed.insertNextCell(NORM_POLYHED, connPolyhed.getValues())
        meshSinglePolyhed.setCoords(skin1.getCoords())

        rem = MEDCouplingRemapper()
        rem.prepare(meshSinglePolyhed, m, "P0P0")
        res = rem.getCrudeMatrix()
        self.assertTrue(all([len(elt) == 1 for elt in res]))
        self.assertTrue(all([elt[0] > 0.99 and elt[0] < 1.01 for elt in res]))

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(), "requires numpy")
    def testShapeFuncAndDerivative0(self):
        """
        Test values returned by MEDCoupling on HEXA27 element of shape function and its derivatives.
        See https://www.code-aster.org/V2/doc/v10/fr/man_r/r3/r3.01.01.pdf
        """
        import numpy as np

        ref_coords_hexa27_med = [
            [-1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, -1.0, 1.0],
            [-1.0, 0.0, -1.0],
            [0.0, 1.0, -1.0],
            [1.0, 0.0, -1.0],
            [0.0, -1.0, -1.0],
            [-1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, -1.0, 1.0],
            [-1.0, -1.0, 0.0],
            [-1.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
            [0.0, 0.0, -1.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
        ]

        def coor2index(coor):
            zeMap = {-1.0: 0, 0.0: 2, 1.0: 1}
            return zeMap[coor]

        vcoor2index = np.vectorize(coor2index)
        node2ijk_hexa27_med = vcoor2index(np.array(ref_coords_hexa27_med))

        def N_1d_quad(x):
            return np.array([-0.5 * x * (1 - x), 0.5 * x * (x + 1), 1.0 - x * x])

        def N_3d_hexa27(x, i, j, k):
            return N_1d_quad(x[0])[i] * N_1d_quad(x[1])[j] * N_1d_quad(x[2])[k]

        def N_hexa27(x):
            return np.array(
                [N_3d_hexa27(x, *node2ijk_hexa27_med[node, :]) for node in range(27)]
            )

        # Implementing shape function derivatives
        def diff_N_1d_quad(x):
            return np.array([x - 0.5, x + 0.5, -2.0 * x])

        def diff_N_3d_hexa27(x, i, j, k):
            return np.array(
                [
                    diff_N_1d_quad(x[0])[i] * N_1d_quad(x[1])[j] * N_1d_quad(x[2])[k],
                    N_1d_quad(x[0])[i] * diff_N_1d_quad(x[1])[j] * N_1d_quad(x[2])[k],
                    N_1d_quad(x[0])[i] * N_1d_quad(x[1])[j] * diff_N_1d_quad(x[2])[k],
                ]
            )

        def diff_N_hexa27(x):
            return np.array(
                [
                    diff_N_3d_hexa27(x, *node2ijk_hexa27_med[node, :])
                    for node in range(27)
                ]
            )

        # computation of ref values
        posInRefCoord = [-0.85685375, -0.90643355, -0.90796825]
        ref = N_hexa27(np.array(posInRefCoord))
        ref2 = diff_N_hexa27(np.array(posInRefCoord))
        # computation using MEDCoupling
        gl = MEDCouplingGaussLocalization(
            NORM_HEXA27, sum(ref_coords_hexa27_med, []), posInRefCoord, [1]
        )
        mcShapeFunc = gl.getShapeFunctionValues()
        mcShapeFunc.rearrange(1)
        self.assertTrue(mcShapeFunc.isEqual(DataArrayDouble(ref), 1e-12))

        mvDevOfShapeFunc = gl.getDerivativeOfShapeFunctionValues()
        mvDevOfShapeFunc.rearrange(1)
        ref2_mc = DataArrayDouble(ref2)
        ref2_mc.rearrange(1)
        self.assertTrue(mvDevOfShapeFunc.isEqual(ref2_mc, 1e-12))

    def testShapeFuncAndDerivative1(self):
        """
        This test focus
        """

        def GetShapeFunc(ref_coord, vec):
            gl3 = MEDCouplingGaussLocalization(gt, sum(ref_coord, []), vec, [1])
            funVal = gl3.getShapeFunctionValues()
            funVal.rearrange(1)
            return funVal

        def GetDerivative(ref_coord, vec):
            gl3 = MEDCouplingGaussLocalization(gt, sum(ref_coord, []), vec, [1])
            funVal = gl3.getDerivativeOfShapeFunctionValues()
            return funVal

        vec = [-0.85685375, -0.90643355, -0.90796825]
        eps = 1e-6
        # 3D cells
        for gt in [
            NORM_TETRA4,
            NORM_TETRA10,
            NORM_HEXA8,
            NORM_PENTA6,
            NORM_PYRA5,
            NORM_PYRA13,
            NORM_PENTA15,
            NORM_PENTA6,
            NORM_PENTA18,
            NORM_HEXA20,
            NORM_HEXA27,
        ]:  # type of cell for which derivatives are implemented
            ref_coord = [
                list(elt)
                for elt in MEDCouplingGaussLocalization.GetDefaultReferenceCoordinatesOf(
                    gt
                ).getValuesAsTuple()
            ]

            der_computed = GetDerivative(ref_coord, vec)
            der_computed.rearrange(3)

            der_deduced = (
                GetShapeFunc(ref_coord, [vec[0] + eps, vec[1], vec[2]])
                - GetShapeFunc(ref_coord, vec)
            ) / eps
            delta_X = der_computed[:, 0] - der_deduced
            delta_X.abs()
            self.assertTrue(delta_X.findIdsNotInRange(-1e-4, +1e-4).empty())

            der_deduced = (
                GetShapeFunc(ref_coord, [vec[0], vec[1] + eps, vec[2]])
                - GetShapeFunc(ref_coord, vec)
            ) / eps
            delta_Y = der_computed[:, 1] - der_deduced
            delta_Y.abs()
            self.assertTrue(delta_Y.findIdsNotInRange(-1e-5, +1e-5).empty())

            der_deduced = (
                GetShapeFunc(ref_coord, [vec[0], vec[1], vec[2] + eps])
                - GetShapeFunc(ref_coord, vec)
            ) / eps
            delta_Z = der_computed[:, 2] - der_deduced
            delta_Z.abs()
            self.assertTrue(delta_Z.findIdsNotInRange(-1e-5, +1e-5).empty())

        for gt, ref_coord in [
            (
                NORM_TETRA4,
                [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
            ),
            (
                NORM_TETRA10,
                [
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 0.5, 0.0],
                    [0.0, 0.0, 0.5],
                    [0.0, 0.5, 0.5],
                    [0.5, 0.5, 0.0],
                    [0.5, 0.0, 0.0],
                    [0.5, 0.0, 0.5],
                ],
            ),
            (
                NORM_HEXA8,
                [
                    [-1.0, -1.0, -1.0],
                    [-1.0, 1.0, -1.0],
                    [1.0, 1.0, -1.0],
                    [1.0, -1.0, -1.0],
                    [-1.0, -1.0, 1.0],
                    [-1.0, 1.0, 1.0],
                    [1.0, 1.0, 1.0],
                    [1.0, -1.0, 1.0],
                ],
            ),
            (
                NORM_HEXA8,
                [
                    [-1.0, 1.0, 0.0],
                    [-1.0, -1.0, 0.0],
                    [1.0, -1.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                ],
            ),
            (
                NORM_HEXA8,
                [
                    [-1.0, -1.0, 0.0],
                    [-1.0, 1.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [1.0, -1.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                ],
            ),
            (
                NORM_PENTA6,
                [
                    [-1.0, 1.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [-1.0, -0.0, 1.0],
                    [1.0, 1.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [1.0, 0.0, 1.0],
                ],
            ),
            (
                NORM_PENTA6,
                [
                    [-1.0, 1.0, 0.0],
                    [-1.0, -1.0, 0.0],
                    [1.0, -1.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                ],
            ),
            (
                NORM_PENTA6,
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                ],
            ),
            (
                NORM_PYRA5,
                [
                    [1.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
            ),
            (
                NORM_PYRA13,
                [
                    [1.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [0.5, -0.5, 0.0],
                    [-0.5, -0.5, 0.0],
                    [-0.5, 0.5, 0.0],
                    [0.5, 0.5, 0.0],
                    [0.5, 0.0, 0.5],
                    [0.0, -0.5, 0.5],
                    [-0.5, 0.0, 0.5],
                    [0.0, 0.5, 0.5],
                ],
            ),
            (
                NORM_PENTA15,
                [
                    [-1.0, 1.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [-1.0, -0.0, 1.0],
                    [1.0, 1.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [1.0, 0.0, 1.0],
                    [-1.0, 0.5, 0.0],
                    [-1.0, 0.0, 0.5],
                    [-1.0, 0.5, 0.5],
                    [1.0, 0.5, 0.0],
                    [1.0, 0.0, 0.5],
                    [1.0, 0.5, 0.5],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
            ),
            (
                NORM_PENTA18,
                [
                    [-1.0, 1.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [-1.0, -0.0, 1.0],
                    [1.0, 1.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [1.0, 0.0, 1.0],
                    [-1.0, 0.5, 0.0],
                    [-1.0, 0.0, 0.5],
                    [-1.0, 0.5, 0.5],
                    [1.0, 0.5, 0.0],
                    [1.0, 0.0, 0.5],
                    [1.0, 0.5, 0.5],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [0.0, 0.5, 0.0],
                    [0.0, 0.0, 0.5],
                    [0.0, 0.5, 0.5],
                ],
            ),
            (
                NORM_HEXA20,
                [
                    [-1.0, -1.0, -1.0],
                    [-1.0, 1.0, -1.0],
                    [1.0, 1.0, -1.0],
                    [1.0, -1.0, -1.0],
                    [-1.0, -1.0, 1.0],
                    [-1.0, 1.0, 1.0],
                    [1.0, 1.0, 1.0],
                    [1.0, -1.0, 1.0],
                    [-1.0, 0.0, -1.0],
                    [0.0, 1.0, -1.0],
                    [1.0, 0.0, -1.0],
                    [0.0, -1.0, -1.0],
                    [-1.0, 0.0, 1.0],
                    [0.0, 1.0, 1.0],
                    [1.0, 0.0, 1.0],
                    [0.0, -1.0, 1.0],
                    [-1.0, -1.0, 0.0],
                    [-1.0, 1.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [1.0, -1.0, 0.0],
                ],
            ),
        ]:  # type of cell for which derivatives are implemented
            der_computed = GetDerivative(ref_coord, vec)
            der_computed.rearrange(3)

            der_deduced = (
                GetShapeFunc(ref_coord, [vec[0] + eps, vec[1], vec[2]])
                - GetShapeFunc(ref_coord, vec)
            ) / eps
            delta_X = der_computed[:, 0] - der_deduced
            delta_X.abs()
            self.assertTrue(delta_X.findIdsNotInRange(-1e-4, +1e-4).empty())

            der_deduced = (
                GetShapeFunc(ref_coord, [vec[0], vec[1] + eps, vec[2]])
                - GetShapeFunc(ref_coord, vec)
            ) / eps
            delta_Y = der_computed[:, 1] - der_deduced
            delta_Y.abs()
            self.assertTrue(delta_Y.findIdsNotInRange(-1e-5, +1e-5).empty())

            der_deduced = (
                GetShapeFunc(ref_coord, [vec[0], vec[1], vec[2] + eps])
                - GetShapeFunc(ref_coord, vec)
            ) / eps
            delta_Z = der_computed[:, 2] - der_deduced
            delta_Z.abs()
            self.assertTrue(delta_Z.findIdsNotInRange(-1e-5, +1e-5).empty())

        # 2D cells
        vec = [0.64, 0.2]

        for gt in [NORM_QUAD4, NORM_QUAD8, NORM_QUAD9, NORM_TRI3, NORM_TRI6, NORM_TRI7]:
            ref_coord = [
                list(elt)
                for elt in MEDCouplingGaussLocalization.GetDefaultReferenceCoordinatesOf(
                    gt
                ).getValuesAsTuple()
            ]

            der_computed = GetDerivative(ref_coord, vec)
            der_computed.rearrange(2)

            der_deduced = (
                GetShapeFunc(ref_coord, [vec[0] + eps, vec[1]])
                - GetShapeFunc(ref_coord, vec)
            ) / eps
            delta_X = der_computed[:, 0] - der_deduced
            delta_X.abs()
            self.assertTrue(delta_X.findIdsNotInRange(-1e-5, +1e-5).empty())

            der_deduced = (
                GetShapeFunc(ref_coord, [vec[0], vec[1] + eps])
                - GetShapeFunc(ref_coord, vec)
            ) / eps
            delta_Y = der_computed[:, 1] - der_deduced
            delta_Y.abs()
            self.assertTrue(delta_Y.findIdsNotInRange(-1e-4, +1e-4).empty())

        # B version of TRI6, QUAD4 and QUAD8
        for gt, ref_coord in [
            (NORM_TRI3, [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
            (
                NORM_TRI6,
                [
                    [0.0, 0.0],
                    [1.0, 0.0],
                    [0.0, 1.0],
                    [0.5, 0.0],
                    [0.5, 0.5],
                    [0.0, 0.5],
                ],
            ),
            (NORM_QUAD4, [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]]),
            (NORM_QUAD4, [[-1.0, -1.0], [-1.0, 1.0], [1.0, 1.0], [1.0, -1.0]]),
            (NORM_QUAD4, [[-1.0, 0.0], [1.0, 0.0], [0.0, 0.0], [0.0, 0.0]]),
            (
                NORM_QUAD8,
                [
                    [-1.0, -1.0],
                    [1.0, -1.0],
                    [1.0, 1.0],
                    [-1.0, 1.0],
                    [0.0, -1.0],
                    [1.0, 0.0],
                    [0.0, 1.0],
                    [-1.0, 0.0],
                ],
            ),
        ]:
            der_computed = GetDerivative(ref_coord, vec)
            der_computed.rearrange(2)

            der_deduced = (
                GetShapeFunc(ref_coord, [vec[0] + eps, vec[1]])
                - GetShapeFunc(ref_coord, vec)
            ) / eps
            delta_X = der_computed[:, 0] - der_deduced
            delta_X.abs()
            self.assertTrue(delta_X.findIdsNotInRange(-1e-5, +1e-5).empty())

            der_deduced = (
                GetShapeFunc(ref_coord, [vec[0], vec[1] + eps])
                - GetShapeFunc(ref_coord, vec)
            ) / eps
            delta_Y = der_computed[:, 1] - der_deduced
            delta_Y.abs()
            self.assertTrue(delta_Y.findIdsNotInRange(-1e-5, +1e-5).empty())

        # 1D cells
        vec = [0.64]

        for gt in [NORM_SEG2, NORM_SEG3, NORM_SEG4]:
            ref_coord = [
                list(elt)
                for elt in MEDCouplingGaussLocalization.GetDefaultReferenceCoordinatesOf(
                    gt
                ).getValuesAsTuple()
            ]

            der_computed = GetDerivative(ref_coord, vec)
            der_computed.rearrange(1)

            der_deduced = (
                GetShapeFunc(ref_coord, [vec[0] + eps]) - GetShapeFunc(ref_coord, vec)
            ) / eps
            delta_X = der_computed[:, 0] - der_deduced
            delta_X.abs()
            self.assertTrue(delta_X.findIdsNotInRange(-1e-5, +1e-5).empty())

        # B version of SEG2
        for gt, ref_coord in [(NORM_SEG2, [[0.0], [1.0]])]:
            der_computed = GetDerivative(ref_coord, vec)
            der_computed.rearrange(1)

            der_deduced = (
                GetShapeFunc(ref_coord, [vec[0] + eps]) - GetShapeFunc(ref_coord, vec)
            ) / eps
            delta_X = der_computed[:, 0] - der_deduced
            delta_X.abs()
            self.assertTrue(delta_X.findIdsNotInRange(-1e-5, +1e-5).empty())

    def testComputeTriangleHeight0(self):
        arr = DataArrayDouble([0, 1])
        m = MEDCouplingCMesh()
        m.setCoords(arr, arr)
        m = m.buildUnstructured()
        m.simplexize(0)
        m = MEDCoupling1SGTUMesh(m)
        res = m.computeTriangleHeight()
        expected = DataArrayDouble(
            [(1.0, 1.0, sqrt(2) / 2.0), (sqrt(2) / 2.0, 1.0, 1.0)]
        )
        self.assertTrue(res.isEqual(expected, 1e-12))
        m.changeSpaceDimension(3, 100)
        res2 = m.computeTriangleHeight()
        self.assertTrue(res2.isEqual(expected, 1e-12))
        expected2 = DataArrayDouble([sqrt(2) / 2.0, sqrt(2) / 2.0])
        self.assertTrue(res2.minPerTuple().isEqual(expected2, 1e-12))

    def testComputeTriangleHeight1(self):
        m = MEDCouplingUMesh("mesh", 2)
        m.setCoords(DataArrayDouble([(0, 0, 0), (0, 0, 0), (10, 0, 0)]))
        m.allocateCells()
        m.insertNextCell(NORM_TRI3, [0, 1, 2])
        m = MEDCoupling1SGTUMesh(m)
        res = m.computeTriangleHeight()
        expected = DataArrayDouble([(10, 0, 0)])
        self.assertTrue(res.isEqual(expected, 1e-12))

    def testDAILocateComponentId0(self):
        arr = DataArrayInt(
            [
                (0, 1, 2),
                (3, 4, 5),
                (6, 2, 3),
                (7, 8, 9),
                (9, 0, 10),
                (11, 12, 13),
                (14, 5, 11),
                (15, 16, 17),
            ]
        )
        valToSearchIntoTuples = DataArrayInt([1, 4, 6, 8, 10, 12, 14, 16, 17])
        tupleIdHint = DataArrayInt([0, 1, 2, 3, 4, 5, 6, 7, 7])
        ret = arr.locateComponentId(valToSearchIntoTuples, tupleIdHint)
        self.assertTrue(ret.isEqual(DataArrayInt([1, 1, 0, 1, 2, 1, 0, 1, 2])))
        pass

    def testMeasureOnGaussPtMeshDimNotEqualSpaceDim0(self):
        """
        [EDF26877] : This test focuses on computation of measure field on field on Gauss Point in the special case where SpaceDim
        are not eqaul to the meshDim.
        """
        seg2 = MEDCouplingUMesh("mesh", 1)
        seg2.setCoords(DataArrayDouble([(3, 3), (4, 4)]))
        seg2.allocateCells()
        seg2.insertNextCell(NORM_SEG2, [0, 1])
        fff = MEDCouplingFieldDouble.New(ON_GAUSS_PT)
        fff.setName("CH1RB")
        fff.setNature(IntensiveMaximum)
        fff.setMesh(seg2)
        fff.setGaussLocalizationOnCells([0], [0.0, 1.0], [0.333333333333333], [1.0])
        disc = fff.getDiscretization()
        # spaceDim = 2 meshDim = 1
        self.assertTrue(
            disc.getMeasureField(seg2, True)
            .getArray()
            .isEqual(DataArrayDouble([sqrt(2.0)]), 1e-10)
        )
        # spaceDim = 3 meshDim = 1
        seg2.setCoords(DataArrayDouble([(3, 3, 3), (4, 4, 4)]))
        self.assertTrue(
            disc.getMeasureField(seg2, True)
            .getArray()
            .isEqual(DataArrayDouble([sqrt(3.0)]), 1e-10)
        )
        # spaceDim = 3 meshDim = 2
        tri = MEDCouplingUMesh("mesh", 2)
        tri.setCoords(DataArrayDouble([(0, 0, 0), (1, 1, 0), (2, 2, 2)]))
        tri.allocateCells()
        tri.insertNextCell(NORM_TRI3, [0, 1, 2])
        fff = MEDCouplingFieldDouble.New(ON_GAUSS_PT)
        fff.setName("CH1RB")
        fff.setNature(IntensiveMaximum)
        fff.setMesh(tri)
        fff.setGaussLocalizationOnCells(
            list(range(0, 1)),
            [0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
            [0.3333333333333333, 0.3333333333333333],
            [0.5],
        )
        disc = fff.getDiscretization()
        self.assertTrue(
            disc.getMeasureField(tri, True)
            .getArray()
            .isEqual(tri.getMeasureField(True).getArray(), 1e-10)
        )
        pass

    def testUMeshExplodeMeshTo(self):
        """
        [EDF27988] : implementation of reduceToCells implies implementation of MEDCouplingUMesh.explodeMeshTo
        """
        arr = DataArrayDouble(5)
        arr.iota()
        m = MEDCouplingCMesh()
        m.setCoords(arr, arr, arr)
        m = m.buildUnstructured()
        m1 = m[::2]
        m2 = m[1::2]
        m1.simplexize(PLANAR_FACE_5)
        m = MEDCouplingUMesh.MergeUMeshesOnSameCoords([m1, m2])
        mE1 = m.explodeMeshTo(-1)
        ref1 = m.buildDescendingConnectivity()
        mE2 = m.explodeMeshTo(-2)
        ref2 = m.explode3DMeshTo1D()
        mE3 = m.explodeMeshTo(-3)
        self.assertTrue(len(mE1) == 5)
        self.assertTrue(
            mE1[0].getNodalConnectivity().isEqual(ref1[0].getNodalConnectivity())
        )
        self.assertTrue(
            mE1[0]
            .getNodalConnectivityIndex()
            .isEqual(ref1[0].getNodalConnectivityIndex())
        )
        self.assertTrue(
            mE1[0].getCoords().getHiddenCppPointer()
            == m.getCoords().getHiddenCppPointer()
        )
        for i in range(1, 5):
            self.assertTrue(mE1[i].isEqual(ref1[i]))
        #
        self.assertTrue(len(mE2) == 5)
        self.assertTrue(
            mE2[0].getNodalConnectivity().isEqual(ref2[0].getNodalConnectivity())
        )
        self.assertTrue(
            mE2[0]
            .getNodalConnectivityIndex()
            .isEqual(ref2[0].getNodalConnectivityIndex())
        )
        self.assertTrue(
            mE2[0].getCoords().getHiddenCppPointer()
            == m.getCoords().getHiddenCppPointer()
        )
        for i in range(1, 5):
            self.assertTrue(mE2[i].isEqual(ref2[i]))
        #
        self.assertTrue(mE3[0].getMeshDimension() == 0)
        self.assertTrue(mE3[0].getNumberOfCells() == mE3[0].getNumberOfNodes())
        a, b = m.getReverseNodalConnectivity()
        self.assertTrue(mE3[3].isEqual(a) and mE3[4].isEqual(b))
        ref3_2 = m.getNodalConnectivityIndex().deltaShiftIndex() - 1
        ref3_2.computeOffsetsFull()
        self.assertTrue(ref3_2.isEqual(mE3[2]))
        tmp = m.getNodalConnectivityIndex().deepCopy()
        tmp.popBackSilent()
        tmp = tmp.buildComplement(len(m.getNodalConnectivity()))
        ref3_1 = m.getNodalConnectivity()[tmp]
        self.assertTrue(ref3_1.isEqual(mE3[1]))
        #
        cellsInPolyh = [37, 160]
        polyh = m[cellsInPolyh]
        polyh.convertAllToPoly()
        m[cellsInPolyh] = polyh
        pE3 = m.explodeMeshTo(-3)
        self.assertTrue(pE3[0].getMeshDimension() == 0)
        self.assertTrue(pE3[0].getNumberOfCells() == pE3[0].getNumberOfNodes())
        a, b = m.getReverseNodalConnectivity()
        self.assertTrue(pE3[3].isEqual(a) and pE3[4].isEqual(b))
        self.assertTrue(pE3[2].isEqual(mE3[2]))  # indexed arrays are the same

        ref_a, ref_b = DataArrayInt.ExtractFromIndexedArrays(
            DataArrayInt(cellsInPolyh).buildComplement(m.getNumberOfCells()),
            mE3[1],
            mE3[2],
        )
        a, b = DataArrayInt.ExtractFromIndexedArrays(
            DataArrayInt(cellsInPolyh).buildComplement(m.getNumberOfCells()),
            pE3[1],
            pE3[2],
        )
        self.assertTrue(ref_a.isEqual(a))
        self.assertTrue(ref_b.isEqual(b))
        for cell in cellsInPolyh:
            ref_c, ref_d = DataArrayInt.ExtractFromIndexedArrays(cell, mE3[1], mE3[2])
            ref_c.sort()
            c, d = DataArrayInt.ExtractFromIndexedArrays(cell, pE3[1], pE3[2])
            self.assertTrue(ref_c.isEqual(c))
            self.assertTrue(ref_d.isEqual(d))

    def testGetCellContainingPointRelativeEps(self):
        """
        See EDF27860 : This test checks that detection of point inside a cell works by normalizing cell around origin with factor equal to the max delta of bbox along axis X, Y or Z.
        """
        # in this test cell is vuluntary far from origin {15260.775604514516, 11197.646906189217, 14187.820484060947}
        # and caracteritic size is ~ 1500
        coo = DataArrayDouble(
            [
                (14724.199858870656, 11928.888084722483, 14442.32726944039),
                (14788.407409534622, 11992.60694822231, 14453.86181555231),
                (15572.175148726046, 10798.586790270576, 14471.54225356788),
                (15643.898717334796, 10853.094666047728, 14477.233802854305),
                (15005.31495255754, 11573.261110174888, 13933.313698681504),
                (15070.29423166349, 11636.377758513776, 13946.650959030132),
                (15797.351350158377, 10466.40572765595, 13965.524190108257),
                (15869.808770928525, 10519.99285973948, 13972.419352086607),
                (15273.866774426142, 11216.458197414971, 13433.169979717744),
                (15340.421031616577, 11277.882145177837, 13446.53598386297),
                (16013.382514001762, 10132.795887638129, 13465.184281842226),
                (16086.979064572806, 10184.802292369684, 13472.147425473782),
            ]
        )
        m = MEDCouplingUMesh("", 3)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_TETRA4, [0, 5, 4, 6])
        m.insertNextCell(NORM_TETRA4, [4, 5, 9, 7])

        ##### See EDF2760 pt is outside cell 0 (6e-4) and 1 (8e-4)
        pt = DataArrayDouble([(15263.41200205526, 11314.957094727113, 13950.0)])
        a, b = m.getCellsContainingPoints(pt, 1e-3)
        self.assertTrue(a.isEqual(DataArrayInt([0, 1])))
        self.assertTrue(b.isEqual(DataArrayInt([0, 2])))

        # by shifting pt by 10 along Z pt in only inside cell # 0
        pt += [0, 0, 10]
        a1, b1 = m.getCellsContainingPoints(pt, 1e-3)
        self.assertTrue(a1.isEqual(DataArrayInt([0])))
        self.assertTrue(b1.isEqual(DataArrayInt([0, 1])))

    def testGetCellContainingPointOnPolyhedronWithPlanarFace(self):
        """
        See CEA spns #40783
        In case of polyhedron with a face defined by several colinear points,
        we must use other non colinear points to be able to define a face from these three points
        to define the relative position of the test point to this face
        """
        eps = 1.0e-12
        coo = DataArrayDouble(
            [
                (0.176, 0.1125, 1.05),
                (0.176, 0.120375, 1.05),
                (0.176, 0.120375, 1.0),
                (0.176, 0.1125, 1.0),
                (0.176000000000000018, 0.12825, 1.05),
                (0.176000000000000018, 0.12825, 1.0),
                (0.207, 0.1125, 1.05),
                (0.207, 0.1125, 1.0),
                (0.207, 0.12825, 1.05),
                (0.207, 0.12825, 1.0),
            ]
        )

        m = MEDCouplingUMesh("Mesh", 3)
        m.setCoords(coo)
        m.allocateCells()
        # put -1 to separate faces connectivity
        # substract -1 from mdump table ids
        m.insertNextCell(
            NORM_POLYHED,
            [
                0,
                1,
                2,
                3,
                -1,
                1,
                4,
                5,
                2,
                -1,
                6,
                7,
                9,
                8,
                -1,
                3,
                7,
                6,
                0,
                -1,
                9,
                5,
                4,
                8,
                -1,
                3,
                2,
                5,
                9,
                7,
                -1,  # PB in this order
                # 7, 3, 2, 5, 9, -1, # OK in this order
                1,
                0,
                6,
                8,
                4,
            ],
        )

        # test point inside the box
        pt_above = (0.2, 0.12, 1.07)
        pt_below = (0.2, 0.12, 0.9)
        pt_inside = (0.2, 0.12, 1.025)
        pts = DataArrayDouble([pt_above, pt_below, pt_inside])
        a, b = m.getCellsContainingPoints(pts, eps)
        self.assertTrue(a.isEqual(DataArrayInt([0])))
        # only the third point is inside
        self.assertTrue(b.isEqual(DataArrayInt([0, 0, 0, 1])))

        # rotate the mesh to see if getCellsContainingPoints works
        # even if point is not inside bounding box
        center = coo[0]
        vector = [1.0, 0.0, 0.0]
        m.rotate(center, vector, -pi / 4.0)

        # test 3 points: above, below and inside
        pt_above = (0.19, 0.09, 1.04)
        pt_below = (0.19, 0.11, 1.02)
        pt_inside = (0.19, 0.10, 1.02)
        pts_rotated = DataArrayDouble([pt_above, pt_below, pt_inside])

        a, b = m.getCellsContainingPoints(pts_rotated, eps)
        self.assertTrue(a.isEqual(DataArrayInt([0])))
        # only the third point is inside
        self.assertTrue(b.isEqual(DataArrayInt([0, 0, 0, 1])))

    def testGetCellContainingPointOnPolyhedronWithPlanarFaceWithManyNodes(self):
        """
        Similar test with many colinear nodes on the planar face
        """
        eps = 1.0e-12
        coo = DataArrayDouble(
            [
                (0.176000000000000018, 0.120375, 1.0),
                (0.176000000000000018, 0.128250, 1.0),
                (0.176000000000000018, 0.136125, 1.0),
                (0.176000000000000018, 0.144, 1.0),
                (0.176000000000000018, 0.151875, 1.0),
                (0.176000000000000018, 0.159750, 1.0),
                (0.176000000000000018, 0.167625, 1.0),
                (0.176000000000000018, 0.1755, 1.0),
                (0.176000000000000018, 0.183375, 1.0),
                (0.176000000000000018, 0.191250, 1.0),
                (0.176000000000000018, 0.199125, 1.0),
                (0.176, 0.207, 1.0),
                (0.207, 0.207, 1.0),
                (0.176, 0.1125, 1.0),
                (0.207, 0.1125, 1.0),
                (0.176, 0.120375, 1.05),
                (0.176000000000000018, 0.128250, 1.05),
                (0.176000000000000018, 0.136125, 1.05),
                (0.176000000000000018, 0.144, 1.05),
                (0.176000000000000018, 0.151875, 1.05),
                (0.176000000000000018, 0.159750, 1.05),
                (0.176000000000000018, 0.167625, 1.05),
                (0.176000000000000018, 0.1755, 1.05),
                (0.176000000000000018, 0.183375, 1.05),
                (0.176000000000000018, 0.191250, 1.05),
                (0.176000000000000018, 0.199125, 1.05),
                (0.176, 0.207, 1.05),
                (0.207, 0.207, 1.05),
                (0.176, 0.1125, 1.05),
                (0.207, 0.1125, 1.05),
            ]
        )

        m = MEDCouplingUMesh("Mesh", 3)
        m.setCoords(coo)
        m.allocateCells()
        # put -1 to separate faces connectivity
        # substract -1 from mdump table ids
        m.insertNextCell(
            NORM_POLYHED,
            [
                13,
                0,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                14,
                -1,  # 1
                29,
                27,
                26,
                25,
                24,
                23,
                22,
                21,
                20,
                19,
                18,
                17,
                16,
                15,
                28,
                -1,  # 2
                14,
                29,
                28,
                13,
                -1,  # 3
                11,
                26,
                27,
                12,
                -1,  # 4
                12,
                27,
                29,
                14,
                -1,  # 5
                13,
                28,
                15,
                0,
                -1,  # 6
                0,
                15,
                16,
                1,
                -1,  # 7
                1,
                16,
                17,
                2,
                -1,  # 8
                2,
                17,
                18,
                3,
                -1,  # 9
                3,
                18,
                19,
                4,
                -1,  # 10
                4,
                19,
                20,
                5,
                -1,  # 11
                5,
                20,
                21,
                6,
                -1,  # 12
                6,
                21,
                22,
                7,
                -1,  # 13
                7,
                22,
                23,
                8,
                -1,  # 14
                8,
                23,
                24,
                9,
                -1,  # 15
                9,
                24,
                25,
                10,
                -1,  # 16
                10,
                25,
                26,
                11,
            ],
        )

        ##### See CEA 40783: error with polyhedrons (box split by on edge on its face)
        pt_above = (0.1915, 0.15975, 1.07)
        pt_below = (0.1915, 0.15975, 0.9)
        pt_inside = (0.1915, 0.15975, 1.025)
        pts = DataArrayDouble([pt_above, pt_below, pt_inside])
        a, b = m.getCellsContainingPoints(pts, eps)
        self.assertTrue(a.isEqual(DataArrayInt([0])))
        # only the third point is inside
        self.assertTrue(b.isEqual(DataArrayInt([0, 0, 0, 1])))

    def testSGTUMeshGetCellsContainingPts1(self):
        """
        EDF29571 : Fix problem of perfomance of GTUMesh::getCellsContainingPts
        """
        arr = DataArrayDouble(31)
        arr.iota()
        m = MEDCouplingCMesh()
        m.setCoords(arr, arr, arr)
        m = m.buildUnstructured()
        m2 = MEDCoupling1SGTUMesh(m)
        pts = m[:100].computeCellCenterOfMass()
        a = datetime.now()
        a1, b1 = m.getCellsContainingPoints(pts, 1e-5)
        b = datetime.now()
        a2, b2 = m2.getCellsContainingPoints(pts, 1e-5)
        c = datetime.now()
        ref = b - a
        toTest = c - a
        self.assertTrue(toTest < 10 * ref)

    def testFuseOfFamilyField0(self):
        """
        EDF30179 : Core algo for family field linked to fusion of entities
        """
        d = DataArrayInt(
            [2, 2, 2, 2, 2, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1]
        )  # 5 x 2 , 4 x 3, 6 x 1

        c = DataArrayInt([])
        ci = DataArrayInt([0])
        #### Case 0 : simplest
        self.assertTrue(
            ci.deltaShiftIndex().empty()
        )  # EDF30179 : extension of deltaShiftIndex to single elt
        a, b, f = d.forThisAsPartitionBuildReduction(c, ci)
        self.assertTrue(a.isEqual(d))
        self.assertTrue(b.empty())
        self.assertTrue(f.isEqual(DataArrayInt([0])))
        #### Case 1 : single fusion
        c = DataArrayInt([3, 6])
        ci = DataArrayInt([0, 2])
        a, b, f = d.forThisAsPartitionBuildReduction(c, ci)
        self.assertTrue(
            a.isEqual(DataArrayInt([2, 2, 2, 4, 2, 3, 3, 3, 1, 1, 1, 1, 1, 1]))
        )
        self.assertTrue(b.isEqual(DataArrayInt([4, 2, 3])))
        self.assertTrue(f.isEqual(DataArrayInt([0, 3])))
        #### Case 2 : single fusion - same partition id
        c = DataArrayInt([6, 7])
        ci = DataArrayInt([0, 2])
        a, b, f = d.forThisAsPartitionBuildReduction(c, ci)
        self.assertTrue(
            a.isEqual(DataArrayInt([2, 2, 2, 2, 2, 3, 4, 3, 1, 1, 1, 1, 1, 1]))
        )
        self.assertTrue(b.isEqual(DataArrayInt([4, 3])))
        self.assertTrue(f.isEqual(DataArrayInt([0, 2])))
        #### Case 3 : multi fusion single tuple
        c = DataArrayInt([2, 7, 3, 6])
        ci = DataArrayInt(
            [0, 2, 4]
        )  # elts (2,7) and (3,6) to merge. These 2 couples refers to partitionIDs (2,3)
        a, b, f = d.forThisAsPartitionBuildReduction(c, ci)
        self.assertTrue(
            a.isEqual(DataArrayInt([2, 2, 4, 4, 2, 3, 3, 1, 1, 1, 1, 1, 1]))
        )
        self.assertTrue(
            b.isEqual(DataArrayInt([4, 2, 3]))
        )  # Fuse element can be located with ID 4
        self.assertTrue(f.isEqual(DataArrayInt([0, 3])))

        #### Case 4 : multi fusion
        c = DataArrayInt([2, 7, 3, 10])
        ci = DataArrayInt([0, 2, 4])
        a, b, f = d.forThisAsPartitionBuildReduction(c, ci)
        self.assertTrue(
            a.isEqual(DataArrayInt([2, 2, 4, 5, 2, 3, 3, 3, 1, 1, 1, 1, 1]))
        )
        self.assertTrue(b.isEqual(DataArrayInt([4, 2, 3, 5, 1, 2])))
        self.assertTrue(f.isEqual(DataArrayInt([0, 3, 6])))

    def testDASortPerTuple1(self):
        """
        EDF30178 : useful method DataArrayInt.sortPerTuple
        """
        arr = DataArrayInt([(1, 2, 3), (5, 4, 6), (9, 8, 7)])
        arr.sortPerTuple(True)
        self.assertTrue(arr.isEqual(DataArrayInt([(1, 2, 3), (4, 5, 6), (7, 8, 9)])))

    def testDAIFindCommonTuples(self):
        """
        EDF30178 : useful method DataArrayInt.findCommonTuples
        """
        arr = DataArrayInt(
            [
                (1, 2, 3),
                (4, 5, 6),
                (-1, 2, 3),
                (3, 2, 1),
                (1, 2, 3),
                (6, 5, 4),
                (4, 5, 6),
                (4, 5, 6),
            ]
        )
        c, ci = arr.findCommonTuples(2)
        self.assertTrue(ci.isEqual(DataArrayInt([0, 2, 5])))
        self.assertTrue(c.isEqual(DataArrayInt([0, 4, 1, 6, 7])))

    def testDAIFromVTKInternalReprOfPolyedra(self):
        """
        EDF31315 : VTK internal representation of polyedra data structure
        """
        faces = DataArrayInt64(
            [
                11,
                4,
                7199,
                6757,
                2950,
                6758,
                5,
                6455,
                1794,
                2400,
                6620,
                7200,
                6,
                6620,
                2400,
                5864,
                2950,
                6757,
                7244,
                6,
                6758,
                2950,
                5864,
                3223,
                6818,
                7245,
                5,
                7246,
                6818,
                3223,
                1794,
                6455,
                4,
                1794,
                3223,
                5864,
                2400,
                4,
                0,
                7246,
                6455,
                7200,
                4,
                0,
                7200,
                6620,
                7244,
                4,
                0,
                7244,
                6757,
                7199,
                4,
                0,
                7199,
                6758,
                7245,
                4,
                0,
                7245,
                6818,
                7246,
                12,
                5,
                6408,
                978,
                1721,
                6441,
                7203,
                4,
                7204,
                7007,
                3987,
                7008,
                6,
                6441,
                1721,
                5813,
                3987,
                7007,
                7247,
                5,
                7248,
                7130,
                4480,
                978,
                6408,
                6,
                7008,
                3987,
                5813,
                4480,
                7130,
                7249,
                4,
                978,
                4480,
                5813,
                1721,
                4,
                1,
                7248,
                6408,
                7203,
                4,
                1,
                7203,
                6441,
                7247,
                4,
                1,
                7247,
                7007,
                7204,
                4,
                1,
                7204,
                7008,
                7249,
                3,
                1,
                7249,
                7130,
                5,
                1,
                7,
                6,
                5,
                4,
            ]
        )
        facesIndex = DataArrayInt([0, 62, 129])
        facesOut, facesIndexOut = DataArrayInt64.FromVTKInternalReprOfPolyedra(
            faces, facesIndex
        )
        m = MEDCoupling1DGTUMesh("mesh", NORM_POLYHED)
        m.setCoords(DataArrayDouble(10000, 3))
        m.setNodalConnectivity(facesOut, facesIndexOut)
        m = m.buildUnstructured()
        self.assertTrue(m.computeNbOfFacesPerCell().isEqual(DataArrayInt([11, 12])))
        connExp = DataArrayInt64(
            [
                31,
                7199,
                6757,
                2950,
                6758,
                -1,
                6455,
                1794,
                2400,
                6620,
                7200,
                -1,
                6620,
                2400,
                5864,
                2950,
                6757,
                7244,
                -1,
                6758,
                2950,
                5864,
                3223,
                6818,
                7245,
                -1,
                7246,
                6818,
                3223,
                1794,
                6455,
                -1,
                1794,
                3223,
                5864,
                2400,
                -1,
                0,
                7246,
                6455,
                7200,
                -1,
                0,
                7200,
                6620,
                7244,
                -1,
                0,
                7244,
                6757,
                7199,
                -1,
                0,
                7199,
                6758,
                7245,
                -1,
                0,
                7245,
                6818,
                7246,
                31,
                6408,
                978,
                1721,
                6441,
                7203,
                -1,
                7204,
                7007,
                3987,
                7008,
                -1,
                6441,
                1721,
                5813,
                3987,
                7007,
                7247,
                -1,
                7248,
                7130,
                4480,
                978,
                6408,
                -1,
                7008,
                3987,
                5813,
                4480,
                7130,
                7249,
                -1,
                978,
                4480,
                5813,
                1721,
                -1,
                1,
                7248,
                6408,
                7203,
                -1,
                1,
                7203,
                6441,
                7247,
                -1,
                1,
                7247,
                7007,
                7204,
                -1,
                1,
                7204,
                7008,
                7249,
                -1,
                1,
                7249,
                7130,
                -1,
                1,
                7,
                6,
                5,
                4,
            ]
        )
        connIndexExp = DataArrayInt([0, 61, 127])
        self.assertTrue(m.getNodalConnectivity().isEqual(connExp))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(connIndexExp))

    def testDAIFromVTK93To94FacesInternaReprOfPolyedra(self):
        """
        EDF31746 : To manage VTK9.3 -> VTK9.4 internal representation of polyhedra datastructure
        """
        faceLocations93 = DataArrayInt([0, 751, 1502, 2753])
        faces93 = DataArrayInt64(
            [
                150,
                4,
                11,
                19,
                67,
                27,
                4,
                27,
                67,
                68,
                28,
                4,
                28,
                68,
                69,
                29,
                4,
                29,
                69,
                70,
                30,
                4,
                30,
                70,
                23,
                13,
                4,
                19,
                20,
                71,
                67,
                4,
                67,
                71,
                72,
                68,
                4,
                68,
                72,
                73,
                69,
                4,
                69,
                73,
                74,
                70,
                4,
                70,
                74,
                24,
                23,
                4,
                20,
                21,
                75,
                71,
                4,
                71,
                75,
                76,
                72,
                4,
                72,
                76,
                77,
                73,
                4,
                73,
                77,
                78,
                74,
                4,
                74,
                78,
                25,
                24,
                4,
                21,
                22,
                79,
                75,
                4,
                75,
                79,
                80,
                76,
                4,
                76,
                80,
                81,
                77,
                4,
                77,
                81,
                82,
                78,
                4,
                78,
                82,
                26,
                25,
                4,
                22,
                12,
                31,
                79,
                4,
                79,
                31,
                32,
                80,
                4,
                80,
                32,
                33,
                81,
                4,
                81,
                33,
                34,
                82,
                4,
                82,
                34,
                14,
                26,
                4,
                15,
                35,
                83,
                39,
                4,
                39,
                83,
                84,
                40,
                4,
                40,
                84,
                85,
                41,
                4,
                41,
                85,
                86,
                42,
                4,
                42,
                86,
                19,
                11,
                4,
                35,
                36,
                87,
                83,
                4,
                83,
                87,
                88,
                84,
                4,
                84,
                88,
                89,
                85,
                4,
                85,
                89,
                90,
                86,
                4,
                86,
                90,
                20,
                19,
                4,
                36,
                37,
                91,
                87,
                4,
                87,
                91,
                92,
                88,
                4,
                88,
                92,
                93,
                89,
                4,
                89,
                93,
                94,
                90,
                4,
                90,
                94,
                21,
                20,
                4,
                37,
                38,
                95,
                91,
                4,
                91,
                95,
                96,
                92,
                4,
                92,
                96,
                97,
                93,
                4,
                93,
                97,
                98,
                94,
                4,
                94,
                98,
                22,
                21,
                4,
                38,
                16,
                43,
                95,
                4,
                95,
                43,
                44,
                96,
                4,
                96,
                44,
                45,
                97,
                4,
                97,
                45,
                46,
                98,
                4,
                98,
                46,
                12,
                22,
                4,
                17,
                47,
                99,
                51,
                4,
                51,
                99,
                100,
                52,
                4,
                52,
                100,
                101,
                53,
                4,
                53,
                101,
                102,
                54,
                4,
                54,
                102,
                35,
                15,
                4,
                47,
                48,
                103,
                99,
                4,
                99,
                103,
                104,
                100,
                4,
                100,
                104,
                105,
                101,
                4,
                101,
                105,
                106,
                102,
                4,
                102,
                106,
                36,
                35,
                4,
                48,
                49,
                107,
                103,
                4,
                103,
                107,
                108,
                104,
                4,
                104,
                108,
                109,
                105,
                4,
                105,
                109,
                110,
                106,
                4,
                106,
                110,
                37,
                36,
                4,
                49,
                50,
                111,
                107,
                4,
                107,
                111,
                112,
                108,
                4,
                108,
                112,
                113,
                109,
                4,
                109,
                113,
                114,
                110,
                4,
                110,
                114,
                38,
                37,
                4,
                50,
                18,
                55,
                111,
                4,
                111,
                55,
                56,
                112,
                4,
                112,
                56,
                57,
                113,
                4,
                113,
                57,
                58,
                114,
                4,
                114,
                58,
                16,
                38,
                4,
                13,
                23,
                115,
                59,
                4,
                59,
                115,
                116,
                60,
                4,
                60,
                116,
                117,
                61,
                4,
                61,
                117,
                118,
                62,
                4,
                62,
                118,
                47,
                17,
                4,
                23,
                24,
                119,
                115,
                4,
                115,
                119,
                120,
                116,
                4,
                116,
                120,
                121,
                117,
                4,
                117,
                121,
                122,
                118,
                4,
                118,
                122,
                48,
                47,
                4,
                24,
                25,
                123,
                119,
                4,
                119,
                123,
                124,
                120,
                4,
                120,
                124,
                125,
                121,
                4,
                121,
                125,
                126,
                122,
                4,
                122,
                126,
                49,
                48,
                4,
                25,
                26,
                127,
                123,
                4,
                123,
                127,
                128,
                124,
                4,
                124,
                128,
                129,
                125,
                4,
                125,
                129,
                130,
                126,
                4,
                126,
                130,
                50,
                49,
                4,
                26,
                14,
                63,
                127,
                4,
                127,
                63,
                64,
                128,
                4,
                128,
                64,
                65,
                129,
                4,
                129,
                65,
                66,
                130,
                4,
                130,
                66,
                18,
                50,
                4,
                11,
                27,
                131,
                42,
                4,
                42,
                131,
                132,
                41,
                4,
                41,
                132,
                133,
                40,
                4,
                40,
                133,
                134,
                39,
                4,
                39,
                134,
                54,
                15,
                4,
                27,
                28,
                135,
                131,
                4,
                131,
                135,
                136,
                132,
                4,
                132,
                136,
                137,
                133,
                4,
                133,
                137,
                138,
                134,
                4,
                134,
                138,
                53,
                54,
                4,
                28,
                29,
                139,
                135,
                4,
                135,
                139,
                140,
                136,
                4,
                136,
                140,
                141,
                137,
                4,
                137,
                141,
                142,
                138,
                4,
                138,
                142,
                52,
                53,
                4,
                29,
                30,
                143,
                139,
                4,
                139,
                143,
                144,
                140,
                4,
                140,
                144,
                145,
                141,
                4,
                141,
                145,
                146,
                142,
                4,
                142,
                146,
                51,
                52,
                4,
                30,
                13,
                59,
                143,
                4,
                143,
                59,
                60,
                144,
                4,
                144,
                60,
                61,
                145,
                4,
                145,
                61,
                62,
                146,
                4,
                146,
                62,
                17,
                51,
                4,
                14,
                34,
                147,
                63,
                4,
                63,
                147,
                148,
                64,
                4,
                64,
                148,
                149,
                65,
                4,
                65,
                149,
                150,
                66,
                4,
                66,
                150,
                55,
                18,
                4,
                34,
                33,
                151,
                147,
                4,
                147,
                151,
                152,
                148,
                4,
                148,
                152,
                153,
                149,
                4,
                149,
                153,
                154,
                150,
                4,
                150,
                154,
                56,
                55,
                4,
                33,
                32,
                155,
                151,
                4,
                151,
                155,
                156,
                152,
                4,
                152,
                156,
                157,
                153,
                4,
                153,
                157,
                158,
                154,
                4,
                154,
                158,
                57,
                56,
                4,
                32,
                31,
                159,
                155,
                4,
                155,
                159,
                160,
                156,
                4,
                156,
                160,
                161,
                157,
                4,
                157,
                161,
                162,
                158,
                4,
                158,
                162,
                58,
                57,
                4,
                31,
                12,
                46,
                159,
                4,
                159,
                46,
                45,
                160,
                4,
                160,
                45,
                44,
                161,
                4,
                161,
                44,
                43,
                162,
                4,
                162,
                43,
                16,
                58,
                150,
                4,
                163,
                167,
                199,
                175,
                4,
                175,
                199,
                200,
                176,
                4,
                176,
                200,
                201,
                177,
                4,
                177,
                201,
                202,
                178,
                4,
                178,
                202,
                171,
                165,
                4,
                167,
                168,
                203,
                199,
                4,
                199,
                203,
                204,
                200,
                4,
                200,
                204,
                205,
                201,
                4,
                201,
                205,
                206,
                202,
                4,
                202,
                206,
                172,
                171,
                4,
                168,
                169,
                207,
                203,
                4,
                203,
                207,
                208,
                204,
                4,
                204,
                208,
                209,
                205,
                4,
                205,
                209,
                210,
                206,
                4,
                206,
                210,
                173,
                172,
                4,
                169,
                170,
                211,
                207,
                4,
                207,
                211,
                212,
                208,
                4,
                208,
                212,
                213,
                209,
                4,
                209,
                213,
                214,
                210,
                4,
                210,
                214,
                174,
                173,
                4,
                170,
                164,
                179,
                211,
                4,
                211,
                179,
                180,
                212,
                4,
                212,
                180,
                181,
                213,
                4,
                213,
                181,
                182,
                214,
                4,
                214,
                182,
                166,
                174,
                4,
                11,
                19,
                215,
                183,
                4,
                183,
                215,
                216,
                184,
                4,
                184,
                216,
                217,
                185,
                4,
                185,
                217,
                218,
                186,
                4,
                186,
                218,
                167,
                163,
                4,
                19,
                20,
                219,
                215,
                4,
                215,
                219,
                220,
                216,
                4,
                216,
                220,
                221,
                217,
                4,
                217,
                221,
                222,
                218,
                4,
                218,
                222,
                168,
                167,
                4,
                20,
                21,
                223,
                219,
                4,
                219,
                223,
                224,
                220,
                4,
                220,
                224,
                225,
                221,
                4,
                221,
                225,
                226,
                222,
                4,
                222,
                226,
                169,
                168,
                4,
                21,
                22,
                227,
                223,
                4,
                223,
                227,
                228,
                224,
                4,
                224,
                228,
                229,
                225,
                4,
                225,
                229,
                230,
                226,
                4,
                226,
                230,
                170,
                169,
                4,
                22,
                12,
                187,
                227,
                4,
                227,
                187,
                188,
                228,
                4,
                228,
                188,
                189,
                229,
                4,
                229,
                189,
                190,
                230,
                4,
                230,
                190,
                164,
                170,
                4,
                27,
                67,
                19,
                11,
                4,
                28,
                68,
                67,
                27,
                4,
                29,
                69,
                68,
                28,
                4,
                30,
                70,
                69,
                29,
                4,
                13,
                23,
                70,
                30,
                4,
                67,
                71,
                20,
                19,
                4,
                68,
                72,
                71,
                67,
                4,
                69,
                73,
                72,
                68,
                4,
                70,
                74,
                73,
                69,
                4,
                23,
                24,
                74,
                70,
                4,
                71,
                75,
                21,
                20,
                4,
                72,
                76,
                75,
                71,
                4,
                73,
                77,
                76,
                72,
                4,
                74,
                78,
                77,
                73,
                4,
                24,
                25,
                78,
                74,
                4,
                75,
                79,
                22,
                21,
                4,
                76,
                80,
                79,
                75,
                4,
                77,
                81,
                80,
                76,
                4,
                78,
                82,
                81,
                77,
                4,
                25,
                26,
                82,
                78,
                4,
                79,
                31,
                12,
                22,
                4,
                80,
                32,
                31,
                79,
                4,
                81,
                33,
                32,
                80,
                4,
                82,
                34,
                33,
                81,
                4,
                26,
                14,
                34,
                82,
                4,
                165,
                171,
                231,
                191,
                4,
                191,
                231,
                232,
                192,
                4,
                192,
                232,
                233,
                193,
                4,
                193,
                233,
                234,
                194,
                4,
                194,
                234,
                23,
                13,
                4,
                171,
                172,
                235,
                231,
                4,
                231,
                235,
                236,
                232,
                4,
                232,
                236,
                237,
                233,
                4,
                233,
                237,
                238,
                234,
                4,
                234,
                238,
                24,
                23,
                4,
                172,
                173,
                239,
                235,
                4,
                235,
                239,
                240,
                236,
                4,
                236,
                240,
                241,
                237,
                4,
                237,
                241,
                242,
                238,
                4,
                238,
                242,
                25,
                24,
                4,
                173,
                174,
                243,
                239,
                4,
                239,
                243,
                244,
                240,
                4,
                240,
                244,
                245,
                241,
                4,
                241,
                245,
                246,
                242,
                4,
                242,
                246,
                26,
                25,
                4,
                174,
                166,
                195,
                243,
                4,
                243,
                195,
                196,
                244,
                4,
                244,
                196,
                197,
                245,
                4,
                245,
                197,
                198,
                246,
                4,
                246,
                198,
                14,
                26,
                4,
                163,
                175,
                247,
                186,
                4,
                186,
                247,
                248,
                185,
                4,
                185,
                248,
                249,
                184,
                4,
                184,
                249,
                250,
                183,
                4,
                183,
                250,
                27,
                11,
                4,
                175,
                176,
                251,
                247,
                4,
                247,
                251,
                252,
                248,
                4,
                248,
                252,
                253,
                249,
                4,
                249,
                253,
                254,
                250,
                4,
                250,
                254,
                28,
                27,
                4,
                176,
                177,
                255,
                251,
                4,
                251,
                255,
                256,
                252,
                4,
                252,
                256,
                257,
                253,
                4,
                253,
                257,
                258,
                254,
                4,
                254,
                258,
                29,
                28,
                4,
                177,
                178,
                259,
                255,
                4,
                255,
                259,
                260,
                256,
                4,
                256,
                260,
                261,
                257,
                4,
                257,
                261,
                262,
                258,
                4,
                258,
                262,
                30,
                29,
                4,
                178,
                165,
                191,
                259,
                4,
                259,
                191,
                192,
                260,
                4,
                260,
                192,
                193,
                261,
                4,
                261,
                193,
                194,
                262,
                4,
                262,
                194,
                13,
                30,
                4,
                166,
                182,
                263,
                195,
                4,
                195,
                263,
                264,
                196,
                4,
                196,
                264,
                265,
                197,
                4,
                197,
                265,
                266,
                198,
                4,
                198,
                266,
                34,
                14,
                4,
                182,
                181,
                267,
                263,
                4,
                263,
                267,
                268,
                264,
                4,
                264,
                268,
                269,
                265,
                4,
                265,
                269,
                270,
                266,
                4,
                266,
                270,
                33,
                34,
                4,
                181,
                180,
                271,
                267,
                4,
                267,
                271,
                272,
                268,
                4,
                268,
                272,
                273,
                269,
                4,
                269,
                273,
                274,
                270,
                4,
                270,
                274,
                32,
                33,
                4,
                180,
                179,
                275,
                271,
                4,
                271,
                275,
                276,
                272,
                4,
                272,
                276,
                277,
                273,
                4,
                273,
                277,
                278,
                274,
                4,
                274,
                278,
                31,
                32,
                4,
                179,
                164,
                190,
                275,
                4,
                275,
                190,
                189,
                276,
                4,
                276,
                189,
                188,
                277,
                4,
                277,
                188,
                187,
                278,
                4,
                278,
                187,
                12,
                31,
                250,
                4,
                39,
                83,
                35,
                15,
                4,
                40,
                84,
                83,
                39,
                4,
                41,
                85,
                84,
                40,
                4,
                42,
                86,
                85,
                41,
                4,
                11,
                19,
                86,
                42,
                4,
                83,
                87,
                36,
                35,
                4,
                84,
                88,
                87,
                83,
                4,
                85,
                89,
                88,
                84,
                4,
                86,
                90,
                89,
                85,
                4,
                19,
                20,
                90,
                86,
                4,
                87,
                91,
                37,
                36,
                4,
                88,
                92,
                91,
                87,
                4,
                89,
                93,
                92,
                88,
                4,
                90,
                94,
                93,
                89,
                4,
                20,
                21,
                94,
                90,
                4,
                91,
                95,
                38,
                37,
                4,
                92,
                96,
                95,
                91,
                4,
                93,
                97,
                96,
                92,
                4,
                94,
                98,
                97,
                93,
                4,
                21,
                22,
                98,
                94,
                4,
                95,
                43,
                16,
                38,
                4,
                96,
                44,
                43,
                95,
                4,
                97,
                45,
                44,
                96,
                4,
                98,
                46,
                45,
                97,
                4,
                22,
                12,
                46,
                98,
                4,
                183,
                215,
                19,
                11,
                4,
                184,
                216,
                215,
                183,
                4,
                185,
                217,
                216,
                184,
                4,
                186,
                218,
                217,
                185,
                4,
                163,
                167,
                218,
                186,
                4,
                215,
                219,
                20,
                19,
                4,
                216,
                220,
                219,
                215,
                4,
                217,
                221,
                220,
                216,
                4,
                218,
                222,
                221,
                217,
                4,
                167,
                168,
                222,
                218,
                4,
                219,
                223,
                21,
                20,
                4,
                220,
                224,
                223,
                219,
                4,
                221,
                225,
                224,
                220,
                4,
                222,
                226,
                225,
                221,
                4,
                168,
                169,
                226,
                222,
                4,
                223,
                227,
                22,
                21,
                4,
                224,
                228,
                227,
                223,
                4,
                225,
                229,
                228,
                224,
                4,
                226,
                230,
                229,
                225,
                4,
                169,
                170,
                230,
                226,
                4,
                227,
                187,
                12,
                22,
                4,
                228,
                188,
                187,
                227,
                4,
                229,
                189,
                188,
                228,
                4,
                230,
                190,
                189,
                229,
                4,
                170,
                164,
                190,
                230,
                4,
                279,
                281,
                314,
                285,
                4,
                285,
                314,
                315,
                286,
                4,
                286,
                315,
                316,
                287,
                4,
                287,
                316,
                317,
                288,
                4,
                288,
                317,
                167,
                163,
                4,
                281,
                282,
                318,
                314,
                4,
                314,
                318,
                319,
                315,
                4,
                315,
                319,
                320,
                316,
                4,
                316,
                320,
                321,
                317,
                4,
                317,
                321,
                168,
                167,
                4,
                282,
                283,
                322,
                318,
                4,
                318,
                322,
                323,
                319,
                4,
                319,
                323,
                324,
                320,
                4,
                320,
                324,
                325,
                321,
                4,
                321,
                325,
                169,
                168,
                4,
                283,
                284,
                326,
                322,
                4,
                322,
                326,
                327,
                323,
                4,
                323,
                327,
                328,
                324,
                4,
                324,
                328,
                329,
                325,
                4,
                325,
                329,
                170,
                169,
                4,
                284,
                0,
                289,
                326,
                4,
                326,
                289,
                290,
                327,
                4,
                327,
                290,
                291,
                328,
                4,
                328,
                291,
                292,
                329,
                4,
                329,
                292,
                164,
                170,
                4,
                280,
                293,
                330,
                297,
                4,
                297,
                330,
                331,
                298,
                4,
                298,
                331,
                332,
                299,
                4,
                299,
                332,
                333,
                300,
                4,
                300,
                333,
                334,
                301,
                4,
                301,
                334,
                335,
                302,
                4,
                302,
                335,
                336,
                303,
                4,
                303,
                336,
                337,
                304,
                4,
                304,
                337,
                338,
                305,
                4,
                305,
                338,
                281,
                279,
                4,
                293,
                294,
                339,
                330,
                4,
                330,
                339,
                340,
                331,
                4,
                331,
                340,
                341,
                332,
                4,
                332,
                341,
                342,
                333,
                4,
                333,
                342,
                343,
                334,
                4,
                334,
                343,
                344,
                335,
                4,
                335,
                344,
                345,
                336,
                4,
                336,
                345,
                346,
                337,
                4,
                337,
                346,
                347,
                338,
                4,
                338,
                347,
                282,
                281,
                4,
                294,
                295,
                348,
                339,
                4,
                339,
                348,
                349,
                340,
                4,
                340,
                349,
                350,
                341,
                4,
                341,
                350,
                351,
                342,
                4,
                342,
                351,
                352,
                343,
                4,
                343,
                352,
                353,
                344,
                4,
                344,
                353,
                354,
                345,
                4,
                345,
                354,
                355,
                346,
                4,
                346,
                355,
                356,
                347,
                4,
                347,
                356,
                283,
                282,
                4,
                295,
                296,
                357,
                348,
                4,
                348,
                357,
                358,
                349,
                4,
                349,
                358,
                359,
                350,
                4,
                350,
                359,
                360,
                351,
                4,
                351,
                360,
                361,
                352,
                4,
                352,
                361,
                362,
                353,
                4,
                353,
                362,
                363,
                354,
                4,
                354,
                363,
                364,
                355,
                4,
                355,
                364,
                365,
                356,
                4,
                356,
                365,
                284,
                283,
                4,
                296,
                1,
                2,
                357,
                4,
                357,
                2,
                3,
                358,
                4,
                358,
                3,
                4,
                359,
                4,
                359,
                4,
                5,
                360,
                4,
                360,
                5,
                6,
                361,
                4,
                361,
                6,
                7,
                362,
                4,
                362,
                7,
                8,
                363,
                4,
                363,
                8,
                9,
                364,
                4,
                364,
                9,
                10,
                365,
                4,
                365,
                10,
                0,
                284,
                4,
                15,
                35,
                366,
                306,
                4,
                306,
                366,
                367,
                307,
                4,
                307,
                367,
                368,
                308,
                4,
                308,
                368,
                369,
                309,
                4,
                309,
                369,
                293,
                280,
                4,
                35,
                36,
                370,
                366,
                4,
                366,
                370,
                371,
                367,
                4,
                367,
                371,
                372,
                368,
                4,
                368,
                372,
                373,
                369,
                4,
                369,
                373,
                294,
                293,
                4,
                36,
                37,
                374,
                370,
                4,
                370,
                374,
                375,
                371,
                4,
                371,
                375,
                376,
                372,
                4,
                372,
                376,
                377,
                373,
                4,
                373,
                377,
                295,
                294,
                4,
                37,
                38,
                378,
                374,
                4,
                374,
                378,
                379,
                375,
                4,
                375,
                379,
                380,
                376,
                4,
                376,
                380,
                381,
                377,
                4,
                377,
                381,
                296,
                295,
                4,
                38,
                16,
                310,
                378,
                4,
                378,
                310,
                311,
                379,
                4,
                379,
                311,
                312,
                380,
                4,
                380,
                312,
                313,
                381,
                4,
                381,
                313,
                1,
                296,
                4,
                15,
                306,
                382,
                39,
                4,
                39,
                382,
                383,
                40,
                4,
                40,
                383,
                384,
                41,
                4,
                41,
                384,
                385,
                42,
                4,
                42,
                385,
                386,
                11,
                4,
                11,
                386,
                387,
                183,
                4,
                183,
                387,
                388,
                184,
                4,
                184,
                388,
                389,
                185,
                4,
                185,
                389,
                390,
                186,
                4,
                186,
                390,
                288,
                163,
                4,
                306,
                307,
                391,
                382,
                4,
                382,
                391,
                392,
                383,
                4,
                383,
                392,
                393,
                384,
                4,
                384,
                393,
                394,
                385,
                4,
                385,
                394,
                395,
                386,
                4,
                386,
                395,
                396,
                387,
                4,
                387,
                396,
                397,
                388,
                4,
                388,
                397,
                398,
                389,
                4,
                389,
                398,
                399,
                390,
                4,
                390,
                399,
                287,
                288,
                4,
                307,
                308,
                400,
                391,
                4,
                391,
                400,
                401,
                392,
                4,
                392,
                401,
                402,
                393,
                4,
                393,
                402,
                403,
                394,
                4,
                394,
                403,
                404,
                395,
                4,
                395,
                404,
                405,
                396,
                4,
                396,
                405,
                406,
                397,
                4,
                397,
                406,
                407,
                398,
                4,
                398,
                407,
                408,
                399,
                4,
                399,
                408,
                286,
                287,
                4,
                308,
                309,
                409,
                400,
                4,
                400,
                409,
                410,
                401,
                4,
                401,
                410,
                411,
                402,
                4,
                402,
                411,
                412,
                403,
                4,
                403,
                412,
                413,
                404,
                4,
                404,
                413,
                414,
                405,
                4,
                405,
                414,
                415,
                406,
                4,
                406,
                415,
                416,
                407,
                4,
                407,
                416,
                417,
                408,
                4,
                408,
                417,
                285,
                286,
                4,
                309,
                280,
                297,
                409,
                4,
                409,
                297,
                298,
                410,
                4,
                410,
                298,
                299,
                411,
                4,
                411,
                299,
                300,
                412,
                4,
                412,
                300,
                301,
                413,
                4,
                413,
                301,
                302,
                414,
                4,
                414,
                302,
                303,
                415,
                4,
                415,
                303,
                304,
                416,
                4,
                416,
                304,
                305,
                417,
                4,
                417,
                305,
                279,
                285,
                4,
                16,
                43,
                418,
                310,
                4,
                310,
                418,
                419,
                311,
                4,
                311,
                419,
                420,
                312,
                4,
                312,
                420,
                421,
                313,
                4,
                313,
                421,
                2,
                1,
                4,
                43,
                44,
                422,
                418,
                4,
                418,
                422,
                423,
                419,
                4,
                419,
                423,
                424,
                420,
                4,
                420,
                424,
                425,
                421,
                4,
                421,
                425,
                3,
                2,
                4,
                44,
                45,
                426,
                422,
                4,
                422,
                426,
                427,
                423,
                4,
                423,
                427,
                428,
                424,
                4,
                424,
                428,
                429,
                425,
                4,
                425,
                429,
                4,
                3,
                4,
                45,
                46,
                430,
                426,
                4,
                426,
                430,
                431,
                427,
                4,
                427,
                431,
                432,
                428,
                4,
                428,
                432,
                433,
                429,
                4,
                429,
                433,
                5,
                4,
                4,
                46,
                12,
                434,
                430,
                4,
                430,
                434,
                435,
                431,
                4,
                431,
                435,
                436,
                432,
                4,
                432,
                436,
                437,
                433,
                4,
                433,
                437,
                6,
                5,
                4,
                12,
                187,
                438,
                434,
                4,
                434,
                438,
                439,
                435,
                4,
                435,
                439,
                440,
                436,
                4,
                436,
                440,
                441,
                437,
                4,
                437,
                441,
                7,
                6,
                4,
                187,
                188,
                442,
                438,
                4,
                438,
                442,
                443,
                439,
                4,
                439,
                443,
                444,
                440,
                4,
                440,
                444,
                445,
                441,
                4,
                441,
                445,
                8,
                7,
                4,
                188,
                189,
                446,
                442,
                4,
                442,
                446,
                447,
                443,
                4,
                443,
                447,
                448,
                444,
                4,
                444,
                448,
                449,
                445,
                4,
                445,
                449,
                9,
                8,
                4,
                189,
                190,
                450,
                446,
                4,
                446,
                450,
                451,
                447,
                4,
                447,
                451,
                452,
                448,
                4,
                448,
                452,
                453,
                449,
                4,
                449,
                453,
                10,
                9,
                4,
                190,
                164,
                292,
                450,
                4,
                450,
                292,
                291,
                451,
                4,
                451,
                291,
                290,
                452,
                4,
                452,
                290,
                289,
                453,
                4,
                453,
                289,
                0,
                10,
            ]
        )
        facesOut, facesOutIndex = (
            DataArrayInt64.FromVTK93To94FacesInternaReprOfPolyedra(
                faces93, faceLocations93
            )
        )
        facesOutExp = DataArrayInt64(
            [
                11,
                19,
                67,
                27,
                27,
                67,
                68,
                28,
                28,
                68,
                69,
                29,
                29,
                69,
                70,
                30,
                30,
                70,
                23,
                13,
                19,
                20,
                71,
                67,
                67,
                71,
                72,
                68,
                68,
                72,
                73,
                69,
                69,
                73,
                74,
                70,
                70,
                74,
                24,
                23,
                20,
                21,
                75,
                71,
                71,
                75,
                76,
                72,
                72,
                76,
                77,
                73,
                73,
                77,
                78,
                74,
                74,
                78,
                25,
                24,
                21,
                22,
                79,
                75,
                75,
                79,
                80,
                76,
                76,
                80,
                81,
                77,
                77,
                81,
                82,
                78,
                78,
                82,
                26,
                25,
                22,
                12,
                31,
                79,
                79,
                31,
                32,
                80,
                80,
                32,
                33,
                81,
                81,
                33,
                34,
                82,
                82,
                34,
                14,
                26,
                15,
                35,
                83,
                39,
                39,
                83,
                84,
                40,
                40,
                84,
                85,
                41,
                41,
                85,
                86,
                42,
                42,
                86,
                19,
                11,
                35,
                36,
                87,
                83,
                83,
                87,
                88,
                84,
                84,
                88,
                89,
                85,
                85,
                89,
                90,
                86,
                86,
                90,
                20,
                19,
                36,
                37,
                91,
                87,
                87,
                91,
                92,
                88,
                88,
                92,
                93,
                89,
                89,
                93,
                94,
                90,
                90,
                94,
                21,
                20,
                37,
                38,
                95,
                91,
                91,
                95,
                96,
                92,
                92,
                96,
                97,
                93,
                93,
                97,
                98,
                94,
                94,
                98,
                22,
                21,
                38,
                16,
                43,
                95,
                95,
                43,
                44,
                96,
                96,
                44,
                45,
                97,
                97,
                45,
                46,
                98,
                98,
                46,
                12,
                22,
                17,
                47,
                99,
                51,
                51,
                99,
                100,
                52,
                52,
                100,
                101,
                53,
                53,
                101,
                102,
                54,
                54,
                102,
                35,
                15,
                47,
                48,
                103,
                99,
                99,
                103,
                104,
                100,
                100,
                104,
                105,
                101,
                101,
                105,
                106,
                102,
                102,
                106,
                36,
                35,
                48,
                49,
                107,
                103,
                103,
                107,
                108,
                104,
                104,
                108,
                109,
                105,
                105,
                109,
                110,
                106,
                106,
                110,
                37,
                36,
                49,
                50,
                111,
                107,
                107,
                111,
                112,
                108,
                108,
                112,
                113,
                109,
                109,
                113,
                114,
                110,
                110,
                114,
                38,
                37,
                50,
                18,
                55,
                111,
                111,
                55,
                56,
                112,
                112,
                56,
                57,
                113,
                113,
                57,
                58,
                114,
                114,
                58,
                16,
                38,
                13,
                23,
                115,
                59,
                59,
                115,
                116,
                60,
                60,
                116,
                117,
                61,
                61,
                117,
                118,
                62,
                62,
                118,
                47,
                17,
                23,
                24,
                119,
                115,
                115,
                119,
                120,
                116,
                116,
                120,
                121,
                117,
                117,
                121,
                122,
                118,
                118,
                122,
                48,
                47,
                24,
                25,
                123,
                119,
                119,
                123,
                124,
                120,
                120,
                124,
                125,
                121,
                121,
                125,
                126,
                122,
                122,
                126,
                49,
                48,
                25,
                26,
                127,
                123,
                123,
                127,
                128,
                124,
                124,
                128,
                129,
                125,
                125,
                129,
                130,
                126,
                126,
                130,
                50,
                49,
                26,
                14,
                63,
                127,
                127,
                63,
                64,
                128,
                128,
                64,
                65,
                129,
                129,
                65,
                66,
                130,
                130,
                66,
                18,
                50,
                11,
                27,
                131,
                42,
                42,
                131,
                132,
                41,
                41,
                132,
                133,
                40,
                40,
                133,
                134,
                39,
                39,
                134,
                54,
                15,
                27,
                28,
                135,
                131,
                131,
                135,
                136,
                132,
                132,
                136,
                137,
                133,
                133,
                137,
                138,
                134,
                134,
                138,
                53,
                54,
                28,
                29,
                139,
                135,
                135,
                139,
                140,
                136,
                136,
                140,
                141,
                137,
                137,
                141,
                142,
                138,
                138,
                142,
                52,
                53,
                29,
                30,
                143,
                139,
                139,
                143,
                144,
                140,
                140,
                144,
                145,
                141,
                141,
                145,
                146,
                142,
                142,
                146,
                51,
                52,
                30,
                13,
                59,
                143,
                143,
                59,
                60,
                144,
                144,
                60,
                61,
                145,
                145,
                61,
                62,
                146,
                146,
                62,
                17,
                51,
                14,
                34,
                147,
                63,
                63,
                147,
                148,
                64,
                64,
                148,
                149,
                65,
                65,
                149,
                150,
                66,
                66,
                150,
                55,
                18,
                34,
                33,
                151,
                147,
                147,
                151,
                152,
                148,
                148,
                152,
                153,
                149,
                149,
                153,
                154,
                150,
                150,
                154,
                56,
                55,
                33,
                32,
                155,
                151,
                151,
                155,
                156,
                152,
                152,
                156,
                157,
                153,
                153,
                157,
                158,
                154,
                154,
                158,
                57,
                56,
                32,
                31,
                159,
                155,
                155,
                159,
                160,
                156,
                156,
                160,
                161,
                157,
                157,
                161,
                162,
                158,
                158,
                162,
                58,
                57,
                31,
                12,
                46,
                159,
                159,
                46,
                45,
                160,
                160,
                45,
                44,
                161,
                161,
                44,
                43,
                162,
                162,
                43,
                16,
                58,
                163,
                167,
                199,
                175,
                175,
                199,
                200,
                176,
                176,
                200,
                201,
                177,
                177,
                201,
                202,
                178,
                178,
                202,
                171,
                165,
                167,
                168,
                203,
                199,
                199,
                203,
                204,
                200,
                200,
                204,
                205,
                201,
                201,
                205,
                206,
                202,
                202,
                206,
                172,
                171,
                168,
                169,
                207,
                203,
                203,
                207,
                208,
                204,
                204,
                208,
                209,
                205,
                205,
                209,
                210,
                206,
                206,
                210,
                173,
                172,
                169,
                170,
                211,
                207,
                207,
                211,
                212,
                208,
                208,
                212,
                213,
                209,
                209,
                213,
                214,
                210,
                210,
                214,
                174,
                173,
                170,
                164,
                179,
                211,
                211,
                179,
                180,
                212,
                212,
                180,
                181,
                213,
                213,
                181,
                182,
                214,
                214,
                182,
                166,
                174,
                11,
                19,
                215,
                183,
                183,
                215,
                216,
                184,
                184,
                216,
                217,
                185,
                185,
                217,
                218,
                186,
                186,
                218,
                167,
                163,
                19,
                20,
                219,
                215,
                215,
                219,
                220,
                216,
                216,
                220,
                221,
                217,
                217,
                221,
                222,
                218,
                218,
                222,
                168,
                167,
                20,
                21,
                223,
                219,
                219,
                223,
                224,
                220,
                220,
                224,
                225,
                221,
                221,
                225,
                226,
                222,
                222,
                226,
                169,
                168,
                21,
                22,
                227,
                223,
                223,
                227,
                228,
                224,
                224,
                228,
                229,
                225,
                225,
                229,
                230,
                226,
                226,
                230,
                170,
                169,
                22,
                12,
                187,
                227,
                227,
                187,
                188,
                228,
                228,
                188,
                189,
                229,
                229,
                189,
                190,
                230,
                230,
                190,
                164,
                170,
                27,
                67,
                19,
                11,
                28,
                68,
                67,
                27,
                29,
                69,
                68,
                28,
                30,
                70,
                69,
                29,
                13,
                23,
                70,
                30,
                67,
                71,
                20,
                19,
                68,
                72,
                71,
                67,
                69,
                73,
                72,
                68,
                70,
                74,
                73,
                69,
                23,
                24,
                74,
                70,
                71,
                75,
                21,
                20,
                72,
                76,
                75,
                71,
                73,
                77,
                76,
                72,
                74,
                78,
                77,
                73,
                24,
                25,
                78,
                74,
                75,
                79,
                22,
                21,
                76,
                80,
                79,
                75,
                77,
                81,
                80,
                76,
                78,
                82,
                81,
                77,
                25,
                26,
                82,
                78,
                79,
                31,
                12,
                22,
                80,
                32,
                31,
                79,
                81,
                33,
                32,
                80,
                82,
                34,
                33,
                81,
                26,
                14,
                34,
                82,
                165,
                171,
                231,
                191,
                191,
                231,
                232,
                192,
                192,
                232,
                233,
                193,
                193,
                233,
                234,
                194,
                194,
                234,
                23,
                13,
                171,
                172,
                235,
                231,
                231,
                235,
                236,
                232,
                232,
                236,
                237,
                233,
                233,
                237,
                238,
                234,
                234,
                238,
                24,
                23,
                172,
                173,
                239,
                235,
                235,
                239,
                240,
                236,
                236,
                240,
                241,
                237,
                237,
                241,
                242,
                238,
                238,
                242,
                25,
                24,
                173,
                174,
                243,
                239,
                239,
                243,
                244,
                240,
                240,
                244,
                245,
                241,
                241,
                245,
                246,
                242,
                242,
                246,
                26,
                25,
                174,
                166,
                195,
                243,
                243,
                195,
                196,
                244,
                244,
                196,
                197,
                245,
                245,
                197,
                198,
                246,
                246,
                198,
                14,
                26,
                163,
                175,
                247,
                186,
                186,
                247,
                248,
                185,
                185,
                248,
                249,
                184,
                184,
                249,
                250,
                183,
                183,
                250,
                27,
                11,
                175,
                176,
                251,
                247,
                247,
                251,
                252,
                248,
                248,
                252,
                253,
                249,
                249,
                253,
                254,
                250,
                250,
                254,
                28,
                27,
                176,
                177,
                255,
                251,
                251,
                255,
                256,
                252,
                252,
                256,
                257,
                253,
                253,
                257,
                258,
                254,
                254,
                258,
                29,
                28,
                177,
                178,
                259,
                255,
                255,
                259,
                260,
                256,
                256,
                260,
                261,
                257,
                257,
                261,
                262,
                258,
                258,
                262,
                30,
                29,
                178,
                165,
                191,
                259,
                259,
                191,
                192,
                260,
                260,
                192,
                193,
                261,
                261,
                193,
                194,
                262,
                262,
                194,
                13,
                30,
                166,
                182,
                263,
                195,
                195,
                263,
                264,
                196,
                196,
                264,
                265,
                197,
                197,
                265,
                266,
                198,
                198,
                266,
                34,
                14,
                182,
                181,
                267,
                263,
                263,
                267,
                268,
                264,
                264,
                268,
                269,
                265,
                265,
                269,
                270,
                266,
                266,
                270,
                33,
                34,
                181,
                180,
                271,
                267,
                267,
                271,
                272,
                268,
                268,
                272,
                273,
                269,
                269,
                273,
                274,
                270,
                270,
                274,
                32,
                33,
                180,
                179,
                275,
                271,
                271,
                275,
                276,
                272,
                272,
                276,
                277,
                273,
                273,
                277,
                278,
                274,
                274,
                278,
                31,
                32,
                179,
                164,
                190,
                275,
                275,
                190,
                189,
                276,
                276,
                189,
                188,
                277,
                277,
                188,
                187,
                278,
                278,
                187,
                12,
                31,
                39,
                83,
                35,
                15,
                40,
                84,
                83,
                39,
                41,
                85,
                84,
                40,
                42,
                86,
                85,
                41,
                11,
                19,
                86,
                42,
                83,
                87,
                36,
                35,
                84,
                88,
                87,
                83,
                85,
                89,
                88,
                84,
                86,
                90,
                89,
                85,
                19,
                20,
                90,
                86,
                87,
                91,
                37,
                36,
                88,
                92,
                91,
                87,
                89,
                93,
                92,
                88,
                90,
                94,
                93,
                89,
                20,
                21,
                94,
                90,
                91,
                95,
                38,
                37,
                92,
                96,
                95,
                91,
                93,
                97,
                96,
                92,
                94,
                98,
                97,
                93,
                21,
                22,
                98,
                94,
                95,
                43,
                16,
                38,
                96,
                44,
                43,
                95,
                97,
                45,
                44,
                96,
                98,
                46,
                45,
                97,
                22,
                12,
                46,
                98,
                183,
                215,
                19,
                11,
                184,
                216,
                215,
                183,
                185,
                217,
                216,
                184,
                186,
                218,
                217,
                185,
                163,
                167,
                218,
                186,
                215,
                219,
                20,
                19,
                216,
                220,
                219,
                215,
                217,
                221,
                220,
                216,
                218,
                222,
                221,
                217,
                167,
                168,
                222,
                218,
                219,
                223,
                21,
                20,
                220,
                224,
                223,
                219,
                221,
                225,
                224,
                220,
                222,
                226,
                225,
                221,
                168,
                169,
                226,
                222,
                223,
                227,
                22,
                21,
                224,
                228,
                227,
                223,
                225,
                229,
                228,
                224,
                226,
                230,
                229,
                225,
                169,
                170,
                230,
                226,
                227,
                187,
                12,
                22,
                228,
                188,
                187,
                227,
                229,
                189,
                188,
                228,
                230,
                190,
                189,
                229,
                170,
                164,
                190,
                230,
                279,
                281,
                314,
                285,
                285,
                314,
                315,
                286,
                286,
                315,
                316,
                287,
                287,
                316,
                317,
                288,
                288,
                317,
                167,
                163,
                281,
                282,
                318,
                314,
                314,
                318,
                319,
                315,
                315,
                319,
                320,
                316,
                316,
                320,
                321,
                317,
                317,
                321,
                168,
                167,
                282,
                283,
                322,
                318,
                318,
                322,
                323,
                319,
                319,
                323,
                324,
                320,
                320,
                324,
                325,
                321,
                321,
                325,
                169,
                168,
                283,
                284,
                326,
                322,
                322,
                326,
                327,
                323,
                323,
                327,
                328,
                324,
                324,
                328,
                329,
                325,
                325,
                329,
                170,
                169,
                284,
                0,
                289,
                326,
                326,
                289,
                290,
                327,
                327,
                290,
                291,
                328,
                328,
                291,
                292,
                329,
                329,
                292,
                164,
                170,
                280,
                293,
                330,
                297,
                297,
                330,
                331,
                298,
                298,
                331,
                332,
                299,
                299,
                332,
                333,
                300,
                300,
                333,
                334,
                301,
                301,
                334,
                335,
                302,
                302,
                335,
                336,
                303,
                303,
                336,
                337,
                304,
                304,
                337,
                338,
                305,
                305,
                338,
                281,
                279,
                293,
                294,
                339,
                330,
                330,
                339,
                340,
                331,
                331,
                340,
                341,
                332,
                332,
                341,
                342,
                333,
                333,
                342,
                343,
                334,
                334,
                343,
                344,
                335,
                335,
                344,
                345,
                336,
                336,
                345,
                346,
                337,
                337,
                346,
                347,
                338,
                338,
                347,
                282,
                281,
                294,
                295,
                348,
                339,
                339,
                348,
                349,
                340,
                340,
                349,
                350,
                341,
                341,
                350,
                351,
                342,
                342,
                351,
                352,
                343,
                343,
                352,
                353,
                344,
                344,
                353,
                354,
                345,
                345,
                354,
                355,
                346,
                346,
                355,
                356,
                347,
                347,
                356,
                283,
                282,
                295,
                296,
                357,
                348,
                348,
                357,
                358,
                349,
                349,
                358,
                359,
                350,
                350,
                359,
                360,
                351,
                351,
                360,
                361,
                352,
                352,
                361,
                362,
                353,
                353,
                362,
                363,
                354,
                354,
                363,
                364,
                355,
                355,
                364,
                365,
                356,
                356,
                365,
                284,
                283,
                296,
                1,
                2,
                357,
                357,
                2,
                3,
                358,
                358,
                3,
                4,
                359,
                359,
                4,
                5,
                360,
                360,
                5,
                6,
                361,
                361,
                6,
                7,
                362,
                362,
                7,
                8,
                363,
                363,
                8,
                9,
                364,
                364,
                9,
                10,
                365,
                365,
                10,
                0,
                284,
                15,
                35,
                366,
                306,
                306,
                366,
                367,
                307,
                307,
                367,
                368,
                308,
                308,
                368,
                369,
                309,
                309,
                369,
                293,
                280,
                35,
                36,
                370,
                366,
                366,
                370,
                371,
                367,
                367,
                371,
                372,
                368,
                368,
                372,
                373,
                369,
                369,
                373,
                294,
                293,
                36,
                37,
                374,
                370,
                370,
                374,
                375,
                371,
                371,
                375,
                376,
                372,
                372,
                376,
                377,
                373,
                373,
                377,
                295,
                294,
                37,
                38,
                378,
                374,
                374,
                378,
                379,
                375,
                375,
                379,
                380,
                376,
                376,
                380,
                381,
                377,
                377,
                381,
                296,
                295,
                38,
                16,
                310,
                378,
                378,
                310,
                311,
                379,
                379,
                311,
                312,
                380,
                380,
                312,
                313,
                381,
                381,
                313,
                1,
                296,
                15,
                306,
                382,
                39,
                39,
                382,
                383,
                40,
                40,
                383,
                384,
                41,
                41,
                384,
                385,
                42,
                42,
                385,
                386,
                11,
                11,
                386,
                387,
                183,
                183,
                387,
                388,
                184,
                184,
                388,
                389,
                185,
                185,
                389,
                390,
                186,
                186,
                390,
                288,
                163,
                306,
                307,
                391,
                382,
                382,
                391,
                392,
                383,
                383,
                392,
                393,
                384,
                384,
                393,
                394,
                385,
                385,
                394,
                395,
                386,
                386,
                395,
                396,
                387,
                387,
                396,
                397,
                388,
                388,
                397,
                398,
                389,
                389,
                398,
                399,
                390,
                390,
                399,
                287,
                288,
                307,
                308,
                400,
                391,
                391,
                400,
                401,
                392,
                392,
                401,
                402,
                393,
                393,
                402,
                403,
                394,
                394,
                403,
                404,
                395,
                395,
                404,
                405,
                396,
                396,
                405,
                406,
                397,
                397,
                406,
                407,
                398,
                398,
                407,
                408,
                399,
                399,
                408,
                286,
                287,
                308,
                309,
                409,
                400,
                400,
                409,
                410,
                401,
                401,
                410,
                411,
                402,
                402,
                411,
                412,
                403,
                403,
                412,
                413,
                404,
                404,
                413,
                414,
                405,
                405,
                414,
                415,
                406,
                406,
                415,
                416,
                407,
                407,
                416,
                417,
                408,
                408,
                417,
                285,
                286,
                309,
                280,
                297,
                409,
                409,
                297,
                298,
                410,
                410,
                298,
                299,
                411,
                411,
                299,
                300,
                412,
                412,
                300,
                301,
                413,
                413,
                301,
                302,
                414,
                414,
                302,
                303,
                415,
                415,
                303,
                304,
                416,
                416,
                304,
                305,
                417,
                417,
                305,
                279,
                285,
                16,
                43,
                418,
                310,
                310,
                418,
                419,
                311,
                311,
                419,
                420,
                312,
                312,
                420,
                421,
                313,
                313,
                421,
                2,
                1,
                43,
                44,
                422,
                418,
                418,
                422,
                423,
                419,
                419,
                423,
                424,
                420,
                420,
                424,
                425,
                421,
                421,
                425,
                3,
                2,
                44,
                45,
                426,
                422,
                422,
                426,
                427,
                423,
                423,
                427,
                428,
                424,
                424,
                428,
                429,
                425,
                425,
                429,
                4,
                3,
                45,
                46,
                430,
                426,
                426,
                430,
                431,
                427,
                427,
                431,
                432,
                428,
                428,
                432,
                433,
                429,
                429,
                433,
                5,
                4,
                46,
                12,
                434,
                430,
                430,
                434,
                435,
                431,
                431,
                435,
                436,
                432,
                432,
                436,
                437,
                433,
                433,
                437,
                6,
                5,
                12,
                187,
                438,
                434,
                434,
                438,
                439,
                435,
                435,
                439,
                440,
                436,
                436,
                440,
                441,
                437,
                437,
                441,
                7,
                6,
                187,
                188,
                442,
                438,
                438,
                442,
                443,
                439,
                439,
                443,
                444,
                440,
                440,
                444,
                445,
                441,
                441,
                445,
                8,
                7,
                188,
                189,
                446,
                442,
                442,
                446,
                447,
                443,
                443,
                447,
                448,
                444,
                444,
                448,
                449,
                445,
                445,
                449,
                9,
                8,
                189,
                190,
                450,
                446,
                446,
                450,
                451,
                447,
                447,
                451,
                452,
                448,
                448,
                452,
                453,
                449,
                449,
                453,
                10,
                9,
                190,
                164,
                292,
                450,
                450,
                292,
                291,
                451,
                451,
                291,
                290,
                452,
                452,
                290,
                289,
                453,
                453,
                289,
                0,
                10,
            ]
        )
        facesOutIndexExp = DataArrayInt(550)
        facesOutIndexExp[:] = 4
        facesOutIndexExp.computeOffsetsFull()
        self.assertTrue(facesOut.isEqual(facesOutExp))
        self.assertTrue(facesOutIndex.isEqual(facesOutIndexExp))
        facesOffset = DataArrayInt([0, 150, 300, 550])
        aa, bb = DataArrayInt64.FromVTK94InternalReprOfPolyedra(
            facesOut, facesOutIndex, facesOffset
        )
        m = MEDCoupling1DGTUMesh("mesh", NORM_POLYHED)
        m.setCoords(DataArrayDouble(10000, 3))
        m.setNodalConnectivity(aa, bb)
        m2 = m.buildUnstructured()
        self.assertTrue(
            m2.computeNbOfFacesPerCell().isEqual(facesOffset.deltaShiftIndex())
        )
        self.assertTrue(aa[aa.findIdsNotEqual(-1)].isEqual(facesOut))

    def testDAWriteLoadForDbg0(self):
        """
        EDF31746: To help debugging in C++
        """
        import tempfile
        from pathlib import Path

        arr = DataArrayDouble(5)
        arr.iota()
        arr = DataArrayDouble.Meld([arr, 2 * arr, 3 * arr])
        arrFloat = arr.convertToFloatArr()
        arrI32 = arr.convertToIntArr()
        arrI64 = arr.convertToInt64Arr()
        with tempfile.TemporaryDirectory() as dir:
            p = Path(dir) / "toto.bin"
            arr.writeForDbg(f"{p}")
            arr2 = DataArrayDouble.LoadForDbg(f"{p}")
            self.assertTrue(arr.isEqual(arr2, 1e-15))
            self.assertRaises(InterpKernelException, DataArrayInt32.LoadForDbg, f"{p}")
            p = Path(dir) / "toto2.bin"
            arrI32.writeForDbg(f"{p}")
            arrI32_2 = DataArrayInt32.LoadForDbg(f"{p}")
            self.assertTrue(arrI32.isEqual(arrI32_2))
            self.assertRaises(InterpKernelException, DataArrayInt64.LoadForDbg, f"{p}")
            p = Path(dir) / "toto3.bin"
            arrI64.writeForDbg(f"{p}")
            arrI64_2 = DataArrayInt64.LoadForDbg(f"{p}")
            self.assertTrue(arrI64.isEqual(arrI64_2))
            self.assertRaises(InterpKernelException, DataArrayInt32.LoadForDbg, f"{p}")
            p = Path(dir) / "toto4.bin"
            arrFloat.writeForDbg(f"{p}")
            arrFloat_2 = DataArrayFloat.LoadForDbg(f"{p}")
            self.assertTrue(arrFloat.isEqual(arrFloat_2, 1e-8))
            self.assertRaises(InterpKernelException, DataArrayDouble.LoadForDbg, f"{p}")
            pass

    def testUMeshConvertToQuadraticBasedOnSeg3_0(self):
        """
        [EDF32060] : pilot nodes in middle position during linear to quadratic transformation
        """
        coo = DataArrayDouble(
            [
                (0, 0),
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 1),
                (1, 2),
                (1, 3),
                (2, 0),
                (2, 1),
                (2, 2),
                (2, 3),
                (3, 0),
                (3, 1),
                (3, 2),
            ]
        )
        m = MEDCouplingUMesh("mesh", 2)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_QUAD4, [0, 1, 4, 3])
        m.insertNextCell(NORM_TRI3, [3, 4, 7])
        m.insertNextCell(NORM_QUAD4, [4, 5, 9, 8])
        m.insertNextCell(NORM_QUAD4, [1, 2, 5, 4])
        m.insertNextCell(NORM_QUAD4, [5, 6, 10, 9])
        m.insertNextCell(NORM_TRI3, [4, 8, 7])
        m.insertNextCell(NORM_QUAD4, [7, 8, 12, 11])
        m.insertNextCell(NORM_QUAD4, [8, 9, 13, 12])
        m.changeSpaceDimension(3, 0.0)
        coords1D = DataArrayDouble([(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 0, 3)])
        m1D = MEDCouplingUMesh.Build1DMeshFromCoords(coords1D)
        m3D = m.buildExtrudedMesh(m1D, 0)
        #
        shufflingArr = DataArrayInt(
            [
                110,
                119,
                54,
                100,
                82,
                48,
                17,
                95,
                24,
                12,
                49,
                125,
                34,
                35,
                92,
                75,
                96,
                5,
                51,
                55,
                68,
                99,
                80,
                39,
                57,
                0,
                46,
                124,
                27,
                109,
                117,
                42,
                103,
                83,
                15,
                86,
                63,
                76,
                98,
                45,
                52,
                115,
                18,
                87,
                44,
                11,
                104,
                28,
                93,
                79,
                13,
                3,
                114,
                65,
                38,
                31,
                94,
                7,
                20,
                120,
                72,
                62,
                32,
                56,
                37,
                122,
                9,
                47,
                107,
                102,
                10,
                22,
                58,
                25,
                50,
                33,
                19,
                64,
                106,
                43,
                71,
                26,
                116,
                111,
                88,
                21,
                74,
                91,
                30,
                1,
                121,
                78,
                16,
                14,
                2,
                40,
                8,
                41,
                61,
                36,
                69,
                70,
                123,
                60,
                118,
                108,
                101,
                84,
                67,
                112,
                89,
                97,
                105,
                66,
                90,
                77,
                85,
                81,
                113,
                6,
                59,
                23,
                53,
                29,
                73,
                4,
            ]
        )
        edges, _, _, _, _ = m3D.explode3DMeshTo1D()
        edges.renumberCells(shufflingArr)
        coo_of_middles = edges.computeCellCenterOfMass()
        edges_1gt = MEDCoupling1GTUMesh.New(edges)
        conn = edges_1gt.getNodalConnectivity()
        conn_seg3 = DataArrayInt(edges_1gt.getNumberOfCells(), 3)
        conn.rearrange(2)
        conn_seg3[:, (0, 1)] = conn
        conn.rearrange(1)
        mid_conn = DataArrayInt(len(coo_of_middles))
        mid_conn.iota()
        mid_conn += m3D.getNumberOfNodes()
        #
        conn_seg3[:, 2] = mid_conn
        m1D_seg3 = MEDCoupling1SGTUMesh("seg3", NORM_SEG3)
        conn_seg3.rearrange(1)
        coo_with_middles = DataArrayDouble.Aggregate([m3D.getCoords(), coo_of_middles])
        m1D_seg3.setCoords(coo_with_middles)
        m1D_seg3.setNodalConnectivity(conn_seg3)
        ref = m1D_seg3.getCoords().getHiddenCppPointer()
        #
        zeRet = m3D.convertToQuadraticBasedOnSeg3(
            m1D_seg3
        )  # <- interesting call is here
        self.assertEqual(zeRet.getCoords().getHiddenCppPointer(), ref)
        self.assertEqual(zeRet.getAllGeoTypes(), [NORM_PENTA15, NORM_HEXA20])
        candidate = zeRet.explode3DMeshTo1D()[0]
        candidate.checkGeoEquivalWith(
            m1D_seg3.buildUnstructured(), 22, 0.0
        )  # <- important test is here
        zeRet.convertQuadraticCellsToLinear()
        zeRet.zipCoords()
        self.assertTrue(zeRet.isEqual(m3D, 1e-12))  # <- less but important test is here

    def testUMeshExtrudeConnectivity_quadratic0(self):
        """
        [EDF32060] : test quadratic case
        """
        coo = DataArrayDouble(
            [
                (0, 0),
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 1),
                (1, 2),
                (1, 3),
                (2, 0),
                (2, 1),
                (2, 2),
                (2, 3),
                (3, 0),
                (3, 1),
                (3, 2),
            ]
        )
        m = MEDCouplingUMesh("mesh", 2)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_QUAD4, [0, 1, 4, 3])
        m.insertNextCell(NORM_TRI3, [3, 4, 7])
        m.insertNextCell(NORM_QUAD4, [4, 5, 9, 8])
        m.insertNextCell(NORM_QUAD4, [1, 2, 5, 4])
        m.insertNextCell(NORM_QUAD4, [5, 6, 10, 9])
        m.insertNextCell(NORM_TRI3, [4, 8, 7])
        m.insertNextCell(NORM_QUAD4, [7, 8, 12, 11])
        m.insertNextCell(NORM_QUAD4, [8, 9, 13, 12])
        m.changeSpaceDimension(3, 0.0)
        m.convertLinearCellsToQuadratic()
        #
        nbExtrusions = 3
        ret = m.extrudeConnectivity(nbExtrusions)
        connExp = DataArrayInt(
            [
                30,
                0,
                1,
                4,
                3,
                70,
                71,
                74,
                73,
                14,
                15,
                16,
                17,
                84,
                85,
                86,
                87,
                35,
                36,
                39,
                38,
                25,
                3,
                4,
                7,
                73,
                74,
                77,
                16,
                18,
                19,
                86,
                88,
                89,
                38,
                39,
                42,
                30,
                4,
                5,
                9,
                8,
                74,
                75,
                79,
                78,
                20,
                21,
                22,
                23,
                90,
                91,
                92,
                93,
                39,
                40,
                44,
                43,
                30,
                1,
                2,
                5,
                4,
                71,
                72,
                75,
                74,
                24,
                25,
                20,
                15,
                94,
                95,
                90,
                85,
                36,
                37,
                40,
                39,
                30,
                5,
                6,
                10,
                9,
                75,
                76,
                80,
                79,
                26,
                27,
                28,
                21,
                96,
                97,
                98,
                91,
                40,
                41,
                45,
                44,
                25,
                4,
                8,
                7,
                74,
                78,
                77,
                23,
                29,
                18,
                93,
                99,
                88,
                39,
                43,
                42,
                30,
                7,
                8,
                12,
                11,
                77,
                78,
                82,
                81,
                29,
                30,
                31,
                32,
                99,
                100,
                101,
                102,
                42,
                43,
                47,
                46,
                30,
                8,
                9,
                13,
                12,
                78,
                79,
                83,
                82,
                22,
                33,
                34,
                30,
                92,
                103,
                104,
                100,
                43,
                44,
                48,
                47,
                30,
                70,
                71,
                74,
                73,
                140,
                141,
                144,
                143,
                84,
                85,
                86,
                87,
                154,
                155,
                156,
                157,
                105,
                106,
                109,
                108,
                25,
                73,
                74,
                77,
                143,
                144,
                147,
                86,
                88,
                89,
                156,
                158,
                159,
                108,
                109,
                112,
                30,
                74,
                75,
                79,
                78,
                144,
                145,
                149,
                148,
                90,
                91,
                92,
                93,
                160,
                161,
                162,
                163,
                109,
                110,
                114,
                113,
                30,
                71,
                72,
                75,
                74,
                141,
                142,
                145,
                144,
                94,
                95,
                90,
                85,
                164,
                165,
                160,
                155,
                106,
                107,
                110,
                109,
                30,
                75,
                76,
                80,
                79,
                145,
                146,
                150,
                149,
                96,
                97,
                98,
                91,
                166,
                167,
                168,
                161,
                110,
                111,
                115,
                114,
                25,
                74,
                78,
                77,
                144,
                148,
                147,
                93,
                99,
                88,
                163,
                169,
                158,
                109,
                113,
                112,
                30,
                77,
                78,
                82,
                81,
                147,
                148,
                152,
                151,
                99,
                100,
                101,
                102,
                169,
                170,
                171,
                172,
                112,
                113,
                117,
                116,
                30,
                78,
                79,
                83,
                82,
                148,
                149,
                153,
                152,
                92,
                103,
                104,
                100,
                162,
                173,
                174,
                170,
                113,
                114,
                118,
                117,
                30,
                140,
                141,
                144,
                143,
                210,
                211,
                214,
                213,
                154,
                155,
                156,
                157,
                224,
                225,
                226,
                227,
                175,
                176,
                179,
                178,
                25,
                143,
                144,
                147,
                213,
                214,
                217,
                156,
                158,
                159,
                226,
                228,
                229,
                178,
                179,
                182,
                30,
                144,
                145,
                149,
                148,
                214,
                215,
                219,
                218,
                160,
                161,
                162,
                163,
                230,
                231,
                232,
                233,
                179,
                180,
                184,
                183,
                30,
                141,
                142,
                145,
                144,
                211,
                212,
                215,
                214,
                164,
                165,
                160,
                155,
                234,
                235,
                230,
                225,
                176,
                177,
                180,
                179,
                30,
                145,
                146,
                150,
                149,
                215,
                216,
                220,
                219,
                166,
                167,
                168,
                161,
                236,
                237,
                238,
                231,
                180,
                181,
                185,
                184,
                25,
                144,
                148,
                147,
                214,
                218,
                217,
                163,
                169,
                158,
                233,
                239,
                228,
                179,
                183,
                182,
                30,
                147,
                148,
                152,
                151,
                217,
                218,
                222,
                221,
                169,
                170,
                171,
                172,
                239,
                240,
                241,
                242,
                182,
                183,
                187,
                186,
                30,
                148,
                149,
                153,
                152,
                218,
                219,
                223,
                222,
                162,
                173,
                174,
                170,
                232,
                243,
                244,
                240,
                183,
                184,
                188,
                187,
            ]
        )
        connIExp = DataArrayInt(
            [
                0,
                21,
                37,
                58,
                79,
                100,
                116,
                137,
                158,
                179,
                195,
                216,
                237,
                258,
                274,
                295,
                316,
                337,
                353,
                374,
                395,
                416,
                432,
                453,
                474,
            ]
        )
        ret.checkConsistency()
        self.assertTrue(ret.getNodalConnectivity().isEqual(connExp))
        self.assertTrue(ret.getNodalConnectivityIndex().isEqual(connIExp))
        self.assertEqual(
            ret.getNumberOfNodes(), (2 * nbExtrusions + 1) * m.getNumberOfNodes()
        )
        for i in range((2 * nbExtrusions + 1)):
            self.assertTrue(
                ret.getCoords()[
                    i * m.getNumberOfNodes() : (i + 1) * m.getNumberOfNodes()
                ].isEqual(m.getCoords(), 1e-12)
            )
        pass

    def testUMeshExtrudeConnectivity_linear0(self):
        """
        [EDF32060] : test linear case
        """
        coo = DataArrayDouble(
            [
                (0, 0),
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 1),
                (1, 2),
                (1, 3),
                (2, 0),
                (2, 1),
                (2, 2),
                (2, 3),
                (3, 0),
                (3, 1),
                (3, 2),
            ]
        )
        m = MEDCouplingUMesh("mesh", 2)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_QUAD4, [0, 1, 4, 3])
        m.insertNextCell(NORM_TRI3, [3, 4, 7])
        m.insertNextCell(NORM_QUAD4, [4, 5, 9, 8])
        m.insertNextCell(NORM_QUAD4, [1, 2, 5, 4])
        m.insertNextCell(NORM_QUAD4, [5, 6, 10, 9])
        m.insertNextCell(NORM_TRI3, [4, 8, 7])
        m.insertNextCell(NORM_QUAD4, [7, 8, 12, 11])
        m.insertNextCell(NORM_QUAD4, [8, 9, 13, 12])
        m.changeSpaceDimension(3, 0.0)
        connExp = DataArrayInt(
            [
                18,
                0,
                1,
                4,
                3,
                14,
                15,
                18,
                17,
                16,
                3,
                4,
                7,
                17,
                18,
                21,
                18,
                4,
                5,
                9,
                8,
                18,
                19,
                23,
                22,
                18,
                1,
                2,
                5,
                4,
                15,
                16,
                19,
                18,
                18,
                5,
                6,
                10,
                9,
                19,
                20,
                24,
                23,
                16,
                4,
                8,
                7,
                18,
                22,
                21,
                18,
                7,
                8,
                12,
                11,
                21,
                22,
                26,
                25,
                18,
                8,
                9,
                13,
                12,
                22,
                23,
                27,
                26,
                18,
                14,
                15,
                18,
                17,
                28,
                29,
                32,
                31,
                16,
                17,
                18,
                21,
                31,
                32,
                35,
                18,
                18,
                19,
                23,
                22,
                32,
                33,
                37,
                36,
                18,
                15,
                16,
                19,
                18,
                29,
                30,
                33,
                32,
                18,
                19,
                20,
                24,
                23,
                33,
                34,
                38,
                37,
                16,
                18,
                22,
                21,
                32,
                36,
                35,
                18,
                21,
                22,
                26,
                25,
                35,
                36,
                40,
                39,
                18,
                22,
                23,
                27,
                26,
                36,
                37,
                41,
                40,
                18,
                28,
                29,
                32,
                31,
                42,
                43,
                46,
                45,
                16,
                31,
                32,
                35,
                45,
                46,
                49,
                18,
                32,
                33,
                37,
                36,
                46,
                47,
                51,
                50,
                18,
                29,
                30,
                33,
                32,
                43,
                44,
                47,
                46,
                18,
                33,
                34,
                38,
                37,
                47,
                48,
                52,
                51,
                16,
                32,
                36,
                35,
                46,
                50,
                49,
                18,
                35,
                36,
                40,
                39,
                49,
                50,
                54,
                53,
                18,
                36,
                37,
                41,
                40,
                50,
                51,
                55,
                54,
            ]
        )
        connIExp = DataArrayInt(
            [
                0,
                9,
                16,
                25,
                34,
                43,
                50,
                59,
                68,
                77,
                84,
                93,
                102,
                111,
                118,
                127,
                136,
                145,
                152,
                161,
                170,
                179,
                186,
                195,
                204,
            ]
        )
        nbExtrusions = 3
        ret = m.extrudeConnectivity(nbExtrusions)
        ret.checkConsistency()
        self.assertEqual(
            ret.getNumberOfNodes(), (nbExtrusions + 1) * m.getNumberOfNodes()
        )
        self.assertTrue(ret.getNodalConnectivity().isEqual(connExp))
        self.assertTrue(ret.getNodalConnectivityIndex().isEqual(connIExp))
        for i in range((nbExtrusions + 1)):
            self.assertEqual(
                ret.getNumberOfNodes(), (nbExtrusions + 1) * m.getNumberOfNodes()
            )

    def testUMeshExtrudeConnectivity_biquadratic0(self):
        """
        [EDF32060] : test biquadratic case
        """
        coo = DataArrayDouble(
            [
                (0, 0),
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 1),
                (1, 2),
                (1, 3),
                (2, 0),
                (2, 1),
                (2, 2),
                (2, 3),
                (3, 0),
                (3, 1),
                (3, 2),
            ]
        )
        m = MEDCouplingUMesh("mesh", 2)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_QUAD4, [0, 1, 4, 3])
        m.insertNextCell(NORM_TRI3, [3, 4, 7])
        m.insertNextCell(NORM_QUAD4, [4, 5, 9, 8])
        m.insertNextCell(NORM_QUAD4, [1, 2, 5, 4])
        m.insertNextCell(NORM_QUAD4, [5, 6, 10, 9])
        m.insertNextCell(NORM_TRI3, [4, 8, 7])
        m.insertNextCell(NORM_QUAD4, [7, 8, 12, 11])
        m.insertNextCell(NORM_QUAD4, [8, 9, 13, 12])

        m.changeSpaceDimension(3, 0.0)
        # conversion quadratic to bi-quadratic
        mapTrad = {NORM_TRI6: NORM_TRI7, NORM_QUAD8: NORM_QUAD9}
        m.convertLinearCellsToQuadratic()
        offset = m.getNumberOfNodes()
        m.setCoords(
            DataArrayDouble.Aggregate([m.getCoords(), m.computeCellCenterOfMass()])
        )
        o2n = m.sortCellsInMEDFileFrmt()
        n2o = o2n.invertArrayO2N2N2O(len(o2n))
        types = m.splitByType()
        ret = []
        offsetCell = 0
        for type in types:
            mt = MEDCoupling1SGTUMesh(type)
            conn = mt.getNodalConnectivity()
            conn.rearrange(mt.getNumberOfNodesPerCell())
            connOut = DataArrayInt.Meld(
                conn, n2o[offsetCell : offsetCell + type.getNumberOfCells()] + offset
            )
            connOut.rearrange(1)
            mtOut = MEDCoupling1SGTUMesh("", mapTrad[mt.getCellModelEnum()])
            mtOut.setNodalConnectivity(connOut)
            mtOut.setCoords(m.getCoords())
            ret.append(mtOut.buildUnstructured())
            offsetCell += type.getNumberOfCells()
        m = MEDCouplingUMesh.MergeUMeshesOnSameCoords(ret)
        m.renumberCells(n2o)
        # conversion quadratic to bi-quadratic
        connExp = DataArrayInt(
            [
                27,
                0,
                1,
                4,
                3,
                86,
                87,
                90,
                89,
                14,
                15,
                16,
                17,
                100,
                101,
                102,
                103,
                43,
                44,
                47,
                46,
                35,
                57,
                58,
                59,
                60,
                121,
                78,
                28,
                3,
                4,
                7,
                89,
                90,
                93,
                16,
                18,
                19,
                102,
                104,
                105,
                46,
                47,
                50,
                59,
                61,
                62,
                27,
                4,
                5,
                9,
                8,
                90,
                91,
                95,
                94,
                20,
                21,
                22,
                23,
                106,
                107,
                108,
                109,
                47,
                48,
                52,
                51,
                37,
                63,
                64,
                65,
                66,
                123,
                80,
                27,
                1,
                2,
                5,
                4,
                87,
                88,
                91,
                90,
                24,
                25,
                20,
                15,
                110,
                111,
                106,
                101,
                44,
                45,
                48,
                47,
                38,
                67,
                68,
                63,
                58,
                124,
                81,
                27,
                5,
                6,
                10,
                9,
                91,
                92,
                96,
                95,
                26,
                27,
                28,
                21,
                112,
                113,
                114,
                107,
                48,
                49,
                53,
                52,
                39,
                69,
                70,
                71,
                64,
                125,
                82,
                28,
                4,
                8,
                7,
                90,
                94,
                93,
                23,
                29,
                18,
                109,
                115,
                104,
                47,
                51,
                50,
                66,
                72,
                61,
                27,
                7,
                8,
                12,
                11,
                93,
                94,
                98,
                97,
                29,
                30,
                31,
                32,
                115,
                116,
                117,
                118,
                50,
                51,
                55,
                54,
                41,
                72,
                73,
                74,
                75,
                127,
                84,
                27,
                8,
                9,
                13,
                12,
                94,
                95,
                99,
                98,
                22,
                33,
                34,
                30,
                108,
                119,
                120,
                116,
                51,
                52,
                56,
                55,
                42,
                65,
                76,
                77,
                73,
                128,
                85,
                27,
                86,
                87,
                90,
                89,
                172,
                173,
                176,
                175,
                100,
                101,
                102,
                103,
                186,
                187,
                188,
                189,
                129,
                130,
                133,
                132,
                121,
                143,
                144,
                145,
                146,
                207,
                164,
                28,
                89,
                90,
                93,
                175,
                176,
                179,
                102,
                104,
                105,
                188,
                190,
                191,
                132,
                133,
                136,
                145,
                147,
                148,
                27,
                90,
                91,
                95,
                94,
                176,
                177,
                181,
                180,
                106,
                107,
                108,
                109,
                192,
                193,
                194,
                195,
                133,
                134,
                138,
                137,
                123,
                149,
                150,
                151,
                152,
                209,
                166,
                27,
                87,
                88,
                91,
                90,
                173,
                174,
                177,
                176,
                110,
                111,
                106,
                101,
                196,
                197,
                192,
                187,
                130,
                131,
                134,
                133,
                124,
                153,
                154,
                149,
                144,
                210,
                167,
                27,
                91,
                92,
                96,
                95,
                177,
                178,
                182,
                181,
                112,
                113,
                114,
                107,
                198,
                199,
                200,
                193,
                134,
                135,
                139,
                138,
                125,
                155,
                156,
                157,
                150,
                211,
                168,
                28,
                90,
                94,
                93,
                176,
                180,
                179,
                109,
                115,
                104,
                195,
                201,
                190,
                133,
                137,
                136,
                152,
                158,
                147,
                27,
                93,
                94,
                98,
                97,
                179,
                180,
                184,
                183,
                115,
                116,
                117,
                118,
                201,
                202,
                203,
                204,
                136,
                137,
                141,
                140,
                127,
                158,
                159,
                160,
                161,
                213,
                170,
                27,
                94,
                95,
                99,
                98,
                180,
                181,
                185,
                184,
                108,
                119,
                120,
                116,
                194,
                205,
                206,
                202,
                137,
                138,
                142,
                141,
                128,
                151,
                162,
                163,
                159,
                214,
                171,
            ]
        )
        connIExp = DataArrayInt(
            [
                0,
                28,
                47,
                75,
                103,
                131,
                150,
                178,
                206,
                234,
                253,
                281,
                309,
                337,
                356,
                384,
                412,
            ]
        )
        #
        offset = m.getNumberOfNodes()
        nbExtrusions = 2
        ret = m.extrudeConnectivity(nbExtrusions)
        ret.checkConsistency()
        self.assertEqual(
            ret.getNumberOfNodes(), (2 * nbExtrusions + 1) * m.getNumberOfNodes()
        )
        self.assertTrue(ret.getNodalConnectivity().isEqual(connExp))
        self.assertTrue(ret.getNodalConnectivityIndex().isEqual(connIExp))
        for i in range((2 * nbExtrusions + 1)):
            self.assertTrue(
                ret.getCoords()[
                    i * m.getNumberOfNodes() : (i + 1) * m.getNumberOfNodes()
                ].isEqual(m.getCoords(), 1e-12)
            )

    def testUMesh_orientCorrectly3DCells(self):
        """
        [EDF32603] : MEDCouplingUMesh.orientCorrectly3DCells
        """
        connInit = [
            0,
            2,
            4,
            1,
            3,
            6,
            7,
            5,
            9,
            12,
            11,
            8,
            16,
            19,
            18,
            14,
            10,
            15,
            17,
            13,
        ]  # with negative volume
        connInverted = [
            0,
            1,
            4,
            2,
            3,
            5,
            7,
            6,
            8,
            11,
            12,
            9,
            14,
            18,
            19,
            16,
            10,
            13,
            17,
            15,
        ]  # with positive volume
        cooRef = DataArrayDouble(
            [
                -14.1150768258072,
                13.715,
                0.0,
                -14.115076825807153,
                13.715,
                0.25,
                -14.135716805787013,
                13.084513834221452,
                0.0,
                -14.6600768258072,
                13.715,
                0.0,
                -14.135716805786988,
                13.08451383422145,
                0.25,
                -14.660076825807154,
                13.715,
                0.25,
                -14.679633820990539,
                13.118881371469758,
                0.0,
                -14.679633820990517,
                13.11888137146976,
                0.25,
                -14.115076825807165,
                13.714999999999998,
                0.125,
                -14.120247868747509,
                13.399419437291247,
                0.0,
                -14.3875768258072,
                13.714999999999998,
                0.0,
                -14.120247868747473,
                13.399419437291245,
                0.25,
                -14.135716805786995,
                13.08451383422145,
                0.125,
                -14.387576825807155,
                13.714999999999998,
                0.25,
                -14.660076825807167,
                13.714999999999998,
                0.125,
                -14.407675313388776,
                13.101697602845606,
                0.0,
                -14.66497658818184,
                13.41662022366878,
                0.0,
                -14.407675313388753,
                13.101697602845604,
                0.25,
                -14.669855323398835,
                13.416940685734879,
                0.25,
                -14.679633820990528,
                13.118881371469758,
                0.125,
            ],
            20,
            3,
        )

        # -14.6796 to -14.1151
        def func(i):
            cell = MEDCouplingUMesh("", 3)
            cell.allocateCells()
            cell.insertNextCell(NORM_HEXA20, connInit if i % 2 == 0 else connInverted)
            cell.setCoords(cooRef)
            cell.translate([i * 0.8, 0, 0])
            return cell

        cells = MEDCouplingUMesh.MergeUMeshes([func(i) for i in range(6)])
        self.assertTrue(
            cells.getMeasureField(False)
            .getArray()
            .findIdsLowerThan(0.0)
            .isEqual(DataArrayInt([0, 2, 4]))
        )  # False is very important !
        refCoords0 = cells.getCoords().getHiddenCppPointer()
        refCoords1 = cells.getCoords().deepCopy()
        cells.orientCorrectly3DCells()  # key point is here !
        cellsCopy = cells.deepCopy()
        self.assertTrue(
            cells.getMeasureField(False).getArray().findIdsLowerThan(0.0).empty()
        )  # the aim of previous line !
        connToTest = MEDCoupling1SGTUMesh(cells).getNodalConnectivity()
        self.assertTrue(
            connToTest.isEqual(
                DataArrayInt.Aggregate(
                    [DataArrayInt(connInverted) + 20 * i for i in range(6)]
                )
            )
        )
        self.assertEqual(cells.getCoords().getHiddenCppPointer(), refCoords0)
        self.assertTrue(cells.getCoords().isEqual(refCoords1, 1e-12))
        # check that additionnal call does nothing
        cells.orientCorrectly3DCells()
        self.assertTrue(cells.isEqual(cellsCopy, 1e-12))

    def test_DAI_fromListOfPairsToIndexArray_1(self):
        """
        [EDF32671] : test of DataArrayInt.fromListOfPairsToIndexArray useful for joints management
        """
        arr = DataArrayInt([0, 3, 5, 4, 0, 7], 3, 2)
        c, ci = arr.fromListOfPairsToIndexArray()
        self.assertTrue(c.isEqual(DataArrayInt([0, 3, 7, 4, 5])))
        self.assertTrue(ci.isEqual(DataArrayInt([0, 3, 5])))

    def test_DAI_fromListOfPairsToIndexArray_2(self):
        """
        [EDF32671] : Deal with case of where P2 sends to P0 and sends to P1
        """
        arr = DataArrayInt([0, 7, 2, 7], 2, 2)
        c, ci = arr.fromListOfPairsToIndexArray()
        self.assertTrue(c.isEqual(DataArrayInt([0, 2, 7])))
        self.assertTrue(ci.isEqual(DataArrayInt([0, 3])))
        # even harder case
        arr = DataArrayInt([0, 3, 10, 17, 5, 4, 12, 17, 0, 107], 5, 2)
        c, ci = arr.fromListOfPairsToIndexArray()
        self.assertTrue(c.isEqual(DataArrayInt([0, 3, 107, 4, 5, 10, 12, 17])))
        self.assertTrue(ci.isEqual(DataArrayInt([0, 3, 5, 8])))


if __name__ == "__main__":
    unittest.main()
