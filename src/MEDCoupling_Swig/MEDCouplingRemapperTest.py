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

from MEDCouplingDataForTest import MEDCouplingDataForTest
from MEDCouplingRemapper import *

from math import *
import unittest


class MEDCouplingBasicsTest(unittest.TestCase):
    def testRemapper1(self):
        sourceMesh = self.build2DSourceMesh_1()
        targetMesh = self.build2DTargetMesh_1()
        remapper = MEDCouplingRemapper()
        remapper.setPrecision(1e-12)
        remapper.setIntersectionType(Triangulation)
        self.assertTrue(remapper.prepare(sourceMesh, targetMesh, "P0P0") == 1)
        srcField = MEDCouplingFieldDouble.New(ON_CELLS)
        srcField.setNature(IntensiveMaximum)
        srcField.setMesh(sourceMesh)
        array = DataArrayDouble.New()
        ptr = sourceMesh.getNumberOfCells() * [None]
        for i in range(sourceMesh.getNumberOfCells()):
            ptr[i] = float(i + 7)
            pass
        array.setValues(ptr, sourceMesh.getNumberOfCells(), 1)
        srcField.setArray(array)
        srcField.setName("abc")
        srcField.setDescription("def")
        srcField.setTime(7.7, 9, 10)
        trgfield = remapper.transferField(srcField, 4.57)
        self.assertEqual("abc", trgfield.getName())
        self.assertEqual("def", trgfield.getDescription())
        a, b, c = trgfield.getTime()
        self.assertAlmostEqual(7.7, a, 14)
        self.assertEqual(b, 9)
        self.assertEqual(c, 10)
        values = trgfield.getArray().getValues()
        valuesExpected = [7.5, 7.0, 7.0, 8.0, 7.5]
        for i in range(targetMesh.getNumberOfCells()):
            self.assertTrue(abs(values[i] - valuesExpected[i]) < 1e-12)
            pass
        self.assertTrue(1 == trgfield.getArray().getNumberOfComponents())
        pass

    def testPrepareEx1(self):
        sourceMesh = self.build2DSourceMesh_1()
        targetMesh = self.build2DTargetMesh_3()
        #
        remapper = MEDCouplingRemapper()
        remapper.setPrecision(1e-12)
        remapper.setIntersectionType(Triangulation)
        srcFt = MEDCouplingFieldTemplate.New(ON_CELLS)
        trgFt = MEDCouplingFieldTemplate.New(ON_CELLS)
        srcFt.setMesh(sourceMesh)
        trgFt.setMesh(targetMesh)
        self.assertEqual(1, remapper.prepareEx(srcFt, trgFt))
        srcField = MEDCouplingFieldDouble.New(ON_CELLS)
        srcField.setNature(IntensiveMaximum)
        srcField.setMesh(sourceMesh)
        array = DataArrayDouble.New()
        ptr = sourceMesh.getNumberOfCells() * [None]
        for i in range(sourceMesh.getNumberOfCells()):
            ptr[i] = float(i + 7)
            pass
        array.setValues(ptr, sourceMesh.getNumberOfCells(), 1)
        srcField.setArray(array)
        trgfield = remapper.transferField(srcField, 4.220173)
        values = trgfield.getArray().getValues()
        valuesExpected = [7.75, 7.0625, 4.220173, 8.0]
        self.assertEqual(4, trgfield.getArray().getNumberOfTuples())
        self.assertEqual(1, trgfield.getArray().getNumberOfComponents())
        for i0 in range(4):
            self.assertAlmostEqual(valuesExpected[i0], values[i0], 12)
            pass
        pass

    def testPartialTransfer1(self):
        sourceMesh = self.build2DSourceMesh_1()
        targetMesh = self.build2DTargetMesh_3()
        #
        remapper = MEDCouplingRemapper()
        remapper.setPrecision(1e-12)
        remapper.setIntersectionType(Triangulation)
        srcFt = MEDCouplingFieldTemplate.New(ON_CELLS)
        trgFt = MEDCouplingFieldTemplate.New(ON_CELLS)
        srcFt.setMesh(sourceMesh)
        trgFt.setMesh(targetMesh)
        self.assertEqual(1, remapper.prepareEx(srcFt, trgFt))
        srcField = MEDCouplingFieldDouble.New(ON_CELLS)
        srcField.setNature(IntensiveMaximum)
        srcField.setMesh(sourceMesh)
        array = DataArrayDouble.New()
        ptr = sourceMesh.getNumberOfCells() * [None]
        for i in range(sourceMesh.getNumberOfCells()):
            ptr[i] = float(i + 7)
            pass
        array.setValues(ptr, sourceMesh.getNumberOfCells(), 1)
        srcField.setArray(array)
        trgfield = MEDCouplingFieldDouble.New(ON_CELLS)
        trgfield.setNature(IntensiveMaximum)
        trgfield.setMesh(targetMesh)
        array = DataArrayDouble.New()
        ptr = targetMesh.getNumberOfCells() * [None]
        for i in range(targetMesh.getNumberOfCells()):
            ptr[i] = 4.220173
            pass
        array.setValues(ptr, targetMesh.getNumberOfCells(), 1)
        trgfield.setArray(array)
        remapper.partialTransfer(srcField, trgfield)
        values = trgfield.getArray().getValues()
        valuesExpected = [7.75, 7.0625, 4.220173, 8.0]
        self.assertEqual(4, trgfield.getArray().getNumberOfTuples())
        self.assertEqual(1, trgfield.getArray().getNumberOfComponents())
        for i0 in range(4):
            self.assertAlmostEqual(valuesExpected[i0], values[i0], 12)
            pass
        pass

    def testPrepareUC(self):
        # 1D
        coords = DataArrayDouble([0.0, 0.5, 0.7])
        src = MEDCouplingUMesh("", 1)
        src.setCoords(coords)
        src.allocateCells(2)
        src.insertNextCell(NORM_SEG2, [0, 1])
        src.insertNextCell(NORM_SEG2, [1, 2])
        src.finishInsertingCells()
        trg = MEDCouplingCMesh()
        arr = DataArrayDouble([-0.7, -0.1, 0.2, 0.7, 2.0, 2.3])
        trg.setCoordsAt(0, arr)
        fieldSrc = MEDCouplingFieldDouble(ON_CELLS, NO_TIME)
        fieldSrc.setMesh(src)
        arrSrc = DataArrayDouble([10.0, 30.0])
        fieldSrc.setNature(ExtensiveMaximum)
        fieldSrc.setArray(arrSrc)
        rem = MEDCouplingRemapper()
        rem.prepare(src, trg, "P0P0")
        trgField = rem.transferField(fieldSrc, -7.0)
        expected1 = [-7.0, 4.0, 36.0, -7.0, -7.0]
        self.assertEqual(5, trgField.getArray().getNumberOfTuples())
        self.assertEqual(5, len(expected1))
        for i, val in enumerate(expected1):
            self.assertAlmostEqual(expected1[i], trgField.getArray().getIJ(i, 0), 12)
            pass
        # 2D
        coords = DataArrayDouble(
            [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.5, -0.2], 5, 2
        )
        src = MEDCouplingUMesh("", 2)
        src.setCoords(coords)
        src.allocateCells(2)
        src.insertNextCell(NORM_TRI3, [0, 1, 2])
        src.insertNextCell(NORM_TRI3, [3, 4, 0])
        src.finishInsertingCells()
        trg = MEDCouplingCMesh()
        arr = DataArrayDouble([-0.7, -0.1, 0.2, 0.7, 2.0, 2.3])
        trg.setCoordsAt(0, arr)
        trg.setCoordsAt(1, arr)
        fieldSrc = MEDCouplingFieldDouble(ON_CELLS, NO_TIME)
        fieldSrc.setMesh(src)
        arrSrc = DataArrayDouble([10.0, 30.0])
        fieldSrc.setNature(ExtensiveMaximum)
        fieldSrc.setArray(arrSrc)
        rem = MEDCouplingRemapper()
        rem.prepare(src, trg, "P0P0")
        trgField = rem.transferField(fieldSrc, -7.0)
        expected2 = [
            -7.0,
            -7.0,
            7.35,
            0.15,
            -7.0,
            -7.0,
            2.8,
            14.85,
            5.25,
            -7.0,
            -7.0,
            2.0,
            2.5,
            -7.0,
            -7.0,
            -7.0,
            1.2,
            3.0,
            0.9,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
        ]
        self.assertEqual(25, trgField.getArray().getNumberOfTuples())
        self.assertEqual(25, len(expected2))
        for i, val in enumerate(expected2):
            self.assertAlmostEqual(expected2[i], trgField.getArray().getIJ(i, 0), 12)
            pass
        # 3D
        coords = DataArrayDouble(
            [
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                1.0,
                1.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.5,
                -0.2,
                0.0,
                0.1,
                0.8,
                1.0,
                0.5,
                0.0,
                1.0,
            ],
            7,
            3,
        )
        src = MEDCouplingUMesh("", 3)
        src.setCoords(coords)
        src.allocateCells(2)
        src.insertNextCell(NORM_TETRA4, [0, 1, 2, 5])
        src.insertNextCell(NORM_TETRA4, [3, 4, 0, 6])
        src.finishInsertingCells()
        trg = MEDCouplingCMesh()
        arr = DataArrayDouble([-0.7, -0.1, 0.2, 0.7, 2.0, 2.3])
        arr2 = DataArrayDouble([-0.7, 0.2, 0.6, 1.2, 2.0])
        trg.setCoordsAt(0, arr)
        trg.setCoordsAt(1, arr)
        trg.setCoordsAt(2, arr2)
        src.checkConsistency(1e-10)
        trg.checkConsistencyLight()
        fieldSrc = MEDCouplingFieldDouble(ON_CELLS, NO_TIME)
        fieldSrc.setMesh(src)
        arrSrc = DataArrayDouble([10.0, 30.0])
        fieldSrc.setNature(ExtensiveMaximum)
        fieldSrc.setArray(arrSrc)
        rem = MEDCouplingRemapper()
        rem.prepare(src, trg, "P0P0")
        trgField = rem.transferField(fieldSrc, -7.0)
        expected3 = [
            -7.0,
            -7.0,
            2.925,
            0.015,
            -7.0,
            -7.0,
            0.9392,
            8.595,
            2.265,
            -7.0,
            -7.0,
            1.1008,
            1.1192,
            -7.0,
            -7.0,
            -7.0,
            0.6392,
            1.6408,
            0.2808,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            0.81,
            -7.0,
            -7.0,
            -7.0,
            0.1208,
            11.55,
            0.96,
            -7.0,
            -7.0,
            1.1752,
            0.6592,
            -7.0,
            -7.0,
            -7.0,
            0.8512,
            1.7744,
            0.0192,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            1.92,
            -7.0,
            -7.0,
            -7.0,
            0.12578571428571422,
            0.007314285714285673,
            -7.0,
            -7.0,
            -7.0,
            0.3189253968253971,
            0.1879746031746033,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
            -7.0,
        ]
        self.assertEqual(100, trgField.getArray().getNumberOfTuples())
        self.assertEqual(100, len(expected3))
        for i, val in enumerate(expected3):
            self.assertAlmostEqual(expected3[i], trgField.getArray().getIJ(i, 0), 12)
            pass
        pass

    def testPrepareCU(self):
        # 1D
        coords = DataArrayDouble([0.0, 0.5, 0.7])
        trg = MEDCouplingUMesh("", 1)
        trg.setCoords(coords)
        trg.allocateCells(2)
        trg.insertNextCell(NORM_SEG2, [0, 1])
        trg.insertNextCell(NORM_SEG2, [1, 2])
        trg.finishInsertingCells()
        src = MEDCouplingCMesh()
        arr = DataArrayDouble([-0.7, -0.1, 0.2, 0.7, 2.0, 2.3])
        src.setCoordsAt(0, arr)
        fieldSrc = MEDCouplingFieldDouble(ON_CELLS, NO_TIME)
        fieldSrc.setMesh(src)
        arrTrg = DataArrayDouble([10.0, 30.0, 40.0, 70.0, 80.0])
        fieldSrc.setNature(ExtensiveMaximum)
        fieldSrc.setArray(arrTrg)
        rem = MEDCouplingRemapper()
        rem.prepare(src, trg, "P0P0")
        trgField = rem.transferField(fieldSrc, -7.0)
        expected1 = [44.0, 16.0]
        self.assertEqual(2.0, trgField.getArray().getNumberOfTuples())
        self.assertEqual(2, len(expected1))
        for i, val in enumerate(expected1):
            self.assertAlmostEqual(expected1[i], trgField.getArray().getIJ(i, 0), 12)
            pass
        # 2D
        coords = DataArrayDouble(
            [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.5, -0.2], 5, 2
        )
        trg = MEDCouplingUMesh("", 2)
        trg.setCoords(coords)
        trg.allocateCells(2)
        trg.insertNextCell(NORM_TRI3, [0, 1, 2])
        trg.insertNextCell(NORM_TRI3, [3, 4, 0])
        trg.finishInsertingCells()
        src = MEDCouplingCMesh()
        arr = DataArrayDouble([-0.7, -0.1, 0.2, 0.7, 2.0, 2.3])
        src.setCoordsAt(0, arr)
        src.setCoordsAt(1, arr)
        fieldSrc = MEDCouplingFieldDouble(ON_CELLS, NO_TIME)
        fieldSrc.setMesh(src)
        arrSrc = DataArrayDouble(
            [
                10.0,
                30.0,
                40.0,
                70.0,
                80.0,
                110.0,
                130.0,
                140.0,
                170.0,
                180.0,
                210.0,
                230.0,
                240.0,
                270.0,
                280.0,
                310.0,
                330.0,
                340.0,
                370.0,
                380.0,
                410.0,
                430.0,
                440.0,
                470.0,
                480.0,
            ]
        )
        fieldSrc.setNature(ExtensiveMaximum)
        fieldSrc.setArray(arrSrc)
        rem = MEDCouplingRemapper()
        rem.prepare(src, trg, "P0P0")
        trgField = rem.transferField(fieldSrc, -7.0)
        expected2 = [441.3050624589086, 68.69529914529915]
        self.assertEqual(2, trgField.getArray().getNumberOfTuples())
        self.assertEqual(2, len(expected2))
        for i, val in enumerate(expected2):
            self.assertAlmostEqual(expected2[i], trgField.getArray().getIJ(i, 0), 12)
            pass
        # 3D
        coords = DataArrayDouble(
            [
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                1.0,
                1.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.5,
                -0.2,
                0.0,
                0.1,
                0.8,
                1.0,
                0.5,
                0.0,
                1.0,
            ],
            7,
            3,
        )
        trg = MEDCouplingUMesh("", 3)
        trg.setCoords(coords)
        trg.allocateCells(2)
        trg.insertNextCell(NORM_TETRA4, [0, 1, 2, 5])
        trg.insertNextCell(NORM_TETRA4, [3, 4, 0, 6])
        trg.finishInsertingCells()
        src = MEDCouplingCMesh()
        arr = DataArrayDouble([-0.7, -0.1, 0.2, 0.7, 2.0, 2.3])
        arr2 = DataArrayDouble([-0.7, 0.2, 0.6, 1.2, 2.0])
        src.setCoordsAt(0, arr)
        src.setCoordsAt(1, arr)
        src.setCoordsAt(2, arr2)
        trg.checkConsistency(1e-10)
        src.checkConsistencyLight()
        fieldSrc = MEDCouplingFieldDouble(ON_CELLS, NO_TIME)
        fieldSrc.setMesh(src)
        arrSrc = DataArrayDouble(100)
        arrSrc.iota(7.7)
        fieldSrc.setNature(ExtensiveMaximum)
        fieldSrc.setArray(arrSrc)
        rem = MEDCouplingRemapper()
        rem.prepare(src, trg, "P0P0")
        trgField = rem.transferField(fieldSrc, -7.0)
        expected3 = [39.635196634558845, 12.13422356758468]
        self.assertEqual(2, trgField.getArray().getNumberOfTuples())
        self.assertEqual(2, len(expected3))
        for i, val in enumerate(expected3):
            self.assertAlmostEqual(expected3[i], trgField.getArray().getIJ(i, 0), 12)
            pass
        pass

    def testPrepareCC(self):
        # 1D
        src = MEDCouplingCMesh()
        arr = DataArrayDouble([-0.7, -0.1, 0.2, 0.7, 2.0, 2.3])
        src.setCoordsAt(0, arr)
        trg = MEDCouplingCMesh()
        arr = DataArrayDouble([-0.9, -0.1, 0.15])
        trg.setCoordsAt(0, arr)
        fieldSrc = MEDCouplingFieldDouble(ON_CELLS, NO_TIME)
        fieldSrc.setMesh(src)
        arrTrg = DataArrayDouble([10.0, 30.0, 40.0, 70.0, 80.0])
        fieldSrc.setNature(ExtensiveMaximum)
        fieldSrc.setArray(arrTrg)
        rem = MEDCouplingRemapper()
        rem.prepare(src, trg, "P0P0")
        trgField = rem.transferField(fieldSrc, -7.0)
        expected1 = [10.0, 25.0]
        self.assertEqual(2.0, trgField.getArray().getNumberOfTuples())
        self.assertEqual(2, len(expected1))
        for i, val in enumerate(expected1):
            self.assertAlmostEqual(expected1[i], trgField.getArray().getIJ(i, 0), 12)
            pass
        # 2D
        src = MEDCouplingCMesh()
        arr = DataArrayDouble([-0.7, -0.1, 0.2, 0.7, 2.0, 2.3])
        src.setCoordsAt(0, arr)
        src.setCoordsAt(1, arr)
        trg = MEDCouplingCMesh()
        arr = DataArrayDouble([-0.9, -0.1, 0.15])
        trg.setCoordsAt(0, arr)
        trg.setCoordsAt(1, arr)
        fieldSrc = MEDCouplingFieldDouble(ON_CELLS, NO_TIME)
        fieldSrc.setMesh(src)
        arrSrc = DataArrayDouble(
            [
                10.0,
                30.0,
                40.0,
                70.0,
                80.0,
                110.0,
                130.0,
                140.0,
                170.0,
                180.0,
                210.0,
                230.0,
                240.0,
                270.0,
                280.0,
                310.0,
                330.0,
                340.0,
                370.0,
                380.0,
                410.0,
                430.0,
                440.0,
                470.0,
                480.0,
            ]
        )
        fieldSrc.setNature(ExtensiveMaximum)
        fieldSrc.setArray(arrSrc)
        rem = MEDCouplingRemapper()
        rem.prepare(src, trg, "P0P0")
        trgField = rem.transferField(fieldSrc, -7.0)
        expected2 = [10.0, 25.0, 91.66666666666666, 90.27777777777777]
        self.assertEqual(4, trgField.getArray().getNumberOfTuples())
        self.assertEqual(4, len(expected2))
        for i, val in enumerate(expected2):
            self.assertAlmostEqual(expected2[i], trgField.getArray().getIJ(i, 0), 12)
            pass
        # 3D
        src = MEDCouplingCMesh()
        arr = DataArrayDouble([-0.7, -0.1, 0.2, 0.7, 2.0, 2.3])
        src.setCoordsAt(0, arr)
        src.setCoordsAt(1, arr)
        src.setCoordsAt(2, arr)
        trg = MEDCouplingCMesh()
        arr = DataArrayDouble([-0.9, -0.1, 0.15])
        trg.setCoordsAt(0, arr)
        trg.setCoordsAt(1, arr)
        trg.setCoordsAt(2, arr)
        fieldSrc = MEDCouplingFieldDouble(ON_CELLS, NO_TIME)
        fieldSrc.setMesh(src)
        arrSrc = DataArrayDouble(125)
        arrSrc.iota(7.7)
        fieldSrc.setNature(ExtensiveMaximum)
        fieldSrc.setArray(arrSrc)
        fieldSrc.checkConsistencyLight()
        rem = MEDCouplingRemapper()
        rem.prepare(src, trg, "P0P0")
        trgField = rem.transferField(fieldSrc, -7.0)
        expected3 = [
            7.7,
            7.249999999999999,
            10.583333333333332,
            9.513888888888886,
            27.25,
            23.40277777777777,
            26.180555555555546,
            22.39583333333333,
        ]
        self.assertEqual(8, trgField.getArray().getNumberOfTuples())
        self.assertEqual(8, len(expected3))
        for i, val in enumerate(expected3):
            self.assertAlmostEqual(expected3[i], trgField.getArray().getIJ(i, 0), 12)
            pass
        pass

    # Bug when source mesh is not homogeneously oriented in source mesh
    def testNonRegressionNonHomegenousOrriented3DCells(self):
        csrc = DataArrayDouble(
            [
                -0.15240000188350677,
                0,
                0,
                -0.1086929515004158,
                0,
                0,
                -0.15240000188350677,
                0.018142856657505035,
                0,
                -0.13054648041725159,
                0.0090714283287525177,
                0.019050000235438347,
                -0.13054648041725159,
                0.0090714283287525177,
                0,
            ],
            5,
            3,
        )
        src1 = MEDCouplingUMesh("src", 3)
        src1.allocateCells(0)
        src1.insertNextCell(NORM_TETRA4, [0, 1, 4, 3])
        src1.insertNextCell(NORM_TETRA4, [2, 0, 4, 3])
        src2 = MEDCouplingUMesh("src", 3)
        src2.allocateCells(0)
        src2.insertNextCell(NORM_TETRA4, [0, 4, 1, 3])
        src2.insertNextCell(NORM_TETRA4, [2, 0, 4, 3])
        src1.setCoords(csrc)
        src2.setCoords(csrc)
        ctrg = DataArrayDouble(
            [
                -0.15240000188350677,
                -0.038100000470876694,
                0,
                0.32379999756813049,
                -0.038100000470876694,
                0,
                -0.15240000188350677,
                0.076200000941753387,
                0,
                0.32379999756813049,
                0.076200000941753387,
                0,
                -0.15240000188350677,
                -0.038100000470876694,
                0.076200000941753387,
                0.32379999756813049,
                -0.038100000470876694,
                0.076200000941753387,
                -0.15240000188350677,
                0.076200000941753387,
                0.076200000941753387,
                0.32379999756813049,
                0.076200000941753387,
                0.076200000941753387,
            ],
            8,
            3,
        )
        trg = MEDCouplingUMesh("trg", 3)
        trg.allocateCells(0)
        trg.insertNextCell(NORM_HEXA8, [0, 1, 3, 2, 4, 5, 7, 6])
        trg.setCoords(ctrg)
        rem1 = MEDCouplingRemapper()
        rem1.setSplittingPolicy(PLANAR_FACE_5)
        rem1.prepare(src1, trg, "P0P0")
        rem2 = MEDCouplingRemapper()
        rem2.setSplittingPolicy(PLANAR_FACE_5)
        rem2.prepare(src1, trg, "P0P0")
        mat1 = rem1.getCrudeMatrix()
        mat2 = rem2.getCrudeMatrix()
        self.assertEqual(1, len(mat1))
        self.assertEqual(1, len(mat2))
        self.assertEqual(list(mat1[0].keys()), list(mat2[0].keys()))
        self.assertEqual([0, 1], list(mat1[0].keys()))
        self.assertAlmostEqual(1.25884108122e-06, mat1[0][0], 16)
        self.assertAlmostEqual(1.25884108122e-06, mat2[0][0], 16)
        self.assertAlmostEqual(1.25884086663e-06, mat1[0][1], 16)
        self.assertAlmostEqual(1.25884086663e-06, mat2[0][1], 16)
        #
        d = DataArrayDouble([13.45, 27.67], 2, 1)
        f1 = MEDCouplingFieldDouble(ON_CELLS, ONE_TIME)
        f1.setMesh(src1)
        f1.setArray(d)
        f1.setNature(IntensiveConservation)
        f2 = MEDCouplingFieldDouble(ON_CELLS, ONE_TIME)
        f2.setMesh(src2)
        f2.setArray(d)
        f2.setNature(IntensiveConservation)
        f11 = rem1.transferField(f1, 1e300)
        f22 = rem2.transferField(f2, 1e300)
        expected1 = DataArrayDouble([0.012480539537637884])
        self.assertTrue(f11.getArray().isEqual(expected1, 1e-15))
        self.assertTrue(f22.getArray().isEqual(expected1, 1e-15))
        #
        f1.setNature(ExtensiveMaximum)
        f2.setNature(ExtensiveMaximum)
        f11 = rem1.transferField(f1, 1e300)
        f22 = rem2.transferField(f2, 1e300)
        #
        expected2 = DataArrayDouble([41.12])
        self.assertTrue(f11.getArray().isEqual(expected2, 1e-13))
        self.assertTrue(f22.getArray().isEqual(expected2, 1e-13))
        pass

    def testCellToNodeReverse3D(self):
        c = DataArrayDouble([0.0, 1.0, 2.5])
        cc = MEDCouplingCMesh()
        cc.setCoords(c, c, c)
        um = cc.buildUnstructured()
        f = um.getMeasureField(False)
        #
        n2o = um.simplexize(PLANAR_FACE_5)
        f.setArray(f.getArray()[n2o])
        f.checkConsistencyLight()
        f.setNature(IntensiveMaximum)
        f.setTime(5.6, 7, 8)
        f.setName("toto")
        f.setDescription("aDescription")
        p = MEDCouplingRemapper()
        p.setIntersectionType(Barycentric)
        p.prepare(um, um, "P1P0")
        fNode = p.reverseTransferField(f, 1e300)
        self.assertEqual("toto", fNode.getName())
        self.assertEqual("aDescription", fNode.getDescription())
        a, b, c = fNode.getTime()
        self.assertAlmostEqual(5.6, a, 14)
        self.assertEqual(7, b)
        self.assertEqual(8, c)
        #
        integExpected = 34.328125
        self.assertAlmostEqual(fNode.integral(False)[0], integExpected, 14)
        self.assertAlmostEqual(f.integral(False)[0], integExpected, 14)
        pass

    def testGauss2Gauss2DValidated(self):
        srcFt = MEDCouplingDataForTest.buildFieldOnGauss_1()
        trgFt = MEDCouplingDataForTest.buildFieldOnGauss_2()
        src = MEDCouplingFieldDouble(srcFt)
        self.assertEqual(
            srcFt.getMesh().getHiddenCppPointer(), src.getMesh().getHiddenCppPointer()
        )
        self.assertEqual(
            srcFt.getDiscretization().getHiddenCppPointer(),
            src.getDiscretization().getHiddenCppPointer(),
        )
        # values given by ASTER usecase
        src.setArray(
            DataArrayDouble(
                [
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                ]
            )
        )
        src.getArray().setInfoOnComponents(["DOMA"])
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepareEx(srcFt, trgFt)
        trg = rem.transferField(src, 1e300)
        self.assertEqual(
            trg.getMesh().getHiddenCppPointer(), trgFt.getMesh().getHiddenCppPointer()
        )
        self.assertEqual(
            trg.getDiscretization().getHiddenCppPointer(),
            trgFt.getDiscretization().getHiddenCppPointer(),
        )
        # values given after interpolation in ASTER
        arrExpected = DataArrayDouble(
            [
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                1.0,
                1.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                1.0,
                1.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                1.0,
                1.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                1.0,
                1.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                1.0,
                1.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
            ]
        )
        arrExpected.setInfoOnComponents(["DOMA"])
        self.assertTrue(trg.getArray().isEqual(arrExpected, 1e-12))
        #
        # second part of the test : reverse source and target
        #
        rem.prepareEx(
            trgFt, srcFt
        )  # sorry trgFt is in the place of source and srcFt in the place of target it is not a bug
        trg = MEDCouplingFieldDouble(trgFt)
        # values given after interpolation in ASTER
        trg.setArray(
            DataArrayDouble(
                [
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                ]
            )
        )
        trg.getArray().setInfoOnComponents(["DOMA"])
        src = rem.transferField(trg, 1e300)
        # values given after interpolation in ASTER
        arrExpected2 = DataArrayDouble(
            [
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
                1.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                0.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
            ]
        )
        arrExpected2.setInfoOnComponents(["DOMA"])
        # modification of values in ASTER due to modification of algorithm
        # target PG 82 in target cell 32(C)/36 PG 1(C)/9 is in source cell 58(C)/120 source Gauss point 113 (1(C)/4). Values must be 1. and not 0.
        arrExpected2.setIJ(82, 0, 1.0)
        self.assertTrue(src.getArray().isEqual(arrExpected2, 1e-12))
        pass

    def testGauss2Gauss3DValidated(self):
        srcFt = MEDCouplingDataForTest.buildFieldOnGauss_3()
        trgFt = MEDCouplingDataForTest.buildFieldOnGauss_4()
        src = MEDCouplingFieldDouble(srcFt)
        self.assertEqual(
            srcFt.getMesh().getHiddenCppPointer(), src.getMesh().getHiddenCppPointer()
        )
        self.assertEqual(
            srcFt.getDiscretization().getHiddenCppPointer(),
            src.getDiscretization().getHiddenCppPointer(),
        )
        # values given by ASTER usecase
        src.setArray(
            DataArrayDouble(
                [
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                ]
            )
        )
        src.getArray().setInfoOnComponents(["DOMA"])
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepareEx(srcFt, trgFt)
        trg = rem.transferField(src, 1e300)
        self.assertEqual(
            trg.getMesh().getHiddenCppPointer(), trgFt.getMesh().getHiddenCppPointer()
        )
        self.assertEqual(
            trg.getDiscretization().getHiddenCppPointer(),
            trgFt.getDiscretization().getHiddenCppPointer(),
        )
        # values given after interpolation in ASTER
        arrExpected = DataArrayDouble(
            [
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
            ]
        )
        arrExpected.setInfoOnComponents(["DOMA"])
        self.assertTrue(trg.getArray().isEqual(arrExpected, 1e-12))
        #
        # second part of the test : reverse source and target
        #
        rem.prepareEx(
            trgFt, srcFt
        )  # sorry trgFt is in the place of source and srcFt in the place of target it is not a bug
        trg = MEDCouplingFieldDouble(trgFt)
        # values given after interpolation in ASTER
        trg.setArray(
            DataArrayDouble(
                [
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                ]
            )
        )
        trg.getArray().setInfoOnComponents(["DOMA"])
        src = rem.transferField(trg, 1e300)
        # values given after interpolation in ASTER
        arrExpected2 = DataArrayDouble(
            [
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
                0.0,
                1.0,
                0.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
                1.0,
                1.0,
                0.0,
                1.0,
            ]
        )
        arrExpected2.setInfoOnComponents(["DOMA"])
        self.assertTrue(src.getArray().isEqual(arrExpected2, 1e-12))
        pass

    def testSwig2MixOfUMesh(self):
        arr0 = DataArrayDouble([0, 1, 1.5])
        arr1 = DataArrayDouble([0, 1])
        sc = MEDCouplingCMesh()
        sc.setCoords(arr0, arr1, arr1)
        tc = sc.deepCopy()
        tc.translate([0.4, 0.3, 0.3])
        # umesh-umesh
        # 90 (umesh-1sgtumesh)
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(Triangulation)
        s = sc.buildUnstructured()
        t = tc.build1SGTUnstructured()
        self.assertTrue(isinstance(s, MEDCouplingUMesh))
        self.assertTrue(isinstance(t, MEDCoupling1SGTUMesh))
        rem.prepare(s, t, "P0P0")
        mat = rem.getCrudeMatrix()
        self.assertEqual(2, len(mat))
        self.assertEqual(2, len(mat[0]))
        self.assertEqual(1, len(mat[1]))
        self.assertAlmostEqual(0.294, mat[0][0], 14)
        self.assertAlmostEqual(0.196, mat[0][1], 14)
        self.assertAlmostEqual(0.049, mat[1][1], 14)
        del s, t
        # 91 (umesh-1dgtumesh)
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(Triangulation)
        s = sc.buildUnstructured()
        t = tc.buildUnstructured()
        t.convertAllToPoly()
        t = MEDCoupling1DGTUMesh(t)
        self.assertTrue(isinstance(s, MEDCouplingUMesh))
        self.assertTrue(isinstance(t, MEDCoupling1DGTUMesh))
        rem.prepare(s, t, "P0P0")
        mat = rem.getCrudeMatrix()
        self.assertEqual(2, len(mat))
        self.assertEqual(2, len(mat[0]))
        self.assertEqual(1, len(mat[1]))
        self.assertAlmostEqual(0.294, mat[0][0], 14)
        self.assertAlmostEqual(0.196, mat[0][1], 14)
        self.assertAlmostEqual(0.049, mat[1][1], 14)
        del s, t
        # 165 (1sgtumesh-umesh)
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(Triangulation)
        s = sc.build1SGTUnstructured()
        t = tc.buildUnstructured()
        self.assertTrue(isinstance(s, MEDCoupling1SGTUMesh))
        self.assertTrue(isinstance(t, MEDCouplingUMesh))
        rem.prepare(s, t, "P0P0")
        mat = rem.getCrudeMatrix()
        self.assertEqual(2, len(mat))
        self.assertEqual(2, len(mat[0]))
        self.assertEqual(1, len(mat[1]))
        self.assertAlmostEqual(0.294, mat[0][0], 14)
        self.assertAlmostEqual(0.196, mat[0][1], 14)
        self.assertAlmostEqual(0.049, mat[1][1], 14)
        del s, t
        # 181 (1dgtumesh-umesh
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(Triangulation)
        s = sc.buildUnstructured()
        s.convertAllToPoly()
        s = MEDCoupling1DGTUMesh(s)
        t = tc.buildUnstructured()
        self.assertTrue(isinstance(s, MEDCoupling1DGTUMesh))
        self.assertTrue(isinstance(t, MEDCouplingUMesh))
        rem.prepare(s, t, "P0P0")
        mat = rem.getCrudeMatrix()
        self.assertEqual(2, len(mat))
        self.assertEqual(2, len(mat[0]))
        self.assertEqual(1, len(mat[1]))
        self.assertAlmostEqual(0.294, mat[0][0], 14)
        self.assertAlmostEqual(0.196, mat[0][1], 14)
        self.assertAlmostEqual(0.049, mat[1][1], 14)
        del s, t
        # 170 (1sgtumesh-1sgtumesh)
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(Triangulation)
        s = sc.build1SGTUnstructured()
        t = tc.build1SGTUnstructured()
        self.assertTrue(isinstance(s, MEDCoupling1SGTUMesh))
        self.assertTrue(isinstance(t, MEDCoupling1SGTUMesh))
        rem.prepare(s, t, "P0P0")
        mat = rem.getCrudeMatrix()
        self.assertEqual(2, len(mat))
        self.assertEqual(2, len(mat[0]))
        self.assertEqual(1, len(mat[1]))
        self.assertAlmostEqual(0.294, mat[0][0], 14)
        self.assertAlmostEqual(0.196, mat[0][1], 14)
        self.assertAlmostEqual(0.049, mat[1][1], 14)
        del s, t
        # 171 (1sgtumesh-1dgtumesh)
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(Triangulation)
        s = sc.build1SGTUnstructured()
        t = tc.buildUnstructured()
        t.convertAllToPoly()
        t = MEDCoupling1DGTUMesh(t)
        self.assertTrue(isinstance(s, MEDCoupling1SGTUMesh))
        self.assertTrue(isinstance(t, MEDCoupling1DGTUMesh))
        rem.prepare(s, t, "P0P0")
        mat = rem.getCrudeMatrix()
        self.assertEqual(2, len(mat))
        self.assertEqual(2, len(mat[0]))
        self.assertEqual(1, len(mat[1]))
        self.assertAlmostEqual(0.294, mat[0][0], 14)
        self.assertAlmostEqual(0.196, mat[0][1], 14)
        self.assertAlmostEqual(0.049, mat[1][1], 14)
        del s, t
        # 186 (1dgtumesh-1sgtumesh)
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(Triangulation)
        s = sc.buildUnstructured()
        s.convertAllToPoly()
        s = MEDCoupling1DGTUMesh(s)
        t = tc.build1SGTUnstructured()
        self.assertTrue(isinstance(s, MEDCoupling1DGTUMesh))
        self.assertTrue(isinstance(t, MEDCoupling1SGTUMesh))
        rem.prepare(s, t, "P0P0")
        mat = rem.getCrudeMatrix()
        self.assertEqual(2, len(mat))
        self.assertEqual(2, len(mat[0]))
        self.assertEqual(1, len(mat[1]))
        self.assertAlmostEqual(0.294, mat[0][0], 14)
        self.assertAlmostEqual(0.196, mat[0][1], 14)
        self.assertAlmostEqual(0.049, mat[1][1], 14)
        del s, t
        # 187 (1dgtumesh-1dgtumesh)
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(Triangulation)
        s = sc.buildUnstructured()
        s.convertAllToPoly()
        s = MEDCoupling1DGTUMesh(s)
        t = tc.buildUnstructured()
        t.convertAllToPoly()
        t = MEDCoupling1DGTUMesh(t)
        self.assertTrue(isinstance(s, MEDCoupling1DGTUMesh))
        self.assertTrue(isinstance(t, MEDCoupling1DGTUMesh))
        rem.prepare(s, t, "P0P0")
        mat = rem.getCrudeMatrix()
        self.assertEqual(2, len(mat))
        self.assertEqual(2, len(mat[0]))
        self.assertEqual(1, len(mat[1]))
        self.assertAlmostEqual(0.294, mat[0][0], 14)
        self.assertAlmostEqual(0.196, mat[0][1], 14)
        self.assertAlmostEqual(0.049, mat[1][1], 14)
        del s, t
        # (umesh-cmesh)
        # 167 (1sgtumesh-cmesh)
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(Triangulation)
        s = sc.build1SGTUnstructured()
        t = tc.deepCopy()
        self.assertTrue(isinstance(s, MEDCoupling1SGTUMesh))
        self.assertTrue(isinstance(t, MEDCouplingCMesh))
        rem.prepare(s, t, "P0P0")
        mat = rem.getCrudeMatrix()
        self.assertEqual(2, len(mat))
        self.assertEqual(2, len(mat[0]))
        self.assertEqual(1, len(mat[1]))
        self.assertAlmostEqual(0.294, mat[0][0], 14)
        self.assertAlmostEqual(0.196, mat[0][1], 14)
        self.assertAlmostEqual(0.049, mat[1][1], 14)
        del s, t
        # 183 (1dgtumesh-cmesh)
        # rem=MEDCouplingRemapper() ; rem.setIntersectionType(Triangulation)
        # s=sc.buildUnstructured() ; s.convertAllToPoly() ; s=MEDCoupling1DGTUMesh(s) ; t=tc.deepCopy()
        # self.assertTrue(isinstance(s,MEDCoupling1DGTUMesh))
        # self.assertTrue(isinstance(t,MEDCouplingCMesh))
        # rem.prepare(s,t,"P0P0")
        # mat=rem.getCrudeMatrix()
        # self.assertEqual(2,len(mat)) ; self.assertEqual(2,len(mat[0])) ; self.assertEqual(1,len(mat[1]))
        # self.assertAlmostEqual(0.294,mat[0][0],14) ; self.assertAlmostEqual(0.196,mat[0][1],14) ; self.assertAlmostEqual(0.049,mat[1][1],14)
        # del s,t
        # (cmesh-umesh)
        # 122 (cmesh-1sgtumesh)
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(Triangulation)
        s = sc.deepCopy()
        t = tc.build1SGTUnstructured()
        self.assertTrue(isinstance(s, MEDCouplingCMesh))
        self.assertTrue(isinstance(t, MEDCoupling1SGTUMesh))
        rem.prepare(s, t, "P0P0")
        mat = rem.getCrudeMatrix()
        self.assertEqual(2, len(mat))
        self.assertEqual(2, len(mat[0]))
        self.assertEqual(1, len(mat[1]))
        self.assertAlmostEqual(0.294, mat[0][0], 14)
        self.assertAlmostEqual(0.196, mat[0][1], 14)
        self.assertAlmostEqual(0.049, mat[1][1], 14)
        del s, t
        # 123 (cmesh-1dgtumesh)
        # rem=MEDCouplingRemapper() ; rem.setIntersectionType(Triangulation)
        # s=sc.deepCopy() ; t=tc.buildUnstructured() ; t.convertAllToPoly() ; t=MEDCoupling1DGTUMesh(t)
        # self.assertTrue(isinstance(s,MEDCouplingCMesh))
        # self.assertTrue(isinstance(t,MEDCoupling1DGTUMesh))
        # rem.prepare(s,t,"P0P0")
        # mat=rem.getCrudeMatrix()
        # self.assertEqual(2,len(mat)) ; self.assertEqual(2,len(mat[0])) ; self.assertEqual(1,len(mat[1]))
        # self.assertAlmostEqual(0.294,mat[0][0],14) ; self.assertAlmostEqual(0.196,mat[0][1],14) ; self.assertAlmostEqual(0.049,mat[1][1],14)
        # del s,t
        pass

    def testSwig2BarycentricP1P13D_1(self):
        sCoo = DataArrayDouble(
            [
                0.313,
                0.00218,
                6.90489,
                0.313,
                0.10692667,
                6.90489,
                0.313,
                0.10692667,
                6.96790167,
                0.313,
                0.00218,
                6.9773125,
                0.313,
                0.21167333,
                6.90489,
                0.313,
                0.21167333,
                6.95849083,
                0.313,
                0.31642,
                6.90489,
                0.313,
                0.31642,
                6.94908,
                0.313,
                0.09383333,
                7.04891667,
                0.313,
                0.00218,
                7.049735,
                0.313,
                0.18548667,
                7.04809833,
                0.313,
                0.27714,
                7.04728,
                0.313,
                0.05782667,
                7.133205,
                0.313,
                0.00218,
                7.1221575,
                0.313,
                0.11347333,
                7.1442525,
                0.313,
                0.16912,
                7.1553,
                0.313,
                0.02509333,
                7.19458,
                0.313,
                0.00218,
                7.19458,
                0.313,
                0.04800667,
                7.19458,
                0.313,
                0.07092,
                7.19458,
                0.31005609,
                0.00218,
                6.90460005,
                0.31005609,
                0.10692667,
                6.90460005,
                0.29776312,
                0.10692667,
                6.96640097,
                0.29592716,
                0.00218,
                6.97563097,
                0.31005609,
                0.21167333,
                6.90460005,
                0.29959908,
                0.21167333,
                6.95717096,
                0.31005609,
                0.31642,
                6.90460005,
                0.30143505,
                0.31642,
                6.94794095,
                0.28195788,
                0.09383333,
                7.04585928,
                0.28179823,
                0.00218,
                7.04666189,
                0.28211753,
                0.18548667,
                7.04505668,
                0.28227718,
                0.27714,
                7.04425407,
                0.26551404,
                0.05782667,
                7.12852804,
                0.2676693,
                0.00218,
                7.11769282,
                0.26335878,
                0.11347333,
                7.13936327,
                0.26120352,
                0.16912,
                7.15019849,
                0.25354037,
                0.02509333,
                7.18872374,
                0.25354037,
                0.00218,
                7.18872374,
                0.25354037,
                0.04800667,
                7.18872374,
                0.25354037,
                0.07092,
                7.18872374,
                0.30722531,
                0.00218,
                6.90374134,
                0.30722531,
                0.10692667,
                6.90374134,
                0.28311179,
                0.10692667,
                6.96195653,
                0.27951042,
                0.00218,
                6.97065101,
                0.30722531,
                0.21167333,
                6.90374134,
                0.28671316,
                0.21167333,
                6.95326205,
                0.30722531,
                0.31642,
                6.90374134,
                0.29031453,
                0.31642,
                6.94456758,
                0.25210869,
                0.09383333,
                7.03680463,
                0.25179553,
                0.00218,
                7.03756067,
                0.25242185,
                0.18548667,
                7.03604859,
                0.25273501,
                0.27714,
                7.03529255,
                0.21985294,
                0.05782667,
                7.1146769,
                0.22408063,
                0.00218,
                7.10447034,
                0.21562524,
                0.11347333,
                7.12488346,
                0.21139755,
                0.16912,
                7.13509002,
                0.19636574,
                0.02509333,
                7.17138,
                0.19636574,
                0.00218,
                7.17138,
                0.19636574,
                0.04800667,
                7.17138,
                0.19636574,
                0.07092,
                7.17138,
                0.30461645,
                0.00218,
                6.90234688,
                0.30461645,
                0.10692667,
                6.90234688,
                0.26960904,
                0.10692667,
                6.95473916,
                0.26438066,
                0.00218,
                6.96256398,
                0.30461645,
                0.21167333,
                6.90234688,
                0.27483742,
                0.21167333,
                6.94691434,
                0.30461645,
                0.31642,
                6.90234688,
                0.2800658,
                0.31642,
                6.93908952,
                0.22459952,
                0.09383333,
                7.02210067,
                0.22414487,
                0.00218,
                7.02278109,
                0.22505416,
                0.18548667,
                7.02142025,
                0.2255088,
                0.27714,
                7.02073983,
                0.17777143,
                0.05782667,
                7.09218386,
                0.18390909,
                0.00218,
                7.0829982,
                0.17163377,
                0.11347333,
                7.10136952,
                0.1654961,
                0.16912,
                7.11055518,
                0.1436733,
                0.02509333,
                7.14321531,
                0.1436733,
                0.00218,
                7.14321531,
                0.1436733,
                0.04800667,
                7.14321531,
                0.1436733,
                0.07092,
                7.14321531,
                0.30232976,
                0.00218,
                6.90047024,
                0.30232976,
                0.10692667,
                6.90047024,
                0.25777378,
                0.10692667,
                6.94502622,
                0.25111932,
                0.00218,
                6.95168068,
                0.30232976,
                0.21167333,
                6.90047024,
                0.26442825,
                0.21167333,
                6.93837175,
                0.30232976,
                0.31642,
                6.90047024,
                0.27108271,
                0.31642,
                6.93171729,
                0.20048753,
                0.09383333,
                7.00231247,
                0.19990888,
                0.00218,
                7.00289112,
                0.20106618,
                0.18548667,
                7.00173382,
                0.20164482,
                0.27714,
                7.00115518,
                0.14088667,
                0.05782667,
                7.06191333,
                0.14869844,
                0.00218,
                7.05410156,
                0.13307491,
                0.11347333,
                7.06972509,
                0.12526315,
                0.16912,
                7.07753685,
                0.097488,
                0.02509333,
                7.105312,
                0.097488,
                0.00218,
                7.105312,
                0.097488,
                0.04800667,
                7.105312,
                0.097488,
                0.07092,
                7.105312,
                0.30045312,
                0.00218,
                6.89818355,
                0.30045312,
                0.10692667,
                6.89818355,
                0.24806084,
                0.10692667,
                6.93319096,
                0.24023602,
                0.00218,
                6.93841934,
                0.30045312,
                0.21167333,
                6.89818355,
                0.25588566,
                0.21167333,
                6.92796258,
                0.30045312,
                0.31642,
                6.89818355,
                0.26371048,
                0.31642,
                6.9227342,
                0.18069933,
                0.09383333,
                6.97820048,
                0.18001891,
                0.00218,
                6.97865513,
                0.18137975,
                0.18548667,
                6.97774584,
                0.18206017,
                0.27714,
                6.9772912,
                0.11061614,
                0.05782667,
                7.02502857,
                0.1198018,
                0.00218,
                7.01889091,
                0.10143048,
                0.11347333,
                7.03116623,
                0.09224482,
                0.16912,
                7.0373039,
                0.05958469,
                0.02509333,
                7.0591267,
                0.05958469,
                0.00218,
                7.0591267,
                0.05958469,
                0.04800667,
                7.0591267,
                0.05958469,
                0.07092,
                7.0591267,
                0.29905866,
                0.00218,
                6.89557469,
                0.29905866,
                0.10692667,
                6.89557469,
                0.24084347,
                0.10692667,
                6.91968821,
                0.23214899,
                0.00218,
                6.92328958,
                0.29905866,
                0.21167333,
                6.89557469,
                0.24953795,
                0.21167333,
                6.91608684,
                0.29905866,
                0.31642,
                6.89557469,
                0.25823242,
                0.31642,
                6.91248547,
                0.16599537,
                0.09383333,
                6.95069131,
                0.16523933,
                0.00218,
                6.95100447,
                0.16675141,
                0.18548667,
                6.95037815,
                0.16750745,
                0.27714,
                6.95006499,
                0.0881231,
                0.05782667,
                6.98294706,
                0.09832966,
                0.00218,
                6.97871937,
                0.07791654,
                0.11347333,
                6.98717476,
                0.06770998,
                0.16912,
                6.99140245,
                0.03142,
                0.02509333,
                7.00643426,
                0.03142,
                0.00218,
                7.00643426,
                0.03142,
                0.04800667,
                7.00643426,
                0.03142,
                0.07092,
                7.00643426,
                0.29819995,
                0.00218,
                6.89274391,
                0.29819995,
                0.10692667,
                6.89274391,
                0.23639903,
                0.10692667,
                6.90503688,
                0.22716903,
                0.00218,
                6.90687284,
                0.29819995,
                0.21167333,
                6.89274391,
                0.24562904,
                0.21167333,
                6.90320092,
                0.29819995,
                0.31642,
                6.89274391,
                0.25485905,
                0.31642,
                6.90136495,
                0.15694072,
                0.09383333,
                6.92084212,
                0.15613811,
                0.00218,
                6.92100177,
                0.15774332,
                0.18548667,
                6.92068247,
                0.15854593,
                0.27714,
                6.92052282,
                0.07427196,
                0.05782667,
                6.93728596,
                0.08510718,
                0.00218,
                6.9351307,
                0.06343673,
                0.11347333,
                6.93944122,
                0.05260151,
                0.16912,
                6.94159648,
                0.01407626,
                0.02509333,
                6.94925963,
                0.01407626,
                0.00218,
                6.94925963,
                0.01407626,
                0.04800667,
                6.94925963,
                0.01407626,
                0.07092,
                6.94925963,
                0.29792818,
                0.00218,
                6.89054043,
                0.29792818,
                0.10692667,
                6.89054043,
                0.23499241,
                0.10692667,
                6.89363227,
                0.22559291,
                0.00218,
                6.89409403,
                0.29792818,
                0.21167333,
                6.89054043,
                0.24439191,
                0.21167333,
                6.8931705,
                0.29792818,
                0.31642,
                6.89054043,
                0.25379141,
                0.31642,
                6.89270873,
                0.154075,
                0.09383333,
                6.89760748,
                0.15325765,
                0.00218,
                6.89764764,
                0.15489234,
                0.18548667,
                6.89756733,
                0.15570969,
                0.27714,
                6.89752718,
                0.06988819,
                0.05782667,
                6.90174332,
                0.08092238,
                0.00218,
                6.90120124,
                0.058854,
                0.11347333,
                6.90228539,
                0.04781981,
                0.16912,
                6.90282747,
                0.00858712,
                0.02509333,
                6.90475485,
                0.00858712,
                0.00218,
                6.90475485,
                0.00858712,
                0.04800667,
                6.90475485,
                0.00858712,
                0.07092,
                6.90475485,
                0.29791,
                0.00218,
                6.820902,
                0.29791,
                0.10692667,
                6.820902,
                0.23489833,
                0.10692667,
                6.820902,
                0.2254875,
                0.00218,
                6.820902,
                0.29791,
                0.21167333,
                6.820902,
                0.24430917,
                0.21167333,
                6.820902,
                0.29791,
                0.31642,
                6.820902,
                0.25372,
                0.31642,
                6.820902,
                0.15388333,
                0.09383333,
                6.820902,
                0.153065,
                0.00218,
                6.820902,
                0.15470167,
                0.18548667,
                6.820902,
                0.15552,
                0.27714,
                6.820902,
                0.069595,
                0.05782667,
                6.820902,
                0.0806425,
                0.00218,
                6.820902,
                0.0585475,
                0.11347333,
                6.820902,
                0.0475,
                0.16912,
                6.820902,
                0.00822,
                0.02509333,
                6.820902,
                0.00822,
                0.00218,
                6.820902,
                0.00822,
                0.04800667,
                6.820902,
                0.00822,
                0.07092,
                6.820902,
            ],
            200,
            3,
        )
        sConn = DataArrayInt(
            [
                0,
                1,
                2,
                3,
                20,
                21,
                22,
                23,
                1,
                4,
                5,
                2,
                21,
                24,
                25,
                22,
                4,
                6,
                7,
                5,
                24,
                26,
                27,
                25,
                3,
                2,
                8,
                9,
                23,
                22,
                28,
                29,
                2,
                5,
                10,
                8,
                22,
                25,
                30,
                28,
                5,
                7,
                11,
                10,
                25,
                27,
                31,
                30,
                9,
                8,
                12,
                13,
                29,
                28,
                32,
                33,
                8,
                10,
                14,
                12,
                28,
                30,
                34,
                32,
                10,
                11,
                15,
                14,
                30,
                31,
                35,
                34,
                13,
                12,
                16,
                17,
                33,
                32,
                36,
                37,
                12,
                14,
                18,
                16,
                32,
                34,
                38,
                36,
                14,
                15,
                19,
                18,
                34,
                35,
                39,
                38,
                20,
                21,
                22,
                23,
                40,
                41,
                42,
                43,
                21,
                24,
                25,
                22,
                41,
                44,
                45,
                42,
                24,
                26,
                27,
                25,
                44,
                46,
                47,
                45,
                23,
                22,
                28,
                29,
                43,
                42,
                48,
                49,
                22,
                25,
                30,
                28,
                42,
                45,
                50,
                48,
                25,
                27,
                31,
                30,
                45,
                47,
                51,
                50,
                29,
                28,
                32,
                33,
                49,
                48,
                52,
                53,
                28,
                30,
                34,
                32,
                48,
                50,
                54,
                52,
                30,
                31,
                35,
                34,
                50,
                51,
                55,
                54,
                33,
                32,
                36,
                37,
                53,
                52,
                56,
                57,
                32,
                34,
                38,
                36,
                52,
                54,
                58,
                56,
                34,
                35,
                39,
                38,
                54,
                55,
                59,
                58,
                40,
                41,
                42,
                43,
                60,
                61,
                62,
                63,
                41,
                44,
                45,
                42,
                61,
                64,
                65,
                62,
                44,
                46,
                47,
                45,
                64,
                66,
                67,
                65,
                43,
                42,
                48,
                49,
                63,
                62,
                68,
                69,
                42,
                45,
                50,
                48,
                62,
                65,
                70,
                68,
                45,
                47,
                51,
                50,
                65,
                67,
                71,
                70,
                49,
                48,
                52,
                53,
                69,
                68,
                72,
                73,
                48,
                50,
                54,
                52,
                68,
                70,
                74,
                72,
                50,
                51,
                55,
                54,
                70,
                71,
                75,
                74,
                53,
                52,
                56,
                57,
                73,
                72,
                76,
                77,
                52,
                54,
                58,
                56,
                72,
                74,
                78,
                76,
                54,
                55,
                59,
                58,
                74,
                75,
                79,
                78,
                60,
                61,
                62,
                63,
                80,
                81,
                82,
                83,
                61,
                64,
                65,
                62,
                81,
                84,
                85,
                82,
                64,
                66,
                67,
                65,
                84,
                86,
                87,
                85,
                63,
                62,
                68,
                69,
                83,
                82,
                88,
                89,
                62,
                65,
                70,
                68,
                82,
                85,
                90,
                88,
                65,
                67,
                71,
                70,
                85,
                87,
                91,
                90,
                69,
                68,
                72,
                73,
                89,
                88,
                92,
                93,
                68,
                70,
                74,
                72,
                88,
                90,
                94,
                92,
                70,
                71,
                75,
                74,
                90,
                91,
                95,
                94,
                73,
                72,
                76,
                77,
                93,
                92,
                96,
                97,
                72,
                74,
                78,
                76,
                92,
                94,
                98,
                96,
                74,
                75,
                79,
                78,
                94,
                95,
                99,
                98,
                80,
                81,
                82,
                83,
                100,
                101,
                102,
                103,
                81,
                84,
                85,
                82,
                101,
                104,
                105,
                102,
                84,
                86,
                87,
                85,
                104,
                106,
                107,
                105,
                83,
                82,
                88,
                89,
                103,
                102,
                108,
                109,
                82,
                85,
                90,
                88,
                102,
                105,
                110,
                108,
                85,
                87,
                91,
                90,
                105,
                107,
                111,
                110,
                89,
                88,
                92,
                93,
                109,
                108,
                112,
                113,
                88,
                90,
                94,
                92,
                108,
                110,
                114,
                112,
                90,
                91,
                95,
                94,
                110,
                111,
                115,
                114,
                93,
                92,
                96,
                97,
                113,
                112,
                116,
                117,
                92,
                94,
                98,
                96,
                112,
                114,
                118,
                116,
                94,
                95,
                99,
                98,
                114,
                115,
                119,
                118,
                100,
                101,
                102,
                103,
                120,
                121,
                122,
                123,
                101,
                104,
                105,
                102,
                121,
                124,
                125,
                122,
                104,
                106,
                107,
                105,
                124,
                126,
                127,
                125,
                103,
                102,
                108,
                109,
                123,
                122,
                128,
                129,
                102,
                105,
                110,
                108,
                122,
                125,
                130,
                128,
                105,
                107,
                111,
                110,
                125,
                127,
                131,
                130,
                109,
                108,
                112,
                113,
                129,
                128,
                132,
                133,
                108,
                110,
                114,
                112,
                128,
                130,
                134,
                132,
                110,
                111,
                115,
                114,
                130,
                131,
                135,
                134,
                113,
                112,
                116,
                117,
                133,
                132,
                136,
                137,
                112,
                114,
                118,
                116,
                132,
                134,
                138,
                136,
                114,
                115,
                119,
                118,
                134,
                135,
                139,
                138,
                120,
                121,
                122,
                123,
                140,
                141,
                142,
                143,
                121,
                124,
                125,
                122,
                141,
                144,
                145,
                142,
                124,
                126,
                127,
                125,
                144,
                146,
                147,
                145,
                123,
                122,
                128,
                129,
                143,
                142,
                148,
                149,
                122,
                125,
                130,
                128,
                142,
                145,
                150,
                148,
                125,
                127,
                131,
                130,
                145,
                147,
                151,
                150,
                129,
                128,
                132,
                133,
                149,
                148,
                152,
                153,
                128,
                130,
                134,
                132,
                148,
                150,
                154,
                152,
                130,
                131,
                135,
                134,
                150,
                151,
                155,
                154,
                133,
                132,
                136,
                137,
                153,
                152,
                156,
                157,
                132,
                134,
                138,
                136,
                152,
                154,
                158,
                156,
                134,
                135,
                139,
                138,
                154,
                155,
                159,
                158,
                140,
                141,
                142,
                143,
                160,
                161,
                162,
                163,
                141,
                144,
                145,
                142,
                161,
                164,
                165,
                162,
                144,
                146,
                147,
                145,
                164,
                166,
                167,
                165,
                143,
                142,
                148,
                149,
                163,
                162,
                168,
                169,
                142,
                145,
                150,
                148,
                162,
                165,
                170,
                168,
                145,
                147,
                151,
                150,
                165,
                167,
                171,
                170,
                149,
                148,
                152,
                153,
                169,
                168,
                172,
                173,
                148,
                150,
                154,
                152,
                168,
                170,
                174,
                172,
                150,
                151,
                155,
                154,
                170,
                171,
                175,
                174,
                153,
                152,
                156,
                157,
                173,
                172,
                176,
                177,
                152,
                154,
                158,
                156,
                172,
                174,
                178,
                176,
                154,
                155,
                159,
                158,
                174,
                175,
                179,
                178,
                160,
                161,
                162,
                163,
                180,
                181,
                182,
                183,
                161,
                164,
                165,
                162,
                181,
                184,
                185,
                182,
                164,
                166,
                167,
                165,
                184,
                186,
                187,
                185,
                163,
                162,
                168,
                169,
                183,
                182,
                188,
                189,
                162,
                165,
                170,
                168,
                182,
                185,
                190,
                188,
                165,
                167,
                171,
                170,
                185,
                187,
                191,
                190,
                169,
                168,
                172,
                173,
                189,
                188,
                192,
                193,
                168,
                170,
                174,
                172,
                188,
                190,
                194,
                192,
                170,
                171,
                175,
                174,
                190,
                191,
                195,
                194,
                173,
                172,
                176,
                177,
                193,
                192,
                196,
                197,
                172,
                174,
                178,
                176,
                192,
                194,
                198,
                196,
                174,
                175,
                179,
                178,
                194,
                195,
                199,
                198,
            ]
        )
        s = MEDCoupling1SGTUMesh("source", NORM_HEXA8)
        s.setCoords(sCoo)
        s.setNodalConnectivity(sConn)
        #
        tCoo = DataArrayDouble(
            [
                0.328,
                0.012,
                6.8598,
                0.328,
                0.168320184237353,
                6.8598,
                0.328,
                0.324640368474706,
                6.8598,
                0.328,
                0.0,
                6.8598,
                0.298,
                0.012,
                6.8598,
                0.1565,
                0.012,
                6.8598,
                0.180205346493166,
                0.144794653506834,
                6.8598,
                0.298,
                0.168320184237353,
                6.8598,
                0.0,
                0.012,
                6.8598,
                0.0916755774886107,
                0.233324422511389,
                6.8598,
                0.298,
                0.324640368474706,
                6.8598,
                0.298,
                0.0,
                6.8598,
                0.1565,
                0.0,
                6.8598,
                0.0,
                0.0,
                6.8598,
                0.328,
                0.012,
                7.2298,
                0.328,
                0.168320184237353,
                7.2298,
                0.328,
                0.324640368474706,
                7.2298,
                0.328,
                0.0,
                7.2298,
                0.298,
                0.012,
                7.2298,
                0.1565,
                0.012,
                7.2298,
                0.180205346493166,
                0.144794653506834,
                7.2298,
                0.298,
                0.168320184237353,
                7.2298,
                0.0,
                0.012,
                7.2298,
                0.0916755774886107,
                0.233324422511389,
                7.2298,
                0.298,
                0.324640368474706,
                7.2298,
                0.298,
                0.0,
                7.2298,
                0.1565,
                0.0,
                7.2298,
                0.0,
                0.0,
                7.2298,
            ],
            28,
            3,
        )
        tConn = DataArrayInt(
            [
                4,
                5,
                6,
                7,
                18,
                19,
                20,
                21,
                5,
                8,
                9,
                6,
                19,
                22,
                23,
                20,
                6,
                9,
                10,
                7,
                20,
                23,
                24,
                21,
                11,
                12,
                5,
                4,
                25,
                26,
                19,
                18,
                12,
                13,
                8,
                5,
                26,
                27,
                22,
                19,
                3,
                11,
                4,
                0,
                17,
                25,
                18,
                14,
                0,
                4,
                7,
                1,
                14,
                18,
                21,
                15,
                1,
                7,
                10,
                2,
                15,
                21,
                24,
                16,
            ]
        )
        t = MEDCoupling1SGTUMesh("target", NORM_HEXA8)
        t.setCoords(tCoo)
        t.setNodalConnectivity(tConn)
        #
        s.simplexize(PLANAR_FACE_5)
        aRemapper = MEDCouplingRemapper()
        aRemapper.setPrecision(1e-12)
        aRemapper.setIntersectionType(Barycentric)
        self.assertEqual(aRemapper.prepare(s, t, "P1P1"), 1)
        m = aRemapper.getCrudeMatrix()
        self.assertEqual(len(m), 28)
        for i in range(28):
            if i not in [5, 6]:
                self.assertEqual(len(m[i]), 0)
                pass
            pass
        self.assertEqual(len(m[5]), 4)
        self.assertEqual(len(m[6]), 4)
        self.assertAlmostEqual(0.10714286103952797, m[5][168], 12)
        self.assertAlmostEqual(0.35691534416938014, m[5][169], 12)
        self.assertAlmostEqual(0.04492099619713096, m[5][163], 12)
        self.assertAlmostEqual(0.49102079859396097, m[5][189], 12)
        self.assertAlmostEqual(0.14039089397104254, m[6][185], 12)
        self.assertAlmostEqual(0.16362822318261033, m[6][162], 12)
        self.assertAlmostEqual(0.3438363717836785, m[6][188], 12)
        self.assertAlmostEqual(0.3521445110626687, m[6][170], 12)
        pass

    def testSwig2MappedBarycentricP1P12D_1(self):
        """Testing mapped barycentric P1P1 projection
        (uses analytical mapping from square to arbitrary convex quadrangle)
        """
        n = 5
        sCoo = DataArrayDouble(n, 1)
        sCoo.iota(0.0)
        sCoo /= float(n - 1)
        m = MEDCouplingCMesh("target")
        m.setCoordsAt(0, sCoo)
        m.setCoordsAt(1, sCoo)
        tgt = m.buildUnstructured()
        coo = tgt.getCoords()
        orig = coo.deepCopy()
        orig[:, 0] = 10.0
        orig[:, 1] = 15.0
        pt_a = coo.deepCopy()
        pt_a[:, 0] = -0.3
        pt_a[:, 1] = 1.0
        pt_b = coo.deepCopy()
        pt_b[:, 0] = 2.0
        pt_b[:, 1] = 3.0
        pt_c = coo.deepCopy()
        pt_c[:, 0] = 1.0
        pt_c[:, 1] = 0.0
        # P = x*C+y*A + xy(B-A-C) + ORIGIN
        coo2 = (
            coo[:, 0] * pt_c
            + coo[:, 1] * pt_a
            + coo[:, 0] * coo[:, 1] * (pt_b - pt_a - pt_c)
            + orig
        )

        tgt.setCoords(coo2)

        sCoo = DataArrayDouble([0.0, 0.0, -0.3, 1.0, 2.0, 3.0, 1.0, 0.0], 4, 2)
        sCoo[:, 0] += 10.0
        sCoo[:, 1] += 15.0
        sConn = DataArrayInt([0, 1, 2, 3])
        s = MEDCoupling1SGTUMesh("source", NORM_QUAD4)
        s.setCoords(sCoo)
        s.setNodalConnectivity(sConn)
        #
        aRemapper = MEDCouplingRemapper()
        aRemapper.setPrecision(1e-12)
        aRemapper.setIntersectionType(MappedBarycentric)
        self.assertEqual(aRemapper.prepare(s, tgt, "P1P1"), 1)
        srcField = MEDCouplingFieldDouble(ON_NODES, ONE_TIME)
        srcField.setNature(IntensiveMaximum)
        srcField.setMesh(s)
        srcField.setName("field")
        srcField.setArray(DataArrayDouble([1.0, 2.0, 3.0, 4.0]))
        tgtF = aRemapper.transferField(srcField, 1e300)
        ref = [
            1.0,
            1.75,
            2.5,
            3.25,
            4.0,
            1.25,
            1.875,
            2.5,
            3.125,
            3.75,
            1.5,
            2.0,
            2.5,
            3.0,
            3.5,
            1.75,
            2.125,
            2.5,
            2.875,
            3.25,
            2.0,
            2.25,
            2.5,
            2.75,
            3.0,
        ]
        val = tgtF.getArray().getValues()
        for i, ref_v in enumerate(ref):
            self.assertAlmostEqual(ref_v, val[i])
        pass

    def testSwig2MappedBarycentricP1P13_1(self):
        """Testing mapped barycentric P1P1 projection in 3D (uses orthogonal distances to
        HEXA8 faces).
        Convention:
              0 ------ 3
             /|       /|
            / |      / |
           1 ------ 2  |
           |  |     |  |
           |  |     |  |
           |  4-----|- 7
           | /      | /
           5 ------ 6
        """
        n = 5
        sCoo = DataArrayDouble(n, 1)
        sCoo.iota(0.0)
        sCoo /= float(n - 1)
        m = MEDCouplingCMesh("target")
        m.setCoordsAt(0, sCoo)
        m.setCoordsAt(1, sCoo)
        m.setCoordsAt(2, sCoo)
        tgt = m.buildUnstructured()
        coo = tgt.getCoords()
        pt_0 = coo.deepCopy()
        pt_0[:, 0] = -0.3
        pt_0[:, 1] = 1.0
        pt_0[:, 2] = 1.0
        pt_1 = coo.deepCopy()
        pt_1[:, 0] = 0.0
        pt_1[:, 1] = 0.0
        pt_1[:, 2] = 1.0
        pt_2 = coo.deepCopy()
        pt_2[:, 0] = 1.0
        pt_2[:, 1] = 0.0
        pt_2[:, 2] = 1.0
        pt_3 = coo.deepCopy()
        pt_3[:, 0] = 2.0
        pt_3[:, 1] = 3.0
        pt_3[:, 2] = 1.0

        pt_4 = coo.deepCopy()
        pt_4[:, 0] = -0.3
        pt_4[:, 1] = 1.0
        pt_4[:, 2] = 0.0
        orig = coo.deepCopy()
        orig[:, 0] = 10.0
        orig[:, 1] = 15.0
        orig[:, 2] = 20.0
        pt_6 = coo.deepCopy()
        pt_6[:, 0] = 1.0
        pt_6[:, 1] = 0.0
        pt_6[:, 2] = 0.0
        pt_7 = coo.deepCopy()
        pt_7[:, 0] = 2.0
        pt_7[:, 1] = 3.0
        pt_7[:, 2] = 0.0
        # P = x*p6 + y*p4 + z*p1 + xy*(p7-p6-p4) + xz*(p2-p1-p6) + yz*(p0-p4-p1) + xyz(p3-p7-p2-p0+p1+p6+p4)
        x, y, z = coo[:, 0], coo[:, 1], coo[:, 2]
        coo2 = (
            x * pt_6
            + y * pt_4
            + z * pt_1
            + x * y * (pt_7 - pt_6 - pt_4)
            + x * z * (pt_2 - pt_1 - pt_6)
            + y * z * (pt_0 - pt_4 - pt_1)
            + x * y * z * (pt_3 - pt_7 - pt_2 - pt_0 + pt_6 + pt_1 + pt_4)
            + orig
        )
        tgt.setCoords(coo2)

        sCoo = DataArrayDouble(
            [
                -0.3,
                1.0,
                1.0,
                0.0,
                0.0,
                1.0,
                1.0,
                0.0,
                1.0,
                2.0,
                3.0,
                1.0,
                -0.3,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                2.0,
                3.0,
                0.0,
            ],
            8,
            3,
        )
        sCoo[:, 0] += 10.0
        sCoo[:, 1] += 15.0
        sCoo[:, 2] += 20.0
        sConn = DataArrayInt([0, 1, 2, 3, 4, 5, 6, 7])
        s = MEDCoupling1SGTUMesh("source", NORM_HEXA8)
        s.setCoords(sCoo)
        s.setNodalConnectivity(sConn)
        #
        aRemapper = MEDCouplingRemapper()
        aRemapper.setPrecision(1e-12)
        aRemapper.setIntersectionType(MappedBarycentric)
        self.assertEqual(aRemapper.prepare(s, tgt, "P1P1"), 1)
        srcField = MEDCouplingFieldDouble(ON_NODES, ONE_TIME)
        srcField.setNature(IntensiveMaximum)
        srcField.setMesh(s)
        srcField.setName("field")
        srcField.setArray(DataArrayDouble([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]))
        tgtF = aRemapper.transferField(srcField, 1e300)
        #        print tgtF.getArray().getValues()
        ref = [
            6.0,
            6.251802698104413,
            6.502397834044702,
            6.7517940736426665,
            7.0,
            5.740554726834594,
            6.1761835575796935,
            6.6052985689637564,
            7.009392769824465,
            7.383488834310164,
            5.487562931129931,
            6.140664596972973,
            6.720290674177548,
            7.220534970454015,
            7.651092836860121,
            5.2407867837524345,
            6.125759809889516,
            6.82853486793175,
            7.390880823876876,
            7.848445254819061,
            5.0,
            6.12211344611157,
            6.925740671133115,
            7.529623182840827,
            8.0,
            5.0,
            5.251802698104413,
            5.502397834044702,
            5.751794073642667,
            6.0,
            4.740554726834594,
            5.1761835575796935,
            5.6052985689637564,
            6.009392769824465,
            6.383488834310163,
            4.487562931129931,
            5.140664596972973,
            5.720290674177548,
            6.220534970454015,
            6.651092836860121,
            4.2407867837524345,
            5.125759809889516,
            5.828534867931749,
            6.390880823876876,
            6.848445254819061,
            4.0,
            5.122113446111569,
            5.925740671133115,
            6.529623182840827,
            7.0,
            4.0,
            4.251802698104413,
            4.502397834044702,
            4.751794073642667,
            5.0,
            3.740554726834594,
            4.176183557579693,
            4.6052985689637564,
            5.009392769824464,
            5.383488834310164,
            3.487562931129931,
            4.140664596972973,
            4.720290674177548,
            5.220534970454015,
            5.651092836860121,
            3.240786783752434,
            4.125759809889516,
            4.82853486793175,
            5.390880823876876,
            5.848445254819061,
            3.0,
            4.122113446111569,
            4.925740671133115,
            5.529623182840827,
            6.0,
            3.0,
            3.2518026981044135,
            3.502397834044702,
            3.7517940736426674,
            4.0,
            2.7405547268345933,
            3.176183557579693,
            3.6052985689637564,
            4.009392769824465,
            4.383488834310164,
            2.487562931129931,
            3.140664596972973,
            3.7202906741775474,
            4.220534970454015,
            4.65109283686012,
            2.2407867837524345,
            3.1257598098895154,
            3.828534867931749,
            4.390880823876876,
            4.848445254819061,
            2.0,
            3.1221134461115687,
            3.9257406711331146,
            4.529623182840826,
            5.0,
            2.0,
            2.2518026981044135,
            2.502397834044702,
            2.7517940736426674,
            3.0,
            1.7405547268345936,
            2.176183557579693,
            2.6052985689637564,
            3.0093927698244642,
            3.3834888343101635,
            1.4875629311299305,
            2.1406645969729734,
            2.720290674177548,
            3.2205349704540143,
            3.6510928368601205,
            1.2407867837524345,
            2.125759809889516,
            2.8285348679317495,
            3.390880823876876,
            3.848445254819061,
            1.0,
            2.1221134461115687,
            2.9257406711331146,
            3.529623182840827,
            4.0,
        ]

        val = tgtF.getArray().getValues()
        for i, ref_v in enumerate(ref):
            self.assertAlmostEqual(ref_v, val[i])
        pass

    @unittest.skipUnless(
        MEDCouplingHasNumPyBindings() and MEDCouplingHasSciPyBindings(),
        "requires numpy AND scipy",
    )
    def testGetCrudeCSRMatrix1(self):
        """testing CSR matrix output using numpy/scipy."""
        from scipy.sparse import spdiags  # diags
        import scipy
        from numpy import array

        arr = DataArrayDouble(3)
        arr.iota()
        m = MEDCouplingCMesh()
        m.setCoords(arr, arr)
        src = m.buildUnstructured()
        trg = src.deepCopy()
        trg = trg[[0, 1, 3]]
        trg.getCoords()[:] *= 0.5
        trg.getCoords()[:] += [0.3, 0.25]
        # Let's interpolate.
        rem = MEDCouplingRemapper()
        rem.prepare(src, trg, "P0P0")
        # Internal crude sparse matrix computed. Let's manipulate it using CSR matrix in scipy.
        for i in range(10):
            m = rem.getCrudeCSRMatrix()
            pass
        m2 = rem.getCrudeCSRMatrix()
        diff = m - m2
        self.assertTrue(isinstance(m, scipy.sparse.csr_matrix))
        self.assertEqual(m.getnnz(), 7)
        self.assertAlmostEqual(m[0, 0], 0.25, 12)
        self.assertAlmostEqual(m[1, 0], 0.1, 12)
        self.assertAlmostEqual(m[1, 1], 0.15, 12)
        self.assertAlmostEqual(m[2, 0], 0.05, 12)
        self.assertAlmostEqual(m[2, 1], 0.075, 12)
        self.assertAlmostEqual(m[2, 2], 0.05, 12)
        self.assertAlmostEqual(m[2, 3], 0.075, 12)
        self.assertEqual(diff.getnnz(), 0)
        # ExtensiveConservation (division by sum of cols)
        colSum = m.sum(axis=0)
        # version 0.12.0 # m_0=m*diags(array(1/colSum),[0])
        m_0 = m * spdiags(array(1 / colSum), [0], colSum.shape[1], colSum.shape[1])
        del colSum
        self.assertAlmostEqual(m_0[0, 0], 0.625, 12)
        self.assertAlmostEqual(m_0[1, 0], 0.25, 12)
        self.assertAlmostEqual(m_0[1, 1], 0.6666666666666667, 12)
        self.assertAlmostEqual(m_0[2, 0], 0.125, 12)
        self.assertAlmostEqual(m_0[2, 1], 0.3333333333333333, 12)
        self.assertAlmostEqual(m_0[2, 2], 1.0, 12)
        self.assertAlmostEqual(m_0[2, 3], 1.0, 12)
        self.assertEqual(m_0.getnnz(), 7)
        # IntensiveMaximum (division by sum of rows)
        rowSum = m.sum(axis=1)
        # version 0.12.0 # m_1=diags(array(1/rowSum.transpose()),[0])*m
        m_1 = (
            spdiags(
                array(1 / rowSum.transpose()), [0], rowSum.shape[0], rowSum.shape[0]
            )
            * m
        )
        del rowSum
        self.assertAlmostEqual(m_1[0, 0], 1.0, 12)
        self.assertAlmostEqual(m_1[1, 0], 0.4, 12)
        self.assertAlmostEqual(m_1[1, 1], 0.6, 12)
        self.assertAlmostEqual(m_1[2, 0], 0.2, 12)
        self.assertAlmostEqual(m_1[2, 1], 0.3, 12)
        self.assertAlmostEqual(m_1[2, 2], 0.2, 12)
        self.assertAlmostEqual(m_1[2, 3], 0.3, 12)
        self.assertEqual(m_1.getnnz(), 7)
        pass

    @unittest.skipUnless(
        MEDCouplingHasNumPyBindings() and MEDCouplingHasSciPyBindings(),
        "requires numpy AND scipy",
    )
    def testP0P1Bary_1(self):
        a = MEDCouplingUMesh("a", 2)
        a.allocateCells()
        conna = [
            0,
            1,
            3,
            2,
            1,
            4,
            5,
            3,
            4,
            6,
            7,
            5,
            6,
            8,
            9,
            7,
            8,
            10,
            11,
            9,
            10,
            12,
            13,
            11,
            12,
            14,
            15,
            13,
            14,
            16,
            17,
            15,
            16,
            18,
            19,
            17,
            18,
            20,
            21,
            19,
            20,
            22,
            23,
            21,
            22,
            24,
            25,
            23,
            24,
            26,
            27,
            25,
        ]
        a.setCoords(
            DataArrayDouble(
                [
                    1.54,
                    0,
                    -0.01,
                    1.54,
                    0.02,
                    -0.01,
                    1.54,
                    0,
                    0.01,
                    1.54,
                    0.02,
                    0.01,
                    1.54,
                    0.04,
                    -0.01,
                    1.54,
                    0.04,
                    0.01,
                    1.54,
                    0.06,
                    -0.01,
                    1.54,
                    0.06,
                    0.01,
                    1.54,
                    0.08,
                    -0.01,
                    1.54,
                    0.08,
                    0.01,
                    1.54,
                    0.1,
                    -0.01,
                    1.54,
                    0.1,
                    0.01,
                    1.54,
                    0.12,
                    -0.01,
                    1.54,
                    0.12,
                    0.01,
                    1.54,
                    0.14,
                    -0.01,
                    1.54,
                    0.14,
                    0.01,
                    1.54,
                    0.16,
                    -0.01,
                    1.54,
                    0.16,
                    0.01,
                    1.54,
                    0.18,
                    -0.01,
                    1.54,
                    0.18,
                    0.01,
                    1.54,
                    0.2,
                    -0.01,
                    1.54,
                    0.2,
                    0.01,
                    1.54,
                    0.22,
                    -0.01,
                    1.54,
                    0.22,
                    0.01,
                    1.54,
                    0.24,
                    -0.01,
                    1.54,
                    0.24,
                    0.01,
                    1.54,
                    0.26,
                    -0.01,
                    1.54,
                    0.26,
                    0.01,
                ],
                28,
                3,
            )
        )
        for i in range(13):
            a.insertNextCell(NORM_QUAD4, conna[4 * i : 4 * (i + 1)])
            pass
        a.finishInsertingCells()
        a.simplexize(0)
        #
        connb = [
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
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
            22,
            23,
            24,
            25,
            26,
            27,
            28,
            29,
            30,
            31,
            32,
            33,
            34,
            35,
            36,
            37,
            38,
            0,
            2,
            39,
            3,
            5,
            40,
            6,
            8,
            41,
            9,
            11,
            42,
            12,
            14,
            43,
            15,
            17,
            44,
            18,
            20,
            45,
            21,
            23,
            46,
            24,
            26,
            47,
            27,
            29,
            48,
            30,
            32,
            49,
            33,
            35,
            50,
            36,
            38,
            51,
            52,
            2,
            39,
            53,
            5,
            40,
            54,
            8,
            41,
            55,
            11,
            42,
            56,
            14,
            43,
            57,
            17,
            44,
            58,
            20,
            45,
            59,
            23,
            46,
            60,
            26,
            47,
            61,
            29,
            48,
            62,
            32,
            49,
            63,
            35,
            50,
            64,
            38,
            51,
            52,
            2,
            65,
            53,
            5,
            66,
            54,
            8,
            67,
            55,
            11,
            68,
            56,
            14,
            69,
            57,
            17,
            70,
            58,
            20,
            71,
            59,
            23,
            72,
            60,
            26,
            73,
            61,
            29,
            74,
            62,
            32,
            75,
            63,
            35,
            76,
            64,
            38,
            77,
            53,
            2,
            65,
            54,
            5,
            66,
            55,
            8,
            67,
            56,
            11,
            68,
            57,
            14,
            69,
            58,
            17,
            70,
            59,
            20,
            71,
            60,
            23,
            72,
            61,
            26,
            73,
            62,
            29,
            74,
            63,
            32,
            75,
            64,
            35,
            76,
            78,
            38,
            77,
            53,
            2,
            40,
            54,
            5,
            41,
            55,
            8,
            42,
            56,
            11,
            43,
            57,
            14,
            44,
            58,
            17,
            45,
            59,
            20,
            46,
            60,
            23,
            47,
            61,
            26,
            48,
            62,
            29,
            49,
            63,
            32,
            50,
            64,
            35,
            51,
            78,
            38,
            79,
            3,
            2,
            40,
            6,
            5,
            41,
            9,
            8,
            42,
            12,
            11,
            43,
            15,
            14,
            44,
            18,
            17,
            45,
            21,
            20,
            46,
            24,
            23,
            47,
            27,
            26,
            48,
            30,
            29,
            49,
            33,
            32,
            50,
            36,
            35,
            51,
            80,
            38,
            79,
            3,
            2,
            1,
            6,
            5,
            4,
            9,
            8,
            7,
            12,
            11,
            10,
            15,
            14,
            13,
            18,
            17,
            16,
            21,
            20,
            19,
            24,
            23,
            22,
            27,
            26,
            25,
            30,
            29,
            28,
            33,
            32,
            31,
            36,
            35,
            34,
            80,
            38,
            37,
        ]
        b = MEDCouplingUMesh("b", 2)
        b.allocateCells()
        for i in range(104):
            b.insertNextCell(NORM_TRI3, connb[3 * i : 3 * (i + 1)])
            pass
        b.setCoords(
            DataArrayDouble(
                [
                    1.54,
                    0,
                    -0.01,
                    1.54,
                    0.01,
                    -0.01,
                    1.54,
                    0.01,
                    0,
                    1.54,
                    0.02,
                    -0.01,
                    1.54,
                    0.03,
                    -0.01,
                    1.54,
                    0.03,
                    0,
                    1.54,
                    0.04,
                    -0.01,
                    1.54,
                    0.05,
                    -0.01,
                    1.54,
                    0.05,
                    0,
                    1.54,
                    0.06,
                    -0.01,
                    1.54,
                    0.07,
                    -0.01,
                    1.54,
                    0.07,
                    0,
                    1.54,
                    0.08,
                    -0.01,
                    1.54,
                    0.09,
                    -0.01,
                    1.54,
                    0.09,
                    0,
                    1.54,
                    0.1,
                    -0.01,
                    1.54,
                    0.11,
                    -0.01,
                    1.54,
                    0.11,
                    0,
                    1.54,
                    0.12,
                    -0.01,
                    1.54,
                    0.13,
                    -0.01,
                    1.54,
                    0.13,
                    0,
                    1.54,
                    0.14,
                    -0.01,
                    1.54,
                    0.15,
                    -0.01,
                    1.54,
                    0.15,
                    0,
                    1.54,
                    0.16,
                    -0.01,
                    1.54,
                    0.17,
                    -0.01,
                    1.54,
                    0.17,
                    0,
                    1.54,
                    0.18,
                    -0.01,
                    1.54,
                    0.19,
                    -0.01,
                    1.54,
                    0.19,
                    0,
                    1.54,
                    0.2,
                    -0.01,
                    1.54,
                    0.21,
                    -0.01,
                    1.54,
                    0.21,
                    0,
                    1.54,
                    0.22,
                    -0.01,
                    1.54,
                    0.23,
                    -0.01,
                    1.54,
                    0.23,
                    0,
                    1.54,
                    0.24,
                    -0.01,
                    1.54,
                    0.25,
                    -0.01,
                    1.54,
                    0.25,
                    0,
                    1.54,
                    0,
                    0,
                    1.54,
                    0.02,
                    0,
                    1.54,
                    0.04,
                    0,
                    1.54,
                    0.06,
                    0,
                    1.54,
                    0.08,
                    0,
                    1.54,
                    0.1,
                    0,
                    1.54,
                    0.12,
                    0,
                    1.54,
                    0.14,
                    0,
                    1.54,
                    0.16,
                    0,
                    1.54,
                    0.18,
                    0,
                    1.54,
                    0.2,
                    0,
                    1.54,
                    0.22,
                    0,
                    1.54,
                    0.24,
                    0,
                    1.54,
                    0,
                    0.01,
                    1.54,
                    0.02,
                    0.01,
                    1.54,
                    0.04,
                    0.01,
                    1.54,
                    0.06,
                    0.01,
                    1.54,
                    0.08,
                    0.01,
                    1.54,
                    0.1,
                    0.01,
                    1.54,
                    0.12,
                    0.01,
                    1.54,
                    0.14,
                    0.01,
                    1.54,
                    0.16,
                    0.01,
                    1.54,
                    0.18,
                    0.01,
                    1.54,
                    0.2,
                    0.01,
                    1.54,
                    0.22,
                    0.01,
                    1.54,
                    0.24,
                    0.01,
                    1.54,
                    0.01,
                    0.01,
                    1.54,
                    0.03,
                    0.01,
                    1.54,
                    0.05,
                    0.01,
                    1.54,
                    0.07,
                    0.01,
                    1.54,
                    0.09,
                    0.01,
                    1.54,
                    0.11,
                    0.01,
                    1.54,
                    0.13,
                    0.01,
                    1.54,
                    0.15,
                    0.01,
                    1.54,
                    0.17,
                    0.01,
                    1.54,
                    0.19,
                    0.01,
                    1.54,
                    0.21,
                    0.01,
                    1.54,
                    0.23,
                    0.01,
                    1.54,
                    0.25,
                    0.01,
                    1.54,
                    0.26,
                    0.01,
                    1.54,
                    0.26,
                    0,
                    1.54,
                    0.26,
                    -0.01,
                ],
                81,
                3,
            )
        )
        #
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(Barycentric)
        rem.prepare(a, b, "P1P0")
        m0 = rem.getCrudeCSRMatrix()
        self.assertEqual(m0.nnz, 312)
        #
        ids = 4 * [None]
        vs = 4 * [None]
        ids[0] = DataArrayInt(
            [
                0,
                3,
                6,
                9,
                12,
                15,
                18,
                21,
                24,
                27,
                30,
                33,
                36,
                39,
                42,
                45,
                48,
                51,
                54,
                57,
                60,
                63,
                66,
                69,
                72,
                75,
                158,
                161,
                164,
                167,
                170,
                173,
                176,
                179,
                182,
                185,
                188,
                191,
                194,
                197,
                200,
                203,
                206,
                209,
                212,
                215,
                218,
                221,
                224,
                227,
                230,
                233,
            ]
        )
        vs[0] = 10.0 / 3.0
        ids[1] = DataArrayInt(
            [
                1,
                2,
                4,
                5,
                7,
                8,
                10,
                11,
                13,
                14,
                16,
                17,
                19,
                20,
                22,
                23,
                25,
                26,
                28,
                29,
                31,
                32,
                34,
                35,
                37,
                38,
                40,
                41,
                43,
                44,
                46,
                47,
                49,
                50,
                52,
                53,
                55,
                56,
                58,
                59,
                61,
                62,
                64,
                65,
                67,
                68,
                70,
                71,
                73,
                74,
                76,
                77,
                80,
                83,
                86,
                89,
                92,
                95,
                98,
                101,
                104,
                107,
                110,
                113,
                116,
                117,
                120,
                123,
                126,
                129,
                132,
                135,
                138,
                141,
                144,
                147,
                150,
                153,
                156,
                157,
                159,
                160,
                162,
                163,
                165,
                166,
                168,
                169,
                171,
                172,
                174,
                175,
                177,
                178,
                180,
                181,
                183,
                184,
                186,
                187,
                189,
                190,
                192,
                193,
                195,
                196,
                198,
                199,
                201,
                202,
                204,
                205,
                207,
                208,
                210,
                211,
                213,
                214,
                216,
                217,
                219,
                220,
                222,
                223,
                225,
                226,
                228,
                229,
                231,
                232,
                234,
                237,
                240,
                243,
                246,
                249,
                252,
                255,
                258,
                261,
                264,
                267,
                270,
                275,
                278,
                281,
                284,
                287,
                290,
                293,
                296,
                299,
                302,
                305,
                308,
                311,
            ]
        )
        vs[1] = 5.0 / 6.0
        ids[2] = DataArrayInt(
            [
                78,
                81,
                84,
                87,
                90,
                93,
                96,
                99,
                102,
                105,
                108,
                111,
                114,
                119,
                122,
                125,
                128,
                131,
                134,
                137,
                140,
                143,
                146,
                149,
                152,
                155,
                236,
                239,
                242,
                245,
                248,
                251,
                254,
                257,
                260,
                263,
                266,
                269,
                272,
                273,
                276,
                279,
                282,
                285,
                288,
                291,
                294,
                297,
                300,
                303,
                306,
                309,
            ]
        )
        vs[2] = 5.0 / 3.0
        ids[3] = DataArrayInt(
            [
                79,
                82,
                85,
                88,
                91,
                94,
                97,
                100,
                103,
                106,
                109,
                112,
                115,
                118,
                121,
                124,
                127,
                130,
                133,
                136,
                139,
                142,
                145,
                148,
                151,
                154,
                235,
                238,
                241,
                244,
                247,
                250,
                253,
                256,
                259,
                262,
                265,
                268,
                271,
                274,
                277,
                280,
                283,
                286,
                289,
                292,
                295,
                298,
                301,
                304,
                307,
                310,
            ]
        )
        vs[3] = 2.5
        vals = DataArrayDouble(312, 1)
        for idd, v in zip(ids, vs):
            vals[idd] = v
            pass
        vals *= 1e-5
        eps0 = DataArrayDouble(m0.data) - vals
        eps0.abs()
        self.assertTrue(eps0.findIdsInRange(1e-17, 1e300).empty())
        self.assertTrue(
            DataArrayInt32(m0.indices).isEqual(
                DataArrayInt32(
                    [
                        0,
                        1,
                        3,
                        1,
                        4,
                        5,
                        4,
                        6,
                        7,
                        6,
                        8,
                        9,
                        8,
                        10,
                        11,
                        10,
                        12,
                        13,
                        12,
                        14,
                        15,
                        14,
                        16,
                        17,
                        16,
                        18,
                        19,
                        18,
                        20,
                        21,
                        20,
                        22,
                        23,
                        22,
                        24,
                        25,
                        24,
                        26,
                        27,
                        0,
                        2,
                        3,
                        1,
                        3,
                        5,
                        4,
                        5,
                        7,
                        6,
                        7,
                        9,
                        8,
                        9,
                        11,
                        10,
                        11,
                        13,
                        12,
                        13,
                        15,
                        14,
                        15,
                        17,
                        16,
                        17,
                        19,
                        18,
                        19,
                        21,
                        20,
                        21,
                        23,
                        22,
                        23,
                        25,
                        24,
                        25,
                        27,
                        0,
                        2,
                        3,
                        1,
                        3,
                        5,
                        4,
                        5,
                        7,
                        6,
                        7,
                        9,
                        8,
                        9,
                        11,
                        10,
                        11,
                        13,
                        12,
                        13,
                        15,
                        14,
                        15,
                        17,
                        16,
                        17,
                        19,
                        18,
                        19,
                        21,
                        20,
                        21,
                        23,
                        22,
                        23,
                        25,
                        24,
                        25,
                        27,
                        0,
                        2,
                        3,
                        1,
                        3,
                        5,
                        4,
                        5,
                        7,
                        6,
                        7,
                        9,
                        8,
                        9,
                        11,
                        10,
                        11,
                        13,
                        12,
                        13,
                        15,
                        14,
                        15,
                        17,
                        16,
                        17,
                        19,
                        18,
                        19,
                        21,
                        20,
                        21,
                        23,
                        22,
                        23,
                        25,
                        24,
                        25,
                        27,
                        0,
                        2,
                        3,
                        1,
                        3,
                        5,
                        4,
                        5,
                        7,
                        6,
                        7,
                        9,
                        8,
                        9,
                        11,
                        10,
                        11,
                        13,
                        12,
                        13,
                        15,
                        14,
                        15,
                        17,
                        16,
                        17,
                        19,
                        18,
                        19,
                        21,
                        20,
                        21,
                        23,
                        22,
                        23,
                        25,
                        24,
                        25,
                        27,
                        0,
                        1,
                        3,
                        1,
                        4,
                        5,
                        4,
                        6,
                        7,
                        6,
                        8,
                        9,
                        8,
                        10,
                        11,
                        10,
                        12,
                        13,
                        12,
                        14,
                        15,
                        14,
                        16,
                        17,
                        16,
                        18,
                        19,
                        18,
                        20,
                        21,
                        20,
                        22,
                        23,
                        22,
                        24,
                        25,
                        24,
                        26,
                        27,
                        0,
                        1,
                        3,
                        1,
                        4,
                        5,
                        4,
                        6,
                        7,
                        6,
                        8,
                        9,
                        8,
                        10,
                        11,
                        10,
                        12,
                        13,
                        12,
                        14,
                        15,
                        14,
                        16,
                        17,
                        16,
                        18,
                        19,
                        18,
                        20,
                        21,
                        20,
                        22,
                        23,
                        22,
                        24,
                        25,
                        24,
                        26,
                        27,
                        0,
                        1,
                        3,
                        1,
                        4,
                        5,
                        4,
                        6,
                        7,
                        6,
                        8,
                        9,
                        8,
                        10,
                        11,
                        10,
                        12,
                        13,
                        12,
                        14,
                        15,
                        14,
                        16,
                        17,
                        16,
                        18,
                        19,
                        18,
                        20,
                        21,
                        20,
                        22,
                        23,
                        22,
                        24,
                        25,
                        24,
                        26,
                        27,
                    ]
                )
            )
        )
        self.assertTrue(
            DataArrayInt32(m0.indptr).isEqual(
                DataArrayInt32(
                    [
                        0,
                        3,
                        6,
                        9,
                        12,
                        15,
                        18,
                        21,
                        24,
                        27,
                        30,
                        33,
                        36,
                        39,
                        42,
                        45,
                        48,
                        51,
                        54,
                        57,
                        60,
                        63,
                        66,
                        69,
                        72,
                        75,
                        78,
                        81,
                        84,
                        87,
                        90,
                        93,
                        96,
                        99,
                        102,
                        105,
                        108,
                        111,
                        114,
                        117,
                        120,
                        123,
                        126,
                        129,
                        132,
                        135,
                        138,
                        141,
                        144,
                        147,
                        150,
                        153,
                        156,
                        159,
                        162,
                        165,
                        168,
                        171,
                        174,
                        177,
                        180,
                        183,
                        186,
                        189,
                        192,
                        195,
                        198,
                        201,
                        204,
                        207,
                        210,
                        213,
                        216,
                        219,
                        222,
                        225,
                        228,
                        231,
                        234,
                        237,
                        240,
                        243,
                        246,
                        249,
                        252,
                        255,
                        258,
                        261,
                        264,
                        267,
                        270,
                        273,
                        276,
                        279,
                        282,
                        285,
                        288,
                        291,
                        294,
                        297,
                        300,
                        303,
                        306,
                        309,
                        312,
                    ]
                )
            )
        )
        #
        rem2 = MEDCouplingRemapper()
        rem2.setIntersectionType(Barycentric)
        rem2.prepare(b, a, "P0P1")
        m1 = rem2.getCrudeCSRMatrix()
        self.assertEqual(m1.nnz, 312)
        #
        m1 = rem2.getCrudeCSRMatrix()
        m1t = m1.transpose()
        delta = m0 - m1t
        self.assertTrue(DataArrayDouble(delta.data).isUniform(0.0, 1e-17))
        pass

    @unittest.skipUnless(
        MEDCouplingHasNumPyBindings() and MEDCouplingHasSciPyBindings(),
        "requires numpy AND scipy",
    )
    def testNonConformWithRemapper_1(self):
        coo = DataArrayDouble(
            [
                -0.396700000780411,
                -0.134843245350081,
                -0.0361311386958691,
                -0.407550009429364,
                -0.13484324535008,
                -0.0361311386958923,
                -0.396700000780411,
                -0.132191446077668,
                -0.0448729493559049,
                -0.407550009429364,
                -0.132191446077666,
                -0.0448729493559254,
                -0.396700000780411,
                -0.128973582738749,
                -0.0534226071577727,
                -0.407550009429364,
                -0.128973582738747,
                -0.0534226071577904,
                -0.396700000780411,
                -0.128348829636458,
                -0.0346583696473619,
                -0.407550009429364,
                -0.128348829636457,
                -0.0346583696473822,
                -0.396700000780411,
                -0.125874740261886,
                -0.0430683597970123,
                -0.407550009429364,
                -0.125874740261885,
                -0.0430683597970302,
                -0.396700000780411,
                -0.122905344829122,
                -0.051310216195766,
                -0.407550009429364,
                -0.12290534482912,
                -0.0513102161957814,
            ],
            12,
            3,
        )
        conn = DataArrayInt(
            [
                2,
                9,
                3,
                11,
                2,
                3,
                5,
                11,
                2,
                8,
                9,
                11,
                2,
                10,
                8,
                11,
                2,
                5,
                4,
                11,
                2,
                4,
                10,
                11,
                3,
                0,
                1,
                6,
                3,
                1,
                7,
                6,
                3,
                2,
                0,
                6,
                3,
                8,
                2,
                6,
                3,
                7,
                9,
                6,
                3,
                9,
                8,
                6,
            ]
        )
        m = MEDCoupling1SGTUMesh("mesh", NORM_TETRA4)
        m.setNodalConnectivity(conn)
        m.setCoords(coo)
        # m is ready
        m1, d, di, rd, rdi = m.buildUnstructured().buildDescendingConnectivity()
        rdi2 = rdi.deltaShiftIndex()
        cellIds = rdi2.findIdsEqual(1)
        skinAndNonConformCells = m1[cellIds]
        skinAndNonConformCells.zipCoords()  # at this point skinAndNonConformCells contains non conform cells and skin cells. Now trying to split them in two parts.
        #
        rem = MEDCouplingRemapper()
        rem.setMaxDistance3DSurfIntersect(1e-12)
        rem.setMinDotBtwPlane3DSurfIntersect(
            0.99
        )  # this line is important it is to tell to remapper to select only cells with very close orientation
        rem.prepare(skinAndNonConformCells, skinAndNonConformCells, "P0P0")
        mat = rem.getCrudeCSRMatrix()
        indptr = DataArrayInt32(mat.indptr)  # not depend on MEDCouplingUse64BitIDs()
        indptr2 = indptr.deltaShiftIndex()
        cellIdsOfNonConformCells = indptr2.findIdsNotEqual(1)
        cellIdsOfSkin = indptr2.findIdsEqual(1)
        self.assertTrue(
            cellIdsOfSkin.isEqual(
                DataArrayInt(
                    [
                        1,
                        2,
                        3,
                        5,
                        6,
                        7,
                        8,
                        9,
                        10,
                        11,
                        12,
                        13,
                        14,
                        15,
                        16,
                        17,
                        19,
                        20,
                        21,
                        23,
                    ]
                )
            )
        )
        self.assertTrue(cellIdsOfNonConformCells.isEqual(DataArrayInt([0, 4, 18, 22])))
        pass

    def test3D1DOnP1P0_1(self):
        """This test focused on P1P0 interpolation with a source with meshDim=1 spaceDim=3 and a target with meshDim=3.
        This test has revealed a bug in remapper. A reverse matrix is computed so a reverse method should be given in input.
        """
        target = MEDCouplingCMesh()
        arrX = DataArrayDouble([0, 1])
        arrY = DataArrayDouble([0, 1])
        arrZ = DataArrayDouble(11)
        arrZ.iota()
        target.setCoords(arrX, arrY, arrZ)
        target = target.buildUnstructured()
        target.setName("TargetSecondaire")
        #
        sourceCoo = DataArrayDouble(
            [
                (0.5, 0.5, 0.1),
                (0.5, 0.5, 1.2),
                (0.5, 0.5, 1.6),
                (0.5, 0.5, 1.8),
                (0.5, 0.5, 2.43),
                (0.5, 0.5, 2.55),
                (0.5, 0.5, 4.1),
                (0.5, 0.5, 4.4),
                (0.5, 0.5, 4.9),
                (0.5, 0.5, 5.1),
                (0.5, 0.5, 7.6),
                (0.5, 0.5, 7.7),
                (0.5, 0.5, 8.2),
                (0.5, 0.5, 8.4),
                (0.5, 0.5, 8.6),
                (0.5, 0.5, 8.8),
                (0.5, 0.5, 9.2),
                (0.5, 0.5, 9.6),
                (0.5, 0.5, 11.5),
            ]
        )
        source = MEDCoupling1SGTUMesh("SourcePrimaire", NORM_SEG2)
        source.setCoords(sourceCoo)
        source.allocateCells()
        for i in range(len(sourceCoo) - 1):
            source.insertNextCell([i, i + 1])
            pass
        source = source.buildUnstructured()
        fsource = MEDCouplingFieldDouble(ON_NODES)
        fsource.setName("field")
        fsource.setMesh(source)
        arr = DataArrayDouble(len(sourceCoo))
        arr.iota(0.7)
        arr *= arr
        fsource.setArray(arr)
        fsource.setNature(IntensiveMaximum)
        #
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(source, target, "P1P0")
        f2Test = rem.transferField(fsource, -27)
        self.assertEqual(f2Test.getName(), fsource.getName())
        self.assertEqual(
            f2Test.getMesh().getHiddenCppPointer(), target.getHiddenCppPointer()
        )
        expArr = DataArrayDouble(
            [
                0.49,
                7.956666666666667,
                27.29,
                -27,
                59.95666666666667,
                94.09,
                -27,
                125.69,
                202.89,
                296.09,
            ]
        )
        self.assertTrue(f2Test.getArray().isEqual(expArr, 1e-12))
        f2Test = rem.reverseTransferField(f2Test, -36)
        self.assertEqual(f2Test.getName(), fsource.getName())
        self.assertEqual(
            f2Test.getMesh().getHiddenCppPointer(), source.getHiddenCppPointer()
        )
        expArr2 = DataArrayDouble(
            [
                0.49,
                7.956666666666667,
                7.956666666666667,
                7.956666666666667,
                27.29,
                27.29,
                59.95666666666667,
                59.95666666666667,
                59.95666666666667,
                94.09,
                125.69,
                125.69,
                202.89,
                202.89,
                202.89,
                202.89,
                296.09,
                296.09,
                -36.0,
            ]
        )
        self.assertTrue(f2Test.getArray().isEqual(expArr2, 1e-12))
        pass

    def testRemapperAMR1(self):
        """This test is the origin of the ref values for MEDCouplingBasicsTest.testAMR2"""
        coarse = DataArrayDouble(35)
        coarse.iota(0)  # X=5,Y=7
        fine = DataArrayDouble(3 * 2 * 4 * 4)
        fine.iota(0)  # X=3,Y=2 refined by 4
        MEDCouplingIMesh.CondenseFineToCoarse(
            [5, 7], fine, [(1, 4), (2, 4)], [4, 4], coarse
        )
        #
        m = MEDCouplingCartesianAMRMesh("mesh", 2, [6, 8], [0.0, 0.0], [1.0, 1.0])
        trgMesh = m.buildUnstructured()
        m.addPatch([(1, 4), (2, 4)], [4, 4])
        srcMesh = m[0].getMesh().buildUnstructured()
        srcField = MEDCouplingFieldDouble(ON_CELLS)
        fine2 = DataArrayDouble(3 * 2 * 4 * 4)
        fine2.iota(0)
        srcField.setArray(fine2)
        srcField.setMesh(srcMesh)
        srcField.setNature(ExtensiveMaximum)
        #
        trgField = MEDCouplingFieldDouble(ON_CELLS)
        coarse2 = DataArrayDouble(35)
        coarse2.iota(0)
        trgField.setArray(coarse2)
        trgField.setMesh(trgMesh)
        trgField.setNature(ExtensiveMaximum)
        #
        rem = MEDCouplingRemapper()
        rem.prepare(srcMesh, trgMesh, "P0P0")
        rem.partialTransfer(srcField, trgField)
        #
        self.assertTrue(coarse.isEqual(trgField.getArray(), 1e-12))
        pass

    @unittest.skipUnless(
        MEDCouplingHasNumPyBindings() and MEDCouplingHasSciPyBindings(),
        "requires numpy AND scipy",
    )
    def test1DPointLocator1(self):
        """This test focuses on PointLocator for P1P1 in 1D and 2DCurve."""
        from numpy import array
        from scipy.sparse import diags, csr_matrix, identity

        ## basic case 1D
        arrS = DataArrayInt.Range(0, 11, 1).convertToDblArr()
        arrT = DataArrayDouble([0.1, 1.7, 5.5, 9.6])
        mS = MEDCouplingCMesh()
        mS.setCoords(arrS)
        mT = MEDCouplingCMesh()
        mT.setCoords(arrT)
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        self.assertEqual(
            rem.prepare(mS.buildUnstructured(), mT.buildUnstructured(), "P1P1"), 1
        )
        m = rem.getCrudeCSRMatrix()
        rowSum = m.sum(axis=1)
        m = diags(array(1 / rowSum.transpose()), [0]) * m
        # expected matrix
        row = array([0, 0, 1, 1, 2, 2, 3, 3])
        col = array([0, 1, 1, 2, 5, 6, 9, 10])
        data = array([0.9, 0.1, 0.3, 0.7, 0.5, 0.5, 0.4, 0.6])
        mExp0 = csr_matrix((data, (row, col)), shape=(4, 11))
        # compute diff and check
        diff = abs(m - mExp0)
        self.assertAlmostEqual(diff.sum(), 0.0, 14)
        ## full specific case 1D where target=source
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        self.assertEqual(
            rem.prepare(mS.buildUnstructured(), mS.buildUnstructured(), "P1P1"), 1
        )
        m = rem.getCrudeCSRMatrix()
        rowSum = m.sum(axis=1)
        m = diags(array(1 / rowSum.transpose()), [0]) * m
        # expected matrix
        mExp1 = identity(11)
        diff = abs(m - mExp1)
        self.assertAlmostEqual(diff.sum(), 0.0, 14)
        ## case where some points in target are not in source
        arrT = DataArrayDouble([-0.2, 0.1, 1.7, 5.5, 10.3])
        mT = MEDCouplingCMesh()
        mT.setCoords(arrT)
        mT = mT.buildUnstructured()
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        self.assertEqual(rem.prepare(mS.buildUnstructured(), mT, "P1P1"), 1)
        m = rem.getCrudeCSRMatrix()
        row = array([1, 1, 2, 2, 3, 3])
        col = array([0, 1, 1, 2, 5, 6])
        data = array([1.8, 0.2, 0.6, 1.4, 1.0, 1.0])
        mExp2 = csr_matrix((data, (row, col)), shape=(5, 11))
        diff = abs(m - mExp2)
        self.assertAlmostEqual(diff.sum(), 0.0, 14)
        ## basic case 2D Curve
        arrS = DataArrayInt.Range(0, 11, 1).convertToDblArr()
        arrT = DataArrayDouble([0.1, 1.7, 5.5, 9.6])
        mS = MEDCouplingCMesh()
        mS.setCoords(arrS)
        mT = MEDCouplingCMesh()
        mT.setCoords(arrT)
        mS = mS.buildUnstructured()
        mS.changeSpaceDimension(2)
        mT = mT.buildUnstructured()
        mT.changeSpaceDimension(2)
        mS.rotate([-1.0, -1.0], 1.2)
        mT.rotate([-1.0, -1.0], 1.2)
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        self.assertEqual(rem.prepare(mS, mT, "P1P1"), 1)
        m = rem.getCrudeCSRMatrix()
        rowSum = m.sum(axis=1)
        m = diags(array(1 / rowSum.transpose()), [0]) * m
        diff = abs(m - mExp0)
        self.assertAlmostEqual(diff.sum(), 0.0, 14)
        pass

    def test3D2Dand2D3DPointLocator1(self):
        """Non regression test solving SIGSEGV when using 3D<->3Dsurf pointlocator."""
        arrX = DataArrayDouble([0, 1, 2])
        arrY = DataArrayDouble([0, 1])
        arrZ = DataArrayDouble([0, 1])
        ms = MEDCouplingCMesh()
        ms.setCoords(arrX, arrY, arrZ)
        ms = ms.buildUnstructured()
        ms.setName("source")
        #
        mt = MEDCouplingUMesh("target", 2)
        mt.allocateCells()
        mt.insertNextCell(NORM_TRI3, [0, 4, 6])
        mt.insertNextCell(NORM_TRI3, [1, 5, 7])
        mt.setCoords(ms.getCoords()[:])
        mt.zipCoords()
        #
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(ms, mt, "P0P0")
        self.assertEqual(rem.getCrudeMatrix(), [{0: 1.0}, {1: 1.0}])
        rem2 = MEDCouplingRemapper()
        rem2.setIntersectionType(PointLocator)
        ##
        # 2D to 3D with point locator does not make sense:
        ##
        self.assertRaises(InterpKernelException, rem2.prepare, mt, ms, "P0P0")
        pass

    def test2D1Dand1D2DPointLocator1(self):
        arrX = DataArrayDouble([0, 1, 2])
        arrY = DataArrayDouble([0, 1])
        ms = MEDCouplingCMesh()
        ms.setCoords(arrX, arrY)
        ms = ms.buildUnstructured()
        mt = MEDCouplingUMesh("target", 1)
        mt.setCoords(ms.getCoords()[:])
        mt.allocateCells()
        mt.insertNextCell(NORM_SEG2, [0, 4])
        mt.insertNextCell(NORM_SEG2, [1, 5])
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(ms, mt, "P0P0")
        self.assertEqual(rem.getCrudeMatrix(), [{0: 1.0}, {1: 1.0}])
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(mt, ms, "P0P0")
        self.assertEqual(rem.getCrudeMatrix(), [{0: 1.0}, {1: 1.0}])
        pass

    def test3D1DPointLocatorBBoxAdjusted(self):
        """In case a 1D segment lies exactly on the interface between two 2D (or 3D) faces, the default
        bounding box logic will make it non-intersecting with the surrounding 2D (or 3D) faces.
        Test bounding box adjustment allowing to widen the BB to capture this.
        """
        m = MEDCouplingCMesh("source")
        di, dd = DataArrayInt, DataArrayDouble
        m.setCoordsAt(0, dd([0.0, 1.0, 2.0]))
        m.setCoordsAt(1, dd([0.0, 1.0]))
        m.setCoordsAt(2, dd([0.0, 1.0]))
        m3d = m.buildUnstructured()
        m1d = MEDCouplingUMesh("target", 1)
        m1d.setCoords(dd([1.0, 0.5, 0.2, 1.0, 0.5, 0.8], 2, 3))
        m1d.setConnectivity(di([NORM_SEG2, 0, 1]), di([0, 3]))

        rem = MEDCouplingRemapper()
        rem.setPrecision(1e-12)
        rem.setIntersectionType(PointLocator)
        rem.prepare(m3d, m1d, "P0P1")
        self.assertEqual(rem.getCrudeMatrix(), [{0: 1.0, 1: 1.0}, {0: 1.0, 1: 1.0}])

        rem = MEDCouplingRemapper()
        rem.setPrecision(1e-12)
        rem.setIntersectionType(PointLocator)
        rem.setBoundingBoxAdjustment(0.0)
        rem.setBoundingBoxAdjustmentAbs(0.0)
        rem.prepare(m3d, m1d, "P0P1")
        self.assertEqual(rem.getCrudeMatrix(), [{}, {}])
        pass

    def testPointLocator3DTo2D(self):
        """Target mesh has spaceDim==3 and meshDim==2. Source has spaceDim==3 and meshDim==3. Here we are on pointlocator alg.
        The test evaluates on each nodes of target mesh the bary coor into source mesh."""
        src = MEDCouplingCMesh()
        arr = DataArrayDouble([0, 1, 2])
        src.setCoords(arr, arr, arr)
        src = src.buildUnstructured()
        src.simplexize(PLANAR_FACE_5)
        fsrc = MEDCouplingFieldDouble(ON_NODES)
        fsrc.setMesh(src)
        fsrc.setArray(
            DataArrayDouble(
                [
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
                    13,
                    14,
                    15,
                    16,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                ]
            )
        )
        #
        trg = MEDCouplingCMesh()
        arr = DataArrayDouble([0, 1])
        trg.setCoords(arr, arr)
        trg = trg.buildUnstructured()
        trg.changeSpaceDimension(3, 0.0)
        trg.translate([0.5, 0.5, 0.5])
        #
        arrTrg = fsrc.getValueOnMulti(trg.getCoords())
        ftrg = MEDCouplingFieldDouble(ON_NODES)
        ftrg.setMesh(trg)
        ftrg.setArray(arrTrg)
        ftrg.checkConsistencyLight()
        ftrg.setNature(IntensiveMaximum)
        #
        fsrc.setNature(IntensiveMaximum)
        remap = MEDCouplingRemapper()
        remap.setIntersectionType(PointLocator)
        self.assertEqual(remap.prepare(src, trg, "P1P1"), 1)
        ftrg2 = remap.transferField(fsrc, 1e300)
        self.assertTrue(ftrg.isEqual(ftrg2, 1e-12, 1e-12))
        pass

    def testPointLocator2D2DNonConvexPolygons(self):
        """PointLocator remapper now correclty support non-convex polygons"""
        src = MEDCouplingUMesh("src", 2)
        coo = DataArrayDouble(
            [
                (6, 1),
                (6, 2),
                (4, 2),
                (4, 3),
                (3, 3),
                (3, 4),
                (2, 4),
                (2, 6),
                (1, 6),
                (1, 8),
                (2, 8),
                (2, 9),
                (3, 9),
                (3, 8),
                (4, 8),
                (4, 9),
                (5, 9),
                (5, 8),
                (6, 8),
                (6, 9),
                (7, 9),
                (7, 8),
                (8, 8),
                (8, 9),
                (9, 9),
                (9, 8),
                (10, 8),
                (10, 9),
                (11, 9),
                (11, 8),
                (12, 8),
                (12, 9),
                (13, 9),
                (13, 8),
                (14, 8),
                (14, 9),
                (15, 9),
                (15, 8),
                (16, 8),
                (16, 6),
                (15, 6),
                (15, 4),
                (14, 4),
                (14, 3),
                (13, 3),
                (13, 2),
                (11, 2),
                (11, 1),
                (16, 11),
                (15, 11),
                (15, 13),
                (14, 13),
                (14, 14),
                (13, 14),
                (13, 15),
                (11, 15),
                (11, 16),
                (6, 16),
                (6, 15),
                (4, 15),
                (4, 14),
                (3, 14),
                (3, 13),
                (2, 13),
                (2, 11),
                (1, 11),
            ]
        )
        src.setCoords(coo)
        c = DataArrayInt(
            [
                5,
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
                13,
                14,
                15,
                16,
                17,
                18,
                19,
                20,
                21,
                22,
                23,
                24,
                25,
                26,
                27,
                28,
                29,
                30,
                31,
                32,
                33,
                34,
                35,
                36,
                37,
                38,
                39,
                40,
                41,
                42,
                43,
                44,
                45,
                46,
                47,
                5,
                9,
                10,
                11,
                12,
                13,
                14,
                15,
                16,
                17,
                18,
                19,
                20,
                21,
                22,
                23,
                24,
                25,
                26,
                27,
                28,
                29,
                30,
                31,
                32,
                33,
                34,
                35,
                36,
                37,
                38,
                48,
                49,
                50,
                51,
                52,
                53,
                54,
                55,
                56,
                57,
                58,
                59,
                60,
                61,
                62,
                63,
                64,
                65,
            ]
        )
        cI = DataArrayInt([0, 49, 98])
        src.setConnectivity(c, cI)
        src.checkConsistency()
        tgt = MEDCouplingCMesh("tgt")
        da = DataArrayDouble(18, 1)
        da.iota()
        tgt.setCoords(da, da)
        tgt = tgt.buildUnstructured()
        srcF = MEDCouplingFieldDouble(ON_CELLS, ONE_TIME)
        srcF.setArray(DataArrayDouble([25.0, 50.0]))
        srcF.setMesh(src)
        srcF.setNature(IntensiveConservation)
        remap = MEDCouplingRemapper()
        remap.setIntersectionType(PointLocator)
        remap.prepare(src, tgt, "P0P0")
        tgtF = remap.transferField(srcF, 0.0)
        ids1 = [
            137,
            139,
            141,
            143,
            145,
            147,
            149,
            151,
            154,
            155,
            156,
            157,
            158,
            159,
            160,
            161,
            162,
            163,
            164,
            165,
            166,
            167,
            168,
            171,
            172,
            173,
            174,
            175,
            176,
            177,
            178,
            179,
            180,
            181,
            182,
            183,
            184,
            185,
            189,
            190,
            191,
            192,
            193,
            194,
            195,
            196,
            197,
            198,
            199,
            200,
            201,
            206,
            207,
            208,
            209,
            210,
            211,
            212,
            213,
            214,
            215,
            216,
            217,
            218,
            224,
            225,
            226,
            227,
            228,
            229,
            230,
            231,
            232,
            233,
            234,
            242,
            243,
            244,
            245,
            246,
            247,
            248,
            249,
            250,
            261,
            262,
            263,
            264,
            265,
        ]
        ids2 = [
            23,
            24,
            25,
            26,
            27,
            38,
            39,
            40,
            41,
            42,
            43,
            44,
            45,
            46,
            54,
            55,
            56,
            57,
            58,
            59,
            60,
            61,
            62,
            63,
            64,
            70,
            71,
            72,
            73,
            74,
            75,
            76,
            77,
            78,
            79,
            80,
            81,
            82,
            87,
            88,
            89,
            90,
            91,
            92,
            93,
            94,
            95,
            96,
            97,
            98,
            99,
            103,
            104,
            105,
            106,
            107,
            108,
            109,
            110,
            111,
            112,
            113,
            114,
            115,
            116,
            117,
            120,
            121,
            122,
            123,
            124,
            125,
            126,
            127,
            128,
            129,
            130,
            131,
            132,
            133,
            134,
            138,
            140,
            142,
            144,
            146,
            148,
            150,
        ]
        ids3 = [
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
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
            22,
            28,
            29,
            30,
            31,
            32,
            33,
            34,
            35,
            36,
            37,
            47,
            48,
            49,
            50,
            51,
            52,
            53,
            65,
            66,
            67,
            68,
            69,
            83,
            84,
            85,
            86,
            100,
            101,
            102,
            118,
            119,
            135,
            136,
            152,
            153,
            169,
            170,
            186,
            187,
            188,
            202,
            203,
            204,
            205,
            219,
            220,
            221,
            222,
            223,
            235,
            236,
            237,
            238,
            239,
            240,
            241,
            251,
            252,
            253,
            254,
            255,
            256,
            257,
            258,
            259,
            260,
            266,
            267,
            268,
            269,
            270,
            271,
            272,
            273,
            274,
            275,
            276,
            277,
            278,
            279,
            280,
            281,
            282,
            283,
            284,
            285,
            286,
            287,
            288,
        ]
        a = tgtF.getArray()
        self.assertTrue(a[ids1].isUniform(50.0, 1e-12))
        self.assertTrue(a[ids2].isUniform(25.0, 1e-12))
        self.assertTrue(a[ids3].isUniform(0.0, 1e-12))
        pass

    def testExtrudedOnDiffZLev1(self):
        """Non regression bug : This test is base on P0P0 ExtrudedExtruded. This test checks that if the input meshes are not based on a same plane // OXY the interpolation works"""
        arrX = DataArrayDouble([0, 1])
        arrY = DataArrayDouble([0, 1])
        arrZ = DataArrayDouble([0, 1, 2])
        src = MEDCouplingCMesh()
        src.setCoords(arrX, arrY, arrZ)
        arrX = DataArrayDouble([0.5, 1.5])
        arrY = DataArrayDouble([0.5, 1.5])
        arrZ = DataArrayDouble([0.5, 2])
        trg = MEDCouplingCMesh()
        trg.setCoords(arrX, arrY, arrZ)
        #
        src = MEDCouplingMappedExtrudedMesh(src)
        trg = MEDCouplingMappedExtrudedMesh(trg)
        pt1 = src.getMesh2D().getCoords().getHiddenCppPointer()
        pt2 = trg.getMesh2D().getCoords().getHiddenCppPointer()
        #
        rem = MEDCouplingRemapper()
        rem.prepare(src, trg, "P0P0")
        self.checkMatrix(
            rem.getCrudeMatrix(), [{0: 0.125, 1: 0.25}], src.getNumberOfCells(), 1e-12
        )
        #
        self.assertEqual(src.getMesh2D().getSpaceDimension(), 3)
        self.assertEqual(trg.getMesh2D().getSpaceDimension(), 3)
        self.assertEqual(src.getMesh2D().getCoords().getHiddenCppPointer(), pt1)
        self.assertEqual(trg.getMesh2D().getCoords().getHiddenCppPointer(), pt2)
        #
        rem2 = MEDCouplingRemapper()
        rem2.setIntersectionType(Geometric2D)
        rem2.prepare(src, trg, "P0P0")
        self.checkMatrix(
            rem2.getCrudeMatrix(), [{0: 0.125, 1: 0.25}], src.getNumberOfCells(), 1e-12
        )
        pass

    def testP0P0WithHEXGP12(self):
        """Test that HEXGP12 are correctly remapped (elements with polygonal faces were not properly handled)"""
        # From Astrid, two disjoint hexagonal prisms:
        coo1 = [
            -4.991193077144312,
            8.644999999999998,
            0.0,
            -9.982386154288623,
            6.112246755425186e-16,
            0.0,
            -4.991193077144315,
            -8.644999999999998,
            0.0,
            4.991193077144309,
            -8.645000000000005,
            0.0,
            9.982386154288626,
            1.1651321638577316e-15,
            0.0,
            4.991193077144314,
            8.645,
            0.0,
            -4.991193077144312,
            8.644999999999998,
            7.561799999999991,
            -9.982386154288623,
            6.112246755425186e-16,
            7.561799999999991,
            -4.991193077144315,
            -8.644999999999998,
            7.561799999999991,
            4.991193077144309,
            -8.645000000000005,
            7.561799999999991,
            9.982386154288626,
            1.1651321638577316e-15,
            7.561799999999991,
            4.991193077144314,
            8.645,
            7.561799999999991,
        ]
        coo2 = [
            -4.991193077144313,
            -8.645,
            0.0,
            -9.982386154288626,
            -1.3992140779350848e-15,
            0.0,
            -19.964772308577256,
            0.0,
            0.0,
            -24.95596538572157,
            -8.644999999999998,
            0.0,
            -19.96477230857726,
            -17.289999999999996,
            0.0,
            -9.982386154288626,
            -17.289999999999996,
            0.0,
            -4.991193077144313,
            -8.645,
            5.041200000000004,
            -9.982386154288626,
            -1.3992140779350848e-15,
            5.041200000000004,
            -19.964772308577256,
            0.0,
            5.041200000000004,
            -24.95596538572157,
            -8.644999999999998,
            5.041200000000004,
            -19.96477230857726,
            -17.289999999999996,
            5.041200000000004,
            -9.982386154288626,
            -17.289999999999996,
            5.041200000000004,
        ]
        conn1 = [
            31,
            0,
            5,
            4,
            3,
            2,
            1,
            -1,
            11,
            6,
            7,
            8,
            9,
            10,
            -1,
            1,
            7,
            6,
            0,
            -1,
            2,
            8,
            7,
            1,
            -1,
            3,
            9,
            8,
            2,
            -1,
            4,
            10,
            9,
            3,
            -1,
            5,
            11,
            10,
            4,
            -1,
            0,
            6,
            11,
            5,
        ]
        cI1 = [0, 44]
        conn2 = [
            31,
            0,
            5,
            4,
            3,
            2,
            1,
            -1,
            6,
            7,
            8,
            9,
            10,
            11,
            -1,
            0,
            1,
            7,
            6,
            -1,
            1,
            2,
            8,
            7,
            -1,
            2,
            3,
            9,
            8,
            -1,
            3,
            4,
            10,
            9,
            -1,
            4,
            5,
            11,
            10,
            -1,
            5,
            0,
            6,
            11,
        ]
        cI2 = [0, 44]
        mTgt = MEDCouplingUMesh("target", 3)
        mSrc = MEDCouplingUMesh("src", 3)
        mTgt.setCoords(DataArrayDouble(coo1, len(coo1) // 3, 3))
        mSrc.setCoords(DataArrayDouble(coo2, len(coo2) // 3, 3))
        mTgt.setConnectivity(DataArrayInt(conn1), DataArrayInt(cI1))
        mSrc.setConnectivity(DataArrayInt(conn2), DataArrayInt(cI2))

        # Recognize the HEXGP12:
        mTgt.unPolyze()
        mSrc.unPolyze()

        rmp = MEDCouplingRemapper()
        rmp.setIntersectionType(Triangulation)
        rmp.prepare(mSrc, mTgt, "P0P0")
        mat = rmp.getCrudeMatrix()
        self.assertEqual(len(mat[0]), 0)
        self.assertEqual(len(mat), 1)
        pass

    def testP0P0KillerTet(self):
        """The killer tetrahedron detected by LMEC!"""
        mesh = MEDCouplingUMesh("SupportOf_ECHIA1_Tin", 3)
        #         # was OK:
        #         coo = DataArrayDouble([(-4.50135,1.95352,4.59608),(-4.50409,1.86642,4.54551), (-4.55175,1.92167,4.64844),(-4.58813,1.94795,4.5283)])
        # was KO:
        coo = DataArrayDouble(
            [
                (-4.501352938826142847, 1.953517433537110159, 4.596082552008083688),
                (-4.504092113061189728, 1.866415526007169978, 4.545507396150389567),
                (-4.551750368181751050, 1.921669328035479962, 4.648439577911889664),
                (-4.588131417812300050, 1.947948377683889953, 4.528298931319220344),
            ]
        )
        mesh.setCoords(coo)
        c = DataArrayInt([14, 2, 0, 3, 1])
        cI = DataArrayInt([0, 5])
        mesh.setConnectivity(c, cI)
        mesh_src, mesh_tgt = mesh.deepCopy(), mesh.deepCopy()
        field_src = mesh_src.fillFromAnalytic(ON_CELLS, 1, "1")
        field_src.setNature(IntensiveMaximum)
        rmp = MEDCouplingRemapper()
        rmp.setIntersectionType(Triangulation)
        rmp.prepare(mesh_src, mesh_tgt, "P0P0")
        self.assertEqual(1, len(rmp.getCrudeMatrix()))
        self.assertEqual(1, len(rmp.getCrudeMatrix()[0]))
        pass

    @unittest.skipUnless(
        MEDCouplingHasNumPyBindings() and MEDCouplingHasSciPyBindings(),
        "requires numpy AND scipy AND C++11",
    )
    def testP1P1PL3DSpaceFrom1DTo0D(self):
        from scipy.sparse import csr_matrix
        from numpy import array

        def generateTrg(eps):
            trgArr = DataArrayDouble(
                [
                    (0.5, 0.5, 0.5),
                    (0.2, 0.2, 0.2),
                    (0.9, 0.9, 0.9),
                    (0.7 + eps * sqrt(3), 0.7 - eps * sqrt(3), 0.7),
                ]
            )
            trg = MEDCouplingUMesh("trg", 0)
            trg.setCoords(trgArr)
            trg.allocateCells()
            RenumTrg = [2, 3, 0, 1]
            for rt in RenumTrg:
                trg.insertNextCell(NORM_POINT1, [rt])
            return trg

        srcArr = DataArrayDouble([(0.0, 0.0, 1.0), (0.0, 0.0, 0.0), (1.0, 1.0, 1.0)])
        src = MEDCouplingUMesh("src", 1)
        src.setCoords(srcArr)
        src.allocateCells()
        src.insertNextCell(NORM_SEG2, [1, 2])
        #
        trg = generateTrg(
            1e-7
        )  # trg point 3 of trg cell 1 is NOT closer enough to source edge #1 -> not intercepted
        #
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        self.assertEqual(rem.prepare(src, trg, "P1P1"), 1)
        mat = rem.getCrudeCSRMatrix()
        row = array([2, 2, 0, 0, 1, 1])  # here no ref to point 3 !
        col = array([1, 2, 1, 2, 1, 2])
        data = array([0.1, 0.9, 0.5, 0.5, 0.8, 0.2])
        mExp = csr_matrix((data, (row, col)), shape=(4, 3))
        delta = abs(mExp - mat)
        self.assertAlmostEqual(delta.sum(), 0.0, 14)
        #
        trg = generateTrg(
            1e-14
        )  # trg point 3 of trg cell 1 is closer enough to source edge #1 -> intercepted
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        self.assertEqual(rem.prepare(src, trg, "P1P1"), 1)
        mat = rem.getCrudeCSRMatrix()
        row = array([2, 2, 3, 3, 0, 0, 1, 1])  # here ref to target point 3
        col = array([1, 2, 1, 2, 1, 2, 1, 2])
        data = array([0.1, 0.9, 0.3, 0.7, 0.5, 0.5, 0.8, 0.2])
        mExp2 = csr_matrix((data, (row, col)), shape=(4, 3))
        delta2 = abs(mExp2 - mat)
        self.assertAlmostEqual(delta2.sum(), 0.0, 14)
        pass

    def testSetMatrix1(self):
        """Remapper has now setCrudeMatrix method to reload matrix to skip prepare phase"""
        cooS = DataArrayDouble([1, 1, 7, 1, 7, 2, 1, 2], 4, 2)
        cooT = DataArrayDouble([0, 0, 3, 0, 3, 3, 0, 3, 6, 0, 12, 0, 12, 3, 6, 3], 8, 2)
        ms = MEDCouplingUMesh("source", 2)
        ms.allocateCells(1)
        ms.insertNextCell(NORM_QUAD4, [0, 1, 2, 3])
        ms.setCoords(cooS)
        mt = MEDCouplingUMesh("target", 2)
        mt.allocateCells(2)
        mt.insertNextCell(NORM_QUAD4, [0, 1, 2, 3])
        mt.insertNextCell(NORM_QUAD4, [4, 5, 6, 7])
        mt.setCoords(cooT)
        rem = MEDCouplingRemapper()
        self.assertEqual(rem.prepare(ms, mt, "P0P0"), 1)  # [{0: 2.0}, {0: 1.0}]
        fs = MEDCouplingFieldDouble(ON_CELLS)
        fs.setMesh(ms)
        fs.setArray(DataArrayDouble([10]))
        fs.checkConsistencyLight()
        #
        fs.setNature(ExtensiveConservation)
        self.assertTrue(
            rem.transferField(fs, 1e300)
            .getArray()
            .isEqual(DataArrayDouble([20.0 / 3, 10.0 / 3.0]), 1e-12)
        )  # sum is equal to 10. First value is twice than second value
        #
        fs.setNature(ExtensiveMaximum)
        self.assertTrue(
            rem.transferField(fs, 1e300)
            .getArray()
            .isEqual(DataArrayDouble([20.0 / 6.0, 10.0 / 6.0]), 1e-12)
        )  # sum is equal to 5 (10/2. because only half part on input cell is intercepted by the target cells). First value is twice than second value
        #
        fs.setNature(IntensiveConservation)
        self.assertTrue(
            rem.transferField(fs, 1e300)
            .getArray()
            .isEqual(DataArrayDouble([2.0 / 9.0 * 10.0, 1.0 / 18.0 * 10.0]), 1e-12)
        )  #
        #
        fs.setNature(IntensiveMaximum)
        self.assertTrue(
            rem.transferField(fs, 1e300)
            .getArray()
            .isEqual(DataArrayDouble([10.0, 10.0]), 1e-12)
        )  #
        ####
        rem2 = MEDCouplingRemapper()
        rem2.setCrudeMatrix(ms, mt, "P0P0", rem.getCrudeMatrix())
        fs.setNature(ExtensiveConservation)
        self.assertTrue(
            rem2.transferField(fs, 1e300)
            .getArray()
            .isEqual(DataArrayDouble([20.0 / 3, 10.0 / 3.0]), 1e-12)
        )
        #
        fs.setNature(ExtensiveMaximum)
        self.assertTrue(
            rem2.transferField(fs, 1e300)
            .getArray()
            .isEqual(DataArrayDouble([20.0 / 6.0, 10.0 / 6.0]), 1e-12)
        )
        #
        fs.setNature(IntensiveConservation)
        self.assertTrue(
            rem2.transferField(fs, 1e300)
            .getArray()
            .isEqual(DataArrayDouble([2.0 / 9.0 * 10.0, 1.0 / 18.0 * 10.0]), 1e-12)
        )
        #
        fs.setNature(IntensiveMaximum)
        self.assertTrue(
            rem2.transferField(fs, 1e300)
            .getArray()
            .isEqual(DataArrayDouble([10.0, 10.0]), 1e-12)
        )
        #
        srcFt = MEDCouplingFieldTemplate.New(ON_CELLS)
        trgFt = MEDCouplingFieldTemplate.New(ON_CELLS)
        srcFt.setMesh(ms)
        trgFt.setMesh(mt)
        rem3 = MEDCouplingRemapper()
        rem3.setCrudeMatrixEx(srcFt, trgFt, rem.getCrudeMatrix())
        fs.setNature(ExtensiveConservation)
        self.assertTrue(
            rem3.transferField(fs, 1e300)
            .getArray()
            .isEqual(DataArrayDouble([20.0 / 3, 10.0 / 3.0]), 1e-12)
        )
        pass

    @unittest.skipUnless(
        MEDCouplingHasNumPyBindings() and MEDCouplingHasSciPyBindings(),
        "requires numpy AND scipy",
    )
    def testSetMatrix2(self):
        """Remapper has now setCrudeMatrix method to reload matrix to skip prepare phase. Same as testSetMatrix1 but with CSR scipy matrix"""
        arrx_s = DataArrayDouble(6)
        arrx_s.iota()
        arry_s = DataArrayDouble(6)
        arry_s.iota()
        ms = MEDCouplingCMesh()
        ms.setCoords(arrx_s, arry_s)
        ms = ms.buildUnstructured()
        #
        arrx_t = DataArrayDouble([2.5, 4.5, 5.5])
        arry_t = DataArrayDouble([2.5, 3.5, 5.5])
        mt = MEDCouplingCMesh()
        mt.setCoords(arrx_t, arry_t)
        mt = mt.buildUnstructured()
        #
        rem = MEDCouplingRemapper()
        self.assertEqual(rem.prepare(ms, mt, "P0P0"), 1)
        #
        fs = MEDCouplingFieldDouble(ON_CELLS)
        fs.setMesh(ms)
        arr = DataArrayDouble(25)
        arr.iota()
        fs.setArray(arr)
        fs.checkConsistencyLight()
        #
        fs.setNature(ExtensiveConservation)
        self.assertTrue(
            rem.transferField(fs, 1e300)
            .getArray()
            .isEqual(DataArrayDouble([54.25, 11.75, 79.25, 16.75]), 1e-12)
        )
        mat = rem.getCrudeCSRMatrix()
        rem2 = MEDCouplingRemapper()
        rem2.setCrudeMatrix(ms, mt, "P0P0", mat)
        self.assertTrue(
            rem2.transferField(fs, 1e300)
            .getArray()
            .isEqual(DataArrayDouble([54.25, 11.75, 79.25, 16.75]), 1e-12)
        )
        pass

    def testSmallTetraCell(self):
        """This test is a non regression test. When using tetra/tetra P0P0 interpolation on very small cells the
        3x3 matrix in the TetraAffine contains very small values and so the determinant is small (cubic).
        So the tetra was detected as flat. Now the infinite norm of matrix is considered to establish if matrix is inversible or not."""
        coords = [
            (-0.019866666666666668, 0.02, 0.002),
            (-0.020000073463967143, 0.019999926535763005, 0.0018666666666666673),
            (-0.020000073463967143, 0.019999926535763005, 0.002),
            (-0.020000072974206463, 0.019866593202430387, 0.002),
        ]
        m = MEDCouplingUMesh("mesh", 3)
        m.allocateCells()
        m.insertNextCell(NORM_TETRA4, [0, 1, 2, 3])
        m.setCoords(DataArrayDouble(coords))
        rem = MEDCouplingRemapper()
        rem.setPrecision(1e-12)
        rem.prepare(m, m, "P0P0")
        mat = rem.getCrudeMatrix()
        self.assertTrue(len(mat) == 1)
        self.assertTrue(len(mat[0]) == 1)
        self.assertTrue(list(mat[0].keys()) == [0])
        res = list(mat[0].values())[0]
        ref = float(m.getMeasureField(True).getArray())
        self.assertTrue(abs(res - ref) / ref < 1e-12)
        pass

    def test3D0DPointLocator(self):
        """
        For pointlocator fans, Remapper support following intersection
        IntersectionType == PointLocator
        - source == 3D
        - target == 0D
        """
        src = MEDCouplingUMesh("src", 3)
        src.allocateCells()
        src.setCoords(DataArrayDouble([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]))
        src.insertNextCell(NORM_TETRA4, [0, 1, 2, 3])
        trg = MEDCouplingUMesh.Build0DMeshFromCoords(
            DataArrayDouble([(0.4, 0.3, 0.07)])
        )
        # P1P1
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(src, trg, "P1P1")
        self.checkMatrix(
            rem.getCrudeMatrix(),
            [{0: 0.23, 1: 0.4, 2: 0.3, 3: 0.07}],
            src.getNumberOfNodes(),
            1e-12,
        )
        # P1P0
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(src, trg, "P1P0")
        self.checkMatrix(
            rem.getCrudeMatrix(),
            [{0: 0.23, 1: 0.4, 2: 0.3, 3: 0.07}],
            src.getNumberOfNodes(),
            1e-12,
        )
        # P0P1
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(src, trg, "P0P1")
        self.checkMatrix(
            rem.getCrudeMatrix(), [{0: 1.0}], src.getNumberOfCells(), 1e-12
        )
        # P0P0
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(src, trg, "P0P0")
        self.checkMatrix(
            rem.getCrudeMatrix(), [{0: 1.0}], src.getNumberOfCells(), 1e-12
        )
        pass

    def test2D0DPointLocator(self):
        """
        For pointlocator fans, Remapper support following intersection
        IntersectionType == PointLocator
        - source == 2D
        - target == 0D
        """
        src = MEDCouplingUMesh("src", 2)
        src.allocateCells()
        src.setCoords(DataArrayDouble([(0, 0), (1, 0), (0, 1)]))
        src.insertNextCell(NORM_TRI3, [0, 1, 2])
        trg = MEDCouplingUMesh.Build0DMeshFromCoords(DataArrayDouble([(0.4, 0.3)]))
        # P1P1
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(src, trg, "P1P1")
        self.checkMatrix(
            rem.getCrudeMatrix(),
            [{0: 0.3, 1: 0.4, 2: 0.3}],
            src.getNumberOfNodes(),
            1e-12,
        )
        # P1P0
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(src, trg, "P1P0")
        self.checkMatrix(
            rem.getCrudeMatrix(),
            [{0: 0.3, 1: 0.4, 2: 0.3}],
            src.getNumberOfNodes(),
            1e-12,
        )
        # P0P1
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(src, trg, "P0P1")
        self.checkMatrix(
            rem.getCrudeMatrix(), [{0: 1.0}], src.getNumberOfNodes(), 1e-12
        )
        # P0P0
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(src, trg, "P0P0")
        self.checkMatrix(
            rem.getCrudeMatrix(), [{0: 1.0}], src.getNumberOfNodes(), 1e-12
        )
        pass

    def test1D0DPointLocator(self):
        """
        For pointlocator fans, Remapper support following intersection
        IntersectionType == PointLocator
        - source == 1D
        - target == 0D
        """
        # P1P1 - 0
        src = MEDCouplingUMesh("src", 1)
        src.allocateCells()
        src.setCoords(DataArrayDouble([0, 1]))
        src.insertNextCell(NORM_SEG2, [0, 1])
        trg = MEDCouplingUMesh.Build0DMeshFromCoords(DataArrayDouble([0.4]))
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(src, trg, "P1P1")
        self.checkMatrix(
            rem.getCrudeMatrix(), [{0: 0.6, 1: 0.4}], src.getNumberOfNodes(), 1e-12
        )
        # P1P1 - 1
        src = MEDCouplingUMesh("src", 1)
        src.allocateCells()
        src.setCoords(DataArrayDouble([0, 1]))
        src.insertNextCell(NORM_SEG2, [1, 0])  # permutation
        trg = MEDCouplingUMesh.Build0DMeshFromCoords(DataArrayDouble([0.4]))
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(src, trg, "P1P1")
        self.checkMatrix(
            rem.getCrudeMatrix(), [{0: 0.6, 1: 0.4}], src.getNumberOfNodes(), 1e-12
        )
        # P1P1 - 2
        src = MEDCouplingUMesh("src", 1)
        src.allocateCells()
        src.setCoords(DataArrayDouble([1, 0]))
        src.insertNextCell(NORM_SEG2, [0, 1])
        trg = MEDCouplingUMesh.Build0DMeshFromCoords(DataArrayDouble([0.4]))
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(src, trg, "P1P1")
        self.checkMatrix(
            rem.getCrudeMatrix(), [{0: 0.4, 1: 0.6}], src.getNumberOfNodes(), 1e-12
        )
        # P1P1 - 3 - 2DCurve
        src = MEDCouplingUMesh("src", 1)
        src.allocateCells()
        src.setCoords(DataArrayDouble([0, 1]))
        src.insertNextCell(NORM_SEG2, [0, 1])
        trg = MEDCouplingUMesh.Build0DMeshFromCoords(DataArrayDouble([0.4]))
        src.changeSpaceDimension(2)
        trg.changeSpaceDimension(2)
        src.rotate([-1.0, -1.0], 1.2)
        trg.rotate([-1.0, -1.0], 1.2)
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(src, trg, "P1P1")
        self.checkMatrix(
            rem.getCrudeMatrix(), [{0: 0.6, 1: 0.4}], src.getNumberOfNodes(), 1e-12
        )
        # P1P1 - 4
        src = MEDCouplingUMesh("src", 1)
        src.allocateCells()
        src.setCoords(DataArrayDouble([1.1, 7.6, 2.3, 5.4]))
        src.insertNextCell(NORM_SEG2, [0, 2])
        src.insertNextCell(NORM_SEG2, [2, 3])
        src.insertNextCell(NORM_SEG2, [3, 1])
        for eps in [0, 1e-13, -1e-13]:
            trg = MEDCouplingUMesh.Build0DMeshFromCoords(
                DataArrayDouble([0.4, 2.3 + eps, 4.0, 7.0])
            )
            rem = MEDCouplingRemapper()
            rem.setIntersectionType(PointLocator)
            rem.prepare(src, trg, "P1P1")
            rem.nullifiedTinyCoeffInCrudeMatrixAbs(1e-12)
            self.checkMatrix(
                rem.getCrudeMatrix(),
                [
                    {},
                    {2: 2.0},
                    {2: 0.4516129032258065, 3: 0.5483870967741935},
                    {1: 0.7272727272727273, 3: 0.27272727272727265},
                ],
                src.getNumberOfNodes(),
                1e-12,
            )
        # P1P1 - 5 - descending order of coords in source mesh
        src = MEDCouplingUMesh("src", 1)
        src.allocateCells()
        src.setCoords(DataArrayDouble([3.0, 1.0]))
        src.insertNextCell(NORM_SEG2, [0, 1])
        trg = MEDCouplingUMesh.Build0DMeshFromCoords(DataArrayDouble([2.3]))
        rem = MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(src, trg, "P1P1")
        self.checkMatrix(
            rem.getCrudeMatrix(), [{0: 0.65, 1: 0.35}], src.getNumberOfNodes(), 1e-12
        )
        pass

    @unittest.skipUnless(
        MEDCouplingHasNumPyBindings() and MEDCouplingHasSciPyBindings(),
        "requires numpy AND scipy",
    )
    def testRemToCSRMatrix(self):
        import scipy

        mPy = [{0: 1.0, 1: 3.0, 3: 7.0, 6: 10.0}, {1: 12.0, 2: 23.0}]
        m = MEDCouplingRemapper.ToCSRMatrix(mPy, 8)
        self.assertTrue(isinstance(m, scipy.sparse.csr_matrix))
        self.assertEqual(m.getnnz(), 6)
        self.assertAlmostEqual(m[0, 0], 1.0, 12)
        self.assertAlmostEqual(m[0, 1], 3.0, 12)
        self.assertAlmostEqual(m[0, 3], 7.0, 12)
        self.assertAlmostEqual(m[0, 6], 10.0, 12)
        self.assertAlmostEqual(m[1, 1], 12.0, 12)
        self.assertAlmostEqual(m[1, 2], 23.0, 12)
        self.assertEqual(m.shape, (2, 8))

    def test_Interpolation2D3D_bbox_adjustment_1(self):
        """Interpolation 2D <-> 3D was not using bounding box adjustment.
        In case of a 2D mesh perfectly aligned with the axis, the bounding box intersection was not working properly (flat bounding box).
        """
        ## Source
        meshS = MEDCouplingUMesh("SupportOf_TEMPERATURE_OUT", 2)
        coo = cooS = DataArrayDouble(
            [
                (
                    -0.00074999999999877595,
                    0.00000000000000000000,
                    0.00032540000000000005,
                ),
                (
                    -0.00049999999999755579,
                    0.00025000000000140708,
                    0.00032540000000000005,
                ),
                (
                    -0.00049999999999755600,
                    0.00000000000000000000,
                    0.00032540000000000005,
                ),
                (
                    -0.00100000000000000002,
                    0.00000000000000000000,
                    0.00032540000000000005,
                ),
                (
                    -0.00100000000000000002,
                    0.00025000000000000543,
                    0.00032540000000000005,
                ),
                (
                    -0.00075651925565617829,
                    0.00034416831541328637,
                    0.00032540000000000005,
                ),
            ]
        )  # the extra 5e-20 on Z is the true culprit :-)
        meshS.setCoords(coo)
        c = DataArrayInt([3, 0, 1, 2, 3, 0, 3, 4, 3, 5, 0, 4, 3, 5, 1, 0])
        cI = DataArrayInt([0, 4, 8, 12, 16])
        meshS.setConnectivity(c, cI)
        meshS.checkConsistency()
        ## Target
        meshT = MEDCouplingUMesh("IJK_mesh", 3)
        coo = DataArrayDouble(
            [
                (-0.001, 0, 0.000303602),
                (-0.0009, 0, 0.000303602),
                (-0.0008, 0, 0.000303602),
                (-0.0007, 0, 0.000303602),
                (-0.0006, 0, 0.000303602),
                (-0.0005, 0, 0.000303602),
                (-0.001, 0.0005, 0.000303602),
                (-0.0009, 0.0005, 0.000303602),
                (-0.0008, 0.0005, 0.000303602),
                (-0.0007, 0.0005, 0.000303602),
                (-0.0006, 0.0005, 0.000303602),
                (-0.0005, 0.0005, 0.000303602),
                (-0.001, 0, 0.0003254),
                (-0.0009, 0, 0.0003254),
                (-0.0008, 0, 0.0003254),
                (-0.0007, 0, 0.0003254),
                (-0.0006, 0, 0.0003254),
                (-0.0005, 0, 0.0003254),
                (-0.001, 0.0005, 0.0003254),
                (-0.0009, 0.0005, 0.0003254),
                (-0.0008, 0.0005, 0.0003254),
                (-0.0007, 0.0005, 0.0003254),
                (-0.0006, 0.0005, 0.0003254),
                (-0.0005, 0.0005, 0.0003254),
            ]
        )
        meshT.setCoords(coo)
        c = DataArrayInt(
            [
                18,
                1,
                0,
                6,
                7,
                13,
                12,
                18,
                19,
                18,
                2,
                1,
                7,
                8,
                14,
                13,
                19,
                20,
                18,
                3,
                2,
                8,
                9,
                15,
                14,
                20,
                21,
                18,
                4,
                3,
                9,
                10,
                16,
                15,
                21,
                22,
                18,
                5,
                4,
                10,
                11,
                17,
                16,
                22,
                23,
            ]
        )
        cI = DataArrayInt([0, 9, 18, 27, 36, 45])
        meshT.setConnectivity(c, cI)
        meshT.checkConsistency()
        ## Dummy field
        fldSrc = MEDCouplingFieldDouble(ON_CELLS, ONE_TIME)
        fldSrc.setMesh(meshS)
        da = DataArrayDouble(meshS.getNumberOfCells())
        da[:] = 50.0
        fldSrc.setArray(da)
        remap = MEDCouplingRemapper()
        # remap.setBoundingBoxAdjustmentAbs(1.0e-5)  # was not taken into account for 2D/3D - but we don't even need it! Default value is OK.
        remap.prepare(meshS, meshT, "P0P0")
        fldSrc.setNature(IntensiveMaximum)
        fldTgt = remap.transferField(fldSrc, -1.0)
        self.assertTrue(fldTgt.getArray().isUniform(50.0, 1e-12))

    def testGrandyBug1(self):
        """
        Non regression test relative to test tuleap26461
        """
        rem = MEDCouplingRemapper()
        src_final = MEDCouplingUMesh("src_final", 3)
        src_final.setCoords(
            DataArrayDouble(
                [
                    0.74763179385813627,
                    2.0528797000000716,
                    0.42830000000000013,
                    0.77426950622837643,
                    2.0528797000000001,
                    0.40000000000000036,
                    0.77426950622837643,
                    2.0528797000000001,
                    0.42830000000000013,
                    0.77426950622830537,
                    2.0262419876297604,
                    0.42830000000000013,
                ],
                4,
                3,
            )
        )
        src_final.allocateCells()
        src_final.insertNextCell(NORM_TETRA4, [0, 3, 2, 1])
        trg_final = MEDCouplingUMesh("trg_final", 3)
        trg_final.setCoords(
            DataArrayDouble(
                [
                    0.81034725000000007,
                    1.9988565499999984,
                    0.40000000000000002,
                    0.75632410000000005,
                    2.0528796999999983,
                    0.40000000000000002,
                    0.75632410000000005,
                    1.9988565499999984,
                    0.41800000000000004,
                    0.81034725000000007,
                    2.0528796999999983,
                    0.41800000000000004,
                ],
                4,
                3,
            )
        )
        trg_final.allocateCells()
        trg_final.insertNextCell(NORM_TETRA4, [0, 2, 1, 3])

        ref_values = [  # ref values coming from geom2medcoupling.py
            (0.0, 1.671615506097834e-08),
            (1e-12, 1.671615506712106e-08),
            (1e-11, 1.671615512239666e-08),
            (1e-10, 1.6716155675164925e-08),
            (1e-9, 1.671616120285316e-08),
            (1e-8, 1.6716216479802182e-08),
            (1e-7, 1.671676925650806e-08),
            (1e-6, 1.672014459316238e-08),
            (1e-5, 1.6805275475457618e-08),
            (1e-4, 1.7608769838220544e-08),
            (1e-3, 2.5791583779126835e-08),
        ]

        for ty, ref_value in ref_values:
            trg_final2 = trg_final.deepCopy()
            trg_final2.translate([0, ty, 0])
            rem.setPrecision(1e-12)
            rem.prepare(src_final, trg_final2, "P0P0")
            mat_mc = rem.getCrudeMatrix()
            csr_new = MEDCouplingRemapper.ToCSRMatrix(
                mat_mc, src_final.getNumberOfCells()
            )
            delta = abs(csr_new[0, 0] - ref_value) / ref_value
            self.assertTrue(delta < 1e-3)

    def testFEFE_QUAD(self):
        """
        [EDF31187] : Test to stress localisation of target points into a source mesh using standard reference FE elements.

        Test for QUAD4, QUAD8 and QUAD9 with non-planar faces
        """

        coo = DataArrayDouble(
            [
                (9.0, 18.0),
                (11.0, 18.0),
                (11.0, 22.0),
                (9.0, 22.0),
                (10, 18.1),
                (11.0, 20.0),
                (10, 21.9),
                (9.0, 20.0),
                (10.0, 20.0),
            ]
        )

        eps = 1e-8

        ref_val0_all = {
            NORM_QUAD4: DataArrayDouble(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    1.5,
                ]
            ),
            NORM_QUAD8: DataArrayDouble(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    7.0,
                    9.5,
                ]
            ),
            NORM_QUAD9: DataArrayDouble(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    7.0,
                    8.0,
                ]
            ),
        }

        inPts_all = {
            NORM_QUAD4: DataArrayDouble(
                [
                    (9.0, 18.0),
                    (11.0, 18.0),
                    (11.0, 22.0),
                    (9.0, 22.0),
                    (10.0, 20.0),
                ]
            ),
            NORM_QUAD8: DataArrayDouble(
                [
                    (9.0, 18.0),
                    (11.0, 18.0),
                    (11.0, 22.0),
                    (9.0, 22.0),
                    (10, 18.1),
                    (11.0, 20.0),
                    (10, 21.9),
                    (9.0, 20.0),
                    (10.0, 20.0),
                ]
            ),
            NORM_QUAD9: DataArrayDouble(
                [
                    (9.0, 18.0),
                    (11.0, 18.0),
                    (11.0, 22.0),
                    (9.0, 22.0),
                    (10, 18.1),
                    (11.0, 20.0),
                    (10, 21.9),
                    (9.0, 20.0),
                    (10.0, 20.0),
                ]
            ),
        }

        outPts_all = {
            NORM_QUAD4: DataArrayDouble(
                [
                    (10.0, 17.99999),
                    (0.0, 0.0),
                    (8.999999, 20.0),
                ]
            ),
            NORM_QUAD8: DataArrayDouble(
                [
                    (10.0, 17.99999),
                    (10.0, 18.09999),
                    (0.0, 0.0),
                    (8.999999, 20.0),
                ]
            ),
            NORM_QUAD9: DataArrayDouble(
                [
                    (10.0, 17.99999),
                    (10.0, 18.09999),
                    (0.0, 0.0),
                    (8.999999, 20.0),
                ]
            ),
        }

        for gt in [NORM_QUAD4, NORM_QUAD8, NORM_QUAD9]:
            nbPtsInCell = MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(gt)
            m = MEDCouplingUMesh("mesh", 2)
            m.setCoords(coo)
            m.allocateCells()
            m.insertNextCell(gt, list(range(nbPtsInCell)))
            m.finishInsertingCells()

            srcField = MEDCouplingFieldDouble(ON_NODES_FE)
            srcField.setMesh(m)
            arr = DataArrayDouble(nbPtsInCell)
            arr.iota()
            srcField.setArray(arr)

            inPts = inPts_all[gt]
            outPts = outPts_all[gt]
            ref_val0 = ref_val0_all[gt]
            nbInPts = len(inPts)

            for i in range(nbInPts):
                self.assertTrue(
                    abs(srcField.getValueOn(inPts[i])[0] - ref_val0[i]) < eps
                )

            self.assertTrue(srcField.getValueOnMulti(inPts).isEqual(ref_val0, eps))

            srcFt = MEDCouplingFieldTemplate(srcField)
            trgFt = MEDCouplingFieldTemplate(ON_NODES_FE)
            trgMesh = MEDCouplingUMesh.Build0DMeshFromCoords(inPts)
            trgFt.setMesh(trgMesh)
            rem = MEDCouplingRemapper()
            rem.setIntersectionType(PointLocator)
            rem.prepareEx(srcFt, trgFt)
            # scan content of matrix computed by remapper
            mat = rem.getCrudeMatrix()
            self.assertEqual(len(mat), nbInPts)
            for irow, row in enumerate(mat):
                self.assertTrue(abs(sum([y for x, y in row.items()]) - 1.0) < eps)

            # ask for MEDCouplingFieldDiscretizationOnNodesFE instance to compute coordination into ref element
            sd = srcField.getDiscretization()
            coosInRefElemFoundByNewton = sd.getCooInRefElement(
                srcField.getMesh(), inPts
            )

            for zePt, cooInRefElemFoundByNewton in zip(
                inPts, coosInRefElemFoundByNewton
            ):
                # now check by performing refCoo -> realCoo
                ftest = MEDCouplingFieldDouble(ON_GAUSS_PT)
                ftest.setMesh(srcField.getMesh())
                ftest.setGaussLocalizationOnType(
                    gt,
                    sum(
                        [
                            list(elt)
                            for elt in MEDCouplingGaussLocalization.GetDefaultReferenceCoordinatesOf(
                                gt
                            ).getValuesAsTuple()
                        ],
                        [],
                    ),
                    list(cooInRefElemFoundByNewton),
                    [1],
                )
                self.assertTrue(
                    float((ftest.getLocalizationOfDiscr() - zePt).magnitude()) < 1e-4
                )

            # testing with outside point
            for zePt in outPts:
                # safer than
                # self.assertRaises(InterpKernelException,sd.getCooInRefElement,srcField.getMesh(),zePt.buildDADouble())
                try:
                    sd.getCooInRefElement(srcField.getMesh(), zePt.buildDADouble())
                except InterpKernelException as e:
                    self.assertTrue("fail to locate point" in e.what())
                else:
                    self.assertTrue(False, f"{zePt}")

    def testFEFE_HEXA(self):
        """
        [EDF31887] : Test to stress localisation of target points into a source mesh using standard reference FE elements.

        Test for HEXA8, HEXA20 and HEXA27 with non-planar faces
        """

        coo = DataArrayDouble(
            [
                (9.0, 18.0, 27.0),
                (9.0, 22.0, 27.0),
                (11.0, 22.0, 27.0),
                (11.0, 18.0, 27.0),
                (9.0, 18.0, 33.0),
                (9.0, 22.0, 33.0),
                (11.0, 22.0, 33.0),
                (11.0, 18.0, 33.0),
                (8.8, 20.0, 26.4),
                (10.0, 21.6, 27.6),
                (11.2, 20.0, 26.4),
                (10.0, 18.4, 27.6),
                (8.8, 20.0, 33.6),
                (10.0, 21.6, 32.4),
                (11.2, 20.0, 33.6),
                (10.0, 18.4, 32.4),
                (8.8, 17.6, 30.0),
                (9.2, 21.6, 30.0),
                (11.2, 22.4, 30.0),
                (10.8, 18.4, 30.0),
                (10.0, 20.0, 26.4),
                (9.2, 20.0, 30.0),
                (10.0, 22.4, 30.0),
                (10.8, 20.0, 30.0),
                (10.0, 17.6, 30.0),
                (10.0, 20.0, 32.4),
                (10.0, 20.0, 30.0),
            ]
        )

        eps = 1e-8

        ref_val0_all = {
            NORM_HEXA8: DataArrayDouble(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    7.0,
                    0.39499999999999968,
                    1.2050000000000001,
                    2.1950000000000003,
                    3.0050000000000008,
                    3.9950000000000019,
                    4.8050000000000015,
                    5.7950000000000026,
                    6.6050000000000022,
                    2.0999999999999992,
                    2.0999999999999992,
                    4.9000000000000004,
                    4.9000000000000012,
                    3.0449999999999999,
                    4.4449999999999985,
                    2.8000000000000007,
                    4.1999999999999993,
                    4.9000000000000004,
                    3.5,
                ]
            ),
            NORM_HEXA20: DataArrayDouble(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    7.0,
                    8.0,
                    9.0,
                    10.0,
                    11.0,
                    12.0,
                    13.0,
                    14.0,
                    15.0,
                    16.0,
                    17.0,
                    18.0,
                    19.0,
                    2.9837271818496278,
                    3.5626912044542562,
                    4.3303322656262804,
                    5.6397881200484505,
                    7.723465224770071,
                    7.9394713811744184,
                    11.719570381471632,
                    9.3911004685421293,
                    10.049844847063252,
                    9.8511656389488813,
                    11.064826449849123,
                    11.719570381471632,
                    13.032378656477301,
                    18.051315791314238,
                    18.906074591293937,
                    19.95383395247174,
                    20.856771727594506,
                    21.928765459430675,
                    25.417967079131408,
                    33.5,
                ]
            ),
            NORM_HEXA27: DataArrayDouble(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    7.0,
                    8.0,
                    9.0,
                    10.0,
                    11.0,
                    12.0,
                    13.0,
                    14.0,
                    15.0,
                    16.0,
                    17.0,
                    18.0,
                    19.0,
                    20.0,
                    21.0,
                    22.0,
                    23.0,
                    24.0,
                    25.0,
                    26.0,
                    6.5782430384766535,
                    5.5057603626436311,
                    7.8099625624768398,
                    7.643290943690746,
                    9.758765408487054,
                    9.06408508454036,
                    11.003779543627997,
                    11.205026363340515,
                    10.56416007071563,
                    18.44359561721225,
                    12.499588353132655,
                    19.85351355074463,
                    14.186041573114885,
                    19.339743214851023,
                    16.084629460041207,
                    20.892007336402663,
                    17.269258227200577,
                    19.549962126638338,
                    19.039562190136063,
                    21.648928068870756,
                    20.667409503475078,
                    22.062499999998867,
                    22.562500000009678,
                    23.812499999505995,
                    24.395833333387696,
                    25.62468592706991,
                    26.0,
                ]
            ),
        }

        inPts_all = {
            NORM_HEXA8: DataArrayDouble(
                [
                    (9.0, 18.0, 27.0),
                    (9.0, 22.0, 27.0),
                    (11.0, 22.0, 27.0),
                    (11.0, 18.0, 27.0),
                    (9.0, 18.0, 33.0),
                    (9.0, 22.0, 33.0),
                    (11.0, 22.0, 33.0),
                    (11.0, 18.0, 33.0),
                    (9.1, 18.2, 27.3),
                    (9.1, 21.8, 27.3),
                    (10.9, 21.8, 27.3),
                    (10.9, 18.2, 27.3),
                    (9.1, 18.2, 32.7),
                    (9.1, 21.8, 32.7),
                    (10.9, 21.8, 32.7),
                    (10.9, 18.2, 32.7),
                    (10.0, 21.4, 27.9),
                    (10.0, 18.6, 27.9),
                    (10.0, 21.4, 32.1),
                    (10.0, 18.6, 32.1),
                    (9.3, 21.4, 30.0),
                    (10.7, 18.6, 30.0),
                    (9.3, 20.0, 30.0),
                    (10.7, 20.0, 30.0),
                    (10.0, 20.0, 32.1),
                    (10.0, 20.0, 30.0),
                ]
            ),
            NORM_HEXA20: DataArrayDouble(
                [
                    (9.0, 18.0, 27.0),
                    (9.0, 22.0, 27.0),
                    (11.0, 22.0, 27.0),
                    (11.0, 18.0, 27.0),
                    (9.0, 18.0, 33.0),
                    (9.0, 22.0, 33.0),
                    (11.0, 22.0, 33.0),
                    (11.0, 18.0, 33.0),
                    (8.8, 20.0, 26.4),
                    (10.0, 21.6, 27.6),
                    (11.2, 20.0, 26.4),
                    (10.0, 18.4, 27.6),
                    (8.8, 20.0, 33.6),
                    (10.0, 21.6, 32.4),
                    (11.2, 20.0, 33.6),
                    (10.0, 18.4, 32.4),
                    (8.8, 17.6, 30.0),
                    (9.2, 21.6, 30.0),
                    (11.2, 22.4, 30.0),
                    (10.8, 18.4, 30.0),
                    (9.1, 18.2, 27.3),
                    (9.1, 21.8, 27.3),
                    (10.9, 21.8, 27.3),
                    (10.9, 18.2, 27.3),
                    (9.1, 18.2, 32.7),
                    (9.1, 21.8, 32.7),
                    (11.1, 20.0, 26.7),
                    (10.9, 21.8, 32.7),
                    (10.9, 18.2, 32.7),
                    (8.9, 20.0, 26.7),
                    (10.0, 21.4, 27.9),
                    (11.1, 20.0, 26.7),
                    (10.0, 18.6, 27.9),
                    (8.9, 17.8, 30.0),
                    (9.3, 21.4, 30.0),
                    (11.1, 22.2, 30.0),
                    (10.7, 18.6, 30.0),
                    (9.3, 20.0, 30.0),
                    (10.7, 20.0, 30.0),
                    (10.0, 20.0, 30.0),
                ]
            ),
            NORM_HEXA27: DataArrayDouble(
                [
                    (9.0, 18.0, 27.0),
                    (9.0, 22.0, 27.0),
                    (11.0, 22.0, 27.0),
                    (11.0, 18.0, 27.0),
                    (9.0, 18.0, 33.0),
                    (9.0, 22.0, 33.0),
                    (11.0, 22.0, 33.0),
                    (11.0, 18.0, 33.0),
                    (8.8, 20.0, 26.4),
                    (10.0, 21.6, 27.6),
                    (11.2, 20.0, 26.4),
                    (10.0, 18.4, 27.6),
                    (8.8, 20.0, 33.6),
                    (10.0, 21.6, 32.4),
                    (11.2, 20.0, 33.6),
                    (10.0, 18.4, 32.4),
                    (8.8, 17.6, 30.0),
                    (9.2, 21.6, 30.0),
                    (11.2, 22.4, 30.0),
                    (10.8, 18.4, 30.0),
                    (10.0, 20.0, 26.4),
                    (9.2, 20.0, 30.0),
                    (10.0, 22.4, 30.0),
                    (10.8, 20.0, 30.0),
                    (10.0, 17.6, 30.0),
                    (10.0, 20.0, 32.4),
                    (10.0, 20.0, 30.0),
                    (9.1, 18.2, 27.3),
                    (9.1, 21.8, 27.3),
                    (10.9, 21.8, 27.3),
                    (10.9, 18.2, 27.3),
                    (9.1, 18.2, 32.7),
                    (9.1, 21.8, 32.7),
                    (10.9, 21.8, 32.7),
                    (10.9, 18.2, 32.7),
                    (8.9, 20.0, 26.7),
                    (10.0, 21.4, 27.9),
                    (11.1, 20.0, 26.7),
                    (10.0, 18.6, 27.9),
                    (8.9, 20.0, 33.3),
                    (10.0, 21.4, 32.1),
                    (11.1, 20.0, 33.3),
                    (10.0, 18.6, 32.1),
                    (8.9, 17.8, 30.0),
                    (9.3, 21.4, 30.0),
                    (11.1, 22.2, 30.0),
                    (10.7, 18.6, 30.0),
                    (10.0, 20.0, 26.7),
                    (9.3, 20.0, 30.0),
                    (10.0, 22.2, 30.0),
                    (10.7, 20.0, 30.0),
                    (10.0, 17.8, 30.0),
                    (10.0, 20.0, 32.1),
                    (10.0, 20.0, 30.0),
                ]
            ),
        }

        outPts_all = {
            NORM_HEXA8: DataArrayDouble(
                [
                    (8.9, 20.0, 26.7),
                    (8.9, 20.0, 33.3),
                    (8.9, 17.8, 30.0),
                    (11.1, 22.2, 30.0),
                    (11.1, 20.0, 33.3),
                    (10.0, 22.2, 30.0),
                    (8.9, 17.8, 26.7),
                    (8.9, 22.2, 26.7),
                    (11.1, 22.2, 26.7),
                    (11.1, 17.8, 26.7),
                    (10.0, 17.8, 30.0),
                    (8.9, 17.8, 33.3),
                    (8.9, 22.2, 33.3),
                    (11.1, 22.2, 33.3),
                    (11.1, 17.8, 33.3),
                    (8.7, 20.0, 26.1),
                    (11.3, 20.0, 26.1),
                    (8.7, 20.0, 33.9),
                    (11.3, 20.0, 33.9),
                    (8.7, 17.4, 30.0),
                    (11.3, 22.6, 30.0),
                    (10.0, 20.0, 26.1),
                    (10.0, 22.6, 30.0),
                    (10.0, 17.4, 30.0),
                    (14.0, 20.0, 30.0),
                ]
            ),
            NORM_HEXA20: DataArrayDouble(
                [
                    (8.9, 20.0, 33.3),
                    (10.0, 21.4, 32.1),
                    (11.1, 20.0, 33.3),
                    (10.0, 22.2, 30.0),
                    (10.0, 20.0, 32.1),
                    (10.0, 20.0, 26.7),
                    (10.0, 17.8, 30.0),
                    (10.0, 18.6, 32.1),
                    (8.9, 17.8, 26.7),
                    (8.9, 22.2, 26.7),
                    (11.1, 22.2, 26.7),
                    (11.1, 17.8, 26.7),
                    (8.9, 17.8, 33.3),
                    (8.9, 22.2, 33.3),
                    (11.1, 22.2, 33.3),
                    (11.1, 17.8, 33.3),
                    (8.7, 20.0, 26.1),
                    (10.0, 21.8, 27.3),
                    (11.3, 20.0, 26.1),
                    (10.0, 18.2, 27.3),
                    (8.7, 20.0, 33.9),
                    (10.0, 21.8, 32.7),
                    (11.3, 20.0, 33.9),
                    (10.0, 18.2, 32.7),
                    (8.7, 17.4, 30.0),
                    (9.1, 21.8, 30.0),
                    (11.3, 22.6, 30.0),
                    (10.9, 18.2, 30.0),
                    (10.0, 20.0, 26.1),
                    (10.0, 22.6, 30.0),
                    (10.0, 17.4, 30.0),
                    (10.0, 20.0, 32.7),
                    (14.0, 20.0, 30.0),
                ]
            ),
            NORM_HEXA27: DataArrayDouble(
                [
                    (8.9, 17.8, 26.7),
                    (8.9, 22.2, 26.7),
                    (11.1, 22.2, 26.7),
                    (11.1, 17.8, 26.7),
                    (8.9, 17.8, 33.3),
                    (8.9, 22.2, 33.3),
                    (11.1, 22.2, 33.3),
                    (11.1, 17.8, 33.3),
                    (8.7, 20.0, 26.1),
                    (10.0, 21.8, 27.3),
                    (11.3, 20.0, 26.1),
                    (10.0, 18.2, 27.3),
                    (8.7, 20.0, 33.9),
                    (10.0, 21.8, 32.7),
                    (11.3, 20.0, 33.9),
                    (10.0, 18.2, 32.7),
                    (8.7, 17.4, 30.0),
                    (9.1, 21.8, 30.0),
                    (11.3, 22.6, 30.0),
                    (10.9, 18.2, 30.0),
                    (10.0, 20.0, 26.1),
                    (9.1, 20.0, 30.0),
                    (10.0, 22.6, 30.0),
                    (10.9, 20.0, 30.0),
                    (10.0, 17.4, 30.0),
                    (10.0, 20.0, 32.7),
                    (14.0, 20.0, 30.0),
                ]
            ),
        }

        for gt in [NORM_HEXA8, NORM_HEXA20, NORM_HEXA27]:
            nbPtsInCell = MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(gt)
            m = MEDCouplingUMesh("mesh", 3)
            m.setCoords(coo)
            m.allocateCells()
            m.insertNextCell(gt, list(range(nbPtsInCell)))
            m.finishInsertingCells()

            srcField = MEDCouplingFieldDouble(ON_NODES_FE)
            srcField.setMesh(m)
            arr = DataArrayDouble(nbPtsInCell)
            arr.iota()
            srcField.setArray(arr)

            inPts = inPts_all[gt]
            outPts = outPts_all[gt]
            ref_val0 = ref_val0_all[gt]
            nbInPts = len(inPts)

            for i in range(nbInPts):
                self.assertTrue(
                    abs(srcField.getValueOn(inPts[i])[0] - ref_val0[i]) < eps
                )

            self.assertTrue(srcField.getValueOnMulti(inPts).isEqual(ref_val0, eps))

            srcFt = MEDCouplingFieldTemplate(srcField)
            trgFt = MEDCouplingFieldTemplate(ON_NODES_FE)
            trgMesh = MEDCouplingUMesh.Build0DMeshFromCoords(inPts)
            trgFt.setMesh(trgMesh)
            rem = MEDCouplingRemapper()
            rem.setIntersectionType(PointLocator)
            rem.prepareEx(srcFt, trgFt)
            # scan content of matrix computed by remapper
            mat = rem.getCrudeMatrix()
            self.assertEqual(len(mat), nbInPts)
            for irow, row in enumerate(mat):
                self.assertTrue(abs(sum([y for x, y in row.items()]) - 1.0) < eps)

            # ask for MEDCouplingFieldDiscretizationOnNodesFE instance to compute coordination into ref element
            sd = srcField.getDiscretization()
            coosInRefElemFoundByNewton = sd.getCooInRefElement(
                srcField.getMesh(), inPts
            )

            for zePt, cooInRefElemFoundByNewton in zip(
                inPts, coosInRefElemFoundByNewton
            ):
                # now check by performing refCoo -> realCoo
                ftest = MEDCouplingFieldDouble(ON_GAUSS_PT)
                ftest.setMesh(srcField.getMesh())
                ftest.setGaussLocalizationOnType(
                    gt,
                    sum(
                        [
                            list(elt)
                            for elt in MEDCouplingGaussLocalization.GetDefaultReferenceCoordinatesOf(
                                gt
                            ).getValuesAsTuple()
                        ],
                        [],
                    ),
                    list(cooInRefElemFoundByNewton),
                    [1],
                )
                self.assertTrue(
                    float((ftest.getLocalizationOfDiscr() - zePt).magnitude()) < 1e-4
                )

            # testing with outside point
            for zePt in outPts:
                # safer than
                # self.assertRaises(InterpKernelException,sd.getCooInRefElement,srcField.getMesh(),zePt.buildDADouble())
                try:
                    sd.getCooInRefElement(srcField.getMesh(), zePt.buildDADouble())
                except InterpKernelException as e:
                    self.assertTrue("fail to locate point" in e.what())
                else:
                    self.assertTrue(False, "")

    def testFEFE_TETRA(self):
        """
        [EDF31187] : Test to stress localisation of target points into a source mesh using standard reference FE elements.

        Test for TETRA4 and TETRAb10 with non-planar faces
        """

        coo = DataArrayDouble(
            [
                (1.0, 5.0, 8.0),
                (1.0, 3.0, 11.0),
                (1.1, 4.0, 8.0),
                (7.0, 3.0, 7.0),
                (1.0, 4.0, 9.5),
                (1.0, 3.5, 9.0),
                (1.1, 4.5, 8),
                (3.0, 4.5, 8.0),
                (3.5, 3.5, 9.0),
                (3.5, 3.0, 7.8),
            ]
        )

        eps = 1e-8

        ref_val0_all = {
            NORM_TETRA4: DataArrayDouble(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    1.4069767441860472,
                ]
            ),
            NORM_TETRA10: DataArrayDouble(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    7.0,
                    8.0,
                    9.0,
                    9.3606435183767083,
                ]
            ),
        }

        inPts_all = {
            NORM_TETRA4: DataArrayDouble(
                [
                    (1.0, 5.0, 8.0),
                    (1.0, 3.0, 11.0),
                    (1.1, 4.0, 8.0),
                    (7.0, 3.0, 7.0),
                    (2.5, 3.8, 8.5),
                ]
            ),
            NORM_TETRA10: DataArrayDouble(
                [
                    (1.0, 5.0, 8.0),
                    (1.0, 3.0, 11.0),
                    (1.1, 4.0, 8.0),
                    (7.0, 3.0, 7.0),
                    (1.0, 4.0, 9.5),
                    (1.0, 3.5, 9.0),
                    (1.1, 4.5, 8),
                    (3.0, 4.5, 8.0),
                    (3.5, 3.5, 9.0),
                    (3.5, 3.0, 7.8),
                    (2.5, 3.8, 8.5),
                ]
            ),
        }

        outPts_all = {
            NORM_TETRA4: DataArrayDouble(
                [
                    (0.0, 0.0, 0.0),
                    (1.0, 3.5, 9.00001),
                    (7.00001, 3.0, 7.0),
                    (6.99, 3.0, 7.0),
                ]
            ),
            NORM_TETRA10: DataArrayDouble(
                [
                    (0.0, 0.0, 0.0),
                    (1.0, 3.48, 9.0),
                    (7.00001, 3.0, 7.0),
                    (6.99, 3.0, 7.0),
                ]
            ),
        }

        for gt in [NORM_TETRA4, NORM_TETRA10]:
            nbPtsInCell = MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(gt)
            m = MEDCouplingUMesh("mesh", 3)
            m.setCoords(coo)
            m.allocateCells()
            m.insertNextCell(gt, list(range(nbPtsInCell)))
            m.finishInsertingCells()

            srcField = MEDCouplingFieldDouble(ON_NODES_FE)
            srcField.setMesh(m)
            arr = DataArrayDouble(nbPtsInCell)
            arr.iota()
            srcField.setArray(arr)

            inPts = inPts_all[gt]
            outPts = outPts_all[gt]
            ref_val0 = ref_val0_all[gt]
            nbInPts = len(inPts)

            for i in range(nbInPts):
                self.assertTrue(
                    abs(srcField.getValueOn(inPts[i])[0] - ref_val0[i]) < eps
                )

            self.assertTrue(srcField.getValueOnMulti(inPts).isEqual(ref_val0, eps))

            srcFt = MEDCouplingFieldTemplate(srcField)
            trgFt = MEDCouplingFieldTemplate(ON_NODES_FE)
            trgMesh = MEDCouplingUMesh.Build0DMeshFromCoords(inPts)
            trgFt.setMesh(trgMesh)
            rem = MEDCouplingRemapper()
            rem.setIntersectionType(PointLocator)
            rem.prepareEx(srcFt, trgFt)
            # scan content of matrix computed by remapper
            mat = rem.getCrudeMatrix()
            self.assertEqual(len(mat), nbInPts)
            for irow, row in enumerate(mat):
                self.assertTrue(abs(sum([y for x, y in row.items()]) - 1.0) < eps)

            # ask for MEDCouplingFieldDiscretizationOnNodesFE instance to compute coordination into ref element
            sd = srcField.getDiscretization()
            coosInRefElemFoundByNewton = sd.getCooInRefElement(
                srcField.getMesh(), inPts
            )

            for zePt, cooInRefElemFoundByNewton in zip(
                inPts, coosInRefElemFoundByNewton
            ):
                # now check by performing refCoo -> realCoo
                ftest = MEDCouplingFieldDouble(ON_GAUSS_PT)
                ftest.setMesh(srcField.getMesh())
                ftest.setGaussLocalizationOnType(
                    gt,
                    sum(
                        [
                            list(elt)
                            for elt in MEDCouplingGaussLocalization.GetDefaultReferenceCoordinatesOf(
                                gt
                            ).getValuesAsTuple()
                        ],
                        [],
                    ),
                    list(cooInRefElemFoundByNewton),
                    [1],
                )
                self.assertTrue(
                    float((ftest.getLocalizationOfDiscr() - zePt).magnitude()) < 1e-4
                )

            # testing with outside point
            for zePt in outPts:
                # safer than
                # self.assertRaises(InterpKernelException,sd.getCooInRefElement,srcField.getMesh(),zePt.buildDADouble())
                try:
                    sd.getCooInRefElement(srcField.getMesh(), zePt.buildDADouble())
                except InterpKernelException as e:
                    self.assertTrue("fail to locate point" in e.what())
                else:
                    self.assertTrue(False, "")

    def testFEFE_PENTA(self):
        """
        [EDF31187] : Test to stress localisation of target points into a source mesh using standard reference FE elements.

        Test for NORM_PENTA6, NORM_PENTA15 and NORM_PENTA18 with non-planar faces
        """

        coo = DataArrayDouble(
            [
                (-3.0, 1.0, 8.0),
                (-3.0, -2.0, 11.0),
                (-2.0, -2.0, 8.0),
                (2.5, 1.3, 7.0),
                (2.8, -0.1, 9.5),
                (2.8, -0.2, 9.0),
                (-3.1, -0.5, 9.5),
                (-3.1, -2.0, 9.5),
                (-2.7, -0.5, 8.0),
                (0.0, 1.15, 8.0),
                (0.0, -1.0, 10.5),
                (0.0, -1.1, 7.80),
                (2.8, 0.7, 8.0),
                (2.8, -0.2, 9.5),
                (2.6, 0.8, 8.0),
                (0.0, 0.0, 9.0),
                (0.0, -1.0, 9.2),
                (0.1, 0.0, 7.9),
            ]
        )

        eps = 1e-8

        ref_val0_all = {
            NORM_PENTA6: DataArrayDouble(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    2.6065933924123641,
                ]
            ),
            NORM_PENTA15: DataArrayDouble(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    7.0,
                    8.0,
                    9.0,
                    10.0,
                    11.0,
                    12.0,
                    13.0,
                    14.0,
                    19.917147597384314,
                ]
            ),
            NORM_PENTA18: DataArrayDouble(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    7.0,
                    8.0,
                    9.0,
                    10.0,
                    11.0,
                    12.0,
                    13.0,
                    14.0,
                    15.0,
                    16.0,
                    17.0,
                    17.868885100760988,
                ]
            ),
        }

        inPts_all = {
            NORM_PENTA6: DataArrayDouble(
                [
                    (-3.0, 1.0, 8.0),
                    (-3.0, -2.0, 11.0),
                    (-2.0, -2.0, 8.0),
                    (2.5, 1.3, 7.0),
                    (2.8, -0.1, 9.5),
                    (2.8, -0.2, 9.0),
                    (0.0, -0.5, 8.8),
                ]
            ),
            NORM_PENTA15: DataArrayDouble(
                [
                    (-3.0, 1.0, 8.0),
                    (-3.0, -2.0, 11.0),
                    (-2.0, -2.0, 8.0),
                    (2.5, 1.3, 7.0),
                    (2.8, -0.1, 9.5),
                    (2.8, -0.2, 9.0),
                    (-3.1, -0.5, 9.5),
                    (-3.1, -2.0, 9.5),
                    (-2.7, -0.5, 8.0),
                    (0.0, 1.15, 8.0),
                    (0.0, -1.0, 10.5),
                    (0.0, -1.1, 7.80),
                    (2.8, 0.7, 8.0),
                    (2.8, -0.2, 9.5),
                    (2.6, 0.8, 8.0),
                    (0.0, -0.5, 8.8),
                ]
            ),
            NORM_PENTA18: DataArrayDouble(
                [
                    (-3.0, 1.0, 8.0),
                    (-3.0, -2.0, 11.0),
                    (-2.0, -2.0, 8.0),
                    (2.5, 1.3, 7.0),
                    (2.8, -0.1, 9.5),
                    (2.8, -0.2, 9.0),
                    (-3.1, -0.5, 9.5),
                    (-3.1, -2.0, 9.5),
                    (-2.7, -0.5, 8.0),
                    (0.0, 1.15, 8.0),
                    (0.0, -1.0, 10.5),
                    (0.0, -1.1, 7.80),
                    (2.8, 0.7, 8.0),
                    (2.8, -0.2, 9.5),
                    (2.6, 0.8, 8.0),
                    (0.0, 0.0, 9.0),
                    (0.0, -1.0, 9.2),
                    (0.1, 0.0, 7.9),
                    (0.0, -0.5, 8.8),
                ]
            ),
        }

        outPts_all = {
            NORM_PENTA6: DataArrayDouble(
                [
                    (-7.0, 10e50, 0.0),
                    (-3.00001, 1.0, 8.0),
                    (-3.0, -2.0, 11.00001),
                ]
            ),
            NORM_PENTA15: DataArrayDouble(
                [
                    (-7.0, 10e50, 0.0),
                    (-3.00001, 1.0, 8.0),
                    (-3.0, -2.0, 11.00001),
                ]
            ),
            NORM_PENTA18: DataArrayDouble(
                [
                    (-7.0, 10e50, 0.0),
                    (-3.00001, 1.0, 8.0),
                    (-3.0, -2.0, 11.00001),
                ]
            ),
        }

        for gt in [NORM_PENTA6, NORM_PENTA15, NORM_PENTA18]:
            nbPtsInCell = MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(gt)
            m = MEDCouplingUMesh("mesh", 3)
            m.setCoords(coo)
            m.allocateCells(1)
            m.insertNextCell(gt, list(range(nbPtsInCell)))
            m.finishInsertingCells()

            srcField = MEDCouplingFieldDouble(ON_NODES_FE)
            srcField.setMesh(m)
            arr = DataArrayDouble(nbPtsInCell)
            arr.iota()
            srcField.setArray(arr)

            inPts = inPts_all[gt]
            outPts = outPts_all[gt]
            ref_val0 = ref_val0_all[gt]
            nbInPts = len(inPts)

            for i in range(nbInPts):
                self.assertTrue(
                    abs(srcField.getValueOn(inPts[i])[0] - ref_val0[i]) < eps
                )

            self.assertTrue(srcField.getValueOnMulti(inPts).isEqual(ref_val0, eps))

            srcFt = MEDCouplingFieldTemplate(srcField)
            trgFt = MEDCouplingFieldTemplate(ON_NODES_FE)
            trgMesh = MEDCouplingUMesh.Build0DMeshFromCoords(inPts)
            trgFt.setMesh(trgMesh)
            rem = MEDCouplingRemapper()
            rem.setIntersectionType(PointLocator)
            rem.prepareEx(srcFt, trgFt)
            # scan content of matrix computed by remapper
            mat = rem.getCrudeMatrix()
            self.assertEqual(len(mat), nbInPts)
            for irow, row in enumerate(mat):
                self.assertTrue(abs(sum([y for x, y in row.items()]) - 1.0) < eps)

            # ask for MEDCouplingFieldDiscretizationOnNodesFE instance to compute coordination into ref element
            sd = srcField.getDiscretization()
            coosInRefElemFoundByNewton = sd.getCooInRefElement(
                srcField.getMesh(), inPts
            )

            for zePt, cooInRefElemFoundByNewton in zip(
                inPts, coosInRefElemFoundByNewton
            ):
                # now check by performing refCoo -> realCoo
                ftest = MEDCouplingFieldDouble(ON_GAUSS_PT)
                ftest.setMesh(srcField.getMesh())
                ftest.setGaussLocalizationOnType(
                    gt,
                    sum(
                        [
                            list(elt)
                            for elt in MEDCouplingGaussLocalization.GetDefaultReferenceCoordinatesOf(
                                gt
                            ).getValuesAsTuple()
                        ],
                        [],
                    ),
                    list(cooInRefElemFoundByNewton),
                    [1],
                )
                self.assertTrue(
                    float((ftest.getLocalizationOfDiscr() - zePt).magnitude()) < 1e-4
                )

            # testing with outside point
            for zePt in outPts:
                # safer than
                # self.assertRaises(InterpKernelException,sd.getCooInRefElement,srcField.getMesh(),zePt.buildDADouble())
                try:
                    sd.getCooInRefElement(srcField.getMesh(), zePt.buildDADouble())
                except InterpKernelException as e:
                    self.assertTrue("fail to locate point" in e.what())
                else:
                    self.assertTrue(False, "")

    def testFEFE_PYRAM(self):
        """
        [EDF31187] : Test to stress localisation of target points into a source mesh using standard reference FE elements.

        Test for NORM_PYRA5 and NORM_PYRA13 with non-planar faces
        """

        coo = DataArrayDouble(
            [
                (-3.0, 3.0, 4.0),
                (-3.0, -2.0, 4.1),
                (2.5, -2.1, 4.0),
                (2.6, 3.0, 4.1),
                (-0.3, -0.1, 9.5),
                (-3.0, 0.5, 4.05),
                (-0.25, -2.05, 4.05),
                (2.55, 0.45, 4.05),
                (-0.2, 3.0, 4.05),
                (-1.5, 1.55, 6.78),
                (-1.65, -1.05, 6.8),
                (1.1, -1.3, 6.72),
                (1.3, 1.45, 6.8),
            ]
        )

        eps = 1e-8

        ref_val0_all = {
            NORM_PYRA5: DataArrayDouble(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    2.8248558088942426,
                ]
            ),
            NORM_PYRA13: DataArrayDouble(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    7.0,
                    8.0,
                    9.0,
                    10.0,
                    11.0,
                    12.0,
                    12.928985087457612,
                ]
            ),
        }

        inPts_all = {
            NORM_PYRA5: DataArrayDouble(
                [
                    (-3.0, 3.0, 4.0),
                    (-3.0, -2.0, 4.1),
                    (2.5, -2.1, 4.0),
                    (2.6, 3.0, 4.1),
                    (-0.3, -0.1, 9.5),
                    (0.5, 1.1, 6.2),
                ]
            ),
            NORM_PYRA13: DataArrayDouble(
                [
                    (-3.0, 3.0, 4.0),
                    (-3.0, -2.0, 4.1),
                    (2.5, -2.1, 4.0),
                    (2.6, 3.0, 4.1),
                    (-0.3, -0.1, 9.5),
                    (-3.0, 0.5, 4.05),
                    (-0.25, -2.05, 4.05),
                    (2.55, 0.45, 4.05),
                    (-0.2, 3.0, 4.05),
                    (-1.5, 1.55, 6.78),
                    (-1.65, -1.05, 6.8),
                    (1.1, -1.3, 6.72),
                    (1.3, 1.45, 6.8),
                    (0.5, 1.1, 6.76),
                ]
            ),
        }

        outPts_all = {
            NORM_PYRA5: DataArrayDouble(
                [
                    (-7.0, 10e50, 0.0),
                    (-3.00001, 3.0, 4.0),
                    (-0.3, -0.1, 9.500001),
                ]
            ),
            NORM_PYRA13: DataArrayDouble(
                [
                    (-7.0, 10e50, 0.0),
                    (-3.00001, 1.0, 8.0),
                    (-0.3, -0.1, 9.500001),
                    (1.10001, -1.3, 6.72),
                ]
            ),
        }

        for gt in [NORM_PYRA5, NORM_PYRA13]:
            nbPtsInCell = MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(gt)
            m = MEDCouplingUMesh("mesh", 3)
            m.setCoords(coo)
            m.allocateCells(1)
            m.insertNextCell(gt, list(range(nbPtsInCell)))
            m.finishInsertingCells()

            srcField = MEDCouplingFieldDouble(ON_NODES_FE)
            srcField.setMesh(m)
            arr = DataArrayDouble(nbPtsInCell)
            arr.iota()
            srcField.setArray(arr)

            inPts = inPts_all[gt]
            outPts = outPts_all[gt]
            ref_val0 = ref_val0_all[gt]
            nbInPts = len(inPts)

            for i in range(nbInPts):
                self.assertTrue(
                    abs(srcField.getValueOn(inPts[i])[0] - ref_val0[i]) < eps
                )

            self.assertTrue(srcField.getValueOnMulti(inPts).isEqual(ref_val0, eps))

            srcFt = MEDCouplingFieldTemplate(srcField)
            trgFt = MEDCouplingFieldTemplate(ON_NODES_FE)
            trgMesh = MEDCouplingUMesh.Build0DMeshFromCoords(inPts)
            trgFt.setMesh(trgMesh)
            rem = MEDCouplingRemapper()
            rem.setIntersectionType(PointLocator)
            rem.prepareEx(srcFt, trgFt)
            # scan content of matrix computed by remapper
            mat = rem.getCrudeMatrix()
            self.assertEqual(len(mat), nbInPts)
            for irow, row in enumerate(mat):
                self.assertTrue(abs(sum([y for x, y in row.items()]) - 1.0) < eps)

            # ask for MEDCouplingFieldDiscretizationOnNodesFE instance to compute coordination into ref element
            sd = srcField.getDiscretization()
            coosInRefElemFoundByNewton = sd.getCooInRefElement(
                srcField.getMesh(), inPts
            )

            for zePt, cooInRefElemFoundByNewton in zip(
                inPts, coosInRefElemFoundByNewton
            ):
                # now check by performing refCoo -> realCoo
                ftest = MEDCouplingFieldDouble(ON_GAUSS_PT)
                ftest.setMesh(srcField.getMesh())
                ftest.setGaussLocalizationOnType(
                    gt,
                    sum(
                        [
                            list(elt)
                            for elt in MEDCouplingGaussLocalization.GetDefaultReferenceCoordinatesOf(
                                gt
                            ).getValuesAsTuple()
                        ],
                        [],
                    ),
                    list(cooInRefElemFoundByNewton),
                    [1],
                )
                self.assertTrue(
                    float((ftest.getLocalizationOfDiscr() - zePt).magnitude()) < 1e-4
                )

            # testing with outside point
            for zePt in outPts:
                # safer than
                # self.assertRaises(InterpKernelException,sd.getCooInRefElement,srcField.getMesh(),zePt.buildDADouble())
                try:
                    sd.getCooInRefElement(srcField.getMesh(), zePt.buildDADouble())
                except InterpKernelException as e:
                    self.assertTrue("fail to locate point" in e.what())
                else:
                    self.assertTrue(False, "")

    def testFEFE_PROJ2D(self):
        """
        [EDF31187] : Test to project a field using standard reference FE elements.
        """

        eps = 1e-8

        # Source mesh
        srcCoo = DataArrayDouble(
            [
                (9.0, 18.0),
                (11.0, 18.0),
                (11.0, 22.0),
                (9.0, 22.0),
                (10, 18.1),
                (11.0, 20.0),
                (10, 21.9),
                (9.0, 20.0),
                (10.0, 20.0),
            ]
        )

        srcMesh = MEDCouplingUMesh("srcMesh", 2)
        srcMesh.setCoords(srcCoo)
        srcMesh.allocateCells()
        srcMesh.insertNextCell(NORM_QUAD4, (0, 4, 8, 7))
        srcMesh.insertNextCell(NORM_QUAD4, (4, 1, 5, 8))
        srcMesh.insertNextCell(NORM_QUAD4, (8, 5, 2, 6))
        srcMesh.insertNextCell(NORM_QUAD4, (7, 8, 6, 3))
        srcMesh.finishInsertingCells()

        # Target mesh
        trgCoo = DataArrayDouble(
            [
                (10.0, 20.0),
                (11.0, 20.0),
                (10.0, 22.5),
                (10.5, 20.0),
                (10.55, 21.26),
                (9.99, 21.25),
            ]
        )

        trgMesh = MEDCouplingUMesh("trgMesh", 2)
        trgMesh.setCoords(trgCoo)
        trgMesh.allocateCells()
        trgMesh.insertNextCell(NORM_TRI6, list(range(6)))
        trgMesh.finishInsertingCells()

        # Source Field
        srcField = MEDCouplingFieldDouble(ON_NODES_FE)
        srcField.setMesh(srcMesh)
        srcField.fillFromAnalytic(1, "x+y")
        srcField.setNature(IntensiveMaximum)

        # Target Field
        trgField = MEDCouplingFieldDouble(ON_NODES_FE)
        trgField.setMesh(trgMesh)
        trgField.fillFromAnalytic(1, "x+y")
        trgField.setNature(IntensiveMaximum)

        # Interpolate nodes to nodes using FEM interpolation
        remap = MEDCouplingRemapper()
        remap.setIntersectionType(PointLocator)

        # Workaround since "FEFE" interpolation does not exist yet
        srcFt = MEDCouplingFieldTemplate(srcField)
        trgFt = MEDCouplingFieldTemplate(ON_NODES_FE)
        trgFt.setMesh(trgMesh)

        remap.prepareEx(srcFt, trgFt)
        # remap.prepare(srcMesh, trgMesh, "FEFE")

        # Check matrix
        myMatrix = remap.getCrudeMatrix()
        self.assertEqual(len(myMatrix), 6)
        # tested values are constructed for this test
        for irow, row in enumerate(myMatrix):
            if irow != 2:
                self.assertTrue(abs(sum([y for x, y in row.items()]) - 1.0) < eps)
            else:
                self.assertTrue(len(row) == 0)

        # Transfer field src -> trg
        trgFieldCV = remap.transferField(srcField, 1e300)
        trgFieldAN = trgField.deepCopy()

        # Projection of a linear field hence no interpolation error
        trgFieldAN.getArray()[2] = 1e300
        self.assertTrue(trgFieldCV.getArray().isEqual(trgFieldAN.getArray(), eps))

        # Transfer field trg -> src (should not work since need two matrices)
        try:
            srcFieldCV = remap.transferField(trgField, 2e300)
        except InterpKernelException as e:
            self.assertTrue("MEDCouplingRemapper::transferUnderground" in e.what())
        else:
            self.assertTrue(False, "")

    def testFEFE_PROJ3D(self):
        """
        [EDF31187] : Test to project a 3D field using standard reference FE elements.
        """

        eps = 1e-8

        # Source mesh
        srcCoo = DataArrayDouble(
            [
                (9.0, 18.0, 0.0),
                (11.0, 18.0, 0.0),
                (11.0, 22.0, 0.0),
                (9.0, 22.0, 0.0),
                (10, 18.1, 0.0),
                (11.0, 20.0, 0.0),
                (10, 21.9, 0.0),
                (9.0, 20.0, 0.0),
                (10.0, 20.0, 0.0),
                (9.0, 18.0, 2.0),
                (11.0, 18.0, 2.0),
                (11.0, 22.0, 2.0),
                (9.0, 22.0, 2.0),
                (10, 18.1, 2.0),
                (11.0, 20.0, 2.0),
                (10, 21.9, 2.0),
                (9.0, 20.0, 2.0),
                (10.0, 20.0, 2.0),
            ]
        )

        srcMesh = MEDCouplingUMesh("srcMesh", 3)
        srcMesh.setCoords(srcCoo)
        srcMesh.allocateCells(4)
        srcMesh.insertNextCell(NORM_HEXA8, (0, 7, 8, 4, 9, 16, 17, 13))
        srcMesh.insertNextCell(NORM_HEXA8, (4, 8, 5, 1, 13, 17, 14, 10))
        srcMesh.insertNextCell(NORM_HEXA8, (8, 6, 2, 5, 17, 15, 11, 14))
        srcMesh.insertNextCell(NORM_HEXA8, (7, 3, 6, 8, 16, 12, 15, 17))
        srcMesh.finishInsertingCells()

        # Target mesh
        trgCoo = DataArrayDouble(
            [
                (10.0, 20.0, 0.0),
                (11.0, 20.0, 0.0),
                (10.0, 22.5, 0.0),
                (10.0, 20.0, 1.0),
                (10.9, 19.9, 1.0),
                (10.0, 21.9, 1.0),
            ]
        )

        trgMesh = MEDCouplingUMesh("trgMesh", 3)
        trgMesh.setCoords(trgCoo)
        trgMesh.allocateCells(1)
        trgMesh.insertNextCell(NORM_PENTA6, list(range(6)))
        trgMesh.finishInsertingCells()

        # Source Field
        srcField = MEDCouplingFieldDouble(ON_NODES_FE)
        srcField.setMesh(srcMesh)
        srcField.fillFromAnalytic(1, "x+y+z")
        srcField.setNature(IntensiveMaximum)

        # Target Field
        trgField = MEDCouplingFieldDouble(ON_NODES_FE)
        trgField.setMesh(trgMesh)
        trgField.fillFromAnalytic(1, "x+y+z")
        trgField.setNature(IntensiveMaximum)

        # Interpolate nodes to nodes using FEM interpolation
        remap = MEDCouplingRemapper()
        remap.setIntersectionType(PointLocator)

        # Workaround since "FEFE" interpolation does not exist yet
        srcFt = MEDCouplingFieldTemplate(srcField)
        trgFt = MEDCouplingFieldTemplate(ON_NODES_FE)
        trgFt.setMesh(trgMesh)

        remap.prepareEx(srcFt, trgFt)
        # remap.prepare(srcMesh, trgMesh, "FEFE")

        # Check matrix
        myMatrix = remap.getCrudeMatrix()
        self.assertEqual(len(myMatrix), 6)
        # tested values are constructed for this test
        for irow, row in enumerate(myMatrix):
            if irow not in (2,):
                self.assertTrue(abs(sum([abs(y) for x, y in row.items()]) - 1.0) < eps)
            else:
                self.assertTrue(len(row) == 0)

        # test specific values
        self.assertTrue(abs(myMatrix[0][8] - 1.0) < eps)
        self.assertTrue(abs(myMatrix[1][5] - 1.0) < eps)
        self.assertTrue(abs(myMatrix[3][8] - 0.5) < eps)
        self.assertTrue(abs(myMatrix[3][17] - 0.5) < eps)
        for i in (4, 1, 5, 8, 15, 10, 14, 17):
            self.assertTrue(myMatrix[4][17] > eps)
        self.assertTrue(abs(myMatrix[5][6] - 0.5) < eps)
        self.assertTrue(abs(myMatrix[5][15] - 0.5) < eps)

        # Transfer field src -> trg
        trgFieldCV = remap.transferField(srcField, 1e300)
        trgFieldAN = trgField.deepCopy()

        # Projection of a linear field hence no interpolation error
        trgFieldAN.getArray()[2] = 1e300
        for i in range(6):
            self.assertTrue((trgFieldCV.getArray()[i] - trgFieldAN.getArray()[i]) < eps)

        # Transfer field trg -> src (should not work since need two matrices)
        try:
            srcFieldCV = remap.transferField(trgField, 2e300)
        except InterpKernelException as e:
            self.assertTrue("MEDCouplingRemapper::transferUnderground" in e.what())
        else:
            self.assertTrue(False, "")

    def testFEFE_HEXA_SURF(self):
        """
        [EDF31887] : Test to stress localisation of target points into a source mesh using standard reference FE elements.

        Test for HEXA8, HEXA20 and HEXA27 with non-planar faces
        """

        coo = DataArrayDouble(
            [
                (9.0, 18.0, 27.0),
                (9.0, 22.0, 27.0),
                (11.0, 22.0, 27.0),
                (11.0, 18.0, 27.0),
                (9.0, 18.0, 33.0),
                (9.0, 22.0, 33.0),
                (11.0, 22.0, 33.0),
                (11.0, 18.0, 33.0),
                (8.8, 20.0, 26.4),
                (10.0, 21.6, 27.6),
                (11.2, 20.0, 26.4),
                (10.0, 18.4, 27.6),
                (8.8, 20.0, 33.6),
                (10.0, 21.6, 32.4),
                (11.2, 20.0, 33.6),
                (10.0, 18.4, 32.4),
                (8.8, 17.6, 30.0),
                (9.2, 21.6, 30.0),
                (11.2, 22.4, 30.0),
                (10.8, 18.4, 30.0),
                (10.0, 20.0, 26.4),
                (9.2, 20.0, 30.0),
                (10.0, 22.4, 30.0),
                (10.8, 20.0, 30.0),
                (10.0, 17.6, 30.0),
                (10.0, 20.0, 32.4),
                (10.0, 20.0, 30.0),
            ]
        )

        eps = 1e-8

        ref_val0_all = {
            NORM_HEXA8: DataArrayDouble(
                [
                    (0.0, 0.0),
                    (0.0, 0.0),
                    (-1.0, 1.0),
                    (-2.0 / 3.0, 1.0),
                ]
            ),
        }

        outPts_all = {
            NORM_HEXA8: DataArrayDouble(
                [
                    (12.0, 20.0, 30.0),
                    (10.8, 20.0, 30.0),
                    (8.0, 18.0, 27.0),
                    (8.0, 17.0, 28.0),
                ]
            ),
        }

        for gt in [
            NORM_HEXA8,
        ]:
            m = MEDCouplingUMesh("mesh", 2)
            m.setCoords(coo)
            m.allocateCells(6)
            m.insertNextCell(NORM_QUAD4, (0, 1, 2, 3))
            m.insertNextCell(NORM_QUAD4, (4, 5, 6, 7))
            m.insertNextCell(NORM_QUAD4, (0, 1, 5, 4))
            m.insertNextCell(NORM_QUAD4, (1, 2, 6, 5))
            m.insertNextCell(NORM_QUAD4, (2, 3, 7, 6))
            m.insertNextCell(NORM_QUAD4, (0, 3, 7, 4))
            m.finishInsertingCells()

            srcField = MEDCouplingFieldDouble(ON_NODES_FE)
            srcField.setMesh(m)
            sd = srcField.getDiscretization()

            # testing with outside point
            coosInRefElemFoundByNewton = sd.getClosestCooInRefElement(
                srcField.getMesh(), outPts_all[gt]
            )

            self.assertTrue(coosInRefElemFoundByNewton.isEqual(ref_val0_all[gt], eps))

    def testFEFE_PROJ3DSurf(self):
        """
        [EDF31187] : Test to project a 3D field using standard reference FE elements.
        """

        def cylinder(nbPt, quadratic=False, simplexe=False):
            # base
            baseCoo = DataArrayDouble(
                [
                    (1.0, 0.0),
                    (2.0, 0.0),
                ]
            )

            baseMesh = MEDCouplingUMesh("baseMesh", 1)
            baseMesh.setCoords(baseCoo)
            baseMesh.allocateCells(1)
            baseMesh.insertNextCell(NORM_SEG2, (0, 1))
            baseMesh.finishInsertingCells()

            lineCoo = []
            for i in range(nbPt + 1):
                angle = i * (pi / nbPt)
                x = cos(angle)
                y = sin(angle)
                lineCoo.append((x, y))

            lineCoo = DataArrayDouble(lineCoo)
            lineMesh = MEDCouplingUMesh("lineMesh", 1)
            lineMesh.setCoords(lineCoo)
            lineMesh.allocateCells(nbPt)
            for i in range(nbPt):
                lineMesh.insertNextCell(NORM_SEG2, (i, i + 1))
            lineMesh.finishInsertingCells()

            extrudedMesh = baseMesh.buildExtrudedMesh(lineMesh, 1)
            if simplexe:
                extrudedMesh.simplexize(0)

            skin = extrudedMesh.computeSkin()
            if quadratic:
                skin.convertLinearCellsToQuadratic(1)

            return skin

        for simplex in [False, True]:
            for quadratic in [False, True]:
                # Source mesh
                srcSkinMesh = cylinder(9, quadratic, simplex)

                # Target mesh
                trgSkinMesh = cylinder(5)

                # Source Field
                srcField = MEDCouplingFieldDouble(ON_NODES_FE)
                srcField.setMesh(srcSkinMesh)
                srcField.fillFromAnalytic(1, "x+y")
                srcField.setNature(IntensiveMaximum)

                # Target Field
                trgField = MEDCouplingFieldDouble(ON_NODES_FE)
                trgField.setMesh(trgSkinMesh)
                trgField.fillFromAnalytic(1, "x+y")
                trgField.setNature(IntensiveMaximum)

                # Interpolate nodes to nodes using FEM interpolation
                remap = MEDCouplingRemapper()
                remap.setIntersectionType(PointLocator)

                # Workaround since "FEFE" interpolation does not exist yet
                srcFt = MEDCouplingFieldTemplate(srcField)
                trgFt = MEDCouplingFieldTemplate(ON_NODES_FE)
                trgFt.setMesh(trgSkinMesh)

                remap.prepareEx(srcFt, trgFt)
                # remap.prepare(srcMesh, trgMesh, "FEFE")

                # Check matrix
                myMatrix = remap.getCrudeMatrix()
                self.assertEqual(len(myMatrix), trgSkinMesh.getNumberOfNodes())
                # tested values are constructed for this test
                for irow, row in enumerate(myMatrix):
                    self.assertTrue(
                        abs(sum([y for x, y in row.items()]) - 1.0) < 1e-8,
                        f"{irow}: {row}",
                    )

                # Transfer field src -> trg
                trgFieldCV = remap.transferField(srcField, 1e300)
                self.assertTrue(
                    trgFieldCV.getArray().isEqual(trgField.getArray(), 0.04)
                )

                # Transfer field trg -> src (should not work since need two matrices)
                try:
                    srcFieldCV = remap.transferField(trgField, 2e300)
                except InterpKernelException as e:
                    self.assertTrue(
                        "MEDCouplingRemapper::transferUnderground" in e.what()
                    )
                else:
                    self.assertTrue(False, "")

    def testFEFE_PROJ3DSurf2(self):
        """
        [EDF35052] : Test to project a 3D field using standard reference FE elements.
                     Comparaison with PROJ_CHAMP of code_aster
        """

        # fmt: off
        def getSrcField():
            coo = DataArrayDouble( [10,0,10,10,0,0,10,0,5.0000000000849765,8.8545602565321282,4.6472317204376292,0,5.6806474673116636,8.229838658936492,0,1.2053668025534483,9.9270887409805137,0,-3.5460488704250497,9.3501624268542649,0,-7.485107481710779,6.6312265824082139,0,-9.7094181742604402,2.3931566428759043,0,-9.7094181742605858,-2.3931566428753115,0,-7.4851074817111876,-6.631226582407753,0,-3.546048870425575,-9.350162426854066,0,1.2053668025530646,-9.9270887409805599,0,5.6806474673114602,-8.2298386589366306,0,8.8545602565320749,-4.6472317204377314,0,8.8545602565321282,4.6472317204376292,10,5.6806474673116636,8.229838658936492,10,1.2053668025534483,9.9270887409805137,10,-3.5460488704250497,9.3501624268542649,10,-7.485107481710779,6.6312265824082139,10,-9.7094181742604402,2.3931566428759043,10,-9.7094181742605858,-2.3931566428753115,10,-7.4851074817111876,-6.631226582407753,10,-3.546048870425575,-9.350162426854066,10,1.2053668025530646,-9.9270887409805599,10,5.6806474673114602,-8.2298386589366306,10,8.8545602565320749,-4.6472317204377314,10,7.7811633976032013,6.2812018101474978,4.9999999981719911,3.6604911852851965,9.305955312725791,4.9999999869607041,-1.1708274595790373,9.9312216297843072,4.9999998730442723,-5.6724724317290667,8.2354754878679426,4.9999997339800499,-8.8533840633895409,4.6494720803678486,4.9999988289974473,-9.9999999919939455,-0.00040015132076599684,4.9999960998526829,-8.8531781926891302,-4.6498640720448616,4.999995531092754,-5.6722873304897901,-8.235602979768089,4.9999980638359167,-1.1707087909851739,-9.9312356193330764,4.9999993701841738,3.6606551656336506,-9.3058908094990915,4.9999988475357213,7.7813069052198003,-6.2810240285146701,4.9999987767825651,5.5340823710774334,1.3487841907658067,0,3.2250863175388838,5.127626938559553,0,-2.1803139901376309,5.2950976220419328,0,-5.5332569132108107,1.2328717980453225,0,-4.3776376627149798,-3.780656132560102,0,0.69809907950161454,-5.808735198191747,0,5.1916325368070906,-3.1217986645435594,0,0.20716076262543112,0.1522685834023709,0,5.5292816307011892,1.3558915370920468,10,3.1640828217357702,5.1077914499421224,10,-2.1109563050631692,5.2681532859939226,10,-5.5254983067041481,1.2430021172291135,10,-4.3725138131205572,-3.7677495656097046,10,0.69678623562337327,-5.7916330222385186,10,5.1863237999917118,-3.1141374383478566,10,0.2148336697050659,0.2364258291455468,10], 54, 3 )
            conn = DataArrayInt( [3,1,2,3,3,2,0,15,3,37,2,14,3,0,2,26,3,2,15,27,3,3,2,27,3,15,16,27,3,4,3,27,3,5,4,28,3,4,27,28,3,6,5,29,3,5,28,29,3,7,6,30,3,6,29,30,3,8,7,31,3,7,30,31,3,9,8,32,3,8,31,32,3,10,9,33,3,9,32,33,3,11,10,34,3,10,33,34,3,12,11,35,3,11,34,35,3,13,12,36,3,12,35,36,3,14,13,37,3,13,36,37,3,26,2,37,3,1,14,2,3,16,17,28,3,17,18,29,3,18,19,30,3,19,20,31,3,20,21,32,3,21,22,33,3,22,23,34,3,23,24,35,3,24,25,36,3,25,26,37,3,35,24,36,3,34,23,35,3,33,22,34,3,27,16,28,3,32,21,33,3,30,19,31,3,29,18,30,3,28,17,29,3,20,32,31,3,25,37,36,3,1,3,38,3,3,4,39,3,38,3,39,3,4,5,39,3,5,6,40,3,39,5,40,3,6,7,40,3,7,8,41,3,40,7,41,3,8,9,41,3,9,10,42,3,41,9,42,3,10,11,42,3,11,12,43,3,42,11,43,3,12,13,43,3,13,14,44,3,43,13,44,3,14,1,44,3,1,38,44,3,38,39,45,3,39,40,45,3,40,41,45,3,41,42,45,3,42,43,45,3,44,38,45,3,45,43,44,3,15,0,46,3,16,15,47,3,15,46,47,3,17,16,47,3,18,17,48,3,17,47,48,3,19,18,48,3,20,19,49,3,21,20,49,3,19,48,49,3,22,21,50,3,21,49,50,3,23,22,50,3,24,23,51,3,25,24,51,3,23,50,51,3,26,25,52,3,0,26,52,3,25,51,52,3,0,52,46,3,48,47,53,3,47,46,53,3,49,48,53,3,50,49,53,3,51,50,53,3,46,52,53,3,52,51,53] )
            connI = DataArrayInt( [0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188,192,196,200,204,208,212,216,220,224,228,232,236,240,244,248,252,256,260,264,268,272,276,280,284,288,292,296,300,304,308,312,316,320,324,328,332,336,340,344,348,352,356,360,364,368,372,376,380,384,388,392,396,400,404,408,412,416] )
            src_mesh = MEDCouplingUMesh("src",2) ; src_mesh.setCoords( coo )
            src_mesh.setConnectivity( conn, connI, True )
            arr = DataArrayDouble( [20.0, 10.0, 15.000000000084977, 13.501791976969757, 13.910486126248156, 11.132455543533961, 5.804113556429215, -0.853880899302565, -7.316261531384535, -12.102574817135897, -14.116334064118941, -12.896211297279642, -8.721721938427496, -2.5491911916251704, 4.2073285360943435, 23.501791976969756, 23.910486126248156, 21.13245554353396, 15.804113556429215, 9.146119100697435, 2.6837384686154646, -2.1025748171358973, -4.116334064118941, -2.896211297279642, 1.2782780615725038, 7.45080880837483, 14.207328536094344, 19.06236520592269, 17.96644648497169, 13.760394043249542, 7.563002790118926, 0.7960868459757551, -5.000404043462029, -8.503046733641238, -8.907892246421962, -6.101945040134076, -0.6452367963297192, 6.500281653487695, 6.88286656184324, 8.352713256098436, 3.114783631904302, -4.300385115165488, -8.158293795275082, -5.110636118690133, 2.069833872263531, 0.359429346027802, 16.885173167793237, 18.271874271677895, 13.157196980930753, 5.717503810524965, 1.8597366212697377, 4.905153213384855, 12.072186361643855, 10.451259498850613] )
            arr.setInfoOnComponents( ["X1"] )
            src_field = MEDCouplingFieldDouble( ON_NODES ) ; src_field.setMesh( src_mesh ) ; src_field.setArray( arr )
            src_field.setNature(IntensiveMaximum)
            return src_field

        def getTrgField():
            coo = DataArrayDouble( [10,0,10,10,0,0,10,0,2.5000000000000306,10,0,5.000000000000048,10,0,7.5000000000000089,9.6858316112863143,2.4868988716485303,0,8.7630668004386543,4.817536741017121,0,7.2896862742141533,6.8454710592868464,0,5.3582679497900294,8.4432792550201103,0,3.0901699437495549,9.5105651629515098,0,0.62790519529321731,9.9802672842827107,0,-1.8738131458571674,9.8228725072869025,0,-4.2577929156506551,9.0482705246602286,0,-6.3742398974868415,7.7051324277579383,0,-8.0901699437494319,5.8778525229247887,0,-9.2977648588824895,3.6812455268468436,0,-9.9211470131447719,1.2533323356431025,0,-9.9211470131447861,-1.2533323356429857,0,-9.2977648588825375,-3.6812455268467215,0,-8.0901699437495012,-5.8778525229246936,0,-6.3742398974869365,-7.7051324277578592,0,-4.2577929156507777,-9.0482705246601718,0,-1.8738131458572813,-9.8228725072868812,0,0.62790519529310185,-9.9802672842827178,0,3.0901699437494474,-9.5105651629515435,0,5.3582679497899521,-8.4432792550201601,0,7.2896862742141071,-6.8454710592868961,0,8.7630668004386312,-4.817536741017161,0,9.6858316112863072,-2.4868988716485623,0,9.6858316112863143,2.4868988716485303,10,8.7630668004386543,4.817536741017121,10,7.2896862742141533,6.8454710592868464,10,5.3582679497900294,8.4432792550201103,10,3.0901699437495549,9.5105651629515098,10,0.62790519529321731,9.9802672842827107,10,-1.8738131458571674,9.8228725072869025,10,-4.2577929156506551,9.0482705246602286,10,-6.3742398974868415,7.7051324277579383,10,-8.0901699437494319,5.8778525229247887,10,-9.2977648588824895,3.6812455268468436,10,-9.9211470131447719,1.2533323356431025,10,-9.9211470131447861,-1.2533323356429857,10,-9.2977648588825375,-3.6812455268467215,10,-8.0901699437495012,-5.8778525229246936,10,-6.3742398974869365,-7.7051324277578592,10,-4.2577929156507777,-9.0482705246601718,10,-1.8738131458572813,-9.8228725072868812,10,0.62790519529310185,-9.9802672842827178,10,3.0901699437494474,-9.5105651629515435,10,5.3582679497899521,-8.4432792550201601,10,7.2896862742141071,-6.8454710592868961,10,8.7630668004386312,-4.817536741017161,10,9.6858316112863072,-2.4868988716485623,10,9.6858316112863143,2.4868988716485307,2.5000000000000138,9.6858316112863143,2.4868988716485307,5.0000000000000195,9.6858316112863161,2.4868988716485294,7.500000000000008,9.6858316112863072,-2.4868988716485623,2.5000000000000187,9.6858316112863072,-2.4868988716485623,5.0000000000000338,9.6858316112863072,-2.4868988716485623,7.5000000000000089,8.7630668004386525,4.8175367410171219,2.5000000000000062,7.2896862742141515,6.8454710592868491,2.5000000000000031,5.3582679497900241,8.4432792550201139,2.5000000000000018,3.0901699437495491,9.5105651629515116,2.5000000000000018,0.62790519529321287,9.9802672842827107,2.5000000000000009,-1.8738131458571674,9.8228725072869025,2.5000000000000013,-4.2577929156506551,9.0482705246602286,2.5000000000000013,-6.3742398974868388,7.7051324277579401,2.5000000000000013,-8.0901699437494319,5.8778525229247887,2.5000000000000018,-9.2977648588824895,3.6812455268468436,2.5000000000000018,-9.9211470131447701,1.253332335643107,2.5000000000000027,-9.9211470131447861,-1.2533323356429769,2.5000000000000036,-9.2977648588825375,-3.6812455268467215,2.5000000000000036,-8.0901699437495065,-5.8778525229246874,2.5000000000000044,-6.3742398974869365,-7.7051324277578592,2.5000000000000031,-4.2577929156507697,-9.0482705246601753,2.5000000000000018,-1.87381314585729,-9.8228725072868777,2.5000000000000027,0.62790519529309297,-9.9802672842827178,2.5000000000000031,3.0901699437494474,-9.5105651629515435,2.5000000000000031,5.3582679497899459,-8.4432792550201636,2.5000000000000036,7.2896862742141,-6.8454710592869024,2.5000000000000058,8.7630668004386276,-4.817536741017169,2.5000000000000107,8.7630668004386525,4.8175367410171219,7.5000000000000053,7.2896862742141515,6.8454710592868491,7.5000000000000027,5.3582679497900223,8.4432792550201157,7.5000000000000018,3.0901699437495465,9.5105651629515116,7.5000000000000009,0.62790519529321287,9.9802672842827107,7.5,-1.8738131458571696,9.8228725072869025,7.5000000000000009,-4.2577929156506586,9.0482705246602269,7.5000000000000009,-6.3742398974868388,7.7051324277579401,7.5000000000000018,-8.0901699437494319,5.8778525229247887,7.5000000000000009,-9.2977648588824895,3.6812455268468436,7.5000000000000053,-9.9211470131447701,1.253332335643107,7.5000000000000036,-9.9211470131447861,-1.2533323356429769,7.5000000000000027,-9.2977648588825375,-3.6812455268467215,7.5000000000000044,-8.0901699437495065,-5.8778525229246874,7.5000000000000053,-6.3742398974869365,-7.7051324277578592,7.5000000000000036,-4.2577929156507697,-9.0482705246601753,7.5000000000000027,-1.87381314585729,-9.8228725072868777,7.5000000000000044,0.62790519529309297,-9.9802672842827178,7.5000000000000044,3.0901699437494474,-9.5105651629515435,7.5000000000000036,5.3582679497899521,-8.4432792550201601,7.5000000000000044,7.2896862742141071,-6.8454710592868961,7.500000000000008,8.7630668004386312,-4.817536741017161,7.500000000000016,7.2896862742141497,6.8454710592868491,5.0000000000000044,8.7630668004386525,4.8175367410171219,5.000000000000008,8.7630668004386276,-4.817536741017169,5.0000000000000178,5.3582679497900223,8.4432792550201157,5.0000000000000018,3.0901699437495447,9.5105651629515133,5.0000000000000009,0.62790519529321065,9.9802672842827107,5.0000000000000009,-1.8738131458571696,9.8228725072869025,5,-4.2577929156506586,9.0482705246602269,5.0000000000000009,-6.3742398974868388,7.7051324277579401,5.0000000000000009,-8.0901699437494319,5.8778525229247887,5.0000000000000018,-9.2977648588824895,3.6812455268468436,5.0000000000000062,-9.9211470131447701,1.253332335643107,5.0000000000000044,-9.9211470131447861,-1.2533323356429769,5.0000000000000036,-9.2977648588825375,-3.6812455268467215,5.0000000000000053,-8.0901699437495065,-5.8778525229246874,5.0000000000000133,-6.3742398974869365,-7.7051324277578592,5.0000000000000053,-4.2577929156507697,-9.0482705246601753,5.0000000000000044,-1.87381314585729,-9.8228725072868777,5.0000000000000053,0.62790519529309297,-9.9802672842827178,5.0000000000000027,3.0901699437494385,-9.5105651629515471,5.0000000000000027,5.3582679497899379,-8.443279255020169,5.0000000000000053,7.2896862742141,-6.8454710592869024,5.0000000000000098,7.2095270388828405,1.6808728117093512,0,7.9763081628945178,-0.27810726502905703,0,6.1086302157897796,3.9964997206538495,0,4.0539722994132834,6.3201009564727997,0,2.4521123970423151,7.5144410149360743,0,0.27303102321078415,7.7231513893245261,0,-1.7880987132448045,7.9331976930326533,0,-3.3227204034723452,7.6939160933452175,0,-4.4513748898922616,6.5231665015157931,0,-6.4457171910796989,2.183550875279674,0,-4.9064971728924593,4.074013802060934,0,-7.3566886740400275,-0.0096034496729340418,0,-6.9644432453634169,-2.8701233774687616,0,-6.4406079582105642,-4.7064475741415972,0,-4.8261667575663578,-6.1516191757087251,0,-3.2184954564163037,-7.4415115519258235,0,-1.8328247568853897,-8.1279643493735083,0,-0.20514259866412757,-7.8356596446710665,0,4.0085965872335194,-5.5742776741746409,0,1.6741845163786164,-6.1068021841980284,0,6.0413511505972881,-4.3763745942353198,0,7.5136822221347002,-1.9987275262158615,0,-3.2637789840303704,5.1255778472766247,0,-0.45474739796773173,-5.9681218968928285,0,3.5494531022935583,2.4800721871736795,0,5.152376792496292,0.47451130406962427,0,0.11819581857869121,5.1484974110223778,0,2.581436880340211,5.0145208684583951,0,-1.8595680568275261,6.3210879792869319,0,-4.9660183199407095,-2.7729446336920458,0,-5.1054778182523126,-0.73603776829594225,0,-4.3701648120050942,0.93605282550997815,0,-3.0579468853463783,2.3173890659791994,0,-3.3793805162051473,-4.7278558532409569,0,-2.3076274324330983,-6.2664860721772397,0,6.0484022490484852,-1.2315926311161065,0,4.7489098239411627,-2.4971358482901804,0,3.0750009815501365,-3.4607900782459704,0,1.1847714333729324,-3.5702814418981892,0,-1.6486730049443152,3.6096730959599497,0,-1.5001799518335719,-3.5802188795623366,0,-3.2600184643838261,-2.3035217803814581,0,-3.0269974635482728,-0.051805425777902675,0,1.4014994623832955,2.674493675946279,0,3.7879538294403994,-0.91140602824714978,0,2.6262086039668935,-1.9557853082416179,0,2.496055247626094,0.41206241371634572,0,-1.4586080227731983,0.51758123580053939,0,0.44453036466952589,-1.4799281721395645,0,-0.38651905760367866,1.898585578962819,0,-1.1970079620378002,-1.4532724043458576,0,0.21490069809469703,-0.059950956309141511,0,1.6792943125839284,-0.94912714466517378,0,0.98150709044035256,1.2359068998856566,0,-2.4955302526328418,-1.2678798353267142,0,7.6569064607699051,-0.18470542218585201,10,7.6409983731915281,1.9059435579894441,10,6.646918880245078,3.4687667926584327,10,5.4883254608184098,5.0891485566898611,10,3.921681421651344,6.2429305650829292,10,2.1195048377507399,6.9876462867604907,10,0.22921688816265734,7.3118553780245472,10,-1.6552552386328752,7.1435622049320022,10,-3.4506422016708824,6.5542577169772729,10,-5.056626749474268,5.5690135021684251,10,-6.3200251764295627,4.1589338094415256,10,-7.1969699769132625,2.4621640263188969,10,-7.649042600888702,0.58302900533194324,10,-7.6020725180550368,-1.3592569447173222,10,-7.0754733432354939,-3.2558112521058216,10,-6.1512015249435796,-5.0036892687111285,10,-4.9355815362984199,-6.4816563960190239,10,-3.375608769251349,-7.4378584406804897,10,-1.4783808261309446,-7.9484180324595588,10,0.57047614642603928,-7.9734725499616523,10,2.6548028034180531,-7.5656056145338875,10,4.5685563308463326,-6.7687843022710297,10,6.1234631966850053,-5.5455541278600169,10,7.2650777563727313,-4.0591659230319879,10,7.6294847934738943,-2.2556990391621055,10,5.0108331686110787,2.1323710464892347,10,4.1924016319898669,3.8212027610254053,10,2.9150908111077491,4.5482511350952723,10,1.4817838361927602,5.0977088851546197,10,0.044417999675292648,5.4327915002796905,10,-1.3877766212144671,5.2069623240742686,10,-2.7953705223308027,4.7113938496820733,10,-4.0418839227342991,4.0802024233453498,10,-4.8982205027140404,2.969192020888654,10,-5.5125039721719187,1.6237752575282207,10,-5.8760835542537322,0.15963171829784051,10,-5.7689778546245991,-1.3403315900343864,10,-5.2831424705481433,-2.8442187863411479,10,-4.5691421095719793,-4.2919791663481561,10,-3.8881242737258881,-5.5656652162713272,10,-2.798087179035698,-6.1078215658552395,10,-1.2434674097378302,-6.3923326579460991,10,0.41535513410622527,-6.2342199219408476,10,2.2388890700215569,-5.8613445426910635,10,3.8719084951175993,-5.3778926211103188,10,5.0901072707632586,-4.4381343717106398,10,6.0665622740080529,-3.597241451087128,10,5.7936694483533069,-2.2817593354536472,10,5.2273164684384872,-0.39924506458064235,10,2.2891075656865496,1.3655255661765742,10,2.3443836996398422,3.210546844049174,10,0.94415190662883219,3.5982539597571592,10,-0.020002473393825545,4.23403796272083,10,-1.0517759948923793,3.685473458843413,10,-3.3332119455567955,3.1693626307111828,10,-3.6692369467848036,2.2251518897046347,10,-4.0995141325682924,1.0531271605189785,10,-4.5704987089897804,-0.096584337336179535,10,-4.3146150972212656,-1.1864013402474334,10,-3.7132339532714815,-2.3887037225788244,10,-2.940072272515073,-3.7375497621147336,10,-2.4840833458471119,-5.0339994623994198,10,-1.2126941927181307,-5.1035483857587547,10,0.063872815239428588,-4.624375757006268,10,1.8857215732990711,-4.1521204187518226,10,3.3192708968909539,-4.2308889729896988,10,4.1481080839155346,-3.1431911079312371,10,2.0227232073463695,-1.7040021641585059,10,-0.20725475852660161,2.0405094819967169,10,-0.82836907341101329,-3.2224356055663987,10,-2.5985542919983851,0.44370428306557347,10,-3.2428460064560412,-0.73178432008425986,10,-2.2196834926971478,-1.7061392359930285,10,-2.2316379891803186,1.9398763330546598,10,-2.3152715370388997,3.3193264272542566,10,-0.84605190699064248,-0.27020027104450139,10,10,0,1.2500000000000153,10,0,3.7500000000000391,10,0,6.2500000000000284,10,0,8.7500000000000036,9.921147013144779,1.2533323356430335,0,9.2977648588825232,3.6812455268467548,0,8.0901699437494994,5.8778525229246945,0,6.3742398974869472,7.7051324277578503,0,4.257792915650799,9.0482705246601611,0,1.8738131458573302,9.8228725072868706,0,-0.627905195293052,9.9802672842827214,0,-3.0901699437493972,9.5105651629515613,0,-5.3582679497899015,8.4432792550201921,0,-7.2896862742140645,6.8454710592869406,0,-8.7630668004386028,4.8175367410172143,0,-9.6858316112862948,2.4868988716486125,0,-10,1.2246467991473533e-15,0,-9.6858316112863267,-2.4868988716484899,0,-8.7630668004386632,-4.8175367410171033,0,-7.2896862742141515,-6.8454710592868482,0,-5.3582679497900152,-8.4432792550201192,0,-3.0901699437495185,-9.5105651629515222,0,-0.62790519529316746,-9.9802672842827143,0,1.8738131458572167,-9.8228725072868919,0,4.2577929156507022,-9.0482705246602073,0,6.3742398974868859,-7.705132427757901,0,8.0901699437494674,-5.8778525229247407,0,9.2977648588825073,-3.6812455268467947,0,9.9211470131447772,-1.2533323356430555,0,9.921147013144779,1.2533323356430335,10,9.2977648588825232,3.6812455268467548,10,8.0901699437494994,5.8778525229246945,10,6.3742398974869472,7.7051324277578503,10,4.257792915650799,9.0482705246601611,10,1.8738131458573302,9.8228725072868706,10,-0.627905195293052,9.9802672842827214,10,-3.0901699437493972,9.5105651629515613,10,-5.3582679497899015,8.4432792550201921,10,-7.2896862742140645,6.8454710592869406,10,-8.7630668004386028,4.8175367410172143,10,-9.6858316112862948,2.4868988716486125,10,-10,1.2246467991473533e-15,10,-9.6858316112863267,-2.4868988716484899,10,-8.7630668004386632,-4.8175367410171033,10,-7.2896862742141515,-6.8454710592868482,10,-5.3582679497900152,-8.4432792550201192,10,-3.0901699437495185,-9.5105651629515222,10,-0.62790519529316746,-9.9802672842827143,10,1.8738131458572167,-9.8228725072868919,10,4.2577929156507022,-9.0482705246602073,10,6.3742398974868859,-7.705132427757901,10,8.0901699437494674,-5.8778525229247407,10,9.2977648588825073,-3.6812455268467947,10,9.9211470131447772,-1.2533323356430555,10,9.6858316112863143,2.4868988716485307,1.2500000000000069,9.921147013144779,1.2533323356430337,2.5000000000000222,9.6858316112863143,2.4868988716485307,3.7500000000000169,9.921147013144779,1.2533323356430337,5.0000000000000338,9.6858316112863143,2.4868988716485307,6.2500000000000142,9.921147013144779,1.2533323356430328,7.5000000000000089,9.6858316112863143,2.4868988716485303,8.7500000000000036,9.6858316112863072,-2.4868988716485623,1.2500000000000093,9.9211470131447772,-1.2533323356430555,2.5000000000000249,9.6858316112863072,-2.4868988716485623,3.7500000000000262,9.9211470131447772,-1.2533323356430555,5.0000000000000409,9.6858316112863072,-2.4868988716485623,6.2500000000000213,9.9211470131447772,-1.2533323356430555,7.5000000000000089,9.6858316112863072,-2.4868988716485623,8.7500000000000036,9.2977648588825232,3.6812455268467548,2.5000000000000098,8.7630668004386543,4.817536741017121,1.2500000000000031,8.0901699437494994,5.8778525229246945,2.5000000000000044,7.2896862742141524,6.8454710592868473,1.2500000000000016,6.3742398974869428,7.7051324277578539,2.5000000000000027,5.3582679497900276,8.4432792550201121,1.2500000000000009,4.2577929156507928,9.0482705246601629,2.5000000000000018,3.0901699437495505,9.5105651629515116,1.2500000000000009,1.8738131458573237,9.8228725072868723,2.5000000000000013,0.62790519529321509,9.9802672842827107,1.2500000000000004,-0.62790519529305422,9.9802672842827214,2.5000000000000009,-1.8738131458571674,9.8228725072869025,1.2500000000000007,-3.0901699437493972,9.5105651629515613,2.5000000000000013,-4.2577929156506551,9.0482705246602286,1.2500000000000007,-5.3582679497899015,8.4432792550201921,2.5000000000000013,-6.3742398974868388,7.7051324277579401,1.2500000000000007,-7.2896862742140645,6.8454710592869406,2.5000000000000018,-8.0901699437494319,5.8778525229247887,1.2500000000000009,-8.7630668004386028,4.8175367410172143,2.5000000000000018,-9.2977648588824895,3.6812455268468436,1.2500000000000009,-9.6858316112862948,2.4868988716486125,2.5000000000000022,-9.9211470131447719,1.2533323356431025,1.2500000000000013,-10,6.3397136178156118e-14,2.5000000000000031,-9.9211470131447861,-1.2533323356429769,1.2500000000000018,-9.6858316112863285,-2.4868988716484814,2.5000000000000036,-9.2977648588825392,-3.681245526846717,1.2500000000000018,-8.7630668004386649,-4.8175367410170997,2.500000000000004,-8.0901699437495065,-5.8778525229246874,1.2500000000000022,-7.2896862742141515,-6.8454710592868482,2.5000000000000036,-6.3742398974869365,-7.7051324277578592,1.2500000000000016,-5.3582679497900081,-8.4432792550201246,2.5000000000000027,-4.2577929156507697,-9.0482705246601753,1.2500000000000009,-3.0901699437495185,-9.5105651629515222,2.5000000000000022,-1.87381314585729,-9.8228725072868777,1.2500000000000013,-0.62790519529317634,-9.9802672842827143,2.5000000000000027,0.62790519529309297,-9.9802672842827178,1.2500000000000016,1.8738131458572078,-9.8228725072868937,2.5000000000000031,3.0901699437494474,-9.5105651629515435,1.2500000000000016,4.2577929156507022,-9.0482705246602073,2.5000000000000036,5.3582679497899459,-8.4432792550201636,1.2500000000000018,6.3742398974868788,-7.7051324277579072,2.5000000000000044,7.2896862742141,-6.8454710592869024,1.2500000000000029,8.0901699437494621,-5.8778525229247469,2.500000000000008,8.7630668004386276,-4.817536741017169,1.2500000000000053,9.2977648588825073,-3.6812455268467947,2.5000000000000147,9.2977648588825232,3.6812455268467543,7.5000000000000071,8.7630668004386543,4.817536741017121,8.7500000000000036,8.0901699437494994,5.8778525229246945,7.5000000000000036,7.2896862742141524,6.8454710592868473,8.7500000000000018,6.3742398974869428,7.7051324277578548,7.5000000000000018,5.3582679497900259,8.4432792550201121,8.75,4.257792915650791,9.0482705246601647,7.5000000000000018,3.0901699437495505,9.5105651629515116,8.75,1.8738131458573237,9.8228725072868723,7.5,0.62790519529321509,9.9802672842827107,8.75,-0.62790519529305422,9.9802672842827214,7.5,-1.8738131458571674,9.8228725072869025,8.75,-3.0901699437494017,9.5105651629515595,7.5000000000000009,-4.2577929156506551,9.0482705246602286,8.75,-5.3582679497899015,8.4432792550201921,7.5000000000000018,-6.3742398974868388,7.7051324277579401,8.75,-7.2896862742140645,6.8454710592869406,7.5000000000000018,-8.0901699437494319,5.8778525229247887,8.75,-8.7630668004386028,4.8175367410172143,7.5000000000000036,-9.2977648588824895,3.6812455268468436,8.7500000000000036,-9.6858316112862948,2.4868988716486125,7.5000000000000044,-9.9211470131447719,1.2533323356431025,8.7500000000000018,-10,6.3397136178156118e-14,7.5000000000000036,-9.9211470131447861,-1.2533323356429769,8.7500000000000018,-9.6858316112863285,-2.4868988716484814,7.5000000000000036,-9.2977648588825392,-3.681245526846717,8.7500000000000018,-8.7630668004386649,-4.8175367410170997,7.5000000000000053,-8.0901699437495065,-5.8778525229246874,8.7500000000000036,-7.2896862742141515,-6.8454710592868482,7.5000000000000044,-6.3742398974869365,-7.7051324277578592,8.7500000000000018,-5.3582679497900081,-8.4432792550201246,7.5000000000000036,-4.2577929156507697,-9.0482705246601753,8.7500000000000018,-3.0901699437495185,-9.5105651629515222,7.5000000000000036,-1.87381314585729,-9.8228725072868777,8.7500000000000018,-0.62790519529317634,-9.9802672842827143,7.5000000000000044,0.62790519529309297,-9.9802672842827178,8.7500000000000018,1.8738131458572078,-9.8228725072868937,7.5000000000000036,3.0901699437494474,-9.5105651629515435,8.7500000000000018,4.2577929156507022,-9.0482705246602073,7.5000000000000036,5.3582679497899521,-8.4432792550201601,8.7500000000000018,6.3742398974868859,-7.705132427757901,7.5000000000000062,7.2896862742141071,-6.8454710592868961,8.7500000000000036,8.0901699437494674,-5.8778525229247407,7.5000000000000124,8.7630668004386312,-4.817536741017161,8.7500000000000071,9.2977648588825073,-3.6812455268467947,7.5000000000000124,8.0901699437494994,5.8778525229246963,5.0000000000000062,8.7630668004386525,4.8175367410171219,6.2500000000000071,7.2896862742141515,6.8454710592868491,6.2500000000000036,9.2977648588825232,3.6812455268467548,5.0000000000000142,8.7630668004386525,4.8175367410171219,3.7500000000000071,8.7630668004386276,-4.817536741017169,3.7500000000000142,9.2977648588825073,-3.6812455268467947,5.0000000000000258,8.7630668004386276,-4.817536741017169,6.2500000000000169,7.2896862742141515,6.8454710592868491,3.7500000000000036,6.3742398974869428,7.7051324277578548,5.0000000000000036,5.3582679497900223,8.4432792550201157,3.7500000000000018,4.2577929156507892,9.0482705246601665,5.0000000000000018,3.0901699437495465,9.5105651629515116,3.7500000000000013,1.8738131458573193,9.8228725072868723,5.0000000000000009,0.62790519529321287,9.9802672842827107,3.7500000000000009,-0.62790519529305644,9.9802672842827196,5,-1.8738131458571674,9.8228725072869025,3.7500000000000009,-3.0901699437494017,9.5105651629515595,5,-4.2577929156506551,9.0482705246602286,3.7500000000000009,-5.3582679497899015,8.4432792550201921,5.0000000000000009,-6.3742398974868388,7.7051324277579401,3.7500000000000009,-7.2896862742140645,6.8454710592869406,5.0000000000000018,-8.0901699437494319,5.8778525229247887,3.7500000000000018,-8.7630668004386028,4.8175367410172143,5.0000000000000036,-9.2977648588824895,3.6812455268468436,3.750000000000004,-9.6858316112862948,2.4868988716486125,5.0000000000000053,-9.9211470131447701,1.253332335643107,3.7500000000000036,-10,6.3397136178156118e-14,5.0000000000000036,-9.9211470131447879,-1.2533323356429724,3.7500000000000036,-9.6858316112863285,-2.4868988716484814,5.0000000000000044,-9.2977648588825392,-3.681245526846717,3.7500000000000044,-8.7630668004386649,-4.8175367410170997,5.0000000000000089,-8.0901699437495065,-5.8778525229246874,3.7500000000000089,-7.2896862742141515,-6.8454710592868482,5.0000000000000089,-6.3742398974869365,-7.7051324277578592,3.7500000000000044,-5.3582679497900081,-8.4432792550201246,5.0000000000000053,-4.2577929156507697,-9.0482705246601753,3.7500000000000031,-3.0901699437495185,-9.5105651629515222,5.0000000000000053,-1.87381314585729,-9.8228725072868777,3.750000000000004,-0.62790519529317634,-9.9802672842827143,5.0000000000000036,0.62790519529309297,-9.9802672842827178,3.7500000000000027,1.8738131458572078,-9.8228725072868937,5.0000000000000027,3.0901699437494385,-9.5105651629515471,3.7500000000000027,4.2577929156506942,-9.0482705246602109,5.0000000000000036,5.3582679497899379,-8.443279255020169,3.7500000000000044,6.3742398974868717,-7.7051324277579125,5.0000000000000071,7.2896862742141,-6.8454710592869024,3.750000000000008,8.0901699437494621,-5.8778525229247469,5.0000000000000142,5.3582679497900223,8.4432792550201157,6.2500000000000018,3.0901699437495465,9.5105651629515116,6.2500000000000009,0.62790519529321287,9.9802672842827107,6.25,-1.8738131458571696,9.8228725072869025,6.25,3.0901699437494385,-9.5105651629515471,6.2500000000000036,5.3582679497899459,-8.4432792550201636,6.2500000000000053,-6.3742398974868388,7.7051324277579401,6.2500000000000018,-8.0901699437494319,5.8778525229247887,6.2500000000000018,-9.2977648588824895,3.6812455268468436,6.2500000000000053,-9.9211470131447701,1.253332335643107,6.2500000000000036,-9.9211470131447879,-1.2533323356429724,6.2500000000000036,0.62790519529309297,-9.9802672842827178,6.2500000000000036,-8.0901699437495065,-5.8778525229246874,6.2500000000000089,-6.3742398974869365,-7.7051324277578592,6.2500000000000044,-4.2577929156507697,-9.0482705246601753,6.2500000000000036,-1.87381314585729,-9.8228725072868777,6.2500000000000053,-9.2977648588825392,-3.681245526846717,6.2500000000000053,7.2896862742141,-6.8454710592869024,6.2500000000000089,-4.2577929156506586,9.0482705246602269,6.2500000000000009,-3.1861038097865002,-5.2998323393353735,10,-2.641085262441405,-5.5709105141273296,10,-3.3431057263807933,-5.8367433910632833,10,-2.0207772943867641,-6.2500771119006693,10,-1.8483887692826213,-5.0687739240790872,10,-1.2280808012279805,-5.7479405218524269,10,-0.52665333275862203,0.8851546054761078,10,2.1559153865164595,-0.16923829899096576,10,0.5883356501778636,-0.9871012176015036,10,1.0409264035799741,1.7030175240866456,10,7.5929176008886792,0.70138277334014709,0,8.9881540814472594,-0.13905363251452851,0,8.447679325084577,2.0838858416789408,0,6.6590786273363101,2.8386862661816004,0,7.4358485081142174,4.4070182308354848,0,0.50749020238980846,2.2865396274545491,0,0.29749401641833695,1.5672462394242377,0,1.1915032764118241,1.9552002879159678,0,0.9470975053393127,-0.50453905048715764,0,0.32971553138211146,-0.76993956422435306,0,1.0619123386267271,-1.2145276584023692,0,1.3625717101265495,7.6187962021302997,0,2.7711411703959348,8.5125030889437916,0,0.45046810925200076,8.8517093368036193,0,-0.75753384501701015,7.8281745411785897,0,-1.8309559295509859,8.8780351001597779,0,-2.555409558358575,7.8135568931889354,0,-3.7902566595614999,8.371093309002724,0,-3.8870476466823032,7.1085412974305058,0,-5.4128073936895511,7.1141494646368653,0,-5.6761071819860796,3.128782338670304,0,-6.4983335583209456,4.9759331624928613,0,-7.8717410249810946,2.9323982010632585,0,-6.9012029325598636,1.0869737128033699,0,-8.6389178435923988,0.62186444298508425,0,1.7387811690332233,0.82398465680100119,0,2.4754762823384269,2.5772829315599792,0,3.0227541749598261,1.4460673004450126,0,2.0876747801050111,-0.26853236547441406,0,0.59820389426752474,0.58797797178825761,0,-5.633387357888461,-5.4290333749251616,0,-7.2653889509800322,-5.292150048533145,0,-5.6002033275266472,-6.9283758017332921,0,-4.0223311069913308,-6.7965653638172743,0,-3.7381441860335407,-8.2448910382929981,0,-2.5256601066508466,-7.7847379506496654,0,-1.8533189513713355,-8.9754184283301939,0,-1.0189836777747585,-7.9818119970222874,0,0.21138129831448715,-8.907963464476893,0,2.841390551806068,-5.8405399291863347,0,2.3821772300640318,-7.8086836735747855,0,4.6834322685117353,-7.0087784645974001,0,5.0249738689154038,-4.9753261342049804,0,6.6655187124056976,-5.610922826761108,0,1.9054900186699131,-2.7630333750699037,0,2.1298862074615346,-3.5155357600720798,0,2.850604792758515,-2.7082876932437943,0,7.744995192514609,-1.1384173956224592,0,8.5997569167105041,-2.242813198932212,0,0.73452095885724444,-6.9712309144345479,0,0.60971855920544238,-6.037462040545428,0,-0.32994499831592966,-6.9018907707819475,0,4.3509149473949247,1.4772917456216519,0,6.1809519156895663,1.0776920578894877,0,4.8290416590416694,3.2382859539137643,0,1.3498163494594511,5.0815091397403869,0,2.5167746386912633,6.2644809416972347,0,0.19561342089473768,6.4358244001734519,0,-0.8706861191244174,5.7347926951546544,0,-1.8238333850361652,7.127142836159793,0,-7.1605659597017226,-1.4398634135708477,0,-5.0357480690965115,-1.7544912009939939,0,-6.23108324614617,-0.37282060898443814,0,-5.9652307826520632,-2.8215340055804035,0,-3.714055848675736,1.6267209457445888,0,-3.982222029119419,3.1957014340200667,0,-5.4079410015423965,1.5598018503948261,0,-4.7378213151287039,0.100007528607018,0,-4.1726994180729289,-3.7504002434665011,0,-5.7033131390756369,-3.7396961039168213,0,-4.1027736368857521,-5.4397375144748406,0,-2.8435039743191228,-5.4971709627090988,0,-2.763061444424701,-6.8539988120515316,0,6.7775166863659937,-3.1875510602255908,0,5.398656036494824,-1.8643642397031435,0,5.3951304872692258,-3.4367552212627501,0,6.7810422355915927,-1.6151600786659839,0,1.4294779748757744,-4.838541813048109,0,3.5417987843918279,-4.5175338762103054,0,3.9119554027456496,-2.9789629632680752,0,2.1527514582754108,-1.4524562264533958,0,0.81465089902122911,-2.5251048070188769,0,5.0813012576015311,5.1583003385633246,0,3.0654449913168849,3.7472965278160375,0,3.3177045898767474,5.6673109124655969,0,-4.6789360313923609,5.298590151788364,0,-4.0851380784614149,4.5997958246687798,0,-3.857576936961316,5.8243721743962089,0,3.1420045385332465,-0.249671807265402,0,4.4701653109683459,-0.2184473620887627,0,-2.456225994487343,4.367625471618287,0,-2.3533099451453467,2.9635310809695747,0,-0.97746367490065178,-4.7741703882275823,0,-0.15770425923031972,-3.5752501607302629,0,-4.1827481413180694,-1.5197797743387,0,-4.1130183921622674,-2.5382332070367521,0,-2.8777743585083337,-1.7857008078540861,0,-4.0662376409002929,-0.39392159703692248,0,-2.7612638580905573,-0.65984263055230841,0,0.75984764048099329,3.9114955434843282,0,-1.0175960312739969,2.7541293374613844,0,-0.76523859318281195,4.379085253491164,0,-2.0702260946592439,-7.1972252107753736,0,-1.3811874152004151,-6.1173039845350345,0,-2.3800992081086991,-2.9418703299718976,0,-2.4397802340193597,-4.1540373664016466,0,-2.2428027431607358,0.23288790501131834,0,-2.2582774540597885,1.4174851508898694,0,-3.6985811377766833,0.44212369986603772,0,-2.5911442301499354,7.0075020363160743,0,-0.92256354018843845,1.2080834073816793,0,-1.348593956935686,-2.5167456419540972,0,-0.37623879868413723,-1.4666002882427112,0,-2.5616735204289482,5.7233329132817783,0,7.0123552059715015,-0.75484994807258177,0,5.6003895207723886,-0.37854066352324112,0,3.2070812167036467,-1.4335956682443838,0,7.4022089755179596,-4.59695566762624,0,-1.846269107335321,-1.3605761198362858,0,-0.62185366233925066,0.22881513974569895,0,-8.6389178435924059,-0.63146789265795988,0,-8.1311040521229767,-3.2756844521577415,0,-6.7025256017869905,-3.7882854758051794,0,6.699158245001966,5.4209853899703475,0,3.2530423482277993,6.9172709857044374,0,4.7061201246016564,7.381690105746455,0,4.2684318266907813,-1.7042709382686652,0,-2.7120778091810926,-4.3857746122570767,10,-1.020531633064572,-4.1629919956625763,10,-1.8842206729630431,-3.4799926838405661,10,-1.3278079924054993,-0.4678455842726591,0,-1.5240262830540805,-2.4642874207797139,10,-1.5328676998438953,-0.9881697535187649,10,0.59717706696767814,-2.4632188848624521,10,-2.4150961405893518,1.1917903080601167,10,-1.7223030994945137,0.086752006010536065,10,-1.2194463738534602,1.9901929075256883,10,-0.38224812908579231,-3.9234056812863334,10,-0.57441068873935108,-4.8639620713825114,10,-2.5553210296848512,4.0153601384681652,10,-3.687547934145547,3.6247825270282661,10,-3.4186272225325509,4.3957981365137115,10,-2.8242417412978478,3.2443445289827197,10,-3.9066723577229108,-0.41418432871021971,10,-3.7787305518386534,-0.95909283016584657,10,-4.4425569031055225,-0.64149283879180641,10,3.2683926658148543,3.5158748025372897,10,2.6297372553737954,3.8793989895722234,10,3.553746221548808,4.1847269480603391,10,7.6489524169807162,0.86061906790179599,10,8.6634149922389216,2.1964212148189874,10,8.8284532303849517,-0.092352711092926004,10,7.1439586267183035,2.6873551753239386,10,7.7049928403418662,4.1431517668377769,10,6.0676221705317435,4.2789576746741469,10,6.389005867516282,5.9673098079883538,10,4.7050034412348767,5.6660395608863947,10,4.6399746857206869,7.3431049100515198,10,3.0205931297010418,6.61528842592171,10,2.6048373907501476,8.2491057248560011,10,1.1743608629566986,7.1497508323925185,10,0.4285610417279373,8.6460613311536285,10,-0.71301917523510894,7.2277087914782747,10,-1.7645341922450213,8.4832173561094528,10,-2.5529487201518788,6.8489099609546376,10,-3.8542175586607685,7.8012641208187503,10,-4.253634475572575,6.061635609572849,10,-5.7154333234805552,6.6370729649631812,10,-5.6883259629519154,4.8639736558049753,10,-7.2050975600894969,5.0183931661831576,10,-6.758497576671413,3.3105489178802112,10,-8.2473674178978769,3.0717047765828704,10,-7.4230062889009822,1.5225965158254202,10,-8.785094807016737,0.91818067048752283,10,-7.6255575594718694,-0.38811396969268941,10,-8.7616097655999106,-1.3062946401801541,10,-7.3387729306452654,-2.3075340984115718,10,-8.1866191010590157,-3.4685283894762717,10,-6.6133374340895372,-4.1297502604084748,10,-7.1206857343465408,-5.4407708958179111,10,-5.5433915306209993,-5.7426728323650762,10,-5.6549107168926778,-7.0933944118884416,10,-4.1555951527748842,-6.9597574183497564,10,-3.8167008424510636,-8.2430644826703308,10,-2.4269947976911466,-7.6931382365700243,10,-1.6760969859941128,-8.8856452698732191,10,-0.45395233985245276,-7.9609452912106056,10,0.59919067085957056,-8.9768699171221851,10,1.6126394749220463,-7.7695390822477695,10,2.87248637358375,-8.5380853887427151,10,3.6116795671321928,-7.1671949584024581,10,4.9634121403181428,-7.6060317786455949,10,5.3460097637656689,-6.1571692150655233,10,6.7065747354495562,-6.1955125935734561,10,6.6942704765288683,-4.8023600254460028,10,8.0140722784056813,-4.438351332024574,10,7.4472812749233128,-3.1574324810970467,10,8.6576582023801016,-2.3712989554053339,10,7.6431956271218997,-1.2202022306739786,10,4.6016174003004728,2.9767869037573202,10,4.8403635464041379,4.455175658857633,10,5.8288760244280784,2.8005689195738337,10,3.4183861163795468,5.3955908500891008,10,2.1984373236502548,4.8229800101249456,10,1.8006443369717502,6.0426775859575557,10,0.76310091793402646,5.2652501927171551,10,0.13681744391897499,6.3723234391521189,10,-0.6716793107695872,5.31987691217698,10,-1.5215159299236711,6.1752622645031359,10,-2.0915735717726349,4.959178086878171,10,-3.1230063620008428,5.6328257833296735,10,-4.5492553361042836,4.8246079627568879,10,-4.4700522127241697,3.5246972221170019,10,-5.6091228395718016,3.5640629151650898,10,-5.20536223744298,2.2964836392084376,10,-6.3547369745425906,2.042969641923559,10,-5.6942937632128254,0.89170348791303067,10,-6.7625630775712171,0.37133036181489187,10,-5.8225307044391652,-0.59034993586827289,10,-6.6855251863398184,-1.3497942673758543,10,-5.5260601625863712,-2.092275188187767,10,-6.1793079068918182,-3.0500150192234847,10,-4.9261422900600618,-3.5680989763446522,10,-5.360171817257779,-4.6478342175296419,10,-4.2286331916489335,-4.9288221913097416,10,-4.4118529050121538,-6.0236608061451751,10,-3.0868479741435237,-6.7728400032678646,10,-1.3609241179343874,-7.1703753452028289,10,-0.41405613781580242,-6.3132762899434738,10,0.49291564026613227,-7.1038462359512504,10,1.327122102063891,-6.0477822323159556,10,2.446845936719805,-6.7134750786124755,10,3.0553987825695783,-5.6196185819006912,10,4.2202324129819662,-6.0733384616906747,10,4.4810078829404292,-4.9080134964104793,10,5.606785233724132,-4.9918442497853288,10,5.5783347723856558,-4.0176879113988839,10,6.6658200151903921,-3.8282036870595579,10,5.9301158611806795,-2.9395003932703876,10,6.711577120913601,-2.2687291873078763,10,5.1190748185247834,0.86656299095429623,10,6.3259157709013039,2.0191573022393392,10,6.4421114646041957,-0.29197524338324721,10,-1.2197763080534232,4.4462178914588408,10,-1.6835237659656395,3.5023999430488351,10,2.3167456326631957,2.2880362051128742,10,3.6499703671488142,1.7489483063329043,10,1.6442678031343372,3.4044004019031666,10,1.2129678714107963,4.3479814224558897,10,0.46207471661750332,3.9161459612389944,10,0.012207763140733553,4.8334147315002607,10,-0.53588923414310241,3.9597557107821215,10,-0.62951537670949043,2.8629914704200647,10,-2.2734547631096094,2.6296013801544582,10,-2.9207001492272129,-0.14404001850934323,10,-2.7312647495765945,-1.2189617780386441,10,-3.5012244461707995,2.697257260207909,10,-4.2837287247494222,2.5971719552966444,10,-3.8843755396765483,1.6391395251118066,10,-4.806009052370106,1.3384512090235996,10,-4.3350064207790364,0.47827141159139941,10,-5.2232911316217567,0.03152369048083048,10,-5.0417964759229328,-1.2633664651409098,10,-4.0139245252463738,-1.7875525314131289,10,-4.4981882119098122,-2.6164612544599861,10,-3.3266531128932773,-3.0631267423467792,10,-3.7546071910435259,-4.0147644642314448,10,1.9542223903227203,-2.9280612914551645,10,0.97479719426924982,-4.3882480878790453,10,0.23961397467282691,-5.4292978394735574,10,2.0623053216603138,-5.0067324807214426,10,2.6024962350950123,-4.1915046958707602,10,3.5955896960042768,-4.8043907970500088,10,3.7336894904032443,-3.687040040460468,10,4.6191076773393966,-3.7906627398209385,10,5.510492958395897,-1.3405022000171447,10,4.9708887661344203,-2.7124752216924422,10,3.6250198378924283,-1.0516236143695741,10,3.0854156456309521,-2.4235966360448717,10,-2.9504374679825611,2.0825141113796475,10,-3.349034212283339,0.748415721792276,10,-2.9664587229843145,-2.0474214792859264,10,0.36844857405111531,2.8193817208769381,10], 774,3 )
            conn = DataArrayInt( [6,220,219,241,483,481,482,6,178,168,174,498,496,497,6,173,177,176,501,499,500,6,162,170,163,537,535,536,6,148,142,144,542,540,541,6,147,133,135,578,576,577,6,154,166,155,586,585,552,6,140,141,159,526,593,563,6,131,132,153,507,600,550,6,146,126,160,538,605,567,6,26,27,145,282,608,534,6,16,17,136,272,611,515,6,6,7,127,262,614,495,6,128,129,152,615,547,575,6,137,138,154,613,560,554,6,156,155,167,558,588,599,6,243,242,249,629,619,628,6,238,237,251,636,634,635,6,207,206,230,639,637,638,6,182,181,205,643,732,692,8,1,2,53,5,256,311,310,260,8,2,3,54,53,257,313,312,311,8,3,4,55,54,258,315,314,313,8,4,0,29,55,259,285,316,315,8,2,1,28,56,256,284,317,318,8,3,2,56,57,257,318,319,320,8,4,3,57,58,258,320,321,322,8,0,4,58,52,259,322,323,309,8,6,5,53,59,261,310,324,325,8,7,6,59,60,262,325,326,327,8,8,7,60,61,263,327,328,329,8,9,8,61,62,264,329,330,331,8,10,9,62,63,265,331,332,333,8,11,10,63,64,266,333,334,335,8,12,11,64,65,267,335,336,337,8,13,12,65,66,268,337,338,339,8,14,13,66,67,269,339,340,341,8,15,14,67,68,270,341,342,343,8,16,15,68,69,271,343,344,345,8,17,16,69,70,272,345,346,347,8,18,17,70,71,273,347,348,349,8,19,18,71,72,274,349,350,351,8,20,19,72,73,275,351,352,353,8,21,20,73,74,276,353,354,355,8,22,21,74,75,277,355,356,357,8,23,22,75,76,278,357,358,359,8,24,23,76,77,279,359,360,361,8,25,24,77,78,280,361,362,363,8,26,25,78,79,281,363,364,365,8,27,26,79,80,282,365,366,367,8,28,27,80,56,283,367,368,317,8,29,30,81,55,286,370,369,316,8,30,31,82,81,287,372,371,370,8,31,32,83,82,288,374,373,372,8,32,33,84,83,289,376,375,374,8,33,34,85,84,290,378,377,376,8,34,35,86,85,291,380,379,378,8,35,36,87,86,292,382,381,380,8,36,37,88,87,293,384,383,382,8,37,38,89,88,294,386,385,384,8,38,39,90,89,295,388,387,386,8,39,40,91,90,296,390,389,388,8,40,41,92,91,297,392,391,390,8,41,42,93,92,298,394,393,392,8,42,43,94,93,299,396,395,394,8,43,44,95,94,300,398,397,396,8,44,45,96,95,301,400,399,398,8,45,46,97,96,302,402,401,400,8,46,47,98,97,303,404,403,402,8,47,48,99,98,304,406,405,404,8,48,49,100,99,305,408,407,406,8,49,50,101,100,306,410,409,408,8,50,51,102,101,307,412,411,410,8,51,52,58,102,308,323,413,412,8,81,82,103,104,371,416,414,415,8,59,53,54,104,324,312,417,418,8,57,56,80,105,319,368,419,420,8,58,57,105,102,321,420,421,413,8,60,59,104,103,326,418,414,422,8,61,60,103,106,328,422,423,424,8,62,61,106,107,330,424,425,426,8,63,62,107,108,332,426,427,428,8,64,63,108,109,334,428,429,430,8,65,64,109,110,336,430,431,432,8,66,65,110,111,338,432,433,434,8,67,66,111,112,340,434,435,436,8,68,67,112,113,342,436,437,438,8,69,68,113,114,344,438,439,440,8,70,69,114,115,346,440,441,442,8,71,70,115,116,348,442,443,444,8,72,71,116,117,350,444,445,446,8,73,72,117,118,352,446,447,448,8,74,73,118,119,354,448,449,450,8,75,74,119,120,356,450,451,452,8,76,75,120,121,358,452,453,454,8,77,76,121,122,360,454,455,456,8,78,77,122,123,362,456,457,458,8,79,78,123,124,364,458,459,460,8,80,79,124,105,366,460,461,419,8,55,81,104,54,369,415,417,314,8,83,84,107,106,375,463,425,462,8,84,85,108,107,377,464,427,463,8,85,86,109,108,379,465,429,464,8,122,99,100,123,466,407,467,457,8,88,89,112,111,385,469,435,468,8,89,90,113,112,387,470,437,469,8,90,91,114,113,389,471,439,470,8,91,92,115,114,391,472,441,471,8,98,99,122,121,405,466,455,473,8,94,95,118,117,397,475,447,474,8,95,96,119,118,399,476,449,475,8,96,97,120,119,401,477,451,476,8,97,98,121,120,403,473,453,477,8,116,93,94,117,478,395,474,445,8,100,101,124,123,409,479,459,467,8,101,102,105,124,411,421,461,479,8,82,83,106,103,373,462,423,416,8,92,93,116,115,393,478,443,472,8,110,87,88,111,480,383,468,433,8,86,87,110,109,381,480,431,465,8,221,220,241,242,484,482,485,486,8,255,248,229,247,487,490,488,489,8,1,5,125,126,260,493,491,492,8,5,6,127,125,261,495,494,493,8,9,10,130,129,265,504,502,503,8,10,11,131,130,266,506,505,504,8,11,12,132,131,267,508,507,506,8,12,13,133,132,268,510,509,508,8,14,15,134,135,270,513,511,512,8,15,16,136,134,271,515,514,513,8,178,171,149,168,516,518,517,498,8,177,171,178,176,519,516,520,499,8,19,20,139,138,275,523,521,522,8,20,21,140,139,276,525,524,523,8,21,22,141,140,277,527,526,525,8,22,23,142,141,278,529,528,527,8,24,25,143,144,280,532,530,531,8,25,26,145,143,281,534,533,532,8,28,1,126,146,284,492,538,539,8,142,23,24,144,529,279,531,540,8,125,127,149,150,494,545,543,544,8,129,130,151,152,502,548,546,547,8,130,131,153,151,505,550,549,548,8,136,137,154,155,551,554,552,553,8,135,134,156,157,511,557,555,556,8,134,136,155,156,514,553,558,557,8,138,139,158,154,521,561,559,560,8,139,140,159,158,524,563,562,561,8,145,146,160,161,564,567,565,566,8,144,143,162,163,530,569,536,568,8,143,145,161,162,533,566,570,569,8,170,177,173,163,571,501,572,535,8,127,128,152,149,573,575,574,545,8,171,169,150,149,579,580,543,518,8,135,157,164,147,556,582,581,577,8,144,163,165,148,568,584,583,541,8,166,179,167,155,587,589,588,585,8,168,151,164,174,590,592,591,496,8,141,142,148,159,528,542,594,593,8,154,158,165,166,559,596,595,586,8,157,156,167,172,555,599,597,598,8,157,172,174,164,598,601,591,582,8,163,173,175,165,572,603,602,584,8,132,133,147,153,509,578,604,600,8,126,125,150,160,491,544,606,605,8,170,169,171,177,607,579,519,571,8,27,28,146,145,283,539,564,608,8,166,165,175,179,595,602,609,587,8,133,13,14,135,510,269,512,576,8,172,176,178,174,610,520,497,601,8,17,18,137,136,273,612,551,611,8,18,19,138,137,274,522,613,612,8,8,9,129,128,264,503,615,616,8,7,8,128,127,263,616,573,614,8,149,152,151,168,574,546,590,517,8,161,160,150,169,565,606,580,617,8,161,169,170,162,617,607,537,570,8,151,153,147,164,549,604,581,592,8,241,240,249,242,618,620,619,485,8,158,159,148,165,562,594,583,596,8,172,167,179,175,597,589,609,621,8,172,175,173,176,621,603,500,610,8,255,247,249,252,489,624,622,623,8,250,253,248,255,625,627,487,626,8,211,254,234,212,630,633,631,632,8,29,0,180,181,285,642,640,641,8,30,29,181,182,286,641,643,644,8,31,30,182,183,287,644,645,646,8,32,31,183,184,288,646,647,648,8,33,32,184,185,289,648,649,650,8,34,33,185,186,290,650,651,652,8,35,34,186,187,291,652,653,654,8,36,35,187,188,292,654,655,656,8,37,36,188,189,293,656,657,658,8,38,37,189,190,294,658,659,660,8,39,38,190,191,295,660,661,662,8,40,39,191,192,296,662,663,664,8,41,40,192,193,297,664,665,666,8,42,41,193,194,298,666,667,668,8,43,42,194,195,299,668,669,670,8,44,43,195,196,300,670,671,672,8,45,44,196,197,301,672,673,674,8,46,45,197,198,302,674,675,676,8,47,46,198,199,303,676,677,678,8,48,47,199,200,304,678,679,680,8,49,48,200,201,305,680,681,682,8,50,49,201,202,306,682,683,684,8,51,50,202,203,307,684,685,686,8,52,51,203,204,308,686,687,688,8,0,52,204,180,309,688,689,642,8,183,182,205,206,645,692,690,691,8,184,183,206,207,647,691,639,693,8,185,184,207,208,649,693,694,695,8,186,185,208,209,651,695,696,697,8,187,186,209,210,653,697,698,699,8,188,187,210,211,655,699,700,701,8,189,188,211,212,657,701,632,702,8,190,189,212,213,659,702,703,704,8,191,190,213,214,661,704,705,706,8,192,191,214,215,663,706,707,708,8,193,192,215,216,665,708,709,710,8,194,193,216,217,667,710,711,712,8,195,194,217,218,669,712,713,714,8,196,195,218,219,671,714,715,716,8,197,196,219,220,673,716,483,717,8,198,197,220,221,675,717,484,718,8,199,198,221,222,677,718,719,720,8,200,199,222,223,679,720,721,722,8,201,200,223,224,681,722,723,724,8,202,201,224,225,683,724,725,726,8,203,202,225,226,685,726,727,728,8,204,203,226,227,687,728,729,730,8,181,180,228,205,640,733,731,732,8,210,233,254,211,734,735,630,700,8,206,205,229,230,690,737,736,637,8,208,207,230,231,694,638,738,739,8,209,208,231,232,696,739,740,741,8,210,209,232,233,698,741,742,734,8,233,248,253,254,743,627,744,735,8,251,250,255,252,745,626,623,746,8,213,212,234,235,703,631,747,748,8,214,213,235,236,705,748,749,750,8,215,214,236,237,707,750,751,752,8,216,215,237,238,709,752,636,753,8,217,216,238,239,711,753,754,755,8,218,217,239,240,713,755,756,757,8,219,218,240,241,715,757,618,481,8,247,244,243,249,758,759,628,624,8,222,221,242,243,719,486,629,760,8,223,222,243,244,721,760,759,761,8,224,223,244,245,723,761,762,763,8,225,224,245,246,725,763,764,765,8,204,227,228,180,730,766,733,689,8,246,227,226,225,767,729,727,765,8,205,228,247,229,731,768,488,737,8,247,228,227,246,768,766,767,769,8,254,253,235,234,744,770,747,633,8,247,246,245,244,769,764,762,758,8,237,236,250,251,751,771,745,634,8,239,238,251,252,754,635,746,772,8,240,239,252,249,756,772,622,620,8,248,233,232,231,743,742,740,773,8,250,236,235,253,771,749,770,625,8,248,231,230,229,773,738,736,490] )
            connI = DataArrayInt( [0,7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126,133,140,149,158,167,176,185,194,203,212,221,230,239,248,257,266,275,284,293,302,311,320,329,338,347,356,365,374,383,392,401,410,419,428,437,446,455,464,473,482,491,500,509,518,527,536,545,554,563,572,581,590,599,608,617,626,635,644,653,662,671,680,689,698,707,716,725,734,743,752,761,770,779,788,797,806,815,824,833,842,851,860,869,878,887,896,905,914,923,932,941,950,959,968,977,986,995,1004,1013,1022,1031,1040,1049,1058,1067,1076,1085,1094,1103,1112,1121,1130,1139,1148,1157,1166,1175,1184,1193,1202,1211,1220,1229,1238,1247,1256,1265,1274,1283,1292,1301,1310,1319,1328,1337,1346,1355,1364,1373,1382,1391,1400,1409,1418,1427,1436,1445,1454,1463,1472,1481,1490,1499,1508,1517,1526,1535,1544,1553,1562,1571,1580,1589,1598,1607,1616,1625,1634,1643,1652,1661,1670,1679,1688,1697,1706,1715,1724,1733,1742,1751,1760,1769,1778,1787,1796,1805,1814,1823,1832,1841,1850,1859,1868,1877,1886,1895,1904,1913,1922,1931,1940,1949,1958,1967,1976,1985,1994,2003,2012,2021,2030,2039,2048,2057,2066,2075,2084,2093,2102,2111,2120,2129,2138,2147,2156,2165,2174,2183,2192,2201,2210,2219,2228,2237,2246,2255,2264,2273,2282,2291,2300,2309,2318,2327,2336] )
            trg_mesh = MEDCouplingUMesh("src",2) ; trg_mesh.setCoords( coo )
            trg_mesh.setConnectivity( conn, connI, True )
            arr = DataArrayDouble( [20.0, 10.0, 12.500000000000032, 15.000000000000048, 17.50000000000001, 11.821617082428425, 13.52053573395015, 13.748153165959753, 13.69689643976962, 12.25837176185719, 10.508868985380936, 7.73186846917319, 4.759906137870237, 1.2813153731211318, -2.12372031146362, -5.504796246688495, -8.44453504958008, -10.95119972086625, -12.650432763719927, -13.70745834583674, -13.714311015177795, -13.080950465982154, -11.38417303660258, -9.209475826678371, -6.242945548738437, -3.0307699734700244, 0.43582905734183236, 3.9442284246025223, 6.986677289529576, 21.821617082428425, 23.515571682858024, 23.716170835651464, 23.68704600373678, 22.226146214083432, 20.494693370964278, 17.700864380461127, 14.741850864473692, 11.252035489189922, 7.854774226792365, 4.468084848923037, 1.5309379724511543, -0.9757266988350075, -2.677551653332913, -3.728963795863692, -3.7435908746453546, -3.099005724293316, -1.4151770805962018, 0.7763485793510357, 3.724829142924064, 6.959379663686013, 10.40384791553673, 13.939264557980716, 16.986677289529577, 14.26791049138478, 16.549026954334494, 19.323227986751434, 9.491018929366323, 11.926619561138022, 14.546329946168406, 15.803114141634904, 16.43739582951134, 16.08988955544173, 14.936258923250186, 12.956393947877872, 10.349722629907557, 7.216801168552728, 3.809812438176911, 0.2902546958542346, -3.07413442900073, -6.115378576862884, -8.584855220341465, -10.35435642428173, -11.335294956174486, -11.424149011029938, -10.664977356817912, -9.053751340175067, -6.749078480517647, -3.835145225555128, -0.5634935959428433, 2.939929684592694, 6.303032787658637, 20.841032275552514, 21.421052662353674, 21.108859558035718, 19.920777314208756, 17.971871168276184, 15.335542264009586, 12.230275843137601, 8.79736269469869, 5.301693519631359, 1.9155773715500164, -1.1063158607582635, -3.5758047646756608, -5.364644304235594, -6.323862767923233, -6.436598597256603, -5.651505734070854, -4.067931662717365, -1.733601693413526, 1.1493731950683967, 4.455475467194579, 7.923586534236313, 11.340941179696829, 18.917530826094314, 18.066831243135564, 8.721230794666406, 18.391361790022067, 17.463743547323595, 15.31299314997292, 12.87251436441844, 9.653547271070948, 6.355522209322096, 2.8224339676858143, -0.4593952533386838, -3.4918311902408004, -5.891922004660152, -7.726552249548187, -8.610383407227127, -8.825801888037466, -7.9544931050639835, -6.5001131607802325, -4.095161285377584, -1.300996733606524, 2.0011776905362684, 5.489956526199346, 8.890399850592193, 7.698200897865459, 10.10512993644363, 10.374073255886083, 9.966553411978389, 7.996182412535308, 6.145098979787848, 4.371195689872872, 2.0717916116235315, -4.262166315800025, -0.8324833708315246, -7.366292123712961, -9.834566622832176, -11.147055532352162, -10.977785933275085, -10.660007008342127, -9.960789106258899, -8.040802243335195, -1.565681086941122, -4.432617667819412, 1.6649765563619687, 5.514954695918839, 1.8617988632462539, -6.42286929486056, 6.029525289467237, 5.626888096565917, 5.266693229601069, 7.595957748798606, 4.461519922459406, -7.738962953632756, -5.841515586548255, -3.434111986495116, -0.740557819367179, -8.107236369446104, -8.57411350461034, 4.816809617932378, 2.251773975650982, -0.3857890966958337, -2.385510008525257, 1.961000091015634, -5.080398831395909, -5.563540244765284, -3.078802889326175, 4.075993138329575, 2.876547801193249, 0.6704232957252753, 2.9081176613424393, -0.9410267869726593, -1.0353978074700387, 1.5120665213591398, -2.6502803663836576, 0.15494974178555554, 0.7301671679187542, 2.217413990326009, -3.763410087959556, 17.472201038584053, 19.546941931180974, 20.11568567290351, 20.577474017508273, 20.16461198673428, 19.107151124511233, 17.541072266187204, 15.488306966299128, 13.103615515306391, 10.512386752694155, 7.838908633011962, 5.265194049405634, 2.9339864044432407, 1.0386705372276401, -0.33128459534131527, -1.1548907936547081, -1.4172379323174442, -0.8134672099318394, 0.5732011414094953, 2.597003596464387, 5.089197188884166, 7.799772028575303, 10.57790906882499, 13.205911833340743, 15.37378575431179, 17.143204215100315, 18.013604393015274, 17.463341946203023, 16.57949272134738, 15.477209499954984, 13.8191857028598, 11.91602332735127, 10.03831850061105, 8.070971518174613, 6.1112712853563025, 4.283548164044108, 2.890690555341014, 1.8726387431107094, 1.1388787240798643, 0.5462105100027844, 1.0940912551090611, 2.3641999323160707, 4.181135212165377, 6.377544527330494, 8.494015874007282, 10.65197289905262, 12.469320822920924, 13.51191011289966, 14.828071403857843, 13.654633131863125, 15.554930543689018, 14.542405866385991, 14.214035489327005, 12.633697463951034, 9.836150685154388, 8.55591494291983, 6.953613027950686, 5.33291695367404, 4.4989835625313015, 3.8980623241496932, 3.322377965370193, 2.481917191753469, 3.6837574215231137, 5.439497058233162, 7.733601154547248, 9.088381923901256, 11.004916975984296, 10.318721043187862, 11.833254723470116, 5.949195321022589, 7.845149991067188, 6.025369673459699, 6.074177271309823, 9.708238343874342, 11.004054890215357, 8.883747821964855, 11.250000000000014, 13.750000000000039, 16.25000000000003, 18.750000000000004, 10.904131550706047, 12.737987300252575, 13.638423705644284, 13.847993615758183, 12.986051200346736, 11.525334063740589, 9.130246768782884, 6.335787347687405, 3.029482224489188, -0.457024780683586, -3.8077454081751942, -7.188109321511137, -9.69356558948805, -12.140898807511608, -13.175378847605122, -14.0408341184303, -13.39265891577735, -12.478994148632335, -10.288574179800998, -7.860246456861957, -4.624612695661728, -1.3359063737445478, 2.205175205360994, 5.470815797449648, 8.504383736145362, 20.904131550706047, 22.737987300252573, 23.611185403505488, 23.82887229631576, 22.95592948935616, 21.509205296121046, 19.098539715966574, 16.323709539419625, 12.996687180715533, 9.535409989276621, 6.158809679596204, 2.809263475342636, 0.27277298736528527, -2.1435260092264437, -3.20882374161157, -4.048399342149275, -3.4254539321503885, -2.4910719394798155, -0.32028118688494467, 2.1236248953554595, 5.345265817152419, 8.644973017385047, 12.177937915431398, 15.470815797449648, 18.504383736145364, 13.071617082428432, 13.404131550706069, 15.359110395681942, 15.772307090645674, 17.936127470542964, 18.404131550706058, 20.57161708242843, 8.236677289529586, 11.004383736145387, 10.659469141747858, 13.492221374775552, 13.236474753653216, 16.00438373614537, 15.736677289529581, 15.030874845851185, 14.711914237337743, 16.327666369195867, 15.092774497735544, 16.358762449421747, 14.959344712627267, 15.663938361739731, 13.597315342553685, 14.040417135778226, 11.794853375332089, 11.748100929517246, 9.040795549540373, 8.833273108994312, 6.024154670398124, 5.557979289544963, 2.545563905649024, 2.05948716774053, -0.908389402619738, -1.3770835904874303, -4.289465337844613, -4.601625778120556, -7.271143709854373, -7.346782910036573, -9.777808381140531, -9.506231466445069, -11.502394594000826, -10.879302508166926, -12.559420176117637, -11.444022023137151, -12.569230013103867, -11.10249691162949, -11.935869463908222, -9.933870466141485, -10.218962188388828, -7.958152483373482, -8.059277946944844, -5.3427053661307, -5.039045387146783, -2.2168123724784192, -1.8256383909155802, 1.191539098320098, 1.6878793709672562, 4.709275832611859, 5.134582575277105, 7.886368769446436, 20.092111287233813, 22.228132791760988, 21.316067230207693, 22.56861174900257, 21.378534802653817, 22.467608442042547, 20.650560589481486, 21.073461764146096, 19.056352206850526, 19.30131017742782, 16.733217599515037, 16.518203322235358, 13.849206856783656, 13.514514467228073, 10.542014386224299, 10.024699091944306, 7.075952565098446, 6.578520488105851, 3.606302202223182, 3.1918311102365258, 0.4151923435603342, 0.2127447806813776, -2.363613396174199, -2.2939198906047733, -4.489417778714383, -4.021097978784251, -5.895916392514247, -5.072510121315025, -6.427558737185626, -5.090094735950979, -6.118461654761635, -4.445509585598938, -4.917937296225226, -2.7415543716567905, -2.9730357690061098, -0.5528218974315005, -0.32676850301101884, 2.4371011689962288, 2.769809869296753, 5.6826243555237195, 6.211313832668555, 9.163717224886522, 9.697676534130986, 12.650796372212024, 12.947595419960034, 18.783978987921085, 19.45393175934404, 20.276279710101136, 17.31791025481687, 16.894314045932063, 7.471483000040164, 10.327885034929633, 10.031085987181612, 17.76556076270995, 18.661037034640167, 17.220434398256188, 18.112758031943823, 16.275202503946684, 16.39747418854726, 14.11793452042366, 14.227403345416453, 11.65233284213692, 11.272478284717003, 8.419849989145973, 8.041252832667988, 5.083510268123398, 4.596693013152901, 1.5421584521955234, 1.0607261979588376, -1.781356725083497, -1.9703229859222, -4.81657179585984, -4.9999247730074305, -7.2513468285451, -6.805535018698874, -9.083216247678081, -8.516349327630522, -9.987020092452513, -8.714079376489524, -10.19418566586005, -8.689249661839751, -9.325057122433483, -7.220924667218355, -7.857167166600399, -5.499309704999985, -5.438879014090452, -2.688328094975076, -2.6312450639634895, 0.25013204968792013, 0.6986511990298793, 3.7570160560102464, 4.206380796219569, 7.131994737486016, 19.750110674028893, 18.768092864271416, 16.642432159124553, 14.160878650451766, -0.13835477885945024, 3.228326578865424, 7.595151985294872, 4.062063743658586, 0.7333624114025823, -2.299073525499531, -4.733863384667902, -2.9143814893955535, -7.46712308757518, -7.682541568385516, -6.802999419567419, -5.348619475283668, -6.568493629555939, 6.717105414528487, 10.941911557104277, 1.514063850878126, 1.7880042234312654, 0.8201508825559235, 1.7291455937125662, 3.0828373066382917, 3.023978676919593, 10.358501272717486, 11.986677087525495, 9.60123443257636, 12.74394392766662, 8.294300374228825, 8.84910044893273, 10.531565166763514, 9.49776489351791, 11.842866738949704, 2.794029829844357, 1.8647402558425743, 3.1467035643277916, 0.442558454852155, -0.44022403284224154, -0.15261531977564213, 8.981367912256848, 11.283644259339724, 9.302177446055618, 7.070640696161578, 7.047079170608793, 5.25814733483036, 4.580836649441224, 3.2214936507482026, 1.7013420709473137, -2.5473248433157747, -1.5224003958280838, -4.939342823917836, -5.814229219756494, -8.017053400607315, 2.5627658258342243, 5.0527592138984065, 4.468821475404837, 1.8191424146305977, 1.1861818660557821, -11.062420732813624, -12.55753899951318, -12.528579129259938, -10.818896470808605, -11.983035224326539, -10.310398057300512, -10.82873737970153, -9.000795674797047, -8.696582166162408, -2.999149377380267, -5.426506443510755, -2.325346196085666, 0.04964773471042383, 1.0545958856445896, -0.8575433563999901, -1.3856495526105457, 0.14231709951472127, 6.60657779689215, 6.356943717778291, -6.236709955577304, -5.427743481339984, -7.231835769097877, 5.828206693016576, 7.258643973579054, 8.067327612955433, 6.431325489199838, 8.781255580388496, 6.6314378210681895, 4.864106576030237, 5.303309451123628, -8.600429373272569, -6.790239270090506, -6.603903855130609, -8.786764788232466, -2.0873349029311465, -0.7865205950993515, -3.8481391511475698, -4.6378137865216855, -7.923099661539431, -9.443009242992456, -9.542511151360593, -8.340674937028222, -9.617060256476233, 3.589965626140403, 3.5342917967916803, 1.958375266006476, 5.165882156925608, -3.409063838172335, -0.9757350918184768, 0.9329924394775748, 0.7002952318220149, -1.7104539079976473, 10.239601596164855, 6.812741519132922, 8.985015502342344, 0.6196541203960033, 0.5146577462073654, 1.9667952374348925, 2.892332731267844, 4.251717948879583, 1.911399477130944, 0.6102211358242285, -5.751634063128233, -3.732954419960582, -5.70252791565677, -6.6512515991990195, -4.66347516636242, -4.460159237937216, -3.4211064886428657, 4.671343183965321, 1.736533306187387, 3.6138466603083526, -9.26745130543462, -7.498491399735451, -5.321969538080598, -6.593817600421005, -2.0099148381494167, -0.8407923031699194, -3.2564574379106457, 4.416357806166138, 0.285519867193241, -3.865339598889784, -1.8428390869268492, 3.16165939285283, 6.25750525789892, 5.221848857249146, 1.773485548459263, 2.8052533078917197, -3.2068452271716072, -0.39303852259355204, -9.270385736250365, -11.406788504280719, -10.49081107759217, 12.120143634972315, 10.170313333932235, 12.087810230348111, 2.5641608884221165, 2.902147578561831, 4.816476371272852, 4.635786643196392, -1.7956535766781587, 6.011686296166205, 7.478962546637339, 8.133958182105227, 8.776694167470765, 8.364448906516023, 10.770746533672227, 5.6943461896278755, 4.561627239878138, 11.460039108783313, 9.937234592882717, 10.97717091398116, 10.42010278768487, 5.679143313566868, 5.2621766179955, 4.91595025810267, 16.784267468352148, 16.50913624494602, 17.73847316960915, 18.509571484882514, 20.85983620705791, 18.736100519292027, 19.831313802042242, 21.848144607179645, 20.34657984520589, 22.35631567550464, 20.371043002121272, 21.983079595772207, 19.635881555622753, 20.85394311560615, 18.32411169534922, 19.074622372881567, 16.514689616243164, 16.718683163864434, 14.295961240802757, 13.947046562157979, 11.808001134000273, 10.921639641482624, 9.17564769285306, 7.813295606093661, 6.552051341208797, 4.8243373586849945, 4.099590226924438, 2.133085863470786, 1.9863284708354412, -0.06790440578006529, 0.35369297094316077, -1.655147490535287, -0.7430876944980116, -2.561456630164453, -1.2860643629860748, -2.7483051287811198, -1.1153525711246415, -2.059765325121395, -0.12013303426117039, -0.5617422558673335, 1.5851023689369421, 1.622320753737385, 3.843100392674277, 4.334400984841035, 6.444484608729735, 7.35738036167255, 9.188840548700146, 10.5110621418761, 11.891910451082866, 13.575720946381107, 14.289848793826266, 16.286359246974772, 16.42299339644792, 17.578404304057795, 19.295539205261772, 18.629444944001914, 18.813976966468648, 17.021417333775197, 17.843321922929306, 16.028351110651183, 16.509140883071094, 14.648197601407393, 14.653746334579466, 12.867604515105535, 12.50981942132883, 10.275352626652605, 9.054645009392832, 7.954940075593289, 7.091121401765458, 5.6882326673809676, 5.197409724700205, 3.6087672842436733, 3.5871193596925606, 1.9646805462843266, 2.3816646492258604, 0.7706770738846958, 1.5057587335952862, -0.008006034787419791, 0.8425446170413243, -0.4355137111573286, 0.14031202258861164, 1.468700536862784, 3.2726675722407244, 3.3890694043148812, 5.279339869747936, 5.73337085810733, 7.435780200668887, 8.146893951291291, 9.572994386529952, 10.614940983938803, 11.560646860986772, 12.837616328130835, 12.990615467910294, 14.442847933605725, 15.985637809479082, 18.345073073140647, 16.150136221220947, 13.226441583405418, 11.818876177083196, 14.604781837776072, 15.39891867348172, 15.048668205037504, 15.560949293866685, 14.378220677856497, 14.845622494640995, 13.423866476639018, 12.233476093710575, 10.356146617044848, 6.935259832263444, 6.049773472384761, 9.196032814037109, 8.313443230547222, 7.754763985435259, 6.532442156653495, 6.143264990812362, 4.808232558859072, 3.694837058936156, 4.198522943340497, 2.8853505336302017, 3.610220144759942, 2.230628344725029, 9.026161098867554, 6.586549106390205, 4.810316135199271, 7.055572840938872, 8.410991539224252, 8.791198898954267, 10.046649449942775, 10.828444937518459, 14.16999075837875, 12.258413544441977, 12.573396223522852, 10.661819009586079, 9.132076643397086, 7.399381509508937, 4.986119797729758, 13.187830294928053] )
            arr.setInfoOnComponents( ["X1"] )
            trg_field = MEDCouplingFieldDouble( ON_NODES ) ; trg_field.setMesh( trg_mesh ) ; trg_field.setArray( arr )
            trg_field.setNature(IntensiveMaximum)
            return trg_field

        # Source field
        srcField = getSrcField()

        # Target field
        trgField = getTrgField()
        # Interpolate nodes to nodes using FEM interpolation
        remap = MEDCouplingRemapper()
        remap.setIntersectionType(PointLocator)

        # Workaround since "FEFE" interpolation does not exist yet
        srcField2 = MEDCouplingFieldDouble(ON_NODES_FE)
        srcField2.setMesh(srcField.getMesh())
        srcField2.setArray(srcField.getArray())
        srcField2.setNature(IntensiveMaximum)
        srcField2.checkConsistencyLight()
        srcFt = MEDCouplingFieldTemplate(srcField2)
        trgFt = MEDCouplingFieldTemplate(ON_NODES_FE)
        trgFt.setMesh(trgField.getMesh())

        remap.prepareEx(srcFt, trgFt)

        trgFieldMC = remap.transferField(srcField2, 1e300)

        # Transfer field src -> trg
        self.assertTrue(trgFieldMC.getArray().isEqual(trgField.getArray(), 1e-8))
        # fmt: on

    def testP1P0OnHexa_1(self):
        """
        See EDF27859 : This test focuses on P1P0 interpolation with source containing HEXA. So P1P0 intersector is going to split into tetras
        the source cell.
        """
        trgMesh = MEDCouplingUMesh("mesh", 3)
        trgMesh.setCoords(
            DataArrayDouble(
                [
                    18500.0,
                    0.0,
                    0.0,
                    18544.0,
                    0.0,
                    0.0,
                    18544.0,
                    0.0,
                    200.0,
                    18500.0,
                    0.0,
                    200.0,
                    18497.96424104365,
                    274.44295777156043,
                    0.0,
                    18541.959399238567,
                    275.0956869684225,
                    0.0,
                    18541.959399238567,
                    275.0956869684225,
                    200.0,
                    18497.96424104365,
                    274.44295777156043,
                    200.0,
                ],
                8,
                3,
            )
        )
        firstPts = DataArrayDouble(3 * 10)
        firstPts[:] = 0.0
        firstPts.rearrange(3)
        trgMesh.setCoords(
            DataArrayDouble.Aggregate(firstPts, trgMesh.getCoords())
        )  # this line is important to check that correct ids are taken into account
        trgMesh.allocateCells(1)
        trgMesh.insertNextCell(NORM_HEXA8, [10, 11, 12, 13, 14, 15, 16, 17])

        srcMesh = trgMesh.deepCopy()
        cc = trgMesh.computeCellCenterOfMass()[0]
        trgMesh.scale(
            cc, 1.01
        )  # This line is to workaround the EDF28414 bug inside 3D intersector

        expectedMatrix0 = [
            {
                10: 503624.09065889206,
                11: 100868.41855508549,
                12: 503863.42469772906,
                13: 100629.0845162416,
                14: 100629.08451623631,
                15: 503863.4246977626,
                16: 100868.418555101,
                17: 503624.0906588909,
            }
        ]
        expectedMatrix1 = [
            {
                10: 604492.509213978,
                11: 201736.8371101737,
                12: 201736.83711016813,
                13: 201497.50307132734,
                14: 201258.16903247262,
                15: 201497.50307133005,
                16: 604492.5092140044,
                17: 201258.16903247265,
            }
        ]
        expectedMatrix2 = [
            {
                10: 302066.754077835,
                11: 302425.7551361466,
                12: 302425.7551361466,
                13: 302066.754077835,
                14: 302066.7540778395,
                15: 302425.7551361595,
                16: 302425.7551361595,
                17: 302066.75407783955,
            }
        ]
        for sp, expectedMatrix in [
            (PLANAR_FACE_5, expectedMatrix0),
            (PLANAR_FACE_6, expectedMatrix1),
            (GENERAL_24, expectedMatrix2),
        ]:
            remap = MEDCouplingRemapper()
            remap.setSplittingPolicy(sp)
            remap.prepare(srcMesh, trgMesh, "P1P0")
            mat = remap.getCrudeMatrix()
            self.checkMatrix(expectedMatrix, mat, 18, 1.0)

    def testP0P0OnMeshDim1SpaceDim3_0(self):
        """
        See EDF31137 : Management of P0P0 on meshes with meshdim == 1 and spacedim == 3
        """
        mS = MEDCouplingUMesh("", 1)
        mS.setCoords(DataArrayDouble([[0.5, 0.5, 0.5], [1.7, 1.7, 1.7]]))
        mS.allocateCells()
        mS.insertNextCell(NORM_SEG2, [0, 1])
        mT = MEDCouplingUMesh("", 1)
        mT.setCoords(DataArrayDouble([[1, 1, 1], [2, 2, 2 + 2e-13]]))
        mT.allocateCells()
        mT.insertNextCell(NORM_SEG2, [0, 1])
        rem = MEDCouplingRemapper()
        rem.prepare(mS, mT, "P0P0")
        mat = rem.getCrudeMatrix()
        expectedMatrix = [{0: 1.212435565298214}]
        self.checkMatrix(mat, expectedMatrix, 1, 1e-13)
        pass

    def checkMatrix(self, mat1, mat2, nbCols, eps):
        self.assertEqual(len(mat1), len(mat2))
        for i in range(len(mat1)):
            if len(mat2[i].keys()) > 0:
                self.assertTrue(max(mat2[i].keys()) < nbCols)
            if len(mat1[i].keys()) > 0:
                self.assertTrue(max(mat1[i].keys()) < nbCols)
            if len(mat2[i].keys()) > 0:
                self.assertTrue(min(mat2[i].keys()) >= 0)
            if len(mat1[i].keys()) > 0:
                self.assertTrue(min(mat1[i].keys()) >= 0)
            s1 = set(mat1[i].keys())
            s2 = set(mat2[i].keys())
            for elt in s1.intersection(s2):
                self.assertTrue(abs(mat1[i][elt] - mat2[i][elt]) < eps)
                pass
            for elt in s1.difference(s2):
                self.assertTrue(abs(mat1[i][elt]) < eps)
                pass
            for elt in s2.difference(s1):
                self.assertTrue(abs(mat2[i][elt]) < eps)
                pass
            pass
        pass

    def build2DSourceMesh_1(self):
        sourceCoords = [-0.3, -0.3, 0.7, -0.3, -0.3, 0.7, 0.7, 0.7]
        sourceConn = [0, 3, 1, 0, 2, 3]
        sourceMesh = MEDCouplingUMesh.New("my name of mesh 2D", 2)
        sourceMesh.allocateCells(2)
        sourceMesh.insertNextCell(NORM_TRI3, 3, sourceConn[0:3])
        sourceMesh.insertNextCell(NORM_TRI3, 3, sourceConn[3:6])
        sourceMesh.finishInsertingCells()
        myCoords = DataArrayDouble.New()
        myCoords.setValues(sourceCoords, 4, 2)
        sourceMesh.setCoords(myCoords)
        return sourceMesh

    def build2DTargetMesh_1(self):
        targetCoords = [
            -0.3,
            -0.3,
            0.2,
            -0.3,
            0.7,
            -0.3,
            -0.3,
            0.2,
            0.2,
            0.2,
            0.7,
            0.2,
            -0.3,
            0.7,
            0.2,
            0.7,
            0.7,
            0.7,
        ]
        targetConn = [0, 3, 4, 1, 1, 4, 2, 4, 5, 2, 6, 7, 4, 3, 7, 8, 5, 4]
        targetMesh = MEDCouplingUMesh.New()
        targetMesh.setMeshDimension(2)
        targetMesh.allocateCells(5)
        targetMesh.insertNextCell(NORM_QUAD4, 4, targetConn[0:4])
        targetMesh.insertNextCell(NORM_TRI3, 3, targetConn[4:7])
        targetMesh.insertNextCell(NORM_TRI3, 3, targetConn[7:10])
        targetMesh.insertNextCell(NORM_QUAD4, 4, targetConn[10:14])
        targetMesh.insertNextCell(NORM_QUAD4, 4, targetConn[14:18])
        targetMesh.finishInsertingCells()
        myCoords = DataArrayDouble.New()
        myCoords.setValues(targetCoords, 9, 2)
        targetMesh.setCoords(myCoords)
        return targetMesh

    def build2DTargetMesh_3(self):
        targetCoords = [
            -0.6,
            -0.4,
            -0.1,
            -0.4,
            1.1,
            -0.4,
            2.1,
            -0.4,
            -0.6,
            0.1,
            -0.1,
            0.1,
            1.1,
            0.1,
            2.1,
            0.1,
            -0.6,
            1.1,
            -0.1,
            1.1,
        ]
        targetConn = [0, 4, 5, 1, 1, 5, 6, 2, 2, 6, 7, 3, 4, 8, 9, 5]
        targetMesh = MEDCouplingUMesh.New()
        targetMesh.setMeshDimension(2)
        targetMesh.allocateCells(4)
        for i in range(4):
            targetMesh.insertNextCell(NORM_QUAD4, 4, targetConn[4 * i : 4 * (i + 1)])
            pass
        targetMesh.finishInsertingCells()
        myCoords = DataArrayDouble.New()
        myCoords.setValues(targetCoords, 10, 2)
        targetMesh.setCoords(myCoords)
        return targetMesh
        pass

    def setUp(self):
        pass

    pass


if __name__ == "__main__":
    unittest.main()
