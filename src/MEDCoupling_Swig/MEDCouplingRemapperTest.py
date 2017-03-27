#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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
        sourceMesh=self.build2DSourceMesh_1();
        targetMesh=self.build2DTargetMesh_1();
        remapper=MEDCouplingRemapper()
        remapper.setPrecision(1e-12);
        remapper.setIntersectionType(Triangulation);
        self.assertTrue(remapper.prepare(sourceMesh,targetMesh,"P0P0")==1);
        srcField=MEDCouplingFieldDouble.New(ON_CELLS);
        srcField.setNature(IntensiveMaximum);
        srcField.setMesh(sourceMesh);
        array=DataArrayDouble.New();
        ptr=sourceMesh.getNumberOfCells()*[None]
        for i in range(sourceMesh.getNumberOfCells()):
            ptr[i]=float(i+7)
            pass
        array.setValues(ptr,sourceMesh.getNumberOfCells(),1);
        srcField.setArray(array);
        srcField.setName("abc") ; srcField.setDescription("def")
        srcField.setTime(7.7,9,10)
        trgfield=remapper.transferField(srcField,4.57);
        self.assertEqual("abc",trgfield.getName())
        self.assertEqual("def",trgfield.getDescription())
        a,b,c=trgfield.getTime()
        self.assertAlmostEqual(7.7,a,14)
        self.assertEqual(b,9)
        self.assertEqual(c,10)
        values=trgfield.getArray().getValues();
        valuesExpected=[7.5 ,7. ,7.,8.,7.5];
        for i in range(targetMesh.getNumberOfCells()):
            self.assertTrue(abs(values[i]-valuesExpected[i])<1e-12);
            pass
        self.assertTrue(1==trgfield.getArray().getNumberOfComponents());
        pass

    def testPrepareEx1(self):
        sourceMesh=self.build2DSourceMesh_1();
        targetMesh=self.build2DTargetMesh_3();
        #
        remapper=MEDCouplingRemapper();
        remapper.setPrecision(1e-12);
        remapper.setIntersectionType(Triangulation);
        srcFt=MEDCouplingFieldTemplate.New(ON_CELLS);
        trgFt=MEDCouplingFieldTemplate.New(ON_CELLS);
        srcFt.setMesh(sourceMesh);
        trgFt.setMesh(targetMesh);
        self.assertEqual(1,remapper.prepareEx(srcFt,trgFt));
        srcField=MEDCouplingFieldDouble.New(ON_CELLS);
        srcField.setNature(IntensiveMaximum);
        srcField.setMesh(sourceMesh);
        array=DataArrayDouble.New();
        ptr=sourceMesh.getNumberOfCells()*[None]
        for i in range(sourceMesh.getNumberOfCells()):
            ptr[i]=float(i+7);
            pass
        array.setValues(ptr,sourceMesh.getNumberOfCells(),1);
        srcField.setArray(array);
        trgfield=remapper.transferField(srcField,4.220173);
        values=trgfield.getArray().getValues();
        valuesExpected=[7.75, 7.0625, 4.220173,8.0]
        self.assertEqual(4,trgfield.getArray().getNumberOfTuples());
        self.assertEqual(1,trgfield.getArray().getNumberOfComponents());
        for i0 in range(4):
            self.assertAlmostEqual(valuesExpected[i0],values[i0],12);
            pass
        pass

    def testPartialTransfer1(self):
        sourceMesh=self.build2DSourceMesh_1();
        targetMesh=self.build2DTargetMesh_3();
        #
        remapper=MEDCouplingRemapper();
        remapper.setPrecision(1e-12);
        remapper.setIntersectionType(Triangulation);
        srcFt=MEDCouplingFieldTemplate.New(ON_CELLS);
        trgFt=MEDCouplingFieldTemplate.New(ON_CELLS);
        srcFt.setMesh(sourceMesh);
        trgFt.setMesh(targetMesh);
        self.assertEqual(1,remapper.prepareEx(srcFt,trgFt));
        srcField=MEDCouplingFieldDouble.New(ON_CELLS);
        srcField.setNature(IntensiveMaximum);
        srcField.setMesh(sourceMesh);
        array=DataArrayDouble.New();
        ptr=sourceMesh.getNumberOfCells()*[None]
        for i in range(sourceMesh.getNumberOfCells()):
            ptr[i]=float(i+7);
            pass
        array.setValues(ptr,sourceMesh.getNumberOfCells(),1);
        srcField.setArray(array);
        trgfield=MEDCouplingFieldDouble.New(ON_CELLS);
        trgfield.setNature(IntensiveMaximum);
        trgfield.setMesh(targetMesh);
        array=DataArrayDouble.New();
        ptr=targetMesh.getNumberOfCells()*[None]
        for i in range(targetMesh.getNumberOfCells()):
            ptr[i]=4.220173;
            pass
        array.setValues(ptr,targetMesh.getNumberOfCells(),1);
        trgfield.setArray(array);
        remapper.partialTransfer(srcField,trgfield);
        values=trgfield.getArray().getValues();
        valuesExpected=[7.75, 7.0625, 4.220173,8.0]
        self.assertEqual(4,trgfield.getArray().getNumberOfTuples());
        self.assertEqual(1,trgfield.getArray().getNumberOfComponents());
        for i0 in range(4):
            self.assertAlmostEqual(valuesExpected[i0],values[i0],12);
            pass
        pass

    def testPrepareUC(self):
        # 1D
        coords=DataArrayDouble([0.,0.5,0.7])
        src=MEDCouplingUMesh("",1) ; src.setCoords(coords)
        src.allocateCells(2) ; src.insertNextCell(NORM_SEG2,[0,1]) ; src.insertNextCell(NORM_SEG2,[1,2]) ; src.finishInsertingCells()
        trg=MEDCouplingCMesh() ; arr=DataArrayDouble([-0.7,-0.1,0.2,0.7,2.,2.3])
        trg.setCoordsAt(0,arr)
        fieldSrc=MEDCouplingFieldDouble(ON_CELLS,NO_TIME) ; fieldSrc.setMesh(src) ; arrSrc=DataArrayDouble([10.,30.])
        fieldSrc.setNature(ExtensiveMaximum) ;  fieldSrc.setArray(arrSrc)
        rem=MEDCouplingRemapper()
        rem.prepare(src,trg,"P0P0")
        trgField=rem.transferField(fieldSrc,-7.)
        expected1=[-7.,4.,36.,-7.,-7.]
        self.assertEqual(5,trgField.getArray().getNumberOfTuples())
        self.assertEqual(5,len(expected1))
        for i,val in enumerate(expected1):
            self.assertAlmostEqual(expected1[i],trgField.getArray().getIJ(i,0),12);
            pass
        # 2D
        coords=DataArrayDouble([0.,0.,0.,1.,1.,1.,1.,0.,0.5,-0.2],5,2)
        src=MEDCouplingUMesh("",2) ; src.setCoords(coords)
        src.allocateCells(2) ; src.insertNextCell(NORM_TRI3,[0,1,2]) ; src.insertNextCell(NORM_TRI3,[3,4,0]) ; src.finishInsertingCells()
        trg=MEDCouplingCMesh() ; arr=DataArrayDouble([-0.7,-0.1,0.2,0.7,2.,2.3])
        trg.setCoordsAt(0,arr) ; trg.setCoordsAt(1,arr)
        fieldSrc=MEDCouplingFieldDouble(ON_CELLS,NO_TIME) ; fieldSrc.setMesh(src) ; arrSrc=DataArrayDouble([10.,30.])
        fieldSrc.setNature(ExtensiveMaximum) ;  fieldSrc.setArray(arrSrc)
        rem=MEDCouplingRemapper()
        rem.prepare(src,trg,"P0P0")
        trgField=rem.transferField(fieldSrc,-7.)
        expected2=[-7.,-7.,7.35,0.15,-7.,-7.,2.8,14.85,5.25,-7.,-7.,2.,2.5,-7.,-7.,-7.,1.2,3.,0.9,-7.,-7.,-7.,-7.,-7.,-7.]
        self.assertEqual(25,trgField.getArray().getNumberOfTuples())
        self.assertEqual(25,len(expected2))
        for i,val in enumerate(expected2):
            self.assertAlmostEqual(expected2[i],trgField.getArray().getIJ(i,0),12);
            pass
        # 3D
        coords=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,1.,0.,1.,0.,0.,0.5,-0.2,0.,0.1,0.8,1.,0.5,0.,1.],7,3)
        src=MEDCouplingUMesh("",3) ; src.setCoords(coords)
        src.allocateCells(2) ; src.insertNextCell(NORM_TETRA4,[0,1,2,5]) ; src.insertNextCell(NORM_TETRA4,[3,4,0,6]) ; src.finishInsertingCells()
        trg=MEDCouplingCMesh() ; arr=DataArrayDouble([-0.7,-0.1,0.2,0.7,2.,2.3]) ; arr2=DataArrayDouble([-0.7,0.2,0.6,1.2,2.])
        trg.setCoordsAt(0,arr) ; trg.setCoordsAt(1,arr) ; trg.setCoordsAt(2,arr2)
        src.checkConsistency(1e-10)
        trg.checkConsistencyLight()
        fieldSrc=MEDCouplingFieldDouble(ON_CELLS,NO_TIME) ; fieldSrc.setMesh(src) ; arrSrc=DataArrayDouble([10.,30.])
        fieldSrc.setNature(ExtensiveMaximum) ;  fieldSrc.setArray(arrSrc)
        rem=MEDCouplingRemapper()
        rem.prepare(src,trg,"P0P0")
        trgField=rem.transferField(fieldSrc,-7.)
        expected3=[-7.,-7.,2.925,0.015,-7.,-7.,0.9392,8.595,2.265,-7.,-7.,1.1008,1.1192,-7.,-7.,-7.,0.6392,1.6408,0.2808,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,0.81,-7.,-7.,-7.,0.1208,11.55,0.96,-7.,-7.,1.1752,0.6592,-7.,-7.,-7.,0.8512,1.7744,0.0192,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,1.92,-7.,-7.,-7.,0.12578571428571422,0.007314285714285673,-7.,-7.,-7.,0.3189253968253971,0.1879746031746033,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.,-7.]
        self.assertEqual(100,trgField.getArray().getNumberOfTuples())
        self.assertEqual(100,len(expected3))
        for i,val in enumerate(expected3):
            self.assertAlmostEqual(expected3[i],trgField.getArray().getIJ(i,0),12);
            pass
        pass

    def testPrepareCU(self):
        # 1D
        coords=DataArrayDouble([0.,0.5,0.7])
        trg=MEDCouplingUMesh("",1) ; trg.setCoords(coords)
        trg.allocateCells(2) ; trg.insertNextCell(NORM_SEG2,[0,1]) ; trg.insertNextCell(NORM_SEG2,[1,2]) ; trg.finishInsertingCells()
        src=MEDCouplingCMesh() ; arr=DataArrayDouble([-0.7,-0.1,0.2,0.7,2.,2.3])
        src.setCoordsAt(0,arr)
        fieldSrc=MEDCouplingFieldDouble(ON_CELLS,NO_TIME) ; fieldSrc.setMesh(src) ; arrTrg=DataArrayDouble([10.,30.,40.,70.,80.])
        fieldSrc.setNature(ExtensiveMaximum) ;  fieldSrc.setArray(arrTrg)
        rem=MEDCouplingRemapper()
        rem.prepare(src,trg,"P0P0")
        trgField=rem.transferField(fieldSrc,-7.)
        expected1=[44.,16.]
        self.assertEqual(2.,trgField.getArray().getNumberOfTuples())
        self.assertEqual(2,len(expected1))
        for i,val in enumerate(expected1):
            self.assertAlmostEqual(expected1[i],trgField.getArray().getIJ(i,0),12);
            pass
        # 2D
        coords=DataArrayDouble([0.,0.,0.,1.,1.,1.,1.,0.,0.5,-0.2],5,2)
        trg=MEDCouplingUMesh("",2) ; trg.setCoords(coords)
        trg.allocateCells(2) ; trg.insertNextCell(NORM_TRI3,[0,1,2]) ; trg.insertNextCell(NORM_TRI3,[3,4,0]) ; trg.finishInsertingCells()
        src=MEDCouplingCMesh() ; arr=DataArrayDouble([-0.7,-0.1,0.2,0.7,2.,2.3])
        src.setCoordsAt(0,arr) ; src.setCoordsAt(1,arr)
        fieldSrc=MEDCouplingFieldDouble(ON_CELLS,NO_TIME) ; fieldSrc.setMesh(src) ; arrSrc=DataArrayDouble([10.,30.,40.,70.,80.,110.,130.,140.,170.,180.,210.,230.,240.,270.,280.,310.,330.,340.,370.,380.,410.,430.,440.,470.,480.])
        fieldSrc.setNature(ExtensiveMaximum) ;  fieldSrc.setArray(arrSrc)
        rem=MEDCouplingRemapper()
        rem.prepare(src,trg,"P0P0")
        trgField=rem.transferField(fieldSrc,-7.)
        expected2=[441.3050624589086,68.69529914529915]
        self.assertEqual(2,trgField.getArray().getNumberOfTuples())
        self.assertEqual(2,len(expected2))
        for i,val in enumerate(expected2):
            self.assertAlmostEqual(expected2[i],trgField.getArray().getIJ(i,0),12);
            pass
        # 3D
        coords=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,1.,0.,1.,0.,0.,0.5,-0.2,0.,0.1,0.8,1.,0.5,0.,1.],7,3)
        trg=MEDCouplingUMesh("",3) ; trg.setCoords(coords)
        trg.allocateCells(2) ; trg.insertNextCell(NORM_TETRA4,[0,1,2,5]) ; trg.insertNextCell(NORM_TETRA4,[3,4,0,6]) ; trg.finishInsertingCells()
        src=MEDCouplingCMesh() ; arr=DataArrayDouble([-0.7,-0.1,0.2,0.7,2.,2.3]) ; arr2=DataArrayDouble([-0.7,0.2,0.6,1.2,2.])
        src.setCoordsAt(0,arr) ; src.setCoordsAt(1,arr) ; src.setCoordsAt(2,arr2)
        trg.checkConsistency(1e-10)
        src.checkConsistencyLight()
        fieldSrc=MEDCouplingFieldDouble(ON_CELLS,NO_TIME) ; fieldSrc.setMesh(src) ; arrSrc=DataArrayDouble(100) ; arrSrc.iota(7.7)
        fieldSrc.setNature(ExtensiveMaximum) ;  fieldSrc.setArray(arrSrc)
        rem=MEDCouplingRemapper()
        rem.prepare(src,trg,"P0P0")
        trgField=rem.transferField(fieldSrc,-7.)
        expected3=[39.635196634558845,12.13422356758468]
        self.assertEqual(2,trgField.getArray().getNumberOfTuples())
        self.assertEqual(2,len(expected3))
        for i,val in enumerate(expected3):
            self.assertAlmostEqual(expected3[i],trgField.getArray().getIJ(i,0),12);
            pass
        pass

    def testPrepareCC(self):
        # 1D
        src=MEDCouplingCMesh() ; arr=DataArrayDouble([-0.7,-0.1,0.2,0.7,2.,2.3])
        src.setCoordsAt(0,arr)
        trg=MEDCouplingCMesh() ; arr=DataArrayDouble([-0.9,-0.1,0.15])
        trg.setCoordsAt(0,arr)
        fieldSrc=MEDCouplingFieldDouble(ON_CELLS,NO_TIME) ; fieldSrc.setMesh(src) ; arrTrg=DataArrayDouble([10.,30.,40.,70.,80.])
        fieldSrc.setNature(ExtensiveMaximum) ;  fieldSrc.setArray(arrTrg)
        rem=MEDCouplingRemapper()
        rem.prepare(src,trg,"P0P0")
        trgField=rem.transferField(fieldSrc,-7.)
        expected1=[10.,25.]
        self.assertEqual(2.,trgField.getArray().getNumberOfTuples())
        self.assertEqual(2,len(expected1))
        for i,val in enumerate(expected1):
            self.assertAlmostEqual(expected1[i],trgField.getArray().getIJ(i,0),12);
            pass
        # 2D
        src=MEDCouplingCMesh() ; arr=DataArrayDouble([-0.7,-0.1,0.2,0.7,2.,2.3])
        src.setCoordsAt(0,arr) ; src.setCoordsAt(1,arr)
        trg=MEDCouplingCMesh() ; arr=DataArrayDouble([-0.9,-0.1,0.15])
        trg.setCoordsAt(0,arr) ; trg.setCoordsAt(1,arr)
        fieldSrc=MEDCouplingFieldDouble(ON_CELLS,NO_TIME) ; fieldSrc.setMesh(src) ; arrSrc=DataArrayDouble([10.,30.,40.,70.,80.,110.,130.,140.,170.,180.,210.,230.,240.,270.,280.,310.,330.,340.,370.,380.,410.,430.,440.,470.,480.])
        fieldSrc.setNature(ExtensiveMaximum) ;  fieldSrc.setArray(arrSrc)
        rem=MEDCouplingRemapper()
        rem.prepare(src,trg,"P0P0")
        trgField=rem.transferField(fieldSrc,-7.)
        expected2=[10.,25.,91.66666666666666,90.27777777777777]
        self.assertEqual(4,trgField.getArray().getNumberOfTuples())
        self.assertEqual(4,len(expected2))
        for i,val in enumerate(expected2):
            self.assertAlmostEqual(expected2[i],trgField.getArray().getIJ(i,0),12);
            pass
        # 3D
        src=MEDCouplingCMesh() ; arr=DataArrayDouble([-0.7,-0.1,0.2,0.7,2.,2.3])
        src.setCoordsAt(0,arr) ; src.setCoordsAt(1,arr) ; src.setCoordsAt(2,arr)
        trg=MEDCouplingCMesh() ; arr=DataArrayDouble([-0.9,-0.1,0.15])
        trg.setCoordsAt(0,arr) ; trg.setCoordsAt(1,arr) ; trg.setCoordsAt(2,arr)
        fieldSrc=MEDCouplingFieldDouble(ON_CELLS,NO_TIME) ; fieldSrc.setMesh(src) ; arrSrc=DataArrayDouble(125) ; arrSrc.iota(7.7)
        fieldSrc.setNature(ExtensiveMaximum) ;  fieldSrc.setArray(arrSrc) ; fieldSrc.checkConsistencyLight()
        rem=MEDCouplingRemapper()
        rem.prepare(src,trg,"P0P0")
        trgField=rem.transferField(fieldSrc,-7.)
        expected3=[7.7, 7.249999999999999, 10.583333333333332, 9.513888888888886, 27.25, 23.40277777777777, 26.180555555555546, 22.39583333333333]
        self.assertEqual(8,trgField.getArray().getNumberOfTuples())
        self.assertEqual(8,len(expected3))
        for i,val in enumerate(expected3):
            self.assertAlmostEqual(expected3[i],trgField.getArray().getIJ(i,0),12);
            pass
        pass

    # Bug when source mesh is not homogeneously oriented in source mesh
    def testNonRegressionNonHomegenousOrriented3DCells(self):
        csrc=DataArrayDouble([-0.15240000188350677,0,0,-0.1086929515004158,0,0,-0.15240000188350677,0.018142856657505035,0,-0.13054648041725159,0.0090714283287525177,0.019050000235438347,-0.13054648041725159,0.0090714283287525177,0],5,3)
        src1=MEDCouplingUMesh("src",3) ; src1.allocateCells(0) ; src1.insertNextCell(NORM_TETRA4,[0,1,4,3]) ; src1.insertNextCell(NORM_TETRA4,[2,0,4,3])
        src2=MEDCouplingUMesh("src",3) ; src2.allocateCells(0) ; src2.insertNextCell(NORM_TETRA4,[0,4,1,3]) ; src2.insertNextCell(NORM_TETRA4,[2,0,4,3])
        src1.setCoords(csrc) ; src2.setCoords(csrc)
        ctrg=DataArrayDouble([-0.15240000188350677,-0.038100000470876694,0,0.32379999756813049,-0.038100000470876694,0,-0.15240000188350677,0.076200000941753387,0,0.32379999756813049,0.076200000941753387,0,-0.15240000188350677,-0.038100000470876694,0.076200000941753387,0.32379999756813049,-0.038100000470876694,0.076200000941753387,-0.15240000188350677,0.076200000941753387,0.076200000941753387,0.32379999756813049,0.076200000941753387,0.076200000941753387],8,3)
        trg=MEDCouplingUMesh("trg",3) ; trg.allocateCells(0) ; trg.insertNextCell(NORM_HEXA8,[0,1,3,2,4,5,7,6])
        trg.setCoords(ctrg)
        rem1=MEDCouplingRemapper() ; rem1.setSplittingPolicy(PLANAR_FACE_5) ; rem1.prepare(src1,trg,"P0P0")
        rem2=MEDCouplingRemapper() ; rem2.setSplittingPolicy(PLANAR_FACE_5) ; rem2.prepare(src1,trg,"P0P0")
        mat1=rem1.getCrudeMatrix() ; mat2=rem2.getCrudeMatrix()
        self.assertEqual(1,len(mat1)) ; self.assertEqual(1,len(mat2))
        self.assertEqual(list(mat1[0].keys()),list(mat2[0].keys())) ; self.assertEqual([0,1],list(mat1[0].keys()))
        self.assertAlmostEqual(1.25884108122e-06,mat1[0][0],16) ; self.assertAlmostEqual(1.25884108122e-06,mat2[0][0],16)
        self.assertAlmostEqual(1.25884086663e-06,mat1[0][1],16) ; self.assertAlmostEqual(1.25884086663e-06,mat2[0][1],16)
        #
        d=DataArrayDouble([13.45,27.67],2,1)
        f1=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f1.setMesh(src1) ; f1.setArray(d) ; f1.setNature(IntensiveConservation)
        f2=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f2.setMesh(src2) ; f2.setArray(d) ; f2.setNature(IntensiveConservation)
        f11=rem1.transferField(f1,1e300) ; f22=rem2.transferField(f2,1e300)
        expected1=DataArrayDouble([0.012480539537637884])
        self.assertTrue(f11.getArray().isEqual(expected1,1e-15))
        self.assertTrue(f22.getArray().isEqual(expected1,1e-15))
        #
        f1.setNature(ExtensiveMaximum) ; f2.setNature(ExtensiveMaximum)
        f11=rem1.transferField(f1,1e300) ; f22=rem2.transferField(f2,1e300)
        #
        expected2=DataArrayDouble([41.12])
        self.assertTrue(f11.getArray().isEqual(expected2,1e-13))
        self.assertTrue(f22.getArray().isEqual(expected2,1e-13))
        pass

    def testCellToNodeReverse3D(self):
        c=DataArrayDouble([0.,1.,2.5])
        cc=MEDCouplingCMesh()
        cc.setCoords(c,c,c)
        um=cc.buildUnstructured()
        f=um.getMeasureField(False)
        #
        n2o=um.simplexize(PLANAR_FACE_5)
        f.setArray(f.getArray()[n2o])
        f.checkConsistencyLight()
        f.setNature(IntensiveMaximum)
        f.setTime(5.6,7,8)
        f.setName("toto") ; f.setDescription("aDescription")
        p=MEDCouplingRemapper()
        p.setIntersectionType(Barycentric)
        p.prepare(um,um,"P1P0")
        fNode=p.reverseTransferField(f,1e300)
        self.assertEqual("toto",fNode.getName())
        self.assertEqual("aDescription",fNode.getDescription())
        a,b,c=fNode.getTime()
        self.assertAlmostEqual(5.6,a,14)
        self.assertEqual(7,b) ; self.assertEqual(8,c)
        #
        integExpected=34.328125
        self.assertAlmostEqual(fNode.integral(False)[0],integExpected,14)
        self.assertAlmostEqual(f.integral(False)[0],integExpected,14)
        pass

    def testGauss2Gauss2DValidated(self):
        srcFt=MEDCouplingDataForTest.buildFieldOnGauss_1()
        trgFt=MEDCouplingDataForTest.buildFieldOnGauss_2()
        src=MEDCouplingFieldDouble(srcFt)
        self.assertEqual(srcFt.getMesh().getHiddenCppPointer(),src.getMesh().getHiddenCppPointer())
        self.assertEqual(srcFt.getDiscretization().getHiddenCppPointer(),src.getDiscretization().getHiddenCppPointer())
        #values given by ASTER usecase
        src.setArray(DataArrayDouble([1.,1.,0.,0.,1.,1.,1.,1.,0.,0.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.]))
        src.getArray().setInfoOnComponents(["DOMA"])
        rem=MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepareEx(srcFt,trgFt)
        trg=rem.transferField(src,1e300)
        self.assertEqual(trg.getMesh().getHiddenCppPointer(),trgFt.getMesh().getHiddenCppPointer())
        self.assertEqual(trg.getDiscretization().getHiddenCppPointer(),trgFt.getDiscretization().getHiddenCppPointer())
        #values given after interpolation in ASTER
        arrExpected=DataArrayDouble([1.,1.,1.,0.,0.,0.,1.,1.,1.,1.,1.,0.,0.,0.,1.,1.,1.,1.,1.,0.,0.,0.,1.,1.,1.,1.,1.,0.,0.,0.,1.,1.,1.,1.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,1.,1.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,1.,1.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,1.,1.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,1.,1.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,1.,1.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,0.,0.,0.,0.,1.,1.,1.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,0.,0.,0.,0.,1.,1.,1.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,0.,0.,0.,0.,1.,1.,1.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,0.,0.,0.,0.,1.,1.,1.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,0.,0.,0.,0.,1.,1.,1.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.]) ; arrExpected.setInfoOnComponents(["DOMA"])
        self.assertTrue(trg.getArray().isEqual(arrExpected,1e-12))
        #
        # second part of the test : reverse source and target
        #
        rem.prepareEx(trgFt,srcFt)# sorry trgFt is in the place of source and srcFt in the place of target it is not a bug
        trg=MEDCouplingFieldDouble(trgFt)
        #values given after interpolation in ASTER
        trg.setArray(DataArrayDouble([1.,1.,0.,0.,1.,0.,1.,1.,1.,1.,0.,0.,1.,0.,1.,1.,1.,1.,0.,0.,1.,0.,1.,1.,1.,1.,0.,0.,1.,0.,1.,1.,1.,1.,0.,0.,1.,0.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,1.,1.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,1.,1.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,1.,1.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,1.,1.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,1.,1.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,0.,0.,0.,0.,1.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.]))
        trg.getArray().setInfoOnComponents(["DOMA"])
        src=rem.transferField(trg,1e300)
        #values given after interpolation in ASTER
        arrExpected2=DataArrayDouble([1.,1.,0.,0.,1.,1.,1.,1.,0.,0.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,1.,1.,1., 1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,1.,0.,1.,1.,1.,1.,0.,0.,0.,0.,1.,0.,0.,0.,1.,1.,1.,0.,1.,1.,1.,1.,1.,1.,0.,0.,1.,1.,0.,1.,1.,1.,0.,1.,0.,0.,0.,1.,0.,0.,0.,1.,1.,1.,0.,1.,1.,1.,1.,1.]) ; arrExpected2.setInfoOnComponents(["DOMA"])
        # modification of values in ASTER due to modification of algorithm
        # target PG 82 in target cell 32(C)/36 PG 1(C)/9 is in source cell 58(C)/120 source Gauss point 113 (1(C)/4). Values must be 1. and not 0.
        arrExpected2.setIJ(82,0,1.)
        self.assertTrue(src.getArray().isEqual(arrExpected2,1e-12))
        pass

    def testGauss2Gauss3DValidated(self):
        srcFt=MEDCouplingDataForTest.buildFieldOnGauss_3()
        trgFt=MEDCouplingDataForTest.buildFieldOnGauss_4()
        src=MEDCouplingFieldDouble(srcFt)
        self.assertEqual(srcFt.getMesh().getHiddenCppPointer(),src.getMesh().getHiddenCppPointer())
        self.assertEqual(srcFt.getDiscretization().getHiddenCppPointer(),src.getDiscretization().getHiddenCppPointer())
        #values given by ASTER usecase
        src.setArray(DataArrayDouble([0.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,1.,1.,0.,0.,1.,1.,1.,1.,0.,0.,1.,1.,0.,0.]))
        src.getArray().setInfoOnComponents(["DOMA"])
        rem=MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepareEx(srcFt,trgFt)
        trg=rem.transferField(src,1e300)
        self.assertEqual(trg.getMesh().getHiddenCppPointer(),trgFt.getMesh().getHiddenCppPointer())
        self.assertEqual(trg.getDiscretization().getHiddenCppPointer(),trgFt.getDiscretization().getHiddenCppPointer())
        #values given after interpolation in ASTER
        arrExpected=DataArrayDouble([0.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,0.,0.,0.,1.,1.,0.,0.,1.,1.,1.,1.,0.,1.,1.,1.,0.,1.]) ; arrExpected.setInfoOnComponents(["DOMA"])
        self.assertTrue(trg.getArray().isEqual(arrExpected,1e-12))
        #
        # second part of the test : reverse source and target
        #
        rem.prepareEx(trgFt,srcFt)# sorry trgFt is in the place of source and srcFt in the place of target it is not a bug
        trg=MEDCouplingFieldDouble(trgFt)
        #values given after interpolation in ASTER
        trg.setArray(DataArrayDouble([0.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,1.,1.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]))
        trg.getArray().setInfoOnComponents(["DOMA"])
        src=rem.transferField(trg,1e300)
        #values given after interpolation in ASTER
        arrExpected2=DataArrayDouble([0.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,0.,1.,0.,1.,1.,1.,0.,1.,1.,1.,1.,0.,1.,1.,1.,0.,1.]) ; arrExpected2.setInfoOnComponents(["DOMA"])
        self.assertTrue(src.getArray().isEqual(arrExpected2,1e-12))
        pass

    def testSwig2MixOfUMesh(self):
        arr0=DataArrayDouble([0,1,1.5]) ; arr1=DataArrayDouble([0,1])
        sc=MEDCouplingCMesh() ; sc.setCoords(arr0,arr1,arr1)
        tc=sc.deepCopy() ; tc.translate([0.4,0.3,0.3])
        # umesh-umesh
        # 90 (umesh-1sgtumesh)
        rem=MEDCouplingRemapper() ; rem.setIntersectionType(Triangulation)
        s=sc.buildUnstructured() ; t=tc.build1SGTUnstructured()
        self.assertTrue(isinstance(s,MEDCouplingUMesh))
        self.assertTrue(isinstance(t,MEDCoupling1SGTUMesh))
        rem.prepare(s,t,"P0P0")
        mat=rem.getCrudeMatrix()
        self.assertEqual(2,len(mat)) ; self.assertEqual(2,len(mat[0])) ; self.assertEqual(1,len(mat[1]))
        self.assertAlmostEqual(0.294,mat[0][0],14) ; self.assertAlmostEqual(0.196,mat[0][1],14) ; self.assertAlmostEqual(0.049,mat[1][1],14)
        del s,t
        # 91 (umesh-1dgtumesh)
        rem=MEDCouplingRemapper() ; rem.setIntersectionType(Triangulation)
        s=sc.buildUnstructured() ; t=tc.buildUnstructured() ; t.convertAllToPoly() ; t=MEDCoupling1DGTUMesh(t)
        self.assertTrue(isinstance(s,MEDCouplingUMesh))
        self.assertTrue(isinstance(t,MEDCoupling1DGTUMesh))
        rem.prepare(s,t,"P0P0")
        mat=rem.getCrudeMatrix()
        self.assertEqual(2,len(mat)) ; self.assertEqual(2,len(mat[0])) ; self.assertEqual(1,len(mat[1]))
        self.assertAlmostEqual(0.294,mat[0][0],14) ; self.assertAlmostEqual(0.196,mat[0][1],14) ; self.assertAlmostEqual(0.049,mat[1][1],14)
        del s,t
        # 165 (1sgtumesh-umesh)
        rem=MEDCouplingRemapper() ; rem.setIntersectionType(Triangulation)
        s=sc.build1SGTUnstructured() ; t=tc.buildUnstructured()
        self.assertTrue(isinstance(s,MEDCoupling1SGTUMesh))
        self.assertTrue(isinstance(t,MEDCouplingUMesh))
        rem.prepare(s,t,"P0P0")
        mat=rem.getCrudeMatrix()
        self.assertEqual(2,len(mat)) ; self.assertEqual(2,len(mat[0])) ; self.assertEqual(1,len(mat[1]))
        self.assertAlmostEqual(0.294,mat[0][0],14) ; self.assertAlmostEqual(0.196,mat[0][1],14) ; self.assertAlmostEqual(0.049,mat[1][1],14)
        del s,t
        # 181 (1dgtumesh-umesh
        rem=MEDCouplingRemapper() ; rem.setIntersectionType(Triangulation)
        s=sc.buildUnstructured() ; s.convertAllToPoly() ; s=MEDCoupling1DGTUMesh(s) ; t=tc.buildUnstructured()
        self.assertTrue(isinstance(s,MEDCoupling1DGTUMesh))
        self.assertTrue(isinstance(t,MEDCouplingUMesh))
        rem.prepare(s,t,"P0P0")
        mat=rem.getCrudeMatrix()
        self.assertEqual(2,len(mat)) ; self.assertEqual(2,len(mat[0])) ; self.assertEqual(1,len(mat[1]))
        self.assertAlmostEqual(0.294,mat[0][0],14) ; self.assertAlmostEqual(0.196,mat[0][1],14) ; self.assertAlmostEqual(0.049,mat[1][1],14)
        del s,t
        # 170 (1sgtumesh-1sgtumesh)
        rem=MEDCouplingRemapper() ; rem.setIntersectionType(Triangulation)
        s=sc.build1SGTUnstructured() ; t=tc.build1SGTUnstructured()
        self.assertTrue(isinstance(s,MEDCoupling1SGTUMesh))
        self.assertTrue(isinstance(t,MEDCoupling1SGTUMesh))
        rem.prepare(s,t,"P0P0")
        mat=rem.getCrudeMatrix()
        self.assertEqual(2,len(mat)) ; self.assertEqual(2,len(mat[0])) ; self.assertEqual(1,len(mat[1]))
        self.assertAlmostEqual(0.294,mat[0][0],14) ; self.assertAlmostEqual(0.196,mat[0][1],14) ; self.assertAlmostEqual(0.049,mat[1][1],14)
        del s,t
        # 171 (1sgtumesh-1dgtumesh)
        rem=MEDCouplingRemapper() ; rem.setIntersectionType(Triangulation)
        s=sc.build1SGTUnstructured() ; t=tc.buildUnstructured() ; t.convertAllToPoly() ; t=MEDCoupling1DGTUMesh(t)
        self.assertTrue(isinstance(s,MEDCoupling1SGTUMesh))
        self.assertTrue(isinstance(t,MEDCoupling1DGTUMesh))
        rem.prepare(s,t,"P0P0")
        mat=rem.getCrudeMatrix()
        self.assertEqual(2,len(mat)) ; self.assertEqual(2,len(mat[0])) ; self.assertEqual(1,len(mat[1]))
        self.assertAlmostEqual(0.294,mat[0][0],14) ; self.assertAlmostEqual(0.196,mat[0][1],14) ; self.assertAlmostEqual(0.049,mat[1][1],14)
        del s,t
        # 186 (1dgtumesh-1sgtumesh)
        rem=MEDCouplingRemapper() ; rem.setIntersectionType(Triangulation)
        s=sc.buildUnstructured() ; s.convertAllToPoly() ; s=MEDCoupling1DGTUMesh(s) ; t=tc.build1SGTUnstructured()
        self.assertTrue(isinstance(s,MEDCoupling1DGTUMesh))
        self.assertTrue(isinstance(t,MEDCoupling1SGTUMesh))
        rem.prepare(s,t,"P0P0")
        mat=rem.getCrudeMatrix()
        self.assertEqual(2,len(mat)) ; self.assertEqual(2,len(mat[0])) ; self.assertEqual(1,len(mat[1]))
        self.assertAlmostEqual(0.294,mat[0][0],14) ; self.assertAlmostEqual(0.196,mat[0][1],14) ; self.assertAlmostEqual(0.049,mat[1][1],14)
        del s,t
        # 187 (1dgtumesh-1dgtumesh)
        rem=MEDCouplingRemapper() ; rem.setIntersectionType(Triangulation)
        s=sc.buildUnstructured() ; s.convertAllToPoly() ; s=MEDCoupling1DGTUMesh(s) ; t=tc.buildUnstructured() ; t.convertAllToPoly() ; t=MEDCoupling1DGTUMesh(t)
        self.assertTrue(isinstance(s,MEDCoupling1DGTUMesh))
        self.assertTrue(isinstance(t,MEDCoupling1DGTUMesh))
        rem.prepare(s,t,"P0P0")
        mat=rem.getCrudeMatrix()
        self.assertEqual(2,len(mat)) ; self.assertEqual(2,len(mat[0])) ; self.assertEqual(1,len(mat[1]))
        self.assertAlmostEqual(0.294,mat[0][0],14) ; self.assertAlmostEqual(0.196,mat[0][1],14) ; self.assertAlmostEqual(0.049,mat[1][1],14)
        del s,t
        # (umesh-cmesh)
        # 167 (1sgtumesh-cmesh)
        rem=MEDCouplingRemapper() ; rem.setIntersectionType(Triangulation)
        s=sc.build1SGTUnstructured() ; t=tc.deepCopy()
        self.assertTrue(isinstance(s,MEDCoupling1SGTUMesh))
        self.assertTrue(isinstance(t,MEDCouplingCMesh))
        rem.prepare(s,t,"P0P0")
        mat=rem.getCrudeMatrix()
        self.assertEqual(2,len(mat)) ; self.assertEqual(2,len(mat[0])) ; self.assertEqual(1,len(mat[1]))
        self.assertAlmostEqual(0.294,mat[0][0],14) ; self.assertAlmostEqual(0.196,mat[0][1],14) ; self.assertAlmostEqual(0.049,mat[1][1],14)
        del s,t
        # 183 (1dgtumesh-cmesh)
        #rem=MEDCouplingRemapper() ; rem.setIntersectionType(Triangulation)
        #s=sc.buildUnstructured() ; s.convertAllToPoly() ; s=MEDCoupling1DGTUMesh(s) ; t=tc.deepCopy()
        #self.assertTrue(isinstance(s,MEDCoupling1DGTUMesh))
        #self.assertTrue(isinstance(t,MEDCouplingCMesh))
        #rem.prepare(s,t,"P0P0")
        #mat=rem.getCrudeMatrix()
        #self.assertEqual(2,len(mat)) ; self.assertEqual(2,len(mat[0])) ; self.assertEqual(1,len(mat[1]))
        #self.assertAlmostEqual(0.294,mat[0][0],14) ; self.assertAlmostEqual(0.196,mat[0][1],14) ; self.assertAlmostEqual(0.049,mat[1][1],14)
        #del s,t
        # (cmesh-umesh)
        # 122 (cmesh-1sgtumesh)
        rem=MEDCouplingRemapper() ; rem.setIntersectionType(Triangulation)
        s=sc.deepCopy() ; t=tc.build1SGTUnstructured()
        self.assertTrue(isinstance(s,MEDCouplingCMesh))
        self.assertTrue(isinstance(t,MEDCoupling1SGTUMesh))
        rem.prepare(s,t,"P0P0")
        mat=rem.getCrudeMatrix()
        self.assertEqual(2,len(mat)) ; self.assertEqual(2,len(mat[0])) ; self.assertEqual(1,len(mat[1]))
        self.assertAlmostEqual(0.294,mat[0][0],14) ; self.assertAlmostEqual(0.196,mat[0][1],14) ; self.assertAlmostEqual(0.049,mat[1][1],14)
        del s,t
        # 123 (cmesh-1dgtumesh)
        #rem=MEDCouplingRemapper() ; rem.setIntersectionType(Triangulation)
        #s=sc.deepCopy() ; t=tc.buildUnstructured() ; t.convertAllToPoly() ; t=MEDCoupling1DGTUMesh(t)
        #self.assertTrue(isinstance(s,MEDCouplingCMesh))
        #self.assertTrue(isinstance(t,MEDCoupling1DGTUMesh))
        #rem.prepare(s,t,"P0P0")
        #mat=rem.getCrudeMatrix()
        #self.assertEqual(2,len(mat)) ; self.assertEqual(2,len(mat[0])) ; self.assertEqual(1,len(mat[1]))
        #self.assertAlmostEqual(0.294,mat[0][0],14) ; self.assertAlmostEqual(0.196,mat[0][1],14) ; self.assertAlmostEqual(0.049,mat[1][1],14)
        #del s,t
        pass

    def testSwig2BarycentricP1P13D_1(self):
        sCoo=DataArrayDouble([0.313,0.00218,6.90489,0.313,0.10692667,6.90489,0.313,0.10692667,6.96790167,0.313,0.00218,6.9773125,0.313,0.21167333,6.90489,0.313,0.21167333,6.95849083,0.313,0.31642,6.90489,0.313,0.31642,6.94908,0.313,0.09383333,7.04891667,0.313,0.00218,7.049735,0.313,0.18548667,7.04809833,0.313,0.27714,7.04728,0.313,0.05782667,7.133205,0.313,0.00218,7.1221575,0.313,0.11347333,7.1442525,0.313,0.16912,7.1553,0.313,0.02509333,7.19458,0.313,0.00218,7.19458,0.313,0.04800667,7.19458,0.313,0.07092,7.19458,0.31005609,0.00218,6.90460005,0.31005609,0.10692667,6.90460005,0.29776312,0.10692667,6.96640097,0.29592716,0.00218,6.97563097,0.31005609,0.21167333,6.90460005,0.29959908,0.21167333,6.95717096,0.31005609,0.31642,6.90460005,0.30143505,0.31642,6.94794095,0.28195788,0.09383333,7.04585928,0.28179823,0.00218,7.04666189,0.28211753,0.18548667,7.04505668,0.28227718,0.27714,7.04425407,0.26551404,0.05782667,7.12852804,0.2676693,0.00218,7.11769282,0.26335878,0.11347333,7.13936327,0.26120352,0.16912,7.15019849,0.25354037,0.02509333,7.18872374,0.25354037,0.00218,7.18872374,0.25354037,0.04800667,7.18872374,0.25354037,0.07092,7.18872374,0.30722531,0.00218,6.90374134,0.30722531,0.10692667,6.90374134,0.28311179,0.10692667,6.96195653,0.27951042,0.00218,6.97065101,0.30722531,0.21167333,6.90374134,0.28671316,0.21167333,6.95326205,0.30722531,0.31642,6.90374134,0.29031453,0.31642,6.94456758,0.25210869,0.09383333,7.03680463,0.25179553,0.00218,7.03756067,0.25242185,0.18548667,7.03604859,0.25273501,0.27714,7.03529255,0.21985294,0.05782667,7.1146769,0.22408063,0.00218,7.10447034,0.21562524,0.11347333,7.12488346,0.21139755,0.16912,7.13509002,0.19636574,0.02509333,7.17138,0.19636574,0.00218,7.17138,0.19636574,0.04800667,7.17138,0.19636574,0.07092,7.17138,0.30461645,0.00218,6.90234688,0.30461645,0.10692667,6.90234688,0.26960904,0.10692667,6.95473916,0.26438066,0.00218,6.96256398,0.30461645,0.21167333,6.90234688,0.27483742,0.21167333,6.94691434,0.30461645,0.31642,6.90234688,0.2800658,0.31642,6.93908952,0.22459952,0.09383333,7.02210067,0.22414487,0.00218,7.02278109,0.22505416,0.18548667,7.02142025,0.2255088,0.27714,7.02073983,0.17777143,0.05782667,7.09218386,0.18390909,0.00218,7.0829982,0.17163377,0.11347333,7.10136952,0.1654961,0.16912,7.11055518,0.1436733,0.02509333,7.14321531,0.1436733,0.00218,7.14321531,0.1436733,0.04800667,7.14321531,0.1436733,0.07092,7.14321531,0.30232976,0.00218,6.90047024,0.30232976,0.10692667,6.90047024,0.25777378,0.10692667,6.94502622,0.25111932,0.00218,6.95168068,0.30232976,0.21167333,6.90047024,0.26442825,0.21167333,6.93837175,0.30232976,0.31642,6.90047024,0.27108271,0.31642,6.93171729,0.20048753,0.09383333,7.00231247,0.19990888,0.00218,7.00289112,0.20106618,0.18548667,7.00173382,0.20164482,0.27714,7.00115518,0.14088667,0.05782667,7.06191333,0.14869844,0.00218,7.05410156,0.13307491,0.11347333,7.06972509,0.12526315,0.16912,7.07753685,0.097488,0.02509333,7.105312,0.097488,0.00218,7.105312,0.097488,0.04800667,7.105312,0.097488,0.07092,7.105312,0.30045312,0.00218,6.89818355,0.30045312,0.10692667,6.89818355,0.24806084,0.10692667,6.93319096,0.24023602,0.00218,6.93841934,0.30045312,0.21167333,6.89818355,0.25588566,0.21167333,6.92796258,0.30045312,0.31642,6.89818355,0.26371048,0.31642,6.9227342,0.18069933,0.09383333,6.97820048,0.18001891,0.00218,6.97865513,0.18137975,0.18548667,6.97774584,0.18206017,0.27714,6.9772912,0.11061614,0.05782667,7.02502857,0.1198018,0.00218,7.01889091,0.10143048,0.11347333,7.03116623,0.09224482,0.16912,7.0373039,0.05958469,0.02509333,7.0591267,0.05958469,0.00218,7.0591267,0.05958469,0.04800667,7.0591267,0.05958469,0.07092,7.0591267,0.29905866,0.00218,6.89557469,0.29905866,0.10692667,6.89557469,0.24084347,0.10692667,6.91968821,0.23214899,0.00218,6.92328958,0.29905866,0.21167333,6.89557469,0.24953795,0.21167333,6.91608684,0.29905866,0.31642,6.89557469,0.25823242,0.31642,6.91248547,0.16599537,0.09383333,6.95069131,0.16523933,0.00218,6.95100447,0.16675141,0.18548667,6.95037815,0.16750745,0.27714,6.95006499,0.0881231,0.05782667,6.98294706,0.09832966,0.00218,6.97871937,0.07791654,0.11347333,6.98717476,0.06770998,0.16912,6.99140245,0.03142,0.02509333,7.00643426,0.03142,0.00218,7.00643426,0.03142,0.04800667,7.00643426,0.03142,0.07092,7.00643426,0.29819995,0.00218,6.89274391,0.29819995,0.10692667,6.89274391,0.23639903,0.10692667,6.90503688,0.22716903,0.00218,6.90687284,0.29819995,0.21167333,6.89274391,0.24562904,0.21167333,6.90320092,0.29819995,0.31642,6.89274391,0.25485905,0.31642,6.90136495,0.15694072,0.09383333,6.92084212,0.15613811,0.00218,6.92100177,0.15774332,0.18548667,6.92068247,0.15854593,0.27714,6.92052282,0.07427196,0.05782667,6.93728596,0.08510718,0.00218,6.9351307,0.06343673,0.11347333,6.93944122,0.05260151,0.16912,6.94159648,0.01407626,0.02509333,6.94925963,0.01407626,0.00218,6.94925963,0.01407626,0.04800667,6.94925963,0.01407626,0.07092,6.94925963,0.29792818,0.00218,6.89054043,0.29792818,0.10692667,6.89054043,0.23499241,0.10692667,6.89363227,0.22559291,0.00218,6.89409403,0.29792818,0.21167333,6.89054043,0.24439191,0.21167333,6.8931705,0.29792818,0.31642,6.89054043,0.25379141,0.31642,6.89270873,0.154075,0.09383333,6.89760748,0.15325765,0.00218,6.89764764,0.15489234,0.18548667,6.89756733,0.15570969,0.27714,6.89752718,0.06988819,0.05782667,6.90174332,0.08092238,0.00218,6.90120124,0.058854,0.11347333,6.90228539,0.04781981,0.16912,6.90282747,0.00858712,0.02509333,6.90475485,0.00858712,0.00218,6.90475485,0.00858712,0.04800667,6.90475485,0.00858712,0.07092,6.90475485,0.29791,0.00218,6.820902,0.29791,0.10692667,6.820902,0.23489833,0.10692667,6.820902,0.2254875,0.00218,6.820902,0.29791,0.21167333,6.820902,0.24430917,0.21167333,6.820902,0.29791,0.31642,6.820902,0.25372,0.31642,6.820902,0.15388333,0.09383333,6.820902,0.153065,0.00218,6.820902,0.15470167,0.18548667,6.820902,0.15552,0.27714,6.820902,0.069595,0.05782667,6.820902,0.0806425,0.00218,6.820902,0.0585475,0.11347333,6.820902,0.0475,0.16912,6.820902,0.00822,0.02509333,6.820902,0.00822,0.00218,6.820902,0.00822,0.04800667,6.820902,0.00822,0.07092,6.820902],200,3)
        sConn=DataArrayInt([0,1,2,3,20,21,22,23,1,4,5,2,21,24,25,22,4,6,7,5,24,26,27,25,3,2,8,9,23,22,28,29,2,5,10,8,22,25,30,28,5,7,11,10,25,27,31,30,9,8,12,13,29,28,32,33,8,10,14,12,28,30,34,32,10,11,15,14,30,31,35,34,13,12,16,17,33,32,36,37,12,14,18,16,32,34,38,36,14,15,19,18,34,35,39,38,20,21,22,23,40,41,42,43,21,24,25,22,41,44,45,42,24,26,27,25,44,46,47,45,23,22,28,29,43,42,48,49,22,25,30,28,42,45,50,48,25,27,31,30,45,47,51,50,29,28,32,33,49,48,52,53,28,30,34,32,48,50,54,52,30,31,35,34,50,51,55,54,33,32,36,37,53,52,56,57,32,34,38,36,52,54,58,56,34,35,39,38,54,55,59,58,40,41,42,43,60,61,62,63,41,44,45,42,61,64,65,62,44,46,47,45,64,66,67,65,43,42,48,49,63,62,68,69,42,45,50,48,62,65,70,68,45,47,51,50,65,67,71,70,49,48,52,53,69,68,72,73,48,50,54,52,68,70,74,72,50,51,55,54,70,71,75,74,53,52,56,57,73,72,76,77,52,54,58,56,72,74,78,76,54,55,59,58,74,75,79,78,60,61,62,63,80,81,82,83,61,64,65,62,81,84,85,82,64,66,67,65,84,86,87,85,63,62,68,69,83,82,88,89,62,65,70,68,82,85,90,88,65,67,71,70,85,87,91,90,69,68,72,73,89,88,92,93,68,70,74,72,88,90,94,92,70,71,75,74,90,91,95,94,73,72,76,77,93,92,96,97,72,74,78,76,92,94,98,96,74,75,79,78,94,95,99,98,80,81,82,83,100,101,102,103,81,84,85,82,101,104,105,102,84,86,87,85,104,106,107,105,83,82,88,89,103,102,108,109,82,85,90,88,102,105,110,108,85,87,91,90,105,107,111,110,89,88,92,93,109,108,112,113,88,90,94,92,108,110,114,112,90,91,95,94,110,111,115,114,93,92,96,97,113,112,116,117,92,94,98,96,112,114,118,116,94,95,99,98,114,115,119,118,100,101,102,103,120,121,122,123,101,104,105,102,121,124,125,122,104,106,107,105,124,126,127,125,103,102,108,109,123,122,128,129,102,105,110,108,122,125,130,128,105,107,111,110,125,127,131,130,109,108,112,113,129,128,132,133,108,110,114,112,128,130,134,132,110,111,115,114,130,131,135,134,113,112,116,117,133,132,136,137,112,114,118,116,132,134,138,136,114,115,119,118,134,135,139,138,120,121,122,123,140,141,142,143,121,124,125,122,141,144,145,142,124,126,127,125,144,146,147,145,123,122,128,129,143,142,148,149,122,125,130,128,142,145,150,148,125,127,131,130,145,147,151,150,129,128,132,133,149,148,152,153,128,130,134,132,148,150,154,152,130,131,135,134,150,151,155,154,133,132,136,137,153,152,156,157,132,134,138,136,152,154,158,156,134,135,139,138,154,155,159,158,140,141,142,143,160,161,162,163,141,144,145,142,161,164,165,162,144,146,147,145,164,166,167,165,143,142,148,149,163,162,168,169,142,145,150,148,162,165,170,168,145,147,151,150,165,167,171,170,149,148,152,153,169,168,172,173,148,150,154,152,168,170,174,172,150,151,155,154,170,171,175,174,153,152,156,157,173,172,176,177,152,154,158,156,172,174,178,176,154,155,159,158,174,175,179,178,160,161,162,163,180,181,182,183,161,164,165,162,181,184,185,182,164,166,167,165,184,186,187,185,163,162,168,169,183,182,188,189,162,165,170,168,182,185,190,188,165,167,171,170,185,187,191,190,169,168,172,173,189,188,192,193,168,170,174,172,188,190,194,192,170,171,175,174,190,191,195,194,173,172,176,177,193,192,196,197,172,174,178,176,192,194,198,196,174,175,179,178,194,195,199,198])
        s=MEDCoupling1SGTUMesh("source",NORM_HEXA8) ; s.setCoords(sCoo)
        s.setNodalConnectivity(sConn)
        #
        tCoo=DataArrayDouble([0.328,0.012,6.8598,0.328,0.168320184237353,6.8598,0.328,0.324640368474706,6.8598,0.328,0.0,6.8598,0.298,0.012,6.8598,0.1565,0.012,6.8598,0.180205346493166,0.144794653506834,6.8598,0.298,0.168320184237353,6.8598,0.0,0.012,6.8598,0.0916755774886107,0.233324422511389,6.8598,0.298,0.324640368474706,6.8598,0.298,0.0,6.8598,0.1565,0.0,6.8598,0.0,0.0,6.8598,0.328,0.012,7.2298,0.328,0.168320184237353,7.2298,0.328,0.324640368474706,7.2298,0.328,0.0,7.2298,0.298,0.012,7.2298,0.1565,0.012,7.2298,0.180205346493166,0.144794653506834,7.2298,0.298,0.168320184237353,7.2298,0.0,0.012,7.2298,0.0916755774886107,0.233324422511389,7.2298,0.298,0.324640368474706,7.2298,0.298,0.0,7.2298,0.1565,0.0,7.2298,0.0,0.0,7.2298],28,3)
        tConn=DataArrayInt([4,5,6,7,18,19,20,21,5,8,9,6,19,22,23,20,6,9,10,7,20,23,24,21,11,12,5,4,25,26,19,18,12,13,8,5,26,27,22,19,3,11,4,0,17,25,18,14,0,4,7,1,14,18,21,15,1,7,10,2,15,21,24,16])
        t=MEDCoupling1SGTUMesh("target",NORM_HEXA8) ; t.setCoords(tCoo)
        t.setNodalConnectivity(tConn)
        #
        s.simplexize(PLANAR_FACE_5)
        aRemapper=MEDCouplingRemapper()
        aRemapper.setPrecision(1e-12)
        aRemapper.setIntersectionType(Barycentric)
        self.assertEqual(aRemapper.prepare(s,t,'P1P1'),1)
        m=aRemapper.getCrudeMatrix()
        self.assertEqual(len(m),28)
        for i in range(28):
            if i not in [5,6]:
                self.assertEqual(len(m[i]),0)
                pass
            pass
        self.assertEqual(len(m[5]),4)
        self.assertEqual(len(m[6]),4)
        self.assertAlmostEqual(0.10714286103952797,m[5][168],12)
        self.assertAlmostEqual(0.35691534416938014,m[5][169],12)
        self.assertAlmostEqual(0.04492099619713096,m[5][163],12)
        self.assertAlmostEqual(0.49102079859396097,m[5][189],12)
        self.assertAlmostEqual(0.14039089397104254,m[6][185],12)
        self.assertAlmostEqual(0.16362822318261033,m[6][162],12)
        self.assertAlmostEqual(0.3438363717836785 ,m[6][188],12)
        self.assertAlmostEqual(0.3521445110626687 ,m[6][170],12)
        pass

    def testSwig2MappedBarycentricP1P12D_1(self):
        """ Testing mapped barycentric P1P1 projection
        (uses analytical mapping from square to arbitrary convex quadrangle)
        """
        n = 5
        sCoo = DataArrayDouble(n,1)
        sCoo.iota(0.0);     sCoo /= float(n-1)
        m = MEDCouplingCMesh("target")
        m.setCoordsAt(0, sCoo)
        m.setCoordsAt(1, sCoo)
        tgt = m.buildUnstructured()
        coo = tgt.getCoords()
        orig = coo.deepCopy();   orig[:,0] = 10.0; orig[:,1] = 15.0
        pt_a = coo.deepCopy();   pt_a[:,0] = -0.3; pt_a[:,1] = 1.0
        pt_b = coo.deepCopy();   pt_b[:,0] = 2.0;  pt_b[:,1] = 3.0
        pt_c = coo.deepCopy();   pt_c[:,0] = 1.0;  pt_c[:,1] = 0.0
        # P = x*C+y*A + xy(B-A-C) + ORIGIN
        coo2 = coo[:,0]*pt_c + coo[:, 1]*pt_a + coo[:, 0]*coo[:, 1]*(pt_b - pt_a - pt_c) + orig

        tgt.setCoords(coo2)

        sCoo = DataArrayDouble([0.0,0.0,  -0.3,1.0,  2.0,3.0,  1.0,0.0],4,2)
        sCoo[:,0] += 10.0;  sCoo[:,1] += 15.0;
        sConn = DataArrayInt([0,1,2,3])
        s = MEDCoupling1SGTUMesh("source",NORM_QUAD4) ; s.setCoords(sCoo)
        s.setNodalConnectivity(sConn)
        #
        aRemapper=MEDCouplingRemapper()
        aRemapper.setPrecision(1e-12)
        aRemapper.setIntersectionType(MappedBarycentric)
        self.assertEqual(aRemapper.prepare(s,tgt,'P1P1'),1)
        srcField = MEDCouplingFieldDouble(ON_NODES, ONE_TIME)
        srcField.setNature(IntensiveMaximum)
        srcField.setMesh(s); srcField.setName("field")
        srcField.setArray(DataArrayDouble([1.0,2.0,3.0,4.0]))
        tgtF = aRemapper.transferField(srcField, 1e+300)
        ref = [1.0, 1.75, 2.5, 3.25, 4.0, 1.25, 1.875, 2.5, 3.125, 3.75, 1.5, 2.0, 2.5, 3.0, 3.5, 1.75,
         2.125, 2.5, 2.875, 3.25, 2.0, 2.25, 2.5, 2.75, 3.0]
        val = tgtF.getArray().getValues()
        for i, ref_v in enumerate(ref):
            self.assertAlmostEqual(ref_v, val[i])
        pass

    def testSwig2MappedBarycentricP1P13_1(self):
        """ Testing mapped barycentric P1P1 projection in 3D (uses orthogonal distances to 
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
        sCoo = DataArrayDouble(n,1)
        sCoo.iota(0.0)
        sCoo /= float(n-1)
        m = MEDCouplingCMesh("target")
        m.setCoordsAt(0, sCoo)
        m.setCoordsAt(1, sCoo)
        m.setCoordsAt(2, sCoo)
        tgt = m.buildUnstructured()
        coo = tgt.getCoords()
        pt_0 = coo.deepCopy(); pt_0[:,0] = -0.3; pt_0[:,1] = 1.0; pt_0[:,2] = 1.0
        pt_1 = coo.deepCopy(); pt_1[:,0] = 0.0; pt_1[:,1] = 0.0; pt_1[:,2] = 1.0
        pt_2 = coo.deepCopy(); pt_2[:,0] = 1.0; pt_2[:,1] = 0.0; pt_2[:,2] = 1.0
        pt_3 = coo.deepCopy(); pt_3[:,0] = 2.0; pt_3[:,1] = 3.0; pt_3[:,2] = 1.0

        pt_4 = coo.deepCopy(); pt_4[:,0] = -0.3; pt_4[:,1] = 1.0; pt_4[:,2] = 0.0
        orig = coo.deepCopy(); orig[:,0] = 10.0; orig[:,1] = 15.0; orig[:,2] = 20.0
        pt_6 = coo.deepCopy(); pt_6[:,0] = 1.0; pt_6[:,1] = 0.0; pt_6[:,2] = 0.0
        pt_7 = coo.deepCopy(); pt_7[:,0] = 2.0; pt_7[:,1] = 3.0; pt_7[:,2] = 0.0
        # P = x*p6 + y*p4 + z*p1 + xy*(p7-p6-p4) + xz*(p2-p1-p6) + yz*(p0-p4-p1) + xyz(p3-p7-p2-p0+p1+p6+p4)
        x,y,z = coo[:,0],coo[:,1],coo[:,2]
        coo2 = x*pt_6 + y*pt_4 + z*pt_1 + \
               x*y*(pt_7 - pt_6 - pt_4) + x*z*(pt_2 - pt_1 - pt_6) + y*z*(pt_0 - pt_4 - pt_1) + \
               x*y*z*(pt_3 - pt_7 - pt_2 - pt_0 + pt_6 + pt_1 + pt_4) + orig
        tgt.setCoords(coo2)

        sCoo = DataArrayDouble([-0.3,1.0,1.0,  0.0,0.0,1.0,  1.0,0.0,1.0,  2.0,3.0,1.0,
                                -0.3,1.0,0.0,  0.0,0.0,0.0,  1.0,0.0,0.0,  2.0,3.0,0.0,],8,3)
        sCoo[:, 0] += 10.0; sCoo[:, 1] += 15.0; sCoo[:, 2] += 20.0;
        sConn = DataArrayInt([0,1,2,3,4, 5,6,7])
        s = MEDCoupling1SGTUMesh("source",NORM_HEXA8) ; s.setCoords(sCoo)
        s.setNodalConnectivity(sConn)
        #
        aRemapper=MEDCouplingRemapper()
        aRemapper.setPrecision(1e-12)
        aRemapper.setIntersectionType(MappedBarycentric)
        self.assertEqual(aRemapper.prepare(s,tgt,'P1P1'),1)
        srcField = MEDCouplingFieldDouble(ON_NODES, ONE_TIME)
        srcField.setNature(IntensiveMaximum)
        srcField.setMesh(s); srcField.setName("field")
        srcField.setArray(DataArrayDouble([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0]))
        tgtF = aRemapper.transferField(srcField, 1e+300)
#        print tgtF.getArray().getValues()
        ref = [6.0, 6.251802698104413, 6.502397834044702, 6.7517940736426665, 7.0, 5.740554726834594,
               6.1761835575796935, 6.6052985689637564, 7.009392769824465, 7.383488834310164,
               5.487562931129931, 6.140664596972973, 6.720290674177548, 7.220534970454015, 7.651092836860121,
               5.2407867837524345, 6.125759809889516, 6.82853486793175, 7.390880823876876, 7.848445254819061,
               5.0, 6.12211344611157, 6.925740671133115, 7.529623182840827, 8.0, 5.0, 5.251802698104413,
               5.502397834044702, 5.751794073642667, 6.0, 4.740554726834594, 5.1761835575796935,
               5.6052985689637564, 6.009392769824465, 6.383488834310163, 4.487562931129931, 5.140664596972973,
                5.720290674177548, 6.220534970454015, 6.651092836860121, 4.2407867837524345, 5.125759809889516,
                5.828534867931749, 6.390880823876876, 6.848445254819061, 4.0, 5.122113446111569, 5.925740671133115,
                6.529623182840827, 7.0, 4.0, 4.251802698104413, 4.502397834044702, 4.751794073642667, 5.0, 3.740554726834594,
                4.176183557579693, 4.6052985689637564, 5.009392769824464, 5.383488834310164, 3.487562931129931,
                4.140664596972973, 4.720290674177548, 5.220534970454015, 5.651092836860121, 3.240786783752434, 4.125759809889516, 4.82853486793175,
                5.390880823876876, 5.848445254819061, 3.0, 4.122113446111569, 4.925740671133115, 5.529623182840827, 6.0, 3.0,
                3.2518026981044135, 3.502397834044702, 3.7517940736426674, 4.0, 2.7405547268345933, 3.176183557579693,
                3.6052985689637564, 4.009392769824465, 4.383488834310164, 2.487562931129931, 3.140664596972973, 3.7202906741775474, 4.220534970454015, 4.65109283686012, 2.2407867837524345, 3.1257598098895154, 3.828534867931749,
                4.390880823876876, 4.848445254819061, 2.0, 3.1221134461115687, 3.9257406711331146, 4.529623182840826, 5.0, 2.0, 2.2518026981044135, 2.502397834044702, 2.7517940736426674, 3.0, 1.7405547268345936, 2.176183557579693, 2.6052985689637564,
                3.0093927698244642, 3.3834888343101635, 1.4875629311299305, 2.1406645969729734, 2.720290674177548,
                3.2205349704540143, 3.6510928368601205, 1.2407867837524345, 2.125759809889516, 2.8285348679317495, 3.390880823876876, 3.848445254819061, 1.0, 2.1221134461115687, 2.9257406711331146, 3.529623182840827, 4.0]

        val = tgtF.getArray().getValues()
        for i, ref_v in enumerate(ref):
            self.assertAlmostEqual(ref_v, val[i])
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings() and MEDCouplingHasSciPyBindings(),"requires numpy AND scipy")
    def testGetCrudeCSRMatrix1(self):
        """ testing CSR matrix output using numpy/scipy.
        """
        from scipy.sparse import spdiags #diags
        import scipy
        from numpy import array
        arr=DataArrayDouble(3) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr)
        src=m.buildUnstructured()
        trg=src.deepCopy() ; trg=trg[[0,1,3]]
        trg.getCoords()[:]*=0.5 ; trg.getCoords()[:]+=[0.3,0.25]
        # Let's interpolate.
        rem=MEDCouplingRemapper()
        rem.prepare(src,trg,"P0P0")
        # Internal crude sparse matrix computed. Let's manipulate it using CSR matrix in scipy.
        for i in range(10):
            m=rem.getCrudeCSRMatrix()
            pass
        m2=rem.getCrudeCSRMatrix()
        diff=m-m2
        self.assertTrue(isinstance(m,scipy.sparse.csr.csr_matrix))
        self.assertEqual(m.getnnz(),7)
        self.assertAlmostEqual(m[0,0],0.25,12)
        self.assertAlmostEqual(m[1,0],0.1,12)
        self.assertAlmostEqual(m[1,1],0.15,12)
        self.assertAlmostEqual(m[2,0],0.05,12)
        self.assertAlmostEqual(m[2,1],0.075,12)
        self.assertAlmostEqual(m[2,2],0.05,12)
        self.assertAlmostEqual(m[2,3],0.075,12)
        self.assertEqual(diff.getnnz(),0)
        # ExtensiveConservation (division by sum of cols)
        colSum=m.sum(axis=0)
        # version 0.12.0 # m_0=m*diags(array(1/colSum),[0])
        m_0=m*spdiags(array(1/colSum),[0],colSum.shape[1],colSum.shape[1])
        del colSum
        self.assertAlmostEqual(m_0[0,0],0.625,12)
        self.assertAlmostEqual(m_0[1,0],0.25,12)
        self.assertAlmostEqual(m_0[1,1],0.6666666666666667,12)
        self.assertAlmostEqual(m_0[2,0],0.125,12)
        self.assertAlmostEqual(m_0[2,1],0.3333333333333333,12)
        self.assertAlmostEqual(m_0[2,2],1.,12)
        self.assertAlmostEqual(m_0[2,3],1.,12)
        self.assertEqual(m_0.getnnz(),7)
        # IntensiveMaximum (division by sum of rows)
        rowSum=m.sum(axis=1)
        # version 0.12.0 # m_1=diags(array(1/rowSum.transpose()),[0])*m
        m_1=spdiags(array(1/rowSum.transpose()),[0],rowSum.shape[0],rowSum.shape[0])*m
        del rowSum
        self.assertAlmostEqual(m_1[0,0],1.,12)
        self.assertAlmostEqual(m_1[1,0],0.4,12)
        self.assertAlmostEqual(m_1[1,1],0.6,12)
        self.assertAlmostEqual(m_1[2,0],0.2,12)
        self.assertAlmostEqual(m_1[2,1],0.3,12)
        self.assertAlmostEqual(m_1[2,2],0.2,12)
        self.assertAlmostEqual(m_1[2,3],0.3,12)
        self.assertEqual(m_1.getnnz(),7)
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings() and MEDCouplingHasSciPyBindings(),"requires numpy AND scipy")
    def testP0P1Bary_1(self):
        a=MEDCouplingUMesh("a",2)
        a.allocateCells()
        conna=[0,1,3,2,1,4,5,3,4,6,7,5,6,8,9,7,8,10,11,9,10,12,13,11,12,14,15,13,14,16,17,15,16,18,19,17,18,20,21,19,20,22,23,21,22,24,25,23,24,26,27,25]
        a.setCoords(DataArrayDouble([1.54,0,-0.01,1.54,0.02,-0.01,1.54,0,0.01,1.54,0.02,0.01,1.54,0.04,-0.01,1.54,0.04,0.01,1.54,0.06,-0.01,1.54,0.06,0.01,1.54,0.08,-0.01,1.54,0.08,0.01,1.54,0.1,-0.01,1.54,0.1,0.01,1.54,0.12,-0.01,1.54,0.12,0.01,1.54,0.14,-0.01,1.54,0.14,0.01,1.54,0.16,-0.01,1.54,0.16,0.01,1.54,0.18,-0.01,1.54,0.18,0.01,1.54,0.2,-0.01,1.54,0.2,0.01,1.54,0.22,-0.01,1.54,0.22,0.01,1.54,0.24,-0.01,1.54,0.24,0.01,1.54,0.26,-0.01,1.54,0.26,0.01],28,3))
        for i in range(13):
            a.insertNextCell(NORM_QUAD4,conna[4*i:4*(i+1)])
            pass
        a.finishInsertingCells() ; a.simplexize(0)
        #
        connb=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,0,2,39,3,5,40,6,8,41,9,11,42,12,14,43,15,17,44,18,20,45,21,23,46,24,26,47,27,29,48,30,32,49,33,35,50,36,38,51,52,2,39,53,5,40,54,8,41,55,11,42,56,14,43,57,17,44,58,20,45,59,23,46,60,26,47,61,29,48,62,32,49,63,35,50,64,38,51,52,2,65,53,5,66,54,8,67,55,11,68,56,14,69,57,17,70,58,20,71,59,23,72,60,26,73,61,29,74,62,32,75,63,35,76,64,38,77,53,2,65,54,5,66,55,8,67,56,11,68,57,14,69,58,17,70,59,20,71,60,23,72,61,26,73,62,29,74,63,32,75,64,35,76,78,38,77,53,2,40,54,5,41,55,8,42,56,11,43,57,14,44,58,17,45,59,20,46,60,23,47,61,26,48,62,29,49,63,32,50,64,35,51,78,38,79,3,2,40,6,5,41,9,8,42,12,11,43,15,14,44,18,17,45,21,20,46,24,23,47,27,26,48,30,29,49,33,32,50,36,35,51,80,38,79,3,2,1,6,5,4,9,8,7,12,11,10,15,14,13,18,17,16,21,20,19,24,23,22,27,26,25,30,29,28,33,32,31,36,35,34,80,38,37]
        b=MEDCouplingUMesh("b",2)
        b.allocateCells()
        for i in range(104):
            b.insertNextCell(NORM_TRI3,connb[3*i:3*(i+1)])
            pass
        b.setCoords(DataArrayDouble([1.54,0,-0.01,1.54,0.01,-0.01,1.54,0.01,0,1.54,0.02,-0.01,1.54,0.03,-0.01,1.54,0.03,0,1.54,0.04,-0.01,1.54,0.05,-0.01,1.54,0.05,0,1.54,0.06,-0.01,1.54,0.07,-0.01,1.54,0.07,0,1.54,0.08,-0.01,1.54,0.09,-0.01,1.54,0.09,0,1.54,0.1,-0.01,1.54,0.11,-0.01,1.54,0.11,0,1.54,0.12,-0.01,1.54,0.13,-0.01,1.54,0.13,0,1.54,0.14,-0.01,1.54,0.15,-0.01,1.54,0.15,0,1.54,0.16,-0.01,1.54,0.17,-0.01,1.54,0.17,0,1.54,0.18,-0.01,1.54,0.19,-0.01,1.54,0.19,0,1.54,0.2,-0.01,1.54,0.21,-0.01,1.54,0.21,0,1.54,0.22,-0.01,1.54,0.23,-0.01,1.54,0.23,0,1.54,0.24,-0.01,1.54,0.25,-0.01,1.54,0.25,0,1.54,0,0,1.54,0.02,0,1.54,0.04,0,1.54,0.06,0,1.54,0.08,0,1.54,0.1,0,1.54,0.12,0,1.54,0.14,0,1.54,0.16,0,1.54,0.18,0,1.54,0.2,0,1.54,0.22,0,1.54,0.24,0,1.54,0,0.01,1.54,0.02,0.01,1.54,0.04,0.01,1.54,0.06,0.01,1.54,0.08,0.01,1.54,0.1,0.01,1.54,0.12,0.01,1.54,0.14,0.01,1.54,0.16,0.01,1.54,0.18,0.01,1.54,0.2,0.01,1.54,0.22,0.01,1.54,0.24,0.01,1.54,0.01,0.01,1.54,0.03,0.01,1.54,0.05,0.01,1.54,0.07,0.01,1.54,0.09,0.01,1.54,0.11,0.01,1.54,0.13,0.01,1.54,0.15,0.01,1.54,0.17,0.01,1.54,0.19,0.01,1.54,0.21,0.01,1.54,0.23,0.01,1.54,0.25,0.01,1.54,0.26,0.01,1.54,0.26,0,1.54,0.26,-0.01],81,3))
        #
        rem=MEDCouplingRemapper() ; rem.setIntersectionType(Barycentric)
        rem.prepare(a,b,"P1P0")
        m0=rem.getCrudeCSRMatrix()
        self.assertEqual(m0.nnz,312)
        #
        ids=4*[None] ; vs=4*[None]
        ids[0]=DataArrayInt([0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,158,161,164,167,170,173,176,179,182,185,188,191,194,197,200,203,206,209,212,215,218,221,224,227,230,233]) ; vs[0]=10./3.
        ids[1]=DataArrayInt([1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35,37,38,40,41,43,44,46,47,49,50,52,53,55,56,58,59,61,62,64,65,67,68,70,71,73,74,76,77,80,83,86,89,92,95,98,101,104,107,110,113,116,117,120,123,126,129,132,135,138,141,144,147,150,153,156,157,159,160,162,163,165,166,168,169,171,172,174,175,177,178,180,181,183,184,186,187,189,190,192,193,195,196,198,199,201,202,204,205,207,208,210,211,213,214,216,217,219,220,222,223,225,226,228,229,231,232,234,237,240,243,246,249,252,255,258,261,264,267,270,275,278,281,284,287,290,293,296,299,302,305,308,311]) ; vs[1]=5./6.
        ids[2]=DataArrayInt([78,81,84,87,90,93,96,99,102,105,108,111,114,119,122,125,128,131,134,137,140,143,146,149,152,155,236,239,242,245,248,251,254,257,260,263,266,269,272,273,276,279,282,285,288,291,294,297,300,303,306,309]) ; vs[2]=5./3.
        ids[3]=DataArrayInt([79,82,85,88,91,94,97,100,103,106,109,112,115,118,121,124,127,130,133,136,139,142,145,148,151,154,235,238,241,244,247,250,253,256,259,262,265,268,271,274,277,280,283,286,289,292,295,298,301,304,307,310]) ; vs[3]=2.5
        vals=DataArrayDouble(312,1)
        for idd,v in zip(ids,vs):
            vals[idd]=v
            pass
        vals*=1e-5
        eps0=DataArrayDouble(m0.data)-vals ; eps0.abs()
        self.assertTrue(eps0.findIdsInRange(1e-17,1e300).empty())
        self.assertTrue(DataArrayInt(m0.indices).isEqual(DataArrayInt([0,1,3,1,4,5,4,6,7,6,8,9,8,10,11,10,12,13,12,14,15,14,16,17,16,18,19,18,20,21,20,22,23,22,24,25,24,26,27,0,2,3,1,3,5,4,5,7,6,7,9,8,9,11,10,11,13,12,13,15,14,15,17,16,17,19,18,19,21,20,21,23,22,23,25,24,25,27,0,2,3,1,3,5,4,5,7,6,7,9,8,9,11,10,11,13,12,13,15,14,15,17,16,17,19,18,19,21,20,21,23,22,23,25,24,25,27,0,2,3,1,3,5,4,5,7,6,7,9,8,9,11,10,11,13,12,13,15,14,15,17,16,17,19,18,19,21,20,21,23,22,23,25,24,25,27,0,2,3,1,3,5,4,5,7,6,7,9,8,9,11,10,11,13,12,13,15,14,15,17,16,17,19,18,19,21,20,21,23,22,23,25,24,25,27,0,1,3,1,4,5,4,6,7,6,8,9,8,10,11,10,12,13,12,14,15,14,16,17,16,18,19,18,20,21,20,22,23,22,24,25,24,26,27,0,1,3,1,4,5,4,6,7,6,8,9,8,10,11,10,12,13,12,14,15,14,16,17,16,18,19,18,20,21,20,22,23,22,24,25,24,26,27,0,1,3,1,4,5,4,6,7,6,8,9,8,10,11,10,12,13,12,14,15,14,16,17,16,18,19,18,20,21,20,22,23,22,24,25,24,26,27])))
        self.assertTrue(DataArrayInt(m0.indptr).isEqual(DataArrayInt([0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102,105,108,111,114,117,120,123,126,129,132,135,138,141,144,147,150,153,156,159,162,165,168,171,174,177,180,183,186,189,192,195,198,201,204,207,210,213,216,219,222,225,228,231,234,237,240,243,246,249,252,255,258,261,264,267,270,273,276,279,282,285,288,291,294,297,300,303,306,309,312])))
        #
        rem2=MEDCouplingRemapper() ; rem2.setIntersectionType(Barycentric)
        rem2.prepare(b,a,"P0P1")
        m1=rem2.getCrudeCSRMatrix()
        self.assertEqual(m1.nnz,312)
        #
        m1=rem2.getCrudeCSRMatrix()
        m1t=m1.transpose()
        delta=m0-m1t
        self.assertTrue(DataArrayDouble(delta.data).isUniform(0.,1e-17))
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings() and MEDCouplingHasSciPyBindings(),"requires numpy AND scipy")
    def testNonConformWithRemapper_1(self):
        coo=DataArrayDouble([-0.396700000780411,-0.134843245350081,-0.0361311386958691,-0.407550009429364,-0.13484324535008,-0.0361311386958923,-0.396700000780411,-0.132191446077668,-0.0448729493559049,-0.407550009429364,-0.132191446077666,-0.0448729493559254,-0.396700000780411,-0.128973582738749,-0.0534226071577727,-0.407550009429364,-0.128973582738747,-0.0534226071577904,-0.396700000780411,-0.128348829636458,-0.0346583696473619,-0.407550009429364,-0.128348829636457,-0.0346583696473822,-0.396700000780411,-0.125874740261886,-0.0430683597970123,-0.407550009429364,-0.125874740261885,-0.0430683597970302,-0.396700000780411,-0.122905344829122,-0.051310216195766,-0.407550009429364,-0.12290534482912,-0.0513102161957814],12,3)
        conn=DataArrayInt([2,9,3,11,2,3,5,11,2,8,9,11,2,10,8,11,2,5,4,11,2,4,10,11,3,0,1,6,3,1,7,6,3,2,0,6,3,8,2,6,3,7,9,6,3,9,8,6])
        m=MEDCoupling1SGTUMesh("mesh",NORM_TETRA4)
        m.setNodalConnectivity(conn)
        m.setCoords(coo)
        # m is ready
        m1,d,di,rd,rdi=m.buildUnstructured().buildDescendingConnectivity()
        rdi2=rdi.deltaShiftIndex()
        cellIds=rdi2.findIdsEqual(1)
        skinAndNonConformCells=m1[cellIds]
        skinAndNonConformCells.zipCoords() # at this point skinAndNonConformCells contains non conform cells and skin cells. Now trying to split them in two parts.
        #
        rem=MEDCouplingRemapper()
        rem.setMaxDistance3DSurfIntersect(1e-12)
        rem.setMinDotBtwPlane3DSurfIntersect(0.99)# this line is important it is to tell to remapper to select only cells with very close orientation
        rem.prepare(skinAndNonConformCells,skinAndNonConformCells,"P0P0")
        mat=rem.getCrudeCSRMatrix()
        indptr=DataArrayInt(mat.indptr)
        indptr2=indptr.deltaShiftIndex()
        cellIdsOfNonConformCells=indptr2.findIdsNotEqual(1)
        cellIdsOfSkin=indptr2.findIdsEqual(1)
        self.assertTrue(cellIdsOfSkin.isEqual(DataArrayInt([1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,23])))
        self.assertTrue(cellIdsOfNonConformCells.isEqual(DataArrayInt([0,4,18,22])))
        pass

    def test3D1DOnP1P0_1(self):
        """ This test focused on P1P0 interpolation with a source with meshDim=1 spaceDim=3 and a target with meshDim=3.
        This test has revealed a bug in remapper. A reverse matrix is computed so a reverse method should be given in input.
        """
        target=MEDCouplingCMesh()
        arrX=DataArrayDouble([0,1]) ; arrY=DataArrayDouble([0,1]) ; arrZ=DataArrayDouble(11) ; arrZ.iota()
        target.setCoords(arrX,arrY,arrZ)
        target=target.buildUnstructured() ; target.setName("TargetSecondaire")
        #
        sourceCoo=DataArrayDouble([(0.5,0.5,0.1),(0.5,0.5,1.2),(0.5,0.5,1.6),(0.5,0.5,1.8),(0.5,0.5,2.43),(0.5,0.5,2.55),(0.5,0.5,4.1),(0.5,0.5,4.4),(0.5,0.5,4.9),(0.5,0.5,5.1),(0.5,0.5,7.6),(0.5,0.5,7.7),(0.5,0.5,8.2),(0.5,0.5,8.4),(0.5,0.5,8.6),(0.5,0.5,8.8),(0.5,0.5,9.2),(0.5,0.5,9.6),(0.5,0.5,11.5)])
        source=MEDCoupling1SGTUMesh("SourcePrimaire",NORM_SEG2)
        source.setCoords(sourceCoo)
        source.allocateCells()
        for i in range(len(sourceCoo) - 1):
            source.insertNextCell([i,i+1])
            pass
        source=source.buildUnstructured()
        fsource=MEDCouplingFieldDouble(ON_NODES) ; fsource.setName("field")
        fsource.setMesh(source)
        arr=DataArrayDouble(len(sourceCoo)) ; arr.iota(0.7) ; arr*=arr
        fsource.setArray(arr)
        fsource.setNature(IntensiveMaximum)
        #
        rem=MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(source,target,"P1P0")
        f2Test=rem.transferField(fsource,-27)
        self.assertEqual(f2Test.getName(),fsource.getName())
        self.assertEqual(f2Test.getMesh().getHiddenCppPointer(),target.getHiddenCppPointer())
        expArr=DataArrayDouble([0.49,7.956666666666667,27.29,-27,59.95666666666667,94.09,-27,125.69,202.89,296.09])
        self.assertTrue(f2Test.getArray().isEqual(expArr,1e-12))
        f2Test=rem.reverseTransferField(f2Test,-36)
        self.assertEqual(f2Test.getName(),fsource.getName())
        self.assertEqual(f2Test.getMesh().getHiddenCppPointer(),source.getHiddenCppPointer())
        expArr2=DataArrayDouble([0.49,7.956666666666667,7.956666666666667,7.956666666666667,27.29,27.29,59.95666666666667,59.95666666666667,59.95666666666667,94.09,125.69,125.69,202.89,202.89,202.89,202.89,296.09,296.09,-36.])
        self.assertTrue(f2Test.getArray().isEqual(expArr2,1e-12))
        pass

    def testRemapperAMR1(self):
        """ This test is the origin of the ref values for MEDCouplingBasicsTest.testAMR2"""
        coarse=DataArrayDouble(35) ; coarse.iota(0) #X=5,Y=7
        fine=DataArrayDouble(3*2*4*4) ; fine.iota(0) #X=3,Y=2 refined by 4
        MEDCouplingIMesh.CondenseFineToCoarse([5,7],fine,[(1,4),(2,4)],[4,4],coarse)
        #
        m=MEDCouplingCartesianAMRMesh("mesh",2,[6,8],[0.,0.],[1.,1.])
        trgMesh=m.buildUnstructured()
        m.addPatch([(1,4),(2,4)],[4,4])
        srcMesh=m[0].getMesh().buildUnstructured()
        srcField=MEDCouplingFieldDouble(ON_CELLS)
        fine2=DataArrayDouble(3*2*4*4) ; fine2.iota(0) ; srcField.setArray(fine2)
        srcField.setMesh(srcMesh) ; srcField.setNature(ExtensiveMaximum)
        #
        trgField=MEDCouplingFieldDouble(ON_CELLS)
        coarse2=DataArrayDouble(35) ; coarse2.iota(0) ; trgField.setArray(coarse2)
        trgField.setMesh(trgMesh) ; trgField.setNature(ExtensiveMaximum)
        #
        rem=MEDCouplingRemapper()
        rem.prepare(srcMesh,trgMesh,"P0P0")
        rem.partialTransfer(srcField,trgField)
        #
        self.assertTrue(coarse.isEqual(trgField.getArray(),1e-12))
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings() and MEDCouplingHasSciPyBindings(),"requires numpy AND scipy")
    def test1DPointLocator1(self):
        """This test focuses on PointLocator for P1P1 in 1D and 2DCurve."""
        from numpy import array
        from scipy.sparse import diags,csr_matrix,identity
        ## basic case 1D
        arrS=DataArrayInt.Range(0,11,1).convertToDblArr()
        arrT=DataArrayDouble([0.1,1.7,5.5,9.6])
        mS=MEDCouplingCMesh() ; mS.setCoords(arrS)
        mT=MEDCouplingCMesh() ; mT.setCoords(arrT)
        rem=MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        self.assertEqual(rem.prepare(mS.buildUnstructured(),mT.buildUnstructured(),"P1P1"),1)
        m=rem.getCrudeCSRMatrix()
        rowSum=m.sum(axis=1)
        m=diags(array(1/rowSum.transpose()),[0])*m
        # expected matrix
        row=array([0,0,1,1,2,2,3,3])
        col=array([0,1,1,2,5,6,9,10])
        data=array([0.9,0.1,0.3,0.7,0.5,0.5,0.4,0.6])
        mExp0=csr_matrix((data,(row,col)),shape=(4,11))
        # compute diff and check
        diff=abs(m-mExp0)
        self.assertAlmostEqual(diff.sum(),0.,14)
        ## full specific case 1D where target=source
        rem=MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        self.assertEqual(rem.prepare(mS.buildUnstructured(),mS.buildUnstructured(),"P1P1"),1)
        m=rem.getCrudeCSRMatrix()
        rowSum=m.sum(axis=1)
        m=diags(array(1/rowSum.transpose()),[0])*m
        # expected matrix
        mExp1=identity(11)
        diff=abs(m-mExp1)
        self.assertAlmostEqual(diff.sum(),0.,14)
        ## case where some points in target are not in source
        arrT=DataArrayDouble([-0.2,0.1,1.7,5.5,10.3])
        mT=MEDCouplingCMesh() ; mT.setCoords(arrT)
        mT=mT.buildUnstructured()
        rem=MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        self.assertEqual(rem.prepare(mS.buildUnstructured(),mT,"P1P1"),1)
        m=rem.getCrudeCSRMatrix()
        row=array([1,1,2,2,3,3])
        col=array([0,1,1,2,5,6])
        data=array([0.9,0.1,0.3,0.7,0.5,0.5])
        mExp2=csr_matrix((data,(row,col)),shape=(5,11))
        diff=abs(m-mExp2)
        self.assertAlmostEqual(diff.sum(),0.,14)
        ## basic case 2D Curve
        arrS=DataArrayInt.Range(0,11,1).convertToDblArr()
        arrT=DataArrayDouble([0.1,1.7,5.5,9.6])
        mS=MEDCouplingCMesh() ; mS.setCoords(arrS)
        mT=MEDCouplingCMesh() ; mT.setCoords(arrT)
        mS=mS.buildUnstructured() ; mS.changeSpaceDimension(2)
        mT=mT.buildUnstructured() ; mT.changeSpaceDimension(2)
        mS.rotate([-1.,-1.],1.2)
        mT.rotate([-1.,-1.],1.2)
        rem=MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        self.assertEqual(rem.prepare(mS,mT,"P1P1"),1)
        m=rem.getCrudeCSRMatrix()
        rowSum=m.sum(axis=1)
        m=diags(array(1/rowSum.transpose()),[0])*m
        diff=abs(m-mExp0)
        self.assertAlmostEqual(diff.sum(),0.,14)
        pass

    def test3D2Dand2D3DPointLocator1(self):
        """ Non regression test solving SIGSEGV when using 3D<->3Dsurf pointlocator."""
        arrX=DataArrayDouble([0,1,2])
        arrY=DataArrayDouble([0,1])
        arrZ=DataArrayDouble([0,1])
        ms=MEDCouplingCMesh() ; ms.setCoords(arrX,arrY,arrZ)
        ms=ms.buildUnstructured() ; ms.setName("source")
        #
        mt=MEDCouplingUMesh("target",2) ; mt.allocateCells()
        mt.insertNextCell(NORM_TRI3,[0,4,6])
        mt.insertNextCell(NORM_TRI3,[1,5,7])
        mt.setCoords(ms.getCoords()[:])
        mt.zipCoords()
        #
        rem=MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(ms,mt,"P0P0")
        self.assertEqual(rem.getCrudeMatrix(),[{0: 1.0}, {1: 1.0}])
        rem2=MEDCouplingRemapper()
        rem2.setIntersectionType(PointLocator)
        rem2.prepare(mt,ms,"P0P0") # reverse mt<->ms
        self.assertEqual(rem2.getCrudeMatrix(),[{0: 1.0}, {1: 1.0}])
        pass

    def test2D1Dand1D2DPointLocator1(self):
        arrX=DataArrayDouble([0,1,2])
        arrY=DataArrayDouble([0,1])
        ms=MEDCouplingCMesh() ; ms.setCoords(arrX,arrY) ; ms=ms.buildUnstructured()
        mt=MEDCouplingUMesh("target",1) ; mt.setCoords(ms.getCoords()[:])
        mt.allocateCells()
        mt.insertNextCell(NORM_SEG2,[0,4]) ; mt.insertNextCell(NORM_SEG2,[1,5])
        rem=MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(ms,mt,"P0P0")
        self.assertEqual(rem.getCrudeMatrix(),[{0:1.},{1:1.}])
        rem=MEDCouplingRemapper()
        rem.setIntersectionType(PointLocator)
        rem.prepare(mt,ms,"P0P0")
        self.assertEqual(rem.getCrudeMatrix(),[{0:1.},{1:1.}])
        pass

    def test3D1DPointLocatorBBoxAdjusted(self):
        """ In case a 1D segment lies exactly on the interface between two 2D (or 3D) faces, the default
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
        m1d.setCoords(dd([1.0,0.5,0.2  ,  1.0,0.5,0.8], 2,3))
        m1d.setConnectivity(di([NORM_SEG2, 0, 1]), di([0,3]))

        rem = MEDCouplingRemapper()
        rem.setPrecision(1e-12)
        rem.setIntersectionType(PointLocator)
        rem.prepare(m3d, m1d,"P0P1")
        self.assertEqual(rem.getCrudeMatrix(), [{0: 1.0, 1: 1.0}, {0: 1.0, 1: 1.0}])

        rem = MEDCouplingRemapper()
        rem.setPrecision(1e-12)
        rem.setIntersectionType(PointLocator)
        rem.setBoundingBoxAdjustment(0.0)
        rem.setBoundingBoxAdjustmentAbs(0.0)
        rem.prepare(m3d, m1d,"P0P1")
        self.assertEqual(rem.getCrudeMatrix(), [{}, {}])
        pass

    def testPointLocator3DTo2D(self):
        """Target mesh has spaceDim==3 and meshDim==2. Source has spaceDim==3 and meshDim==3. Here we are on pointlocator alg.
        The test evaluates on each nodes of target mesh the bary coor into source mesh."""
        src=MEDCouplingCMesh()
        arr=DataArrayDouble([0,1,2])
        src.setCoords(arr,arr,arr)
        src=src.buildUnstructured()
        src.simplexize(PLANAR_FACE_5)
        fsrc=MEDCouplingFieldDouble(ON_NODES) ; fsrc.setMesh(src)
        fsrc.setArray(DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]))
        #
        trg=MEDCouplingCMesh()
        arr=DataArrayDouble([0,1])
        trg.setCoords(arr,arr)
        trg=trg.buildUnstructured()
        trg.changeSpaceDimension(3,0.)
        trg.translate([0.5,0.5,0.5])
        #
        arrTrg=fsrc.getValueOnMulti(trg.getCoords())
        ftrg=MEDCouplingFieldDouble(ON_NODES)
        ftrg.setMesh(trg)
        ftrg.setArray(arrTrg)
        ftrg.checkConsistencyLight()
        ftrg.setNature(IntensiveMaximum)
        #
        fsrc.setNature(IntensiveMaximum)
        remap=MEDCouplingRemapper()
        remap.setIntersectionType(PointLocator)
        self.assertEqual(remap.prepare(src,trg,"P1P1"),1)
        ftrg2=remap.transferField(fsrc,1e300)
        self.assertTrue(ftrg.isEqual(ftrg2,1e-12,1e-12))
        pass

    def testExtrudedOnDiffZLev1(self):
        """Non regression bug : This test is base on P0P0 ExtrudedExtruded. This test checks that if the input meshes are not based on a same plane // OXY the interpolation works"""
        arrX=DataArrayDouble([0,1]) ; arrY=DataArrayDouble([0,1]) ; arrZ=DataArrayDouble([0,1,2])
        src=MEDCouplingCMesh() ; src.setCoords(arrX,arrY,arrZ)
        arrX=DataArrayDouble([0.5,1.5]) ; arrY=DataArrayDouble([0.5,1.5]) ; arrZ=DataArrayDouble([0.5,2])
        trg=MEDCouplingCMesh() ; trg.setCoords(arrX,arrY,arrZ)
        #
        src=MEDCouplingMappedExtrudedMesh(src) ; trg=MEDCouplingMappedExtrudedMesh(trg)
        pt1=src.getMesh2D().getCoords().getHiddenCppPointer() ; pt2=trg.getMesh2D().getCoords().getHiddenCppPointer()
        #
        rem=MEDCouplingRemapper()
        rem.prepare(src,trg,"P0P0")
        self.checkMatrix(rem.getCrudeMatrix(),[{0:0.125,1:0.25}],src.getNumberOfCells(),1e-12)
        #
        self.assertEqual(src.getMesh2D().getSpaceDimension(),3)
        self.assertEqual(trg.getMesh2D().getSpaceDimension(),3)
        self.assertEqual(src.getMesh2D().getCoords().getHiddenCppPointer(),pt1)
        self.assertEqual(trg.getMesh2D().getCoords().getHiddenCppPointer(),pt2)
        #
        rem2=MEDCouplingRemapper()
        rem2.setIntersectionType(Geometric2D)
        rem2.prepare(src,trg,"P0P0")
        self.checkMatrix(rem2.getCrudeMatrix(),[{0:0.125,1:0.25}],src.getNumberOfCells(),1e-12)
        pass

    def testP0P0WithHEXGP12(self):
        """ Test that HEXGP12 are correclty remapped (elements with polygonal faces were not properly handled) """
        # From Astrid, two disjoint hexagonal prisms:
        coo1 = [-4.991193077144312, 8.644999999999998, 0.0, -9.982386154288623, 6.112246755425186e-16, 0.0, -4.991193077144315, -8.644999999999998, 0.0, 4.991193077144309, -8.645000000000005, 0.0, 9.982386154288626, 1.1651321638577316e-15, 0.0, 4.991193077144314, 8.645, 0.0, -4.991193077144312, 8.644999999999998, 7.561799999999991, -9.982386154288623, 6.112246755425186e-16, 7.561799999999991, -4.991193077144315, -8.644999999999998, 7.561799999999991, 4.991193077144309, -8.645000000000005, 7.561799999999991, 9.982386154288626, 1.1651321638577316e-15, 7.561799999999991, 4.991193077144314, 8.645, 7.561799999999991]
        coo2 = [-4.991193077144313, -8.645, 0.0, -9.982386154288626, -1.3992140779350848e-15, 0.0, -19.964772308577256, 0.0, 0.0, -24.95596538572157, -8.644999999999998, 0.0, -19.96477230857726, -17.289999999999996, 0.0, -9.982386154288626, -17.289999999999996, 0.0, -4.991193077144313, -8.645, 5.041200000000004, -9.982386154288626, -1.3992140779350848e-15, 5.041200000000004, -19.964772308577256, 0.0, 5.041200000000004, -24.95596538572157, -8.644999999999998, 5.041200000000004, -19.96477230857726, -17.289999999999996, 5.041200000000004, -9.982386154288626, -17.289999999999996, 5.041200000000004]
        conn1 = [31, 0, 5, 4, 3, 2, 1, -1, 11, 6, 7, 8, 9, 10, -1, 1, 7, 6, 0, -1, 2, 8, 7, 1, -1, 3, 9, 8, 2, -1, 4, 10, 9, 3, -1, 5, 11, 10, 4, -1, 0, 6, 11, 5]
        cI1 = [0, 44]
        conn2 = [31, 0, 5, 4, 3, 2, 1, -1, 6, 7, 8, 9, 10, 11, -1, 0, 1, 7, 6, -1, 1, 2, 8, 7, -1, 2, 3, 9, 8, -1, 3, 4, 10, 9, -1, 4, 5, 11, 10, -1, 5, 0, 6, 11]
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
        """ The killer tetrahedron detected by LMEC!"""
        mesh = MEDCouplingUMesh('SupportOf_ECHIA1_Tin', 3)
#         # was OK:
#         coo = DataArrayDouble([(-4.50135,1.95352,4.59608),(-4.50409,1.86642,4.54551), (-4.55175,1.92167,4.64844),(-4.58813,1.94795,4.5283)])
        # was KO:
        coo = DataArrayDouble([(-4.501352938826142847,1.953517433537110159,4.596082552008083688),(-4.504092113061189728,1.866415526007169978,4.545507396150389567),(-4.551750368181751050,1.921669328035479962,4.648439577911889664),(-4.588131417812300050,1.947948377683889953,4.528298931319220344)])
        mesh.setCoords(coo)
        c = DataArrayInt([14, 2, 0, 3, 1]); cI = DataArrayInt([0, 5])
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

    def checkMatrix(self,mat1,mat2,nbCols,eps):
        self.assertEqual(len(mat1),len(mat2))
        for i in range(len(mat1)):
            self.assertTrue(max(mat2[i].keys())<nbCols)
            self.assertTrue(max(mat1[i].keys())<nbCols)
            self.assertTrue(min(mat2[i].keys())>=0)
            self.assertTrue(min(mat1[i].keys())>=0)
            s1=set(mat1[i].keys()) ; s2=set(mat2[i].keys())
            for elt in s1.intersection(s2):
                self.assertTrue(abs(mat1[i][elt]-mat2[i][elt])<eps)
                pass
            for elt in s1.difference(s2):
                self.assertTrue(abs(mat1[i][elt])<eps)
                pass
            for elt in s2.difference(s1):
                self.assertTrue(abs(mat2[i][elt])<eps)
                pass
            pass
        pass

    def build2DSourceMesh_1(self):
        sourceCoords=[-0.3,-0.3, 0.7,-0.3, -0.3,0.7, 0.7,0.7]
        sourceConn=[0,3,1,0,2,3]
        sourceMesh=MEDCouplingUMesh.New("my name of mesh 2D",2)
        sourceMesh.allocateCells(2);
        sourceMesh.insertNextCell(NORM_TRI3,3,sourceConn[0:3]);
        sourceMesh.insertNextCell(NORM_TRI3,3,sourceConn[3:6]);
        sourceMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(sourceCoords,4,2);
        sourceMesh.setCoords(myCoords);
        return sourceMesh;

    def build2DTargetMesh_1(self):
        targetCoords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        targetConn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.setMeshDimension(2);
        targetMesh.allocateCells(5);
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[0:4]);
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[4:7]);
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[7:10]);
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[10:14]);
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[14:18]);
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,9,2);
        targetMesh.setCoords(myCoords);
        return targetMesh;

    def build2DTargetMesh_3(self):
        targetCoords=[-0.6,-0.4, -0.1,-0.4, 1.1,-0.4, 2.1,-0.4, -0.6,0.1,  -0.1,0.1,  1.1,0.1,  2.1,0.1, -0.6,1.1,  -0.1,1.1]
        targetConn=[0,4,5,1, 1,5,6,2, 2,6,7,3, 4,8,9,5]
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.setMeshDimension(2);
        targetMesh.allocateCells(4);
        for i in range(4):
            targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[4*i:4*(i+1)])
            pass
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,10,2);
        targetMesh.setCoords(myCoords);
        return targetMesh;
        pass

    def setUp(self):
        pass
    pass

if __name__ == "__main__":
  unittest.main()
