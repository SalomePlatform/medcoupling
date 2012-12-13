#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2012  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License.
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
        self.failUnless(remapper.prepare(sourceMesh,targetMesh,"P0P0")==1);
        srcField=MEDCouplingFieldDouble.New(ON_CELLS);
        srcField.setNature(ConservativeVolumic);
        srcField.setMesh(sourceMesh);
        array=DataArrayDouble.New();
        ptr=sourceMesh.getNumberOfCells()*[None]
        for i in xrange(sourceMesh.getNumberOfCells()):
            ptr[i]=float(i+7)
            pass
        array.setValues(ptr,sourceMesh.getNumberOfCells(),1);
        srcField.setArray(array);
        trgfield=remapper.transferField(srcField,4.57);
        values=trgfield.getArray().getValues();
        valuesExpected=[7.5 ,7. ,7.,8.,7.5];
        for i in xrange(targetMesh.getNumberOfCells()):
            self.failUnless(abs(values[i]-valuesExpected[i])<1e-12);
            pass
        self.failUnless(1==trgfield.getArray().getNumberOfComponents());
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
        srcField.setNature(ConservativeVolumic);
        srcField.setMesh(sourceMesh);
        array=DataArrayDouble.New();
        ptr=sourceMesh.getNumberOfCells()*[None]
        for i in xrange(sourceMesh.getNumberOfCells()):
            ptr[i]=float(i+7);
            pass
        array.setValues(ptr,sourceMesh.getNumberOfCells(),1);
        srcField.setArray(array);
        trgfield=remapper.transferField(srcField,4.220173);
        values=trgfield.getArray().getValues();
        valuesExpected=[7.75, 7.0625, 4.220173,8.0]
        self.assertEqual(4,trgfield.getArray().getNumberOfTuples());
        self.assertEqual(1,trgfield.getArray().getNumberOfComponents());
        for i0 in xrange(4):
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
        srcField.setNature(ConservativeVolumic);
        srcField.setMesh(sourceMesh);
        array=DataArrayDouble.New();
        ptr=sourceMesh.getNumberOfCells()*[None]
        for i in xrange(sourceMesh.getNumberOfCells()):
            ptr[i]=float(i+7);
            pass
        array.setValues(ptr,sourceMesh.getNumberOfCells(),1);
        srcField.setArray(array);
        trgfield=MEDCouplingFieldDouble.New(ON_CELLS);
        trgfield.setNature(ConservativeVolumic);
        trgfield.setMesh(targetMesh);
        array=DataArrayDouble.New();
        ptr=targetMesh.getNumberOfCells()*[None]
        for i in xrange(targetMesh.getNumberOfCells()):
            ptr[i]=4.220173;
            pass
        array.setValues(ptr,targetMesh.getNumberOfCells(),1);
        trgfield.setArray(array);
        remapper.partialTransfer(srcField,trgfield);
        values=trgfield.getArray().getValues();
        valuesExpected=[7.75, 7.0625, 4.220173,8.0]
        self.assertEqual(4,trgfield.getArray().getNumberOfTuples());
        self.assertEqual(1,trgfield.getArray().getNumberOfComponents());
        for i0 in xrange(4):
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
        fieldSrc.setNature(Integral) ;  fieldSrc.setArray(arrSrc)
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
        fieldSrc.setNature(Integral) ;  fieldSrc.setArray(arrSrc)
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
        src.checkCoherency2(1e-10)
        trg.checkCoherency()
        fieldSrc=MEDCouplingFieldDouble(ON_CELLS,NO_TIME) ; fieldSrc.setMesh(src) ; arrSrc=DataArrayDouble([10.,30.])
        fieldSrc.setNature(Integral) ;  fieldSrc.setArray(arrSrc)
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
        fieldSrc.setNature(Integral) ;  fieldSrc.setArray(arrTrg)
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
        fieldSrc.setNature(Integral) ;  fieldSrc.setArray(arrSrc)
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
        trg.checkCoherency2(1e-10)
        src.checkCoherency()
        fieldSrc=MEDCouplingFieldDouble(ON_CELLS,NO_TIME) ; fieldSrc.setMesh(src) ; arrSrc=DataArrayDouble(100) ; arrSrc.iota(7.7)
        fieldSrc.setNature(Integral) ;  fieldSrc.setArray(arrSrc)
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
        fieldSrc.setNature(Integral) ;  fieldSrc.setArray(arrTrg)
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
        fieldSrc.setNature(Integral) ;  fieldSrc.setArray(arrSrc)
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
        fieldSrc.setNature(Integral) ;  fieldSrc.setArray(arrSrc) ; fieldSrc.checkCoherency()
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
        for i in xrange(4):
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

unittest.main()
