#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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
