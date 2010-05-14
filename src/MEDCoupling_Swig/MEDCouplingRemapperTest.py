#  -*- coding: iso-8859-1 -*-
#  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

from libMEDCouplingRemapper_Swig import *
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
    
    def setUp(self):
        pass
    pass

unittest.main()
