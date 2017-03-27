#  -*- coding: utf-8 -*-
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

from MEDCoupling import *
import unittest
from math import pi,e,sqrt,cos,sin
from datetime import datetime
from MEDCouplingDataForTest import MEDCouplingDataForTest
import rlcompleter,readline # this line has to be here, to ensure a usability of MEDCoupling/MEDLoader. B4 removing it please notify to anthony.geay@cea.fr

class MEDCouplingBasicsTest2(unittest.TestCase):
    def testMinMaxFields1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_3();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setTime(2.3,5,6);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[7.,107.,207.,8.,108.,208.,9.,109.,209.,10.,110.,210.,11.,111.,211.,12.,112.,212.,13.,113.,213.,14.,114.,214.,15.,115.,215.,16.,116.,216.]
        array.setValues(arr1,mesh1.getNumberOfCells(),3);
        f1.setArray(array);
        f2=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f2.setTime(7.8,4,5);
        f2.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr2=[6.,108.,206.,9.,107.,209.,8.,110.,208.,11.,109.,211.,10.,112.,210.,13.,111.,213.,12.,114.,212.,15.,113.,215.,14.,116.,214.,17.,115.,217.]
        array.setValues(arr2,mesh1.getNumberOfCells(),3);
        f2.setArray(array);
        #
        f3=f1.max(f2);
        expected1=[7.,108.,207.,9.,108.,209.,9.,110.,209.,11.,110.,211.,11.,112.,211.,13.,112.,213.,13.,114.,213.,15.,114.,215.,15.,116.,215.,17.,116.,217.]
        for i in range(30):
            self.assertAlmostEqual(expected1[i],f3.getIJ(0,i),9);
            pass
        #
        f4=f1.min(f2);
        expected2=[6.,107.,206.,8.,107.,208.,8.,109.,208.,10.,109.,210.,10.,111.,210.,12.,111.,212.,12.,113.,212.,14.,113.,214.,14.,115.,214.,16.,115.,216.]
        for i in range(30):
            self.assertAlmostEqual(expected2[i],f4.getIJ(0,i),9);
            pass
        #
        pass

    def testApplyLin1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_3();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,LINEAR_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr=[7.,107.,8.,108.,9.,109.,10.,110.,11.,111.,12.,112.,13.,113.,14.,114.,15.,115.,16.,116.]
        array.setValues(arr,mesh1.getNumberOfCells(),2);
        f1.setArray(array);
        #
        f1.applyLin(2.,3.,0);
        expected1=[17.,107.,19.,108.,21.,109.,23.,110.,25.,111.,27.,112.,29.,113.,31.,114.,33.,115.,35.,116.]
        for i in range(20):
            self.assertAlmostEqual(expected1[i],f1.getIJ(0,i),9);
            pass
        #
        arr2=[2.,102.,3.,103.,4.,104.,5.,105.,6.,106.,7.,107.,8.,108.,9.,109.,10.,110.,11.,111.]
        array=DataArrayDouble.New();
        array.setValues(arr2,mesh1.getNumberOfCells(),2);
        f1.setEndArray(array);
        #
        f1.applyLin(4.,5.,1);
        #
        expected2=[17.,433.,19.,437.,21.,441.,23.,445.,25.,449.,27.,453.,29.,457.,31.,461.,33.,465.,35.,469.]
        for i in range(20):
            self.assertAlmostEqual(expected2[i],f1.getIJ(0,i),9);
            pass
        expected3=[2.,413.,3.,417.,4.,421.,5.,425.,6.,429.,7.,433.,8.,437.,9.,441.,10.,445.,11.,449.]
        for i in range(20):
            self.assertAlmostEqual(expected3[i],f1.getEndArray().getIJ(0,i),9);
            pass
        #
        pass

    def testGetIdsInRange1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_3();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setTime(2.3,5,6);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[2.,8.,6.,5.,11.,7.,9.,3.,10.,4.]
        array.setValues(arr1,mesh1.getNumberOfCells(),1);
        f1.setArray(array);
        #
        f1.checkConsistencyLight();
        da=f1.findIdsInRange(2.9,7.1);
        self.assertEqual(5,da.getNbOfElems());
        expected1=[2,3,5,7,9]
        self.assertEqual(expected1,list(da.getValues()));
        da=f1.findIdsInRange(8.,12.);
        self.assertEqual(4,da.getNbOfElems());
        expected2=[1,4,6,8]
        self.assertEqual(expected2,list(da.getValues()));
        #
        pass

    def testBuildSubPart1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setTime(2.3,5,6);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[3.,103.,4.,104.,5.,105.,6.,106.,7.,107.]
        array.setValues(arr1,mesh1.getNumberOfCells(),2);
        f1.setArray(array);
        #
        part1=[2,1,4]
        f2=f1[part1];
        f2.zipCoords()
        self.assertEqual(3,f2.getNumberOfTuples());
        self.assertEqual(2,f2.getNumberOfComponents());
        expected1=[5.,105.,4.,104.,7.,107.]
        for i in range(6):
            self.assertAlmostEqual(f2.getIJ(0,i),expected1[i],12);
            pass
        self.assertEqual(3,f2.getMesh().getNumberOfCells());
        self.assertEqual(6,f2.getMesh().getNumberOfNodes());
        self.assertEqual(2,f2.getMesh().getSpaceDimension());
        self.assertEqual(2,f2.getMesh().getMeshDimension());
        m2C=f2.getMesh();
        self.assertEqual(13,m2C.getNodalConnectivityArrayLen());
        expected2=[0.2, -0.3, 0.7, -0.3, 0.2, 0.2, 0.7, 0.2, 0.2, 0.7, 0.7, 0.7]
        for i in range(12):
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12);
            pass
        expected3=[3,2,3,1,3,0,2,1,4,4,5,3,2]
        self.assertEqual(expected3,list(m2C.getNodalConnectivity().getValues()));
        expected4=[0,4,8,13]
        self.assertEqual(expected4,list(m2C.getNodalConnectivityIndex().getValues()));
        # Test with field on nodes.
        f1=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME);
        f1.setTime(2.3,5,6);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr2=[3.,103.,4.,104.,5.,105.,6.,106.,7.,107.,8.,108.,9.,109.,10.,110.,11.,111.]
        array.setValues(arr2,mesh1.getNumberOfNodes(),2);
        f1.setArray(array);
        part2=[1,2]
        f2=f1.buildSubPart(part2);
        self.assertEqual(4,f2.getNumberOfTuples());
        self.assertEqual(2,f2.getNumberOfComponents());
        expected5=[4.,104.,5.,105.,7.,107.,8.,108.]
        for i in range(8):
            self.assertAlmostEqual(f2.getIJ(0,i),expected5[i],12);
            pass
        self.assertEqual(2,f2.getMesh().getNumberOfCells());
        self.assertEqual(4,f2.getMesh().getNumberOfNodes());
        self.assertEqual(2,f2.getMesh().getSpaceDimension());
        self.assertEqual(2,f2.getMesh().getMeshDimension());
        m2C=f2.getMesh();
        self.assertEqual(8,m2C.getNodalConnectivityArrayLen());
        for i in range(8):  # 8 is not an error
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12);
            pass
        self.assertEqual(expected3[:4],list(m2C.getNodalConnectivity().getValues())[4:]);
        self.assertEqual(expected3[4:8],list(m2C.getNodalConnectivity().getValues())[:4]);
        self.assertEqual(expected4[:3],list(m2C.getNodalConnectivityIndex().getValues()));
        #idem previous because nodes of cell#4 are not fully present in part3
        part3=[1,2]
        arrr=DataArrayInt.New();
        arrr.setValues(part3,2,1);
        f2=f1.buildSubPart(arrr);
        self.assertEqual(4,f2.getNumberOfTuples());
        self.assertEqual(2,f2.getNumberOfComponents());
        for i in range(8):
            self.assertAlmostEqual(f2.getIJ(0,i),expected5[i],12);
            pass
        self.assertEqual(2,f2.getMesh().getNumberOfCells());
        self.assertEqual(4,f2.getMesh().getNumberOfNodes());
        self.assertEqual(2,f2.getMesh().getSpaceDimension());
        self.assertEqual(2,f2.getMesh().getMeshDimension());
        m2C=f2.getMesh();
        self.assertEqual(8,m2C.getNodalConnectivityArrayLen());
        for i in range(8):  # 8 is not an error
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12);
            pass
        self.assertEqual(expected3[:4],list(m2C.getNodalConnectivity().getValues())[4:8]);
        self.assertEqual(expected3[4:8],list(m2C.getNodalConnectivity().getValues())[:4]);
        self.assertEqual(expected4[:3],list(m2C.getNodalConnectivityIndex().getValues()));
        #
        part4=[1,2,4]
        f2=f1.buildSubPart(part4);
        self.assertEqual(6,f2.getNumberOfTuples());
        self.assertEqual(2,f2.getNumberOfComponents());
        expected6=[4.,104.,5.,105.,7.,107.,8.,108.,10.,110.,11.,111.]
        for i in range(12):
            self.assertAlmostEqual(f2.getIJ(0,i),expected6[i],12);
            pass
        self.assertEqual(3,f2.getMesh().getNumberOfCells());
        self.assertEqual(6,f2.getMesh().getNumberOfNodes());
        self.assertEqual(2,f2.getMesh().getSpaceDimension());
        self.assertEqual(2,f2.getMesh().getMeshDimension());
        m2C=f2.getMesh();
        self.assertEqual(13,m2C.getNodalConnectivityArrayLen());
        for i in range(12):
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12);
            pass
        self.assertEqual(expected3[0:4],list(m2C.getNodalConnectivity().getValues())[4:8]);
        self.assertEqual(expected3[4:8],list(m2C.getNodalConnectivity().getValues())[0:4]);
        self.assertEqual(expected3[8:13],list(m2C.getNodalConnectivity().getValues())[8:13]);
        self.assertEqual(expected4,list(m2C.getNodalConnectivityIndex().getValues()));
        pass

    def testDoublyContractedProduct1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5]
        array.setValues(arr1,mesh1.getNumberOfCells(),6);
        f1.setArray(array);
        f1.checkConsistencyLight();
        #
        f2=f1.doublyContractedProduct();
        f2.checkConsistencyLight();
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        for i in range(5):
            self.assertAlmostEqual(3906.56,f2.getIJ(i,0),9);
            pass
        #
        pass

    def testDeterminant1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,CONST_ON_TIME_INTERVAL);
        f1.setTime(2.3,5,6);
        f1.setEndTime(3.8,7,3);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5]
        array.setValues(arr1,mesh1.getNumberOfCells(),4);
        f1.setArray(array);
        #4 components
        f1.checkConsistencyLight();
        f2=f1.determinant();
        f2.checkConsistencyLight();
        self.assertEqual(CONST_ON_TIME_INTERVAL,f2.getTimeDiscretization());
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfValues());
        for i in range(5):
            self.assertAlmostEqual(-2.42,f2.getIJ(i,0),13);
            pass
        #6 components multi arrays with end array not defined
        f1=MEDCouplingFieldDouble.New(ON_NODES,LINEAR_TIME);
        f1.setTime(2.3,5,6);
        f1.setEndTime(3.8,7,3);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr2=[1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7,
              1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7]
        array.setValues(arr2,mesh1.getNumberOfNodes(),6);
        f1.setArray(array);
        self.assertRaises(InterpKernelException,f1.checkConsistencyLight);#no end array specified !
        #
        f2=f1.determinant();
        self.assertEqual(LINEAR_TIME,f2.getTimeDiscretization());
        self.assertEqual(1,f2.getArray().getNumberOfComponents());
        self.assertEqual(9,f2.getNumberOfTuples());
        for i in range(9):
            self.assertAlmostEqual(137.335,f2.getIJ(i,0),10);
            pass
        #6 components multi arrays with end array defined
        array=DataArrayDouble.New();
        arr3=[7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5,
              7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5]
        array.setValues(arr3,mesh1.getNumberOfNodes(),6);
        f1.setEndArray(array);
        f1.checkConsistencyLight();
        f2=f1.determinant();
        f2.checkConsistencyLight();
        self.assertEqual(LINEAR_TIME,f2.getTimeDiscretization());
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(9,f2.getNumberOfTuples());
        time2,it,order=f2.getTime()
        self.assertAlmostEqual(2.3,time2,12);
        self.assertEqual(5,it);
        self.assertEqual(6,order);
        time2,it,order=f2.getEndTime()
        self.assertAlmostEqual(3.8,time2,12);
        self.assertEqual(7,it);
        self.assertEqual(3,order);
        for i in range(9):
            self.assertAlmostEqual(137.335,f2.getIJ(i,0),10);
            self.assertAlmostEqual(1289.685,f2.getEndArray().getIJ(i,0),9);
            pass
        #9 components
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setTime(7.8,10,2);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr4=[1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1]
        array.setValues(arr4,mesh1.getNumberOfCells(),9);
        f1.setArray(array);
        #
        f1.checkConsistencyLight();
        f2=f1.determinant();
        f2.checkConsistencyLight();
        self.assertEqual(ONE_TIME,f2.getTimeDiscretization());
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        time2,it,order=f2.getTime()
        self.assertAlmostEqual(7.8,time2,12);
        self.assertEqual(10,it);
        self.assertEqual(2,order);
        for i in range(5):
            self.assertAlmostEqual(3.267,f2.getIJ(i,0),13);
            pass
        pass

    def testEigenValues1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7]
        array.setValues(arr1,mesh1.getNumberOfCells(),6);
        f1.setArray(array);
        f1.checkConsistencyLight();
        #
        f2=f1.eigenValues();
        f2.checkConsistencyLight();
        self.assertEqual(3,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        expected1=[13.638813677891717,-4.502313844635971,-2.2364998332557486]
        for i in range(5):
            self.assertAlmostEqual(expected1[0],f2.getIJ(i,0),13);
            self.assertAlmostEqual(expected1[1],f2.getIJ(i,1),13);
            self.assertAlmostEqual(expected1[2],f2.getIJ(i,2),13);
            pass
        pass

    def testEigenVectors1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7]
        array.setValues(arr1,mesh1.getNumberOfCells(),6);
        f1.setArray(array);
        f1.checkConsistencyLight();
        #
        f2=f1.eigenVectors();
        f2.checkConsistencyLight();
        self.assertEqual(9,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        expected1=[0.5424262364180696, 0.5351201064614425, 0.6476266283176001,#eigenvect 0
                   0.7381111277307373, 0.06458838384003074, -0.6715804522117897,#eigenvect 1
                   -0.4012053603397987, 0.8423032781211455, -0.3599436712889738#eigenvect 2
                   ]
        for i in range(5):
            self.assertAlmostEqual(expected1[0],f2.getIJ(i,0),13);
            self.assertAlmostEqual(expected1[1],f2.getIJ(i,1),13);
            self.assertAlmostEqual(expected1[2],f2.getIJ(i,2),13);
            self.assertAlmostEqual(expected1[3],f2.getIJ(i,3),13);
            self.assertAlmostEqual(expected1[4],f2.getIJ(i,4),13);
            self.assertAlmostEqual(expected1[5],f2.getIJ(i,5),13);
            self.assertAlmostEqual(expected1[6],f2.getIJ(i,6),13);
            self.assertAlmostEqual(expected1[7],f2.getIJ(i,7),13);
            self.assertAlmostEqual(expected1[8],f2.getIJ(i,8),13);
            pass
        #
        pass

    def testInverse1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1]
        array.setValues(arr1,mesh1.getNumberOfCells(),9);
        f1.setArray(array);
        f1.checkConsistencyLight();
        #
        f2=f1.inverse();
        f2.checkConsistencyLight();
        self.assertEqual(9,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        expected1=[-2.6538108356290113, 2.855831037649208, -1.1111111111111067, 3.461891643709813, -4.775022956841121, 2.2222222222222143, -1.1111111111111054, 2.222222222222214, -1.1111111111111072]
        for i in range(5):
            self.assertAlmostEqual(expected1[0],f2.getIJ(i,0),13);
            self.assertAlmostEqual(expected1[1],f2.getIJ(i,1),13);
            self.assertAlmostEqual(expected1[2],f2.getIJ(i,2),13);
            self.assertAlmostEqual(expected1[3],f2.getIJ(i,3),13);
            self.assertAlmostEqual(expected1[4],f2.getIJ(i,4),13);
            self.assertAlmostEqual(expected1[5],f2.getIJ(i,5),13);
            self.assertAlmostEqual(expected1[6],f2.getIJ(i,6),13);
            self.assertAlmostEqual(expected1[7],f2.getIJ(i,7),13);
            self.assertAlmostEqual(expected1[8],f2.getIJ(i,8),13);
            pass
        #
        array=DataArrayDouble.New();
        arr3=[7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5]
        array.setValues(arr3,mesh1.getNumberOfCells(),6);
        f1.setArray(array);
        f1.checkConsistencyLight();
        #
        f2=f1.inverse();
        f2.checkConsistencyLight();
        self.assertEqual(6,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        expected3=[-0.3617705098531818, -0.8678630828458127, -0.026843764174972983, 0.5539957431465833, 0.13133439560823013, -0.05301294502145887]
        for i in range(5):
            self.assertAlmostEqual(expected3[0],f2.getIJ(i,0),13);
            self.assertAlmostEqual(expected3[1],f2.getIJ(i,1),13);
            self.assertAlmostEqual(expected3[2],f2.getIJ(i,2),13);
            self.assertAlmostEqual(expected3[3],f2.getIJ(i,3),13);
            self.assertAlmostEqual(expected3[4],f2.getIJ(i,4),13);
            self.assertAlmostEqual(expected3[5],f2.getIJ(i,5),13);
            pass
        #
        array=DataArrayDouble.New();
        arr2=[1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5]
        array.setValues(arr2,mesh1.getNumberOfCells(),4);
        f1.setArray(array);
        f1.checkConsistencyLight();
        #
        f2=f1.inverse();
        f2.checkConsistencyLight();
        self.assertEqual(4,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        expected2=[-1.8595041322314059, 0.9504132231404963, 1.404958677685951, -0.49586776859504156]
        for i in range(5):
            self.assertAlmostEqual(expected2[0],f2.getIJ(i,0),13);
            self.assertAlmostEqual(expected2[1],f2.getIJ(i,1),13);
            self.assertAlmostEqual(expected2[2],f2.getIJ(i,2),13);
            self.assertAlmostEqual(expected2[3],f2.getIJ(i,3),13);
            pass
        #
        pass

    def testTrace1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1]
        array.setValues(arr1,mesh1.getNumberOfCells(),9);
        f1.setArray(array);
        f1.checkConsistencyLight();
        #
        f2=f1.trace();
        f2.checkConsistencyLight();
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        for i in range(5):
            self.assertAlmostEqual(15.9,f2.getIJ(i,0),13);
            pass
        #
        array=DataArrayDouble.New();
        arr3=[7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5]
        array.setValues(arr3,mesh1.getNumberOfCells(),6);
        f1.setArray(array);
        f1.checkConsistencyLight();
        #
        f2=f1.trace();
        f2.checkConsistencyLight();
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        for i in range(5):
            self.assertAlmostEqual(25.8,f2.getIJ(i,0),13);
            pass
        #
        array=DataArrayDouble.New();
        arr2=[1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5]
        array.setValues(arr2,mesh1.getNumberOfCells(),4);
        f1.setArray(array);
        f1.checkConsistencyLight();
        #
        f2=f1.trace();
        f2.checkConsistencyLight();
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        for i in range(5):
            self.assertAlmostEqual(5.7,f2.getIJ(i,0),13);
            pass
        #
        pass

    def testDeviator1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7]
        array.setValues(arr1,mesh1.getNumberOfCells(),6);
        f1.setArray(array);
        f1.checkConsistencyLight();
        #
        f2=f1.deviator();
        f2.checkConsistencyLight();
        self.assertEqual(6,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        expected1=[-1.1,0.,1.1,4.5,5.6,6.7]
        for i in range(5):
            self.assertAlmostEqual(expected1[0],f2.getIJ(i,0),13);
            self.assertAlmostEqual(expected1[1],f2.getIJ(i,1),13);
            self.assertAlmostEqual(expected1[2],f2.getIJ(i,2),13);
            self.assertAlmostEqual(expected1[3],f2.getIJ(i,3),13);
            self.assertAlmostEqual(expected1[4],f2.getIJ(i,4),13);
            self.assertAlmostEqual(expected1[5],f2.getIJ(i,5),13);
            pass
        #
        pass

    def testMagnitude1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6, 1.2,2.3,3.4,4.5,5.6, 1.2,2.3,3.4,4.5,5.6, 1.2,2.3,3.4,4.5,5.6, 1.2,2.3,3.4,4.5,5.6]
        array.setValues(arr1,mesh1.getNumberOfCells(),5);
        f1.setArray(array);
        f1.checkConsistencyLight();
        #
        f2=f1.magnitude();
        f2.checkConsistencyLight();
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        for i in range(5):
            self.assertAlmostEqual(8.3606219864313918,f2.getIJ(i,0),13);
            pass
        #
        pass

    def testMaxPerTuple1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6, 1.2,3.4,4.5,5.6,2.3, 3.4,4.5,5.6,1.2,2.3, 5.6,1.2,2.3,3.4,4.5, 4.5,5.6,1.2,2.3,3.4]
        array.setValues(arr1,mesh1.getNumberOfCells(),5);
        f1.setArray(array);
        f1.checkConsistencyLight();
        #
        f2=f1.maxPerTuple();
        f2.checkConsistencyLight();
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        for i in range(5):
            self.assertAlmostEqual(5.6,f2.getIJ(i,0),13);
            pass
        #
        d2,d2I=array.maxPerTupleWithCompoId()
        self.assertEqual(1,d2.getNumberOfComponents());
        self.assertEqual(5,d2.getNumberOfTuples());
        for i in range(5):
            self.assertAlmostEqual(5.6,d2.getIJ(i,0),13);
            pass
        self.assertTrue(d2I.isEqual(DataArrayInt([4,3,2,0,1])))
        pass

    def testChangeNbOfComponents(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6, 1.2,3.4,4.5,5.6,2.3, 3.4,4.5,5.6,1.2,2.3, 5.6,1.2,2.3,3.4,4.5, 4.5,5.6,1.2,2.3,3.4]
        array.setValues(arr1,mesh1.getNumberOfCells(),5);
        f1.setArray(array);
        f1.checkConsistencyLight();
        #
        f1.changeNbOfComponents(3,7.77);
        f1.checkConsistencyLight();
        self.assertEqual(3,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        expected1=[1.2,2.3,3.4, 1.2,3.4,4.5, 3.4,4.5,5.6, 5.6,1.2,2.3, 4.5,5.6,1.2]
        for i in range(15):
            self.assertAlmostEqual(expected1[i],f1.getIJ(0,i),13);
            pass
        f1.changeNbOfComponents(4,7.77);
        f1.checkConsistencyLight();
        self.assertEqual(4,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        expected2=[1.2,2.3,3.4,7.77, 1.2,3.4,4.5,7.77, 3.4,4.5,5.6,7.77, 5.6,1.2,2.3,7.77, 4.5,5.6,1.2,7.77]
        for i in range(20):
            self.assertAlmostEqual(expected2[i],f1.getIJ(0,i),13);
            pass
        #
        pass

    def testSortPerTuple1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6, 1.2,3.4,4.5,5.6,2.3, 3.4,4.5,5.6,1.2,2.3, 5.6,1.2,2.3,3.4,4.5, 4.5,5.6,1.2,2.3,3.4]
        array.setValues(arr1,mesh1.getNumberOfCells(),5);
        f1.setArray(array);
        f1.checkConsistencyLight();
        #
        f1.sortPerTuple(True);
        f1.checkConsistencyLight();
        self.assertEqual(5,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        for i in range(5):
            self.assertAlmostEqual(arr1[0],f1.getIJ(i,0),13);
            self.assertAlmostEqual(arr1[1],f1.getIJ(i,1),13);
            self.assertAlmostEqual(arr1[2],f1.getIJ(i,2),13);
            self.assertAlmostEqual(arr1[3],f1.getIJ(i,3),13);
            self.assertAlmostEqual(arr1[4],f1.getIJ(i,4),13);
            pass
        #
        f1.sortPerTuple(False);
        f1.checkConsistencyLight();
        self.assertEqual(5,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        for i in range(5):
            self.assertAlmostEqual(arr1[4],f1.getIJ(i,0),13);
            self.assertAlmostEqual(arr1[3],f1.getIJ(i,1),13);
            self.assertAlmostEqual(arr1[2],f1.getIJ(i,2),13);
            self.assertAlmostEqual(arr1[1],f1.getIJ(i,3),13);
            self.assertAlmostEqual(arr1[0],f1.getIJ(i,4),13);
            pass
        #
        pass

    def testIsEqualWithoutConsideringStr1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        mesh2=MEDCouplingDataForTest.build2DTargetMesh_1();
        #
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        mesh2.setName("rr");
        self.assertTrue(not mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        da1,da2=mesh1.checkGeoEquivalWith(mesh2,2,1e-12);
        self.assertRaises(InterpKernelException,mesh1.checkGeoEquivalWith,mesh2,0,1e-12);
        mesh2.setName("");
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        mesh2.getCoords().setInfoOnComponent(0,"tty");
        self.assertTrue(not mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        mesh2.getCoords().setInfoOnComponent(0,"");
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        mesh2.getCoords().setInfoOnComponent(1,"tty");
        self.assertTrue(not mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        mesh2.getCoords().setInfoOnComponent(1,"");
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        tmp=mesh2.getCoords().getIJ(0,3);
        mesh2.getCoords().setIJ(0,3,9999.);
        self.assertTrue(not mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(not mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        mesh2.getCoords().setIJ(0,3,tmp);
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        tmp2=mesh2.getNodalConnectivity().getIJ(0,4);
        mesh2.getNodalConnectivity().setIJ(0,4,0);
        self.assertTrue(not mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(not mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        mesh2.getNodalConnectivity().setIJ(0,4,tmp2);
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        #
        f1=mesh1.getMeasureField(True);
        f2=mesh2.getMeasureField(True);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        self.assertTrue(f1.isEqualWithoutConsideringStr(f2,1e-12,1e-12));
        f2.setName("ftest");
        self.assertTrue(not f1.isEqual(f2,1e-12,1e-12));
        self.assertTrue(f1.isEqualWithoutConsideringStr(f2,1e-12,1e-12));
        f1.setName("ftest");
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        self.assertTrue(f1.isEqualWithoutConsideringStr(f2,1e-12,1e-12));
        #
        f2.getArray().setInfoOnComponent(0,"eee");
        self.assertTrue(not f1.isEqual(f2,1e-12,1e-12));
        self.assertTrue(f1.isEqualWithoutConsideringStr(f2,1e-12,1e-12));
        f2.getArray().setInfoOnComponent(0,"");
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        self.assertTrue(f1.isEqualWithoutConsideringStr(f2,1e-12,1e-12));
        #
        f2.getArray().setIJ(1,0,0.123);
        self.assertTrue(not f1.isEqual(f2,1e-12,1e-12));
        self.assertTrue(not f1.isEqualWithoutConsideringStr(f2,1e-12,1e-12));
        f2.getArray().setIJ(1,0,0.125);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        self.assertTrue(f1.isEqualWithoutConsideringStr(f2,1e-12,1e-12));
        #
        pass
    
    def testGetNodeIdsOfCell1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        li=mesh1.getNodeIdsOfCell(1)
        expected1=[1, 4, 2]
        self.assertEqual(expected1,list(li))
        li=mesh1.getCoordinatesOfNode(4)
        self.assertEqual(2,len(li))
        self.assertAlmostEqual(0.2,li[0],13);
        self.assertAlmostEqual(0.2,li[1],13);
        li=mesh1.getCoords().getValuesAsTuple()
        self.assertEqual(9,len(li))
        li2=mesh1.getNodalConnectivityIndex().getValuesAsTuple()
        self.assertEqual(6,len(li2))
        pass

    def testGetEdgeRatioField1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m1.setTime(3.4,5,6); m1.setTimeUnit("us");
        f1=m1.getEdgeRatioField();
        self.assertAlmostEqual(3.4,f1.getTime()[0],12) ; self.assertEqual(5,f1.getTime()[1]) ; self.assertEqual(6,f1.getTime()[2])
        self.assertEqual("us",f1.getTimeUnit())
        self.assertEqual(m1.getNumberOfCells(),f1.getNumberOfTuples());
        self.assertEqual(5,f1.getNumberOfTuples());
        self.assertEqual(1,f1.getNumberOfComponents());
        expected1=[1.,1.4142135623730951, 1.4142135623730951,1.,1.]
        for i in range(5):
            self.assertAlmostEqual(expected1[i],f1.getIJ(i,0),14);
            pass
        #
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        f1=m1.getEdgeRatioField();
        self.assertEqual(m1.getNumberOfCells(),f1.getNumberOfTuples());
        self.assertEqual(5,f1.getNumberOfTuples());
        self.assertEqual(1,f1.getNumberOfComponents());
        expected2=[1.4142135623730951, 1.7320508075688772, 1.7320508075688772, 1.4142135623730951, 1.4142135623730951]
        for i in range(5):
            self.assertAlmostEqual(expected2[i],f1.getIJ(i,0),14);
            pass
        pass

    def testFillFromAnalytic3(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
        self.assertRaises(InterpKernelException,f1.fillFromAnalytic,1,"y+x");
        f1.setMesh(m)
        f1.setName("myField");
        f1.fillFromAnalytic(1,"y+x");
        f1.checkConsistencyLight();
        self.assertEqual(f1.getName(),"myField");
        self.assertEqual(f1.getTypeOfField(),ON_CELLS);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        values1=[-0.1,0.23333333333333336,0.56666666666666665,0.4,0.9]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values1),len(tmp))
        for i in range(len(values1)):
            self.assertTrue(abs(values1[i]-tmp[i])<1.e-12);
            pass
        #
        f1=MEDCouplingFieldDouble.New(ON_NODES,CONST_ON_TIME_INTERVAL)
        f1.setMesh(m)
        f1.fillFromAnalytic(1,"y+2*x");
        f1.setEndTime(1.2,3,4);
        f1.checkConsistencyLight();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),CONST_ON_TIME_INTERVAL);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        values2=[-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values2),len(tmp))
        for i in range(len(values2)):
            self.assertTrue(abs(values2[i]-tmp[i])<1.e-12);
            pass
        f1=MEDCouplingFieldDouble.New(ON_NODES,LINEAR_TIME);
        f1.setMesh(m)
        f1.fillFromAnalytic(1,"2.*x+y");
        f1.setEndTime(1.2,3,4);
        f1.checkConsistencyLight();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),LINEAR_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        tmp=f1.getArray().getValues();
        values2Bis=[-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1]
        self.assertEqual(len(values2Bis),len(tmp))
        for i in range(len(values2Bis)):
            self.assertTrue(abs(values2Bis[i]-tmp[i])<1.e-12);
            pass
        tmp=f1.getEndArray().getValues();
        self.assertEqual(len(values2Bis),len(tmp))
        for i in range(len(values2Bis)):
            self.assertTrue(abs(values2Bis[i]-tmp[i])<1.e-12);
            pass
        #
        f1=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME);
        f1.setMesh(m)
        f1.fillFromAnalytic(2,"(x+y)*IVec+2*(x+y)*JVec");
        f1.checkConsistencyLight();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(2,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        values3=[-0.6,-1.2,-0.1,-0.2,0.4,0.8,-0.1,-0.2,0.4,0.8,0.9,1.8,0.4,0.8,0.9,1.8,1.4,2.8]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values3),len(tmp))
        for i in range(len(values3)):
            self.assertTrue(abs(values3[i]-tmp[i])<1.e-12);
            pass
        values4=f1.accumulate();
        self.assertTrue(abs(3.6-values4[0])<1.e-12);
        self.assertTrue(abs(7.2-values4[1])<1.e-12);
        values4=f1.integral(True);
        self.assertTrue(abs(0.5-values4[0])<1.e-12);
        self.assertTrue(abs(1.-values4[1])<1.e-12);
        #
        f1=MEDCouplingFieldDouble.New(ON_NODES,NO_TIME);
        f1.setMesh(m);
        self.assertRaises(InterpKernelException,f1.fillFromAnalytic,1,"1./(x-0.2)");
        pass

    def testFieldDoubleOpEqual1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        self.assertRaises(InterpKernelException,f1.assign,0.07);
        f1.setMesh(m);
        f1.assign(0.07);
        f1.checkConsistencyLight();
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        for i in range(5):
            self.assertAlmostEqual(0.07,f1.getIJ(i,0),16);
            pass
        f1.assign(0.09);
        f1.checkConsistencyLight();
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        for i in range(5):
            self.assertAlmostEqual(0.09,f1.getIJ(i,0),16);
            pass
        #
        f1=MEDCouplingFieldDouble.New(ON_NODES,LINEAR_TIME);
        f1.setEndTime(4.5,2,3);
        f1.setMesh(m);
        f1.assign(0.08);
        f1.checkConsistencyLight();
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        for i in range(9):
            self.assertAlmostEqual(0.08,f1.getIJ(i,0),16);
            pass
        self.assertEqual(1,f1.getEndArray().getNumberOfComponents());
        self.assertEqual(9,f1.getEndArray().getNumberOfTuples());
        for i in range(9):
            self.assertAlmostEqual(0.08,f1.getEndArray().getIJ(i,0),16);
            pass
        pass

    def testAreaBary3D2(self):
        coordsForHexa8=[-75.45749305371, 180.95495078401, 39.515472018008,
                        -9.755591679144, 23.394927935279, 5.108794294848,
                        14.337630157832, 61.705351002702, 160.42422501908,
                        -27.273893776752, 167.567731083961, 192.830034145464,
                        99.857193154796,264.499264735586,-8.287335493412,
                        144.939882761126,156.38626563134,-31.896173894226,
                        161.34096835726,182.4654895809,73.832387065572,
                        132.680430393685,255.37973247196,96.15235602819];
        volHexa8=3258520.29637466;
        baryHexa8=[43.925705821778, 155.31893955289, 65.874418109644]
        
        coordsForPenta6=[-68.199829618726,178.938498373416,62.608505919588,
                         8.461744647847,76.653979804423,165.00018874933,
                         -27.273893776752,167.567731083961,192.830034145464,
                         106.586501038965,262.629609408327,13.124533008813,
                         155.465082847275,197.414118382622,78.408350795821,
                         132.680430393685,255.37973247196,96.15235602819];
        volPenta6=944849.868507338;
        baryPenta6=[39.631002313543,182.692711783428,106.98540473964]
        
        coordsForPyra5=[132.680430393685,255.37973247196,96.15235602819,
                        -27.273893776752,167.567731083961,192.830034145464,
                        8.461744647847,76.653979804423,165.00018874933,
                        155.465082847275,197.414118382622,78.408350795821,
                        -68.199829618726,178.938498373416,62.608505919588];
        volPyra5=756943.92980254;
        baryPyra5=[29.204294116618,172.540129749156,118.01035951483]
        mesh=MEDCouplingUMesh.New("Bary3D2",3);
        coo=DataArrayDouble.New();
        tmp=coordsForHexa8+coordsForPenta6+coordsForPyra5
        coo.setValues(tmp,19,3);
        mesh.setCoords(coo);
        #
        tmpConn=[0,1,2,3,4,5,6,7]
        mesh.allocateCells(3);
        self.assertRaises(InterpKernelException,mesh.insertNextCell,NORM_HEXA8,9,tmpConn[0:8])
        mesh.insertNextCell(NORM_HEXA8,tmpConn[0:8])
        mesh.insertNextCell(NORM_PENTA6,6,[i+8 for i in tmpConn])
        mesh.insertNextCell(NORM_PYRA5,5,[i+14 for i in tmpConn])
        mesh.finishInsertingCells();
        mesh.checkConsistencyLight();
        mesh.mergeNodes(1e-7)
        self.assertEqual(12,mesh.getNumberOfNodes());
        vols=mesh.getMeasureField(True);
        self.assertEqual(3,vols.getNumberOfTuples());
        self.assertEqual(1,vols.getNumberOfComponents());
        self.assertAlmostEqual(volHexa8,vols.getIJ(0,0),6);
        self.assertAlmostEqual(volPenta6,vols.getIJ(1,0),7);
        self.assertAlmostEqual(volPyra5,vols.getIJ(2,0),7);
        bary=mesh.computeCellCenterOfMass();
        self.assertEqual(3,bary.getNumberOfTuples());
        self.assertEqual(3,bary.getNumberOfComponents());
        self.assertAlmostEqual(baryHexa8[0],bary.getIJ(0,0),11);
        self.assertAlmostEqual(baryHexa8[1],bary.getIJ(0,1),11);
        self.assertAlmostEqual(baryHexa8[2],bary.getIJ(0,2),11);
        self.assertAlmostEqual(baryPenta6[0],bary.getIJ(1,0),11);
        self.assertAlmostEqual(baryPenta6[1],bary.getIJ(1,1),11);
        self.assertAlmostEqual(baryPenta6[2],bary.getIJ(1,2),11);
        self.assertAlmostEqual(baryPyra5[0],bary.getIJ(2,0),11);
        self.assertAlmostEqual(baryPyra5[1],bary.getIJ(2,1),11);
        self.assertAlmostEqual(baryPyra5[2],bary.getIJ(2,2),11);
        pass

    def testGetMeasureFieldCMesh1(self):
        m=MEDCouplingCMesh.New();
        da=DataArrayDouble.New();
        discX=[2.3,3.4,5.8,10.2]
        discY=[12.3,23.4,45.8]
        discZ=[-0.7,1.2,1.25,2.13,2.67]
        da.setValues(discX,4,1);
        m.setCoordsAt(0,da);
        m.checkConsistencyLight();
        self.assertEqual(4,m.getNumberOfNodes());
        self.assertEqual(3,m.getNumberOfCells());
        self.assertEqual(1,m.getSpaceDimension());
        f=m.getMeasureField(True);
        self.assertEqual(3,f.getNumberOfTuples());
        self.assertEqual(1,f.getNumberOfComponents());
        expected1=[1.1,2.4,4.4]
        for i in range(3):
            self.assertAlmostEqual(expected1[i],f.getIJ(i,0),12);
            pass
        coords=m.getCoordinatesAndOwner();
        self.assertEqual(4,coords.getNumberOfTuples());
        self.assertEqual(1,coords.getNumberOfComponents());
        for i in range(4):
            self.assertAlmostEqual(discX[i],coords.getIJ(i,0),12);
            pass
        coords=m.computeCellCenterOfMass();
        self.assertEqual(3,coords.getNumberOfTuples());
        self.assertEqual(1,coords.getNumberOfComponents());
        expected1_3=[2.85,4.6,8.]
        for i in range(3):
            self.assertAlmostEqual(expected1_3[i],coords.getIJ(i,0),12);
            pass
        #
        da=DataArrayDouble.New();
        da.setValues(discY,3,1);
        m.setCoordsAt(1,da);
        m.checkConsistencyLight();
        self.assertEqual(12,m.getNumberOfNodes());
        self.assertEqual(6,m.getNumberOfCells());
        self.assertEqual(2,m.getSpaceDimension());
        f=m.getMeasureField(True);
        self.assertEqual(6,f.getNumberOfTuples());
        self.assertEqual(1,f.getNumberOfComponents());
        expected2=[12.21,26.64,48.84,24.64,53.76,98.56]
        for i in range(6):
            self.assertAlmostEqual(expected2[i],f.getIJ(i,0),12);
            pass
        coords=m.getCoordinatesAndOwner();
        self.assertEqual(12,coords.getNumberOfTuples());
        self.assertEqual(2,coords.getNumberOfComponents());
        expected2_2=[2.3,12.3,3.4,12.3,5.8,12.3,10.2,12.3, 2.3,23.4,3.4,23.4,5.8,23.4,10.2,23.4, 2.3,45.8,3.4,45.8,5.8,45.8,10.2,45.8]
        for i in range(24):
            self.assertAlmostEqual(expected2_2[i],coords.getIJ(0,i),12);
            pass
        coords=m.computeCellCenterOfMass();
        self.assertEqual(6,coords.getNumberOfTuples());
        self.assertEqual(2,coords.getNumberOfComponents());
        expected2_3=[2.85,17.85,4.6,17.85,8.,17.85, 2.85,34.6,4.6,34.6,8.,34.6]
        for i in range(12):
            self.assertAlmostEqual(expected2_3[i],coords.getIJ(0,i),12);
            pass
        #
        da=DataArrayDouble.New();
        da.setValues(discZ,5,1);
        m.setCoordsAt(2,da);
        m.checkConsistencyLight();
        self.assertEqual(60,m.getNumberOfNodes());
        self.assertEqual(24,m.getNumberOfCells());
        self.assertEqual(3,m.getSpaceDimension());
        f=m.getMeasureField(True);
        self.assertEqual(24,f.getNumberOfTuples());
        self.assertEqual(1,f.getNumberOfComponents());
        expected3=[23.199, 50.616, 92.796, 46.816, 102.144, 187.264, 0.6105, 1.332, 2.442, 1.232, 2.688, 4.928, 10.7448, 23.4432, 42.9792, 21.6832, 47.3088, 86.7328, 6.5934, 14.3856, 26.3736, 13.3056, 29.0304, 53.2224]
        for i in range(24):
            self.assertAlmostEqual(expected3[i],f.getIJ(i,0),12);
            pass
        coords=m.getCoordinatesAndOwner();
        self.assertEqual(60,coords.getNumberOfTuples());
        self.assertEqual(3,coords.getNumberOfComponents());
        expected3_2=[
            2.3,12.3,-0.7, 3.4,12.3,-0.7, 5.8,12.3,-0.7, 10.2,12.3,-0.7, 2.3,23.4,-0.7, 3.4,23.4,-0.7, 5.8,23.4,-0.7, 10.2,23.4,-0.7, 2.3,45.8,-0.7, 3.4,45.8,-0.7, 5.8,45.8,-0.7, 10.2,45.8,-0.7,
            2.3,12.3,1.2, 3.4,12.3,1.2, 5.8,12.3,1.2, 10.2,12.3,1.2, 2.3,23.4,1.2, 3.4,23.4,1.2, 5.8,23.4,1.2, 10.2,23.4,1.2, 2.3,45.8,1.2, 3.4,45.8,1.2, 5.8,45.8,1.2, 10.2,45.8,1.2,
            2.3,12.3,1.25, 3.4,12.3,1.25, 5.8,12.3,1.25, 10.2,12.3,1.25, 2.3,23.4,1.25, 3.4,23.4,1.25, 5.8,23.4,1.25, 10.2,23.4,1.25, 2.3,45.8,1.25, 3.4,45.8,1.25, 5.8,45.8,1.25, 10.2,45.8,1.25,
            2.3,12.3,2.13, 3.4,12.3,2.13, 5.8,12.3,2.13, 10.2,12.3,2.13, 2.3,23.4,2.13, 3.4,23.4,2.13, 5.8,23.4,2.13, 10.2,23.4,2.13, 2.3,45.8,2.13, 3.4,45.8,2.13, 5.8,45.8,2.13, 10.2,45.8,2.13,
            2.3,12.3,2.67, 3.4,12.3,2.67, 5.8,12.3,2.67, 10.2,12.3,2.67, 2.3,23.4,2.67, 3.4,23.4,2.67, 5.8,23.4,2.67, 10.2,23.4,2.67, 2.3,45.8,2.67, 3.4,45.8,2.67, 5.8,45.8,2.67, 10.2,45.8,2.67];
        for i in range(180):
            self.assertAlmostEqual(expected3_2[i],coords.getIJ(0,i),12);
            pass
        coords=m.computeCellCenterOfMass();
        self.assertEqual(24,coords.getNumberOfTuples());
        self.assertEqual(3,coords.getNumberOfComponents());
        expected3_3=[
            2.85,17.85,0.25,4.6,17.85,0.25,8.,17.85,0.25, 2.85,34.6,0.25,4.6,34.6,0.25,8.,34.6,0.25,
            2.85,17.85,1.225,4.6,17.85,1.225,8.,17.85,1.225, 2.85,34.6,1.225,4.6,34.6,1.225,8.,34.6,1.225,
            2.85,17.85,1.69,4.6,17.85,1.69,8.,17.85,1.69, 2.85,34.6,1.69,4.6,34.6,1.69,8.,34.6,1.69,
            2.85,17.85,2.4,4.6,17.85,2.4,8.,17.85,2.4, 2.85,34.6,2.4,4.6,34.6,2.4,8.,34.6,2.4];
        for i in range(72):
            self.assertAlmostEqual(expected3_3[i],coords.getIJ(0,i),12);
            pass
        pass

    def testFieldDoubleZipCoords1(self):
        m=MEDCouplingDataForTest.build2DTargetMeshMergeNode_1();
        f=m.fillFromAnalytic(ON_NODES,2,"x*2.");
        f.getArray().setInfoOnComponent(0,"titi");
        f.getArray().setInfoOnComponent(1,"tutu");
        f.checkConsistencyLight();
        self.assertEqual(18,f.getNumberOfTuples());
        self.assertEqual(2,f.getNumberOfComponents());
        expected1=[-0.6, -0.6, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, 0.4, 0.4]
        for i in range(36):
            self.assertAlmostEqual(expected1[i],f.getIJ(0,i),12);
            pass
        self.assertTrue(f.zipCoords());
        f.checkConsistencyLight();
        expected2=[-0.6, -0.6, 1.4, 1.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, 0.4, 0.4]
        for i in range(30):
            self.assertAlmostEqual(expected2[i],f.getIJ(0,i),12);
            pass
        self.assertTrue(not f.zipCoords());
        f.checkConsistencyLight();
        for i in range(30):
            self.assertAlmostEqual(expected2[i],f.getIJ(0,i),12);
            pass
        self.assertTrue(f.getArray().getInfoOnComponent(0)=="titi");
        self.assertTrue(f.getArray().getInfoOnComponent(1)=="tutu");
        pass

    def testFieldDoubleZipConnectivity1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        cells1=[2,3,4]
        m3_1=m2.buildPartOfMySelf(cells1,True);
        m3=m3_1;
        m4=MEDCouplingDataForTest.build2DSourceMesh_1();
        m5=MEDCouplingUMesh.MergeUMeshes(m1,m3);
        m6=MEDCouplingUMesh.MergeUMeshes(m5,m4);
        #
        self.assertEqual(10,m6.getNumberOfCells());
        self.assertEqual(22,m6.getNumberOfNodes());
        arr,areNodesMerged,newNbOfNodes=m6.mergeNodes(1e-13);
        self.assertEqual(9,m6.getNumberOfNodes());
        f=m6.fillFromAnalytic(ON_CELLS,2,"x");
        f2=m6.fillFromAnalytic(ON_NODES,2,"x");
        self.assertEqual(10,f.getNumberOfTuples());
        self.assertEqual(2,f.getNumberOfComponents());
        expected1=[-0.05, -0.05, 0.3666666666666667, 0.3666666666666667, 0.53333333333333321, 0.53333333333333321,
                   -0.05, -0.05, 0.45, 0.45, 0.53333333333333321, 0.53333333333333321, -0.05, -0.05, 0.45, 0.45,
                   0.36666666666666659, 0.36666666666666659, 0.033333333333333326, 0.033333333333333326];
        for i in range(20):
            self.assertAlmostEqual(expected1[i],f.getIJ(0,i),12);
            pass
        f.getArray().setInfoOnComponent(0,"titi");
        f.getArray().setInfoOnComponent(1,"tutu");
        f.checkConsistencyLight();
        self.assertTrue(f.zipConnectivity(0));
        expected2=[-0.05, -0.05, 0.3666666666666667, 0.3666666666666667, 0.53333333333333321, 0.53333333333333321,
                   -0.05, -0.05, 0.45, 0.45, 0.36666666666666659, 0.36666666666666659, 0.033333333333333326, 0.033333333333333326];
        self.assertEqual(7,f.getNumberOfTuples());
        self.assertEqual(2,f.getNumberOfComponents());
        for i in range(14):
            self.assertAlmostEqual(expected2[i],f.getIJ(0,i),12);
            pass
        self.assertTrue(f.getArray().getInfoOnComponent(0)=="titi");
        self.assertTrue(f.getArray().getInfoOnComponent(1)=="tutu");
        self.assertTrue(not f.zipConnectivity(0));
        #
        expected3=[-0.3, -0.3, 0.2, 0.2, 0.7, 0.7, -0.3, -0.3, 0.2, 0.2, 0.7, 0.7,
                   -0.3, -0.3, 0.2, 0.2, 0.7, 0.7];
        self.assertEqual(9,f2.getNumberOfTuples());
        self.assertEqual(2,f2.getNumberOfComponents());
        for i in range(18):
            self.assertAlmostEqual(expected3[i],f2.getIJ(0,i),12);
            pass
        self.assertTrue(f2.zipConnectivity(0));
        self.assertEqual(9,f2.getNumberOfTuples());
        self.assertEqual(2,f2.getNumberOfComponents());
        for i in range(18):
            self.assertAlmostEqual(expected3[i],f2.getIJ(0,i),12);
            pass
        pass

    def testDaDoubleRenumber1(self):
        a=DataArrayDouble.New();
        arr1=[1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1]
        a.setValues(arr1,7,2);
        a.setInfoOnComponent(0,"toto");
        a.setInfoOnComponent(1,"tata");
        #
        arr2=[3,1,0,6,5,4,2]
        b=a.renumber(arr2);
        self.assertEqual(7,b.getNumberOfTuples());
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertTrue(b.getInfoOnComponent(0)=="toto");
        self.assertTrue(b.getInfoOnComponent(1)=="tata");
        expected1=[3.1, 13.1, 2.1, 12.1, 7.1, 17.1, 1.1, 11.1, 6.1, 16.1, 5.1, 15.1, 4.1, 14.1]
        for i in range(14):
            self.assertAlmostEqual(expected1[i],b.getIJ(0,i),14);
            pass
        #
        c=DataArrayInt.New();
        arr3=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        c.setValues(arr3,7,2);
        c.setInfoOnComponent(0,"toto");
        c.setInfoOnComponent(1,"tata");
        d=c.renumber(arr2);
        self.assertEqual(7,d.getNumberOfTuples());
        self.assertEqual(2,d.getNumberOfComponents());
        self.assertTrue(d.getInfoOnComponent(0)=="toto");
        self.assertTrue(d.getInfoOnComponent(1)=="tata");
        expected2=[3, 13, 2, 12, 7, 17, 1, 11, 6, 16, 5, 15, 4, 14]
        for i in range(14):
            self.assertEqual(expected2[i],d.getIJ(0,i));
            pass
        pass

    def testDaDoubleRenumberAndReduce1(self):
        a=DataArrayDouble.New();
        arr1=[1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1]
        a.setValues(arr1,7,2);
        a.setInfoOnComponent(0,"toto");
        a.setInfoOnComponent(1,"tata");
        #
        arr2=[2,-1,1,-1,0,4,3]
        b=a.renumberAndReduce(arr2,5);
        self.assertEqual(5,b.getNumberOfTuples());
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertTrue(b.getInfoOnComponent(0)=="toto");
        self.assertTrue(b.getInfoOnComponent(1)=="tata");
        expected1=[5.1,15.1,3.1,13.1,1.1,11.1,7.1,17.1,6.1,16.1]
        for i in range(10):
            self.assertAlmostEqual(expected1[i],b.getIJ(0,i),14);
            pass
        #
        c=DataArrayInt.New();
        arr3=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        c.setValues(arr3,7,2);
        c.setInfoOnComponent(0,"toto");
        c.setInfoOnComponent(1,"tata");
        d=c.renumberAndReduce(arr2,5);
        self.assertEqual(5,d.getNumberOfTuples());
        self.assertEqual(2,d.getNumberOfComponents());
        self.assertTrue(d.getInfoOnComponent(0)=="toto");
        self.assertTrue(d.getInfoOnComponent(1)=="tata");
        expected2=[5,15,3,13,1,11,7,17,6,16]
        for i in range(10):
            self.assertEqual(expected2[i],d.getIJ(0,i));
            pass
        pass

    def testDaDoubleRenumberInPlace1(self):
        a=DataArrayDouble.New();
        arr1=[1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1]
        a.setValues(arr1,7,2);
        #
        arr2=[3,1,0,6,5,4,2]
        a.renumberInPlace(arr2);
        self.assertEqual(7,a.getNumberOfTuples());
        self.assertEqual(2,a.getNumberOfComponents());
        expected1=[3.1, 13.1, 2.1, 12.1, 7.1, 17.1, 1.1, 11.1, 6.1, 16.1, 5.1, 15.1, 4.1, 14.1]
        for i in range(14):
            self.assertAlmostEqual(expected1[i],a.getIJ(0,i),14);
            pass
        #
        c=DataArrayInt.New();
        arr3=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        c.setValues(arr3,7,2);
        c.renumberInPlace(arr2);
        self.assertEqual(7,c.getNumberOfTuples());
        self.assertEqual(2,c.getNumberOfComponents());
        expected2=[3, 13, 2, 12, 7, 17, 1, 11, 6, 16, 5, 15, 4, 14]
        for i in range(14):
            self.assertEqual(expected2[i],c.getIJ(0,i));
            pass
        pass

    def testDaDoubleRenumberR1(self):
        a=DataArrayDouble.New();
        arr1=[1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1]
        a.setValues(arr1,7,2);
        a.setInfoOnComponent(0,"toto");
        a.setInfoOnComponent(1,"tata");
        #
        arr2=[3,1,0,6,5,4,2]
        b=a.renumberR(arr2);
        self.assertEqual(7,b.getNumberOfTuples());
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertTrue(b.getInfoOnComponent(0)=="toto");
        self.assertTrue(b.getInfoOnComponent(1)=="tata");
        expected1=[4.1, 14.1, 2.1, 12.1, 1.1, 11.1, 7.1, 17.1, 6.1, 16.1, 5.1, 15.1, 3.1, 13.1]
        for i in range(14):
            self.assertAlmostEqual(expected1[i],b.getIJ(0,i),14);
            pass
        #
        c=DataArrayInt.New();
        arr3=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        c.setValues(arr3,7,2);
        c.setInfoOnComponent(0,"toto");
        c.setInfoOnComponent(1,"tata");
        d=c.renumberR(arr2);
        self.assertEqual(7,d.getNumberOfTuples());
        self.assertEqual(2,d.getNumberOfComponents());
        self.assertTrue(d.getInfoOnComponent(0)=="toto");
        self.assertTrue(d.getInfoOnComponent(1)=="tata");
        expected2=[4, 14, 2, 12, 1, 11, 7, 17, 6, 16, 5, 15, 3, 13]
        for i in range(14):
            self.assertEqual(expected2[i],d.getIJ(0,i));
            pass
        pass

    def testDaDoubleRenumberInPlaceR1(self):
        a=DataArrayDouble.New();
        arr1=[1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1]
        a.setValues(arr1,7,2);
        #
        arr2=[3,1,0,6,5,4,2]
        a.renumberInPlaceR(arr2);
        self.assertEqual(7,a.getNumberOfTuples());
        self.assertEqual(2,a.getNumberOfComponents());
        expected1=[4.1, 14.1, 2.1, 12.1, 1.1, 11.1, 7.1, 17.1, 6.1, 16.1, 5.1, 15.1, 3.1, 13.1]
        for i in range(14):
            self.assertAlmostEqual(expected1[i],a.getIJ(0,i),14);
            pass
        #
        c=DataArrayInt.New();
        arr3=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        c.setValues(arr3,7,2);
        c.renumberInPlaceR(arr2);
        self.assertEqual(7,c.getNumberOfTuples());
        self.assertEqual(2,c.getNumberOfComponents());
        expected2=[4, 14, 2, 12, 1, 11, 7, 17, 6, 16, 5, 15, 3, 13]
        for i in range(14):
            self.assertEqual(expected2[i],c.getIJ(0,i));
            pass
        pass

    def testDaDoubleSelectByTupleId1(self):
        a=DataArrayDouble.New();
        arr1=[1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1]
        a.setValues(arr1,7,2);
        a.setInfoOnComponent(0,"toto");
        a.setInfoOnComponent(1,"tata");
        #
        arr2=[4,2,0,6,5]
        b=a.selectByTupleId(arr2);
        self.assertEqual(5,b.getNumberOfTuples());
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertTrue(b.getInfoOnComponent(0)=="toto");
        self.assertTrue(b.getInfoOnComponent(1)=="tata");
        expected1=[5.1,15.1,3.1,13.1,1.1,11.1,7.1,17.1,6.1,16.1]
        for i in range(10):
            self.assertAlmostEqual(expected1[i],b.getIJ(0,i),14);
            pass
        #
        c=DataArrayInt.New();
        arr3=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        c.setValues(arr3,7,2);
        c.setInfoOnComponent(0,"toto");
        c.setInfoOnComponent(1,"tata");
        d=c.selectByTupleId(arr2);
        self.assertEqual(5,d.getNumberOfTuples());
        self.assertEqual(2,d.getNumberOfComponents());
        self.assertTrue(d.getInfoOnComponent(0)=="toto");
        self.assertTrue(d.getInfoOnComponent(1)=="tata");
        expected2=[5,15,3,13,1,11,7,17,6,16]
        for i in range(10):
            self.assertEqual(expected2[i],d.getIJ(0,i));
            pass
        pass

    def testDaDoubleGetMinMaxValues1(self):
        a=DataArrayDouble.New();
        arr1=[2.34,4.56,-6.77,4.55,4.56,2.24,2.34,1.02,4.56]
        a.setValues(arr1,9,1);
        m,where=a.getMaxValue();
        self.assertEqual(1,where);
        self.assertAlmostEqual(4.56,m,12);
        m,ws=a.getMaxValue2();
        self.assertAlmostEqual(4.56,m,12);
        self.assertEqual(3,ws.getNumberOfTuples());
        self.assertEqual(1,ws.getNumberOfComponents());
        expected1=[1,4,8]
        for i in range(3):
            self.assertEqual(expected1[i],ws.getIJ(i,0));
            pass
        a=DataArrayDouble.New();
        arr2=[-2.34,-4.56,6.77,-4.55,-4.56,-2.24,-2.34,-1.02,-4.56]
        a.setValues(arr2,9,1);
        m,where=a.getMinValue();
        self.assertEqual(1,where);
        self.assertAlmostEqual(-4.56,m,12);
        m,ws=a.getMinValue2();
        self.assertAlmostEqual(-4.56,m,12);
        self.assertEqual(3,ws.getNumberOfTuples());
        self.assertEqual(1,ws.getNumberOfComponents());
        for i in range(3):
            self.assertEqual(expected1[i],ws.getIJ(i,0));
            pass
        pass

    def testFieldDoubleGetMinMaxValues2(self):
        m2,m1=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        self.assertEqual(18,m2.getNumberOfCells());
        arr1=[8.71,4.53,-12.41,8.71,-8.71,8.7099,4.55,8.71,5.55,6.77,-1e-200,4.55,8.7099,0.,1.23,0.,2.22,8.71]
        f=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        a=DataArrayDouble.New();
        a.setValues(arr1,18,1);
        f.setArray(a);
        f.setMesh(m2);
        #
        f.checkConsistencyLight();
        m=f.getMaxValue();
        self.assertAlmostEqual(8.71,m,12);
        m,ws=f.getMaxValue2();
        self.assertAlmostEqual(8.71,m,12);
        self.assertEqual(4,ws.getNumberOfTuples());
        self.assertEqual(1,ws.getNumberOfComponents());
        expected1=[0,3,7,17]
        for i in range(4):
            self.assertEqual(expected1[i],ws.getIJ(i,0));
            pass
        #
        arr2=[-8.71,-4.53,12.41,-8.71,8.71,-8.7099,-4.55,-8.71,-5.55,-6.77,1e-200,-4.55,-8.7099,0.,-1.23,0.,-2.22,-8.71]
        a.setValues(arr2,18,1);
        f.checkConsistencyLight();
        m=f.getMinValue();
        self.assertAlmostEqual(-8.71,m,12);
        m,ws=f.getMinValue2();
        self.assertAlmostEqual(-8.71,m,12);
        self.assertEqual(4,ws.getNumberOfTuples());
        self.assertEqual(1,ws.getNumberOfComponents());
        for i in range(4):
            self.assertEqual(expected1[i],ws.getIJ(i,0));
            pass
        pass

    def testBuildUnstructuredCMesh1(self):
        m=MEDCouplingCMesh.New();
        da=DataArrayDouble.New();
        discX=[2.3,3.4,5.8,10.2]
        discY=[12.3,23.4,45.8]
        discZ=[-0.7,1.2,1.25,2.13,2.67]
        da.setValues(discX,4,1);
        m.setCoordsAt(0,da);
        m.checkConsistencyLight();
        self.assertEqual(0,m.getCellContainingPoint([2.4],1e-12));
        self.assertEqual(1,m.getCellContainingPoint([3.7],1e-12));
        self.assertEqual(2,m.getCellContainingPoint([5.9],1e-12));
        self.assertEqual(-1,m.getCellContainingPoint([10.3],1e-12));
        self.assertEqual(-1,m.getCellContainingPoint([1.3],1e-12));
        #
        m2=m.buildUnstructured();
        m2.checkConsistencyLight();
        f1=m.getMeasureField(False);
        f2=m2.getMeasureField(False);
        self.assertTrue(isinstance(f1.getMesh(),MEDCouplingCMesh))
        self.assertEqual(f1.getNumberOfTuples(),3);
        self.assertEqual(f2.getNumberOfTuples(),3);
        self.assertEqual(1,m2.getMeshDimension());
        self.assertEqual(1,m2.getSpaceDimension());
        for i in range(3):
            self.assertAlmostEqual(f1.getIJ(i,0),f2.getIJ(i,0),10);
            pass
        da=DataArrayDouble.New();
        da.setValues(discY,3,1);
        m.setCoordsAt(1,da);
        #
        m2=m.buildUnstructured();
        m2.checkConsistencyLight();
        f1=m.getMeasureField(False);
        f2=m2.getMeasureField(False);
        self.assertEqual(f1.getNumberOfTuples(),6);
        self.assertEqual(f2.getNumberOfTuples(),6);
        self.assertEqual(2,m2.getMeshDimension());
        self.assertEqual(2,m2.getSpaceDimension());
        for i in range(6):
            self.assertAlmostEqual(f1.getIJ(i,0),f2.getIJ(i,0),10);
            pass
        #
        da=DataArrayDouble.New();
        da.setValues(discZ,5,1);
        m.setCoordsAt(2,da);
        m2=m.buildUnstructured();
        m2.checkConsistencyLight();
        f1=m.getMeasureField(False);
        f2=m2.getMeasureField(False);
        self.assertEqual(f1.getNumberOfTuples(),24);
        self.assertEqual(f2.getNumberOfTuples(),24);
        self.assertEqual(3,m2.getMeshDimension());
        self.assertEqual(3,m2.getSpaceDimension());
        for i in range(24):
            self.assertAlmostEqual(f1.getIJ(i,0),f2.getIJ(i,0),10);
            pass
        #
        pos1=[5.,30.,2.]
        self.assertEqual(16,m.getCellContainingPoint(pos1,1e-12));
        #
        elems=m2.getCellsInBoundingBox([3.5,6.,12.2,25.,0.,1.5],1e-7)
        self.assertEqual([1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17],elems.getValues())
        #
        pt=[2.4,12.7,-3.4]
        m.scale(pt,3.7);
        m3=m.buildUnstructured();
        m2.scale(pt,3.7);
        self.assertTrue(m3.isEqual(m2,1e-12));
        pass

    def testDataArrayIntInvertO2NNO21(self):
        arr1=[2,0,4,1,5,3]
        da=DataArrayInt.New();
        da.setValues(arr1,6,1);
        da2=da.invertArrayO2N2N2O(6);
        self.assertEqual(6,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        expected1=[1,3,0,5,2,4]
        for i in range(6):
            self.assertEqual(expected1[i],da2.getIJ(i,0));
            pass
        da3=da2.invertArrayN2O2O2N(6);
        for i in range(6):
            self.assertEqual(arr1[i],da3.getIJ(i,0));
            pass
        #
        arr2=[3,-1,5,4,-1,0,-1,1,2,-1]
        da=DataArrayInt.New();
        da.setValues(arr2,10,1);
        da2=da.invertArrayO2N2N2O(6);
        self.assertEqual(6,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        expected2=[5,7,8,0,3,2]
        for i in range(6):
            self.assertEqual(expected2[i],da2.getIJ(i,0));
            pass
        da3=da2.invertArrayN2O2O2N(10);
        for i in range(10):
            self.assertEqual(arr2[i],da3.getIJ(i,0));
            pass
        pass
    
    def testKeepSetSelectedComponent1(self):
        arr1=[1.,2.,3.,4., 11.,12.,13.,14., 21.,22.,23.,24., 31.,32.,33.,34., 41.,42.,43.,44.]
        a1=DataArrayDouble.New();
        a1.setValues(arr1,5,4);
        expp=[21.,22.,23.,24.]
        self.assertEqual(4,len(a1.getTuple(2)));
        for i in range(4):
            self.assertAlmostEqual(expp[i],a1.getTuple(2)[i],12)
            pass
        a1.setInfoOnComponent(0,"aaaa");
        a1.setInfoOnComponent(1,"bbbb");
        a1.setInfoOnComponent(2,"cccc");
        a1.setInfoOnComponent(3,"dddd");
        arr2V=[1,2,1,2,0,0]
        a2=a1.keepSelectedComponents(arr2V);
        self.assertEqual(6,a2.getNumberOfComponents());
        self.assertEqual(5,a2.getNumberOfTuples());
        self.assertTrue(a2.getInfoOnComponent(0)=="bbbb");
        self.assertTrue(a2.getInfoOnComponent(1)=="cccc");
        self.assertTrue(a2.getInfoOnComponent(2)=="bbbb");
        self.assertTrue(a2.getInfoOnComponent(3)=="cccc");
        self.assertTrue(a2.getInfoOnComponent(4)=="aaaa");
        self.assertTrue(a2.getInfoOnComponent(5)=="aaaa");
        expected1=[2.,3.,2.,3.,1.,1., 12.,13.,12.,13.,11.,11., 22.,23.,22.,23.,21.,21., 32.,33.,32.,33.,31.,31., 42.,43.,42.,43.,41.,41.]
        for i in range(30):
            self.assertAlmostEqual(expected1[i],a2.getIJ(0,i),14);
            pass
        a3=a1.convertToIntArr();
        self.assertEqual([21,22,23,24],a3.getTuple(2))
        a4=a3.keepSelectedComponents(arr2V);
        self.assertEqual(6,a4.getNumberOfComponents());
        self.assertEqual(5,a4.getNumberOfTuples());
        self.assertTrue(a4.getInfoOnComponent(0)=="bbbb");
        self.assertTrue(a4.getInfoOnComponent(1)=="cccc");
        self.assertTrue(a4.getInfoOnComponent(2)=="bbbb");
        self.assertTrue(a4.getInfoOnComponent(3)=="cccc");
        self.assertTrue(a4.getInfoOnComponent(4)=="aaaa");
        self.assertTrue(a4.getInfoOnComponent(5)=="aaaa");
        for i in range(30):
            self.assertEqual(int(expected1[i]),a4.getIJ(0,i));
            pass
        # setSelectedComponents
        arr3V=[3,2]
        a5=a1.keepSelectedComponents(arr3V);
        a5.setInfoOnComponent(0,"eeee");
        a5.setInfoOnComponent(1,"ffff");
        arr4V=[1,2]
        a2.setSelectedComponents(a5,arr4V);
        self.assertEqual(6,a2.getNumberOfComponents());
        self.assertEqual(5,a2.getNumberOfTuples());
        self.assertTrue(a2.getInfoOnComponent(0)=="bbbb");
        self.assertTrue(a2.getInfoOnComponent(1)=="eeee");
        self.assertTrue(a2.getInfoOnComponent(2)=="ffff");
        self.assertTrue(a2.getInfoOnComponent(3)=="cccc");
        self.assertTrue(a2.getInfoOnComponent(4)=="aaaa");
        self.assertTrue(a2.getInfoOnComponent(5)=="aaaa");
        expected2=[2.,4.,3.,3.,1.,1., 12.,14.,13.,13.,11.,11., 22.,24.,23.,23.,21.,21., 32.,34.,33.,33.,31.,31., 42.,44.,43.,43.,41.,41.]
        for i in range(30):
            self.assertAlmostEqual(expected2[i],a2.getIJ(0,i),14);
            pass
        a6=a5.convertToIntArr();
        a6.setInfoOnComponent(0,"eeee");
        a6.setInfoOnComponent(1,"ffff");
        a4.setSelectedComponents(a6,arr4V);
        self.assertEqual(6,a4.getNumberOfComponents());
        self.assertEqual(5,a4.getNumberOfTuples());
        self.assertTrue(a4.getInfoOnComponent(0)=="bbbb");
        self.assertTrue(a4.getInfoOnComponent(1)=="eeee");
        self.assertTrue(a4.getInfoOnComponent(2)=="ffff");
        self.assertTrue(a4.getInfoOnComponent(3)=="cccc");
        self.assertTrue(a4.getInfoOnComponent(4)=="aaaa");
        self.assertTrue(a4.getInfoOnComponent(5)=="aaaa");
        for i in range(30):
            self.assertEqual(int(expected2[i]),a4.getIJ(0,i));
            pass
        # test of throw
        arr5V=[2,3,6]
        arr6V=[2,7,5]
        arr7V=[2,1,4,6]
        self.assertRaises(InterpKernelException,a2.keepSelectedComponents,arr5V);
        self.assertRaises(InterpKernelException,a2.keepSelectedComponents,arr6V);
        self.assertRaises(InterpKernelException,a2.setSelectedComponents,a1,arr7V);
        arr7V=arr7V[0:3]
        self.assertRaises(InterpKernelException,a2.setSelectedComponents,a1,arr7V);
        #
        pass

    def testKeepSetSelectedComponent2(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        arr1=[1.,2.,3.,4., 11.,12.,13.,14., 21.,22.,23.,24., 31.,32.,33.,34., 41.,42.,43.,44.]
        a1=DataArrayDouble.New();
        a1.setValues(arr1,5,4);
        a1.setInfoOnComponent(0,"aaaa");
        a1.setInfoOnComponent(1,"bbbb");
        a1.setInfoOnComponent(2,"cccc");
        a1.setInfoOnComponent(3,"dddd");
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setTime(2.3,4,5);
        f1.setMesh(m1);
        f1.setName("f1");
        f1.setArray(a1);
        f1.checkConsistencyLight();
        #
        arr2V=[1,2,1,2,0,0]
        f2=f1.keepSelectedComponents(arr2V);
        self.assertTrue(f2.getTimeDiscretization()==ONE_TIME);
        t,dt,it=f2.getTime()
        self.assertAlmostEqual(2.3,t,13);
        self.assertEqual(4,dt);
        self.assertEqual(5,it);
        f2.checkConsistencyLight();
        self.assertEqual(6,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        self.assertTrue(f2.getArray().getInfoOnComponent(0)=="bbbb");
        self.assertTrue(f2.getArray().getInfoOnComponent(1)=="cccc");
        self.assertTrue(f2.getArray().getInfoOnComponent(2)=="bbbb");
        self.assertTrue(f2.getArray().getInfoOnComponent(3)=="cccc");
        self.assertTrue(f2.getArray().getInfoOnComponent(4)=="aaaa");
        self.assertTrue(f2.getArray().getInfoOnComponent(5)=="aaaa");
        expected1=[2.,3.,2.,3.,1.,1., 12.,13.,12.,13.,11.,11., 22.,23.,22.,23.,21.,21., 32.,33.,32.,33.,31.,31., 42.,43.,42.,43.,41.,41.]
        for i in range(30):
            self.assertAlmostEqual(expected1[i],f2.getIJ(0,i),14);
            pass
        #setSelectedComponents
        arr3V=[3,2]
        f5=f1.keepSelectedComponents(arr3V);
        f5.setTime(6.7,8,9);
        f5.getArray().setInfoOnComponent(0,"eeee");
        f5.getArray().setInfoOnComponent(1,"ffff");
        f5.checkConsistencyLight();
        arr4V=[1,2]
        f2.setSelectedComponents(f5,arr4V);
        self.assertEqual(6,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        f2.checkConsistencyLight();
        t,dt,it=f2.getTime()
        self.assertAlmostEqual(2.3,t,13);
        self.assertEqual(4,dt);
        self.assertEqual(5,it);
        self.assertTrue(f2.getArray().getInfoOnComponent(0)=="bbbb");
        self.assertTrue(f2.getArray().getInfoOnComponent(1)=="eeee");
        self.assertTrue(f2.getArray().getInfoOnComponent(2)=="ffff");
        self.assertTrue(f2.getArray().getInfoOnComponent(3)=="cccc");
        self.assertTrue(f2.getArray().getInfoOnComponent(4)=="aaaa");
        self.assertTrue(f2.getArray().getInfoOnComponent(5)=="aaaa");
        expected2=[2.,4.,3.,3.,1.,1., 12.,14.,13.,13.,11.,11., 22.,24.,23.,23.,21.,21., 32.,34.,33.,33.,31.,31., 42.,44.,43.,43.,41.,41.]
        for i in range(30):
            self.assertAlmostEqual(expected2[i],f2.getIJ(0,i),14);
            pass
        #
        pass
    
    def testElementaryDAThrowAndSpecialCases(self):
        da=DataArrayInt.New();
        self.assertRaises(InterpKernelException, da.checkAllocated);
        self.assertRaises(InterpKernelException, da.fillWithValue, 1);
        self.assertRaises(InterpKernelException, da.iota, 1);
        da.alloc(7,1);
        da.fillWithValue(11); #11,11,11,11...
        da.iota(10); #10,11,12,13...
        
        db=DataArrayInt.New();
        db.alloc(7,2);
        
        dbl2=DataArrayDouble.New();
        dbl2.alloc(7,2);
        self.assertRaises(InterpKernelException, dbl2.isUniform, 10., 1e-15);
        self.assertRaises(InterpKernelException, dbl2.sort);
        self.assertRaises(InterpKernelException, dbl2.iota, 10.);
        
        dbl=DataArrayDouble.New();
        #DataArrayDouble not allocated yet
        self.assertRaises(InterpKernelException, dbl.iota, 10.);
        self.assertRaises(InterpKernelException, dbl.isUniform, 10., 1e-15);
        self.assertRaises(InterpKernelException, dbl.sort);
        self.assertRaises(InterpKernelException, dbl.fromNoInterlace);
        self.assertRaises(InterpKernelException, dbl.toNoInterlace);
        
        dbl.alloc(7,1);
        dbl.iota(10.);
        self.assertTrue(not dbl.isUniform(10.,1e-15));
        dbl.sort();
        self.assertTrue(dbl.isMonotonic(True, .99));
        self.assertTrue(dbl.isMonotonic(True, -.99));
        self.assertTrue(not dbl.isMonotonic(True, 1.1));
        self.assertTrue(not dbl.isMonotonic(True, -1.1));
        dbl.reverse();
        self.assertTrue(dbl.isMonotonic(False, .99));
        self.assertTrue(not dbl.isMonotonic(False, 1.1));
        self.assertTrue(not dbl.isMonotonic(False, -1.1));
        
        dc=DataArrayInt.New();
        dc.alloc(14,1);
        
        dd=DataArrayDouble.New();
        self.assertRaises(InterpKernelException, dd.checkAllocated);
        self.assertRaises(InterpKernelException, dd.fillWithValue, 1.);
        self.assertRaises(InterpKernelException, dd.iota, 1.);
        self.assertTrue(not ((dd.repr().find("No data"))==-1));
        
        dd.alloc(0,1); #Allocated but nbOfElements==0!
        self.assertTrue(not ((dd.repr().find("Number of tuples : 0"))==-1));
        self.assertTrue(not ((dd.repr().find("Empty Data"))==-1));
        dd.fillWithValue(11); #?!...ok
        dd.iota(10); #?!...ok
        self.assertTrue(dd.isMonotonic(True, 1.));  #nothing is monotonic
        self.assertTrue(dd.isMonotonic(False, 1.));
        
        self.assertRaises(InterpKernelException, db.copyStringInfoFrom, da);
        self.assertRaises(InterpKernelException, db.copyStringInfoFrom, da);
        cIds=[2,2]
        self.assertRaises(InterpKernelException, da.copyPartOfStringInfoFrom, db, cIds);
        cIds[0]=1;
        cIds[0]=-1;
        self.assertRaises(InterpKernelException, da.copyPartOfStringInfoFrom, db, cIds);
        
        info=["infoOfOneComponent"]*2;
        self.assertRaises(InterpKernelException, da.setInfoOnComponents, info);
        self.assertRaises(InterpKernelException, da.setInfoOnComponent, 1, info[0]);
        db.setInfoOnComponents(info);
        
        self.assertRaises(InterpKernelException, da.getInfoOnComponent, -1);
        self.assertRaises(InterpKernelException, da.getInfoOnComponent, 2);
        self.assertTrue(db.getInfoOnComponent(1)==db.getInfoOnComponent(0));
        self.assertRaises(InterpKernelException, db.getVarOnComponent, -1);
        self.assertRaises(InterpKernelException, db.getVarOnComponent, 2);
        self.assertRaises(InterpKernelException, db.getUnitOnComponent, -1);
        self.assertRaises(InterpKernelException, db.getUnitOnComponent, 2);
        
        self.assertTrue(da.GetVarNameFromInfo("varname unit ")=="varname unit ");
        self.assertTrue(da.GetVarNameFromInfo("varname]unit[")=="varname]unit[");
        self.assertTrue(da.GetVarNameFromInfo("[unit]")=="");
        self.assertTrue(da.GetVarNameFromInfo("varname [unit]")=="varname");
        
        self.assertTrue(da.GetUnitFromInfo("varname unit ")=="");
        self.assertTrue(da.GetUnitFromInfo("varname]unit[")=="");
        self.assertTrue(da.GetUnitFromInfo("[unit]")=="unit");
        self.assertTrue(da.GetUnitFromInfo("varname [unit]")=="unit");
        
        self.assertRaises(InterpKernelException, da.checkNbOfTuplesAndComp, db, "theMessageInThrow");
        self.assertRaises(InterpKernelException, da.checkNbOfTuplesAndComp, dc, "theMessageInThrow");
        self.assertRaises(InterpKernelException, db.checkNbOfTuplesAndComp, dc, "theMessageInThrow");
        
        self.assertRaises(InterpKernelException, da.checkNbOfTuplesAndComp, 7, 2, "theMessageInThrow");
        da.checkNbOfTuplesAndComp(7,1,"theMessageInThrow");
        
        self.assertRaises(InterpKernelException, db.checkNbOfElems, 7*2+1, "theMessageInThrow");
        db.checkNbOfElems(7*2,"theMessageInThrow");
        
        self.assertRaises(InterpKernelException, db.GetNumberOfItemGivenBES, 10, 9, 1, "theMessageInThrow");
        self.assertRaises(InterpKernelException, db.GetNumberOfItemGivenBES, 0, 1, -1, "theMessageInThrow");
        self.assertEqual(10,db.GetNumberOfItemGivenBES(0,10,1,"theMessageInThrow"));
        self.assertEqual(5,db.GetNumberOfItemGivenBES(0,10,2,"theMessageInThrow"));
        self.assertEqual(6,db.GetNumberOfItemGivenBES(0,11,2,"theMessageInThrow"));
        
        self.assertTrue(not ((da.repr().find("Number of components : 1"))==-1));
        self.assertTrue(not ((dd.repr().find("Number of components : 1"))==-1));
        self.assertTrue(not ((dbl.repr().find("Number of components : 1"))==-1));
        
        self.assertTrue(not ((da.reprZip().find("Number of components : 1"))==-1));
        self.assertTrue(not ((dd.reprZip().find("Number of components : 1"))==-1));
        self.assertTrue(not ((dbl.reprZip().find("Number of components : 1"))==-1));
        
        self.assertRaises(InterpKernelException, dbl.selectByTupleIdSafeSlice, 0, 1, -1);
        self.assertRaises(InterpKernelException, dbl.subArray, -1, 1);
        self.assertRaises(InterpKernelException, dbl.subArray, 8, 1);
        self.assertRaises(InterpKernelException, dbl.subArray, 0, 8);
        self.assertRaises(InterpKernelException, dbl.meldWith, dd);
        
        self.assertRaises(InterpKernelException, dbl.setPartOfValuesAdv, dbl2, da); #dbl dbl2 not have the same number of components
        self.assertRaises(InterpKernelException, dbl.setPartOfValuesAdv, dd, da);  #da tuple selector DataArrayInt instance not have exactly 2 components
        
        dbl3=DataArrayDouble.New();
        dbl3.alloc(6,2);
        dbl3.fillWithValue(11.);
        #bad number of components
        self.assertRaises(InterpKernelException, dbl3.getMaxValue);
        self.assertRaises(InterpKernelException, dd.getMaxValue);
        self.assertRaises(InterpKernelException, dbl3.getMinValue);
        self.assertRaises(InterpKernelException, dd.getMinValue);
        self.assertRaises(InterpKernelException, dbl3.getAverageValue);
        self.assertRaises(InterpKernelException, dd.getAverageValue);
        self.assertRaises(InterpKernelException, dd.accumulate, 100);
        self.assertRaises(InterpKernelException, dbl.fromPolarToCart);
        self.assertRaises(InterpKernelException, dbl3.fromCylToCart);
        self.assertRaises(InterpKernelException, dbl3.fromSpherToCart);
        self.assertRaises(InterpKernelException, dbl3.doublyContractedProduct);
        self.assertRaises(InterpKernelException, dbl3.determinant);
        self.assertRaises(InterpKernelException, dbl3.eigenValues);
        self.assertRaises(InterpKernelException, dbl3.eigenVectors);
        self.assertRaises(InterpKernelException, dbl3.inverse);
        self.assertRaises(InterpKernelException, dbl3.trace);
        self.assertRaises(InterpKernelException, dbl3.deviator);
        
        dbl3.setIJ(5,1,12.);
        self.assertTrue(dbl3.getMaxValueInArray()==12.);
        self.assertTrue(dbl3.getMinValueInArray()==11.);
        
        db.fillWithValue(100); #bad Ids
        self.assertRaises(InterpKernelException, dbl3.setPartOfValuesAdv, dbl2, db);
        db.fillWithValue(-1); #bad Ids
        self.assertRaises(InterpKernelException, dbl3.setPartOfValuesAdv, dbl2, db);
        db.fillWithValue(6); #bad Ids for dbl3
        self.assertRaises(InterpKernelException, dbl3.setPartOfValuesAdv, dbl2, db);
        
        dbl3.checkNoNullValues();
        dbl3.setIJ(5,0,0.);
        self.assertRaises(InterpKernelException, dbl3.checkNoNullValues);
        self.assertRaises(InterpKernelException, dbl3.applyInv, 1.);  #div by zero
        self.assertRaises(InterpKernelException, dbl2.findIdsInRange, 1., 2.);
        a=[]
        self.assertRaises(InterpKernelException, DataArrayDouble_Aggregate, a);
        self.assertRaises(InterpKernelException, DataArrayDouble_Meld, a);
        
        a=[dbl2,dbl]; #Nb of components mismatch
        self.assertRaises(InterpKernelException, DataArrayDouble_Aggregate, a);
        
        self.assertRaises(InterpKernelException, DataArrayDouble_Dot, dbl2, dbl);
        
        self.assertRaises(InterpKernelException, DataArrayDouble_CrossProduct, dbl2, dbl); #Nb of components mismatch
        self.assertRaises(InterpKernelException, DataArrayDouble_CrossProduct, dbl2, dbl2); #Nb of components must be equal to 3
        dbl4=DataArrayDouble.New();
        dbl4.alloc(6,3);
        dbl5=DataArrayDouble.New();
        dbl5.alloc(7,3);
        self.assertRaises(InterpKernelException, DataArrayDouble_CrossProduct, dbl4, dbl5); #Nb of tuples mismatch
        
        a[0]=dbl4; #Nb of tuple mismatch
        a[1]=dbl5; #Nb of tuple mismatch
        self.assertRaises(InterpKernelException, DataArrayDouble_Meld, a);
        self.assertRaises(InterpKernelException, DataArrayDouble_Dot, dbl4, dbl5);
        pass

    def testDAIGetIdsEqual1(self):
        tab1=[5,-2,-4,-2,3,2,-2];
        da=DataArrayInt.New();
        da.setValues(tab1,7,1);
        da2=da.findIdsEqual(-2);
        self.assertEqual(3,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        expected1=[1,3,6];
        self.assertEqual(expected1,da2.getValues());
        pass

    def testDAIGetIdsEqualList1(self):
        tab1=[5,-2,-4,-2,3,2,-2];
        da=DataArrayInt.New();
        da.setValues(tab1,7,1);
        da2=da.findIdsEqualList([3,-2,0]);
        self.assertEqual(4,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        expected1=[1,3,4,6];
        self.assertEqual(expected1,da2.getValues());
        pass

    def testDAFromNoInterlace1(self):
        tab1=[1,11,21,31,41,2,12,22,32,42,3,13,23,33,43]
        da=DataArrayInt.New();
        da.setValues(tab1,5,3);
        da2=da.fromNoInterlace();
        expected1=[1,2,3,11,12,13,21,22,23,31,32,33,41,42,43]
        self.assertEqual(5,da2.getNumberOfTuples());
        self.assertEqual(3,da2.getNumberOfComponents());# it's not a bug. Avoid to have 1 million components !
        self.assertEqual(expected1,da2.getValues());
        da3=da.convertToDblArr();
        da4=da3.fromNoInterlace();
        self.assertEqual(5,da4.getNumberOfTuples());
        self.assertEqual(3,da4.getNumberOfComponents());# it's not a bug. Avoid to have 1 million components !
        for i in range(15):
            self.assertAlmostEqual(expected1[i],da4.getIJ(0,i),14);
            pass
        pass
    
    def testDAToNoInterlace1(self):
        tab1=[1,2,3,11,12,13,21,22,23,31,32,33,41,42,43]
        da=DataArrayInt.New();
        da.setValues(tab1,5,3);
        da2=da.toNoInterlace();
        expected1=[1,11,21,31,41,2,12,22,32,42,3,13,23,33,43]
        self.assertEqual(5,da2.getNumberOfTuples());
        self.assertEqual(3,da2.getNumberOfComponents());# it's not a bug. Avoid to have 1 million components !
        self.assertEqual(expected1,da2.getValues());
        da3=da.convertToDblArr();
        da4=da3.toNoInterlace();
        self.assertEqual(5,da4.getNumberOfTuples());
        self.assertEqual(3,da4.getNumberOfComponents());# it's not a bug. Avoid to have 1 million components !
        for i in range(15):
            self.assertAlmostEqual(expected1[i],da4.getIJ(0,i),14);
            pass
        pass
    
    def testDAIsUniform1(self):
        tab1=[1,1,1,1,1]
        da=DataArrayInt.New();
        da.setValues(tab1,5,1);
        self.assertTrue(da.isUniform(1));
        da.setIJ(2,0,2);
        self.assertTrue(not da.isUniform(1));
        da.setIJ(2,0,1);
        self.assertTrue(da.isUniform(1));
        da2=da.convertToDblArr();
        self.assertTrue(da2.isUniform(1.,1.e-12));
        da2.setIJ(1,0,1.+1.e-13);
        self.assertTrue(da2.isUniform(1.,1.e-12));
        da2.setIJ(1,0,1.+1.e-11);
        self.assertTrue(not da2.isUniform(1.,1.e-12));
        pass

    def testDAHasUniqueValues1(self):
        da=DataArrayInt([1,2,3,4,5])
        self.assertTrue(da.hasUniqueValues())
        da[1,0] = 5
        self.assertFalse(da.hasUniqueValues())
        da=DataArrayInt([])
        self.assertTrue(da.hasUniqueValues())
        da=DataArrayInt([(1,2), (2,3)]) # wrong num of compo
        self.assertRaises(InterpKernelException, da.hasUniqueValues)
        da=DataArrayInt()  # non allocated array
        self.assertRaises(InterpKernelException, da.hasUniqueValues)
        pass
    
    def testDADFromPolarToCart1(self):
        tab1=[2.,0.2,2.5,0.7]
        da=DataArrayDouble.New();
        da.setValues(tab1,2,2);
        da2=da.fromPolarToCart();
        expected1=[1.9601331556824833,0.39733866159012243, 1.9121054682112213,1.6105442180942275]
        for i in range(4):
            self.assertAlmostEqual(expected1[i],da2.getIJ(0,i),13);
            pass
        pass
    
    def testDADFromCylToCart1(self):
        tab1=[2.,0.2,4.,2.5,0.7,9.]
        da=DataArrayDouble.New();
        da.setValues(tab1,2,3);
        da2=da.fromCylToCart();
        expected1=[1.9601331556824833,0.39733866159012243,4., 1.9121054682112213,1.6105442180942275,9.]
        for i in range(6):
            self.assertAlmostEqual(expected1[i],da2.getIJ(0,i),13);
            pass
        pass
    
    def testDADFromSpherToCart1(self):
        tab1=[2.,0.2,0.3,2.5,0.7,0.8]
        da=DataArrayDouble.New();
        da.setValues(tab1,2,3);
        da2=da.fromSpherToCart();
        expected1=[0.37959212195737485,0.11742160338765303,1.9601331556824833, 1.1220769624465328,1.1553337045129035,1.9121054682112213]
        for i in range(6):
            self.assertAlmostEqual(expected1[i],da2.getIJ(0,i),13);
            pass
        pass

    def testUnPolyze1(self):
        elts=[0,1,2,3,4,5,6,7]
        eltsV=elts;
        mesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        mesh.convertToPolyTypes(eltsV);
        mesh.unPolyze();
        mesh2=MEDCouplingDataForTest.build3DTargetMesh_1();
        mesh.checkConsistencyLight();
        self.assertTrue(mesh.isEqual(mesh2,1e-12));
        mesh.convertToPolyTypes(eltsV);
        self.assertTrue(not mesh.isEqual(mesh2,1e-12));
        mesh.getNodalConnectivity().setIJ(0,6,10);
        mesh.getNodalConnectivity().setIJ(0,7,9);
        mesh.getNodalConnectivity().setIJ(0,8,12);
        mesh.getNodalConnectivity().setIJ(0,9,13);
        mesh.unPolyze();
        self.assertTrue(mesh.isEqual(mesh2,1e-12));
        mesh.convertToPolyTypes(eltsV);
        mesh.getNodalConnectivity().setIJ(0,6,12);
        mesh.getNodalConnectivity().setIJ(0,7,13);
        mesh.getNodalConnectivity().setIJ(0,8,10);
        mesh.getNodalConnectivity().setIJ(0,9,9);
        mesh.unPolyze();
        self.assertTrue(mesh.isEqual(mesh2,1e-12));
        mesh.convertToPolyTypes(eltsV);
        mesh.getNodalConnectivity().setIJ(0,6,12);
        mesh.getNodalConnectivity().setIJ(0,7,10);
        mesh.getNodalConnectivity().setIJ(0,8,13);
        mesh.getNodalConnectivity().setIJ(0,9,9);
        mesh.unPolyze();
        self.assertTrue(not mesh.isEqual(mesh2,1e-12));
        # Test for 2D mesh
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        mesh2=MEDCouplingDataForTest.build2DTargetMesh_1();
        eltsV=eltsV[:5];
        mesh.convertToPolyTypes(eltsV);
        self.assertTrue(not mesh.isEqual(mesh2,1e-12));
        mesh.unPolyze();
        self.assertTrue(mesh.isEqual(mesh2,1e-12));
        pass

    def testConvertDegeneratedCells1(self):
        mesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        conn=[0,1,3,3,9,10,12,12, 0,1,3,4,9,9,9,9, 1,1,1,1,10,12,9,10, 10,11,12,9,1,1,1,1]
        mesh.allocateCells(4);
        mesh.insertNextCell(NORM_HEXA8,8,conn[0:8])
        mesh.insertNextCell(NORM_HEXA8,8,conn[8:16])
        mesh.insertNextCell(NORM_HEXA8,8,conn[16:24])
        mesh.insertNextCell(NORM_HEXA8,8,conn[24:32])
        mesh.finishInsertingCells();
        mesh.checkConsistencyLight();
        self.assertEqual(4,mesh.getNumberOfCells());
        self.assertEqual(NORM_HEXA8,mesh.getTypeOfCell(0));
        self.assertEqual(NORM_HEXA8,mesh.getTypeOfCell(1));
        self.assertEqual(NORM_HEXA8,mesh.getTypeOfCell(2));
        self.assertEqual(NORM_HEXA8,mesh.getTypeOfCell(3));
        f1=mesh.getMeasureField(True);
        mesh.convertDegeneratedCells();
        mesh.checkConsistencyLight();
        f2=mesh.getMeasureField(True);
        self.assertEqual(4,mesh.getNumberOfCells());
        self.assertEqual(NORM_PENTA6,mesh.getTypeOfCell(0));
        self.assertEqual(NORM_PYRA5,mesh.getTypeOfCell(1));
        self.assertEqual(NORM_TETRA4,mesh.getTypeOfCell(2));
        self.assertEqual(NORM_PYRA5,mesh.getTypeOfCell(3));
        for i in range(4):
            self.assertAlmostEqual(f1.getArray().getIJ(0,i),f2.getArray().getIJ(0,i),5);
            pass
        pass

    def testGetNodeIdsNearPoints1(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        coords=mesh.getCoords();
        tmp=DataArrayDouble.New();
        vals=[0.2,0.2,0.1,0.2,0.2,0.2]
        tmp.setValues(vals,3,2);
        tmp2=DataArrayDouble.Aggregate(coords,tmp);
        mesh.setCoords(tmp2);
        pts=[0.2,0.2,0.1,0.3,-0.3,0.7]
        c=mesh.getNodeIdsNearPoint(pts[:2],1e-7);
        self.assertEqual([4,9,11],c.getValues());
        c,cI=mesh.getNodeIdsNearPoints(pts,3,1e-7);
        self.assertEqual([0,3,3,4],cI.getValues());
        self.assertEqual([4,9,11,6],c.getValues());
        c,cI=mesh.getNodeIdsNearPoints(pts,1e-7);
        self.assertEqual([0,3,3,4],cI.getValues());
        self.assertEqual([4,9,11,6],c.getValues());
        c,cI=mesh.getNodeIdsNearPoints(DataArrayDouble.New(pts,3,2),1e-7);
        self.assertEqual([0,3,3,4],cI.getValues());
        self.assertEqual([4,9,11,6],c.getValues());
        self.assertRaises(InterpKernelException,mesh.getNodeIdsNearPoints,DataArrayDouble.New(pts,2,3),1e-7);
        pass

    def testFieldCopyTinyAttrFrom1(self):
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setName("f1");
        f1.setTimeTolerance(1.e-5);
        f1.setDescription("f1Desc");
        f1.setTime(1.23,4,5);
        f2=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f2.setName("f2");
        f2.setDescription("f2Desc");
        f2.setTime(6.78,9,10);
        f2.setTimeTolerance(4.556e-12);
        #
        f1.copyTinyAttrFrom(f2);
        self.assertAlmostEqual(4.556e-12,f1.getTimeTolerance(),24);
        t,dt,it=f1.getTime()
        self.assertAlmostEqual(6.78,t,12);
        self.assertEqual(9,dt);
        self.assertEqual(10,it);
        self.assertTrue(f1.getName()=="f1");#name unchanged
        self.assertTrue(f1.getDescription()=="f1Desc");#description unchanged
        #
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setName("f1");
        f1.setTimeTolerance(1.e-5);
        f1.setDescription("f1Desc");
        f2=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f2.setName("f2");
        f2.setDescription("f2Desc");
        f2.setTimeTolerance(4.556e-12);
        #
        f1.copyTinyAttrFrom(f2);
        self.assertAlmostEqual(4.556e-12,f1.getTimeTolerance(),24);
        self.assertTrue(f1.getName()=="f1");#name unchanged
        self.assertTrue(f1.getDescription()=="f1Desc");#description unchanged
        #
        f1=MEDCouplingFieldDouble.New(ON_CELLS,CONST_ON_TIME_INTERVAL);
        f1.setName("f1");
        f1.setTimeTolerance(1.e-5);
        f1.setDescription("f1Desc");
        f1.setTime(1.23,4,5);
        f1.setEndTime(5.43,2,1);
        f2=MEDCouplingFieldDouble.New(ON_CELLS,CONST_ON_TIME_INTERVAL);
        f2.setName("f2");
        f2.setDescription("f2Desc");
        f2.setTimeTolerance(4.556e-12);
        f2.setTime(6.78,9,10);
        f2.setEndTime(10.98,7,6);
        #
        f1.copyTinyAttrFrom(f2);
        self.assertAlmostEqual(4.556e-12,f1.getTimeTolerance(),24);
        self.assertTrue(f1.getName()=="f1");#name unchanged
        self.assertTrue(f1.getDescription()=="f1Desc");#description unchanged
        t,dt,it=f1.getTime()
        self.assertAlmostEqual(6.78,t,12);
        self.assertEqual(9,dt);
        self.assertEqual(10,it);
        t,dt,it=f1.getEndTime()
        self.assertAlmostEqual(10.98,t,12);
        self.assertEqual(7,dt);
        self.assertEqual(6,it);
        #
        f1=MEDCouplingFieldDouble.New(ON_CELLS,LINEAR_TIME);
        f1.setName("f1");
        f1.setTimeTolerance(1.e-5);
        f1.setDescription("f1Desc");
        f1.setTime(1.23,4,5);
        f1.setEndTime(5.43,2,1);
        f2=MEDCouplingFieldDouble.New(ON_CELLS,LINEAR_TIME);
        f2.setName("f2");
        f2.setDescription("f2Desc");
        f2.setTimeTolerance(4.556e-12);
        f2.setTime(6.78,9,10);
        f2.setEndTime(10.98,7,6);
        #
        f1.copyTinyAttrFrom(f2);
        self.assertAlmostEqual(4.556e-12,f1.getTimeTolerance(),24);
        self.assertTrue(f1.getName()=="f1");#name unchanged
        self.assertTrue(f1.getDescription()=="f1Desc");#description unchanged
        t,dt,it=f1.getTime()
        self.assertAlmostEqual(6.78,t,12);
        self.assertEqual(9,dt);
        self.assertEqual(10,it);
        t,dt,it=f1.getEndTime()
        self.assertAlmostEqual(10.98,t,12);
        self.assertEqual(7,dt);
        self.assertEqual(6,it);
        pass

    def testExtrudedMesh5(self):
        coo1=[0.,1.,2.,3.5]
        a=DataArrayDouble.New();
        a.setValues(coo1,4,1);
        b=MEDCouplingCMesh.New();
        b.setCoordsAt(0,a);
        c=b.buildUnstructured();
        self.assertEqual(1,c.getSpaceDimension());
        c.changeSpaceDimension(2);
        #
        d=DataArrayDouble.New();
        d.alloc(13,1);
        d.iota();
        e=MEDCouplingCMesh.New();
        e.setCoordsAt(0,d);
        f=e.buildUnstructured();
        g=f.getCoords().applyFunc(2,"3.5*IVec+x/6*3.14159265359*JVec");
        self.assertRaises(InterpKernelException,f.getCoords().applyFunc,2,"3.5*IVec+x/6*3.14159265359*KVec"); # KVec refers to component #2 and there is only 2 components !
        h=g.fromPolarToCart();
        f.setCoords(h);
        i=c.buildExtrudedMesh(f,1);
        self.assertEqual(52,i.getNumberOfNodes());
        tmp,tmp2,tmp3=i.mergeNodes(1e-9);
        self.assertTrue(tmp2);
        self.assertEqual(37,tmp3);
        i.convertDegeneratedCells();
        i.checkConsistencyLight();
        self.assertEqual(36,i.getNumberOfCells());
        self.assertEqual(37,i.getNumberOfNodes());
        self.assertEqual(12,i.getNumberOfCellsWithType(NORM_TRI3));
        self.assertEqual(24,i.getNumberOfCellsWithType(NORM_QUAD4));
        expected1=[0.25,0.75,2.0625]
        j=i.getMeasureField(True);
        for ii in range(12):
            for k in range(3):
                self.assertAlmostEqual(expected1[k],j.getIJ(0,ii*3+k),10);
                pass
            pass
        expected2=[0.62200846792814113, 0.16666666666681595, 1.4513530918323276, 0.38888888888923495, 2.6293994326053212, 0.7045454545460802, 0.45534180126145435, 0.45534180126150181, 1.0624642029433926, 1.0624642029435025, 1.9248539780597826, 1.9248539780599816, 0.16666666666661334, 0.62200846792815856, 0.38888888888876294, 1.4513530918323678, 0.70454545454522521, 2.629399432605394, -0.16666666666674007, 0.62200846792812436, -0.38888888888906142, 1.4513530918322881, -0.70454545454576778, 2.6293994326052488, -0.45534180126154766, 0.45534180126140844, -1.0624642029436118, 1.0624642029432834, -1.9248539780601803, 1.9248539780595841, -0.62200846792817499, 0.1666666666665495, -1.451353091832408, 0.388888888888613, -2.6293994326054668, 0.70454545454495332, -0.62200846792810593, -0.16666666666680507, -1.451353091832247, -0.38888888888921297, -2.6293994326051746, -0.70454545454604123, -0.45534180126135926, -0.45534180126159562, -1.0624642029431723, -1.0624642029437235, -1.9248539780593836, -1.9248539780603811, -0.1666666666664828, -0.62200846792819242, -0.38888888888846079, -1.4513530918324489, -0.70454545454467987, -2.6293994326055397, 0.16666666666687083, -0.62200846792808862, 0.38888888888936374, -1.4513530918322073, 0.70454545454631357, -2.6293994326051022, 0.45534180126164348, -0.45534180126131207, 1.0624642029438327, -1.0624642029430627, 1.9248539780605791, -1.9248539780591853, 0.62200846792821063, -0.16666666666641802, 1.4513530918324888, -0.38888888888831086, 2.6293994326056125, -0.70454545454440853]
        m=i.computeCellCenterOfMass();
        for i in range(72):
            self.assertAlmostEqual(expected2[i],m.getIJ(0,i),10);
            pass
        #
        pass

    def testExtrudedMesh6(self):
        coo1=[0.,1.,2.,3.5]
        a=DataArrayDouble.New();
        a.setValues(coo1,4,1);
        b=MEDCouplingCMesh.New();
        b.setCoordsAt(0,a);
        c=b.buildUnstructured();
        self.assertEqual(1,c.getSpaceDimension());
        c.changeSpaceDimension(2);
        #
        d=DataArrayDouble.New();
        d.alloc(5);
        d.iota();
        e=MEDCouplingCMesh.New();
        e.setCoordsAt(0,d);
        f=e.buildUnstructured();
        d2=f.getCoords().applyFunc("x*x/2");
        f.setCoords(d2);
        f.changeSpaceDimension(2);
        #
        center=[0.,0.]
        f.rotate(center,None,pi/3);
        g=c.buildExtrudedMesh(f,0);
        g.checkConsistencyLight();
        expected1=[ 0.4330127018922193, 0.4330127018922193, 0.649519052838329, 1.2990381056766578, 1.299038105676658, 1.948557158514987, 2.1650635094610955, 2.1650635094610964, 3.2475952641916446, 3.031088913245533, 3.0310889132455352, 4.546633369868303 ]
        f1=g.getMeasureField(True);
        for i in range(12):
            self.assertAlmostEqual(expected1[i],f1.getIJ(0,i),12);
            pass
        expected2=[0.625, 0.21650635094610962, 1.625, 0.21650635094610959, 2.8750000000000004, 0.21650635094610965, 1.1250000000000002, 1.0825317547305482, 2.125, 1.0825317547305482, 3.3750000000000004, 1.0825317547305484, 2.125, 2.8145825622994254, 3.125, 2.8145825622994254, 4.375, 2.8145825622994254, 3.6250000000000009, 5.4126587736527414, 4.625, 5.4126587736527414, 5.875, 5.4126587736527414]
        f2=g.computeCellCenterOfMass();
        for i in range(24):
            self.assertAlmostEqual(expected2[i],f2.getIJ(0,i),12);
            pass
        pass

    def testExtrudedMesh7(self):
        coo1=[0.,1.,2.,3.5]
        a=DataArrayDouble.New();
        a.setValues(coo1,4,1);
        b=MEDCouplingCMesh.New();
        b.setCoordsAt(0,a);
        c=b.buildUnstructured();
        self.assertEqual(1,c.getSpaceDimension());
        c.changeSpaceDimension(2);
        #
        d=DataArrayDouble.New();
        d.alloc(13,1);
        d.iota();
        e=MEDCouplingCMesh.New();
        e.setCoordsAt(0,d);
        f=e.buildUnstructured();
        g=f.getCoords().applyFunc(2,"3.5*IVec+x/6*3.14159265359*JVec");
        h=g.fromPolarToCart();
        f.setCoords(h);
        i=c.buildExtrudedMesh(f,1);
        self.assertEqual(52,i.getNumberOfNodes());
        tmp,tmp2,tmp3=i.mergeNodes(1e-9);
        self.assertTrue(tmp2);
        self.assertEqual(37,tmp3);
        i.convertDegeneratedCells();
        vec1=[10.,0]
        i.translate(vec1);
        g2=h.applyFunc(3,"13.5/3.5*x*IVec+0*JVec+13.5/3.5*y*KVec");
        f.setCoords(g2);
        i.changeSpaceDimension(3);
        i3=i.buildExtrudedMesh(f,1);
        f2=i3.getMeasureField(True);
        tmp,tmp2,tmp3=i.mergeNodes(1e-9);
        self.assertTrue(tmp2);
        self.assertEqual(444,tmp3);
        expected1=[1.327751058489274, 4.2942574094314701, 13.024068164857139, 1.3069177251569044, 4.1484240761012954, 12.297505664866796, 1.270833333332571, 3.8958333333309674, 11.039062499993179, 1.2291666666659207, 3.6041666666644425, 9.585937499993932, 1.1930822748415895, 3.3515759238941376, 8.3274943351204556, 1.1722489415082769, 3.2057425905609289, 7.6009318351210622, 1.1722489415082862, 3.2057425905609884, 7.6009318351213713, 1.1930822748416161, 3.3515759238943001, 8.3274943351212727, 1.2291666666659564, 3.6041666666646734, 9.5859374999950777, 1.2708333333326081, 3.8958333333311868, 11.039062499994293, 1.3069177251569224, 4.1484240761014384, 12.297505664867627, 1.3277510584902354, 4.2942574094346071, 13.024068164866796]
        for ii in range(12):
            for jj in range(36):
                self.assertAlmostEqual(expected1[jj],f2.getIJ(0,ii*36+jj),9);
                pass
        #
        pass

    def testSimplexize1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m.convertToPolyTypes([3]);
        da=m.simplexize(0);
        self.assertEqual(7,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        expected2=[0,0,1,2,3,4,4]
        for i in range(7):
            self.assertEqual(expected2[i],da.getIJ(i,0));
            pass
        m.checkConsistencyLight();
        self.assertEqual(7,m.getNumberOfCells());
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(0));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(1));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(2));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(3));
        self.assertEqual(NORM_POLYGON,m.getTypeOfCell(4));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(5));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(6));
        expected1=[0.125,0.125,0.125,0.125,0.25,0.125,0.125]
        f=m.getMeasureField(False);
        for i in range(7):
            self.assertAlmostEqual(expected1[i]*sqrt(2.),f.getIJ(i,0),10);
            pass
        types=m.getAllGeoTypes();
        self.assertEqual([NORM_TRI3,NORM_POLYGON],types);
        #
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m.convertToPolyTypes([3]);
        da=m.simplexize(1);
        self.assertEqual(7,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        for i in range(7):
            self.assertEqual(expected2[i],da.getIJ(i,0));
            pass
        m.checkConsistencyLight();
        types=m.getAllGeoTypes();
        self.assertEqual([NORM_TRI3,NORM_POLYGON],types);
        self.assertEqual(7,m.getNumberOfCells());
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(0));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(1));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(2));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(3));
        self.assertEqual(NORM_POLYGON,m.getTypeOfCell(4));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(5));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(6));
        f=m.getMeasureField(False);
        for i in range(7):
            self.assertAlmostEqual(expected1[i]*sqrt(2.),f.getIJ(i,0),10);
            pass
        pass

    def testSimplexize2(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m.convertToPolyTypes([3]);
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setMesh(m);
        arr=DataArrayDouble.New();
        arr1=[10.,110.,20.,120.,30.,130.,40.,140.,50.,150.]
        arr.setValues(arr1,5,2);
        f1.setArray(arr);
        #
        f1.checkConsistencyLight();
        self.assertTrue(f1.simplexize(0));
        f1.checkConsistencyLight();
        expected1=[10.,110.,10.,110.,20.,120.,30.,130.,40.,140.,50.,150.,50.,150.]
        for i in range(14):
            self.assertAlmostEqual(expected1[i],f1.getIJ(0,i),10);
            pass
        self.assertTrue(not f1.simplexize(0));
        for i in range(14):
            self.assertAlmostEqual(expected1[i],f1.getIJ(0,i),10);
            pass
        #
        pass

    def testDAMeld1(self):
        da1=DataArrayDouble.New();
        da1.alloc(7,2);
        da2=DataArrayDouble.New();
        da2.alloc(7,1);
        #
        da1.fillWithValue(7.);
        da2.iota(0.);
        da3=da2.applyFunc(3,"10*x*IVec+100*x*JVec+1000*x*KVec");
        #
        da1.setInfoOnComponent(0,"c0da1");
        da1.setInfoOnComponent(1,"c1da1");
        da3.setInfoOnComponent(0,"c0da3");
        da3.setInfoOnComponent(1,"c1da3");
        da3.setInfoOnComponent(2,"c2da3");
        #
        da1C=da1.deepCopy();
        da1.meldWith(da3);
        self.assertEqual(5,da1.getNumberOfComponents());
        self.assertEqual(7,da1.getNumberOfTuples());
        self.assertTrue(da1.getInfoOnComponent(0)=="c0da1");
        self.assertTrue(da1.getInfoOnComponent(1)=="c1da1");
        self.assertTrue(da1.getInfoOnComponent(2)=="c0da3");
        self.assertTrue(da1.getInfoOnComponent(3)=="c1da3");
        self.assertTrue(da1.getInfoOnComponent(4)=="c2da3");
        #
        expected1=[7.,7.,0.,0.,0., 7.,7.,10.,100.,1000., 7.,7.,20.,200.,2000., 7.,7.,30.,300.,3000., 7.,7.,40.,400.,4000.,7.,7.,50.,500.,5000.,7.,7.,60.,600.,6000.]
        for i in range(35):
            self.assertAlmostEqual(expected1[i],da1.getIJ(0,i),10);
            pass
        #
        dai1=da1C.convertToIntArr();
        dai3=da3.convertToIntArr();
        dai1.meldWith(dai3);
        self.assertEqual(5,dai1.getNumberOfComponents());
        self.assertEqual(7,dai1.getNumberOfTuples());
        self.assertTrue(dai1.getInfoOnComponent(0)=="c0da1");
        self.assertTrue(dai1.getInfoOnComponent(1)=="c1da1");
        self.assertTrue(dai1.getInfoOnComponent(2)=="c0da3");
        self.assertTrue(dai1.getInfoOnComponent(3)=="c1da3");
        self.assertTrue(dai1.getInfoOnComponent(4)=="c2da3");
        for i in range(35):
            self.assertEqual(int(expected1[i]),dai1.getIJ(0,i));
            pass
        # test of static method DataArrayDouble::meld
        da4=DataArrayDouble.Meld(da1C,da3);
        tmp=DataArrayDouble.Meld([da1C,da3]);
        self.assertTrue(da4.isEqual(tmp,1e-10))
        self.assertEqual(5,da4.getNumberOfComponents());
        self.assertEqual(7,da4.getNumberOfTuples());
        self.assertTrue(da4.getInfoOnComponent(0)=="c0da1");
        self.assertTrue(da4.getInfoOnComponent(1)=="c1da1");
        self.assertTrue(da4.getInfoOnComponent(2)=="c0da3");
        self.assertTrue(da4.getInfoOnComponent(3)=="c1da3");
        self.assertTrue(da4.getInfoOnComponent(4)=="c2da3");
        for i in range(35):
            self.assertAlmostEqual(expected1[i],da4.getIJ(0,i),10);
            pass
        # test of static method DataArrayInt::meld
        dai1=da1C.convertToIntArr();
        dai4=DataArrayInt.Meld(dai1,dai3);
        tmp=DataArrayInt.Meld([dai1,dai3]);
        self.assertTrue(dai4.isEqual(tmp))
        self.assertEqual(5,dai4.getNumberOfComponents());
        self.assertEqual(7,dai4.getNumberOfTuples());
        self.assertTrue(dai4.getInfoOnComponent(0)=="c0da1");
        self.assertTrue(dai4.getInfoOnComponent(1)=="c1da1");
        self.assertTrue(dai4.getInfoOnComponent(2)=="c0da3");
        self.assertTrue(dai4.getInfoOnComponent(3)=="c1da3");
        self.assertTrue(dai4.getInfoOnComponent(4)=="c2da3");
        for i in range(35):
            self.assertEqual(int(expected1[i]),dai4.getIJ(0,i));
            pass
        pass

    def testFieldMeld1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setMesh(m);
        da1=DataArrayDouble.New();
        arr1=[12.,23.,34.,45.,56.]
        da1.setValues(arr1,5,1);
        da1.setInfoOnComponent(0,"aaa");
        f1.setArray(da1);
        f1.setTime(3.4,2,1);
        f1.checkConsistencyLight();
        #
        f2=f1.deepCopy();
        f2.setMesh(f1.getMesh());
        f2.checkConsistencyLight();
        f2.changeNbOfComponents(2,5.);
        f2.assign(5.);
        f2.getArray().setInfoOnComponent(0,"bbb");
        f2.getArray().setInfoOnComponent(1,"ccc");
        f2.checkConsistencyLight();
        #
        f3=MEDCouplingFieldDouble.MeldFields(f2,f1);
        f3.checkConsistencyLight();
        self.assertEqual(5,f3.getNumberOfTuples());
        self.assertEqual(3,f3.getNumberOfComponents());
        self.assertTrue(f3.getArray().getInfoOnComponent(0)=="bbb");
        self.assertTrue(f3.getArray().getInfoOnComponent(1)=="ccc");
        self.assertTrue(f3.getArray().getInfoOnComponent(2)=="aaa");
        expected1=[5.,5.,12.,5.,5.,23.,5.,5.,34.,5.,5.,45.,5.,5.,56.]
        for i in range(15):
            self.assertAlmostEqual(expected1[i],f3.getIJ(0,i),12);
            pass
        time,dt,it=f3.getTime();
        self.assertAlmostEqual(3.4,time,14);
        self.assertEqual(2,dt);
        self.assertEqual(1,it);
        #
        f4=f2.buildNewTimeReprFromThis(NO_TIME,False);
        f5=f1.buildNewTimeReprFromThis(NO_TIME,False);
        f6=MEDCouplingFieldDouble.MeldFields(f4,f5);
        f6.checkConsistencyLight();
        self.assertEqual(5,f6.getNumberOfTuples());
        self.assertEqual(3,f6.getNumberOfComponents());
        self.assertTrue(f6.getArray().getInfoOnComponent(0)=="bbb");
        self.assertTrue(f6.getArray().getInfoOnComponent(1)=="ccc");
        self.assertTrue(f6.getArray().getInfoOnComponent(2)=="aaa");
        for i in range(15):
            self.assertAlmostEqual(expected1[i],f6.getIJ(0,i),12);
            pass
        #
        pass

    def testMergeNodes2(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        vec=[0.002,0.]
        m2.translate(vec);
        #
        m3=MEDCouplingUMesh.MergeUMeshes([m1,m2]);
        da,b,newNbOfNodes=m3.mergeNodesCenter(0.01);
        self.assertEqual(9,m3.getNumberOfNodes());
        expected1=[-0.299,-0.3, 0.201,-0.3, 0.701,-0.3, -0.299,0.2, 0.201,0.2, 0.701,0.2, -0.299,0.7, 0.201,0.7, 0.701,0.7]
        for i in range(18):
            self.assertAlmostEqual(expected1[i],m3.getCoords().getIJ(0,i),13);
            pass
        #
        pass

    def testMergeField2(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setMesh(m);
        arr=DataArrayDouble.New();
        arr.alloc(5,2);
        arr.fillWithValue(2.);
        f1.setArray(arr);
        f2=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f2.setMesh(m);
        arr=DataArrayDouble.New();
        arr.alloc(5,2);
        arr.fillWithValue(5.);
        f2.setArray(arr);
        f3=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f3.setMesh(m);
        arr=DataArrayDouble.New();
        arr.alloc(5,2);
        arr.fillWithValue(7.);
        f3.setArray(arr);
        #
        f4=MEDCouplingFieldDouble.MergeFields([f1,f2,f3]);
        self.assertEqual(15,f4.getMesh().getNumberOfCells());
        expected1=[2.,2.,2.,2.,2.,2.,2.,2.,2.,2., 5.,5.,5.,5.,5.,5.,5.,5.,5.,5., 7.,7.,7.,7.,7.,7.,7.,7.,7.,7.]
        for i in range(30):
            self.assertAlmostEqual(expected1[i],f4.getIJ(0,i),13);
            pass
        #
        pass

    def testDAIBuildComplement1(self):
        a=DataArrayInt.New();
        tab=[3,1,7,8]
        a.setValues(tab,4,1);
        b=a.buildComplement(12);
        self.assertEqual(8,b.getNumberOfTuples());
        self.assertEqual(1,b.getNumberOfComponents());
        expected1=[0,2,4,5,6,9,10,11]
        for i in range(8):
            self.assertEqual(expected1[i],b.getIJ(0,i));
            pass
        pass

    def testDAIBuildUnion1(self):
        a=DataArrayInt.New();
        tab1=[3,1,7,8]
        a.setValues(tab1,4,1);
        c=DataArrayInt.New();
        tab2=[5,3,0,18,8]
        c.setValues(tab2,5,1);
        b=a.buildUnion(c);
        self.assertEqual(7,b.getNumberOfTuples());
        self.assertEqual(1,b.getNumberOfComponents());
        expected1=[0,1,3,5,7,8,18]
        for i in range(7):
            self.assertEqual(expected1[i],b.getIJ(0,i));
            pass
        b=DataArrayInt.BuildUnion([a,c]);
        self.assertEqual(7,b.getNumberOfTuples());
        self.assertEqual(1,b.getNumberOfComponents());
        expected1=[0,1,3,5,7,8,18]
        for i in range(7):
            self.assertEqual(expected1[i],b.getIJ(0,i));
            pass
        pass

    def testDAIBuildIntersection1(self):
        a=DataArrayInt.New();
        tab1=[3,1,7,8]
        a.setValues(tab1,4,1);
        c=DataArrayInt.New();
        tab2=[5,3,0,18,8]
        c.setValues(tab2,5,1);
        b=a.buildIntersection(c);
        self.assertEqual(2,b.getNumberOfTuples());
        self.assertEqual(1,b.getNumberOfComponents());
        expected1=[3,8]
        for i in range(2):
            self.assertEqual(expected1[i],b.getIJ(0,i));
            pass
        b=DataArrayInt.BuildIntersection([a,c]);
        self.assertEqual(2,b.getNumberOfTuples());
        self.assertEqual(1,b.getNumberOfComponents());
        expected1=[3,8]
        for i in range(2):
            self.assertEqual(expected1[i],b.getIJ(0,i));
            pass
        pass

    def testDAIDeltaShiftIndex1(self):
        a=DataArrayInt.New();
        tab=[1,3,6,7,7,9,15]
        a.setValues(tab,7,1);
        b=a.deltaShiftIndex();
        self.assertEqual(6,b.getNumberOfTuples());
        self.assertEqual(1,b.getNumberOfComponents());
        expected1=[2,3,1,0,2,6]
        for i in range(6):
            self.assertEqual(expected1[i],b.getIJ(0,i));
            pass
        pass

    def testDaDoubleSelectByTupleIdSafe1(self):
        a=DataArrayDouble.New();
        arr1=[1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1]
        a.setValues(arr1,7,2);
        a.setInfoOnComponent(0,"toto");
        a.setInfoOnComponent(1,"tata");
        #
        arr2=[4,2,0,6,5]
        b=a.selectByTupleIdSafe(arr2);
        self.assertEqual(5,b.getNumberOfTuples());
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertTrue(b.getInfoOnComponent(0)=="toto");
        self.assertTrue(b.getInfoOnComponent(1)=="tata");
        expected1=[5.1,15.1,3.1,13.1,1.1,11.1,7.1,17.1,6.1,16.1]
        for i in range(10):
            self.assertAlmostEqual(expected1[i],b.getIJ(0,i),14);
            pass
        arr4=[4,-1,0,6,5]
        self.assertRaises(InterpKernelException,a.selectByTupleIdSafe,arr4);
        arr5=[4,2,0,6,7]
        self.assertRaises(InterpKernelException,a.selectByTupleIdSafe,arr5);
        #
        c=DataArrayInt.New();
        arr3=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        c.setValues(arr3,7,2);
        c.setInfoOnComponent(0,"toto");
        c.setInfoOnComponent(1,"tata");
        d=c.selectByTupleIdSafe(arr2);
        self.assertEqual(5,d.getNumberOfTuples());
        self.assertEqual(2,d.getNumberOfComponents());
        self.assertTrue(d.getInfoOnComponent(0)=="toto");
        self.assertTrue(d.getInfoOnComponent(1)=="tata");
        expected2=[5,15,3,13,1,11,7,17,6,16]
        for i in range(10):
            self.assertEqual(expected2[i],d.getIJ(0,i));
            pass
        self.assertRaises(InterpKernelException,c.selectByTupleIdSafe,arr4);
        self.assertRaises(InterpKernelException,c.selectByTupleIdSafe,arr5);
        pass

    def testAreCellsIncludedIn1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        pt=[1,3]
        m2=m.buildPartOfMySelf(pt,True);
        ret,tmp=m.areCellsIncludedIn(m2,0)
        self.assertTrue(ret);
        self.assertEqual(2,tmp.getNumberOfTuples());
        self.assertEqual(1,tmp.getNumberOfComponents());
        self.assertEqual(pt[0],tmp.getIJ(0,0));
        self.assertEqual(pt[1],tmp.getIJ(0,1));
        ret,tmp=m2.areCellsIncludedIn(m,0)
        self.assertTrue(not ret);
        m3=MEDCouplingUMesh.MergeUMeshesOnSameCoords(m,m2)
        c,cI=m3.findCommonCells(2,m.getNumberOfCells())
        self.assertTrue(c.isEqual(DataArrayInt([1,5,3,6])))
        self.assertTrue(cI.isEqual(DataArrayInt([0,2,4])))
        pass

    def testSwigErrorProtection1(self):
        m=MEDCouplingDataForTest.build3DTargetMesh_1();
        m.rotate([0.,0.,0.],[0.3,0.6,1.2],0.37)
        m.rotate([0.,0.,0.],[0.3,6,1.2],0.37)
        self.assertRaises(InterpKernelException,m.rotate,[0.,0.,0.],(0.3,6,"1.2"),0.37)
        self.assertRaises(InterpKernelException,m.rotate,[0.,"0.",0.],[0.3,0.6,1.2],0.37)
        self.assertRaises(InterpKernelException,m.rotate,[0.,0.,0.],[0.3,'0.6',1.2],0.37)
        m2=m.buildPartOfMySelf([2,5],True)
        m3=m.buildPartOfMySelf((2,5),True)
        self.assertTrue(m2.isEqual(m3,1e-12))
        self.assertRaises(InterpKernelException,m.buildPartOfMySelf,[2,5.],True)
        da1=m.getCoords().keepSelectedComponents([1])
        da2=m.getCoords().keepSelectedComponents((1,))
        self.assertTrue(da1.isEqual(da2,1e-12))
        self.assertRaises(InterpKernelException,m.getCoords().keepSelectedComponents,["1"])
        pass

    def testDAIBuildSubstraction1(self):
        a=DataArrayInt.New()
        aa=[2,3,6,8,9]
        a.setValues(aa,5,1)
        b=DataArrayInt.New()
        bb=[1,3,5,9,11]
        b.setValues(bb,5,1)
        self.assertEqual([2,6,8],a.buildSubstraction(b).getValues())
        pass

    def testBuildOrthogonalField2(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        d1=DataArrayInt.New();
        d2=DataArrayInt.New();
        d3=DataArrayInt.New();
        d4=DataArrayInt.New();
        m1=m.buildDescendingConnectivity(d1,d2,d3,d4);
        #
        f1=m1.buildOrthogonalField();
        da1=f1.getArray();
        self.assertEqual(2,da1.getNumberOfComponents());
        self.assertEqual(13,da1.getNumberOfTuples());
        #
        expected1=[-1.,0.,0.,1.,1.,0.,0.,-1.,0.707106781186548,0.707106781186548,0.,-1.,0.,1.,1.,0.,0.,1.,1.,0.,-1.,0.,0.,1.,1.,0.];
        for i in range(26):
            self.assertAlmostEqual(expected1[i],da1.getIJ(0,i),14);
            pass
        pass

    def testSwigErrorProtection2(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        coo=m.getCoords()
        c=m.getNodalConnectivity()
        ci=m.getNodalConnectivityIndex()
        del m
        self.assertEqual(2,coo.getNumberOfComponents());
        self.assertEqual(6,ci.getNumberOfTuples());
        self.assertEqual(23,c.getNumberOfTuples());
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f=m.getMeasureField(True)
        c=f.getArray()
        del f
        self.assertEqual(1,c.getNumberOfComponents());
        m=MEDCouplingCMesh.New()
        x=DataArrayDouble.New()
        x.setValues([1.,2.,4.],3,1)
        m.setCoordsAt(0,x)
        del x
        xx=m.getCoordsAt(0)
        del m
        self.assertEqual(3,xx.getNumberOfTuples());
        #
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f=m.getMeasureField(True)
        m2=f.getMesh()
        del m
        del f
        self.assertEqual(5,m2.getNumberOfCells());
        pass

    def testUMInsertNextCell1(self):
        targetCoords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        targetConn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.allocateCells(5);
        self.assertRaises(InterpKernelException,targetMesh.insertNextCell,NORM_QUAD4,4,targetConn[0:4])
        targetMesh.setMeshDimension(2);
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[0:4])
        self.assertRaises(InterpKernelException,targetMesh.insertNextCell,NORM_TETRA4,4,targetConn[0:4])
        self.assertRaises(InterpKernelException,targetMesh.insertNextCell,NORM_SEG2,2,targetConn[0:2])
        self.assertRaises(InterpKernelException,targetMesh.insertNextCell,NORM_POINT1,1,targetConn[0:1])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[4:7])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[7:10])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[10:14])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[14:18])
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,9,2);
        targetMesh.setCoords(myCoords);
        targetMesh.checkConsistencyLight();
        pass

    def testFieldOperatorDivDiffComp1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        m1,d0,d1,d2,d3=m.buildDescendingConnectivity();
        #
        f1=m1.buildOrthogonalField();
        arr1=[2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.]
        arr=DataArrayDouble.New();
        arr.setValues(arr1,13,1);
        f2=MEDCouplingFieldDouble.New(ON_CELLS);
        f2.setArray(arr);
        f2.setMesh(m1);
        f2.checkConsistencyLight();
        #
        f3=f1/f2;
        self.assertRaises(InterpKernelException,f2.__div__,f1)
        f3.checkConsistencyLight();
        f1/=f2;
        #self.assertRaises(InterpKernelException,f2.__idiv__,f1) # mem leaks
        self.assertTrue(f1.isEqual(f3,1e-10,1e-10));
        expected1=[-0.5, 0.0, 0.0, 0.33333333333333331, 0.25, 0.0, 0.0, -0.20000000000000001, 0.117851130197758, 0.117851130197758, 0.0, -0.14285714285714285, 0.0, 0.125, 0.1111111111111111, 0.0, 0.0, 0.10000000000000001, 0.090909090909090912, 0.0, -0.083333333333333329, 0.0, 0.0, 0.076923076923076927, 0.071428571428571425, 0.0]
        for i in range(26):
            self.assertAlmostEqual(expected1[i],f3.getIJ(0,i),10);
            pass
        pass

    def testDARearrange1(self):
        da1=DataArrayInt.New();
        da1.alloc(12,1);
        da1.iota(0);
        #
        self.assertEqual(12,da1.getNbOfElems());
        self.assertEqual(1,da1.getNumberOfComponents());
        self.assertEqual(12,da1.getNumberOfTuples());
        da1.rearrange(4);
        self.assertEqual(12,da1.getNbOfElems());
        self.assertEqual(4,da1.getNumberOfComponents());
        self.assertEqual(3,da1.getNumberOfTuples());
        for i in range(12):
            self.assertEqual(i,da1.getIJ(0,i));
        #
        da1.rearrange(6);
        self.assertEqual(12,da1.getNbOfElems());
        self.assertEqual(6,da1.getNumberOfComponents());
        self.assertEqual(2,da1.getNumberOfTuples());
        for i in range(12):
            self.assertEqual(i,da1.getIJ(0,i));
        #
        self.assertRaises(InterpKernelException,da1.rearrange,7);
        #
        da1.rearrange(12);
        self.assertEqual(12,da1.getNbOfElems());
        self.assertEqual(12,da1.getNumberOfComponents());
        self.assertEqual(1,da1.getNumberOfTuples());
        for i in range(12):
            self.assertEqual(i,da1.getIJ(0,i));
        #
        da1.rearrange(3);
        self.assertEqual(12,da1.getNbOfElems());
        self.assertEqual(3,da1.getNumberOfComponents());
        self.assertEqual(4,da1.getNumberOfTuples());
        for i in range(12):
            self.assertEqual(i,da1.getIJ(0,i));
        #double
        da2=da1.convertToDblArr();
        st=da2.getHiddenCppPointer()
        #
        self.assertEqual(12,da2.getNbOfElems());
        self.assertEqual(3,da2.getNumberOfComponents());
        self.assertEqual(4,da2.getNumberOfTuples());
        da2.rearrange(4);
        self.assertEqual(12,da2.getNbOfElems());
        self.assertEqual(4,da2.getNumberOfComponents());
        self.assertEqual(3,da2.getNumberOfTuples());
        for i in range(12):
            self.assertAlmostEqual(float(i),da2.getIJ(0,i),14);
        #
        da2.rearrange(6);
        self.assertEqual(12,da2.getNbOfElems());
        self.assertEqual(6,da2.getNumberOfComponents());
        self.assertEqual(2,da2.getNumberOfTuples());
        for i in range(12):
            self.assertAlmostEqual(float(i),da2.getIJ(0,i),14);
        #
        self.assertRaises(InterpKernelException,da2.rearrange,7);
        #
        da2.rearrange(1);
        self.assertEqual(st,da2.getHiddenCppPointer())
        self.assertEqual(12,da2.getNbOfElems());
        self.assertEqual(1,da2.getNumberOfComponents());
        self.assertEqual(12,da2.getNumberOfTuples());
        for i in range(12):
            self.assertAlmostEqual(float(i),da2.getIJ(0,i),14);
        #
        da2.rearrange(3);
        self.assertEqual(12,da2.getNbOfElems());
        self.assertEqual(3,da2.getNumberOfComponents());
        self.assertEqual(4,da2.getNumberOfTuples());
        for i in range(12):
            self.assertAlmostEqual(float(i),da2.getIJ(0,i),14);
        pass

    def testDARearrange2(self):
        da1=DataArrayInt.New();
        arr=[1,2,3,2,2,3,5,1,5,5,2,2]
        da1.setValues(arr,4,3);
        s=da1.getDifferentValues();
        expected1=DataArrayInt([1,2,3,5])
        self.assertTrue(expected1.isEqual(s));
        pass

    def testSwigErrorProtection3(self):
        da=DataArrayInt.New()
        da.setValues([1,2,3,4,0,0,0,0,0,0,0,0],4,3)
        self.assertEqual([1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        da=DataArrayInt.New()
        da.setValues(((1,2,3),(4,4,3),(0,0,0),(0,0,0)),4,3)
        self.assertEqual([1, 2, 3, 4, 4, 3, 0, 0, 0, 0, 0, 0],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        da.setValues((10*[1]+290*[2])[:12],4,3)
        self.assertEqual(10*[1]+[2,2],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        #
        da=DataArrayDouble.New()
        da.setValues([1,2,3.,4,0,0,0,0,0,0,0,0],4,3)
        self.assertEqual([1., 2., 3., 4., 0., 0., 0., 0., 0., 0., 0., 0.],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        da=DataArrayDouble.New()
        da.setValues(((1,2,3),(4.,4,3),(0,0,0),(0,0,0)),4,3)
        self.assertEqual([1., 2., 3., 4., 4., 3., 0., 0., 0., 0., 0., 0.],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        da.setValues((10*[1]+290*[2])[:12],4,3)
        self.assertEqual(10*[1.]+[2.,2.],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        pass

    def testDAIBuildPermutationArr1(self):
        a=DataArrayInt.New()
        a.setValues([4,5,6,7,8],5,1)
        b=DataArrayInt.New()
        b.setValues([5,4,8,6,7],5,1)
        c=a.buildPermutationArr(b)
        self.assertEqual([1,0,4,2,3],c.getValues())
        self.assertTrue(a.isEqualWithoutConsideringStrAndOrder(b))
        b.setIJ(0,0,9)
        self.assertTrue(not a.isEqualWithoutConsideringStrAndOrder(b))
        self.assertRaises(InterpKernelException,a.buildPermutationArr,b)
        a.setIJ(3,0,4)
        b.setIJ(0,0,5)
        b.setIJ(4,0,4)#a==[4,5,6,4,8] and b==[5,4,8,6,4]
        self.assertTrue(a.isEqualWithoutConsideringStrAndOrder(b))
        c=a.buildPermutationArr(b)
        self.assertEqual([1,3,4,2,3],c.getValues())
        d=b.convertToDblArr()
        expect3=[4,4,5,6,8]
        b.sort()
        self.assertEqual(expect3,b.getValues())
        d.sort()
        self.assertEqual(5,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        for i in range(5):
            self.assertAlmostEqual(float(expect3[i]),d.getIJ(i,0),14);
            pass
        pass

    def testAreCellsIncludedIn2(self):
        myName="Vitoo";
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m2=m.buildPartOfMySelf([],True);
        self.assertEqual(0,m2.getNumberOfCells());
        self.assertEqual(3,m2.getSpaceDimension());
        self.assertEqual(2,m2.getMeshDimension());
        m2.setName(myName);
        test,tmp=m.areCellsIncludedIn(m2,0)
        self.assertTrue(test);
        self.assertEqual(myName,tmp.getName());
        self.assertEqual(0,tmp.getNumberOfTuples())
        self.assertEqual(1,tmp.getNumberOfComponents())
        pass

    def testUMeshGetPartBarycenterAndOwner1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        part1=[1,0,4];
        part=DataArrayInt.New();
        part.setValues(part1,3,1);
        b=m1.getPartBarycenterAndOwner(part);
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertEqual(3,b.getNumberOfTuples());
        expected1=[0.36666666666666665,-0.13333333333333333,-0.05,-0.05,0.45,0.45];
        for i in range(6):
            self.assertAlmostEqual(expected1[i],b.getIJ(0,i),14);
            pass
        pass

    def testUMeshGetPartMeasureField1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        part1=[1,0,4];
        part=DataArrayInt.New();
        part.setValues(part1,3,1);
        b=m1.getPartMeasureField(True,part);
        self.assertEqual(1,b.getNumberOfComponents());
        self.assertEqual(3,b.getNumberOfTuples());
        expected1=[0.125,0.25,0.25];
        for i in range(3):
            self.assertAlmostEqual(expected1[i],b.getIJ(0,i),14);
            pass
        pass

    def testUMeshBuildPartOrthogonalField1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m1.changeSpaceDimension(3);
        part1=[1,0,4];
        part=DataArrayInt.New();
        part.setValues(part1,3,1);
        b=m1.buildPartOrthogonalField(part);
        self.assertEqual(3,b.getArray().getNumberOfComponents());
        self.assertEqual(3,b.getArray().getNumberOfTuples());
        expected1=[0.,0.,-1.,0.,0.,-1.,0.,0.,-1.];
        for i in range(9):
            self.assertAlmostEqual(expected1[i],b.getArray().getIJ(0,i),14);
            pass
        pass

    def testUMeshGetTypesOfPart1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        part1=[0,3,4];
        p1=DataArrayInt.New()
        p1.setValues(part1,3,1)
        s=m1.getTypesOfPart(p1);
        self.assertEqual([NORM_QUAD4],s);
        part2=[2,2,2,1];
        p2=DataArrayInt.New()
        p2.setValues(part2,4,1)
        s=m1.getTypesOfPart(p2);
        self.assertEqual([NORM_TRI3],s);
        part3=[3,2,1];
        p3=DataArrayInt.New()
        p3.setValues(part3,3,1)
        s=m1.getTypesOfPart(p3);
        self.assertEqual(s,[NORM_TRI3,NORM_QUAD4]);
        pass

    def testUMeshKeepCellIdsByType1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        part1=[0,3,4]
        p1=DataArrayInt.New()
        p1.setValues(part1,3,1)
        p1.setName("p1")
        a=m1.keepCellIdsByType(NORM_TRI3,p1);
        self.assertEqual("p1",a.getName())
        self.assertEqual(1,a.getNumberOfComponents());
        self.assertEqual(0,a.getNumberOfTuples());
        #
        part2=[3,2,0,2,4]
        p2=DataArrayInt.New()
        p2.setValues(part2,5,1)
        p2.setName("p2")
        a=m1.keepCellIdsByType(NORM_TRI3,p2);
        self.assertEqual("p2",a.getName())
        self.assertEqual(1,a.getNumberOfComponents());
        self.assertEqual(2,a.getNumberOfTuples());
        self.assertEqual(2,a.getIJ(0,0));
        self.assertEqual(2,a.getIJ(1,0));
        #
        a=m1.keepCellIdsByType(NORM_QUAD4,p2);
        self.assertEqual("p2",a.getName())
        self.assertEqual(1,a.getNumberOfComponents());
        self.assertEqual(3,a.getNumberOfTuples());
        self.assertEqual(3,a.getIJ(0,0));
        self.assertEqual(0,a.getIJ(1,0));
        self.assertEqual(4,a.getIJ(2,0));
        pass
    
    def testSwigErrorDaIntSelectByTupleId1(self):
        a=DataArrayInt.New();
        arr1=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        a.setValues(arr1,7,2);
        a.setInfoOnComponent(0,"toto");
        a.setInfoOnComponent(1,"tata");
        #
        arr2=[4,2,0,6,5]
        b=a.selectByTupleId(arr2);
        self.assertEqual(5,b.getNumberOfTuples());
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertTrue(b.getInfoOnComponent(0)=="toto");
        self.assertTrue(b.getInfoOnComponent(1)=="tata");
        expected1=[5,15,3,13,1,11,7,17,6,16]
        self.assertEqual(expected1,b.getValues())
        #
        a2=DataArrayInt.New()
        a2.setValues(arr2,5,1)
        b=a.selectByTupleId(a2);
        self.assertEqual(5,b.getNumberOfTuples());
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertTrue(b.getInfoOnComponent(0)=="toto");
        self.assertTrue(b.getInfoOnComponent(1)=="tata");
        expected1=[5,15,3,13,1,11,7,17,6,16]
        self.assertEqual(expected1,b.getValues())
        pass

    def testSwigErrorRenum(self):
        da=DataArrayDouble.New()
        da.setValues([7.,107.,8.,108.,9.,109.,10.,110.,11.,111.,12.,112.,13.,113.,14.,114.,15.,115.,16.,116.],10,2)
        d=DataArrayInt.New()
        d.setValues([0,2,3,1,4,5,6,8,7,9],10,1)
        da.renumberInPlace(d)
        da.renumber(d)
        pass
      
    pass

if __name__ == '__main__':
    unittest.main()
