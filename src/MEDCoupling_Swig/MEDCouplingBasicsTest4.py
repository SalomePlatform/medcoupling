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
import rlcompleter,readline # this line has to be here, to ensure a usability of MEDCoupling/MEDLoader. B4 removing it please notify to anthony.geay@edf.fr
from sys import platform

def checkFreeMemory(size):
    """
    Get node total memory and memory usage
    """
    ret = True
    dic = {}
    if platform not in ["win32"]:
        with open('/proc/meminfo', 'r') as mem:
            tmp = 0
            for i in mem:
                sline = i.split()
                if str(sline[0]) == 'MemTotal:':
                    dic['total'] = int(sline[1])
                elif str(sline[0]) in ('MemFree:', 'Buffers:', 'Cached:'):
                    tmp += int(sline[1])
            dic['free'] = tmp
            dic['used'] = int(dic['total']) - int(dic['free'])
            ret = dic['free'] > size
    #TODO: extend this method for Windows OS            
    return ret


class MEDCouplingBasicsTest4(unittest.TestCase):
    def testSwigDADOp4(self):
        da = DataArrayDouble.New(list(range(6, 30)), 12, 2)
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in range(24):
            self.assertAlmostEqual(da.getIJ(0,i),float(i+6),13)
            pass
        # operator transpose
        da.transpose()
        self.assertEqual(2,da.getNumberOfTuples());
        self.assertEqual(12,da.getNumberOfComponents());
        for i in range(24):
            self.assertAlmostEqual(da.getIJ(0,i),float(i+6),13)
            pass
        da.transpose()
        # operator __neg__
        da2=DataArrayDouble.New(12,1)
        da2.iota(0.)
        dabis=-da
        for i in range(24):
            self.assertAlmostEqual(dabis.getIJ(0,i),-float(i+6),13)
            pass
        # operator+=
        da+=da2
        expected1=[6.,7.,9.,10.,12.,13.,15.,16.,18.,19.,21.,22.,24.,25.,27.,28.,30.,31.,33.,34.,36.,37.,39.,40.]
        for i in range(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected1[i],13)
            pass
        da=-dabis
        da+=[100.,101.]
        expected2=[106.,108.,108.,110.,110.,112.,112.,114.,114.,116.,116.,118.,118.,120.,120.,122.,122.,124.,124.,126.,126.,128.,128.,130.]
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in range(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected2[i],13)
            pass
        for pos,elt in enumerate(dabis):
            da[pos]+=elt
            pass
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for elt in da:
            li=elt[:]
            self.assertAlmostEqual(li[0],100.,13) ; self.assertAlmostEqual(li[1],101.,13)
            pass
        # operator-=
        da = DataArrayDouble.New(list(range(6, 30)), 12, 2)
        da2 = DataArrayDouble.New(list(range(12)), 12, 1)
        dabis=-da
        da-=da2
        expected1=[6.,7.,7.,8.,8.,9.,9.,10.,10.,11.,11.,12.,12.,13.,13.,14.,14.,15.,15.,16.,16.,17.,17.,18.]
        for i in range(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected1[i],13)
            pass
        da=-dabis
        da-=[100.,101.]
        expected2=[-94.,-94.,-92.,-92.,-90.,-90.,-88.,-88.,-86.,-86.,-84.,-84.,-82.,-82.,-80.,-80.,-78.,-78.,-76.,-76.,-74.,-74.,-72.,-72.]
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in range(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected2[i],13)
            pass
        for pos,elt in enumerate(dabis):
            da[pos]-=elt
            pass
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        expected3=[-88.,-87.,-84.,-83.,-80.,-79.,-76.,-75.,-72.,-71.,-68.,-67.,-64.,-63.,-60.,-59.,-56.,-55.,-52.,-51.,-48.,-47.,-44.,-43.]
        for i in range(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected3[i],13)
            pass
        # operator*=
        da = DataArrayDouble.New(list(range(6, 30)), 12, 2)
        da2 = DataArrayDouble.New(list(range(12)), 12, 1)
        dabis=-da
        da*=da2
        expected1=[0.,0.,8.,9.,20.,22.,36.,39.,56.,60.,80.,85.,108.,114.,140.,147.,176.,184.,216.,225.,260.,270.,308.,319.]
        for i in range(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected1[i],13)
            pass
        da=-dabis
        da*=[100.,101.]
        expected2=[600.,707.,800.,909.,1000.,1111.,1200.,1313.,1400.,1515.,1600.,1717.,1800.,1919.,2000.,2121.,2200.,2323.,2400.,2525.,2600.,2727.,2800.,2929.]
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in range(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected2[i],13)
            pass
        for pos,elt in enumerate(dabis):
            da[pos]*=elt
            pass
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        expected3=[-3600.,-4949.,-6400.,-8181.,-10000.,-12221.,-14400.,-17069.,-19600.,-22725.,-25600.,-29189.,-32400.,-36461.,-40000.,-44541.,-48400.,-53429.,-57600.,-63125.,-67600.,-73629.,-78400.,-84941.0]
        for i in range(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected3[i],13)
            pass
        # operator/=
        da = DataArrayDouble.New(list(range(6, 30)), 12, 2)
        da2 = DataArrayDouble.New(list(range(1, 13)), 12, 1)
        dabis=-da
        da/=da2
        expected1=[6.0,7.0,4.0,4.5,3.3333333333333335,3.6666666666666665,3.0,3.25,2.8,3.0,2.6666666666666665,2.8333333333333335,2.5714285714285716,2.7142857142857144,2.5,2.625,2.4444444444444446,2.5555555555555554,2.4,2.5,2.3636363636363638,2.4545454545454546,2.3333333333333335,2.4166666666666665]
        for i in range(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected1[i],13)
            pass
        da=-dabis
        da/=[100.,101.]
        expected2=[0.06,0.06930693069306931,0.08,0.0891089108910891,0.1,0.10891089108910891,0.12,0.12871287128712872,0.14,0.1485148514851485,0.16,0.16831683168316833,0.18,0.18811881188118812,0.2,0.2079207920792079,0.22,0.22772277227722773,0.24,0.24752475247524752,0.26,0.26732673267326734,0.28,0.2871287128712871]
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in range(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected2[i],13)
            pass
        for pos,elt in enumerate(dabis):
            da[pos]/=elt
            pass
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        expected3=[-0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.0099009900990099]
        for i in range(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected3[i],13)
            pass
        pass

    def testSwigDAIOp4(self):
        da = DataArrayInt.New(list(range(6, 30)), 12, 2)
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in range(24):
            self.assertEqual(da.getIJ(0,i),i+6)
            pass
        # operator transpose
        da.transpose()
        self.assertEqual(2,da.getNumberOfTuples());
        self.assertEqual(12,da.getNumberOfComponents());
        for i in range(24):
            self.assertEqual(da.getIJ(0,i),i+6)
            pass
        da.transpose()
        # operator __neg__
        da2=DataArrayInt.New(12,1)
        da2.iota(0)
        dabis=-da
        for i in range(24):
            self.assertEqual(dabis.getIJ(0,i),-(i+6))
            pass
        # operator+=
        da+=da2
        expected1=[6,7,9,10,12,13,15,16,18,19,21,22,24,25,27,28,30,31,33,34,36,37,39,40]
        for i in range(24):
            self.assertEqual(da.getIJ(0,i),expected1[i])
            pass
        da=-dabis
        da+=[100,101]
        expected2=[106,108,108,110,110,112,112,114,114,116,116,118,118,120,120,122,122,124,124,126,126,128,128,130]
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in range(24):
            self.assertEqual(da.getIJ(0,i),expected2[i])
            pass
        for pos,elt in enumerate(dabis):
            da[pos]+=elt
            pass
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for elt in da:
            li=elt[:]
            self.assertEqual(li[0],100) ; self.assertEqual(li[1],101)
            pass
        # operator-=
        da = DataArrayInt.New(list(range(6, 30)), 12, 2)
        da2 = DataArrayInt.New(list(range(12)), 12, 1)
        dabis=-da
        da-=da2
        expected1=[6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18]
        for i in range(24):
            self.assertEqual(da.getIJ(0,i),expected1[i])
            pass
        da=-dabis
        da-=[100,101]
        expected2=[-94,-94,-92,-92,-90,-90,-88,-88,-86,-86,-84,-84,-82,-82,-80,-80,-78,-78,-76,-76,-74,-74,-72,-72]
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in range(24):
            self.assertEqual(da.getIJ(0,i),expected2[i])
            pass
        for pos,elt in enumerate(dabis):
            da[pos]-=elt
            pass
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        expected3=[-88,-87,-84,-83,-80,-79,-76,-75,-72,-71,-68,-67,-64,-63,-60,-59,-56,-55,-52,-51,-48,-47,-44,-43]
        for i in range(24):
            self.assertEqual(da.getIJ(0,i),expected3[i])
            pass
        # operator*=
        da = DataArrayInt.New(list(range(6, 30)), 12, 2)
        da2 = DataArrayInt.New(list(range(12)), 12, 1)
        dabis=-da
        da*=da2
        expected1=[0,0,8,9,20,22,36,39,56,60,80,85,108,114,140,147,176,184,216,225,260,270,308,319]
        for i in range(24):
            self.assertEqual(da.getIJ(0,i),expected1[i])
            pass
        da=-dabis
        da*=[100,101]
        expected2=[600,707,800,909,1000,1111,1200,1313,1400,1515,1600,1717,1800,1919,2000,2121,2200,2323,2400,2525,2600,2727,2800,2929]
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in range(24):
            self.assertEqual(da.getIJ(0,i),expected2[i])
            pass
        for pos,elt in enumerate(dabis):
            da[pos]*=elt
            pass
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        expected3=[-3600,-4949,-6400,-8181,-10000,-12221,-14400,-17069,-19600,-22725,-25600,-29189,-32400,-36461,-40000,-44541,-48400,-53429,-57600,-63125,-67600,-73629,-78400,-84941.0]
        for i in range(24):
            self.assertEqual(da.getIJ(0,i),expected3[i])
            pass
        # operator/=
        da = DataArrayInt.New(list(range(6, 30)), 12, 2)
        da2 = DataArrayInt.New(list(range(1, 13)), 12, 1)
        dabis=-da
        da/=da2
        expected1=[6,7,4,4,3,3,3,3,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
        for i in range(24):
            self.assertEqual(da.getIJ(0,i),expected1[i])
            pass
        da=-dabis
        da/=DataArrayInt.New([2,3],1,2)
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        expected2=[3,2,4,3,5,3,6,4,7,5,8,5,9,6,10,7,11,7,12,8,13,9,14,9]
        for i in range(24):
            self.assertEqual(da.getIJ(0,i),expected2[i])
            pass
        pass

    def testSwigDADOp5(self):
        da=DataArrayDouble.New([5,6,7,8,9,6,7,-2,3,9,8,10])
        da.rearrange(3)
        da2=DataArrayDouble.New([5.,8.,10.,12])
        self.assertEqual(4,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        da3=da+da2
        self.assertEqual(4,da3.getNumberOfTuples());
        self.assertEqual(3,da3.getNumberOfComponents());
        expected1=[10.,11.,12.,16.,17.,14.,17.,8.,13.,21.,20.,22.]
        for i in range(12):
            self.assertAlmostEqual(da3.getIJ(0,i),expected1[i],13)
            pass
        da3=da2+da
        self.assertEqual(4,da3.getNumberOfTuples());
        self.assertEqual(3,da3.getNumberOfComponents());
        for i in range(12):
            self.assertAlmostEqual(da3.getIJ(0,i),expected1[i],13)
            pass
        # Test new API of classmethod DataArrayDouble.New
        vals=[5,6,7,8,9,6,7,-2,3,9,8,10]
        da=DataArrayDouble.New(vals)
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        for i in range(12):
            self.assertAlmostEqual(da.getIJ(0,i),vals[i],13)
            pass
        da=DataArrayDouble.New(vals,12)
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        for i in range(12):
            self.assertAlmostEqual(da.getIJ(0,i),vals[i],13)
            pass
        da=DataArrayDouble.New(vals,1,12)
        self.assertEqual(1,da.getNumberOfTuples());
        self.assertEqual(12,da.getNumberOfComponents());
        for i in range(12):
            self.assertAlmostEqual(da.getIJ(0,i),vals[i],13)
            pass
        da=DataArrayDouble.New(vals,6,2)
        self.assertEqual(6,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in range(12):
            self.assertAlmostEqual(da.getIJ(0,i),vals[i],13)
            pass
        da=DataArrayDouble.New(vals,4,3)
        self.assertEqual(4,da.getNumberOfTuples());
        self.assertEqual(3,da.getNumberOfComponents());
        for i in range(12):
            self.assertAlmostEqual(da.getIJ(0,i),vals[i],13)
            pass
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,11);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,13);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,5,2);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,7,2);
        pass

    def testSwigDADOp6(self):
        da=DataArrayInt.New([5,6,7,8,9,6,7,-2,3,9,8,10])
        da.rearrange(3)
        da2=DataArrayInt.New([5,8,10,12])
        self.assertEqual(4,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        da3=da+da2
        self.assertEqual(4,da3.getNumberOfTuples());
        self.assertEqual(3,da3.getNumberOfComponents());
        expected1=[10,11,12,16,17,14,17,8,13,21,20,22]
        for i in range(12):
            self.assertEqual(da3.getIJ(0,i),expected1[i])
            pass
        da3=da2+da
        self.assertEqual(4,da3.getNumberOfTuples());
        self.assertEqual(3,da3.getNumberOfComponents());
        for i in range(12):
            self.assertEqual(da3.getIJ(0,i),expected1[i])
            pass
        da3=da+DataArrayInt.New(da2.getValues())
        # Test new API of classmethod DataArrayInt.New
        vals=[5,6,7,8,9,6,7,-2,3,9,8,10]
        da=DataArrayDouble.New(vals)
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        for i in range(12):
            self.assertEqual(da.getIJ(0,i),vals[i])
            pass
        da=DataArrayDouble.New(vals,12)
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        for i in range(12):
            self.assertEqual(da.getIJ(0,i),vals[i])
            pass
        da=DataArrayDouble.New(vals,1,12)
        self.assertEqual(1,da.getNumberOfTuples());
        self.assertEqual(12,da.getNumberOfComponents());
        for i in range(12):
            self.assertEqual(da.getIJ(0,i),vals[i])
            pass
        da=DataArrayDouble.New(vals,6,2)
        self.assertEqual(6,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in range(12):
            self.assertEqual(da.getIJ(0,i),vals[i])
            pass
        da=DataArrayDouble.New(vals,4,3)
        self.assertEqual(4,da.getNumberOfTuples());
        self.assertEqual(3,da.getNumberOfComponents());
        for i in range(12):
            self.assertEqual(da.getIJ(0,i),vals[i])
            pass
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,11);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,13);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,5,2);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,7,2);
        pass

    def testSwigDADOp9(self):
        l1=[(1.,2.,3),(4.,5.,6.),(7.,8.,9.),[10.,11.,12.]]
        da1=DataArrayDouble(l1,4,3)
        self.assertEqual(4,da1.getNumberOfTuples());
        self.assertEqual(3,da1.getNumberOfComponents());
        da2=DataArrayDouble(12) ; da2.iota(1.) ; da2.rearrange(3)
        self.assertTrue(da2.isEqual(da1,1e-12))
        self.assertRaises(InterpKernelException,DataArrayDouble.New,l1,3,4);
        da3=DataArrayDouble(l1,4)
        self.assertTrue(da3.isEqual(da1,1e-12))
        self.assertRaises(InterpKernelException,DataArrayDouble.New,l1,3);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,l1,5);
        l1=[(1.,2.,3),(4.,(5.),((6.))),(7.,8.,9.),[10.,11.,12.]]
        da1=DataArrayDouble(l1,4,3)
        self.assertEqual(4,da1.getNumberOfTuples());
        self.assertEqual(3,da1.getNumberOfComponents());
        da2=DataArrayDouble(12) ; da2.iota(1.) ; da2.rearrange(3)
        self.assertTrue(da2.isEqual(da1,1e-12))
        self.assertRaises(InterpKernelException,DataArrayDouble.New,l1,3,4);
        da3=DataArrayDouble(l1,4)
        self.assertTrue(da3.isEqual(da1,1e-12))
        self.assertRaises(InterpKernelException,DataArrayDouble.New,l1,3);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,l1,5);
        #
        l1=[(1,2,3),(4,5,6),(7,8,9),[10,11,12]]
        da1=DataArrayInt(l1,4,3)
        self.assertEqual(4,da1.getNumberOfTuples());
        self.assertEqual(3,da1.getNumberOfComponents());
        da2=DataArrayInt(12) ; da2.iota(1) ; da2.rearrange(3)
        self.assertTrue(da2.isEqual(da1))
        self.assertRaises(InterpKernelException,DataArrayInt.New,l1,3,4);
        da3=DataArrayInt(l1,4)
        self.assertTrue(da3.isEqual(da1))
        self.assertRaises(InterpKernelException,DataArrayInt.New,l1,3);
        self.assertRaises(InterpKernelException,DataArrayInt.New,l1,5);
        l1=[(1,[2],3),(4,[(5)],6),((([7])),8,9),[10,11,12]]
        da1=DataArrayInt(l1,4,3)
        self.assertEqual(4,da1.getNumberOfTuples());
        self.assertEqual(3,da1.getNumberOfComponents());
        da2=DataArrayInt(12) ; da2.iota(1) ; da2.rearrange(3)
        self.assertTrue(da2.isEqual(da1))
        self.assertRaises(InterpKernelException,DataArrayInt.New,l1,3,4);
        da3=DataArrayInt(l1,4)
        self.assertTrue(da3.isEqual(da1))
        self.assertRaises(InterpKernelException,DataArrayInt.New,l1,3);
        self.assertRaises(InterpKernelException,DataArrayInt.New,l1,5);
        pass

    def testRenumberNodesInConn1(self):
        mesh2DCoords=[-0.3,-0.3,0., 0.2,-0.3,0., 0.7,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.7,0.2,0., -0.3,0.7,0., 0.2,0.7,0., 0.7,0.7,0. ]
        mesh2DConn=[1,4,2, 4,5,2, 0,3,4,1, 6,7,4,3, 7,8,5,4]
        mesh2D=MEDCouplingUMesh.New("mesh",2);
        mesh2D.allocateCells(5);
        mesh2D.insertNextCell(NORM_TRI3,3,mesh2DConn[0:3])
        mesh2D.insertNextCell(NORM_TRI3,3,mesh2DConn[3:6])
        mesh2D.insertNextCell(NORM_QUAD4,4,mesh2DConn[6:10])
        mesh2D.insertNextCell(NORM_QUAD4,4,mesh2DConn[10:14])
        mesh2D.insertNextCell(NORM_QUAD4,4,mesh2DConn[14:18])
        mesh2D.finishInsertingCells();
        myCoords=DataArrayDouble.New(mesh2DCoords,9,3);
        mesh2D.setCoords(myCoords);
        mesh2D.checkConsistencyLight();
        #
        mesh3DCoords=[-0.3,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.2,-0.3,0., -0.3,-0.3,1., -0.3,0.2,1., 0.2,0.2,1., 0.2,-0.3,1. ]
        mesh3DConn=[0,1,2,3,4,5,6,7]
        mesh3D=MEDCouplingUMesh.New("mesh",3);
        mesh3D.allocateCells(1);
        mesh3D.insertNextCell(NORM_HEXA8,8,mesh3DConn[:])
        mesh3D.finishInsertingCells();
        myCoords3D=DataArrayDouble.New(mesh3DCoords,8,3);
        mesh3D.setCoords(myCoords3D);
        mesh3D.checkConsistencyLight();
        #
        mesh3D_2=mesh3D.deepCopy();
        mesh2D_2=mesh2D.deepCopy();
        mesh3D_4=mesh3D.deepCopy();
        mesh2D_4=mesh2D.deepCopy();
        oldNbOf3DNodes=mesh3D.getNumberOfNodes();
        renumNodes=DataArrayInt.New();
        renumNodes.alloc(mesh2D.getNumberOfNodes(),1);
        renumNodes.iota(oldNbOf3DNodes);
        coo=DataArrayDouble.Aggregate(mesh3D.getCoords(),mesh2D.getCoords());
        mesh3D.setCoords(coo);
        mesh2D.setCoords(coo);
        mesh2DCpy=mesh2D.deepCopy()
        mesh2D_3=mesh2D.deepCopy();
        mesh2D_3.shiftNodeNumbersInConn(oldNbOf3DNodes);
        mesh2D.renumberNodesInConn(renumNodes);
        mesh2DCpy.renumberNodesInConn(renumNodes.getValues());
        self.assertTrue(mesh2D.isEqual(mesh2DCpy,1e-12))
        self.assertTrue(mesh2D.isEqual(mesh2D_3,1e-12))
        #
        da1,da2=mesh3D.checkGeoEquivalWith(mesh3D_2,10,1e-12);
        self.assertTrue(da1==None);
        self.assertEqual(8,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        expected1=[8,11,12,9,4,5,6,7]
        for i in range(8):
            self.assertEqual(expected1[i],da2.getIJ(i,0));
            pass
        #
        da1,da2=mesh2D.checkGeoEquivalWith(mesh2D_2,10,1e-12);
        self.assertTrue(da1==None);
        self.assertEqual(9,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        for i in range(9):
            self.assertEqual(8+i,da2.getIJ(i,0));
            pass
        #
        mesh2D_5=mesh2D_4.deepCopy();
        mesh2D_5.translate([1.,0.,0.]);
        meshes=[mesh3D_4,mesh2D_4,mesh2D_5];
        MEDCouplingUMesh.PutUMeshesOnSameAggregatedCoords(meshes);
        self.assertTrue(mesh3D_4.getCoords().getHiddenCppPointer()==mesh2D_4.getCoords().getHiddenCppPointer());
        self.assertTrue(mesh2D_4.getCoords().getHiddenCppPointer()==mesh2D_5.getCoords().getHiddenCppPointer());
        mesh3D_4.checkConsistencyLight(); mesh2D_4.checkConsistencyLight(); mesh2D_5.checkConsistencyLight();
        self.assertEqual(26,mesh3D_4.getNumberOfNodes());
        self.assertEqual(3,mesh3D_4.getSpaceDimension());
        self.assertEqual(9,mesh3D_4.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(23,mesh2D_4.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(23,mesh2D_5.getNodalConnectivity().getNumberOfTuples());
        expected2=[18,0,1,2,3,4,5,6,7]
        expected3=[3,9,12,10, 3,12,13,10, 4,8,11,12,9, 4,14,15,12,11, 4,15,16,13,12]
        expected4=[3,18,21,19, 3,21,22,19, 4,17,20,21,18, 4,23,24,21,20, 4,24,25,22,21]
        expected5=[-0.3,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.2,-0.3,0., -0.3,-0.3,1., -0.3,0.2,1., 0.2,0.2,1., 0.2,-0.3,1., -0.3,-0.3,0., 0.2,-0.3,0., 0.7,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.7,0.2,0., -0.3,0.7,0., 0.2,0.7,0., 0.7,0.7,0., 0.7, -0.3, 0.0, 1.2, -0.3, 0.0, 1.7, -0.3, 0.0, 0.7, 0.2, 0.0, 1.2, 0.2, 0.0, 1.7, 0.2, 0.0, 0.7, 0.7, 0.0, 1.2, 0.7, 0.0, 1.7, 0.7, 0.0]
        self.assertEqual(expected2,mesh3D_4.getNodalConnectivity().getValues());
        self.assertEqual(expected3,mesh2D_4.getNodalConnectivity().getValues());
        self.assertEqual(expected4,mesh2D_5.getNodalConnectivity().getValues());
        for i in range(78):
            self.assertAlmostEqual(expected5[i],mesh3D_4.getCoords().getIJ(0,i),12);
            pass
        #
        MEDCouplingUMesh.MergeNodesOnUMeshesSharingSameCoords(meshes,1e-12);
        mesh3D_4.checkConsistencyLight(); mesh2D_4.checkConsistencyLight(); mesh2D_5.checkConsistencyLight();
        self.assertTrue(mesh3D_4.getCoords().getHiddenCppPointer()==mesh2D_4.getCoords().getHiddenCppPointer());
        self.assertTrue(mesh2D_4.getCoords().getHiddenCppPointer()==mesh2D_5.getCoords().getHiddenCppPointer());
        self.assertEqual(19,mesh3D_4.getNumberOfNodes());
        self.assertEqual(3,mesh3D_4.getSpaceDimension());
        self.assertEqual(9,mesh3D_4.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(23,mesh2D_4.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(23,mesh2D_5.getNodalConnectivity().getNumberOfTuples());
        expected6=[18,0,1,2,3,4,5,6,7]
        expected7=[3,3,2,8, 3,2,9,8, 4,0,1,2,3, 4,10,11,2,1, 4,11,12,9,2]
        expected8=[3,13,15,14, 3,15,16,14, 4,8,9,15,13, 4,12,17,15,9, 4,17,18,16,15]
        expected9=[-0.3, -0.3, 0., -0.3, 0.2, 0., 0.2, 0.2, 0., 0.2, -0.3, 0., -0.3, -0.3, 1., -0.3, 0.2, 1.,
                    0.2, 0.2, 1., 0.2, -0.3, 1., 0.7, -0.3, 0., 0.7, 0.2, 0., -0.3, 0.7, 0., 0.2, 0.7, 0.,
                    0.7, 0.7, 0., 1.2, -0.3, 0., 1.7, -0.3, 0., 1.2, 0.2, 0., 1.7, 0.2, 0., 1.2, 0.7, 0., 1.7, 0.7, 0.]
        self.assertEqual(expected6,mesh3D_4.getNodalConnectivity().getValues());
        self.assertEqual(expected7,mesh2D_4.getNodalConnectivity().getValues());
        self.assertEqual(expected8,mesh2D_5.getNodalConnectivity().getValues());
        for i in range(57):
            self.assertAlmostEqual(expected9[i],mesh3D_4.getCoords().getIJ(0,i),12);
            pass
        #
        pass
    
    def testComputeNeighborsOfCells1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        d1,d2=m.computeNeighborsOfCells();
        self.assertEqual(6,d2.getNumberOfTuples());
        self.assertEqual(10,d1.getNumberOfTuples());
        expected1=[0,2,4,6,8,10]
        expected2=[3,1,0,2,4,1,4,0,2,3]
        self.assertEqual(expected1,d2.getValues());
        self.assertEqual(expected2,d1.getValues());
        pass

    def testCheckButterflyCellsBug1(self):
        mesh2DCoords=[323.85,120.983748908684,317.5,131.982271536747,336.55,120.983748908686,330.2,131.982271536751,323.85,142.98079416481]
        mesh2DConn=[4,1,0,2,3]
        mesh2D=MEDCouplingUMesh.New("mesh",2);
        mesh2D.allocateCells(1);
        mesh2D.insertNextCell(NORM_POLYGON,5,mesh2DConn[0:5])
        mesh2D.finishInsertingCells();
        myCoords=DataArrayDouble.New(mesh2DCoords,5,2);
        mesh2D.setCoords(myCoords);
        mesh2D.checkConsistencyLight();
        #
        v=mesh2D.checkButterflyCells();
        self.assertTrue(v.empty());
        pass

    def testDataArrayIntRange1(self):
        d=DataArrayInt.Range(2,17,7);
        expected1=[2,9,16]
        self.assertEqual(3,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        self.assertEqual(expected1,d.getValues());
        #
        d=DataArrayInt.Range(2,23,7);
        self.assertEqual(3,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        self.assertEqual(expected1,d.getValues());
        #
        d=DataArrayInt.Range(2,24,7);
        expected2=[2,9,16,23]
        self.assertEqual(4,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        self.assertEqual(expected2,d.getValues());
        #
        d=DataArrayInt.Range(24,2,-7);
        expected3=[24,17,10,3]
        self.assertEqual(4,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        self.assertEqual(expected3,d.getValues());
        #
        d=DataArrayInt.Range(23,2,-7);
        expected4=[23,16,9]
        self.assertEqual(3,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        self.assertEqual(expected4,d.getValues());
        #
        d=DataArrayInt.Range(23,22,-7);
        self.assertEqual(1,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        self.assertEqual(23,d.getIJ(0,0));
        #
        d=DataArrayInt.Range(22,23,7);
        self.assertEqual(1,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        self.assertEqual(22,d.getIJ(0,0));
        #
        d=DataArrayInt.Range(22,22,7);
        self.assertEqual(0,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        #
        d=DataArrayInt.Range(22,22,-7);
        self.assertEqual(0,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        #
        self.assertRaises(InterpKernelException,DataArrayInt.Range,22,23,-7);
        self.assertRaises(InterpKernelException,DataArrayInt.Range,23,22,7);
        self.assertRaises(InterpKernelException,DataArrayInt.Range,23,22,0);
        self.assertRaises(InterpKernelException,DataArrayInt.Range,22,23,0);
        pass

    def testSwigUMeshGetItem1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        subMesh=m.buildPartOfMySelf([1,3],True);
        self.assertTrue(isinstance(subMesh,MEDCouplingUMesh))
        m1=m[[1,3]]
        self.assertTrue(isinstance(m1,MEDCouplingUMesh))
        m2=m[(1,3)]
        self.assertTrue(isinstance(m2,MEDCouplingUMesh))
        m3=m[1::2]
        self.assertTrue(isinstance(m3,MEDCouplingUMesh))
        m4=m[DataArrayInt.New([1,3])]
        m5_1=m[1]
        self.assertTrue(isinstance(m5_1,MEDCouplingUMesh))
        m5_2=m[3]
        self.assertTrue(isinstance(m5_2,MEDCouplingUMesh))
        m5=MEDCouplingUMesh.MergeUMeshesOnSameCoords([m5_1,m5_2]);
        m5.setName(subMesh.getName())
        self.assertTrue(isinstance(m4,MEDCouplingUMesh))
        self.assertTrue(subMesh.isEqual(m1,1e-12))
        self.assertTrue(subMesh.isEqual(m2,1e-12))
        self.assertTrue(subMesh.isEqual(m3,1e-12))
        self.assertTrue(subMesh.isEqual(m4,1e-12))
        self.assertTrue(subMesh.isEqual(m5,1e-12))
        self.assertRaises(InterpKernelException,m.buildPartOfMySelf,[1,5],True);
        pass
    
    def testSwigGetItem3(self):
        da=DataArrayInt.New([4,5,6])
        self.assertEqual(5,da[1])
        self.assertEqual(6,da[-1])
        self.assertRaises(InterpKernelException,da.__getitem__,3)
        da=DataArrayInt.New([4,5,6,7,8,9],2,3)
        self.assertEqual(9,da[1,2])
        da=DataArrayDouble.New([4.1,5.2,6.3])
        self.assertAlmostEqual(5.2,da[1],12)
        self.assertAlmostEqual(6.3,da[-1],12)
        self.assertRaises(InterpKernelException,da.__getitem__,3)
        da=DataArrayDouble.New([4.12,5.12,6.12,7.12,8.12,9.12],2,3)
        self.assertAlmostEqual(9.12,da[1,2],12)
        pass

    def testSwigDADISub1(self):
        mesh3D,mesh2D=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        bary=mesh3D.computeCellCenterOfMass()
        bary=bary[:,:2]
        pts=bary.getDifferentValues(1e-12)
        expected=[[0,6,12],[1,7,13],[2,8,14],[3,9,15],[4,10,16],[5,11,17]]
        for pos,pt in enumerate(pts):
            bary2=bary[:,:2]
            bary2[:]-=pt
            norm=bary2.magnitude()
            self.assertEqual(expected[pos],norm.findIdsInRange(-1.,1e-5).getValues())
            pass
        expected2=[[3.,54.],[-141.,180.],[21.,54.],[39.,72.],[-15.,90.],[21.,90.]]
        for pos,pt in enumerate(pts):
            bary2=bary[:,:2]
            bary2[:]+=pt
            self.assertAlmostEqual(expected2[pos][0],bary2.accumulate()[0],12);
            self.assertAlmostEqual(expected2[pos][1],bary2.accumulate()[1],12);
            pass
        expected3=[[-3.,22.5],[45.,337.5],[-9., 22.5],[-15.,67.5],[3.,112.5],[-9.,112.5]]
        for pos,pt in enumerate(pts):
            bary2=bary[:,:2]
            bary2[:]*=pt
            self.assertAlmostEqual(expected3[pos][0],bary2.accumulate()[0],12);
            self.assertAlmostEqual(expected3[pos][1],bary2.accumulate()[1],12);
            pass
        expected4=[[-12.,90.],[0.8,6.],[-4,90.],[-2.4,30.],[12.,18],[-4,18.]]
        for pos,pt in enumerate(pts):
            bary2=bary[:,:2]
            bary2[:]/=pt
            self.assertAlmostEqual(expected4[pos][0],bary2.accumulate()[0],12);
            self.assertAlmostEqual(expected4[pos][1],bary2.accumulate()[1],12);
            pass
        #
        d=DataArrayInt.New([1,2,0,1,0,2],3,2)
        e=DataArrayInt.New([1,11,101,2,12,102,3,13,103,4,14,104],4,3)
        expected5=[[1,11,101,77,77,77,77,77,77,4,14,104],[77,77,77,77,77,77,3,13,103,4,14,104],[77,77,77,2,12,102,77,77,77,4,14,104]]
        expected6=[[1,77,77,2,77,77,3,77,77,4,77,77],[77,77,101,77,77,102,77,77,103,77,77,104],[77,11,77,77,12,77,77,13,77,77,14,77]]
        for pos,tup in enumerate(d):
            f=e[:]
            self.assertTrue(isinstance(f,DataArrayInt))
            f[tup]=77
            self.assertEqual(expected5[pos],f.getValues())
            self.assertEqual(6*[77],f[tup].getValues())
            f=e[:]
            f[:,tup]=77
            self.assertEqual(expected6[pos],f.getValues())
            self.assertEqual(8*[77],f[:,tup].getValues())
            pass
        #
        e=e.convertToDblArr()
        for pos,tup in enumerate(d):
            f=e[:]
            self.assertTrue(isinstance(f,DataArrayDouble))
            f[tup]=77.
            self.assertEqual(expected5[pos],f.convertToIntArr().getValues())
            self.assertEqual(6*[77],f[tup].convertToIntArr().getValues())
            f=e[:]
            f[:,tup]=77.
            self.assertEqual(expected6[pos],f.convertToIntArr().getValues())
            self.assertEqual(8*[77],f[:,tup].convertToIntArr().getValues())
            pass
        pass

    def testDataArrayDoubleGetMinMaxPerComponent1(self):
        values1=[1.,2.,3.,-0.9,2.1,3.,1.3,1.7,3.,1.,1.8,3.]
        d1=DataArrayDouble.New();
        self.assertRaises(InterpKernelException,d1.getMinMaxPerComponent)
        d1=DataArrayDouble.New(values1,4,3);
        res=d1.getMinMaxPerComponent();
        self.assertTrue(isinstance(res,list))
        self.assertEqual(3,len(res))
        for i in range(3):
            self.assertTrue(isinstance(res[i],tuple))
            self.assertEqual(2,len(res[i]))
            pass
        expected1=[-0.9,1.3,1.7,2.1,3.,3.]
        for i in range(6):
            self.assertAlmostEqual(expected1[i], res[i // 2][i % 2], 14)
            pass
        #
        d1.rearrange(2);
        res=d1.getMinMaxPerComponent();
        self.assertTrue(isinstance(res,list))
        self.assertEqual(2,len(res))
        for i in range(2):
            self.assertTrue(isinstance(res[i],tuple))
            self.assertEqual(2,len(res[i]))
            pass
        expected2=[1.,3.,-0.9,3.]
        for i in range(4):
            self.assertAlmostEqual(expected2[i], res[i // 2][i % 2], 14)
            pass
        #
        d1.rearrange(1);
        res=d1.getMinMaxPerComponent();
        self.assertTrue(isinstance(res,list))
        self.assertEqual(1,len(res))
        for i in range(1):
            self.assertTrue(isinstance(res[i],tuple))
            self.assertEqual(2,len(res[i]))
            pass
        expected3=[-0.9,3.]
        for i in range(2):
            self.assertAlmostEqual(expected3[i], res[i // 2][i % 2], 14)
            pass
        pass

    def testDataArrayIntGetHashCode1(self):
        d1 = DataArrayInt.New(list(range(3545)))
        d2 = DataArrayInt.New(list(range(3545)))
        self.assertEqual(d2.getHashCode(),d1.getHashCode())
        self.assertEqual(232341068,d1.getHashCode())
        d1[886]=6
        self.assertEqual(232340188,d1.getHashCode())
        pass

    def testZipConnectivityPol1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        cells1=[2,3,4]
        m2_1=m1.buildPartOfMySelf(cells1,True);
        m2=m2_1
        self.assertTrue(isinstance(m2,MEDCouplingUMesh))
        # no permutation policy 0
        isOk,arr=m1.areCellsIncludedIn(m2,0)
        self.assertTrue(isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells1,arr.getValues())
        # no permutation policy 1
        isOk,arr=m1.areCellsIncludedIn(m2,1)
        self.assertTrue(isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells1,arr.getValues())
        # no permutation policy 2
        isOk,arr=m1.areCellsIncludedIn(m2,2)
        self.assertTrue(isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells1,arr.getValues())
        # some modification into m2
        modif1=[2,4,5]
        m2.getNodalConnectivity()[1:4]=modif1
        #policy 0 fails because cell0 in m2 has same orientation be not same connectivity
        expected1=[5,3,4]
        isOk,arr=m1.areCellsIncludedIn(m2,0)
        self.assertTrue(not isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(expected1,arr.getValues())
        #policy 1 succeeds because cell0 in m2 has not exactly the same conn
        isOk,arr=m1.areCellsIncludedIn(m2,1)
        self.assertTrue(isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells1,arr.getValues())
        #policy 2 succeeds because cell0 in m2 has same nodes in connectivity
        isOk,arr=m1.areCellsIncludedIn(m2,2)
        self.assertTrue(isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells1,arr.getValues())
        #some new modification into m2
        modif2=[2,5,4]
        m2.getNodalConnectivity()[1:4]=modif2
        #policy 0 fails because cell0 in m2 has not exactly the same conn
        isOk,arr=m1.areCellsIncludedIn(m2,0)
        self.assertTrue(not isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(expected1,arr.getValues())
        #policy 1 fails too because cell0 in m2 has not same orientation
        isOk,arr=m1.areCellsIncludedIn(m2,1)
        self.assertTrue(not isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(expected1,arr.getValues())
        #policy 2 succeeds because cell0 in m2 has same nodes in connectivity
        isOk,arr=m1.areCellsIncludedIn(m2,2)
        self.assertTrue(isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells1,arr.getValues())
        # Now 1D
        cells2=[3,2]
        m1=MEDCouplingDataForTest.build1DSourceMesh_2();
        m2_1=m1.buildPartOfMySelf(cells2,True);
        m2=m2_1
        self.assertTrue(isinstance(m2,MEDCouplingUMesh))
        # no permutation policy 0
        isOk,arr=m1.areCellsIncludedIn(m2,0)
        self.assertTrue(isOk);
        self.assertEqual(2,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells2,arr.getValues())
        # no permutation policy 1
        isOk,arr=m1.areCellsIncludedIn(m2,1)
        self.assertTrue(isOk);
        self.assertEqual(2,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells2,arr.getValues())
        # no permutation policy 2
        isOk,arr=m1.areCellsIncludedIn(m2,2)
        self.assertTrue(isOk);
        self.assertEqual(2,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells2,arr.getValues())
        # some modification into m2
        modif3=[4,3]
        m2.getNodalConnectivity()[1:3]=modif3
        #policy 0 fails because cell0 in m2 has not exactly the same conn
        expected2=[4,2]
        isOk,arr=m1.areCellsIncludedIn(m2,0)
        self.assertTrue(not isOk);
        self.assertEqual(2,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(expected2,arr.getValues())
        #policy 1 fails too because cell0 in m2 has not same orientation
        isOk,arr=m1.areCellsIncludedIn(m2,1)
        self.assertTrue(not isOk);
        self.assertEqual(2,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(expected2,arr.getValues())
        #policy 2 succeeds because cell0 in m2 has same nodes in connectivity
        isOk,arr=m1.areCellsIncludedIn(m2,2)
        self.assertTrue(isOk);
        self.assertEqual(2,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells2,arr.getValues())
        pass

    def toSeeIfDaIIopsAreOK(self,d):
        d+=5
        d*=6
        d/=3
        d-=2
        d%=7
        pass
        
    def testSwigDAIOp5(self):
        d=DataArrayInt.New([4,5,6,10,3,-1],2,3)
        self.toSeeIfDaIIopsAreOK(d)
        dExp=DataArrayInt.New([2,4,6,0,0,6],2,3)
        self.assertTrue(d.isEqual(dExp));
        pass
    
    def toSeeIfDaDIopsAreOK(self,d):
        d+=5
        d*=6
        d/=3
        d-=2
        pass

    def testSwigDADOp7(self):
        d=DataArrayDouble.New([4.,5.,6.,10.,3.,-1.],2,3)
        self.toSeeIfDaDIopsAreOK(d)
        dExp=DataArrayDouble.New([16.,18.,20.,28.,14.,6.],2,3)
        self.assertTrue(d.isEqual(dExp,1e-14));
        pass

    def testConvexEnvelop2D1(self):
        coords=[7.54758495819e-14,-1.12270326253e-12,8.43143594193,-1.02835845055e-12,4.21571797096,7.30183771609,-4.21571797097,7.30183771609,-8.43143594193,-1.09439981894e-12,-4.21571797097,-7.30183771609,4.21571797097,-7.30183771609,16.8628718839,-1.02835845055e-12,12.6471539129,7.30183771609,8.43143594193,14.6036754322,2.26427548746e-13,14.6036754322,-8.43143594193,14.6036754322,-12.6471539129,7.30183771609,-16.8628718839,-1.39630321727e-12,-12.6471539129,-7.30183771609,-8.43143594193,-14.6036754322,3.7737924791e-14,-14.6036754322,8.43143594193,-14.6036754322,12.6471539129,-7.30183771609,25.2943078258,-1.07553085654e-12,21.0785898548,7.30183771609,16.8628718839,14.6036754322,12.6471539129,21.9055131483,4.21571797096,21.9055131483,-4.21571797097,21.9055131483,-12.6471539129,21.9055131483,-16.8628718839,14.6036754322,-21.0785898548,7.30183771609,-25.2943078258,-1.02835845055e-12,-21.0785898548,-7.30183771609,-16.8628718839,-14.6036754322,-12.6471539129,-21.9055131483,-4.21571797097,-21.9055131483,4.21571797097,-21.9055131483,12.6471539129,-21.9055131483,16.8628718839,-14.6036754322,21.0785898548,-7.30183771609,33.7257437677,-7.45324014622e-13,29.5100257968,7.30183771609,25.2943078258,14.6036754322,21.0785898548,21.9055131483,16.8628718839,29.2073508644,8.43143594193,29.2073508644,-1.20761359331e-12,29.2073508644,-8.43143594193,29.2073508644,-16.8628718839,29.2073508644,-21.0785898548,21.9055131483,-25.2943078258,14.6036754322,-29.5100257968,7.30183771609,-33.7257437677,-7.26455052226e-13,-29.5100257968,-7.30183771609,-25.2943078258,-14.6036754322,-21.0785898548,-21.9055131483,-16.8628718839,-29.2073508644,-8.43143594193,-29.2073508644,4.15117172701e-13,-29.2073508644,8.43143594193,-29.2073508644,16.8628718839,-29.2073508644,21.0785898548,-21.9055131483,25.2943078258,-14.6036754322,29.5100257968,-7.30183771609,42.1571797097,-1.86802727715e-12,37.9414617387,7.30183771609,33.7257437677,14.6036754322,29.5100257968,21.9055131483,25.2943078258,29.2073508644,21.0785898548,36.5091885805,12.6471539129,36.5091885805,4.21571797096,36.5091885805,-4.21571797096,36.5091885805,-12.6471539129,36.5091885805,-21.0785898548,36.5091885805,-25.2943078258,29.2073508644,-29.5100257968,21.9055131483,-33.7257437677,14.6036754322,-37.9414617387,7.30183771609,-42.1571797097,-9.81186044565e-13,-37.9414617387,-7.30183771609,-33.7257437677,-14.6036754322,-29.5100257968,-21.9055131483,-25.2943078258,-29.2073508644,-21.0785898548,-36.5091885805,-12.6471539129,-36.5091885805,-4.21571797097,-36.5091885805,4.21571797097,-36.5091885805,12.6471539129,-36.5091885805,21.0785898548,-36.5091885805,25.2943078258,-29.2073508644,29.5100257968,-21.9055131483,33.7257437677,-14.6036754322,37.9414617387,-7.30183771609,50.5886156516,-6.98151608633e-13,46.3728976806,7.30183771609,42.1571797097,14.6036754322,37.9414617387,21.9055131483,33.7257437677,29.2073508644,29.5100257968,36.5091885805,25.2943078258,43.8110262966,16.8628718839,43.8110262966,8.43143594193,43.8110262966,-1.84915831476e-12,43.8110262966,-8.43143594193,43.8110262966,-16.8628718839,43.8110262966,-25.2943078258,43.8110262966,-29.5100257968,36.5091885805,-33.7257437677,29.2073508644,-37.9414617387,21.9055131483,-42.1571797097,14.6036754322,-46.3728976806,7.30183771609,-50.5886156516,-1.47177906685e-12,-46.3728976806,-7.30183771609,-42.1571797097,-14.6036754322,-37.9414617387,-21.9055131483,-33.7257437677,-29.2073508644,-29.5100257968,-36.5091885805,-25.2943078258,-43.8110262966,-16.8628718839,-43.8110262966,-8.43143594193,-43.8110262966,7.54758495819e-14,-43.8110262966,8.43143594193,-43.8110262966,16.8628718839,-43.8110262966,25.2943078258,-43.8110262966,29.5100257968,-36.5091885805,33.7257437677,-29.2073508644,37.9414617387,-21.9055131483,42.1571797097,-14.6036754322,46.3728976806,-7.30183771609,59.0200515935,-7.9249642061e-13,54.8043336225,7.30183771609,50.5886156516,14.6036754322,46.3728976806,21.9055131483,42.1571797097,29.2073508644,37.9414617387,36.5091885805,33.7257437677,43.8110262966,29.5100257968,51.1128640127,21.0785898548,51.1128640127,12.6471539129,51.1128640127,4.21571797096,51.1128640127,-4.21571797096,51.1128640127,-12.6471539129,51.1128640127,-21.0785898548,51.1128640127,-29.5100257968,51.1128640127,-33.7257437677,43.8110262966,-37.9414617387,36.5091885805,-42.1571797097,29.2073508644,-46.3728976806,21.9055131483,-50.5886156516,14.6036754322,-54.8043336226,7.30183771609,-59.0200515935,-1.31139288649e-12,-54.8043336226,-7.30183771609,-50.5886156516,-14.6036754322,-46.3728976806,-21.9055131483,-42.1571797097,-29.2073508644,-37.9414617387,-36.5091885805,-33.7257437677,-43.8110262966,-29.5100257968,-51.1128640127,-21.0785898548,-51.1128640127,-12.6471539129,-51.1128640127,-4.21571797097,-51.1128640127,4.21571797097,-51.1128640127,12.6471539129,-51.1128640127,21.0785898548,-51.1128640127,29.5100257968,-51.1128640127,33.7257437677,-43.8110262966,37.9414617387,-36.5091885805,42.1571797097,-29.2073508644,46.3728976806,-21.9055131483,50.5886156516,-14.6036754322,54.8043336225,-7.30183771609,67.4514875354,-2.14162723189e-12,63.2357695645,7.30183771609,59.0200515935,14.6036754322,54.8043336226,21.9055131483,50.5886156516,29.2073508644,46.3728976806,36.5091885805,42.1571797097,43.8110262966,37.9414617387,51.1128640127,33.7257437677,58.4147017287,25.2943078258,58.4147017287,16.8628718839,58.4147017287,8.43143594193,58.4147017287,6.79282646237e-13,58.4147017287,-8.43143594193,58.4147017287,-16.8628718839,58.4147017287,-25.2943078258,58.4147017287,-33.7257437677,58.4147017287,-37.9414617387,51.1128640127,-42.1571797097,43.8110262966,-46.3728976806,36.5091885805,-50.5886156516,29.2073508644,-54.8043336226,21.9055131483,-59.0200515935,14.6036754322,-63.2357695645,7.30183771609,-67.4514875354,-1.16044118732e-12,-63.2357695645,-7.30183771609,-59.0200515935,-14.6036754322,-54.8043336226,-21.9055131483,-50.5886156516,-29.2073508644,-46.3728976806,-36.5091885805,-42.1571797097,-43.8110262966,-37.9414617387,-51.1128640127,-33.7257437677,-58.4147017287,-25.2943078258,-58.4147017287,-16.8628718839,-58.4147017287,-8.43143594193,-58.4147017287,-5.66068871864e-14,-58.4147017287,8.43143594193,-58.4147017287,16.8628718839,-58.4147017287,25.2943078258,-58.4147017287,33.7257437677,-58.4147017287,37.9414617387,-51.1128640127,42.1571797097,-43.8110262966,46.3728976806,-36.5091885805,50.5886156516,-29.2073508644,54.8043336226,-21.9055131483,59.0200515935,-14.6036754322,63.2357695645,-7.30183771609,75.8829234774,-2.29257893105e-12,71.6672055064,7.30183771609,67.4514875354,14.6036754322,63.2357695645,21.9055131483,59.0200515935,29.2073508644,54.8043336226,36.5091885805,50.5886156516,43.8110262966,46.3728976806,51.1128640127,42.1571797097,58.4147017287,37.9414617387,65.7165394448,29.5100257968,65.7165394448,21.0785898548,65.7165394448,12.6471539129,65.7165394448,4.21571797097,65.7165394448,-4.21571797096,65.7165394448,-12.6471539129,65.7165394448,-21.0785898548,65.7165394448,-29.5100257968,65.7165394448,-37.9414617387,65.7165394448,-42.1571797097,58.4147017287,-46.3728976806,51.1128640127,-50.5886156516,43.8110262966,-54.8043336226,36.5091885805,-59.0200515935,29.2073508644,-63.2357695645,21.9055131483,-67.4514875354,14.6036754322,-71.6672055064,7.30183771609,-75.8829234774,-1.31139288649e-12,-71.6672055064,-7.30183771609,-67.4514875354,-14.6036754322,-63.2357695645,-21.9055131483,-59.0200515935,-29.2073508644,-54.8043336226,-36.5091885805,-50.5886156516,-43.8110262966,-46.3728976806,-51.1128640127,-42.1571797097,-58.4147017287,-37.9414617387,-65.7165394448,-29.5100257968,-65.7165394448,-21.0785898548,-65.7165394448,-12.6471539129,-65.7165394448,-4.21571797097,-65.7165394448,4.21571797097,-65.7165394448,12.6471539129,-65.7165394448,21.0785898548,-65.7165394448,29.5100257968,-65.7165394448,37.9414617387,-65.7165394448,42.1571797097,-58.4147017287,46.3728976806,-51.1128640127,50.5886156516,-43.8110262966,54.8043336226,-36.5091885805,59.0200515935,-29.2073508644,63.2357695645,-21.9055131483,67.4514875354,-14.6036754322,71.6672055064,-7.30183771609,84.3143594193,-1.49064802924e-12,80.0986414483,7.30183771609,75.8829234774,14.6036754322,71.6672055064,21.9055131483,67.4514875354,29.2073508644,63.2357695645,36.5091885805,59.0200515935,43.8110262966,54.8043336226,51.1128640127,50.5886156516,58.4147017287,46.3728976806,65.7165394448,42.1571797097,73.0183771609,33.7257437677,73.0183771609,25.2943078258,73.0183771609,16.8628718839,73.0183771609,8.43143594193,73.0183771609,2.0755858635e-12,73.0183771609,-8.43143594193,73.0183771609,-16.8628718839,73.0183771609,-25.2943078258,73.0183771609,-33.7257437677,73.0183771609,-42.1571797097,73.0183771609,-46.3728976806,65.7165394448,-50.5886156516,58.4147017287,-54.8043336226,51.1128640127,-59.0200515935,43.8110262966,-63.2357695645,36.5091885805,-67.4514875354,29.2073508644,-71.6672055064,21.9055131483,-75.8829234774,14.6036754322,-80.0986414483,7.30183771609,-84.3143594193,-1.11326878133e-12,-80.0986414483,-7.30183771609,-75.8829234774,-14.6036754322,-71.6672055064,-21.9055131483,-67.4514875354,-29.2073508644,-63.2357695645,-36.5091885805,-59.0200515935,-43.8110262966,-54.8043336226,-51.1128640127,-50.5886156516,-58.4147017287,-46.3728976806,-65.7165394448,-42.1571797097,-73.0183771609,-33.7257437677,-73.0183771609,-25.2943078258,-73.0183771609,-16.8628718839,-73.0183771609,-8.43143594193,-73.0183771609,-5.66068871864e-14,-73.0183771609,8.43143594193,-73.0183771609,16.8628718839,-73.0183771609,25.2943078258,-73.0183771609,33.7257437677,-73.0183771609,42.1571797097,-73.0183771609,46.3728976806,-65.7165394448,50.5886156516,-58.4147017287,54.8043336226,-51.1128640127,59.0200515935,-43.8110262966,63.2357695645,-36.5091885805,67.4514875354,-29.2073508644,71.6672055064,-21.9055131483,75.8829234774,-14.6036754322,80.0986414483,-7.3018377161]
        conn=[0,2,3,4,5,6,1,1,8,2,0,6,18,7,2,9,10,3,0,1,8,3,10,11,12,4,0,2,4,3,12,13,14,5,0,5,0,4,14,15,16,6,6,1,0,5,16,17,18,7,20,8,1,18,36,19,8,21,9,2,1,7,20,9,22,23,10,2,8,21,10,23,24,11,3,2,9,11,24,25,26,12,3,10,12,11,26,27,13,4,3,13,12,27,28,29,14,4,14,4,13,29,30,15,5,15,5,14,30,31,32,16,16,6,5,15,32,33,17,17,18,6,16,33,34,35,18,7,1,6,17,35,36,19,38,20,7,36,60,37,20,39,21,8,7,19,38,21,40,22,9,8,20,39,22,41,42,23,9,21,40,23,42,43,24,10,9,22,24,43,44,25,11,10,23,25,44,45,46,26,11,24,26,25,46,47,27,12,11,27,26,47,48,28,13,12,28,27,48,49,50,29,13,29,13,28,50,51,30,14,30,14,29,51,52,31,15,31,15,30,52,53,54,32,32,16,15,31,54,55,33,33,17,16,32,55,56,34,34,35,17,33,56,57,58,35,36,18,17,34,58,59,36,19,7,18,35,59,60,37,62,38,19,60,90,61,38,63,39,20,19,37,62,39,64,40,21,20,38,63,40,65,41,22,21,39,64,41,66,67,42,22,40,65,42,67,68,43,23,22,41,43,68,69,44,24,23,42,44,69,70,45,25,24,43,45,70,71,72,46,25,44,46,45,72,73,47,26,25,47,46,73,74,48,27,26,48,47,74,75,49,28,27,49,48,75,76,77,50,28,50,28,49,77,78,51,29,51,29,50,78,79,52,30,52,30,51,79,80,53,31,53,31,52,80,81,82,54,54,32,31,53,82,83,55,55,33,32,54,83,84,56,56,34,33,55,84,85,57,57,58,34,56,85,86,87,58,59,35,34,57,87,88,59,60,36,35,58,88,89,60,37,19,36,59,89,90,61,92,62,37,90,126,91,62,93,63,38,37,61,92,63,94,64,39,38,62,93,64,95,65,40,39,63,94,65,96,66,41,40,64,95,66,97,98,67,41,65,96,67,98,99,68,42,41,66,68,99,100,69,43,42,67,69,100,101,70,44,43,68,70,101,102,71,45,44,69,71,102,103,104,72,45,70,72,71,104,105,73,46,45,73,72,105,106,74,47,46,74,73,106,107,75,48,47,75,74,107,108,76,49,48,76,75,108,109,110,77,49,77,49,76,110,111,78,50,78,50,77,111,112,79,51,79,51,78,112,113,80,52,80,52,79,113,114,81,53,81,53,80,114,115,116,82,82,54,53,81,116,117,83,83,55,54,82,117,118,84,84,56,55,83,118,119,85,85,57,56,84,119,120,86,86,87,57,85,120,121,122,87,88,58,57,86,122,123,88,89,59,58,87,123,124,89,90,60,59,88,124,125,90,61,37,60,89,125,126,91,128,92,61,126,168,127,92,129,93,62,61,91,128,93,130,94,63,62,92,129,94,131,95,64,63,93,130,95,132,96,65,64,94,131,96,133,97,66,65,95,132,97,134,135,98,66,96,133,98,135,136,99,67,66,97,99,136,137,100,68,67,98,100,137,138,101,69,68,99,101,138,139,102,70,69,100,102,139,140,103,71,70,101,103,140,141,142,104,71,102,104,103,142,143,105,72,71,105,104,143,144,106,73,72,106,105,144,145,107,74,73,107,106,145,146,108,75,74,108,107,146,147,109,76,75,109,108,147,148,149,110,76,110,76,109,149,150,111,77,111,77,110,150,151,112,78,112,78,111,151,152,113,79,113,79,112,152,153,114,80,114,80,113,153,154,115,81,115,81,114,154,155,156,116,116,82,81,115,156,157,117,117,83,82,116,157,158,118,118,84,83,117,158,159,119,119,85,84,118,159,160,120,120,86,85,119,160,161,121,121,122,86,120,161,162,163,122,123,87,86,121,163,164,123,124,88,87,122,164,165,124,125,89,88,123,165,166,125,126,90,89,124,166,167,126,91,61,90,125,167,168,127,170,128,91,168,216,169,128,171,129,92,91,127,170,129,172,130,93,92,128,171,130,173,131,94,93,129,172,131,174,132,95,94,130,173,132,175,133,96,95,131,174,133,176,134,97,96,132,175,134,177,178,135,97,133,176,135,178,179,136,98,97,134,136,179,180,137,99,98,135,137,180,181,138,100,99,136,138,181,182,139,101,100,137,139,182,183,140,102,101,138,140,183,184,141,103,102,139,141,184,185,186,142,103,140,142,141,186,187,143,104,103,143,142,187,188,144,105,104,144,143,188,189,145,106,105,145,144,189,190,146,107,106,146,145,190,191,147,108,107,147,146,191,192,148,109,108,148,147,192,193,194,149,109,149,109,148,194,195,150,110,150,110,149,195,196,151,111,151,111,150,196,197,152,112,152,112,151,197,198,153,113,153,113,152,198,199,154,114,154,114,153,199,200,155,115,155,115,154,200,201,202,156,156,116,115,155,202,203,157,157,117,116,156,203,204,158,158,118,117,157,204,205,159,159,119,118,158,205,206,160,160,120,119,159,206,207,161,161,121,120,160,207,208,162,162,163,121,161,208,209,210,163,164,122,121,162,210,211,164,165,123,122,163,211,212,165,166,124,123,164,212,213,166,167,125,124,165,213,214,167,168,126,125,166,214,215,168,127,91,126,167,215,216,169,218,170,127,216,270,217,170,219,171,128,127,169,218,171,220,172,129,128,170,219,172,221,173,130,129,171,220,173,222,174,131,130,172,221,174,223,175,132,131,173,222,175,224,176,133,132,174,223,176,225,177,134,133,175,224,177,226,227,178,134,176,225,178,227,228,179,135,134,177,179,228,229,180,136,135,178,180,229,230,181,137,136,179,181,230,231,182,138,137,180,182,231,232,183,139,138,181,183,232,233,184,140,139,182,184,233,234,185,141,140,183,185,234,235,236,186,141,184,186,185,236,237,187,142,141,187,186,237,238,188,143,142,188,187,238,239,189,144,143,189,188,239,240,190,145,144,190,189,240,241,191,146,145,191,190,241,242,192,147,146,192,191,242,243,193,148,147,193,192,243,244,245,194,148,194,148,193,245,246,195,149,195,149,194,246,247,196,150,196,150,195,247,248,197,151,197,151,196,248,249,198,152,198,152,197,249,250,199,153,199,153,198,250,251,200,154,200,154,199,251,252,201,155,201,155,200,252,253,254,202,202,156,155,201,254,255,203,203,157,156,202,255,256,204,204,158,157,203,256,257,205,205,159,158,204,257,258,206,206,160,159,205,258,259,207,207,161,160,206,259,260,208,208,162,161,207,260,261,209,209,210,162,208,261,262,263,210,211,163,162,209,263,264,211,212,164,163,210,264,265,212,213,165,164,211,265,266,213,214,166,165,212,266,267,214,215,167,166,213,267,268,215,216,168,167,214,268,269,216,169,127,168,215,269,270,217,272,218,169,270,330,271,218,273,219,170,169,217,272,219,274,220,171,170,218,273,220,275,221,172,171,219,274,221,276,222,173,172,220,275,222,277,223,174,173,221,276,223,278,224,175,174,222,277,224,279,225,176,175,223,278,225,280,226,177,176,224,279,226,281,282,227,177,225,280,227,282,283,228,178,177,226,228,283,284,229,179,178,227,229,284,285,230,180,179,228,230,285,286,231,181,180,229,231,286,287,232,182,181,230,232,287,288,233,183,182,231,233,288,289,234,184,183,232,234,289,290,235,185,184,233,235,290,291,292,236,185,234,236,235,292,293,237,186,185,237,236,293,294,238,187,186,238,237,294,295,239,188,187,239,238,295,296,240,189,188,240,239,296,297,241,190,189,241,240,297,298,242,191,190,242,241,298,299,243,192,191,243,242,299,300,244,193,192,244,243,300,301,302,245,193,245,193,244,302,303,246,194,246,194,245,303,304,247,195,247,195,246,304,305,248,196,248,196,247,305,306,249,197,249,197,248,306,307,250,198,250,198,249,307,308,251,199,251,199,250,308,309,252,200,252,200,251,309,310,253,201,253,201,252,310,311,312,254,254,202,201,253,312,313,255,255,203,202,254,313,314,256,256,204,203,255,314,315,257,257,205,204,256,315,316,258,258,206,205,257,316,317,259,259,207,206,258,317,318,260,260,208,207,259,318,319,261,261,209,208,260,319,320,262,262,263,209,261,320,321,322,263,264,210,209,262,322,323,264,265,211,210,263,323,324,265,266,212,211,264,324,325,266,267,213,212,265,325,326,267,268,214,213,266,326,327,268,269,215,214,267,327,328,269,270,216,215,268,328,329,270,217,169,216,269,329,330,271,272,217,330,273,218,217,271,274,219,218,272,275,220,219,273,276,221,220,274,277,222,221,275,278,223,222,276,279,224,223,277,280,225,224,278,281,226,225,279,281,282,226,280,283,227,226,281,284,228,227,282,285,229,228,283,286,230,229,284,287,231,230,285,288,232,231,286,289,233,232,287,290,234,233,288,291,235,234,289,291,292,235,290,291,293,236,235,292,294,237,236,293,295,238,237,294,296,239,238,295,297,240,239,296,298,241,240,297,299,242,241,298,300,243,242,299,301,244,243,301,300,302,244,244,301,303,245,245,302,304,246,246,303,305,247,247,304,306,248,248,305,307,249,249,306,308,250,250,307,309,251,251,308,310,252,252,309,311,253,311,253,310,312,254,253,311,313,255,254,312,314,256,255,313,315,257,256,314,316,258,257,315,317,259,258,316,318,260,259,317,319,261,260,318,320,262,261,319,321,321,322,262,320,323,263,262,321,324,264,263,322,325,265,264,323,326,266,265,324,327,267,266,325,328,268,267,326,329,269,268,327,330,270,269,328,271,217,270,329]
        connI=[0,7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126,133,140,147,154,161,168,175,182,189,196,203,210,217,224,231,238,245,252,259,266,273,280,287,294,301,308,315,322,329,336,343,350,357,364,371,378,385,392,399,406,413,420,427,434,441,448,455,462,469,476,483,490,497,504,511,518,525,532,539,546,553,560,567,574,581,588,595,602,609,616,623,630,637,644,651,658,665,672,679,686,693,700,707,714,721,728,735,742,749,756,763,770,777,784,791,798,805,812,819,826,833,840,847,854,861,868,875,882,889,896,903,910,917,924,931,938,945,952,959,966,973,980,987,994,1001,1008,1015,1022,1029,1036,1043,1050,1057,1064,1071,1078,1085,1092,1099,1106,1113,1120,1127,1134,1141,1148,1155,1162,1169,1176,1183,1190,1197,1204,1211,1218,1225,1232,1239,1246,1253,1260,1267,1274,1281,1288,1295,1302,1309,1316,1323,1330,1337,1344,1351,1358,1365,1372,1379,1386,1393,1400,1407,1414,1421,1428,1435,1442,1449,1456,1463,1470,1477,1484,1491,1498,1505,1512,1519,1526,1533,1540,1547,1554,1561,1568,1575,1582,1589,1596,1603,1610,1617,1624,1631,1638,1645,1652,1659,1666,1673,1680,1687,1694,1701,1708,1715,1722,1729,1736,1743,1750,1757,1764,1771,1778,1785,1792,1799,1806,1813,1820,1827,1834,1841,1848,1855,1862,1869,1876,1883,1890,1897,1901,1905,1909,1913,1917,1921,1925,1929,1933,1937,1941,1945,1949,1953,1957,1961,1965,1969,1973,1977,1981,1985,1989,1993,1997,2001,2005,2009,2013,2017,2021,2025,2029,2033,2037,2041,2045,2049,2053,2057,2061,2065,2069,2073,2077,2081,2085,2089,2093,2097,2101,2105,2109,2113,2117,2121,2125,2129,2133,2137]
        #
        m=MEDCouplingUMesh.New("convexhull",2);
        m.allocateCells(331);
        for i in range(331):
            m.insertNextCell(NORM_POLYGON,conn[connI[i]:connI[i+1]]);
            pass
        m.finishInsertingCells();
        coordsDa=DataArrayDouble.New(coords,331,2);
        m.setCoords(coordsDa);
        m.checkConsistencyLight();
        #
        da=m.convexEnvelop2D();
        m.checkConsistencyLight()
        self.assertEqual(coordsDa.getHiddenCppPointer(),m.getCoords().getHiddenCppPointer())
        daC=da.buildComplement(m.getNumberOfCells());
        expected2=DataArrayInt.New([271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,302,303,304,305,306,307,308,309,310,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330]);
        self.assertTrue(expected2.isEqual(daC));
        #
        vals=m.getMeasureField(False).getArray()
        ref=271*[184.69493088478035]+3*[-61.564976961404426,-92.34746544254946,-92.34746544259811,-92.34746544253488,-92.3474654425349,-92.34746544180479,-92.34746544253493,-92.3474654419026,-92.34746544190256,-92.34746544253491]+2*[61.564976961404426,-92.34746544254946,-92.34746544259811,-92.34746544253488,-92.3474654425349,-92.34746544180479,-92.34746544253493,-92.3474654419026,-92.34746544190256,-92.34746544253491]+[-61.564976961404426,-92.34746544254946,-92.34746544259811,-92.34746544253488,-92.3474654425349,-92.34746544180479,-92.34746544253493,-92.3474654419026,-92.34746544190256,-92.34746544253491]
        vals-=DataArrayDouble.New(ref)
        vals.abs()
        theTest=vals.findIdsInRange(-1.,1e-7)
        self.assertTrue(theTest.isIota(331))
        pass

    def testSwigDAIOp8(self):
        da=DataArrayInt.New([7,5,6,7,8,9,9,10,12,13,47,15])
        self.assertTrue(7 in da)
        self.assertTrue(47 in da)
        self.assertTrue(15 in da)
        self.assertEqual(0,da.index(7))
        self.assertEqual(10,da.index(47))
        self.assertTrue(14 not in da)
        self.assertEqual(5,da.findIdSequence([9,9]))
        self.assertEqual(-1,da.findIdSequence([5,8]))
        da.rearrange(2)
        self.assertTrue([47,16] not in da)
        self.assertTrue([5,6] not in da)
        self.assertTrue([6,7] in da)
        self.assertEqual(4,da.index([12,13]))
        pass

    def testDataArraySort1(self):
        arr=DataArrayInt.New();
        self.assertRaises(InterpKernelException,arr.sort,True)
        self.assertRaises(InterpKernelException,arr.sort,False)
        values=[2,1,6,5,4,7]
        arr.alloc(3,2);
        self.assertRaises(InterpKernelException,arr.sort,True)
        self.assertRaises(InterpKernelException,arr.sort,False)
        arr.rearrange(1);
        arr.setValues(values,6,1)
        arr1=arr.deepCopy();
        arr2=arr.deepCopy();
        arr1.sort(True);
        expected1=[1,2,4,5,6,7]
        self.assertEqual(6,arr1.getNumberOfTuples());
        self.assertEqual(1,arr1.getNumberOfComponents());
        self.assertEqual(expected1,arr1.getValues());
        arr2.sort(False);
        expected2=[7,6,5,4,2,1]
        self.assertEqual(6,arr2.getNumberOfTuples());
        self.assertEqual(1,arr2.getNumberOfComponents());
        self.assertTrue(expected2,arr2.getValues());
        #
        ard=DataArrayDouble.New();
        self.assertRaises(InterpKernelException,ard.sort,True)
        self.assertRaises(InterpKernelException,ard.sort,False)
        valuesD=[2.,1.,6.,5.,4.,7.]
        ard.alloc(3,2);
        self.assertRaises(InterpKernelException,ard.sort,True)
        self.assertRaises(InterpKernelException,ard.sort,False)
        ard.rearrange(1);
        ard.setValues(valuesD,6,1)
        ard1=ard.deepCopy();
        ard2=ard.deepCopy();
        ard1.sort(True);
        expected3=[1.,2.,4.,5.,6.,7.]
        self.assertEqual(6,ard1.getNumberOfTuples());
        self.assertEqual(1,ard1.getNumberOfComponents());
        for i in range(6):
            self.assertAlmostEqual(expected3[i],ard1.getIJ(i,0),12)
            pass
        ard2.sort(False);
        expected4=[7.,6.,5.,4.,2.,1.]
        self.assertEqual(6,ard2.getNumberOfTuples());
        self.assertEqual(1,ard2.getNumberOfComponents());
        for i in range(6):
            self.assertAlmostEqual(expected4[i],ard2.getIJ(i,0),12)
            pass
        pass
    
    def testPartitionBySpreadZone1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        m4=MEDCouplingUMesh.MergeUMeshes([m,m[-3:],m[0:2]]);
        m4.renumberCells([5,2,9,6,4,7,0,1,3,8]);
        #
        v2=m4.partitionBySpreadZone();
        self.assertTrue(3,len(v2));
        self.assertTrue(v2[0].isEqual(DataArrayInt.New([0,1,7])))
        self.assertTrue(v2[1].isEqual(DataArrayInt.New([2,4,5,6,9])))
        self.assertTrue(v2[2].isEqual(DataArrayInt.New([3,8])))
        #
        m5=m4.buildSpreadZonesWithPoly();
        self.assertEqual(3,m5.getNumberOfCells());
        self.assertTrue(m5.getCoords().getHiddenCppPointer()==m4.getCoords().getHiddenCppPointer());
        self.assertEqual([5,15,16,17,14,11,13,12,5,2,1,0,3,6,7,8,5,5,18,21,22,20,19],m5.getNodalConnectivity().getValues())
        self.assertEqual([0,8,17,23],m5.getNodalConnectivityIndex().getValues())
        #
        pass

    def testGiveCellsWithType1(self):
        expected0=[1,2]
        expected1=[0,3,4]
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        da=m.giveCellsWithType(NORM_TRI3);
        self.assertEqual(2,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        self.assertEqual(expected0,da.getValues())
        #
        da=m.giveCellsWithType(NORM_QUAD4);
        self.assertEqual(3,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        self.assertEqual(expected1,da.getValues())
        #
        da=m.giveCellsWithType(NORM_TRI6);
        self.assertEqual(0,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        #
        self.assertRaises(InterpKernelException,m.giveCellsWithType,NORM_SEG2)
        self.assertRaises(InterpKernelException,m.giveCellsWithType,NORM_HEXA8)
        pass

    def testSwigDAOp1(self):
        d=DataArrayDouble.New(5,2)
        d.rearrange(1) ; d.iota(2.) ; d.rearrange(2)
        d.setInfoOnComponents(["X [m]","Y [m]"])
        d.setName("AName")
        #
        d1=d+[8,9]
        self.assertTrue(d1.isEqualWithoutConsideringStr(DataArrayDouble.New([10.0,12.0,12.0,14.0,14.0,16.0,16.0,18.0,18.0,20.0]),1e-12))
        d1bis=DataArrayDouble.New([8,9],1,2)+d
        self.assertTrue(d1bis.isEqual(d1,1e-12))
        d1ter=[8,9]+d
        self.assertTrue(d1ter.isEqual(d1,1e-12))
        #
        d2=d1-[8,9]
        self.assertTrue(d2.isEqual(d,1e-12))
        self.assertRaises(InterpKernelException,d1.__rsub__,[8,9])#[8,9]-d1
        #
        d3=d*[8,9]
        self.assertTrue(d3.isEqualWithoutConsideringStr(DataArrayDouble.New([16.0,27.0,32.0,45.0,48.0,63.0,64.0,81.0,80.0,99.0]),1e-12))
        d3bis=DataArrayDouble.New([8,9],1,2)*d
        self.assertTrue(d3bis.isEqual(d3,1e-12))
        d3ter=[8,9]*d
        self.assertTrue(d3ter.isEqual(d3,1e-12))
        #
        d4=d3/[8,9]
        self.assertTrue(d4.isEqual(d,1e-12))
        #
        d=DataArrayInt.New(5,2)
        d.rearrange(1) ; d.iota(2) ; d.rearrange(2)
        d.setInfoOnComponents(["X [m]","Y [m]"])
        d.setName("AName")
        #
        d1=d+[8,9]
        self.assertEqual(d1.getValues(),[10,12,12,14,14,16,16,18,18,20])
        d1bis=DataArrayInt.New([8,9],1,2)+d
        self.assertTrue(d1bis.isEqual(d1))
        d1ter=[8,9]+d
        self.assertTrue(d1ter.isEqual(d1))
        #
        d2=d1-[8,9]
        self.assertTrue(d2.isEqual(d))
        self.assertRaises(InterpKernelException,d1.__rsub__,[8,9])
        #
        d3=d*[8,9]
        self.assertEqual(d3.getValues(),[16,27,32,45,48,63,64,81,80,99])
        d3bis=DataArrayInt.New([8,9],1,2)*d
        self.assertTrue(d3bis.isEqual(d3))
        d3ter=[8,9]*d
        self.assertTrue(d3ter.isEqual(d3))
        #
        d4=d3/[8,9]
        self.assertTrue(d4.isEqual(d))
        #
        d5=d%[4,5]
        self.assertEqual(d5.getValues(),[2,3,0,0,2,2,0,4,2,1])
        pass

    def testSwigSelectTupleId2DAIBug1(self):
        da=DataArrayInt.New([0,1,2,3,12,13,4,5,6,7,14,15,8,9,10,11,16,17])
        self.assertEqual([2,6,10],da[2::6].getValues())
        self.assertEqual([0,4,8],da[::6].getValues())
        self.assertEqual([5,9],da[7::6].getValues())
        self.assertEqual([5],da[7:-5:6].getValues())
        pass

    def testSwigCpp5Safe1(self):
        m=MEDCouplingUMesh.New("toto",2)
        coords=DataArrayDouble.New([0.,0.,1.,0.,1.,1.,0.,1.],4,2)
        m.setCoords(coords)
        vecs=DataArrayDouble.New([2.,3.,4.,5.,6.,7.],3,2)
        expected1=[[2.,3.,3.,3.,3.,4.,2.,4.0],[4.,5.,5.,5.,5.,6.,4.,6.0],[6.,7.,7.,7.,7.,8.,6.,8.0]]
        for pos,vec in enumerate(vecs):
            m2=m.deepCopy()
            m2.translate(vec)
            self.assertTrue(m2.getCoords().isEqual(DataArrayDouble.New(expected1[pos],4,2),1e-12))
            pass
        for pos,vec in enumerate(vecs):
            m2=m.deepCopy()
            m2.translate(vec.buildDADouble())
            self.assertTrue(m2.getCoords().isEqual(DataArrayDouble.New(expected1[pos],4,2),1e-12))
            pass
        pass

    def testSwigBugNonRegressionZipDA(self):
        angles = [pi / 3 * x for x in range(6)]
        radius=3
        #
        dad=DataArrayDouble.New(6, 2)
        dad[:,0]=radius
        dad[:,1]=angles
        #
        dad2=dad.fromPolarToCart()
        dads=[dad2.deepCopy() for elt in 7*[None]]
        #
        translationToPerform=[[0.01,0.02],[3./2.*radius,-radius*sqrt(3.)/2],[3./2.*radius,radius*sqrt(3.)/2],[0.,radius*sqrt(3.)],[-3./2.*radius,radius*sqrt(3.)/2],[-3./2.*radius,-radius*sqrt(3.)/2],[0.,-radius*sqrt(3.)]]
        for d,t in zip(dads,translationToPerform):
            d+=t
            pass
        for elt in dads:
            self.assertTrue(not dad2.isEqual(elt,1e-12))
            pass
        for d,t in zip(dads,translationToPerform):
            d-=t
            pass
        for elt in dads:
            self.assertTrue(dad2.isEqual(elt,1e-12))
            pass
        pass

    def testBuildSlice3D2(self):
        mesh3D,mesh2D=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        vec1=[-0.07,1.,0.07]
        origin1=[1.524,1.4552,1.74768]
        slice1,ids=mesh3D.buildSlice3D(origin1,vec1,1e-10);
        f=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME)
        f.setTime(4.5,6,7) ; f.setMesh(mesh3D)
        arr=DataArrayDouble(mesh3D.getNumberOfCells(),2)
        arr.rearrange(1) ; arr.iota(2.) ; arr.rearrange(2)
        f.setArray(arr)
        f.checkConsistencyLight()
        expected1=DataArrayInt([1,3,4,7,9,10,13,15,16])
        self.assertTrue(expected1.isEqual(ids))
        arr2=arr[expected1]
        #
        f2=f.extractSlice3D(origin1,vec1,1e-10)
        self.assertTrue(f2.getArray().isEqual(arr2,1e-12));
        self.assertTrue(slice1.isEqual(f2.getMesh(),1e-12))
        self.assertEqual(6,f2.getTime()[1]) ; self.assertEqual(7,f2.getTime()[2])
        self.assertAlmostEqual(4.5,f2.getTime()[0],12);
        pass

    def testComputeTupleIdsToSelectFromCellIds1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_3()
        f=MEDCouplingFieldDouble.New(ON_GAUSS_NE,NO_TIME);
        f.setMesh(m);
        arr=DataArrayDouble(52,2) ; arr.rearrange(1) ; arr.iota(7.) ; arr.rearrange(2)
        f.setArray(arr)
        #
        f2=f.buildSubPart([1,5,9])
        f2.checkConsistencyLight()
        cI=m.computeNbOfNodesPerCell()
        cI.computeOffsetsFull()
        sel=DataArrayInt([1,5,9])
        res=sel.buildExplicitArrByRanges(cI)
        arr2=arr[res]
        self.assertTrue(arr2.isEqual(DataArrayDouble([13,14,15,16,17,18,19,20,59,60,61,62,63,64,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110],15,2),1e-12))
        self.assertTrue(arr2.isEqual(f2.getArray(),1e-12))
        pass

    def testComputeSkin1(self):
        arrX=DataArrayDouble([2.,3.4,5.6,7.7,8.0]) ; arrY=DataArrayDouble([2.,3.4,5.6,7.7,9.0,14.2])
        cmesh=MEDCouplingCMesh() ; cmesh.setCoordsAt(0,arrX) ; cmesh.setCoordsAt(1,arrY)
        umesh=cmesh.buildUnstructured()
        #
        skin=umesh.computeSkin()
        self.assertEqual(18,skin.getNumberOfCells())
        self.assertEqual(1,skin.getMeshDimension())
        self.assertTrue(skin.getCoords().getHiddenCppPointer()==umesh.getCoords().getHiddenCppPointer())
        self.assertEqual([0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54],skin.getNodalConnectivityIndex().getValues())
        self.assertEqual([1,1,0,1,0,5,1,2,1,1,3,2,1,4,3,1,9,4,1,5,10,1,14,9,1,10,15,1,19,14,1,15,20,1,24,19,1,20,25,1,25,26,1,26,27,1,27,28,1,28,29,1,29,24],skin.getNodalConnectivity().getValues())
        ids=skin.computeFetchedNodeIds()
        self.assertEqual([0,1,2,3,4,5,9,10,14,15,19,20,24,25,26,27,28,29],ids.getValues())
        part=umesh.buildFacePartOfMySelfNode(ids,True)
        part.setName(skin.getName());
        self.assertTrue(part.isEqual(skin,1e-12))
        part2=part[1::2]
        part[::2]=part2
        self.assertTrue(not part.isEqual(skin,1e-12))
        trad=part.zipConnectivityTraducer(0)
        self.assertEqual(9,part.getNumberOfCells())
        self.assertEqual([0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8],trad.getValues())
        pass

    def testUMeshSetPartOfMySelf2(self):
        # resize with explicit ids list
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        self.assertEqual([3,4],m.getAllGeoTypes())
        part=m[[0,3,4]]
        part.simplexize(0)
        part2=part[[1,2,5]]
        m[[0,3,4]]=part2
        self.assertEqual([3,0,4,1,3,1,4,2,3,4,5,2,3,6,7,4,3,7,5,4],m.getNodalConnectivity().getValues())
        self.assertEqual([0,4,8,12,16,20],m.getNodalConnectivityIndex().getValues())
        self.assertEqual([3],m.getAllGeoTypes())
        # no resize with explicit ids list
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        part=m[[0,3]]
        part.convertAllToPoly()
        m[[3,4]]=part
        self.assertEqual([4,0,3,4,1,3,1,4,2,3,4,5,2,5,0,3,4,1,5,6,7,4,3],m.getNodalConnectivity().getValues())
        self.assertEqual([0,5,9,13,18,23],m.getNodalConnectivityIndex().getValues())
        self.assertEqual([3,4,5],m.getAllGeoTypes())
        # resize with range ids
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        part=m[3:]
        m[1:3]=part
        self.assertEqual([4,0,3,4,1,4,6,7,4,3,4,7,8,5,4,4,6,7,4,3,4,7,8,5,4],m.getNodalConnectivity().getValues())
        self.assertEqual([0,5,10,15,20,25],m.getNodalConnectivityIndex().getValues())
        self.assertEqual([4],m.getAllGeoTypes())
        # no resize with range ids
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        part=m[0::3]
        part.convertAllToPoly()
        m[3:]=part
        self.assertEqual([4,0,3,4,1,3,1,4,2,3,4,5,2,5,0,3,4,1,5,6,7,4,3],m.getNodalConnectivity().getValues())
        self.assertEqual([0,5,9,13,18,23],m.getNodalConnectivityIndex().getValues())
        self.assertEqual([3,4,5],m.getAllGeoTypes())
        # no resize with range ids negative direction
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        part=m[3::-3]
        part.convertAllToPoly()
        m[:-3:-1]=part
        self.assertEqual([4,0,3,4,1,3,1,4,2,3,4,5,2,5,0,3,4,1,5,6,7,4,3],m.getNodalConnectivity().getValues())
        self.assertEqual([0,5,9,13,18,23],m.getNodalConnectivityIndex().getValues())
        self.assertEqual([3,4,5],m.getAllGeoTypes())
        pass

    def testUnPolyze3(self):
        coord=[0.0,0.5,-0.5,-0.5,-0.5,-0.5,0.5,-0.5,-0.5,0.0,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5,0.5]
        conn=[1,2,5,4,-1,4,3,0,1,-1,2,0,3,5,-1,0,2,1,-1,4,5,3]
        m=MEDCouplingUMesh.New("a mesh",3);
        m.allocateCells(1);
        m.insertNextCell(NORM_POLYHED,22,conn[0:22])
        m.finishInsertingCells();
        coords=DataArrayDouble(coord,6,3);
        m.setCoords(coords);
        m.checkConsistencyLight();
        #
        vol=m.getMeasureField(False);
        self.assertEqual(1,vol.getArray().getNumberOfTuples());
        self.assertAlmostEqual(0.5,vol.getArray().getIJ(0,0),12)
        #
        m.unPolyze();
        #
        self.assertEqual([NORM_PENTA6],m.getAllGeoTypes())
        self.assertTrue(DataArrayInt([0,7]).isEqual(m.getNodalConnectivityIndex()))
        self.assertTrue(DataArrayInt([16,0,2,1,3,5,4]).isEqual(m.getNodalConnectivity()))
        #
        vol=m.getMeasureField(False);
        self.assertEqual(1,vol.getArray().getNumberOfTuples());
        self.assertAlmostEqual(0.5,vol.getArray().getIJ(0,0),12)
        pass

    def testKrSpatialDiscretization1(self):
        srcPointCoordsX=[0.8401877171547095, 0.7830992237586059, 0.9116473579367843, 0.335222755714889, 0.2777747108031878, 0.4773970518621602, 0.3647844727918433, 0.9522297251747128, 0.6357117279599009, 0.1416025553558034]
        srcFieldValsOnPoints=[2.129892434968836, 2.295320474540621, 1.931948594981134, 2.728013590937196, 2.715603240418478, 2.661778472822935, 2.695696990104364, 1.893710234970982, 2.529628016549284, 2.728432341300668]
        targetPointCoordsX=[-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,-6.93889390391e-17,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45]
        targetFieldValsExpected=[2.975379475824351, 2.95613491917003, 2.936890362515361, 2.917645805861018, 2.898401249206574, 2.879156692552137, 2.859912135897732, 2.840667579243201, 2.821423022588731, 2.802178465934342, 2.78293390927989, 2.763689352625457, 2.744444795971001, 2.725209522098197, 2.709077577124666, 2.706677252549218, 2.727467797847971, 2.713338094723676, 2.671342424824244, 2.664877370146978, 2.653840141412181, 2.619607861392791, 2.569777214476479, 2.513263929794591, 2.450732752808528, 2.368313560985155, 2.250909795670307, 2.098194272085416, 1.954257891732065, 1.895040660973802, 1.865256788315972, 1.835475248687992, 1.80569370905998, 1.775912169431971, 1.746130629803976, 1.716349090175918, 1.686567550547855, 1.656786010919941, 1.627004471291988, 1.597222931663817]
        coeffsExpected=DataArrayDouble.New([52.238272642008695, 26.186513281350948, -173.42106377948534, 324.56733663875184, -104.64968873410248, 34.375030568158316, -256.12372208190425, 105.2292032463934, -16.239907618144965, 7.838025836978943, 2.621910745077291, -0.4902609628247241])
        #
        nbOfInputPoints=10;
        f=MEDCouplingFieldDouble.New(ON_NODES_KR,ONE_TIME);
        srcArrX=DataArrayDouble.New(srcPointCoordsX,nbOfInputPoints,1);
        cmesh=MEDCouplingCMesh.New("aMesh");
        cmesh.setCoordsAt(0,srcArrX);
        umesh=cmesh.buildUnstructured();
        f.setMesh(umesh);
        srcVals=DataArrayDouble.New(srcFieldValsOnPoints,nbOfInputPoints,1);
        f.setArray(srcVals);
        f.checkConsistencyLight();
        #
        res0=f.getValueOn(targetPointCoordsX[:1]);
        self.assertAlmostEqual(targetFieldValsExpected[0],res0[0],10)
        #
        valuesToTest=f.getValueOnMulti(targetPointCoordsX);
        self.assertEqual(40,valuesToTest.getNumberOfTuples());
        self.assertEqual(1,valuesToTest.getNumberOfComponents());
        for i in range(40):
            self.assertAlmostEqual(targetFieldValsExpected[i],valuesToTest.getIJ(i,0),10)
            pass
        fd=f.getDiscretization()
        del f
        self.assertTrue(isinstance(fd,MEDCouplingFieldDiscretizationKriging))
        coeffs,isDrift=fd.computeVectorOfCoefficients(umesh,srcVals)
        self.assertEqual(2,isDrift)
        self.assertTrue(coeffsExpected.isEqual(coeffs,1e-8))
        #
        pass

    def testDuplicateEachTupleNTimes1(self):
        d=DataArrayDouble.New([9.,8.,7.,6.],4,1) ; d.setInfoOnComponents(["mass [kg]"]) ; d.setName("aname")
        d2=d.duplicateEachTupleNTimes(3)
        self.assertTrue(d2.isEqualWithoutConsideringStr(DataArrayDouble.New([9.,9.,9.,8.,8.,8.,7.,7.,7.,6.,6.,6.],4*3,1),1e-14))
        self.assertEqual("aname",d2.getName())
        self.assertEqual(["mass [kg]"],d2.getInfoOnComponents())
        #
        d=DataArrayInt.New([9,8,7,6],4,1) ; d.setInfoOnComponents(["mass [kg]"]) ; d.setName("aname")
        d2=d.duplicateEachTupleNTimes(3)
        self.assertTrue(d2.isEqualWithoutConsideringStr(DataArrayInt.New([9,9,9,8,8,8,7,7,7,6,6,6],4*3,1)))
        self.assertEqual("aname",d2.getName())
        self.assertEqual(["mass [kg]"],d2.getInfoOnComponents())
        pass

    def testSwigComputeTupleIdsNearTuples1(self):
        da=DataArrayDouble([5.,6.,-5.,-6.,5.,-6.,-5.,6.,5.,6.],5,2)
        arr,arrI=da.computeTupleIdsNearTuples(DataArrayDouble([5.,-6.,5.,6.,-5.,-6.],3,2),1e-10)
        self.assertEqual([2,0,4,1],arr.getValues())
        self.assertEqual([0,1,3,4],arrI.getValues())
        arr,arrI=da.computeTupleIdsNearTuples([5.,-6.,5.,6.,-5.,-6.],1e-10)
        self.assertEqual([2,0,4,1],arr.getValues())
        self.assertEqual([0,1,3,4],arrI.getValues())
        expected0=[[2],[0,4],[1]]
        expected1=[[0,1],[0,2],[0,1]]
        for pos,it in enumerate(DataArrayDouble([5.,-6.,5.,6.,-5.,-6.],3,2)):
            arr,arrI=da.computeTupleIdsNearTuples(it,1e-10)
            self.assertEqual(expected0[pos],arr.getValues())
            self.assertEqual(expected1[pos],arrI.getValues())
            pass
        pass

    def testSwigDataTupleIOp1(self):
        d=DataArrayDouble(10,1)
        d.iota(7.)
        for elt in d:
            elt+=2.
            pass
        toTest=DataArrayDouble([9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0])
        self.assertTrue(toTest.isEqual(d,1e-12))
        for elt in d:
            elt-=2.
            pass
        toTest=DataArrayDouble([7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0])
        self.assertTrue(toTest.isEqual(d,1e-12))
        for elt in d:
            elt*=2.
            pass
        toTest=DataArrayDouble([14.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0])
        self.assertTrue(toTest.isEqual(d,1e-12))
        for elt in d:
            elt/=2.
            pass
        toTest=DataArrayDouble([7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0])
        self.assertTrue(toTest.isEqual(d,1e-12))
        #
        d=DataArrayInt(10,1)
        d.iota(7)
        for elt in d:
            elt+=2
            pass
        self.assertEqual(d.getValues(),[9,10,11,12,13,14,15,16,17,18])
        for elt in d:
            elt-=2
            pass
        self.assertEqual(d.getValues(),[7,8,9,10,11,12,13,14,15,16])
        for elt in d:
            elt*=2
            pass
        self.assertEqual(d.getValues(),[14,16,18,20,22,24,26,28,30,32])
        for elt in d:
            elt/=2
            pass
        self.assertEqual(d.getValues(),[7,8,9,10,11,12,13,14,15,16])
        for elt in d:
            elt%=3
            pass
        self.assertEqual(d.getValues(),[1,2,0,1,2,0,1,2,0,1])
        pass

    def testDAIBuildUnique1(self):
        d=DataArrayInt([1,2,2,3,3,3,3,4,5,5,7,7,7,19])
        e=d.buildUnique()
        self.assertTrue(e.isEqual(DataArrayInt([1,2,3,4,5,7,19])))
        pass

    def testDAIPartitionByDifferentValues1(self):
        d=DataArrayInt([1,0,1,2,0,2,2,-3,2])
        expected=[[-3,[7]],[0,[1,4]],[1,[0,2]],[2,[3,5,6,8]]]
        for i,elt in enumerate(zip(*d.partitionByDifferentValues())):
            self.assertEqual(expected[i][0],elt[1])
            self.assertEqual(expected[i][1],elt[0].getValues())
            pass
        pass

    def testFieldGaussMultiDiscPerType1(self):
        coords=DataArrayDouble([0.,0.,0.,1.,1.,1.,1.,0.,0.,0.5,0.5,1.,1.,0.5,0.5,0.],8,2)
        mQ8=MEDCouplingUMesh("",2) ; mQ8.setCoords(coords)
        mQ8.allocateCells(1)
        mQ8.insertNextCell(NORM_QUAD8, list(range(8)))
        mQ8.finishInsertingCells()
        mQ4=MEDCouplingUMesh("",2) ; mQ4.setCoords(coords)
        mQ4.allocateCells(1)
        mQ4.insertNextCell(NORM_QUAD4, list(range(4)))
        mQ4.finishInsertingCells()
        mT3=MEDCouplingUMesh("",2) ; mT3.setCoords(coords)
        mT3.allocateCells(1)
        mT3.insertNextCell(NORM_TRI3, list(range(3)))
        mT3.finishInsertingCells()

        tr=[[0.,0.],[2.,0.], [0.,2.],[2.,2.],[4.,2.],[6.,2.],[8.,2.],[10.,2.],[12.,2.],[0.,4.],[2.,4.],[4.,4.],[6.,4.],[8.,4.],[10.,4.],[12.,4.],[14.,4.],[16.,4.],[18.,4.],[20.,4.],[22.,4.]]
        ms=2*[mQ4]+7*[mQ8]+11*[mT3]
        ms[:]=(elt.deepCopy() for elt in ms)
        for m,t in zip(ms,tr):
            d=m.getCoords() ; d+= t
            pass
        m=MEDCouplingUMesh.MergeUMeshes(ms)
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,NO_TIME)
        f.setMesh(m)
        # throw because cell 0,1 are QUAD4 and cell 3 is QUAD8
        self.assertRaises(InterpKernelException,f.setGaussLocalizationOnCells,[0,1,3],[0.,0.,1.,0.,1.,1.,0.,1.],[0.3,0.3,0.7,0.7],[0.8,0.2])
        f.setGaussLocalizationOnCells([0,1],[0.,0.,1.,0.,1.,1.,0.,1.],[0.3,0.3,0.7,0.7],[0.8,0.2])
        f.setGaussLocalizationOnCells([3,2,5],[0.,0.,1.,0.,1.,1.,0.,1.,0.5,0.,1.,0.5,0.5,1.,0.,0.5],[0.3,0.3,0.7,0.7,0.9,0.9],[0.8,0.05,0.15])
        f.setGaussLocalizationOnCells([4,6,8,7],[0.,0.,1.,0.,1.,1.,0.,1.,0.5,0.,1.,0.5,0.5,1.,0.,0.5],[0.3,0.3,0.7,0.7,0.9,0.9,-0.1,0.3],[0.7,0.05,0.15,0.1])
        f.setGaussLocalizationOnCells([9,10,11,12,13],[0.,0.,1.,0.,1.,1.],[0.4,0.4],[1.])
        f.setGaussLocalizationOnCells([14,15,16,17,18,19],[0.,0.,1.,0.,1.,1.],[0.4,0.4,0.14,0.16],[0.22,0.78])
        self.assertEqual(46,f.getNumberOfTuplesExpected())
        vals=DataArrayDouble.New(46*3,1) ; vals.iota(7.7) ; vals.rearrange(3)
        f.setArray(vals)
        f.checkConsistencyLight()
        #f.getLocalizationOfDiscr()
        self.assertRaises(InterpKernelException,f.getGaussLocalizationIdOfOneType,NORM_QUAD8) #throw because several loc
        self.assertEqual([1,2],f.getGaussLocalizationIdsOfOneType(NORM_QUAD8))
        self.assertEqual([0,0,1,1,2,1,2,2,2,3,3,3,3,3,4,4,4,4,4,4],f.getDiscretization().getArrayOfDiscIds().getValues())
        fc=f[[1,2,3,8]]
        fc.checkConsistencyLight()
        self.assertTrue(DataArrayDouble([13.7,14.7,15.7,16.7,17.7,18.7,19.7,20.7,21.7,22.7,23.7,24.7,25.7,26.7,27.7,28.7,29.7,30.7,31.7,32.7,33.7,34.7,35.7,36.7,82.7,83.7,84.7,85.7,86.7,87.7,88.7,89.7,90.7,91.7,92.7,93.7],12,3).isEqual(fc.getArray(),1e-10))
        fc.renumberCells([3,2,0,1])
        self.assertTrue(DataArrayDouble([28.7, 29.7, 30.7, 31.7, 32.7, 33.7, 34.7, 35.7, 36.7, 82.7, 83.7, 84.7, 85.7, 86.7, 87.7, 88.7, 89.7, 90.7, 91.7, 92.7, 93.7, 19.7, 20.7, 21.7, 22.7, 23.7, 24.7, 25.7, 26.7, 27.7, 13.7, 14.7, 15.7, 16.7, 17.7, 18.7],12,3).isEqual(fc.getArray(),1e-10))
        fc.getArray()
        pass

    def testSwigRotate(self):
        d=DataArrayDouble([1.,2.,3.,4.,6.,5.],2,3)
        MEDCouplingPointSet.Rotate3DAlg([0.,0.,0.],[0.,1.,0.],1.5707963267948966,d)
        self.assertTrue(d.isEqual(DataArrayDouble([3.,2.,-1.,5.,6.,-4.],2,3),1e-12))
        d=DataArrayDouble([1.,2.,3.,4.,6.,5.],3,2)
        MEDCouplingPointSet.Rotate2DAlg([0.,0.],1.5707963267948966,d)
        self.assertTrue(d.isEqual(DataArrayDouble([-2.,1.,-4.,3.,-5.,6.],3,2),1e-12))
        pass

    def testSwigCMeshProtection(self):
        cm=MEDCouplingCMesh()
        self.assertRaises(InterpKernelException,cm.setCoordsAt,0,DataArrayDouble([4.,4.5,6.,7.],2,2))
        self.assertRaises(InterpKernelException,cm.setCoords,DataArrayDouble([4.,4.5,6.,7.],2,2))
        pass

    def testSwigCellsInBoundingBox1(self):
        m3D=MEDCouplingDataForTest.build3DExtrudedUMesh_1()[0]
        self.assertTrue(m3D.getCellsInBoundingBox([(0,3),(0,3),(0,1)],-1e-12).isEqual(DataArrayInt([0,1,2,3,4,5])))
        self.assertRaises(InterpKernelException,m3D.getCellsInBoundingBox,[(0,3,0),(3,0,1)],-1e-12)
        pass

    def testDAICheckMonotonic1(self):
        data1=[-1,0,2,2,4,5]
        data2=[6,2,0,-8,-9,-56]
        data3=[-1,0,3,2,4,6]
        data4=[7,5,2,3,0,-6]
        d=DataArrayInt.New(data1);
        self.assertTrue(d.isMonotonic(True));
        self.assertTrue(not d.isMonotonic(False));
        d.checkMonotonic(True);
        self.assertRaises(InterpKernelException,d.checkMonotonic,False)
        d=DataArrayInt.New(data2);
        self.assertTrue(d.isMonotonic(False));
        self.assertTrue(not d.isMonotonic(True));
        d.checkMonotonic(False);
        self.assertRaises(InterpKernelException,d.checkMonotonic,True)
        d=DataArrayInt.New(data3);
        self.assertTrue(not d.isMonotonic(False));
        self.assertTrue(not d.isMonotonic(True));
        self.assertRaises(InterpKernelException,d.checkMonotonic,True)
        self.assertRaises(InterpKernelException,d.checkMonotonic,False)
        d=DataArrayInt.New(data4);
        self.assertTrue(not d.isMonotonic(False));
        self.assertTrue(not d.isMonotonic(True));
        self.assertRaises(InterpKernelException,d.checkMonotonic,True)
        self.assertRaises(InterpKernelException,d.checkMonotonic,False)
        d=DataArrayInt.New(0,1)
        self.assertTrue(d.isMonotonic(True));
        self.assertTrue(d.isMonotonic(False));
        d.checkMonotonic(True);
        d.checkMonotonic(False);
        d=DataArrayInt.New(data4,3,2);#throw because nbComp!=1
        self.assertRaises(InterpKernelException,d.isMonotonic,True)
        self.assertRaises(InterpKernelException,d.isMonotonic,False)
        self.assertRaises(InterpKernelException,d.checkMonotonic,True)
        self.assertRaises(InterpKernelException,d.checkMonotonic,False)
        pass

    def testSwigDASetItemOnEmpty1(self):
        d=DataArrayInt(0,1)
        isThrow=False
        try:
            d[0:1000:2]=4
        except InterpKernelException as e:
            isThrow=True
            pass
        self.assertTrue(isThrow)
        d[:]=4
        d[::2]=5
        #
        d=DataArrayDouble(0,1)
        isThrow=False
        try:
            d[0:1000:2]=4
        except InterpKernelException as e:
            isThrow=True
            pass
        self.assertTrue(isThrow)
        d[:]=4
        d[::2]=5
        d=DataArrayInt([],0,1)
        d2=DataArrayInt(0)
        self.assertTrue(d2.isEqual(d))
        d=DataArrayDouble([],0,1)
        d2=DataArrayDouble(0)
        self.assertTrue(d2.isEqual(d,1e-12))
        pass

    def testSwigDAITransformWithIndArr1(self):
        arr=DataArrayInt([0,4,5,1])
        d=DataArrayInt([7,8,9,10])
        self.assertRaises(InterpKernelException,arr.transformWithIndArr,d)
        pass

    def testIntersect2DMeshesTmp6(self):
        # coordinates
        coords=DataArrayDouble.New([2.7554552980815448e-15,45,-45,5.5109105961630896e-15,-31.819805153394636,31.81980515339464,2.8779199779962799e-15,47,2.8166876380389124e-15,46,-47,5.7558399559925599e-15,-33.234018715767732,33.234018715767739,-46,5.6333752760778247e-15],8,2);
        # connectivity
        conn=DataArrayInt.New([8,0,3,5,1,4,6,7,2])
        connI=DataArrayInt.New([0,9]);
        m1=MEDCouplingUMesh.New("Fixe",2);
        m1.setCoords(coords);
        m1.setConnectivity(conn,connI,True);
        #
        coords=DataArrayDouble.New([-7.3800475508445391,41.854329503018846,-3.7041190667754655,42.338274668899189,-3.7041190667754655,45.338274668899189,-7.3800475508445382,44.854329503018839,-5.5473631693521845,42.136406608386956,-3.7041190667754655,43.838274668899189,-5.5420833088100014,45.09630208595901,-7.3800475508445382,43.354329503018839,-3.7041190667754651,52.338274668899189,-7.3800475508445382,51.854329503018839,-3.7041190667754655,48.838274668899189,-5.5420833088100014,52.09630208595901,-7.3800475508445382,48.354329503018839],13,2);
        # connectivity
        conn=DataArrayInt.New([8,0,1,2,3,4,5,6,7,8,3,2,8,9,6,10,11,12]);
        connI=DataArrayInt.New([0,9,18]);
        #
        m2=MEDCouplingUMesh.New("Mobile",2);
        m2.setCoords(coords);
        m2.setConnectivity(conn,connI,True);
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m1,m2,1e-10);
        self.assertTrue(d1.isEqual(DataArrayInt([0,0,0,0])));
        self.assertTrue(d2.isEqual(DataArrayInt([0,1,-1,-1])));
        self.assertEqual(4,m3.getNumberOfCells());
        self.assertEqual(4,d1.getNumberOfTuples());
        self.assertEqual(4,d2.getNumberOfTuples());
        self.assertEqual(43,m3.getNumberOfNodes());
        dI,areMerged,newNbOfNodes=m3.mergeNodes(1e-12)
        self.assertEqual(35,m3.getNumberOfNodes());
        m3.zipCoords();
        self.assertEqual(23,m3.getNumberOfNodes());
        #
        f=m3.getMeasureField(True);
        valuesExpected=DataArrayDouble([1.6603638692585716,5.747555728471923,129.68907101754394,7.4162714498559694])
        self.assertTrue(f.getArray().isEqual(valuesExpected,1e-12))
        pass

    def testDAPushBack(self):
        d=DataArrayDouble(0,1)
        for i in range(8):
            d.pushBackSilent(i)
            pass
        self.assertEqual(d.getNumberOfTuples(),8)
        self.assertEqual(d.getNbOfElemAllocated(),8)
        d.pushBackSilent(4.44)
        self.assertEqual(d.getNumberOfTuples(),9)
        self.assertEqual(d.getNbOfElemAllocated(),16)
        self.assertTrue(d.isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,4.44]),1e-12))
        e=d.deepCopy()
        self.assertEqual(e.getNumberOfTuples(),9)
        self.assertEqual(e.getNbOfElemAllocated(),9)
        self.assertTrue(e.isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,4.44]),1e-12))
        self.assertAlmostEqual(d.popBackSilent(),4.44,12)
        self.assertEqual(d.getNumberOfTuples(),8)
        self.assertEqual(d.getNbOfElemAllocated(),16)
        self.assertTrue(d.isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.]),1e-12))
        f=DataArrayDouble()
        f.reserve(1000)
        f.pushBackSilent(4.)
        self.assertTrue(f.isEqual(DataArrayDouble([4.]),1e-12))
        self.assertEqual(f.getNumberOfTuples(),1)
        self.assertEqual(f.getNbOfElemAllocated(),1000)
        ff=f[:]
        self.assertTrue(ff.isEqual(DataArrayDouble([4.]),1e-12))
        self.assertEqual(ff.getNumberOfTuples(),1)
        self.assertEqual(ff.getNbOfElemAllocated(),1)
        d=DataArrayDouble()
        d.pushBackSilent(4.44)
        d.pushBackSilent(5.55)
        d.pushBackSilent(6.66)
        self.assertTrue(d.isEqual(DataArrayDouble([4.44,5.55,6.66]),1e-12))
        #
        d=DataArrayInt(0,1)
        for i in range(8):
            d.pushBackSilent(i)
            pass
        self.assertEqual(d.getNumberOfTuples(),8)
        self.assertEqual(d.getNbOfElemAllocated(),8)
        d.pushBackSilent(444)
        self.assertEqual(d.getNumberOfTuples(),9)
        self.assertEqual(d.getNbOfElemAllocated(),16)
        self.assertTrue(d.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,444])))
        e=d.deepCopy()
        self.assertEqual(e.getNumberOfTuples(),9)
        self.assertEqual(e.getNbOfElemAllocated(),9)
        self.assertTrue(e.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,444])))
        self.assertEqual(d.popBackSilent(),444)
        self.assertEqual(d.getNumberOfTuples(),8)
        self.assertEqual(d.getNbOfElemAllocated(),16)
        self.assertTrue(d.isEqual(DataArrayInt([0,1,2,3,4,5,6,7])))
        f=DataArrayInt()
        f.reserve(1000)
        f.pushBackSilent(4)
        self.assertTrue(f.isEqual(DataArrayInt([4])))
        self.assertEqual(f.getNumberOfTuples(),1)
        self.assertEqual(f.getNbOfElemAllocated(),1000)
        ff=f[:]
        self.assertTrue(ff.isEqual(DataArrayInt([4])))
        self.assertEqual(ff.getNumberOfTuples(),1)
        self.assertEqual(ff.getNbOfElemAllocated(),1)
        d=DataArrayInt()
        d.pushBackSilent(444)
        d.pushBackSilent(555)
        d.pushBackSilent(666)
        self.assertTrue(d.isEqual(DataArrayInt([444,555,666])))
        #
        d=DataArrayInt()
        d.alloc(10,1)
        d.setInfoOnComponent(0,"ABC")
        d.setName("dEf")
        d.iota(7)
        e=DataArrayInt([7,8,9,10,11,12,13,14,15,16]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e))
        self.assertEqual(10,d.getNbOfElemAllocated())
        d.pushBackSilent(55)
        e=DataArrayInt([7,8,9,10,11,12,13,14,15,16,55]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e))
        self.assertEqual(20,d.getNbOfElemAllocated())
        d.reserve(4)
        e=DataArrayInt([7,8,9,10]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e))
        self.assertEqual(4,d.getNbOfElemAllocated())
        d.pushBackSilent(5)
        e=DataArrayInt([7,8,9,10,5]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e))
        self.assertEqual(8,d.getNbOfElemAllocated())
        self.assertEqual(5,d.popBackSilent())
        e=DataArrayInt([7,8,9,10]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e))
        self.assertEqual(8,d.getNbOfElemAllocated())
        self.assertRaises(OverflowError,d.reserve,-1)
        e=DataArrayInt([7,8,9,10]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e))
        self.assertEqual(8,d.getNbOfElemAllocated())
        d.reserve(0)
        e=DataArrayInt([]) ; e.setInfoOnComponent(0,"ABC") ; e.setName("dEf") ; self.assertTrue(d.isEqual(e))
        self.assertEqual(0,d.getNbOfElemAllocated())
        #
        d=DataArrayDouble()
        d.alloc(10,1)
        d.setInfoOnComponent(0,"ABC")
        d.setName("dEf")
        d.iota(7)
        e=DataArrayDouble([7,8,9,10,11,12,13,14,15,16]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e,1e-14))
        self.assertEqual(10,d.getNbOfElemAllocated())
        d.pushBackSilent(55)
        e=DataArrayDouble([7,8,9,10,11,12,13,14,15,16,55]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e,1e-14))
        self.assertEqual(20,d.getNbOfElemAllocated())
        d.reserve(4)
        e=DataArrayDouble([7,8,9,10]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e,1e-14))
        self.assertEqual(4,d.getNbOfElemAllocated())
        d.pushBackSilent(5)
        e=DataArrayDouble([7,8,9,10,5]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e,1e-14))
        self.assertEqual(8,d.getNbOfElemAllocated())
        self.assertEqual(5.,d.popBackSilent())
        e=DataArrayDouble([7,8,9,10]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e,1e-14))
        self.assertEqual(8,d.getNbOfElemAllocated())
        self.assertRaises(OverflowError,d.reserve,-1)
        e=DataArrayDouble([7,8,9,10]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e,1e-14))
        self.assertEqual(8,d.getNbOfElemAllocated())
        d.reserve(0)
        e=DataArrayDouble([]) ; e.setInfoOnComponent(0,"ABC") ; e.setName("dEf") ; self.assertTrue(d.isEqual(e,1e-14))
        self.assertEqual(0,d.getNbOfElemAllocated())
        pass

    def testDAIBuildSubstractionOptimized1(self):
        da1=DataArrayInt.New([1,3,5,6,7,9,13])
        da2=DataArrayInt.New([3,5,9])
        da3=DataArrayInt.New([1,3,5])
        da4=DataArrayInt.New([1,3,5,6,7,9,13])
        #
        a=da1.buildSubstractionOptimized(da2);
        self.assertTrue(a.isEqual(DataArrayInt([1,6,7,13])));
        #
        a=da1.buildSubstractionOptimized(da3);
        self.assertTrue(a.isEqual(DataArrayInt([6,7,9,13])));
        #
        a=da1.buildSubstractionOptimized(da4);
        self.assertTrue(a.isEqual(DataArrayInt([])));
        pass

    def testDAIIsStrictlyMonotonic1(self):
        da1=DataArrayInt.New([1,3,5,6,7,9,13])
        self.assertTrue(da1.isStrictlyMonotonic(True));
        da1.checkStrictlyMonotonic(True);
        self.assertTrue(da1.isMonotonic(True));
        da1.checkMonotonic(True);
        self.assertTrue(not da1.isStrictlyMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,False)
        self.assertTrue(not da1.isMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,False)
        #
        da1=DataArrayInt.New([1,3,5,6,6,9,13])
        self.assertTrue(not da1.isStrictlyMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,True)
        self.assertTrue(da1.isMonotonic(True));
        da1.checkMonotonic(True);
        self.assertTrue(not da1.isStrictlyMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,False)
        self.assertTrue(not da1.isMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,False)
        #
        da1=DataArrayInt.New([1,3,5,6,5,9,13])
        self.assertTrue(not da1.isStrictlyMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,True)
        self.assertTrue(not da1.isMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,True)
        self.assertTrue(not da1.isStrictlyMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,False)
        self.assertTrue(not da1.isMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,False)
        #
        da1=DataArrayInt.New([13,9,7,6,5,3,1])
        self.assertTrue(not da1.isStrictlyMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,True)
        self.assertTrue(not da1.isMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,True)
        self.assertTrue(da1.isStrictlyMonotonic(False));
        da1.checkStrictlyMonotonic(False);
        self.assertTrue(da1.isMonotonic(False));
        da1.checkMonotonic(False);
        #
        da1=DataArrayInt.New([13,9,6,6,5,3,1])
        self.assertTrue(not da1.isStrictlyMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,True)
        self.assertTrue(not da1.isMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,True)
        self.assertTrue(not da1.isStrictlyMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,False)
        self.assertTrue(da1.isMonotonic(False));
        da1.checkMonotonic(False);
        #
        da1=DataArrayInt.New([13,9,5,6,5,3,1])
        self.assertTrue(not da1.isStrictlyMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,True)
        self.assertTrue(not da1.isMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,True)
        self.assertTrue(not da1.isStrictlyMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,False)
        self.assertTrue(not da1.isMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,False)
        #
        da1=DataArrayInt.New([])
        self.assertTrue(da1.isStrictlyMonotonic(True));
        da1.checkStrictlyMonotonic(True);
        self.assertTrue(da1.isMonotonic(True));
        da1.checkMonotonic(True);
        self.assertTrue(da1.isStrictlyMonotonic(False));
        da1.checkStrictlyMonotonic(False);
        self.assertTrue(da1.isMonotonic(False));
        da1.checkMonotonic(False);
        #
        da1=DataArrayInt.New([13])
        self.assertTrue(da1.isStrictlyMonotonic(True));
        da1.checkStrictlyMonotonic(True);
        self.assertTrue(da1.isMonotonic(True));
        da1.checkMonotonic(True);
        self.assertTrue(da1.isStrictlyMonotonic(False));
        da1.checkStrictlyMonotonic(False);
        self.assertTrue(da1.isMonotonic(False));
        da1.checkMonotonic(False);
        pass

    def testFindAndCorrectBadOriented3DCells1(self):
        nbOfDisc=20
        vects=([0,0,-1],[0.3,0.7,0.2],[-0.3,0.7,0.2],[-0.3,-0.7,0.2])
        #
        m0=MEDCouplingUMesh("m",3) ; m0.allocateCells(0); m0.insertNextCell(NORM_TETRA4,[0,1,2,3]); #Well oriented
        m1=MEDCouplingUMesh("m",3) ; m1.allocateCells(0); m1.insertNextCell(NORM_PYRA5,[0,1,2,3,4]); #Well oriented
        m2=MEDCouplingUMesh("m",3) ; m2.allocateCells(0); m2.insertNextCell(NORM_PENTA6,[0,1,2,3,4,5]); #Well oriented 
        m3=MEDCouplingUMesh("m",3) ; m3.allocateCells(0); m3.insertNextCell(NORM_HEXA8,[0,1,2,3,4,5,6,7]); #Well oriented
        m4=MEDCouplingUMesh("m",3) ; m4.allocateCells(0)
        self.assertRaises(InterpKernelException,m4.insertNextCell,NORM_HEXGP12,[0,1,2,3,4,5,6,7,8,9,10,11,12]);
        m4.insertNextCell(NORM_HEXGP12,[0,1,2,3,4,5,6,7,8,9,10,11]); #Well oriented
        c0=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,0.,0.,0.,0.,1.],4,3) ; m0.setCoords(c0)
        c1=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,1.,0.,1.,0.,0.,0.,0.,1.],5,3) ; m1.setCoords(c1)
        c2=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,0.,0., 0.,0.,1.,0.,1.,1.,1.,0.,1.],6,3) ; m2.setCoords(c2)
        c3=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,1.,0.,1.,0.,0.,0.,0.,1.,0.,1.,1.,1.,1.,1.,1.,0.,1.],8,3) ; m3.setCoords(c3)
        c4=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,1.,0.,1.,0.,0.,0.8,0.,0.,0.45,0.,0.,   0.,0.,1.,0.,1.,1.,1.,1.,1.,1.,0.,1.,0.8,0.,1.,0.45,0.,1.],12,3) ; m4.setCoords(c4)
        m=MEDCouplingMesh.MergeMeshes([m0,m1,m2,m3,m4])
        expected1=DataArrayDouble([0.16666666666666666,0.3333333333333333,0.5,1.,1.])
        for v in vects:
            for i in range(nbOfDisc):
                mm=m.deepCopy()
                mm.rotate([0.,0.,0.],[0.3,0.7,0.2],float(i)/float(nbOfDisc)*2*pi)
                mm2=mm.deepCopy()
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected1,1e-14))
                self.assertTrue(mm.findAndCorrectBadOriented3DCells().empty())
                self.assertTrue(mm.isEqual(mm2,1e-14))
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected1,1e-14))
                mm.convertAllToPoly()
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected1,1e-14))
                pass
            pass
        #
        mOK=m.deepCopy()
        m0=MEDCouplingUMesh("m",3) ; m0.allocateCells(0); m0.insertNextCell(NORM_TETRA4,[0,2,1,3]); #Not well oriented
        m1=MEDCouplingUMesh("m",3) ; m1.allocateCells(0); m1.insertNextCell(NORM_PYRA5,[0,1,2,3,4]); #Well oriented 
        m2=MEDCouplingUMesh("m",3) ; m2.allocateCells(0); m2.insertNextCell(NORM_PENTA6,[0,1,2,3,4,5]); #Well oriented 
        m3=MEDCouplingUMesh("m",3) ; m3.allocateCells(0); m3.insertNextCell(NORM_HEXA8,[0,3,2,1,4,7,6,5]); #Not well oriented
        m4=MEDCouplingUMesh("m",3) ; m4.allocateCells(0); m4.insertNextCell(NORM_HEXGP12,[0,5,4,3,2,1,6,11,10,9,8,7]); #Not well oriented
        m0.setCoords(c0) ; m1.setCoords(c1) ; m2.setCoords(c2) ; m3.setCoords(c3) ; m4.setCoords(c4)
        m=MEDCouplingMesh.MergeMeshes([m0,m1,m2,m3,m4])
        expected2=DataArrayDouble([-0.16666666666666666,0.3333333333333333,0.5,-1.,-1.])
        for v in vects:
            for i in range(nbOfDisc):
                mm=m.deepCopy()
                mm.rotate([0.,0.,0.],[0.3,0.7,0.2],float(i)/float(nbOfDisc)*2*pi)
                mm2=mm.deepCopy() ; mm3=mm.deepCopy() ; mm3.convertAllToPoly()
                self.assertTrue(mm3.getMeasureField(False).getArray().isEqual(expected2,1e-14))
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected2,1e-14))
                self.assertTrue(mm.findAndCorrectBadOriented3DCells().isEqual(DataArrayInt([0,3,4])))
                mOK.setCoords(mm.getCoords())
                self.assertTrue(mm.isEqual(mOK,1e-14))
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected1,1e-14))
                mmm=mm.deepCopy()
                self.assertTrue(mmm.findAndCorrectBadOriented3DCells().empty())
                mm.convertAllToPoly()
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected1,1e-14))
                pass
            pass
        #
        m0=MEDCouplingUMesh("m",3) ; m0.allocateCells(0); m0.insertNextCell(NORM_TETRA4,[0,1,2,3]); #Well oriented
        m1=MEDCouplingUMesh("m",3) ; m1.allocateCells(0); m1.insertNextCell(NORM_PYRA5,[0,3,2,1,4]); #Not well oriented 
        m2=MEDCouplingUMesh("m",3) ; m2.allocateCells(0); m2.insertNextCell(NORM_PENTA6,[0,2,1,3,5,4]); #Not well oriented 
        m3=MEDCouplingUMesh("m",3) ; m3.allocateCells(0); m3.insertNextCell(NORM_HEXA8,[0,1,2,3,4,5,6,7]); #Well oriented
        m4 = MEDCouplingUMesh("m", 3) ; m4.allocateCells(0); m4.insertNextCell(NORM_HEXGP12, list(range(12)));  # Well oriented
        m0.setCoords(c0) ; m1.setCoords(c1) ; m2.setCoords(c2) ; m3.setCoords(c3) ; m4.setCoords(c4)
        m=MEDCouplingMesh.MergeMeshes([m0,m1,m2,m3,m4])
        expected3=DataArrayDouble([0.16666666666666666,-0.3333333333333333,-0.5,1.,1.])
        for v in vects:
            for i in range(nbOfDisc):
                mm=m.deepCopy()
                mm.rotate([0.,0.,0.],[0.3,0.7,0.2],float(i)/float(nbOfDisc)*2*pi)
                mm2=mm.deepCopy() ; mm3=mm.deepCopy() ; mm3.convertAllToPoly()
                self.assertTrue(mm3.getMeasureField(False).getArray().isEqual(expected3,1e-14))
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected3,1e-14))
                self.assertTrue(mm.findAndCorrectBadOriented3DCells().isEqual(DataArrayInt([1,2])))
                mOK.setCoords(mm.getCoords())
                self.assertTrue(mm.isEqual(mOK,1e-14))
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected1,1e-14))
                mmm=mm.deepCopy()
                self.assertTrue(mmm.findAndCorrectBadOriented3DCells().empty())
                mm.convertAllToPoly()
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected1,1e-14))
                pass
            pass
        pass

    def testSwig2CellOrientation1(self):
        coords=DataArrayDouble([-0.21606,-0.10803,0.29999999999999999,-0.21606,-0.10803,0.37700000000000006,0,-0.10803,0.29999999999999999,0,-0.10803,0.37700000000000006,0,0.10803,0.29999999999999999,0,0.10803,0.37700000000000006,-0.21606,0.10803,0.29999999999999999,-0.21606,0.10803,0.37700000000000006,0,0.03601,0.29999999999999999,0,0.03601,0.37700000000000006,0,-0.03601,0.29999999999999999,0,-0.03601,0.37700000000000006],12,3)
        conn=[[0,2,10,8,4,6],[1,3,11,9,5,7],[0,1,3,2],[2,3,11,10],[10,11,9,8],[8,9,5,4],[4,5,7,6],[6,7,1,0]]
        for i in range(256):
            mesh=MEDCouplingUMesh("FluidMesh_1",3);
            mesh.allocateCells(0)
            conn2=[elt[:] for elt in conn]
            code=bin(i)[2:] ; code='0'*(8-len(code))+code
            for face,rev in zip(conn2,code):
                if bool(int(rev)):
                    face.reverse()
                    pass
                pass
            conn3=[elt+[-1] for elt in conn2]
            conn3=sum(conn3,[])[:-1]
            mesh.insertNextCell(NORM_POLYHED,conn3)
            mesh.setCoords(coords)
            mesh.orientCorrectlyPolyhedrons()
            self.assertTrue(mesh.computeCellCenterOfMass().isEqual(DataArrayDouble([-0.10803,0.,0.3385],1,3),1e-12))
            pass
        pass

    def testSwig2CheckConsecutiveCellTypesForMEDFileFrmt1(self):
        m1=MEDCouplingUMesh("",2) ; m1.allocateCells(0)
        m1.insertNextCell(NORM_QUAD4,[0,1,2,3])
        m1.insertNextCell(NORM_TRI3,[0,1,2])
        d=DataArrayDouble(4,3) ; d[:]=0.
        m1.setCoords(d)
        self.assertTrue(m1.checkConsecutiveCellTypes())
        self.assertTrue(not m1.checkConsecutiveCellTypesForMEDFileFrmt())
        m1.renumberCells([1,0])
        self.assertTrue(m1.checkConsecutiveCellTypes())
        self.assertTrue(m1.checkConsecutiveCellTypesForMEDFileFrmt())
        pass

    def testSwig2DAAccumulate1(self):
        d=DataArrayInt(10) ; d.iota(0)
        self.assertEqual([45],d.accumulate())
        self.assertEqual(45,d.accumulate(0))
        d=DataArrayInt(30) ; d.iota(0) ; d.rearrange(3)
        self.assertEqual([135,145,155],d.accumulate())
        self.assertEqual(135,d.accumulate(0))
        self.assertEqual(145,d.accumulate(1))
        self.assertEqual(155,d.accumulate(2))
        d=DataArrayDouble(10) ; d.iota(0.)
        self.assertEqual([45.],d.accumulate())
        self.assertEqual(45.,d.accumulate(0))
        d=DataArrayDouble(30) ; d.iota(0) ; d.rearrange(3)
        self.assertEqual([135.,145.,155.],d.accumulate())
        self.assertEqual(135.,d.accumulate(0))
        self.assertEqual(145.,d.accumulate(1))
        self.assertEqual(155.,d.accumulate(2))
        pass

    def testSwig2UMeshDistanceToMesh1(self):
        m=MEDCouplingUMesh("toto",2)
        coords=DataArrayDouble([2.3,3.4,5.6,6.5,-4.3,3.2,-9.8,7.6,-5.4],3,3)
        m.setCoords(coords)
        m.allocateCells(0)
        m.insertNextCell(NORM_TRI3,[0,1,2])
        a,b=m.distanceToPoint([-0.335,2.27,1.21])
        self.assertEqual(0,b)
        self.assertAlmostEqual(0.0223609881003,a,12);
        a,b=m.distanceToPoint(DataArrayDouble([-0.335,2.27,1.21],1,3))
        self.assertEqual(0,b)
        self.assertAlmostEqual(0.0223609881003,a,12);
        a,b=coords.distanceToTuple([-0.335,2.27,1.21])
        self.assertAlmostEqual(5.243302871282566,a,14)
        self.assertEqual(0,b)
        #
        m=MEDCouplingUMesh("toto",2)
        coords=DataArrayDouble([0.,0.,0., 8.,0.,0., 8.,8.,0., 0.,8.,0.],4,3)
        m.setCoords(coords)
        m.allocateCells(0)
        m.insertNextCell(NORM_QUAD4,[0,1,2,3])
        m.checkConsistency()
        self.assertEqual([4,0,1,2,3],m.getNodalConnectivity().getValues())
        a,b=m.distanceToPoint([5.,2.,0.1])
        self.assertAlmostEqual(0.1,a,14) ; self.assertEqual(0,b)
        a,b=m.distanceToPoint([5.,-2.,4.])
        self.assertAlmostEqual(sqrt(2*2+4*4),a,14) ; self.assertEqual(0,b)
        m.allocateCells(0)
        m.insertNextCell(NORM_POLYGON,[0,1,2,3])
        m.checkConsistency()
        self.assertEqual([5,0,1,2,3],m.getNodalConnectivity().getValues())
        a,b=m.distanceToPoint([11.,3.,4.])
        self.assertAlmostEqual(sqrt(3*3+4*4),a,14) ; self.assertEqual(0,b)
        a,b=m.distanceToPoint([4.,12.,5.])
        self.assertAlmostEqual(sqrt(4*4+5*5),a,14) ; self.assertEqual(0,b)
        d=DataArrayDouble([-1.2,3.,2.],1,3)
        for elt in d:
            a,b=m.distanceToPoint(d)
            self.assertAlmostEqual(sqrt(1.2*1.2+2*2),a,14) ; self.assertEqual(0,b)
            pass
        #
        m=MEDCouplingUMesh("toto",1)
        coords=DataArrayDouble([0.,0.,4.,0.,0.,4.],3,2) ; m.setCoords(coords)
        m.allocateCells(0) ; m.insertNextCell(NORM_SEG2,[0,1]) ; m.insertNextCell(NORM_SEG2,[1,2])
        a,b=m.distanceToPoint([-0.1,4.1])
        self.assertAlmostEqual(0.14142135623730925,a,14)  # b==1 self.assertEqual(2,c)
        a,b=m.distanceToPoint([0.,3.9])
        self.assertAlmostEqual(0.07071067811865482,a,14) ; self.assertEqual(1,b) # self.assertEqual(2,c)
        pass

    def testSwig2UMeshDistanceToMesh2(self):
        mesh = MEDCouplingUMesh('Solid_3', 2)
        coo = DataArrayDouble([(99.75,-1.42109e-14,102.75),(99.75,200,102.75),(2.5,0,200),(2.5,200,200),(197,0,200),(197,200,200)])
        mesh.setCoords(coo)
        c = DataArrayInt([3, 4, 0, 1, 3, 4, 1, 5, 3, 1, 0, 3, 3, 3, 0, 2])
        cI = DataArrayInt([0, 4, 8, 12, 16])
        mesh.setConnectivity(c, cI)
        mesh.checkConsistency()
        pt = [125.0, 175.0, 175.0]
        # Values computed from GEOM:
        exp1, exp2, exp3, exp4 = 54.0633707597, 33.2340187158, 68.9429111657, 99.5221476482
        d1, _ = mesh[0].distanceToPoint(pt)
        d2, _ = mesh[1].distanceToPoint(pt)
        d3, _ = mesh[2].distanceToPoint(pt)
        d4, _ = mesh[3].distanceToPoint(pt)
        self.assertAlmostEqual(exp1,d1,10)
        self.assertAlmostEqual(exp2,d2,10)
        self.assertAlmostEqual(exp3,d3,10)
        self.assertAlmostEqual(exp4,d4,10)
        pass

    def testSwig2NonRegressionPartitionBySpreadZone1(self):
        m=MEDCouplingCMesh()
        arr=DataArrayDouble(6) ; arr.iota(0.)
        m.setCoords(arr,arr,arr)
        m=m.buildUnstructured()
        mPart=m[50,80,85,87,92,122]
        zones=mPart.partitionBySpreadZone()
        self.assertEqual(4,len(zones))
        self.assertTrue(zones[0].isEqual(DataArrayInt([0])))
        self.assertTrue(zones[1].isEqual(DataArrayInt([1,2])))
        self.assertTrue(zones[2].isEqual(DataArrayInt([3,4])))
        self.assertTrue(zones[3].isEqual(DataArrayInt([5])))
        #
        n,ni=m.computeNeighborsOfCells()
        a,b=MEDCouplingUMesh.ComputeSpreadZoneGraduallyFromSeed(0,n,ni)
        self.assertEqual(13,b) ; self.assertTrue(a.isIota(125))
        a,b=MEDCouplingUMesh.ComputeSpreadZoneGraduallyFromSeed([1],n,ni)
        self.assertEqual(12,b) ; self.assertTrue(a.isIota(125))
        a,b=MEDCouplingUMesh.ComputeSpreadZoneGraduallyFromSeed((2,),n,ni)
        self.assertEqual(11,b) ; self.assertTrue(a.isIota(125))
        a,b=MEDCouplingUMesh.ComputeSpreadZoneGraduallyFromSeed(DataArrayInt([3]),n,ni)
        self.assertEqual(12,b) ; self.assertTrue(a.isIota(125))
        pass

    def testSwigUMeshInsertNextCell1(self):
        m=MEDCouplingUMesh("toto",2)
        #
        coords=DataArrayDouble([0.,0.,1.,1.,1.,0.]) ; m.setCoords(coords)
        da=DataArrayInt([0,1,2])
        m.allocateCells(0)
        for i in range(5):
            m.insertNextCell(NORM_TRI3,da)
            pass
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,8,12,16,20])))
        #
        da=DataArrayInt([0,1,2,3])
        m.allocateCells(0)
        for i in range(5):
            m.insertNextCell(NORM_TRI3,3,da)
            pass
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,8,12,16,20])))
        #
        da=DataArrayInt([0,1])
        m.allocateCells(0)
        self.assertRaises(InterpKernelException,m.insertNextCell,NORM_TRI3,3,da)
        #
        da=DataArrayInt([0,1,2,0,1,3,0,1,4,0,1,5,0,1,6],5,3)
        m.allocateCells(0)
        for t in da:
            m.insertNextCell(NORM_TRI3,t)
            pass
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([3,0,1,2,3,0,1,3,3,0,1,4,3,0,1,5,3,0,1,6])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,8,12,16,20])))
        self.assertRaises(InterpKernelException,m.insertNextCell,NORM_TRI3,None)
        pass

    def testSwigCurveLinearMesh1(self):
        m=MEDCouplingCurveLinearMesh("toto")
        m.setNodeGridStructure([2,3])
        coords=DataArrayDouble([0.,0., 2.,0., 0.,1., 1.9,1.1, 0.3,1.9, 2.2,2.1],6,2)
        m.setCoords(coords)
        m.checkConsistencyLight()
        m0=m.deepCopy()
        self.assertTrue(m0.isEqual(m,1e-12))
        m.getCoords().setInfoOnComponents(["X [m]","Y [m]"])
        self.assertTrue(not m0.isEqual(m,1e-12))
        m0=m.deepCopy()
        self.assertTrue(m0.isEqual(m,1e-12))
        self.assertEqual(m.getNodeGridStructure(),(2,3))
        pass

    def testSimplexize3(self):
        m=MEDCouplingUMesh("toto",3)
        m.allocateCells(0)
        m.insertNextCell(NORM_TETRA4,[0,1,2,3])
        self.assertEqual([NORM_TETRA4],m.getAllGeoTypesSorted())
        m.insertNextCell(NORM_HEXA8,[4,5,6,7,8,9,10,11])
        self.assertEqual([NORM_TETRA4,NORM_HEXA8],m.getAllGeoTypesSorted())
        m.insertNextCell(NORM_HEXA8,[12,13,14,15,16,17,18,19])
        self.assertEqual([NORM_TETRA4,NORM_HEXA8],m.getAllGeoTypesSorted())
        m.insertNextCell(NORM_TETRA4,[20,21,22,23])
        self.assertEqual([NORM_TETRA4,NORM_HEXA8,NORM_TETRA4],m.getAllGeoTypesSorted())
        c1=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,0.,0.,0.,0.,1.],4,3)
        c2=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,1.,0.,1.,0.,0., 0.,0.,1.,0.,1.,1.,1.,1.,1.,1.,0.,1.],8,3) ; c2+=[2.,0.,0.]
        c3=c2+[2.,0.,0.]
        c4=c1+[6.,0.,0.]
        c=DataArrayDouble.Aggregate([c1,c2,c3,c4])
        m.setCoords(c)
        m.checkConsistency()
        #
        m1=m.deepCopy()
        d1=m1.simplexize(PLANAR_FACE_5)
        m1.checkConsistency()
        vol1=m1.getMeasureField(False).getArray()
        self.assertTrue(vol1.isEqual(DataArrayDouble([1./6, 1./6, 1./6,1./6, 1./6, 1./3,1./6, 1./6, 1./6, 1./6, 1./3, 1./6]),1e-12))
        self.assertEqual(m1.getNodalConnectivity().getValues(),[14,0,1,2,3,14,4,9,5,6,14,4,8,9,11,14,4,7,11,6,14,9,11,10,6,14,4,9,6,11,14,12,17,13,14,14,12,16,17,19,14,12,15,19,14,14,17,19,18,14,14,12,17,14,19,14,20,21,22,23])
        self.assertEqual(m1.getNodalConnectivityIndex().getValues(),[0,5,10,15,20,25,30,35,40,45,50,55,60])
        self.assertTrue(d1.isEqual(DataArrayInt([0,1,1,1,1,1,2,2,2,2,2,3])))
        #
        m2=m.deepCopy()
        d2=m2.simplexize(PLANAR_FACE_6)
        m2.checkConsistency()
        vol2=m2.getMeasureField(False).getArray()
        self.assertTrue(vol2.isEqual(DataArrayDouble([1./6, 1./6, 1./6,1./6, 1./6, 1./6,1./6,1./6, 1./6, 1./6, 1./6, 1./6,1./6,1./6]),1e-12))
        self.assertEqual(m2.getNodalConnectivity().getValues(),[14,0,1,2,3,14,4,9,5,10,14,4,5,6,10,14,4,8,9,10,14,4,11,8,10,14,4,6,7,10,14,4,7,11,10,14,12,17,13,18,14,12,13,14,18,14,12,16,17,18,14,12,19,16,18,14,12,14,15,18,14,12,15,19,18,14,20,21,22,23])
        self.assertEqual(m2.getNodalConnectivityIndex().getValues(),[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70])
        self.assertTrue(d2.isEqual(DataArrayInt([0,1,1,1,1,1,1,2,2,2,2,2,2,3])))
        pass

    def testSwig2CurveLinearMesh2(self):
        c=MEDCouplingCMesh()
        #2D
        arr1=DataArrayDouble([0,1,3,7])
        arr2=DataArrayDouble([0,1,1.5])
        c.setCoords(arr1,arr2)
        u=c.buildUnstructured()
        coo=u.getCoords()
        cl=MEDCouplingCurveLinearMesh()
        cl.setCoords(coo)
        cl.setNodeGridStructure([4,3])
        cl.checkConsistency()
        li1=[1.,2.,4.,0.5,1.,2.]
        self.assertTrue(cl.getMeasureField(False).getArray().isEqual(DataArrayDouble(li1),1e-14))
        self.assertTrue(u.getMeasureField(False).getArray().isEqual(DataArrayDouble(li1),1e-14))
        li1_1=[0.5,0.5,2.,0.5,5.,0.5,0.5,1.25,2.,1.25,5.,1.25]
        self.assertTrue(cl.computeCellCenterOfMass().isEqual(DataArrayDouble(li1_1,6,2),1e-14))
        self.assertTrue(u.computeCellCenterOfMass().isEqual(DataArrayDouble(li1_1,6,2),1e-14))
        #3D
        c.setCoords(arr1,arr2,arr2)
        u=c.buildUnstructured()
        coo=u.getCoords()
        cl=MEDCouplingCurveLinearMesh()
        cl.setCoords(coo)
        cl.setNodeGridStructure([4,3,3])
        cl.checkConsistency()
        li2=[1.,2.,4.,0.5, 1.,2.,0.5,1.,2.,0.25,0.5,1.]
        li2_1=[0.5,0.5,0.5,2.,0.5,0.5,5.,0.5,0.5,0.5,1.25,0.5,2.,1.25,0.5,5.,1.25,0.5,0.5,0.5,1.25,2.,0.5,1.25,5.,0.5,1.25,0.5,1.25,1.25,2.,1.25,1.25,5.,1.25,1.25]
        self.assertTrue(cl.getMeasureField(False).getArray().isEqual(DataArrayDouble(li2),1e-14))
        self.assertTrue(u.getMeasureField(False).getArray().isEqual(DataArrayDouble(li2),1e-14))
        self.assertTrue(cl.computeCellCenterOfMass().isEqual(DataArrayDouble(li2_1,12,3),1e-14))
        self.assertTrue(u.computeCellCenterOfMass().isEqual(DataArrayDouble(li2_1,12,3),1e-14))
        #1D spaceDim 1
        coo=DataArrayDouble(5) ; coo.iota(0.)
        coo=coo*coo
        cl.setCoords(coo)
        cl.setNodeGridStructure([5])
        cl.checkConsistency()
        li3=[1.,3.,5.,7.]
        li3_1=[0.5,2.5,6.5,12.5]
        self.assertTrue(cl.getMeasureField(False).getArray().isEqual(DataArrayDouble(li3),1e-14))
        self.assertTrue(cl.buildUnstructured().getMeasureField(False).getArray().isEqual(DataArrayDouble(li3),1e-14))
        self.assertTrue(cl.computeCellCenterOfMass().isEqual(DataArrayDouble(li3_1),1e-14))
        self.assertTrue(cl.buildUnstructured().computeCellCenterOfMass().isEqual(DataArrayDouble(li3_1),1e-14))
        #1D spaceDim 2
        coo=DataArrayDouble.Meld(coo,coo)
        cl.setCoords(coo)
        cl.checkConsistency()
        li4=[sqrt(2.)*elt for elt in [1.,3.,5.,7.]]
        li4_1=[0.5,0.5,2.5,2.5,6.5,6.5,12.5,12.5]
        self.assertEqual(2,cl.getSpaceDimension())
        self.assertEqual(1,cl.getMeshDimension())
        self.assertEqual(4,cl.getNumberOfCells())
        self.assertEqual(5,cl.getNumberOfNodes())
        self.assertTrue(cl.getMeasureField(False).getArray().isEqual(DataArrayDouble(li4),1e-14))
        self.assertTrue(cl.buildUnstructured().getMeasureField(False).getArray().isEqual(DataArrayDouble(li4),1e-14))
        self.assertTrue(cl.computeCellCenterOfMass().isEqual(DataArrayDouble(li4_1,4,2),1e-14))
        self.assertTrue(cl.buildUnstructured().computeCellCenterOfMass().isEqual(DataArrayDouble(li4_1,4,2),1e-14))
        pass

    def testSwig2CurveLinearMeshNonRegression1(self):
        coords=DataArrayDouble([0.0, 0.0, 0.10000000149011612, 0.6000000238418579, 0.10000000149011612, 0.30000001192092896, 1.100000023841858, 0.10000000149011612, 0.20000000298023224, 0.10000000149011612, 0.6000000238418579, 0.20000000298023224, 0.699999988079071, 0.6000000238418579, 0.10000000149011612, 1.2000000476837158, 0.6000000238418579, 0.30000001192092896, 0.10000000149011612, 1.100000023841858, 0.30000001192092896, 0.5, 1.100000023841858, 0.20000000298023224, 1.0, 1.2000000476837158, 0.10000000149011612, 0.0, 0.10000000149011612, 0.5, 0.5, 0.10000000149011612, 0.6000000238418579, 1.2000000476837158, 0.10000000149011612, 0.699999988079071, 0.10000000149011612, 0.6000000238418579, 0.699999988079071, 0.6000000238418579, 0.6000000238418579, 0.5, 1.100000023841858, 0.6000000238418579, 0.6000000238418579, 0.10000000149011612, 1.0, 0.6000000238418579, 0.699999988079071, 1.2000000476837158, 0.699999988079071, 0.8999999761581421, 1.0, 0.5, 0.10000000149011612, 0.10000000149011612, 1.2000000476837158, 0.699999988079071, 0.10000000149011612, 1.0, 1.0, 0.10000000149011612, 1.100000023841858, 0.10000000149011612, 0.6000000238418579, 1.100000023841858, 0.6000000238418579, 0.6000000238418579, 1.100000023841858, 1.100000023841858, 0.6000000238418579, 1.2000000476837158, 0.10000000149011612, 1.2000000476837158, 1.0, 0.5, 1.100000023841858, 1.2000000476837158, 1.2000000476837158, 1.100000023841858, 1.0],27,3)
        m=MEDCouplingCurveLinearMesh("toto")
        m.setCoords(coords)
        m.setNodeGridStructure([3,3,3])
        #
        vol=m.getMeasureField(False).getArray()
        self.assertTrue(vol.isEqual(DataArrayDouble([0.11450000709295281, 0.10583334351579375,0.11149999939029423,0.08866666863113633, 0.1404166805123294,0.1250000135352219,0.1270833433481557,0.13258334288001067]),1e-12))
        self.assertTrue(vol.isEqual(m.buildUnstructured().getMeasureField(False).getArray(),1e-12))
        #
        self.assertTrue(m.computeCellCenterOfMass().isEqual(m.buildUnstructured().computeCellCenterOfMass(),1e-12))
        pass

    def testSwig2NonRegressionDASetSelectedComponents1(self):
        da=DataArrayDouble.New([1.,2.,3.,4.,5.,6.],3,2)
        dv=DataArrayDouble.New();
        dv.alloc(4,4)
        dv.fillWithZero()
        # da has less tuples than dv
        dv.setSelectedComponents(da,[1,0])
        #
        self.assertTrue(dv.isEqual(DataArrayDouble([2.,1.,0.,0.,4.,3.,0.,0.,6.,5.,0.,0.,0.,0.,0.,0.],4,4),1e-14))
        #
        da=DataArrayInt.New([1,2,3,4,5,6],3,2)
        dv=DataArrayInt.New();
        dv.alloc(4,4)
        dv.fillWithZero()
        # da has less tuples than dv
        dv.setSelectedComponents(da,[1,0])
        #
        self.assertTrue(dv.isEqual(DataArrayInt([2,1,0,0,4,3,0,0,6,5,0,0,0,0,0,0],4,4)))
        pass

    def testSwigSetItem3(self):
        # 1-2 
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[3]=[1,2]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,0,0,0,0,1,2,0,0,0,0],6,2),1e-14))
        # 2-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[5,3,2]]=[1,2]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,0,0,1,2,1,2,0,0,1,2],6,2),1e-14))
        # 3-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[:]=[1,2]
        self.assertTrue(d.isEqual(DataArrayDouble([1,2,1,2,1,2,1,2,1,2,1,2],6,2),1e-14))
        # 4-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[DataArrayInt([0,3,4])]=[1,2]
        self.assertTrue(d.isEqual(DataArrayDouble([1,2,0,0,0,0,1,2,1,2,0,0],6,2),1e-14))
        # 5-2
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[5,1]=[7]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,7],6,2),1e-14))
        # 6-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[3,5],1]=[7]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,0,0,0,0,0,7,0,0,0,7],6,2),1e-14))
        # 7-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[:-1:2,1]=[7]
        self.assertTrue(d.isEqual(DataArrayDouble([0,7,0,0,0,7,0,0,0,7,0,0],6,2),1e-14))
        # 8-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[DataArrayInt([0,3,4]),1]=[7]
        self.assertTrue(d.isEqual(DataArrayDouble([0,7,0,0,0,0,0,7,0,7,0,0],6,2),1e-14))
        # 9-2
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[3,[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,0,0,0,0,8,7,0,0,0,0],6,2),1e-14))
        # 10-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[1,3,4],[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,8,7,0,0,8,7,8,7,0,0],6,2),1e-14))
        # 11-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[1::2,[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,8,7,0,0,8,7,0,0,8,7],6,2),1e-14))
        # 12-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[DataArrayInt([1,4]),[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,8,7,0,0,0,0,8,7,0,0],6,2),1e-14))
        # 13-2
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[1,:-1]=[9]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,9,0,0,0,0,0,0,0,0,0],6,2),1e-14))
        # 14-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[1,4,5],:]=[7,8]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,7,8,0,0,0,0,7,8,7,8],6,2),1e-14))
        # 15-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[1::2,:]=[3,9]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,3,9,0,0,3,9,0,0,3,9],6,2),1e-14))
        # 1-2 
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[3]=[1,2]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,0,0,0,0,1,2,0,0,0,0],6,2)))
        # 2-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[5,3,2]]=[1,2]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,0,0,1,2,1,2,0,0,1,2],6,2)))
        # 3-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[:]=[1,2]
        self.assertTrue(d.isEqual(DataArrayInt([1,2,1,2,1,2,1,2,1,2,1,2],6,2)))
        # 4-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[DataArrayInt([0,3,4])]=[1,2]
        self.assertTrue(d.isEqual(DataArrayInt([1,2,0,0,0,0,1,2,1,2,0,0],6,2)))
        # 5-2
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[5,1]=[7]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,7],6,2)))
        # 6-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[3,5],1]=[7]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,0,0,0,0,0,7,0,0,0,7],6,2)))
        # 7-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[:-1:2,1]=[7]
        self.assertTrue(d.isEqual(DataArrayInt([0,7,0,0,0,7,0,0,0,7,0,0],6,2)))
        # 8-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[DataArrayInt([0,3,4]),1]=[7]
        self.assertTrue(d.isEqual(DataArrayInt([0,7,0,0,0,0,0,7,0,7,0,0],6,2)))
        # 9-2
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[3,[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,0,0,0,0,8,7,0,0,0,0],6,2)))
        # 10-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[1,3,4],[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,8,7,0,0,8,7,8,7,0,0],6,2)))
        # 11-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[1::2,[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,8,7,0,0,8,7,0,0,8,7],6,2)))
        # 12-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[DataArrayInt([1,4]),[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,8,7,0,0,0,0,8,7,0,0],6,2)))
        # 13-2
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[1,:-1]=[9]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,9,0,0,0,0,0,0,0,0,0],6,2)))
        # 14-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[1,4,5],:]=[7,8]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,7,8,0,0,0,0,7,8,7,8],6,2)))
        # 15-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[1::2,:]=[3,9]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,3,9,0,0,3,9,0,0,3,9],6,2)))
        pass

    def testSwig2ConvertLinearCellsToQuadratic1(self):
        coordsExp=DataArrayDouble([-0.3,-0.3,0.2,-0.3,0.7,-0.3,-0.3,0.2,0.2,0.2,0.7,0.2,-0.3,0.7,0.2,0.7,0.7,0.7,-0.3,-0.05,-0.05,0.2,0.2,-0.05,-0.05,-0.3,0.45,-0.05,0.45,-0.3,0.45,0.2,0.7,-0.05,-0.05,0.7,0.2,0.45,-0.3,0.45,0.45,0.7,0.7,0.45],22,2)
        # 2D
        m2D=MEDCouplingDataForTest.build2DTargetMesh_1()
        m2D.convertLinearCellsToQuadratic(0)
        m2D.checkConsistency()
        self.assertEqual(m2D.getNodalConnectivity().getValues(),[8,0,3,4,1,9,10,11,12,6,1,4,2,11,13,14,6,4,5,2,15,16,13,8,6,7,4,3,17,18,10,19,8,7,8,5,4,20,21,15,18])
        self.assertEqual(m2D.getNodalConnectivityIndex().getValues(),[0,9,16,23,32,41])
        self.assertTrue(m2D.getCoords().isEqual(coordsExp,1e-14))
        # 1D
        m1D=MEDCouplingDataForTest.build2DTargetMesh_1().buildDescendingConnectivity()[0]
        m1D.convertLinearCellsToQuadratic(0)
        m1D.checkConsistency()
        self.assertEqual(m1D.getNodalConnectivity().getValues(),[2,0,3,9,2,3,4,10,2,4,1,11,2,1,0,12,2,4,2,13,2,2,1,14,2,4,5,15,2,5,2,16,2,6,7,17,2,7,4,18,2,3,6,19,2,7,8,20,2,8,5,21])
        self.assertEqual(m1D.getNodalConnectivityIndex().getValues(),[0,4,8,12,16,20,24,28,32,36,40,44,48,52])
        self.assertTrue(m1D.getCoords().isEqual(coordsExp,1e-14))
        # 3D
        m2D=MEDCouplingDataForTest.build2DTargetMesh_1()
        m2D.changeSpaceDimension(3)
        arr=DataArrayDouble(4);  arr.iota(0) ; z=MEDCouplingCMesh() ; z.setCoords(arr)
        m1D=z.buildUnstructured() ; m1D.setCoords(arr.changeNbOfComponents(3,0.))
        m1D.getCoords()[:]=m1D.getCoords()[:,[1,2,0]]
        cooTmp=m2D.getCoords()[:]
        m3D=m2D.buildExtrudedMesh(m1D,0)
        m3D.convertLinearCellsToQuadratic(0)
        m3D.checkConsistency()
        # check of new m3D content
        coordsExp2 = [coordsExp.changeNbOfComponents(3, i) for i in range(4)]
        coordsExp3 = [DataArrayDouble.Meld(cooTmp[:, [0, 1]], cooTmp[:, 2] + (0.5 + float(i))) for i in range(3)]
        coordsExp4=DataArrayDouble.Aggregate([coordsExp2[0],coordsExp3[0],coordsExp2[1],coordsExp3[1],coordsExp2[2],coordsExp3[2],coordsExp2[3]])
        c=DataArrayDouble.Aggregate(m3D.getCoords(),coordsExp4)
        self.assertEqual(len(coordsExp4),115)
        self.assertEqual(len(m3D.getCoords()),115)
        a,b=c.findCommonTuples(1e-14)
        self.assertEqual(len(b),len(coordsExp4)+1)
        e,f=DataArrayInt.ConvertIndexArrayToO2N(2*115,a,b)
        self.assertEqual(f,115)
        self.assertTrue(e.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,0,1,2,3,4,5,6,7,8,36,37,38,39,48,49,53,54,58,59,60,66,67,44,47,52,45,46,57,64,65,70,9,10,11,12,13,14,15,16,17,40,41,42,43,50,51,55,56,61,62,63,68,69,75,78,81,76,77,84,88,89,92,18,19,20,21,22,23,24,25,26,71,72,73,74,79,80,82,83,85,86,87,90,91,97,100,103,98,99,106,110,111,114,27,28,29,30,31,32,33,34,35,93,94,95,96,101,102,104,105,107,108,109,112,113])))
        self.assertTrue(DataArrayInt([30,0,3,4,1,9,12,13,10,36,37,38,39,40,41,42,43,44,45,46,47,25,1,4,2,10,13,11,38,48,49,42,50,51,47,46,52,25,4,5,2,13,14,11,53,54,48,55,56,50,46,57,52,30,6,7,4,3,15,16,13,12,58,59,37,60,61,62,41,63,64,65,46,45,30,7,8,5,4,16,17,14,13,66,67,53,59,68,69,55,62,65,70,57,46,30,9,12,13,10,18,21,22,19,40,41,42,43,71,72,73,74,75,76,77,78,25,10,13,11,19,22,20,42,50,51,73,79,80,78,77,81,25,13,14,11,22,23,20,55,56,50,82,83,79,77,84,81,30,15,16,13,12,24,25,22,21,61,62,41,63,85,86,72,87,88,89,77,76,30,16,17,14,13,25,26,23,22,68,69,55,62,90,91,82,86,89,92,84,77,30,18,21,22,19,27,30,31,28,71,72,73,74,93,94,95,96,97,98,99,100,25,19,22,20,28,31,29,73,79,80,95,101,102,100,99,103,25,22,23,20,31,32,29,82,83,79,104,105,101,99,106,103,30,24,25,22,21,33,34,31,30,85,86,72,87,107,108,94,109,110,111,99,98,30,25,26,23,22,34,35,32,31,90,91,82,86,112,113,104,108,111,114,106,99]).isEqual(m3D.getNodalConnectivity()))
        self.assertTrue(DataArrayInt([0,21,37,53,74,95,116,132,148,169,190,211,227,243,264,285]).isEqual(m3D.getNodalConnectivityIndex()))
        # testing explode3DMeshTo1D
        m3DSlice0=m3D[:5]
        m3DSlice0.zipCoords()
        a,b,c,d,e=m3DSlice0.explode3DMeshTo1D()
        self.assertTrue(b.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,2,12,13,6,14,15,11,10,16,17,18,12,19,20,14,10,21,16,22,23,1,24,25,26,5,27,28,29,10,9,30,31,17,23,32,33,19,26,29,34,21,10])))
        self.assertTrue(c.isEqual(DataArrayInt([0,12,21,30,42,54])))
        self.assertTrue(d.isEqual(DataArrayInt([0,0,3,0,1,0,0,0,3,0,1,0,0,0,3,0,1,2,3,4,0,1,1,2,1,1,2,1,1,2,2,4,2,2,4,2,2,4,3,3,4,3,3,3,4,3,3,3,4,4,4,4,4,4])))
        self.assertTrue(e.isEqual(DataArrayInt([0,1,3,5,6,7,9,11,12,13,15,20,22,24,25,27,28,30,32,33,35,36,38,39,41,42,43,45,46,47,49,50,51,52,53,54])))
        self.assertTrue(a.getNodalConnectivity().isEqual(DataArrayInt([2,0,3,18,2,3,4,19,2,4,1,20,2,1,0,21,2,9,12,22,2,12,13,23,2,13,10,24,2,10,9,25,2,0,9,26,2,3,12,27,2,4,13,28,2,1,10,29,2,4,2,30,2,2,1,31,2,13,11,32,2,11,10,33,2,2,11,34,2,4,5,35,2,5,2,36,2,13,14,37,2,14,11,38,2,5,14,39,2,6,7,40,2,7,4,41,2,3,6,42,2,15,16,43,2,16,13,44,2,12,15,45,2,6,15,46,2,7,16,47,2,7,8,48,2,8,5,49,2,16,17,50,2,17,14,51,2,8,17,52])))
        self.assertTrue(a.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,132,136,140])))
        self.assertTrue(a.getCoords().isEqual(DataArrayDouble([-0.3,-0.3,0.0,0.2,-0.3,0.0,0.7,-0.3,0.0,-0.3,0.2,0.0,0.2,0.2,0.0,0.7,0.2,0.0,-0.3,0.7,0.0,0.2,0.7,0.0,0.7,0.7,0.0,-0.3,-0.3,1.0,0.2,-0.3,1.0,0.7,-0.3,1.0,-0.3,0.2,1.0,0.2,0.2,1.0,0.7,0.2,1.0,-0.3,0.7,1.0,0.2,0.7,1.0,0.7,0.7,1.0,-0.3,-0.05,0.0,-0.05,0.2,0.0,0.2,-0.05,0.0,-0.05,-0.3,0.0,-0.3,-0.05,1.0,-0.05,0.2,1.0,0.2,-0.05,1.0,-0.05,-0.3,1.0,-0.3,-0.3,0.5,-0.3,0.2,0.5,0.2,0.2,0.5,0.2,-0.3,0.5,0.45,-0.05,0.0,0.45,-0.3,0.0,0.45,-0.05,1.0,0.45,-0.3,1.0,0.7,-0.3,0.5,0.45,0.2,0.0,0.7,-0.05,0.0,0.45,0.2,1.0,0.7,-0.05,1.0,0.7,0.2,0.5,-0.05,0.7,0.0,0.2,0.45,0.0,-0.3,0.45,0.0,-0.05,0.7,1.0,0.2,0.45,1.0,-0.3,0.45,1.0,-0.3,0.7,0.5,0.2,0.7,0.5,0.45,0.7,0.0,0.7,0.45,0.0,0.45,0.7,1.0,0.7,0.45,1.0,0.7,0.7,0.5],53,3),1e-14))
        pass

    def testSwig2DataArrayPushBackValsSilent1(self):
        d=DataArrayDouble()
        d.pushBackValsSilent([4,5,6])
        self.assertTrue(d.isEqual(DataArrayDouble([4.,5.,6.]),1e-14))
        e=DataArrayDouble([1,2,3],1,3)
        for t in e: d.pushBackValsSilent(t)
        self.assertTrue(d.isEqual(DataArrayDouble([4.,5.,6.,1.,2.,3.]),1e-14))
        d.pushBackValsSilent(DataArrayDouble([9,10.]))
        self.assertTrue(d.isEqual(DataArrayDouble([4.,5.,6.,1.,2.,3.,9.,10.]),1e-14))
        d.pushBackValsSilent(DataArrayDouble(0,1))
        self.assertTrue(d.isEqual(DataArrayDouble([4.,5.,6.,1.,2.,3.,9.,10.]),1e-14))
        e=DataArrayDouble([1,2,3],3,1)
        for t in e: d.pushBackValsSilent(t)
        self.assertTrue(d.isEqual(DataArrayDouble([4.,5.,6.,1.,2.,3.,9.,10.,1.,2.,3.]),1e-14))
        d.pushBackValsSilent(77)
        self.assertTrue(d.isEqual(DataArrayDouble([4.,5.,6.,1.,2.,3.,9.,10.,1.,2.,3.,77.]),1e-14))
        #
        d=DataArrayInt()
        d.pushBackValsSilent([4,5,6])
        self.assertTrue(d.isEqual(DataArrayInt([4,5,6])))
        e=DataArrayInt([1,2,3],1,3)
        for t in e: d.pushBackValsSilent(t)
        self.assertTrue(d.isEqual(DataArrayInt([4,5,6,1,2,3])))
        d.pushBackValsSilent(DataArrayInt([9,10]))
        self.assertTrue(d.isEqual(DataArrayInt([4,5,6,1,2,3,9,10])))
        d.pushBackValsSilent(DataArrayInt(0,1))
        self.assertTrue(d.isEqual(DataArrayInt([4,5,6,1,2,3,9,10])))
        e=DataArrayInt([1,2,3],3,1)
        for t in e: d.pushBackValsSilent(t)
        self.assertTrue(d.isEqual(DataArrayInt([4,5,6,1,2,3,9,10,1,2,3])))
        d.pushBackValsSilent(77)
        self.assertTrue(d.isEqual(DataArrayInt([4,5,6,1,2,3,9,10,1,2,3,77])))
        pass

    def testSwig2ConvertLinearCellsToQuadratic2(self):
        m2D=MEDCouplingDataForTest.build2DTargetMesh_1()
        ret=m2D.convertLinearCellsToQuadratic(1)
        self.assertTrue(ret.isIota(5))
        m2D.checkConsistency()
        coordsExp=DataArrayDouble([-0.3,-0.3,0.2,-0.3,0.7,-0.3,-0.3,0.2,0.2,0.2,0.7,0.2,-0.3,0.7,0.2,0.7,0.7,0.7,-0.3,-0.05,-0.05,0.2,0.2,-0.05,-0.05,-0.3,0.45,-0.05,0.45,-0.3,0.45,0.2,0.7,-0.05,-0.05,0.7,0.2,0.45,-0.3,0.45,0.45,0.7,0.7,0.45,-0.05,-0.05,0.3666666666666667,-0.1333333333333333,0.5333333333333332,0.03333333333333334,-0.05,0.45,0.45,0.45],27,2)
        self.assertTrue(m2D.getCoords().isEqual(coordsExp,1e-14))
        self.assertTrue(m2D.getNodalConnectivity().isEqual(DataArrayInt([9,0,3,4,1,9,10,11,12,22,7,1,4,2,11,13,14,23,7,4,5,2,15,16,13,24,9,6,7,4,3,17,18,10,19,25,9,7,8,5,4,20,21,15,18,26])))
        self.assertTrue(m2D.getNodalConnectivityIndex().isEqual(DataArrayInt([0,10,18,26,36,46])))
        #
        m2D=MEDCouplingDataForTest.build2DTargetMesh_1()[(0,3)] ; m2D.zipCoords()
        m2D.changeSpaceDimension(3)
        arr=DataArrayDouble(3);  arr.iota(0) ; z=MEDCouplingCMesh() ; z.setCoords(arr)
        m1D=z.buildUnstructured() ; m1D.setCoords(arr.changeNbOfComponents(3,0.))
        m1D.getCoords()[:]=m1D.getCoords()[:,[1,2,0]]
        cooTmp=m2D.getCoords()[:]
        m3D=m2D.buildExtrudedMesh(m1D,0)
        ret=m3D.convertLinearCellsToQuadratic(1)
        self.assertTrue(ret.isIota(4))
        m3D.checkConsistency()
        coordsExp2=DataArrayDouble([-0.3,-0.3,0.0,0.2,-0.3,0.0,-0.3,0.2,0.0,0.2,0.2,0.0,-0.3,0.7,0.0,0.2,0.7,0.0,-0.3,-0.3,1.0,0.2,-0.3,1.0,-0.3,0.2,1.0,0.2,0.2,1.0,-0.3,0.7,1.0,0.2,0.7,1.0,-0.3,-0.3,2.0,0.2,-0.3,2.0,-0.3,0.2,2.0,0.2,0.2,2.0,-0.3,0.7,2.0,0.2,0.7,2.0,-0.3,-0.05,0.0,-0.05,0.2,0.0,0.2,-0.05,0.0,-0.05,-0.3,0.0,-0.3,-0.05,1.0,-0.05,0.2,1.0,0.2,-0.05,1.0,-0.05,-0.3,1.0,-0.3,-0.3,0.5,-0.3,0.2,0.5,0.2,0.2,0.5,0.2,-0.3,0.5,-0.05,0.7,0.0,0.2,0.45,0.0,-0.3,0.45,0.0,-0.05,0.7,1.0,0.2,0.45,1.0,-0.3,0.45,1.0,-0.3,0.7,0.5,0.2,0.7,0.5,-0.3,-0.05,2.0,-0.05,0.2,2.0,0.2,-0.05,2.0,-0.05,-0.3,2.0,-0.3,-0.3,1.5,-0.3,0.2,1.5,0.2,0.2,1.5,0.2,-0.3,1.5,-0.05,0.7,2.0,0.2,0.45,2.0,-0.3,0.45,2.0,-0.3,0.7,1.5,0.2,0.7,1.5,-0.05,-0.05,0.0,-0.3,-0.05,0.5,-0.05,0.2,0.5,0.2,-0.05,0.5,-0.05,-0.3,0.5,-0.05,-0.05,1.0,-0.05,0.45,0.0,-0.05,0.7,0.5,0.2,0.45,0.5,-0.3,0.45,0.5,-0.05,0.45,1.0,-0.3,-0.05,1.5,-0.05,0.2,1.5,0.2,-0.05,1.5,-0.05,-0.3,1.5,-0.05,-0.05,2.0,-0.05,0.7,1.5,0.2,0.45,1.5,-0.3,0.45,1.5,-0.05,0.45,2.0,-0.05,-0.05,0.5,-0.05,0.45,0.5,-0.05,-0.05,1.5,-0.05,0.45,1.5],75,3)
        self.assertTrue(m3D.getCoords().isEqual(coordsExp2,1e-14))
        self.assertTrue(m3D.getNodalConnectivity().isEqual(DataArrayInt([27,0,2,3,1,6,8,9,7,18,19,20,21,22,23,24,25,26,27,28,29,51,52,53,54,55,56,71,27,4,5,3,2,10,11,9,8,30,31,19,32,33,34,23,35,36,37,28,27,57,58,59,53,60,61,72,27,6,8,9,7,12,14,15,13,22,23,24,25,38,39,40,41,42,43,44,45,56,62,63,64,65,66,73,27,10,11,9,8,16,17,15,14,33,34,23,35,46,47,39,48,49,50,44,43,61,67,68,63,69,70,74])))
        self.assertTrue(m3D.getNodalConnectivityIndex().isEqual(DataArrayInt([0,28,56,84,112])))
        pass

    def testSwig2GaussNEIntegral1(self):
        m2D=MEDCouplingDataForTest.build2DTargetMesh_1()
        m0=m2D[0] ; m0.zipCoords()
        m1=m2D[[1,2]] ; m1.zipCoords()
        m2=m2D[[3,4]] ; m2.zipCoords()
        m0.convertLinearCellsToQuadratic(1)
        m1.convertLinearCellsToQuadratic(0)
        m2.convertLinearCellsToQuadratic(1)
        m=MEDCouplingUMesh.MergeUMeshes([m0,m1,m2])
        m.mergeNodes(1e-12)
        f=MEDCouplingFieldDouble(ON_GAUSS_NE)
        f.setMesh(m)
        arr=DataArrayDouble([1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9,
                             11.1,12.2,13.3,14.4,15.5,16.6,
                             21.1,22.2,23.3,24.4,25.5,26.6,
                             31.1,32.2,33.3,34.4,35.5,36.6,37.7,38.8,39.9,
                             41.1,42.2,43.3,44.4,45.5,46.6,47.7,48.8,49.9])
        arr2=DataArrayDouble(len(arr),2)
        arr2[:,0]=arr ; arr2[:,1]=arr+100
        f.setArray(arr2)
        f.checkConsistencyLight()
        res=f.integral(False)
        # a=25./81 ; b=40./81 ; c=64./81
        # p1=0.11169079483905 ; p2=0.0549758718227661
        # 1st compo
        # c0=(a*(1.1+2.2+3.3+4.4)+b*(5.5+6.6+7.7+8.8)+c*9.9)*0.25/3.9999999999999978 ; c0=1.5837962962962973
        # c1=(p2*(11.1+12.2+13.3)+p1*(14.4+15.5+16.6))*0.125/0.4999999999854482 ; c1=1.8014347172346943
        # c2=(p2*(21.1+22.2+23.3)+p1*(24.4+25.5+26.6))*0.125/0.4999999999854482 ; c2=3.0514347172346943
        # c3=(a*(31.1+32.2+33.3+34.4)+b*(35.5+36.6+37.7+38.8)+c*39.9)*0.25/3.9999999999999978 ; c3=9.0837962962963
        # c4=(a*(41.1+42.2+43.3+44.4)+b*(45.5+46.6+47.7+48.8)+c*49.9)*0.25/3.9999999999999978 ; c4=11.583796296296303
        # c0+c1+c2+c3+c4=27.104258323358287
        integExp0=27.104258323358287
        self.assertAlmostEqual(res[0],integExp0,13)
        # 2nd compo
        # c0=(a*(101.1+102.2+103.3+104.4)+b*(105.5+106.6+107.7+108.8)+c*109.9)*0.25/3.9999999999999978 ; c0=26.58379629629631
        # c1=(p2*(111.1+112.2+113.3)+p1*(114.4+115.5+116.6))*0.125/0.4999999999854482 ; c1=14.301434717234699
        # c2=(p2*(121.1+122.2+123.3)+p1*(124.4+125.5+126.6))*0.125/0.4999999999854482 ; c2=15.5514347172347
        # c3=(a*(131.1+132.2+133.3+134.4)+b*(135.5+136.6+137.7+138.8)+c*139.9)*0.25/3.9999999999999978 ; c3=34.08379629629631
        # c4=(a*(141.1+142.2+143.3+144.4)+b*(145.5+146.6+147.7+148.8)+c*149.9)*0.25/3.9999999999999978 ; c4=36.58379629629632
        # c0+c1+c2+c3+c4=127.10425832335835
        integExp1=127.10425832335835
        self.assertAlmostEqual(res[1],integExp1,12)
        meas=f.getDiscretization().getMeasureField(f.getMesh(),False)
        intPerTuple=meas*f
        res2=intPerTuple.accumulate()
        self.assertAlmostEqual(res2[0],integExp0,13)
        self.assertAlmostEqual(res2[1],integExp1,12)
        #
        meas2=f.buildMeasureField(False)
        intPerTuple=meas2*f
        res3=intPerTuple.accumulate()
        self.assertAlmostEqual(res3[0],integExp0,13)
        self.assertAlmostEqual(res3[1],integExp1,12)
        #
        res4=f.getWeightedAverageValue(False) # res4==res2 because sum of area of mesh is equal to 1
        self.assertAlmostEqual(res4[0],integExp0,13)
        self.assertAlmostEqual(res4[1],integExp1,12)
        #
        m.scale([0,0],2.)
        #
        res5=f.getWeightedAverageValue() # res4==res4 because weighted average is not sensitive to the scaling
        self.assertAlmostEqual(res5[0],integExp0,13)
        self.assertAlmostEqual(res5[1],integExp1,12)
        meas3=f.buildMeasureField(False)
        delta=4*meas2.getArray()-meas3.getArray()
        delta.abs()
        self.assertTrue(delta.isUniform(0.,1e-16))
        res6=f.integral(False)
        self.assertAlmostEqual(res6[0],4.*integExp0,12)
        self.assertAlmostEqual(res6[1],4.*integExp1,11)
        pass

    def testSwig2SlowDADFindClosestTupleId(self):
        nbPts=[10,]
        for nbPt in nbPts:
            d=DataArrayDouble(nbPt) ; d.iota() ; d*=1./(nbPt-1)
            c=MEDCouplingCMesh() ; c.setCoords(d,d) ; m=c.buildUnstructured() ; pts=m.getCoords() ; del m
            #
            d0=DataArrayDouble((nbPt-1)*(nbPt-1)) ; d0.iota() ; d0*=(3./((nbPt-1)*(nbPt-1))) ; d0=d0.applyFunc("exp(x)-1")
            d1=DataArrayDouble((nbPt-1)*(nbPt-1)) ; d1.iota()
            d2=DataArrayDouble.Meld(d0,d1) ; d2=d2.fromPolarToCart() ; d2+=[0.32,0.73]
            ids=pts.findClosestTupleId(d2)
            #print "Start of costly computation"
            idsExpected=DataArrayInt(len(d2))
            tmp=1e300
            for i,elt in enumerate(d2):
                l,m=(pts-elt).magnitude().getMinValue()
                idsExpected.setIJSilent(i,0,m)
                if l<tmp:
                    tmp=l ; tmp1=m ; tmp2=i
                    pass
                pass
            #print "End of costly computation"
            self.assertTrue(idsExpected.isEqual(ids))
            a,b,c=pts.minimalDistanceTo(d2)
            self.assertEqual(tmp,a)
            self.assertEqual(tmp1,b)
            self.assertEqual(tmp2,c)
            #
            l=[d2[:,i] for i in [0,1]]
            for elt in l: elt.reverse()
            d2i=DataArrayDouble.Meld(l)
            ids1=pts.findClosestTupleId(d2i)
            idsExpectedI=idsExpected.deepCopy() ; idsExpectedI.reverse()
            self.assertTrue(idsExpectedI.isEqual(ids1))
            #
            l=[pts[:,i] for i in [0,1]]
            for elt in l: elt.reverse()
            ptsi=DataArrayDouble.Meld(l)
            ids2=ptsi.findClosestTupleId(d2)
            idsExpected2=nbPt*nbPt-1-ids
            self.assertTrue(idsExpected2.isEqual(ids2))
            #
            ids3=ptsi.findClosestTupleId(d2i)
            idsExpected3=idsExpected2.deepCopy() ; idsExpected3.reverse()
            self.assertTrue(idsExpected3.isEqual(ids3))
            pass

    def testSwig2DataArrayAsciiChar1(self):
        alpha=DataArrayInt(26) ; alpha.iota(ord("A"))
        d=DataArrayAsciiChar(alpha.getValues(),2,13)
        d.setInfoOnComponents(["c%i" % (v) for v in range(13)])
        self.assertEqual('ABCDEFGHIJKLM',d.getTuple(0))
        self.assertEqual('NOPQRSTUVWXYZ',d.getTuple(1))
        self.assertEqual(2,d.getNumberOfTuples())
        self.assertEqual(26,d.getNbOfElems())
        self.assertEqual(13,d.getNumberOfComponents())
        dd=d.deepCopy()
        self.assertTrue(d.isEqual(dd))
        dd.setIJ(0,3,'d')
        self.assertTrue(not d.isEqual(dd))
        d.setIJ(0,3,ord('d'))
        self.assertTrue(d.isEqual(dd))
        d.rearrange(1)
        d.reserve(20)
        self.assertEqual(20,d.getNumberOfTuples())
        self.assertEqual(20,d.getNbOfElems())
        self.assertEqual(1,d.getNumberOfComponents())
        #
        d0=DataArrayAsciiChar([ord('a')],1,1)
        self.assertEqual('a',d0.asciiCharValue())
        self.assertTrue(not d0.empty())
        d0=DataArrayAsciiChar(0,3)
        self.assertTrue(d0.empty())
        d.pushBackSilent("U") ; d.pushBackSilent("V") ; d.pushBackSilent("W")
        self.assertEqual("W",d.popBackSilent())
        d.rearrange(2)
        self.assertEqual(['AB','Cd','EF','GH','IJ','KL','MN','OP','QR','ST','UV'],d.toStrList())
        d.fillWithZero()
        self.assertEqual(11*[''],d.toStrList())
        d.fillWithValue('T')
        self.assertEqual(11*["TT"],d.toStrList())
        d.rearrange(1)
        self.assertTrue(d.isUniform("T"))
        d.rearrange(2)
        #
        dd.rearrange(2)
        dd2=dd.deepCopy()
        dd.renumberInPlace([3,1,2,4,0,11,10,9,8,7,5,12,6])
        self.assertEqual(dd.toStrList(),['IJ','Cd','EF','AB','GH','UV','YZ','ST','QR','OP','MN','KL','WX'])
        dd.renumberInPlaceR([3,1,2,4,0,11,10,9,8,7,5,12,6])
        self.assertEqual(['AB','Cd','EF','GH','IJ','KL','MN','OP','QR','ST','UV','WX','YZ'],dd.toStrList())
        e=dd.renumber([3,1,2,4,0,11,10,9,8,7,5,12,6])
        self.assertEqual(e.toStrList(),['IJ','Cd','EF','AB','GH','UV','YZ','ST','QR','OP','MN','KL','WX'])
        e=dd.renumberR([3,1,2,4,0,11,10,9,8,7,5,12,6])
        self.assertEqual(e.toStrList(),['GH','Cd','EF','IJ','AB','WX','UV','ST','QR','OP','KL','YZ','MN'])
        e=dd.renumberAndReduce([1,1,1,1,1,1,1,2,0,0,0,0,0],3)
        self.assertEqual(['YZ','MN','OP'],e.toStrList())
        self.assertEqual(['GH','IJ'],dd.selectByTupleIdSafe([3,4]).toStrList())
        self.assertEqual(['AB','GH','MN','ST','YZ'],dd.selectByTupleIdSafeSlice(0,13,3).toStrList())
        dd3=dd.changeNbOfComponents(3,"G")
        self.assertEqual(['ABG','CdG','EFG','GHG','IJG','KLG','MNG','OPG','QRG','STG','UVG','WXG','YZG'],dd3.toStrList())
        dd3.rearrange(1) ; self.assertEqual("G",dd3.back()) ; dd3.rearrange(3)
        self.assertTrue(dd3.changeNbOfComponents(2,"\0").isEqual(dd))
        self.assertEqual(len(dd),13)
        d=DataArrayAsciiChar(13,2) ; d.fillWithValue('Y')
        dd3.meldWith(d)
        self.assertEqual(['ABGYY','CdGYY','EFGYY','GHGYY','IJGYY','KLGYY','MNGYY','OPGYY','QRGYY','STGYY','UVGYY','WXGYY','YZGYY'],dd3.toStrList())
        self.assertEqual("d",dd3.getIJ(0,6))
        self.assertRaises(InterpKernelException,dd3.getIJSafe,0,6)
        self.assertEqual("d",dd3.getIJSafe(1,1))
        dd3.rearrange(1)
        e=dd3.findIdsEqual("Y")
        self.assertTrue(e.isEqual(DataArrayInt([3,4,8,9,13,14,18,19,23,24,28,29,33,34,38,39,43,44,48,49,53,54,58,59,60,63,64])))
        e=dd3.findIdsNotEqual("Y")
        self.assertTrue(e.isEqual(DataArrayInt([0,1,2,5,6,7,10,11,12,15,16,17,20,21,22,25,26,27,30,31,32,35,36,37,40,41,42,45,46,47,50,51,52,55,56,57,61,62])))
        self.assertEqual(("d",6),dd3.getMaxValue())
        self.assertEqual(("A",0),dd3.getMinValue())
        self.assertEqual(26,dd3.findIdSequence("LGYYM"))
        self.assertEqual(-1,dd3.findIdSequence("LGYYN"))
        dd3.rearrange(5)
        self.assertEqual(7,dd3.findIdFirstEqualTuple("OPGYY"))
        self.assertTrue("OPGYY" in dd3)
        self.assertEqual(7,dd3.index("OPGYY"))
        self.assertEqual(-1,dd3.findIdFirstEqualTuple("OPGYP"))
        dd3.rearrange(1)
        self.assertEqual(2,dd3.findIdFirstEqual("OPGYY"))
        self.assertTrue(dd3.presenceOfValue("OPGYY"))
        self.assertTrue("O" in dd3)
        self.assertTrue(not dd3.presenceOfValue("z"))
        self.assertTrue("z" not in dd3)
        dd3.rearrange(5)
        l=list(dd3)
        self.assertEqual([e.buildDAAsciiChar().toStrList()[0] for e in list(dd3)],dd3.toStrList())
        dd3.reAlloc(5)
        dd4=DataArrayChar.Aggregate(dd3,dd3)
        self.assertEqual(['ABGYY','CdGYY','EFGYY','GHGYY','IJGYY','ABGYY','CdGYY','EFGYY','GHGYY','IJGYY'],dd4.toStrList())
        dd5=DataArrayChar.Aggregate([dd4,dd3,dd4])
        self.assertEqual(['ABGYY','CdGYY','EFGYY','GHGYY','IJGYY','ABGYY','CdGYY','EFGYY','GHGYY','IJGYY','ABGYY','CdGYY','EFGYY','GHGYY','IJGYY','ABGYY','CdGYY','EFGYY','GHGYY','IJGYY','ABGYY','CdGYY','EFGYY','GHGYY','IJGYY'],dd5.toStrList())
        # getitem,__iter__,__setitem__
        a=list(dd3)
        self.assertEqual("ABGYY",str(a[0]))
        dd4=dd3[::2]
        self.assertEqual(['ABGYY','EFGYY','IJGYY'],dd4.toStrList())
        dd4=dd3[(3,2,1)]
        self.assertEqual(['GHGYY','EFGYY','CdGYY'],dd4.toStrList())
        dd4=dd3[:]
        dd4[::2]=["12","345","67890"]
        self.assertEqual(['12   ','CdGYY','345  ','GHGYY','67890'],dd4.toStrList())
        dd4=dd3[:]
        dd4[[1,2]]=" "
        self.assertEqual(['ABGYY','     ','     ','GHGYY','IJGYY'],dd4.toStrList())
        dd4=dd3[:]
        dd4[4]='12345'
        self.assertEqual(['ABGYY','CdGYY','EFGYY','GHGYY','12345'],dd4.toStrList())
        dd4[0]=dd4[1]
        self.assertEqual(['CdGYY','CdGYY','EFGYY','GHGYY','12345'],dd4.toStrList())
        dd4=DataArrayAsciiChar(["abc","de","fghi"])
        self.assertEqual(['abc ','de  ','fghi'],dd4.toStrList())
        dd4=DataArrayAsciiChar(["abc","de","fghi"],"t")
        self.assertEqual(['abct','dett','fghi'],dd4.toStrList())
        pass

    def testSwig2GaussNELocalizationOfDiscValues(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_3()[[0,1,3,4,5,6,8,9]] # suppression of polygons
        f=MEDCouplingFieldDouble(ON_GAUSS_NE)
        f.setMesh(m)
        loc=f.getLocalizationOfDiscr()
        self.assertEqual(42,len(loc))
        self.assertTrue(loc.isEqual(DataArrayDouble([0.,0.,1.,0.,0.5,1.,0.,0.,1.,0.,1.,1.,0.,1.,0.,0.,1.,0.,0.5,1.,0.5,0.,0.75,0.5,0.25,0.5,0.,0.,1.,0.,1.,1.,0.,1.,0.5,0.,1.,0.5,0.5,1.,0.,0.5,0.,0.,0.5,1.,1.,0.,0.,0.,0.,1.,1.,1.,1.,0.,0.,0.,0.5,1.,1.,0.,0.25,0.5,0.75,0.5,0.5,0.,0.,0.,0.,1.,1.,1.,1.,0.,0.,0.5,0.5,1.,1.,0.5,0.5,0.],42,2),1e-13))
        m.changeSpaceDimension(3)
        m.getCoords()[:,2]=7.
        loc=f.getLocalizationOfDiscr()
        self.assertEqual(42,len(loc))
        self.assertTrue(loc.isEqual(DataArrayDouble([0.,0.,7.,1.,0.,7.,0.5,1.,7.,0.,0.,7.,1.,0.,7.,1.,1.,7.,0.,1.,7.,0.,0.,7.,1.,0.,7.,0.5,1.,7.,0.5,0.,7.,0.75,0.5,7.,0.25,0.5,7.,0.,0.,7.,1.,0.,7.,1.,1.,7.,0.,1.,7.,0.5,0.,7.,1.,0.5,7.,0.5,1.,7.,0.,0.5,7.,0.,0.,7.,0.5,1.,7.,1.,0.,7.,0.,0.,7.,0.,1.,7.,1.,1.,7.,1.,0.,7.,0.,0.,7.,0.5,1.,7.,1.,0.,7.,0.25,0.5,7.,0.75,0.5,7.,0.5,0.,7.,0.,0.,7.,0.,1.,7.,1.,1.,7.,1.,0.,7.,0.,0.5,7.,0.5,1.,7.,1.,0.5,7.,0.5,0.,7.],42,3),1e-13))
        pass

    def testSwig2GaussMeasureAndIntegral(self):
        ft=MEDCouplingDataForTest.buildFieldOnGauss_1()
        mea=ft.buildMeasureField(False)
        mea.checkConsistencyLight()
        self.assertTrue(mea.getArray().isEqual(DataArrayDouble([-0.08504076274779823,-0.06378057206084897,-0.08504076274779869,-0.10630095343474463,-0.12756114412169625,-0.10630095343474734,-0.0637805720608491,-0.0850407627477968,-0.1063009534347449,-0.0850407627477994,-0.10630095343474809,-0.1275611441216954,-0.037205333702161475,-0.037205333702161475,-0.037205333702161475,-0.037205333702161475,-0.047835429045636084,-0.047835429045636084,-0.047835429045636084,-0.047835429045636084,-0.05846552438911087,-0.05846552438911087,-0.05846552438911087,-0.05846552438911087,-0.037205333702161725,-0.037205333702161725,-0.037205333702161725,-0.037205333702161725,-0.047835429045635834,-0.047835429045635834,-0.047835429045635834,-0.047835429045635834,-0.05846552438911058,-0.05846552438911058,-0.05846552438911058,-0.05846552438911058,-0.03879154890291829,-0.03879154890291829,-0.03879154890291829,-0.04120270848015563,-0.04120270848015563,-0.04120270848015563,-0.03393028948486933,-0.03393028948486933,-0.03393028948486933,-0.03151955746491709,-0.03151955746491709,-0.03151955746491709,-0.02424752187358276,-0.02424752187358276,-0.02424752187358276,-0.026657914642918758,-0.026657914642918758,-0.026657914642918758,-0.04120270848015456,-0.04120270848015456,-0.04120270848015456,-0.03879154890291757,-0.03879154890291757,-0.03879154890291757,-0.031519557464916595,-0.031519557464916595,-0.031519557464916595,-0.03393028948487046,-0.03393028948487046,-0.03393028948487046,-0.0266579146429191,-0.0266579146429191,-0.0266579146429191,-0.024247521873582645,-0.024247521873582645,-0.024247521873582645,-0.01851718920904466,-0.01851718920904466,-0.01851718920904466,-0.01851718920904466,-0.029627502734471456,-0.029627502734471456,-0.029627502734471456,-0.029627502734471456,-0.04740400437515433,-0.015150427534672922,-0.015150427534672922,-0.015150427534672922,-0.015150427534672922,-0.024240684055476674,-0.024240684055476674,-0.024240684055476674,-0.024240684055476674,-0.038785094488762675,-0.011783665860301345,-0.011783665860301345,-0.011783665860301345,-0.011783665860301345,-0.018853865376482152,-0.018853865376482152,-0.018853865376482152,-0.018853865376482152,-0.030166184602371443,-0.018517189209044892,-0.018517189209044892,-0.018517189209044892,-0.018517189209044892,-0.029627502734471827,-0.029627502734471827,-0.029627502734471827,-0.029627502734471827,-0.04740400437515492,-0.015150427534672776,-0.015150427534672776,-0.015150427534672776,-0.015150427534672776,-0.02424068405547644,-0.02424068405547644,-0.02424068405547644,-0.02424068405547644,-0.03878509448876231,-0.011783665860301277,-0.011783665860301277,-0.011783665860301277,-0.011783665860301277,-0.01885386537648204,-0.01885386537648204,-0.01885386537648204,-0.01885386537648204,-0.030166184602371266]),1e-14))
        f=MEDCouplingFieldDouble(ft)
        arr=DataArrayDouble(126,2)
        arr[:, 0] = list(range(126))
        arr[:, 1] = list(range(126))
        arr[:,1]+=1000
        f.setArray(arr)
        f.checkConsistencyLight()
        self.assertTrue(DataArrayDouble(f.integral(False)).isEqual(DataArrayDouble([-211.66121638700983,-4863.9563007698835]),1e-11))
        self.assertTrue(DataArrayDouble(f.getWeightedAverageValue()).isEqual(DataArrayDouble([45.4960858131136,1045.496085813114]),1e-11))
        self.assertTrue(DataArrayDouble(f.normL1()).isEqual(DataArrayDouble([45.49608581311362,1045.496085813114]),1e-11))
        self.assertTrue(DataArrayDouble(f.normL2()).isEqual(DataArrayDouble([58.16846378340898,1046.1241521947334]),1e-11))
        pass

    def testSwig2FieldDiscretizationComputeMeshRestrictionFromTupleIds1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_GAUSS_NE)
        f.setMesh(m)
        a,b=f.getDiscretization().computeMeshRestrictionFromTupleIds(f.getMesh(),[3,4,5,6,8,9,10,14,15,16,17])
        self.assertTrue(a.isEqual(DataArrayInt([1,4])))
        self.assertTrue(b.isEqual(DataArrayInt([4,5,6,14,15,16,17])))
        a,b=f.getDiscretization().computeMeshRestrictionFromTupleIds(f.getMesh(),DataArrayInt([0,1,2,3,5,7,8,9,10,11,12,18]))
        self.assertTrue(a.isEqual(DataArrayInt([0,2])))
        self.assertTrue(b.isEqual(DataArrayInt([0,1,2,3,7,8,9])))
        #
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        a,b=f.getDiscretization().computeMeshRestrictionFromTupleIds(f.getMesh(),[3,4])
        self.assertTrue(a.isEqual(DataArrayInt([3,4])))
        self.assertTrue(b.isEqual(DataArrayInt([3,4])))
        #
        f=MEDCouplingFieldDouble(ON_NODES)
        f.setMesh(m)
        a,b=f.getDiscretization().computeMeshRestrictionFromTupleIds(f.getMesh(),[1,2,3,4])
        self.assertTrue(a.isEqual(DataArrayInt([1])))
        self.assertTrue(b.isEqual(DataArrayInt([1,2,4])))
        #
        f=MEDCouplingDataForTest.buildFieldOnGauss_1()
        a,b=f.getDiscretization().computeMeshRestrictionFromTupleIds(f.getMesh(),[0,11,12,13,14,15,17,18,19,36,37,38,115,117,118,119,120,121,122,123,124,125])
        self.assertTrue(a.isEqual(DataArrayInt([0,11,12,18,35])))
        self.assertTrue(b.isEqual(DataArrayInt([0,11,12,13,14,15,36,37,38,117,118,119,120,121,122,123,124,125])))
        #
        d=DataArrayInt([0,3,7,9,15,18])
        e=DataArrayInt([0,1,2,3,7,8,15,16,17])
        a,b=d.findIdsRangesInListOfIds(e)
        self.assertTrue(a.isEqual(DataArrayInt([0,2,4])))
        self.assertTrue(b.isEqual(DataArrayInt([0,1,2,7,8,15,16,17])))
        pass

    @unittest.skipUnless(checkFreeMemory((223456789*16)/(1024)), "Not enough memory")
    def testSwig2BigMem(self):
        if MEDCouplingSizeOfVoidStar()==64:
            d=DataArrayAsciiChar(223456789,16)
            self.assertTrue(d.getNumberOfTuples(),223456789)
            self.assertTrue(d.getNumberOfComponents(),16)
            d.setIJ(223456788,5,"r")
            self.assertTrue(d.getIJ(223456788,5),'r')
            d[223456787]="1234567890123456"
            self.assertTrue(d[223456787],'1234567890123456')
            self.assertRaises(InterpKernelException,d.rearrange,1)# fails because it would lead to nb of tuples > 2147483647
            pass
        pass

    def testSwig2DAReverseMultiCompo1(self):
        d=DataArrayDouble(6,2)
        d[:, 0] = list(range(6))
        d[:, 1] = list(range(10, 16))
        d.reverse()
        self.assertTrue(d.isEqual(DataArrayDouble([5.,15.,4.,14.,3.,13.,2.,12.,1.,11.,0.,10.],6,2),1e-14))
        d=DataArrayDouble(7,2)
        d[:, 0] = list(range(7))
        d[:, 1] = list(range(10, 17))
        d.reverse()
        self.assertTrue(d.isEqual(DataArrayDouble([6.,16.,5.,15.,4.,14.,3.,13.,2.,12.,1.,11.,0.,10.],7,2),1e-14))
        #
        d=DataArrayInt(6,2)
        d[:, 0] = list(range(6))
        d[:, 1] = list(range(10, 16))
        d.reverse()
        self.assertTrue(d.isEqual(DataArrayInt([5,15,4,14,3,13,2,12,1,11,0,10],6,2)))
        d=DataArrayInt(7,2)
        d[:, 0] = list(range(7))
        d[:, 1] = list(range(10, 17))
        d.reverse()
        self.assertTrue(d.isEqual(DataArrayInt([6,16,5,15,4,14,3,13,2,12,1,11,0,10],7,2)))
        pass

    def testSwigDAPow1(self):
        d=DataArrayInt(10)
        d.iota(0)
        d1=d.deepCopy()
        d.setIJ(2,0,-2)
        self.assertTrue((d**2).isEqual(DataArrayInt([0,1,4,9,16,25,36,49,64,81])))
        self.assertTrue((d**3).isEqual(DataArrayInt([0,1,-8,27,64,125,216,343,512,729])))
        for elt in [d]:
            elt**=2
            pass
        self.assertTrue(d.isEqual(DataArrayInt([0,1,4,9,16,25,36,49,64,81])))
        self.assertTrue((d1[:4]**d1[:4]).isEqual(DataArrayInt([1,1,4,27])))
        self.assertTrue((3**d1[:4]).isEqual(DataArrayInt([1,3,9,27])))
        d2=d1[:4]
        d2**=d2
        self.assertTrue(d2.isEqual(DataArrayInt([1,1,4,27])))
        self.assertRaises(InterpKernelException,d2.__pow__,-1)#non supporting negative pow in DataArrayInt.__pow__
        self.assertRaises(InterpKernelException,d2.__ipow__,-1)#non supporting negative pow in DataArrayInt.__pow__
        #
        d=DataArrayDouble(10)
        d.iota(0)
        d1=d.deepCopy()
        d.setIJ(2,0,-2.)
        self.assertTrue((d**2).isEqual(DataArrayDouble([0,1,4,9,16,25,36,49,64,81]),1e-12))
        self.assertTrue((d**3).isEqual(DataArrayDouble([0,1,-8,27,64,125,216,343,512,729]),1e-12))
        self.assertRaises(InterpKernelException,d.__pow__,3.1)#3.1 is double not integer -> not supporting negative values in d
        for elt in [d]:
            elt**=2
            pass
        self.assertTrue(d.isEqual(DataArrayDouble([0,1,4,9,16,25,36,49,64,81]),1e-12))
        self.assertTrue((d1[:4]**d1[:4]).isEqual(DataArrayDouble([1,1,4,27]),1e-12))
        self.assertTrue((3**d1[:4]).isEqual(DataArrayDouble([1,3,9,27]),1e-12))
        d2=d1[:4]
        d2**=d2
        self.assertTrue(d2.isEqual(DataArrayDouble([1,1,4,27]),1e-12))
        d2**=-0.5
        self.assertTrue(d2.isEqual(DataArrayDouble([1,1,1./2,1./sqrt(27.)]),1e-14))
        d3=-1./d1[1:5]
        self.assertTrue((3**d3).isEqual(DataArrayDouble([0.3333333333333333,0.5773502691896257,0.6933612743506348,0.7598356856515925]),1e-14))
        d4=d3.deepCopy() ; d4.abs()
        self.assertTrue((d4**d3).isEqual(DataArrayDouble([1.,sqrt(2.),1.4422495703074083,sqrt(2.)]),1e-14))
        d4**=d3
        self.assertTrue(d4.isEqual(DataArrayDouble([1.,sqrt(2.),1.4422495703074083,sqrt(2.)]),1e-14))
        pass
    
    def testSwig2Baryenter3DForCellsWithVolumeZero1(self):
        coo=DataArrayDouble([0.,0.,0.,1.,0.,0.,0.,1.,0.],3,3)
        m2=MEDCouplingUMesh("mesh",2)
        m2.allocateCells(0)
        m2.insertNextCell(NORM_POLYGON,[0,1,2])
        m2.setCoords(coo)
        m2.checkConsistency()
        #
        coo2=DataArrayDouble([0.,0.,0.,0.,0.,0.,0.,0.,2.],3,3)
        m1=MEDCouplingUMesh("mesh",1)
        m1.allocateCells(0)
        m1.insertNextCell(NORM_SEG2,[0,1])
        m1.insertNextCell(NORM_SEG2,[1,2])
        m1.setCoords(coo2)
        m1.checkConsistency()
        #
        m3=m2.buildExtrudedMesh(m1,0)
        m3.insertNextCell(NORM_POLYHED,[3,4,5,-1,8,7,6,-1,4,3,6,7,-1,5,4,7,8,-1,5,4,-1,3,5,8,6])# addition of face #4 with null surface
        self.assertTrue(m3.computeCellCenterOfMass().isEqual(DataArrayDouble([0.3333333333333333,0.3333333333333333,0.,0.3333333333333333,0.3333333333333333,1.,0.3333333333333333,0.3333333333333333,1.],3,3),1e-13))
        m4,a,b,c,d=m3.buildDescendingConnectivity()
        self.assertTrue(m4.computeCellCenterOfMass().isEqual(DataArrayDouble([0.3333333333333333,0.3333333333333333,0.,0.3333333333333333,0.3333333333333333,0.,0.5,0.,0.,0.5,0.5,0.,0.,0.5,0.,0.3333333333333333,0.3333333333333333,2.,0.5,0.,1.,0.5,0.5,1.,0.,0.5,1.,0.5,0.5,0.],10,3),1e-13))
        pass

    def testSwigRepr1(self):
        d=DataArrayDouble()
        self.assertTrue(len(d.__repr__())<120)
        d.alloc(1000,0) ; self.assertTrue(len(d.__repr__())<100)
        for i in range(100):
            d.alloc(i,1) ; d.iota(1.1234567890123456) ; d*=1e123
            self.assertTrue(len(d.__repr__())<500)
            pass
        for i in range(50):
            d.alloc(i,2) ; d.rearrange(1) ; d.iota(1.1234567890123456) ; d.rearrange(2) ; d*=1e123
            self.assertTrue(len(d.__repr__())<500)
            pass
        d.alloc(4000,1) ; d.iota() ; self.assertTrue(len(d.__repr__())<500)
        for i in range(2, 4):
            d.alloc(362880,1) ; d.iota() ; d.rearrange(i) ; self.assertTrue(len(d.__repr__())<500)
            pass
        d.alloc(0,9)
        self.assertTrue(len(d.__repr__())<120)
        #
        d=DataArrayInt()
        self.assertTrue(len(d.__repr__())<100)
        d.alloc(1000,0) ; self.assertTrue(len(d.__repr__())<100)
        for i in range(100):
            d.alloc(i,1) ; d.iota(123456789)
            self.assertTrue(len(d.__repr__())<500)
            pass
        for i in range(50):
            d.alloc(i,2) ; d.rearrange(1) ; d.iota(123456789) ; d.rearrange(2)
            self.assertTrue(len(d.__repr__())<500)
            pass
        d.alloc(4000,1) ; d.iota() ; self.assertTrue(len(d.__repr__())<500)
        for i in range(2, 10):
            d.alloc(362880,1) ; d.iota() ; d.rearrange(i) ; self.assertTrue(len(d.__repr__())<500)
            pass
        d.alloc(0,9)
        self.assertTrue(len(d.__repr__())<100)
        #
        d=DataArrayAsciiChar()
        d.alloc(1000,0) ; self.assertTrue(len(d.__repr__())<100)
        d.alloc(2,16) ; d[:]='1234567890ABCDEF'
        self.assertTrue(len(d.__repr__())<500)
        d.alloc(2000,16) ; d[:]='1234567890ABCDEF'
        self.assertTrue(len(d.__repr__())<500)
        d.alloc(0,16) ; d[:]='1234567890ABCDEF'
        self.assertTrue(len(d.__repr__())<120)
        #
        d=DataArrayByte()
        self.assertTrue(len(d.__repr__())<100)
        d.alloc(1000,0) ; self.assertTrue(len(d.__repr__())<100)
        d.alloc(0,16) ; self.assertTrue(len(d.__repr__())<100)
        d.alloc(5,1) ; d.fillWithValue(127)
        self.assertTrue(len(d.__repr__())<200)
        d.alloc(1000,1) ; d.fillWithValue(127)
        self.assertTrue(len(d.__repr__())<500)
        d.alloc(1000,3) ; d.fillWithValue(127)
        self.assertTrue(len(d.__repr__())<500)
        pass
    
    def testSwig2MeshComputeIsoBarycenterOfNodesPerCell1(self):
        coo=DataArrayDouble([26.17509821414239,5.0374,200.,26.175098214142388,-5.0374,200.,17.450065476094927,20.1496,200.,8.725032738047464,25.187,200.,43.62516369023732,5.0374,200.,34.90013095218986,10.0748,200.,34.900130952189855,-10.0748,200.,43.625163690237315,-5.0374,200.,26.175098214142402,25.187,200.,26.175098214142395,35.2618,200.,17.45006547609493,40.2992,200.,8.725032738047469,35.2618,200.,26.17509821414239,5.0374,200.,26.175098214142388,-5.0374,200.,17.450065476094927,20.1496,200.,8.725032738047464,25.187,200.,43.62516369023732,5.0374,200.,34.90013095218986,10.0748,200.,34.900130952189855,-10.0748,200.,43.625163690237315,-5.0374,200.,26.175098214142402,25.187,200.,26.175098214142395,35.2618,200.,17.45006547609493,40.2992,200.,8.725032738047469,35.2618,200.],24,3)
        m=MEDCouplingUMesh.New("toto",3)
        m.allocateCells(0)
        m.insertNextCell(NORM_POLYHED,[4,5,0,1,6,7,-1,19,18,13,12,17,16,-1,5,4,16,17,-1,0,5,17,12,-1,1,0,12,13,-1,6,1,13,18,-1,7,6,18,19,-1,4,7,19,16])
        m.insertNextCell(NORM_POLYHED,[9,10,11,3,2,8,-1,20,14,15,23,22,21,-1,10,9,21,22,-1,11,10,22,23,-1,3,11,23,15,-1,2,3,15,14,-1,8,2,14,20,-1,9,8,20,21])
        m.setCoords(coo)
        m.checkConsistency()
        #
        dReference=DataArrayDouble([(34.900130952189848,0.,200),(17.450065476094931,30.2244,200.)])
        self.assertTrue(m.computeIsoBarycenterOfNodesPerCell().isEqual(dReference,1e-12))
        m.getNodalConnectivity().setIJ(87,0,24)
        self.assertRaises(InterpKernelException,m.computeIsoBarycenterOfNodesPerCell)
        m.getNodalConnectivity().setIJ(87,0,-2)
        self.assertRaises(InterpKernelException,m.computeIsoBarycenterOfNodesPerCell)
        m.getNodalConnectivity().setIJ(87,0,21)# put again 21 as at the beginning
        #
        self.assertTrue(m.unPolyze())
        self.assertEqual([NORM_HEXGP12],m.getAllGeoTypes())
        self.assertTrue(m.computeIsoBarycenterOfNodesPerCell().isEqual(dReference,1e-12))
        m.getNodalConnectivity().setIJ(25,0,24)
        self.assertRaises(InterpKernelException,m.computeIsoBarycenterOfNodesPerCell)
        m.getNodalConnectivity().setIJ(25,0,-1)
        self.assertRaises(InterpKernelException,m.computeIsoBarycenterOfNodesPerCell)
        pass

    def testSwig2NonRegressionBugDescHexa20(self):
        coo=DataArrayDouble([0.,0.,0.,1.23,0.,0.,0.615,0.,0.,0.,2.1,0.,0.615,2.1,0.,1.23,2.1,0.,1.23,1.05,0.,0.,1.05,0.,0.,0.,2.16,1.23,0.,2.16,1.23,2.1,2.16,0.,2.1,2.16,0.,0.,4.32,0.615,0.,4.32,1.23,0.,4.32,1.23,1.05,4.32,1.23,2.1,4.32,0.615,2.1,4.32,0.,2.1,4.32,0.,1.05,4.32],20,3)
        m=MEDCouplingUMesh('mesh',3)
        m.allocateCells(0)
        m.insertNextCell(NORM_HEXA20,[0,3,5,1,12,18,16,14,7,4,6,2,19,17,15,13,8,11,10,9])
        m.setCoords(coo)
        m.checkConsistency()
        #
        a,b,c,d,e=m.buildDescendingConnectivity()
        m2=MEDCouplingUMesh('mesh',2)
        m2.allocateCells(0)
        m2.setCoords(coo)
        conn2=[[0,3,5,1,7,4,6,2],[12,14,16,18,13,15,17,19],[0,12,18,3,8,19,11,7],[3,18,16,5,11,17,10,4],[5,16,14,1,10,15,9,6],[1,14,12,0,9,13,8,2]]
        for i in range(6):
            m2.insertNextCell(NORM_QUAD8,conn2[i])
            pass
        self.assertTrue(m2.isEqual(a,1e-12))
        self.assertTrue(b.isEqual(DataArrayInt([0,1,2,3,4,5])))
        self.assertTrue(c.isEqual(DataArrayInt([0,6])))
        self.assertTrue(d.isEqual(DataArrayInt([0,0,0,0,0,0])))
        self.assertTrue(e.isEqual(DataArrayInt([0,1,2,3,4,5,6])))
        #
        m.convertQuadraticCellsToLinear() ; m.zipCoords()
        m.convertLinearCellsToQuadratic(1)
        #
        coo2=DataArrayDouble([0.,0.,0.,1.23,0.,0.,0.,2.1,0.,1.23,2.1,0.,0.,0.,4.32,1.23,0.,4.32,1.23,2.1,4.32,0.,2.1,4.32,0.,1.05,0.,0.615,2.1,0.,1.23,1.05,0.,0.615,0.,0.,0.,1.05,4.32,0.615,2.1,4.32,1.23,1.05,4.32,0.615,0.,4.32,0.,0.,2.16,0.,2.1,2.16,1.23,2.1,2.16,1.23,0.,2.16,0.615,1.05,0.,0.,1.05,2.16,0.615,2.1,2.16,1.23,1.05,2.16,0.615,0.,2.16,0.615,1.05,4.32,0.615,1.05,2.16],27,3)
        m3=MEDCouplingUMesh("mesh",3)
        m3.allocateCells(1)
        m3.insertNextCell(NORM_HEXA27,[0,2,3,1,4,7,6,5,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26])
        m3.setCoords(coo2)
        self.assertTrue(m3.isEqual(m,1e-12))
        #
        a,b,c,d,e=m.buildDescendingConnectivity()
        conn4=[[0,2,3,1,8,9,10,11,20],[4,5,6,7,15,14,13,12,25],[0,4,7,2,16,12,17,8,21],[2,7,6,3,17,13,18,9,22],[3,6,5,1,18,14,19,10,23],[1,5,4,0,19,15,16,11,24]]
        m4=MEDCouplingUMesh("mesh",2)
        m4.allocateCells(0)
        for i in range(6):
            m4.insertNextCell(NORM_QUAD9,conn4[i])
            pass
        m4.setCoords(coo2)
        self.assertTrue(m4.isEqual(a,1e-12))
        self.assertTrue(b.isEqual(DataArrayInt([0,1,2,3,4,5])))
        self.assertTrue(c.isEqual(DataArrayInt([0,6])))
        self.assertTrue(d.isEqual(DataArrayInt([0,0,0,0,0,0])))
        self.assertTrue(e.isEqual(DataArrayInt([0,1,2,3,4,5,6])))
        pass
    
    def testSwigAdvGauss(self):
        f=MEDCouplingFieldTemplate(ON_GAUSS_PT)
        f.setDiscretization(None)
        f.__repr__() ; f.__str__()
        #
        f=MEDCouplingFieldTemplate(ON_GAUSS_PT)
        d=f.getDiscretization()
        i=DataArrayInt() ; i.alloc(10,1) ; i.iota(1)
        d.setArrayOfDiscIds(i)
        f.__repr__() ; f.__str__()
        i2=d.getArrayOfDiscIds()
        self.assertEqual(i.__repr__(),i2.__repr__())
        #
        f=MEDCouplingFieldDouble(ON_GAUSS_PT)
        f.setDiscretization(None)
        f.__repr__() ; f.__str__()
        #
        f=MEDCouplingFieldDouble(ON_GAUSS_PT)
        d=f.getDiscretization()
        i=DataArrayInt() ; i.alloc(10,1) ; i.iota(1)
        d.setArrayOfDiscIds(i)
        f.__repr__() ; f.__str__()
        #
        gl=MEDCouplingGaussLocalization(NORM_SEG2,[0,1],[0.5],[1.])
        gl.setWeights([3.])
        gl.__repr__() ; gl.__str__()
        gl=MEDCouplingGaussLocalization(NORM_ERROR)
        gl.setWeights([3.])
        gl.__repr__() ; gl.__str__()
        pass

    def testSwig2NonRegressionBugSubstractInPlaceDM(self):
        m0=MEDCouplingCMesh()
        arr=DataArrayDouble(5,1) ; arr.iota(0.)
        m0.setCoords(arr,arr)
        m0=m0.buildUnstructured()
        m00=m0[::2] ; m00.simplexize(0) ; m01=m0[1::2]
        m0=MEDCouplingUMesh.MergeUMeshes([m00,m01])
        m0.getCoords()[:]*=1/4.
        m0.setName("mesh")
        #
        NodeField=MEDCouplingFieldDouble(ON_NODES,ONE_TIME) ; NodeField.setTime(5.6,5,6) ; NodeField.setMesh(m0)
        NodeField.setName("NodeField")
        NodeField.fillFromAnalytic(1,"exp(-((x-1)*(x-1)+(y-1)*(y-1)))") ; NodeField.getArray().setInfoOnComponent(0,"powernode [W]")
        proc0=m0.getCellsInBoundingBox([(0.,0.4),(0.,0.4)],1e-10)
        proc1=proc0.buildComplement(m0.getNumberOfCells())
        #
        NodeField0=NodeField[proc0] ; NodeField0.getMesh().setName(m0.getName())
        NodeField1=NodeField[proc1] ; NodeField1.getMesh().setName(m0.getName())
        #
        NodeField_read=MEDCouplingFieldDouble.MergeFields([NodeField0,NodeField1])
        NodeField_read.mergeNodes(1e-10)
        NodeFieldCpy=NodeField.deepCopy()
        NodeFieldCpy.mergeNodes(1e-10)
        NodeField.checkConsistencyLight()
        self.assertTrue(not NodeField.getArray().isUniform(0.,1e-12))
        NodeField.substractInPlaceDM(NodeField_read,10,1e-12)
        self.assertTrue(NodeField.getArray().isUniform(0.,1e-12))
        pass

    def testSwigFieldOperationOpen1(self):
        ## MEDCouplingFieldDouble.__add__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:, 0] = list(range(5)) ; arr[:, 1] = 2 * arr[:, 0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__add__,2)
        self.assertRaises(InterpKernelException, f.__add__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__add__,arr)
        self.assertRaises(InterpKernelException,f.__add__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__add__,2)
        self.assertRaises(InterpKernelException, f.__add__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__add__,arr)
        self.assertRaises(InterpKernelException,f.__add__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:, 0] = list(range(5)) ; f.getArray()[:, 1] = f.getArray()[:, 0] + 7
        ff=f+2
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(2,9),(3,10),(4,11),(5,12),(6,13)]),1e-12))
        ff=f+arr
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,7),(2,10),(4,13),(6,16),(8,19)]),1e-12))
        self.assertRaises(InterpKernelException,f.__add__,f2)
        f2.setArray(arr)
        ff=f+f2
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,7),(2,10),(4,13),(6,16),(8,19)]),1e-12))
        ff=f+[5,8]
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(5,15),(6,16),(7,17),(8,18),(9,19)]),1e-12))
        ### MEDCouplingFieldDouble.__sub__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:, 0] = list(range(5)) ; arr[:, 1] = 2 * arr[:, 0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__sub__,2)
        self.assertRaises(InterpKernelException, f.__sub__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__sub__,arr)
        self.assertRaises(InterpKernelException,f.__sub__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__sub__,2)
        self.assertRaises(InterpKernelException, f.__sub__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__sub__,arr)
        self.assertRaises(InterpKernelException,f.__sub__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:, 0] = list(range(5)) ; f.getArray()[:, 1] = f.getArray()[:, 0] + 7
        ff=f-2
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(-2,5),(-1,6),(0,7),(1,8),(2,9)]),1e-12))
        ff=f-arr
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,7),(0,6),(0,5),(0,4),(0,3)]),1e-12))
        self.assertRaises(InterpKernelException,f.__sub__,f2)
        f2.setArray(arr)
        ff=f-f2
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,7),(0,6),(0,5),(0,4),(0,3)]),1e-12))
        ff=f-[5,8]
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(-5,-1),(-4,0),(-3,1),(-2,2),(-1,3)]),1e-12))
        ### MEDCouplingFieldDouble.__mul__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:, 0] = list(range(5)) ; arr[:, 1] = 2 * arr[:, 0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__mul__,2)
        self.assertRaises(InterpKernelException, f.__mul__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__mul__,arr)
        self.assertRaises(InterpKernelException,f.__mul__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__mul__,2)
        self.assertRaises(InterpKernelException, f.__mul__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__mul__,arr)
        self.assertRaises(InterpKernelException,f.__mul__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:, 0] = list(range(5)) ; f.getArray()[:, 1] = f.getArray()[:, 0] + 7
        ff=f*2
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,14),(2,16),(4,18),(6,20),(8,22)]),1e-12))
        ff=f*arr
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,0),(1,16),(4,36),(9,60),(16,88)]),1e-12))
        self.assertRaises(InterpKernelException,f.__mul__,f2)
        f2.setArray(arr)
        ff=f*f2
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,0),(1,16),(4,36),(9,60),(16,88)]),1e-12))
        ff=f*[5,8]
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,56),(5,64),(10,72),(15,80),(20,88)]),1e-12))
        ### MEDCouplingFieldDouble.__div__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:, 0] = list(range(1, 6)) ; arr[:, 1] = 2 * arr[:, 0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__div__,2)
        self.assertRaises(InterpKernelException, f.__div__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__div__,arr)
        self.assertRaises(InterpKernelException,f.__div__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__div__,2)
        self.assertRaises(InterpKernelException, f.__div__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__div__,arr)
        self.assertRaises(InterpKernelException,f.__div__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:, 0] = list(range(5)) ; f.getArray()[:, 1] = f.getArray()[:, 0] + 7
        self.assertRaises(InterpKernelException,f.__div__,0)
        ff=f/2
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,3.5),(0.5,4),(1,4.5),(1.5,5),(2,5.5)]),1e-12))
        ff=f/arr
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,3.5),(0.5,2),(0.6666666666666666,1.5),(0.75,1.25),(0.8,1.1)]),1e-12))
        self.assertRaises(InterpKernelException,f.__div__,f2)
        f2.setArray(arr)
        ff=f/f2
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,3.5),(0.5,2),(0.6666666666666666,1.5),(0.75,1.25),(0.8,1.1)]),1e-12))
        ff=f/[5,8]
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,0.875),(0.2,1),(0.4,1.125),(0.6,1.25),(0.8,1.375)]),1e-12))
        ### MEDCouplingFieldDouble.__pow__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5)
        arr[:]=[1,1,3,2,0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__div__,2)
        self.assertRaises(InterpKernelException, f.__div__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__div__,arr)
        self.assertRaises(InterpKernelException,f.__div__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__div__,2)
        self.assertRaises(InterpKernelException, f.__div__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__div__,arr)
        self.assertRaises(InterpKernelException,f.__div__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,1)
        f.getArray()[:] = list(range(2, 7))
        ff=f**2
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([4,9,16,25,36]),1e-12))
        ff=f**arr
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([2,3,64,25,1]),1e-12))
        f2.setArray(arr)
        ff=f**f2
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([2,3,64,25,1]),1e-12))
        ## MEDCouplingFieldDouble.__iadd__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:, 0] = list(range(5)) ; arr[:, 1] = 2 * arr[:, 0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__iadd__,2)
        self.assertRaises(InterpKernelException, f.__iadd__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__iadd__,arr)
        self.assertRaises(InterpKernelException,f.__iadd__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__iadd__,2)
        self.assertRaises(InterpKernelException, f.__iadd__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__iadd__,arr)
        self.assertRaises(InterpKernelException,f.__iadd__,f2)
        f.getArray().alloc(5,2)
        f.getArray()[:, 0] = list(range(5)) ; f.getArray()[:, 1] = f.getArray()[:, 0] + 7
        f.checkConsistencyLight()
        f+=2
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(2,9),(3,10),(4,11),(5,12),(6,13)]),1e-12))
        f+=arr
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(2,9),(4,12),(6,15),(8,18),(10,21)]),1e-12))
        f2.setArray(arr)
        f+=f2
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(2,9),(5,14),(8,19),(11,24),(14,29)]),1e-12))
        f+=[0.1,0.2]
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(2.1,9.2),(5.1,14.2),(8.1,19.2),(11.1,24.2),(14.1,29.2)]),1e-12))
        ## MEDCouplingFieldDouble.__isub__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:, 0] = list(range(5)) ; arr[:, 1] = 2 * arr[:, 0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__isub__,2)
        self.assertRaises(InterpKernelException, f.__isub__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__isub__,arr)
        self.assertRaises(InterpKernelException,f.__isub__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__isub__,2)
        self.assertRaises(InterpKernelException, f.__isub__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__isub__,arr)
        self.assertRaises(InterpKernelException,f.__isub__,f2)
        f.getArray().alloc(5,2)
        f.getArray()[:, 0] = list(range(5)) ; f.getArray()[:, 1] = f.getArray()[:, 0] + 7
        f.checkConsistencyLight()
        f-=2
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(-2,5),(-1,6),(0,7),(1,8),(2,9)]),1e-12))
        f-=arr
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(-2,5),(-2,4),(-2,3),(-2,2),(-2,1)]),1e-12))
        f2.setArray(arr)
        f-=f2
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(-2,5),(-3,2),(-4,-1),(-5,-4),(-6,-7)]),1e-12))
        f-=[0.1,0.2]
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(-2.1,4.8),(-3.1,1.8),(-4.1,-1.2),(-5.1,-4.2),(-6.1,-7.2)]),1e-12))
        ## MEDCouplingFieldDouble.__imul__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:, 0] = list(range(5)) ; arr[:, 1] = 2 * arr[:, 0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__imul__,2)
        self.assertRaises(InterpKernelException, f.__imul__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__imul__,arr)
        self.assertRaises(InterpKernelException,f.__imul__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__imul__,2)
        self.assertRaises(InterpKernelException, f.__imul__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__imul__,arr)
        self.assertRaises(InterpKernelException,f.__imul__,f2)
        f.getArray().alloc(5,2)
        f.getArray()[:, 0] = list(range(5)) ; f.getArray()[:, 1] = f.getArray()[:, 0] + 7
        f.checkConsistencyLight()
        f*=2
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,14),(2,16),(4,18),(6,20),(8,22)]),1e-12))
        f*=arr
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,0),(2,32),(8,72),(18,120),(32,176)]),1e-12))
        f2.setArray(arr)
        f*=f2
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,0),(2,64),(16,288),(54,720),(128,1408)]),1e-12))
        f*=[0.1,0.2]
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,0),(0.2,12.8),(1.6,57.6),(5.4,144),(12.8,281.6)]),1e-12))
        ## MEDCouplingFieldDouble.__idiv__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:, 0] = list(range(1, 6)) ; arr[:, 1] = 2 * arr[:, 0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__idiv__,2)
        self.assertRaises(InterpKernelException, f.__idiv__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__idiv__,arr)
        self.assertRaises(InterpKernelException,f.__idiv__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__idiv__,2)
        self.assertRaises(InterpKernelException, f.__idiv__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__idiv__,arr)
        self.assertRaises(InterpKernelException,f.__idiv__,f2)
        f.getArray().alloc(5,2)
        f.getArray()[:, 0] = list(range(5)) ; f.getArray()[:, 1] = f.getArray()[:, 0] + 7
        f.checkConsistencyLight()
        f/=2
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,3.5),(0.5,4),(1,4.5),(1.5,5),(2,5.5)]),1e-12))
        f/=arr
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,1.75),(0.25,1),(0.3333333333333333,0.75),(0.375,0.625),(0.4,0.55)]),1e-12))
        f2.setArray(arr)
        f/=f2
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,0.875),(0.125,0.25),(0.1111111111111111,0.125),(0.09375,0.078125),(0.08,0.055)]),1e-12))
        f/=[0.1,0.2]
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,4.375),(1.25,1.25),(1.1111111111111111,0.625),(0.9375,0.390625),(0.8,0.275)]),1e-12))
        ## MEDCouplingFieldDouble.__ipow__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:, 0] = list(range(1, 6)) ; arr[:, 1] = 2 * arr[:, 0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__ipow__,2)
        self.assertRaises(InterpKernelException, f.__ipow__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__ipow__,arr)
        self.assertRaises(InterpKernelException,f.__ipow__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__ipow__,2)
        self.assertRaises(InterpKernelException, f.__ipow__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__ipow__,arr)
        self.assertRaises(InterpKernelException,f.__ipow__,f2)
        f.getArray().alloc(5,2)
        f.getArray()[:, 0] = list(range(5)) ; f.getArray()[:, 1] = f.getArray()[:, 0] + 7
        f.checkConsistencyLight()
        f**=2
        f.checkConsistencyLight()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,49),(1,64),(4,81),(9,100),(16,121)]),1e-12))
         ## MEDCouplingFieldDouble.__radd__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:, 0] = list(range(5)) ; arr[:, 1] = 2 * arr[:, 0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__radd__,2)
        self.assertRaises(InterpKernelException, f.__radd__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__radd__,arr)
        self.assertRaises(InterpKernelException,f.__radd__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__radd__,2)
        self.assertRaises(InterpKernelException, f.__radd__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__radd__,arr)
        self.assertRaises(InterpKernelException,f.__radd__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:, 0] = list(range(5)) ; f.getArray()[:, 1] = f.getArray()[:, 0] + 7
        ff=2+f
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(2,9),(3,10),(4,11),(5,12),(6,13)]),1e-12))
        ff=arr+f
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,7),(2,10),(4,13),(6,16),(8,19)]),1e-12))
        self.assertRaises(InterpKernelException,f.__radd__,f2)
        ff=[5,8]+f
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(5,15),(6,16),(7,17),(8,18),(9,19)]),1e-12))
        ### MEDCouplingFieldDouble.__rsub__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:, 0] = list(range(5)) ; arr[:, 1] = 2 * arr[:, 0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__rsub__,2)
        self.assertRaises(InterpKernelException, f.__rsub__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__rsub__,arr)
        self.assertRaises(InterpKernelException,f.__rsub__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__rsub__,2)
        self.assertRaises(InterpKernelException, f.__rsub__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__rsub__,arr)
        self.assertRaises(InterpKernelException,f.__rsub__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:, 0] = list(range(5)) ; f.getArray()[:, 1] = f.getArray()[:, 0] + 7
        ff=2-f
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(2,-5),(1,-6),(0,-7),(-1,-8),(-2,-9)]),1e-12))
        ff=arr-f
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,-7),(0,-6),(0,-5),(0,-4),(0,-3)]),1e-12))
        self.assertRaises(InterpKernelException,f.__rsub__,f2)
        ### MEDCouplingFieldDouble.__rmul__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:, 0] = list(range(5)) ; arr[:, 1] = 2 * arr[:, 0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__rmul__,2)
        self.assertRaises(InterpKernelException, f.__rmul__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__rmul__,arr)
        self.assertRaises(InterpKernelException,f.__rmul__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__rmul__,2)
        self.assertRaises(InterpKernelException, f.__rmul__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__rmul__,arr)
        self.assertRaises(InterpKernelException,f.__rmul__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:, 0] = list(range(5)) ; f.getArray()[:, 1] = f.getArray()[:, 0] + 7
        ff=2*f
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,14),(2,16),(4,18),(6,20),(8,22)]),1e-12))
        ff=arr*f
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,0),(1,16),(4,36),(9,60),(16,88)]),1e-12))
        self.assertRaises(InterpKernelException,f.__rmul__,f2)
        ff=f*[5,8]
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,56),(5,64),(10,72),(15,80),(20,88)]),1e-12))
        ### MEDCouplingFieldDouble.__rdiv__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:, 0] = list(range(1, 6)) ; arr[:, 1] = 2 * arr[:, 0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__rdiv__,2)
        self.assertRaises(InterpKernelException, f.__rdiv__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__rdiv__,arr)
        self.assertRaises(InterpKernelException,f.__rdiv__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__rdiv__,2)
        self.assertRaises(InterpKernelException, f.__rdiv__, list(range(5)))
        self.assertRaises(InterpKernelException,f.__rdiv__,arr)
        self.assertRaises(InterpKernelException,f.__rdiv__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:, 0] = list(range(1, 6)) ; f.getArray()[:, 1] = f.getArray()[:, 0] + 7
        ff=2/f
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(2,0.25),(1,0.22222222222222221),(0.66666666666666663,0.20000000000000001),(0.5,0.18181818181818182),(0.40000000000000002,0.16666666666666666)]),1e-12))
        ff=arr/f
        ff.checkConsistencyLight()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(1,0.25),(1,0.44444444444444442),(1,0.59999999999999998),(1,0.72727272727272729),(1,0.83333333333333337)]),1e-12))
        self.assertRaises(InterpKernelException,f.__rdiv__,f2)
        pass

    pass

if __name__ == '__main__':
    unittest.main()

