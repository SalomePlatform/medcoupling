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

class MEDCouplingBasicsTest3(unittest.TestCase):
    def testSwigGetItem1(self):
        da=DataArrayInt.New()
        da.alloc(16,3)
        da.rearrange(1)
        da.iota(7)
        da.rearrange(3)
        da.setInfoOnComponent(0,"X [m]")
        da.setInfoOnComponent(1,"Y [m]")
        da.setInfoOnComponent(2,"Z [km]")
        da2=da[5:-1]
        self.assertEqual([22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51],da2.getValues())
        da2=da[4]
        self.assertEqual([19, 20, 21],da2.getValues())
        try:
            da2=da[4:17]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        da2=da[5:-2,2]
        self.assertEqual([24, 27, 30, 33, 36, 39, 42, 45, 48],da2.getValues())
        da2=da[5:8,:]
        self.assertEqual([22, 23, 24, 25, 26, 27, 28, 29, 30],da2.getValues())
        da2=da[:]
        self.assertTrue(da2.isEqual(da))
        da2=da[:,:]
        self.assertTrue(da2.isEqual(da))
        try:
            da2=da[:,:,:]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        self.assertTrue(da[5:8,-2].isEqualWithoutConsideringStr(DataArrayInt([23,26,29])))
        da2=da[5:8,:-2]
        self.assertEqual([22, 25, 28],da2.getValues())
        try:
            da2=da[5:-18,2]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        da2=da[5:5,2]
        self.assertEqual([],da2.getValues())
        pass

    def testSwigGetItem2(self):
        da=DataArrayDouble.New()
        da.alloc(16,3)
        da.rearrange(1)
        da.iota(7)
        da.rearrange(3)
        da.setInfoOnComponent(0,"X [m]")
        da.setInfoOnComponent(1,"Y [m]")
        da.setInfoOnComponent(2,"Z [km]")
        da2=da[5:-1]
        self.assertEqual([22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 51.],da2.getValues())
        da2=da[4]
        self.assertEqual([19., 20., 21],da2.getValues())
        try:
            da2=da[4:17]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        da2=da[5:-2,2]
        self.assertEqual([24., 27., 30., 33., 36., 39., 42., 45., 48.],da2.getValues())
        da2=da[5:8,:]
        self.assertEqual([22., 23., 24., 25., 26., 27., 28., 29., 30.],da2.getValues())
        da2=da[:]
        self.assertTrue(da2.isEqual(da,1e-12))
        da2=da[:,:]
        self.assertTrue(da2.isEqual(da,1e-12))
        try:
            da2=da[:,:,:]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        self.assertTrue(da[5:8,-2].isEqualWithoutConsideringStr(DataArrayDouble([23.,26.,29.]),1e-12))
        da2=da[5:8,:-2]
        self.assertEqual([22., 25., 28.],da2.getValues())
        try:
            da2=da[5:-18,2]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        da2=da[5:5,2]
        self.assertEqual([],da2.getValues())
        pass

    def testSwigSetItem1(self):
        da=DataArrayInt.New()
        da.alloc(20,1)
        da.iota(7)
        da.rearrange(5)
        da.setInfoOnComponent(0,"X [m]") ; da.setInfoOnComponent(1,"Y [km]") ; da.setInfoOnComponent(2,"Y [m]")
        da.setInfoOnComponent(3,"Z [W]") ; da.setInfoOnComponent(4,"ZZ [km]") ; 
        da[:,2]=3
        self.assertEqual([7, 8, 3, 10, 11, 12, 13, 3, 15, 16, 17, 18, 3, 20, 21, 22, 23, 3, 25, 26],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[2]=3
        self.assertEqual([7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 3, 3, 3, 3, 3, 22, 23, 24, 25, 26],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[[0,3]]=-1
        self.assertEqual([-1, -1, -1, -1, -1, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, -1, -1, -1, -1, -1],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[:,[1,3,4]]=-3
        self.assertEqual([7, -3, 9, -3, -3, 12, -3, 14, -3, -3, 17, -3, 19, -3, -3, 22, -3, 24, -3, -3],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da2=DataArrayInt.New() ; da2.setValues([0,2,3],3,1)
        da[da2]=-7
        self.assertEqual([-7, -7, -7, -7, -7, 12, 13, 14, 15, 16, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,-2:]=-7
        self.assertEqual([7, 8, 9, -7, -7, 12, 13, 14, 15, 16, 17, 18, 19, -7, -7, 22, 23, 24, -7, -7],da.getValues())
        # Let's test with DAI right hand side
        da1=DataArrayInt.New()
        da1.setValues([25,26,27,125,126,127],2,3)
        #
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[-2:,1:4]=da1
        self.assertEqual([7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 25, 26, 27, 21, 22, 125, 126, 127, 26],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[1:,3]=[225,226,227]
        self.assertEqual([7, 8, 9, 10, 11, 12, 13, 14, 225, 16, 17, 18, 19, 226, 21, 22, 23, 24, 227, 26],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[1,2:]=[225,226,227]
        self.assertEqual([7, 8, 9, 10, 11, 12, 13, 225, 226, 227, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,-2:]=[88,99,1010,1111,1212,1313]
        self.assertEqual([7, 8, 9, 88, 99, 12, 13, 14, 15, 16, 17, 18, 19, 1010, 1111, 22, 23, 24, 1212, 1313],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da3=DataArrayInt.New(); da3.setValues([88,99,1010,1111,1212,1313],3,2)
        da[da2,-2:]=da3
        self.assertEqual([7, 8, 9, 88, 99, 12, 13, 14, 15, 16, 17, 18, 19, 1010, 1111, 22, 23, 24, 1212, 1313],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,[0,2]]=da3
        self.assertEqual([88, 8, 99, 10, 11, 12, 13, 14, 15, 16, 1010, 18, 1111, 20, 21, 1212, 23, 1313, 25, 26],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,0:3:2]=da3
        self.assertEqual([88, 8, 99, 10, 11, 12, 13, 14, 15, 16, 1010, 18, 1111, 20, 21, 1212, 23, 1313, 25, 26],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,0:3:2]=-8
        self.assertEqual([-8, 8, -8, 10, 11, 12, 13, 14, 15, 16, -8, 18, -8, 20, 21, -8, 23, -8, 25, 26],da.getValues())
        pass

    def testSwigSetItem2(self):
        da=DataArrayDouble.New()
        da.alloc(20,1)
        da.iota(7)
        da.rearrange(5)
        da.setInfoOnComponent(0,"X [m]") ; da.setInfoOnComponent(1,"Y [km]") ; da.setInfoOnComponent(2,"Y [m]")
        da.setInfoOnComponent(3,"Z [W]") ; da.setInfoOnComponent(4,"ZZ [km]") ; 
        da[:,2]=3.
        self.assertEqual([7., 8., 3., 10., 11., 12., 13., 3., 15., 16., 17., 18., 3., 20., 21., 22., 23., 3., 25., 26.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[2]=3.
        self.assertEqual([7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 3., 3., 3., 3., 3., 22., 23., 24., 25., 26.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[[0,3]]=-1.
        self.assertEqual([-1., -1., -1., -1., -1., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., -1., -1., -1., -1., -1.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[:,[1,3,4]]=-3.
        self.assertEqual([7., -3., 9., -3., -3., 12., -3., 14., -3., -3., 17., -3., 19., -3., -3., 22., -3., 24., -3., -3.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da2=DataArrayInt.New() ; da2.setValues([0,2,3],3,1)
        da[da2]=-7.
        self.assertEqual([-7., -7., -7., -7., -7., 12., 13., 14., 15., 16., -7., -7., -7., -7., -7., -7., -7., -7., -7., -7.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,-2:]=-7
        self.assertEqual([7., 8., 9., -7., -7., 12., 13., 14., 15., 16., 17., 18., 19., -7., -7., 22., 23., 24., -7., -7.],da.getValues())
        # Let's test with DAI right hand side
        da1=DataArrayDouble.New()
        da1.setValues([25,26,27,125,126,127],2,3)
        #
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[-2:,1:4]=da1
        self.assertEqual([7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 25., 26., 27., 21., 22., 125., 126., 127., 26.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[1:,3]=[225.,226.,227.]
        self.assertEqual([7., 8., 9., 10., 11., 12., 13., 14., 225., 16., 17., 18., 19., 226., 21., 22., 23., 24., 227., 26.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[1,2:]=[225,226,227]
        self.assertEqual([7., 8., 9., 10., 11., 12., 13., 225., 226., 227., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,-2:]=[88,99,1010,1111,1212,1313]
        self.assertEqual([7., 8., 9., 88., 99., 12., 13., 14., 15., 16., 17., 18., 19., 1010., 1111., 22., 23., 24., 1212., 1313.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da3=DataArrayDouble.New(); da3.setValues([88,99,1010,1111,1212,1313],3,2)
        da[da2,-2:]=da3
        self.assertEqual([7., 8., 9., 88., 99., 12., 13., 14., 15., 16., 17., 18., 19., 1010., 1111., 22., 23., 24., 1212., 1313.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,[0,2]]=da3
        self.assertEqual([88., 8., 99., 10., 11., 12., 13., 14., 15., 16., 1010., 18., 1111., 20., 21., 1212., 23., 1313., 25., 26.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,0:3:2]=da3
        self.assertEqual([88., 8., 99., 10., 11., 12., 13., 14., 15., 16., 1010., 18., 1111., 20., 21., 1212., 23., 1313., 25., 26.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,0:3:2]=-8.
        self.assertEqual([-8., 8., -8., 10., 11., 12., 13., 14., 15., 16., -8., 18., -8., 20., 21., -8., 23., -8., 25., 26.],da.getValues())
        pass

    def testSwigDADOp(self):
        da=DataArrayDouble.New()
        da.alloc(12,1)
        da.iota(7.)
        da1=DataArrayDouble.New()
        da1.alloc(12,1)
        da1.iota(8.)
        da2=da+da1
        self.assertEqual([15., 17., 19., 21., 23., 25., 27., 29., 31., 33., 35., 37.],da2.getValues())
        da2=da+3
        da3=3+da
        self.assertTrue(da2.isEqual(da3,1e-12))
        da2=da-1.
        self.assertEqual([6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0],da2.getValues())
        da2=1-da
        self.assertEqual([-6.0, -7.0, -8.0, -9.0, -10.0, -11.0, -12.0, -13.0, -14.0, -15.0, -16.0, -17.0],da2.getValues())
        da2=da*3
        self.assertEqual([21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0],da2.getValues())
        da2=3.*da
        self.assertEqual([21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0],da2.getValues())
        da2=da*da1
        self.assertEqual([56.0, 72.0, 90.0, 110.0, 132.0, 156.0, 182.0, 210.0, 240.0, 272.0, 306.0, 342.0],da2.getValues())
        da2=da/4.
        self.assertEqual([1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5],da2.getValues())
        da3=4./da
        da4=da3*da2
        self.assertTrue(da4.isUniform(1.,1e-12))
        st1=da.getHiddenCppPointer()
        da+=1
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertTrue(da.isEqual(da1,1e-12))
        da-=8
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual(list(range(12)), da.getValues())
        da+=da1
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0],da.getValues())
        da*=0.5
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0],da.getValues())
        da*=da1
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([32.0, 45.0, 60.0, 77.0, 96.0, 117.0, 140.0, 165.0, 192.0, 221.0, 252.0, 285.0],da.getValues())
        da/=da1
        self.assertEqual(st1,st2)
        self.assertEqual([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0],da.getValues())
        da/=2
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5],da.getValues())
        da.rearrange(3)
        da5=DataArrayDouble.New()
        da5.setValues([5.,4.,3.,2.],4,1)
        da*=da5 # it works with unmathing number of compo
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([10.0, 12.5, 15.0, 14.0, 16.0, 18.0, 15.0, 16.5, 18.0, 13.0, 14.0, 15.0],da.getValues())
        #
        da.alloc(30,1)
        da.iota(7.)
        da.rearrange(3)
        ids=DataArrayInt.New()
        ids.setValues([3,4,7],3,1)
        da[ids,:]=[5.,8.,9.]
        self.assertEqual([7.,8.,9.,10.,11.,12.,13.,14.,15.,5.,8.,9.,5.,8.,9.,22.,23.,24.,25.,26.,27.,5.,8.,9.,31.,32.,33.,34.,35.,36.0],da.getValues())
        #
        da.rearrange(1) ; da.iota(7) ; da.rearrange(3)
        da[ids,[1,2]]=[5,8]
        self.assertEqual([7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,5.,8.,19.,5.,8.,22.,23.,24.,25.,26.,27.,28.,5.,8.,31.,32.,33.,34.,35.,36.],da.getValues())
        pass

    def testSwigDAIOp(self):
        da=DataArrayInt.New()
        da.alloc(12,1)
        da.iota(7)
        da1=DataArrayInt.New()
        da1.alloc(12,1)
        da1.iota(8)
        da2=da+da1
        self.assertEqual([15,17,19,21,23,25,27,29,31,33,35,37],da2.getValues())
        da2=da+3
        da3=3+da
        self.assertTrue(da2.isEqual(da3))
        da2=da-1
        self.assertEqual([6,7,8,9,10,11,12,13,14,15,16,17],da2.getValues())
        da2=1-da
        self.assertEqual([-6,-7,-8,-9,-10,-11,-12,-13,-14,-15,-16,-17],da2.getValues())
        da2=da*3
        self.assertEqual([21,24,27,30,33,36,39,42,45,48,51,54.0],da2.getValues())
        da2=3*da
        self.assertEqual([21,24,27,30,33,36,39,42,45,48,51,54.0],da2.getValues())
        da2=da*da1
        self.assertEqual([56,72,90,110,132,156,182,210,240,272,306,342.0],da2.getValues())
        da2=da/4
        self.assertEqual([1,2,2,2,2,3,3,3,3,4,4,4],da2.getValues())
        da3=4/da
        da4=da3*da2
        self.assertTrue(da4.isUniform(0))
        st1=da.getHiddenCppPointer()
        da+=1
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertTrue(da.isEqual(da1))
        da-=8
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual(list(range(12)), da.getValues())
        da+=da1
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([8,10,12,14,16,18,20,22,24,26,28,30],da.getValues())
        da/=2
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([4,5,6,7,8,9,10,11,12,13,14,15],da.getValues())
        da*=da1
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([32,45,60,77,96,117,140,165,192,221,252,285],da.getValues())
        da/=da1
        self.assertEqual(st1,st2)
        self.assertEqual([4,5,6,7,8,9,10,11,12,13,14,15],da.getValues())
        da/=2
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([2,2, 3,3, 4,4, 5,5, 6,6, 7,7],da.getValues())
        da.rearrange(3)
        da5=DataArrayInt.New()
        da5.setValues([5,4,3,2],4,1)
        da*=da5 # it works with unmathing number of compo
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([10,10, 15,12,16,16,15,15, 18,12,14,14],da.getValues())
        da%=6
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([4,4,3,0,4,4,3,3,0,0,2,2],da.getValues())
        #
        da.alloc(30,1)
        da.iota(7)
        da.rearrange(3)
        ids=DataArrayInt.New()
        ids.setValues([3,4,7],3,1)
        da[ids,:]=[5,8,9]
        self.assertEqual([7,8,9,10,11,12,13,14,15,5,8,9,5,8,9,22,23,24,25,26,27,5,8,9,31,32,33,34,35,36],da.getValues())
        #
        da.rearrange(1) ; da.iota(7) ; da.rearrange(3)
        da[ids,[1,2]]=[5,8]
        self.assertEqual([7,8,9,10,11,12,13,14,15,16,5,8,19,5,8,22,23,24,25,26,27,28,5,8,31,32,33,34,35,36],da.getValues())
        pass

    def testSwigDAIOp2(self):
        da=DataArrayInt.New()
        st=da.getHiddenCppPointer()
        da.alloc(10,3)
        da.rearrange(1)
        da.iota(0)
        da.rearrange(3)
        da[:,1]+=4
        da[-2:,2]+=10
        da[-2:,2]+=10
        da[:,2]+=da[:,0]
        da[da[0],:]=7
        self.assertEqual(st,da.getHiddenCppPointer())
        self.assertEqual(da.getValues(),[7,7,7,3,8,8,7,7,7,9,14,20,12,17,26,7,7,7,18,23,38,21,26,44,24,29,70,27,32,76])
        pass

    def testSwigDAIOp3(self):
        da=DataArrayInt.New()
        self.assertRaises(InterpKernelException,da.__len__)
        self.assertRaises(InterpKernelException,da.__int__)
        for elt in da:
            self.assertTrue(False)
            pass
        da.alloc(12,3)
        da.rearrange(1) ; da.fillWithZero()
        l1=list(da)
        self.assertEqual(36,len(da));
        da.rearrange(3)
        tmp=da[0]
        self.assertRaises(InterpKernelException,tmp.__int__)
        self.assertEqual(12,len(da));
        l=list(da)
        for elt in enumerate(l):
            elt[1][2]=elt[0]
            pass
        ref=[0,0,0,0,0,1,0,0,2,0,0,3,0,0,4,0,0,5,0,0,6,0,0,7,0,0,8,0,0,9,0,0,10,0,0,11]
        self.assertEqual(ref,da.getValues());
        da.rearrange(1)
        l=[int(elt) for elt in l1]
        self.assertEqual(ref,da.getValues());
        self.assertEqual(11,int(da[-1:]))
        pass

    def testSwigDADOp3(self):
        da=DataArrayDouble.New()
        self.assertRaises(InterpKernelException,da.__len__)
        self.assertRaises(InterpKernelException,da.__float__)
        for elt in da:
            self.assertTrue(False)
            pass
        da.alloc(12,3)
        da.rearrange(1) ; da.fillWithZero()
        l1=list(da)
        self.assertEqual(36,len(da));
        da.rearrange(3)
        tmp=da[0]
        self.assertRaises(InterpKernelException,tmp.__float__)
        self.assertEqual(12,len(da));
        l=list(da)
        for elt in enumerate(l):
            elt[1][2]=elt[0]
            pass
        ref=[0.,0.,0.,0.,0.,1.,0.,0.,2.,0.,0.,3.,0.,0.,4.,0.,0.,5.,0.,0.,6.,0.,0.,7.,0.,0.,8.,0.,0.,9.,0.,0.,10.,0.,0.,11.]
        self.assertEqual(ref,da.getValues());
        da.rearrange(1)
        l=[float(elt) for elt in l1]
        self.assertEqual(ref,da.getValues());
        self.assertEqual(11.,float(da[-1:]))
        pass

    def testSwigDataArrayIntIterator1(self):
        da=DataArrayInt.New()
        da.alloc(12,1)
        da.iota(2)
        da.rearrange(3)
        # __getitem__ testing
        li=[]
        for it in da:
            li+=it[1:]
            pass
        self.assertEqual([3, 4, 6, 7, 9, 10, 12, 13],li)
        li=[]
        for it in da:
            li+=[it[-1]]
            pass
        self.assertEqual([4, 7, 10, 13],li)
        li=[]
        for it in da:
            li+=it[[2,1,0]]
            pass
        self.assertEqual([4, 3, 2, 7, 6, 5, 10, 9, 8, 13, 12, 11],li)
        # __setitem__ testing
        da3=da.deepCopy()
        da2=DataArrayInt.New()
        da2.alloc(12,1)
        da2.iota(2002)
        da2.rearrange(3)
        it2=da2.__iter__()
        i=0
        for it in da:
            pt = next(it2)
            it[:]=pt
            pass
        self.assertTrue(da.isEqual(da2))
        da=da3
        da3=da.deepCopy()
        #
        for it in da:
            it[:]=5
            pass
        da.rearrange(1)
        self.assertTrue(da.isUniform(5))
        da=da3
        da3=da.deepCopy()
        #
        for it in da:
            it[:]=[8,9,12]
            pass
        self.assertEqual([8, 9, 12, 8, 9, 12, 8, 9, 12, 8, 9, 12],da.getValues())
        da=da3
        da3=da.deepCopy()
        #
        for it in da:
            it[2]=[7]
            pass
        self.assertEqual([2, 3, 7, 5, 6, 7, 8, 9, 7, 11, 12, 7],da.getValues())
        pass

    def testSwigDataArrayDoubleIterator1(self):
        da=DataArrayDouble.New()
        da.alloc(12,1)
        da.iota(2)
        da.rearrange(3)
        # __getitem__ testing
        li=[]
        for it in da:
            li+=it[1:]
            pass
        self.assertEqual([3, 4, 6, 7, 9, 10, 12, 13],li)
        li=[]
        for it in da:
            li+=[it[-1]]
            pass
        self.assertEqual([4, 7, 10, 13],li)
        li=[]
        for it in da:
            li+=it[[2,1,0]]
            pass
        self.assertEqual([4, 3, 2, 7, 6, 5, 10, 9, 8, 13, 12, 11],li)
        # __setitem__ testing
        da3=da.deepCopy()
        da2=DataArrayDouble.New()
        da2.alloc(12,1)
        da2.iota(2002)
        da2.rearrange(3)
        it2=da2.__iter__()
        i=0
        for it in da:
            pt = next(it2)
            it[:]=pt
            pass
        self.assertTrue(da.isEqual(da2,1e-12))
        da=da3
        da3=da.deepCopy()
        #
        for it in da:
            it[:]=5
            pass
        da.rearrange(1)
        self.assertTrue(da.isUniform(5,1e-12))
        da=da3
        da3=da.deepCopy()
        #
        for it in da:
            it[:]=[8,9,12]
            pass
        self.assertEqual([8, 9, 12, 8, 9, 12, 8, 9, 12, 8, 9, 12],da.getValues())
        da=da3
        da3=da.deepCopy()
        #
        for it in da:
            it[2]=[7]
            pass
        self.assertEqual([2, 3, 7, 5, 6, 7, 8, 9, 7, 11, 12, 7],da.getValues())
        pass

    def testSwigUMeshIterator1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        li1=[]
        li2=[]
        for cell in m:
            li1+=cell.getAllConn()[1:]
            li2+=[cell.getType()]
            pass
        self.assertEqual(li1,[0, 3, 4, 1, 1, 4, 2, 4, 5, 2, 6, 7, 4, 3, 7, 8, 5, 4])
        self.assertEqual(li2,[4, 3, 3, 4, 4])
        pass

    def testSwigUMeshIterator2(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        self.assertRaises(InterpKernelException,m.cellsByType);
        m.rearrange2ConsecutiveCellTypes()
        li1=[]
        li2=[]
        li3=[]
        for cellsByType in m.cellsByType():
            li1.append(cellsByType.getType())
            li2.append(cellsByType.getNumberOfElems())
            temp=[]
            for cell in cellsByType:
                t=[None,None]
                t[0]=cell.getType()
                t[1]=cell.getAllConn()[1:]
                temp.append(t)
                pass
            li3.append(temp)
            pass
        self.assertEqual(li1,[4, 3])
        self.assertEqual(li2,[3, 2])
        self.assertEqual(li3,[[[4, (0, 3, 4, 1)], [4, (6, 7, 4, 3)], [4, (7, 8, 5, 4)]], [[3, (1, 4, 2)], [3, (4, 5, 2)]]])
        pass

    def testDAIAggregateMulti1(self):
        a=DataArrayInt.New()
        a.setValues(list(range(4)), 2, 2)
        a.setName("aa")
        b=DataArrayInt.New()
        b.setValues(list(range(6)), 3, 2)
        c=DataArrayInt.Aggregate([a,b])
        self.assertEqual(list(range(4)) + list(range(6)), c.getValues())
        self.assertEqual("aa",c.getName())
        self.assertEqual(5,c.getNumberOfTuples())
        self.assertEqual(2,c.getNumberOfComponents())
        pass

    def testMergeUMeshes2(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m2=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m3=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        #
        vec1=[0,2,3]
        m2_2=m2.buildPartOfMySelf(vec1,False);
        vec2=[1,1]
        m3_2=m3.buildPartOfMySelf(vec2,False);
        #
        ms=[m1,m2_2,m3_2];
        #
        self.assertRaises(InterpKernelException,MEDCouplingUMesh.MergeUMeshes,ms+[None]);
        self.assertRaises(InterpKernelException,MEDCouplingUMesh.MergeUMeshes,ms+[3.4])
        m4=MEDCouplingUMesh.MergeUMeshes(ms);
        m4.checkConsistencyLight();
        self.assertEqual(10,m4.getNumberOfCells());
        self.assertEqual(20,m4.getNumberOfNodes());
        self.assertEqual(45,m4.getNodalConnectivityArrayLen());
        m4bis=MEDCouplingMesh.MergeMeshes(ms);
        self.assertTrue(m4.isEqual(m4bis,1e-12))
        del m4bis
        #
        vec3=[0,1,2,3,4]
        m4_1=m4.buildPartOfMySelf(vec3,False);
        m4_1.setName(m1.getName());
        self.assertTrue(m4_1.isEqual(m1,1e-12));
        #
        vec4=[5,6,7]
        m4_2=m4.buildPartOfMySelf(vec4,False);
        cellCor,nodeCor=m4_2.checkGeoEquivalWith(m2_2,10,1e-12);
        #
        vec5=[8,9]
        m4_3=m4.buildPartOfMySelf(vec5,False);
        self.assertEqual(2,m4_3.getNumberOfCells());
        self.assertEqual(3,m4_3.getNumberOfNodes());
        m3_2.zipCoords();
        m4_3.setName(m3_2.getName());
        self.assertTrue(m4_3.isEqual(m3_2,1e-12));
        #
        pass

    def testBuild0DMeshFromCoords1(self):
        sourceCoords=[-0.3,-0.3,0., 0.7,-0.3,0., -0.3,0.7,0., 0.7,0.7,0.]
        coo=DataArrayDouble.New();
        coo.setValues(sourceCoords,4,3);
        coo.setName("My0D");
        m=MEDCouplingUMesh.Build0DMeshFromCoords(coo);
        m.checkConsistencyLight();
        self.assertEqual(4,m.getNumberOfNodes());
        self.assertEqual(4,m.getNumberOfCells());
        self.assertEqual(3,m.getSpaceDimension());
        self.assertEqual(0,m.getMeshDimension());
        types1=m.getAllGeoTypes();
        self.assertEqual([NORM_POINT1],types1);
        for i in range(4):
            conn=m.getNodeIdsOfCell(i);
            self.assertEqual([i],conn);
            self.assertTrue(NORM_POINT1==m.getTypeOfCell(i));
            pass
        self.assertEqual(m.getName(),"My0D");
        pass

    def testDescriptionInMeshTimeUnit1(self):
        text1="totoTTEDD";
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        m.setDescription(text1);
        self.assertEqual(m.getDescription(),text1);
        m2=m.deepCopy();
        self.assertTrue(m.isEqual(m2,1e-12));
        self.assertEqual(m2.getDescription(),text1);
        m2.setDescription("ggg");
        self.assertTrue(not m.isEqual(m2,1e-12));
        #
        f=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f.setTimeUnit(text1);
        self.assertEqual(f.getTimeUnit(),text1);
        f2=f.deepCopy();
        self.assertEqual(f2.getTimeUnit(),text1);
        #
        pass

    def testMultiFields1(self):
        mfs=MEDCouplingDataForTest.buildMultiFields_1();
        ms=mfs.getMeshes();
        dms,refs=mfs.getDifferentMeshes()
        das=mfs.getArrays();
        das2,refs2=mfs.getDifferentArrays()
        self.assertEqual(5,len(mfs.getFields()))
        self.assertEqual(1,len(mfs.getFields()[0].getArrays()));
        self.assertEqual(2,len(mfs.getFields()[1].getArrays()));
        self.assertEqual(1,len(mfs.getFields()[2].getArrays()));
        self.assertEqual(1,len(mfs.getFields()[3].getArrays()));
        self.assertEqual(1,len(mfs.getFields()[4].getArrays()));
        self.assertEqual(5,len(ms));
        self.assertEqual(2,len(dms));
        self.assertEqual(6,len(das));
        self.assertEqual(5,len(das2));
        mfs2=mfs.deepCopy();
        self.assertTrue(mfs.isEqual(mfs2,1e-12,1e-12))
        pass

    def testFieldOverTime1(self):
        fs=MEDCouplingDataForTest.buildMultiFields_2();
        self.assertRaises(InterpKernelException,MEDCouplingFieldOverTime.New,fs);
        f4bis=fs[4].buildNewTimeReprFromThis(ONE_TIME,False);
        fs[4]=f4bis;
        self.assertRaises(InterpKernelException,MEDCouplingFieldOverTime.New,fs);
        f4bis.setTime(2.7,20,21);
        fot=MEDCouplingFieldOverTime.New(fs);
        dt=fot.getDefinitionTimeZone();
        hs=dt.getHotSpotsTime();
        self.assertEqual(6,len(hs));
        expected1=[0.2,0.7,1.2,1.35,1.7,2.7]
        for i in range(6):
            self.assertAlmostEqual(expected1[i],hs[i],12);
            pass
        meshId,arrId,arrIdInField,fieldId=dt.getIdsOnTimeRight(0.2);
        self.assertEqual(0,meshId);
        self.assertEqual(0,arrId);
        self.assertEqual(0,arrIdInField);
        self.assertEqual(0,fieldId);
        #
        meshId,arrId,arrIdInField,fieldId=dt.getIdsOnTimeRight(0.7);
        self.assertEqual(0,meshId);
        self.assertEqual(1,arrId);
        self.assertEqual(0,arrIdInField);
        self.assertEqual(1,fieldId);
        #
        meshId,arrId,arrIdInField,fieldId=dt.getIdsOnTimeLeft(1.2);#**** WARNING left here
        self.assertEqual(0,meshId);
        self.assertEqual(2,arrId);
        self.assertEqual(1,arrIdInField);
        self.assertEqual(1,fieldId);
        #
        meshId,arrId,arrIdInField,fieldId=dt.getIdsOnTimeRight(1.2);#**** WARNING right again here
        self.assertEqual(1,meshId);
        self.assertEqual(3,arrId);
        self.assertEqual(0,arrIdInField);
        self.assertEqual(2,fieldId);
        #
        meshId,arrId,arrIdInField,fieldId=dt.getIdsOnTimeRight(1.35);
        self.assertEqual(1,meshId);
        self.assertEqual(3,arrId);
        self.assertEqual(0,arrIdInField);
        self.assertEqual(2,fieldId);
        #
        meshId,arrId,arrIdInField,fieldId=dt.getIdsOnTimeRight(1.7);
        self.assertEqual(0,meshId);
        self.assertEqual(3,arrId);
        self.assertEqual(0,arrIdInField);
        self.assertEqual(3,fieldId);
        #
        meshId,arrId,arrIdInField,fieldId=dt.getIdsOnTimeRight(2.7);
        self.assertEqual(1,meshId);
        self.assertEqual(4,arrId);
        self.assertEqual(0,arrIdInField);
        self.assertEqual(4,fieldId);
        #
        dt2=MEDCouplingDefinitionTime();
        self.assertTrue(not dt2.isEqual(dt));
        dt2.assign(dt);
        dt2.assign(dt);#to check memory management
        self.assertTrue(dt2.isEqual(dt));
        #
        dt3=MEDCouplingDefinitionTime();
        #
        pass

    def testDAICheckAndPreparePermutation1(self):
        vals1=[9,10,0,6,4,11,3,7];
        expect1=[5,6,0,3,2,7,1,4];
        vals2=[9,10,0,6,10,11,3,7];
        da=DataArrayInt.New();
        da.setValues(vals1,8,1);
        da2=da.checkAndPreparePermutation();
        self.assertEqual(8,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        for i in range(8):
            self.assertEqual(expect1[i],da2.getIJ(i,0));
            pass
        #
        da=DataArrayInt.New();
        da.alloc(8,1);
        da.iota(0);
        da2=da.checkAndPreparePermutation();
        self.assertEqual(1,da2.getNumberOfComponents());
        self.assertTrue(da2.isIota(8));
        #
        da=DataArrayInt.New();
        da.alloc(8,1);
        da.setValues(vals2,8,1);
        self.assertRaises(InterpKernelException,da.checkAndPreparePermutation);
        pass

    def testDAIChangeSurjectiveFormat1(self):
        vals1=[0,3,2,3,2,2,1,2]
        expected1=[0,1,2,6,8]
        expected2=[0,  6,  2,4,5,7,  1,3]
        da=DataArrayInt.New();
        da.setValues(vals1,8,1);
        #
        da2,da2I=da.changeSurjectiveFormat(4);
        self.assertEqual(5,da2I.getNumberOfTuples());
        self.assertEqual(8,da2.getNumberOfTuples());
        self.assertEqual(expected1,da2I.getValues());
        self.assertEqual(expected2,da2.getValues());
        #
        self.assertRaises(InterpKernelException,da.changeSurjectiveFormat,3);
        #
        pass

    def testUMeshGetCellIdsLyingOnNodes1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        nodeIds1=[1,2,3,4,6]
        nodeIds2=[6,7]
        da=m.getCellIdsLyingOnNodes(nodeIds1,True);
        self.assertEqual(1,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        self.assertEqual(1,da.getIJ(0,0));
        da2=DataArrayInt.New()
        da2.setValues(nodeIds2,2,1)
        da=m.getCellIdsLyingOnNodes(da2,False);
        self.assertEqual(2,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        self.assertEqual(3,da.getIJ(0,0));
        self.assertEqual(4,da.getIJ(1,0));
        pass

    def testUMeshFindCellIdsOnBoundary1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        da5=m.findCellIdsOnBoundary();
        self.assertTrue(da5.isIota(5));
        pass

    def testMeshSetTime1(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m2=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        #
        self.assertTrue(m1.isEqual(m2,1e-12));
        m1.setTime(3.14,6,7);
        tmp3,tmp1,tmp2=m1.getTime();
        self.assertEqual(6,tmp1);
        self.assertEqual(7,tmp2);
        self.assertAlmostEqual(3.14,tmp3,12);
        self.assertTrue(not m1.isEqual(m2,1e-12));
        m2.setTime(3.14,6,7);
        self.assertTrue(m1.isEqual(m2,1e-12));
        m1.setTimeUnit("ms");
        self.assertTrue(m1.getTimeUnit()=="ms");
        m1.setTimeUnit("us");
        self.assertTrue(m1.getTimeUnit()=="us");
        self.assertTrue(not m1.isEqual(m2,1e-12));
        m2.setTimeUnit("us");
        self.assertTrue(m1.isEqual(m2,1e-12));
        m2.setTime(3.14,6,8);
        self.assertTrue(not m1.isEqual(m2,1e-12));
        m2.setTime(3.14,7,7);
        self.assertTrue(not m1.isEqual(m2,1e-12));
        m2.setTime(3.15,6,7);
        self.assertTrue(not m1.isEqual(m2,1e-12));
        #
        m1.setTime(10.34,55,12);
        m3=m1.deepCopy();
        self.assertTrue(m1.isEqual(m3,1e-12));
        tmp3,tmp1,tmp2=m3.getTime();
        self.assertEqual(55,tmp1);
        self.assertEqual(12,tmp2);
        self.assertAlmostEqual(10.34,tmp3,12);
        #
        # testing CMesh
        coo1=[0.,1.,2.,3.5]
        a=DataArrayDouble.New();
        a.setValues(coo1,4,1);
        b=MEDCouplingCMesh.New();
        b.setCoordsAt(0,a);
        #
        b.setTime(5.67,8,100);
        tmp3,tmp1,tmp2=b.getTime();
        self.assertEqual(8,tmp1);
        self.assertEqual(100,tmp2);
        self.assertAlmostEqual(5.67,tmp3,12);
        c=b.deepCopy();
        self.assertTrue(c.isEqual(b,1e-12));
        tmp3,tmp1,tmp2=c.getTime();
        self.assertEqual(8,tmp1);
        self.assertEqual(100,tmp2);
        self.assertAlmostEqual(5.67,tmp3,12);
        pass

    def testApplyFuncTwo1(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setMesh(m1);
        #
        vals=[1.,11.,21.,2.,12.,22.,3.,13.,23.,4.,14.,24.,5.,15.,25.]
        da=DataArrayDouble.New();
        da.setValues(vals,5,3);
        f1.setArray(da);
        #
        self.assertRaises(InterpKernelException,da.applyFuncCompo,1,"y+z");
        da.setInfoOnComponent(0,"x [m]");
        da.setInfoOnComponent(1,"y [mm]");
        da.setInfoOnComponent(2,"z [km]");

        self.assertRaises(InterpKernelException, da.applyFuncCompo, 1, "x+y+zz+zzz");
        self.assertRaises(InterpKernelException, da.applyFuncCompo, 1, "toto(x+y)");
        self.assertRaises(InterpKernelException, da.applyFuncCompo, 1, "x/0");

        da2=da.applyFuncCompo(1,"y+z");
        self.assertEqual(1,da2.getNumberOfComponents());
        self.assertEqual(5,da2.getNumberOfTuples());
        expected1=[32.,34.,36.,38.,40.]
        for i in range(5):
            self.assertAlmostEqual(expected1[i],da2.getIJ(0,i),12);
            pass
        da2=da.applyFunc(1,"y+z");
        expected2=[12.,14.,16.,18.,20.]
        for i in range(5):
            self.assertAlmostEqual(expected2[i],da2.getIJ(0,i),12);
            pass
        #
        self.assertEqual(3,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        f1.applyFuncCompo(1,"y+z");
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        for i in range(5):
            self.assertAlmostEqual(expected1[i],f1.getArray().getIJ(0,i),12);
            pass
        #
        pass

    def testApplyFuncThree1(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setMesh(m1);
        #
        vals=[1.,11.,21.,2.,12.,22.,3.,13.,23.,4.,14.,24.,5.,15.,25.]
        da=DataArrayDouble.New();
        da.setValues(vals,5,3);
        f1.setArray(da);
        #
        vs=3*[None];
        vs[0]="x"; vs[1]="Y"; vs[2]="z";
        self.assertRaises(InterpKernelException, da.applyFuncNamedCompo, 1, vs, "y+z");
        self.assertRaises(InterpKernelException, da.applyFuncNamedCompo, 1, vs, "x+Y+z+zz+zzz");
        self.assertRaises(InterpKernelException, da.applyFuncNamedCompo, 1, vs, "x/0");
        vs[1]="y";
        da2=da.applyFuncNamedCompo(1,vs,"y+z");
        expected1=[32.,34.,36.,38.,40.]
        for i in range(5):
            self.assertAlmostEqual(expected1[i],da2.getIJ(0,i),12);
            pass
        self.assertRaises(InterpKernelException, da.applyFuncNamedCompo, 1, ["x","y","z","a"],"x+a")
        f1.setArray(da);
        self.assertEqual(3,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        f1.applyFuncNamedCompo(1,vs,"y+z");
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        for i in range(5):
            self.assertAlmostEqual(expected1[i],f1.getArray().getIJ(0,i),12);
            pass
        pass

    def testFillFromAnalyticTwo1(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m1.setTime(3.4,5,6); m1.setTimeUnit("us");
        self.assertRaises(InterpKernelException,m1.fillFromAnalyticCompo,ON_NODES,1,"y+z");
        m1.getCoords().setInfoOnComponent(0,"x [m]");
        m1.getCoords().setInfoOnComponent(1,"y");
        m1.getCoords().setInfoOnComponent(2,"z");
        f1=m1.fillFromAnalyticCompo(ON_NODES,1,"y+z");
        self.assertAlmostEqual(3.4,f1.getTime()[0],12) ; self.assertEqual(5,f1.getTime()[1]) ; self.assertEqual(6,f1.getTime()[2])
        self.assertEqual("us",f1.getTimeUnit())
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        expected1=[0.2, 0.7, 1.2, 0.7, 1.2, 1.7, 1.2, 1.7, 2.2]
        for i in range(9):
            self.assertAlmostEqual(expected1[i],f1.getArray().getIJ(0,i),12);
            pass
        pass

    def testFillFromAnalyticThree1(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m1.setTime(3.4,5,6); m1.setTimeUnit("us");
        vs=3*[None];
        vs[0]="x"; vs[1]="Y"; vs[2]="z";
        self.assertRaises(InterpKernelException,m1.fillFromAnalyticNamedCompo,ON_NODES,1,vs,"y+z");
        vs[1]="y";
        f1=m1.fillFromAnalyticNamedCompo(ON_NODES,1,vs,"y+z");
        self.assertAlmostEqual(3.4,f1.getTime()[0],12) ; self.assertEqual(5,f1.getTime()[1]) ; self.assertEqual(6,f1.getTime()[2])
        self.assertEqual("us",f1.getTimeUnit())
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        expected1=[0.2, 0.7, 1.2, 0.7, 1.2, 1.7, 1.2, 1.7, 2.2]
        for i in range(9):
            self.assertAlmostEqual(expected1[i],f1.getArray().getIJ(0,i),12);
            pass
        pass

    def testDAUnitVar1(self):
        da=DataArrayDouble.New();
        da.alloc(1,3);
        da.setInfoOnComponent(0,"XPS [m]");
        st1=da.getVarOnComponent(0);
        self.assertTrue(st1=="XPS");
        st2=da.getUnitOnComponent(0);
        self.assertTrue(st2=="m");
        #
        da.setInfoOnComponent(0,"XPS         [m]");
        st1=da.getVarOnComponent(0);
        self.assertTrue(st1=="XPS");
        st2=da.getUnitOnComponent(0);
        self.assertTrue(st2=="m");
        #
        da.setInfoOnComponent(0,"XPP         [m]");
        st1=da.getVarOnComponent(0);
        self.assertTrue(st1=="XPP");
        st2=da.getUnitOnComponent(0);
        self.assertTrue(st2=="m");
        #
        da.setInfoOnComponent(0,"XPP kdep  kefer   [ m  ]");
        st1=da.getVarOnComponent(0);
        self.assertTrue(st1=="XPP kdep  kefer");
        st2=da.getUnitOnComponent(0);
        self.assertTrue(st2==" m  ");
        #
        da.setInfoOnComponent(0,"     XPP k[  dep  k]efer   [ m^ 2/s^3*kJ  ]");
        st1=da.getVarOnComponent(0);
        self.assertTrue(st1=="     XPP k[  dep  k]efer");
        st2=da.getUnitOnComponent(0);
        self.assertTrue(st2==" m^ 2/s^3*kJ  ");
        #
        da.setInfoOnComponent(0,"     XPP kefer   ");
        st1=da.getVarOnComponent(0);
        self.assertTrue(st1=="     XPP kefer   ");
        st2=da.getUnitOnComponent(0);
        self.assertTrue(st2=="");
        #
        da.setInfoOnComponent(0,"temperature( bof)");
        st1=da.getVarOnComponent(0);
        self.assertTrue(st1=="temperature( bof)");
        st2=da.getUnitOnComponent(0);
        self.assertTrue(st2=="");
        #
        da.setInfoOnComponent(0,"kkk [m]");
        da.setInfoOnComponent(1,"ppp   [m^2/kJ]");
        da.setInfoOnComponent(2,"abcde   [MW/s]");
        #
        vs=da.getVarsOnComponent();
        self.assertEqual(3,len(vs));
        self.assertTrue(vs[0]=="kkk");
        self.assertTrue(vs[1]=="ppp");
        self.assertTrue(vs[2]=="abcde");
        vs=da.getUnitsOnComponent();
        self.assertEqual(3,len(vs));
        self.assertTrue(vs[0]=="m");
        self.assertTrue(vs[1]=="m^2/kJ");
        self.assertTrue(vs[2]=="MW/s");
        pass

    def testGaussCoordinates1(self):
        #Testing 1D cell types
        m1=MEDCouplingDataForTest.build1DMultiTypes_1();
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,ONE_TIME);
        f.setMesh(m1);
        wg1=[0.3];
        gsCoo1=[0.2];
        refCoo1=[-1.0,1.0];
        f.setGaussLocalizationOnType(NORM_SEG2,refCoo1,gsCoo1,wg1);
        wg2=wg1;
        gsCoo2=[0.2];
        refCoo2=[-1.0,1.0,0.0];
        f.setGaussLocalizationOnType(NORM_SEG3,refCoo2,gsCoo2,wg2);
        #
        resToTest=f.getLocalizationOfDiscr();
        self.assertEqual(3,resToTest.getNumberOfComponents());
        self.assertEqual(2,resToTest.getNumberOfTuples());
        expected1=[0.6,0.6,0.6, 0.6,0.6,0.6]
        for i in range(6):
            self.assertAlmostEqual(expected1[i],resToTest.getIJ(0,i),14);
            pass
        #
        #Testing 2D cell types
        m2=MEDCouplingDataForTest.build2DMultiTypes_1();
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,ONE_TIME);
        f.setMesh(m2);
        wg3=[0.3,0.3];
        tria3CooGauss=[ 0.1, 0.8, 0.2, 0.7 ]
        gsCoo3=tria3CooGauss
        tria3CooRef=[ 0.0, 0.0, 1.0 , 0.0, 0.0, 1.0 ]
        refCoo3=tria3CooRef;
        f.setGaussLocalizationOnType(NORM_TRI3,refCoo3,gsCoo3,wg3);
        wg4=[0.3,0.3,0.3];
        tria6CooGauss=[ 0.3, 0.2, 0.2, 0.1, 0.2, 0.4 ]
        gsCoo4=tria6CooGauss;
        tria6CooRef=[0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5]
        refCoo4=tria6CooRef;
        f.setGaussLocalizationOnType(NORM_TRI6,refCoo4,gsCoo4,wg4);
        wg5=[0.3,0.3,0.3,0.3];
        quad4CooGauss=[ 0.3, 0.2, 0.2, 0.1, 0.2, 0.4, 0.15, 0.27 ]
        gsCoo5=quad4CooGauss;
        quad4CooRef=[-1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0]
        refCoo5=quad4CooRef;
        f.setGaussLocalizationOnType(NORM_QUAD4,refCoo5,gsCoo5,wg5);
        wg6=[0.3,0.3,0.3,0.3];
        quad8CooGauss=[ 0.34, 0.16, 0.21, 0.3, 0.23, 0.4, 0.14, 0.37 ]
        gsCoo6=quad8CooGauss;
        quad8CooRef=[ -1.0, -1.0, 1.0, -1.0, 1.0,  1.0, -1.0,  1.0, 0.0, -1.0, 1.0,  0.0, 0.0,  1.0, -1.0,  0.0]
        refCoo6=quad8CooRef;
        f.setGaussLocalizationOnType(NORM_QUAD8,refCoo6,gsCoo6,wg6);
        #
        resToTest=f.getLocalizationOfDiscr();
        self.assertEqual(3,resToTest.getNumberOfComponents());
        self.assertEqual(13,resToTest.getNumberOfTuples());#2+3+4+4 gauss points for resp TRI3,TRI6,QUAD4,QUAD8
        expected2=[5.1,1.55,0.0, 4.7,1.65,0.0,
                   2.32,1.52,0.0, 1.6,1.32,0.0, 3.52,1.26,0.0,#TRI6
                   2.6,1.6,0.0, 2.4,1.8,0.0, 2.4,1.2,0.0, 2.3,1.46,0.0,#QUAD4
                   2.32,2.68,0.0, 2.6,2.42,0.0, 2.8,2.46,0.0, 2.74,2.28,0.0 ];#QUAD8
        for i in range(39):
            self.assertAlmostEqual(expected2[i],resToTest.getIJ(0,i),14);
            pass
        #
        #Testing 3D cell types
        m3=MEDCouplingDataForTest.build3DMultiTypes_1();
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,ONE_TIME);
        f.setMesh(m3);
        #
        wg7=[0.3];
        tetra4CooGauss=[0.34, 0.16, 0.21]
        gsCoo7=tetra4CooGauss;
        tetra4CooRef=[0.0,1.0,0.0, 0.0,0.0,1.0, 0.0,0.0,0.0, 1.0,0.0,0.0]
        refCoo7=tetra4CooRef;
        f.setGaussLocalizationOnType(NORM_TETRA4,refCoo7,gsCoo7,wg7);
        wg8=[0.3];
        tetra10CooGauss=[0.2, 0.3, 0.1]
        gsCoo8=tetra10CooGauss;
        tetra10CooRef=[0.0,1.0,0.0, 0.0,0.0,0.0, 0.0,0.0,1.0, 1.0,0.0,0.0, 0.0,0.5,0.0, 0.0,0.0,0.5, 0.0,0.5,0.5, 0.5,0.5,0.0, 0.5,0.0,0.0, 0.5,0.0,0.5]
        refCoo8=tetra10CooRef;
        f.setGaussLocalizationOnType(NORM_TETRA10,refCoo8,gsCoo8,wg8);
        wg9=[0.3];
        pyra5CooGauss=[0.2, 0.3, 0.1]
        gsCoo9=pyra5CooGauss;
        pyra5CooRef=[1.0,0.0,0.0, 0.0,1.0,0.0, -1.0,0.0,0.0, 0.0,-1.0,0.0, 0.0,0.0,1.0]
        refCoo9=pyra5CooRef;
        f.setGaussLocalizationOnType(NORM_PYRA5,refCoo9,gsCoo9,wg9);
        wg10=[0.3];
        pyra13CooGauss=[0.1, 0.2, 0.7]
        gsCoo10=pyra13CooGauss;
        pyra13CooRef=[1.0,0.0,0.0, 0.0,1.0,0.0,-1.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.5,0.5,0.0,-0.5,0.5,0.0,-0.5,-0.5,0.0,0.5,-0.5,0.0,0.5,0.0,0.5,0.0,0.5,0.5,-0.5,0.0,0.5,0.0,-0.5,0.5]
        refCoo10=pyra13CooRef;
        f.setGaussLocalizationOnType(NORM_PYRA13,refCoo10,gsCoo10,wg10);
        wg11=[0.3];
        penta6CooGauss=[0.2, 0.3, 0.1]
        gsCoo11=penta6CooGauss;
        penta6CooRef=[-1.0,1.0,0.0,-1.0,-0.0,1.0,-1.0,0.0,0.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,0.0,0.0]
        refCoo11=penta6CooRef;
        f.setGaussLocalizationOnType(NORM_PENTA6,refCoo11,gsCoo11,wg11);
        wg12=[0.3];
        penta15CooGauss=[0.2, 0.3,0.15]
        gsCoo12=penta15CooGauss;
        penta15CooRef=[-1.0,1.0,0.0,-1.0,0.0,1.0,-1.0,0.0,0.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,0.0,0.0,-1.0,0.5,0.5,-1.0,0.0,0.5,-1.0,0.5,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.5,0.5,1.0,0.0, 0.5,1.0,0.5,0.0]
        refCoo12=penta15CooRef;
        f.setGaussLocalizationOnType(NORM_PENTA15,refCoo12,gsCoo12,wg12);
        wg13=[0.3];
        hexa8CooGauss=[0.2,0.3,0.15]
        gsCoo13=hexa8CooGauss;
        hexa8CooRef=[-1.0,-1.0,-1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,1.0,1.0,1.0,1.0,-1.0,1.0,1.0]
        refCoo13=hexa8CooRef;
        f.setGaussLocalizationOnType(NORM_HEXA8,refCoo13,gsCoo13,wg13);
        wg14=[0.3];
        hexa20CooGauss=[0.11,0.3,0.55]
        gsCoo14=hexa20CooGauss;
        hexa20CooRef=[-1.0,-1.0,-1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,1.0,1.0,1.0,1.0,-1.0,1.0,1.0,0.0,-1.0,-1.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,0.0,-1.0,-1.0,-1.0,0.0,1.0,-1.0,0.0,1.0,1.0,0.0,-1.0,1.0,0.0,0.0,-1.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,-1.0,0.0,1.0]
        refCoo14=hexa20CooRef;
        f.setGaussLocalizationOnType(NORM_HEXA20,refCoo14,gsCoo14,wg14);
        #
        resToTest=f.getLocalizationOfDiscr();
        self.assertEqual(3,resToTest.getNumberOfComponents());
        self.assertEqual(8,resToTest.getNumberOfTuples());#2+3+4+4 gauss points for resp TRI3,TRI6,QUAD4,QUAD8
        expected3=[1.312,3.15,1.02, 0.56,3.3,0.6, 2.18,1.1,0.2, 1.18,1.54,0.98, 1.56,0.3,3.6, 1.613,0.801,4.374, 2.6,2.4,2.3, 2.31232,2.3933985,1.553255]
        for i in range(24):
            self.assertAlmostEqual(expected3[i],resToTest.getIJ(0,i),14);
            pass
        #
        pass

    def testP2Localization1(self):
        m=MEDCouplingUMesh.New("testP2",2);
        coords=[0.,2.,3.5,0.,4.5,1.5,1.2,0.32,3.4,1.,2.1,2.4]
        conn=[0,1,2,3,4,5]
        coo=DataArrayDouble.New();
        coo.setValues(coords,6,2);
        m.setCoords(coo);
        m.allocateCells(1);
        m.insertNextCell(NORM_TRI6,6,conn[0:6])
        m.finishInsertingCells();
        #
        f=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME);
        f.setMesh(m);
        da=DataArrayDouble.New();
        vals1=[1.2,2.3,3.4, 2.2,3.3,4.4, 3.2,4.3,5.4, 4.2,5.3,6.4, 5.2,6.3,7.4, 6.2,7.3,8.4]
        da.setValues(vals1,6,3);
        f.setArray(da);
        #
        loc=[2.27,1.3]
        locs=f.getValueOnMulti(loc);
        expected1=[6.0921164547752236, 7.1921164547752232, 8.2921164547752255]
        for i in range(3):
            self.assertAlmostEqual(expected1[i],locs.getIJ(0,i),12);
            pass
        pass

    def testP2Localization2(self):
        m=MEDCouplingUMesh.New("testP2_2",3);
        coords=[0.33312787792955395, -0.35155740179580952, -0.03567564825034563, 1.307146326477638, -0.57234557776250305, -0.08608044208272235, 0.5551834466499993, 0.62324964668794192, -0.014638951108536295, 0.37761817224442129, -0.38324019806913578, 0.96283164472856886, 0.79494856035658679, -0.40628057809270046, 0.0021004190225864614, 1.023740446371799, 0.07665912970471335, -0.072889657161871096, 0.54564584619517376, 0.11132872093429744, 0.039647326652013051, 0.27164784387819052, -0.42018012100866675, 0.46563376500745146, 0.89501965094896418, -0.56148455362735061, 0.43337469695473035, 0.49118025152924394, 0.093884938060727313, 0.47216346905220891]
        conn=[0,1,2,3,4,5,6,7,8,9]
        coo=DataArrayDouble.New();
        coo.setValues(coords,10,3);
        m.setCoords(coo);
        m.allocateCells(1);
        m.insertNextCell(NORM_TETRA10,10,conn[0:10])
        m.finishInsertingCells();
        #
        f=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME);
        f.setMesh(m);
        da=DataArrayDouble.New();
        vals1=[1.1,2.1,3.1,4.1,5.2,6.2,7.2,8.2,9.2,10.2]
        da.setValues(vals1,10,1);
        f.setArray(da);
        #
        loc=[0.64637931739890486, -0.16185896817550552, 0.22678966365273748]
        locs=f.getValueOnMulti(loc);
        expected1=[10.0844021968047]
        for i in range(1):
            self.assertAlmostEqual(expected1[i],locs.getIJ(0,i),12);
            pass
        pass

    def testGetValueOn2(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f.setMesh(m);
        arr=DataArrayDouble.New();
        nbOfCells=m.getNumberOfCells();
        f.setArray(arr);
        values1=[7.,107.,10007.,8.,108.,10008.,9.,109.,10009.,10.,110.,10010.,11.,111.,10011.]
        arr.setValues(values1,nbOfCells,3);
        loc=[-0.05,-0.05, 0.55,-0.25, 0.55,0.15, -0.05,0.45, 0.45,0.45]
        f.checkConsistencyLight();
        locs=f.getValueOnMulti(loc);
        self.assertEqual(5,locs.getNumberOfTuples());
        self.assertEqual(3,locs.getNumberOfComponents());
        for j in range(15):
            self.assertAlmostEqual(values1[j],locs.getIJ(0,j),12);
            pass
        # Testing ON_NODES
        f=MEDCouplingFieldDouble.New(ON_NODES,NO_TIME);
        f.setMesh(m);
        arr=DataArrayDouble.New();
        nbOfNodes=m.getNumberOfNodes();
        f.setArray(arr);
        values2=[7.,107.,10007.,8.,108.,10008.,9.,109.,10009.,10.,110.,10010.,11.,111.,10011.,12.,112.,10012.,13.,113.,10013.,14.,114.,10014.,15.,115.,10015.]
        arr.setValues(values2,nbOfNodes,3);
        loc2=[0.5432,-0.2432, 0.5478,0.1528, 0.5432,-0.2432, 0.5432,-0.2432]
        expected2=[9.0272, 109.0272, 10009.0272, 11.4124,111.4124,10011.4124, 9.0272, 109.0272, 10009.0272, 9.0272, 109.0272, 10009.0272]
        f.checkConsistencyLight();
        loc3=DataArrayDouble.New()
        loc3.setValues(loc2,4,2);
        locs=f.getValueOnMulti(loc3);
        self.assertEqual(4,locs.getNumberOfTuples());
        self.assertEqual(3,locs.getNumberOfComponents());
        for i in range(12):
            self.assertAlmostEqual(expected2[i],locs.getIJ(0,i),12);
            pass
        #
        pass

    def testDAIGetIdsNotEqual1(self):
        d=DataArrayInt.New();
        vals1=[2,3,5,6,8,5,5,6,1,-5]
        d.setValues(vals1,10,1);
        d2=d.findIdsNotEqual(5);
        self.assertEqual(7,d2.getNumberOfTuples());
        self.assertEqual(1,d2.getNumberOfComponents());
        expected1=[0,1,3,4,7,8,9]
        for i in range(7):
            self.assertEqual(expected1[i],d2.getIJ(0,i));
            pass
        d.rearrange(2);
        self.assertRaises(InterpKernelException,d.findIdsNotEqual,5);
        vals2=[-4,5,6]
        vals3=vals2;
        d.rearrange(1);
        d3=d.findIdsNotEqualList(vals3);
        self.assertEqual(5,d3.getNumberOfTuples());
        self.assertEqual(1,d3.getNumberOfComponents());
        expected2=[0,1,4,8,9]
        for i in range(5):
            self.assertEqual(expected2[i],d3.getIJ(0,i));
            pass
        pass

    def testDAIComputeOffsets1(self):
        d=DataArrayInt.New();
        vals1=[3,5,1,2,0,8]
        expected1=[0,3,8,9,11,11]
        d.setValues(vals1,6,1);
        d.computeOffsets();
        self.assertEqual(6,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        for i in range(6):
            self.assertEqual(expected1[i],d.getIJ(0,i));
            pass
        pass

    def testUMeshHexagonPrism1(self):
        coords=[0.8660254037844386, 0.5, 0.0, 0.0, 1.0, 0.0, -0.8660254037844386, 0.5, 0.0, -0.8660254037844386, -0.5, 0.0, 0.0, -1.0, 0.0, 0.8660254037844386, -0.5, 0.0,
                0.8660254037844386, 0.5, 2.0, 0.0, 1.0, 2.0, -0.8660254037844386, 0.5, 2.0, -0.8660254037844386, -0.5, 2.0, 0.0, -1.0, 2.0, 0.8660254037844386, -0.5, 2.0];
        conn=[1,2,3,4,5,0,7,8,9,10,11,6]
        mesh=MEDCouplingUMesh.New("MyFirstHexagonalPrism",3);
        coo=DataArrayDouble.New();
        coo.setValues(coords,12,3);
        mesh.setCoords(coo);
        mesh.allocateCells(1);
        mesh.insertNextCell(NORM_HEXGP12,12,conn[0:12])
        mesh.finishInsertingCells();
        #
        mesh.checkConsistencyLight();
        vols=mesh.getMeasureField(False);
        self.assertEqual(1,vols.getNumberOfTuples());
        self.assertEqual(1,vols.getNumberOfComponents());
        self.assertAlmostEqual(-5.196152422706632,vols.getIJ(0,0),12);
        bary=mesh.computeCellCenterOfMass();
        self.assertEqual(1,bary.getNumberOfTuples());
        self.assertEqual(3,bary.getNumberOfComponents());
        expected1=[0.,0.,1.]
        for i in range(3):
            self.assertAlmostEqual(expected1[i],bary.getIJ(0,i),12);
            pass
        d1=DataArrayInt.New();
        d2=DataArrayInt.New();
        d3=DataArrayInt.New();
        d4=DataArrayInt.New();
        m2=mesh.buildDescendingConnectivity(d1,d2,d3,d4);
        self.assertEqual(8,m2.getNumberOfCells());
        expected4=[[1,2,3,4,5,0],[7,6,11,10,9,8],[1,7,8,2],[2,8,9,3],[3,9,10,4],[4,10,11,5],[5,11,6,0],[0,6,7,1]];
        expected2=[NORM_POLYGON, NORM_POLYGON, NORM_QUAD4, NORM_QUAD4, NORM_QUAD4, NORM_QUAD4, NORM_QUAD4, NORM_QUAD4];
        expected3=[6,6,4,4,4,4,4,4]
        for i in range(8):
            self.assertTrue(m2.getTypeOfCell(i)==expected2[i]);
            v=m2.getNodeIdsOfCell(i);
            self.assertTrue(len(v)==expected3[i]);
            self.assertEqual(expected4[i],v);
        #
        mesh.convertAllToPoly();
        self.assertTrue(NORM_POLYHED==mesh.getTypeOfCell(0));
        mesh.unPolyze();
        self.assertTrue(NORM_HEXGP12==mesh.getTypeOfCell(0));
        self.assertEqual(13,mesh.getNodalConnectivityArrayLen());
        #
        pass

    def testDADCheckIsMonotonic(self):
        da=DataArrayDouble.New();
        da.setValues([-1.,1.01,2.03,6.],2,2);
        self.assertRaises(InterpKernelException,da.isMonotonic,True,1e-12);
        da.rearrange(1);
        self.assertTrue(da.isMonotonic(True,1e-12));
        da.checkMonotonic(True,1e-12);
        da.setIJ(2,0,6.1);
        self.assertTrue(not da.isMonotonic(True,1e-12));
        self.assertRaises(InterpKernelException,da.checkMonotonic,True,1e-12);
        da.setIJ(2,0,5.99);
        self.assertTrue(da.isMonotonic(True,1e-12));
        self.assertTrue(not da.isMonotonic(True,1e-1));
        pass

    def testCheckCoherencyDeeper1(self):
        m=MEDCouplingDataForTest.build3DSourceMesh_1();
        m.checkConsistencyLight();
        m.checkConsistency();
        m.getNodalConnectivity().setIJ(8,0,-1);
        m.checkConsistencyLight();
        self.assertRaises(InterpKernelException,m.checkConsistency);
        m.getNodalConnectivity().setIJ(8,0,-6);
        m.checkConsistencyLight();
        self.assertRaises(InterpKernelException,m.checkConsistency);
        m.getNodalConnectivity().setIJ(8,0,9);#9>=NbOfNodes
        m.checkConsistencyLight();
        self.assertRaises(InterpKernelException,m.checkConsistency);
        m.getNodalConnectivity().setIJ(8,0,8);#OK
        m.checkConsistencyLight();
        m.checkConsistency();
        elts=[1,5]
        m.convertToPolyTypes(elts);
        m.checkConsistencyLight();
        m.checkConsistency();
        m.getNodalConnectivity().setIJ(2,0,9);#9>=NbOfNodes
        m.checkConsistencyLight();
        self.assertRaises(InterpKernelException,m.checkConsistency);
        m.getNodalConnectivity().setIJ(2,0,-3);
        m.checkConsistencyLight();
        self.assertRaises(InterpKernelException,m.checkConsistency);
        m.getNodalConnectivity().setIJ(2,0,-1);
        m.checkConsistencyLight();
        self.assertRaises(InterpKernelException,m.checkConsistency);#Throw because cell#0 is not a polyhedron
        m.getNodalConnectivity().setIJ(2,0,4);
        m.checkConsistencyLight();
        m.checkConsistency();
        m.getNodalConnectivity().setIJ(7,0,-1);
        m.checkConsistencyLight();
        m.checkConsistency();#OK because we are in polyhedron connec
        m.getNodalConnectivity().setIJ(36,0,14);
        m.checkConsistencyLight();
        self.assertRaises(InterpKernelException,m.checkConsistency);#Throw beacause now cell 5 is a TETRA4 (14) so mimatch of number index and static type.
        pass

    def testUnPolyze2(self):
        m=MEDCouplingUMesh.New("jjj",3);
        coo=DataArrayDouble.New();
        coo.alloc(4,3);
        coo.rearrange(1);
        coo.iota(0);
        coo.rearrange(3);
        m.setCoords(coo);
        m.allocateCells(2);
        m.insertNextCell(NORM_TETRA4,4,[0,1,2,3]);
        m.insertNextCell(NORM_TETRA4,4,[0,1,2,3]);
        m.finishInsertingCells();
        m2=MEDCouplingUMesh.MergeUMeshesOnSameCoords(4*[m]);
        m2.convertToPolyTypes([2]);
        m2.unPolyze();
        self.assertEqual(NORM_TETRA4,m2.getTypeOfCell(2));
        self.assertEqual(40,m2.getNodalConnectivityArrayLen());
        temp2=m2.getNodeIdsOfCell(2);
        self.assertEqual(temp2,[0,1,2,3]);
        m2.checkConsistency();
        m3=m2.deepCopy();
        m2.unPolyze();
        self.assertTrue(m3.isEqual(m2,1e-12));
        pass

    def testDACpyFrom1(self):
        d=DataArrayDouble.New();
        d.alloc(12,1);
        d.iota(14.);
        d.rearrange(3);
        d.setName("Toto");
        d.setInfoOnComponent(0,"X [m]");
        d.setInfoOnComponent(1,"Y [m]");
        d.setInfoOnComponent(2,"Z [m]");
        #
        d1=DataArrayDouble.New();
        self.assertTrue(not d.isEqual(d1,1e-12));
        d1.deepCopyFrom(d);
        self.assertTrue(d.isEqual(d1,1e-12));
        d1.deepCopyFrom(d);
        self.assertTrue(d.isEqual(d1,1e-12));
        d1.rearrange(2);
        self.assertTrue(not d.isEqual(d1,1e-12));
        d1.deepCopyFrom(d);
        self.assertTrue(d.isEqual(d1,1e-12));
        #
        d2=d.convertToIntArr();
        d4=DataArrayInt.New();
        self.assertTrue(not d2.isEqual(d4));
        d4.deepCopyFrom(d2);
        self.assertTrue(d2.isEqual(d4));
        d4.deepCopyFrom(d2);
        self.assertTrue(d2.isEqual(d4));
        d4.rearrange(2);
        self.assertTrue(not d2.isEqual(d4));
        d4.deepCopyFrom(d2);
        self.assertTrue(d2.isEqual(d4));
        pass

    def testDAITransformWithIndArr1(self):
        tab1=[17,18,22,19]
        tab2=[0,1,1,3,3,0,1,3,2,2,3,0]
        expected=[17,18,18,19,19,17,18,19,22,22,19,17]
        d=DataArrayInt.New();
        d.setValues(tab1,4,1);
        d1=DataArrayInt.New();
        d1.setValues(tab2,12,1);
        d2=d1[:]
        #
        d1.transformWithIndArr(d);
        self.assertEqual(12,d1.getNumberOfTuples());
        self.assertEqual(1,d1.getNumberOfComponents());
        for i in range(12):
            self.assertEqual(expected[i],d1.getIJ(i,0));
            pass
        #
        d1=d2
        d1.transformWithIndArr(tab1)
        self.assertEqual(12,d1.getNumberOfTuples());
        self.assertEqual(1,d1.getNumberOfComponents());
        for i in range(12):
            self.assertEqual(expected[i],d1.getIJ(i,0));
            pass
        pass

    def testDAIBuildPermArrPerLevel1(self):
        arr=[2,0,1,1,0,1,2,0,1,1,0,0]
        expected1=[10,0,5,6,1,7,11,2,8,9,3,4]
        da=DataArrayInt.New();
        da.setValues(arr,12,1);
        da2=da.buildPermArrPerLevel();
        self.assertEqual(12,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        for i in range(12):
            self.assertEqual(expected1[i],da2.getIJ(i,0));
            pass
        pass

    def testDAIOperations1(self):
        arr1=[-1,-2,4,7,3,2,6,6,4,3,0,1]
        da=DataArrayInt.New();
        da.setValues(arr1,4,3);
        da1=DataArrayInt.New();
        da1.alloc(12,1);
        da1.iota(2);
        self.assertRaises(InterpKernelException,DataArrayInt.Add,da,da1);#not same number of tuples/Components
        da1.rearrange(3);
        da2=DataArrayInt.Add(da,da1);
        self.assertEqual(4,da2.getNumberOfTuples());
        self.assertEqual(3,da2.getNumberOfComponents());
        expected1=[1,1,8,12,9,9,14,15,14,14,12,14]
        for i in range(12):
            self.assertEqual(expected1[i],da2.getIJ(0,i));
            pass
        da1.substractEqual(da);
        expected2=[3,5,0,-2,3,5,2,3,6,8,12,12]
        for i in range(12):
            self.assertEqual(expected2[i],da1.getIJ(0,i));
            pass
        da1.rearrange(1); da1.iota(2); da1.rearrange(3);
        da1.addEqual(da);
        for i in range(12):
            self.assertEqual(expected1[i],da1.getIJ(0,i));
            pass
        da1.rearrange(1); da1.iota(2); da1.rearrange(3);
        da2=DataArrayInt.Multiply(da,da1);
        self.assertEqual(4,da2.getNumberOfTuples());
        self.assertEqual(3,da2.getNumberOfComponents());
        expected3=[-2,-6,16,35,18,14,48,54,40,33,0,13]
        for i in range(12):
            self.assertEqual(expected3[i],da2.getIJ(0,i));
            pass
        da.divideEqual(da1);
        self.assertEqual(4,da.getNumberOfTuples());
        self.assertEqual(3,da.getNumberOfComponents());
        expected4=[0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        for i in range(12):
            self.assertEqual(expected4[i],da.getIJ(0,i));
            pass
        da.setValues(arr1,4,3);
        da1.multiplyEqual(da);
        self.assertEqual(4,da1.getNumberOfTuples());
        self.assertEqual(3,da1.getNumberOfComponents());
        for i in range(12):
            self.assertEqual(expected3[i],da1.getIJ(0,i));
            pass
        da1.rearrange(1); da1.iota(2); da1.rearrange(3);
        da2=DataArrayInt.Divide(da,da1);
        self.assertEqual(4,da2.getNumberOfTuples());
        self.assertEqual(3,da2.getNumberOfComponents());
        for i in range(12):
            self.assertEqual(expected4[i],da2.getIJ(0,i));
            pass
        da1.applyInv(321);
        self.assertEqual(4,da1.getNumberOfTuples());
        self.assertEqual(3,da1.getNumberOfComponents());
        expected5=[160,107,80,64,53,45,40,35,32,29,26,24]
        for i in range(12):
            self.assertEqual(expected5[i],da1.getIJ(0,i));
            pass
        da1.applyDivideBy(2);
        self.assertEqual(4,da1.getNumberOfTuples());
        self.assertEqual(3,da1.getNumberOfComponents());
        expected6=[80,53,40,32,26,22,20,17,16,14,13,12]
        for i in range(12):
            self.assertEqual(expected6[i],da1.getIJ(0,i));
            pass
        expected7=[3,4,5,4,5,1,6,3,2,0,6,5]
        da1.applyModulus(7);
        for i in range(12):
            self.assertEqual(expected7[i],da1.getIJ(0,i));
            pass
        da1.applyLin(1,1);
        expected8=[3,3,3,3,3,1,3,3,0,0,3,3]
        da1.applyRModulus(3);
        for i in range(12):
            self.assertEqual(expected8[i],da1.getIJ(0,i));
            pass
        pass

    def testEmulateMEDMEMBDC1(self):
        m,m1=MEDCouplingDataForTest.buildPointe_1();
        m2,da1,da2,da3,da4,da5,da0=m.emulateMEDMEMBDC(m1)
        expected0=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,36,37,32,33,34,35,38,39,40,41,42,43,44,45,46]
        expected1=[1,32,29,23,41,36]
        self.assertEqual(47,da0.getNumberOfTuples());
        self.assertEqual(1,da0.getNumberOfComponents());
        for i in range(47):
            self.assertEqual(expected0[i],da0.getIJ(0,i));
            pass
        self.assertEqual(6,da5.getNumberOfTuples());
        self.assertEqual(1,da5.getNumberOfComponents());
        for i in range(6):
            self.assertEqual(expected1[i],da5.getIJ(0,i));
            pass
        expected2=[0,1,2,3,4,0,5,6,7,4,8,9,1,7,10,11,12,13,14,5,15,16,17,8,18,19,20,10,21,22,23,2,13,24,25,21,16,26,27,12,19,28,29,15,22,30,31,18,36,26,28,30,24,37,32,33,34,35,38,36,39,40,41,42,37,38,43,44,45,46]
        self.assertEqual(70,da1.getNumberOfTuples());
        self.assertEqual(1,da1.getNumberOfComponents());
        for i in range(70):
            self.assertEqual(expected2[i],da1.getIJ(0,i));
            pass
        expected3=[0,4,8,12,16,20,24,28,32,36,40,44,48,53,58,64,70]
        self.assertEqual(17,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        for i in range(17):
            self.assertEqual(expected3[i],da2.getIJ(0,i));
            pass
        expected4=[0,2,4,6,7,9,11,12,14,16,17,19,20,22,24,25,27,29,30,32,34,35,37,39,40,42,43,45,46,48,49,51,52,53,54,55,56,58,60,62,63,64,65,66,67,68,69,70]
        #expected4=[0,2,4,6,7,9,11,12,14,16,17,19,20,22,24,25,27,29,30,32,34,35,37,39,40,42,43,45,46,48,49,51,52,54,56,57,58,59,60,62,63,64,65,66,67,68,69,70];
        self.assertEqual(48,da4.getNumberOfTuples());
        self.assertEqual(1,da4.getNumberOfComponents());
        for i in range(48):
            self.assertEqual(expected4[i],da4.getIJ(0,i));
            pass
        expected5=[0,1,0,3,0,7,0,1,2,1,4,1,2,3,2,5,2,3,6,3,4,9,4,8,4,5,10,5,9,5,6,11,6,10,6,7,8,7,11,7,8,12,8,9,12,9,10,12,10,11,12,11,13,13,13,13,12,14,13,15,14,15,14,14,14,14,15,15,15,15]
        self.assertEqual(70,da3.getNumberOfTuples());
        self.assertEqual(1,da3.getNumberOfComponents());
        for i in range(70):
            self.assertEqual(expected5[i],da3.getIJ(0,i));
            pass
        pass

    def testGetLevArrPerCellTypes1(self):
        m,m1=MEDCouplingDataForTest.buildPointe_1();
        m1,d0,d1,d2,d3=m.buildDescendingConnectivity();
        order=[NORM_TRI3,NORM_QUAD4];
        da0,da1=m1.getLevArrPerCellTypes(order);
        expected0=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1]
        expected1=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,36,37,32,33,34,35,38,39,40,41,42,43,44,45,46]
        self.assertEqual(47,da0.getNumberOfTuples());
        self.assertEqual(1,da0.getNumberOfComponents());
        for i in range(47):
            self.assertEqual(expected0[i],da0.getIJ(0,i));
            pass
        self.assertEqual(2,da1.getNumberOfTuples());
        self.assertEqual(1,da1.getNumberOfComponents());
        self.assertEqual(36,da1.getIJ(0,0));#36 TRI3
        self.assertEqual(11,da1.getIJ(1,0));#11 QUAD4
        #
        da2=da0.buildPermArrPerLevel();
        #
        self.assertEqual(47,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        for i in range(47):
            self.assertEqual(expected1[i],da2.getIJ(0,i));
            pass
        pass

    def testSortCellsInMEDFileFrmt1(self):
        m,m1=MEDCouplingDataForTest.buildPointe_1();
        m2=m.deepCopy()
        da=DataArrayInt.New()
        da.setValues([0,1,2,14,3,12,4,5,15,6,7,8,9,10,11,13],16,1)
        daa=da.invertArrayN2O2O2N(16)
        m.renumberCells(daa,False)
        da2=m.sortCellsInMEDFileFrmt()
        self.assertEqual(da2.getValues(),[0,1,2,14,3,12,4,5,15,6,7,8,9,10,11,13])
        self.assertTrue(m.isEqual(m2,1e-12))
        self.assertTrue(da.isEqual(da2))
        pass

    def testBuildPartAndReduceNodes1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        arr=[1,0]
        m2,da=m.buildPartAndReduceNodes(arr);
        self.assertEqual(5,m2.getNumberOfNodes());
        self.assertEqual(2,m2.getNumberOfCells());
        f=m2.getMeasureField(True);
        self.assertAlmostEqual(0.125,f.getArray().getIJ(0,0),12);
        self.assertAlmostEqual(0.25,f.getArray().getIJ(1,0),12);
        #
        arr2=DataArrayInt.New()
        arr2.setValues(arr,2,1)
        m2,da=m.buildPartAndReduceNodes(arr2);
        self.assertEqual(5,m2.getNumberOfNodes());
        self.assertEqual(2,m2.getNumberOfCells());
        f=m2.getMeasureField(True);
        self.assertAlmostEqual(0.125,f.getArray().getIJ(0,0),12);
        self.assertAlmostEqual(0.25,f.getArray().getIJ(1,0),12);
        pass

    def testDAITransformWithIndArrR1(self):
        tab1=[2,4,5,3,6,7]
        tab2=[-1,-1,0,1,2,3,4,5,-1,-1,-1,-1]
        expected=[0,3,1,2,4,5]
        d=DataArrayInt.New();
        d.setValues(tab1,6,1);
        d1=DataArrayInt.New();
        d1.setValues(tab2,12,1);
        d2=d1[:]
        #
        d3=d.transformWithIndArrR(d1);
        self.assertEqual(6,d3.getNumberOfTuples());
        self.assertEqual(1,d3.getNumberOfComponents());
        for i in range(6):
            self.assertEqual(expected[i],d3.getIJ(i,0));
            pass
        #
        d1=d2
        d3=d.transformWithIndArrR(tab2)
        self.assertEqual(6,d3.getNumberOfTuples());
        self.assertEqual(1,d3.getNumberOfComponents());
        for i in range(6):
            self.assertEqual(expected[i],d3.getIJ(i,0));
            pass
        pass

    def testDAISplitByValueRange1(self):
        val1=[6,5,0,3,2,7,8,1,4]
        val2=[0,4,9]
        d=DataArrayInt.New();
        d.setValues(val1,9,1);
        e,f,g=d.splitByValueRange(val2);
        self.assertEqual(9,e.getNumberOfTuples());
        self.assertEqual(1,e.getNumberOfComponents());
        self.assertEqual(9,f.getNumberOfTuples());
        self.assertEqual(1,f.getNumberOfComponents());
        self.assertEqual(2,g.getNumberOfTuples());
        self.assertEqual(1,g.getNumberOfComponents());
        #
        expected1=[1,1,0,0,0,1,1,0,1]
        expected2=[2,1,0,3,2,3,4,1,0]
        for i in range(9):
            self.assertEqual(expected1[i],e.getIJ(i,0));
            self.assertEqual(expected2[i],f.getIJ(i,0));
            pass
        self.assertEqual(0,g.getIJ(0,0));
        self.assertEqual(1,g.getIJ(1,0));
        #
        d.setIJ(6,0,9);
        self.assertRaises(InterpKernelException,d.splitByValueRange,val2);
        # non regression test in python wrapping
        rg=DataArrayInt([0,10,29,56,75,102,121,148,167,194,213,240,259,286,305,332,351,378,397,424,443,470,489,516])
        a,b,c=DataArrayInt([75]).splitByValueRange(rg)
        assert(a.isEqual(DataArrayInt([4])))
        assert(b.isEqual(DataArrayInt([0])))
        assert(c.isEqual(DataArrayInt([4])))
        pass

    def testUMeshSplitProfilePerType1(self):
        val0=[2,0,1,3,4]
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        m.renumberCells(val0,False);
        #
        val1=[0,2,3]
        d=DataArrayInt.New();
        d.setValues(val1,3,1);
        d.setName("sup")
        code,idsInPflPerType,pfls=m.splitProfilePerType(d);
        self.assertEqual(2,len(code));
        self.assertEqual(2,len(idsInPflPerType));
        expected1=[[3,1,0], [4,2,1]]
        self.assertEqual(expected1,code)
        self.assertEqual(2,len(idsInPflPerType));
        self.assertEqual(1,idsInPflPerType[0].getNumberOfTuples());
        self.assertEqual(0,idsInPflPerType[0].getIJ(0,0));
        self.assertEqual(2,idsInPflPerType[1].getNumberOfTuples());
        self.assertEqual(1,idsInPflPerType[1].getIJ(0,0));
        self.assertEqual(2,idsInPflPerType[1].getIJ(1,0));
        #
        self.assertEqual(2,len(pfls));
        self.assertEqual("sup",pfls[0].getName())
        self.assertEqual(1,pfls[0].getNumberOfTuples());
        self.assertEqual(0,pfls[0].getIJ(0,0));
        self.assertEqual("sup",pfls[1].getName())
        self.assertEqual(2,pfls[1].getNumberOfTuples());
        self.assertEqual(0,pfls[1].getIJ(0,0));
        self.assertEqual(1,pfls[1].getIJ(1,0));
        #
        val2=[0,2,3,4]
        d=DataArrayInt.New();
        d.setValues(val2,4,1);
        code,idsInPflPerType,pfls=m.splitProfilePerType(d);
        self.assertEqual(2,len(code));
        self.assertEqual(2,len(idsInPflPerType));
        expected2=[[3,1,0], [4,3,-1]]
        self.assertEqual(expected2,code);
        self.assertEqual(2,len(idsInPflPerType));
        self.assertEqual(1,idsInPflPerType[0].getNumberOfTuples());
        self.assertEqual(0,idsInPflPerType[0].getIJ(0,0));
        self.assertEqual(3,idsInPflPerType[1].getNumberOfTuples());
        self.assertEqual(1,idsInPflPerType[1].getIJ(0,0));
        self.assertEqual(2,idsInPflPerType[1].getIJ(1,0));
        self.assertEqual(3,idsInPflPerType[1].getIJ(2,0));
        #
        self.assertEqual(1,len(pfls));
        self.assertEqual(1,pfls[0].getNumberOfTuples());
        self.assertEqual(0,pfls[0].getIJ(0,0));
        #
        val3=[1,0,2]
        d=DataArrayInt.New();
        d.setValues(val3,3,1);
        code,idsInPflPerType,pfls=m.splitProfilePerType(d);
        self.assertEqual(2,len(code));
        self.assertEqual(2,len(idsInPflPerType));
        expected3=[[3,2,0], [4,1,1]]
        self.assertEqual(expected3,code);
        self.assertEqual(2,len(idsInPflPerType));
        self.assertEqual(2,idsInPflPerType[0].getNumberOfTuples());
        self.assertEqual(0,idsInPflPerType[0].getIJ(0,0));
        self.assertEqual(1,idsInPflPerType[0].getIJ(1,0));
        self.assertEqual(1,idsInPflPerType[1].getNumberOfTuples());
        self.assertEqual(2,idsInPflPerType[1].getIJ(0,0));
        #
        self.assertEqual(2,len(pfls));
        self.assertEqual(2,pfls[0].getNumberOfTuples());
        self.assertEqual(1,pfls[0].getIJ(0,0));
        self.assertEqual(0,pfls[0].getIJ(1,0));
        self.assertEqual(0,pfls[1].getIJ(0,0));
        #
        val4=[3,4]
        d=DataArrayInt.New();
        d.setValues(val4,2,1);
        code,idsInPflPerType,pfls=m.splitProfilePerType(d);
        self.assertEqual(1,len(code));
        self.assertEqual(1,len(idsInPflPerType));
        expected4=[[4,2,0]]
        self.assertEqual(expected4,code);
        self.assertEqual(1,len(idsInPflPerType));
        self.assertEqual(2,idsInPflPerType[0].getNumberOfTuples());
        self.assertEqual(0,idsInPflPerType[0].getIJ(0,0));
        self.assertEqual(1,idsInPflPerType[0].getIJ(1,0));
        #
        self.assertEqual(1,len(pfls));
        self.assertEqual(2,pfls[0].getNumberOfTuples());
        self.assertEqual(1,pfls[0].getIJ(0,0));
        self.assertEqual(2,pfls[0].getIJ(1,0));
        pass

    def testDAIBuildExplicitArrByRanges1(self):
        d=DataArrayInt.New();
        vals1=[0,2,3]
        d.setValues(vals1,3,1);
        e=DataArrayInt.New();
        vals2=[0,3,6,10,14,20]
        e.setValues(vals2,6,1);
        #
        f=d.buildExplicitArrByRanges(e);
        self.assertEqual(11,f.getNumberOfTuples());
        self.assertEqual(1,f.getNumberOfComponents());
        expected1=[0,1,2,6,7,8,9,10,11,12,13]
        for i in range(11):
            self.assertEqual(expected1[i],f.getIJ(i,0));
            pass
        pass

    def testDAIComputeOffsets2(self):
        d=DataArrayInt.New();
        vals1=[3,5,1,2,0,8]
        expected1=[0,3,8,9,11,11,19]
        d.setValues(vals1,6,1);
        d.computeOffsetsFull();
        self.assertEqual(7,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        for i in range(7):
            self.assertEqual(expected1[i],d.getIJ(0,i));
            pass
        pass

    def testMergeField3(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        m.getCoords().setInfoOnComponent(0,"x [m]");
        m.getCoords().setInfoOnComponent(1,"z [km]");
        m.setName("m");
        m.setDescription("desc");
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setName("f1");
        f1.setMesh(m);
        arr=DataArrayDouble.New();
        arr.alloc(5,2);
        arr.setInfoOnComponent(0,"X [m]");
        arr.setInfoOnComponent(1,"YY [mm]");
        arr.fillWithValue(2.);
        f1.setArray(arr);
        #
        f2=MEDCouplingFieldDouble.MergeFields([f1]);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        #
        pass
    
    def testGetDistributionOfTypes1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        tab1=[2,0,1,3,4]
        self.assertRaises(InterpKernelException,m.getDistributionOfTypes);
        m.renumberCells(tab1,False);
        code=m.getDistributionOfTypes();
        self.assertEqual(2,len(code));
        self.assertEqual(3,code[0][0]);
        self.assertEqual(2,code[0][1]);
        self.assertEqual(-1,code[0][2]);
        self.assertEqual(4,code[1][0]);
        self.assertEqual(3,code[1][1]);
        self.assertEqual(-1,code[1][2]);
        pass

    def testNorm2_1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f.setMesh(m);
        #
        d=DataArrayDouble.New();
        tab=[1.2,1.3,2.2,2.3,3.2,3.3,4.2,4.3,5.2,5.3]
        d.setValues(tab,5,2);
        f.setArray(d);
        f.checkConsistencyLight();
        #
        self.assertAlmostEqual(11.209371079592289,f.norm2(),14);
        #
        pass

    def testNormMax1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f.setMesh(m);
        #
        d=DataArrayDouble.New();
        tab=[2.3,-1.2,6.3,-7.8,2.9,7.7,2.1,0.,3.6,-7.6]
        d.setValues(tab,5,2);
        f.setArray(d);
        f.checkConsistencyLight();
        #
        self.assertAlmostEqual(7.8,f.normMax(),14);
        #
        pass

    def testFindAndCorrectBadOriented3DExtrudedCells1(self):
        coords=[0.0011180339887498999, -0.0011755705045849499, 0.0, -0.0012331070204200001, -0.0011755705045849499, 0.0, -0.00067557050458494599, -0.00145964954842536, 0.0, -0.00050000000000000001, -0.00086602540378443902, 0.0, 0.00140211303259031, -0.00061803398874989504, 0.0, 0.00086602540378443902, -0.00050000000000000001, 0.0, 0.001, 0.0, 0.0, 0.00034561537182258202, 0.000269164072574575, 0.0, 0.0, 0.001, 0.0, -0.00050000000000000001, 0.00086602540378443902, 0.0, -0.000269164072574575, 0.00034561537182258202, 0.0, -0.001, 0.0, 0.0, -0.00086602540378443902, -0.00050000000000000001, 0.0, -0.00034561537182258202, -0.000269164072574575, 0.0, 0.0, -0.001, 0.0, 0.00050000000000000001, -0.00086602540378443902, 0.0, 0.000269164072574575, -0.00034561537182258202, 0.0, 0.0015, -6.01853107621011e-36, 0.0, 0.00056049747291484397, -0.00145964954842536, 0.0, 0.0011180339887498999, -0.0011755705045849499, 0.00050000000000000001, -0.0012331070204200001, -0.0011755705045849499, 0.00050000000000000001, -0.00067557050458494599, -0.00145964954842536, 0.00050000000000000001, -0.00050000000000000001, -0.00086602540378443902, 0.00050000000000000001, 0.00140211303259031, -0.00061803398874989504, 0.00050000000000000001, 0.00086602540378443902, -0.00050000000000000001, 0.00050000000000000001, 0.001, 0.0, 0.00050000000000000001, 0.00034561537182258202, 0.000269164072574575, 0.00050000000000000001, 0.0, 0.001, 0.00050000000000000001, -0.00050000000000000001, 0.00086602540378443902, 0.00050000000000000001, -0.000269164072574575, 0.00034561537182258202, 0.00050000000000000001, -0.001, 0.0, 0.00050000000000000001, -0.00086602540378443902, -0.00050000000000000001, 0.00050000000000000001, -0.00034561537182258202, -0.000269164072574575, 0.00050000000000000001, 0.0, -0.001, 0.00050000000000000001, 0.00050000000000000001, -0.00086602540378443902, 0.00050000000000000001, 0.000269164072574575, -0.00034561537182258202, 0.00050000000000000001, 0.0015, -6.01853107621011e-36, 0.00050000000000000001, 0.00056049747291484397, -0.00145964954842536, 0.00050000000000000001];
        conn=[2, 1, 3, 21, 20, 22, 4, 0, 5, 23, 19, 24, 8, 9, 10, 27, 28, 29, 11, 12, 13, 30, 31, 32, 0, 18, 15, 5, 19, 37, 34, 24, 6, 17, 4, 5, 25, 36, 23, 24, 3, 14, 16, 13, 22, 33, 35, 32, 13, 16, 7, 10, 32, 35, 26, 29]
        connExp=[16, 2, 1, 3, 21, 20, 22, 16, 4, 0, 5, 23, 19, 24, 16, 8, 10, 9, 27, 29, 28, 16, 11, 13, 12, 30, 32, 31, 18, 0, 18, 15, 5, 19, 37, 34, 24,18, 6, 17, 4, 5, 25, 36, 23, 24, 18, 3, 13, 16, 14, 22, 32, 35, 33, 18, 13, 10, 7, 16, 32, 29, 26, 35]
        invalidCells=[2,3,6,7]
        m=MEDCouplingUMesh.New("Example",3);
        coo=DataArrayDouble.New();
        coo.setValues(coords,38,3);
        m.setCoords(coo);
        m.allocateCells(8);
        m.insertNextCell(NORM_PENTA6,6,conn[0:6])
        m.insertNextCell(NORM_PENTA6,6,conn[6:12])
        m.insertNextCell(NORM_PENTA6,6,conn[12:18])
        m.insertNextCell(NORM_PENTA6,6,conn[18:24])
        m.insertNextCell(NORM_HEXA8,8,conn[24:32])
        m.insertNextCell(NORM_HEXA8,8,conn[32:40])
        m.insertNextCell(NORM_HEXA8,8,conn[40:48])
        m.insertNextCell(NORM_HEXA8,8,conn[48:56])
        m.finishInsertingCells();
        #
        v=m.findAndCorrectBadOriented3DExtrudedCells();
        self.assertEqual(4,len(v));
        self.assertEqual(v.getValues(),invalidCells);
        self.assertEqual(connExp,m.getNodalConnectivity().getValues());
        self.assertTrue(m.findAndCorrectBadOriented3DExtrudedCells().empty())
        #
        pass

    def testConvertExtrudedPolyhedra1(self):
        conn=[1,2,3,4, 5,6,7,8,9,10,11,12, 13,14,15,16, 17,18,19,20,21,22, 23,24,25,26,27,28, 29,30,31,32,33,34,35,36,37,38, 39,40,41,42,43,44,45,46, 47,48,49,50,51,52,53,54,55,56,57,58, 59,60,61,62,63,64,65,66,67,68,69,70,71,72]
        m=MEDCouplingUMesh.New("Example",3);
        coo=DataArrayDouble.New();
        coo.alloc(73,3);
        coo.rearrange(1); coo.iota(0); coo.rearrange(3);
        m.setCoords(coo);
        m.allocateCells(9);
        m.insertNextCell(NORM_TETRA4,4,conn[0:4])
        m.insertNextCell(NORM_HEXA8,8,conn[4:12])
        m.insertNextCell(NORM_TETRA4,4,conn[12:16])
        m.insertNextCell(NORM_POLYHED,6,conn[16:22])
        m.insertNextCell(NORM_PENTA6,6,conn[22:28])
        m.insertNextCell(NORM_POLYHED,10,conn[28:38])
        m.insertNextCell(NORM_HEXA8,8,conn[38:46])
        m.insertNextCell(NORM_HEXGP12,12,conn[46:58])
        m.insertNextCell(NORM_POLYHED,14,conn[58:72])
        m.finishInsertingCells();
        #
        m.convertExtrudedPolyhedra();
        da=m.getNodalConnectivity();
        dai=m.getNodalConnectivityIndex();
        self.assertEqual(10,dai.getNbOfElems());
        self.assertEqual(159,da.getNbOfElems());
        #
        expected1=[14,1,2,3,4,18,5,6,7,8,9,10,11,12,14,13,14,15,16,31,17,18,19,-1,20,22,21,-1,17,20,21,18,-1,18,21,22,19,-1,19,22,20,17,16,23,24,25,26,27,28,31,29,30,31,32,33,-1,34,38,37,36,35,-1,29,34,35,30,-1,30,35,36,31,-1,31,36,37,32,-1,32,37,38,33,-1,33,38,34,29,18,39,40,41,42,43,44,45,46,22,47,48,49,50,51,52,53,54,55,56,57,58,31,59,60,61,62,63,64,65,-1,66,72,71,70,69,68,67,-1,59,66,67,60,-1,60,67,68,61,-1,61,68,69,62,-1,62,69,70,63,-1,63,70,71,64,-1,64,71,72,65,-1,65,72,66,59];
        expected2=[0,5,14,19,42,49,86,95,108,159]
        self.assertEqual(expected1,da.getValues());
        self.assertEqual(expected2,dai.getValues());
        m.checkConsistency()
        pass

    def testNonRegressionCopyTinyStrings(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f1=m.getMeasureField(True)
        f1.getArray().setInfoOnComponent(0,"P [N/m^2]")
        bary=m.computeCellCenterOfMass()
        f2=f1.buildNewTimeReprFromThis(NO_TIME,False)
        f2.setArray(bary)
        self.assertRaises(InterpKernelException,f1.copyTinyAttrFrom,f2)
        pass

    def testDaDSetPartOfValuesAdv1(self):
        tab1=[3.,4.,5., 13.,14.,15., 23.,24.,25., 33.,34.,35., 43.,44.,45., 53.,54.,55.]
        tab2=[6.,7.,8., 16.,17.,18., 26.,27.,28.]
        tab3=[4,1, 2,2, 3,0]
        a=DataArrayDouble.New();
        a.setValues(tab1,6,3);
        b=DataArrayDouble.New();
        b.setValues(tab2,3,3);
        c=DataArrayInt.New();
        c.setValues(tab3,3,2);
        #
        a.setPartOfValuesAdv(b,c);
        expected1=[3.,4.,5., 13.,14.,15., 26.,27.,28., 6.,7.,8., 16.,17.,18., 53.,54.,55.]
        self.assertEqual(expected1,a.getValues());
        pass

    def testUMeshBuildSetInstanceFromThis1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m2=m.buildSetInstanceFromThis(3);
        self.assertTrue(m.isEqual(m2,1e-12));
        #
        m=MEDCouplingUMesh.New("toto",2);
        m2=m.buildSetInstanceFromThis(3);
        self.assertEqual(0,m2.getNumberOfNodes());
        self.assertEqual(0,m2.getNumberOfCells());
        pass

    def testUMeshMergeMeshesCVW1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m2=MEDCouplingUMesh.New("toto",2);
        m3=MEDCouplingUMesh.MergeUMeshes([m,m2]);
        m3.setName(m.getName());
        self.assertTrue(m.isEqual(m3,1e-12));
        pass
    
    def testChangeUnderlyingMeshWithCMesh1(self):
        mesh=MEDCouplingCMesh.New();
        coordsX=DataArrayDouble.New();
        arrX=[ -1., 1., 2., 4. ]
        coordsX.setValues(arrX,4,1);
        coordsY=DataArrayDouble.New();
        arrY=[ -2., 2., 4., 8. ]
        coordsY.setValues(arrY,4,1);
        coordsZ=DataArrayDouble.New();
        arrZ=[ -3., 3., 6., 12. ]
        coordsZ.setValues(arrZ,4,1);
        mesh.setCoords(coordsX,coordsY,coordsZ);
        f=mesh.getMeasureField(True)
        mesh2=mesh.deepCopy()
        for myId in [0,1,2,10,11,12,20,21,22]:
            f=mesh.getMeasureField(True)
            f.changeUnderlyingMesh(mesh2,myId,1e-12);
            pass
        mesh2.setName("uuuu")
        for myId in [1,2,10,11,12,20,21,22]:
            f=mesh.getMeasureField(True)
            f.changeUnderlyingMesh(mesh2,myId,1e-12);
            pass
        pass

    def testDADFindCommonTuples1(self):
        da=DataArrayDouble.New();
        # nbOftuples=1
        array1=[2.3,1.2,1.3,2.3,2.301,0.8]
        da.setValues(array1,6,1)
        c,cI=da.findCommonTuples(1e-2);
        expected1=[0,3,4]
        expected2=[0,3]
        self.assertEqual(3,c.getNbOfElems());
        self.assertEqual(2,cI.getNbOfElems());
        self.assertEqual(expected1,c.getValues())
        self.assertEqual(expected2,cI.getValues())
        c,cI=da.findCommonTuples(2e-1)
        expected3=[0,3,4,1,2]
        expected4=[0,3,5]
        self.assertEqual(5,c.getNbOfElems());
        self.assertEqual(3,cI.getNbOfElems());
        self.assertEqual(expected3,c.getValues())
        self.assertEqual(expected4,cI.getValues())
        # nbOftuples=2
        array2=[2.3,2.3,1.2,1.2,1.3,1.3,2.3,2.3,2.301,2.301,0.8,0.8]
        da.setValues(array2,6,2)
        c,cI=da.findCommonTuples(1e-2);
        self.assertEqual(3,c.getNbOfElems());
        self.assertEqual(2,cI.getNbOfElems());
        self.assertEqual(expected1,c.getValues())
        self.assertEqual(expected2,cI.getValues())
        c,cI=da.findCommonTuples(2e-1)
        self.assertEqual(5,c.getNbOfElems());
        self.assertEqual(3,cI.getNbOfElems());
        self.assertEqual(expected3,c.getValues())
        self.assertEqual(expected4,cI.getValues())
        # nbOftuples=3
        array3=[2.3,2.3,2.3,1.2,1.2,1.2,1.3,1.3,1.3,2.3,2.3,2.3,2.301,2.301,2.301,0.8,0.8,0.8]
        da.setValues(array3,6,3)
        c,cI=da.findCommonTuples(1e-2);
        self.assertEqual(3,c.getNbOfElems());
        self.assertEqual(2,cI.getNbOfElems());
        self.assertEqual(expected1,c.getValues())
        self.assertEqual(expected2,cI.getValues())
        c,cI=da.findCommonTuples(2e-1)
        self.assertEqual(5,c.getNbOfElems());
        self.assertEqual(3,cI.getNbOfElems());
        self.assertEqual(expected3,c.getValues())
        self.assertEqual(expected4,cI.getValues())
        # nbOftuples=1, no common groups
        array11=[2.3,1.2,1.3,2.4,2.5,0.8]
        da.setValues(array11,6,1)
        c,cI=da.findCommonTuples(1e-2);
        self.assertEqual(0,c.getNbOfElems());
        self.assertEqual(1,cI.getNbOfElems());
        self.assertEqual([0],cI.getValues())
        
        array12=[0.]*(6*5)
        da.setValues(array12,6,5) #bad NumberOfComponents
        self.assertRaises(InterpKernelException, da.findCommonTuples, 1e-2);
        pass

    def testDABack1(self):
        da=DataArrayDouble.New();
        array1=[2.3,1.2,1.3,2.3,2.301,0.8]
        da.setValues(array1,6,1);
        self.assertAlmostEqual(0.8,da.back(),14);
        da.rearrange(2);
        self.assertRaises(InterpKernelException,da.back);
        da.alloc(0,1);
        self.assertRaises(InterpKernelException,da.back);
        #
        da=DataArrayInt.New();
        array2=[4,7,8,2]
        da.setValues(array2,4,1);
        self.assertEqual(2,da.back());
        da.rearrange(2);
        self.assertRaises(InterpKernelException,da.back);
        da.alloc(0,1);
        self.assertRaises(InterpKernelException,da.back);
        pass

    def testDADGetDifferentValues1(self):
        da=DataArrayDouble.New();
        array1=[2.3,1.2,1.3,2.3,2.301,0.8]
        da.setValues(array1,6,1)
        #
        expected1=[2.301,1.2,1.3,0.8]
        dv=da.getDifferentValues(1e-2);
        self.assertEqual(4,dv.getNbOfElems());
        for i in range(4):
            self.assertAlmostEqual(expected1[i],dv.getIJ(i,0),14);
            pass
        #
        dv=da.getDifferentValues(2e-1);
        expected2=[2.301,1.3,0.8]
        self.assertEqual(3,dv.getNbOfElems());
        for i in range(3):
            self.assertAlmostEqual(expected2[i],dv.getIJ(i,0),14);
            pass
        pass

    def testDAIBuildOld2NewArrayFromSurjectiveFormat2(self):
        arr=[0,3, 5,7,9]
        arrI=[0,2,5]
        a=DataArrayInt.New();
        a.setValues(arr,5,1);
        b=DataArrayInt.New();
        b.setValues(arrI,3,1);
        ret,newNbTuple=DataArrayInt.ConvertIndexArrayToO2N(10,a,b);
        expected=[0,1,2,0,3,4,5,4,6,4]
        self.assertEqual(10,ret.getNbOfElems());
        self.assertEqual(7,newNbTuple);
        self.assertEqual(1,ret.getNumberOfComponents());
        self.assertEqual(expected,ret.getValues());
        self.assertRaises(InterpKernelException,DataArrayInt.ConvertIndexArrayToO2N,9,a,b);
        pass

    def testDADIReverse1(self):
        arr=[0,3,5,7,9,2]
        a=DataArrayInt.New();
        a.setValues(arr,6,1);
        self.assertEqual(2,a.back());
        a.reverse();
        for i in range(6):
            self.assertEqual(arr[5-i],a.getIJ(i,0));
            pass
        a.setValues(arr[:-1],5,1);
        a.reverse();
        for i in range(5):
            self.assertEqual(arr[4-i],a.getIJ(i,0));
            pass
        #
        arr2=[0.,3.,5.,7.,9.,2.]
        b=DataArrayDouble.New();
        b.setValues(arr2,6,1);
        b.reverse();
        for i in range(6):
            self.assertAlmostEqual(arr2[5-i],b.getIJ(i,0),14);
            pass
        b.setValues(arr2[:5],5,1);
        self.assertAlmostEqual(9.,b.back(),14)
        b.reverse();
        for i in range(5):
            self.assertAlmostEqual(arr2[4-i],b.getIJ(i,0),14);
            pass
        pass

    def testGetNodeIdsInUse1(self):
        m0=MEDCouplingDataForTest.build2DTargetMesh_1();
        CellIds=[1,2]
        m1=m0.buildPartOfMySelf(CellIds,True);
        arr,newNbOfNodes=m1.getNodeIdsInUse();
        expected=[-1,0,1,-1,2,3,-1,-1,-1]
        self.assertEqual(4,newNbOfNodes);
        self.assertEqual(9,arr.getNbOfElems());
        self.assertEqual(expected,arr.getValues());
        arr2=arr.invertArrayO2N2N2O(newNbOfNodes);
        self.assertEqual(4,arr2.getNbOfElems());
        expected2=[1,2,4,5]
        self.assertEqual(expected2,arr2.getValues());
        pass

    def testBuildDescendingConnec2(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        #
        mesh2,desc,descIndx,revDesc,revDescIndx=mesh.buildDescendingConnectivity2();
        mesh2.checkConsistencyLight();
        self.assertEqual(1,mesh2.getMeshDimension());
        self.assertEqual(13,mesh2.getNumberOfCells());
        self.assertEqual(14,revDescIndx.getNbOfElems()); self.assertEqual(14,revDescIndx.getNumberOfTuples());
        self.assertEqual(6,descIndx.getNbOfElems()); self.assertEqual(6,descIndx.getNumberOfTuples());
        self.assertEqual(18,desc.getNbOfElems()); self.assertEqual(18,desc.getNumberOfTuples());
        self.assertEqual(18,revDesc.getNbOfElems()); self.assertEqual(18,revDesc.getNumberOfTuples());
        expected1=[1,2,3,4,-3,5,6, 7,8,-5,9,10,-2,11, 12,13,-7,-10]
        self.assertEqual(expected1,desc.getValues());
        expected2=[0,4,7,10,14,18]
        self.assertEqual(expected2,descIndx.getValues());
        expected3=[0,1,3,5,6,8,9,11,12,13,15,16,17,18]
        self.assertEqual(expected3,revDescIndx.getValues());
        expected4=[0, 0,3, 0,1, 0, 1,2, 1, 2,4, 2, 3, 3,4, 3, 4, 4]
        self.assertEqual(expected4,revDesc.getValues());
        conn=mesh2.getNodalConnectivity();
        connIndex=mesh2.getNodalConnectivityIndex();
        expected5=[0,3,6,9,12,15,18,21,24,27,30,33,36,39]
        self.assertEqual(expected5,connIndex.getValues());
        expected6=[1, 0, 3, 1, 3, 4, 1, 4, 1, 1, 1, 0, 1, 4, 2, 1, 2, 1, 1, 4, 5, 1, 5, 2, 1, 6, 7, 1, 7, 4, 1, 3, 6, 1, 7, 8, 1, 8, 5]
        self.assertEqual(expected6,conn.getValues());
        pass

    def testIntersect2DMeshesTmp1(self):
        m1c=MEDCouplingCMesh.New();
        coordsX=DataArrayDouble.New();
        arrX=[ -1., 1., 2., 4. ]
        coordsX.setValues(arrX,4,1);
        m1c.setCoordsAt(0,coordsX);
        coordsY=DataArrayDouble.New();
        arrY=[ -2., 2., 4., 8. ]
        coordsY.setValues(arrY,4,1);
        m1c.setCoordsAt(1,coordsY);
        m1=m1c.buildUnstructured()
        m1bis=m1.buildPartOfMySelf([3,4,5],False)
        m2=m1.deepCopy()
        m2=m2.buildPartOfMySelf([0,1,2],False)
        m2.translate([0.5,0.5])
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m1bis,m2,1e-10)
        expected1=[0,0,1,1,1,2,2,2]
        expected2=[0,-1,0,1,-1,1,2,-1]
        self.assertEqual(8,d1.getNumberOfTuples());
        self.assertEqual(8,d2.getNumberOfTuples());
        self.assertEqual(8,m3.getNumberOfCells());
        self.assertEqual(22,m3.getNumberOfNodes());
        self.assertEqual(2,m3.getSpaceDimension());
        self.assertEqual(expected1,d1.getValues());
        self.assertEqual(expected2,d2.getValues());
        expected3=[5,17,1,16,12,5,16,0,4,5,17,12,5,18,1,17,13,5,19,2,18,13,5,17,5,6,19,13,5,20,2,19,14,5,21,3,20,14,5,19,6,7,21,14]
        expected4=[0,5,12,17,22,28,33,38,44]
        expected5=[-1.0,2.0,1.0,2.0,2.0,2.0,4.0,2.0,-1.0,4.0,1.0,4.0,2.0,4.0,4.0,4.0,-0.5,-1.5,1.5,-1.5,2.5,-1.5,4.5,-1.5,-0.5,2.5,1.5,2.5,2.5,2.5,4.5,2.5,-0.5,2.0,1.0,2.5,1.5,2.0,2.0,2.5,2.5,2.0,4.0,2.5]
        self.assertEqual(44,m3.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(9,m3.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m3.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m3.getNodalConnectivityIndex().getValues());
        for i in range(44):
            self.assertAlmostEqual(expected5[i],m3.getCoords().getIJ(0,i),12);
            pass
        pass

    def testFindNodesOnLine1(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        pt=[-0.3,-0.3]
        pt2=[0.,0.,0.]
        pt3=[-0.3,0.,0.]
        vec=[0.,1.]
        vec2=[1.,0.,0.]
        vec3=[0.,1.,1.]
        expected1=[0,3,6]
        res=mesh.findNodesOnLine(pt,vec,1e-12);
        self.assertEqual(3,len(res));
        self.assertEqual(expected1,res.getValues());
        #
        mesh.changeSpaceDimension(3);
        mesh.rotate(pt2,vec2,pi/4.);
        res=mesh.findNodesOnLine(pt3,vec3,1e-12);
        self.assertEqual(3,len(res));
        self.assertEqual(expected1,res.getValues());
        pass

    def testIntersect2DMeshesTmp2(self):
        m1c=MEDCouplingCMesh.New();
        coordsX1=DataArrayDouble.New();
        arrX1=[ 0., 1., 1.5, 2. ]
        coordsX1.setValues(arrX1,4,1);
        m1c.setCoordsAt(0,coordsX1);
        coordsY1=DataArrayDouble.New();
        arrY1=[ 0., 1.5, 3.]
        coordsY1.setValues(arrY1,3,1);
        m1c.setCoordsAt(1,coordsY1);
        m1=m1c.buildUnstructured();
        m2c=MEDCouplingCMesh.New();
        coordsX2=DataArrayDouble.New();
        arrX2=[ 0., 1., 2. ]
        coordsX2.setValues(arrX2,3,1);
        m2c.setCoordsAt(0,coordsX2);
        coordsY2=DataArrayDouble.New();
        arrY2=[ 0., 1., 3.]
        coordsY2.setValues(arrY2,3,1);
        m2c.setCoordsAt(1,coordsY2);
        m2=m2c.buildUnstructured();
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m1,m2,1e-10)
        #
        expected1=[0,0,1,1,2,2,3,4,5]
        expected2=[0,2,1,3,1,3,2,3,3]
        self.assertEqual(9,d1.getNumberOfTuples());
        self.assertEqual(9,d2.getNumberOfTuples());
        self.assertEqual(9,m3.getNumberOfCells());
        self.assertEqual(22,m3.getNumberOfNodes());
        self.assertEqual(2,m3.getSpaceDimension());
        self.assertEqual(expected1,d1.getValues());
        self.assertEqual(expected2,d2.getValues());
        expected3=[5,16,13,12,15,5,15,4,5,16,5,21,2,13,16,5,16,5,6,21,5,17,14,2,21,5,21,6,7,17,5,4,18,19,5,5,5,19,10,6,5,6,10,20,7]
        expected4=[0,5,10,15,20,25,30,35,40,45]
        expected5=[0.0,0.0,1.0,0.0,1.5,0.0,2.0,0.0,0.0,1.5,1.0,1.5,1.5,1.5,2.0,1.5,0.0,3.0,1.0,3.0,1.5,3.0,2.0,3.0,0.0,0.0,1.0,0.0,2.0,0.0,0.0,1.0,1.0,1.0,2.0,1.0,0.0,3.0,1.0,3.0,2.0,3.0,1.5,1.0]
        self.assertEqual(45,m3.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(10,m3.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m3.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m3.getNodalConnectivityIndex().getValues());
        for i in range(44):
            self.assertAlmostEqual(expected5[i],m3.getCoords().getIJ(0,i),12);
            pass
        pass
    
    def testBuildPartOfMySelfSafe1(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1()
        self.assertRaises(InterpKernelException,mesh.buildPartOfMySelf,[0,-1,4,2],True)
        self.assertRaises(InterpKernelException,mesh.buildPartOfMySelf,[0,4,5,4],True)
        pass

    def testIntersect2DMeshesTmp3(self):
        m1Coords=[0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1,0.,-1.5,0.5,0.,1.25,0.,0.70710678118654757,0.70710678118654757,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.70710678118654757,0.70710678118654757,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.70710678118654757,-0.70710678118654757,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.70710678118654757,-0.70710678118654757,1.0606601717798214,-1.0606601717798214];
        m1Conn=[0,3,1,13,11,9, 3,4,2,1,14,12,10,11, 5,3,0,15,13,17, 6,4,3,5,16,14,15,18, 5,0,7,17,21,19, 6,5,7,8,18,19,22,20, 0,1,7,9,23,21, 1,2,8,7,10,24,22,23];
        m1=MEDCouplingUMesh.New();
        m1.setMeshDimension(2);
        m1.allocateCells(8);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[0:6]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[6:14]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[14:20]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[20:28]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[28:34]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[34:42]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[42:48]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[48:56]);
        m1.finishInsertingCells();
        myCoords1=DataArrayDouble.New();
        myCoords1.setValues(m1Coords,25,2);
        m1.setCoords(myCoords1);
        #
        m2Coords=[0.,0.,1.1,0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1,-1.1,-1.,0.,-1.,1.1,-1,1.7,-1.]
        m2Conn=[0,3,2,1, 1,2,5,4, 7,6,3,0, 8,9,6,7, 7,0,12,11, 8,7,11,10, 0,1,13,12, 1,4,14,13]
        m2=MEDCouplingUMesh.New();
        m2.setMeshDimension(2);
        m2.allocateCells(8);
        for i in range(8):
            m2.insertNextCell(NORM_QUAD4,4,m2Conn[4*i:4*(i+1)])
            pass
        m2.finishInsertingCells();
        myCoords2=DataArrayDouble.New();
        myCoords2.setValues(m2Coords,15,2);
        m2.setCoords(myCoords2);
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m1,m2,1e-10)
        m3.unPolyze()
        #
        expected1=[0,1,1,1,2,3,3,3,4,5,5,5,6,7,7,7]
        expected2=[0,0,1,-1,2,2,3,-1,4,4,5,-1,6,6,7,-1]
        self.assertEqual(16,d1.getNumberOfTuples());
        self.assertEqual(16,d2.getNumberOfTuples());
        self.assertEqual(16,m3.getNumberOfCells());
        self.assertEqual(104,m3.getNumberOfNodes());
        self.assertEqual(2,m3.getSpaceDimension());
        self.assertEqual(expected1,d1.getValues());
        self.assertEqual(expected2,d2.getValues());
        expected3=[6,28,1,25,44,45,46,8,26,1,28,27,47,48,49,50,8,40,2,26,27,51,52,53,54,8,28,4,40,27,55,56,57,58,6,28,25,5,59,60,61,8,28,5,32,31,62,63,64,65,8,32,6,41,31,66,67,68,69,8,41,4,28,31,70,71,72,73,6,25,37,5,74,75,76,8,32,5,37,36,77,78,79,80,8,42,6,32,36,81,82,83,84,8,37,8,42,36,85,86,87,88,6,1,37,25,89,90,91,8,37,1,26,38,92,93,94,95,8,26,2,43,38,96,97,98,99,8,43,8,37,38,100,101,102,103]
        expected4=[0,7,16,25,34,41,50,59,68,75,84,93,102,109,118,127,136]
        expected5=[0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1.,0.,-1.5,0.5,0.,1.25,0.,0.7071067811865476,0.7071067811865476,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.7071067811865476,0.7071067811865476,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.7071067811865476,-0.7071067811865476,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.7071067811865476,-0.7071067811865476,1.0606601717798214,-1.0606601717798214,0.,0.,1.1,0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1.,-1.1,-1.,0.,-1.,1.1,-1.,1.7,-1.,1.118033988749895,1.,-1.118033988749895,1.,-1.118033988749895,-1.,1.118033988749895,-1.,0.7071067811865477,0.7071067811865476,0.5,0.,0.,0.5,1.05,0.,0.7071067811865475,0.7071067811865477,0.55,1.,1.1,0.5,1.4012585384440737,0.535233134659635,1.3,0.,1.1,0.5,1.1090169943749475,1.,0.,1.25,0.6123724356957946,1.369306393762915,1.1090169943749475,1.,0.55,1.,0.,0.5,-0.5,0.,-0.7071067811865477,0.7071067811865476,-0.7071067811865475,0.7071067811865477,-1.05,0.,-1.1,0.5,-0.55,1.,-1.3,0.,-1.4012585384440737,0.5352331346596344,-1.1090169943749475,1.,-1.1,0.5,-0.6123724356957941,1.3693063937629155,0.,1.25,-0.55,1.,-1.1090169943749475,1.,0.,-0.5,-0.7071067811865475,-0.7071067811865477,-0.5,0.,-1.05,0.,-0.7071067811865478,-0.7071067811865475,-0.55,-1.,-1.1,-0.5,-1.4012585384440734,-0.5352331346596354,-1.3,0.,-1.1,-0.5,-1.1090169943749475,-1.,0.,-1.25,-0.6123724356957945,-1.369306393762915,-1.1090169943749475,-1.,-0.55,-1.,0.7071067811865475,-0.7071067811865477,0.,-0.5,0.5,0.,0.7071067811865477,-0.7071067811865475,1.05,0.,1.1,-0.5,0.55,-1.,1.3,0.,1.4012585384440737,-0.535233134659635,1.1090169943749475,-1.,1.1,-0.5,0.6123724356957946,-1.369306393762915,0.,-1.25,0.55,-1.,1.1090169943749475,-1.0]
        self.assertEqual(136,m3.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(17,m3.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m3.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m3.getNodalConnectivityIndex().getValues());
        for i in range(208):
            self.assertAlmostEqual(expected5[i],m3.getCoords().getIJ(0,i),12);
            pass
        pass

    def testUMeshTessellate2D1(self):
        m1Coords=[0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1,0.,-1.5,0.5,0.,1.25,0.,0.70710678118654757,0.70710678118654757,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.70710678118654757,0.70710678118654757,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.70710678118654757,-0.70710678118654757,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.70710678118654757,-0.70710678118654757,1.0606601717798214,-1.0606601717798214];
        m1Conn=[0,3,1,13,11,9, 3,4,2,1,14,12,10,11, 5,3,0,15,13,17, 6,4,3,5,16,14,15,18, 5,0,7,17,21,19, 6,5,7,8,18,19,22,20, 0,1,7,9,23,21, 1,2,8,7,10,24,22,23];
        m1=MEDCouplingUMesh.New();
        m1.setMeshDimension(2);
        m1.allocateCells(8);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[0:6]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[6:14]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[14:20]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[20:28]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[28:34]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[34:42]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[42:48]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[48:56]);
        m1.finishInsertingCells();
        myCoords1=DataArrayDouble.New();
        myCoords1.setValues(m1Coords,25,2);
        m1.setCoords(myCoords1);
        #
        m11=m1.deepCopy();
        m11.tessellate2D(1.);
        self.assertTrue(m11.getCoords().isEqual(m11.getCoords(),1e-12));
        expected1=[5,0,3,11,1,5,3,4,12,2,1,11,5,5,15,3,0,5,6,16,4,3,15,5,5,5,0,7,19,5,6,5,19,7,8,20,5,0,1,23,7,5,1,2,24,8,7,23]
        expected2=[0,5,12,17,24,29,36,41,48]
        self.assertEqual(48,m11.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(9,m11.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected1,m11.getNodalConnectivity().getValues());
        self.assertEqual(expected2,m11.getNodalConnectivityIndex().getValues());
        #
        m12=m1.deepCopy();
        m12.tessellate2D(0.5);
        self.assertEqual(41,m12.getNumberOfNodes());
        expected3=[5,0,3,25,26,1,5,3,4,27,28,2,1,26,25,5,5,29,30,3,0,5,6,31,32,4,3,30,29,5,5,5,0,7,33,34,5,6,5,34,33,7,8,35,36,5,0,1,37,38,7,5,1,2,39,40,8,7,38,37]
        expected4=[0,6,15,21,30,36,45,51,60]
        expected5=[0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1.,0.,-1.5,0.5,0.,1.25,0.,0.7071067811865476,0.7071067811865476,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.7071067811865476,0.7071067811865476,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.7071067811865476,-0.7071067811865476,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.7071067811865476,-0.7071067811865476,1.0606601717798214,-1.0606601717798214,0.479425538604203,0.8775825618903728,0.8414709848078964,0.54030230586814,0.7191383079063044,1.3163738428355591,1.2622064772118446,0.8104534588022099,-0.877582561890373,0.4794255386042027,-0.5403023058681399,0.8414709848078964,-1.3163738428355596,0.7191383079063038,-0.8104534588022098,1.2622064772118446,-0.4794255386042031,-0.8775825618903728,-0.8414709848078965,-0.5403023058681399,-0.7191383079063045,-1.3163738428355591,-1.2622064772118449,-0.8104534588022098,0.8775825618903729,-0.47942553860420295,0.54030230586814,-0.8414709848078964,1.3163738428355594,-0.7191383079063043,0.8104534588022099,-1.2622064772118446]
        for i in range(82):
            self.assertAlmostEqual(expected5[i],m12.getCoords().getIJ(0,i),12);
            pass
        self.assertEqual(60,m12.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(9,m12.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m12.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m12.getNodalConnectivityIndex().getValues());
        pass

    def testUMeshTessellate2DCurve1(self):
        # A quarter of circle:
        mcoords = [0.4,0.0,   0.0,-0.4,   0.283,-0.283]
        mconnec = [0,1,2]

        m1 = MEDCouplingUMesh.New()
        m1.setMeshDimension(1)
        m1.allocateCells(1)
        m1.insertNextCell(NORM_SEG3, mconnec)

        myCoords = DataArrayDouble.New(mcoords, 3, 2)
        m1.setCoords(myCoords)
        
        m2 = m1.deepCopy()
        m2.tessellate2D(0.1)
        # If the following raises, the test will fail automatically:
        m2.checkConsistency(0.0) # eps param not used

    def testIntersect2DMeshesTmp4(self):
        m1Coords=[0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1,0.,-1.5,0.5,0.,1.25,0.,0.70710678118654757,0.70710678118654757,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.70710678118654757,0.70710678118654757,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.70710678118654757,-0.70710678118654757,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.70710678118654757,-0.70710678118654757,1.0606601717798214,-1.0606601717798214];
        m1Conn=[0,3,1,13,11,9, 3,4,2,1,14,12,10,11, 5,3,0,15,13,17, 6,4,3,5,16,14,15,18, 5,0,7,17,21,19, 6,5,7,8,18,19,22,20, 0,1,7,9,23,21, 1,2,8,7,10,24,22,23];
        m1=MEDCouplingUMesh.New();
        m1.setMeshDimension(2);
        m1.allocateCells(8);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[0:6]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[6:14]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[14:20]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[20:28]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[28:34]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[34:42]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[42:48]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[48:56]);
        m1.finishInsertingCells();
        myCoords1=DataArrayDouble.New();
        myCoords1.setValues(m1Coords,25,2);
        m1.setCoords(myCoords1);
        #
        m2Coords=[0.,0.,1.1,0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1,-1.1,-1.,0.,-1.,1.1,-1,1.7,-1.]
        m2Conn=[0,3,2,1, 1,2,5,4, 7,6,3,0, 8,9,6,7, 7,0,12,11, 8,7,11,10, 0,1,13,12, 1,4,14,13]
        m2=MEDCouplingUMesh.New();
        m2.setMeshDimension(2);
        m2.allocateCells(8);
        for i in range(8):
            m2.insertNextCell(NORM_QUAD4,4,m2Conn[4*i:4*(i+1)])
            pass
        m2.finishInsertingCells();
        myCoords2=DataArrayDouble.New();
        myCoords2.setValues(m2Coords,15,2);
        m2.setCoords(myCoords2);
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m2,m1,1e-10)
        m3.unPolyze()
        #
        expected1=[0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7]
        expected2=[0,1,1,-1,2,3,3,-1,4,5,5,-1,6,7,7,-1]
        self.assertEqual(16,d1.getNumberOfTuples());
        self.assertEqual(16,d2.getNumberOfTuples());
        self.assertEqual(16,m3.getNumberOfCells());
        self.assertEqual(104,m3.getNumberOfNodes());
        self.assertEqual(2,m3.getSpaceDimension());
        self.assertEqual(expected1,d1.getValues());
        self.assertEqual(expected2,d2.getValues());
        expected3=[6,16,15,18,44,45,46,8,18,2,1,16,47,48,49,50,8,17,1,2,40,51,52,53,54,8,40,5,4,17,55,56,57,58,6,18,15,20,59,60,61,8,20,7,6,18,62,63,64,65,8,41,6,7,21,66,67,68,69,8,21,8,9,41,70,71,72,73,6,20,15,22,74,75,76,8,22,11,7,20,77,78,79,80,8,21,7,11,42,81,82,83,84,8,42,10,8,21,85,86,87,88,6,22,15,16,89,90,91,8,16,1,13,22,92,93,94,95,8,43,13,1,17,96,97,98,99,8,17,4,14,43,100,101,102,103]
        expected4=[0,7,16,25,34,41,50,59,68,75,84,93,102,109,118,127,136]
        expected5=[0.,0.,1.1, 0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1.,-1.1,-1.,0.,-1.,1.1,-1.,1.7,-1.,0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1.,0.,-1.5,0.5,0.,1.25,0.,0.7071067811865476,0.7071067811865476,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.7071067811865476,0.7071067811865476,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.7071067811865476,-0.7071067811865476,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.7071067811865476,-0.7071067811865476,1.0606601717798214,-1.0606601717798214,1.1180339887498951,1.,-1.1180339887498951,1.,-1.1180339887498951,-1.,1.1180339887498951,-1.,0.5,0.,0.,0.5,0.7071067811865477,0.7071067811865476,0.55,1.,1.1,0.5,1.05,0.,0.7071067811865477,0.7071067811865475,1.3,0.,1.1,0.5,1.1090169943749475,1.,1.4012585384440737,0.535233134659635,1.4090169943749475,1.,1.7,0.5,1.6,0.,1.4012585384440737,0.535233134659635,0.,0.5,-0.5,0.,-0.7071067811865477,0.7071067811865476,-1.05,0.,-1.1,0.5,-0.55,1.,-0.7071067811865478,0.7071067811865475,-1.1090169943749475,1.,-1.1,0.5,-1.3,0.,-1.4012585384440737,0.5352331346596344,-1.6,0.,-1.7,0.5,-1.4090169943749475,1.,-1.4012585384440737,0.5352331346596344,-0.5,0.,0.,-0.5,-0.7071067811865475,-0.7071067811865477,-0.55,-1.,-1.1,-0.5,-1.05,0.,-0.7071067811865475,-0.7071067811865477,-1.3,0.,-1.1,-0.5,-1.1090169943749475,-1.,-1.4012585384440734,-0.5352331346596354,-1.4090169943749475,-1.,-1.7,-0.5,-1.6,0.,-1.4012585384440732,-0.5352331346596354,0.,-0.5,0.5,0.,0.7071067811865475,-0.7071067811865477,1.05,0.,1.1,-0.5,0.55,-1.,0.7071067811865475,-0.7071067811865477,1.1090169943749475,-1.,1.1,-0.5,1.3,0.,1.4012585384440737,-0.535233134659635,1.6,0.,1.7,-0.5,1.4090169943749475,-1.,1.4012585384440737,-0.535233134659635]
        self.assertEqual(136,m3.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(17,m3.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m3.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m3.getNodalConnectivityIndex().getValues());
        for i in range(208):
            self.assertAlmostEqual(expected5[i],m3.getCoords().getIJ(0,i),12);
            pass
        pass

    def testGetCellIdsCrossingPlane1(self):
        mesh3D,mesh2D=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        vec=[-0.07,1.,0.07]
        origin=[1.524,1.4552,1.74768]
        ids1=mesh3D.getCellIdsCrossingPlane(origin,vec,1e-10)
        self.assertEqual([1,3,4,7,9,10,13,15,16],ids1.getValues())
        vec2=[0.,0.,1.]
        ids2=mesh3D.getCellIdsCrossingPlane(origin,vec2,1e-10)
        self.assertEqual([6,7,8,9,10,11],ids2.getValues())
        pass

    def testBuildSlice3D1(self):
        mesh3D,mesh2D=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        vec1=[-0.07,1.,0.07]
        origin1=[1.524,1.4552,1.74768]
        slice1,ids=mesh3D.buildSlice3D(origin1,vec1,1e-10);
        expected1=[1,3,4,7,9,10,13,15,16]
        expected2=[5,42,41,40,43,44,5,42,46,45,41,5,44,43,40,47,48,5,49,42,44,50,5,49,51,46,42,5,50,44,48,52,5,53,49,50,54,5,53,55,51,49,5,54,50,52,56]
        expected3=[0,6,11,17,22,27,32,37,42,47]
        expected4=[1.,1.,0.,1.,1.25,0.,1.,1.5,0.,2.,1.,0.,1.,2.,0.,0.,2.,0.,3.,1.,0.,3.,2.,0.,0.,1.,0.,2.,2.,0.,1.,1.,1.,1.,1.25,1.,1.,1.5,1.,2.,1.,1.,1.,2.,1.,0.,2.,1.,3.,1.,1.,3.,2.,1.,0.,1.,1.,2.,2.,1.,1.,1.,2.,1.,1.25,2.,1.,1.5,2.,2.,1.,2.,1.,2.,2.,0.,2.,2.,3.,1.,2.,3.,2.,2.,0.,1.,2.,2.,2.,2.,1.,1.,3.,1.,1.25,3.,1.,1.5,3.,2.,1.,3.,1.,2.,3.,0.,2.,3.,3.,1.,3.,3.,2.,3.,0.,1.,3.,2.,2.,3.,1.,1.5408576,0.,2.,1.6108576000000001,0.,2.,1.5408576,1.,1.,1.5,0.5836800000000008,1.,1.4708576,1.,3.,1.6808576,0.,3.,1.6108576000000001,1.,0.,1.4708576,0.,0.,1.4008576,1.,2.,1.4708576,2.,1.,1.4008576000000001,2.,3.,1.5408575999999998,2.,0.,1.3308575999999999,2.,2.,1.4008576,3.,1.,1.3308576,3.,3.,1.4708576,3.,0.,1.2608576,3.]
        self.assertEqual(2,slice1.getMeshDimension());
        self.assertEqual(3,slice1.getSpaceDimension());
        self.assertEqual(57,slice1.getNumberOfNodes());
        self.assertEqual(9,slice1.getNumberOfCells());
        self.assertEqual(9,ids.getNumberOfTuples());
        self.assertEqual(47,slice1.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(10,slice1.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected1,ids.getValues());
        self.assertEqual(expected2,slice1.getNodalConnectivity().getValues());
        self.assertEqual(expected3,slice1.getNodalConnectivityIndex().getValues());
        for i in range(171):
            self.assertAlmostEqual(expected4[i],slice1.getCoords().getIJ(0,i),12);
            pass
        # 2nd slice based on already existing nodes of mesh3D.
        vec2=[0.,3.,1.]
        origin2=[2.5,1.,3.]
        slice1,ids=mesh3D.buildSlice3D(origin2,vec2,1e-10);
        expected5=[5,50,10,4,51,5,50,52,7,10,5,51,4,5,53,5,54,50,51,55,56,5,54,57,52,50,5,56,55,51,53,58,5,38,59,56,54,43,5,54,57,46,43,5,38,59,56,58,48]
        expected6=[0,5,10,15,21,26,32,38,43,49]
        expected7=[1.,1.,0.,1.,1.25,0.,1.,1.5,0.,2.,1.,0.,1.,2.,0.,0.,2.,0.,3.,1.,0.,3.,2.,0.,0.,1.,0.,1.,3.,0.,2.,2.,0.,2.,3.,0.,1.,1.,1.,1.,1.25,1.,1.,1.5,1.,2.,1.,1.,1.,2.,1.,0.,2.,1.,3.,1.,1.,3.,2.,1.,0.,1.,1.,1.,3.,1.,2.,2.,1.,2.,3.,1.,0.,0.,2.,1.,1.,2.,1.,1.25,2.,1.,0.,2.,1.,1.5,2.,2.,0.,2.,2.,1.,2.,1.,2.,2.,0.,2.,2.,3.,1.,2.,3.,2.,2.,0.,1.,2.,2.,2.,2.,0.,0.,3.,1.,1.,3.,1.,1.25,3.,1.,0.,3.,1.,1.5,3.,2.,0.,3.,2.,1.,3.,1.,2.,3.,0.,2.,3.,3.,1.,3.,3.,2.,3.,0.,1.,3.,2.,2.,3.,2.,1.6666666666666667,1.,1.,1.6666666666666667,1.,3.,1.6666666666666667,1.,0.,1.6666666666666667,1.,2.,1.3333333333333335,2.,1.,1.5,1.5,1.,1.3333333333333333,2.,3.,1.3333333333333335,2.,0.,1.3333333333333335,2.,1.,1.25,2.25]
        self.assertEqual(2,slice1.getMeshDimension());
        self.assertEqual(3,slice1.getSpaceDimension());
        self.assertEqual(60,slice1.getNumberOfNodes());
        self.assertEqual(9,slice1.getNumberOfCells());
        self.assertEqual(9,ids.getNumberOfTuples());
        self.assertEqual(49,slice1.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(10,slice1.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected1,ids.getValues());
        self.assertEqual(expected5,slice1.getNodalConnectivity().getValues());
        self.assertEqual(expected6,slice1.getNodalConnectivityIndex().getValues());
        for i in range(180):
            self.assertAlmostEqual(expected7[i],slice1.getCoords().getIJ(0,i),12);
            pass
        # 3rd slice based on shared face of mesh3D.
        vec3=[0.,0.,1.]
        origin3=[2.5,1.,2.]
        slice1,ids=mesh3D.buildSlice3D(origin3,vec3,1e-10);
        expected8=[6,7,8,9,10,11,12,13,14,15,16,17]
        expected9=[5,15,26,16,18,5,16,21,28,22,19,17,5,18,20,21,16,5,21,24,25,28,5,26,16,17,19,22,23,5,22,27,29,28,5,15,26,16,18,5,16,21,28,22,19,17,5,18,20,21,16,5,21,24,25,28,5,26,16,17,19,22,23,5,22,27,29,28]
        expected10=[0,5,12,17,22,29,34,39,46,51,56,63,68]
        expected11=[0.,0.,1.,1.,1.,1.,1.,1.25,1.,1.,0.,1.,1.,1.5,1.,2.,0.,1.,2.,1.,1.,1.,2.,1.,0.,2.,1.,3.,1.,1.,3.,2.,1.,0.,1.,1.,1.,3.,1.,2.,2.,1.,2.,3.,1.,0.,0.,2.,1.,1.,2.,1.,1.25,2.,1.,0.,2.,1.,1.5,2.,2.,0.,2.,2.,1.,2.,1.,2.,2.,0.,2.,2.,3.,1.,2.,3.,2.,2.,0.,1.,2.,1.,3.,2.,2.,2.,2.,2.,3.,2.,0.,0.,3.,1.,1.,3.,1.,1.25,3.,1.,0.,3.,1.,1.5,3.,2.,0.,3.,2.,1.,3.,1.,2.,3.,0.,2.,3.,3.,1.,3.,3.,2.,3.,0.,1.,3.,1.,3.,3.,2.,2.,3.,2.,3.,3.]
        self.assertEqual(2,slice1.getMeshDimension());
        self.assertEqual(3,slice1.getSpaceDimension());
        self.assertEqual(45,slice1.getNumberOfNodes());
        self.assertEqual(12,slice1.getNumberOfCells());
        self.assertEqual(12,ids.getNumberOfTuples());
        self.assertEqual(68,slice1.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(13,slice1.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected8,ids.getValues());
        self.assertEqual(expected9,slice1.getNodalConnectivity().getValues());
        self.assertEqual(expected10,slice1.getNodalConnectivityIndex().getValues());
        for i in range(135):
            self.assertAlmostEqual(expected11[i],slice1.getCoords().getIJ(0,i),12);
            pass
        pass

    def testBuildSlice3DSurf1(self):
        mesh3D,mesh2D=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        mesh2D=mesh3D.buildDescendingConnectivity()[0];
        vec1=[-0.07,1.,0.07]
        origin1=[1.524,1.4552,1.74768]
        slice1,ids=mesh2D.buildSlice3DSurf(origin1,vec1,1e-10);
        expected1=[6,8,10,11,13,18,19,21,23,25,26,38,41,43,47,49,52,53,64,67,69,73,75,78,79]
        expected2=[1,40,41,1,42,41,1,40,43,1,44,43,1,42,44,1,45,41,1,42,46,1,46,45,1,47,40,1,47,48,1,44,48,1,49,42,1,44,50,1,49,50,1,49,51,1,51,46,1,48,52,1,50,52,1,53,49,1,50,54,1,53,54,1,53,55,1,55,51,1,52,56,1,54,56]
        expected3=[0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75];
        expected4=[1.,1.,0.,1.,1.25,0.,1.,1.5,0.,2.,1.,0.,1.,2.,0.,0.,2.,0.,3.,1.,0.,3.,2.,0.,0.,1.,0.,2.,2.,0.,1.,1.,1.,1.,1.25,1.,1.,1.5,1.,2.,1.,1.,1.,2.,1.,0.,2.,1.,3.,1.,1.,3.,2.,1.,0.,1.,1.,2.,2.,1.,1.,1.,2.,1.,1.25,2.,1.,1.5,2.,2.,1.,2.,1.,2.,2.,0.,2.,2.,3.,1.,2.,3.,2.,2.,0.,1.,2.,2.,2.,2.,1.,1.,3.,1.,1.25,3.,1.,1.5,3.,2.,1.,3.,1.,2.,3.,0.,2.,3.,3.,1.,3.,3.,2.,3.,0.,1.,3.,2.,2.,3.,1.,1.5408576,0.,2.,1.6108576000000001,0.,2.,1.5408576,1.,1.,1.5,0.5836800000000008,1.,1.4708576,1.,3.,1.6808576,0.,3.,1.6108576000000001,1.,0.,1.4708576,0.,0.,1.4008576,1.,2.,1.4708576,2.,1.,1.4008576000000001,2.,3.,1.5408575999999998,2.,0.,1.3308575999999999,2.,2.,1.4008576,3.,1.,1.3308576,3.,3.,1.4708576,3.,0.,1.2608576,3.]
        self.assertEqual(1,slice1.getMeshDimension());
        self.assertEqual(3,slice1.getSpaceDimension());
        self.assertEqual(57,slice1.getNumberOfNodes());
        self.assertEqual(25,slice1.getNumberOfCells());
        self.assertEqual(25,ids.getNumberOfTuples());
        self.assertEqual(75,slice1.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(26,slice1.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected1,ids.getValues());
        self.assertEqual(expected2,slice1.getNodalConnectivity().getValues());
        self.assertEqual(expected3,slice1.getNodalConnectivityIndex().getValues());
        for i in range(171):
            self.assertAlmostEqual(expected4[i],slice1.getCoords().getIJ(0,i),12);
            pass
        #
        vec2=[0.,0.,1.]
        origin2=[2.5,1.,2.]
        slice1,ids=mesh2D.buildSlice3DSurf(origin2,vec2,1e-10);
        expected5=[32,32,32,32,33,34,35,36,37,38,39,40,41,42,43,43,43,43,43,43,44,44,44,44,45,46,47,47,47,47,48,49,50,51,52,53,53,53,53,53,53,54,54,54,54,55,56,57,59,60,61,62,63,64,65,66,67,68,71,72,74,75,76,77,78,81,82,83]
        expected6=[1,15,18,1,18,16,1,16,26,1,26,15,1,26,15,1,16,26,1,18,16,1,15,18,1,16,21,1,21,28,1,22,28,1,19,22,1,17,19,1,16,17,1,16,21,1,21,28,1,28,22,1,22,19,1,19,17,1,17,16,1,16,18,1,18,20,1,20,21,1,21,16,1,20,21,1,18,20,1,28,21,1,21,24,1,24,25,1,25,28,1,25,28,1,24,25,1,21,24,1,23,22,1,26,23,1,26,16,1,16,17,1,17,19,1,19,22,1,22,23,1,23,26,1,22,28,1,28,29,1,29,27,1,27,22,1,27,22,1,29,27,1,28,29,1,26,15,1,16,26,1,18,16,1,15,18,1,16,21,1,21,28,1,22,28,1,19,22,1,17,19,1,16,17,1,20,21,1,18,20,1,25,28,1,24,25,1,21,24,1,23,22,1,26,23,1,27,22,1,29,27,1,28,29]
        expected7=[0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102,105,108,111,114,117,120,123,126,129,132,135,138,141,144,147,150,153,156,159,162,165,168,171,174,177,180,183,186,189,192,195,198,201,204];
        expected8=[0.,0.,1.,1.,1.,1.,1.,1.25, 1.,1.,0.,1.,1.,1.5, 1.,2.,0.,1.,2.,1.,1.,1.,2.,1.,0.,2.,1.,3.,1.,1.,3.,2.,1.,0.,1.,1.,1.,3.,1.,2.,2.,1.,2.,3.,1.,0.,0.,2.,1.,1.,2.,1.,1.25, 2.,1.,0.,2.,1.,1.5, 2.,2.,0.,2.,2.,1.,2.,1.,2.,2.,0.,2.,2.,3.,1.,2.,3.,2.,2.,0.,1.,2.,1.,3.,2.,2.,2.,2.,2.,3.,2.,0.,0.,3.,1.,1.,3.,1.,1.25, 3.,1.,0.,3.,1.,1.5, 3.,2.,0.,3.,2.,1.,3.,1.,2.,3.,0.,2.,3.,3.,1.,3.,3.,2.,3.,0.,1.,3.,1.,3.,3.,2.,2.,3.,2.,3.,3.]
        self.assertEqual(1,slice1.getMeshDimension());
        self.assertEqual(3,slice1.getSpaceDimension());
        self.assertEqual(45,slice1.getNumberOfNodes());
        self.assertEqual(68,slice1.getNumberOfCells());
        self.assertEqual(68,ids.getNumberOfTuples());
        self.assertEqual(204,slice1.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(69,slice1.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected5,ids.getValues());
        self.assertEqual(expected6,slice1.getNodalConnectivity().getValues());
        self.assertEqual(expected7,slice1.getNodalConnectivityIndex().getValues());
        for i in range(135):
            self.assertAlmostEqual(expected8[i],slice1.getCoords().getIJ(0,i),12);
            pass
        pass

    def testDataArrayDoubleAdvSetting1(self):
        data1=[1.,11.,2.,12.,3.,13.,4.,14.,5.,15.,6.,16.,7.,17.]
        data2=[8.,38.,9.,39.,0.,30.,11.,41.,12.,42.]
        compsCpp=["comp1","comp2"]
        da=DataArrayDouble.New();
        da.setInfoAndChangeNbOfCompo(compsCpp);
        da.setName("da");
        da.alloc(7,2);
        compsCpp=compsCpp[:-1]
        self.assertRaises(InterpKernelException,da.setInfoAndChangeNbOfCompo,compsCpp);
        da.setValues(data1,7,2)
        #
        p=[(0,3),(3,5),(5,7)]
        tmp=da.selectByTupleRanges(p);
        self.assertTrue(tmp.isEqual(da,1e-14));
        p=[(0,2),(3,4),(5,7)]
        tmp=da.selectByTupleRanges(p);
        expected1=[1.,11.,2.,12.,4.,14.,6.,16.,7.,17.]
        self.assertEqual(5,tmp.getNumberOfTuples());
        self.assertEqual(2,tmp.getNumberOfComponents());
        for i in range(10):
            self.assertAlmostEqual(expected1[i],tmp.getIJ(0,i),14);
            pass
        p=[(0,2),(0,2),(5,6)]
        tmp=da.selectByTupleRanges(p);
        expected2=[1.,11.,2.,12.,1.,11.,2.,12.,6.,16.]
        self.assertEqual(5,tmp.getNumberOfTuples());
        self.assertEqual(2,tmp.getNumberOfComponents());
        for i in range(10):
            self.assertAlmostEqual(expected2[i],tmp.getIJ(0,i),14);
            pass
        p=[(0,2),(-1,2),(5,6)]
        self.assertRaises(InterpKernelException,da.selectByTupleRanges,p);
        p=[(0,2),(0,2),(5,8)]
        self.assertRaises(InterpKernelException,da.selectByTupleRanges,p);
        #
        da2=DataArrayDouble.New();
        da2.setValues(data2,5,2);
        #
        dac=da.deepCopy();
        dac.setContigPartOfSelectedValuesSlice(1,da2,2,4,1);
        expected3=[1.,11.,0.,30.,11.,41.,4.,14.,5.,15.,6.,16.,7.,17.]
        for i in range(14):
            self.assertAlmostEqual(expected3[i],dac.getIJ(0,i),14);
            pass
        #
        dac=da.deepCopy();
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValuesSlice,3,da2,0,5,1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValuesSlice,0,da2,4,6,1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValuesSlice,3,da2,5,0,1);
        dac.setContigPartOfSelectedValuesSlice(3,da2,1,5,1);
        expected4=[1.,11.,2.,12.,3.,13.,9.,39.,0.,30.,11.,41.,12.,42.]
        for i in range(14):
            self.assertAlmostEqual(expected4[i],dac.getIJ(0,i),14);
            pass
        #
        ids=DataArrayInt.New();
        ids.alloc(3,1);
        dac=da.deepCopy();
        ids.setIJ(0,0,2); ids.setIJ(1,0,0); ids.setIJ(2,0,4);
        dac.setContigPartOfSelectedValues(2,da2,ids);
        expected5=[1.,11.,2.,12.,0.,30.,8.,38.,12.,42.,6.,16.,7.,17.]
        for i in range(14):
            self.assertAlmostEqual(expected5[i],dac.getIJ(0,i),14);
            pass
        #
        dac=da.deepCopy();
        ids.setIJ(0,0,2); ids.setIJ(1,0,5); ids.setIJ(2,0,4);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues,1,da2,ids);
        ids.setIJ(0,0,2); ids.setIJ(1,0,2); ids.setIJ(2,0,-1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues,1,da2,ids);
        ids.setIJ(0,0,2); ids.setIJ(1,0,2); ids.setIJ(2,0,1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues,5,da2,ids);
        #
        ids.setIJ(0,0,2); ids.setIJ(1,0,2); ids.setIJ(2,0,1);
        dac=da.deepCopy();
        dac.setContigPartOfSelectedValues(4,da2,ids);
        expected6=[1.,11.,2.,12.,3.,13.,4.,14.,0.,30.,0.,30.,9.,39.]
        for i in range(14):
            self.assertAlmostEqual(expected6[i],dac.getIJ(0,i),14);
            pass
        pass

    def testDataArrayIntAdvSetting1(self):
        data1=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        data2=[8,38,9,39,0,30,11,41,12,42]
        compsCpp=["comp1","comp2"]
        da=DataArrayInt.New();
        da.setInfoAndChangeNbOfCompo(compsCpp);
        da.setName("da");
        da.alloc(7,2);
        compsCpp=compsCpp[:-1]
        self.assertRaises(InterpKernelException,da.setInfoAndChangeNbOfCompo,compsCpp);
        da.setValues(data1,7,2)
        #
        p=[(0,3),(3,5),(5,7)]
        tmp=da.selectByTupleRanges(p);
        self.assertTrue(tmp.isEqual(da));
        p=[(0,2),(3,4),(5,7)]
        tmp=da.selectByTupleRanges(p);
        expected1=[1,11,2,12,4,14,6,16,7,17]
        self.assertEqual(5,tmp.getNumberOfTuples());
        self.assertEqual(2,tmp.getNumberOfComponents());
        for i in range(10):
            self.assertEqual(expected1[i],tmp.getIJ(0,i));
            pass
        p=[(0,2),(0,2),(5,6)]
        tmp=da.selectByTupleRanges(p);
        expected2=[1,11,2,12,1,11,2,12,6,16]
        self.assertEqual(5,tmp.getNumberOfTuples());
        self.assertEqual(2,tmp.getNumberOfComponents());
        for i in range(10):
            self.assertEqual(expected2[i],tmp.getIJ(0,i));
            pass
        p=[(0,2),(-1,2),(5,6)]
        self.assertRaises(InterpKernelException,da.selectByTupleRanges,p);
        p=[(0,2),(0,2),(5,8)]
        self.assertRaises(InterpKernelException,da.selectByTupleRanges,p);
        #
        da2=DataArrayInt.New();
        da2.setValues(data2,5,2);
        #
        dac=da.deepCopy();
        dac.setContigPartOfSelectedValuesSlice(1,da2,2,4,1);
        expected3=[1,11,0,30,11,41,4,14,5,15,6,16,7,17]
        for i in range(14):
            self.assertEqual(expected3[i],dac.getIJ(0,i));
            pass
        #
        dac=da.deepCopy();
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValuesSlice,3,da2,0,5,1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValuesSlice,0,da2,4,6,1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValuesSlice,3,da2,5,0,1);
        dac.setContigPartOfSelectedValuesSlice(3,da2,1,5,1);
        expected4=[1,11,2,12,3,13,9,39,0,30,11,41,12,42]
        for i in range(14):
            self.assertEqual(expected4[i],dac.getIJ(0,i));
            pass
        #
        ids=DataArrayInt.New();
        ids.alloc(3,1);
        dac=da.deepCopy();
        ids.setIJ(0,0,2); ids.setIJ(1,0,0); ids.setIJ(2,0,4);
        dac.setContigPartOfSelectedValues(2,da2,ids);
        expected5=[1,11,2,12,0,30,8,38,12,42,6,16,7,17]
        for i in range(14):
            self.assertEqual(expected5[i],dac.getIJ(0,i));
            pass
        #
        dac=da.deepCopy();
        ids.setIJ(0,0,2); ids.setIJ(1,0,5); ids.setIJ(2,0,4);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues,1,da2,ids);
        ids.setIJ(0,0,2); ids.setIJ(1,0,2); ids.setIJ(2,0,-1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues,1,da2,ids);
        ids.setIJ(0,0,2); ids.setIJ(1,0,2); ids.setIJ(2,0,1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues,5,da2,ids);
        #
        ids.setIJ(0,0,2); ids.setIJ(1,0,2); ids.setIJ(2,0,1);
        dac=da.deepCopy();
        dac.setContigPartOfSelectedValues(4,da2,ids);
        expected6=[1,11,2,12,3,13,4,14,0,30,0,30,9,39]
        for i in range(14):
            self.assertEqual(expected6[i],dac.getIJ(0,i));
            pass
        pass

    def testBuildDescendingConnec2Of3DMesh1(self):
        mesh=MEDCouplingDataForTest.build3DSourceMesh_1();
        #
        mesh2,desc,descIndx,revDesc,revDescIndx=mesh.buildDescendingConnectivity2();
        mesh2.checkConsistencyLight();
        self.assertEqual(2,mesh2.getMeshDimension());
        self.assertEqual(30,mesh2.getNumberOfCells());
        self.assertEqual(31,revDescIndx.getNbOfElems()); self.assertEqual(31,revDescIndx.getNumberOfTuples());
        self.assertEqual(13,descIndx.getNbOfElems()); self.assertEqual(13,descIndx.getNumberOfTuples());
        self.assertEqual(48,desc.getNbOfElems()); self.assertEqual(48,desc.getNumberOfTuples());
        self.assertEqual(48,revDesc.getNbOfElems()); self.assertEqual(48,revDesc.getNumberOfTuples());
        expected1=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,-10,15,-5,-13,16,17,-14,18,-4,19,-2,20,21,22,23,24,25,-11,26,-1,-12,-25,-22,27,28,-7,-20,-24,29,-16,-18,30,-8,-28]
        self.assertEqual(expected1,desc.getValues());
        expected2=[0,4,8,12,16,20,24,28,32,36,40,44,48]
        self.assertEqual(expected2,descIndx.getValues());
        expected3=[0,2,4,5,7,9,10,12,14,15,17,19,21,23,25,26,28,29,31,32,34,35,37,38,40,42,43,44,46,47,48]
        self.assertEqual(expected3,revDescIndx.getValues());
        expected4=[0,8,0,6,0,0,5,1,4,1,1,9,1,11,2,2,3,2,7,2,8,3,4,3,5,3,4,10,4,5,11,5,6,10,6,6,9,7,7,10,7,8,8,9,9,11,10,11]
        self.assertEqual(expected4,revDesc.getValues());
        conn=mesh2.getNodalConnectivity();
        connIndex=mesh2.getNodalConnectivityIndex();
        expected5=[0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120]
        self.assertEqual(expected5,connIndex.getValues());
        expected6=[3,8,1,7,3,8,3,1,3,1,3,7,3,7,3,8,3,6,0,8,3,6,2,0,3,0,2,8,3,8,2,6,3,7,4,5,3,7,8,4,3,4,8,5,3,5,8,7,3,6,8,4,3,6,7,8,3,4,7,6,3,8,4,0,3,0,4,6,3,6,3,8,3,7,3,6,3,8,0,1,3,1,0,3,3,3,0,8,3,4,1,5,3,4,8,1,3,1,8,5,3,1,7,5,3,0,2,3,3,3,2,8,3,1,4,0,3,3,2,6]
        self.assertEqual(expected6,conn.getValues());
        pass

    def testAre2DCellsNotCorrectlyOriented1(self):
        m1Coords=[1.,1.,-1.,-1.,-1.,-1.,1.,-1.]
        m1Conn=[0,3,1,2]
        m1=MEDCouplingUMesh.New();
        m1.setMeshDimension(2);
        m1.allocateCells(1);
        m1.insertNextCell(NORM_QUAD4,4,m1Conn[0:4])
        m1.finishInsertingCells();
        myCoords1=DataArrayDouble.New();
        myCoords1.setValues(m1Coords,4,2);
        m1.setCoords(myCoords1);
        #
        vec1=[0.,0.,1.]
        for i in range(18):
            vec2=[3.*cos(pi/9.*i),3.*sin(pi/9.*i)];
            m1Cpy=m1.deepCopy();
            m1Cpy.translate(vec2);
            self.assertRaises(InterpKernelException,m1Cpy.are2DCellsNotCorrectlyOriented,vec1,False);
            m1Cpy.changeSpaceDimension(3);
            res=m1Cpy.are2DCellsNotCorrectlyOriented(vec1,False)
            self.assertEqual([0],res.getValues());
            pass
        pass

    def testDataArrayAbs1(self):
        d1=DataArrayDouble.New();
        val1=[2.,-3.,-5.,6.,-7.,-8.,9.,10.,-11.,-12.,-13.,-15.]
        expected1=[2.,3.,5.,6.,7.,8.,9.,10.,11.,12.,13.,15.]
        d1.setValues(val1,6,2);
        d2=d1.convertToIntArr();
        #
        d1.abs();
        for i in range(12):
            self.assertAlmostEqual(expected1[i],d1.getIJ(0,i),14);
            pass
        #
        expected2=[2,3,5,6,7,8,9,10,11,12,13,15]
        d2.abs();
        for i in range(12):
            self.assertEqual(expected2[i],d2.getIJ(0,i));
            pass
        #
        pass

    # test on 1D
    def testGetValueOn3(self):
        v=[0.,1.,1.5,2.]
        v2=[0.7,1.25,0.,2.,1.5]
        disp=[5.,50.,500.,6.,60.,600.,7.,70.,700.,8.,80.,800.]
        m=MEDCouplingUMesh.New("myMesh",1)
        nbNodes=len(v)
        nbCells=nbNodes-1
        m.allocateCells(nbCells)
        coords=DataArrayDouble.New() ; coords.setValues(v,nbNodes,1)
        m.setCoords(coords)
        m.insertNextCell(NORM_SEG2,2,[0,1])
        m.insertNextCell(NORM_SEG2,2,[2,1])
        m.insertNextCell(NORM_SEG2,2,[2,3])
        m.finishInsertingCells()
        f=MEDCouplingFieldDouble.New(ON_NODES)
        f.setMesh(m)
        array=DataArrayDouble.New(); array.setValues(disp,m.getNumberOfNodes(),3)
        f.setArray(array)
        arr1=f.getValueOnMulti(v2)
        self.assertEqual(5,arr1.getNumberOfTuples());
        self.assertEqual(3,arr1.getNumberOfComponents());
        expected1=[5.7,57.,570.,6.5,65.,650.,5.,50.,500.,8.,80.,800.,7.,70.,700.]
        for i in range(15):
            self.assertAlmostEqual(expected1[i],arr1.getIJ(0,i),14);
            pass
        pass

    def testGetNodeIdsOfCell2(self):
        m1c=MEDCouplingCMesh.New();
        coordsX=DataArrayDouble.New();
        arrX=[ -1., 1., 2., 4., 4.5 ]
        coordsX.setValues(arrX,5,1);
        coordsY=DataArrayDouble.New();
        arrY=[ -2., 2., 4., 8.]
        coordsY.setValues(arrY,4,1);
        coordsZ=DataArrayDouble.New();
        arrZ=[ -2., 2., 4.]
        coordsZ.setValues(arrZ,3,1);
        # test in 1D
        m1c.setCoordsAt(0,coordsX);
        expected1=[[0,1],[1,2],[2,3],[3,4]]
        self.assertEqual(4,m1c.getNumberOfCells())
        for i in range(m1c.getNumberOfCells()):
            self.assertEqual(expected1[i],m1c.getNodeIdsOfCell(i))
            pass
        # test in 2D
        m1c.setCoordsAt(1,coordsY);
        self.assertEqual(12,m1c.getNumberOfCells())
        self.assertEqual(20,m1c.getNumberOfNodes())
        expected2=[[0,1,6,5],[1,2,7,6],[2,3,8,7],[3,4,9,8],[5,6,11,10],[6,7,12,11],[7,8,13,12],[8,9,14,13],[10,11,16,15],[11,12,17,16],[12,13,18,17],[13,14,19,18]]
        for i in range(m1c.getNumberOfCells()):
            self.assertEqual(expected2[i],m1c.getNodeIdsOfCell(i))
            pass
        # test in 3D
        m1c.setCoordsAt(2,coordsZ);
        self.assertEqual(24,m1c.getNumberOfCells())
        self.assertEqual(60,m1c.getNumberOfNodes())
        expected3=[[0,1,6,5,20,21,26,25],[1,2,7,6,21,22,27,26],[2,3,8,7,22,23,28,27],[3,4,9,8,23,24,29,28],[5,6,11,10,25,26,31,30],[6,7,12,11,26,27,32,31],[7,8,13,12,27,28,33,32],[8,9,14,13,28,29,34,33],[10,11,16,15,30,31,36,35],[11,12,17,16,31,32,37,36],[12,13,18,17,32,33,38,37],[13,14,19,18,33,34,39,38],[20,21,26,25,40,41,46,45],[21,22,27,26,41,42,47,46],[22,23,28,27,42,43,48,47],[23,24,29,28,43,44,49,48],[25,26,31,30,45,46,51,50],[26,27,32,31,46,47,52,51],[27,28,33,32,47,48,53,52],[28,29,34,33,48,49,54,53],[30,31,36,35,50,51,56,55],[31,32,37,36,51,52,57,56],[32,33,38,37,52,53,58,57],[33,34,39,38,53,54,59,58]]
        self.assertEqual(24,m1c.getNumberOfCells())
        for i in range(m1c.getNumberOfCells()):
            self.assertEqual(expected3[i],m1c.getNodeIdsOfCell(i))
            pass
        pass
          
    pass

if __name__ == '__main__':
    unittest.main()
