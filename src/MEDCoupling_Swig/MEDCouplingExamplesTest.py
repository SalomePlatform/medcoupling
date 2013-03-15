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

from MEDCoupling import *
import unittest

class MEDCouplingBasicsTest(unittest.TestCase):
    def testExample_DataArrayInt_(self):
#! [PySnippet_DataArrayInt__1]
        pass

    def testExample_DataArrayInt_getTuple(self):
#! [Snippet_DataArrayInt_getTuple_1]
        dv=DataArrayInt.New();
        dv.alloc( 6, 1 )
        dv.iota(7)
        dv.rearrange( 2 )
        print dv.getTuple( 1 )
#! [Snippet_DataArrayInt_getTuple_1]
#! [Snippet_DataArrayInt_getTuple_2]
        for tpl in dv:
            print tpl
#! [Snippet_DataArrayInt_getTuple_2]
        pass

    def testExample_DataArrayInt_buildPermutationArr(self):
#! [PySnippet_DataArrayInt_buildPermutationArr_1]
        a=DataArrayInt.New()
        a.setValues([4,5,6,7,8],5,1)
        b=DataArrayInt.New()
        b.setValues([5,4,8,6,7],5,1)
        c=a.buildPermutationArr(b)
#! [PySnippet_DataArrayInt_buildPermutationArr_1]
        self.assertEqual([1,0,4,2,3],c.getValues())
        pass

    def testExample_DataArrayInt_invertArrayO2N2N2O(self):
#! [PySnippet_DataArrayInt_invertArrayO2N2N2O_1]
        arr1=[2,0,4,1,5,3]
        da=DataArrayInt.New();
        da.setValues(arr1,6,1);
        da2=da.invertArrayO2N2N2O(6);
        expected1=[1,3,0,5,2,4]
        for i in xrange(6):
            self.assertEqual(expected1[i],da2.getIJ(i,0));
            pass
#! [PySnippet_DataArrayInt_invertArrayO2N2N2O_1]
        pass

    def testExample_DataArrayInt_invertArrayN2O2O2N(self):
#! [PySnippet_DataArrayInt_invertArrayN2O2O2N_1]
        arr1=[2,0,4,1,5,3]
        da=DataArrayInt.New();
        da.setValues(arr1,6,1);
        da2=da.invertArrayN2O2O2N(7);
        expected1=[1,3,0,5,2,4,-1]
        for i in xrange(6):
            self.assertEqual(expected1[i],da2.getIJ(i,0));
            pass
#! [PySnippet_DataArrayInt_invertArrayN2O2O2N_1]
        pass


    def testExample_DataArrayDouble_getIdsInRange(self):
#! [PySnippet_DataArrayDouble_getIdsInRange_1]
        da=DataArrayDouble.New()
        da.alloc( 10, 1 )
        da[ :, :] = range(10)
        da2 = da.getIdsInRange( 2.5, 6 )
#! [PySnippet_DataArrayDouble_getIdsInRange_1]
        pass

    def testExample_DataArrayDouble_setPartOfValues2(self):
#! [Snippet_DataArrayDouble_setPartOfValues2_1]
        da=DataArrayDouble.New()
        da.alloc( 4, 7 )
        #
        dv=DataArrayDouble.New();
        dv.alloc( 6, 1 )
        dv.iota(7)
        dv.rearrange( 2 )
#! [Snippet_DataArrayDouble_setPartOfValues2_1]
#! [Snippet_DataArrayDouble_setPartOfValues2_2]
        da.fillWithZero()
        da[ [0,1,2], [1,3] ] = dv
#! [Snippet_DataArrayDouble_setPartOfValues2_2]
#! [Snippet_DataArrayDouble_setPartOfValues2_3]
        da.fillWithZero()
        dv.rearrange( 6 )
        da[ [0,2,3], [0,2,3,4,5,6]] = dv
#! [Snippet_DataArrayDouble_setPartOfValues2_3]
        pass

    def testExample_DataArrayInt_setPartOfValues2(self):
#! [Snippet_DataArrayInt_setPartOfValues2_1]
        da=DataArrayInt.New()
        da.alloc( 4, 7 )
        #
        dv=DataArrayInt.New();
        dv.alloc( 6, 1 )
        dv.iota(7)
        dv.rearrange( 2 )
#! [Snippet_DataArrayInt_setPartOfValues2_1]
#! [Snippet_DataArrayInt_setPartOfValues2_2]
        da.fillWithZero()
        da[ [0,1,2], [1,3] ] = dv
#! [Snippet_DataArrayInt_setPartOfValues2_2]
#! [Snippet_DataArrayInt_setPartOfValues2_3]
        da.fillWithZero()
        dv.rearrange( 6 )
        da[ [0,2,3], [0,2,3,4,5,6]] = dv
#! [Snippet_DataArrayInt_setPartOfValues2_3]
        pass

    def testExample_DataArrayDouble_setPartOfValues3(self):
#! [Snippet_DataArrayDouble_setPartOfValues3_1]
        da=DataArrayDouble.New()
        da.alloc( 4, 7 )
        #
        dv=DataArrayDouble.New();
        dv.alloc( 6, 1 )
        dv.iota(7)
        dv.rearrange( 2 )
#! [Snippet_DataArrayDouble_setPartOfValues3_1]
#! [Snippet_DataArrayDouble_setPartOfValues3_2]
        da.fillWithZero()
        da[ 0:3, [1,3] ] = dv
#! [Snippet_DataArrayDouble_setPartOfValues3_2]
#! [Snippet_DataArrayDouble_setPartOfValues3_3]
        da.fillWithZero()
        dv.rearrange( 6 )
        da[ 0:4:2, [0,2,3,4,5,6]] = dv
#! [Snippet_DataArrayDouble_setPartOfValues3_3]
        pass

    def testExample_DataArrayInt_setPartOfValues3(self):
#! [Snippet_DataArrayInt_setPartOfValues3_1]
        da=DataArrayInt.New()
        da.alloc( 4, 7 )
        #
        dv=DataArrayInt.New();
        dv.alloc( 6, 1 )
        dv.iota(7)
        dv.rearrange( 2 )
#! [Snippet_DataArrayInt_setPartOfValues3_1]
#! [Snippet_DataArrayInt_setPartOfValues3_2]
        da.fillWithZero()
        da[ 0:3, [1,3] ] = dv
#! [Snippet_DataArrayInt_setPartOfValues3_2]
#! [Snippet_DataArrayInt_setPartOfValues3_3]
        da.fillWithZero()
        dv.rearrange( 6 )
        da[ 0:4:2, [0,2,3,4,5,6]] = dv
#! [Snippet_DataArrayInt_setPartOfValues3_3]
        pass

    def testExample_DataArrayDouble_setPartOfValues1(self):
#! [Snippet_DataArrayDouble_setPartOfValues1_1]
        da=DataArrayDouble.New()
        da.alloc( 4, 4 )
        da.setInfoOnComponents( ["v1","v2","v3","v4"])
        #
        dv=DataArrayDouble.New();
        dv.alloc( 4, 1 )
        dv.iota(7)
        dv.rearrange( 2 )
        dv.setInfoOnComponents( ["a1","a2"])
#! [Snippet_DataArrayDouble_setPartOfValues1_1]
#! [Snippet_DataArrayDouble_setPartOfValues1_2]
        da.fillWithZero()
        da.setPartOfValues1( dv, 1,3,1, 1,3,1, True )
#! [Snippet_DataArrayDouble_setPartOfValues1_2]
#! [Snippet_DataArrayDouble_setPartOfValues1_3]
        da.fillWithZero()
        da.setPartOfValues1( dv, 0,4,1, 1,2,1, False )
#! [Snippet_DataArrayDouble_setPartOfValues1_3]
#! [Snippet_DataArrayDouble_setPartOfValues1_4]
        da.fillWithZero()
        da.setPartOfValues1( dv, 1,2,1, 0,4,1, False )
#! [Snippet_DataArrayDouble_setPartOfValues1_4]
#! [Snippet_DataArrayDouble_setPartOfValues1_5]
        da.fillWithZero()
        da.setPartOfValues1( dv, 0,3,2, 1,4,2, True )
#! [Snippet_DataArrayDouble_setPartOfValues1_5]
#! [Snippet_DataArrayDouble_setPartOfValues1_6]
        da2 = da.deepCpy()
        da2.fillWithZero()
        da2[ 0:3:2, 1:4:2 ] = dv
        self.assertTrue( da.isEqual( da2, 1e-20 ))
#! [Snippet_DataArrayDouble_setPartOfValues1_6]
        pass

    def testExample_DataArrayInt_setPartOfValues1(self):
#! [Snippet_DataArrayInt_setPartOfValues1_1]
        da=DataArrayInt.New()
        da.alloc( 4, 4 )
        da.setInfoOnComponents( ["v1","v2","v3","v4"])
        #
        dv=DataArrayInt.New();
        dv.alloc( 4, 1 )
        dv.iota(7)
        dv.rearrange( 2 )
        dv.setInfoOnComponents( ["a1","a2"])
#! [Snippet_DataArrayInt_setPartOfValues1_1]
#! [Snippet_DataArrayInt_setPartOfValues1_2]
        da.fillWithZero()
        da.setPartOfValues1( dv, 1,3,1, 1,3,1, True )
#! [Snippet_DataArrayInt_setPartOfValues1_2]
#! [Snippet_DataArrayInt_setPartOfValues1_3]
        da.fillWithZero()
        da.setPartOfValues1( dv, 0,4,1, 1,2,1, False )
#! [Snippet_DataArrayInt_setPartOfValues1_3]
#! [Snippet_DataArrayInt_setPartOfValues1_4]
        da.fillWithZero()
        da.setPartOfValues1( dv, 1,2,1, 0,4,1, False )
#! [Snippet_DataArrayInt_setPartOfValues1_4]
#! [Snippet_DataArrayInt_setPartOfValues1_5]
        da.fillWithZero()
        da.setPartOfValues1( dv, 0,3,2, 1,4,2, True )
#! [Snippet_DataArrayInt_setPartOfValues1_5]
#! [Snippet_DataArrayInt_setPartOfValues1_6]
        da2 = da.deepCpy()
        da2.fillWithZero()
        da2[ 0:3:2, 1:4:2 ] = dv
        self.assertTrue( da.isEqual( da2 ))
#! [Snippet_DataArrayInt_setPartOfValues1_6]
        pass

    def testExample_DataArrayDouble_setPartOfValuesSimple1(self):
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_1]
        da=DataArrayDouble.New()
        da.alloc( 4, 4 )
        dv = 7
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_1]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_2]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 1,3,1, 1,3,1 )
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_2]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_3]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 0,4,1, 1,2,1 )
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_3]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_4]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 1,2,1, 0,4,1 )
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_4]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_5]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 0,3,2, 1,4,2 )
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_5]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_6]
        da2 = da.deepCpy()
        da2.fillWithZero()
        da2[ 0:3:2, 1:4:2 ] = dv
        self.assertTrue( da.isEqual( da2, 1e-20 ))
#! [Snippet_DataArrayDouble_setPartOfValuesSimple1_6]
        pass

    def testExample_DataArrayInt_setPartOfValuesSimple1(self):
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_1]
        da=DataArrayInt.New()
        da.alloc( 4, 4 )
        dv = 7
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_1]
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_2]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 1,3,1, 1,3,1 )
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_2]
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_3]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 0,4,1, 1,2,1 )
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_3]
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_4]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 1,2,1, 0,4,1 )
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_4]
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_5]
        da.fillWithZero()
        da.setPartOfValuesSimple1( dv, 0,3,2, 1,4,2 )
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_5]
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_6]
        da2 = da.deepCpy()
        da2.fillWithZero()
        da2[ 0:3:2, 1:4:2 ] = dv
        self.assertTrue( da.isEqual( da2 ))
#! [Snippet_DataArrayInt_setPartOfValuesSimple1_6]
        pass

    def testExample_DataArrayDouble_setPartOfValuesSimple2(self):
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_1]
        da=DataArrayDouble.New()
        da.alloc( 4, 4 )
        dv = 7
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_1]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_2]
        da.fillWithZero()
        da[[1,2], [1,2]] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_2]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_3]
        da.fillWithZero()
        da[[0,1,2,3], [1]] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_3]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_4]
        da.fillWithZero()
        da[[1], [0,1,2,3]] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_4]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_5]
        da.fillWithZero()
        da[[0,2], [1,3]] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple2_5]
        pass

    def testExample_DataArrayInt_setPartOfValuesSimple2(self):
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_1]
        da=DataArrayInt.New()
        da.alloc( 4, 4 )
        dv = 7
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_1]
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_2]
        da.fillWithZero()
        da[[1,2], [1,2]] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_2]
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_3]
        da.fillWithZero()
        da[[0,1,2,3], [1]] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_3]
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_4]
        da.fillWithZero()
        da[[1], [0,1,2,3]] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_4]
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_5]
        da.fillWithZero()
        da[[0,2], [1,3]] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple2_5]
        pass

    def testExample_DataArrayDouble_setPartOfValuesSimple3(self):
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_1]
        da=DataArrayDouble.New()
        da.alloc( 4, 4 )
        dv = 7
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_1]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_2]
        da.fillWithZero()
        da[[1,2], 1:3] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_2]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_3]
        da.fillWithZero()
        da[[0,1,2,3], 1:2] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_3]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_4]
        da.fillWithZero()
        da[[1], 0:4] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_4]
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_5]
        da.fillWithZero()
        da[[0,2], 1:4:2] = dv
#! [Snippet_DataArrayDouble_setPartOfValuesSimple3_5]
        pass

    def testExample_DataArrayInt_setPartOfValuesSimple3(self):
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_1]
        da=DataArrayInt.New()
        da.alloc( 4, 4 )
        dv = 7
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_1]
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_2]
        da.fillWithZero()
        da[[1,2], 1:3] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_2]
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_3]
        da.fillWithZero()
        da[[0,1,2,3], 1:2] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_3]
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_4]
        da.fillWithZero()
        da[[1], 0:4] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_4]
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_5]
        da.fillWithZero()
        da[[0,2], 1:4:2] = dv
#! [Snippet_DataArrayInt_setPartOfValuesSimple3_5]
        pass

    def testExample_DataArrayDouble_setSelectedComponents(self):
#! [Snippet_DataArrayDouble_setSelectedComponents1]
        da=DataArrayDouble.New();
        array1=[1.,2., 3.,4., 5.,6.]
        da.setValues(array1,3,2)
        da.setInfoOnComponents( ["a1","a2"])
#! [Snippet_DataArrayDouble_setSelectedComponents1]
#! [Snippet_DataArrayDouble_setSelectedComponents2]
        dv=DataArrayDouble.New();
        dv.alloc( 4, 4 )
        dv.fillWithZero()
        dv.setInfoOnComponents( ["v1","v2","v3","v4"])
        dv2 = dv.deepCpy()
        dv.setSelectedComponents( da, [1,0] )
#! [Snippet_DataArrayDouble_setSelectedComponents2]
#! [Snippet_DataArrayDouble_setSelectedComponents3]
        dv2[:3,[1,0]] = da
        self.assertTrue( dv.isEqualWithoutConsideringStr( dv2, 1e-20 ))
#! [Snippet_DataArrayDouble_setSelectedComponents3]
        pass

    def testExample_DataArrayInt_setSelectedComponents(self):
#! [Snippet_DataArrayInt_setSelectedComponents1]
        da=DataArrayInt.New();
        array1=[1,2, 3,4, 5,6]
        da.setValues(array1,3,2)
        da.setInfoOnComponents( ["a1","a2"])
#! [Snippet_DataArrayInt_setSelectedComponents1]
#! [Snippet_DataArrayInt_setSelectedComponents2]
        dv=DataArrayInt.New();
        dv.alloc( 4, 4 )
        dv.fillWithZero()
        dv.setInfoOnComponents( ["v1","v2","v3","v4"])
        dv2 = dv.deepCpy()
        dv.setSelectedComponents( da, [1,0] )
#! [Snippet_DataArrayInt_setSelectedComponents2]
#! [Snippet_DataArrayInt_setSelectedComponents3]
        dv2[:3,[1,0]] = da
        self.assertTrue( dv.isEqualWithoutConsideringStr( dv2 ))
#! [Snippet_DataArrayInt_setSelectedComponents3]
        pass

    def testExample_DataArrayDouble_getDifferentValues(self):
#! [Snippet_DataArrayDouble_getDifferentValues1]
        da=DataArrayDouble.New();
        array1=[2.3,1.2,1.3,2.3,2.301,0.8]
        da.setValues(array1,6,1)
        #
        dv=da.getDifferentValues(2e-1);
        expected2=[2.301,1.3,0.8]
        self.assertEqual(3,dv.getNbOfElems());
        for i in xrange(3):
            self.assertAlmostEqual(expected2[i],dv.getIJ(i,0),14);
            pass
#! [Snippet_DataArrayDouble_getDifferentValues1]
        pass

    def testExample_DataArrayDouble_findCommonTuples1(self):
#! [PySnippet_DataArrayDouble_findCommonTuples1]
        da=DataArrayDouble.New();
        array2=[2.3,2.3, 1.2,1.2, 1.3,1.3, 2.3,2.3, 2.301,2.301, 0.8,0.8]
        da.setValues(array2,6,2)
#! [PySnippet_DataArrayDouble_findCommonTuples1]
#! [PySnippet_DataArrayDouble_findCommonTuples2]
        c,cI=da.findCommonTuples(1e-1);
        expected3=[0,3,4,1,2]
        expected4=[0,3,5]
        self.assertEqual(expected3,c.getValues())
        self.assertEqual(expected4,cI.getValues())
#! [PySnippet_DataArrayDouble_findCommonTuples2]
        pass

    def testExampleDataArrayDoubleMeldWith(self):
#! [PySnippet_DataArrayDouble_Meld1_1]
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
        da1C=da1.deepCpy();
        da1.meldWith(da3);
#! [PySnippet_DataArrayDouble_Meld1_1]

    def testExampleDataArrayIntMeldWith(self):
#! [PySnippet_DataArrayInt_Meld1_1]
        da1=DataArrayInt.New();
        da1.alloc(7,2);
        da2=DataArrayInt.New();
        da2.alloc(7,1);
        #
        da1.fillWithValue(7);
        da2.iota(0);
        #
        da1.setInfoOnComponent(0,"c0da1");
        da1.setInfoOnComponent(1,"c1da1");
        da2.setInfoOnComponent(0,"c0da2");
        #
        da1.meldWith(da2);
#! [PySnippet_DataArrayInt_Meld1_1]

    def testExampleDataArrayDoubleKeepSelectedComponents1(self):
#! [SnippeDataArrayDoubleKeepSelectedComponents1_1]
        arr1=[1.,2.,3.,4.,     # tuple 0
              11.,12.,13.,14., # tuple 1
              21.,22.,23.,24., # ...
              31.,32.,33.,34.,
              41.,42.,43.,44.]
        a1=DataArrayDouble.New()
        a1.setValues(arr1,5,4)
        a1.setInfoOnComponent(0,"a");
        a1.setInfoOnComponent(1,"b");
        a1.setInfoOnComponent(2,"c");
        a1.setInfoOnComponent(3,"d");
#! [SnippeDataArrayDoubleKeepSelectedComponents1_1]
#! [SnippeDataArrayDoubleKeepSelectedComponents1_2]
        arr2V=[1,2,1,2,0,0]
        a2=a1.keepSelectedComponents(arr2V)
#! [SnippeDataArrayDoubleKeepSelectedComponents1_2]
        pass

    def testExampleDataArrayIntKeepSelectedComponents1(self):
#! [SnippeDataArrayIntKeepSelectedComponents1_1]
        arr1=[1,2,3,4,     # tuple 0
              11,12,13,14, # tuple 1
              21,22,23,24, # 
              31,32,33,34,
              41,42,43,44]
        a1=DataArrayInt.New()
        a1.setValues(arr1,5,4)
        a1.setInfoOnComponent(0,"a");
        a1.setInfoOnComponent(1,"b");
        a1.setInfoOnComponent(2,"c");
        a1.setInfoOnComponent(3,"d");
#! [SnippeDataArrayIntKeepSelectedComponents1_1]
#! [SnippeDataArrayIntKeepSelectedComponents1_2]
        arr2V=[1,2,1,2,0,0]
        a2=a1.keepSelectedComponents(arr2V)
#! [SnippeDataArrayIntKeepSelectedComponents1_2]
#! [SnippeDataArrayIntKeepSelectedComponents1_3]
        a3=a1[:,arr2V ]
#! [SnippeDataArrayIntKeepSelectedComponents1_3]
        pass

    def testExampleFieldDoubleBuildSubPart1(self):
        from MEDCouplingDataForTest import MEDCouplingDataForTest
#! [PySnippetFieldDoubleBuildSubPart1_1]
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1()
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
        f1.setTime(2.3,5,6)
        f1.setMesh(mesh1)
        array=DataArrayDouble.New()
        arr1=[3.,103.,4.,104.,5.,105.,6.,106.,7.,107.]
        array.setValues(arr1,mesh1.getNumberOfCells(),2)
        f1.setArray(array)
# ! [PySnippetFieldDoubleBuildSubPart1_1]
# ! [PySnippetFieldDoubleBuildSubPart1_2]
        part1=[2,1,4]
        f2=f1.buildSubPart(part1)
# ! [PySnippetFieldDoubleBuildSubPart1_2]
        f2.zipCoords()
        self.assertEqual(3,f2.getNumberOfTuples())
        self.assertEqual(2,f2.getNumberOfComponents())
        expected1=[5.,105.,4.,104.,7.,107.]
        for i in xrange(6):
            self.assertAlmostEqual(f2.getIJ(0,i),expected1[i],12)
            pass
        self.assertEqual(3,f2.getMesh().getNumberOfCells())
        self.assertEqual(6,f2.getMesh().getNumberOfNodes())
        self.assertEqual(2,f2.getMesh().getSpaceDimension())
        self.assertEqual(2,f2.getMesh().getMeshDimension())
        m2C=f2.getMesh()
        self.assertEqual(13,m2C.getMeshLength())
        expected2=[0.2, -0.3, 0.7, -0.3, 0.2, 0.2, 0.7, 0.2, 0.2, 0.7, 0.7, 0.7]
        for i in xrange(12):
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12)
            pass
        expected3=[3,2,3,1,3,0,2,1,4,4,5,3,2]
        self.assertEqual(expected3,list(m2C.getNodalConnectivity().getValues()))
        expected4=[0,4,8,13]
        self.assertEqual(expected4,list(m2C.getNodalConnectivityIndex().getValues()))
        # Test with field on nodes.
# ! [PySnippetFieldDoubleBuildSubPart1_3]
        f1=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME)
        f1.setTime(2.3,5,6)
        f1.setMesh(mesh1)
        array=DataArrayDouble.New()
        arr2=[3.,103.,4.,104.,5.,105.,6.,106.,7.,107.,8.,108.,9.,109.,10.,110.,11.,111.]
        array.setValues(arr2,mesh1.getNumberOfNodes(),2)
        f1.setArray(array)
# ! [PySnippetFieldDoubleBuildSubPart1_3]
# ! [PySnippetFieldDoubleBuildSubPart1_4]
        part2=[1,2]
        f2=f1.buildSubPart(part2)
# ! [PySnippetFieldDoubleBuildSubPart1_4]
        self.assertEqual(4,f2.getNumberOfTuples())
        self.assertEqual(2,f2.getNumberOfComponents())
        expected5=[4.,104.,5.,105.,7.,107.,8.,108.]
        for i in xrange(8):
            self.assertAlmostEqual(f2.getIJ(0,i),expected5[i],12)
            pass
        self.assertEqual(2,f2.getMesh().getNumberOfCells())
        self.assertEqual(4,f2.getMesh().getNumberOfNodes())
        self.assertEqual(2,f2.getMesh().getSpaceDimension())
        self.assertEqual(2,f2.getMesh().getMeshDimension())
        m2C=f2.getMesh()
        self.assertEqual(8,m2C.getMeshLength())
        for i in xrange(8):#8 is not an error
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12)
            pass
        self.assertEqual(expected3[:4],list(m2C.getNodalConnectivity().getValues())[4:])
        self.assertEqual(expected3[4:8],list(m2C.getNodalConnectivity().getValues())[:4])
        self.assertEqual(expected4[:3],list(m2C.getNodalConnectivityIndex().getValues()))
        #idem previous because nodes of cell#4 are not fully present in part3
        part3=[1,2]
        arrr=DataArrayInt.New()
        arrr.setValues(part3,2,1)
        f2=f1.buildSubPart(arrr)
        self.assertEqual(4,f2.getNumberOfTuples())
        self.assertEqual(2,f2.getNumberOfComponents())
        for i in xrange(8):
            self.assertAlmostEqual(f2.getIJ(0,i),expected5[i],12)
            pass
        self.assertEqual(2,f2.getMesh().getNumberOfCells())
        self.assertEqual(4,f2.getMesh().getNumberOfNodes())
        self.assertEqual(2,f2.getMesh().getSpaceDimension())
        self.assertEqual(2,f2.getMesh().getMeshDimension())
        m2C=f2.getMesh()
        self.assertEqual(8,m2C.getMeshLength())
        for i in xrange(8):#8 is not an error
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12)
            pass
        self.assertEqual(expected3[:4],list(m2C.getNodalConnectivity().getValues())[4:8])
        self.assertEqual(expected3[4:8],list(m2C.getNodalConnectivity().getValues())[:4])
        self.assertEqual(expected4[:3],list(m2C.getNodalConnectivityIndex().getValues()))
        part4=[1,2,4]
        f2=f1.buildSubPart(part4)
        self.assertEqual(6,f2.getNumberOfTuples())
        self.assertEqual(2,f2.getNumberOfComponents())
        expected6=[4.,104.,5.,105.,7.,107.,8.,108.,10.,110.,11.,111.]
        for i in xrange(12):
            self.assertAlmostEqual(f2.getIJ(0,i),expected6[i],12)
            pass
        self.assertEqual(3,f2.getMesh().getNumberOfCells())
        self.assertEqual(6,f2.getMesh().getNumberOfNodes())
        self.assertEqual(2,f2.getMesh().getSpaceDimension())
        self.assertEqual(2,f2.getMesh().getMeshDimension())
        m2C=f2.getMesh()
        self.assertEqual(13,m2C.getMeshLength())
        for i in xrange(12):
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12)
            pass
        self.assertEqual(expected3[0:4],list(m2C.getNodalConnectivity().getValues())[4:8])
        self.assertEqual(expected3[4:8],list(m2C.getNodalConnectivity().getValues())[0:4])
        self.assertEqual(expected3[8:13],list(m2C.getNodalConnectivity().getValues())[8:13])
        self.assertEqual(expected4,list(m2C.getNodalConnectivityIndex().getValues()))
        pass

    def testExampleUMeshStdBuild1(self):
# ! [PySnippetUMeshStdBuild1_1]
        coords=[-0.3,-0.3,0.,   0.2,-0.3,0.,   0.7,-0.3,0.,   -0.3,0.2,0.,   0.2,0.2,0., 
                 0.7,0.2,0.,    -0.3,0.7,0.,    0.2,0.7,0.,     0.7,0.7,0. ]
        nodalConnPerCell=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
# ! [PySnippetUMeshStdBuild1_1]
# ! [PySnippetUMeshStdBuild1_2]
        mesh=MEDCouplingUMesh.New("My2DMesh",2)
# ! [PySnippetUMeshStdBuild1_2]
# ! [PySnippetUMeshStdBuild1_3]
        mesh.allocateCells(5)#You can put more than 5 if you want but not less.
        mesh.insertNextCell(NORM_QUAD4,nodalConnPerCell[:4])
        mesh.insertNextCell(NORM_TRI3,nodalConnPerCell[4:7])
        mesh.insertNextCell(NORM_TRI3,nodalConnPerCell[7:10])
        mesh.insertNextCell(NORM_QUAD4,nodalConnPerCell[10:14])
        mesh.insertNextCell(NORM_QUAD4,nodalConnPerCell[14:])
        mesh.finishInsertingCells()
# ! [PySnippetUMeshStdBuild1_3]
# ! [PySnippetUMeshStdBuild1_4]
        myCoords=DataArrayDouble.New(coords,9,3)#here myCoords are declared to have 3 components, mesh will deduce that its spaceDim==3. 
        mesh.setCoords(myCoords)#myCorrds contains 9 tuples, that is to say mesh contains 9 nodes.
# ! [PySnippetUMeshStdBuild1_4]
# ! [PySnippetUMeshStdBuild1_5]
# ! [PySnippetUMeshStdBuild1_5]
        mesh.checkCoherency()
        pass

    def testExampleCMeshStdBuild1(self):
# ! [PySnippetCMeshStdBuild1_1]
        XCoords=[-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22] # 9 values along X
        YCoords=[0.,0.1,0.37,0.45,0.47,0.49,1.007] # 7 values along Y
        arrX=DataArrayDouble.New(XCoords)
        arrX.setInfoOnComponent(0,"X [m]")
        arrY=DataArrayDouble.New(YCoords)
        arrY.setInfoOnComponent(0,"Y [m]")
# ! [PySnippetCMeshStdBuild1_1]
# ! [PySnippetCMeshStdBuild1_2]
        mesh=MEDCouplingCMesh.New("My2D_CMesh")
        mesh.setCoords(arrX,arrY)
# ! [PySnippetCMeshStdBuild1_2]
# ! [PySnippetCMeshStdBuild1_3]
        self.assertEqual(8*6,mesh.getNumberOfCells())
        self.assertEqual(9*7,mesh.getNumberOfNodes())
        self.assertEqual(2,mesh.getSpaceDimension())
        self.assertEqual(2,mesh.getMeshDimension())
# ! [PySnippetCMeshStdBuild1_3]
        mesh=MEDCouplingCMesh.New("My2D_CMesh")
# ! [PySnippetCMeshStdBuild1_2bis]
        mesh.setCoordsAt(0,arrX)
        mesh.setCoordsAt(1,arrY)
# ! [PySnippetCMeshStdBuild1_2bis]
        self.assertEqual(8*6,mesh.getNumberOfCells())
        self.assertEqual(9*7,mesh.getNumberOfNodes())
        self.assertEqual(2,mesh.getSpaceDimension())
        self.assertEqual(2,mesh.getMeshDimension())
        pass

    def testExampleUMeshAdvBuild1(self):
# ! [PySnippetUMeshAdvBuild1_1]
        coords=[-0.3,-0.3,0.,   0.2,-0.3,0.,   0.7,-0.3,0.,   -0.3,0.2,0.,   0.2,0.2,0., 
                 0.7,0.2,0.,    -0.3,0.7,0.,    0.2,0.7,0.,     0.7,0.7,0. ]
        nodalConnPerCell=[4,0,3,4,1, 3,1,4,2, 3,4,5,2, 4,6,7,4,3, 4,7,8,5,4]
        nodalConnPerCellIndex=[0,5,9,13,18,23]
# ! [PySnippetUMeshAdvBuild1_1]
# ! [PySnippetUMeshAdvBuild1_2]
        mesh=MEDCouplingUMesh.New("My2DMesh",2)
# ! [PySnippetUMeshAdvBuild1_2]
# ! [PySnippetUMeshAdvBuild1_3]
        nodalConn=DataArrayInt.New(nodalConnPerCell,23,1)
        nodalConnI=DataArrayInt.New(nodalConnPerCellIndex,6,1)
        mesh.setConnectivity(nodalConn,nodalConnI,True)
# ! [PySnippetUMeshAdvBuild1_3]
# ! [PySnippetUMeshAdvBuild1_4]
        myCoords=DataArrayDouble.New(coords,9,3)#here myCoords are declared to have 3 components, mesh will deduce that its spaceDim==3. 
        mesh.setCoords(myCoords)#myCorrds contains 9 tuples, that is to say mesh contains 9 nodes.
# ! [PySnippetUMeshAdvBuild1_4]
# ! [PySnippetUMeshAdvBuild1_5]
# ! [PySnippetUMeshAdvBuild1_5]
        mesh.checkCoherency()
        pass

    def testExampleDataArrayBuild1(self):
# ! [PySnippetDataArrayBuild1_0]
        dataDouble=[0.,10.,20.,1.,11.,21.,2.,12.,22.,3.,13.,23.,4.,14.,24.]
# ! [PySnippetDataArrayBuild1_0]
# ! [PySnippetDataArrayBuild1_1]
        arrayDouble=DataArrayDouble.New()
        arrayDouble.setValues(dataDouble,5,3)# 5 tuples containing each 3 components
# ! [PySnippetDataArrayBuild1_1]
# ! [PySnippetDataArrayBuild1_1bis]
        arrayDouble=DataArrayDouble.New(dataDouble,5,3)
# ! [PySnippetDataArrayBuild1_1bis]
# ! [PySnippetDataArrayBuild1_2]
        dataInt=[0, 10, 20, 1, 11, 21, 2, 12, 22, 3, 13, 23, 4, 14, 24]
# ! [PySnippetDataArrayBuild1_2]
# ! [PySnippetDataArrayBuild1_3]
        arrayInt=DataArrayInt.New()
        arrayInt.setValues(dataInt,5,3)# 5 tuples containing each 3 components
# ! [PySnippetDataArrayBuild1_3]
# ! [PySnippetDataArrayBuild1_3bis]
        arrayInt=DataArrayInt.New(dataInt,5,3)
# ! [PySnippetDataArrayBuild1_3bis]
        pass

    def testExampleFieldDoubleBuild1(self):
        XCoords=[-0.3,0.07,0.1,0.3,0.45,0.47,0.49,1.,1.22] ; arrX=DataArrayDouble.New(XCoords)
        YCoords=[0.07,0.1,0.37,0.45,0.47,0.49,1.007] ; arrY=DataArrayDouble.New(YCoords)
        mesh=MEDCouplingCMesh.New("My2D_CMesh")
        mesh.setCoords(arrX,arrY)
# ! [PySnippetFieldDoubleBuild1_1]
        fieldOnCells=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME)
        fieldOnCells.setName("MyTensorFieldOnCellNoTime")
        fieldOnCells.setMesh(mesh)
        array=DataArrayDouble.New()
        array.alloc(fieldOnCells.getMesh().getNumberOfCells(),9) # Implicitely fieldOnCells will be a 9 components field.
        array.fillWithValue(7.)
        fieldOnCells.setArray(array)
        # fieldOnCells is now usable
        # ...
# ! [PySnippetFieldDoubleBuild1_1]
# ! [PySnippetFieldDoubleBuild1_2]
        f1=mesh.fillFromAnalytic(ON_CELLS,1,"x*x+y*y*3+2.*x") # f1 is scalar
        f2=mesh.fillFromAnalytic(ON_CELLS,1,"cos(x+y/x)") # f2 is scalar too
        f2bis=mesh.fillFromAnalytic(ON_CELLS,2,"x*x*IVec+3*y*JVec") # f2bis is a vectors field
        f3=f1+f2 # f3 scalar
        f4=f3/f2 # f4 scalar
        f2bis.applyFunc(1,"sqrt(x*x+y*y)") # f2bis becomes scalar
        f5=f2bis*f4 # f5 scalar
        pos1=[0.48,0.38]
        res=f4.getValueOn(pos1) # f4 is scalar so the returned value is of size 1.
        # ...
# ! [PySnippetFieldDoubleBuild1_2]
# ! [PySnippetFieldDoubleBuild1_3]
# ! [PySnippetFieldDoubleBuild1_3]
        pass

    def testExampleFieldDoubleBuild2(self):
        XCoords=[-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22] ; arrX=DataArrayDouble.New(XCoords)
        YCoords=[0.,0.1,0.37,0.45,0.47,0.49,1.007] ; arrY=DataArrayDouble.New(YCoords)
        mesh=MEDCouplingCMesh.New("My2D_CMesh")
        mesh.setCoords(arrX,arrY)
# ! [PySnippetFieldDoubleBuild2_1]
        fieldOnNodes=MEDCouplingFieldDouble.New(ON_NODES,NO_TIME)
        fieldOnNodes.setName("MyScalarFieldOnNodeNoTime")
        fieldOnNodes.setMesh(mesh)
        array=DataArrayDouble.New()
        array.alloc(fieldOnNodes.getMesh().getNumberOfNodes(),1) # Implicitely fieldOnNodes will be a 1 component field.
        array.fillWithValue(7.)
        fieldOnNodes.setArray(array)
        # fieldOnNodes is now usable
        # ...
# ! [PySnippetFieldDoubleBuild2_1]
        pass

    def testExampleFieldDoubleBuild3(self):
        XCoords=[-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22] ; arrX=DataArrayDouble.New(XCoords)
        YCoords=[0.,0.1,0.37,0.45,0.47,0.49,1.007] ; arrY=DataArrayDouble.New(YCoords)
        mesh=MEDCouplingCMesh.New("My2D_CMesh")
        mesh.setCoords(arrX,arrY)
# ! [PySnippetFieldDoubleBuild3_1]
        fieldOnCells=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
        fieldOnCells.setName("MyTensorFieldOnCellNoTime")
        fieldOnCells.setTimeUnit("ms") # Time unit is ms.
        fieldOnCells.setTime(4.22,2,-1) # Time attached is 4.22 ms, iteration id is 2 and order id (or sub iteration id) is -1
        fieldOnCells.setMesh(mesh)
        array=DataArrayDouble.New()
        array.alloc(fieldOnCells.getMesh().getNumberOfCells(),2) # Implicitely fieldOnCells will be a 2 components field.
        array.fillWithValue(7.)
        fieldOnCells.setArray(array)
        # fieldOnCells is now usable
        # ...
# ! [PySnippetFieldDoubleBuild3_1]
        pass

    def testExampleFieldDoubleBuild4(self):
        XCoords=[-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22] ; arrX=DataArrayDouble.New(XCoords)
        YCoords=[0.,0.1,0.37,0.45,0.47,0.49,1.007] ; arrY=DataArrayDouble.New(YCoords)
        mesh=MEDCouplingCMesh.New("My2D_CMesh")
        mesh.setCoords(arrX,arrY)
# ! [PySnippetFieldDoubleBuild4_1]
        fieldOnNodes=MEDCouplingFieldDouble.New(ON_NODES,CONST_ON_TIME_INTERVAL)
        fieldOnNodes.setName("MyVecFieldOnNodeWithConstTime")
        fieldOnNodes.setTimeUnit("ms") # Time unit is ms.
        fieldOnNodes.setStartTime(4.22,2,-1)
        fieldOnNodes.setEndTime(6.44,4,-1)# fieldOnNodes is defined in interval [4.22 ms,6.44 ms]
        fieldOnNodes.setMesh(mesh)
        array=DataArrayDouble.New()
        array.alloc(fieldOnNodes.getMesh().getNumberOfNodes(),3) # Implicitely fieldOnNodes will be a 3 components field.
        array.fillWithValue(7.)
        fieldOnNodes.setArray(array)
        # fieldOnNodes is now usable
        # ...
# ! [PySnippetFieldDoubleBuild4_1]
        pass

    def testExampleDataArrayApplyFunc1(self):
# ! [PySnippetDataArrayApplyFunc1_1]
        d=DataArrayDouble.New([1.,2.,11.,12.,21.,22.,31.,41.],4,2)
        self.assertRaises(InterpKernelException,d.applyFunc,"x*y")
# ! [PySnippetDataArrayApplyFunc1_1]
# ! [PySnippetDataArrayApplyFunc1_2]
        d=DataArrayDouble.New([1.,2.,11.,12.,21.,22.,31.,41.],4,2)
        d1=d.applyFunc("smth*smth")
        self.assertTrue(d1.isEqual(DataArrayDouble([1.,4.,121.,144.,441.,484.,961.,1681.],4,2),1e-12))
# ! [PySnippetDataArrayApplyFunc1_2]
# ! [PySnippetDataArrayApplyFunc1_3]
        d2=d.applyFunc("smth*IVec+2*smth*JVec")
        self.assertTrue(d2.isEqual(DataArrayDouble([1.,4.,11.,24.,21.,44.,31.,82.],4,2),1e-12))
# ! [PySnippetDataArrayApplyFunc1_3]
# ! [PySnippetDataArrayApplyFunc1_4]
        dd=DataArrayDouble.New([1.,4.,3.,11.,144.,13.,21.,484.,23.,31.,1024.,33.],4,3)
# ! [PySnippetDataArrayApplyFunc1_4]
# ! [PySnippetDataArrayApplyFunc1_5]
        dd1=dd.applyFunc(1,"f+sqrt(g)+h")
        self.assertTrue(dd1.isEqual(DataArrayDouble([6.,36.,66.,96.],4,1),1e-12))
# ! [PySnippetDataArrayApplyFunc1_5]
# ! [PySnippetDataArrayApplyFunc1_6]
        dd2=dd.applyFunc(1,"a+0.*b+c")
        self.assertTrue(dd2.isEqual(DataArrayDouble([4.,24.,44.,64.],4,1),1e-12))
# ! [PySnippetDataArrayApplyFunc1_6]
# ! [PySnippetDataArrayApplyFunc1_7]
        ddd=DataArrayDouble.New([1.,4.,3.,11.,144.,13.,21.,484.,23.,31.,1024.,33.],4,3)
        ddd.setInfoOnComponents(["Y [m]","AA [m/s]","GG [MW]"])
# ! [PySnippetDataArrayApplyFunc1_7]
# ! [PySnippetDataArrayApplyFunc1_8]
        ddd1=ddd.applyFunc2(1,"Y+GG")
        self.assertTrue(ddd1.isEqual(DataArrayDouble([4.,24.,44.,64.],4,1),1e-12))
# ! [PySnippetDataArrayApplyFunc1_8]
# ! [PySnippetDataArrayApplyFunc1_9]
        ddd1=ddd.applyFunc3(1,["X","Y","Z"],"X+Z")
        self.assertTrue(ddd1.isEqual(DataArrayDouble([4.,24.,44.,64.],4,1),1e-12))
# ! [PySnippetDataArrayApplyFunc1_9]
        pass

    pass

unittest.main()
