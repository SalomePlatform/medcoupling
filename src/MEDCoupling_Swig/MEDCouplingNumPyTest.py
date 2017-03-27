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

from MEDCoupling import *

if MEDCouplingHasNumPyBindings():
    from numpy import *
    pass

from platform import architecture
from sys import getrefcount

import os,gc,weakref,unittest

class MEDCouplingNumPyTest(unittest.TestCase):
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test1(self):
        sz=20
        a=array(0,dtype=int32)
        a.resize(sz)
        a[:]=4
        self.assertEqual(getrefcount(a),2)
        a=a.cumsum(dtype=int32)
        a=array(a,dtype=int64) ; a=array(a,dtype=int32)
        self.assertEqual(getrefcount(a),2)
        d=DataArrayInt(a)
        d[:]=2
        #
        e=DataArrayInt(sz) ; e.fillWithValue(2)
        self.assertTrue(d.isEqual(e))
        #
        a[:]=4 ; e.fillWithValue(4)
        self.assertTrue(d.isEqual(e))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test2(self):
        sz=20
        a=array(0,dtype=int32)
        a.resize(sz,2)
        self.assertEqual(getrefcount(a),2)
        b=a.reshape(2*sz)
        self.assertEqual(getrefcount(a),3)
        self.assertEqual(getrefcount(b),2)
        b[:]=5
        d=DataArrayInt(b)
        #
        e=DataArrayInt(sz*2) ; e.fillWithValue(5)
        self.assertTrue(d.isEqual(e))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test3(self):
        sz=10
        a=array(0,dtype=int32)
        a.resize(sz,2)
        b=a.reshape(2*sz)
        c=a.reshape(2,sz)
        b[:]=6
        b[7:17]=7
        d=DataArrayInt(b)
        self.assertTrue(d.isEqual(DataArrayInt([6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,6,6,6])))
        #
        a=zeros((10,2),dtype=int32)
        b=a.T
        c=b.view()
        a.shape=20
        a[3:]=10.
        d=DataArrayInt(a)
        self.assertTrue(d.isEqual(DataArrayInt([0,0,0,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10])))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test4(self):
        a=zeros(20,dtype=int32)
        b = a[::-1]
        self.assertRaises(InterpKernelException,DataArrayInt.New,b) # b is not contiguous in memory
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test5(self):
        a=arange(20,dtype=int32)
        self.assertEqual(weakref.getweakrefcount(a),0)
        d=DataArrayInt(a)
        self.assertEqual(weakref.getweakrefcount(a),1)
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(d.isIota(20))
        a[:]=2 # modifying a and d because a and d share the same chunk of data
        self.assertTrue(d.isUniform(2))
        del d # d is destroyed, a retrieves its ownership of its initial chunk of data
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called 
        import gc
        gc.collect()
        self.assertTrue(a.flags["OWNDATA"])
        a[:]=4 # a can be used has usual
        self.assertTrue(DataArrayInt(a).isUniform(4))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test6(self):
        a=arange(20,dtype=int32)
        d=DataArrayInt(a) # d owns data of a
        e=DataArrayInt(a) # a not owned -> e only an access to chunk of a 
        self.assertTrue(d.isIota(d.getNumberOfTuples()))
        self.assertTrue(e.isIota(e.getNumberOfTuples()))
        a[:]=6
        self.assertTrue(d.isUniform(6))
        self.assertTrue(e.isUniform(6))
        del a # a destroyed -> d no change because owned and e array is has no more data set
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called 
        import gc
        gc.collect()
        self.assertTrue(d.isUniform(6))
        self.assertTrue(not e.isAllocated())
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test7(self):
        a=array(0,dtype=int32) ; a.resize(10,2)
        b=a.reshape(20)
        c=a.reshape(2,10)
        d=DataArrayInt(b) # d owns data of a
        e=DataArrayInt(b) # a not owned -> e only an access to chunk of a
        f=DataArrayInt(b) # a not owned -> e only an access to chunk of a
        del d # d removed -> a ownes again data
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called 
        import gc
        gc.collect()
        self.assertTrue(e.isUniform(0))
        e[:]=6
        self.assertTrue(e.isUniform(6))
        self.assertTrue(f.isUniform(6))
        self.assertEqual(b.tolist(),[6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6])
        self.assertEqual(a.tolist(),[[6,6],[6,6],[6,6],[6,6],[6,6],[6,6],[6,6],[6,6],[6,6],[6,6]])
        b[:]=arange(20)
        del b # no impact on e and f because a is the base of a.
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        gc.collect()
        self.assertTrue(f.isIota(f.getNumberOfTuples()))
        self.assertTrue(e.isIota(e.getNumberOfTuples()))
        del a # a destroyed, but as c has its base set to a, a exists -> e and f not allocated
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        gc.collect()
        self.assertTrue(f.isIota(f.getNumberOfTuples()))
        self.assertTrue(e.isIota(e.getNumberOfTuples()))
        del c # c killed -> a killed -> e and d are put into not allocated state
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        gc.collect()        
        self.assertTrue(not e.isAllocated())
        self.assertTrue(not f.isAllocated())
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test8(self):
        a=arange(20,dtype=int32)
        self.assertTrue(a.flags["OWNDATA"])
        d=DataArrayInt(a) # d owns data of a
        self.assertTrue(not a.flags["OWNDATA"])
        d.pushBackSilent(20)# d pushBack so release of chunk of data -> a becomes owner of its data again
        self.assertTrue(a.flags["OWNDATA"])
        self.assertTrue(d.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])))
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test9(self):
        sz=20
        a=array(0,dtype=float64)
        a.resize(sz)
        a[:]=4
        self.assertEqual(getrefcount(a),2)
        a=a.cumsum(dtype=float64)
        self.assertEqual(getrefcount(a),2)
        d=DataArrayDouble(a)
        d[:]=2
        #
        e=DataArrayDouble(sz) ; e.fillWithValue(2)
        self.assertTrue(d.isEqual(e,1e-14))
        #
        a[:]=4 ; e.fillWithValue(4)
        self.assertTrue(d.isEqual(e,1e-14))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test10(self):
        sz=20
        a=array(0,dtype=float64)
        a.resize(sz,2)
        self.assertEqual(getrefcount(a),2)
        b=a.reshape(2*sz)
        self.assertEqual(getrefcount(a),3)
        self.assertEqual(getrefcount(b),2)
        b[:]=5
        d=DataArrayDouble(b)
        #
        e=DataArrayDouble(sz*2) ; e.fillWithValue(5)
        self.assertTrue(d.isEqual(e,1e-14))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test11(self):
        sz=10
        a=array(0,dtype=float64)
        a.resize(sz,2)
        b=a.reshape(2*sz)
        c=a.reshape(2,sz)
        b[:]=6
        b[7:17]=7
        d=DataArrayDouble(b)
        self.assertTrue(d.isEqual(DataArrayDouble([6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,6,6,6]),1e-14))
        #
        a=zeros((10,2),dtype=float64)
        b=a.T
        c=b.view()
        a.shape=20
        a[3:]=10.
        d=DataArrayDouble(a)
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,0,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10]),1e-14))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test12(self):
        a=zeros(20,dtype=float64)
        b = a[::-1]
        self.assertRaises(InterpKernelException,DataArrayDouble.New,b) # b is not contiguous in memory
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test13(self):
        a=arange(20,dtype=float64)
        self.assertEqual(weakref.getweakrefcount(a),0)
        d=DataArrayDouble(a)
        self.assertEqual(weakref.getweakrefcount(a),1)
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(d.isEqual(DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]),1e-14))
        self.assertEqual(len(d),20)
        a[:]=2 # modifying a and d because a and d share the same chunk of data
        self.assertTrue(d.isUniform(2,1e-14))
        del d # d is destroyed, a retrieves its ownership of its initial chunk of data
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        import gc
        gc.collect()
        self.assertTrue(a.flags["OWNDATA"])
        a[:]=4 # a can be used has usual
        self.assertTrue(DataArrayDouble(a).isUniform(4,1e-14))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test14(self):
        a=arange(20,dtype=float64)
        d=DataArrayDouble(a) # d owns data of a
        e=DataArrayDouble(a) # a not owned -> e only an access to chunk of a 
        self.assertTrue(d.isEqual(DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]),1e-14))
        self.assertTrue(e.isEqual(DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]),1e-14))
        a[:]=6
        self.assertTrue(d.isUniform(6,1e-14))
        self.assertTrue(e.isUniform(6,1e-14))
        del a # a destroyed -> d no change because owned and e array is has no more data set
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        import gc
        gc.collect()
        self.assertTrue(d.isUniform(6,1e-14))
        self.assertTrue(not e.isAllocated())
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test15(self):
        a=array(0,dtype=float64) ; a.resize(10,2)
        b=a.reshape(20)
        c=a.reshape(2,10)
        d=DataArrayDouble(b) # d owns data of a
        e=DataArrayDouble(b) # a not owned -> e only an access to chunk of a
        f=DataArrayDouble(b) # a not owned -> e only an access to chunk of a
        del d # d removed -> a ownes again data
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        import gc
        gc.collect()
        self.assertTrue(e.isUniform(0,1e-14))
        e[:]=6
        self.assertTrue(e.isUniform(6,1e-14))
        self.assertTrue(f.isUniform(6,1e-14))
        self.assertEqual(b.tolist(),[6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6])
        self.assertEqual(a.tolist(),[[6,6],[6,6],[6,6],[6,6],[6,6],[6,6],[6,6],[6,6],[6,6],[6,6]])
        b[:]=arange(20)
        del b # no impact on e and f because a is the base of a.
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        gc.collect()
        self.assertTrue(f.isEqual(DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]),1e-14))
        self.assertTrue(e.isEqual(DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]),1e-14))
        del a # a destroyed, but as c has its base set to a, a exists -> e and f not allocated
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        gc.collect()
        self.assertTrue(f.isEqual(DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]),1e-14))
        self.assertTrue(e.isEqual(DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]),1e-14))
        del c # c killed -> a killed -> e and d are put into not allocated state
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        gc.collect()
        self.assertTrue(not e.isAllocated())
        self.assertTrue(not f.isAllocated())
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test16(self):
        a=arange(20,dtype=float64)
        self.assertTrue(a.flags["OWNDATA"])
        d=DataArrayDouble(a) # d owns data of a
        self.assertTrue(not a.flags["OWNDATA"])
        d.pushBackSilent(20)# d pushBack so release of chunk of data -> a becomes owner of its data again
        self.assertTrue(a.flags["OWNDATA"])
        self.assertTrue(d.isEqual(DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]),1e-14))
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test17(self):
        d=DataArrayInt.Range(0,20,1)
        a=d.toNumPyArray()
        self.assertTrue(not a.flags["OWNDATA"])
        a[-2:]=100
        self.assertTrue(d.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,100,100])))
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,100,100])
        del a
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        import gc
        gc.collect()
        self.assertTrue(d.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,100,100])))
        #
        d.rearrange(2)
        a=d.toNumPyArray()
        self.assertTrue(not a.flags["OWNDATA"])
        a[-2:]=200
        self.assertTrue(d.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,200,200,200,200],10,2)))
        self.assertEqual(a.tolist(),[[0,1],[2,3],[4,5],[6,7],[8,9],[10,11],[12,13],[14,15],[200,200],[200,200]])
        del a
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        gc.collect()
        self.assertTrue(d.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,200,200,200,200],10,2)))
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test18(self):
        d=DataArrayInt.Range(0,20,1)
        d=d.convertToDblArr()
        a=d.toNumPyArray()
        self.assertTrue(not a.flags["OWNDATA"])
        a[-2:]=100
        self.assertTrue(d.isEqual(DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,100,100]),1e-14))
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,100,100])
        del a
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        import gc
        gc.collect()
        self.assertTrue(d.isEqual(DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,100,100]),1e-14))
        #
        d.rearrange(2)
        a=d.toNumPyArray()
        self.assertTrue(not a.flags["OWNDATA"])
        a[-2:]=200
        self.assertTrue(d.isEqual(DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,200,200,200,200],10,2),1e-14))
        self.assertEqual(a.tolist(),[[0,1],[2,3],[4,5],[6,7],[8,9],[10,11],[12,13],[14,15],[200,200],[200,200]])
        del a
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        gc.collect()
        self.assertTrue(d.isEqual(DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,200,200,200,200],10,2),1e-14))
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test19(self):
        sz=20
        a=array(0,dtype=int32)
        a.resize(sz/2,2)
        a[:]=4
        self.assertEqual(getrefcount(a),2)
        d=DataArrayInt(a)
        self.assertEqual(10,d.getNumberOfTuples())
        self.assertEqual(2,d.getNumberOfComponents())
        self.assertEqual(sz,d.getNbOfElems())
        self.assertTrue(d.isEqual(DataArrayInt([(4,4),(4,4),(4,4),(4,4),(4,4),(4,4),(4,4),(4,4),(4,4),(4,4)])))
        a[:]=7
        self.assertTrue(d.isEqual(DataArrayInt([(7,7),(7,7),(7,7),(7,7),(7,7),(7,7),(7,7),(7,7),(7,7),(7,7)])))
        #
        b=a.reshape((2,5,2))
        self.assertRaises(InterpKernelException,DataArrayInt.New,b) # b has not dimension in [0,1] !
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test20(self):
        sz=20
        a=array(0,dtype=float64)
        a.resize(sz/2,2)
        a[:]=4
        self.assertEqual(getrefcount(a),2)
        d=DataArrayDouble(a)
        self.assertEqual(10,d.getNumberOfTuples())
        self.assertEqual(2,d.getNumberOfComponents())
        self.assertEqual(sz,d.getNbOfElems())
        self.assertTrue(d.isEqual(DataArrayDouble([(4.,4.),(4.,4.),(4.,4.),(4.,4.),(4.,4.),(4.,4.),(4.,4.),(4.,4.),(4.,4.),(4.,4.)]),1e-14))
        a[:]=7
        self.assertTrue(d.isEqual(DataArrayDouble([(7.,7.),(7.,7.),(7.,7.),(7.,7.),(7.,7.),(7.,7.),(7.,7.),(7.,7.),(7.,7.),(7.,7.)]),1e-14))
        #
        b=a.reshape((2,5,2))
        self.assertRaises(InterpKernelException,DataArrayDouble.New,b) # b has not dimension in [0,1] !
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test21(self):
        #tests that only DataArray*(npArray) contructor is available
        a=array(0,dtype=int32)
        a.resize(20)
        DataArrayInt(a)
        self.assertRaises(InterpKernelException,DataArrayInt.New,a,20)
        self.assertRaises(InterpKernelException,DataArrayInt.New,a,20,1)
        a=array(0,dtype=float64)
        a.resize(20)
        DataArrayDouble(a)
        self.assertRaises(InterpKernelException,DataArrayDouble.New,a,20)
        self.assertRaises(InterpKernelException,DataArrayDouble.New,a,20,1)
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test22(self):
        d=DataArrayDouble(10)
        d.iota()
        a=d.toNumPyArray()
        self.assertTrue(not a.flags["OWNDATA"])
        del d
        gc.collect()
        self.assertTrue(a.flags["OWNDATA"])
        self.assertEqual(a.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        #
        d=DataArrayInt(10)
        d.iota()
        a=d.toNumPyArray()
        self.assertTrue(not a.flags["OWNDATA"])
        del d
        gc.collect()
        self.assertTrue(a.flags["OWNDATA"])
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9])
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test23(self):
        d=DataArrayDouble(10)
        d.iota()
        a=d.toNumPyArray()
        b=d.toNumPyArray()
        c=d.toNumPyArray()
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(not b.flags["OWNDATA"])
        self.assertTrue(not c.flags["OWNDATA"])
        self.assertEqual(a.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertEqual(b.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertEqual(c.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        del d
        gc.collect()
        self.assertTrue(a.flags["OWNDATA"])
        self.assertTrue(not b.flags["OWNDATA"])
        self.assertTrue(not c.flags["OWNDATA"])
        self.assertEqual(a.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertEqual(b.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertEqual(c.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        #
        d=DataArrayInt(10)
        d.iota()
        a=d.toNumPyArray()
        b=d.toNumPyArray()
        c=d.toNumPyArray()
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(not b.flags["OWNDATA"])
        self.assertTrue(not c.flags["OWNDATA"])
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertEqual(b.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertEqual(c.tolist(),[0,1,2,3,4,5,6,7,8,9])
        del d
        gc.collect()
        self.assertTrue(a.flags["OWNDATA"])
        self.assertTrue(not b.flags["OWNDATA"])
        self.assertTrue(not c.flags["OWNDATA"])
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertEqual(b.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertEqual(c.tolist(),[0,1,2,3,4,5,6,7,8,9])
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test24(self):
        d=DataArrayDouble(10)
        d.iota()
        a=d.toNumPyArray()
        self.assertEqual(a.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(a.base is None)
        del a
        gc.collect()
        a=d.toNumPyArray()
        self.assertEqual(a.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(a.base is None)
        b=d.toNumPyArray()
        self.assertEqual(a.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertEqual(b.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(not b.flags["OWNDATA"])
        self.assertTrue(b.base is a)
        del a
        gc.collect()
        self.assertEqual(b.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertTrue(not b.flags["OWNDATA"])
        del d
        gc.collect()
        self.assertEqual(b.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertTrue(not b.flags["OWNDATA"])
        #
        d=DataArrayInt(10)
        d.iota()
        a=d.toNumPyArray()
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(a.base is None)
        del a
        gc.collect()
        a=d.toNumPyArray()
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(a.base is None)
        b=d.toNumPyArray()
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertEqual(b.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(not b.flags["OWNDATA"])
        self.assertTrue(b.base is a)
        del a
        gc.collect()
        self.assertEqual(b.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertTrue(not b.flags["OWNDATA"])
        del d
        gc.collect()
        self.assertEqual(b.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertTrue(not b.flags["OWNDATA"])
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test25(self):
        a=arange(10,dtype=int32)
        b=DataArrayInt(a)
        c=DataArrayInt(a)
        d=DataArrayInt(a)
        self.assertTrue(b.isIota(10))
        self.assertTrue(c.isIota(10))
        self.assertTrue(d.isIota(10))
        c.pushBackSilent(10) # c and a,b are dissociated
        self.assertTrue(b.isIota(10))
        self.assertTrue(c.isIota(11))
        self.assertTrue(d.isIota(10))
        del a
        gc.collect()
        self.assertTrue(b.isIota(10))
        self.assertTrue(c.isIota(11))
        self.assertTrue(not d.isAllocated())
        del b
        gc.collect()
        self.assertTrue(c.isIota(11))
        #
        a=arange(10,dtype=int32)
        b=DataArrayInt(a)
        c=DataArrayInt(a)
        self.assertTrue(b.isIota(10))
        self.assertTrue(c.isIota(10))
        b.pushBackSilent(10) # c and a,b are dissociated
        self.assertTrue(b.isIota(11))
        self.assertTrue(c.isIota(10))
        del a
        gc.collect()
        self.assertTrue(b.isIota(11))
        self.assertTrue(not c.isAllocated())
        del b
        gc.collect()
        self.assertTrue(not c.isAllocated())
        #
        a=float64(arange(5,dtype=int32))
        b=DataArrayDouble(a)
        c=DataArrayDouble(a)
        d=DataArrayDouble(a)
        self.assertTrue(b.isEqual(DataArrayDouble([0.,1.,2.,3.,4.]),1e-12))
        self.assertTrue(c.isEqual(DataArrayDouble([0.,1.,2.,3.,4.]),1e-12))
        self.assertTrue(d.isEqual(DataArrayDouble([0.,1.,2.,3.,4.]),1e-12))
        c.pushBackSilent(10.) # c and a,b are dissociated
        self.assertTrue(b.isEqual(DataArrayDouble([0.,1.,2.,3.,4.]),1e-12))
        self.assertTrue(c.isEqual(DataArrayDouble([0.,1.,2.,3.,4.,10.]),1e-12))
        self.assertTrue(d.isEqual(DataArrayDouble([0.,1.,2.,3.,4.]),1e-12))
        del a
        gc.collect()
        self.assertTrue(b.isEqual(DataArrayDouble([0.,1.,2.,3.,4.]),1e-12))
        self.assertTrue(c.isEqual(DataArrayDouble([0.,1.,2.,3.,4.,10.]),1e-12))
        self.assertTrue(not d.isAllocated())
        del b
        gc.collect()
        self.assertTrue(c.isEqual(DataArrayDouble([0.,1.,2.,3.,4.,10.]),1e-12))
        #
        a=float64(arange(5,dtype=int32))
        b=DataArrayDouble(a)
        c=DataArrayDouble(a)
        self.assertTrue(b.isEqual(DataArrayDouble([0.,1.,2.,3.,4.]),1e-12))
        self.assertTrue(c.isEqual(DataArrayDouble([0.,1.,2.,3.,4.]),1e-12))
        b.pushBackSilent(10.) # c and a,b are dissociated
        self.assertTrue(b.isEqual(DataArrayDouble([0.,1.,2.,3.,4.,10.]),1e-12))
        self.assertTrue(c.isEqual(DataArrayDouble([0.,1.,2.,3.,4.]),1e-12))
        del a
        gc.collect()
        self.assertTrue(b.isEqual(DataArrayDouble([0.,1.,2.,3.,4.,10.]),1e-12))
        self.assertTrue(not c.isAllocated())
        del b
        gc.collect()
        self.assertTrue(not c.isAllocated())
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test26(self):
        d=DataArrayInt(15) ; d.iota()
        d.rearrange(3)
        a=d.toNumPyArray()
        self.assertEqual(a.ndim,2)
        self.assertEqual(a.size,15)
        self.assertEqual(a.shape,(5,3))
        self.assertEqual(a.strides,(12,4))
        self.assertEqual(a.nbytes,60)
        self.assertEqual(a.itemsize,4)
        self.assertEqual(a.tolist(),[[0,1,2],[3,4,5],[6,7,8],[9,10,11],[12,13,14]])
        #
        d2=d.convertToDblArr()
        a2=d2.toNumPyArray()
        self.assertEqual(a2.ndim,2)
        self.assertEqual(a2.size,15)
        self.assertEqual(a2.shape,(5,3))
        self.assertEqual(a2.strides,(24,8))
        self.assertEqual(a2.nbytes,120)
        self.assertEqual(a2.itemsize,8)
        self.assertEqual(a2.tolist(),[[0.,1.,2.],[3.,4.,5.],[6.,7.,8.],[9.,10.,11.],[12.,13.,14.]])
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test27(self):
        m0=DenseMatrix(DataArrayDouble([2,3,4,5,1,6]),2,3)
        m0np=m0.toNumPyMatrix()
        self.assertEqual(m0np.shape,(2,3))
        self.assertEqual(m0np.tolist(),[[2.0,3.0,4.0],[5.0,1.0,6.0]])
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test28(self):
        """Test on DataArrayBytes"""
        # use case 1
        d=DataArrayByte(256)
        for i in range(len(d)):
            d[i]=-128+i
            pass
        arr=d.toNumPyArray()
        for i in range(len(d)):
            self.assertEqual(int(arr[i]),-128+i)
            pass
        d[0]=7
        self.assertEqual(int(arr[0]),7)
        arr[0]=8
        self.assertEqual(int(d.getIJ(0,0)),8)
        del arr
        gc.collect()
        del d
        gc.collect()
        # use case 2
        d=DataArrayByte(256)
        for i in range(len(d)):
            d[i]=-128+i
            pass
        arr=d.toNumPyArray()
        for i in range(len(d)):
            self.assertEqual(int(arr[i]),-128+i)
            pass
        del d
        gc.collect()
        del arr
        gc.collect()
        # use case 3
        d=DataArrayByte(256)
        for i in range(len(d)):
            d[i]=-128+i
            pass
        arr1=d.toNumPyArray()
        arr2=d.toNumPyArray()
        arr3=d.toNumPyArray()
        d[0]=10
        self.assertEqual(int(arr1[0]),10) ; self.assertEqual(int(arr2[0]),10) ; self.assertEqual(int(arr3[0]),10)
        arr2[0]=15 ; self.assertEqual(int(d.getIJ(0,0)),15) ; self.assertEqual(int(arr1[0]),15) ; self.assertEqual(int(arr3[0]),15)
        arr1[0]=-128
        for i in range(len(d)):
            self.assertEqual(int(arr1[i]),-128+i)
            self.assertEqual(int(arr2[i]),-128+i)
            self.assertEqual(int(arr3[i]),-128+i)
            pass
        del arr2
        gc.collect()
        for i in range(len(d)):
            self.assertEqual(int(arr1[i]),-128+i)
            self.assertEqual(int(arr3[i]),-128+i)
            pass
        del arr1
        gc.collect()
        for i in range(len(d)):
            self.assertEqual(int(arr3[i]),-128+i)
            pass
        del arr3
        gc.collect()
        # use case 4
        arr=array(0,dtype=int8)
        arr.resize(256)
        for i in range(256):
            arr[i]=-128+i
            pass
        d=DataArrayByte(arr)
        for i in range(256):
            self.assertEqual(int(d.getIJ(i,0)),-128+i)
            pass
        del arr
        gc.collect()
        del d
        gc.collect()
        # use case 5
        arr=array(0,dtype=int8)
        arr.resize(256)
        for i in range(256):
            arr[i]=-128+i
            pass
        d=DataArrayByte(arr)
        for i in range(256):
            self.assertEqual(int(d.getIJ(i,0)),-128+i)
            pass
        del d
        gc.collect()
        del arr
        gc.collect()
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test29(self):
        """Same as test9 with float32"""
        sz=20
        a=array(0,dtype=float32)
        a.resize(sz)
        a[:]=4
        self.assertEqual(getrefcount(a),2)
        a=a.cumsum(dtype=float32)
        self.assertEqual(getrefcount(a),2)
        d=DataArrayFloat(a)
        d[:]=2
        #
        e=DataArrayFloat(sz) ; e.fillWithValue(2)
        self.assertTrue(d.isEqual(e,1e-7))
        #
        a[:]=4 ; e.fillWithValue(4)
        self.assertTrue(d.isEqual(e,1e-7))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test30(self):
        """Same as test10 with float32"""
        sz=20
        a=array(0,dtype=float32)
        a.resize(sz,2)
        self.assertEqual(getrefcount(a),2)
        b=a.reshape(2*sz)
        self.assertEqual(getrefcount(a),3)
        self.assertEqual(getrefcount(b),2)
        b[:]=5
        d=DataArrayFloat(b)
        #
        e=DataArrayFloat(sz*2) ; e.fillWithValue(5)
        self.assertTrue(d.isEqual(e,1e-7))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test31(self):
        """Same as test11 with float32"""
        sz=10
        a=array(0,dtype=float32)
        a.resize(sz,2)
        b=a.reshape(2*sz)
        c=a.reshape(2,sz)
        b[:]=6
        b[7:17]=7
        d=DataArrayFloat(b)
        self.assertTrue(d.isEqual(DataArrayFloat([6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,6,6,6]),1e-7))
        #
        a=zeros((10,2),dtype=float32)
        b=a.T
        c=b.view()
        a.shape=20
        a[3:]=10.
        d=DataArrayFloat(a)
        self.assertTrue(d.isEqual(DataArrayFloat([0,0,0,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10]),1e-7))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test32(self):
        """Same as test12 with float32"""
        a=zeros(20,dtype=float32)
        b = a[::-1]
        self.assertRaises(InterpKernelException,DataArrayFloat.New,b) # b is not contiguous in memory
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test33(self):
        """Same as test13 with float32"""
        a=arange(20,dtype=float32)
        self.assertEqual(weakref.getweakrefcount(a),0)
        d=DataArrayFloat(a)
        self.assertEqual(weakref.getweakrefcount(a),1)
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(d.isEqual(DataArrayFloat([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]),1e-7))
        self.assertEqual(len(d),20)
        a[:]=2 # modifying a and d because a and d share the same chunk of data
        self.assertTrue(d.isUniform(2,1e-7))
        del d # d is destroyed, a retrieves its ownership of its initial chunk of data
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        import gc
        gc.collect()
        self.assertTrue(a.flags["OWNDATA"])
        a[:]=4 # a can be used has usual
        self.assertTrue(DataArrayFloat(a).isUniform(4,1e-7))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test34(self):
        """Same as test14 with float32"""
        a=arange(20,dtype=float32)
        d=DataArrayFloat(a) # d owns data of a
        e=DataArrayFloat(a) # a not owned -> e only an access to chunk of a 
        self.assertTrue(d.isEqual(DataArrayFloat([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]),1e-7))
        self.assertTrue(e.isEqual(DataArrayFloat([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]),1e-7))
        a[:]=6
        self.assertTrue(d.isUniform(6,1e-7))
        self.assertTrue(e.isUniform(6,1e-7))
        del a # a destroyed -> d no change because owned and e array is has no more data set
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        import gc
        gc.collect()
        self.assertTrue(d.isUniform(6,1e-7))
        self.assertTrue(not e.isAllocated())
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test35(self):
        """Same as test15 with float32"""
        a=array(0,dtype=float32) ; a.resize(10,2)
        b=a.reshape(20)
        c=a.reshape(2,10)
        d=DataArrayFloat(b) # d owns data of a
        e=DataArrayFloat(b) # a not owned -> e only an access to chunk of a
        f=DataArrayFloat(b) # a not owned -> e only an access to chunk of a
        del d # d removed -> a ownes again data
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        import gc
        gc.collect()
        self.assertTrue(e.isUniform(0,1e-7))
        e[:]=6
        self.assertTrue(e.isUniform(6,1e-7))
        self.assertTrue(f.isUniform(6,1e-7))
        self.assertEqual(b.tolist(),[6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6])
        self.assertEqual(a.tolist(),[[6,6],[6,6],[6,6],[6,6],[6,6],[6,6],[6,6],[6,6],[6,6],[6,6]])
        b[:]=arange(20)
        del b # no impact on e and f because a is the base of a.
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        gc.collect()
        self.assertTrue(f.isEqual(DataArrayFloat([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]),1e-7))
        self.assertTrue(e.isEqual(DataArrayFloat([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]),1e-7))
        del a # a destroyed, but as c has its base set to a, a exists -> e and f not allocated
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        gc.collect()
        self.assertTrue(f.isEqual(DataArrayFloat([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]),1e-7))
        self.assertTrue(e.isEqual(DataArrayFloat([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]),1e-7))
        del c # c killed -> a killed -> e and d are put into not allocated state
        ##@@ Ensure a pass of the garbage collector so that the de-allocator of d is called
        gc.collect()
        self.assertTrue(not e.isAllocated())
        self.assertTrue(not f.isAllocated())
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test36(self):
        """Same as test16 with float32"""
        a=arange(20,dtype=float32)
        self.assertTrue(a.flags["OWNDATA"])
        d=DataArrayFloat(a) # d owns data of a
        self.assertTrue(not a.flags["OWNDATA"])
        d.pushBackSilent(20)# d pushBack so release of chunk of data -> a becomes owner of its data again
        self.assertTrue(a.flags["OWNDATA"])
        self.assertTrue(d.isEqual(DataArrayFloat([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]),1e-7))
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test37(self):
        """Same as test20 with float32"""
        sz=20
        a=array(0,dtype=float32)
        a.resize(sz/2,2)
        a[:]=4
        self.assertEqual(getrefcount(a),2)
        d=DataArrayFloat(a)
        self.assertEqual(10,d.getNumberOfTuples())
        self.assertEqual(2,d.getNumberOfComponents())
        self.assertEqual(sz,d.getNbOfElems())
        self.assertTrue(d.isEqual(DataArrayFloat([(4.,4.),(4.,4.),(4.,4.),(4.,4.),(4.,4.),(4.,4.),(4.,4.),(4.,4.),(4.,4.),(4.,4.)]),1e-7))
        a[:]=7
        self.assertTrue(d.isEqual(DataArrayFloat([(7.,7.),(7.,7.),(7.,7.),(7.,7.),(7.,7.),(7.,7.),(7.,7.),(7.,7.),(7.,7.),(7.,7.)]),1e-7))
        #
        b=a.reshape((2,5,2))
        self.assertRaises(InterpKernelException,DataArrayFloat.New,b) # b has not dimension in [0,1] !
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test38(self):
        """Same as test22 with float32"""
        d=DataArrayFloat(10)
        d.iota()
        a=d.toNumPyArray()
        self.assertTrue(not a.flags["OWNDATA"])
        del d
        gc.collect()
        self.assertTrue(a.flags["OWNDATA"])
        self.assertEqual(a.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        #
        d=DataArrayInt(10)
        d.iota()
        a=d.toNumPyArray()
        self.assertTrue(not a.flags["OWNDATA"])
        del d
        gc.collect()
        self.assertTrue(a.flags["OWNDATA"])
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9])
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test39(self):
        """Same as test23 with float32"""
        d=DataArrayFloat(10)
        d.iota()
        a=d.toNumPyArray()
        b=d.toNumPyArray()
        c=d.toNumPyArray()
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(not b.flags["OWNDATA"])
        self.assertTrue(not c.flags["OWNDATA"])
        self.assertEqual(a.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertEqual(b.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertEqual(c.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        del d
        gc.collect()
        self.assertTrue(a.flags["OWNDATA"])
        self.assertTrue(not b.flags["OWNDATA"])
        self.assertTrue(not c.flags["OWNDATA"])
        self.assertEqual(a.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertEqual(b.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertEqual(c.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        #
        d=DataArrayInt(10)
        d.iota()
        a=d.toNumPyArray()
        b=d.toNumPyArray()
        c=d.toNumPyArray()
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(not b.flags["OWNDATA"])
        self.assertTrue(not c.flags["OWNDATA"])
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertEqual(b.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertEqual(c.tolist(),[0,1,2,3,4,5,6,7,8,9])
        del d
        gc.collect()
        self.assertTrue(a.flags["OWNDATA"])
        self.assertTrue(not b.flags["OWNDATA"])
        self.assertTrue(not c.flags["OWNDATA"])
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertEqual(b.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertEqual(c.tolist(),[0,1,2,3,4,5,6,7,8,9])
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test40(self):
        """Same as test24 with float32"""
        d=DataArrayFloat(10)
        d.iota()
        a=d.toNumPyArray()
        self.assertEqual(a.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(a.base is None)
        del a
        gc.collect()
        a=d.toNumPyArray()
        self.assertEqual(a.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(a.base is None)
        b=d.toNumPyArray()
        self.assertEqual(a.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertEqual(b.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(not b.flags["OWNDATA"])
        self.assertTrue(b.base is a)
        del a
        gc.collect()
        self.assertEqual(b.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertTrue(not b.flags["OWNDATA"])
        del d
        gc.collect()
        self.assertEqual(b.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertTrue(not b.flags["OWNDATA"])
        #
        d=DataArrayInt(10)
        d.iota()
        a=d.toNumPyArray()
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(a.base is None)
        del a
        gc.collect()
        a=d.toNumPyArray()
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(a.base is None)
        b=d.toNumPyArray()
        self.assertEqual(a.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertEqual(b.tolist(),[0,1,2,3,4,5,6,7,8,9])
        self.assertTrue(not a.flags["OWNDATA"])
        self.assertTrue(not b.flags["OWNDATA"])
        self.assertTrue(b.base is a)
        del a
        gc.collect()
        self.assertEqual(b.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertTrue(not b.flags["OWNDATA"])
        del d
        gc.collect()
        self.assertEqual(b.tolist(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
        self.assertTrue(not b.flags["OWNDATA"])
        pass

    def setUp(self):
        pass
    pass

#gc.set_debug(gc.DEBUG_LEAK)
unittest.main()
