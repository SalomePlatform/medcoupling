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
from MEDCouplingDataForTest import MEDCouplingDataForTest

if MEDCouplingHasNumPyBindings():
    from numpy import *
    pass

from platform import architecture
from sys import getrefcount

import os, gc, weakref, unittest
import sys
if sys.version_info.major < 3:
  import cPickle as pickle
else:
  import pickle

class MEDCouplingPickleTest(unittest.TestCase):
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test1(self):
        """ Test of a simple DataArrayDouble."""
        x=DataArrayDouble(10,1) ; x.iota() ; x.rearrange(2) ; x.setInfoOnComponents(["aa","bbb"])
        x.setName("toto")
        pickled=pickle.dumps(x,pickle.HIGHEST_PROTOCOL)
        xx=pickle.loads(pickled)
        self.assertTrue(xx.isEqual(x,1e-16))
        # Bigger to check that the behavior is OK for large strings.
        x=DataArrayDouble(1200) ; x.iota() ; x.setInfoOnComponents(["aa"])
        x.setName("titi")
        pickled=pickle.dumps(x,pickle.HIGHEST_PROTOCOL)
        xx=pickle.loads(pickled)
        self.assertTrue(xx.isEqual(x,1e-16))
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test2(self):
        """ Test of a simple DataArrayInt."""
        x=DataArrayInt(10) ; x.iota() ; x.rearrange(2) ; x.setInfoOnComponents(["aa","bbb"])
        x.setName("toto")
        pickled=pickle.dumps(x,pickle.HIGHEST_PROTOCOL)
        xx=pickle.loads(pickled)
        self.assertTrue(xx.isEqual(x))
        # Bigger to check that the behavior is OK for large strings.
        x=DataArrayInt(1200) ; x.iota() ; x.setInfoOnComponents(["aa"])
        x.setName("titi")
        pickled=pickle.dumps(x,pickle.HIGHEST_PROTOCOL)
        xx=pickle.loads(pickled)
        self.assertTrue(xx.isEqual(x))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test3(self):
        """ Test of a MEDCouplingUMesh pickeling."""
        arr=DataArrayDouble(10) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr,arr)
        m=m.buildUnstructured()
        m.setName("mesh")
        m.getCoords().setInfoOnComponents(["aa","bbb","ddddd"])
        m.checkConsistencyLight()
        st=pickle.dumps(m,pickle.HIGHEST_PROTOCOL)
        m2=pickle.loads(st)
        self.assertTrue(m2.isEqual(m,1e-16))
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test4(self):
        """ Idem test3 except that here serialization/deserialization is done explicitely."""
        arr=DataArrayDouble(10) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr,arr)
        m=m.buildUnstructured()
        m.setName("mesh")
        m.getCoords().setInfoOnComponents(["aa","bbb","ddddd"])
        m.checkConsistencyLight()
        #
        a0,a1,a2=m.getTinySerializationInformation()
        b0,b1=m.serialize()
        m2=MEDCouplingUMesh()
        m2.unserialization(a0,a1,b0,b1,a2);
        self.assertTrue(m2.isEqual(m,1e-16))
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test5(self):
        """ Test of a MEDCouplingCMesh pickeling."""
        arrX=DataArrayDouble(10) ; arrX.iota() ; arrX.setInfoOnComponents(["aa"])
        arrY=DataArrayDouble(5) ; arrY.iota() ; arrY.setInfoOnComponents(["bbb"])
        arrZ=DataArrayDouble(7) ; arrZ.iota() ; arrZ.setInfoOnComponents(["cccc"])
        #
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY,arrZ)
        m.setName("mesh")
        m.checkConsistencyLight()
        st=pickle.dumps(m,pickle.HIGHEST_PROTOCOL)
        m2=pickle.loads(st)
        self.assertTrue(m2.isEqual(m,1e-16))
        self.assertTrue(m2.getCoordsAt(0).isEqual(arrX,1e-16))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test6(self):
        """ Test of a MEDCoupling1SGTUMesh pickeling."""
        arr=DataArrayDouble(10) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr)
        m=m.build1SGTUnstructured()
        self.assertTrue(isinstance(m,MEDCoupling1SGTUMesh))
        st=pickle.dumps(m,pickle.HIGHEST_PROTOCOL)
        m2=pickle.loads(st)
        self.assertTrue(m2.isEqual(m,1e-16))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test7(self):
        """ Test of a MEDCoupling1DGTUMesh pickeling."""
        arr=DataArrayDouble(10) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr)
        m=m.buildUnstructured() ; m.convertAllToPoly()
        m=MEDCoupling1DGTUMesh(m)
        self.assertTrue(isinstance(m,MEDCoupling1DGTUMesh))
        st=pickle.dumps(m,pickle.HIGHEST_PROTOCOL)
        m2=pickle.loads(st)
        self.assertTrue(m2.isEqual(m,1e-16))
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test8(self):
        """ Test of a MEDCouplingMappedExtrudedMesh pickeling."""
        arrX=DataArrayDouble(10) ; arrX.iota() ; arrX.setInfoOnComponents(["aa"])
        arrY=DataArrayDouble(5) ; arrY.iota() ; arrY.setInfoOnComponents(["bbb"])
        arrZ=DataArrayDouble(7) ; arrZ.iota() ; arrZ.setInfoOnComponents(["cccc"])
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY,arrZ)
        mesh3D=m.buildUnstructured() ; del m
        #
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY)
        mesh2D=m.buildUnstructured() ; del m
        #
        mesh2D.setCoords(mesh3D.getCoords())
        mesh=MEDCouplingMappedExtrudedMesh(mesh3D,mesh2D,0) ; del mesh3D,mesh2D
        self.assertTrue(isinstance(mesh,MEDCouplingMappedExtrudedMesh))
        st=pickle.dumps(mesh,pickle.HIGHEST_PROTOCOL)
        m2=pickle.loads(st)
        self.assertTrue(m2.isEqual(mesh,1e-16))
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test9(self):
        """ Test of a MEDCouplingCurveLinearMesh pickeling."""
        arrX=DataArrayDouble(10) ; arrX.iota() ; arrX.setInfoOnComponents(["aa"])
        arrY=DataArrayDouble(5) ; arrY.iota() ; arrY.setInfoOnComponents(["bbb"])
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY)
        m=m.buildUnstructured()
        #
        mesh=MEDCouplingCurveLinearMesh() ; mesh.setCoords(m.getCoords()) ; del m
        mesh.setNodeGridStructure([10,5])
        st=pickle.dumps(mesh,pickle.HIGHEST_PROTOCOL)
        m2=pickle.loads(st)
        self.assertTrue(m2.isEqual(mesh,1e-16))
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test10(self):
        """ Test of a MEDCouplingIMesh pickeling."""
        m=MEDCouplingIMesh("mesh",3,DataArrayInt([3,1,4]),DataArrayDouble([1.5,2.5,3.5]),DataArrayDouble((0.5,1.,0.25))) ; m.setAxisUnit("km")
        m.checkConsistencyLight()
        st=pickle.dumps(m,pickle.HIGHEST_PROTOCOL)
        m2=pickle.loads(st)
        self.assertTrue(m2.isEqual(m,1e-16))
        self.assertEqual(m2.getName(),m.getName())
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test11(self):
        """  Test of MEDCouplingFieldDouble lying on MEDCouplingCMesh pickeling. """
        arrX=DataArrayDouble(10) ; arrX.iota() ; arrX.setInfoOnComponents(["aa"])
        arrY=DataArrayDouble(5) ; arrY.iota() ; arrY.setInfoOnComponents(["bbb"])
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY)
        f=m.getMeasureField(True)
        f.setName("aname")
        a=f.getArray()
        b=a[:] ; b.iota(7000.)
        f.setArray(DataArrayDouble.Meld(a,b))
        f.getArray().setInfoOnComponents(["u1","vv2"])
        f.checkConsistencyLight();
        #
        st=pickle.dumps(f,pickle.HIGHEST_PROTOCOL)
        f2=pickle.loads(st)
        self.assertTrue(f2.isEqual(f,1e-16,1e-16))
        self.assertTrue(f2.getMesh().isEqual(f.getMesh(),1e-16))
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test12(self):
        """  Test of MEDCouplingFieldDouble on Gauss Points lying on MEDCouplingUMesh pickeling."""
        _a=0.446948490915965;
        _b=0.091576213509771;
        _p1=0.11169079483905;
        _p2=0.0549758718227661;
        refCoo1=[ 0.,0., 1.,0., 0.,1. ]
        gsCoo1=[ 2*_b-1, 1-4*_b, 2*_b-1, 2.07*_b-1, 1-4*_b,
                 2*_b-1, 1-4*_a, 2*_a-1, 2*_a-1, 1-4*_a, 2*_a-1, 2*_a-1 ]
        wg1=[ 4*_p2, 4*_p2, 4*_p2, 4*_p1, 4*_p1, 4*_p1 ]
        _refCoo1=refCoo1
        _gsCoo1=gsCoo1
        _wg1=wg1
        #
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,NO_TIME);
        f.setMesh(m);
        self.assertEqual(5,f.getNumberOfMeshPlacesExpected());
        self.assertEqual(0,f.getNbOfGaussLocalization());
        f.setGaussLocalizationOnType(NORM_TRI3,_refCoo1,_gsCoo1,_wg1);
        f.setGaussLocalizationOnType(NORM_TRI3,_refCoo1,_gsCoo1,_wg1); # not a bug only to check that it works well
        self.assertRaises(InterpKernelException,f.setGaussLocalizationOnType,NORM_QUAD4,_refCoo1,_gsCoo1,_wg1)
        self.assertEqual(1,f.getNbOfGaussLocalization());
        refCoo2=[ 0.,0., 1.,0., 1.,1., 0.,1. ]
        _refCoo2=refCoo2
        _gsCoo1=_gsCoo1[0:4]
        _wg1=_wg1[0:2]
        f.setGaussLocalizationOnType(NORM_QUAD4,_refCoo2,_gsCoo1,_wg1);
        self.assertEqual(2,f.getNbOfGaussLocalization());
        array=DataArrayDouble.New();
        ptr=18*2*[None]
        for i in range(18 * 2):
            ptr[i]=float(i+1)
        array.setValues(ptr,18,2);
        ptr=array.getPointer();
        f.setArray(array);
        f.setName("MyFirstFieldOnGaussPoint");
        f.checkConsistencyLight();
        self.assertAlmostEqual(27.,f.getIJK(2,5,0),14);
        self.assertAlmostEqual(16.,f.getIJK(1,5,1),14);
        #
        f.clearGaussLocalizations();
        self.assertEqual(0,f.getNbOfGaussLocalization());
        self.assertRaises(InterpKernelException,f.checkConsistencyLight);
        ids1=[0,1,3,4]
        self.assertRaises(InterpKernelException,f.setGaussLocalizationOnCells,ids1,_refCoo2,_gsCoo1,_wg1);
        self.assertEqual(0,f.getNbOfGaussLocalization());
        ids2=[0,4]
        f.setGaussLocalizationOnCells(ids2,_refCoo2,_gsCoo1,_wg1);
        self.assertEqual(1,f.getNbOfGaussLocalization());
        self.assertEqual(0,f.getGaussLocalizationIdOfOneCell(0));
        self.assertRaises(InterpKernelException,f.getGaussLocalizationIdOfOneCell,1);
        ids3=[1,2]
        f.setGaussLocalizationOnCells(ids3,_refCoo1,_gsCoo1,_wg1);
        self.assertEqual(2,f.getNbOfGaussLocalization());
        self.assertEqual(0,f.getGaussLocalizationIdOfOneCell(0));
        self.assertEqual(1,f.getGaussLocalizationIdOfOneCell(1));
        self.assertEqual(1,f.getGaussLocalizationIdOfOneCell(2));
        self.assertRaises(InterpKernelException,f.checkConsistencyLight);#<- cell 3 has no localization
        ids4=[3]
        _gsCoo2=_gsCoo1;
        _wg2=_wg1;
        _gsCoo2[0]=0.8888777776666;
        _wg2[0]=0.1234567892377;
        f.setGaussLocalizationOnCells(ids4,_refCoo2,_gsCoo2,_wg2);
        self.assertEqual(3,f.getNbOfGaussLocalization());
        tmpIds=f.getCellIdsHavingGaussLocalization(0);
        self.assertEqual(ids2,list(tmpIds.getValues()));
        self.assertRaises(InterpKernelException,f.checkConsistencyLight);#<- it's always not ok because undelying array not with the good size.
        array2=f.getArray().subArray(0,10);
        f.setArray(array2);
        f.checkConsistencyLight();
        ####
        st=pickle.dumps(f,pickle.HIGHEST_PROTOCOL)
        f2=pickle.loads(st)
        self.assertTrue(f2.isEqual(f,1e-16,1e-16))
        self.assertTrue(f2.getMesh().isEqual(f.getMesh(),1e-16))
        pass

    def test13(self):
        eStr="This is an exception."
        e=InterpKernelException(eStr)
        self.assertEqual(e.what(),eStr)
        st=pickle.dumps(e,pickle.HIGHEST_PROTOCOL)
        e2=pickle.loads(st)
        self.assertTrue(e is not e2)
        self.assertTrue(isinstance(e2,InterpKernelException))
        self.assertEqual(e2.what(),eStr)
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test14(self):
        """Pickelization of DataArrayBytes"""
        x=DataArrayByte(256,1)
        for i in range(256):
            x[i]=-128+i
            pass
        x.rearrange(2) ; x.setInfoOnComponents(["aa","bbb"])
        x.setName("toto")
        st=pickle.dumps(x,pickle.HIGHEST_PROTOCOL)
        x2=pickle.loads(st)
        self.assertTrue(x2.isEqual(x))
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test15(self):
        """Pickelization of DataArrayFloat"""
        x=DataArrayFloat(256) ; x.iota()
        x.rearrange(2) ; x.setInfoOnComponents(["aa","bbb"])
        x.setName("toto")
        st = pickle.dumps(x, pickle.HIGHEST_PROTOCOL)
        x2 = pickle.loads(st)
        self.assertTrue(x2.isEqual(x,1e-7))
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test16(self):
        """  Test of MEDCouplingFieldInt lying on MEDCouplingCMesh pickeling. """
        arrX=DataArrayDouble(10) ; arrX.iota() ; arrX.setInfoOnComponents(["aa"])
        arrY=DataArrayDouble(5) ; arrY.iota() ; arrY.setInfoOnComponents(["bbb"])
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY)
        f=m.getMeasureField(True)
        f=f.convertToIntField()
        self.assertTrue(isinstance(f,MEDCouplingFieldInt))
        f.setName("aname")
        a=f.getArray()
        b=a[:] ; b.iota(7000)
        f.setArray(DataArrayInt.Meld(a,b))
        f.getArray().setInfoOnComponents(["u1","vv2"])
        f.checkConsistencyLight();
        #
        st = pickle.dumps(f, pickle.HIGHEST_PROTOCOL)
        f2 = pickle.loads(st)
        self.assertTrue(f2.isEqual(f,1e-16,0))
        self.assertTrue(f2.getMesh().isEqual(f.getMesh(),1e-16))
        pass
    
    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def test17(self):
        """  Test of MEDCouplingFieldInt lying on MEDCouplingCMesh pickeling. """
        arrX=DataArrayDouble(10) ; arrX.iota() ; arrX.setInfoOnComponents(["aa"])
        arrY=DataArrayDouble(5) ; arrY.iota() ; arrY.setInfoOnComponents(["bbb"])
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY)
        f2=m.getMeasureField(True)
        f=MEDCouplingFieldFloat(ON_CELLS)
        f.setMesh(m) ; f.setArray(f2.getArray().convertToFloatArr())
        self.assertTrue(isinstance(f,MEDCouplingFieldFloat))
        f.setName("aname")
        a=f.getArray()
        b=a[:] ; b.iota(7000.)
        f.setArray(DataArrayFloat.Meld(a,b))
        f.getArray().setInfoOnComponents(["u1","vv2"])
        f.checkConsistencyLight();
        #
        st = pickle.dumps(f, pickle.HIGHEST_PROTOCOL)
        f2 = pickle.loads(st)
        self.assertTrue(f2.isEqual(f,1e-16,0))
        self.assertTrue(f2.getMesh().isEqual(f.getMesh(),1e-16))
        pass

    def setUp(self):
        pass
    pass

if __name__=="__main__":
    unittest.main()
