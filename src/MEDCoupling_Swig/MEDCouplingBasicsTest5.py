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

class MEDCouplingBasicsTest5(unittest.TestCase):
    def testSwig2FieldDoubleBuildSubPartRange1(self):
        #ON_CELLS
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr = DataArrayDouble(5, 2) ; arr[:, 0] = list(range(7, 12)) ; arr[:, 1] = 100 + arr[:, 0]
        f.setArray(arr)
        f.checkConsistencyLight()
        ff=f[1:-1:2]
        ff.checkConsistencyLight()
        self.assertTrue((m.buildPartOfMySelf([1,3],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(2,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[1,3]],1e-12))
        #
        a,b=f.buildSubMeshDataRange(2,5,1)
        self.assertTrue(m.buildPartOfMySelf([2,3,4],True).isEqual(a,1e-12))
        self.assertEqual(b,slice(2,5,1))
        ff=f[2:]
        ff.checkConsistencyLight()
        self.assertTrue((m.buildPartOfMySelf([2,3,4],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[2,3,4]],1e-12))
        #
        ff=f[-2:0:-1]
        ff.checkConsistencyLight()
        self.assertTrue((m.buildPartOfMySelf([3,2,1],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[3,2,1]],1e-12))
        self.assertTrue(f[-2:0:-1,1].getArray().isEqual(arr[[3,2,1],1],1e-12))
        #ON_NODES
        f=MEDCouplingFieldDouble(ON_NODES)
        f.setMesh(m)
        arr = DataArrayDouble(9, 2) ; arr[:, 0] = list(range(7, 16)) ; arr[:, 1] = 100 + arr[:, 0]
        f.setArray(arr)
        f.checkConsistencyLight()
        ff=f[1:-1:2]
        ff.checkConsistencyLight()
        self.assertTrue((m.buildPartOfMySelf([1,3],False)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(6,ff.getMesh().getNumberOfNodes())
        self.assertTrue(2,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[1,2,3,4,6,7]],1e-12))
        #
        m2=m.buildPartRange(2,5,1)
        self.assertTrue(m.buildPartOfMySelf([2,3,4],True).isEqual(m2,1e-12))
        m2,b=m.buildPartRangeAndReduceNodes(2,5,1)
        self.assertTrue(m.buildPartOfMySelf([2,3,4],False).isEqual(m2,1e-12))
        self.assertTrue(b.isEqual(DataArrayInt([-1,-1,0,1,2,3,4,5,6])))
        a,b=f.buildSubMeshDataRange(2,5,1)
        self.assertTrue(m.buildPartOfMySelf([2,3,4],False).isEqual(a,1e-12))
        self.assertTrue(b.isEqual(DataArrayInt([2,3,4,5,6,7,8])))
        ff=f[2:]
        ff.checkConsistencyLight()
        self.assertTrue((m.buildPartOfMySelf([2,3,4],False)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(7,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[2,3,4,5,6,7,8]],1e-12))
        #
        ff=f[-2:0:-1]
        ff.checkConsistencyLight()
        self.assertTrue((m.buildPartOfMySelf([3,2,1],False)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(7,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[1,2,3,4,5,6,7]],1e-12))
        self.assertTrue(f[-2:0:-1,1].getArray().isEqual(arr[[1,2,3,4,5,6,7],1],1e-12))
        #ON_GAUSS_NE
        f=MEDCouplingFieldDouble(ON_GAUSS_NE)
        f.setMesh(m)
        arr = DataArrayDouble(18, 2) ; arr[:, 0] = list(range(7, 25)) ; arr[:, 1] = 100 + arr[:, 0]
        f.setArray(arr)
        f.checkConsistencyLight()
        ff=f[1:-1:2]
        ff.checkConsistencyLight()
        self.assertTrue((m.buildPartOfMySelf([1,3],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(2,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[4,5,6,10,11,12,13]],1e-12))
        #
        a,b=f.buildSubMeshDataRange(2,5,1)
        self.assertTrue(m.buildPartOfMySelf([2,3,4],True).isEqual(a,1e-12))
        self.assertEqual(b,slice(7,18,1))
        ff=f[2:]
        ff.checkConsistencyLight()
        self.assertTrue((m.buildPartOfMySelf([2,3,4],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[7,8,9,10,11,12,13,14,15,16,17]],1e-12))
        #
        ff=f[-2:0:-1]
        ff.checkConsistencyLight()
        self.assertTrue((m.buildPartOfMySelf([3,2,1],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[10,11,12,13,7,8,9,4,5,6]],1e-12))
        self.assertTrue(f[-2:0:-1,1].getArray().isEqual(arr[[10,11,12,13,7,8,9,4,5,6],1],1e-12))
        #ON_GAUSS_PT
        f=MEDCouplingFieldDouble(ON_GAUSS_PT)
        f.setMesh(m)
        f.setGaussLocalizationOnCells([0,4],[0,0,1,0,1,1,1,0],[1.1,1.1,2.2,2.2],[0.2,0.8]);
        f.setGaussLocalizationOnCells([3],[0,0,1,0,1,1,1,0],[1.1,1.1,2.2,2.2,3.,3.],[0.2,0.4,0.4]);
        f.setGaussLocalizationOnCells([1],[0,0,1,0,1,0],[1.1,1.1,2.2,2.2,3.,3.,4.,4.],[0.1,0.1,0.4,0.4]);
        f.setGaussLocalizationOnCells([2],[0,0,1,0,1,0],[1.1,1.1,2.2,2.2,3.,3.,4.,4.,5.,5.],[0.1,0.1,0.4,0.3,0.1]);
        arr = DataArrayDouble(16, 2) ; arr[:, 0] = list(range(7, 23)) ; arr[:, 1] = 100 + arr[:, 0]
        f.setArray(arr)
        f.checkConsistencyLight()
        ff=f[1:-1:2]
        ff.checkConsistencyLight()
        self.assertTrue((m.buildPartOfMySelf([1,3],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(2,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[2,3,4,5,11,12,13]],1e-12))
        #
        a,b=f.buildSubMeshDataRange(2,5,1)
        self.assertTrue(m.buildPartOfMySelf([2,3,4],True).isEqual(a,1e-12))
        self.assertEqual(b,slice(6,16,1))
        ff=f[2:]
        ff.checkConsistencyLight()
        self.assertTrue((m.buildPartOfMySelf([2,3,4],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[6,7,8,9,10,11,12,13,14,15]],1e-12))
        #
        ff=f[-2:0:-1]
        ff.checkConsistencyLight()
        self.assertTrue((m.buildPartOfMySelf([3,2,1],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[11,12,13,6,7,8,9,10,2,3,4,5]],1e-12))
        self.assertTrue(f[-2:0:-1,0].getArray().isEqual(arr[[11,12,13,6,7,8,9,10,2,3,4,5],0],1e-12))
        pass

    def testSwig2FieldDoubleApplyFuncBug1(self):
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(MEDCouplingDataForTest.build2DTargetMesh_1())
        f.applyFunc(3,700.)
        f.checkConsistencyLight()
        self.assertEqual(3,f.getArray().getNumberOfComponents())
        f.getArray().rearrange(1)
        self.assertTrue(f.getArray().isUniform(700.,1e-10))
        f.getArray().rearrange(3)
        f.checkConsistencyLight()
        f.applyFunc(4,800.)
        f.checkConsistencyLight()
        self.assertEqual(4,f.getArray().getNumberOfComponents())
        f.getArray().rearrange(1)
        self.assertTrue(f.getArray().isUniform(800.,1e-10))
        f.getArray().rearrange(4)
        f.checkConsistencyLight()
        pass

    def testSwig2ComputeTupleIdsNearTupleBug1(self):
        coords=[1.1,0.0, 1.1,0.0 ];
        coordsArr=DataArrayDouble(coords,2,2);
        mesh=MEDCouplingUMesh();
        mesh.setCoords(coordsArr);
        points=[1.1, 0.002]
        c,cI=mesh.getNodeIdsNearPoints(points,0.00185);
        self.assertTrue(c.isEqual(DataArrayInt([])))
        self.assertTrue(cI.isEqual(DataArrayInt([0,0])))
        c,cI=mesh.getNodeIdsNearPoints(points,0.00200000000000001);
        self.assertTrue(c.isEqual(DataArrayInt([0,1])))
        self.assertTrue(cI.isEqual(DataArrayInt([0,2])))
        pass

    def testSwig2NonRegressionBugChangeUnderlyingWithZeroCells(self):
        coords1=[0.,1.,2.,3.]
        coords2=[2.,1.,0.,3.] #0 <==> #2
        # mesh 1
        mesh1=MEDCouplingUMesh.New();
        coordsArr=DataArrayDouble.New(coords1,4,1);
        mesh1.setCoords(coordsArr);
        mesh1.setMeshDimension(0);
        mesh1.allocateCells(0);
        mesh1.finishInsertingCells();
        # mesh 2
        mesh2=mesh1.deepCopy();
        coordsArr=DataArrayDouble.New(coords2,4,1);
        mesh2.setCoords(coordsArr);
        field = mesh1.fillFromAnalytic(ON_NODES,1,"x")
        field.checkConsistencyLight()
        levOfCheck = 10
        field.changeUnderlyingMesh( mesh2, levOfCheck, 1e-13, 0 )
        self.assertTrue( field.getArray().getValues() == coords2 )
        pass

    def testSwig2UMeshDistanceToMesh2(self):
        sz=5
        m=MEDCouplingCMesh()
        arr=DataArrayDouble(sz+1) ; arr.iota() ; arr/=sz
        m.setCoords(arr,arr,arr)
        m=m.buildUnstructured()
        m1=m.computeSkin()
        m1.zipCoords()
        c=m1.getCoords()[:]
        d=2*(c-[0.5,0.5,0.5])+[0.5,0.5,0.5]
        time_deb = datetime.now()
        #print "go.."
        a,b=m1.distanceToPoints(d)
        #print 'time spent in distanceToPoints %s ' %str(datetime.now() - time_deb)
        time_deb = datetime.now()
        a1=DataArrayDouble(len(d))
        b1=DataArrayInt(len(d))
        m1s = [m1[i] for i in range(m1.getNumberOfCells())]
        for j,pt in enumerate(d):
            eter=1e308
            fter=-1
            for i,miter in enumerate(m1s):
                e,f=miter.distanceToPoint(pt)
                self.assertEqual(0,f)
                if e<eter:
                    eter=e ; fter=i
                    pass
                pass
            a1[j]=eter
            b1[j]=fter
            pass
        #print 'time spent in naive distanceToPoints  %s ' %str(datetime.now() - time_deb)
        self.assertTrue(a.isEqual(a1,1e-12))
        self.assertTrue(b.isEqual(b1))
        self.assertTrue(a.isEqual(DataArrayDouble([0.8660254037844386,0.714142842854285,0.7071067811865476,0.7071067811865476,0.714142842854285,0.8660254037844386,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706632,0.714142842854285,0.7071067811865475,0.5099019513592785,0.5,0.5,0.5099019513592785,0.7071067811865476,0.7071067811865475,0.5099019513592785,0.5,0.5,0.5099019513592785,0.7071067811865476,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706632,0.714142842854285,0.8660254037844386,0.714142842854285,0.7071067811865476,0.7071067811865476,0.714142842854285,0.8660254037844386,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706631,0.714142842854285,0.5196152422706631,0.5196152422706632,0.5099019513592784,0.5099019513592785,0.5099019513592784,0.5099019513592785,0.5196152422706631,0.5196152422706632,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706632,0.714142842854285,0.7071067811865475,0.5099019513592785,0.5,0.5,0.5099019513592784,0.7071067811865475,0.5099019513592784,0.5099019513592785,0.5,0.5,0.5,0.5,0.5099019513592785,0.5099019513592785,0.7071067811865476,0.5099019513592785,0.5,0.5,0.5099019513592785,0.7071067811865476,0.7071067811865475,0.5099019513592785,0.5,0.5,0.5099019513592784,0.7071067811865475,0.5099019513592784,0.5099019513592785,0.5,0.5,0.5,0.5,0.5099019513592785,0.5099019513592785,0.7071067811865476,0.5099019513592785,0.5,0.5,0.5099019513592785,0.7071067811865476,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706631,0.714142842854285,0.5196152422706631,0.5196152422706632,0.5099019513592784,0.5099019513592785,0.5099019513592784,0.5099019513592785,0.5196152422706631,0.5196152422706632,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706632,0.714142842854285,0.8660254037844386,0.714142842854285,0.7071067811865476,0.7071067811865476,0.714142842854285,0.8660254037844386,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706632,0.714142842854285,0.7071067811865475,0.5099019513592785,0.5,0.5,0.5099019513592785,0.7071067811865476,0.7071067811865475,0.5099019513592785,0.5,0.5,0.5099019513592785,0.7071067811865476,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706632,0.714142842854285,0.8660254037844386,0.714142842854285,0.7071067811865476,0.7071067811865476,0.714142842854285,0.8660254037844386]),1e-12))
        self.assertTrue(b.isEqual(DataArrayInt([0,0,3,7,9,9,0,0,3,7,9,9,12,12,14,16,17,17,26,26,28,30,31,31,33,33,36,40,42,42,33,33,36,40,42,42,0,0,3,7,11,9,0,9,12,17,26,31,33,42,33,33,36,40,42,42,45,45,47,49,51,51,45,50,52,53,56,57,58,63,58,58,60,62,63,63,85,85,87,89,91,91,85,90,92,93,96,97,98,103,98,98,100,102,103,103,105,105,108,112,116,114,105,114,117,122,131,136,138,147,138,138,141,145,147,147,105,105,108,112,114,114,105,105,108,112,114,114,117,117,119,121,122,122,131,131,133,135,136,136,138,138,141,145,147,147,138,138,141,145,147,147])))
        pass

    def testSwig2NonRegressionBugDistance1(self):
        pt=DataArrayDouble([(8.8452994616207476,3.1547005383792515,3.1547005383792515)])
        coo=DataArrayDouble([(8,0,0),(8,0,8),(8,8,8),(8,8,0),(16,0,0),(16,0,8),(16,8,8),(16,8,0),(8,0,4),(8,4,8),(8,8,4),(8,4,0),(16,0,4),(16,4,8),(16,8,4),(16,4,0),(12,0,0),(12,0,8),(12,8,8),(12,8,0),(8,4,4),(16,4,4),(12,0,4),(12,4,8),(12,8,4),(12,4,0)])
        conn=DataArrayInt([4,15,21,12,4,16,25,15,12,22,16,4,0,8,20,11,16,0,11,25,22,8,0,16,15,7,14,21,15,25,19,7,7,19,24,14,11,20,10,3,25,11,3,19,19,3,10,24,12,21,13,5,13,23,17,5,5,17,22,12,8,1,9,20,23,9,1,17,17,1,8,22,21,14,6,13,14,24,18,6 ,6,18,23,13,20,9,2,10,24,10,2,18,18,2,9,23])
        m=MEDCouplingUMesh("mesh",2)
        m.setCoords(coo)
        m.allocateCells()
        for i in range(24):
            m.insertNextCell(NORM_QUAD4,conn[4*i:4*i+4])
            pass
        m.checkConsistency()
        m0=m[3] ; m0.zipCoords()
        expectedDist=0.8452994616207476
        a,b=m0.distanceToPoint(pt)
        self.assertAlmostEqual(expectedDist,a,14)
        self.assertEqual(0,b)
        #
        a,b=m.distanceToPoint(pt)
        self.assertAlmostEqual(expectedDist,a,14)
        self.assertEqual(3,b)
        #
        fd=MEDCouplingFieldDiscretization.New(ON_CELLS)
        self.assertEqual(24,fd.getNumberOfTuples(m))
        fd=MEDCouplingFieldDiscretization.New(ON_NODES)
        self.assertEqual(26,fd.getNumberOfTuples(m))
        pass

    def testSwig2AreaBarySeg3Quad8Tri6QPolyg(self):
        #QUAD8 representing a circle of center zeBary and radius zeRadius
        zeBary=[5,6]
        zeRadius=3
        d=DataArrayDouble(8,2)
        d[:,0]=zeRadius
        d[:,1]=[87,-100,-170,110,5,-130,175,95] # angle in degree
        d[:,1]*=pi/180. # angle in radian
        d=d.fromPolarToCart()
        d+=zeBary
        m = MEDCouplingUMesh("quad8", 2) ; m.allocateCells() ; m.insertNextCell(NORM_QUAD8, list(range(8))) ; m.setCoords(d)
        self.assertTrue(m.computeCellCenterOfMass().isEqual(DataArrayDouble(zeBary,1,2),1e-13))
        self.assertAlmostEqual(float(m.getMeasureField(False).getArray()),pi*zeRadius*zeRadius,12)
        tri32D=m.buildDescendingConnectivity()[0][0] ; tri32D.zipCoords()
        # spaceDim=3 QUAD8 becomes QUAD4 ... for the moment
        m.setCoords(m.getCoords().changeNbOfComponents(3,0.))
        m2=m.deepCopy()
        m2.convertQuadraticCellsToLinear()
        self.assertAlmostEqual(float(m.getMeasureField(False).getArray()),float(m2.getMeasureField(False).getArray()),12)
        self.assertTrue(m.computeCellCenterOfMass().isEqual(m2.computeCellCenterOfMass(),1e-13))
        #TRI6 representing a circle of center zeBary and radius zeRadius
        zeBary=[5,6]
        zeRadius=3
        d=DataArrayDouble(6,2)
        d[:,0]=zeRadius
        d[:,1]=[87,-100,110,5,175,95] # angle in degree
        d[:,1]*=pi/180. # angle in radian
        d=d.fromPolarToCart()
        d+=zeBary
        m = MEDCouplingUMesh("tri6", 2) ; m.allocateCells() ; m.insertNextCell(NORM_TRI6, list(range(6))) ; m.setCoords(d)
        self.assertTrue(m.computeCellCenterOfMass().isEqual(DataArrayDouble(zeBary,1,2),1e-13))
        self.assertAlmostEqual(float(m.getMeasureField(False).getArray()),pi*zeRadius*zeRadius,12)
        # spaceDim=3 TRI6 becomes TRI3 ... for the moment
        m.setCoords(m.getCoords().changeNbOfComponents(3,0.))
        m2=m.deepCopy()
        m2.convertQuadraticCellsToLinear()
        self.assertAlmostEqual(float(m.getMeasureField(False).getArray()),float(m2.getMeasureField(False).getArray()),12)
        self.assertTrue(m.computeCellCenterOfMass().isEqual(m2.computeCellCenterOfMass(),1e-13))
        # QPOLYG representing a circle of center zeBary and radius zeRadius
        zeBary=[5,6]
        zeRadius=3
        d=DataArrayDouble(10,2)
        d[:,0]=zeRadius
        d[:,1]=[87,-80,-100,-170,110,5,-90,-130,175,95] # angle in degree
        d[:,1]*=pi/180. # angle in radian
        d=d.fromPolarToCart()
        d+=zeBary
        m = MEDCouplingUMesh("qpolyg", 2) ; m.allocateCells() ; m.insertNextCell(NORM_QPOLYG, list(range(10))) ; m.setCoords(d)
        self.assertTrue(m.computeCellCenterOfMass().isEqual(DataArrayDouble(zeBary,1,2),1e-13))
        self.assertAlmostEqual(float(m.getMeasureField(False).getArray()),pi*zeRadius*zeRadius,12)
        # spaceDim=3 QPOLYG becomes POLYG ... for the moment
        m.setCoords(m.getCoords().changeNbOfComponents(3,0.))
        m2=m.deepCopy()
        m2.convertQuadraticCellsToLinear() ; m2.checkConsistency()
        self.assertTrue(m2.getAllGeoTypes()==[NORM_POLYGON] and m2.getNodalConnectivity().getValues()==[5,0,1,2,3,4])
        self.assertAlmostEqual(float(m.getMeasureField(False).getArray()),float(m2.getMeasureField(False).getArray()),12)
        self.assertTrue(m.computeCellCenterOfMass().isEqual(m2.computeCellCenterOfMass(),1e-13))
        # TRI3
        self.assertAlmostEqual(float(tri32D.getMeasureField(False).getArray()),(87+100)*pi/180*zeRadius,13)
        exp=DataArrayDouble(1,2) ; exp[:,0]=3 ; exp[:,1]=(87-100)/2. ; exp[:,1]*=pi/180. ;  exp=exp.fromPolarToCart() ; exp+=DataArrayDouble([5,6],1,2)
        self.assertTrue(tri32D.computeCellCenterOfMass().isEqual(exp,1e-12))
        # spaceDim=3 TRI3 becomes TRI2 ... for the moment
        tri32D.changeSpaceDimension(3)
        tri2=tri32D.deepCopy() ; tri2.convertQuadraticCellsToLinear()
        self.assertAlmostEqual(float(tri32D.getMeasureField(False).getArray()),float(tri2.getMeasureField(False).getArray()),13)
        self.assertTrue(tri32D.computeCellCenterOfMass().isEqual(tri2.computeCellCenterOfMass(),1e-12))
        tri32D.changeSpaceDimension(1)
        self.assertAlmostEqual(float(tri32D.getMeasureField(False).getArray()),-0.67795240172962323,12)
        pass

    # this bug 5/6/2013 is swig specific
    def testSwigNonRegressionBugRotate3D1(self):
        m=MEDCouplingUMesh.New()
        dataArray=DataArrayDouble.New(100,3)
        dataArray[:]=0.
        dataArray[0]=[0.,1,3]
        m.setCoords(dataArray[0])
        m1=m.deepCopy()
        m.rotate([0.,0.,3.],[1.,0.,0.],0.5*pi)
        self.assertTrue(m.getCoords().isEqual(DataArrayDouble([0.,0.,4.],1,3),1e-15))
        #
        d1=DataArrayDouble([0.,0.,3.],1,3) ; d2=DataArrayDouble([1.,0.,0.],1,3)
        pts=[[0.,0.,3.],[(0.,0.,3.)],DataArrayDouble([0.,0.,3.],1,3),list(d1)[0]]
        vec=[[1.,0.,0.],[(1.,0.,0.)],DataArrayDouble([1.,0.,0.],1,3),list(d2)[0]]
        for p in pts:
            for v in vec:
                m2=m1.deepCopy()
                m2.rotate(p,v,0.5*pi)
                self.assertTrue(m2.getCoords().isEqual(DataArrayDouble([0.,0.,4.],1,3),1e-15))
                pass
        pass

    def testSwig2DataArrayCount1(self):
        d=DataArrayInt([])
        self.assertEqual(0,d.getNumberOfTuples())
        self.assertEqual(1,d.getNumberOfComponents())
        self.assertEqual(0,d.count(0))
        self.assertEqual(0,d.count(1))
        self.assertEqual(0,d.count(-1))
        d=DataArrayInt([2,1,-2,-3,2,0,0,7,2,-2,3,0])
        self.assertEqual(12,d.getNumberOfTuples())
        self.assertEqual(1,d.getNumberOfComponents())
        self.assertEqual(3,d.count(0))
        self.assertEqual(1,d.count(1))
        self.assertEqual(0,d.count(-1))
        self.assertEqual(2,d.count(-2))
        self.assertEqual(3,d.count(2))
        e=d.getDifferentValues()
        f=DataArrayInt()
        for it in e:
            f.pushBackSilent(d.count(int(it)))
            pass
        self.assertEqual(12,f.accumulate()[0])
        #
        eps=1e-12
        d=DataArrayDouble([])
        self.assertEqual(0,d.getNumberOfTuples())
        self.assertEqual(1,d.getNumberOfComponents())
        self.assertEqual(0,d.count(0,eps))
        self.assertEqual(0,d.count(1,eps))
        self.assertEqual(0,d.count(-1,eps))
        d=DataArrayDouble([2,1,-2,-3,2,0,eps/10,7,2+eps/10,-2,3,0])
        self.assertEqual(12,d.getNumberOfTuples())
        self.assertEqual(1,d.getNumberOfComponents())
        self.assertEqual(3,d.count(0,eps))
        self.assertEqual(1,d.count(1,eps))
        self.assertEqual(0,d.count(-1,eps))
        self.assertEqual(2,d.count(-2,eps))
        self.assertEqual(3,d.count(2,eps))
        self.assertEqual(3,d.count(2,eps))
        self.assertEqual(2,d.count(2,eps/100))
        e=d.getDifferentValues(eps)
        f=DataArrayInt()
        for it in e:
            f.pushBackSilent(d.count(float(it),eps))
            pass
        self.assertEqual(12,f.accumulate()[0])
        pass

    def testSwig2DataArrayGetSlice1(self):
        s=slice(2,18,1)
        self.assertEqual(DataArray.GetNumberOfItemGivenBESRelative(s),16)
        self.assertEqual(DataArray.GetNumberOfItemGivenBES(s),16)
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(2,6,1))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(6,10,1))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(10,14,1))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(14,18,1))
        #
        s=slice(2,18,2)
        self.assertEqual(DataArray.GetNumberOfItemGivenBESRelative(s),8)
        self.assertEqual(DataArray.GetNumberOfItemGivenBES(s),8)
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(2,6,2))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(6,10,2))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(10,14,2))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(14,18,2))
        #
        s=slice(1,18,1)
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(1,5,1))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(5,9,1))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(9,13,1))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(13,18,1))# 18 not 17
        #
        s=slice(1,18,2)
        self.assertEqual(DataArray.GetNumberOfItemGivenBESRelative(s),9)
        self.assertEqual(DataArray.GetNumberOfItemGivenBES(s),9)
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(1,5,2))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(5,9,2))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(9,13,2))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(13,18,2))# 18 not 17
        #
        s=slice(18,2,-1)
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(18,14,-1))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(14,10,-1))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(10,6,-1))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(6,2,-1))
        #
        s=slice(18,2,-2)
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(18,14,-2))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(14,10,-2))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(10,6,-2))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(6,2,-2))
        #
        s=slice(18,1,-1)
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(18,14,-1))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(14,10,-1))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(10,6,-1))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(6,1,-1))# 1 not 2
        #
        s=slice(18,1,-2)
        self.assertEqual(DataArray.GetNumberOfItemGivenBESRelative(s),9)
        self.assertRaises(InterpKernelException,DataArray.GetNumberOfItemGivenBES,s)
        self.assertEqual(sum([DataArray.GetNumberOfItemGivenBESRelative(DataArray.GetSlice(s, i, 4)) for i in range(4)]), DataArray.GetNumberOfItemGivenBESRelative(s))
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(18,14,-2))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(14,10,-2))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(10,6,-2))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(6,1,-2))# 1 not 2
        self.assertRaises(InterpKernelException,DataArray.GetSlice,slice(0,None,2),0,4)
        #
        d=DataArrayInt.Range(0,18,1)
        s=slice(2,None,1)
        self.assertEqual(d.getNumberOfItemGivenBES(s),16)
        self.assertEqual(d.getNumberOfItemGivenBESRelative(s),16)
        self.assertEqual(d.getSlice(s,0,4),slice(2,6,1))
        self.assertEqual(d.getSlice(s,1,4),slice(6,10,1))
        self.assertEqual(d.getSlice(s,2,4),slice(10,14,1))
        self.assertEqual(d.getSlice(s,3,4),slice(14,18,1))
        #
        d=DataArrayInt.Range(0,18,1)
        s=slice(2,-2,1)
        self.assertEqual(d.getSlice(s,0,4),slice(2,5,1))
        self.assertEqual(d.getSlice(s,1,4),slice(5,8,1))
        self.assertEqual(d.getSlice(s,2,4),slice(8,11,1))
        self.assertEqual(d.getSlice(s,3,4),slice(11,16,1))
        #
        d=DataArrayInt.Range(0,18,1)
        s=slice(None,None,1)
        self.assertEqual(d.getSlice(s,0,4),slice(0,4,1))
        self.assertEqual(d.getSlice(s,1,4),slice(4,8,1))
        self.assertEqual(d.getSlice(s,2,4),slice(8,12,1))
        self.assertEqual(d.getSlice(s,3,4),slice(12,18,1))
        #
        d=DataArrayInt.Range(0,18,1)
        s=slice(None,2,-2)
        self.assertRaises(InterpKernelException,d.getNumberOfItemGivenBES,s)
        self.assertEqual(d.getNumberOfItemGivenBESRelative(s),8)
        self.assertEqual(d.getSlice(s,0,4),slice(17,13,-2))
        self.assertEqual(d.getSlice(s,1,4),slice(13,9,-2))
        self.assertEqual(d.getSlice(s,2,4),slice(9,5,-2))
        self.assertEqual(d.getSlice(s,3,4),slice(5,2,-2))
        pass

    def testSwig2AccumulatePerChunk1(self):
        arr=DataArrayDouble(11) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr)
        m=m.buildUnstructured()
        m0=m[::2] ; ids0=m0.simplexize(0) ; m1=m[1::2]
        m=MEDCouplingUMesh.MergeUMeshesOnSameCoords(m0,m1) ; m.setName("mesh")
        m.checkConsecutiveCellTypesForMEDFileFrmt()
        #
        formula="7-sqrt((x-5.)*(x-5.)+(y-5.)*(y-5.))"
        f=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f.setMesh(m)
        f.fillFromAnalytic(1,formula)
        f.setName("Field1") ; f.setTime(1.1,1,-1)
        f.checkConsistencyLight()
        #
        arr=f.getArray()
        arr2=DataArrayDouble(len(arr),2) ; arr2[:,0]=arr
        arr2=DataArrayDouble(len(arr),2) ; arr2[:,0]=arr ; arr2[:,1]=2*arr
        f.setArray(arr2)
        f.checkConsistencyLight()
        # here the compact code to obviously put field on cell to nodes
        rn,rni=f.getMesh().getReverseNodalConnectivity()
        arr2=f.getArray()[rn]
        arr4=arr2.accumulatePerChunck(rni)
        nbOfCellsSharingNodes=rni.deltaShiftIndex()
        arr4/=nbOfCellsSharingNodes.convertToDblArr()
        #
        maxNbCSN=nbOfCellsSharingNodes.getMaxValue()[0]
        arr3=DataArrayDouble(f.getMesh().getNumberOfNodes(),f.getArray().getNumberOfComponents()) ; arr3[:]=0.
        for i in range(1, maxNbCSN + 1):
            ids=nbOfCellsSharingNodes.findIdsEqual(i)
            if len(ids)==0:
                continue
            for j in range(i):
                rni2=rni[ids] ; rni2+=j
                arr3[ids]+=arr2[rni2]
                pass
            arr3[ids]/=i
            pass
        fNode=MEDCouplingFieldDouble(ON_NODES,ONE_TIME) ; fNode.setMesh(m)
        fNode.setName("Field1Node") ; fNode.setTime(1.1,1,-1)
        fNode.setArray(arr3) ; fNode.checkConsistencyLight()
        self.assertTrue(arr3.isEqual(arr4,1e-12))
        #
        d=DataArrayInt.Range(0,20,1)
        self.assertTrue(d.accumulatePerChunck([2,4,12]).isEqual(DataArrayInt([5,60])))
        #
        a=DataArrayDouble(12) ; a.iota() ; a.rearrange(3)
        b=DataArrayDouble(12) ; b.iota(20) ; b.rearrange(3)
        ids=DataArrayInt([])
        self.assertEqual(len(a[ids]),0)
        self.assertEqual(len(b[ids]),0)
        a2=a.deepCopy() ;  a2[ids]+=b[ids] ; self.assertTrue(a2.isEqual(a,1e-15))
        a2=a.deepCopy() ;  a2[ids]*=b[ids] ; self.assertTrue(a2.isEqual(a,1e-15))
        a2=a.deepCopy() ;  a2[ids]/=b[ids] ; self.assertTrue(a2.isEqual(a,1e-15))
        a2=a.deepCopy() ;  a2[ids]-=b[ids] ; self.assertTrue(a2.isEqual(a,1e-15))
        pass

    def testSwig2CheckAndPreparePermutation1(self):
        a=DataArrayInt([10003,9999999,5,67])
        self.assertTrue(a.checkAndPreparePermutation().isEqual(DataArrayInt([2,3,0,1])))
        a=DataArrayInt([10003,-9999999,5,67])
        self.assertTrue(a.checkAndPreparePermutation().isEqual(DataArrayInt([3,0,1,2])))
        a=DataArrayInt([])
        self.assertTrue(a.checkAndPreparePermutation().isEqual(DataArrayInt([])))
        a=DataArrayInt([])
        a.iota();
        self.assertTrue(a.isEqual(DataArrayInt([])))
        pass

    def testSwig21SGTUMesh1(self):
        m=MEDCoupling1GTUMesh.New("m",NORM_PENTA6)
        m.__repr__() ; m.__str__()
        self.assertTrue(isinstance(m,MEDCoupling1SGTUMesh))
        m.setCoords(DataArrayDouble(20,3))
        m.allocateCells()
        m.__repr__() ; m.__str__()
        m.insertNextCell([0,1,2,5,7,2])
        self.assertEqual(1,m.getNumberOfCells())
        self.assertTrue(DataArrayInt([6]).isEqual(m.computeNbOfNodesPerCell()))
        self.assertTrue(DataArrayInt([5]).isEqual(m.computeNbOfFacesPerCell()))
        m.__repr__() ; m.__str__()
        m.checkConsistencyLight()
        m.checkConsistency()
        #
        cm=MEDCouplingCMesh() ; cm.setName("m")
        arr0=DataArrayDouble(6) ; arr0.iota()
        arr1=DataArrayDouble([0,1])
        cm.setCoords(arr0,arr1,arr1) ; um=cm.buildUnstructured()
        #
        m=MEDCoupling1SGTUMesh("m",NORM_QUAD4)
        mem_m=m.getHeapMemorySize()
        m.allocateCells(5)
        self.assertIn(m.getHeapMemorySize() - mem_m, list(range(5 * 4 * 4, 5 * 4 * 4 + 32)))
        self.assertEqual(m.getNodalConnectivity().getNbOfElemAllocated(),20)
        m.setCoords(um.getCoords())
        m.insertNextCell([1,0,6,7])
        self.assertEqual(1,m.getNumberOfCells())
        m.insertNextCell([2,1,7,8])
        m.insertNextCell([3,2,8,9])
        m.insertNextCell([4,3,9,10])
        m.insertNextCell([5,4,10,11])
        self.assertEqual(5,m.getNumberOfCells())
        self.assertRaises(InterpKernelException,m.insertNextCell,[0,6,7])
        self.assertRaises(InterpKernelException,m.insertNextCell,[0,6,7,1,2])
        self.assertEqual(m.getNodalConnectivity().getNbOfElemAllocated(),20)
        f=m.getMeasureField(False)
        self.assertEqual(f.getMesh().getHiddenCppPointer(),m.getHiddenCppPointer())
        self.assertTrue(f.getArray().isUniform(1,1e-14))
        self.assertEqual(m.getType(),10)
        self.assertEqual(m.getCellModelEnum(),NORM_QUAD4)
        mo=MEDCoupling1SGTUMesh("m",NORM_QUAD4) ; mo.setCoords(m.getCoords())
        mo.setNodalConnectivity(DataArrayInt([1,0,6,7,2,1,7,8,3,2,8,9,4,3,9,10,5,4,10,11]))
        self.assertTrue(m.isEqual(mo,1e-12))
        #
        mo2=MEDCoupling1SGTUMesh.Merge1SGTUMeshesOnSameCoords([m[[0,1]],m[[2]],m[[3,4]]])
        mo2.setName(m.getName())
        self.assertTrue(m.isEqual(mo2,1e-12))
        #
        mp0=m[[0]] ; mp0.zipCoords() ; mp1=m[2] ; mp1.zipCoords() ; mp2=m[4] ; mp2.zipCoords()
        mo3=MEDCoupling1SGTUMesh.Merge1SGTUMeshes([mp0,mp1,mp2])
        self.assertTrue(isinstance(mo3,MEDCoupling1SGTUMesh))
        mo3.setName(m.getName())
        m_ref=m[(0,2,4)] ; m_ref.zipCoords()
        m_ref.tryToShareSameCoordsPermute(mo3,1e-12)
        self.assertTrue(m_ref.isEqual(mo3,1e-12))
        #
        m1=um.buildDescendingConnectivity()[0]
        ids=m1.getCellIdsFullyIncludedInNodeIds(DataArrayInt.Range(0,12,1))
        m1=m1[ids]
        m1c=m1.convertIntoSingleGeoTypeMesh()
        self.assertTrue(isinstance(m1c,MEDCoupling1SGTUMesh))
        self.assertEqual(m1c.getCoords().getHiddenCppPointer(),m.getCoords().getHiddenCppPointer())
        m1c.checkConsistency()
        self.assertTrue(m1c.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,2,1,7,8,3,2,8,9,4,3,9,10,5,4,10,11])))
        self.assertEqual(20,m1c.getNodalConnectivityLength())
        self.assertTrue(m.isEqual(m1c,1e-12))
        m.getNodalConnectivity().setIJ(1,0,1)
        self.assertTrue(not m.isEqual(m1c,1e-12))
        m.getNodalConnectivity().setIJ(1,0,0)
        self.assertTrue(m.isEqual(m1c,1e-12))
        m1c.setCoords(m.getCoords().deepCopy())
        self.assertTrue(m.isEqual(m1c,1e-12))
        m1c.getCoords().setIJ(0,1,0.1)
        self.assertTrue(not m.isEqual(m1c,1e-12))
        m1c.getCoords().setIJ(0,1,0)
        self.assertTrue(m.isEqual(m1c,1e-12))
        m1c.getCoords().setInfoOnComponent(1,"X")
        self.assertTrue(not m.isEqual(m1c,1e-12) and m.isEqualWithoutConsideringStr(m1c,1e-12))
        m.getCoords().setInfoOnComponent(1,"X")
        self.assertTrue(m.isEqual(m1c,1e-12) and m.isEqualWithoutConsideringStr(m1c,1e-12))
        m.setName("m2")
        self.assertTrue(not m.isEqual(m1c,1e-12) and m.isEqualWithoutConsideringStr(m1c,1e-12))
        #
        m.checkConsistencyLight() ; m.checkConsistency() ; m.checkConsistency()
        self.assertEqual(m.getMeshDimension(),2)
        self.assertTrue(m.giveCellsWithType(NORM_QUAD4).isEqual(DataArrayInt([0,1,2,3,4])))
        self.assertTrue(m.giveCellsWithType(NORM_TRI3).isEqual(DataArrayInt([])))
        self.assertEqual(m.getNumberOfCellsWithType(NORM_QUAD4),5)
        self.assertEqual(m.getNumberOfCellsWithType(NORM_TRI3),0)
        self.assertEqual(m.getTypeOfCell(3),NORM_QUAD4)
        self.assertRaises(InterpKernelException,m.getTypeOfCell,5)
        self.assertEqual(m.getAllGeoTypes(),[NORM_QUAD4])
        self.assertEqual(m.getDistributionOfTypes(),[[NORM_QUAD4,5,-1]])
        ##
        pfl1=DataArrayInt([1,3,4])
        a,b,c=m.splitProfilePerType(pfl1)
        d,e,f=m.buildUnstructured().splitProfilePerType(pfl1)
        self.assertTrue(a==[[4,3,0]] and len(b)==1 and b[0].isEqual(DataArrayInt([0,1,2])) and len(c)==1 and c[0].getHiddenCppPointer()==pfl1.getHiddenCppPointer())
        self.assertTrue(a==d and len(b)==1 and b[0].isEqual(e[0]) and len(c)==1 and c[0].isEqual(f[0]))
        #
        pfl2=DataArrayInt([0,1,2,3])
        a,b,c=m.splitProfilePerType(pfl2)
        d,e,f=m.buildUnstructured().splitProfilePerType(pfl2)
        self.assertTrue(a==[[4,4,0]] and len(b)==1 and b[0].isEqual(DataArrayInt([0,1,2,3])) and len(c)==1 and c[0].getHiddenCppPointer()==pfl2.getHiddenCppPointer())
        self.assertTrue(a==d and len(b)==1 and b[0].isEqual(e[0]) and len(c)==1 and c[0].isEqual(f[0]))
        #
        pfl3=DataArrayInt([0,1,2,3,4])
        a,b,c=m.splitProfilePerType(pfl3)
        d,e,f=m.buildUnstructured().splitProfilePerType(pfl3)
        self.assertTrue(a==[[4,5,-1]] and len(b)==1 and b[0].isEqual(DataArrayInt([0,1,2,3,4])) and c==[])
        self.assertTrue(a==d and len(b)==1 and b[0].isEqual(e[0]) and c==[])
        #
        invalidPfl=DataArrayInt([1,2,3,4,5])
        self.assertRaises(InterpKernelException,m.splitProfilePerType,invalidPfl)
        self.assertRaises(InterpKernelException,m.buildUnstructured().splitProfilePerType,invalidPfl)
        ##
        pfl1=DataArrayInt([1,2,3])
        a=m.checkTypeConsistencyAndContig([NORM_QUAD4,3,0],[pfl1])
        b=m.buildUnstructured().checkTypeConsistencyAndContig([NORM_QUAD4,3,0],[pfl1])
        self.assertTrue(a.isEqual(b) and pfl1.getHiddenCppPointer(),a.getHiddenCppPointer())
        #
        pfl2=DataArrayInt([0,1,2,3])
        a=m.checkTypeConsistencyAndContig([NORM_QUAD4,4,0],[pfl2])
        b=m.buildUnstructured().checkTypeConsistencyAndContig([NORM_QUAD4,4,0],[pfl2])
        self.assertTrue(a.isEqual(b) and pfl2.getHiddenCppPointer()==a.getHiddenCppPointer())
        #
        pfl3=DataArrayInt([0,1,2,3,4])
        a=m.checkTypeConsistencyAndContig([NORM_QUAD4,4,0],[pfl3])
        b=m.buildUnstructured().checkTypeConsistencyAndContig([NORM_QUAD4,5,0],[pfl3])
        self.assertTrue(a.isEqual(b) and pfl3.getHiddenCppPointer()==a.getHiddenCppPointer())
        #
        invalidPfl=DataArrayInt([1,2,3,4,5])
        self.assertRaises(InterpKernelException,m.checkTypeConsistencyAndContig,[NORM_QUAD4,5,0],[invalidPfl])
        self.assertRaises(InterpKernelException,m.buildUnstructured().checkTypeConsistencyAndContig,[NORM_QUAD4,5,0],[invalidPfl])
        ##
        self.assertTrue(DataArrayInt([4,4,4,4,4]).isEqual(m.computeNbOfNodesPerCell()))
        ##
        self.assertEqual(m.getNodeIdsOfCell(1),[2,1,7,8])
        ##
        self.assertTrue(m.computeIsoBarycenterOfNodesPerCell().isEqual(DataArrayDouble([(0.5,0.5,0),(1.5,0.5,0),(2.5,0.5,0),(3.5,0.5,0),(4.5,0.5,0)]),1e-13))
        ##
        ref=m.getCoords().getHiddenCppPointer()
        mcpy=m.deepCopy() ; mcpy.insertNextCell([1,0,6,7])
        c=m.getNodalConnectivity().deepCopy()
        o2n=DataArrayInt([2,0,1,4,3])
        m.renumberCells(o2n,False)
        c.rearrange(4) ; c.renumberInPlace(o2n) ; c.rearrange(1)
        self.assertTrue(c.isEqual(m.getNodalConnectivity()))
        self.assertEqual(ref,m.getCoords().getHiddenCppPointer())
        m2=mcpy.mergeMyselfWith(m)
        self.assertTrue(isinstance(m2,MEDCoupling1SGTUMesh))
        self.assertEqual(11,m2.getNumberOfCells())
        self.assertEqual(48,m2.getNumberOfNodes())
        self.assertTrue(m2.getCoords().isEqual(DataArrayDouble.Aggregate([m.getCoords(),m.getCoords()]),1e-12))
        self.assertTrue(m2.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,2,1,7,8,3,2,8,9,4,3,9,10,5,4,10,11,1,0,6,7,26,25,31,32,27,26,32,33,25,24,30,31,29,28,34,35,28,27,33,34])))
        ##
        mu=m.buildUnstructured()
        mu.checkConsistency()
        self.assertEqual(mu.getCoords().getHiddenCppPointer(),m.getCoords().getHiddenCppPointer())
        self.assertEqual(2,mu.getMeshDimension())
        self.assertEqual([NORM_QUAD4],mu.getAllGeoTypes())
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([4,2,1,7,8,4,3,2,8,9,4,1,0,6,7,4,5,4,10,11,4,4,3,9,10])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,10,15,20,25])))
        ##
        for typ in [0,1]:
            mcpy2=m.deepCopy() ; umcpy2=mcpy2.buildUnstructured()
            ids=mcpy2.simplexize(typ) ; ids2=umcpy2.simplexize(typ)
            self.assertTrue(ids.isEqual(ids2))
            mcpy3=umcpy2.convertIntoSingleGeoTypeMesh()
            self.assertTrue(mcpy2.isEqual(mcpy3,1e-14))
            pass
        um1=um.convertIntoSingleGeoTypeMesh()
        self.assertEqual(8,um1.getNumberOfNodesPerCell())
        for typ in [PLANAR_FACE_5,PLANAR_FACE_6]:
            mcpy2=um1.deepCopy() ; umcpy2=mcpy2.buildUnstructured()
            ids=mcpy2.simplexize(typ) ; ids2=umcpy2.simplexize(typ)
            self.assertTrue(ids.isEqual(ids2))
            mcpy3=umcpy2.convertIntoSingleGeoTypeMesh()
            self.assertTrue(mcpy2.isEqual(mcpy3,1e-14))
            pass
        ##
        self.assertRaises(InterpKernelException,mcpy.mergeMyselfWithOnSameCoords,m)
        mcpy.tryToShareSameCoords(m,1e-14)
        m3=mcpy.mergeMyselfWithOnSameCoords(m)
        self.assertTrue(isinstance(m3,MEDCoupling1SGTUMesh))
        self.assertEqual(11,m3.getNumberOfCells())
        self.assertEqual(24,m3.getNumberOfNodes())
        self.assertEqual(m3.getCoords().getHiddenCppPointer(),mcpy.getCoords().getHiddenCppPointer())
        self.assertTrue(m3.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,2,1,7,8,3,2,8,9,4,3,9,10,5,4,10,11,1,0,6,7,2,1,7,8,3,2,8,9,1,0,6,7,5,4,10,11,4,3,9,10])))
        ##
        ref=mcpy.getCoords().deepCopy()
        c3=mcpy.getNodalConnectivity()[:]
        mcpy.getNodalConnectivity().setIJ(int(c3.findIdsEqual(11)),0,24)
        c2=DataArrayDouble.Aggregate([mcpy.getCoords(),mcpy.getCoords()[11:]])
        mcpy.setCoords(c2)
        mcpy.checkConsistency()
        a,b=mcpy.getNodeIdsInUse()
        self.assertEqual(12,b)
        self.assertTrue(a.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])))
        ids=mcpy.zipCoordsTraducer()
        self.assertTrue(ids.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])))
        self.assertTrue(mcpy.getCoords().isEqual(ref[:12],1e-12))
        self.assertTrue(mcpy.getNodalConnectivity().isEqual(c3))
        mcpy.checkConsistency()
        ##
        m4=mcpy[DataArrayInt([0,3,4])]
        m5=mcpy.buildPartOfMySelfKeepCoords(DataArrayInt([0,3,4]))
        self.assertTrue(isinstance(m4,MEDCoupling1SGTUMesh))
        self.assertTrue(m4.isEqual(m5,-1e-14))# < 0 not a bug it proves that coordinates pointer are equal
        self.assertTrue(m4.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,4,3,9,10,5,4,10,11])))
        m6=mcpy[::2]
        self.assertTrue(isinstance(m6,MEDCoupling1SGTUMesh))
        self.assertTrue(m6.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,3,2,8,9,5,4,10,11])))
        ##
        mcpy.setCoords(DataArrayDouble.Aggregate([mcpy.getCoords(),mcpy.getCoords()]))
        mcpy.checkConsistency()
        ##
        mcppy=mcpy.deepCopyConnectivityOnly()
        self.assertTrue(mcppy.isEqual(mcpy,1e-12))
        self.assertTrue(mcppy.getCoords().getHiddenCppPointer()==mcpy.getCoords().getHiddenCppPointer())
        self.assertTrue(mcppy.getNodalConnectivity().isEqual(mcpy.getNodalConnectivity()))
        self.assertTrue(mcppy.getNodalConnectivity().getHiddenCppPointer()!=mcpy.getNodalConnectivity().getHiddenCppPointer())
        ##
        a,b=mcpy.getReverseNodalConnectivity()
        self.assertTrue(a.isEqual(DataArrayInt([0,5,0,1,5,1,2,2,3,3,4,4,0,5,0,1,5,1,2,2,3,3,4,4])))
        self.assertTrue(b.isEqual(DataArrayInt([0,2,5,7,9,11,12,14,17,19,21,23,24,24,24,24,24,24,24,24,24,24,24,24,24])))
        self.assertTrue(mcpy.fillCellIdsToKeepFromNodeIds([0,1,6,7],False).isEqual(DataArrayInt([0,1,5])))
        self.assertTrue(mcpy.fillCellIdsToKeepFromNodeIds([0,1,6,7],True).isEqual(DataArrayInt([0,5])))
        self.assertTrue(mcpy.getCellsInBoundingBox([(0,1),(0,1),(0,1)],1e-12).isEqual(DataArrayInt([0,1,5])))
        f=mcpy.buildOrthogonalField()
        self.assertEqual(f.getMesh().getHiddenCppPointer(),mcpy.getHiddenCppPointer())
        self.assertTrue(f.getArray().isEqual(DataArrayDouble(6*[(0,0,-1)]),1e-12))
        mcpy.changeSpaceDimension(2)
        self.assertEqual(1,mcpy.getCellContainingPoint([1.5,0.5],1e-12))
        ##
        self.assertTrue(mcpy.fillCellIdsToKeepFromNodeIds(DataArrayInt([6,7]),False).isEqual(DataArrayInt([0,1,5])))
        ##
        mcpy2=mcpy.deepCopy()
        self.assertEqual([None,None],mcpy.checkGeoEquivalWith(mcpy2,1,1e-12))#fast equal
        mcpy.checkFastEquivalWith(mcpy2,1e-12)
        mcpy2.renumberCells([0,2,4,3,1,5])
        mcpy.checkFastEquivalWith(mcpy2,1e-12)
        self.assertEqual([None,None],mcpy.checkGeoEquivalWith(mcpy2,1,1e-12))#fast equal
        mcpy2.renumberCells([0,2,4,3,1,5])
        mcpy2.renumberCells([1,3,5,0,2,4])
        self.assertRaises(InterpKernelException,mcpy.checkFastEquivalWith,mcpy2,1e-12)
        self.assertRaises(InterpKernelException,mcpy.checkGeoEquivalWith,mcpy2,1,1e-12)#fast equal
        pass

    def testSwig21DGTUMesh1(self):
        a0=DataArrayInt([0,2,3,5,6,8])
        a1=DataArrayInt([0,4,7,11,14,18,21,25])
        a2=DataArrayInt([0,1,4,5])
        self.assertTrue(DataArrayInt.AggregateIndexes([a0,a1,a2]).isEqual(DataArrayInt([0,2,3,5,6,8,12,15,19,22,26,29,33,34,37,38])))
        self.assertEqual(a1[3:].front(),11)
        self.assertEqual(a1[4:].convertToDblArr().front(),14.)
        a1c=DataArrayInt([5,7,1,2, 8,11,0, 5,6,3,12, 1,5,2, 13,12,11,7, 6,1,0, 20,21,19,17])
        d,e=MEDCouplingUMesh.ExtractFromIndexedArraysSlice(1,5,2,a1c,a1)
        self.assertTrue(d.isEqual(DataArrayInt([8,11,0,1,5,2])))
        self.assertTrue(e.isEqual(DataArrayInt([0,3,6])))
        #
        m=MEDCouplingDataForTest.build2DTargetMesh_1()[0,3,4]
        ref=DataArrayInt([0,3,4,1,6,7,4,3,7,8,5,4])
        self.assertTrue(m.convertNodalConnectivityToStaticGeoTypeMesh().isEqual(ref))
        d,e=m.convertNodalConnectivityToDynamicGeoTypeMesh()
        self.assertTrue(d.isEqual(ref))
        self.assertTrue(e.isEqual(DataArrayInt.Range(0,13,4)))
        self.assertTrue(m.fillCellIdsToKeepFromNodeIds(DataArrayInt([6,7]),False).isEqual(DataArrayInt([1,2])))
        #
        m=MEDCoupling1GTUMesh.New("m",NORM_POLYHED)
        self.assertTrue(isinstance(m,MEDCoupling1DGTUMesh))
        m.__repr__() ; m.__str__()
        m.setCoords(DataArrayDouble(20,3))
        m.allocateCells()
        m.__repr__() ; m.__str__()
        m.insertNextCell([0,1,2,5,7,2,-1,1,3])
        self.assertEqual(1,m.getNumberOfCells())
        self.assertTrue(DataArrayInt([8]).isEqual(m.computeNbOfNodesPerCell()))
        self.assertTrue(DataArrayInt([2]).isEqual(m.computeNbOfFacesPerCell()))
        m.__repr__() ; m.__str__()
        m.checkConsistencyLight()
        m.checkConsistency()
        #
        cm=MEDCouplingCMesh() ; cm.setName("m")
        arr0=DataArrayDouble(6) ; arr0.iota()
        arr1=DataArrayDouble([0,1])
        cm.setCoords(arr0,arr1,arr1) ; um=cm.buildUnstructured() ; um.convertAllToPoly()
        um2=um.deepCopyConnectivityOnly()
        self.assertTrue(um2.isEqual(um,1e-12))
        self.assertEqual(um2.getCoords().getHiddenCppPointer(),um.getCoords().getHiddenCppPointer())
        self.assertTrue(um2.getNodalConnectivity().isEqual(um.getNodalConnectivity()))
        self.assertTrue(um2.getNodalConnectivity().getHiddenCppPointer()!=um.getNodalConnectivity().getHiddenCppPointer())
        self.assertTrue(um2.getNodalConnectivityIndex().isEqual(um.getNodalConnectivityIndex()))
        self.assertTrue(um2.getNodalConnectivityIndex().getHiddenCppPointer()!=um.getNodalConnectivityIndex().getHiddenCppPointer())
        #
        self.assertRaises(InterpKernelException,MEDCoupling1SGTUMesh.New,"m",NORM_POLYHED)
        m=MEDCoupling1DGTUMesh("m",NORM_POLYHED)
        m.allocateCells(5)
        self.assertEqual(15,m.getNodalConnectivity().getNbOfElemAllocated())
        self.assertEqual(6,m.getNodalConnectivityIndex().getNbOfElemAllocated())
        m.setCoords(um.getCoords())
        m.insertNextCell([1,0,6,7,-1,7,6,1])
        self.assertEqual(1,m.getNumberOfCells())
        m.insertNextCell([2,1,7,8,-1,2,1,-1,8,-1,7])
        m.insertNextCell([3,2,8,9])
        m.insertNextCell([4,3,9,10,-1,5,3,9])
        m.insertNextCell([5,4,10,11,-1,11,10,-1,5])
        m.checkConsistencyLight()
        m.checkConsistency()
        self.assertEqual(5,m.getNumberOfCells())
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,8,19,23,31,40])))
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,-1,7,6,1,2,1,7,8,-1,2,1,-1,8,-1,7,3,2,8,9,4,3,9,10,-1,5,3,9,5,4,10,11,-1,11,10,-1,5])))
        #
        m4=m.deepCopy()
        self.assertTrue(m.isEqual(m4,1e-12))
        m4.getNodalConnectivity().setIJ(2,0,5)
        self.assertTrue(not m.isEqual(m4,1e-12))
        m4.getNodalConnectivity().setIJ(2,0,6)
        self.assertTrue(m.isEqual(m4,1e-12))
        m4.getNodalConnectivityIndex().setIJ(2,0,21)
        self.assertTrue(not m.isEqual(m4,1e-12))
        m4.getNodalConnectivityIndex().setIJ(2,0,19)
        self.assertTrue(m.isEqual(m4,1e-12))
        m4.getCoords().setIJ(10,1,1.1)
        self.assertTrue(not m.isEqual(m4,1e-12))
        m4.getCoords().setIJ(10,1,1.)
        self.assertTrue(m.isEqual(m4,1e-12))
        m4.getNodalConnectivity().pushBackSilent(7)
        self.assertTrue(not m.isEqual(m4,1e-12))
        self.assertEqual(7,m4.getNodalConnectivity().popBackSilent())
        self.assertTrue(m.isEqual(m4,1e-12))
        m4.setName("m4")
        self.assertTrue(not m.isEqual(m4,1e-12))
        m4.setName("m")
        self.assertTrue(m.isEqual(m4,1e-12))
        #
        self.assertEqual(6,m.getNodalConnectivityIndex().getNbOfElemAllocated())
        self.assertEqual(60,m.getNodalConnectivity().getNbOfElemAllocated())
        self.assertTrue(m.computeNbOfNodesPerCell().isEqual(DataArrayInt([7,8,4,7,7])))
        self.assertTrue(m.computeNbOfFacesPerCell().isEqual(DataArrayInt([2,4,1,2,3])))
        self.assertEqual(m.getNodeIdsOfCell(1),[2,1,7,8,-1,2,1,-1,8,-1,7])
        f=m.computeIsoBarycenterOfNodesPerCell()
        self.assertTrue(DataArrayDouble([(0.5714285714285714,0.5714285714285714,0),(1.5,0.5,0),(2.5,0.5,0),(3.5714285714285712,0.42857142857142855,0),(4.5714285714285712,0.5714285714285714,0)]).isEqual(f,1e-14))
        mu0=m.buildUnstructured()
        o2n=[1,2,0,4,3]
        m2=m.deepCopy()
        m3=m.deepCopyConnectivityOnly()
        self.assertTrue(m3.isEqual(m,1e-12))
        self.assertEqual(m3.getCoords().getHiddenCppPointer(),m.getCoords().getHiddenCppPointer())
        self.assertTrue(m3.getNodalConnectivity().getHiddenCppPointer()!=m.getNodalConnectivity().getHiddenCppPointer())
        self.assertTrue(m3.getNodalConnectivity().isEqual(m.getNodalConnectivity()))
        self.assertTrue(m3.getNodalConnectivityIndex().getHiddenCppPointer()!=m.getNodalConnectivityIndex().getHiddenCppPointer())
        self.assertTrue(m3.getNodalConnectivityIndex().isEqual(m.getNodalConnectivityIndex()))
        m.renumberCells(o2n)
        mu0.renumberCells(o2n)
        self.assertTrue(mu0.isEqual(m.buildUnstructured(),1e-12))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,12,23,32,40])))
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([3,2,8,9,1,0,6,7,-1,7,6,1,2,1,7,8,-1,2,1,-1,8,-1,7,5,4,10,11,-1,11,10,-1,5,4,3,9,10,-1,5,3,9])))
        #
        mcpy0=m.buildUnstructured()
        self.assertTrue(isinstance(mcpy0,MEDCouplingUMesh))
        self.assertTrue(mcpy0.getNodalConnectivity().isEqual(DataArrayInt([31,3,2,8,9,31,1,0,6,7,-1,7,6,1,31,2,1,7,8,-1,2,1,-1,8,-1,7,31,5,4,10,11,-1,11,10,-1,5,31,4,3,9,10,-1,5,3,9])))
        self.assertTrue(mcpy0.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,14,26,36,45])))
        self.assertEqual(mcpy0.getAllGeoTypes(),[NORM_POLYHED])
        mcpy0.checkConsistencyLight()
        mcpy0.checkConsistency()
        mcpy1=mcpy0.convertIntoSingleGeoTypeMesh()
        self.assertTrue(mcpy1.isEqual(m,1e-12))
        #
        m_mrg=MEDCoupling1DGTUMesh.Merge1DGTUMeshes([m2,m,m2])
        self.assertTrue(m_mrg.getNodalConnectivityIndex().isEqual(DataArrayInt([0,8,19,23,31,40,44,52,63,72,80,88,99,103,111,120])))
        self.assertTrue(m_mrg.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,-1,7,6,1,2,1,7,8,-1,2,1,-1,8,-1,7,3,2,8,9,4,3,9,10,-1,5,3,9,5,4,10,11,-1,11,10,-1,5,27,26,32,33,25,24,30,31,-1,31,30,25,26,25,31,32,-1,26,25,-1,32,-1,31,29,28,34,35,-1,35,34,-1,29,28,27,33,34,-1,29,27,33,49,48,54,55,-1,55,54,49,50,49,55,56,-1,50,49,-1,56,-1,55,51,50,56,57,52,51,57,58,-1,53,51,57,53,52,58,59,-1,59,58,-1,53])))
        m_mrg2=MEDCoupling1DGTUMesh.Merge1DGTUMeshesOnSameCoords([m3,m,m3])
        self.assertTrue(m_mrg2.getNodalConnectivityIndex().isEqual(DataArrayInt([0,8,19,23,31,40,44,52,63,72,80,88,99,103,111,120])))
        self.assertTrue(m_mrg2.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,-1,7,6,1,2,1,7,8,-1,2,1,-1,8,-1,7,3,2,8,9,4,3,9,10,-1,5,3,9,5,4,10,11,-1,11,10,-1,5,3,2,8,9,1,0,6,7,-1,7,6,1,2,1,7,8,-1,2,1,-1,8,-1,7,5,4,10,11,-1,11,10,-1,5,4,3,9,10,-1,5,3,9,1,0,6,7,-1,7,6,1,2,1,7,8,-1,2,1,-1,8,-1,7,3,2,8,9,4,3,9,10,-1,5,3,9,5,4,10,11,-1,11,10,-1,5])))
        a,b=m_mrg2.getReverseNodalConnectivity()
        self.assertTrue(b.isEqual(DataArrayInt([0,3,15,24,33,39,48,54,66,75,84,93,99,99,99,99,99,99,99,99,99,99,99,99,99])))
        self.assertTrue(a.isEqual(DataArrayInt([0,6,10,0,0,1,1,6,6,7,7,10,10,11,11,1,1,2,5,7,7,11,11,12,2,3,3,5,9,9,12,13,13,3,4,8,9,13,14,3,4,4,8,8,9,13,14,14,0,0,6,6,10,10,0,0,1,1,6,6,7,7,10,10,11,11,1,1,2,5,7,7,11,11,12,2,3,3,5,9,9,12,13,13,3,4,4,8,8,9,13,14,14,4,4,8,8,14,14])))
        self.assertTrue(m_mrg2.fillCellIdsToKeepFromNodeIds([7],False).isEqual(DataArrayInt([0,1,6,7,10,11])))
        self.assertTrue(m_mrg2.fillCellIdsToKeepFromNodeIds([0,1,6,7],True).isEqual(DataArrayInt([0,6,10])))
        #
        self.assertTrue(m_mrg2.isPacked())
        self.assertEqual(120,m_mrg2.getNodalConnectivityIndex().popBackSilent())
        self.assertEqual(m_mrg2.getNumberOfCells(),14)
        m_mrg2.checkConsistency()
        self.assertTrue(not m_mrg2.isPacked())
        m_mrg4,b=m_mrg2.copyWithNodalConnectivityPacked()
        self.assertTrue(not b)
        m_mrg4.checkConsistency()
        self.assertEqual(m_mrg4.getNumberOfCells(),14)
        self.assertTrue(m_mrg4.getNodalConnectivityIndex().isEqual(m_mrg2.getNodalConnectivityIndex()))
        self.assertEqual(len(m_mrg4.getNodalConnectivity()),111)
        self.assertEqual(len(m_mrg2.getNodalConnectivity()),120)
        self.assertTrue(m_mrg4.getNodalConnectivity().isEqual(m_mrg2.getNodalConnectivity()[:111]))
        #
        m0=m_mrg2[:5]
        m1=m_mrg2[[5,6,7,8,9]]
        m2=m_mrg2[10:]
        self.assertTrue(m1.isEqualWithoutConsideringStr(m,1e-12))
        a,b=m.checkGeoEquivalWith(m0,12,1e-12)
        self.assertTrue(a.isEqual(DataArrayInt(o2n)))
        self.assertTrue(b is None)
        pass

    def testSwig2DADAreIncludedInMe1(self):
        a=DataArrayDouble(30) ; a.iota() ; a.rearrange(3)
        p=DataArrayInt([5,2,1,9])
        b,c=a.areIncludedInMe(a[p],1e-12)
        self.assertTrue(b)
        self.assertTrue(c.isEqual(p))
        d=a[p]
        d.setIJ(3,1,28.1)
        b,c=a.areIncludedInMe(d,1e-12)
        self.assertTrue(not b)
        self.assertTrue(c.isEqual(DataArrayInt([5,2,1,10])))
        pass

    def testSwig2DADesallocate1(self):
        d=DataArrayDouble([(1,2),(6,7),(6,8)]) ; d.setInfoOnComponents(["aa","bbb"])
        self.assertTrue(d.isAllocated())
        d.checkAllocated()
        self.assertEqual(d.getInfoOnComponents(),["aa","bbb"])
        ref=d.getHeapMemorySize()
        d.desallocate()
        self.assertEqual(ref-d.getHeapMemorySize(),6*8)
        self.assertTrue(not d.isAllocated())
        self.assertEqual(d.getInfoOnComponents(),["aa","bbb"])
        self.assertRaises(InterpKernelException,d.checkAllocated)
        #
        d=DataArrayInt([(1,2),(6,7),(6,8)]) ; d.setInfoOnComponents(["aa","bbb"])
        self.assertTrue(d.isAllocated())
        d.checkAllocated()
        self.assertEqual(d.getInfoOnComponents(),["aa","bbb"])
        ref=d.getHeapMemorySize()
        d.desallocate()
        self.assertEqual(ref-d.getHeapMemorySize(),6*4)
        self.assertTrue(not d.isAllocated())
        self.assertEqual(d.getInfoOnComponents(),["aa","bbb"])
        self.assertRaises(InterpKernelException,d.checkAllocated)
        pass

    def testSwig2IsPartStructured1(self):
        #dim 1
        d10=DataArrayInt([2,3,4,5,6,7,8,9,10,11])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d10,[13])
        self.assertTrue(a) ; self.assertEqual(b,[(2,12)])
        d11=DataArrayInt([2,3,4,5,6,7,8,10,9,11])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d11,[13])
        self.assertTrue(not a)
        self.assertRaises(InterpKernelException,MEDCouplingStructuredMesh.IsPartStructured,d10,[11])
        #dim 2
        st=[10,4]
        d20=DataArrayInt([1,2,3,4,11,12,13,14,21,22,23,24])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d20,st)
        self.assertTrue(a) ; self.assertEqual(b,[(1,5),(0,3)])
        self.assertEqual(12,MEDCouplingStructuredMesh.DeduceNumberOfGivenRangeInCompactFrmt(b))
        self.assertEqual(0,MEDCouplingStructuredMesh.DeduceNumberOfGivenRangeInCompactFrmt([(1,5),(1,3),(2,2)]))
        self.assertEqual(0,MEDCouplingStructuredMesh.DeduceNumberOfGivenRangeInCompactFrmt([(5,5),(3,3),(2,2)]))
        self.assertEqual(36,MEDCouplingStructuredMesh.DeduceNumberOfGivenStructure([3,2,6]))
        self.assertEqual(126,MEDCouplingStructuredMesh.DeduceNumberOfGivenStructure((3,7,6)))
        d20=DataArrayInt([1,2,3,4,12,11,13,14,21,22,23,24])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d20,st)
        self.assertTrue(not a)
        d20=DataArrayInt([1,2,3,4,11,12,13,15,21,22,23,24])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d20,st)
        self.assertTrue(not a)
        d21=DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d21,st)
        self.assertTrue(a) ; self.assertEqual(b,[(0,10),(0,4)])
        d22=DataArrayInt([1,2,3,4,11,12,13,14,21,22,23,24,31,32,33,34,41,42,43,44])
        self.assertRaises(InterpKernelException,MEDCouplingStructuredMesh.IsPartStructured,d22,st)
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d22,[10,5])
        self.assertTrue(a) ; self.assertEqual(b,[(1,5),(0,5)])
        #dim 3
        d30=DataArrayInt([11,12,13,14,21,22,23,24,51,52,53,54,61,62,63,64])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d30,[10,4,2])
        self.assertTrue(a) ; self.assertEqual(b,[(1,5),(1,3),(0,2)])
        d31=DataArrayInt([11,12,13,14,21,22,24,23,51,52,53,54,61,62,63,64])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d31,[10,4,2])
        self.assertTrue(not a)
        self.assertRaises(InterpKernelException,MEDCouplingStructuredMesh.IsPartStructured,d30,[10,4,1])
        pass

    def testSwig2PartStructured1(self):
        c=MEDCouplingCMesh() ; c.setName("toto")
        arr0=DataArrayDouble(10); arr0.iota()
        arr1=DataArrayDouble(4) ; arr1.iota(3)
        c.setCoords(arr0,arr1)
        self.assertEqual(c.getNodeGridStructure(),(10,4))
        self.assertEqual(c.getCellGridStructure(),(9,3))
        d20=DataArrayInt([1,2,3,4,10,11,12,13,19,20,21,22])
        self.assertEqual(27,c.getNumberOfCells())
        self.assertEqual(40,c.getNumberOfNodes())
        self.assertEqual(2,c.getMeshDimension())
        c.checkConsistencyLight()
        #
        arr2=MEDCouplingStructuredMesh.BuildExplicitIdsFrom([9,3],[(1,5),(0,3)])
        self.assertTrue(arr2.isEqual(DataArrayInt([1,2,3,4,10,11,12,13,19,20,21,22])))
        # CMesh
        c2=c.buildStructuredSubPart([(1,5),(0,3)])
        c2.checkConsistencyLight()
        self.assertTrue(isinstance(c2,MEDCouplingCMesh))
        self.assertEqual(12,c2.getNumberOfCells())
        self.assertEqual(20,c2.getNumberOfNodes())
        self.assertEqual(2,c2.getMeshDimension())
        self.assertEqual("toto",c2.getName())
        self.assertTrue(c2.getCoordsAt(0).isEqual(DataArrayDouble([1.,2.,3.,4.,5.]),1e-12))
        self.assertTrue(c2.getCoordsAt(1).isEqual(DataArrayDouble([3.,4.,5.,6.]),1e-12))
        #
        a,b=c.buildPartAndReduceNodes(d20)
        a.checkConsistencyLight()
        exp2=DataArrayInt([-1,0,1,2,3,4,-1,-1,-1,-1,-1,5,6,7,8,9,-1,-1,-1,-1,-1,10,11,12,13,14,-1,-1,-1,-1,-1,15,16,17,18,19,-1,-1,-1,-1])
        self.assertTrue(exp2.isEqual(b))
        self.assertTrue(isinstance(a,MEDCouplingCMesh))
        self.assertTrue(a.buildUnstructured().isEqual(c.buildUnstructured().buildPartAndReduceNodes(d20)[0],1e-12))
        # CurveLinearMesh
        c2=MEDCouplingCurveLinearMesh() ; c2.setName("toto")
        c2.setCoords(c.buildUnstructured().getCoords())
        c2.setNodeGridStructure([10,4])
        c2.checkConsistencyLight()
        a,b=c2.buildPartAndReduceNodes(d20)
        a.checkConsistencyLight()
        self.assertTrue(exp2.isEqual(b))
        self.assertTrue(isinstance(a,MEDCouplingCurveLinearMesh))
        self.assertTrue(a.buildUnstructured().isEqual(c2.buildUnstructured().buildPartAndReduceNodes(d20)[0],1e-12))
        pass

    def testSwig2FindPermutationFromFirstToSecond1(self):
        ids1=DataArrayInt([3,1,103,4,6,10,-7,205])
        ids2=DataArrayInt([-7,1,205,10,6,3,103,4])
        ids3=DataArrayInt.FindPermutationFromFirstToSecond(ids1,ids2)
        self.assertTrue(ids3.isEqual(DataArrayInt([5,1,6,7,4,3,0,2])))
        ids2ToTest=ids1.renumber(ids3)
        self.assertTrue(ids2ToTest.isEqual(ids2))
        self.assertRaises(InterpKernelException,DataArrayInt.FindPermutationFromFirstToSecond,DataArrayInt([3,1,103]),DataArrayInt([1,103]))
        self.assertRaises(InterpKernelException,DataArrayInt.FindPermutationFromFirstToSecond,DataArrayInt([3,1,103]),DataArrayInt([1,103,2]))
        self.assertRaises(InterpKernelException,DataArrayInt.FindPermutationFromFirstToSecond,DataArrayInt([3,1,103]),DataArrayInt([1,103,1]))
        self.assertTrue(DataArrayInt.FindPermutationFromFirstToSecond(DataArrayInt([]),DataArrayInt([])).empty())
        pass

    def testSwig2BugStructuredMeshGetNodeIdsOfCell1(self):
        m=MEDCouplingCMesh("mesh")
        coordsX=DataArrayDouble([0,1.1,2.2,3.3,4.4]) ; coordsX.setInfoOnComponents(["XX [m]"])
        coordsY=DataArrayDouble([0,1.7,3.4]) ; coordsY.setInfoOnComponents(["YYY [km]"])
        m.setCoords(coordsX,coordsY)
        self.assertEqual([2,3,8,7],m.getNodeIdsOfCell(2))
        self.assertEqual([3,4,9,8],m.getNodeIdsOfCell(3))
        self.assertEqual([7,8,13,12],m.getNodeIdsOfCell(6))
        self.assertEqual([8,9,14,13],m.getNodeIdsOfCell(7))
        pass

    def testSwig2ThrowOnDAIInvertN2O2ON2(self):
        p1=DataArrayInt([3,5,8])
        p2=DataArrayInt([0,3,4,5,6,7,8,9,10])
        p1.transformWithIndArr(p2.invertArrayN2O2O2N(11))
        self.assertTrue(p1.isEqual(DataArrayInt([1,3,6])))
        self.assertTrue(p2.invertArrayN2O2O2N(11).isEqual(DataArrayInt([0,-1,-1,1,2,3,4,5,6,7,8])))
        self.assertRaises(InterpKernelException,p2.invertArrayN2O2O2N,10)
        pass

    def testSwig2ComputeEffectiveNbOfNodesPerCell1(self):
        coords=DataArrayDouble([ 0.241310763507 , 0.0504777305619 , 0.0682283524903 , 0.252501053866 , -0.0625176732937 , 0.137272639894 ,
                 0.152262663601 , 0.241816569527 , 0.133812556197 , 0.18047750211 , -0.0789949051358 , 0.339098173401 ,
                 0.151741971857 , 0.238885278571 , 0.137715037333 , 0.242532155481 , -0.0928169086456 , 0.0678043417367 ,
                 0.240941965335 , -0.015461491464 , 0.0617186345825 , 0.24127650112 , 0.0499427876717 , 0.0679634099148 ,
                 -0.145828917428 , 0.206291632565 , 0.0310071927543 , 0.0125651775307 , 0.266262085828 , 0.105228430543 ,
                 -0.0994066533286 , 0.233224271238 , 0.0572213839567 , -0.0951345338317 , 0.234819509426 , 0.0592126284538 ,
                 0.136580574205 , -0.205486212579 , 0.0572866072014 , 0.0637270784978 , -0.168886355238 , 0.446614057077 ,
                 0.041337157151 , -0.213402568198 , 0.372407095999 , 0.0411601970268 , -0.202387875756 , 0.411334979491 ,
                 -0.108355701857 , 0.193636239335 , 0.204886756738 , 0.00639779029829 , 0.155296981517 , 0.252585892979 ,
                 0.0262473111702 , -0.112919732543 , 0.424286639249 ,-0.224103052733 , -0.139430015438 , -0.0122352295701 ,
                -0.0312760589481 , -0.274272003594 , 0.0323959636568 , -0.166663422532 , -0.217754445175 , 0.00392109070364 ,
                 -0.30586619777 , -0.0475168041091 , -0.0144585228182 , -0.280881480586 , 0.135571293538 , 0.00623923647986 ,
                 -0.25548538234 , 0.156819217766 , 0.0645277879769 , -0.131567009284 , 0.184133752309 , 0.206021802753 ,
                 -0.196204010965 , 0.151602971681 , 0.212974777736 , -0.183713879463 , 0.0802946639531 , 0.260115662599 ,
                 -0.244241178767 , -0.0738873389604 , 0.144590565817 , -0.155804057829 , -0.164892720025 , 0.210613950558 ,
                 -0.170950800428 , -0.215099334026 , 0.00610122860092 , -0.30552634869 , -0.0490020791904 , -0.0132786533145 ,
                 0.271831011884 , 0.15105657296 , 0.0230534827908 , 0.281919192283 , 0.0898544306288 , -0.0625201489143 ,
                 0.260240727276 , -0.0120688706637 , -0.0532316588626 , 0.244947737722 , 0.0197984684293 , 0.0309341209233 ,
                 0.23439631578 , 0.229825279875 , 0.0508520585381 , 0.160921316875 , 0.265078502128 , 0.121716560626 ,
                 -0.315088694175 , 0.0747700471918 , -0.245836615071 , -0.327728781776 , 0.0857114674649 , -0.239431905957 ,
                 -0.308385460634 , 0.145142997084 , -0.149886828433 , 0.0488236045164 , 0.309462801914 , 0.0849169148265 ,
                -0.0244964803395 , 0.33145611751 , -0.0476415818061 , 0.0060567994229 , 0.32418412014 , 0.0367779543812 ,
                 -0.0950221448063 , 0.236675326003 , 0.0572594453983 , 0.248723023186 , 0.0886648784791 , -0.176629430538 ,
                 0.116796984 , 0.256596599567 , -0.292863523603 , 0.118024552914 , 0.229154257843 , -0.34233232501 ,
                 0.217507892549 , -0.0417822335742 , -0.176771782888 , -0.224429321304 , 0.0125595300114 , -0.362064725588 ,
                 0.0937301100955 , -0.0500824832657 , -0.299713548444 , -0.244162220397 , 0.0383853931293 , -0.389856984411 ,
                 -0.0281989366102 , 0.097392811563 , -0.458244577284 , -0.385010847162 , 0.10122766194 , -0.140052859922 ,
                 -0.377936358012 , 0.110875172128 , -0.176207095463 , 0.244483045556 , -0.0991073977045 , 0.0575134372934 ,
                0.262605120167 , -0.100243191645 , -0.0495620806935 , 0.240306880972 , -0.136153701579 , -0.114745281696 ,
                 0.215763176129 , -0.0836766059189 , -0.183249640616 , 0.237870396603 , -0.132449578286 , -0.121598854639 ,
                 -0.0637683083097 , -0.27921020214 , -0.149112321992 , -0.0856211014977 , -0.2973233473 , -0.0446878139589 ,
                 0.104675342288 , -0.0625908305324 , -0.290346256534 , 0.0248264249186 , -0.247797708548 , -0.165830884019 ,
                 0.0719302438309 , -0.178468260473 , -0.211432157345 , 0.142871843159 , -0.208769948542 , 0.0454101128246 ,
                 0.167803379307 , -0.207851396623 , -0.088802726124 , 0.12868717152 , -0.230920439715 , 0.00760508389036 ,
                 -0.0372812069535 , -0.286740286332 , 0.00963701291166 ], 69, 3)
        connN = [ #polyhedron 0
            0 , 1 , 3 , 4 , 2 , -1 , 1 , 5 , 6 , 7 , 0 , -1 , 0 , 7 , 8 , 10 , 11 , 9 , 2 , -1 , 1 , 5 , 12 , 14 , 15 , 13 , 3 , -1 , 16 , 9 , 2 , 4 , 17 , -1
            , 4 , 3 , 13 , 18 , 17 , -1 , 5 , 6 , 19 , 21 , 20 , 12 , -1 , 6 , 7 , 8 , 23 , 22 , 19 , -1 , 23 , 24 , 10 , 8 , -1 , 25 , 11 , 9 , 16 , -1
            , 24 , 26 , 25 , 11 , 10 , -1 , 12 , 14 , 20 , -1 , 27 , 28 , 29 , 15 , 13 , 18 , -1 , 14 , 15 , 29 , 30 , 21 , 20 , -1 , 26 , 27 , 18 , 17 , 16 , 25 , -1
            , 22 , 19 , 21 , 30 , 31 , -1 , 22 , 31 , 28 , 27 , 26 , 24 , 23 , -1 , 31 , 30 , 29 , 28,
            # polyhedron 1
            0 , 7 , 8 , 10 , 11 , 9 , 2 , -1 , 32 , 0 , 7 , 35 , 34 , 33 , -1 , 32 , 0 , 2 , 37 , 36 , -1 , 35 , 7 , 8 , 40 , 39 , 38 , -1
            , 2 , 37 , 41 , 9 , -1 , 40 , 8 , 10 , 44 , 43 , 42 , -1 , 41 , 9 , 11 , 44 , 43 , -1 , 44 , 11 , 10 , -1 , 32 , 33 , 45 , 47 , 46 , 36 , -1
            , 33 , 34 , 48 , 45 , -1 , 35 , 34 , 48 , 50 , 49 , 38 , -1 , 41 , 43 , 42 , 46 , 36 , 37 , -1 , 38 , 39 , 51 , 49 , -1
            , 39 , 40 , 42 , 46 , 47 , 52 , 51 , -1 , 45 , 47 , 52 , 50 , 48 , -1 , 52 , 51 , 49 , 50,
            # polyhedron 2
            6 , 7 , 8 , 23 , 22 , 19 , -1 , 6 , 35 , 7 , -1 , 6 , 35 , 38 , 19 , -1 , 35 , 7 , 8 , 40 , 39 , 38 , -1 , 53 , 22 , 19 , 38 , 39 , 54 , -1
            , 23 , 53 , 54 , 40 , 8 , -1 , 53 , 22 , 23 , -1 , 39 , 54 , 40,
            # polyhedron 3
            35 , 34 , 48 , 50 , 49 , 38 , -1 , 6 , 35 , 34 , 56 , 55 , 5 , -1 , 6 , 35 , 38 , 19 , -1 , 34 , 56 , 57 , 59 , 58 , 48 , -1
            , 60 , 61 , 21 , 19 , 38 , 49 , -1 , 62 , 50 , 48 , 58 , -1 , 60 , 63 , 64 , 62 , 50 , 49 , -1 , 5 , 6 , 19 , 21 , 20 , 12 , -1
            , 55 , 5 , 12 , 65 , -1 , 66 , 67 , 65 , 55 , 56 , 57 , -1 , 63 , 66 , 57 , 59 , 64 , -1 , 64 , 62 , 58 , 59 , -1
            , 60 , 63 , 66 , 67 , 68 , 61 , -1 , 61 , 68 , 20 , 21 , -1 , 67 , 68 , 20 , 12 , 65]
        meshN=MEDCouplingUMesh.New()
        meshN.setName("ForBary")
        meshN.setMeshDimension(3) ; meshN.setCoords(coords)
        meshN.allocateCells(4)
        meshN.insertNextCell(NORM_POLYHED,113,connN);
        meshN.insertNextCell(NORM_POLYHED,99,connN[113:])
        meshN.insertNextCell(NORM_POLYHED,43,connN[212:])
        meshN.insertNextCell(NORM_POLYHED,92,connN[255:])
        d=meshN.computeEffectiveNbOfNodesPerCell()
        e=meshN.computeNbOfNodesPerCell()
        self.assertTrue(d.isEqual(DataArrayInt([32,28,12,26])))
        self.assertTrue(e.isEqual(DataArrayInt([96,84,36,78])))
        m0=MEDCoupling1DGTUMesh(meshN)
        c=MEDCouplingCMesh()
        arr=DataArrayDouble(3) ; arr.iota(10)
        c.setCoords(arr,arr,arr)
        m10=c.buildUnstructured()
        m11=c.build1SGTUnstructured()
        m12=MEDCoupling1SGTUMesh.New(m10)
        self.assertTrue(m12.isEqual(m11,1e-12))
        m12.setCoords(m0.getCoords()) # m12 is not OK geometrically but the aim of the test is only connectivity values
        m3=MEDCoupling1GTUMesh.AggregateOnSameCoordsToUMesh([m12,m0])
        m3.checkConsistencyLight()
        self.assertEqual(m3.getCoords().getHiddenCppPointer(),m12.getCoords().getHiddenCppPointer())
        self.assertTrue(m3.getNodalConnectivity().isEqual(DataArrayInt([18,1,0,3,4,10,9,12,13,18,2,1,4,5,11,10,13,14,18,4,3,6,7,13,12,15,16,18,5,4,7,8,14,13,16,17,18,10,9,12,13,19,18,21,22,18,11,10,13,14,20,19,22,23,18,13,12,15,16,22,21,24,25,18,14,13,16,17,23,22,25,26,31,0,1,3,4,2,-1,1,5,6,7,0,-1,0,7,8,10,11,9,2,-1,1,5,12,14,15,13,3,-1,16,9,2,4,17,-1,4,3,13,18,17,-1,5,6,19,21,20,12,-1,6,7,8,23,22,19,-1,23,24,10,8,-1,25,11,9,16,-1,24,26,25,11,10,-1,12,14,20,-1,27,28,29,15,13,18,-1,14,15,29,30,21,20,-1,26,27,18,17,16,25,-1,22,19,21,30,31,-1,22,31,28,27,26,24,23,-1,31,30,29,28,31,0,7,8,10,11,9,2,-1,32,0,7,35,34,33,-1,32,0,2,37,36,-1,35,7,8,40,39,38,-1,2,37,41,9,-1,40,8,10,44,43,42,-1,41,9,11,44,43,-1,44,11,10,-1,32,33,45,47,46,36,-1,33,34,48,45,-1,35,34,48,50,49,38,-1,41,43,42,46,36,37,-1,38,39,51,49,-1,39,40,42,46,47,52,51,-1,45,47,52,50,48,-1,52,51,49,50,31,6,7,8,23,22,19,-1,6,35,7,-1,6,35,38,19,-1,35,7,8,40,39,38,-1,53,22,19,38,39,54,-1,23,53,54,40,8,-1,53,22,23,-1,39,54,40,31,35,34,48,50,49,38,-1,6,35,34,56,55,5,-1,6,35,38,19,-1,34,56,57,59,58,48,-1,60,61,21,19,38,49,-1,62,50,48,58,-1,60,63,64,62,50,49,-1,5,6,19,21,20,12,-1,55,5,12,65,-1,66,67,65,55,56,57,-1,63,66,57,59,64,-1,64,62,58,59,-1,60,63,66,67,68,61,-1,61,68,20,21,-1,67,68,20,12,65])))
        self.assertTrue(m3.getNodalConnectivityIndex().isEqual(DataArrayInt([0,9,18,27,36,45,54,63,72,186,286,330,423])))
        pass

    def testSwig2Tetrahedrize1(self):
        d=DataArrayInt([0,3,6,10,14,20])
        d2=d.buildExplicitArrOfSliceOnScaledArr(slice(0,5,2))
        self.assertTrue(d2.isEqual(DataArrayInt([0,0,0, 2,2,2,2, 4,4,4,4,4,4])))
        m=MEDCouplingUMesh("Penta6",3)
        m.setCoords(DataArrayDouble([0,0,0,0,1,0,1,0,0,0,0,2,0,1,2,1,0,2],6,3)) ; m.getCoords().setInfoOnComponents(["X","YY","ZZZ"])
        m.allocateCells()
        m.insertNextCell(NORM_PENTA6,[1,2,0,4,5,3])
        st=m.getCoords().getHiddenCppPointer()
        c,a,b=m.tetrahedrize(PLANAR_FACE_5)
        c.checkConsistency()
        self.assertTrue(a.isEqual(DataArrayInt([0,0,0])))
        self.assertEqual(0,b)
        self.assertEqual(m.getCoords().getHiddenCppPointer(),c.getCoords().getHiddenCppPointer())
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([1,2,0,4,4,3,5,0,5,0,2,4])))
        del m,c
        #
        m2=MEDCouplingUMesh("octa12",3)
        coords=DataArrayDouble([1.,0.,0.,0.5,0.8660254037844386,0.,-0.5,0.8660254037844387,0.,-1.,1.2246467991473532e-16,0.,-0.5,-0.8660254037844384,0.,0.5,-0.866025403784439,0.,1.,0.,2.,0.5,0.8660254037844386,2.,-0.5,0.8660254037844387,2.,-1.,1.2246467991473532e-16,2.,-0.5,-0.8660254037844384,2.,0.5,-0.866025403784439,2.0],12,3)
        m2.setCoords(coords)
        m2.allocateCells()
        m2.insertNextCell(NORM_HEXGP12,[3,2,1,0,5,4,9,8,7,6,11,10])
        c,a,b=m2.tetrahedrize(PLANAR_FACE_5)
        c.checkConsistency()
        self.assertTrue(a.isEqual(DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0])))
        self.assertEqual(0,b)
        self.assertEqual(c.getCoords().getHiddenCppPointer(),coords.getHiddenCppPointer())
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([3,2,4,9,9,10,8,4,8,4,2,9,2,5,4,8,8,10,11,4,11,4,5,8,2,1,5,8,8,11,7,5,7,5,1,8,1,0,5,7,7,11,6,5,6,5,0,7])))
        del m2,coords,c
        #
        coords=DataArrayDouble([0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,2.,1.,0.,2.,1.,1.,2.,0.,1.,2.],8,3) ; coords.setInfoOnComponents(["X","YY","ZZZ"])
        m3=MEDCouplingUMesh("hexa8",3)
        m3.setCoords(coords)
        m3.allocateCells(0)
        m3.insertNextCell(NORM_HEXA8,[3,2,1,0,7,6,5,4])
        st=m3.getCoords().getHiddenCppPointer()
        c,a,b=m3.tetrahedrize(PLANAR_FACE_5)
        c.checkConsistency()
        a.isEqual(DataArrayInt([0,0,0,0,0]))
        self.assertEqual(0,b)
        self.assertEqual(m3.getCoords().getHiddenCppPointer(),coords.getHiddenCppPointer())
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([3,6,2,1,3,7,6,4,3,0,4,1,6,4,5,1,3,6,1,4])))
        #
        m4=MEDCouplingUMesh("hexa8",3)
        m4.setCoords(coords)
        m4.allocateCells(0)
        m4.insertNextCell(NORM_HEXA8,[3,2,1,0,7,6,5,4])
        c,a,b=m4.tetrahedrize(PLANAR_FACE_6)
        c.checkConsistency()
        a.isEqual(DataArrayInt([0,0,0,0,0,0]))
        self.assertEqual(0,b)
        self.assertEqual(c.getCoords().getHiddenCppPointer(),coords.getHiddenCppPointer())
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([3,6,2,5,3,2,1,5,3,7,6,5,3,4,7,5,3,1,0,5,3,0,4,5])))
        #
        m4=MEDCouplingUMesh("hexa8",3)
        m4.setCoords(coords)
        m4.allocateCells(0)
        m4.insertNextCell(NORM_HEXA8,[3,2,1,0,7,6,5,4])
        st=m4.getCoords().getHiddenCppPointer()
        c,a,b=m4.tetrahedrize(GENERAL_24)
        c.checkConsistency()
        a.isEqual(DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))
        self.assertEqual(7,b)
        self.assertTrue(c.getCoords().getHiddenCppPointer()!=coords.getHiddenCppPointer())
        self.assertTrue(c.getCoords()[:8].isEqual(coords,0))
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([3,7,8,14,7,6,8,14,6,2,8,14,2,3,8,14,3,2,9,14,2,1,9,14,1,0,9,14,0,3,9,14,3,0,10,14,0,4,10,14,4,7,10,14,7,3,10,14,2,6,11,14,6,5,11,14,5,1,11,14,1,2,11,14,7,4,12,14,4,5,12,14,5,6,12,14,6,7,12,14,1,5,13,14,5,4,13,14,4,0,13,14,0,1,13,14])))
        m4CoordsExp=DataArrayDouble([0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,2.,1.,0.,2.,1.,1.,2.,0.,1.,2.,0.5,1.,1.,0.5,0.5,0.,0.,0.5,1.,1.,0.5,1.,0.5,0.5,2.,0.5,0.,1.,0.5,0.5,1.],15,3)
        m4CoordsExp.setInfoOnComponents(["X","YY","ZZZ"])
        self.assertTrue(c.getCoords().isEqual(m4CoordsExp,1e-12))
        self.assertAlmostEqual(2.,c.getMeasureField(False).accumulate()[0],12)
        #
        m6=MEDCouplingUMesh("hexa8",3)
        m6.setCoords(coords)
        m6.allocateCells(0)
        m6.insertNextCell(NORM_HEXA8,[3,2,1,0,7,6,5,4])
        st=m6.getCoords().getHiddenCppPointer()
        c,a,b=m6.tetrahedrize(GENERAL_48)
        c.checkConsistency()
        a.isEqual(DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))
        self.assertEqual(19,b)
        self.assertTrue(c.getCoords().getHiddenCppPointer()!=coords.getHiddenCppPointer())
        self.assertTrue(c.getCoords()[:8].isEqual(coords,0))
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([3,20,8,26,3,8,21,26,3,9,20,26,3,22,9,26,3,21,12,26,3,12,22,26,8,10,2,23,8,2,13,23,8,20,10,23,8,26,20,23,8,13,21,23,8,21,26,23,12,26,21,25,12,21,16,25,12,22,26,25,12,17,22,25,12,16,0,25,12,0,17,25,21,23,13,18,21,13,1,18,21,26,23,18,21,25,26,18,21,1,16,18,21,16,25,18,9,11,20,24,9,20,26,24,9,7,11,24,9,14,7,24,9,26,22,24,9,22,14,24,20,6,10,15,20,10,23,15,20,11,6,15,20,24,11,15,20,23,26,15,20,26,24,15,22,24,26,19,22,26,25,19,22,14,24,19,22,4,14,19,22,25,17,19,22,17,4,19,26,15,23,5,26,23,18,5,26,24,15,5,26,19,24,5,26,18,25,5,26,25,19,5])))
        m6CoordsExp=DataArrayDouble([0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,2.,1.,0.,2.,1.,1.,2.,0.,1.,2.,0.5,1.,0.,0.,1.,1.,1.,1.,1.,0.5,1.,2.,0.,0.5,0.,1.,0.5,0.,0.,0.5,2.,1.,0.5,2.,0.5,0.,0.,0.,0.,1.,1.,0.,1.,0.5,0.,2.,0.5,1.,1.,0.5,0.5,0.,0.,0.5,1.,1.,0.5,1.,0.5,0.5,2.,0.5,0.,1.,0.5,0.5,1.],27,3)
        m6CoordsExp.setInfoOnComponents(["X","YY","ZZZ"])
        self.assertTrue(c.getCoords().isEqual(m6CoordsExp,1e-12))
        self.assertAlmostEqual(2.,c.getMeasureField(False).accumulate()[0],12)
        #
        m7=MEDCouplingUMesh("polyhed",3)
        coords=DataArrayDouble([1.,0.,0.,0.5,0.8660254037844386,0.,-0.5,0.8660254037844387,0.,-1.,0.,0.,-0.5,-0.8660254037844384,0.,0.5,-0.866025403784439,0.,1.,0.,2.,0.5,0.8660254037844386,2.,-0.5,0.8660254037844387,2.,-1.,0.,2.,-0.5,-0.8660254037844384,2.,0.5,-0.866025403784439,2.0],12,3) ; coords.setInfoOnComponents(["X","YY","ZZZ"])
        m7.setCoords(coords)
        m7.allocateCells()
        m7.insertNextCell(NORM_POLYHED,[3,2,1,0,5,4,-1,9,10,11,6,7,8,-1,3,9,8,2,-1,2,8,7,1,-1,1,7,6,0,-1,0,6,11,5,-1,5,11,10,4,-1,4,10,9,3])
        c,a,b=m7.tetrahedrize(PLANAR_FACE_5)
        c.checkConsistency()
        self.assertTrue(a.isEqual(DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])))
        self.assertEqual(9,b)
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([3,2,12,20,2,1,12,20,1,0,12,20,0,5,12,20,5,4,12,20,4,3,12,20,9,10,13,20,10,11,13,20,11,6,13,20,6,7,13,20,7,8,13,20,8,9,13,20,3,9,14,20,9,8,14,20,8,2,14,20,2,3,14,20,2,8,15,20,8,7,15,20,7,1,15,20,1,2,15,20,1,7,16,20,7,6,16,20,6,0,16,20,0,1,16,20,0,6,17,20,6,11,17,20,11,5,17,20,5,0,17,20,5,11,18,20,11,10,18,20,10,4,18,20,4,5,18,20,4,10,19,20,10,9,19,20,9,3,19,20,3,4,19,20])))
        self.assertAlmostEqual(5.196152422706635,c.getMeasureField(False).accumulate()[0],12)
        m7CoordsExp=DataArrayDouble([1.0,0.0,0.0,0.5,0.8660254037844386,0.0,-0.5,0.8660254037844387,0.0,-1.0,0.,0.0,-0.5,-0.8660254037844384,0.0,0.5,-0.866025403784439,0.0,1.0,0.0,2.0,0.5,0.8660254037844386,2.0,-0.5,0.8660254037844387,2.0,-1.0,0.,2.0,-0.5,-0.8660254037844384,2.0,0.5,-0.866025403784439,2.0,0.0,0.0,0.0,0.0,0.,2.0,-0.75,0.4330127018922194,1.0,0.0,0.8660254037844386,1.0,0.75,0.4330127018922193,1.0,0.75,-0.4330127018922195,1.0,0.0,-0.8660254037844387,1.0,-0.75,-0.4330127018922191,1.0,0.0,0.,1.0],21,3)
        m7CoordsExp.setInfoOnComponents(["X","YY","ZZZ"])
        self.assertTrue(c.getCoords().isEqual(m7CoordsExp,1e-12))
        del m7,coords,c
        #
        coords=DataArrayDouble([0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,2.,1.,0.,2.,1.,1.,2.,0.,1.,2.],8,3) ; coords.setInfoOnComponents(["X","YY","ZZZ"])
        m8=MEDCouplingUMesh("pyra5",3)
        m8.setCoords(coords)
        m8.allocateCells(0)
        m8.insertNextCell(NORM_PYRA5,[3,2,1,0,7])
        st=m8.getCoords().getHiddenCppPointer()
        c,a,b=m8.tetrahedrize(PLANAR_FACE_5)
        self.assertEqual(m8.getCoords().getHiddenCppPointer(),coords.getHiddenCppPointer())
        c.checkConsistency()
        self.assertTrue(a.isEqual(DataArrayInt([0,0])))
        self.assertEqual(0,b)
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([3,2,1,7,3,1,0,7])))
        self.assertAlmostEqual(0.6666666666666667,c.getMeasureField(False).accumulate()[0],12)
        pass

    def testDualMesh3D1(self):
        arr=DataArrayDouble(2) ; arr.iota()
        c=MEDCouplingCMesh() ; c.setCoords(arr,arr,arr)
        m=c.buildUnstructured()
        t=m.tetrahedrize(PLANAR_FACE_5)[0]
        d=t.computeDualMesh()
        self.assertTrue(d.getNodalConnectivityIndex().isEqual(DataArrayInt([0,29,118,207,236,325,354,383,472])))
        self.assertTrue(d.getNodalConnectivity().isEqual(DataArrayInt([26,11,42,8,-1,25,8,42,10,-1,29,10,42,11,-1,0,26,8,25,-1,0,25,10,29,-1,0,29,11,26,24,9,42,8,-1,26,8,42,11,-1,27,11,42,9,-1,1,24,8,26,-1,1,26,11,27,-1,30,13,43,12,-1,24,12,43,15,-1,32,15,43,13,-1,1,30,12,24,-1,1,32,13,30,-1,35,17,44,16,-1,32,16,44,19,-1,27,19,44,17,-1,1,35,16,32,-1,1,27,17,35,-1,24,15,46,9,-1,27,9,46,19,-1,32,19,46,15,27,9,42,11,-1,29,11,42,10,-1,28,10,42,9,-1,2,29,10,28,-1,2,27,11,29,-1,27,17,44,19,-1,38,19,44,18,-1,37,18,44,17,-1,2,37,17,27,-1,2,38,18,37,-1,28,21,45,23,-1,41,23,45,22,-1,38,22,45,21,-1,2,41,22,38,-1,2,28,23,41,-1,27,19,46,9,-1,28,9,46,21,-1,38,21,46,19,35,16,44,17,-1,36,18,44,16,-1,37,17,44,18,-1,3,36,16,35,-1,3,35,17,37,-1,3,37,18,36,24,8,42,9,-1,25,10,42,8,-1,28,9,42,10,-1,4,25,8,24,-1,4,28,10,25,-1,24,15,43,12,-1,31,12,43,14,-1,34,14,43,15,-1,4,24,12,31,-1,4,31,14,34,-1,34,21,45,20,-1,40,20,45,23,-1,28,23,45,21,-1,4,34,20,40,-1,4,40,23,28,-1,24,9,46,15,-1,28,21,46,9,-1,34,15,46,21,30,12,43,13,-1,31,14,43,12,-1,33,13,43,14,-1,5,31,12,30,-1,5,30,13,33,-1,5,33,14,31,40,23,45,20,-1,39,20,45,22,-1,41,22,45,23,-1,6,40,20,39,-1,6,39,22,41,-1,6,41,23,40,32,13,43,15,-1,34,15,43,14,-1,33,14,43,13,-1,7,33,13,32,-1,7,34,14,33,-1,32,19,44,16,-1,36,16,44,18,-1,38,18,44,19,-1,7,32,16,36,-1,7,36,18,38,-1,34,20,45,21,-1,39,22,45,20,-1,38,21,45,22,-1,7,39,20,34,-1,7,38,22,39,-1,32,15,46,19,-1,38,19,46,21,-1,34,21,46,15])))
        self.assertTrue(d.getCoords().isEqual(DataArrayDouble([0.,0.,0.,1.,0.,0.,0.,1.,0.,1.,1.,0.,0.,0.,1.,1.,0.,1.,0.,1.,1.,1.,1.,1.,0.3333333333333333,0.,0.3333333333333333,0.3333333333333333,0.3333333333333333,0.3333333333333333,0.,0.3333333333333333,0.3333333333333333,0.3333333333333333,0.3333333333333333,0.,0.6666666666666666,0.,0.6666666666666666,1.,0.3333333333333333,0.6666666666666666,0.6666666666666666,0.3333333333333333,1.,0.6666666666666666,0.3333333333333333,0.6666666666666666,1.,0.6666666666666666,0.3333333333333333,0.6666666666666666,0.6666666666666666,0.,0.6666666666666666,1.,0.3333333333333333,0.6666666666666666,0.6666666666666666,0.3333333333333333,0.3333333333333333,0.6666666666666666,1.,0.3333333333333333,0.6666666666666666,0.6666666666666666,0.3333333333333333,1.,0.6666666666666666,0.,0.6666666666666666,0.6666666666666666,0.5,0.,0.5,0.,0.,0.5,0.5,0.,0.,0.5,0.5,0.,0.,0.5,0.5,0.,0.5,0.,1.,0.,0.5,0.5,0.,1.,1.,0.5,0.5,1.,0.5,1.,0.5,0.5,1.,1.,0.5,0.,1.,1.,0.5,0.5,1.,0.,0.5,1.,0.5,0.5,1.,1.,0.,0.5,1.,0.,1.,0.5,0.25,0.25,0.25,0.75,0.25,0.75,0.75,0.75,0.25,0.25,0.75,0.75,0.5,0.5,0.5],47,3),1e-12))
        self.assertAlmostEqual(1.,d.getMeasureField(False).accumulate()[0],1e-13)
        pass

    def testDualMesh2D1(self):
        arr=DataArrayDouble(5) ; arr.iota()
        c=MEDCouplingCMesh() ; c.setCoords(arr,arr)
        m=c.buildUnstructured()
        m.simplexize(0)
        t=MEDCoupling1SGTUMesh(m)
        d=t.computeDualMesh()
        self.assertTrue(d.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,12,20,28,34,42,54,66,78,86,94,106,118,130,138,146,158,170,182,190,196,204,212,220,224])))
        self.assertTrue(d.getNodalConnectivity().isEqual(DataArrayInt([26,81,25,0,25,81,27,82,29,83,30,1,30,83,31,84,33,85,34,2,34,85,35,86,37,87,38,3,38,87,39,88,41,4,27,81,26,5,42,89,28,82,29,82,28,89,43,90,45,91,32,84,31,83,33,84,32,91,46,92,48,93,36,86,35,85,37,86,36,93,49,94,51,95,40,88,39,87,41,88,40,95,52,96,54,9,43,89,42,10,55,97,44,90,45,90,44,97,56,98,58,99,47,92,46,91,48,92,47,99,59,100,61,101,50,94,49,93,51,94,50,101,62,102,64,103,53,96,52,95,54,96,53,103,65,104,67,14,56,97,55,15,68,105,57,98,58,98,57,105,69,106,71,107,60,100,59,99,61,100,60,107,72,108,74,109,63,102,62,101,64,102,63,109,75,110,77,111,66,104,65,103,67,104,66,111,78,112,80,19,69,105,68,20,70,106,71,106,70,21,73,108,72,107,74,108,73,22,76,110,75,109,77,110,76,23,79,112,78,111,80,112,79,24])))
        self.assertTrue(d.getCoords().isEqual(DataArrayDouble([0.,0.,1.,0.,2.,0.,3.,0.,4.,0.,0.,1.,1.,1.,2.,1.,3.,1.,4.,1.,0.,2.,1.,2.,2.,2.,3.,2.,4.,2.,0.,3.,1.,3.,2.,3.,3.,3.,4.,3.,0.,4.,1.,4.,2.,4.,3.,4.,4.,4.,0.5,0.,0.,0.5,0.5,0.5,0.5,1.,1.,0.5,1.5,0.,1.5,0.5,1.5,1.,2.,0.5,2.5,0.,2.5,0.5,2.5,1.,3.,0.5,3.5,0.,3.5,0.5,3.5,1.,4.,0.5,0.,1.5,0.5,1.5,0.5,2.,1.,1.5,1.5,1.5,1.5,2.,2.,1.5,2.5,1.5,2.5,2.,3.,1.5,3.5,1.5,3.5,2.,4.,1.5,0.,2.5,0.5,2.5,0.5,3.,1.,2.5,1.5,2.5,1.5,3.,2.,2.5,2.5,2.5,2.5,3.,3.,2.5,3.5,2.5,3.5,3.,4.,2.5,0.,3.5,0.5,3.5,0.5,4.,1.,3.5,1.5,3.5,1.5,4.,2.,3.5,2.5,3.5,2.5,4.,3.,3.5,3.5,3.5,3.5,4.,4.,3.5,0.3333333333333333,0.3333333333333333,0.6666666666666666,0.6666666666666666,1.3333333333333333,0.3333333333333333,1.6666666666666665,0.6666666666666666,2.333333333333333,0.3333333333333333,2.6666666666666665,0.6666666666666666,3.333333333333333,0.3333333333333333,3.6666666666666665,0.6666666666666666,0.3333333333333333,1.3333333333333333,0.6666666666666666,1.6666666666666665,1.3333333333333333,1.3333333333333333,1.6666666666666665,1.6666666666666665,2.333333333333333,1.3333333333333333,2.6666666666666665,1.6666666666666665,3.333333333333333,1.3333333333333333,3.6666666666666665,1.6666666666666665,0.3333333333333333,2.333333333333333,0.6666666666666666,2.6666666666666665,1.3333333333333333,2.333333333333333,1.6666666666666665,2.6666666666666665,2.333333333333333,2.333333333333333,2.6666666666666665,2.6666666666666665,3.333333333333333,2.333333333333333,3.6666666666666665,2.6666666666666665,0.3333333333333333,3.333333333333333,0.6666666666666666,3.6666666666666665,1.3333333333333333,3.333333333333333,1.6666666666666665,3.6666666666666665,2.333333333333333,3.333333333333333,2.6666666666666665,3.6666666666666665,3.333333333333333,3.333333333333333,3.6666666666666665,3.6666666666666665],113,2),1e-12))
        self.assertAlmostEqual(16.,d.getMeasureField(False).accumulate()[0],1e-13)
        pass

    def testSwig2LoadBalanceBBox1(self):
        arr=DataArrayDouble(5) ; arr.iota()
        t=MEDCouplingCMesh() ; t.setCoords(arr,arr)
        arr=DataArrayDouble(16) ; arr.iota() ; arr*=2./15
        s=MEDCouplingCMesh() ; s.setCoords(arr,arr[:]) ; s.translate([2.,1.])
        #
        s1=s.build1SGTUnstructured()
        t1=t.build1SGTUnstructured()
        w=MEDCouplingPointSet.ComputeNbOfInteractionsWithSrcCells(s1,t1,1e-12)
        wExp=DataArrayInt([0,0,0,0,0,0,64,64,0,0,64,64,0,0,0,0])
        self.assertTrue(w.isEqual(wExp))
        slcs=w.splitInBalancedSlices(4)
        self.assertEqual(len(slcs),4)
        self.assertEqual(slcs,[slice(0,7,1),slice(7,8,1),slice(8,11,1),slice(11,16,1)])
        bbs=s1.getBoundingBoxForBBTree()
        bbt=t1.getBoundingBoxForBBTree()
        self.assertTrue(bbt.computeNbOfInteractionsWith(bbs,1e-12).isEqual(wExp))
        pass

    def testKrSpatialDiscretization2(self):
        srcPointCoordsXY=DataArrayDouble([0.8401877171547095,0.39438292681909304,0.7830992237586059,0.7984400334760733,0.9116473579367843,0.19755136929338396,0.335222755714889,0.768229594811904,0.2777747108031878,0.5539699557954305,0.47739705186216025,0.6288709247619244,0.36478447279184334,0.5134009101956155,0.9522297251747128,0.9161950680037007,0.6357117279599009,0.7172969294326831,0.14160255535580338,0.6069688762570586,0.01630057162432958,0.24288677062973696,0.13723157678601872,0.8041767542269904,0.15667908925408455,0.4009443942461835,0.12979044678145574,0.10880880202576929,0.998924518003559,0.21825690531090688,0.5129323944043984,0.8391122346926072,0.6126398325956612,0.29603161769734304,0.6375522677030192,0.5242871900667843,0.493582986990727,0.9727750238835695,0.29251678441302703,0.7713576977939148,0.5267449792133388,0.7699138362751873,0.4002286220901779,0.8915294520051822,0.2833147460051415,0.3524583472648907,0.8077245200088827,0.9190264739650424,0.06975527623191256,0.9493270753646861,0.5259953502221011,0.08605584785624214,0.19221384599442307,0.6632269270081198,0.8902326025488938,0.3488929352485076,0.06417132078864207,0.02002304886468828,0.4577017372742769,0.06309583832653977,0.23827995417559517,0.9706341316786754,0.9022080734848082,0.8509197867712563,0.2666657493760184,0.5397603407221662,0.3752069763723793,0.7602487363667454,0.5125353641400744,0.6677237607854063,0.5316064341606602,0.039280343353413204,0.4376375965949323,0.9318350562508382,0.9308097953585953,0.7209523430657351,0.28429340305006756,0.7385343149018168,0.6399788165651163,0.3540486797476414,0.687861390266503,0.16597416632155615,0.4401045276038835,0.880075236260926,0.829201093329676,0.3303371296871161,0.22896817104377232,0.8933724145839793,0.35036017855180435,0.6866699083180492,0.9564682529105192,0.5886401331930609,0.6573040395310633,0.8586763259296661,0.4395599194986559,0.9239697889070817,0.39843666665183225,0.8147668963366965,0.6842185252738271,0.9109720307919067,0.4824906566564416,0.21582495896882609,0.9502523741453198,0.9201282537170352,0.14766001475400292,0.8810621695039152,0.641080596317109,0.43195341826973177,0.6195964839400707,0.281059412416564,0.7860020980173732,0.3074578737409124,0.44703357920378145,0.22610662515559543,0.18753310953617705,0.27623467206779617,0.5564437553083728,0.4165012805799494,0.16960708618611428,0.9068039338601771,0.10317118843233734,0.1260753390966334,0.49544406658757667,0.7604752284290619,0.9847516650262995,0.9350039865518939,0.6844450168704823,0.3831883312124705,0.7497708824229291,0.36866354167864823,0.2941603620043771,0.2322615386137094,0.5844885006474743,0.24441273568403568,0.15238979186508328,0.7321485158671385,0.12547490472228962,0.7934703881821923,0.164101933671209,0.7450713891280216,0.07452980059875632,0.9501040316885822,0.05252926240327268,0.5215633798025378,0.1762106563785163,0.24006237240511102,0.797798051870334,0.732654411686889,0.6565636529850605,0.9674051385221095,0.6394583455470663,0.7597348418830591,0.09348047715308166,0.13490241166898162,0.5202100698464597,0.07823214171371988,0.06990639775521419,0.2046550862512808,0.4614204733918516,0.8196772801781433,0.5733186283955903,0.7555808353962288,0.05193881879185271,0.1578071285774033,0.9999935710802644,0.204328610656936,0.8899556444445419,0.12546847580255405,0.9977989993047895,0.054057577650089554,0.8705398649305757,0.07232879943788462,0.004161608873010431,0.9230691273338484,0.5938921792404224,0.180372265717188,0.16313149927329806,0.3916902306450951,0.9130266774040771,0.8196951527240198,0.35909536870154335,0.552485022485482,0.5794299941414176,0.452575845854625,0.687387434620125,0.09964006352221597,0.5308079880340062,0.7572938323753392,0.30429514977349675,0.9922284614258579,0.5769711125534824,0.877613778169087,0.7478092963564253,0.6289099313453351,0.03542090674649035,0.7478028669710285,0.8332385420022712,0.9253765511910322,0.8732713427735824,0.8310375408413995],100,2)
        srcFieldValsOnPoints=DataArrayDouble([0.7643742528498438,-0.023507696856211995,1.1082895131907775,0.6299357452572031,0.8892623544912389,0.72212114810697,0.9196401044320336,-0.759961711221917,0.40801932617748826,0.8441134300809151,0.982483804252809,0.6752368914020778,0.9924403977479798,1.1063334970204484,0.9403055261137516,0.3624481886322733,1.1344772505996308,0.7522965618948239,0.17077741651388564,0.6504551671311436,0.45843479588425423,0.41098905950326753,1.0681420394050904,-0.3483587903820091,0.5620151050607809,1.384969776596035,0.7948875141132845,0.7931192000237167,1.062498042490183,1.3709072529577366,0.44929346605311893,-0.4469683401788374,0.9035857424514101,0.6137249300593463,0.6355610879026966,1.4318174829507697,0.3097567072129551,-0.20515052260807165,0.6922559820922779,1.0341638749443423,1.3072652153341024,0.38511367353000436,0.9160514929274943,0.54513408530581,0.722252267913328,0.06684522818576251,0.10571899758067793,0.3193844999960903,0.5213532270828706,-0.04834998649603944,1.2408805068350615,-0.7632951295676795,0.5980054665011202,0.9064738717547436,1.1541070755096696,1.008234260272265,1.2225806960553827,1.0788560195121106,0.9818990282104452,0.5621951325841853,1.0796757508374188,0.5082872315589883,-0.9153702001062469,0.9560418838920791,0.9251098559152824,1.1603063610984021,1.2122303611181837,0.7379539363312343,0.6877611899207183,0.723966552446608,0.5596025827162566,0.8849725005989729,1.0908363665075547,0.08956512916455672,-0.10247645571248344,0.3236718069555875,1.069478546398975,1.3900071080692746,1.0322398863403262,0.45315515354558034,0.4249870238786733,1.030226761858634,0.974024629584669,1.2838885424020365,1.3451943506525155,1.4029933267831995,0.6025539675442462,1.2947650597767038,1.0006061239483002,-0.4017336259949164,0.8771165113201297,0.9158909024218246,1.403798605551443,0.4742904006425974,0.3671787905896653,0.20646491720419674,0.40739337434288925,0.7341932402033597,-0.4295893651836911,-0.3187777570661546],100,1)
        targetPointCoordsXY=DataArrayDouble([-0.5,-0.5,-0.5,-0.35,-0.5,-0.2,-0.5,-0.05,-0.5,0.1,-0.5,0.25,-0.5,0.4,-0.5,0.55,-0.5,0.7,-0.5,0.85,-0.5,1.0,-0.5,1.15,-0.5,1.3,-0.5,1.45,-0.35,-0.5,-0.35,-0.35,-0.35,-0.2,-0.35,-0.05,-0.35,0.1,-0.35,0.25,-0.35,0.4,-0.35,0.55,-0.35,0.7,-0.35,0.85,-0.35,1.0,-0.35,1.15,-0.35,1.3,-0.35,1.45,-0.2,-0.5,-0.2,-0.35,-0.2,-0.2,-0.2,-0.05,-0.2,0.1,-0.2,0.25,-0.2,0.4,-0.2,0.55,-0.2,0.7,-0.2,0.85,-0.2,1.0,-0.2,1.15,-0.2,1.3,-0.2,1.45,-0.05,-0.5,-0.05,-0.35,-0.05,-0.2,-0.05,-0.05,-0.05,0.1,-0.05,0.25,-0.05,0.4,-0.05,0.55,-0.05,0.7,-0.05,0.85,-0.05,1.0,-0.05,1.15,-0.05,1.3,-0.05,1.45,0.1,-0.5,0.1,-0.35,0.1,-0.2,0.1,-0.05,0.1,0.1,0.1,0.25,0.1,0.4,0.1,0.55,0.1,0.7,0.1,0.85,0.1,1.0,0.1,1.15,0.1,1.3,0.1,1.45,0.25,-0.5,0.25,-0.35,0.25,-0.2,0.25,-0.05,0.25,0.1,0.25,0.25,0.25,0.4,0.25,0.55,0.25,0.7,0.25,0.85,0.25,1.0,0.25,1.15,0.25,1.3,0.25,1.45,0.4,-0.5,0.4,-0.35,0.4,-0.2,0.4,-0.05,0.4,0.1,0.4,0.25,0.4,0.4,0.4,0.55,0.4,0.7,0.4,0.85,0.4,1.0,0.4,1.15,0.4,1.3,0.4,1.45,0.55,-0.5,0.55,-0.35,0.55,-0.2,0.55,-0.05,0.55,0.1,0.55,0.25,0.55,0.4,0.55,0.55,0.55,0.7,0.55,0.85,0.55,1.0,0.55,1.15,0.55,1.3,0.55,1.45,0.7,-0.5,0.7,-0.35,0.7,-0.2,0.7,-0.05,0.7,0.1,0.7,0.25,0.7,0.4,0.7,0.55,0.7,0.7,0.7,0.85,0.7,1.0,0.7,1.15,0.7,1.3,0.7,1.45,0.85,-0.5,0.85,-0.35,0.85,-0.2,0.85,-0.05,0.85,0.1,0.85,0.25,0.85,0.4,0.85,0.55,0.85,0.7,0.85,0.85,0.85,1.0,0.85,1.15,0.85,1.3,0.85,1.45,1.0,-0.5,1.0,-0.35,1.0,-0.2,1.0,-0.05,1.0,0.1,1.0,0.25,1.0,0.4,1.0,0.55,1.0,0.7,1.0,0.85,1.0,1.0,1.0,1.15,1.0,1.3,1.0,1.45,1.15,-0.5,1.15,-0.35,1.15,-0.2,1.15,-0.05,1.15,0.1,1.15,0.25,1.15,0.4,1.15,0.55,1.15,0.7,1.15,0.85,1.15,1.0,1.15,1.15,1.15,1.3,1.15,1.45,1.3,-0.5,1.3,-0.35,1.3,-0.2,1.3,-0.05,1.3,0.1,1.3,0.25,1.3,0.4,1.3,0.55,1.3,0.7,1.3,0.85,1.3,1.0,1.3,1.15,1.3,1.3,1.3,1.45,1.45,-0.5,1.45,-0.35,1.45,-0.2,1.45,-0.05,1.45,0.1,1.45,0.25,1.45,0.4,1.45,0.55,1.45,0.7,1.45,0.85,1.45,1.0,1.45,1.15,1.45,1.3,1.45,1.45],196,2)
        targetFieldValsExpected=DataArrayDouble([1.645976003316459, 1.454458180060204, 1.286087532859835, 1.147305389930914, 1.040143042030752, 0.9592075185603157, 0.8932542207607532, 0.8296417057622609, 0.7572539678257579, 0.6669048311361028, 0.551329882743212, 0.4064445075734602, 0.2323703965460786, 0.03253142054561309, 1.615321686989539, 1.414941300553572, 1.238383118538708, 1.096701655702075, 0.9955792747382535, 0.9271194507282707, 0.8741000712825546, 0.8201879508155141, 0.7537335933761495, 0.6656210809234322, 0.5470285414729397, 0.3927301586610237, 0.2044036897887453, -0.01181672742825013, 1.609602552867195, 1.400625195269133, 1.213287847440801, 1.065318574929208, 0.9717609562002842, 0.9182626517777217, 0.8760698972315855, 0.8258196104516153, 0.7586487405165288, 0.6686168424854784, 0.5434121624038266, 0.3741815029337978, 0.1661376046619205, -0.0704038088420833, 1.635421686625182, 1.422642113482769, 1.225977424080963, 1.066864693789366, 0.9864801043792362, 0.9486639217909161, 0.9075176697327381, 0.8471248730261529, 0.7660983406349626, 0.6675300501188994, 0.5320013361909732, 0.3404583135353376, 0.1074346390951333, -0.1520751802856468, 1.695346918429566, 1.489526279573347, 1.297678617961701, 1.139921240332637, 1.080508463804929, 1.036847769764088, 0.9687840669352359, 0.8790397822170175, 0.76938768351059, 0.6441978169925557, 0.4915328571013788, 0.2742929463574293, 0.0148214290833748, -0.2671755287427691, 1.782761788232491, 1.59423004798623, 1.422317125787222, 1.286999529473285, 1.20500638941831, 1.127058114031519, 1.022332539190471, 0.8945753999401338, 0.7469190939381181, 0.582396906110898, 0.4015920181411496, 0.1584700483835366, -0.1251860255418387, -0.4254052799545267, 1.881794862747652, 1.712890309994015, 1.557517508390291, 1.422727414977963, 1.308048056353061, 1.187569766723152, 1.03942150436647, 0.8677583087532357, 0.6766652050643343, 0.4703897480238999, 0.2497994532908829, -0.02005989176786582, -0.3224387891441491, -0.6331519303649853, 1.973114284621266, 1.820187301531605, 1.673403730111759, 1.528504440482262, 1.379693463484634, 1.207642134784147, 1.008217764780293, 0.7863328498822348, 0.5465383049529959, 0.2944879513187435, 0.03250657765404452, -0.2670900851421072, -0.5806516907976924, -0.8911331026431459, 2.038729888975378, 1.895652364645637, 1.751759791756183, 1.594035761810714, 1.403016809171641, 1.171403152610878, 0.913267035125007, 0.6343281031932027, 0.3434843176189371, 0.04195410032095204, -0.2645533663891493, -0.58577400250975, -0.8958218846257981, -1.192230697656513, 2.064018033720731, 1.922048791644444, 1.773847180028208, 1.600340336378483, 1.361620036333164, 1.060873411411508, 0.7373484802125152, 0.3868966266761109, 0.04316272760227413, -0.3009370030949727, -0.6505233805563486, -0.9669887470696283, -1.250005719852354, -1.519122595631787, 2.039938287785342, 1.887400820799651, 1.722008733683987, 1.523879290022419, 1.23834392230135, 0.8606985727866472, 0.4844892131548788, 0.08077959236877175, -0.3195742594962179, -0.726291368696764, -1.094357645641832, -1.359078900303776, -1.604725656501341, -1.845297168323687, 1.965762248218393, 1.791665198563286, 1.595056719739704, 1.353692777435502, 1.033006623003495, 0.6416349531117889, 0.2290046916364761, -0.1993180965088852, -0.6311618804827295, -1.051489875129883, -1.409404344854132, -1.681249363331096, -1.917859637689007, -2.145034400762945, 1.849053542205925, 1.648479366622312, 1.418493963148431, 1.141939527533839, 0.8042385795619003, 0.4127534639189761, -0.008572116677791453, -0.4428317297963555, -0.8745477268718713, -1.281769237471681, -1.635421857742795, -1.926210204560556, -2.175577364628722, -2.405762639746138, 1.701519686999922, 1.475879908746998, 1.219065416294153, 0.9203732349759972, 0.5740137315474942, 0.1856460506119944, -0.2298288912529738, -0.6558565521653752, -1.075391078040103, -1.469402631469075, -1.820558929095151, -2.123592211415966, -2.388177455227765, -2.628832075944413])
        coeffsExpected=DataArrayDouble([0.3953237723894342,-0.17220705170185724,0.620727139132215,-0.01938292763088709,-0.007524685306185282,0.0016277944443884584,-0.0005209587893117361,-1.8992696595839718,-0.13154330748345855,0.11248800965389728,-0.47310750305033406,0.03685741122098605,0.21362468750754374,0.8082608687799991,-0.6775548200221704,-0.027683208482275873,-0.007806877014495724,-0.013539239795959668,0.3478535665778018,0.005145793726360813,0.03708618549628136,-0.18235332489209385,-0.04517273339177797,-0.081755114492025,0.12791746560435255,0.09659355695676189,-0.024809653129318366,0.08327587452569823,-1.790380673650165,-0.10622983512164165,0.14989029282340274,0.05949513762355707,0.004548072841131278,0.011252095917834793,-0.004848057194721367,-0.2658537133108412,0.016651579133606154,-0.021640915366981317,0.008975511042160175,-0.021052213988815974,-0.09347841701844657,0.03533229488135717,-0.014556185287109863,-0.27228591670520086,0.002989987191209683,-0.5489428537951813,-0.02134456783001304,-0.22462281620064825,0.005230853443767429,-0.1894678262257301,0.0033140729457334884,5.295483062326795,-0.2724500716060311,0.026433905662192683,0.01368706308878908,-0.03014264855048227,0.053679001877659956,0.08109477254132096,-0.005004603067203444,0.016907143132293558,0.2105509502082437,0.003657404455024417,-4.904755847017426,0.01634808163992959,-0.008325515865305198,0.062188432751569676,-0.013114633511406406,0.11020519384963083,-0.008599402366091309,-0.012125149710784723,0.31723729052927313,-0.10298398036815914,-0.07250078775612204,0.39976713701763433,0.45897498107347223,0.01018626210400031,0.20163425809089347,0.19729093298588943,0.42863333455911523,0.015595097081693168,0.06060353651437489,-0.16379444813161725,-0.43290344196574165,-0.5931022701412187,1.1906610004748832,0.44418106894148945,0.06536220001548931,0.010261694323554562,-0.05943099382075491,-0.04939614579484797,0.002234505477641322,-0.011262130967449935,0.09644905007708474,-0.029518792883267808,0.41564004027396634,-0.18459770295961597,0.3100981306103734,-0.2509873737065425,0.5434321443668653,0.3009912967350914,1.9560655796099518,-0.7143435150084513,-1.5123449469879784])
        #
        nbOfInputPoints=100;
        f=MEDCouplingFieldDouble.New(ON_NODES_KR,ONE_TIME);
        mesh=MEDCoupling1SGTUMesh.New("aMesh",NORM_POINT1);
        mesh.setCoords(srcPointCoordsXY);
        f.setMesh(mesh);
        f.setArray(srcFieldValsOnPoints);
        f.checkConsistencyLight();
        #
        res0=f.getValueOn([-0.5,-0.5]);
        self.assertAlmostEqual(targetFieldValsExpected.getIJ(0,0),res0[0],10)
        #
        valuesToTest=f.getValueOnMulti(targetPointCoordsXY);
        self.assertEqual(196,valuesToTest.getNumberOfTuples());
        self.assertEqual(1,valuesToTest.getNumberOfComponents());
        for i in range(40):
            self.assertAlmostEqual(targetFieldValsExpected[i],valuesToTest.getIJ(i,0),10)
            pass
        fd=f.getDiscretization()
        del f
        self.assertTrue(isinstance(fd,MEDCouplingFieldDiscretizationKriging))
        coeffs,isDrift=fd.computeVectorOfCoefficients(mesh,srcFieldValsOnPoints)
        self.assertEqual(3,isDrift)
        self.assertTrue(coeffsExpected.isEqual(coeffs,1e-8))
        # testing matrix
        pts3=[-0.5,-0.5,-0.5,-0.35,-0.35,-0.2]
        mesh.setCoords(srcPointCoordsXY[:4])
        m,nbCols=fd.computeEvaluationMatrixOnGivenPts(mesh,pts3)
        self.assertTrue(m.isEqual(DataArrayDouble([0.05768877688524917,-4.438982030395039,1.9495386255911573,3.431754627918642,0.11803848510231275,-4.138339658420563,1.6630742187104417,3.357226954607818,0.14630203028580618,-3.5156045565871734,1.414680070737206,2.954622455564169]),1e-12))
        if MEDCouplingHasNumPyBindings():
            import numpy as np
            m0=m.toNumPyArray() ; m0=m0.reshape(3,nbCols) ; m0=np.matrix(m0)
            srcFieldValsOnPoints2=DataArrayDouble(4,2) ; srcFieldValsOnPoints2[:,0]=srcFieldValsOnPoints[:4] ; srcFieldValsOnPoints2[:,1]=2*srcFieldValsOnPoints[:4]
            n0=srcFieldValsOnPoints2.toNumPyArray() ; n0=n0.reshape(4,2) ; n0=np.matrix(n0)
            #
            f=MEDCouplingFieldDouble.New(ON_NODES_KR,ONE_TIME) ;  f.setMesh(mesh) ; f.setArray(srcFieldValsOnPoints2) ; f.checkConsistencyLight()
            self.assertTrue(DataArrayDouble(np.array((m0*n0))).isEqual(f.getValueOnMulti(pts3),1e-14))
            pass
        #
        pass

    # test the when input slice is all the same object is return by MEDCouplingMesh.buildPartRange
    def testSwig2MeshPartSlice1(self):
        a=DataArrayDouble(4) ; a.iota()
        c=MEDCouplingCMesh() ; c.setCoords(a,a) ; m=c.buildUnstructured()
        fc0=c.getMeasureField(False) ; fc1=fc0[:] ; fc2=fc0*fc1 ; fc2.setName(fc0.getName())
        self.assertEqual(fc0.getMesh().getHiddenCppPointer(),fc1.getMesh().getHiddenCppPointer())
        self.assertEqual(fc2.getMesh().getHiddenCppPointer(),fc1.getMesh().getHiddenCppPointer())
        self.assertTrue(fc2.isEqual(fc1,1e-12,1e-12))
        #
        fm0=m.getMeasureField(False) ; fm1=fm0[:] ; fm2=fm0*fm1 ; fm2.setName(fm0.getName())
        self.assertEqual(fm0.getMesh().getHiddenCppPointer(),fm1.getMesh().getHiddenCppPointer())
        self.assertEqual(fm2.getMesh().getHiddenCppPointer(),fm1.getMesh().getHiddenCppPointer())
        self.assertTrue(fm2.isEqual(fm1,1e-12,1e-12))
        pass

    # test the correct behaviour when attempting to aggregate two fields whose mesh is null
    def testSwig2MergeFieldsOnFieldsHavingNoMesh(self):
        a=DataArrayDouble(4) ; a.iota() ; a*=1.5
        c=MEDCouplingCMesh() ; c.setCoords(a,a) ; f1=c.getMeasureField(False)
        f1.setMesh(None) ; f2=f1.deepCopy() ; f2*=2
        f3=MEDCouplingFieldDouble.MergeFields(f1,f2)
        daExp=DataArrayDouble([2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5])
        self.assertTrue(f3.getArray().isEqual(daExp,1e-12))
        self.assertEqual(f3.getTypeOfField(),ON_CELLS)
        self.assertEqual(f3.getMesh(),None)
        f4=MEDCouplingFieldDouble.MergeFields([f1,f2])
        self.assertTrue(f4.getArray().isEqual(daExp,1e-12))
        self.assertEqual(f4.getTypeOfField(),ON_CELLS)
        self.assertEqual(f4.getMesh(),None)
        pass

    # test a simple node to cell convertion of a field
    def testSwig2NodeToCellDiscretization1(self):
        f=MEDCouplingFieldDouble(ON_NODES) ; f.setTime(1.1,2,3)
        a1=DataArrayDouble(4) ; a1.iota()
        a2=DataArrayDouble(3) ; a2.iota()
        m=MEDCouplingCMesh() ; m.setCoords(a1,a2)
        f.setMesh(m)
        arr=DataArrayDouble([21.,121.,20.,120.,19.,119.,18.,118.,17.,117.,16.,116.,15.,115.,14.,114.,13.,113.,12.,112.,11.,111.,10.,110.],12,2) ; arr.setInfoOnComponents(["aa [km]","bbb [kJ]"])
        f.setArray(arr) ; f.setName("toto")
        #
        f2=f.nodeToCellDiscretization()
        self.assertEqual(ON_CELLS,f2.getTypeOfField())
        self.assertEqual("toto",f2.getName())
        self.assertEqual([1.1,2,3],f2.getTime())
        self.assertEqual(["aa [km]","bbb [kJ]"],f2.getArray().getInfoOnComponents())
        self.assertEqual(6,f2.getArray().getNumberOfTuples())
        self.assertEqual(f.getMesh().getHiddenCppPointer(),f2.getMesh().getHiddenCppPointer())
        exp=DataArrayDouble([18.5,118.5,17.5,117.5,16.5,116.5,14.5,114.5,13.5,113.5,12.5,112.5],6,2) ; exp.setInfoOnComponents(["aa [km]","bbb [kJ]"])
        self.assertTrue(f2.getArray().isEqual(exp,1e-13))
        pass

    def testSwig2MeshOrientCorrectly2DCells1(self):
        m=MEDCouplingUMesh("mesh",2)
        coo=DataArrayDouble([1.,0.,0.5,-0.1,0.,1.,0.,0.,0.07,0.5,0.59,0.5],6,2)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_TRI6,[3,0,2,1,5,4])
        m.insertNextCell(NORM_QPOLYG,[3,0,2,1,5,4])
        self.assertTrue(DataArrayDouble([-0.58093333350930543,-0.58093333350930543]).isEqual(m.getMeasureField(False).getArray(),1e-12))
        m.changeSpaceDimension(3)
        m.orientCorrectly2DCells([0.,0.,-1.],False)
        #
        m.checkConsistencyLight()
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([6,3,2,0,4,5,1, 32,3,2,0,4,5,1])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,7,14])))
        m.changeSpaceDimension(2)
        self.assertTrue(DataArrayDouble([0.58093333350930543,0.58093333350930543]).isEqual(m.getMeasureField(False).getArray(),1e-12))
        pass

    def testSwig2Hexa8HavingFacesWarped1(self):
        """ This test is bases on a "error" of interpolation detected. After investigation cell #3 of src is warped that leads to the fact that when trg is
        intersected with src the sum of intersection volume is greater than the volume of the trg cell.
        A test that can be done is to split the cell #3 of src into tetrohedrons and by summing all the volumes it does not fit the volume computed of cell#3 unsplitted (expect for
        GENERAL_24).
        """
        srcCoo=DataArrayDouble([0.15694071546650565,0.09383333333333337,6.920842121738133,0.15774332475430292,0.185486666666667,6.920682472824616,0.1585459340420992,0.27713999999999994,6.9205228239111,0.07427195882345167,0.05782666666666668,6.937285959830335,0.06343673343819695,0.11347333333333297,6.939441220162809,0.05260150805294228,0.16911999999999996,6.941596480495282,0.014076262238703396,0.04800666666666667,6.949259628344076,0.014076262238703396,0.07092000000000007,6.949259628344076,0.15407499632681992,0.09383333333333338,6.897607484780063,0.15489234394181514,0.18548666666666702,6.897567331066572,0.15570969155680933,0.27714,6.897527177353081,0.06988819198237989,0.05782666666666669,6.901743317269663,0.05885399917995321,0.11347333333333298,6.9022853924017955,0.047819806377526586,0.16912,6.902827467533927,0.0085871208577874,0.048006666666666684,6.9047548457815076,0.0085871208577874,0.07092000000000008,6.9047548457815076,0.153883333333333,0.09383333333333338,6.820902,0.154701666666667,0.18548666666666702,6.820902,0.15551999999999996,0.27714,6.820902,0.06959499999999999,0.05782666666666669,6.820902,0.058547499999999975,0.11347333333333298,6.820902,0.04749999999999999,0.16912,6.820902],22,3)
        src=MEDCouplingUMesh("TBmesh3D",3) ; src.setCoords(srcCoo)
        src.allocateCells()
        src.insertNextCell(NORM_HEXA8,[0,1,4,3,8,9,12,11])
        src.insertNextCell(NORM_HEXA8,[1,2,5,4,9,10,13,12])
        src.insertNextCell(NORM_HEXA8,[4,5,7,6,12,13,15,14])
        src.insertNextCell(NORM_HEXA8,[8,9,12,11,16,17,20,19])
        src.insertNextCell(NORM_HEXA8,[9,10,13,12,17,18,21,20])
        src.checkConsistency()
        # trg is useless here but I keep it in case of MEDCouplingRemapper were expected to do something about warped NORM_HEXA8
        trgCoo=DataArrayDouble([0.0960891897852753,0.105088620541845,6.8598,0.0599574480546212,0.118434267436059,6.8598,0.113514510609589,0.14874473653263,6.8598,0.0831322609794463,0.167319109733883,6.8598,0.0960891897852753,0.105088620541845,6.92146666666667,0.0599574480546212,0.118434267436059,6.92146666666667,0.113514510609589,0.14874473653263,6.92146666666667,0.0831322609794463,0.167319109733883,6.92146666666667],8,3)
        trg=MEDCouplingUMesh("MESH",3) ; trg.setCoords(trgCoo)
        trg.allocateCells()
        trg.insertNextCell(NORM_HEXA8,[0,1,3,2,4,5,7,6])
        #
        srcFace=src.buildDescendingConnectivity()[0]
        conn=MEDCoupling1SGTUMesh(srcFace).getNodalConnectivity() ; conn.rearrange(4)
        eqFaces=srcFace.computePlaneEquationOf3DFaces()
        nodeIdInCell=3
        e=(srcFace.getCoords()[conn[:,nodeIdInCell]]*eqFaces[:,:-1]).sumPerTuple()+eqFaces[:,3]# e represent the error between the expected 'a*X+b*Y+c*Z+d' in eqFaces and 0. Closer e to 0. is closer the 4th point is to the plane built with the 3 first points
        lambd=-e/(eqFaces[:,:3]**2).sumPerTuple()
        pts=lambd*eqFaces[:,:-1]+srcFace.getCoords()[conn[:,nodeIdInCell]]#pts represent the projection of the last points of each NORM_QUAD4 to the plane defined by the 3 first points of the NORM_QUAD4 cell
        shouldBeZero=(pts*eqFaces[:,:-1]).sumPerTuple()+eqFaces[:,3]# this line is useless only to be sure that pts are on the plane.
        check=(pts-srcFace.getCoords()[conn[:,nodeIdInCell]]).magnitude() # check contains the distance of the last point to its plane
        idsToTest=check.findIdsNotInRange(0.,1e-10)
        self.assertTrue(idsToTest.isEqual(DataArrayInt([17,18,19,20,22,23,24])))
        idsToTest2=idsToTest.findIdsNotInRange(18,22)
        self.assertTrue(idsToTest2.isEqual(DataArrayInt([0,4,5,6])))
        idsToTest2.rearrange(2)
        self.assertTrue(idsToTest2.sumPerTuple().isEqual(DataArrayInt([4,11])))
        pass

    def testSwig2SortHexa8EachOther1(self):
        """
        testing MEDCoupling1SGTUMesh.sortHexa8EachOther method
        """
        coords1=DataArrayDouble([(-0.5,0.5,-0.5),(0.5,-0.5,-0.5),(-0.5,-0.5,0.5),(-0.5,-0.5,-0.5),(0.5,-0.5,0.5),(-0.5,0.5,0.5),(0.5,0.5,0.5),(0.5,0.5,-0.5)])
        m1=MEDCouplingUMesh("m1",3) ; m1.setCoords(coords1)
        m1.allocateCells() ; m1.insertNextCell(NORM_HEXA8,[7,1,3,0,6,4,2,5])
        m1.checkConsistencyLight()
        #
        m2=m1.deepCopy() ; m2.setName("m2")
        #
        trs=[[0.,0.,-1.],[0.,0.,1.],[1.,0.,0.],[0.,-1.,0.],[-1.,0.,0.],[0.,1.,0.]]
        for i,t in enumerate(trs):
            for j in range(64):
                j2=(j//16) ; j1=((j%16)//4) ; j0=(j%4)
                m11=m1.deepCopy()
                m11.rotate([0.,0.,0.],[0.,0.,1.],float(j0)*pi/2)
                m11.rotate([0.,0.,0.],[0.,1.,0.],float(j1)*pi/2)
                m11.rotate([0.,0.,0.],[1.,0.,0.],float(j2)*pi/2)
                m11.translate(t)
                #
                m=MEDCouplingUMesh.MergeUMeshes(m2,m11)
                m.mergeNodes(1e-12)
                self.assertEqual(12,m.getNumberOfNodes())
                m=MEDCoupling1SGTUMesh(m)
                m.sortHexa8EachOther()
                tmp0=m.buildUnstructured().tetrahedrize(PLANAR_FACE_6)[0].buildUnstructured()
                self.assertEqual(20,tmp0.computeSkin().getNumberOfCells())
                pass
            pass
        pass

    def testSwig2normMinComputeAbs1(self):
        d=DataArrayDouble([4,-5,2,6.1,-7.33,1,-1,3e2,0.07,-0.009,-6,-1e30],4,3)
        d.setInfoOnComponents(["XX [m]","YYY [km]","ABSJJ [MW]"])
        d0=d.computeAbs()
        dExp=d.deepCopy() ; dExp.abs()
        self.assertTrue(dExp.isEqual(d0,1e-12))
        e=d0-DataArrayDouble([4,5,2,6.1,7.33,1,1,3e2,0.07,0.009,6,1e30],4,3)
        self.assertAlmostEqual(0.,e.normMin(),13)
        self.assertAlmostEqual(0.009,d.normMin(),13)
        #
        di=DataArrayInt([3,-12,5,6,14,16,-23,100,23,-1,0,-6],4,3)
        di.setInfoOnComponents(["XX [m]","YYY [km]","ABSJJ [MW]"])
        d0i=di.computeAbs()
        diExp=di.deepCopy() ; diExp.abs()
        self.assertTrue(diExp.isEqual(d0i))
        self.assertEqual([3,12,5,6,14,16,23,100,23,1,0,6],d0i.getValues())
        pass

    def testSwig2GetCellsContainingPointsForNonConvexPolygon1(self):
        coo=DataArrayDouble([-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5,-0.5,0.,-0.5,0.,0.,0.5,0.,],7,2)
        m=MEDCouplingUMesh("Intersect2D",2) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_POLYGON,[6,3,4,5])
        m.insertNextCell(NORM_POLYGON,[4,0,1,2,6,5])
        m.checkConsistency()
        #
        self.assertTrue(m.getCellsContainingPoint((0.4,-0.4),1e-12).isEqual(DataArrayInt([0])))
        self.assertTrue(m.getCellsContainingPoint((-0.4,-0.4),1e-12).isEqual(DataArrayInt([1])))
        self.assertTrue(m.getCellsContainingPoint((0.,-0.4),1e-12).isEqual(DataArrayInt([0,1])))
        pass

    def testSwig2GetCellsContainingPointsForNonConvexPolygon2(self):
        coo=DataArrayDouble([-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5,-0.5,-2.0816681711721685e-17,-2.0816681711721685e-17,-0.17677669529663687,0.1767766952966369,0.,0.5,0.5,0.,0.17677669529663684,-0.17677669529663692,0.17677669529663692,0.17677669529663684,-0.17677669529663692,-0.17677669529663687,0.,-0.5,-0.5,0.,0.33838834764831843,-0.3383883476483185,-0.33838834764831843,0.33838834764831843,-0.21213203435596423,0.21213203435596426,0.2121320343559642,-0.2121320343559643,0.21213203435596426,0.2121320343559642,-0.21213203435596423,-0.21213203435596428,0.3560660171779821,-0.35606601717798214,-0.35606601717798214,0.35606601717798214,0.19445436482630052,-0.19445436482630063,-0.19445436482630055,0.19445436482630057,0.,0.27],24,2)
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_QPOLYG,[8,5,4,9])
        m.insertNextCell(NORM_QPOLYG,[5,8,4,10])
        m.insertNextCell(NORM_QPOLYG,[16,8,5,15,21,9,22,17])
        m.insertNextCell(NORM_QPOLYG,[15,1,2,3,16,20,6,7,19,17])
        m.insertNextCell(NORM_QPOLYG,[15,5,8,16,22,10,21,18])
        m.insertNextCell(NORM_QPOLYG,[16,3,0,1,15,19,11,12,20,18])
        m.checkConsistency()
        self.assertTrue(m.getCellsContainingPoint([0.,0.27],1e-12).isEqual(DataArrayInt([2])))
        pass

    def testSwig2DAIGetIdsEqualTuple1(self):
        da=DataArrayInt([0,7,1,2,4,1,2,1,1,2,0,1,2,1,5,1,1,2],9,2)
        self.assertTrue(da.findIdsEqualTuple([1,2]).isEqual(DataArrayInt([1,4,8])))
        self.assertTrue(da.findIdsEqualTuple((1,2)).isEqual(DataArrayInt([1,4,8])))
        self.assertTrue(da.findIdsEqualTuple(DataArrayInt([1,2])).isEqual(DataArrayInt([1,4,8])))
        da.rearrange(3)
        self.assertRaises(InterpKernelException,da.findIdsEqualTuple,[1,2])# mismatch nb of compo (3) and nb of elts in input tuple (2)
        self.assertTrue(da.findIdsEqualTuple([2,0,1]).isEqual(DataArrayInt([3])))
        self.assertTrue(da.findIdsEqualTuple([2,0,7]).isEqual(DataArrayInt([])))
        da.rearrange(1)
        self.assertTrue(da.findIdsEqualTuple(2).isEqual(DataArrayInt([3,6,9,12,17])))
        self.assertTrue(da.findIdsEqualTuple(2).isEqual(da.findIdsEqual(2)))
        pass

    def testSwig2GaussNEStaticInfo1(self):
        self.assertTrue(DataArrayDouble(MEDCouplingFieldDiscretizationGaussNE.GetWeightArrayFromGeometricType(NORM_TRI3)).isEqual(DataArrayDouble([0.16666666666666666,0.16666666666666666,0.16666666666666666]),1e-12))
        self.assertTrue(DataArrayDouble(MEDCouplingFieldDiscretizationGaussNE.GetRefCoordsFromGeometricType(NORM_TRI3)).isEqual(DataArrayDouble([0.,0.,1.,0.,0.,1.]),1e-12))
        self.assertTrue(DataArrayDouble(MEDCouplingFieldDiscretizationGaussNE.GetLocsFromGeometricType(NORM_TRI3)).isEqual(DataArrayDouble([0.16666666666666666,0.16666666666666666,0.6666666666666667,0.16666666666666666,0.16666666666666666,0.6666666666666667]),1e-12))
        pass

    def testSwigReverseNodalConnOnStructuredMesh(self):
        # 1D - standard
        c=MEDCouplingCMesh() ; arr=DataArrayDouble(10) ; arr.iota()
        c.setCoordsAt(0,arr)
        rn,rni=c.getReverseNodalConnectivity()
        rn2,rni2=c.buildUnstructured().getReverseNodalConnectivity()
        self.assertTrue(rn.isEqual(DataArrayInt([0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8])))
        self.assertTrue(rni.isEqual(DataArrayInt([0,1,3,5,7,9,11,13,15,17,18])))
        self.assertTrue(rn.isEqual(rn2)) ; self.assertTrue(rni.isEqual(rni2))
        # 1D - limit
        c=MEDCouplingCMesh() ; arr=DataArrayDouble(1) ; arr.iota()
        c.setCoordsAt(0,arr)
        rn,rni=c.getReverseNodalConnectivity()
        rn2,rni2=c.buildUnstructured().getReverseNodalConnectivity()
        self.assertTrue(rn.isEqual(DataArrayInt([0])))
        self.assertTrue(rni.isEqual(DataArrayInt([0,1])))
        self.assertTrue(rn.isEqual(rn2)) ; self.assertTrue(rni.isEqual(rni2))
        # 1D - limit
        c=MEDCouplingCMesh() ; arr=DataArrayDouble(0) ; arr.iota()
        c.setCoordsAt(0,arr)
        rn,rni=c.getReverseNodalConnectivity()
        rn.isEqual(DataArrayInt([]))
        rni.isEqual(DataArrayInt([0]))
        # 2D - standard
        c=MEDCouplingCMesh() ; arr=DataArrayDouble(5) ; arr.iota() ; arr2=DataArrayDouble(4) ; arr.iota()
        c.setCoords(arr,arr2)
        rn,rni=c.getReverseNodalConnectivity()
        rn2,rni2=c.buildUnstructured().getReverseNodalConnectivity()
        self.assertTrue(rn.isEqual(DataArrayInt([0,0,1,1,2,2,3,3,0,4,0,1,4,5,1,2,5,6,2,3,6,7,3,7,4,8,4,5,8,9,5,6,9,10,6,7,10,11,7,11,8,8,9,9,10,10,11,11])))
        self.assertTrue(rni.isEqual(DataArrayInt([0,1,3,5,7,8,10,14,18,22,24,26,30,34,38,40,41,43,45,47,48])))
        self.assertTrue(rn.isEqual(rn2)) ; self.assertTrue(rni.isEqual(rni2))
        # 2D - limit
        c=MEDCouplingCMesh() ; arr=DataArrayDouble(10) ; arr.iota() ; arr2=DataArrayDouble(1) ; arr.iota()
        c.setCoords(arr,arr2)
        rn,rni=c.getReverseNodalConnectivity()
        self.assertTrue(rn.isEqual(DataArrayInt([0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8])))
        self.assertTrue(rni.isEqual(DataArrayInt([0,1,3,5,7,9,11,13,15,17,18])))
        # 2D - limit
        c=MEDCouplingCMesh() ; arr=DataArrayDouble(10) ; arr.iota() ; arr2=DataArrayDouble(1) ; arr.iota()
        c.setCoords(arr2,arr)
        rn,rni=c.getReverseNodalConnectivity()
        self.assertTrue(rn.isEqual(DataArrayInt([0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8])))
        self.assertTrue(rni.isEqual(DataArrayInt([0,1,3,5,7,9,11,13,15,17,18])))
        # 3D - standard
        c=MEDCouplingCMesh() ; arr0=DataArrayDouble(5) ; arr0.iota() ; arr1=DataArrayDouble(3) ; arr1.iota() ; arr2=DataArrayDouble(4) ; arr2.iota()
        c.setCoords(arr0,arr1,arr2)
        rn,rni=c.getReverseNodalConnectivity()
        self.assertTrue(rn.isEqual(DataArrayInt([0,0,1,1,2,2,3,3,0,4,0,1,4,5,1,2,5,6,2,3,6,7,3,7,4,4,5,5,6,6,7,7,0,8,0,1,8,9,1,2,9,10,2,3,10,11,3,11,0,4,8,12,0,1,4,5,8,9,12,13,1,2,5,6,9,10,13,14,2,3,6,7,10,11,14,15,3,7,11,15,4,12,4,5,12,13,5,6,13,14,6,7,14,15,7,15,8,16,8,9,16,17,9,10,17,18,10,11,18,19,11,19,8,12,16,20,8,9,12,13,16,17,20,21,9,10,13,14,17,18,21,22,10,11,14,15,18,19,22,23,11,15,19,23,12,20,12,13,20,21,13,14,21,22,14,15,22,23,15,23,16,16,17,17,18,18,19,19,16,20,16,17,20,21,17,18,21,22,18,19,22,23,19,23,20,20,21,21,22,22,23,23])))
        self.assertTrue(rni.isEqual(DataArrayInt([0,1,3,5,7,8,10,14,18,22,24,25,27,29,31,32,34,38,42,46,48,52,60,68,76,80,82,86,90,94,96,98,102,106,110,112,116,124,132,140,144,146,150,154,158,160,161,163,165,167,168,170,174,178,182,184,185,187,189,191,192])))
        rn2,rni2=c.buildUnstructured().getReverseNodalConnectivity()
        self.assertTrue(rn.isEqual(rn2)) ; self.assertTrue(rni.isEqual(rni2))
        pass

    def testSwig2CellToNodeDiscretization1(self):
        m=MEDCouplingCMesh() ; arr0=DataArrayDouble(5) ; arr0.iota() ; arr1=DataArrayDouble(4) ; arr1.iota() ; m.setCoords(arr0,arr1)
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(m) ; f.setTime(1.1,5,6)
        arr=DataArrayDouble(12) ; arr.iota()
        arr=DataArrayDouble.Meld(arr,arr+100.) ; arr.setInfoOnComponents(["aaa","bbb"])
        f.setArray(arr)
        f.checkConsistencyLight()
        #
        ref=DataArrayDouble([0.,0.5,1.5,2.5,3.,2.,2.5,3.5,4.5,5.,6.,6.5,7.5,8.5,9.,8.,8.5,9.5,10.5,11.])
        ref=DataArrayDouble.Meld(ref,ref+100.) ; ref.setInfoOnComponents(["aaa","bbb"])
        f2=f.cellToNodeDiscretization()
        f2.checkConsistencyLight()
        self.assertEqual(f2.getTime()[1:],[5,6])
        self.assertAlmostEqual(f2.getTime()[0],1.1,15)
        self.assertEqual(f2.getMesh().getHiddenCppPointer(),m.getHiddenCppPointer())
        self.assertTrue(f2.getArray().isEqual(ref,1e-12))
        rn,rni=m.getReverseNodalConnectivity()
        rni2=(rni.deltaShiftIndex()).convertToDblArr()
        arr2=(f.getArray()[rn]).accumulatePerChunck(rni)/rni2
        self.assertTrue(f2.getArray().isEqual(arr2,1e-12))
        del f2
        #
        u=m.buildUnstructured() ; f.setMesh(u) ; del m
        f3=f.cellToNodeDiscretization()
        f3.checkConsistencyLight()
        self.assertEqual(f3.getTime()[1:],[5,6])
        self.assertAlmostEqual(f3.getTime()[0],1.1,15)
        self.assertEqual(f3.getMesh().getHiddenCppPointer(),u.getHiddenCppPointer())
        self.assertTrue(f3.getArray().isEqual(ref,1e-12))
        pass

    def testSwig2GetMeshSpaceDimensionCMesh1(self):
        c=MEDCouplingCMesh()
        arr0=DataArrayDouble([0,1,2])
        arr1=DataArrayDouble([0])
        c.setCoords(arr0,arr0,arr0)
        self.assertEqual(c.getMeshDimension(),3)
        self.assertEqual(c.getSpaceDimension(),3)
        #
        c.setCoords(arr0,arr0,arr1)
        self.assertEqual(c.getMeshDimension(),2)
        self.assertEqual(c.getSpaceDimension(),3)
        #
        c.setCoords(arr0,arr0)
        self.assertEqual(c.getMeshDimension(),2)
        self.assertEqual(c.getSpaceDimension(),2)
        #
        c.setCoords(arr0,arr1)
        self.assertEqual(c.getMeshDimension(),1)
        self.assertEqual(c.getSpaceDimension(),2)
        #
        c.setCoords(arr0)
        self.assertEqual(c.getMeshDimension(),1)
        self.assertEqual(c.getSpaceDimension(),1)
        #
        c.setCoords(arr1)
        self.assertEqual(c.getMeshDimension(),0)
        self.assertEqual(c.getSpaceDimension(),1)
        pass

    def testSwig2BuildSpreadZonesWithPolyOnQPolyg1(self):
        nx=6
        ny=6
        m=MEDCouplingCMesh()
        arr1=DataArrayDouble(nx) ; arr1.iota()
        arr2=DataArrayDouble(ny) ; arr2.iota()
        m.setCoords(arr1,arr2)
        m=m.buildUnstructured()
        da=DataArrayInt.Range(nx-1,(nx-1)*(ny-1),nx)
        m2=m[da] ; m2.simplexize(0)
        dan=da.buildComplement(m.getNumberOfCells())
        m1=m[dan]
        m=MEDCouplingUMesh.MergeUMeshesOnSameCoords(m1,m2)
        #
        m.convertLinearCellsToQuadratic()
        m1=m[::2] ; m2=m[1::2] ; m2.convertAllToPoly()
        m=MEDCouplingUMesh.MergeUMeshesOnSameCoords(m1,m2)
        p=m.buildSpreadZonesWithPoly()
        self.assertTrue(p.getNodalConnectivity().isEqual(DataArrayInt([32,1,0,6,12,18,24,30,31,32,33,34,35,29,23,17,11,5,4,3,2,36,37,94,62,72,83,84,86,89,99,92,93,82,71,60,51,49,46,43,40])))
        self.assertTrue(p.getNodalConnectivityIndex().isEqual(DataArrayInt([0,41])))
        self.assertTrue(p.getCoords().isEqual(DataArrayDouble([0.,0.,1.,0.,2.,0.,3.,0.,4.,0.,5.,0.,0.,1.,1.,1.,2.,1.,3.,1.,4.,1.,5.,1.,0.,2.,1.,2.,2.,2.,3.,2.,4.,2.,5.,2.,0.,3.,1.,3.,2.,3.,3.,3.,4.,3.,5.,3.,0.,4.,1.,4.,2.,4.,3.,4.,4.,4.,5.,4.,0.,5.,1.,5.,2.,5.,3.,5.,4.,5.,5.,5.,0.5,0.,0.,0.5,0.5,1.,1.,0.5,1.5,0.,1.5,1.,2.,0.5,2.5,0.,2.5,1.,3.,0.5,3.5,0.,3.5,1.,4.,0.5,4.5,0.,4.5,1.,5.,0.5,1.,1.5,1.5,2.,2.,1.5,2.5,2.,3.,1.5,3.5,2.,4.,1.5,4.5,2.,5.,1.5,0.5,2.,0.,2.5,0.5,3.,1.,2.5,2.,2.5,2.5,3.,3.,2.5,3.5,3.,4.,2.5,4.5,3.,5.,2.5,0.,3.5,0.5,4.,1.,3.5,1.5,3.,1.5,4.,2.,3.5,3.,3.5,3.5,4.,4.,3.5,4.5,4.,5.,3.5,0.,4.5,0.5,5.,1.,4.5,1.5,5.,2.,4.5,2.5,4.,2.5,5.,3.,4.5,4.,4.5,4.5,5.,5.,4.5,0.,1.5,0.5,1.5,1.5,2.5,2.5,3.5,3.5,4.5,3.5,5.0],100,2),1e-13))
        pass

    def testSwigExtendedSlice1(self):
        d=DataArrayInt([5,6,7])
        self.assertTrue(d[2:].isEqual(DataArrayInt([7])))
        self.assertTrue(d[3:].isEqual(DataArrayInt([])))
        try:
            d[4:]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        d=DataArrayInt([5,6,7,8])
        self.assertEqual(d[-1],8)
        self.assertEqual(d[-4],5)
        try:
            d[-5]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        self.assertTrue(d[2::-1].isEqual(DataArrayInt([7,6,5])))
        self.assertTrue(d[0::-1].isEqual(DataArrayInt([5])))
        self.assertTrue(d[-1::-1].isEqual(DataArrayInt([8,7,6,5])))
        self.assertTrue(d[-3::-1].isEqual(DataArrayInt([6,5])))
        self.assertTrue(d[-5::-1].isEqual(DataArrayInt([])))
        try:
            d[-6::-1]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        d=DataArrayInt([])
        self.assertTrue(d[0:].isEqual(DataArrayInt([])))
        #
        d=DataArrayDouble([5,6,7])
        self.assertTrue(d[2:].isEqual(DataArrayDouble([7]),1e-12))
        self.assertTrue(d[3:].isEqual(DataArrayDouble([]),1e-12))
        try:
            d[4:]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        d=DataArrayDouble([5,6,7,8])
        self.assertAlmostEqual(d[-1],8.,12)
        self.assertAlmostEqual(d[-4],5.,12)
        try:
            d[-5]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        self.assertTrue(d[2::-1].isEqual(DataArrayDouble([7,6,5]),1e-12))
        self.assertTrue(d[0::-1].isEqual(DataArrayDouble([5]),1e-12))
        self.assertTrue(d[-1::-1].isEqual(DataArrayDouble([8,7,6,5]),1e-12))
        self.assertTrue(d[-3::-1].isEqual(DataArrayDouble([6,5]),1e-12))
        self.assertTrue(d[-5::-1].isEqual(DataArrayDouble([]),1e-12))
        try:
            d[-6::-1]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        d=DataArrayDouble([])
        self.assertTrue(d[0:].isEqual(DataArrayDouble([]),1e-12))
        pass

    def testSwig2Hexa27GP1(self):
        """ This test focused on shape functions of hexa27.
        """
        coo=DataArrayDouble([[0.,2.,2.],[0.,0.,2.],[2.,0.,2.],[2.,2.,2.],[0.,2.,0.],[0.,0.,0.],[2.,0.,0.],[2.,2.,0.], [0.,1.,2.],[1.,0.,2.],[2.,1.,2.],[1.,2.,2.], [0.,1.,0.],[1.,0.,0.],[2.,1.,0.],[1.,2.,0.], [0.,2.,1.],[0.,0.,1.],[2.,0.,1.],[2.,2.,1.], [1.,1.,2.], [0.,1.,1.],[1.,0.,1.],[2.,1.,1.],[1.,2.,1.], [1.,1.,0.], [1.,1.,1.]])
        m=MEDCouplingUMesh("mesh",3) ; m.setCoords(coo)
        m.allocateCells()
        # the cell description is exactly those described in the description of HEXA27 in MED file 3.0.7 documentation
        m.insertNextCell(NORM_HEXA27,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26])
        refCoo=[-1.,-1.,-1.,-1.,1.,-1.,1.,1.,-1.,1.,-1.,-1.,-1.,-1.,1.,-1.,1.,1.,1.,1.,1.,1.,-1.,1.,-1.,0.,-1.,0.,1.,-1.,1.,0.,-1.,0.,-1.,-1.,-1.,0.,1.,0.,1.,1.,1.,0.,1.,0.,-1.,1.,-1.,-1.,0.,-1.,1.,0.,1.,1.,0.,1.,-1.,0.,0.,0.,-1.,-1.,0.,0.,0.,1.,0.,1.,0.,0.,0.,-1.,0.,0.,0.,1.,0.,0.,0.]
        weights=[0.1714677640603571,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.43895747599451346,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.43895747599451346,0.27434842249657115,0.43895747599451346,0.7023319615912209,0.43895747599451346,0.27434842249657115,0.43895747599451346,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.43895747599451346,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.1714677640603571]
        gCoords=[-0.774596669241483,-0.774596669241483,-0.774596669241483,-0.774596669241483,-0.774596669241483,0.0,-0.774596669241483,-0.774596669241483,0.774596669241483,-0.774596669241483,0.0,-0.774596669241483,-0.774596669241483,0.0,0.0,-0.774596669241483,0.0,0.774596669241483,-0.774596669241483,0.774596669241483,-0.774596669241483,-0.774596669241483,0.774596669241483,0.0,-0.774596669241483,0.774596669241483,0.774596669241483,0.0,-0.774596669241483,-0.774596669241483,0.0,-0.774596669241483,0.0,0.0,-0.774596669241483,0.774596669241483,0.0,0.0,-0.774596669241483,0.0,0.0,0.0,0.0,0.0,0.774596669241483,0.0,0.774596669241483,-0.774596669241483,0.0,0.774596669241483,0.0,0.0,0.774596669241483,0.774596669241483,0.774596669241483,-0.774596669241483,-0.774596669241483,0.774596669241483,-0.774596669241483,0.0,0.774596669241483,-0.774596669241483,0.774596669241483,0.774596669241483,0.0,-0.774596669241483,0.774596669241483,0.0,0.0,0.774596669241483,0.0,0.774596669241483,0.774596669241483,0.774596669241483,-0.774596669241483,0.774596669241483,0.774596669241483,0.0,0.774596669241483,0.774596669241483,0.774596669241483]
        fGauss=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fGauss.setName("fGauss")
        fGauss.setMesh(m)
        fGauss.setGaussLocalizationOnType(NORM_HEXA27,refCoo,gCoords,weights)
        arr=DataArrayDouble(fGauss.getNumberOfTuplesExpected()) ; arr.iota()
        fGauss.setArray(arr)
        arrOfDisc=fGauss.getLocalizationOfDiscr()
        # the test is here
        self.assertTrue(arrOfDisc.isEqual(DataArrayDouble([0.2254033307585172,1.7745966692414836,1.7745966692414834,0.22540333075851715,1.7745966692414834,1.,0.22540333075851715,1.7745966692414836,0.22540333075851715,0.22540333075851715,1.,1.7745966692414834,0.2254033307585171,1.,1.,0.22540333075851715,1.0000000000000002,0.2254033307585171,0.22540333075851715,0.22540333075851715,1.7745966692414838,0.22540333075851715,0.22540333075851715,1.,0.22540333075851715,0.22540333075851715,0.22540333075851715,1.,1.7745966692414832,1.7745966692414834,1.,1.774596669241483,1.,1.0000000000000002,1.7745966692414832,0.22540333075851712,1.,1.,1.774596669241483,1.,1.,1.,1.,1.,0.2254033307585171,1.,0.22540333075851715,1.7745966692414834,1.,0.2254033307585171,1.,1.0000000000000002,0.22540333075851715,0.2254033307585171,1.7745966692414834,1.7745966692414834,1.7745966692414836,1.7745966692414832,1.7745966692414834,1.0000000000000002,1.7745966692414834,1.7745966692414836,0.22540333075851712,1.7745966692414832,1.,1.7745966692414834,1.774596669241483,1.,1.,1.7745966692414832,1.0000000000000002,0.22540333075851712,1.7745966692414836,0.22540333075851715,1.7745966692414836,1.7745966692414832,0.22540333075851715,1.,1.7745966692414836,0.22540333075851715,0.22540333075851715],27,3),1e-12))
        #
        weights=27*[1]
        gCoords=refCoo
        fGauss.setGaussLocalizationOnType(NORM_HEXA27,refCoo,gCoords,weights)
        arrOfDisc2=fGauss.getLocalizationOfDiscr()
        self.assertTrue(arrOfDisc2.isEqual(coo,1e-12))
        pass

    def testSwig2Pyra13GP1(self):
        coo=DataArrayDouble([[0.,2.,0.],[2.,2.,0.],[2.,0.,0.],[0.,0.,0.],[1.,1.,2.],[1.,2.,0.],[2.,1.,0.],[1.,0.,0.],[0.,1.,0.],[0.5,1.5,1.],[1.5,1.5,1.],[1.5,0.5,1.],[0.5,0.5,1.]])
        m=MEDCouplingUMesh("mesh",3) ; m.setCoords(coo)
        m.allocateCells()
        # the cell description is exactly those described in the description of PYRA13 in MED file 3.0.7 documentation
        m.insertNextCell(NORM_PYRA13,[0,1,2,3,4,5,6,7,8,9,10,11,12])
        refCoords=[1.,0.,0.,0.,-1.,0.,-1.,0.,0.,0.,1.,0.,0.,0.,1.,0.5,-0.5,0.,-0.5,-0.5,0.,-0.5,0.5,0.,0.5,0.5,0.,0.5,0.,0.5,0.,-0.5,0.5,-0.5,0.,0.5,0.,0.5,0.5]
        gaussCoords=[0.,0.,0.5,0.21210450275,0.21210450275,0.5,-0.21210450275,0.21210450275,0.5,-0.21210450275,-0.21210450275,0.5,0.21210450275,-0.21210450275,0.5,0.,0.,0.07579099449999999,0.,0.,0.9242090055000001,0.5394929090572634,0.,0.17359176399999998,0.,0.5394929090572634,0.17359176399999998,-0.5394929090572634,0.,0.17359176399999998,0.,-0.5394929090572634,0.17359176399999998,0.1133235629427366,0.,0.826408236,0.,0.1133235629427366,0.826408236,-0.1133235629427366,0.,0.826408236,0.,-0.1133235629427366,0.826408236,0.5826406005183961,0.5826406005183961,-0.053206449499999975,-0.5826406005183961,0.5826406005183961,-0.053206449499999975,-0.5826406005183961,-0.5826406005183961,-0.053206449499999975,0.5826406005183961,-0.5826406005183961,-0.053206449499999975,0.5532064495,0.,0.5,0.,0.5532064495,0.5,-0.5532064495,0.,0.5,0.,-0.5532064495,0.5,-0.029434151018396033,-0.029434151018396033,1.0532064495,0.029434151018396033,-0.029434151018396033,1.0532064495,0.029434151018396033,0.029434151018396033,1.0532064495,-0.029434151018396033,0.029434151018396033,1.0532064495]
        weights=[0.0492545926875,0.031210562625,0.031210562625,0.031210562625,0.031210562625,0.10663554205740113,0.0007171281994273535,0.0816994048010844,0.0816994048010844,0.0816994048010844,0.0816994048010844,0.0036048554264914074,0.0036048554264914074,0.0036048554264914074,0.0036048554264914074,0.008958181586640837,0.008958181586640837,0.008958181586640837,0.008958181586640837,0.002018983875,0.002018983875,0.002018983875,0.002018983875,2.286237794882217e-05,2.286237794882217e-05,2.286237794882217e-05,2.286237794882217e-05]
        fGauss=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fGauss.setName("fGauss")
        fGauss.setMesh(m)
        fGauss.setGaussLocalizationOnType(NORM_PYRA13,refCoords,gaussCoords,weights)
        arr=DataArrayDouble(fGauss.getNumberOfTuplesExpected()) ; arr.iota()
        fGauss.setArray(arr)
        arrOfDisc=fGauss.getLocalizationOfDiscr()
        # the test is here
        self.assertTrue(arrOfDisc.isEqual(DataArrayDouble([1.,1.,1.,0.5757909945,1.,1.,1.,0.5757909945,1.,1.4242090055,1.,1.,1.,1.4242090055,1.,1.,1.,0.151581989,1.,1.,1.848418011,0.4605070909427367,1.5394929090572635,0.347183528,0.4605070909427367,0.4605070909427367,0.347183528,1.5394929090572638,0.4605070909427366,0.347183528,1.5394929090572635,1.5394929090572638,0.347183528,0.8866764370572636,1.1133235629427367,1.652816472,0.8866764370572636,0.8866764370572636,1.652816472,1.1133235629427367,0.8866764370572636,1.652816472,1.1133235629427365,1.1133235629427367,1.652816472,-0.16528120103679209,1.,-0.106412899,1.,-0.1652812010367921,-0.106412899,2.1652812010367914,1.,-0.106412899,1.,2.165281201036791,-0.106412899,0.4467935505,1.5532064495,1.,0.4467935505,0.4467935505,1.,1.5532064495,0.4467935505,1.,1.5532064495,1.5532064495,1.,1.0588683020367922,1.,2.106412899,1.,1.0588683020367922,2.106412899,0.9411316979632077,1.,2.106412899,1.,0.9411316979632078,2.106412899],27,3),1e-12))
        #
        weights=13*[1]
        gaussCoords=refCoords[:] ; gaussCoords[14]=0.9999999999999 # change z of point #4 0.999... instead of 1. because with shape function it leads to division by 0. !
        fGauss.setGaussLocalizationOnType(NORM_PYRA13,refCoords,gaussCoords,weights)
        arrOfDisc2=fGauss.getLocalizationOfDiscr()
        self.assertTrue(arrOfDisc2.isEqual(coo,1e-10)) # be less exigent 1e-10 instead of 1e-12 due to shape function sensitivity arount 0.,0.,1. !
        pass

    def testSwig2Tri7GP1(self):
        coo=DataArrayDouble([[0,0],[0,2],[2,0],[0,1],[1,1],[1,0],[0.6666666666666667,0.6666666666666667]])
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo)
        m.allocateCells()
        # the cell description is exactly those described in the description of TRI7 in MED file 3.0.7 documentation
        m.insertNextCell(NORM_TRI7, list(range(7)))
        refCoords=[0.,0.,1.,0.,0.,1.,0.5,0.,0.5,0.5,0.,0.5,0.3333333333333333,0.3333333333333333]
        gaussCoords=[0.3333333333333333,0.3333333333333333,0.470142064105115,0.470142064105115,0.05971587178977,0.470142064105115,0.470142064105115,0.05971587178977,0.101286507323456,0.101286507323456,0.797426985353088,0.101286507323456,0.101286507323456,0.797426985353088]
        weights=[0.062969590272413,0.062969590272413,0.062969590272413,0.066197076394253,0.066197076394253,0.066197076394253,0.1125]
        fGauss=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fGauss.setName("fGauss")
        fGauss.setMesh(m)
        fGauss.setGaussLocalizationOnType(NORM_TRI7,refCoords,gaussCoords,weights)
        arr=DataArrayDouble(fGauss.getNumberOfTuplesExpected()) ; arr.iota()
        fGauss.setArray(arr)
        arrOfDisc=fGauss.getLocalizationOfDiscr()
        self.assertTrue(arrOfDisc.isEqual(DataArrayDouble([0.666666666666667,0.666666666666667,0.9402841282102293,0.9402841282102293,0.9402841282102299,0.11943174357954002,0.11943174357953992,0.9402841282102299,0.20257301464691194,0.20257301464691196,0.20257301464691205,1.5948539707061757,1.5948539707061757,0.20257301464691202],7,2),1e-12))
        #
        weights=7*[1]
        gaussCoords=refCoords
        fGauss.setGaussLocalizationOnType(NORM_TRI7,refCoords,gaussCoords,weights)
        arrOfDisc2=fGauss.getLocalizationOfDiscr()
        self.assertTrue(arrOfDisc2.isEqual(coo,1e-12))
        pass

    def testSwig2StructuredDesc1(self):
        c=MEDCouplingCMesh()
        arr0=DataArrayDouble(3) ; arr0.iota()
        arr1=DataArrayDouble(4) ; arr1.iota()
        arr2=DataArrayDouble(5) ; arr2.iota()
        c.setCoords(arr0,arr1,arr2)
        #
        self.assertEqual(98,c.getNumberOfCellsOfSubLevelMesh())
        m=c.build1SGTSubLevelMesh()
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([0,12,15,3,12,24,27,15,24,36,39,27,36,48,51,39,3,15,18,6,15,27,30,18,27,39,42,30,39,51,54,42,6,18,21,9,18,30,33,21,30,42,45,33,42,54,57,45,1,13,16,4,13,25,28,16,25,37,40,28,37,49,52,40,4,16,19,7,16,28,31,19,28,40,43,31,40,52,55,43,7,19,22,10,19,31,34,22,31,43,46,34,43,55,58,46,2,14,17,5,14,26,29,17,26,38,41,29,38,50,53,41,5,17,20,8,17,29,32,20,29,41,44,32,41,53,56,44,8,20,23,11,20,32,35,23,32,44,47,35,44,56,59,47,0,12,13,1,12,24,25,13,24,36,37,25,36,48,49,37,1,13,14,2,13,25,26,14,25,37,38,26,37,49,50,38,3,15,16,4,15,27,28,16,27,39,40,28,39,51,52,40,4,16,17,5,16,28,29,17,28,40,41,29,40,52,53,41,6,18,19,7,18,30,31,19,30,42,43,31,42,54,55,43,7,19,20,8,19,31,32,20,31,43,44,32,43,55,56,44,9,21,22,10,21,33,34,22,33,45,46,34,45,57,58,46,10,22,23,11,22,34,35,23,34,46,47,35,46,58,59,47,0,1,4,3,3,4,7,6,6,7,10,9,1,2,5,4,4,5,8,7,7,8,11,10,12,13,16,15,15,16,19,18,18,19,22,21,13,14,17,16,16,17,20,19,19,20,23,22,24,25,28,27,27,28,31,30,30,31,34,33,25,26,29,28,28,29,32,31,31,32,35,34,36,37,40,39,39,40,43,42,42,43,46,45,37,38,41,40,40,41,44,43,43,44,47,46,48,49,52,51,51,52,55,54,54,55,58,57,49,50,53,52,52,53,56,55,55,56,59,58])))
        self.assertEqual(NORM_QUAD4,m.getCellModelEnum())
        #
        self.assertTrue(MEDCouplingStructuredMesh.Build1GTNodalConnectivityOfSubLevelMesh([3,7]).isEqual(DataArrayInt([0,3,3,6,6,9,9,12,12,15,15,18,1,4,4,7,7,10,10,13,13,16,16,19,2,5,5,8,8,11,11,14,14,17,17,20,0,1,1,2,3,4,4,5,6,7,7,8,9,10,10,11,12,13,13,14,15,16,16,17,18,19,19,20])))
        pass

    def testSwig2Colinearize2D1(self):
        coo=DataArrayDouble([-5.,0.,-1.,0.,4.,3.,7.,0.,1.,6.,1.,0.,-3.,0.,6.,1.,5.,0.,3.,0.],10,2)
        #
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_POLYGON,[5,9,8,3,7,2,4,0,6,1])
        refPtr=m.getCoords().getHiddenCppPointer()
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([0])))
        self.assertEqual(refPtr,m.getCoords().getHiddenCppPointer())
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([5,0,3,4])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4])))
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([])))
        self.assertEqual(refPtr,m.getCoords().getHiddenCppPointer())
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([5,0,3,4])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4])))
        #
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_POLYGON,[8,3,7,2,4,0,6,1,5,9])
        refPtr=m.getCoords().getHiddenCppPointer()
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([0])))
        self.assertEqual(refPtr,m.getCoords().getHiddenCppPointer())
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([5,0,3,4])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4])))
        #
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_POLYGON,[3,7,2,4,0,6,1,5,9,8])
        refPtr=m.getCoords().getHiddenCppPointer()
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([0])))
        self.assertEqual(refPtr,m.getCoords().getHiddenCppPointer())
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([5,3,4,0])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4])))
        #
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_POLYGON,[4,0,6,1,5,9,8,3,7,2,])
        refPtr=m.getCoords().getHiddenCppPointer()
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([0])))
        self.assertEqual(refPtr,m.getCoords().getHiddenCppPointer())
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([5,4,0,3])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4])))
        ## false quadratic
        coo2=DataArrayDouble([(-5,0),(-1,0),(4,3),(7,0),(1,6),(1,0),(-3,0),(6,1),(5,0),(3,0),(2,0),(4,0),(6,0),(6.5,0.5),(5,2),(2.5,4.5),(-2,3),(-4,0),(-2,0),(0,0)])
        coo2.setInfoOnComponents(["aa","bbbb"])
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo2) ; m.allocateCells()
        m.insertNextCell(NORM_QPOLYG,[5,9,8,3,7,2,4,0,6,1,10,11,12,13,14,15,16,17,18,19])
        refPtr=m.getCoords().getHiddenCppPointer()
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([0])))
        self.assertNotEqual(refPtr,m.getCoords().getHiddenCppPointer())#not same coordinates here
        self.assertEqual(["aa","bbbb"],m.getCoords().getInfoOnComponents())
        refPtr=m.getCoords().getHiddenCppPointer()
        self.assertTrue(coo2.isEqual(m.getCoords()[:20],1e-12))
        self.assertTrue(m.getCoords()[20:].isEqualWithoutConsideringStr(DataArrayDouble([(1.,0.),(4.,3.)]),1e-12))
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([32,0,3,4,20,21,16])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,7])))
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([])))
        self.assertEqual(refPtr,m.getCoords().getHiddenCppPointer())
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([32,0,3,4,20,21,16])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,7])))
        # mix of quadratic and linear inside a QPOLYG cell
        coo2=DataArrayDouble([(-5,0),(-1,0),(7.,6.),(7,0),(1,6),(1,0),(-3,0),(8.2426406871192839,3),(5,0),(3,0),  (2,0),(4,0),(6,0),(7.9196888946291288,1.3764116995614091),(7.9196888946291288,4.6235883004385911),(4,7.2426406871192848),(-2,3),(-4,0),(-2,0),(0,0)])
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo2) ; m.allocateCells()
        m.insertNextCell(NORM_QPOLYG,[5,9,8,3,7,2,4,0,6,1,10,11,12,13,14,15,16,17,18,19])
        refPtr=m.getCoords().getHiddenCppPointer()
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([0])))
        self.assertNotEqual(refPtr,m.getCoords().getHiddenCppPointer())#not same coordinates here
        self.assertTrue(coo2.isEqual(m.getCoords()[:20],1e-12))
        self.assertTrue(m.getCoords()[20:].isEqual(DataArrayDouble([(1.,0.),(7.,6.)]),1e-12))
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([32,0,3,4,20,21,16])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,7])))
        pass

    def testSwig2BoundingBoxForBBTree1(self):
        """ This test appears simple but it checks that bounding box are correctly computed for quadratic polygons. It can help a lot to reduce the amount of intersections !
        """
        coo=DataArrayDouble([-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5,-0.5,0.45,0.,0.3181980515339464,0.31819805153394637,0.,0.45,-0.31819805153394637,0.3181980515339464,-0.45,0.,-0.3181980515339465,-0.31819805153394637,0.,-0.45,0.3181980515339463,-0.3181980515339465,-0.5,0.0,0.0,0.5,0.5,0.0,0.0,-0.5,-0.4090990257669732,-0.4090990257669732,0.40909902576697316,-0.4090990257669732],18,2)
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_QPOLYG,[0,1,2,3,11,5,7,9,12,13,14,17,4,6,8,16])
        m.insertNextCell(NORM_QPOLYG,[3,0,9,11,15,16,10,17])
        self.assertTrue(m.getBoundingBoxForBBTree().isEqual(DataArrayDouble([-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,-0.31819805153394637],2,4),1e-12))
        pass

    def testSwig2CartBuildUnstructuredOnExoticCases1(self):
        """ Test focusing on traduction from cartesian to unstructured mesh when spaceDim greater than meshDim.
        """
        #
        m=MEDCouplingCMesh()
        arrX=DataArrayDouble(3) ; arrX.iota()
        arrY=DataArrayDouble(4) ; arrY.iota()
        arrZ=DataArrayDouble(1) ; arrZ.iota()
        m.setCoords(arrX,arrY,arrZ)
        self.assertEqual(2,m.getMeshDimension())
        self.assertEqual(3,m.getSpaceDimension())
        mu=m.buildUnstructured()
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([4,1,0,3,4,4,2,1,4,5,4,4,3,6,7,4,5,4,7,8,4,7,6,9,10,4,8,7,10,11])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,10,15,20,25,30])))
        coo0=DataArrayDouble([(0,0,0),(1,0,0),(2,0,0),(0,1,0),(1,1,0),(2,1,0),(0,2,0),(1,2,0),(2,2,0),(0,3,0),(1,3,0),(2,3,0)])
        self.assertTrue(mu.getCoords().isEqual(coo0,1e-12))
        #
        m=MEDCouplingCMesh()
        arrX=DataArrayDouble(3) ; arrX.iota()
        arrY=DataArrayDouble(1) ; arrY.iota()
        arrZ=DataArrayDouble(4) ; arrZ.iota()
        m.setCoords(arrX,arrY,arrZ)
        self.assertEqual(2,m.getMeshDimension())
        self.assertEqual(3,m.getSpaceDimension())
        mu=m.buildUnstructured()
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([4,1,0,3,4,4,2,1,4,5,4,4,3,6,7,4,5,4,7,8,4,7,6,9,10,4,8,7,10,11])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,10,15,20,25,30])))
        coo1=DataArrayDouble([(0,0,0),(1,0,0),(2,0,0),(0,0,1),(1,0,1),(2,0,1),(0,0,2),(1,0,2),(2,0,2),(0,0,3),(1,0,3),(2,0,3)])
        self.assertTrue(mu.getCoords().isEqual(coo1,1e-12))
        #
        m=MEDCouplingCMesh()
        arrX=DataArrayDouble(1) ; arrX.iota() ; arrX+=9
        arrY=DataArrayDouble(3) ; arrY.iota()
        arrZ=DataArrayDouble(4) ; arrZ.iota()
        m.setCoords(arrX,arrY,arrZ)
        self.assertEqual(2,m.getMeshDimension())
        self.assertEqual(3,m.getSpaceDimension())
        mu=m.buildUnstructured()
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([4,1,0,3,4,4,2,1,4,5,4,4,3,6,7,4,5,4,7,8,4,7,6,9,10,4,8,7,10,11])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,10,15,20,25,30])))
        coo2=DataArrayDouble([(9,0,0),(9,1,0),(9,2,0),(9,0,1),(9,1,1),(9,2,1),(9,0,2),(9,1,2),(9,2,2),(9,0,3),(9,1,3),(9,2,3)])
        self.assertTrue(mu.getCoords().isEqual(coo2,1e-12))
        #
        m=MEDCouplingCMesh()
        arrX=DataArrayDouble(3) ; arrX.iota()
        arrY=DataArrayDouble(1) ; arrY.iota(7)
        arrZ=DataArrayDouble(1) ; arrZ.iota(8)
        m.setCoords(arrX,arrY,arrZ)
        self.assertEqual(1,m.getMeshDimension())
        self.assertEqual(3,m.getSpaceDimension())
        mu=m.buildUnstructured()
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([1,0,1,1,1,2])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6])))
        coo3=DataArrayDouble([(0,7,8),(1,7,8),(2,7,8)])
        self.assertTrue(mu.getCoords().isEqual(coo3,1e-12))
        #
        m=MEDCouplingCMesh()
        arrX=DataArrayDouble(1) ; arrX.iota(7)
        arrY=DataArrayDouble(1) ; arrY.iota(8)
        arrZ=DataArrayDouble(3) ; arrZ.iota()
        m.setCoords(arrX,arrY,arrZ)
        self.assertEqual(1,m.getMeshDimension())
        self.assertEqual(3,m.getSpaceDimension())
        mu=m.buildUnstructured()
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([1,0,1,1,1,2])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6])))
        coo4=DataArrayDouble([(7,8,0),(7,8,1),(7,8,2)])
        self.assertTrue(mu.getCoords().isEqual(coo4,1e-12))
        #
        m=MEDCouplingCMesh()
        arrX=DataArrayDouble(3) ; arrX.iota()
        arrY=DataArrayDouble(1) ; arrY.iota(7)
        m.setCoords(arrX,arrY)
        self.assertEqual(1,m.getMeshDimension())
        self.assertEqual(2,m.getSpaceDimension())
        mu=m.buildUnstructured()
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([1,0,1,1,1,2])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6])))
        coo5=DataArrayDouble([(0,7),(1,7),(2,7)])
        self.assertTrue(mu.getCoords().isEqual(coo5,1e-12))
        #
        m=MEDCouplingCMesh()
        arrX=DataArrayDouble(1) ; arrX.iota(7)
        arrY=DataArrayDouble(3) ; arrY.iota()
        m.setCoords(arrX,arrY)
        self.assertEqual(1,m.getMeshDimension())
        self.assertEqual(2,m.getSpaceDimension())
        mu=m.buildUnstructured()
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([1,0,1,1,1,2])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6])))
        coo6=DataArrayDouble([(7,0),(7,1),(7,2)])
        self.assertTrue(mu.getCoords().isEqual(coo6,1e-12))
        pass

    def testSwig2Colinearize2D2(self):
        """ simple non regression test but that has revealed a bug"""
        coo=DataArrayDouble([(0,0),(0,0.5),(0,1),(1,1),(1,0),(0.5,0)])
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo)
        m.allocateCells() ; m.insertNextCell(NORM_POLYGON,[0,1,2,3,4,5])
        m.checkConsistency()
        refPtr=m.getCoords().getHiddenCppPointer()
        #
        m.colinearize2D(1e-12)
        m.checkConsistency()
        self.assertEqual(refPtr,m.getCoords().getHiddenCppPointer())
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([NORM_POLYGON,0,2,3,4])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5])))
        pass

    def testSwig2Colinearize2D3(self):
        """ colinearize was too agressive, potentially producing cells with one edge """
        # Flat polygon  with 3 edges - nothing should happen (min number of edges for a linear polyg)
        coo = DataArrayDouble([0.0,0.0,  2.0,0.0,   1.5,0.0,  1.0,0.0,  0.5,0.0], 5,2)
        m = MEDCouplingUMesh("m", 2)
        c, cI = [DataArrayInt(l) for l in [[NORM_POLYGON, 0,1,2], [0,4]] ]
        m.setCoords(coo); m.setConnectivity(c, cI)
        m.colinearize2D(1e-10)
        m.checkConsistency()
        self.assertEqual(c.getValues(), m.getNodalConnectivity().getValues())
        self.assertEqual(cI.getValues(), m.getNodalConnectivityIndex().getValues())

        # Flat quad polygon, 2 edges - nothing should happen (min number of edges for a quad polyg)
        m = MEDCouplingUMesh("m", 2)
        c, cI = [DataArrayInt(l) for l in [[NORM_QPOLYG, 0,1,  2,3], [0,5]] ]
        m.setCoords(coo); m.setConnectivity(c, cI)
        m.colinearize2D(1e-10)
        m.checkConsistency()
        self.assertEqual(c.getValues(), m.getNodalConnectivity().getValues())
        self.assertEqual(cI.getValues(), m.getNodalConnectivityIndex().getValues())

        # Flat polygon, 4 edges - one reduction should happen
        m = MEDCouplingUMesh("m", 2)
        c, cI = [DataArrayInt(l) for l in [[NORM_POLYGON, 0,1,2,3], [0,5]] ]
        m.setCoords(coo); m.setConnectivity(c, cI)
        m.colinearize2D(1e-10)
        m.checkConsistency()
        self.assertEqual([NORM_POLYGON, 3,1,2], m.getNodalConnectivity().getValues())
        self.assertEqual([0,4], m.getNodalConnectivityIndex().getValues())

        # Flat quad polygon, 3 edges - one reduction expected
        m = MEDCouplingUMesh("m", 2)
        c, cI = [DataArrayInt(l) for l in [[NORM_QPOLYG, 0,1,3,  3,2,4], [0,7]] ]
        m.setCoords(coo); m.setConnectivity(c, cI)
        m.colinearize2D(1e-10)
        m.checkConsistency()
        self.assertEqual([NORM_QPOLYG, 3,1, 5,2], m.getNodalConnectivity().getValues())
        self.assertTrue( m.getCoords()[5].isEqual( DataArrayDouble([(1.5,0.0)]), 1.0e-12 ) )
        self.assertEqual([0,5], m.getNodalConnectivityIndex().getValues())

        # Now an actual (neutronic) case: circle made of 4 SEG3. Should be reduced to 2 SEG3
        m = MEDCouplingDataForTest.buildCircle2(0.0, 0.0, 1.0)
        c, cI = [DataArrayInt(l) for l in [[NORM_QPOLYG, 7,5,3,1,  6,4,2,0], [0,9]] ]
        m.colinearize2D(1e-10)
        m.checkConsistency()
        self.assertEqual([NORM_QPOLYG, 3,5,  8,4], m.getNodalConnectivity().getValues())
        self.assertTrue( m.getCoords()[8].isEqual( DataArrayDouble([(1.0,0.0)]), 1.0e-12 ) )
        self.assertEqual([0,5], m.getNodalConnectivityIndex().getValues())

    def testSwig2CheckAndPreparePermutation2(self):
        a=DataArrayInt([10003,9999999,5,67])
        self.assertTrue(DataArrayInt.CheckAndPreparePermutation(a).isEqual(DataArrayInt([2,3,0,1])))
        a=DataArrayInt([10003,-9999999,5,67])
        self.assertTrue(DataArrayInt.CheckAndPreparePermutation(a).isEqual(DataArrayInt([3,0,1,2])))
        a=DataArrayInt([])
        self.assertTrue(DataArrayInt.checkAndPreparePermutation(a).isEqual(DataArrayInt([])))
        pass

    def testSwig2ComputeNeighborsOfNodes1(self):
        arrX=DataArrayDouble(3) ; arrX.iota()
        arrY=DataArrayDouble(4) ; arrY.iota()
        arrZ=DataArrayDouble(5) ; arrZ.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY,arrZ) ; m=m.buildUnstructured()
        # 3D
        a,b=m.computeNeighborsOfNodes()
        self.assertTrue(a.isEqual(DataArrayInt([1,3,12,0,4,13,2,1,5,14,0,4,15,6,3,1,16,5,7,4,2,17,8,3,7,18,9,6,4,19,8,10,7,5,20,11,6,10,21,9,7,22,11,10,8,23,13,15,0,24,12,16,1,14,25,13,17,2,26,12,16,3,18,27,15,13,4,17,19,28,16,14,5,20,29,15,19,6,21,30,18,16,7,20,22,31,19,17,8,23,32,18,22,9,33,21,19,10,23,34,22,20,11,35,25,27,12,36,24,28,13,26,37,25,29,14,38,24,28,15,30,39,27,25,16,29,31,40,28,26,17,32,41,27,31,18,33,42,30,28,19,32,34,43,31,29,20,35,44,30,34,21,45,33,31,22,35,46,34,32,23,47,37,39,24,48,36,40,25,38,49,37,41,26,50,36,40,27,42,51,39,37,28,41,43,52,40,38,29,44,53,39,43,30,45,54,42,40,31,44,46,55,43,41,32,47,56,42,46,33,57,45,43,34,47,58,46,44,35,59,49,51,36,48,52,37,50,49,53,38,48,52,39,54,51,49,40,53,55,52,50,41,56,51,55,42,57,54,52,43,56,58,55,53,44,59,54,58,45,57,55,46,59,58,56,47])))
        self.assertTrue(b.isEqual(DataArrayInt([0,3,7,10,14,19,23,27,32,36,39,43,46,50,55,59,64,70,75,80,86,91,95,100,104,108,113,117,122,128,133,138,144,149,153,158,162,166,171,175,180,186,191,196,202,207,211,216,220,223,227,230,234,239,243,247,252,256,259,263,266])))
        # 2D
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY) ; m=m.buildUnstructured()
        a,b=m.computeNeighborsOfNodes()
        self.assertTrue(a.isEqual(DataArrayInt([1,3,0,4,2,1,5,0,4,6,3,1,5,7,4,2,8,3,7,9,6,4,8,10,7,5,11,6,10,9,7,11,10,8])))
        self.assertTrue(b.isEqual(DataArrayInt([0,2,5,7,10,14,17,20,24,27,29,32,34])))
        # 1D
        m=m.buildDescendingConnectivity()[0]
        a,b=m.computeNeighborsOfNodes()
        self.assertTrue(a.isEqual(DataArrayInt([1,3,0,4,2,1,5,0,4,6,3,1,5,7,4,2,8,3,7,9,6,4,8,10,7,5,11,6,10,9,7,11,10,8])))
        self.assertTrue(b.isEqual(DataArrayInt([0,2,5,7,10,14,17,20,24,27,29,32,34])))
        pass

    def testSwigBugOnUnpackingTuplesInDataArray1(self):
        inp=DataArrayDouble([(1,2,3),(4,5,6),(7,8,9),(10,11,12)])
        it=inp.__iter__()
        r = next(it)
        self.assertRaises(StopIteration,r.__getitem__,4)
        self.assertEqual(len(r),3)
        a,b,c=r
        r = next(it)
        self.assertEqual(len(r),3)
        d,e,f=r
        r = next(it)
        self.assertEqual(len(r),3)
        g,h,i=r
        r = next(it)
        self.assertEqual(len(r),3)
        j,k,l=r
        self.assertTrue(inp.isEqual(DataArrayDouble([a,b,c,d,e,f,g,h,i,j,k,l],4,3),1e-12))
        ########
        inp=DataArrayInt([(1,2,3),(4,5,6),(7,8,9),(10,11,12)])
        it=inp.__iter__()
        r = next(it)
        self.assertRaises(StopIteration,r.__getitem__,4)
        self.assertEqual(len(r),3)
        a,b,c=r
        r = next(it)
        self.assertEqual(len(r),3)
        d,e,f=r
        r = next(it)
        self.assertEqual(len(r),3)
        g,h,i=r
        r = next(it)
        self.assertEqual(len(r),3)
        j,k,l=r
        self.assertTrue(inp.isEqual(DataArrayInt([a,b,c,d,e,f,g,h,i,j,k,l],4,3)))
        pass

    def testSwig2IMesh1(self):
        """ 1st test of image grid mesh.
        """
        m=MEDCouplingIMesh()
        self.assertEqual(m.getSpaceDimension(),-1)
        self.assertEqual(1,len(m.__repr__().split("\n")))
        self.assertEqual(6,len(m.__str__().split("\n")))
        self.assertRaises(InterpKernelException,m.getNodeStruct)
        self.assertRaises(InterpKernelException,m.getOrigin)
        self.assertRaises(InterpKernelException,m.getDXYZ)
        m.setSpaceDimension(3)
        self.assertEqual(9,len(m.__str__().split("\n")))
        self.assertEqual(4,len(m.__repr__().split("\n")))
        self.assertEqual((0,0,0),m.getNodeStruct())
        self.assertEqual((0.,0.,0.),m.getOrigin())
        self.assertEqual((0.,0.,0.),m.getDXYZ())
        self.assertRaises(InterpKernelException,m.setNodeStruct,[3,4])
        m.setNodeStruct([3,4,2])
        self.assertEqual((3,4,2),m.getNodeStruct())
        m.setOrigin(DataArrayDouble([1.5,2.5,3.5]))
        self.assertEqual((1.5,2.5,3.5),m.getOrigin())
        m.setDXYZ((0.5,1.,0.25))
        self.assertEqual((0.5,1.,0.25),m.getDXYZ())
        for it in DataArrayDouble([(1.5,2.5,3.5)]):
            m2=MEDCouplingIMesh("",3,DataArrayInt([3,4,2]),it,DataArrayDouble((0.5,1.,0.25)))
            pass
        self.assertEqual(3,m.getSpaceDimension())
        self.assertEqual((3,4,2),m2.getNodeStruct())
        self.assertEqual((1.5,2.5,3.5),m2.getOrigin())
        self.assertEqual((0.5,1.,0.25),m2.getDXYZ())
        self.assertEqual(24,m2.getNumberOfNodes())
        self.assertEqual(6,m2.getNumberOfCells())
        self.assertTrue(m.isEqual(m2,1e-12)) ; self.assertTrue(m.isEqualWithoutConsideringStr(m2,1e-12))
        m2.setAxisUnit("m")
        self.assertTrue(not m.isEqual(m2,1e-12)) ; self.assertTrue(m.isEqualWithoutConsideringStr(m2,1e-12))
        m.setAxisUnit("m")
        self.assertTrue(m.isEqual(m2,1e-12)) ; self.assertTrue(m.isEqualWithoutConsideringStr(m2,1e-12))
        m.setName("mesh")
        self.assertTrue(not m.isEqual(m2,1e-12)) ; self.assertTrue(m.isEqualWithoutConsideringStr(m2,1e-12))
        m2.setName("mesh")
        self.assertTrue(m.isEqual(m2,1e-12)) ; self.assertTrue(m.isEqualWithoutConsideringStr(m2,1e-12))
        m2.setTime(1.1,0,3)
        self.assertTrue(not m.isEqual(m2,1e-12))
        m.setTime(1.1,0,3)
        self.assertTrue(m.isEqual(m2,1e-12))
        m.setTimeUnit("ms")
        self.assertTrue(not m.isEqual(m2,1e-12)) ; self.assertTrue(m.isEqualWithoutConsideringStr(m2,1e-12))
        m2.setTimeUnit("ms")
        self.assertTrue(m.isEqual(m2,1e-12)) ; self.assertTrue(m.isEqualWithoutConsideringStr(m2,1e-12))
        #
        m2.setNodeStruct([3,2,4])
        self.assertTrue(not m.isEqual(m2,1e-12))
        m.setNodeStruct([3,2,4])
        self.assertTrue(m.isEqual(m2,1e-12))
        m.setOrigin(DataArrayDouble([1.5,3.5,2.5]))
        self.assertTrue(not m.isEqual(m2,1e-12))
        m2.setOrigin([1.5,3.5,2.5])
        self.assertTrue(m.isEqual(m2,1e-12))
        m.setDXYZ((0.5,0.25,1.))
        self.assertTrue(not m.isEqual(m2,1e-12))
        m2.setDXYZ(DataArrayDouble((0.5,0.25,1.)))
        self.assertTrue(m.isEqual(m2,1e-12))
        m2bis=m2.deepCopy()
        self.assertTrue(m2bis.isEqual(m2,1e-12))
        #
        self.assertEqual(6,m2bis.getNumberOfCells())#3,2,4
        m2bis.refineWithFactor([3,3,3])
        self.assertEqual(162,m2bis.getNumberOfCells())
        self.assertEqual((7,4,10),m2bis.getNodeStruct())
        self.assertEqual((1.5,3.5,2.5),m2bis.getOrigin())
        self.assertTrue(DataArrayDouble([0.16666666666666666,0.08333333333333333,0.3333333333333333]).isEqual(DataArrayDouble(m2bis.getDXYZ()),1e-12))
        #
        self.assertEqual(3,m.getMeshDimension())
        self.assertAlmostEqual(0.125,m.getMeasureOfAnyCell(),16);
        mu=MEDCoupling1SGTUMesh(m.buildUnstructured())
        mu.checkConsistency()
        cooExp=DataArrayDouble([(1.5,3.5,2.5),(2,3.5,2.5),(2.5,3.5,2.5),(1.5,3.75,2.5),(2,3.75,2.5),(2.5,3.75,2.5),(1.5,3.5,3.5),(2,3.5,3.5),(2.5,3.5,3.5),(1.5,3.75,3.5),(2,3.75,3.5),(2.5,3.75,3.5),(1.5,3.5,4.5),(2,3.5,4.5),(2.5,3.5,4.5),(1.5,3.75,4.5),(2,3.75,4.5),(2.5,3.75,4.5),(1.5,3.5,5.5),(2,3.5,5.5),(2.5,3.5,5.5),(1.5,3.75,5.5),(2,3.75,5.5),(2.5,3.75,5.5)]) ; cooExp.setInfoOnComponents(["X [m]","Y [m]","Z [m]"])
        self.assertTrue(isinstance(mu,MEDCoupling1SGTUMesh))
        self.assertEqual(NORM_HEXA8,mu.getCellModelEnum())
        self.assertTrue(mu.getCoords().isEqual(cooExp,1e-12))
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([1,0,3,4,7,6,9,10,2,1,4,5,8,7,10,11,7,6,9,10,13,12,15,16,8,7,10,11,14,13,16,17,13,12,15,16,19,18,21,22,14,13,16,17,20,19,22,23])))
        bary=m.computeCellCenterOfMass()
        baryExp=DataArrayDouble([(1.75,3.625,3),(2.25,3.625,3),(1.75,3.625,4),(2.25,3.625,4),(1.75,3.625,5),(2.25,3.625,5)]) ; baryExp.setInfoOnComponents(["X [m]","Y [m]","Z [m]"])
        self.assertTrue(bary.isEqual(baryExp,1e-12))
        #
        c=m.convertToCartesian()
        c.checkConsistencyLight()
        self.assertEqual([1.1,0,3],c.getTime())
        self.assertEqual("ms",c.getTimeUnit())
        self.assertEqual(3,c.getMeshDimension())
        self.assertEqual(3,c.getSpaceDimension())
        arrX=DataArrayDouble([1.5,2.,2.5]) ; arrX.setInfoOnComponents(["X [m]"])
        self.assertTrue(c.getCoordsAt(0).isEqual(arrX,1e-12))
        arrY=DataArrayDouble([3.5,3.75]) ; arrY.setInfoOnComponents(["Y [m]"])
        self.assertTrue(c.getCoordsAt(1).isEqual(arrY,1e-12))
        arrZ=DataArrayDouble([2.5,3.5,4.5,5.5]) ; arrZ.setInfoOnComponents(["Z [m]"])
        self.assertTrue(c.getCoordsAt(2).isEqual(arrZ,1e-12))
        self.assertTrue(c.buildUnstructured().isEqual(m.buildUnstructured(),1e-12))
        #
        a,b=m.getCellsContainingPoints(baryExp,1e-12)
        self.assertTrue(a.isEqual(DataArrayInt([0,1,2,3,4,5])))
        self.assertTrue(b.isEqual(DataArrayInt([0,1,2,3,4,5,6])))
        for a,b in enumerate(baryExp):
            self.assertEqual(a,m.getCellContainingPoint(b,1e-12))
            pass
        #
        m.translate([1.,2.,4.])
        self.assertEqual((3,2,4),m.getNodeStruct())
        self.assertEqual((2.5,5.5,6.5),m.getOrigin())
        self.assertEqual((0.5,0.25,1.),m.getDXYZ())
        m.scale([0.,1.,3.],2.)
        self.assertAlmostEqual(1.,m.getMeasureOfAnyCell(),16);
        self.assertEqual((3,2,4),m.getNodeStruct())
        self.assertEqual((5.,10.,10.),m.getOrigin())
        self.assertEqual((1.,0.5,2.),m.getDXYZ())
        #
        f=m.getMeasureField(False)
        f2=MEDCouplingFieldDouble(ON_CELLS) ; f2.setTime(1.1,0,3) ; f2.setMesh(m) ; arr=DataArrayDouble(6) ; arr[:]=1. ; f2.setArray(arr) ; f2.setTimeUnit("ms")
        f2.setName("MeasureOfMesh_mesh")
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        #
        m3=m.buildStructuredSubPart([(1,2),(0,1),(1,3)])
        self.assertEqual((2,2,3),m3.getNodeStruct())
        self.assertEqual((6.,10.,12.),m3.getOrigin())
        self.assertEqual((1.,0.5,2.),m3.getDXYZ())
        # now playing with 3D surf
        m4=MEDCouplingIMesh("",3,DataArrayInt([3,1,4]),DataArrayDouble([1.5,2.5,3.5]),DataArrayDouble((0.5,1.,0.25))) ; m4.setAxisUnit("km")
        self.assertEqual([(1.5,2.5),(2.5,3.5),(3.5,4.25)],m4.getBoundingBox())
        self.assertEqual(3,m4.getSpaceDimension())
        self.assertEqual(2,m4.getMeshDimension())
        self.assertEqual(12,m4.getNumberOfNodes())
        self.assertEqual(6,m4.getNumberOfCells())
        mu=MEDCoupling1SGTUMesh(m4.buildUnstructured())
        mu.checkConsistency()
        self.assertTrue(isinstance(mu,MEDCoupling1SGTUMesh))
        self.assertEqual(NORM_QUAD4,mu.getCellModelEnum())
        coordsExp=DataArrayDouble([(1.5,2.5,3.5),(2,2.5,3.5),(2.5,2.5,3.5),(1.5,2.5,3.75),(2,2.5,3.75),(2.5,2.5,3.75),(1.5,2.5,4),(2,2.5,4),(2.5,2.5,4),(1.5,2.5,4.25),(2,2.5,4.25),(2.5,2.5,4.25)]) ; coordsExp.setInfoOnComponents(["X [km]","Y [km]","Z [km]"])
        self.assertTrue(mu.getCoords().isEqual(coordsExp,1e-12))
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([1,0,3,4,2,1,4,5,4,3,6,7,5,4,7,8,7,6,9,10,8,7,10,11])))
        pass

    def testSwig1GetValuesAsTuple1(self):
        d=DataArrayDouble()
        self.assertEqual(d.getValues(),[])
        self.assertEqual(d.getValuesAsTuple(),[])
        d=DataArrayDouble(24) ; d.iota() ; d.rearrange(3)
        self.assertEqual(d.getValues(),[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.])
        self.assertEqual(d.getValuesAsTuple(),[(0.,1.,2.0),(3.,4.,5.0),(6.,7.,8.0),(9.,10.,11.0),(12.,13.,14.0),(15.,16.,17.0),(18.,19.,20.0),(21.,22.,23.)])
        d=DataArrayInt()
        self.assertEqual(d.getValues(),[])
        self.assertEqual(d.getValuesAsTuple(),[])
        d=DataArrayInt(24) ; d.iota() ; d.rearrange(3)
        self.assertEqual(d.getValues(),[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23])
        self.assertEqual(d.getValuesAsTuple(),[(0,1,2),(3,4,5),(6,7,8),(9,10,11),(12,13,14),(15,16,17),(18,19,20),(21,22,23)])
        pass

    def testSwig2AMR1(self):
        self.assertEqual((1,3,12),MEDCouplingStructuredMesh.GetSplitVectFromStruct([3,4,5]))
        self.assertEqual((3,2),MEDCouplingStructuredMesh.GetDimensionsFromCompactFrmt([(1,4),(2,4)]))
        #
        amr=MEDCouplingCartesianAMRMesh("",2,[3,3],[0,0],[1,1])
        self.assertEqual(4,amr.getNumberOfCellsAtCurrentLevel())
        self.assertEqual(4,amr.getNumberOfCellsRecursiveWithOverlap())
        self.assertEqual(4,amr.getNumberOfCellsRecursiveWithoutOverlap())
        self.assertEqual(0,amr.getNumberOfPatches())
        self.assertEqual(1,amr.getMaxNumberOfLevelsRelativeToThis())
        self.assertEqual(2,amr.getSpaceDimension())
        amr.addPatch([(1,2),(0,1)],[4,4])
        self.assertEqual(4,amr.getNumberOfCellsAtCurrentLevel())
        self.assertEqual(20,amr.getNumberOfCellsRecursiveWithOverlap())
        self.assertEqual(19,amr.getNumberOfCellsRecursiveWithoutOverlap())
        self.assertEqual(1,amr.getNumberOfPatches())
        self.assertEqual(2,amr.getMaxNumberOfLevelsRelativeToThis())
        self.assertEqual(2,amr.getSpaceDimension())
        amr[0].addPatch([(2,3),(1,3)],[3,2])
        self.assertEqual(amr[0].getBLTRRange(),[(1,2),(0,1)])
        self.assertEqual(4,amr.getNumberOfCellsAtCurrentLevel())
        self.assertEqual(32,amr.getNumberOfCellsRecursiveWithOverlap())
        self.assertEqual(29,amr.getNumberOfCellsRecursiveWithoutOverlap())
        self.assertEqual(1,amr.getNumberOfPatches())
        self.assertEqual(3,amr.getMaxNumberOfLevelsRelativeToThis())
        self.assertEqual(2,amr.getSpaceDimension())
        amr[0].addPatch([(0,2),(3,4)],[3,2])
        self.assertEqual(16,amr[0].getMesh().getNumberOfCellsAtCurrentLevel())
        self.assertEqual(44,amr.getNumberOfCellsRecursiveWithOverlap())
        self.assertEqual(39,amr.getNumberOfCellsRecursiveWithoutOverlap())
        self.assertEqual(2,amr[0].getMesh().getNumberOfPatches())
        self.assertEqual(3,amr.getMaxNumberOfLevelsRelativeToThis())
        self.assertEqual(2,amr.getSpaceDimension())
        del amr[0][1]
        self.assertEqual(amr[0].getBLTRRange(),[(1,2),(0,1)])
        self.assertEqual(4,amr.getNumberOfCellsAtCurrentLevel())
        self.assertEqual(32,amr.getNumberOfCellsRecursiveWithOverlap())
        self.assertEqual(29,amr.getNumberOfCellsRecursiveWithoutOverlap())
        self.assertEqual(1,amr.getNumberOfPatches())
        self.assertEqual(3,amr.getMaxNumberOfLevelsRelativeToThis())
        self.assertEqual(2,amr.getSpaceDimension())
        pass

    def testSwig2NonRegressionTestPAL1164(self):
        """ Test PAL1164 Protection of applyLin against error in compoId ( #CEA22584 ) """
        xarr=DataArrayDouble(3,1)
        xarr.iota(0.)
        cmesh=MEDCouplingCMesh()
        cmesh.setCoords(xarr,xarr,xarr)
        mesh=cmesh.buildUnstructured()
        f=mesh.fillFromAnalytic(ON_CELLS,1,"(x-5.)*(x-5.)+(y-5.)*(y-5.)+(z-5.)*(z-5.)")
        f.setName("MyField")
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([60.75,52.75,52.75,44.75,52.75,44.75,44.75,36.75]),1e-12))
        self.assertRaises(InterpKernelException,f.applyLin,2.,0.,1)# compoId 1 whereas f has only one component !
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([60.75,52.75,52.75,44.75,52.75,44.75,44.75,36.75]),1e-12))
        f.applyLin(2.,0.,0)# here it is OK !
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([121.5,105.5,105.5,89.5,105.5,89.5,89.5,73.5]),1e-12))
        f.applyLin(2.,0.)
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([243.,211.,211.,179.,211.,179.,179.,147.]),1e-12))
        pass

    def testSwig2StructurizeMe1(self):
        arrx=DataArrayDouble(3) ; arrx.iota() ; arrx*=2.
        arry=DataArrayDouble(4) ; arry.iota() ; arry+=3.
        arrz=DataArrayDouble(5) ; arrz.iota() ; arrz*=0.5 ; arrz+=2.
        c=MEDCouplingCMesh() ; c.setCoords(arrx,arry,arrz)
        c.setName("mesh") ; c.setDescription("mesh descr") ; c.setTimeUnit("us") ; c.setTime(1.2,3,4)
        u=c.buildUnstructured()
        cp=DataArrayInt([3,5,6,1,0,9,8,7,12,11,16,10,17,23,22,21,19,20,18,14,13,2,4,15])
        np=DataArrayInt([3,33,5,35,6,36,1,31,0,30,9,39,8,38,7,37,12,42,11,41,16,46,10,40,17,47,23,53,22,52,21,51,19,49,20,50,18,48,14,44,13,43,2,32,4,34,15,45,29,59,28,58,27,57,26,56,25,55,24,54])
        u.renumberCells(cp)
        u.renumberNodes(np,len(np))
        u=MEDCoupling1SGTUMesh(u)
        #
        e,d,f=u.structurizeMe()
        self.assertTrue(c.isEqual(e,1e-12))
        self.assertTrue(d.isEqual(cp))
        self.assertTrue(f.isEqual(np))
        pass

    def testSwig2DenseMatrix1(self):
        m0=DenseMatrix(DataArrayDouble([2,3,4,5,1,6]),2,3)
        self.assertEqual(m0.getNumberOfRows(),2)
        self.assertEqual(m0.getNumberOfCols(),3)
        self.assertEqual(m0.getNbOfElems(),6)
        ref=m0.getData().getHiddenCppPointer()
        m00=m0.deepCopy()
        self.assertTrue(m0.isEqual(m00,1e-12))
        m00.getData().setIJ(0,0,2.1)
        self.assertTrue(not m0.isEqual(m00,1e-12))
        m00.getData().setIJ(0,0,2.)
        self.assertTrue(m0.isEqual(m00,1e-12))
        self.assertTrue(m0.getData().isEqual(DataArrayDouble([2,3,4,5,1,6]),1e-12))
        #
        m000=m0*DataArrayDouble([5,9,3])
        self.assertTrue(m000.getData().isEqual(DataArrayDouble([49.,52.]),1e-12))
        #
        m0.reShape(3,2)
        self.assertTrue(not m0.isEqual(m00,1e-12))
        self.assertEqual(m0.getNumberOfRows(),3)
        self.assertEqual(m0.getNumberOfCols(),2)
        self.assertEqual(ref,m0.getData().getHiddenCppPointer())
        self.assertTrue(m0.getData().isEqual(DataArrayDouble([2,3,4,5,1,6]),1e-12))
        m0.reShape(2,3)
        self.assertTrue(m0.isEqual(m00,1e-12))
        self.assertEqual(ref,m0.getData().getHiddenCppPointer())
        self.assertEqual(m0.getNumberOfRows(),2)
        self.assertEqual(m0.getNumberOfCols(),3)
        self.assertTrue(m0.getData().isEqual(DataArrayDouble([2,3,4,5,1,6]),1e-12))
        #m0np=m0.getData().toNumPyArray() ; m0np=matrix(m0np.reshape(m0.getNumberOfRows(),m0.getNumberOfCols()))
        m1=m0.deepCopy()
        self.assertEqual(m1.getNumberOfRows(),2)
        self.assertEqual(m1.getNumberOfCols(),3)
        self.assertTrue(m1.getData().isEqual(DataArrayDouble([2,3,4,5,1,6]),1e-12))
        m11=m0.deepCopy() ; m11+=m1
        self.assertEqual(m11.getNumberOfRows(),2)
        self.assertEqual(m11.getNumberOfCols(),3)
        self.assertTrue(m11.getData().isEqual(DataArrayDouble([4,6,8,10,2,12]),1e-12))
        m11=m11+m1
        self.assertEqual(m11.getNumberOfRows(),2)
        self.assertEqual(m11.getNumberOfCols(),3)
        self.assertTrue(m11.getData().isEqual(DataArrayDouble([6,9,12,15,3,18]),1e-12))
        m11=m11-m1
        self.assertEqual(m11.getNumberOfRows(),2)
        self.assertEqual(m11.getNumberOfCols(),3)
        self.assertTrue(m11.getData().isEqual(DataArrayDouble([4,6,8,10,2,12]),1e-12))
        m11-=m1
        self.assertEqual(m1.getNumberOfRows(),2)
        self.assertEqual(m1.getNumberOfCols(),3)
        self.assertTrue(m1.getData().isEqual(DataArrayDouble([2,3,4,5,1,6]),1e-12))
        m1.transpose()
        self.assertEqual(m1.getNumberOfRows(),3)
        self.assertEqual(m1.getNumberOfCols(),2)
        self.assertTrue(m1.getData().isEqual(DataArrayDouble([2,5,3,1,4,6]),1e-12))
        #m1np=m0np.transpose()
        m2=m0*m1
        self.assertEqual(m2.getNumberOfRows(),2)
        self.assertEqual(m2.getNumberOfCols(),2)
        self.assertTrue(m2.getData().isEqual(DataArrayDouble([29,37,37,62]),1e-12))
        pass

    def testSwig2AMR2(self):
        """ Test condensation of fine IMesh instance into a coarse one, with a factor. See testRemapperAMR1 in MEDCouplingRemapperTest.py file to see how the expected value is obtained."""
        coarse=DataArrayDouble(35) ; coarse.iota(0) #X=5,Y=7
        fine=DataArrayDouble(3*2*4*4) ; fine.iota(0) #X=3,Y=2 refined by 4
        MEDCouplingIMesh.CondenseFineToCoarse([5,7],fine,[(1,4),(2,4)],[4,4],coarse)
        self.assertTrue(coarse.isEqual(DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,312,376,440,14,15,1080,1144,1208,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34]),1e-12))
        # 3D
        coarse=DataArrayDouble(175) ; coarse.iota(0) #X=5,Y=7,Z=5
        fine=DataArrayDouble(3*2*3*4*4*4) ; fine.iota(0) #X=3,Y=2,Z=3 refined by 4
        MEDCouplingIMesh.CondenseFineToCoarse([5,7,5],fine,[(1,4),(2,4),(1,4)],[4,4,4],coarse)
        self.assertTrue(coarse.isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,35.,36.,37.,38.,39.,40.,41.,42.,43.,44.,45.,10464.,10720.,10976.,49.,50.,13536.,13792.,14048.,54.,55.,56.,57.,58.,59.,60.,61.,62.,63.,64.,65.,66.,67.,68.,69.,70.,71.,72.,73.,74.,75.,76.,77.,78.,79.,80.,35040.,35296.,35552.,84.,85.,38112.,38368.,38624.,89.,90.,91.,92.,93.,94.,95.,96.,97.,98.,99.,100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,110.,111.,112.,113.,114.,115.,59616.,59872.,60128.,119.,120.,62688.,62944.,63200.,124.,125.,126.,127.,128.,129.,130.,131.,132.,133.,134.,135.,136.,137.,138.,139.,140.,141.,142.,143.,144.,145.,146.,147.,148.,149.,150.,151.,152.,153.,154.,155.,156.,157.,158.,159.,160.,161.,162.,163.,164.,165.,166.,167.,168.,169.,170.,171.,172.,173.,174.]),1e-12))
        # 1D
        coarse=DataArrayDouble(5) ; coarse.iota(0) #X=5
        fine=DataArrayDouble(3*4) ; fine.iota(0) #X=3 refined by 4
        MEDCouplingIMesh.CondenseFineToCoarse([5],fine,[(1,4)],[4],coarse)
        self.assertTrue(coarse.isEqual(DataArrayDouble([0,6,22,38,4]),1e-12))
        pass

    def testSwig2AMR3(self):
        """ Test spread of coarse IMesh instance into a fine one, with a factor."""
        coarse=DataArrayDouble(35) ; coarse.iota(0) #X=5,Y=7
        fine=DataArrayDouble(3*2*4*4) ; fine.iota(0) #X=3,Y=2 refined by 4
        MEDCouplingIMesh.SpreadCoarseToFine(coarse,[5,7],fine,[(1,4),(2,4)],[4,4])
        self.assertTrue(fine.isEqual(DataArrayDouble([11.,11.,11.,11.,12.,12.,12.,12.,13.,13.,13.,13.,11.,11.,11.,11.,12.,12.,12.,12.,13.,13.,13.,13.,11.,11.,11.,11.,12.,12.,12.,12.,13.,13.,13.,13.,11.,11.,11.,11.,12.,12.,12.,12.,13.,13.,13.,13.,16.,16.,16.,16.,17.,17.,17.,17.,18.,18.,18.,18.,16.,16.,16.,16.,17.,17.,17.,17.,18.,18.,18.,18.,16.,16.,16.,16.,17.,17.,17.,17.,18.,18.,18.,18.,16.,16.,16.,16.,17.,17.,17.,17.,18.,18.,18.,18.]),1e-12))
        # 3D
        coarse=DataArrayDouble(175) ; coarse.iota(0) #X=5,Y=7,Z=5
        fine=DataArrayDouble(3*2*3*4*4*4) ; fine.iota(0) #X=3,Y=2,Z=3 refined by 4
        MEDCouplingIMesh.SpreadCoarseToFine(coarse,[5,7,5],fine,[(1,4),(2,4),(1,4)],[4,4,4])
        self.assertTrue(fine.isEqual(DataArrayDouble([46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,46.,46.,46.,46.,47.,47.,47.,47.,48.,48.,48.,48.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,51.,51.,51.,51.,52.,52.,52.,52.,53.,53.,53.,53.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,81.,81.,81.,81.,82.,82.,82.,82.,83.,83.,83.,83.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,86.,86.,86.,86.,87.,87.,87.,87.,88.,88.,88.,88.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,116.,116.,116.,116.,117.,117.,117.,117.,118.,118.,118.,118.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.,121.,121.,121.,121.,122.,122.,122.,122.,123.,123.,123.,123.]),1e-12))
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(MEDCouplingIMesh("",3,DataArrayInt([6,8,6]),[0.,0.,0.],DataArrayDouble((1.,1.,1.)))) ; f.setArray(coarse) ; f.setName("tutu") ; f.checkConsistencyLight()
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(MEDCouplingIMesh("",3,DataArrayInt([13,9,13]),[1.,2.,1.],DataArrayDouble((0.25,0.25,0.25)))) ; f.setArray(fine) ; f.setName("tutu") ; f.checkConsistencyLight()
        # 1D
        coarse=DataArrayDouble(5) ; coarse.iota(0) #X=5
        fine=DataArrayDouble(3*4) ; fine.iota(0) #X=3 refined by 4
        MEDCouplingIMesh.SpreadCoarseToFine(coarse,[5],fine,[(1,4)],[4])
        self.assertTrue(fine.isEqual(DataArrayDouble([1.,1.,1.,1.,2.,2.,2.,2.,3.,3.,3.,3.]),1e-12))
        pass

    def testSwig2AMR4(self):
        """This test focuses on MEDCouplingCartesianAMRMesh.createPatchesFromCriterion method. To test it a field containing 0 everywhere except in the annulus (centered on the center of the mesh) value is 1."""
        im=MEDCouplingIMesh("mesh",2,[51,51],[0.,0.],[0.04,0.04])
        b=im.computeCellCenterOfMass() ; b-=[1.,1.] ; b=b.magnitude()
        ids=b.findIdsInRange(0.4,0.7)
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(im) ; f.setName("toto") ; arr=DataArrayDouble(im.getNumberOfCells()) ; arr[:]=0. ; arr[ids]=1. ; f.setArray(arr)
        # f.write("test.vti")
        amr=MEDCouplingCartesianAMRMesh(MEDCouplingIMesh("mesh",2,[51,51],[0.,0.],[0.04,0.04]))
        arr2=DataArrayByte(im.getNumberOfCells()) ; arr2[:]=0 ; arr2[ids]=1
        bso=BoxSplittingOptions() ; bso.setEfficiencyGoal(0.5); bso.setEfficiencyThreshold(0.8) ; bso.setMaximumNbOfCellsInPatch(3000) ; bso.setMinimumPatchLength(6) ; bso.setMaximumPatchLength(11)
        amr.createPatchesFromCriterion(bso,arr2,[2,2])
        m=amr.getImageMesh() ; m=m.buildUnstructured() ; m.changeSpaceDimension(3,1.)
        self.assertEqual(12,amr.getNumberOfPatches())
        exp0=[[(9,19),(9,19)],[(9,19),(31,41)],[(31,41),(9,19)],[(8,17),(19,25)],[(8,17),(25,31)],[(19,25),(8,17)],[(25,31),(8,17)],[(19,25),(33,42)],[(25,31),(33,42)],[(31,41),(31,41)],[(33,42),(19,25)],[(33,42),(25,31)]]
        for i,bltr in enumerate(exp0):
            self.assertEqual(amr[i].getBLTRRange(),bltr)
            pass
        self.assertAlmostEqual(0.666666666667,amr[3].getMesh().getImageMesh().computeSquareness(),12)
        #
        self.assertEqual(MEDCouplingStructuredMesh.ChangeReferenceToGlobalOfCompactFrmt([(8,32),(4,17)],[(0,24),(2,12)]),[(8,32),(6,16)])
        self.assertEqual(MEDCouplingStructuredMesh.ChangeReferenceFromGlobalOfCompactFrmt([(8,32),(4,17)],[(8,32),(6,16)]),[(0,24),(2,12)])
        self.assertTrue(amr.getImageMesh().isEqual(im,1e-12))
        m=amr.getImageMesh().asSingleCell().build1SGTUnstructured()
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([1,0,2,3])))
        self.assertTrue(m.getCoords().isEqualWithoutConsideringStr(DataArrayDouble([(0,0),(2,0),(0,2),(2,2)]),1e-12))
        pass

    def testSwig2AMR5(self):
        """ Idem testAMR3, test spread of coarse IMesh instance into a fine one, with a factor, but here ghost is used !"""
        # 1D
        coarse=DataArrayDouble(5+2) ; coarse.iota(-1) #X=5 with ghostLev=1
        fine=DataArrayDouble(3*4+2) ; fine.iota(1000) #X=3 refined by 4 with ghostLev=1
        MEDCouplingIMesh.SpreadCoarseToFineGhost(coarse,[5],fine,[(1,4)],[4],1)
        self.assertTrue(fine.isEqual(DataArrayDouble([0,1,1,1,1,2,2,2,2,3,3,3,3,4]),1e-12))
        coarse.iota(-1000)
        MEDCouplingIMesh.CondenseFineToCoarseGhost([5],fine,[(1,4)],[4],coarse,1)
        self.assertTrue(coarse.isEqual(DataArrayDouble([-1000.,-999.,4.,8.,12.,-995.,-994.]),1e-12))
        # 2D
        coarse=DataArrayDouble((5+2*1)*(7+2*1)) ; coarse.iota(0) #X=5,Y=7 with ghostLev=1
        fine=DataArrayDouble((3*4+2*1)*(2*4+2*1)) ; fine.iota(1000) #X=3,Y=2 refined by 4
        MEDCouplingIMesh.SpreadCoarseToFineGhost(coarse,[5,7],fine,[(1,4),(2,4)],[4,4],1)
        self.assertTrue(fine.isEqual(DataArrayDouble([15.,16.,16.,16.,16.,17.,17.,17.,17.,18.,18.,18.,18.,19.,22.,23.,23.,23.,23.,24.,24.,24.,24.,25.,25.,25.,25.,26.,22.,23.,23.,23.,23.,24.,24.,24.,24.,25.,25.,25.,25.,26.,22.,23.,23.,23.,23.,24.,24.,24.,24.,25.,25.,25.,25.,26.,22.,23.,23.,23.,23.,24.,24.,24.,24.,25.,25.,25.,25.,26.,29.,30.,30.,30.,30.,31.,31.,31.,31.,32.,32.,32.,32.,33.,29.,30.,30.,30.,30.,31.,31.,31.,31.,32.,32.,32.,32.,33.,29.,30.,30.,30.,30.,31.,31.,31.,31.,32.,32.,32.,32.,33.,29.,30.,30.,30.,30.,31.,31.,31.,31.,32.,32.,32.,32.,33.,36.,37.,37.,37.,37.,38.,38.,38.,38.,39.,39.,39.,39.,40.]),1e-12))
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(MEDCouplingIMesh("",2,DataArrayInt([8,10]),[0.,0.],DataArrayDouble((1.,1.)))) ; f.setArray(coarse) ; f.setName("tutu") ; f.checkConsistencyLight()
        coarse.iota(-1000)
        fine2=DataArrayDouble.Meld(fine,3*fine) ; coarse2=DataArrayDouble.Meld(coarse,3*coarse)
        MEDCouplingIMesh.CondenseFineToCoarseGhost([5,7],fine,[(1,4),(2,4)],[4,4],coarse,1)
        MEDCouplingIMesh.CondenseFineToCoarseGhost([5,7],fine2,[(1,4),(2,4)],[4,4],coarse2,1)
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(MEDCouplingIMesh("",2,DataArrayInt([8,10]),[0.,0.],DataArrayDouble((1.,1.)))) ; f.setArray(coarse) ; f.setName("tutu") ; f.checkConsistencyLight()
        coarseExp=DataArrayDouble([-1000.,-999.,-998.,-997.,-996.,-995.,-994.,-993.,-992.,-991.,-990.,-989.,-988.,-987.,-986.,-985.,-984.,-983.,-982.,-981.,-980.,-979.,-978.,368.,384.,400.,-974.,-973.,-972.,-971.,480.,496.,512.,-967.,-966.,-965.,-964.,-963.,-962.,-961.,-960.,-959.,-958.,-957.,-956.,-955.,-954.,-953.,-952.,-951.,-950.,-949.,-948.,-947.,-946.,-945.,-944.,-943.,-942.,-941.,-940.,-939.,-938.])
        self.assertTrue(coarse.isEqual(coarseExp,1e-12))
        self.assertTrue(coarse2[:,0].isEqual(coarseExp,1e-12))
        self.assertTrue(coarse2[:,1].isEqual(3*coarseExp,1e-12))
        pass

    def testSwig2AMR6(self):
        """ Idem testSwig2AMR5, except that only 2D is considered here, and fine to fine is considered here. At the end of the test some checks about typing with AMR structs."""
        amr=MEDCouplingCartesianAMRMesh("",2,[6,6],[0,0],[1,1])
        da=DataArrayDouble((5+2)*(5+2)) ; da.iota() ; da+=0.9
        amr.addPatch([(1,4),(2,4)],[4,4])
        amr.addPatch([(0,1),(0,1)],[4,4])
        amr.addPatch([(4,5),(3,4)],[4,4])
        amr.addPatch([(4,5),(1,3)],[4,4])
        amr.addPatch([(0,1),(1,4)],[4,4])
        da0=DataArrayDouble((3*4+2)*(2*4+2)) ; da0.iota() ; da0[:]+=0.2
        da1=DataArrayDouble((1*4+2)*(1*4+2)) ; da1.iota() ; da1[:]+=0.4
        da2=DataArrayDouble((1*4+2)*(1*4+2)) ; da2.iota() ; da2[:]+=0.6
        da3=DataArrayDouble((1*4+2)*(2*4+2)) ; da3.iota() ; da3[:]+=0.7
        da4=DataArrayDouble((1*4+2)*(3*4+2)) ; da4.iota() ; da4[:]+=0.8
        self.assertEqual(5,amr.getNumberOfPatches())
        l=[da0,da1,da2,da3,da4]
        lCpy=[elt.deepCopy() for elt in l]
        l2=[DataArrayDouble.Meld(elt,3*elt) for elt in l]
        amr.fillCellFieldOnPatchGhostAdv(0,da,1,l,False)
        amr.fillCellFieldOnPatchGhostAdv(0,DataArrayDouble.Meld(da,3*da),1,l2,False)
        amr.fillCellFieldOnPatchOnlyOnGhostZone(0,da,lCpy[0],1)
        #
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(amr.getImageMesh().buildWithGhost(1)) ; f.setArray(da) ; f.setName("all")
        f0=MEDCouplingFieldDouble(ON_CELLS) ; f0.setMesh(amr[0].getMesh().getImageMesh().buildWithGhost(1)) ; f0.setArray(da0) ; f0.setName("p0") ; f0.checkConsistencyLight()
        f1=MEDCouplingFieldDouble(ON_CELLS) ; f1.setMesh(amr[1].getMesh().getImageMesh().buildWithGhost(1)) ; f1.setArray(da1) ; f1.setName("p1") ; f1.checkConsistencyLight()
        f2=MEDCouplingFieldDouble(ON_CELLS) ; f2.setMesh(amr[2].getMesh().getImageMesh().buildWithGhost(1)) ; f2.setArray(da2) ; f2.setName("p2") ; f2.checkConsistencyLight()
        f3=MEDCouplingFieldDouble(ON_CELLS) ; f3.setMesh(amr[3].getMesh().getImageMesh().buildWithGhost(1)) ; f3.setArray(da3) ; f3.setName("p3") ; f3.checkConsistencyLight()
        f4=MEDCouplingFieldDouble(ON_CELLS) ; f4.setMesh(amr[4].getMesh().getImageMesh().buildWithGhost(1)) ; f4.setArray(da4) ; f4.setName("p4") ; f4.checkConsistencyLight()
        #
        da0Exp=DataArrayDouble([28.8,16.9,16.9,16.9,16.9,17.9,17.9,17.9,17.9,18.9,18.9,18.9,18.9,25.7,34.8,23.9,23.9,23.9,23.9,24.9,24.9,24.9,24.9,25.9,25.9,25.9,25.9,31.7,40.8,23.9,23.9,23.9,23.9,24.9,24.9,24.9,24.9,25.9,25.9,25.9,25.9,37.7,46.8,23.9,23.9,23.9,23.9,24.9,24.9,24.9,24.9,25.9,25.9,25.9,25.9,43.7,52.8,23.9,23.9,23.9,23.9,24.9,24.9,24.9,24.9,25.9,25.9,25.9,25.9,49.7,58.8,30.9,30.9,30.9,30.9,31.9,31.9,31.9,31.9,32.9,32.9,32.9,32.9,7.6,64.8,30.9,30.9,30.9,30.9,31.9,31.9,31.9,31.9,32.9,32.9,32.9,32.9,13.6,70.8,30.9,30.9,30.9,30.9,31.9,31.9,31.9,31.9,32.9,32.9,32.9,32.9,19.6,76.8,30.9,30.9,30.9,30.9,31.9,31.9,31.9,31.9,32.9,32.9,32.9,32.9,25.6,36.9,37.9,37.9,37.9,37.9,38.9,38.9,38.9,38.9,39.9,39.9,39.9,39.9,40.9])
        da0Exp2=DataArrayDouble([15.9,16.9,16.9,16.9,16.9,17.9,17.9,17.9,17.9,18.9,18.9,18.9,18.9,19.9,22.9,15.2,16.2,17.2,18.2,19.2,20.2,21.2,22.2,23.2,24.2,25.2,26.2,26.9,22.9,29.2,30.2,31.2,32.2,33.2,34.2,35.2,36.2,37.2,38.2,39.2,40.2,26.9,22.9,43.2,44.2,45.2,46.2,47.2,48.2,49.2,50.2,51.2,52.2,53.2,54.2,26.9,22.9,57.2,58.2,59.2,60.2,61.2,62.2,63.2,64.2,65.2,66.2,67.2,68.2,26.9,29.9,71.2,72.2,73.2,74.2,75.2,76.2,77.2,78.2,79.2,80.2,81.2,82.2,33.9,29.9,85.2,86.2,87.2,88.2,89.2,90.2,91.2,92.2,93.2,94.2,95.2,96.2,33.9,29.9,99.2,100.2,101.2,102.2,103.2,104.2,105.2,106.2,107.2,108.2,109.2,110.2,33.9,29.9,113.2,114.2,115.2,116.2,117.2,118.2,119.2,120.2,121.2,122.2,123.2,124.2,33.9,36.9,37.9,37.9,37.9,37.9,38.9,38.9,38.9,38.9,39.9,39.9,39.9,39.9,40.9])
        self.assertTrue(da0.isEqual(da0Exp,1e-12))
        self.assertTrue(l2[0][:,0].isEqual(da0Exp,1e-12))
        self.assertTrue(l2[0][:,1].isEqual(3*da0Exp,1e-12))
        self.assertTrue(lCpy[0].isEqual(da0Exp2,1e-12))
        #
        g0=amr.retrieveGridsAt(0)
        self.assertEqual(1,len(g0))
        self.assertTrue(isinstance(g0[0],MEDCouplingCartesianAMRPatchGF))
        g1=amr.retrieveGridsAt(1)
        self.assertEqual(5,len(g1))
        for i in range(5):
            self.assertTrue(isinstance(g1[i],MEDCouplingCartesianAMRPatch))
            pass
        pass

    def testSwig2AMR7(self):
        """Idem testSwig2AMR6 except that we are in 1D"""
        amr=MEDCouplingCartesianAMRMesh("",1,[6],[0],[1])
        da=DataArrayDouble(5+2) ; da.iota() ; da+=0.9
        amr.addPatch([(1,4)],[4])
        amr.addPatch([(0,1)],[4])
        da0=DataArrayDouble(3*4+2) ; da0.iota() ; da0[:]+=0.2
        da1=DataArrayDouble(1*4+2) ; da1.iota() ; da1[:]+=0.4
        self.assertEqual(2,amr.getNumberOfPatches())
        l=[da0,da1]
        lCpy=[elt.deepCopy() for elt in l]
        l2=[DataArrayDouble.Meld(elt,3*elt) for elt in l]
        amr.fillCellFieldOnPatchGhostAdv(0,da,1,l,False)
        amr.fillCellFieldOnPatchGhostAdv(0,DataArrayDouble.Meld(da,3*da),1,l2,False)
        amr.fillCellFieldOnPatchOnlyOnGhostZone(0,da,lCpy[0],1)
        #
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(amr.getImageMesh().buildWithGhost(1)) ; f.setArray(da) ; f.setName("all")
        f0=MEDCouplingFieldDouble(ON_CELLS) ; f0.setMesh(amr[0].getMesh().getImageMesh().buildWithGhost(1)) ; f0.setArray(da0) ; f0.setName("p0") ; f0.checkConsistencyLight()
        f1=MEDCouplingFieldDouble(ON_CELLS) ; f1.setMesh(amr[1].getMesh().getImageMesh().buildWithGhost(1)) ; f1.setArray(da1) ; f1.setName("p1") ; f1.checkConsistencyLight()
        #
        da0Exp=DataArrayDouble([4.4,2.9,2.9,2.9,2.9,3.9,3.9,3.9,3.9,4.9,4.9,4.9,4.9,5.9])
        da0Exp2=DataArrayDouble([1.9,1.2,2.2,3.2,4.2,5.2,6.2,7.2,8.2,9.2,10.2,11.2,12.2,5.9])
        self.assertTrue(da0.isEqual(da0Exp,1e-12))
        self.assertTrue(l2[0][:,0].isEqual(da0Exp,1e-12))
        self.assertTrue(l2[0][:,1].isEqual(3*da0Exp,1e-12))
        self.assertTrue(lCpy[0].isEqual(da0Exp2,1e-12))
        pass

    def testSwig2AMR8(self):
        """This test checks 'basic' operations for ghost update."""
        ghostSz=1
        amr=MEDCouplingCartesianAMRMesh("",2,[6,7],[0,0],[1,1])
        amr.addPatch([(1,4),(2,4)],[4,4])
        amr.addPatch([(4,5),(3,5)],[4,4])
        amr.addPatch([(0,1),(4,6)],[4,4])
        amr[0].addPatch([(10,12),(5,8)],[2,2])
        amr[1].addPatch([(0,1),(0,5)],[2,2])
        amr[2].addPatch([(3,4),(0,3)],[2,2])
        m=amr.buildMeshFromPatchEnvelop()
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([1,0,2,3,5,4,6,7,9,8,10,11])))
        self.assertTrue(m.getCoords().isEqualWithoutConsideringStr(DataArrayDouble([1.,2.,4.,2.,1.,4.,4.,4.,4.,3.,5.,3.,4.,5.,5.,5.,0.,4.,1.,4.,0.,6.,1.,6.],12,2),1e-12))
        self.assertEqual(3,amr.getMaxNumberOfLevelsRelativeToThis())
        att=MEDCouplingAMRAttribute(amr,[("Field",["X"])],ghostSz)
        att.alloc()
        d=att.getFieldOn(amr,"Field")
        self.assertEqual(56,d.getNumberOfTuples())
        self.assertEqual(1,d.getNumberOfComponents())
        d.iota() ; d+=0.1
        d0=att.getFieldOn(amr[0].getMesh(),"Field")
        self.assertEqual(140,d0.getNumberOfTuples())
        self.assertEqual(1,d0.getNumberOfComponents())
        d0.iota() ; d0+=0.2
        d1=att.getFieldOn(amr[1].getMesh(),"Field")
        self.assertEqual(60,d1.getNumberOfTuples())
        self.assertEqual(1,d1.getNumberOfComponents())
        d1.iota() ; d1+=0.3
        d2=att.getFieldOn(amr[2].getMesh(),"Field")
        self.assertEqual(60,d2.getNumberOfTuples())
        self.assertEqual(1,d2.getNumberOfComponents())
        d2.iota() ; d2+=0.4
        d00=att.getFieldOn(amr[0][0].getMesh(),"Field")
        self.assertEqual(48,d00.getNumberOfTuples())
        self.assertEqual(1,d00.getNumberOfComponents())
        d00.iota() ; d00+=0.5
        d10=att.getFieldOn(amr[1][0].getMesh(),"Field")
        self.assertEqual(48,d10.getNumberOfTuples())
        self.assertEqual(1,d10.getNumberOfComponents())
        d10.iota() ; d10+=0.6
        d20=att.getFieldOn(amr[2][0].getMesh(),"Field")
        self.assertEqual(32,d20.getNumberOfTuples())
        self.assertEqual(1,d20.getNumberOfComponents())
        d20.iota() ; d20+=0.7
        f=att.buildCellFieldOnRecurseWithoutOverlapWithoutGhost(amr,"Field")
        arrExp=DataArrayDouble([8.1,9.1,10.1,11.1,12.1,15.1,16.1,17.1,18.1,19.1,22.1,26.1,29.1,37.1,38.1,39.1,44.1,45.1,46.1,47.1,15.2,16.2,17.2,18.2,19.2,20.2,21.2,22.2,23.2,24.2,25.2,26.2,29.2,30.2,31.2,32.2,33.2,34.2,35.2,36.2,37.2,38.2,39.2,40.2,43.2,44.2,45.2,46.2,47.2,48.2,49.2,50.2,51.2,52.2,53.2,54.2,57.2,58.2,59.2,60.2,61.2,62.2,63.2,64.2,65.2,66.2,67.2,68.2,71.2,72.2,73.2,74.2,75.2,76.2,77.2,78.2,79.2,80.2,81.2,82.2,85.2,86.2,87.2,88.2,89.2,90.2,91.2,92.2,93.2,94.2,99.2,100.2,101.2,102.2,103.2,104.2,105.2,106.2,107.2,108.2,113.2,114.2,115.2,116.2,117.2,118.2,119.2,120.2,121.2,122.2,7.5,8.5,9.5,10.5,13.5,14.5,15.5,16.5,19.5,20.5,21.5,22.5,25.5,26.5,27.5,28.5,31.5,32.5,33.5,34.5,37.5,38.5,39.5,40.5,8.3,9.3,10.3,14.3,15.3,16.3,20.3,21.3,22.3,26.3,27.3,28.3,32.3,33.3,34.3,37.3,38.3,39.3,40.3,43.3,44.3,45.3,46.3,49.3,50.3,51.3,52.3,5.6,6.6,9.6,10.6,13.6,14.6,17.6,18.6,21.6,22.6,25.6,26.6,29.6,30.6,33.6,34.6,37.6,38.6,41.6,42.6,7.4,8.4,9.4,13.4,14.4,15.4,19.4,20.4,21.4,25.4,26.4,27.4,28.4,31.4,32.4,33.4,34.4,37.4,38.4,39.4,40.4,43.4,44.4,45.4,46.4,49.4,50.4,51.4,52.4,5.7,6.7,9.7,10.7,13.7,14.7,17.7,18.7,21.7,22.7,25.7,26.7])
        arrExp.setName("Field") ; arrExp.setInfoOnComponents(["X"])
        self.assertTrue(f.getArray().isEqual(arrExp,1e-12))
        m=MEDCoupling1SGTUMesh(f.getMesh())
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,2,1,7,8,3,2,8,9,4,3,9,10,5,4,10,11,7,6,12,13,8,7,13,14,9,8,14,15,10,9,15,16,11,10,16,17,13,12,18,19,17,16,20,21,19,18,22,23,24,23,27,28,25,24,28,29,26,25,29,30,28,27,32,33,29,28,33,34,30,29,34,35,31,30,35,36,38,37,50,51,39,38,51,52,40,39,52,53,41,40,53,54,42,41,54,55,43,42,55,56,44,43,56,57,45,44,57,58,46,45,58,59,47,46,59,60,48,47,60,61,49,48,61,62,51,50,63,64,52,51,64,65,53,52,65,66,54,53,66,67,55,54,67,68,56,55,68,69,57,56,69,70,58,57,70,71,59,58,71,72,60,59,72,73,61,60,73,74,62,61,74,75,64,63,76,77,65,64,77,78,66,65,78,79,67,66,79,80,68,67,80,81,69,68,81,82,70,69,82,83,71,70,83,84,72,71,84,85,73,72,85,86,74,73,86,87,75,74,87,88,77,76,89,90,78,77,90,91,79,78,91,92,80,79,92,93,81,80,93,94,82,81,94,95,83,82,95,96,84,83,96,97,85,84,97,98,86,85,98,99,87,86,99,100,88,87,100,101,90,89,102,103,91,90,103,104,92,91,104,105,93,92,105,106,94,93,106,107,95,94,107,108,96,95,108,109,97,96,109,110,98,97,110,111,99,98,111,112,100,99,112,113,101,100,113,114,103,102,115,116,104,103,116,117,105,104,117,118,106,105,118,119,107,106,119,120,108,107,120,121,109,108,121,122,110,109,122,123,111,110,123,124,112,111,124,125,116,115,126,127,117,116,127,128,118,117,128,129,119,118,129,130,120,119,130,131,121,120,131,132,122,121,132,133,123,122,133,134,124,123,134,135,125,124,135,136,127,126,137,138,128,127,138,139,129,128,139,140,130,129,140,141,131,130,141,142,132,131,142,143,133,132,143,144,134,133,144,145,135,134,145,146,136,135,146,147,149,148,153,154,150,149,154,155,151,150,155,156,152,151,156,157,154,153,158,159,155,154,159,160,156,155,160,161,157,156,161,162,159,158,163,164,160,159,164,165,161,160,165,166,162,161,166,167,164,163,168,169,165,164,169,170,166,165,170,171,167,166,171,172,169,168,173,174,170,169,174,175,171,170,175,176,172,171,176,177,174,173,178,179,175,174,179,180,176,175,180,181,177,176,181,182,184,183,187,188,185,184,188,189,186,185,189,190,188,187,191,192,189,188,192,193,190,189,193,194,192,191,195,196,193,192,196,197,194,193,197,198,196,195,199,200,197,196,200,201,198,197,201,202,200,199,204,205,201,200,205,206,202,201,206,207,204,203,208,209,205,204,209,210,206,205,210,211,207,206,211,212,209,208,213,214,210,209,214,215,211,210,215,216,212,211,216,217,214,213,218,219,215,214,219,220,216,215,220,221,217,216,221,222,224,223,226,227,225,224,227,228,227,226,229,230,228,227,230,231,230,229,232,233,231,230,233,234,233,232,235,236,234,233,236,237,236,235,238,239,237,236,239,240,239,238,241,242,240,239,242,243,242,241,244,245,243,242,245,246,245,244,247,248,246,245,248,249,248,247,250,251,249,248,251,252,251,250,253,254,252,251,254,255,257,256,260,261,258,257,261,262,259,258,262,263,261,260,264,265,262,261,265,266,263,262,266,267,265,264,268,269,266,265,269,270,267,266,270,271,269,268,273,274,270,269,274,275,271,270,275,276,272,271,276,277,274,273,278,279,275,274,279,280,276,275,280,281,277,276,281,282,279,278,283,284,280,279,284,285,281,280,285,286,282,281,286,287,284,283,288,289,285,284,289,290,286,285,290,291,287,286,291,292,289,288,293,294,290,289,294,295,291,290,295,296,292,291,296,297,299,298,301,302,300,299,302,303,302,301,304,305,303,302,305,306,305,304,307,308,306,305,308,309,308,307,310,311,309,308,311,312,311,310,313,314,312,311,314,315,314,313,316,317,315,314,317,318])))
        self.assertTrue(m.getCoords().isEqualWithoutConsideringStr(DataArrayDouble([0.,0.,1.,0.,2.,0.,3.,0.,4.,0.,5.,0.,0.,1.,1.,1.,2.,1.,3.,1.,4.,1.,5.,1.,0.,2.,1.,2.,2.,2.,3.,2.,4.,2.,5.,2.,0.,3.,1.,3.,4.,3.,5.,3.,0.,4.,1.,4.,2.,4.,3.,4.,4.,4.,1.,5.,2.,5.,3.,5.,4.,5.,5.,5.,1.,6.,2.,6.,3.,6.,4.,6.,5.,6.,1.,2.,1.25,2.,1.5,2.,1.75,2.,2.,2.,2.25,2.,2.5,2.,2.75,2.,3.,2.,3.25,2.,3.5,2.,3.75,2.,4.,2.,1.,2.25,1.25,2.25,1.5,2.25,1.75,2.25,2.,2.25,2.25,2.25,2.5,2.25,2.75,2.25,3.,2.25,3.25,2.25,3.5,2.25,3.75,2.25,4.,2.25,1.,2.5,1.25,2.5,1.5,2.5,1.75,2.5,2.,2.5,2.25,2.5,2.5,2.5,2.75,2.5,3.,2.5,3.25,2.5,3.5,2.5,3.75,2.5,4.,2.5,1.,2.75,1.25,2.75,1.5,2.75,1.75,2.75,2.,2.75,2.25,2.75,2.5,2.75,2.75,2.75,3.,2.75,3.25,2.75,3.5,2.75,3.75,2.75,4.,2.75,1.,3.,1.25,3.,1.5,3.,1.75,3.,2.,3.,2.25,3.,2.5,3.,2.75,3.,3.,3.,3.25,3.,3.5,3.,3.75,3.,4.,3.,1.,3.25,1.25,3.25,1.5,3.25,1.75,3.25,2.,3.25,2.25,3.25,2.5,3.25,2.75,3.25,3.,3.25,3.25,3.25,3.5,3.25,3.75,3.25,4.,3.25,1.,3.5,1.25,3.5,1.5,3.5,1.75,3.5,2.,3.5,2.25,3.5,2.5,3.5,2.75,3.5,3.,3.5,3.25,3.5,3.5,3.5,1.,3.75,1.25,3.75,1.5,3.75,1.75,3.75,2.,3.75,2.25,3.75,2.5,3.75,2.75,3.75,3.,3.75,3.25,3.75,3.5,3.75,1.,4.,1.25,4.,1.5,4.,1.75,4.,2.,4.,2.25,4.,2.5,4.,2.75,4.,3.,4.,3.25,4.,3.5,4.,3.5,3.25,3.625,3.25,3.75,3.25,3.875,3.25,4.,3.25,3.5,3.375,3.625,3.375,3.75,3.375,3.875,3.375,4.,3.375,3.5,3.5,3.625,3.5,3.75,3.5,3.875,3.5,4.,3.5,3.5,3.625,3.625,3.625,3.75,3.625,3.875,3.625,4.,3.625,3.5,3.75,3.625,3.75,3.75,3.75,3.875,3.75,4.,3.75,3.5,3.875,3.625,3.875,3.75,3.875,3.875,3.875,4.,3.875,3.5,4.,3.625,4.,3.75,4.,3.875,4.,4.,4.,4.25,3.,4.5,3.,4.75,3.,5.,3.,4.25,3.25,4.5,3.25,4.75,3.25,5.,3.25,4.25,3.5,4.5,3.5,4.75,3.5,5.,3.5,4.25,3.75,4.5,3.75,4.75,3.75,5.,3.75,4.25,4.,4.5,4.,4.75,4.,5.,4.,4.,4.25,4.25,4.25,4.5,4.25,4.75,4.25,5.,4.25,4.,4.5,4.25,4.5,4.5,4.5,4.75,4.5,5.,4.5,4.,4.75,4.25,4.75,4.5,4.75,4.75,4.75,5.,4.75,4.,5.,4.25,5.,4.5,5.,4.75,5.,5.,5.,4.,3.,4.125,3.,4.25,3.,4.,3.125,4.125,3.125,4.25,3.125,4.,3.25,4.125,3.25,4.25,3.25,4.,3.375,4.125,3.375,4.25,3.375,4.,3.5,4.125,3.5,4.25,3.5,4.,3.625,4.125,3.625,4.25,3.625,4.,3.75,4.125,3.75,4.25,3.75,4.,3.875,4.125,3.875,4.25,3.875,4.,4.,4.125,4.,4.25,4.,4.,4.125,4.125,4.125,4.25,4.125,4.,4.25,4.125,4.25,4.25,4.25,0.,4.,0.25,4.,0.5,4.,0.75,4.,0.,4.25,0.25,4.25,0.5,4.25,0.75,4.25,0.,4.5,0.25,4.5,0.5,4.5,0.75,4.5,0.,4.75,0.25,4.75,0.5,4.75,0.75,4.75,1.,4.75,0.,5.,0.25,5.,0.5,5.,0.75,5.,1.,5.,0.,5.25,0.25,5.25,0.5,5.25,0.75,5.25,1.,5.25,0.,5.5,0.25,5.5,0.5,5.5,0.75,5.5,1.,5.5,0.,5.75,0.25,5.75,0.5,5.75,0.75,5.75,1.,5.75,0.,6.,0.25,6.,0.5,6.,0.75,6.,1.,6.,0.75,4.,0.875,4.,1.,4.,0.75,4.125,0.875,4.125,1.,4.125,0.75,4.25,0.875,4.25,1.,4.25,0.75,4.375,0.875,4.375,1.,4.375,0.75,4.5,0.875,4.5,1.,4.5,0.75,4.625,0.875,4.625,1.,4.625,0.75,4.75,0.875,4.75,1.,4.75],319,2),1e-12))
        # the test is here ! To be called after iteration with no remesh
        att.synchronizeAllGhostZones()
        f=att.buildCellFieldOnWithGhost(amr,"Field") ; f.checkConsistencyLight()
        ftmp=att.buildCellFieldOnWithoutGhost(amr,"Field") ; ftmp.checkConsistencyLight() ; self.assertTrue(ftmp.getArray().isEqualWithoutConsideringStr(DataArrayDouble([8.1,9.1,10.1,11.1,12.1,15.1,16.1,17.1,18.1,19.1,22.1,23.1,24.1,25.1,26.1,29.1,30.1,31.1,32.1,33.1,36.1,37.1,38.1,39.1,40.1,43.1,44.1,45.1,46.1,47.1]),1e-12))
        f0=att.buildCellFieldOnWithGhost(amr[0].getMesh(),"Field")
        f1=att.buildCellFieldOnWithGhost(amr[1].getMesh(),"Field")
        f2=att.buildCellFieldOnWithGhost(amr[2].getMesh(),"Field")
        f00=att.buildCellFieldOnWithGhost(amr[0][0].getMesh(),"Field")
        f10=att.buildCellFieldOnWithGhost(amr[1][0].getMesh(),"Field")
        f20=att.buildCellFieldOnWithGhost(amr[2][0].getMesh(),"Field")
        self.assertTrue(f.getArray().isEqualWithoutConsideringStr(DataArrayDouble([0.1,1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1,9.1,10.1,11.1,12.1,13.1,14.1,15.1,16.1,17.1,18.1,19.1,20.1,21.1,22.1,23.1,24.1,25.1,26.1,27.1,28.1,29.1,30.1,31.1,32.1,33.1,34.1,35.1,36.1,37.1,38.1,39.1,40.1,41.1,42.1,43.1,44.1,45.1,46.1,47.1,48.1,49.1,50.1,51.1,52.1,53.1,54.1,55.1]),1e-12))
        self.assertTrue(f0.getArray().isEqualWithoutConsideringStr(DataArrayDouble([15.1,16.1,16.1,16.1,16.1,17.1,17.1,17.1,17.1,18.1,18.1,18.1,18.1,19.1,22.1,15.2,16.2,17.2,18.2,19.2,20.2,21.2,22.2,23.2,24.2,25.2,26.2,26.1,22.1,29.2,30.2,31.2,32.2,33.2,34.2,35.2,36.2,37.2,38.2,39.2,40.2,26.1,22.1,43.2,44.2,45.2,46.2,47.2,48.2,49.2,50.2,51.2,52.2,53.2,54.2,26.1,22.1,57.2,58.2,59.2,60.2,61.2,62.2,63.2,64.2,65.2,66.2,67.2,68.2,26.1,29.1,71.2,72.2,73.2,74.2,75.2,76.2,77.2,78.2,79.2,80.2,81.2,82.2,7.3,29.1,85.2,86.2,87.2,88.2,89.2,90.2,91.2,92.2,93.2,94.2,95.2,96.2,13.3,29.1,99.2,100.2,101.2,102.2,103.2,104.2,105.2,106.2,107.2,108.2,109.2,110.2,19.3,29.1,113.2,114.2,115.2,116.2,117.2,118.2,119.2,120.2,121.2,122.2,123.2,124.2,25.3,10.4,37.1,37.1,37.1,37.1,38.1,38.1,38.1,38.1,39.1,39.1,39.1,39.1,31.3]),1e-12))
        self.assertTrue(f1.getArray().isEqualWithoutConsideringStr(DataArrayDouble([68.2,26.1,26.1,26.1,26.1,27.1,82.2,7.3,8.3,9.3,10.3,34.1,96.2,13.3,14.3,15.3,16.3,34.1,110.2,19.3,20.3,21.3,22.3,34.1,124.2,25.3,26.3,27.3,28.3,34.1,39.1,31.3,32.3,33.3,34.3,41.1,39.1,37.3,38.3,39.3,40.3,41.1,39.1,43.3,44.3,45.3,46.3,41.1,39.1,49.3,50.3,51.3,52.3,41.1,46.1,47.1,47.1,47.1,47.1,48.1]),1e-12))
        self.assertTrue(f2.getArray().isEqualWithoutConsideringStr(DataArrayDouble([28.1,29.1,29.1,29.1,29.1,113.2,35.1,7.4,8.4,9.4,10.4,37.1,35.1,13.4,14.4,15.4,16.4,37.1,35.1,19.4,20.4,21.4,22.4,37.1,35.1,25.4,26.4,27.4,28.4,37.1,42.1,31.4,32.4,33.4,34.4,44.1,42.1,37.4,38.4,39.4,40.4,44.1,42.1,43.4,44.4,45.4,46.4,44.1,42.1,49.4,50.4,51.4,52.4,44.1,49.1,50.1,50.1,50.1,50.1,51.1]),1e-12))
        self.assertTrue(f00.getArray().isEqualWithoutConsideringStr(DataArrayDouble([80.2,81.2,81.2,82.2,82.2,9.6,94.2,7.5,8.5,9.5,10.5,13.6,94.2,13.5,14.5,15.5,16.5,17.6,108.2,19.5,20.5,21.5,22.5,21.6,108.2,25.5,26.5,27.5,28.5,25.6,122.2,31.5,32.5,33.5,34.5,29.6,122.2,37.5,38.5,39.5,40.5,33.6,39.1,39.1,39.1,39.1,39.1,37.6]),1e-12))
        self.assertTrue(f10.getArray().isEqualWithoutConsideringStr(DataArrayDouble([68.2,26.1,26.1,26.1,82.2,5.6,6.6,8.3,82.2,9.6,10.6,8.3,10.5,13.6,14.6,14.3,16.5,17.6,18.6,14.3,22.5,21.6,22.6,20.3,28.5,25.6,26.6,20.3,34.5,29.6,30.6,26.3,40.5,33.6,34.6,26.3,39.1,37.6,38.6,32.3,39.1,41.6,42.6,32.3,39.1,37.3,37.3,38.3]),1e-12))
        self.assertTrue(f20.getArray().isEqualWithoutConsideringStr(DataArrayDouble([29.1,29.1,29.1,113.2,9.4,5.7,6.7,37.1,9.4,9.7,10.7,37.1,15.4,13.7,14.7,37.1,15.4,17.7,18.7,37.1,21.4,21.7,22.7,37.1,21.4,25.7,26.7,37.1,27.4,28.4,28.4,37.1]),1e-12))
        pass

    def testSwig2AMR9(self):
        """ Equivalent to testSwig2AMR8 except that here the ghost level is 2 !"""
        ghostSz=2
        amr=MEDCouplingCartesianAMRMesh("",2,[6,7],[0,0],[1,1])
        amr.addPatch([(1,4),(2,4)],[4,4])
        amr.addPatch([(4,5),(3,5)],[4,4])
        amr.addPatch([(0,1),(4,6)],[4,4])
        amr[0].addPatch([(10,12),(5,8)],[2,2])
        amr[1].addPatch([(0,1),(0,5)],[2,2])
        amr[2].addPatch([(3,4),(0,3)],[2,2])
        self.assertEqual(3,amr.getMaxNumberOfLevelsRelativeToThis())
        att=MEDCouplingAMRAttribute(amr,[("Field",["X"])],ghostSz)
        att.alloc()
        d=att.getFieldOn(amr,"Field")
        self.assertEqual(90,d.getNumberOfTuples())
        self.assertEqual(1,d.getNumberOfComponents())
        d.iota() ; d+=0.1
        d0=att.getFieldOn(amr[0].getMesh(),"Field")
        self.assertEqual(192,d0.getNumberOfTuples())
        self.assertEqual(1,d0.getNumberOfComponents())
        d0.iota() ; d0+=0.2
        d1=att.getFieldOn(amr[1].getMesh(),"Field")
        self.assertEqual(96,d1.getNumberOfTuples())
        self.assertEqual(1,d1.getNumberOfComponents())
        d1.iota() ; d1+=0.3
        d2=att.getFieldOn(amr[2].getMesh(),"Field")
        self.assertEqual(96,d2.getNumberOfTuples())
        self.assertEqual(1,d2.getNumberOfComponents())
        d2.iota() ; d2+=0.4
        d00=att.getFieldOn(amr[0][0].getMesh(),"Field")
        self.assertEqual(80,d00.getNumberOfTuples())
        self.assertEqual(1,d00.getNumberOfComponents())
        d00.iota() ; d00+=0.5
        d10=att.getFieldOn(amr[1][0].getMesh(),"Field")
        self.assertEqual(84,d10.getNumberOfTuples())
        self.assertEqual(1,d10.getNumberOfComponents())
        d10.iota() ; d10+=0.6
        d20=att.getFieldOn(amr[2][0].getMesh(),"Field")
        self.assertEqual(60,d20.getNumberOfTuples())
        self.assertEqual(1,d20.getNumberOfComponents())
        d20.iota() ; d20+=0.7
        # the test is here ! To be called after iteration with no remesh
        att.synchronizeAllGhostZones()
        f=att.buildCellFieldOnWithGhost(amr,"Field")
        f0=att.buildCellFieldOnWithGhost(amr[0].getMesh(),"Field")
        f1=att.buildCellFieldOnWithGhost(amr[1].getMesh(),"Field")
        f2=att.buildCellFieldOnWithGhost(amr[2].getMesh(),"Field")
        f00=att.buildCellFieldOnWithGhost(amr[0][0].getMesh(),"Field")
        f10=att.buildCellFieldOnWithGhost(amr[1][0].getMesh(),"Field")
        f20=att.buildCellFieldOnWithGhost(amr[2][0].getMesh(),"Field")
        self.assertTrue(f0.getArray().isEqualWithoutConsideringStr(DataArrayDouble([29.1,29.1,30.1,30.1,30.1,30.1,31.1,31.1,31.1,31.1,32.1,32.1,32.1,32.1,33.1,33.1,29.1,29.1,30.1,30.1,30.1,30.1,31.1,31.1,31.1,31.1,32.1,32.1,32.1,32.1,33.1,33.1,38.1,38.1,34.2,35.2,36.2,37.2,38.2,39.2,40.2,41.2,42.2,43.2,44.2,45.2,42.1,42.1,38.1,38.1,50.2,51.2,52.2,53.2,54.2,55.2,56.2,57.2,58.2,59.2,60.2,61.2,42.1,42.1,38.1,38.1,66.2,67.2,68.2,69.2,70.2,71.2,72.2,73.2,74.2,75.2,76.2,77.2,42.1,42.1,38.1,38.1,82.2,83.2,84.2,85.2,86.2,87.2,88.2,89.2,90.2,91.2,92.2,93.2,42.1,42.1,47.1,47.1,98.2,99.2,100.2,101.2,102.2,103.2,104.2,105.2,106.2,107.2,108.2,109.2,18.3,19.3,47.1,47.1,114.2,115.2,116.2,117.2,118.2,119.2,120.2,121.2,122.2,123.2,124.2,125.2,26.3,27.3,47.1,47.1,130.2,131.2,132.2,133.2,134.2,135.2,136.2,137.2,138.2,139.2,140.2,141.2,34.3,35.3,47.1,47.1,146.2,147.2,148.2,149.2,150.2,151.2,152.2,153.2,154.2,155.2,156.2,157.2,42.3,43.3,20.4,21.4,57.1,57.1,57.1,57.1,58.1,58.1,58.1,58.1,59.1,59.1,59.1,59.1,50.3,51.3,28.4,29.4,57.1,57.1,57.1,57.1,58.1,58.1,58.1,58.1,59.1,59.1,59.1,59.1,58.3,59.3]),1e-12))
        self.assertTrue(f1.getArray().isEqualWithoutConsideringStr(DataArrayDouble([76.2,77.2,42.1,42.1,42.1,42.1,43.1,43.1,92.2,93.2,42.1,42.1,42.1,42.1,43.1,43.1,108.2,109.2,18.3,19.3,20.3,21.3,52.1,52.1,124.2,125.2,26.3,27.3,28.3,29.3,52.1,52.1,140.2,141.2,34.3,35.3,36.3,37.3,52.1,52.1,156.2,157.2,42.3,43.3,44.3,45.3,52.1,52.1,59.1,59.1,50.3,51.3,52.3,53.3,61.1,61.1,59.1,59.1,58.3,59.3,60.3,61.3,61.1,61.1,59.1,59.1,66.3,67.3,68.3,69.3,61.1,61.1,59.1,59.1,74.3,75.3,76.3,77.3,61.1,61.1,68.1,68.1,69.1,69.1,69.1,69.1,70.1,70.1,68.1,68.1,69.1,69.1,69.1,69.1,70.1,70.1]),1e-12))
        self.assertTrue(f2.getArray().isEqualWithoutConsideringStr(DataArrayDouble([46.1,46.1,47.1,47.1,47.1,47.1,130.2,131.2,46.1,46.1,47.1,47.1,47.1,47.1,146.2,147.2,55.1,55.1,18.4,19.4,20.4,21.4,57.1,57.1,55.1,55.1,26.4,27.4,28.4,29.4,57.1,57.1,55.1,55.1,34.4,35.4,36.4,37.4,57.1,57.1,55.1,55.1,42.4,43.4,44.4,45.4,57.1,57.1,64.1,64.1,50.4,51.4,52.4,53.4,66.1,66.1,64.1,64.1,58.4,59.4,60.4,61.4,66.1,66.1,64.1,64.1,66.4,67.4,68.4,69.4,66.1,66.1,64.1,64.1,74.4,75.4,76.4,77.4,66.1,66.1,73.1,73.1,74.1,74.1,74.1,74.1,75.1,75.1,73.1,73.1,74.1,74.1,74.1,74.1,75.1,75.1]),1e-12))
        self.assertTrue(f00.getArray().isEqualWithoutConsideringStr(DataArrayDouble([107.2,107.2,108.2,108.2,109.2,109.2,14.6,15.6,107.2,107.2,108.2,108.2,109.2,109.2,20.6,21.6,123.2,123.2,18.5,19.5,20.5,21.5,26.6,27.6,123.2,123.2,26.5,27.5,28.5,29.5,32.6,33.6,139.2,139.2,34.5,35.5,36.5,37.5,38.6,39.6,139.2,139.2,42.5,43.5,44.5,45.5,44.6,45.6,155.2,155.2,50.5,51.5,52.5,53.5,50.6,51.6,155.2,155.2,58.5,59.5,60.5,61.5,56.6,57.6,59.1,59.1,59.1,59.1,59.1,59.1,62.6,63.6,59.1,59.1,59.1,59.1,59.1,59.1,68.6,69.6]),1e-12))
        self.assertTrue(f10.getArray().isEqualWithoutConsideringStr(DataArrayDouble([93.2,93.2,42.1,42.1,42.1,42.1,93.2,93.2,42.1,42.1,42.1,42.1,109.2,109.2,14.6,15.6,19.3,19.3,109.2,109.2,20.6,21.6,19.3,19.3,20.5,21.5,26.6,27.6,27.3,27.3,28.5,29.5,32.6,33.6,27.3,27.3,36.5,37.5,38.6,39.6,35.3,35.3,44.5,45.5,44.6,45.6,35.3,35.3,52.5,53.5,50.6,51.6,43.3,43.3,60.5,61.5,56.6,57.6,43.3,43.3,59.1,59.1,62.6,63.6,51.3,51.3,59.1,59.1,68.6,69.6,51.3,51.3,59.1,59.1,58.3,58.3,59.3,59.3,59.1,59.1,58.3,58.3,59.3,59.3]),1e-12))
        self.assertTrue(f20.getArray().isEqualWithoutConsideringStr(DataArrayDouble([47.1,47.1,47.1,47.1,146.2,146.2,47.1,47.1,47.1,47.1,146.2,146.2,20.4,20.4,14.7,15.7,57.1,57.1,20.4,20.4,20.7,21.7,57.1,57.1,28.4,28.4,26.7,27.7,57.1,57.1,28.4,28.4,32.7,33.7,57.1,57.1,36.4,36.4,38.7,39.7,57.1,57.1,36.4,36.4,44.7,45.7,57.1,57.1,44.4,44.4,45.4,45.4,57.1,57.1,44.4,44.4,45.4,45.4,57.1,57.1]),1e-12))
        self.assertTrue(MEDCouplingStructuredMesh.ComputeCornersGhost([3],1).isEqual(DataArrayInt([0,4])))
        self.assertTrue(MEDCouplingStructuredMesh.ComputeCornersGhost([3],2).isEqual(DataArrayInt([0,1,5,6])))
        self.assertTrue(MEDCouplingStructuredMesh.ComputeCornersGhost([5,6],1).isEqual(DataArrayInt([0,6,49,55])))
        self.assertTrue(MEDCouplingStructuredMesh.ComputeCornersGhost([5,6],2).isEqual(DataArrayInt([0,8,10,16,73,79,81,89])))
        self.assertTrue(MEDCouplingStructuredMesh.ComputeCornersGhost([5,6,3],1).isEqual(DataArrayInt([0,6,49,55,224,230,273,279])))
        self.assertTrue(MEDCouplingStructuredMesh.ComputeCornersGhost([5,6,3],2).isEqual(DataArrayInt([0,8,81,89,100,106,163,169,460,466,523,529,540,548,621,629])))
        pass

    def testSwig2AMR10(self):
        """ This test, focuses on basic operations of coarse to fine and fine to coarse and ghost zone update with a ghost size set to 2 and dimension equal to 2."""
        szGhost=2
        amr=MEDCouplingCartesianAMRMesh("",2,[11,11],[0,0],[0.1,0.1])
        amr.addPatch([(3,8),(0,3)],[2,2])
        amr[0].addPatch([(0,10),(3,6)],[3,3])
        amr[0].addPatch([(2,6),(0,3)],[3,3])
        amr[0].addPatch([(6,10),(2,3)],[3,3])
        amr.addPatch([(3,8),(3,6)],[2,2])
        amr[1].addPatch([(0,4),(0,6)],[3,3])
        amr[1].addPatch([(7,10),(0,4)],[3,3])
        amr[1].addPatch([(4,7),(0,3)],[3,3])
        amr[1].addPatch([(4,7),(3,6)],[3,3])
        amr.addPatch([(0,3),(6,10)],[2,2])
        self.assertEqual(([(30,39),(27,36)],[6,6]),amr[1][3].getMesh().positionRelativeToGodFather())
        self.assertEqual(([(6,16),(6,12)],[2,2]),amr[1].getMesh().positionRelativeToGodFather())
        self.assertTrue(not MEDCouplingStructuredMesh.AreRangesIntersect([(30,39),(27,36)],[(6,16),(6,12)]))
        self.assertTrue(MEDCouplingStructuredMesh.AreRangesIntersect([(30,39),(27,36)],[(28,32),(35,37)]))
        da=DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.])
        MEDCouplingStructuredMesh.AssignPartOfFieldOfDoubleUsing([3,4],da,[(1,3),(2,3)],DataArrayDouble([7.7,8.8]))
        self.assertTrue(da.isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.7,8.8,9.,10.,11.]),1e-12))
        att=MEDCouplingAMRAttribute(amr,[("YY",1)],szGhost)
        att.spillNatures([IntensiveMaximum])
        att.alloc()
        yy=att.getFieldOn(amr,"YY") ; yy.iota(0.01)
        yy=att.getFieldOn(amr[0].getMesh(),"YY") ; yy.iota(0.02)
        yy=att.getFieldOn(amr[1].getMesh(),"YY") ; yy.iota(0.03)
        yy=att.getFieldOn(amr[0][0].getMesh(),"YY") ; yy.iota(0.04)
        yy=att.getFieldOn(amr[0][1].getMesh(),"YY") ; yy.iota(0.05)
        yy=att.getFieldOn(amr[0][2].getMesh(),"YY") ; yy.iota(0.06)
        yy=att.getFieldOn(amr[1][0].getMesh(),"YY") ; yy.iota(0.07)
        yy=att.getFieldOn(amr[1][1].getMesh(),"YY") ; yy.iota(0.08)
        yy=att.getFieldOn(amr[1][2].getMesh(),"YY") ; yy.iota(0.09)
        yy=att.getFieldOn(amr[1][3].getMesh(),"YY") ; yy.iota(0.10)
        yy=att.getFieldOn(amr[2].getMesh(),"YY") ; yy.iota(0.11)
        att2=att.deepCopy() ; att3=att2.deepCopy() ; att4=att3.deepCopy() ; att5=att4.deepCopy() ; att6=att5.deepCopy()
        ###
        att.synchronizeFineToCoarseBetween(2,1)
        ###
        for pos in [(),(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(1,3)]:
            self.assertTrue(att.getFieldOn(att.getMyGodFather().getMeshAtPosition(pos),"YY").isEqual(att2.getFieldOn(att2.getMyGodFather().getMeshAtPosition(pos),"YY"),1e-12))
            pass
        for pos in [(0,),(1,)]:
            self.assertTrue(not att.getFieldOn(att.getMyGodFather().getMeshAtPosition(pos),"YY").isEqual(att2.getFieldOn(att2.getMyGodFather().getMeshAtPosition(pos),"YY"),1e-12))
            pass
        self.assertTrue(att.getFieldOn(amr[0].getMesh(),"YY").isEqualWithoutConsideringStr(DataArrayDouble([0.02,1.02,2.02,3.02,4.02,5.02,6.02,7.02,8.02,9.02,10.02,11.02,12.02,13.02,14.02,15.02,16.02,17.02,18.02,19.02,20.02,21.02,22.02,23.02,24.02,25.02,26.02,27.02,28.02,29.02,30.02,31.02,51.05,54.05,57.05,60.05,36.02,37.02,38.02,39.02,40.02,41.02,42.02,43.02,44.02,45.02,99.05,102.05,105.05,108.05,50.02,51.02,52.02,53.02,54.02,55.02,56.02,57.02,58.02,59.02,147.05,150.05,153.05,156.05,51.06,54.06,57.06,60.06,68.02,69.02,70.02,71.02,105.04,108.04,111.04,114.04,117.04,120.04,123.04,126.04,129.04,132.04,82.02,83.02,84.02,85.02,207.04,210.04,213.04,216.04,219.04,222.04,225.04,228.04,231.04,234.04,96.02,97.02,98.02,99.02,309.04,312.04,315.04,318.04,321.04,324.04,327.04,330.04,333.04,336.04,110.02,111.02,112.02,113.02,114.02,115.02,116.02,117.02,118.02,119.02,120.02,121.02,122.02,123.02,124.02,125.02,126.02,127.02,128.02,129.02,130.02,131.02,132.02,133.02,134.02,135.02,136.02,137.02,138.02,139.02]),1e-12))
        self.assertTrue(att.getFieldOn(amr[1].getMesh(),"YY").isEqualWithoutConsideringStr(DataArrayDouble([0.03,1.03,2.03,3.03,4.03,5.03,6.03,7.03,8.03,9.03,10.03,11.03,12.03,13.03,14.03,15.03,16.03,17.03,18.03,19.03,20.03,21.03,22.03,23.03,24.03,25.03,26.03,27.03,28.03,29.03,51.07,54.07,57.07,60.07,42.09,45.09,48.09,42.08,45.08,48.08,40.03,41.03,42.03,43.03,99.07,102.07,105.07,108.07,81.09,84.09,87.09,81.08,84.08,87.08,54.03,55.03,56.03,57.03,147.07,150.07,153.07,156.07,120.09,123.09,126.09,120.08,123.08,126.08,68.03,69.03,70.03,71.03,195.07,198.07,201.07,204.07,42.1,45.1,48.1,159.08,162.08,165.08,82.03,83.03,84.03,85.03,243.07,246.07,249.07,252.07,81.1,84.1,87.1,93.03,94.03,95.03,96.03,97.03,98.03,99.03,291.07,294.07,297.07,300.07,120.1,123.1,126.1,107.03,108.03,109.03,110.03,111.03,112.03,113.03,114.03,115.03,116.03,117.03,118.03,119.03,120.03,121.03,122.03,123.03,124.03,125.03,126.03,127.03,128.03,129.03,130.03,131.03,132.03,133.03,134.03,135.03,136.03,137.03,138.03,139.03]),1e-12))
        del att
        ####
        att2.synchronizeAllGhostZonesOfDirectChidrenOf(att2.getMyGodFather())
        ### Only the 3 (0) (1) and (2) are modified (0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (1,3) are not modified.
        exp2=DataArrayDouble([0.11,1.11,2.11,3.11,4.11,5.11,6.11,7.11,86.03,87.03,10.11,11.11,12.11,13.11,14.11,15.11,16.11,17.11,100.03,101.03,20.11,21.11,22.11,23.11,24.11,25.11,26.11,27.11,28.11,29.11,30.11,31.11,32.11,33.11,34.11,35.11,36.11,37.11,38.11,39.11,40.11,41.11,42.11,43.11,44.11,45.11,46.11,47.11,48.11,49.11,50.11,51.11,52.11,53.11,54.11,55.11,56.11,57.11,58.11,59.11,60.11,61.11,62.11,63.11,64.11,65.11,66.11,67.11,68.11,69.11,70.11,71.11,72.11,73.11,74.11,75.11,76.11,77.11,78.11,79.11,80.11,81.11,82.11,83.11,84.11,85.11,86.11,87.11,88.11,89.11,90.11,91.11,92.11,93.11,94.11,95.11,96.11,97.11,98.11,99.11,100.11,101.11,102.11,103.11,104.11,105.11,106.11,107.11,108.11,109.11,110.11,111.11,112.11,113.11,114.11,115.11,116.11,117.11,118.11,119.11])
        self.assertTrue(att2.getFieldOn(att2.getMyGodFather().getMeshAtPosition((2,)),"YY").isEqualWithoutConsideringStr(exp2,1e-12))
        exp3=DataArrayDouble([0.03,1.03,86.02,87.02,88.02,89.02,90.02,91.02,92.02,93.02,94.02,95.02,12.03,13.03,14.03,15.03,100.02,101.02,102.02,103.02,104.02,105.02,106.02,107.02,108.02,109.02,26.03,27.03,28.03,29.03,30.03,31.03,32.03,33.03,34.03,35.03,36.03,37.03,38.03,39.03,40.03,41.03,42.03,43.03,44.03,45.03,46.03,47.03,48.03,49.03,50.03,51.03,52.03,53.03,54.03,55.03,56.03,57.03,58.03,59.03,60.03,61.03,62.03,63.03,64.03,65.03,66.03,67.03,68.03,69.03,70.03,71.03,72.03,73.03,74.03,75.03,76.03,77.03,78.03,79.03,80.03,81.03,82.03,83.03,84.03,85.03,86.03,87.03,88.03,89.03,90.03,91.03,92.03,93.03,94.03,95.03,96.03,97.03,98.03,99.03,100.03,101.03,102.03,103.03,104.03,105.03,106.03,107.03,108.03,109.03,110.03,111.03,26.11,27.11,114.03,115.03,116.03,117.03,118.03,119.03,120.03,121.03,122.03,123.03,124.03,125.03,36.11,37.11,128.03,129.03,130.03,131.03,132.03,133.03,134.03,135.03,136.03,137.03,138.03,139.03])
        self.assertTrue(att2.getFieldOn(att2.getMyGodFather().getMeshAtPosition((1,)),"YY").isEqualWithoutConsideringStr(exp3,1e-12))
        exp4=DataArrayDouble([0.02,1.02,2.02,3.02,4.02,5.02,6.02,7.02,8.02,9.02,10.02,11.02,12.02,13.02,14.02,15.02,16.02,17.02,18.02,19.02,20.02,21.02,22.02,23.02,24.02,25.02,26.02,27.02,28.02,29.02,30.02,31.02,32.02,33.02,34.02,35.02,36.02,37.02,38.02,39.02,40.02,41.02,42.02,43.02,44.02,45.02,46.02,47.02,48.02,49.02,50.02,51.02,52.02,53.02,54.02,55.02,56.02,57.02,58.02,59.02,60.02,61.02,62.02,63.02,64.02,65.02,66.02,67.02,68.02,69.02,70.02,71.02,72.02,73.02,74.02,75.02,76.02,77.02,78.02,79.02,80.02,81.02,82.02,83.02,84.02,85.02,86.02,87.02,88.02,89.02,90.02,91.02,92.02,93.02,94.02,95.02,96.02,97.02,98.02,99.02,100.02,101.02,102.02,103.02,104.02,105.02,106.02,107.02,108.02,109.02,110.02,111.02,112.02,113.02,30.03,31.03,32.03,33.03,34.03,35.03,36.03,37.03,38.03,39.03,124.02,125.02,126.02,127.02,44.03,45.03,46.03,47.03,48.03,49.03,50.03,51.03,52.03,53.03,138.02,139.02])
        self.assertTrue(att2.getFieldOn(att2.getMyGodFather().getMeshAtPosition((0,)),"YY").isEqualWithoutConsideringStr(exp4,1e-12))
        for pos,iot in [((),0.01),((0,0),0.04),((0,1),0.05),((0,2),0.06),((1,0),0.07),((1,1),0.08),((1,2),0.09),((1,3),0.10)]:
            vals=att2.getFieldOn(att2.getMyGodFather().getMeshAtPosition(pos),"YY")
            l=vals.getNumberOfTuples()
            exps=DataArrayDouble(l) ; exps.iota(iot)
            self.assertTrue(vals.isEqualWithoutConsideringStr(exps,1e-12))
            pass
        del att2
        ###
        att3.synchronizeCoarseToFineBetween(1,2)
        ###
        for pos in [(),(0,),(1,),(2,)]:
            self.assertTrue(att3.getFieldOn(att3.getMyGodFather().getMeshAtPosition(pos),"YY").isEqual(att4.getFieldOn(att4.getMyGodFather().getMeshAtPosition(pos),"YY"),1e-12))
            pass
        exp5=DataArrayDouble([57.02,57.02,58.02,58.02,58.02,59.02,59.02,59.02,60.02,60.02,60.02,61.02,61.02,61.02,62.02,62.02,62.02,63.02,63.02,63.02,64.02,64.02,64.02,65.02,65.02,65.02,66.02,66.02,66.02,67.02,67.02,67.02,68.02,68.02,57.02,57.02,58.02,58.02,58.02,59.02,59.02,59.02,60.02,60.02,60.02,61.02,61.02,61.02,62.02,62.02,62.02,63.02,63.02,63.02,64.02,64.02,64.02,65.02,65.02,65.02,66.02,66.02,66.02,67.02,67.02,67.02,68.02,68.02,71.02,71.02,72.02,72.02,72.02,73.02,73.02,73.02,74.02,74.02,74.02,75.02,75.02,75.02,76.02,76.02,76.02,77.02,77.02,77.02,78.02,78.02,78.02,79.02,79.02,79.02,80.02,80.02,80.02,81.02,81.02,81.02,82.02,82.02,71.02,71.02,72.02,72.02,72.02,73.02,73.02,73.02,74.02,74.02,74.02,75.02,75.02,75.02,76.02,76.02,76.02,77.02,77.02,77.02,78.02,78.02,78.02,79.02,79.02,79.02,80.02,80.02,80.02,81.02,81.02,81.02,82.02,82.02,71.02,71.02,72.02,72.02,72.02,73.02,73.02,73.02,74.02,74.02,74.02,75.02,75.02,75.02,76.02,76.02,76.02,77.02,77.02,77.02,78.02,78.02,78.02,79.02,79.02,79.02,80.02,80.02,80.02,81.02,81.02,81.02,82.02,82.02,85.02,85.02,86.02,86.02,86.02,87.02,87.02,87.02,88.02,88.02,88.02,89.02,89.02,89.02,90.02,90.02,90.02,91.02,91.02,91.02,92.02,92.02,92.02,93.02,93.02,93.02,94.02,94.02,94.02,95.02,95.02,95.02,96.02,96.02,85.02,85.02,86.02,86.02,86.02,87.02,87.02,87.02,88.02,88.02,88.02,89.02,89.02,89.02,90.02,90.02,90.02,91.02,91.02,91.02,92.02,92.02,92.02,93.02,93.02,93.02,94.02,94.02,94.02,95.02,95.02,95.02,96.02,96.02,85.02,85.02,86.02,86.02,86.02,87.02,87.02,87.02,88.02,88.02,88.02,89.02,89.02,89.02,90.02,90.02,90.02,91.02,91.02,91.02,92.02,92.02,92.02,93.02,93.02,93.02,94.02,94.02,94.02,95.02,95.02,95.02,96.02,96.02,99.02,99.02,100.02,100.02,100.02,101.02,101.02,101.02,102.02,102.02,102.02,103.02,103.02,103.02,104.02,104.02,104.02,105.02,105.02,105.02,106.02,106.02,106.02,107.02,107.02,107.02,108.02,108.02,108.02,109.02,109.02,109.02,110.02,110.02,99.02,99.02,100.02,100.02,100.02,101.02,101.02,101.02,102.02,102.02,102.02,103.02,103.02,103.02,104.02,104.02,104.02,105.02,105.02,105.02,106.02,106.02,106.02,107.02,107.02,107.02,108.02,108.02,108.02,109.02,109.02,109.02,110.02,110.02,99.02,99.02,100.02,100.02,100.02,101.02,101.02,101.02,102.02,102.02,102.02,103.02,103.02,103.02,104.02,104.02,104.02,105.02,105.02,105.02,106.02,106.02,106.02,107.02,107.02,107.02,108.02,108.02,108.02,109.02,109.02,109.02,110.02,110.02,113.02,113.02,114.02,114.02,114.02,115.02,115.02,115.02,116.02,116.02,116.02,117.02,117.02,117.02,118.02,118.02,118.02,119.02,119.02,119.02,120.02,120.02,120.02,121.02,121.02,121.02,122.02,122.02,122.02,123.02,123.02,123.02,124.02,124.02,113.02,113.02,114.02,114.02,114.02,115.02,115.02,115.02,116.02,116.02,116.02,117.02,117.02,117.02,118.02,118.02,118.02,119.02,119.02,119.02,120.02,120.02,120.02,121.02,121.02,121.02,122.02,122.02,122.02,123.02,123.02,123.02,124.02,124.02])
        self.assertTrue(att3.getFieldOn(att3.getMyGodFather().getMeshAtPosition((0,0)),"YY").isEqualWithoutConsideringStr(exp5,1e-12))
        exp6=DataArrayDouble([17.02,17.02,18.02,18.02,18.02,19.02,19.02,19.02,20.02,20.02,20.02,21.02,21.02,21.02,22.02,22.02,17.02,17.02,18.02,18.02,18.02,19.02,19.02,19.02,20.02,20.02,20.02,21.02,21.02,21.02,22.02,22.02,31.02,31.02,32.02,32.02,32.02,33.02,33.02,33.02,34.02,34.02,34.02,35.02,35.02,35.02,36.02,36.02,31.02,31.02,32.02,32.02,32.02,33.02,33.02,33.02,34.02,34.02,34.02,35.02,35.02,35.02,36.02,36.02,31.02,31.02,32.02,32.02,32.02,33.02,33.02,33.02,34.02,34.02,34.02,35.02,35.02,35.02,36.02,36.02,45.02,45.02,46.02,46.02,46.02,47.02,47.02,47.02,48.02,48.02,48.02,49.02,49.02,49.02,50.02,50.02,45.02,45.02,46.02,46.02,46.02,47.02,47.02,47.02,48.02,48.02,48.02,49.02,49.02,49.02,50.02,50.02,45.02,45.02,46.02,46.02,46.02,47.02,47.02,47.02,48.02,48.02,48.02,49.02,49.02,49.02,50.02,50.02,59.02,59.02,60.02,60.02,60.02,61.02,61.02,61.02,62.02,62.02,62.02,63.02,63.02,63.02,64.02,64.02,59.02,59.02,60.02,60.02,60.02,61.02,61.02,61.02,62.02,62.02,62.02,63.02,63.02,63.02,64.02,64.02,59.02,59.02,60.02,60.02,60.02,61.02,61.02,61.02,62.02,62.02,62.02,63.02,63.02,63.02,64.02,64.02,73.02,73.02,74.02,74.02,74.02,75.02,75.02,75.02,76.02,76.02,76.02,77.02,77.02,77.02,78.02,78.02,73.02,73.02,74.02,74.02,74.02,75.02,75.02,75.02,76.02,76.02,76.02,77.02,77.02,77.02,78.02,78.02])
        self.assertTrue(att3.getFieldOn(att3.getMyGodFather().getMeshAtPosition((0,1)),"YY").isEqualWithoutConsideringStr(exp6,1e-12))
        exp7=DataArrayDouble([49.02,49.02,50.02,50.02,50.02,51.02,51.02,51.02,52.02,52.02,52.02,53.02,53.02,53.02,54.02,54.02,49.02,49.02,50.02,50.02,50.02,51.02,51.02,51.02,52.02,52.02,52.02,53.02,53.02,53.02,54.02,54.02,63.02,63.02,64.02,64.02,64.02,65.02,65.02,65.02,66.02,66.02,66.02,67.02,67.02,67.02,68.02,68.02,63.02,63.02,64.02,64.02,64.02,65.02,65.02,65.02,66.02,66.02,66.02,67.02,67.02,67.02,68.02,68.02,63.02,63.02,64.02,64.02,64.02,65.02,65.02,65.02,66.02,66.02,66.02,67.02,67.02,67.02,68.02,68.02,77.02,77.02,78.02,78.02,78.02,79.02,79.02,79.02,80.02,80.02,80.02,81.02,81.02,81.02,82.02,82.02,77.02,77.02,78.02,78.02,78.02,79.02,79.02,79.02,80.02,80.02,80.02,81.02,81.02,81.02,82.02,82.02])
        self.assertTrue(att3.getFieldOn(att3.getMyGodFather().getMeshAtPosition((0,2)),"YY").isEqualWithoutConsideringStr(exp7,1e-12))
        exp8=DataArrayDouble([15.03,15.03,16.03,16.03,16.03,17.03,17.03,17.03,18.03,18.03,18.03,19.03,19.03,19.03,20.03,20.03,15.03,15.03,16.03,16.03,16.03,17.03,17.03,17.03,18.03,18.03,18.03,19.03,19.03,19.03,20.03,20.03,29.03,29.03,30.03,30.03,30.03,31.03,31.03,31.03,32.03,32.03,32.03,33.03,33.03,33.03,34.03,34.03,29.03,29.03,30.03,30.03,30.03,31.03,31.03,31.03,32.03,32.03,32.03,33.03,33.03,33.03,34.03,34.03,29.03,29.03,30.03,30.03,30.03,31.03,31.03,31.03,32.03,32.03,32.03,33.03,33.03,33.03,34.03,34.03,43.03,43.03,44.03,44.03,44.03,45.03,45.03,45.03,46.03,46.03,46.03,47.03,47.03,47.03,48.03,48.03,43.03,43.03,44.03,44.03,44.03,45.03,45.03,45.03,46.03,46.03,46.03,47.03,47.03,47.03,48.03,48.03,43.03,43.03,44.03,44.03,44.03,45.03,45.03,45.03,46.03,46.03,46.03,47.03,47.03,47.03,48.03,48.03,57.03,57.03,58.03,58.03,58.03,59.03,59.03,59.03,60.03,60.03,60.03,61.03,61.03,61.03,62.03,62.03,57.03,57.03,58.03,58.03,58.03,59.03,59.03,59.03,60.03,60.03,60.03,61.03,61.03,61.03,62.03,62.03,57.03,57.03,58.03,58.03,58.03,59.03,59.03,59.03,60.03,60.03,60.03,61.03,61.03,61.03,62.03,62.03,71.03,71.03,72.03,72.03,72.03,73.03,73.03,73.03,74.03,74.03,74.03,75.03,75.03,75.03,76.03,76.03,71.03,71.03,72.03,72.03,72.03,73.03,73.03,73.03,74.03,74.03,74.03,75.03,75.03,75.03,76.03,76.03,71.03,71.03,72.03,72.03,72.03,73.03,73.03,73.03,74.03,74.03,74.03,75.03,75.03,75.03,76.03,76.03,85.03,85.03,86.03,86.03,86.03,87.03,87.03,87.03,88.03,88.03,88.03,89.03,89.03,89.03,90.03,90.03,85.03,85.03,86.03,86.03,86.03,87.03,87.03,87.03,88.03,88.03,88.03,89.03,89.03,89.03,90.03,90.03,85.03,85.03,86.03,86.03,86.03,87.03,87.03,87.03,88.03,88.03,88.03,89.03,89.03,89.03,90.03,90.03,99.03,99.03,100.03,100.03,100.03,101.03,101.03,101.03,102.03,102.03,102.03,103.03,103.03,103.03,104.03,104.03,99.03,99.03,100.03,100.03,100.03,101.03,101.03,101.03,102.03,102.03,102.03,103.03,103.03,103.03,104.03,104.03,99.03,99.03,100.03,100.03,100.03,101.03,101.03,101.03,102.03,102.03,102.03,103.03,103.03,103.03,104.03,104.03,113.03,113.03,114.03,114.03,114.03,115.03,115.03,115.03,116.03,116.03,116.03,117.03,117.03,117.03,118.03,118.03,113.03,113.03,114.03,114.03,114.03,115.03,115.03,115.03,116.03,116.03,116.03,117.03,117.03,117.03,118.03,118.03])
        self.assertTrue(att3.getFieldOn(att3.getMyGodFather().getMeshAtPosition((1,0)),"YY").isEqualWithoutConsideringStr(exp8,1e-12))
        exp9=DataArrayDouble([22.03,22.03,23.03,23.03,23.03,24.03,24.03,24.03,25.03,25.03,25.03,26.03,26.03,22.03,22.03,23.03,23.03,23.03,24.03,24.03,24.03,25.03,25.03,25.03,26.03,26.03,36.03,36.03,37.03,37.03,37.03,38.03,38.03,38.03,39.03,39.03,39.03,40.03,40.03,36.03,36.03,37.03,37.03,37.03,38.03,38.03,38.03,39.03,39.03,39.03,40.03,40.03,36.03,36.03,37.03,37.03,37.03,38.03,38.03,38.03,39.03,39.03,39.03,40.03,40.03,50.03,50.03,51.03,51.03,51.03,52.03,52.03,52.03,53.03,53.03,53.03,54.03,54.03,50.03,50.03,51.03,51.03,51.03,52.03,52.03,52.03,53.03,53.03,53.03,54.03,54.03,50.03,50.03,51.03,51.03,51.03,52.03,52.03,52.03,53.03,53.03,53.03,54.03,54.03,64.03,64.03,65.03,65.03,65.03,66.03,66.03,66.03,67.03,67.03,67.03,68.03,68.03,64.03,64.03,65.03,65.03,65.03,66.03,66.03,66.03,67.03,67.03,67.03,68.03,68.03,64.03,64.03,65.03,65.03,65.03,66.03,66.03,66.03,67.03,67.03,67.03,68.03,68.03,78.03,78.03,79.03,79.03,79.03,80.03,80.03,80.03,81.03,81.03,81.03,82.03,82.03,78.03,78.03,79.03,79.03,79.03,80.03,80.03,80.03,81.03,81.03,81.03,82.03,82.03,78.03,78.03,79.03,79.03,79.03,80.03,80.03,80.03,81.03,81.03,81.03,82.03,82.03,92.03,92.03,93.03,93.03,93.03,94.03,94.03,94.03,95.03,95.03,95.03,96.03,96.03,92.03,92.03,93.03,93.03,93.03,94.03,94.03,94.03,95.03,95.03,95.03,96.03,96.03])
        self.assertTrue(att3.getFieldOn(att3.getMyGodFather().getMeshAtPosition((1,1)),"YY").isEqualWithoutConsideringStr(exp9,1e-12))
        exp10=DataArrayDouble([19.03,19.03,20.03,20.03,20.03,21.03,21.03,21.03,22.03,22.03,22.03,23.03,23.03,19.03,19.03,20.03,20.03,20.03,21.03,21.03,21.03,22.03,22.03,22.03,23.03,23.03,33.03,33.03,34.03,34.03,34.03,35.03,35.03,35.03,36.03,36.03,36.03,37.03,37.03,33.03,33.03,34.03,34.03,34.03,35.03,35.03,35.03,36.03,36.03,36.03,37.03,37.03,33.03,33.03,34.03,34.03,34.03,35.03,35.03,35.03,36.03,36.03,36.03,37.03,37.03,47.03,47.03,48.03,48.03,48.03,49.03,49.03,49.03,50.03,50.03,50.03,51.03,51.03,47.03,47.03,48.03,48.03,48.03,49.03,49.03,49.03,50.03,50.03,50.03,51.03,51.03,47.03,47.03,48.03,48.03,48.03,49.03,49.03,49.03,50.03,50.03,50.03,51.03,51.03,61.03,61.03,62.03,62.03,62.03,63.03,63.03,63.03,64.03,64.03,64.03,65.03,65.03,61.03,61.03,62.03,62.03,62.03,63.03,63.03,63.03,64.03,64.03,64.03,65.03,65.03,61.03,61.03,62.03,62.03,62.03,63.03,63.03,63.03,64.03,64.03,64.03,65.03,65.03,75.03,75.03,76.03,76.03,76.03,77.03,77.03,77.03,78.03,78.03,78.03,79.03,79.03,75.03,75.03,76.03,76.03,76.03,77.03,77.03,77.03,78.03,78.03,78.03,79.03,79.03])
        self.assertTrue(att3.getFieldOn(att3.getMyGodFather().getMeshAtPosition((1,2)),"YY").isEqualWithoutConsideringStr(exp10,1e-12))
        exp11=DataArrayDouble([61.03,61.03,62.03,62.03,62.03,63.03,63.03,63.03,64.03,64.03,64.03,65.03,65.03,61.03,61.03,62.03,62.03,62.03,63.03,63.03,63.03,64.03,64.03,64.03,65.03,65.03,75.03,75.03,76.03,76.03,76.03,77.03,77.03,77.03,78.03,78.03,78.03,79.03,79.03,75.03,75.03,76.03,76.03,76.03,77.03,77.03,77.03,78.03,78.03,78.03,79.03,79.03,75.03,75.03,76.03,76.03,76.03,77.03,77.03,77.03,78.03,78.03,78.03,79.03,79.03,89.03,89.03,90.03,90.03,90.03,91.03,91.03,91.03,92.03,92.03,92.03,93.03,93.03,89.03,89.03,90.03,90.03,90.03,91.03,91.03,91.03,92.03,92.03,92.03,93.03,93.03,89.03,89.03,90.03,90.03,90.03,91.03,91.03,91.03,92.03,92.03,92.03,93.03,93.03,103.03,103.03,104.03,104.03,104.03,105.03,105.03,105.03,106.03,106.03,106.03,107.03,107.03,103.03,103.03,104.03,104.03,104.03,105.03,105.03,105.03,106.03,106.03,106.03,107.03,107.03,103.03,103.03,104.03,104.03,104.03,105.03,105.03,105.03,106.03,106.03,106.03,107.03,107.03,117.03,117.03,118.03,118.03,118.03,119.03,119.03,119.03,120.03,120.03,120.03,121.03,121.03,117.03,117.03,118.03,118.03,118.03,119.03,119.03,119.03,120.03,120.03,120.03,121.03,121.03])
        self.assertTrue(att3.getFieldOn(att3.getMyGodFather().getMeshAtPosition((1,3)),"YY").isEqualWithoutConsideringStr(exp11,1e-12))
        del att3
        ###
        att4.synchronizeAllGhostZonesAtASpecifiedLevel(2)
        for pos in [(),(0,),(1,),(2,)]:
            self.assertTrue(att4.getFieldOn(att4.getMyGodFather().getMeshAtPosition(pos),"YY").isEqual(att5.getFieldOn(att5.getMyGodFather().getMeshAtPosition(pos),"YY"),1e-12))
            pass
        exp12=DataArrayDouble([0.04,1.04,2.04,3.04,4.04,5.04,6.04,7.04,146.05,147.05,148.05,149.05,150.05,151.05,152.05,153.05,154.05,155.05,156.05,157.05,50.06,51.06,52.06,53.06,54.06,55.06,56.06,57.06,58.06,59.06,60.06,61.06,32.04,33.04,34.04,35.04,36.04,37.04,38.04,39.04,40.04,41.04,162.05,163.05,164.05,165.05,166.05,167.05,168.05,169.05,170.05,171.05,172.05,173.05,66.06,67.06,68.06,69.06,70.06,71.06,72.06,73.06,74.06,75.06,76.06,77.06,66.04,67.04,68.04,69.04,70.04,71.04,72.04,73.04,74.04,75.04,76.04,77.04,78.04,79.04,80.04,81.04,82.04,83.04,84.04,85.04,86.04,87.04,88.04,89.04,90.04,91.04,92.04,93.04,94.04,95.04,96.04,97.04,98.04,99.04,100.04,101.04,102.04,103.04,104.04,105.04,106.04,107.04,108.04,109.04,110.04,111.04,112.04,113.04,114.04,115.04,116.04,117.04,118.04,119.04,120.04,121.04,122.04,123.04,124.04,125.04,126.04,127.04,128.04,129.04,130.04,131.04,132.04,133.04,134.04,135.04,136.04,137.04,138.04,139.04,140.04,141.04,142.04,143.04,144.04,145.04,146.04,147.04,148.04,149.04,150.04,151.04,152.04,153.04,154.04,155.04,156.04,157.04,158.04,159.04,160.04,161.04,162.04,163.04,164.04,165.04,166.04,167.04,168.04,169.04,170.04,171.04,172.04,173.04,174.04,175.04,176.04,177.04,178.04,179.04,180.04,181.04,182.04,183.04,184.04,185.04,186.04,187.04,188.04,189.04,190.04,191.04,192.04,193.04,194.04,195.04,196.04,197.04,198.04,199.04,200.04,201.04,202.04,203.04,204.04,205.04,206.04,207.04,208.04,209.04,210.04,211.04,212.04,213.04,214.04,215.04,216.04,217.04,218.04,219.04,220.04,221.04,222.04,223.04,224.04,225.04,226.04,227.04,228.04,229.04,230.04,231.04,232.04,233.04,234.04,235.04,236.04,237.04,238.04,239.04,240.04,241.04,242.04,243.04,244.04,245.04,246.04,247.04,248.04,249.04,250.04,251.04,252.04,253.04,254.04,255.04,256.04,257.04,258.04,259.04,260.04,261.04,262.04,263.04,264.04,265.04,266.04,267.04,268.04,269.04,270.04,271.04,272.04,273.04,274.04,275.04,276.04,277.04,278.04,279.04,280.04,281.04,282.04,283.04,284.04,285.04,286.04,287.04,288.04,289.04,290.04,291.04,292.04,293.04,294.04,295.04,296.04,297.04,298.04,299.04,300.04,301.04,302.04,303.04,304.04,305.04,306.04,307.04,308.04,309.04,310.04,311.04,312.04,313.04,314.04,315.04,316.04,317.04,318.04,319.04,320.04,321.04,322.04,323.04,324.04,325.04,326.04,327.04,328.04,329.04,330.04,331.04,332.04,333.04,334.04,335.04,336.04,337.04,338.04,339.04,340.04,341.04,342.04,343.04,344.04,345.04,346.04,347.04,348.04,349.04,350.04,351.04,352.04,353.04,354.04,355.04,356.04,357.04,358.04,359.04,360.04,361.04,362.04,363.04,364.04,365.04,366.04,367.04,368.04,369.04,370.04,371.04,372.04,373.04,374.04,375.04,34.07,35.07,36.07,37.07,38.07,39.07,40.07,41.07,42.07,43.07,44.07,45.07,28.09,29.09,30.09,31.09,32.09,33.09,34.09,35.09,36.09,28.08,29.08,30.08,31.08,32.08,33.08,34.08,35.08,36.08,406.04,407.04,408.04,409.04,50.07,51.07,52.07,53.07,54.07,55.07,56.07,57.07,58.07,59.07,60.07,61.07,41.09,42.09,43.09,44.09,45.09,46.09,47.09,48.09,49.09,41.08,42.08,43.08,44.08,45.08,46.08,47.08,48.08,49.08,440.04,441.04])
        self.assertTrue(att4.getFieldOn(att4.getMyGodFather().getMeshAtPosition((0,0)),"YY").isEqualWithoutConsideringStr(exp12,1e-12))
        exp13=DataArrayDouble([[0.05,1.05,2.05,3.05,4.05,5.05,6.05,7.05,8.05,9.05,10.05,11.05,12.05,13.05,14.05,15.05,16.05,17.05,18.05,19.05,20.05,21.05,22.05,23.05,24.05,25.05,26.05,27.05,28.05,29.05,30.05,31.05,32.05,33.05,34.05,35.05,36.05,37.05,38.05,39.05,40.05,41.05,42.05,43.05,44.05,45.05,46.05,47.05,48.05,49.05,50.05,51.05,52.05,53.05,54.05,55.05,56.05,57.05,58.05,59.05,60.05,61.05,62.05,63.05,64.05,65.05,66.05,67.05,68.05,69.05,70.05,71.05,72.05,73.05,74.05,75.05,76.05,77.05,78.05,79.05,80.05,81.05,82.05,83.05,84.05,85.05,86.05,87.05,88.05,89.05,90.05,91.05,92.05,93.05,94.05,95.05,96.05,97.05,98.05,99.05,100.05,101.05,102.05,103.05,104.05,105.05,106.05,107.05,108.05,109.05,110.05,111.05,112.05,113.05,114.05,115.05,116.05,117.05,118.05,119.05,120.05,121.05,122.05,123.05,124.05,125.05,126.05,127.05,128.05,129.05,130.05,131.05,132.05,133.05,134.05,135.05,136.05,137.05,138.05,139.05,140.05,141.05,34.06,35.06,144.05,145.05,146.05,147.05,148.05,149.05,150.05,151.05,152.05,153.05,154.05,155.05,156.05,157.05,50.06,51.06,160.05,161.05,162.05,163.05,164.05,165.05,166.05,167.05,168.05,169.05,170.05,171.05,172.05,173.05,66.06,67.06,74.04,75.04,76.04,77.04,78.04,79.04,80.04,81.04,82.04,83.04,84.04,85.04,86.04,87.04,88.04,89.04,108.04,109.04,110.04,111.04,112.04,113.04,114.04,115.04,116.04,117.04,118.04,119.04,120.04,121.04,122.04,123.04]])
        self.assertTrue(att4.getFieldOn(att4.getMyGodFather().getMeshAtPosition((0,1)),"YY").isEqualWithoutConsideringStr(exp13,1e-12))
        exp14=DataArrayDouble([108.05,109.05,2.06,3.06,4.06,5.06,6.06,7.06,8.06,9.06,10.06,11.06,12.06,13.06,14.06,15.06,124.05,125.05,18.06,19.06,20.06,21.06,22.06,23.06,24.06,25.06,26.06,27.06,28.06,29.06,30.06,31.06,140.05,141.05,34.06,35.06,36.06,37.06,38.06,39.06,40.06,41.06,42.06,43.06,44.06,45.06,46.06,47.06,156.05,157.05,50.06,51.06,52.06,53.06,54.06,55.06,56.06,57.06,58.06,59.06,60.06,61.06,62.06,63.06,172.05,173.05,66.06,67.06,68.06,69.06,70.06,71.06,72.06,73.06,74.06,75.06,76.06,77.06,78.06,79.06,86.04,87.04,88.04,89.04,90.04,91.04,92.04,93.04,94.04,95.04,96.04,97.04,98.04,99.04,94.06,95.06,120.04,121.04,122.04,123.04,124.04,125.04,126.04,127.04,128.04,129.04,130.04,131.04,132.04,133.04,110.06,111.06])
        self.assertTrue(att4.getFieldOn(att4.getMyGodFather().getMeshAtPosition((0,2)),"YY").isEqualWithoutConsideringStr(exp14,1e-12))
        exp15=DataArrayDouble([0.07,1.07,308.04,309.04,310.04,311.04,312.04,313.04,314.04,315.04,316.04,317.04,318.04,319.04,320.04,321.04,16.07,17.07,342.04,343.04,344.04,345.04,346.04,347.04,348.04,349.04,350.04,351.04,352.04,353.04,354.04,355.04,32.07,33.07,34.07,35.07,36.07,37.07,38.07,39.07,40.07,41.07,42.07,43.07,44.07,45.07,28.09,29.09,48.07,49.07,50.07,51.07,52.07,53.07,54.07,55.07,56.07,57.07,58.07,59.07,60.07,61.07,41.09,42.09,64.07,65.07,66.07,67.07,68.07,69.07,70.07,71.07,72.07,73.07,74.07,75.07,76.07,77.07,54.09,55.09,80.07,81.07,82.07,83.07,84.07,85.07,86.07,87.07,88.07,89.07,90.07,91.07,92.07,93.07,67.09,68.09,96.07,97.07,98.07,99.07,100.07,101.07,102.07,103.07,104.07,105.07,106.07,107.07,108.07,109.07,80.09,81.09,112.07,113.07,114.07,115.07,116.07,117.07,118.07,119.07,120.07,121.07,122.07,123.07,124.07,125.07,93.09,94.09,128.07,129.07,130.07,131.07,132.07,133.07,134.07,135.07,136.07,137.07,138.07,139.07,140.07,141.07,106.09,107.09,144.07,145.07,146.07,147.07,148.07,149.07,150.07,151.07,152.07,153.07,154.07,155.07,156.07,157.07,119.09,120.09,160.07,161.07,162.07,163.07,164.07,165.07,166.07,167.07,168.07,169.07,170.07,171.07,172.07,173.07,132.09,133.09,176.07,177.07,178.07,179.07,180.07,181.07,182.07,183.07,184.07,185.07,186.07,187.07,188.07,189.07,28.1,29.1,192.07,193.07,194.07,195.07,196.07,197.07,198.07,199.07,200.07,201.07,202.07,203.07,204.07,205.07,41.1,42.1,208.07,209.07,210.07,211.07,212.07,213.07,214.07,215.07,216.07,217.07,218.07,219.07,220.07,221.07,54.1,55.1,224.07,225.07,226.07,227.07,228.07,229.07,230.07,231.07,232.07,233.07,234.07,235.07,236.07,237.07,67.1,68.1,240.07,241.07,242.07,243.07,244.07,245.07,246.07,247.07,248.07,249.07,250.07,251.07,252.07,253.07,80.1,81.1,256.07,257.07,258.07,259.07,260.07,261.07,262.07,263.07,264.07,265.07,266.07,267.07,268.07,269.07,93.1,94.1,272.07,273.07,274.07,275.07,276.07,277.07,278.07,279.07,280.07,281.07,282.07,283.07,284.07,285.07,106.1,107.1,288.07,289.07,290.07,291.07,292.07,293.07,294.07,295.07,296.07,297.07,298.07,299.07,300.07,301.07,119.1,120.1,304.07,305.07,306.07,307.07,308.07,309.07,310.07,311.07,312.07,313.07,314.07,315.07,316.07,317.07,132.1,133.1,320.07,321.07,322.07,323.07,324.07,325.07,326.07,327.07,328.07,329.07,330.07,331.07,332.07,333.07,334.07,335.07,336.07,337.07,338.07,339.07,340.07,341.07,342.07,343.07,344.07,345.07,346.07,347.07,348.07,349.07,350.07,351.07])
        self.assertTrue(att4.getFieldOn(att4.getMyGodFather().getMeshAtPosition((1,0)),"YY").isEqualWithoutConsideringStr(exp15,1e-12))
        exp16=DataArrayDouble([327.04,328.04,329.04,330.04,331.04,332.04,333.04,334.04,335.04,336.04,337.04,11.08,12.08,361.04,362.04,363.04,364.04,365.04,366.04,367.04,368.04,369.04,370.04,371.04,24.08,25.08,35.09,36.09,28.08,29.08,30.08,31.08,32.08,33.08,34.08,35.08,36.08,37.08,38.08,48.09,49.09,41.08,42.08,43.08,44.08,45.08,46.08,47.08,48.08,49.08,50.08,51.08,61.09,62.09,54.08,55.08,56.08,57.08,58.08,59.08,60.08,61.08,62.08,63.08,64.08,74.09,75.09,67.08,68.08,69.08,70.08,71.08,72.08,73.08,74.08,75.08,76.08,77.08,87.09,88.09,80.08,81.08,82.08,83.08,84.08,85.08,86.08,87.08,88.08,89.08,90.08,100.09,101.09,93.08,94.08,95.08,96.08,97.08,98.08,99.08,100.08,101.08,102.08,103.08,113.09,114.09,106.08,107.08,108.08,109.08,110.08,111.08,112.08,113.08,114.08,115.08,116.08,126.09,127.09,119.08,120.08,121.08,122.08,123.08,124.08,125.08,126.08,127.08,128.08,129.08,139.09,140.09,132.08,133.08,134.08,135.08,136.08,137.08,138.08,139.08,140.08,141.08,142.08,35.1,36.1,145.08,146.08,147.08,148.08,149.08,150.08,151.08,152.08,153.08,154.08,155.08,48.1,49.1,158.08,159.08,160.08,161.08,162.08,163.08,164.08,165.08,166.08,167.08,168.08,61.1,62.1,171.08,172.08,173.08,174.08,175.08,176.08,177.08,178.08,179.08,180.08,181.08,74.1,75.1,184.08,185.08,186.08,187.08,188.08,189.08,190.08,191.08,192.08,193.08,194.08,87.1,88.1,197.08,198.08,199.08,200.08,201.08,202.08,203.08,204.08,205.08,206.08,207.08])
        self.assertTrue(att4.getFieldOn(att4.getMyGodFather().getMeshAtPosition((1,1)),"YY").isEqualWithoutConsideringStr(exp16,1e-12))
        exp17=DataArrayDouble([318.04,319.04,320.04,321.04,322.04,323.04,324.04,325.04,326.04,327.04,328.04,329.04,330.04,352.04,353.04,354.04,355.04,356.04,357.04,358.04,359.04,360.04,361.04,362.04,363.04,364.04,44.07,45.07,28.09,29.09,30.09,31.09,32.09,33.09,34.09,35.09,36.09,28.08,29.08,60.07,61.07,41.09,42.09,43.09,44.09,45.09,46.09,47.09,48.09,49.09,41.08,42.08,76.07,77.07,54.09,55.09,56.09,57.09,58.09,59.09,60.09,61.09,62.09,54.08,55.08,92.07,93.07,67.09,68.09,69.09,70.09,71.09,72.09,73.09,74.09,75.09,67.08,68.08,108.07,109.07,80.09,81.09,82.09,83.09,84.09,85.09,86.09,87.09,88.09,80.08,81.08,124.07,125.07,93.09,94.09,95.09,96.09,97.09,98.09,99.09,100.09,101.09,93.08,94.08,140.07,141.07,106.09,107.09,108.09,109.09,110.09,111.09,112.09,113.09,114.09,106.08,107.08,156.07,157.07,119.09,120.09,121.09,122.09,123.09,124.09,125.09,126.09,127.09,119.08,120.08,172.07,173.07,132.09,133.09,134.09,135.09,136.09,137.09,138.09,139.09,140.09,132.08,133.08,188.07,189.07,28.1,29.1,30.1,31.1,32.1,33.1,34.1,35.1,36.1,145.08,146.08,204.07,205.07,41.1,42.1,43.1,44.1,45.1,46.1,47.1,48.1,49.1,158.08,159.08])
        self.assertTrue(att4.getFieldOn(att4.getMyGodFather().getMeshAtPosition((1,2)),"YY").isEqualWithoutConsideringStr(exp17,1e-12))
        exp18=DataArrayDouble([156.07,157.07,119.09,120.09,121.09,122.09,123.09,124.09,125.09,126.09,127.09,119.08,120.08,172.07,173.07,132.09,133.09,134.09,135.09,136.09,137.09,138.09,139.09,140.09,132.08,133.08,188.07,189.07,28.1,29.1,30.1,31.1,32.1,33.1,34.1,35.1,36.1,145.08,146.08,204.07,205.07,41.1,42.1,43.1,44.1,45.1,46.1,47.1,48.1,49.1,158.08,159.08,220.07,221.07,54.1,55.1,56.1,57.1,58.1,59.1,60.1,61.1,62.1,171.08,172.08,236.07,237.07,67.1,68.1,69.1,70.1,71.1,72.1,73.1,74.1,75.1,76.1,77.1,252.07,253.07,80.1,81.1,82.1,83.1,84.1,85.1,86.1,87.1,88.1,89.1,90.1,268.07,269.07,93.1,94.1,95.1,96.1,97.1,98.1,99.1,100.1,101.1,102.1,103.1,284.07,285.07,106.1,107.1,108.1,109.1,110.1,111.1,112.1,113.1,114.1,115.1,116.1,300.07,301.07,119.1,120.1,121.1,122.1,123.1,124.1,125.1,126.1,127.1,128.1,129.1,316.07,317.07,132.1,133.1,134.1,135.1,136.1,137.1,138.1,139.1,140.1,141.1,142.1,143.1,144.1,145.1,146.1,147.1,148.1,149.1,150.1,151.1,152.1,153.1,154.1,155.1,156.1,157.1,158.1,159.1,160.1,161.1,162.1,163.1,164.1,165.1,166.1,167.1,168.1])
        self.assertTrue(att4.getFieldOn(att4.getMyGodFather().getMeshAtPosition((1,3)),"YY").isEqualWithoutConsideringStr(exp18,1e-12))
        del att4
        ###
        att5.synchronizeAllGhostZonesAtASpecifiedLevelUsingOnlyFather(2)
        for pos in [(),(0,),(1,),(2,)]:
            self.assertTrue(att5.getFieldOn(att5.getMyGodFather().getMeshAtPosition(pos),"YY").isEqual(att6.getFieldOn(att6.getMyGodFather().getMeshAtPosition(pos),"YY"),1e-12))
            pass
        att5.buildCellFieldOnWithGhost(att5.getMyGodFather().getMeshAtPosition((0,0)),"YY")
        exp19=DataArrayDouble([57.02,57.02,58.02,58.02,58.02,59.02,59.02,59.02,60.02,60.02,60.02,61.02,61.02,61.02,62.02,62.02,62.02,63.02,63.02,63.02,64.02,64.02,64.02,65.02,65.02,65.02,66.02,66.02,66.02,67.02,67.02,67.02,68.02,68.02,57.02,57.02,58.02,58.02,58.02,59.02,59.02,59.02,60.02,60.02,60.02,61.02,61.02,61.02,62.02,62.02,62.02,63.02,63.02,63.02,64.02,64.02,64.02,65.02,65.02,65.02,66.02,66.02,66.02,67.02,67.02,67.02,68.02,68.02,71.02,71.02,70.04,71.04,72.04,73.04,74.04,75.04,76.04,77.04,78.04,79.04,80.04,81.04,82.04,83.04,84.04,85.04,86.04,87.04,88.04,89.04,90.04,91.04,92.04,93.04,94.04,95.04,96.04,97.04,98.04,99.04,82.02,82.02,71.02,71.02,104.04,105.04,106.04,107.04,108.04,109.04,110.04,111.04,112.04,113.04,114.04,115.04,116.04,117.04,118.04,119.04,120.04,121.04,122.04,123.04,124.04,125.04,126.04,127.04,128.04,129.04,130.04,131.04,132.04,133.04,82.02,82.02,71.02,71.02,138.04,139.04,140.04,141.04,142.04,143.04,144.04,145.04,146.04,147.04,148.04,149.04,150.04,151.04,152.04,153.04,154.04,155.04,156.04,157.04,158.04,159.04,160.04,161.04,162.04,163.04,164.04,165.04,166.04,167.04,82.02,82.02,85.02,85.02,172.04,173.04,174.04,175.04,176.04,177.04,178.04,179.04,180.04,181.04,182.04,183.04,184.04,185.04,186.04,187.04,188.04,189.04,190.04,191.04,192.04,193.04,194.04,195.04,196.04,197.04,198.04,199.04,200.04,201.04,96.02,96.02,85.02,85.02,206.04,207.04,208.04,209.04,210.04,211.04,212.04,213.04,214.04,215.04,216.04,217.04,218.04,219.04,220.04,221.04,222.04,223.04,224.04,225.04,226.04,227.04,228.04,229.04,230.04,231.04,232.04,233.04,234.04,235.04,96.02,96.02,85.02,85.02,240.04,241.04,242.04,243.04,244.04,245.04,246.04,247.04,248.04,249.04,250.04,251.04,252.04,253.04,254.04,255.04,256.04,257.04,258.04,259.04,260.04,261.04,262.04,263.04,264.04,265.04,266.04,267.04,268.04,269.04,96.02,96.02,99.02,99.02,274.04,275.04,276.04,277.04,278.04,279.04,280.04,281.04,282.04,283.04,284.04,285.04,286.04,287.04,288.04,289.04,290.04,291.04,292.04,293.04,294.04,295.04,296.04,297.04,298.04,299.04,300.04,301.04,302.04,303.04,110.02,110.02,99.02,99.02,308.04,309.04,310.04,311.04,312.04,313.04,314.04,315.04,316.04,317.04,318.04,319.04,320.04,321.04,322.04,323.04,324.04,325.04,326.04,327.04,328.04,329.04,330.04,331.04,332.04,333.04,334.04,335.04,336.04,337.04,110.02,110.02,99.02,99.02,342.04,343.04,344.04,345.04,346.04,347.04,348.04,349.04,350.04,351.04,352.04,353.04,354.04,355.04,356.04,357.04,358.04,359.04,360.04,361.04,362.04,363.04,364.04,365.04,366.04,367.04,368.04,369.04,370.04,371.04,110.02,110.02,113.02,113.02,114.02,114.02,114.02,115.02,115.02,115.02,116.02,116.02,116.02,117.02,117.02,117.02,118.02,118.02,118.02,119.02,119.02,119.02,120.02,120.02,120.02,121.02,121.02,121.02,122.02,122.02,122.02,123.02,123.02,123.02,124.02,124.02,113.02,113.02,114.02,114.02,114.02,115.02,115.02,115.02,116.02,116.02,116.02,117.02,117.02,117.02,118.02,118.02,118.02,119.02,119.02,119.02,120.02,120.02,120.02,121.02,121.02,121.02,122.02,122.02,122.02,123.02,123.02,123.02,124.02,124.02])
        self.assertTrue(att5.getFieldOn(att5.getMyGodFather().getMeshAtPosition((0,0)),"YY").isEqualWithoutConsideringStr(exp19,1e-12))
        exp20=DataArrayDouble([17.02,17.02,18.02,18.02,18.02,19.02,19.02,19.02,20.02,20.02,20.02,21.02,21.02,21.02,22.02,22.02,17.02,17.02,18.02,18.02,18.02,19.02,19.02,19.02,20.02,20.02,20.02,21.02,21.02,21.02,22.02,22.02,31.02,31.02,34.05,35.05,36.05,37.05,38.05,39.05,40.05,41.05,42.05,43.05,44.05,45.05,36.02,36.02,31.02,31.02,50.05,51.05,52.05,53.05,54.05,55.05,56.05,57.05,58.05,59.05,60.05,61.05,36.02,36.02,31.02,31.02,66.05,67.05,68.05,69.05,70.05,71.05,72.05,73.05,74.05,75.05,76.05,77.05,36.02,36.02,45.02,45.02,82.05,83.05,84.05,85.05,86.05,87.05,88.05,89.05,90.05,91.05,92.05,93.05,50.02,50.02,45.02,45.02,98.05,99.05,100.05,101.05,102.05,103.05,104.05,105.05,106.05,107.05,108.05,109.05,50.02,50.02,45.02,45.02,114.05,115.05,116.05,117.05,118.05,119.05,120.05,121.05,122.05,123.05,124.05,125.05,50.02,50.02,59.02,59.02,130.05,131.05,132.05,133.05,134.05,135.05,136.05,137.05,138.05,139.05,140.05,141.05,64.02,64.02,59.02,59.02,146.05,147.05,148.05,149.05,150.05,151.05,152.05,153.05,154.05,155.05,156.05,157.05,64.02,64.02,59.02,59.02,162.05,163.05,164.05,165.05,166.05,167.05,168.05,169.05,170.05,171.05,172.05,173.05,64.02,64.02,73.02,73.02,74.02,74.02,74.02,75.02,75.02,75.02,76.02,76.02,76.02,77.02,77.02,77.02,78.02,78.02,73.02,73.02,74.02,74.02,74.02,75.02,75.02,75.02,76.02,76.02,76.02,77.02,77.02,77.02,78.02,78.02])
        self.assertTrue(att5.getFieldOn(att5.getMyGodFather().getMeshAtPosition((0,1)),"YY").isEqualWithoutConsideringStr(exp20,1e-12))
        exp21=DataArrayDouble([49.02,49.02,50.02,50.02,50.02,51.02,51.02,51.02,52.02,52.02,52.02,53.02,53.02,53.02,54.02,54.02,49.02,49.02,50.02,50.02,50.02,51.02,51.02,51.02,52.02,52.02,52.02,53.02,53.02,53.02,54.02,54.02,63.02,63.02,34.06,35.06,36.06,37.06,38.06,39.06,40.06,41.06,42.06,43.06,44.06,45.06,68.02,68.02,63.02,63.02,50.06,51.06,52.06,53.06,54.06,55.06,56.06,57.06,58.06,59.06,60.06,61.06,68.02,68.02,63.02,63.02,66.06,67.06,68.06,69.06,70.06,71.06,72.06,73.06,74.06,75.06,76.06,77.06,68.02,68.02,77.02,77.02,78.02,78.02,78.02,79.02,79.02,79.02,80.02,80.02,80.02,81.02,81.02,81.02,82.02,82.02,77.02,77.02,78.02,78.02,78.02,79.02,79.02,79.02,80.02,80.02,80.02,81.02,81.02,81.02,82.02,82.02])
        self.assertTrue(att5.getFieldOn(att5.getMyGodFather().getMeshAtPosition((0,2)),"YY").isEqualWithoutConsideringStr(exp21,1e-12))
        exp22=DataArrayDouble([15.03,15.03,16.03,16.03,16.03,17.03,17.03,17.03,18.03,18.03,18.03,19.03,19.03,19.03,20.03,20.03,15.03,15.03,16.03,16.03,16.03,17.03,17.03,17.03,18.03,18.03,18.03,19.03,19.03,19.03,20.03,20.03,29.03,29.03,34.07,35.07,36.07,37.07,38.07,39.07,40.07,41.07,42.07,43.07,44.07,45.07,34.03,34.03,29.03,29.03,50.07,51.07,52.07,53.07,54.07,55.07,56.07,57.07,58.07,59.07,60.07,61.07,34.03,34.03,29.03,29.03,66.07,67.07,68.07,69.07,70.07,71.07,72.07,73.07,74.07,75.07,76.07,77.07,34.03,34.03,43.03,43.03,82.07,83.07,84.07,85.07,86.07,87.07,88.07,89.07,90.07,91.07,92.07,93.07,48.03,48.03,43.03,43.03,98.07,99.07,100.07,101.07,102.07,103.07,104.07,105.07,106.07,107.07,108.07,109.07,48.03,48.03,43.03,43.03,114.07,115.07,116.07,117.07,118.07,119.07,120.07,121.07,122.07,123.07,124.07,125.07,48.03,48.03,57.03,57.03,130.07,131.07,132.07,133.07,134.07,135.07,136.07,137.07,138.07,139.07,140.07,141.07,62.03,62.03,57.03,57.03,146.07,147.07,148.07,149.07,150.07,151.07,152.07,153.07,154.07,155.07,156.07,157.07,62.03,62.03,57.03,57.03,162.07,163.07,164.07,165.07,166.07,167.07,168.07,169.07,170.07,171.07,172.07,173.07,62.03,62.03,71.03,71.03,178.07,179.07,180.07,181.07,182.07,183.07,184.07,185.07,186.07,187.07,188.07,189.07,76.03,76.03,71.03,71.03,194.07,195.07,196.07,197.07,198.07,199.07,200.07,201.07,202.07,203.07,204.07,205.07,76.03,76.03,71.03,71.03,210.07,211.07,212.07,213.07,214.07,215.07,216.07,217.07,218.07,219.07,220.07,221.07,76.03,76.03,85.03,85.03,226.07,227.07,228.07,229.07,230.07,231.07,232.07,233.07,234.07,235.07,236.07,237.07,90.03,90.03,85.03,85.03,242.07,243.07,244.07,245.07,246.07,247.07,248.07,249.07,250.07,251.07,252.07,253.07,90.03,90.03,85.03,85.03,258.07,259.07,260.07,261.07,262.07,263.07,264.07,265.07,266.07,267.07,268.07,269.07,90.03,90.03,99.03,99.03,274.07,275.07,276.07,277.07,278.07,279.07,280.07,281.07,282.07,283.07,284.07,285.07,104.03,104.03,99.03,99.03,290.07,291.07,292.07,293.07,294.07,295.07,296.07,297.07,298.07,299.07,300.07,301.07,104.03,104.03,99.03,99.03,306.07,307.07,308.07,309.07,310.07,311.07,312.07,313.07,314.07,315.07,316.07,317.07,104.03,104.03,113.03,113.03,114.03,114.03,114.03,115.03,115.03,115.03,116.03,116.03,116.03,117.03,117.03,117.03,118.03,118.03,113.03,113.03,114.03,114.03,114.03,115.03,115.03,115.03,116.03,116.03,116.03,117.03,117.03,117.03,118.03,118.03])
        self.assertTrue(att5.getFieldOn(att5.getMyGodFather().getMeshAtPosition((1,0)),"YY").isEqualWithoutConsideringStr(exp22,1e-12))
        exp23=DataArrayDouble([22.03,22.03,23.03,23.03,23.03,24.03,24.03,24.03,25.03,25.03,25.03,26.03,26.03,22.03,22.03,23.03,23.03,23.03,24.03,24.03,24.03,25.03,25.03,25.03,26.03,26.03,36.03,36.03,28.08,29.08,30.08,31.08,32.08,33.08,34.08,35.08,36.08,40.03,40.03,36.03,36.03,41.08,42.08,43.08,44.08,45.08,46.08,47.08,48.08,49.08,40.03,40.03,36.03,36.03,54.08,55.08,56.08,57.08,58.08,59.08,60.08,61.08,62.08,40.03,40.03,50.03,50.03,67.08,68.08,69.08,70.08,71.08,72.08,73.08,74.08,75.08,54.03,54.03,50.03,50.03,80.08,81.08,82.08,83.08,84.08,85.08,86.08,87.08,88.08,54.03,54.03,50.03,50.03,93.08,94.08,95.08,96.08,97.08,98.08,99.08,100.08,101.08,54.03,54.03,64.03,64.03,106.08,107.08,108.08,109.08,110.08,111.08,112.08,113.08,114.08,68.03,68.03,64.03,64.03,119.08,120.08,121.08,122.08,123.08,124.08,125.08,126.08,127.08,68.03,68.03,64.03,64.03,132.08,133.08,134.08,135.08,136.08,137.08,138.08,139.08,140.08,68.03,68.03,78.03,78.03,145.08,146.08,147.08,148.08,149.08,150.08,151.08,152.08,153.08,82.03,82.03,78.03,78.03,158.08,159.08,160.08,161.08,162.08,163.08,164.08,165.08,166.08,82.03,82.03,78.03,78.03,171.08,172.08,173.08,174.08,175.08,176.08,177.08,178.08,179.08,82.03,82.03,92.03,92.03,93.03,93.03,93.03,94.03,94.03,94.03,95.03,95.03,95.03,96.03,96.03,92.03,92.03,93.03,93.03,93.03,94.03,94.03,94.03,95.03,95.03,95.03,96.03,96.03])
        self.assertTrue(att5.getFieldOn(att5.getMyGodFather().getMeshAtPosition((1,1)),"YY").isEqualWithoutConsideringStr(exp23,1e-12))
        exp24=DataArrayDouble([19.03,19.03,20.03,20.03,20.03,21.03,21.03,21.03,22.03,22.03,22.03,23.03,23.03,19.03,19.03,20.03,20.03,20.03,21.03,21.03,21.03,22.03,22.03,22.03,23.03,23.03,33.03,33.03,28.09,29.09,30.09,31.09,32.09,33.09,34.09,35.09,36.09,37.03,37.03,33.03,33.03,41.09,42.09,43.09,44.09,45.09,46.09,47.09,48.09,49.09,37.03,37.03,33.03,33.03,54.09,55.09,56.09,57.09,58.09,59.09,60.09,61.09,62.09,37.03,37.03,47.03,47.03,67.09,68.09,69.09,70.09,71.09,72.09,73.09,74.09,75.09,51.03,51.03,47.03,47.03,80.09,81.09,82.09,83.09,84.09,85.09,86.09,87.09,88.09,51.03,51.03,47.03,47.03,93.09,94.09,95.09,96.09,97.09,98.09,99.09,100.09,101.09,51.03,51.03,61.03,61.03,106.09,107.09,108.09,109.09,110.09,111.09,112.09,113.09,114.09,65.03,65.03,61.03,61.03,119.09,120.09,121.09,122.09,123.09,124.09,125.09,126.09,127.09,65.03,65.03,61.03,61.03,132.09,133.09,134.09,135.09,136.09,137.09,138.09,139.09,140.09,65.03,65.03,75.03,75.03,76.03,76.03,76.03,77.03,77.03,77.03,78.03,78.03,78.03,79.03,79.03,75.03,75.03,76.03,76.03,76.03,77.03,77.03,77.03,78.03,78.03,78.03,79.03,79.03])
        self.assertTrue(att5.getFieldOn(att5.getMyGodFather().getMeshAtPosition((1,2)),"YY").isEqualWithoutConsideringStr(exp24,1e-12))
        exp25=DataArrayDouble([61.03,61.03,62.03,62.03,62.03,63.03,63.03,63.03,64.03,64.03,64.03,65.03,65.03,61.03,61.03,62.03,62.03,62.03,63.03,63.03,63.03,64.03,64.03,64.03,65.03,65.03,75.03,75.03,28.1,29.1,30.1,31.1,32.1,33.1,34.1,35.1,36.1,79.03,79.03,75.03,75.03,41.1,42.1,43.1,44.1,45.1,46.1,47.1,48.1,49.1,79.03,79.03,75.03,75.03,54.1,55.1,56.1,57.1,58.1,59.1,60.1,61.1,62.1,79.03,79.03,89.03,89.03,67.1,68.1,69.1,70.1,71.1,72.1,73.1,74.1,75.1,93.03,93.03,89.03,89.03,80.1,81.1,82.1,83.1,84.1,85.1,86.1,87.1,88.1,93.03,93.03,89.03,89.03,93.1,94.1,95.1,96.1,97.1,98.1,99.1,100.1,101.1,93.03,93.03,103.03,103.03,106.1,107.1,108.1,109.1,110.1,111.1,112.1,113.1,114.1,107.03,107.03,103.03,103.03,119.1,120.1,121.1,122.1,123.1,124.1,125.1,126.1,127.1,107.03,107.03,103.03,103.03,132.1,133.1,134.1,135.1,136.1,137.1,138.1,139.1,140.1,107.03,107.03,117.03,117.03,118.03,118.03,118.03,119.03,119.03,119.03,120.03,120.03,120.03,121.03,121.03,117.03,117.03,118.03,118.03,118.03,119.03,119.03,119.03,120.03,120.03,120.03,121.03,121.03])
        self.assertTrue(att5.getFieldOn(att5.getMyGodFather().getMeshAtPosition((1,3)),"YY").isEqualWithoutConsideringStr(exp25,1e-12))
        pass

    def testSwig2AMR11(self):
        """ Some tests in 3D with CondenseFineToCoarseGhost and SpreadCoarseToFineGhost"""
        coarse=DataArrayDouble((6+4)*(7+4)*(5+4)) ; coarse.iota()
        fine=DataArrayDouble((4*2+4)*(2*3+4)*(3*4+4))
        MEDCouplingIMesh.SpreadCoarseToFineGhost(coarse,[6,7,5],fine,[(1,5),(2,4),(1,4)],[2,3,4],2)
        exp0=DataArrayDouble([252.,252.,253.,253.,254.,254.,255.,255.,256.,256.,257.,257.,252.,252.,253.,253.,254.,254.,255.,255.,256.,256.,257.,257.,262.,262.,263.,263.,264.,264.,265.,265.,266.,266.,267.,267.,262.,262.,263.,263.,264.,264.,265.,265.,266.,266.,267.,267.,262.,262.,263.,263.,264.,264.,265.,265.,266.,266.,267.,267.,272.,272.,273.,273.,274.,274.,275.,275.,276.,276.,277.,277.,272.,272.,273.,273.,274.,274.,275.,275.,276.,276.,277.,277.,272.,272.,273.,273.,274.,274.,275.,275.,276.,276.,277.,277.,282.,282.,283.,283.,284.,284.,285.,285.,286.,286.,287.,287.,282.,282.,283.,283.,284.,284.,285.,285.,286.,286.,287.,287.])
        exp1=DataArrayDouble([362.,362.,363.,363.,364.,364.,365.,365.,366.,366.,367.,367.,362.,362.,363.,363.,364.,364.,365.,365.,366.,366.,367.,367.,372.,372.,373.,373.,374.,374.,375.,375.,376.,376.,377.,377.,372.,372.,373.,373.,374.,374.,375.,375.,376.,376.,377.,377.,372.,372.,373.,373.,374.,374.,375.,375.,376.,376.,377.,377.,382.,382.,383.,383.,384.,384.,385.,385.,386.,386.,387.,387.,382.,382.,383.,383.,384.,384.,385.,385.,386.,386.,387.,387.,382.,382.,383.,383.,384.,384.,385.,385.,386.,386.,387.,387.,392.,392.,393.,393.,394.,394.,395.,395.,396.,396.,397.,397.,392.,392.,393.,393.,394.,394.,395.,395.,396.,396.,397.,397.])
        exp2=DataArrayDouble([472.,472.,473.,473.,474.,474.,475.,475.,476.,476.,477.,477.,472.,472.,473.,473.,474.,474.,475.,475.,476.,476.,477.,477.,482.,482.,483.,483.,484.,484.,485.,485.,486.,486.,487.,487.,482.,482.,483.,483.,484.,484.,485.,485.,486.,486.,487.,487.,482.,482.,483.,483.,484.,484.,485.,485.,486.,486.,487.,487.,492.,492.,493.,493.,494.,494.,495.,495.,496.,496.,497.,497.,492.,492.,493.,493.,494.,494.,495.,495.,496.,496.,497.,497.,492.,492.,493.,493.,494.,494.,495.,495.,496.,496.,497.,497.,502.,502.,503.,503.,504.,504.,505.,505.,506.,506.,507.,507.,502.,502.,503.,503.,504.,504.,505.,505.,506.,506.,507.,507.])
        exp3=DataArrayDouble([582.,582.,583.,583.,584.,584.,585.,585.,586.,586.,587.,587.,582.,582.,583.,583.,584.,584.,585.,585.,586.,586.,587.,587.,592.,592.,593.,593.,594.,594.,595.,595.,596.,596.,597.,597.,592.,592.,593.,593.,594.,594.,595.,595.,596.,596.,597.,597.,592.,592.,593.,593.,594.,594.,595.,595.,596.,596.,597.,597.,602.,602.,603.,603.,604.,604.,605.,605.,606.,606.,607.,607.,602.,602.,603.,603.,604.,604.,605.,605.,606.,606.,607.,607.,602.,602.,603.,603.,604.,604.,605.,605.,606.,606.,607.,607.,612.,612.,613.,613.,614.,614.,615.,615.,616.,616.,617.,617.,612.,612.,613.,613.,614.,614.,615.,615.,616.,616.,617.,617.])
        exp4=DataArrayDouble([692.,692.,693.,693.,694.,694.,695.,695.,696.,696.,697.,697.,692.,692.,693.,693.,694.,694.,695.,695.,696.,696.,697.,697.,702.,702.,703.,703.,704.,704.,705.,705.,706.,706.,707.,707.,702.,702.,703.,703.,704.,704.,705.,705.,706.,706.,707.,707.,702.,702.,703.,703.,704.,704.,705.,705.,706.,706.,707.,707.,712.,712.,713.,713.,714.,714.,715.,715.,716.,716.,717.,717.,712.,712.,713.,713.,714.,714.,715.,715.,716.,716.,717.,717.,712.,712.,713.,713.,714.,714.,715.,715.,716.,716.,717.,717.,722.,722.,723.,723.,724.,724.,725.,725.,726.,726.,727.,727.,722.,722.,723.,723.,724.,724.,725.,725.,726.,726.,727.,727.])
        exp=DataArrayDouble.Aggregate([exp0,exp0,exp1,exp1,exp1,exp1,exp2,exp2,exp2,exp2,exp3,exp3,exp3,exp3,exp4,exp4])
        self.assertTrue(fine.isEqual(exp,1e-12))
        #
        fine.iota()
        coarse.iota(0.5)
        MEDCouplingIMesh.CondenseFineToCoarseGhost([6,7,5],fine,[(1,5),(2,4),(1,4)],[2,3,4],coarse,2)
        amr=MEDCouplingCartesianAMRMesh("mesh",3,[7,8,6],[0.,0.,0.],[1.,1.,1.])
        amr.addPatch([(1,5),(2,4),(1,4)],[2,3,4])
        att=MEDCouplingAMRAttribute(amr,[("YY",1)],2)
        att.alloc()
        exp1=DataArrayDouble(990) ; exp1.iota(0.5)
        ids=DataArrayInt([373,374,375,376,383,384,385,386,483,484,485,486,493,494,495,496,593,594,595,596,603,604,605,606])
        vals=DataArrayDouble([11004.,11052.,11100.,11148.,11868.,11916.,11964.,12012.,22524.,22572.,22620.,22668.,23388.,23436.,23484.,23532.,34044.,34092.,34140.,34188.,34908.,34956.,35004.,35052.])
        exp1[ids]=vals
        self.assertTrue(coarse.isEqual(exp1,1e-12))
        #
        MEDCouplingStructuredMesh.MultiplyPartOf([10,11,9],[(3,7),(4,6),(3,6)],1/24.,coarse)
        exp2=DataArrayDouble(990) ; exp2.iota(0.5)
        exp2[ids]=vals/24.
        self.assertTrue(coarse.isEqual(exp2,1e-12))
        #
        coarse.iota(0.5) ; fine.iota(0.1)
        MEDCouplingIMesh.SpreadCoarseToFineGhostZone(coarse,[6,7,5],fine,[(1,5),(2,4),(1,4)],[2,3,4],2)
        #
        coarse.iota(0.5) ; fine.iota(0.1)
        MEDCouplingIMesh.SpreadCoarseToFineGhostZone(coarse,[6,7,5],fine,[(1,5),(2,4),(1,4)],[2,3,4],2)
        exp00=DataArrayDouble.Aggregate([exp0,exp0]) ; exp00+=0.5
        self.assertTrue(fine[:240].isEqual(exp00,1e-12))
        exp44=DataArrayDouble.Aggregate([exp4,exp4]) ; exp44+=0.5
        self.assertTrue(fine[-240:].isEqual(exp44,1e-12))
        self.assertTrue(fine[240:-240].isEqual(DataArrayDouble([362.5,362.5,363.5,363.5,364.5,364.5,365.5,365.5,366.5,366.5,367.5,367.5,362.5,362.5,363.5,363.5,364.5,364.5,365.5,365.5,366.5,366.5,367.5,367.5,372.5,372.5,266.1,267.1,268.1,269.1,270.1,271.1,272.1,273.1,377.5,377.5,372.5,372.5,278.1,279.1,280.1,281.1,282.1,283.1,284.1,285.1,377.5,377.5,372.5,372.5,290.1,291.1,292.1,293.1,294.1,295.1,296.1,297.1,377.5,377.5,382.5,382.5,302.1,303.1,304.1,305.1,306.1,307.1,308.1,309.1,387.5,387.5,382.5,382.5,314.1,315.1,316.1,317.1,318.1,319.1,320.1,321.1,387.5,387.5,382.5,382.5,326.1,327.1,328.1,329.1,330.1,331.1,332.1,333.1,387.5,387.5,392.5,392.5,393.5,393.5,394.5,394.5,395.5,395.5,396.5,396.5,397.5,397.5,392.5,392.5,393.5,393.5,394.5,394.5,395.5,395.5,396.5,396.5,397.5,397.5,362.5,362.5,363.5,363.5,364.5,364.5,365.5,365.5,366.5,366.5,367.5,367.5,362.5,362.5,363.5,363.5,364.5,364.5,365.5,365.5,366.5,366.5,367.5,367.5,372.5,372.5,386.1,387.1,388.1,389.1,390.1,391.1,392.1,393.1,377.5,377.5,372.5,372.5,398.1,399.1,400.1,401.1,402.1,403.1,404.1,405.1,377.5,377.5,372.5,372.5,410.1,411.1,412.1,413.1,414.1,415.1,416.1,417.1,377.5,377.5,382.5,382.5,422.1,423.1,424.1,425.1,426.1,427.1,428.1,429.1,387.5,387.5,382.5,382.5,434.1,435.1,436.1,437.1,438.1,439.1,440.1,441.1,387.5,387.5,382.5,382.5,446.1,447.1,448.1,449.1,450.1,451.1,452.1,453.1,387.5,387.5,392.5,392.5,393.5,393.5,394.5,394.5,395.5,395.5,396.5,396.5,397.5,397.5,392.5,392.5,393.5,393.5,394.5,394.5,395.5,395.5,396.5,396.5,397.5,397.5,362.5,362.5,363.5,363.5,364.5,364.5,365.5,365.5,366.5,366.5,367.5,367.5,362.5,362.5,363.5,363.5,364.5,364.5,365.5,365.5,366.5,366.5,367.5,367.5,372.5,372.5,506.1,507.1,508.1,509.1,510.1,511.1,512.1,513.1,377.5,377.5,372.5,372.5,518.1,519.1,520.1,521.1,522.1,523.1,524.1,525.1,377.5,377.5,372.5,372.5,530.1,531.1,532.1,533.1,534.1,535.1,536.1,537.1,377.5,377.5,382.5,382.5,542.1,543.1,544.1,545.1,546.1,547.1,548.1,549.1,387.5,387.5,382.5,382.5,554.1,555.1,556.1,557.1,558.1,559.1,560.1,561.1,387.5,387.5,382.5,382.5,566.1,567.1,568.1,569.1,570.1,571.1,572.1,573.1,387.5,387.5,392.5,392.5,393.5,393.5,394.5,394.5,395.5,395.5,396.5,396.5,397.5,397.5,392.5,392.5,393.5,393.5,394.5,394.5,395.5,395.5,396.5,396.5,397.5,397.5,362.5,362.5,363.5,363.5,364.5,364.5,365.5,365.5,366.5,366.5,367.5,367.5,362.5,362.5,363.5,363.5,364.5,364.5,365.5,365.5,366.5,366.5,367.5,367.5,372.5,372.5,626.1,627.1,628.1,629.1,630.1,631.1,632.1,633.1,377.5,377.5,372.5,372.5,638.1,639.1,640.1,641.1,642.1,643.1,644.1,645.1,377.5,377.5,372.5,372.5,650.1,651.1,652.1,653.1,654.1,655.1,656.1,657.1,377.5,377.5,382.5,382.5,662.1,663.1,664.1,665.1,666.1,667.1,668.1,669.1,387.5,387.5,382.5,382.5,674.1,675.1,676.1,677.1,678.1,679.1,680.1,681.1,387.5,387.5,382.5,382.5,686.1,687.1,688.1,689.1,690.1,691.1,692.1,693.1,387.5,387.5,392.5,392.5,393.5,393.5,394.5,394.5,395.5,395.5,396.5,396.5,397.5,397.5,392.5,392.5,393.5,393.5,394.5,394.5,395.5,395.5,396.5,396.5,397.5,397.5,472.5,472.5,473.5,473.5,474.5,474.5,475.5,475.5,476.5,476.5,477.5,477.5,472.5,472.5,473.5,473.5,474.5,474.5,475.5,475.5,476.5,476.5,477.5,477.5,482.5,482.5,746.1,747.1,748.1,749.1,750.1,751.1,752.1,753.1,487.5,487.5,482.5,482.5,758.1,759.1,760.1,761.1,762.1,763.1,764.1,765.1,487.5,487.5,482.5,482.5,770.1,771.1,772.1,773.1,774.1,775.1,776.1,777.1,487.5,487.5,492.5,492.5,782.1,783.1,784.1,785.1,786.1,787.1,788.1,789.1,497.5,497.5,492.5,492.5,794.1,795.1,796.1,797.1,798.1,799.1,800.1,801.1,497.5,497.5,492.5,492.5,806.1,807.1,808.1,809.1,810.1,811.1,812.1,813.1,497.5,497.5,502.5,502.5,503.5,503.5,504.5,504.5,505.5,505.5,506.5,506.5,507.5,507.5,502.5,502.5,503.5,503.5,504.5,504.5,505.5,505.5,506.5,506.5,507.5,507.5,472.5,472.5,473.5,473.5,474.5,474.5,475.5,475.5,476.5,476.5,477.5,477.5,472.5,472.5,473.5,473.5,474.5,474.5,475.5,475.5,476.5,476.5,477.5,477.5,482.5,482.5,866.1,867.1,868.1,869.1,870.1,871.1,872.1,873.1,487.5,487.5,482.5,482.5,878.1,879.1,880.1,881.1,882.1,883.1,884.1,885.1,487.5,487.5,482.5,482.5,890.1,891.1,892.1,893.1,894.1,895.1,896.1,897.1,487.5,487.5,492.5,492.5,902.1,903.1,904.1,905.1,906.1,907.1,908.1,909.1,497.5,497.5,492.5,492.5,914.1,915.1,916.1,917.1,918.1,919.1,920.1,921.1,497.5,497.5,492.5,492.5,926.1,927.1,928.1,929.1,930.1,931.1,932.1,933.1,497.5,497.5,502.5,502.5,503.5,503.5,504.5,504.5,505.5,505.5,506.5,506.5,507.5,507.5,502.5,502.5,503.5,503.5,504.5,504.5,505.5,505.5,506.5,506.5,507.5,507.5,472.5,472.5,473.5,473.5,474.5,474.5,475.5,475.5,476.5,476.5,477.5,477.5,472.5,472.5,473.5,473.5,474.5,474.5,475.5,475.5,476.5,476.5,477.5,477.5,482.5,482.5,986.1,987.1,988.1,989.1,990.1,991.1,992.1,993.1,487.5,487.5,482.5,482.5,998.1,999.1,1000.1,1001.1,1002.1,1003.1,1004.1,1005.1,487.5,487.5,482.5,482.5,1010.1,1011.1,1012.1,1013.1,1014.1,1015.1,1016.1,1017.1,487.5,487.5,492.5,492.5,1022.1,1023.1,1024.1,1025.1,1026.1,1027.1,1028.1,1029.1,497.5,497.5,492.5,492.5,1034.1,1035.1,1036.1,1037.1,1038.1,1039.1,1040.1,1041.1,497.5,497.5,492.5,492.5,1046.1,1047.1,1048.1,1049.1,1050.1,1051.1,1052.1,1053.1,497.5,497.5,502.5,502.5,503.5,503.5,504.5,504.5,505.5,505.5,506.5,506.5,507.5,507.5,502.5,502.5,503.5,503.5,504.5,504.5,505.5,505.5,506.5,506.5,507.5,507.5,472.5,472.5,473.5,473.5,474.5,474.5,475.5,475.5,476.5,476.5,477.5,477.5,472.5,472.5,473.5,473.5,474.5,474.5,475.5,475.5,476.5,476.5,477.5,477.5,482.5,482.5,1106.1,1107.1,1108.1,1109.1,1110.1,1111.1,1112.1,1113.1,487.5,487.5,482.5,482.5,1118.1,1119.1,1120.1,1121.1,1122.1,1123.1,1124.1,1125.1,487.5,487.5,482.5,482.5,1130.1,1131.1,1132.1,1133.1,1134.1,1135.1,1136.1,1137.1,487.5,487.5,492.5,492.5,1142.1,1143.1,1144.1,1145.1,1146.1,1147.1,1148.1,1149.1,497.5,497.5,492.5,492.5,1154.1,1155.1,1156.1,1157.1,1158.1,1159.1,1160.1,1161.1,497.5,497.5,492.5,492.5,1166.1,1167.1,1168.1,1169.1,1170.1,1171.1,1172.1,1173.1,497.5,497.5,502.5,502.5,503.5,503.5,504.5,504.5,505.5,505.5,506.5,506.5,507.5,507.5,502.5,502.5,503.5,503.5,504.5,504.5,505.5,505.5,506.5,506.5,507.5,507.5,582.5,582.5,583.5,583.5,584.5,584.5,585.5,585.5,586.5,586.5,587.5,587.5,582.5,582.5,583.5,583.5,584.5,584.5,585.5,585.5,586.5,586.5,587.5,587.5,592.5,592.5,1226.1,1227.1,1228.1,1229.1,1230.1,1231.1,1232.1,1233.1,597.5,597.5,592.5,592.5,1238.1,1239.1,1240.1,1241.1,1242.1,1243.1,1244.1,1245.1,597.5,597.5,592.5,592.5,1250.1,1251.1,1252.1,1253.1,1254.1,1255.1,1256.1,1257.1,597.5,597.5,602.5,602.5,1262.1,1263.1,1264.1,1265.1,1266.1,1267.1,1268.1,1269.1,607.5,607.5,602.5,602.5,1274.1,1275.1,1276.1,1277.1,1278.1,1279.1,1280.1,1281.1,607.5,607.5,602.5,602.5,1286.1,1287.1,1288.1,1289.1,1290.1,1291.1,1292.1,1293.1,607.5,607.5,612.5,612.5,613.5,613.5,614.5,614.5,615.5,615.5,616.5,616.5,617.5,617.5,612.5,612.5,613.5,613.5,614.5,614.5,615.5,615.5,616.5,616.5,617.5,617.5,582.5,582.5,583.5,583.5,584.5,584.5,585.5,585.5,586.5,586.5,587.5,587.5,582.5,582.5,583.5,583.5,584.5,584.5,585.5,585.5,586.5,586.5,587.5,587.5,592.5,592.5,1346.1,1347.1,1348.1,1349.1,1350.1,1351.1,1352.1,1353.1,597.5,597.5,592.5,592.5,1358.1,1359.1,1360.1,1361.1,1362.1,1363.1,1364.1,1365.1,597.5,597.5,592.5,592.5,1370.1,1371.1,1372.1,1373.1,1374.1,1375.1,1376.1,1377.1,597.5,597.5,602.5,602.5,1382.1,1383.1,1384.1,1385.1,1386.1,1387.1,1388.1,1389.1,607.5,607.5,602.5,602.5,1394.1,1395.1,1396.1,1397.1,1398.1,1399.1,1400.1,1401.1,607.5,607.5,602.5,602.5,1406.1,1407.1,1408.1,1409.1,1410.1,1411.1,1412.1,1413.1,607.5,607.5,612.5,612.5,613.5,613.5,614.5,614.5,615.5,615.5,616.5,616.5,617.5,617.5,612.5,612.5,613.5,613.5,614.5,614.5,615.5,615.5,616.5,616.5,617.5,617.5,582.5,582.5,583.5,583.5,584.5,584.5,585.5,585.5,586.5,586.5,587.5,587.5,582.5,582.5,583.5,583.5,584.5,584.5,585.5,585.5,586.5,586.5,587.5,587.5,592.5,592.5,1466.1,1467.1,1468.1,1469.1,1470.1,1471.1,1472.1,1473.1,597.5,597.5,592.5,592.5,1478.1,1479.1,1480.1,1481.1,1482.1,1483.1,1484.1,1485.1,597.5,597.5,592.5,592.5,1490.1,1491.1,1492.1,1493.1,1494.1,1495.1,1496.1,1497.1,597.5,597.5,602.5,602.5,1502.1,1503.1,1504.1,1505.1,1506.1,1507.1,1508.1,1509.1,607.5,607.5,602.5,602.5,1514.1,1515.1,1516.1,1517.1,1518.1,1519.1,1520.1,1521.1,607.5,607.5,602.5,602.5,1526.1,1527.1,1528.1,1529.1,1530.1,1531.1,1532.1,1533.1,607.5,607.5,612.5,612.5,613.5,613.5,614.5,614.5,615.5,615.5,616.5,616.5,617.5,617.5,612.5,612.5,613.5,613.5,614.5,614.5,615.5,615.5,616.5,616.5,617.5,617.5,582.5,582.5,583.5,583.5,584.5,584.5,585.5,585.5,586.5,586.5,587.5,587.5,582.5,582.5,583.5,583.5,584.5,584.5,585.5,585.5,586.5,586.5,587.5,587.5,592.5,592.5,1586.1,1587.1,1588.1,1589.1,1590.1,1591.1,1592.1,1593.1,597.5,597.5,592.5,592.5,1598.1,1599.1,1600.1,1601.1,1602.1,1603.1,1604.1,1605.1,597.5,597.5,592.5,592.5,1610.1,1611.1,1612.1,1613.1,1614.1,1615.1,1616.1,1617.1,597.5,597.5,602.5,602.5,1622.1,1623.1,1624.1,1625.1,1626.1,1627.1,1628.1,1629.1,607.5,607.5,602.5,602.5,1634.1,1635.1,1636.1,1637.1,1638.1,1639.1,1640.1,1641.1,607.5,607.5,602.5,602.5,1646.1,1647.1,1648.1,1649.1,1650.1,1651.1,1652.1,1653.1,607.5,607.5,612.5,612.5,613.5,613.5,614.5,614.5,615.5,615.5,616.5,616.5,617.5,617.5,612.5,612.5,613.5,613.5,614.5,614.5,615.5,615.5,616.5,616.5,617.5,617.5]),1e-12))
        pass

    def testSwig2AMR12(self):
        """ This test check the MEDCouplingAMRAttribute.projectTo method."""
        amr0=MEDCouplingCartesianAMRMesh("mesh",2,[11,11],[0.,0.],[1.,1.])
        amr0.addPatch([(3,8),(0,3)],[2,2])
        amr0.addPatch([(3,8),(3,6)],[2,2])
        att0=MEDCouplingAMRAttribute(amr0,[("YY",1)],2)
        att0.alloc()
        att0.getFieldOn(amr0,"YY").iota(0.01)
        att0.getFieldOn(amr0[0].getMesh(),"YY").iota(0.02)
        att0.getFieldOn(amr0[1].getMesh(),"YY").iota(0.03)
        amr1=MEDCouplingCartesianAMRMesh("mesh",2,[11,11],[0.,0.],[1.,1.])
        amr1.addPatch([(2,5),(1,4)],[2,2])
        att1=att0.projectTo(amr1)
        self.assertTrue(att1.getFieldOn(amr1,"YY").isEqualWithoutConsideringStr(att0.getFieldOn(amr0,"YY"),1e-12))
        self.assertTrue(att1.getFieldOn(amr1[0].getMesh(),"YY").isEqualWithoutConsideringStr(DataArrayDouble([31.01,31.01,32.01,32.01,33.01,33.01,34.01,34.01,35.01,35.01,31.01,31.01,32.01,32.01,33.01,33.01,34.01,34.01,35.01,35.01,45.01,45.01,46.01,46.01,58.02,59.02,60.02,61.02,49.01,49.01,45.01,45.01,46.01,46.01,72.02,73.02,74.02,75.02,49.01,49.01,59.01,59.01,60.01,60.01,86.02,87.02,88.02,89.02,63.01,63.01,59.01,59.01,60.01,60.01,100.02,101.02,102.02,103.02,63.01,63.01,73.01,73.01,74.01,74.01,30.03,31.03,32.03,33.03,77.01,77.01,73.01,73.01,74.01,74.01,44.03,45.03,46.03,47.03,77.01,77.01,87.01,87.01,88.01,88.01,89.01,89.01,90.01,90.01,91.01,91.01,87.01,87.01,88.01,88.01,89.01,89.01,90.01,90.01,91.01,91.01]),1e-12))
        #
        amr0=MEDCouplingCartesianAMRMesh("mesh",2,[11,11],[0.,0.],[1.,1.])
        amr0.addPatch([(2,5),(2,7)],[2,2])
        amr0.addPatch([(5,8),(2,7)],[2,2])
        att0=MEDCouplingAMRAttribute(amr0,[("YY",1)],2)
        att0.alloc()
        att0.getFieldOn(amr0,"YY").iota(0.01)
        att0.getFieldOn(amr0[0].getMesh(),"YY").iota(0.02)
        att0.getFieldOn(amr0[1].getMesh(),"YY").iota(0.03)
        amr1=MEDCouplingCartesianAMRMesh("mesh",2,[11,11],[0.,0.],[1.,1.])
        amr1.addPatch([(3,6),(2,7)],[2,2])
        amr1.addPatch([(6,9),(2,7)],[2,2])
        att1=att0.projectTo(amr1)
        self.assertTrue(att1.getFieldOn(amr1,"YY").isEqual(att0.getFieldOn(amr0,"YY"),1e-12))
        self.assertTrue(att1.getFieldOn(amr1[0].getMesh(),"YY").isEqualWithoutConsideringStr(DataArrayDouble([46.01,46.01,47.01,47.01,48.01,48.01,49.01,49.01,50.01,50.01,46.01,46.01,47.01,47.01,48.01,48.01,49.01,49.01,50.01,50.01,60.01,60.01,24.02,25.02,26.02,27.02,22.03,23.03,64.01,64.01,60.01,60.01,34.02,35.02,36.02,37.02,32.03,33.03,64.01,64.01,74.01,74.01,44.02,45.02,46.02,47.02,42.03,43.03,78.01,78.01,74.01,74.01,54.02,55.02,56.02,57.02,52.03,53.03,78.01,78.01,88.01,88.01,64.02,65.02,66.02,67.02,62.03,63.03,92.01,92.01,88.01,88.01,74.02,75.02,76.02,77.02,72.03,73.03,92.01,92.01,102.01,102.01,84.02,85.02,86.02,87.02,82.03,83.03,106.01,106.01,102.01,102.01,94.02,95.02,96.02,97.02,92.03,93.03,106.01,106.01,116.01,116.01,104.02,105.02,106.02,107.02,102.03,103.03,120.01,120.01,116.01,116.01,114.02,115.02,116.02,117.02,112.03,113.03,120.01,120.01,130.01,130.01,131.01,131.01,132.01,132.01,133.01,133.01,134.01,134.01,130.01,130.01,131.01,131.01,132.01,132.01,133.01,133.01,134.01,134.01]),1e-12))
        self.assertTrue(att1.getFieldOn(amr1[1].getMesh(),"YY").isEqualWithoutConsideringStr(DataArrayDouble([49.01,49.01,50.01,50.01,51.01,51.01,52.01,52.01,53.01,53.01,49.01,49.01,50.01,50.01,51.01,51.01,52.01,52.01,53.01,53.01,63.01,63.01,24.03,25.03,26.03,27.03,66.01,66.01,67.01,67.01,63.01,63.01,34.03,35.03,36.03,37.03,66.01,66.01,67.01,67.01,77.01,77.01,44.03,45.03,46.03,47.03,80.01,80.01,81.01,81.01,77.01,77.01,54.03,55.03,56.03,57.03,80.01,80.01,81.01,81.01,91.01,91.01,64.03,65.03,66.03,67.03,94.01,94.01,95.01,95.01,91.01,91.01,74.03,75.03,76.03,77.03,94.01,94.01,95.01,95.01,105.01,105.01,84.03,85.03,86.03,87.03,108.01,108.01,109.01,109.01,105.01,105.01,94.03,95.03,96.03,97.03,108.01,108.01,109.01,109.01,119.01,119.01,104.03,105.03,106.03,107.03,122.01,122.01,123.01,123.01,119.01,119.01,114.03,115.03,116.03,117.03,122.01,122.01,123.01,123.01,133.01,133.01,134.01,134.01,135.01,135.01,136.01,136.01,137.01,137.01,133.01,133.01,134.01,134.01,135.01,135.01,136.01,136.01,137.01,137.01]),1e-12))
        pass

    def testSwig2AMR13(self):
        """ non regression test"""
        for fact,len1,len2 in [([2,2],64,48),([3,3],100,70),([4,4],144,96)]:
            amr=MEDCouplingCartesianAMRMesh("mesh",2,[5,5],[0.,0.],[1.,1.])
            amr.addPatch([(1,3),(0,2)],fact)
            amr.addPatch([(1,3),(3,4)],fact)
            att=MEDCouplingAMRAttribute(amr,[("YY",1)],2)
            att.alloc()
            att.getFieldOn(amr,"YY").iota(0.1)
            att.getFieldOn(amr[0].getMesh(),"YY").iota(0.2)
            att.getFieldOn(amr[1].getMesh(),"YY").iota(0.3)
            att.synchronizeAllGhostZonesOfDirectChidrenOf(amr)
            exp=DataArrayDouble(64) ; exp.iota(0.1)
            self.assertTrue(att.getFieldOn(amr,"YY").isEqualWithoutConsideringStr(exp,1e-12))
            exp0=DataArrayDouble(len1) ; exp0.iota(0.2)
            self.assertTrue(att.getFieldOn(amr[0].getMesh(),"YY").isEqualWithoutConsideringStr(exp0,1e-12))
            exp1=DataArrayDouble(len2) ; exp1.iota(0.3)
            self.assertTrue(att.getFieldOn(amr[1].getMesh(),"YY").isEqualWithoutConsideringStr(exp1,1e-12))
            pass
        pass

    def testSwig2AMR14(self):
        """ non regression linked to VTHB write."""
        fact=[2,2] ; fact2=[3,3]
        amr=MEDCouplingCartesianAMRMesh("mesh",2,[5,5],[0.,0.],[1.,1.])
        amr.addPatch([(1,3),(0,2)],fact)
        amr.addPatch([(1,3),(3,4)],fact)
        amr[0].addPatch([(1,3),(1,3)],fact2)
        amr[1].addPatch([(1,3),(1,2)],fact2)
        att=MEDCouplingAMRAttribute(amr,[("YY",1)],2)
        att.alloc()
        att.getFieldOn(amr,"YY").iota(0.1)
        att.getFieldOn(amr[0].getMesh(),"YY").iota(0.2)
        att.getFieldOn(amr[1].getMesh(),"YY").iota(0.3)
        att.getFieldOn(amr[0][0].getMesh(),"YY").iota(0.4)
        att.getFieldOn(amr[1][0].getMesh(),"YY").iota(0.5)
        self.assertEqual(amr[0].getBLTRRangeRelativeToGF(),[(2,6),(0,4)])
        self.assertEqual(amr[1].getBLTRRangeRelativeToGF(),[(2,6),(6,8)])
        self.assertEqual(amr[0][0].getBLTRRangeRelativeToGF(),[(9,15),(3,9)])
        self.assertEqual(amr[1][0].getBLTRRangeRelativeToGF(),[(9,15),(21,24)])
        pass

    def testOrderConsecutiveCells1D1(self):
        """A line in several unconnected pieces:"""
        m2 = MEDCouplingUMesh.New("bla", 1)
        c = DataArrayInt([NORM_SEG2,0,1,NORM_SEG3,1,3,2, NORM_SEG2,3,4,
                               NORM_SEG3,5,7,6, NORM_SEG3,7,9,8, NORM_SEG2,9,10,
                               NORM_SEG2,11,12,NORM_SEG2,12,13,
                               NORM_SEG2,14,15])
        cI = DataArrayInt([0,3,7,10,14,18,21,24,27,30])
        coords2 = DataArrayDouble([float(i) for i in range(32)], 16, 2)
        m2.setCoords(coords2);
        m2.setConnectivity(c, cI);
        m2.checkConsistency(1.0e-8);

        # Shuffle a bit :-)
        m2.renumberCells(DataArrayInt([0,3,6,8,1,4,7,5,2]), True);
        res = m2.orderConsecutiveCells1D()
        expRes = [0,3,6,8,1,4,2,7,5]
        self.assertEqual(m2.getNumberOfCells(),res.getNumberOfTuples())
        self.assertEqual(expRes, res.getValues())

        # A closed line (should also work)
        m3 = MEDCouplingUMesh.New("bla3", 1)
        conn3A = DataArrayInt([NORM_SEG2,0,1,NORM_SEG3,1,3,2, NORM_SEG2,3,0])
        coord3 = coords2[0:5]
        c.reAlloc(10)
        cI.reAlloc(4)

        m3.setCoords(coord3)
        m3.setConnectivity(conn3A, cI)
        m3.checkConsistency(1.0e-8)
        res2 = m3.orderConsecutiveCells1D()
        expRes2 = [0,1,2]
        self.assertEqual(m3.getNumberOfCells(),res2.getNumberOfTuples())
        self.assertEqual(expRes2, res2.getValues())
        pass

    def testDADApplyFuncOnThis1(self):
        d=DataArrayDouble(5) ; d.iota(0.)
        d.applyFuncOnThis("2*x+1")
        self.assertTrue(d.isEqual(DataArrayDouble([1.,3.,5.,7.,9.]),1e-12))
        d=DataArrayDouble(6) ; d.iota(0.) ; d.rearrange(2)
        d.applyFuncOnThis("2*x+1")
        self.assertTrue(d.isEqual(DataArrayDouble([1.,3.,5.,7.,9.,11.],3,2),1e-12))
        d.applyFuncOnThis("1+2*3")
        self.assertTrue(d.isEqual(DataArrayDouble([(7.,7.),(7.,7.),(7.,7.)]),1e-12))
        pass

    def testSwig2PointSetComputeFetchedNodeIds1(self):
        arr=DataArrayDouble(6) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr,arr)
        m=m.buildUnstructured()
        m0=m[[0,1,5,6,25,26,30,31,124]]
        ref=DataArrayInt([0,1,2,6,7,8,12,13,14,36,37,38,42,43,44,48,49,50,72,73,74,78,79,80,84,85,86,172,173,178,179,208,209,214,215])
        self.assertTrue(m0.computeFetchedNodeIds().isEqual(ref))
        self.assertTrue(MEDCoupling1SGTUMesh(m0).computeFetchedNodeIds().isEqual(ref))
        self.assertEqual(m0.getAllGeoTypes(),[NORM_HEXA8])
        m0.convertAllToPoly()
        self.assertEqual(m0.getAllGeoTypes(),[NORM_POLYHED])
        self.assertTrue(MEDCoupling1DGTUMesh(m0).computeFetchedNodeIds().isEqual(ref))
        pass

    def testSwig2PartDefinition1(self):
        pd=PartDefinition.New(5,22,3)
        self.assertTrue(isinstance(pd,SlicePartDefinition))
        self.assertTrue(pd.toDAI().isEqual(DataArrayInt([5,8,11,14,17,20])))
        self.assertEqual(pd.getNumberOfElems(),6)
        self.assertEqual(pd.getEffectiveStop(),23)
        pd=PartDefinition.New(5,23,3)
        self.assertTrue(isinstance(pd,SlicePartDefinition))
        self.assertTrue(pd.toDAI().isEqual(DataArrayInt([5,8,11,14,17,20])))
        self.assertEqual(pd.getNumberOfElems(),6)
        self.assertEqual(pd.getEffectiveStop(),23)
        self.assertEqual(pd.getSlice(),slice(5,23,3))
        pd=PartDefinition.New(5,22,1)
        self.assertTrue(isinstance(pd,SlicePartDefinition))
        self.assertTrue(pd.toDAI().isEqual(DataArrayInt([5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21])))
        self.assertEqual(pd.getNumberOfElems(),17)
        self.assertEqual(pd.getEffectiveStop(),22)
        pd=PartDefinition.New(5,23,3)+PartDefinition.New(23,27,3)
        self.assertTrue(isinstance(pd,SlicePartDefinition))
        self.assertEqual(pd.getNumberOfElems(),8)
        self.assertTrue(pd.toDAI().isEqual(DataArrayInt([5,8,11,14,17,20,23,26])))
        self.assertEqual(pd.getEffectiveStop(),29)
        pd=SlicePartDefinition(5,22,1)
        self.assertTrue(isinstance(pd,SlicePartDefinition))
        self.assertTrue(pd.toDAI().isEqual(DataArrayInt([5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21])))
        self.assertEqual(pd.getNumberOfElems(),17)
        self.assertEqual(pd.getEffectiveStop(),22)
        d=DataArrayInt([2,4,5,6,10])
        pd=PartDefinition.New(d)
        self.assertTrue(isinstance(pd,DataArrayPartDefinition))
        self.assertEqual(pd.toDAI().getHiddenCppPointer(),d.getHiddenCppPointer())
        pd=DataArrayPartDefinition(d)
        self.assertEqual(pd.toDAI().getHiddenCppPointer(),d.getHiddenCppPointer())
        pd=DataArrayPartDefinition(d)+DataArrayPartDefinition(DataArrayInt([12,14,20]))
        self.assertTrue(isinstance(pd,DataArrayPartDefinition))
        self.assertEqual(pd.getNumberOfElems(),8)
        self.assertTrue(pd.toDAI().isEqual(DataArrayInt([2,4,5,6,10,12,14,20])))
        pass

    def testSwig2SortEachPairToMakeALinkedList1(self):
        d=DataArrayInt([(50,49),(50,51),(51,52),(53,52),(53,54),(55,54),(55,56),(56,57),(58,57),(58,59),(60,59),(60,61),(61,62),(63,62),(63,64),(65,64),(65,66),(66,67)])
        d.sortEachPairToMakeALinkedList()
        self.assertTrue(d.isEqual(DataArrayInt([(49,50),(50,51),(51,52),(52,53),(53,54),(54,55),(55,56),(56,57),(57,58),(58,59),(59,60),(60,61),(61,62),(62,63),(63,64),(64,65),(65,66),(66,67)])))
        d=DataArrayInt([(0,2),(1,2),(1,3)])
        d.sortEachPairToMakeALinkedList()
        self.assertTrue(d.isEqual(DataArrayInt([(0,2),(2,1),(1,3)])))
        d=DataArrayInt([(0,2),(1,2),(3,1)])
        d.sortEachPairToMakeALinkedList()
        self.assertTrue(d.isEqual(DataArrayInt([(0,2),(2,1),(1,3)])))
        d=DataArrayInt([(8,6062),(6062,472),(472,6292),(6292,960)])
        d.sortEachPairToMakeALinkedList()
        self.assertTrue(d.isEqual(DataArrayInt([(8,6062),(6062,472),(472,6292),(6292,960)])))
        pass

    def testSwig2DAIIsRange(self):
        d=DataArrayInt([2,6,10])
        a,b=d.isRange()
        self.assertTrue(a)
        self.assertEqual(b,slice(2,11,4))
        self.assertTrue(DataArrayInt.Range(b.start,b.stop,b.step).isEqual(d))
        #
        d=DataArrayInt([2,7,10])
        a,b=d.isRange()
        self.assertTrue(not a)
        self.assertTrue(b is None)
        #
        d=DataArrayInt([22,17,12])
        a,b=d.isRange()
        self.assertTrue(a)
        self.assertEqual(b,slice(22,11,-5))
        self.assertTrue(DataArrayInt.Range(b.start,b.stop,b.step).isEqual(d))
        #
        d=DataArrayInt([22,16,12])
        a,b=d.isRange()
        self.assertTrue(not a)
        self.assertTrue(b is None)
        #
        d=DataArrayInt([33])
        a,b=d.isRange()
        self.assertTrue(a)
        self.assertEqual(b,slice(33,34,1))
        self.assertTrue(DataArrayInt.Range(b.start,b.stop,b.step).isEqual(d))
        #
        d=DataArrayInt([])
        a,b=d.isRange()
        self.assertTrue(a)
        self.assertEqual(b,slice(0,0,1))
        self.assertTrue(DataArrayInt.Range(b.start,b.stop,b.step).isEqual(d))
        #
        d=DataArrayInt([2,6,10,2])
        a,b=d.isRange()
        self.assertTrue(not a)
        self.assertTrue(b is None)
        pass

    def testSwig2PartDefinitionComposeWith1(self):
        f=PartDefinition.New(DataArrayInt([0,1,2,3,6,7,8,9]))
        g=PartDefinition.New(4,14,1)
        g2=g.deepCopy()
        self.assertTrue(g2.isEqual(g)[0])
        h=f.composeWith(g)
        self.assertTrue(isinstance(h,DataArrayPartDefinition))
        self.assertTrue(h.toDAI().isEqual(DataArrayInt([4,5,6,7,10,11,12,13])))
        f2=f.tryToSimplify()
        g2=g.tryToSimplify()
        self.assertEqual(f2.getHiddenCppPointer(),f.getHiddenCppPointer())# same because no simplification due to content of array
        self.assertEqual(g2.getHiddenCppPointer(),g.getHiddenCppPointer())# same because no simplification linked to type of PartDef
        p=PartDefinition.New(DataArrayInt([2,6,10]))
        p2=p.tryToSimplify()
        self.assertNotEqual(p2.getHiddenCppPointer(),p.getHiddenCppPointer())
        self.assertTrue(isinstance(p2,SlicePartDefinition))
        self.assertEqual(p2.getSlice(),slice(2,11,4))
        self.assertTrue(p2.isEqual(SlicePartDefinition(2,11,4))[0])
        self.assertTrue(p2.isEqual(p2.deepCopy())[0])
        self.assertTrue(not p2.isEqual(SlicePartDefinition(1,11,4))[0])
        self.assertTrue(not p2.isEqual(SlicePartDefinition(2,10,4))[0])
        self.assertTrue(not p2.isEqual(SlicePartDefinition(2,11,3))[0])
        pass

    def testSwig2DAIGetIdsStrictlyNegative1(self):
        d=DataArrayInt([4,-5,-1,0,3,99,-7])
        self.assertTrue(d.findIdsStrictlyNegative().isEqual(DataArrayInt([1,2,6])))
        pass

    def testSwig2DAIReplaceOneValByInThis1(self):
        d=DataArrayInt([4,-5,-1,0,-5,99,-7,5])
        d.changeValue(-5,900)
        self.assertTrue(d.isEqual(DataArrayInt([4,900,-1,0,900,99,-7,5])))
        pass

    def testSwig2DAIGetMinMaxValues1(self):
        d=DataArrayInt([4,-5,-1,0,3,99,-7])
        a,b=d.getMinMaxValues()
        self.assertEqual(a,-7)
        self.assertEqual(b,99)
        pass

    def testSwig2DAIBuildUniqueNotSorted1(self):
        d=DataArrayInt([-5,3,2,-1,2,3,-6,4,2,-5,3,7])
        self.assertTrue(d.buildUniqueNotSorted().isEqual(DataArrayInt([-5,3,2,-1,-6,4,7])))
        pass

    def testSwig2UMeshChangeOrientationOfCells1(self):
        """ Here testing changeOrientationOfCell method on unstructured meshes lying on no coords."""
        m=MEDCouplingUMesh("mesh",1)
        c=DataArrayInt([NORM_SEG2,4,5,NORM_SEG2,10,8,NORM_SEG3,20,7,33,NORM_SEG3,13,15,12,NORM_SEG2,3,2,NORM_SEG4,5,6,8,10,NORM_SEG4,34,33,3,2])
        cI=DataArrayInt([0,3,6,10,14,17,22,27])
        m.setConnectivity(c,cI)
        m.changeOrientationOfCells()
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([NORM_SEG2,5,4,NORM_SEG2,8,10,NORM_SEG3,7,20,33,NORM_SEG3,15,13,12,NORM_SEG2,2,3,NORM_SEG4,6,5,10,8,NORM_SEG4,33,34,2,3])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(cI))
        # testing 2D cells
        m=MEDCouplingUMesh("mesh",2)
        c=DataArrayInt([NORM_TRI3,0,1,2,NORM_QUAD4,3,4,5,6,NORM_POLYGON,7,8,9,10,11,NORM_TRI6,12,13,14,15,16,17,NORM_QUAD8,18,19,20,21,22,23,24,25,NORM_QPOLYG,26,27,28,29,30,31,32,33,34,35])
        cI=DataArrayInt([0,4,9,15,22,31,42])
        m.setConnectivity(c,cI)
        m.changeOrientationOfCells()
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([NORM_TRI3,0,2,1,NORM_QUAD4,3,6,5,4,NORM_POLYGON,7,11,10,9,8,NORM_TRI6,12,14,13,17,16,15,NORM_QUAD8,18,21,20,19,25,24,23,22,NORM_QPOLYG,26,30,29,28,27,35,34,33,32,31])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(cI))
        pass

    def testSwig2StructuredMeshCellLocation1(self):
        # 3D
        arrX=DataArrayDouble(5) ; arrX.iota()
        arrY=DataArrayDouble(4) ; arrY.iota()
        arrZ=DataArrayDouble(3) ; arrZ.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY,arrZ)
        li=[]
        liExp3D=[(0,0,0),(1,0,0),(2,0,0),(3,0,0),(0,1,0),(1,1,0),(2,1,0),(3,1,0),(0,2,0),(1,2,0),(2,2,0),(3,2,0),(0,0,1),(1,0,1),(2,0,1),(3,0,1),(0,1,1),(1,1,1),(2,1,1),(3,1,1),(0,2,1),(1,2,1),(2,2,1),(3,2,1)]
        self.assertEqual(24,m.getNumberOfCells())
        for i in range(m.getNumberOfCells()):
            li.append(m.getLocationFromCellId(i))
            pass
        self.assertEqual(liExp3D,li)
        self.assertRaises(InterpKernelException,m.getLocationFromCellId,24)
        self.assertRaises(InterpKernelException,m.getLocationFromCellId,-1)
        # 2D
        arrX=DataArrayDouble(5) ; arrX.iota()
        arrY=DataArrayDouble(4) ; arrY.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY)
        li=[]
        liExp2D=[(0,0),(1,0),(2,0),(3,0),(0,1),(1,1),(2,1),(3,1),(0,2),(1,2),(2,2),(3,2)]
        self.assertEqual(12,m.getNumberOfCells())
        for i in range(m.getNumberOfCells()):
            li.append(m.getLocationFromCellId(i))
            pass
        self.assertEqual(liExp2D,li)
        self.assertRaises(InterpKernelException,m.getLocationFromCellId,12)
        self.assertRaises(InterpKernelException,m.getLocationFromCellId,-1)
        # 1D
        arrX=DataArrayDouble(5) ; arrX.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arrX)
        self.assertEqual(4,m.getNumberOfCells())
        for i in range(m.getNumberOfCells()):
            self.assertEqual((i,),m.getLocationFromCellId(i))
            pass
        self.assertRaises(InterpKernelException,m.getLocationFromCellId,4)
        self.assertRaises(InterpKernelException,m.getLocationFromCellId,-1)
        pass

    def testSwig2StructuredMeshNodeLocation1(self):
        # 3D
        arrX=DataArrayDouble(5) ; arrX.iota()
        arrY=DataArrayDouble(4) ; arrY.iota()
        arrZ=DataArrayDouble(3) ; arrZ.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY,arrZ)
        li=[]
        liExp3D=[(0,0,0),(1,0,0),(2,0,0),(3,0,0),(4,0,0),(0,1,0),(1,1,0),(2,1,0),(3,1,0),(4,1,0),(0,2,0),(1,2,0),(2,2,0),(3,2,0),(4,2,0),(0,3,0),(1,3,0),(2,3,0),(3,3,0),(4,3,0),(0,0,1),(1,0,1),(2,0,1),(3,0,1),(4,0,1),(0,1,1),(1,1,1),(2,1,1),(3,1,1),(4,1,1),(0,2,1),(1,2,1),(2,2,1),(3,2,1),(4,2,1),(0,3,1),(1,3,1),(2,3,1),(3,3,1),(4,3,1),(0,0,2),(1,0,2),(2,0,2),(3,0,2),(4,0,2),(0,1,2),(1,1,2),(2,1,2),(3,1,2),(4,1,2),(0,2,2),(1,2,2),(2,2,2),(3,2,2),(4,2,2),(0,3,2),(1,3,2),(2,3,2),(3,3,2),(4,3,2)]
        self.assertEqual(60,m.getNumberOfNodes())
        for i in range(m.getNumberOfNodes()):
            li.append(m.getLocationFromNodeId(i))
            pass
        self.assertEqual(liExp3D,li)
        self.assertRaises(InterpKernelException,m.getLocationFromNodeId,60)
        self.assertRaises(InterpKernelException,m.getLocationFromNodeId,-1)
        # 2D
        arrX=DataArrayDouble(5) ; arrX.iota()
        arrY=DataArrayDouble(4) ; arrY.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY)
        li=[]
        liExp2D=[(0,0),(1,0),(2,0),(3,0),(4,0),(0,1),(1,1),(2,1),(3,1),(4,1),(0,2),(1,2),(2,2),(3,2),(4,2),(0,3),(1,3),(2,3),(3,3),(4,3)]
        self.assertEqual(20,m.getNumberOfNodes())
        for i in range(m.getNumberOfNodes()):
            li.append(m.getLocationFromNodeId(i))
            pass
        self.assertEqual(liExp2D,li)
        self.assertRaises(InterpKernelException,m.getLocationFromNodeId,20)
        self.assertRaises(InterpKernelException,m.getLocationFromNodeId,-1)
        # 1D
        arrX=DataArrayDouble(5) ; arrX.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arrX)
        self.assertEqual(5,m.getNumberOfNodes())
        for i in range(m.getNumberOfNodes()):
            self.assertEqual((i,),m.getLocationFromNodeId(i))
            pass
        self.assertRaises(InterpKernelException,m.getLocationFromCellId,5)
        self.assertRaises(InterpKernelException,m.getLocationFromCellId,-1)
        pass

    def testSwig2DataArrayPrintNotTooLong1(self):
        """ Now that DataArrayDouble and DataArrayInt and pickelized they can appear in YACS ports. Avoid to have too heavy string representation of them."""
        d=DataArrayDouble(2000) ; d.iota() ; d.rearrange(2)
        st0=d.repr() ; st1=str(d) ; st2=d.reprNotTooLong()
        self.assertEqual(st0,st1) # 1000 tuples ( >=0 and <= 1000) -> str(d)==d.repr()
        self.assertEqual(st1,st2)
        #
        d=DataArrayDouble(2002) ; d.iota() ; d.rearrange(2)
        st0=d.repr() ; st1=str(d) ; st2=d.reprNotTooLong()
        self.assertNotEqual(st0,st1) # 1001 tuples ( > 1000) -> str(d)==d.reprNotTooLong()
        self.assertEqual(st1,st2)
        self.assertIn(len(st2), list(range(0, 1000)))  # no more than 1000 characters
        ## Now for DataArrayInt
        d=DataArrayInt(2000) ; d.iota() ; d.rearrange(2)
        st0=d.repr() ; st1=str(d) ; st2=d.reprNotTooLong()
        self.assertEqual(st0,st1) # 1000 tuples ( >=0 and <= 1000) -> str(d)==d.repr()
        self.assertEqual(st1,st2)
        #
        d=DataArrayInt(2002) ; d.iota() ; d.rearrange(2)
        st0=d.repr() ; st1=str(d) ; st2=d.reprNotTooLong()
        self.assertNotEqual(st0,st1) # 1001 tuples ( > 1000) -> str(d)==d.reprNotTooLong()
        self.assertEqual(st1,st2)
        self.assertIn(len(st2), list(range(0, 1000)))  # no more than 1000 characters
        pass

    def testExtrudedMeshWithoutZipCoords1(self):
        """This test checks that MEDCouplingUMesh.buildExtrudedMesh do not perform a zipCoords."""
        arr=DataArrayDouble([(0.,0.),(1.,0.),(2.,0.),(3.,0.)])
        m=MEDCouplingUMesh("mesh",1) ; m.setCoords(arr)
        m.allocateCells()
        m.insertNextCell(NORM_SEG2,[1,2])
        arr1D=DataArrayDouble([(0.,0.),(0.,1.5),(0.,2.)])
        m1D=MEDCouplingUMesh("mesh1D",1) ; m1D.setCoords(arr1D)
        m1D.allocateCells()
        m1D.insertNextCell(NORM_SEG2,[0,1])
        m1D.insertNextCell(NORM_SEG2,[1,2])
        m2D=m.buildExtrudedMesh(m1D,0)
        self.assertEqual(m.getCoords().getHiddenCppPointer(),m2D.getCoords().getHiddenCppPointer())
        coo=DataArrayDouble([(0,0),(1,0),(2,0),(3,0),(0,1.5),(1,1.5),(2,1.5),(3,1.5),(0,2),(1,2),(2,2),(3,2)])
        self.assertTrue(m.getCoords().isEqual(coo,1e-12))
        self.assertTrue(m2D.getNodalConnectivity().isEqual(DataArrayInt([4,1,2,6,5,4,5,6,10,9])))
        self.assertTrue(m2D.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,10])))
        pass

    def testPointSetAreAllNodesFetched1(self):
        m=MEDCouplingCMesh() ; arr=DataArrayDouble(10) ; arr.iota()
        m.setCoords(arr,arr)
        m=m.buildUnstructured()
        self.assertTrue(m.areAllNodesFetched())
        m2=m[[0,2,3,4,5]]
        self.assertTrue(not m2.areAllNodesFetched())
        m2.zipCoords()
        self.assertTrue(m2.areAllNodesFetched())
        pass

    def testMEDCouplingPointSetComputeDiameterField1(self):
        arrX=DataArrayDouble([0.,1.1,1.7,2.1])
        arrY=DataArrayDouble([0.,0.7,0.8,1.9])
        arrZ=DataArrayDouble([0.,1.3,2.1,2.4])
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY,arrZ) ; m=m.buildUnstructured()
        f=m.computeDiameterField()
        f.checkConsistencyLight()
        exp=DataArrayDouble([1.8411952639521971,1.5937377450509227,1.5297058540778357,1.705872210923198,1.4352700094407325,1.3638181696985856,2.0273134932713295,1.8055470085267789,1.7492855684535902,1.5297058540778357,1.2206555615733703,1.1357816691600546,1.3638181696985856,1.004987562112089,0.9,1.7492855684535902,1.4866068747318506,1.4177446878757824,1.3379088160259651,0.9695359714832656,0.8602325267042626,1.1445523142259597,0.6782329983125266,0.5099019513592785,1.5842979517754858,1.2884098726725124,1.208304597359457])
        self.assertTrue(exp.isEqual(f.getArray(),1e-12))
        m1=m[::2]
        m2=m[1::2]
        m2.simplexize(PLANAR_FACE_5)
        m3=MEDCouplingUMesh.MergeUMeshesOnSameCoords(m1,m2)
        f=m3.computeDiameterField()
        f.checkConsistencyLight()
        exp2=DataArrayDouble([1.8411952639521971,1.5297058540778357,1.4352700094407325,2.0273134932713295,1.7492855684535902,1.2206555615733703,1.3638181696985856,0.9,1.4866068747318506,1.3379088160259651,0.8602325267042626,0.6782329983125266,1.5842979517754858,1.208304597359457,1.47648230602334,1.47648230602334,1.47648230602334,1.47648230602334,1.47648230602334,1.7029386365926402,1.7029386365926402,1.7029386365926402,1.7029386365926402,1.7029386365926402,1.3601470508735445,1.3601470508735445,1.3601470508735445,1.3601470508735445,1.3601470508735445,1.70293863659264,1.70293863659264,1.70293863659264,1.70293863659264,1.70293863659264,1.3601470508735445,1.3601470508735445,1.3601470508735445,1.3601470508735445,1.3601470508735445,1.063014581273465,1.063014581273465,1.063014581273465,1.063014581273465,1.063014581273465,1.0,1.0,1.0,1.0,1.0,1.5556349186104046,1.5556349186104046,1.5556349186104046,1.5556349186104046,1.5556349186104046,1.3601470508735443,1.3601470508735443,1.3601470508735443,1.3601470508735443,1.3601470508735443,0.9219544457292886,0.9219544457292886,0.9219544457292886,0.9219544457292886,0.9219544457292886,1.140175425099138,1.140175425099138,1.140175425099138,1.140175425099138,1.140175425099138,0.5,0.5,0.5,0.5,0.5,1.2529964086141667,1.2529964086141667,1.2529964086141667,1.2529964086141667,1.2529964086141667])
        self.assertTrue(exp2.isEqual(f.getArray(),1e-12))
        # TRI3 - spacedim = 2
        coo=DataArrayDouble([(1,1),(5,1.9),(2.1,3)])
        m=MEDCoupling1SGTUMesh("mesh",NORM_TRI3) ; m.setCoords(coo)
        for c in [[0,1,2],[0,2,1],[2,1,0]]:
            m.setNodalConnectivity(DataArrayInt(c))
            self.assertAlmostEqual(m.computeDiameterField().getArray()[0],4.1,12)
            m2=m.buildUnstructured() ; m2.convertLinearCellsToQuadratic(0)
            self.assertAlmostEqual(m2.computeDiameterField().getArray()[0],4.1,12)
            m3=m.buildUnstructured() ; m3.convertLinearCellsToQuadratic(1)
            self.assertAlmostEqual(m3.computeDiameterField().getArray()[0],4.1,12)
        # TRI3 - spacedim = 3
        coo=DataArrayDouble([(1.3198537928820775,1.0991902391274959,-0.028645697595823361),(5.2486835106806335,2.2234012799688281,0.30368935050077939),(2.2973688139447361,3.1572023778066649,0.10937756365410012)])
        m=MEDCoupling1SGTUMesh("mesh",NORM_TRI3) ; m.setCoords(coo)
        for c in [[0,1,2],[0,2,1],[2,1,0]]:
            m.setNodalConnectivity(DataArrayInt(c))
            self.assertAlmostEqual(m.computeDiameterField().getArray()[0],4.1,12)
            m2=m.buildUnstructured() ; m2.convertLinearCellsToQuadratic(0)
            self.assertAlmostEqual(m2.computeDiameterField().getArray()[0],4.1,12)
            m3=m.buildUnstructured() ; m3.convertLinearCellsToQuadratic(1)
            self.assertAlmostEqual(m3.computeDiameterField().getArray()[0],4.1,12)
        # QUAD4 - spacedim = 2
        coo=DataArrayDouble([(0,2),(2,0),(6,4),(4,9)])
        m=MEDCoupling1SGTUMesh("mesh",NORM_QUAD4) ; m.setCoords(coo)
        exp3=sqrt(85.)
        for delta in range(4):
            c = [(elt + delta) % 4 for elt in range(4)]
            m.setNodalConnectivity(DataArrayInt(c))
            self.assertAlmostEqual(m.computeDiameterField().getArray()[0],exp3,12)
            m2=m.buildUnstructured() ; m2.convertLinearCellsToQuadratic(0)
            self.assertAlmostEqual(m2.computeDiameterField().getArray()[0],exp3,12)
            m3=m.buildUnstructured() ; m3.convertLinearCellsToQuadratic(1)
            self.assertAlmostEqual(m3.computeDiameterField().getArray()[0],exp3,12)
            c.reverse()
            m.setNodalConnectivity(DataArrayInt(c))
            self.assertAlmostEqual(m.computeDiameterField().getArray()[0],exp3,12)
            m2=m.buildUnstructured() ; m2.convertLinearCellsToQuadratic(0)
            self.assertAlmostEqual(m2.computeDiameterField().getArray()[0],exp3,12)
            m3=m.buildUnstructured() ; m3.convertLinearCellsToQuadratic(1)
            self.assertAlmostEqual(m3.computeDiameterField().getArray()[0],exp3,12)
        # QUAD4 - spacedim = 3
        coo=DataArrayDouble([(0.26570992384234871,2.0405889913271817,-0.079134238105786903),(2.3739976619218064,0.15779148692781009,0.021842842914139737),(6.1207841448393197,4.3755532938679655,0.43666375769970678),(3.8363255342943359,9.2521096041694229,0.41551170895942313)])
        m=MEDCoupling1SGTUMesh("mesh",NORM_QUAD4) ; m.setCoords(coo)
        for delta in range(4):
            c = [(elt + delta) % 4 for elt in range(4)]
            m.setNodalConnectivity(DataArrayInt(c))
            self.assertAlmostEqual(m.computeDiameterField().getArray()[0],exp3,12)
            m2=m.buildUnstructured() ; m2.convertLinearCellsToQuadratic(0)
            self.assertAlmostEqual(m2.computeDiameterField().getArray()[0],exp3,12)
            m3=m.buildUnstructured() ; m3.convertLinearCellsToQuadratic(1)
            self.assertAlmostEqual(m3.computeDiameterField().getArray()[0],exp3,12)
            c.reverse()
            m.setNodalConnectivity(DataArrayInt(c))
            self.assertAlmostEqual(m.computeDiameterField().getArray()[0],exp3,12)
            m2=m.buildUnstructured() ; m2.convertLinearCellsToQuadratic(0)
            self.assertAlmostEqual(m2.computeDiameterField().getArray()[0],exp3,12)
            m3=m.buildUnstructured() ; m3.convertLinearCellsToQuadratic(1)
            self.assertAlmostEqual(m3.computeDiameterField().getArray()[0],exp3,12)
        # PENTA6
        # noise of coo=DataArrayDouble([(0,0,0),(1,0,0),(0,1,0),(0,0,2),(1,0,2),(0,1,2)]) + rotation([0.7,-1.2,0.6],[-4,-1,10],0.3)
        coo=DataArrayDouble([(-0.28594726851554486,-0.23715005500928255,-0.10268080010083136),(0.6167364988633947,-0.008923258436324799,-0.08574087516687756),(-0.6132873463333834,0.6943403970881654,-0.2806118260037991),(-0.40705974936532896,-0.05868487929989308,1.7724055544436323),(0.5505955507861958,0.19145393798144705,1.8788156352163994),(-0.6092686217773406,0.812502961290914,1.685712743757831)])
        m=MEDCoupling1SGTUMesh("mesh",NORM_PENTA6) ; m.setCoords(coo)
        exp4=2.5041256256889888
        self.assertAlmostEqual(exp4,coo.buildEuclidianDistanceDenseMatrix().getMaxValue()[0],12)# <- the definition of diameter
        for delta in range(3):
            c = [(elt + delta) % 3 for elt in range(3)]
            c+=[elt+3 for elt in c]
            m.setNodalConnectivity(DataArrayInt(c))
            self.assertAlmostEqual(m.computeDiameterField().getArray()[0],exp4,12)
            m2=m.buildUnstructured() ; m2.convertLinearCellsToQuadratic(0)
            self.assertAlmostEqual(m2.computeDiameterField().getArray()[0],exp4,12)
            c.reverse()
            m.setNodalConnectivity(DataArrayInt(c))
            self.assertAlmostEqual(m.computeDiameterField().getArray()[0],exp4,12)
            m2=m.buildUnstructured() ; m2.convertLinearCellsToQuadratic(0)
            self.assertAlmostEqual(m2.computeDiameterField().getArray()[0],exp4,12)
        # HEXA8
        # noise of coo=DataArrayDouble([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,2),(1,0,2),(1,1,2),(0,1,2)]) + rotation([0.7,-1.2,0.6],[-4,-1,10],0.3)
        coo=DataArrayDouble([(-0.21266406388867243,-0.3049569460042527,-0.11012394815006032),(0.7641037943272584,-0.06990814759929553,-0.0909613877456491),(0.47406560768559974,0.8681310650341907,-0.2577311403703061),(-0.5136830410871793,0.644390554940524,-0.21319015989794698),(-0.4080167737381202,-0.12853761670628505,1.7869166291979348),(0.5650318811550441,0.20476257733110748,1.8140158890821603),(0.3230844436386215,1.1660778242678538,1.7175073141333406),(-0.6656588358432984,0.918357550969698,1.7566470691880265)])
        m=MEDCoupling1SGTUMesh("mesh",NORM_HEXA8) ; m.setCoords(coo)
        exp5=2.5366409441884215
        self.assertAlmostEqual(exp5,coo.buildEuclidianDistanceDenseMatrix().getMaxValue()[0],12)# <- the definition of diameter
        for delta in range(4):
            c = [(elt + delta) % 4 for elt in range(4)]
            c+=[elt+4 for elt in c]
            m.setNodalConnectivity(DataArrayInt(c))
            self.assertAlmostEqual(m.computeDiameterField().getArray()[0],exp5,12)
            m2=m.buildUnstructured() ; m2.convertLinearCellsToQuadratic(0)
            self.assertAlmostEqual(m2.computeDiameterField().getArray()[0],exp5,12)
            c.reverse()
            m.setNodalConnectivity(DataArrayInt(c))
            self.assertAlmostEqual(m.computeDiameterField().getArray()[0],exp5,12)
            m2=m.buildUnstructured() ; m2.convertLinearCellsToQuadratic(0)
            self.assertAlmostEqual(m2.computeDiameterField().getArray()[0],exp5,12)
        # PYRA5 (1) 5th node is further
        # noise of coo=DataArrayDouble([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0.5,0.5,2)]) + rotation([0.7,-1.2,0.6],[-4,-1,10],0.3)
        coo=DataArrayDouble([(-0.31638393672228626,-0.3157865246451914,-0.12555467233075002),(0.7281379795666488,0.03836511217237115,-0.08431662762197323),(0.4757967840735147,0.8798897996143908,-0.2680890320119049),(-0.5386339871809047,0.5933159894201252,-0.2975311238319419),(0.012042592988768974,0.534282135495012,1.7859521682027926)])
        m=MEDCoupling1SGTUMesh("mesh",NORM_PYRA5) ; m.setCoords(coo)
        exp6=2.1558368027391386
        self.assertAlmostEqual(exp6,coo.buildEuclidianDistanceDenseMatrix().getMaxValue()[0],12)# <- the definition of diameter
        for delta in range(4):
            c = [(elt + delta) % 4 for elt in range(4)]
            c+=[4]
            m.setNodalConnectivity(DataArrayInt(c))
            self.assertAlmostEqual(m.computeDiameterField().getArray()[0],exp6,12)
            m2=m.buildUnstructured() ; m2.convertLinearCellsToQuadratic(0)
            self.assertAlmostEqual(m2.computeDiameterField().getArray()[0],exp6,12)
            pass
        # PYRA5 (2) 5th node is closer
        # noise of coo=DataArrayDouble([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0.5,0.5,0.1)]) + rotation([0.7,-1.2,0.6],[-4,-1,10],0.3)
        coo=DataArrayDouble([(-0.31638393672228626,-0.3157865246451914,-0.12555467233075002),(0.7281379795666488,0.03836511217237115,-0.08431662762197323),(0.4757967840735147,0.8798897996143908,-0.2680890320119049),(-0.5386339871809047,0.5933159894201252,-0.2975311238319419),(0.092964408350795,0.33389670321297005,-0.10171764888060142)])
        m=MEDCoupling1SGTUMesh("mesh",NORM_PYRA5) ; m.setCoords(coo)
        exp7=1.4413563787228953
        self.assertAlmostEqual(exp7,coo.buildEuclidianDistanceDenseMatrix().getMaxValue()[0],12)# <- the definition of diameter
        for delta in range(4):
            c = [(elt + delta) % 4 for elt in range(4)]
            c+=[4]
            m.setNodalConnectivity(DataArrayInt(c))
            self.assertAlmostEqual(m.computeDiameterField().getArray()[0],exp7,12)
            m2=m.buildUnstructured() ; m2.convertLinearCellsToQuadratic(0)
            self.assertAlmostEqual(m2.computeDiameterField().getArray()[0],exp7,12)
            pass
        # TETRA4
        # noise of coo=DataArrayDouble([(0,0,0),(1,0,0),(0,1,0),(1,1,1)]) + rotation([0.7,-1.2,0.6],[-4,-1,10],0.3)
        coo=DataArrayDouble([(-0.2256894071281369,-0.27631691290428106,-0.20266086543995965),(0.655458695100186,-0.08173323565551605,-0.19254662462061933),(-0.49893490718947264,0.5848097154568599,-0.3039928255382145),(0.2988102920828487,1.0582266398878504,0.7347375047372364)])
        m=MEDCoupling1SGTUMesh("mesh",NORM_TETRA4) ; m.setCoords(coo)
        exp8=1.7131322579364157
        self.assertAlmostEqual(exp8,coo.buildEuclidianDistanceDenseMatrix().getMaxValue()[0],12)# <- the definition of diameter
        for c in [[0,1,2,3],[0,3,2,1],[0,1,3,2],[0,2,3,1],[0,3,1,2],[0,2,1,3]]:
            for i in range(4):
                m.setNodalConnectivity(DataArrayInt([(elt+i)%4 for elt in c]))
                self.assertAlmostEqual(m.computeDiameterField().getArray()[0],exp8,12)
                m2=m.buildUnstructured() ; m2.convertLinearCellsToQuadratic(0)
                self.assertAlmostEqual(m2.computeDiameterField().getArray()[0],exp8,12)
                pass
            pass
        pass

    def testMEDCouplingSkyLineArray(self):
        index = DataArrayInt([ 0, 3, 5, 6, 6 ])
        value = DataArrayInt([ 1, 2, 3, 2, 3, 3 ])

        sla0 = MEDCouplingSkyLineArray()
        self.assertEqual( -1, sla0.getNumberOf() )
        self.assertEqual( 0,  sla0.getLength() )
        sla0.set( index, value )
        self.assertTrue( index.isEqual( sla0.getIndexArray() ))
        self.assertTrue( value.isEqual( sla0.getValuesArray() ))
        self.assertEqual( 4, sla0.getNumberOf() )
        self.assertEqual( 6, sla0.getLength() )

        sla1 = MEDCouplingSkyLineArray( index, value )
        self.assertTrue( index.isEqual( sla1.getIndexArray() ))
        self.assertTrue( value.isEqual( sla1.getValuesArray() ))
        self.assertEqual( 4, sla1.getNumberOf() )
        self.assertEqual( 6, sla1.getLength() )

        sla2 = MEDCouplingSkyLineArray( sla1 )
        self.assertTrue( index.isEqual( sla2.getIndexArray() ))
        self.assertTrue( value.isEqual( sla2.getValuesArray() ))
        self.assertEqual( 4, sla2.getNumberOf() )
        self.assertEqual( 6, sla2.getLength() )

        indexVec = ivec(); indexVec.reserve( len( index ))
        for i in index: indexVec.push_back( i[0] )
        valueVec = ivec(); valueVec.reserve( len( value ))
        for i in value: valueVec.push_back( i[0] )
        sla3 = MEDCouplingSkyLineArray( indexVec, valueVec )
        self.assertTrue( index.isEqual( sla3.getIndexArray() ))
        self.assertTrue( value.isEqual( sla3.getValuesArray() ))
        self.assertEqual( 4, sla3.getNumberOf() )
        self.assertEqual( 6, sla3.getLength() )

        pass

    def testMEDCouplingSkyLineArrayThreeLevels(self):
        #  [[28,1,4]] , [[2,35,8], [9,10,1,12]]
        superi = DataArrayInt([ 0,1,3 ])
        index = DataArrayInt ([ 0,3,6,10 ])
        value = DataArrayInt ([ 28,1,4,2,35,8,9,10,1,12 ])

        sla0 = MEDCouplingSkyLineArray()
        self.assertEqual( -1, sla0.getSuperNumberOf() )
        self.assertEqual( -1, sla0.getNumberOf() )
        self.assertEqual( 0,  sla0.getLength() )
        sla0.set3( superi.deepCopy(), index.deepCopy(), value.deepCopy() )
        self.assertTrue( superi.isEqual( sla0.getSuperIndexArray() ))

        pack = sla0.getSimplePackSafe(2)
        self.assertEqual([9,10,1,12], pack)
        ids = sla0.findPackIds([0,1], [9,10,1,12])
        self.assertEqual([-1,1], ids)

        sla0.deletePack(1, 1)
        si, idx, val = sla0.getSuperIndexArray(), sla0.getIndexArray(), sla0.getValuesArray()
        self.assertEqual([28,1,4,2,35,8], val.getValues())
        self.assertEqual([0,3,6], idx.getValues())
        self.assertEqual([0,1,2], si.getValues())

        sla0.pushBackPack(0, [3,2,1,0])
        si, idx, val = sla0.getSuperIndexArray(), sla0.getIndexArray(), sla0.getValuesArray()
        self.assertEqual([0,2,3], si.getValues())
        self.assertEqual([0,3,7,10], idx.getValues())
        self.assertEqual([28,1,4,3,2,1,0,  2,35,8], val.getValues())

        # Build connectivity from POLYHED connectivity
        cI = [0,16,41]
        c = [NORM_POLYHED, 1,2,3,-1,  2,3,4,-1,  3,4,5,-1,  4,5,6,
             NORM_POLYHED, 7,8,9,10,-1,  9,10,11,12,-1,  3,4,5,6,-1,  5,6,7,8,-1,  9,10,11,12]
        sla0 = MEDCouplingSkyLineArray.BuildFromPolyhedronConn(DataArrayInt(c), DataArrayInt(cI))
        si, idx, val = sla0.getSuperIndexArray(), sla0.getIndexArray(), sla0.getValuesArray()
        self.assertEqual([0,4,9], si.getValues())
        self.assertEqual([0,3,6,9,12,16,20,24,28,32], idx.getValues())
        self.assertEqual([1,2,3,  2,3,4,  3,4,5,  4,5,6,
                          7,8,9,10,   9,10,11,12,  3,4,5,6,  5,6,7,8,  9,10,11,12], val.getValues())
        c1, cI1 = sla0.convertToPolyhedronConn()
        self.assertEqual(c1.getValues(), c)
        self.assertEqual(cI1.getValues(), cI)
        pass

    def testMEDCouplingSkyLineArrayThreeLevels2(self):
        si = [0, 9, 15, 21]
        siRef = [0, 9, 16, 22]
        idx = [0,4,8,12,16,20,23,26,29,  32,36,40,44,48,52,  56,60,64,68,72,76,80]
        c = [1,0,2,3,  5,7,6,4,  1,5,4,0,  0,4,6,2,  2,6,7,3,  3,7,8,  7,5,8,  5,1,8,  1,3,8,
             9,1,3,10,  11,12,7,5,  9,11,5,1,  1,5,7,3,  3,7,12,10,  10,12,11,9,
             11,5,7,12,  14,16,15,13,  11,14,13,5,  5,13,15,7,  7,15,16,12,  12,16,14,11]
        idxRef = [0,4,8,12,16,20,23,26,29,32,36,40,44,48,52,55,58, 62, 66, 70, 74, 78, 82 ]
        cRef = [1,0,2,3,  5,7,6,4,  1,5,4,0,  0,4,6,2,  2,6,7,3,  3,7,8,  7,5,8,  5,1,8,  1,3,8,
             9,1,3,10,  11,12,7,5,  9,11,5,1,  3,7,12,10,  10,12,11,9,  3,7,8,  7,5,8,
             11,5,7,12,  14,16,15,13,  11,14,13,5,  5,13,15,7,  7,15,16,12,  12,16,14,11]
        sla0 = MEDCouplingSkyLineArray()
        sla0.set3( DataArrayInt(si), DataArrayInt(idx), DataArrayInt(c) )
        ids = sla0.findPackIds([1], [1,5,7,3])
        sla0.deletePack(1, ids[0])
        sla0.pushBackPack(1, [3,7,8])
        sla0.pushBackPack(1, [7,5,8])
        si, idx, val = sla0.getSuperIndexArray(), sla0.getIndexArray(), sla0.getValuesArray()
        self.assertEqual(siRef, si.getValues())
        self.assertEqual(idxRef, idx.getValues())
        self.assertEqual(cRef, val.getValues())

        idxRef2 = [0,4,8,12,16,20,23,26,29,32,36,40,42,46,50,53,56, 60, 64, 68, 72, 76, 80 ]
        cRef2 = [1,0,2,3,  5,7,6,4,  1,5,4,0,  0,4,6,2,  2,6,7,3,  3,7,8,  7,5,8,  5,1,8,  1,3,8,
             9,1,3,10,  11,12,7,5,  300,300,  3,7,12,10,  10,12,11,9,  3,7,8,  7,5,8,
             11,5,7,12,  14,16,15,13,  11,14,13,5,  5,13,15,7,  7,15,16,12,  12,16,14,11]
        sla0.replacePack(1,2, [300,300])
        si, idx, val = sla0.getSuperIndexArray(), sla0.getIndexArray(), sla0.getValuesArray()
        self.assertEqual(siRef, si.getValues())
        self.assertEqual(idxRef2, idx.getValues())
        self.assertEqual(cRef2, val.getValues())

        sla0.replacePack(1,2, [9,11,5,1])
        si, idx, val = sla0.getSuperIndexArray(), sla0.getIndexArray(), sla0.getValuesArray()
        self.assertEqual(siRef, si.getValues())
        self.assertEqual(idxRef, idx.getValues())
        self.assertEqual(cRef, val.getValues())

        sla0.replaceSimplePack(11, [300,300])  # 11 is the abs index of pack (superIdx=1,idx=2)
        si, idx, val = sla0.getSuperIndexArray(), sla0.getIndexArray(), sla0.getValuesArray()
        self.assertEqual(siRef, si.getValues())
        self.assertEqual(idxRef2, idx.getValues())
        self.assertEqual(cRef2, val.getValues())

        sla0.replaceSimplePack(11, [9,11,5,1])  # 11 is the abs index of pack (superIdx=1,idx=2)
        si, idx, val = sla0.getSuperIndexArray(), sla0.getIndexArray(), sla0.getValuesArray()
        self.assertEqual(siRef, si.getValues())
        self.assertEqual(idxRef, idx.getValues())
        self.assertEqual(cRef, val.getValues())
        pass

    def testMEDCouplingUMeshgenerateGraph(self):
        # cartesian mesh 3x3
        arr=DataArrayDouble(4) ; arr.iota()
        c=MEDCouplingCMesh() ; c.setCoords(arr,arr)
        m=c.buildUnstructured()
        graph = m.generateGraph()
        # 0 1 2
        # 3 4 5
        # 6 7 8
        valRef=[ 0,1,3,
                 0,1,2,4,
                 1,2,5,
                 0,3,4,6,
                 1,3,4,5,7,
                 2,4,5,8,
                 3,6,7,
                 4,6,7,8,
                 5,7,8]
        self.assertEqual(valRef,list(graph.getValuesArray().getValues()));

        indRef=[0, 3, 7, 10, 14, 19, 23, 26, 30, 33]
        self.assertEqual(indRef,list(graph.getIndexArray().getValues()));
        pass

    def testSwig2MEDCouplingCurveLinearReprQuick1(self):
        """Non regression test. Error in m.__str__ when m is a MEDCouplingCurveLinear with spaceDim != meshDim."""
        arr=DataArrayDouble(12) ; arr.iota() ; arr.rearrange(2)
        m=MEDCouplingCurveLinearMesh()
        m.setCoords(arr)
        m.setNodeGridStructure([3,2])
        m.checkConsistencyLight()
        self.assertEqual(m.getMeshDimension(),2)
        self.assertEqual(m.getSpaceDimension(),2)
        self.assertTrue(not "mismatch" in m.__str__())
        self.assertTrue(not "mismatch" in m.__repr__())
        #
        arr=DataArrayDouble(18) ; arr.iota() ; arr.rearrange(3)
        m.setCoords(arr)
        self.assertEqual(m.getMeshDimension(),2)
        self.assertEqual(m.getSpaceDimension(),3)
        self.assertTrue(not "mismatch" in m.__str__())
        self.assertTrue(not "mismatch" in m.__repr__())# bug was here !
        pass

    def testSwig2BugComputeOffsets1(self):
        """Non regression test. computeOffsetsFull on empty array must return 0."""
        d=DataArrayInt([3])
        d.computeOffsetsFull()
        self.assertTrue(d.isEqual(DataArrayInt([0,3])))
        d=DataArrayInt([])
        d.computeOffsets()
        self.assertTrue(d.isEqual(DataArrayInt([])))
        d=DataArrayInt([])
        d.computeOffsetsFull()
        self.assertTrue(d.isEqual(DataArrayInt([0]))) # <- bug was here
        pass

    def testSwig2Cartesianize1(self):
        """Test of engine of cartesianize mechanism in medcoupling"""
        # cyl 2D
        arr=DataArrayDouble([(3,0.2),(2,1.6)]) ; arr.setInfoOnComponents(["A","BB"])
        arr2=arr.cartesianize(AX_CYL)
        arr2_exp=DataArrayDouble([(2.940199733523725,0.5960079923851836),(-0.05839904460257763,1.9991472060830102)]) ; arr2_exp.setInfoOnComponents(["A","BB"])
        self.assertTrue(arr2_exp.isEqual(arr2,1e-14))
        # spher 2D
        arr3=arr.cartesianize(AX_SPHER)
        self.assertTrue(arr2_exp.isEqual(arr3,1e-14))
        # cyl 3D
        arr=DataArrayDouble([(3,0.2,7.1),(2,1.6,12.3)]) ; arr.setInfoOnComponents(["A","BB","CCC"])
        arr4=arr.cartesianize(AX_CYL)
        arr4_exp=DataArrayDouble([(2.940199733523725,0.5960079923851836,7.1),(-0.05839904460257763,1.9991472060830102,12.3)]) ; arr4_exp.setInfoOnComponents(["A","BB","CCC"])
        self.assertTrue(arr4_exp.isEqual(arr4,1e-14))
        # spher 3D
        arr=DataArrayDouble([(3,0.2,0.5),(2,1.3,5.8)]) ; arr.setInfoOnComponents(["A","BB","CCC"])
        arr5=arr.cartesianize(AX_SPHER)
        arr5_exp=DataArrayDouble([(0.5230462208645272,0.2857414527616764,2.940199733523725),(1.706499157790973,-0.8953424658735863,0.5349976572491747)]) ; arr5_exp.setInfoOnComponents(["A","BB","CCC"])
        self.assertTrue(arr5_exp.isEqual(arr5,1e-14))
        #
        m=MEDCouplingCMesh() ; m.setName("aa") ; m.setDescription("bbb") ; m.setTime(4.125,5,6) ; m.setTimeUnit("ms")
        arrX=DataArrayDouble([0,1,2]) ; arrX.setInfoOnComponent(0,"ccc")
        arrY=DataArrayDouble([3,4,5,6]) ; arrY.setInfoOnComponent(0,"dddd")
        m.setCoords(arrX,arrY)
        m2=m.buildCurveLinear()
        #
        self.assertTrue(isinstance(m2,MEDCouplingCurveLinearMesh))
        self.assertEqual(m2.getName(),"aa")
        self.assertEqual(m2.getDescription(),"bbb")
        self.assertEqual(m2.getTime(),[4.125,5,6])
        self.assertEqual(m2.getTimeUnit(),"ms")
        m2c_exp=DataArrayDouble([(0.,3.),(1.,3.),(2.,3.),(0.,4.),(1.,4.),(2.,4.),(0.,5.),(1.,5.),(2.,5.),(0.,6.),(1.,6.),(2.,6.)]) ; m2c_exp.setInfoOnComponents(["ccc","dddd"])
        self.assertTrue(m2.getCoords().isEqual(m2c_exp,1e-14))
        self.assertEqual(m2.getNodeGridStructure(),(3,4))
        pass

    def testRemoveIdsFromIndexedArrays1(self):
        arr=DataArrayInt([101,102,103,201,202,203,204,301,501,502,503,504,505,601,602])
        arrI=DataArrayInt([0,3,7,8,8,13,15])
        # case where all elts in inputs are in
        arr2=arr.deepCopy() ; arrI2=arrI.deepCopy()
        self.assertTrue(MEDCouplingUMesh.RemoveIdsFromIndexedArrays([501,502],arr2,arrI2))
        self.assertTrue(arr2.isEqual(DataArrayInt([101,102,103,201,202,203,204,301,503,504,505,601,602])))
        self.assertTrue(arrI2.isEqual(DataArrayInt([0,3,7,8,8,11,13])))
        # case where part of elts in inputs are in
        arr2=arr.deepCopy() ; arrI2=arrI.deepCopy()
        self.assertTrue(MEDCouplingUMesh.RemoveIdsFromIndexedArrays([504,507],arr2,arrI2))
        self.assertTrue(arr2.isEqual(DataArrayInt([101,102,103,201,202,203,204,301,501,502,503,505,601,602])))
        self.assertTrue(arrI2.isEqual(DataArrayInt([0,3,7,8,8,12,14])))
        # case where no elts in inputs are in
        arr2=arr.deepCopy() ; arrI2=arrI.deepCopy()
        self.assertTrue(not MEDCouplingUMesh.RemoveIdsFromIndexedArrays([1,5,701],arr2,arrI2))
        self.assertTrue(arr2.isEqual(arr))
        self.assertTrue(arrI2.isEqual(arrI))
        pass

    def testFieldIntIsOnStage1(self):
        """ My first test with field int."""
        m=MEDCouplingCMesh()
        m.setName("mesh")
        arrX=DataArrayDouble([0,1,2,3])
        m.setCoords(arrX,arrX)
        f=MEDCouplingFieldInt(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayInt(8) ; arr.iota() ;f.setArray(arr)
        self.assertRaises(InterpKernelException,f.checkConsistencyLight)
        arr=DataArrayInt(9) ; arr.iota() ;f.setArray(arr)
        f.checkConsistencyLight()
        f.setTimeUnit("ms")
        self.assertEqual(f.getTimeUnit(),"ms")
        f.setTime(3.2,5,6)
        a,b,c=f.getTime()
        self.assertEqual(b,5)
        self.assertEqual(c,6)
        self.assertEqual(a,3.2,12)
        pass

    def testNoThrowOn1DGTU2UOnNullCells(self):
        """ Non regression test : no throw when trying to convert 1DGTUMesh to UMesh on an empty mesh"""
        m=MEDCoupling1DGTUMesh("",NORM_POLYGON) ; m.setCoords(DataArrayDouble([],0,3))
        m.setNodalConnectivity(DataArrayInt([]),DataArrayInt([0]))
        m=m.buildUnstructured()
        pass

    def testExplodeMeshIntoMicroEdges1(self):
        """ test for new functionality MEDCouplingUMesh.explodeMeshIntoMicroEdges"""
        m=MEDCouplingUMesh("mesh",2)
        coo=DataArrayDouble([2,0,10,0,12,0,0,3,4,5,10,5,12,7,3,2.5,7,2.5,6,0,10,2.5,11,2.5,11,0,7,5],14,2)
        m.setCoords(coo)
        m.allocateCells()
        # here a mix of quadratic, linear cells. Non conform but conform considering micro edges
        m.insertNextCell(NORM_TRI6,[0,4,1,7,8,9])
        m.insertNextCell(NORM_TRI6,[1,5,2,10,11,12])
        m.insertNextCell(NORM_TRI6,[5,1,4,10,8,13])
        m.insertNextCell(NORM_TRI3,[3,4,7])
        m.insertNextCell(NORM_TRI3,[3,7,0])
        m.insertNextCell(NORM_TRI3,[6,2,11])
        m.insertNextCell(NORM_TRI3,[6,11,5])
        m.insertNextCell(NORM_TRI3,[6,5,13])
        m.insertNextCell(NORM_TRI3,[6,13,4])
        edges,d,di,rd,rdi=m.explodeMeshIntoMicroEdges() # <- new method
        self.assertTrue(MEDCoupling1SGTUMesh(edges).getNodalConnectivity().isEqual(DataArrayInt([0,7,7,4,4,8,8,1,1,9,9,0,1,10,10,5,5,11,11,2,2,12,12,1,4,13,13,5,3,4,7,3,0,3,6,2,11,6,5,6,13,6,4,6])))
        self.assertEqual(edges.getCoords().getHiddenCppPointer(),coo.getHiddenCppPointer())
        self.assertTrue(d.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,7,6,3,2,12,13,14,1,15,15,0,16,17,9,18,18,8,19,19,13,20,20,12,21])))
        self.assertTrue(di.isEqual(DataArrayInt([0,6,12,18,21,24,27,30,33,36])))
        self.assertTrue(rd.isEqual(DataArrayInt([0,4,0,3,0,2,0,2,0,0,1,2,1,2,1,6,1,5,1,1,2,8,2,7,3,3,4,4,5,5,6,6,7,7,8,8])))
        self.assertTrue(rdi.isEqual(DataArrayInt([0,2,4,6,8,9,10,12,14,16,18,19,20,22,24,25,27,28,29,31,33,35,36])))
        pass

    def testFieldIntIsOnStage2(self):
        """ Very important test to check that isEqual of MEDCouplingFieldInt is OK !"""
        m1=MEDCouplingCMesh() ; m1.setCoords(DataArrayDouble([0,1,2,3]),DataArrayDouble([0,1,2,3,4]))
        m1=m1.buildUnstructured() ; m1.setName("mesh")
        f1=MEDCouplingFieldInt(ON_CELLS) ; f1.setMesh(m1)
        arr1=DataArrayInt([(0,1),(2,3),(4,5),(6,7),(8,9),(10,11),(12,13),(14,15),(16,17),(18,19),(20,21),(22,23)]) ; arr1.setInfoOnComponents(["aa","bbb"])
        f1.setArray(arr1) ; f1.setName("f1") ; f1.setTime(2.,3,4)
        #
        m2=MEDCouplingCMesh() ; m2.setCoords(DataArrayDouble([0,1,2,3]),DataArrayDouble([0,1,2,3,4]))
        m2=m2.buildUnstructured() ; m2.setName("mesh")
        f2=MEDCouplingFieldInt(ON_CELLS) ; f2.setMesh(m2)
        arr2=DataArrayInt([(0,1),(2,3),(4,5),(6,7),(8,9),(10,11),(12,13),(14,15),(16,17),(18,19),(20,21),(22,23)]) ; arr2.setInfoOnComponents(["aa","bbb"])
        f2.setArray(arr2) ; f2.setName("f1") ; f2.setTime(2.,3,4)
        #
        self.assertTrue(f1.isEqual(f2,1e-12,0))
        f1.getArray()[:]*=2
        self.assertTrue(not f1.isEqual(f2,1e-12,0))
        self.assertTrue(not f1.isEqualWithoutConsideringStr(f2,1e-12,0))
        f1.getArray()[:]/=2
        self.assertTrue(f1.isEqual(f2,1e-12,0))
        #
        f1.setName("F1")
        self.assertTrue(not f1.isEqual(f2,1e-12,0))
        f1.setName("f1")
        self.assertTrue(f1.isEqual(f2,1e-12,0))
        #
        f1.getArray().setInfoOnComponents(["aa","bbbb"])
        self.assertTrue(not f1.isEqual(f2,1e-12,0))
        self.assertTrue(f1.isEqualWithoutConsideringStr(f2,1e-12,0))
        f1.getArray().setInfoOnComponents(["aa","bbb"])
        self.assertTrue(f1.isEqual(f2,1e-12,0))
        #
        f3=f2.deepCopy()
        self.assertTrue(f1.isEqual(f3,1e-12,0))
        #
        for fd,expected in ((ON_NODES,False),(ON_CELLS,True)):
            f4=MEDCouplingFieldInt(fd) ; f4.setMesh(m2) ; f4.setTime(2.,3,4)
            arr4=DataArrayInt([(0,1),(2,3),(4,5),(6,7),(8,9),(10,11),(12,13),(14,15),(16,17),(18,19),(20,21),(22,23)]) ; arr4.setInfoOnComponents(["aa","bbb"])
            f4.setArray(arr4) ; f4.setName("f1")
            self.assertEqual(f1.isEqual(f4,1e-12,0),expected)
            pass
        pass

    def testDADSymmetry1(self):
        arr=DataArrayDouble([2,3,4],1,3)
        res=arr.symmetry3DPlane([0.,0.,0.],[0.,0.,2.])
        self.assertTrue(res.isEqual(DataArrayDouble([2,3,-4],1,3),1e-14))
        #
        res=arr.symmetry3DPlane([-1000,100,-1],[0.,0.,2.])
        self.assertTrue(res.isEqual(DataArrayDouble([2,3,-6],1,3),1e-14))
        #
        res=arr.symmetry3DPlane([0,0,0],[1.,0.,0.])
        self.assertTrue(res.isEqual(DataArrayDouble([-2,3,4],1,3),1e-14))
        #
        res=arr.symmetry3DPlane([0,0,0],[0.,1.,0.])
        self.assertTrue(res.isEqual(DataArrayDouble([2,-3,4],1,3),1e-14))
        #
        res=arr.symmetry3DPlane([0,0,0],[-1.,1.,0.])
        self.assertTrue(res.isEqual(DataArrayDouble([3,2,4],1,3),1e-14))
        #
        plane=[5.,4.,-7.]
        a=DataArrayDouble(DataArrayDouble.GiveBaseForPlane(plane))
        self.assertAlmostEqual(DataArrayDouble.Dot(a[0],a[1]).magnitude()[0],0.,13)
        self.assertAlmostEqual(DataArrayDouble.Dot(a[0],a[2]).magnitude()[0],0.,13)
        self.assertAlmostEqual(DataArrayDouble.Dot(a[1],a[2]).magnitude()[0],0.,13)
        coo=DataArrayDouble.Aggregate([10*a[0]+10*a[1],-10*a[0]+10*a[1],-10*a[0]-10*a[1],10*a[0]-10*a[1]])
        m=MEDCouplingUMesh("",2) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_QUAD4,[0,1,2,3])
        d,_=m.distanceToPoint(arr)
        res=arr.symmetry3DPlane([0.,0.,0.],plane) #
        d2,_=m.distanceToPoint(res)
        self.assertAlmostEqual(abs(d-d2),0.,12)
        self.assertAlmostEqual(DataArrayDouble.Dot(res-arr,a[0])[0],0.,12)
        self.assertAlmostEqual(DataArrayDouble.Dot(res-arr,a[1])[0],0.,12)
        self.assertAlmostEqual((res-arr).magnitude()[0]-2*d,0.,12)
        self.assertTrue(res.isEqual(DataArrayDouble([2.666666666666667,3.5333333333333333,3.0666666666666666],1,3),1e-12))
        pass

    def testExtrudedMeshBuildUnstructured1(self):
        """ Non reg test. ExtrudedMesh.buildUnstructured used to modify the coordinates of this. It used to lead to an extra amount of memory consumtion. The aim of the test here is to check that buildUnstructured method do not alter the content of the mesh"""
        arr=DataArrayDouble(11) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr,arr)
        m=m.buildUnstructured()
        faces=MEDCouplingCMesh() ; faces.setCoords(arr,arr)
        faces=faces.buildUnstructured()
        faces.setCoords(m.getCoords())
        em=MEDCouplingMappedExtrudedMesh(m,faces,0)
        self.assertTrue(em.buildUnstructured().isEqual(m,1e-12))
        self.assertTrue(em.buildUnstructured().isEqual(m,1e-12)) # the bug was here ... buildUnstructured used to modify em ...
        self.assertTrue(em.buildUnstructured().isEqual(m,1e-12)) # the bug was here ... buildUnstructured used to modify em ...
        pass

    def testExtrudedMeshFromCMesh1(self):
        arrX=DataArrayDouble([0,1,2,3]) ; arrY=DataArrayDouble([0,1,2,3,4]) ; arrZ=DataArrayDouble([0,1,2,3,4,5])
        mesh3D=MEDCouplingCMesh() ; mesh3D.setCoords(arrX,arrY,arrZ)
        ex=MEDCouplingMappedExtrudedMesh(mesh3D)
        self.assertTrue(ex.buildUnstructured().isEqual(mesh3D.buildUnstructured(),1e-12))
        pass

    def testCylSpherPolarCartFiesta(self):
        """Test to check new capabilities from to cyl spher polar cart conversions"""
        da0=DataArrayDouble([(7,13,2.1),(15,2,-4.2),(-6,12,1.4),(-1,10,-3.5),(-2.1,-3.3,2.7),(-1.4,-0.2,-4),(1.2,-1.3,2.8),(2.5,-0.4,-3)])
        self.assertTrue(da0.fromCartToCyl().fromCylToCart().isEqual(da0,1e-12))
        self.assertTrue(da0.fromCartToSpher().fromSpherToCart().isEqual(da0,1e-12))
        da1=da0[:,:2]
        self.assertTrue(da1.fromCartToPolar().fromPolarToCart().isEqual(da1,1e-12))
        #
        da2=da0[::-1]
        pt=[-2.1,0.3,1.1]
        vect=[1.,-0.5,0.7]
        #
        expected=DataArrayDouble([(2.023252607860588,14.699865529518792,1.4934531458504392),(10.91440936818929,7.5640431386495965,8.384564361982669),(-7.1057844983810705,7.853310978767742,-8.354240440239513),(-8.414001990391881,-1.1910713519565301,-6.405928468241733),(-4.35426264858532,1.5616250027467273,1.0916611827536211),(-2.0571195416878396,-2.0266572603615365,-3.1082019786735042),(-1.5714718759210784,0.39735366651452453,2.8883535460356216),(0.8733250236104675,-3.800053532703407,0.45485882614734185)])
        da4=da0.fromCartToCylGiven(da2,pt,vect)
        self.assertTrue(da4.isEqual(expected,1e-12))
        #
        m=MEDCouplingUMesh.Build0DMeshFromCoords(da2)
        self.assertEqual(m.getDirectAccessOfCoordsArrIfInStructure().getHiddenCppPointer(),da2.getHiddenCppPointer())
        f0=MEDCouplingFieldDouble(ON_NODES) ; f0.setMesh(m) ; f0.setArray(da0)
        f=f0.computeVectorFieldCyl(pt,vect)
        f.checkConsistencyLight()
        self.assertEqual(f.getMesh().getHiddenCppPointer(),m.getHiddenCppPointer())
        self.assertTrue(f.getArray().isEqual(expected,1e-12))
        pass

    def testDAIIndicesOfSubPart(self):
        a=DataArrayInt([9,10,0,6,4,11,3,8])
        b=DataArrayInt([6,0,11,8])
        c=a.indicesOfSubPart(b)
        self.assertTrue(c.isEqual(DataArrayInt([3,2,5,7])))
        #
        d=DataArrayInt([9,10,0,6,4,11,0,8])
        self.assertRaises(InterpKernelException,d.indicesOfSubPart,b) # 0 appears twice in the d array
        f=DataArrayInt([6,0,11,8,12])
        self.assertRaises(InterpKernelException,a.indicesOfSubPart,f) # 12 in f does not exist in a
        pass

    def testDACirPermAndRev1(self):
        d=DataArrayInt([1,2,3,4,5,6])
        d2=d.deepCopy() ; d2.circularPermutation(1)
        self.assertTrue(d2.isEqual(DataArrayInt([2,3,4,5,6,1])))
        d2=d.deepCopy() ; d2.circularPermutation()
        self.assertTrue(d2.isEqual(DataArrayInt([2,3,4,5,6,1])))
        d2=d.deepCopy() ; d2.circularPermutation(2)
        self.assertTrue(d2.isEqual(DataArrayInt([3,4,5,6,1,2])))
        d2=d.deepCopy() ; d2.circularPermutation(3)
        self.assertTrue(d2.isEqual(DataArrayInt([4,5,6,1,2,3])))
        d2=d.deepCopy() ; d2.circularPermutation(4)
        self.assertTrue(d2.isEqual(DataArrayInt([5,6,1,2,3,4])))
        d2=d.deepCopy() ; d2.circularPermutation(5)
        self.assertTrue(d2.isEqual(DataArrayInt([6,1,2,3,4,5])))
        d2=d.deepCopy() ; d2.circularPermutation(6)
        self.assertTrue(d2.isEqual(d))
        d2=d.deepCopy() ; d2.circularPermutation(7)
        self.assertTrue(d2.isEqual(DataArrayInt([2,3,4,5,6,1])))
        d2=d.deepCopy() ; d2.circularPermutation(-1)
        self.assertTrue(d2.isEqual(DataArrayInt([6,1,2,3,4,5])))
        d2=d.deepCopy() ; d2.circularPermutation(-2)
        self.assertTrue(d2.isEqual(DataArrayInt([5,6,1,2,3,4])))
        d2=d.deepCopy() ; d2.circularPermutation(-3)
        self.assertTrue(d2.isEqual(DataArrayInt([4,5,6,1,2,3])))
        d2=d.deepCopy() ; d2.circularPermutation(-4)
        self.assertTrue(d2.isEqual(DataArrayInt([3,4,5,6,1,2])))
        d2=d.deepCopy() ; d2.circularPermutation(-5)
        self.assertTrue(d2.isEqual(DataArrayInt([2,3,4,5,6,1])))
        d2=d.deepCopy() ; d2.circularPermutation(-6)
        self.assertTrue(d2.isEqual(d))
        d2=d.deepCopy() ; d2.circularPermutation(-7)
        self.assertTrue(d2.isEqual(DataArrayInt([6,1,2,3,4,5])))
        ####
        d=DataArrayInt([1,2,3,4,5,6],2,3)
        d2=d.deepCopy() ; d2.circularPermutationPerTuple(0)
        self.assertTrue(d2.isEqual(d))
        d2=d.deepCopy() ; d2.circularPermutationPerTuple(1)
        self.assertTrue(d2.isEqual(DataArrayInt([2,3,1,5,6,4],2,3)))
        d2=d.deepCopy() ; d2.circularPermutationPerTuple()
        self.assertTrue(d2.isEqual(DataArrayInt([2,3,1,5,6,4],2,3)))
        d2=d.deepCopy() ; d2.circularPermutationPerTuple(2)
        self.assertTrue(d2.isEqual(DataArrayInt([3,1,2,6,4,5],2,3)))
        d2=d.deepCopy() ; d2.circularPermutationPerTuple(3)
        self.assertTrue(d2.isEqual(d))
        d2=d.deepCopy() ; d2.circularPermutationPerTuple(-1)
        self.assertTrue(d2.isEqual(DataArrayInt([3,1,2,6,4,5],2,3)))
        d2=d.deepCopy() ; d2.circularPermutationPerTuple(-2)
        self.assertTrue(d2.isEqual(DataArrayInt([2,3,1,5,6,4],2,3)))
        d2=d.deepCopy() ; d2.circularPermutationPerTuple(-3)
        self.assertTrue(d2.isEqual(d))
        d.setInfoOnComponents(["a","b","c"])
        d2=d.deepCopy() ; d2.circularPermutationPerTuple(1)
        self.assertEqual(d2.getInfoOnComponents(),["b","c","a"])
        d2=d.deepCopy() ; d2.circularPermutationPerTuple()
        self.assertEqual(d2.getInfoOnComponents(),["b","c","a"])
        d2=d.deepCopy() ; d2.circularPermutationPerTuple(2)
        self.assertEqual(d2.getInfoOnComponents(),["c","a","b"])
        d2=d.deepCopy() ; d2.circularPermutationPerTuple(3)
        self.assertEqual(d2.getInfoOnComponents(),["a","b","c"])
        d2=d.deepCopy() ; d2.circularPermutationPerTuple(4)
        self.assertEqual(d2.getInfoOnComponents(),["b","c","a"])
        d2=d.deepCopy() ; d2.circularPermutationPerTuple(-1)
        self.assertEqual(d2.getInfoOnComponents(),["c","a","b"])
        d2=d.deepCopy() ; d2.circularPermutationPerTuple(-2)
        self.assertEqual(d2.getInfoOnComponents(),["b","c","a"])
        d2=d.deepCopy() ; d2.circularPermutationPerTuple(-3)
        self.assertEqual(d2.getInfoOnComponents(),["a","b","c"])
        d2=d.deepCopy() ; d2.circularPermutationPerTuple(-4)
        self.assertEqual(d2.getInfoOnComponents(),["c","a","b"])
        ####
        d2=d.deepCopy() ; d2.reversePerTuple()
        d3Exp=DataArrayInt([3,2,1,6,5,4],2,3) ; d3Exp.setInfoOnComponents(["c","b","a"])
        self.assertTrue(d3Exp.isEqual(d2))
        pass

    def testDAExplodeComponents1(self):
        d=DataArrayDouble([(1,2),(3,4),(5,6)])
        d.setName("toto")
        d.setInfoOnComponents(["a","b"])
        d2=d.explodeComponents()
        self.assertEqual(len(d2),2)
        #
        d3=DataArrayDouble([1,3,5]) ; d3.setName("toto") ; d3.setInfoOnComponents(["a"])
        self.assertTrue(d3.isEqual(d2[0],1e-14))
        d4=DataArrayDouble([2,4,6]) ; d4.setName("toto") ; d4.setInfoOnComponents(["b"])
        self.assertTrue(d4.isEqual(d2[1],1e-14))
        #
        d=DataArrayInt([(1,2),(3,4),(5,6)])
        d.setName("toto")
        d.setInfoOnComponents(["a","b"])
        d2=d.explodeComponents()
        self.assertEqual(len(d2),2)
        #
        d3=DataArrayInt([1,3,5]) ; d3.setName("toto") ; d3.setInfoOnComponents(["a"])
        self.assertTrue(d3.isEqual(d2[0]))
        d4=DataArrayInt([2,4,6]) ; d4.setName("toto") ; d4.setInfoOnComponents(["b"])
        self.assertTrue(d4.isEqual(d2[1]))
        pass

    def testVoronoi2D_1(self):
        """ Check of voronize on 2D mesh method of MEDCouplingFieldDouble that converts field on Gauss Points to a field on cell"""
        tmp=MEDCouplingCMesh("mesh")
        arr=DataArrayDouble(5) ; arr.iota()
        tmp.setCoords(arr,arr)
        tmp=tmp.build1SGTUnstructured()
        conn=tmp.getNodalConnectivity()
        conn.rearrange(4)
        conn.reversePerTuple()
        conn.circularPermutationPerTuple(2)
        conn.rearrange(1)
        coo=tmp.getCoords().deepCopy()
        coo.circularPermutationPerTuple(2) ; coo*=0.1
        coo.reverse()
        coo2=DataArrayDouble(len(tmp.getCoords())*tmp.getSpaceDimension()) ; coo2.iota() ; coo2.rearrange(tmp.getSpaceDimension())
        coo2*=0.14
        coo2.circularPermutationPerTuple(2)
        tmp.getCoords()[:]+=coo2*coo
        #
        field=MEDCouplingFieldDouble(ON_GAUSS_PT)
        field.setName("MyFieldPG") ; field.setMesh(tmp)
        field.setGaussLocalizationOnType(NORM_QUAD4,[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.8,-0.8, 0.8,0.8, -0.8,0.8, -0.8,-0.8, 0.,0., 0.2,0.2, 0.1,0.3],[0.1,0.1,0.1,0.1,0.1,0.1,0.4])
        arr=DataArrayDouble(field.getNumberOfTuplesExpected()) ; arr.iota()
        field.setArray(arr)
        field.checkConsistencyLight()
        ####
        fieldOnCell=field.voronoize(1e-12) # hot point
        fieldOnCell.checkConsistencyLight()
        self.assertEqual(fieldOnCell.getMesh().getNumberOfCells(),112)
        self.assertEqual(fieldOnCell.getMesh().getNumberOfNodes(),256)
        self.assertTrue(fieldOnCell.getArray().isEqual(field.getArray(),1e-12))
        meaRef=field.getMesh().getMeasureField(True).getArray()
        mea=fieldOnCell.getMesh().getMeasureField(True).getArray()
        self.assertEqual(field.getDiscretization().getNbOfGaussLocalization(),1)
        self.assertEqual(field.getDiscretization().getGaussLocalization(0).getNumberOfGaussPt(),7)
        mea.rearrange(7)
        mea2=mea.sumPerTuple()
        self.assertTrue(mea2.isEqual(meaRef,1e-12))
        pass

    def testVoronoi2D_2(self):
        """More aggressive 2D test. No warping here. To check data"""
        tmp=MEDCouplingCMesh("mesh")
        arr=DataArrayDouble([-1.,1.])
        tmp.setCoords(arr,arr)
        tmp=tmp.buildUnstructured()
        field=MEDCouplingFieldDouble(ON_GAUSS_PT)
        field.setName("MyFieldPG") ; field.setMesh(tmp)
        field.setGaussLocalizationOnType(NORM_QUAD4,[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.8,-0.8, 0.8,0.8, -0.8,0.8, -0.8,-0.8, 0.,0., 0.2,0.2, 0.1,0.3],[0.1,0.1,0.1,0.1,0.1,0.1,0.4])
        arr=DataArrayDouble(field.getNumberOfTuplesExpected()) ; arr.iota()
        field.setArray(arr)
        field.checkConsistencyLight()
        #
        fieldOnCell=field.voronoize(1e-12) # hot point
        fieldOnCell.checkConsistencyLight()
        self.assertEqual(fieldOnCell.getMesh().getNumberOfCells(),7)
        self.assertEqual(fieldOnCell.getMesh().getNumberOfNodes(),16)
        self.assertTrue(fieldOnCell.getArray().isEqual(field.getArray(),1e-12))
        meaRef=DataArrayDouble([0.65,0.4710714285714285,0.59875,0.68,0.73875,0.4,0.46142857142857235])
        mea=fieldOnCell.getMesh().getMeasureField(True).getArray()
        self.assertTrue(mea.isEqual(meaRef,1e-12))# the first important test is here
        self.assertEqual(field.getDiscretization().getNbOfGaussLocalization(),1)
        self.assertEqual(field.getDiscretization().getGaussLocalization(0).getNumberOfGaussPt(),7)
        #
        gsPt=field.getLocalizationOfDiscr()
        a,b=fieldOnCell.getMesh().getCellsContainingPoints(gsPt,1e-12)
        self.assertTrue(a.isIota(7))# the second important test is here ! Check that Gauss points are inside the associated cell in fieldOnCell !
        self.assertTrue(b.isIota(8))
        #
        self.assertEqual(fieldOnCell.getMesh().buildDescendingConnectivity()[0].getNumberOfCells(),22)# last little test to reduce chance of errors. For humans there 21 but last tiny edge is split into 2 subedges due to alg
        pass

    def testVoronoi3D_1(self):
        """ Check of voronize on 3D mesh method of MEDCouplingFieldDouble that converts field on Gauss Points to a field on cell"""
        tmp=MEDCouplingCMesh("mesh")
        arr=DataArrayDouble(5) ; arr.iota()
        tmp.setCoords(arr,arr)
        tmp=tmp.build1SGTUnstructured()
        conn=tmp.getNodalConnectivity()
        conn.rearrange(4)
        conn.reversePerTuple()
        conn.circularPermutationPerTuple(2)
        conn.rearrange(1)
        coo=tmp.getCoords().deepCopy()
        coo.circularPermutationPerTuple(2) ; coo*=0.1
        coo.reverse()
        coo2=DataArrayDouble(len(tmp.getCoords())*tmp.getSpaceDimension()) ; coo2.iota() ; coo2.rearrange(tmp.getSpaceDimension())
        coo2*=0.14
        coo2.circularPermutationPerTuple(2)
        tmp.getCoords()[:]+=coo2*coo
        #
        tmp.changeSpaceDimension(3,0.)
        #
        arrZ=DataArrayDouble(5) ; arrZ.iota()
        mz=MEDCouplingCMesh() ; mz.setCoords(arrZ) ; mz=mz.buildUnstructured()
        mz.changeSpaceDimension(3,0.)
        mz.getCoords().circularPermutationPerTuple(1)
        tmp=tmp.buildUnstructured().buildExtrudedMesh(mz,0)
        #
        field=MEDCouplingFieldDouble(ON_GAUSS_PT)
        field.setName("MyFieldPG") ; field.setMesh(tmp)
        field.setGaussLocalizationOnType(NORM_HEXA8,[-1,-1,-1,  1,-1,-1,  1,1,-1,  -1,1,-1, -1,-1,1, 1,-1,1, 1,1,1, -1,1,1],[0.8,-0.8,0., 0.8,0.8,0., -0.8,0.8,0., -0.8,-0.8,0., 0.,0.,0., 0.2,0.2,0., 0.1,0.3,0.],[0.1,0.1,0.1,0.1,0.1,0.1,0.4])
        arr=DataArrayDouble(field.getNumberOfTuplesExpected()) ; arr.iota()
        field.setArray(arr)
        field.checkConsistencyLight()
        ####
        fieldOnCell=field.voronoize(1e-12) # hot point
        fieldOnCell.checkConsistencyLight()
        self.assertTrue(fieldOnCell.getArray().isEqual(field.getArray(),1e-12))
        meaRef=field.getMesh().getMeasureField(True).getArray()
        mea=fieldOnCell.getMesh().getMeasureField(False).getArray()
        self.assertEqual(field.getDiscretization().getNbOfGaussLocalization(),1)
        self.assertEqual(field.getDiscretization().getGaussLocalization(0).getNumberOfGaussPt(),7)
        mea.rearrange(7)
        mea2=mea.sumPerTuple()
        delta=(meaRef-mea2)
        delta.abs()
        delta/=meaRef
        self.assertEqual(len(delta.findIdsNotInRange(0,1e-2)),0) # 1e-2 because hexa8 are warped !
        pass

    def testVoronoi3D_2(self):
        """More aggressive 3D test. No warping here. To check data"""
        tmp=MEDCouplingCMesh("mesh")
        arr=DataArrayDouble([-1.,1.])
        tmp.setCoords(arr,arr,arr)
        tmp=tmp.buildUnstructured()
        field=MEDCouplingFieldDouble(ON_GAUSS_PT)
        field.setName("MyFieldPG") ; field.setMesh(tmp)
        field.setGaussLocalizationOnType(NORM_HEXA8,[-1,-1,-1,  1,-1,-1,  1,1,-1,  -1,1,-1, -1,-1,1, 1,-1,1, 1,1,1, -1,1,1],[0.8,-0.8,0., 0.8,0.8,0., -0.8,0.8,0., -0.8,-0.8,0., 0.,0.,0., 0.2,0.2,0., 0.1,0.3,0.],[0.1,0.1,0.1,0.1,0.1,0.1,0.4])
        arr=DataArrayDouble(field.getNumberOfTuplesExpected()) ; arr.iota()
        field.setArray(arr)
        field.checkConsistencyLight()
        #
        fieldOnCell=field.voronoize(1e-12) # hot point
        fieldOnCell.checkConsistencyLight()
        self.assertEqual(fieldOnCell.getMesh().getNumberOfCells(),7)
        self.assertEqual(fieldOnCell.getMesh().getNumberOfNodes(),34)
        self.assertTrue(fieldOnCell.getArray().isEqual(field.getArray(),1e-12))
        meaRef=DataArrayDouble([1.3,0.9421428571428572,1.1975,1.36,1.4775,0.8,0.922857142857143])
        mea=fieldOnCell.getMesh().getMeasureField(True).getArray()
        self.assertTrue(mea.isEqual(meaRef,1e-12))# the first important test is here
        self.assertEqual(field.getDiscretization().getNbOfGaussLocalization(),1)
        self.assertEqual(field.getDiscretization().getGaussLocalization(0).getNumberOfGaussPt(),7)
        #
        gsPt=field.getLocalizationOfDiscr()
        a,b=fieldOnCell.getMesh().getCellsContainingPoints(gsPt,1e-12)
        self.assertTrue(a.isIota(7))# the second important test is here ! Check that Gauss points are inside the associated cell in fieldOnCell !
        self.assertTrue(b.isIota(8))
        #
        self.assertEqual(fieldOnCell.getMesh().buildDescendingConnectivity()[0].getNumberOfCells(),2*7+21)
        pass

    def testVoronoi3D_8(self):
        """More aggressive 3D test. Bug EDF 15094"""
        mesh = MEDCouplingUMesh("myMeshForAnthony",3)
        coords = [2.20449946892035, 0.0015302058397972198, -0.014025000000000001, 2.20449522028465, 0.00459061457029268, -0.0109750000232271, 2.20449946892035, 0.0015302058397972198, -0.0125000000116135, 2.20577243296484, 0.00153108944037966, -0.0137555135576553, 2.20517315768831, 0.0045920262990614006, -0.010764118475206199, 2.2054749202977, 0.0015308829283677198, -0.012259816016430801, 2.20449787568164, 0.00306041094231961, -0.0125000000116135, 2.20449787568164, 0.00306041094231961, -0.011737500017420301, 2.20449946892035, 0.0015302058397972198, -0.0132625000058068, 2.20513595094259, 0.0015306476400884401, -0.0138902567788277, 2.20483418898648, 0.0045913204346770395, -0.0108695592492167, 2.20498719460902, 0.00153054438408247, -0.0123799080140222, 2.20547332635401, 0.0030617651191343705, -0.012259816016430801, 2.20532457012796, 0.00306155860717217, -0.0115119672458185, 2.20562367663127, 0.0015309861843736902, -0.013007664787043, 2.20582504233773, 0.0045933837758852306, -0.010139577890770399, 2.20642582267143, 0.004594634833691141, -0.009125379014333041, 2.20612543250458, 0.00459400930478819, -0.00963247845255172, 2.2069524110381, 0.004595731395029229, -0.00776049693994639, 2.20668911685476, 0.004595183114360191, -0.00844293797713971, 2.20832419990944, 0.0076643330146060895, -0.0108392857142857, 2.20832419990944, 0.0076643330146060895, -0.008671428571428571, 2.20704504094678, 0.00765989349423635, -0.008671428571428571, 2.20704504094678, 0.00765989349423635, -0.0108392857142857, 2.2062381754171, 0.00459424407928538, -0.00868052596233734, 2.20832419990944, 0.0076643330146060895, -0.00975535714285714, 2.20768462042811, 0.00766211325442122, -0.008671428571428571, 2.20704504094678, 0.00765989349423635, -0.00975535714285714, 2.20768462042811, 0.00766211325442122, -0.0108392857142857, 2.20737554490036, 0.00612882358882901, -0.009982332364309381, 2.20763883863969, 0.00612955462931014, -0.00821596275568748, 2.2066421405633703, 0.00612678727660696, -0.00867597726688296, 2.20643557437203, 0.006126213741329251, -0.0104894318025281, 2.2065952932276, 0.00459498773715731, -0.00822051145114186, 2.20603160887741, 0.00459381392758531, -0.00941005192655387]
        da = DataArrayDouble.New(coords,35,3)
        mesh.setCoords(da)
        mesh.allocateCells()
        mesh.insertNextCell(NORM_PENTA15, [0, 2, 1, 3, 5, 4, 8, 7, 6, 14, 13, 12, 9, 11, 10])
        mesh.insertNextCell(NORM_HEXA20, [20, 23, 22, 21, 16, 15, 24, 18, 28, 27, 26, 25, 17, 34, 33, 19, 29, 32, 31, 30])
        mesh.zipCoords()
        f=MEDCouplingFieldDouble(ON_GAUSS_PT)
        f.setMesh(mesh)
        f.setName("myFieldForAnthony")
        f.setGaussLocalizationOnCells([0],[-1, 1, 0, -1, 0, 0, -1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, -1, 0.5, 0, -1, 0, 0.5, -1, 0.5, 0.5, 1, 0.5, 0, 1, 0, 0.5, 1, 0.5, 0.5, 0, 1, 0, 0, 0, 0, 0, 0, 1],[-0.774597, 0.333333, 0.333333, -0.774597, 0.470142, 0.470142, -0.774597, 0.0597159, 0.470142, -0.774597, 0.470142, 0.0597159, -0.774597, 0.101287, 0.101287, -0.774597, 0.797427, 0.101287, -0.774597, 0.101287, 0.797427, 0, 0.333333, 0.333333, 0, 0.470142, 0.470142, 0, 0.0597159, 0.470142, 0, 0.470142, 0.0597159, 0, 0.101287, 0.101287, 0, 0.797427, 0.101287, 0, 0.101287, 0.797427, 0.774597, 0.333333, 0.333333, 0.774597, 0.470142, 0.470142, 0.774597, 0.0597159, 0.470142, 0.774597, 0.470142, 0.0597159, 0.774597, 0.101287, 0.101287, 0.774597, 0.797427, 0.101287, 0.774597, 0.101287, 0.797427],[0.0625, 0.0367762, 0.0367762, 0.0367762, 0.0349831, 0.0349831, 0.0349831, 0.1, 0.0588418, 0.0588418, 0.0588418, 0.055973, 0.055973, 0.055973, 0.0625, 0.0367762, 0.0367762, 0.0367762, 0.0349831, 0.0349831, 0.0349831])
        f.setGaussLocalizationOnCells([1],[-1, -1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 1, -1, 0, -1, 0, 1, -1, 1, 0, -1, 0, -1, -1, -1, 0, 1, 0, 1, 1, 1, 0, 1, 0, -1, 1, -1, -1, 0, -1, 1, 0, 1, 1, 0, 1, -1, 0],[-0.774597, -0.774597, -0.774597, -0.774597, -0.774597, 0, -0.774597, -0.774597, 0.774597, -0.774597, 0, -0.774597, -0.774597, 0, 0, -0.774597, 0, 0.774597, -0.774597, 0.774597, -0.774597, -0.774597, 0.774597, 0, -0.774597, 0.774597, 0.774597, 0, -0.774597, -0.774597, 0, -0.774597, 0, 0, -0.774597, 0.774597, 0, 0, -0.774597, 0, 0, 0, 0, 0, 0.774597, 0, 0.774597, -0.774597, 0, 0.774597, 0, 0, 0.774597, 0.774597, 0.774597, -0.774597, -0.774597, 0.774597, -0.774597, 0, 0.774597, -0.774597, 0.774597, 0.774597, 0, -0.774597, 0.774597, 0, 0, 0.774597, 0, 0.774597, 0.774597, 0.774597, -0.774597, 0.774597, 0.774597, 0, 0.774597, 0.774597, 0.774597],[0.171468, 0.274348, 0.171468, 0.274348, 0.438957, 0.274348, 0.171468, 0.274348, 0.171468, 0.274348, 0.438957, 0.274348, 0.438957, 0.702332, 0.438957, 0.274348, 0.438957, 0.274348, 0.171468, 0.274348, 0.171468, 0.274348, 0.438957, 0.274348, 0.171468, 0.274348, 0.171468])
        arr = DataArrayDouble(48, 3)
        arr[:, 0] = list(range(48))
        arr[:, 1] = 100 + arr[:, 0]
        arr[:, 2] = 200 + arr[:, 0]
        f.setArray(arr)
        fieldOnCell=f.voronoize(1e-12) # hot point
        fieldOnCell.checkConsistencyLight()
        self.assertEqual(fieldOnCell.getMesh().getNumberOfCells(),48)
        self.assertEqual(fieldOnCell.getMesh().getNumberOfNodes(),127)
        meaRef=f.getMesh().getMeasureField(True).getArray(); meaRef.rearrange(2); meaRef2 = meaRef.sumPerTuple()
        mea=fieldOnCell.getMesh().getMeasureField(True).getArray(); mea.rearrange(48); mea2 = mea.sumPerTuple()
        self.assertTrue(mea2.isEqual(meaRef2,1e-9))
        pass

    def testVoronoi3DSurf_1(self):
        tmp=MEDCouplingCMesh("mesh")
        arr=DataArrayDouble(5) ; arr.iota()
        tmp.setCoords(arr,arr)
        tmp=tmp.build1SGTUnstructured()
        conn=tmp.getNodalConnectivity()
        conn.rearrange(4)
        conn.reversePerTuple()
        conn.circularPermutationPerTuple(2)
        conn.rearrange(1)
        coo=tmp.getCoords().deepCopy()
        coo.circularPermutationPerTuple(2) ; coo*=0.1
        coo.reverse()
        coo2=DataArrayDouble(len(tmp.getCoords())*tmp.getSpaceDimension()) ; coo2.iota() ; coo2.rearrange(tmp.getSpaceDimension())
        coo2*=0.14
        coo2.circularPermutationPerTuple(2)
        tmp.getCoords()[:]+=coo2*coo
        #
        tmp.changeSpaceDimension(3,0.)    # force 3D surf
        tmp.rotate([0,0,0],[1,0,0],pi/3)  # force 3D surf
        #
        field=MEDCouplingFieldDouble(ON_GAUSS_PT)
        field.setName("MyFieldPG") ; field.setMesh(tmp)
        field.setGaussLocalizationOnType(NORM_QUAD4,[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.8,-0.8, 0.8,0.8, -0.8,0.8, -0.8,-0.8, 0.,0., 0.2,0.2, 0.1,0.3],[0.1,0.1,0.1,0.1,0.1,0.1,0.4])
        arr=DataArrayDouble(field.getNumberOfTuplesExpected()) ; arr.iota()
        field.setArray(arr)
        field.checkConsistencyLight()
        #####
        fieldOnCell=field.voronoize(1e-12);
        fieldOnCell.checkConsistencyLight()
        self.assertEqual(fieldOnCell.getMesh().getSpaceDimension(),3)
        self.assertEqual(fieldOnCell.getMesh().getMeshDimension(),2)
        self.assertEqual(field.getMesh().getSpaceDimension(),fieldOnCell.getMesh().getSpaceDimension())
        self.assertTrue(fieldOnCell.getArray().isEqual(field.getArray(),1e-12))
        meaRef=field.getMesh().getMeasureField(True).getArray()
        mea=fieldOnCell.getMesh().getMeasureField(True).getArray()
        self.assertEqual(field.getDiscretization().getNbOfGaussLocalization(),1)
        self.assertEqual(field.getDiscretization().getGaussLocalization(0).getNumberOfGaussPt(),7)
        mea.rearrange(7)
        mea2=mea.sumPerTuple()
        self.assertTrue(mea2.isEqual(meaRef,1e-12))
        pass

    def testVoronoi1D_1(self):
        tmp=MEDCouplingCMesh("mesh")
        arr=DataArrayDouble(5) ; arr.iota()
        tmp.setCoords(arr)
        tmp=tmp.build1SGTUnstructured()
        tmp1=tmp.deepCopy()
        tmp.changeSpaceDimension(2,0.)
        tmp.getCoords()[:,1]=pi/(len(arr)-1)*tmp.getCoords()[:,0]
        tmp.getCoords()[:,0]=1.
        tmp.setCoords(tmp.getCoords().fromPolarToCart())
        tmp.changeSpaceDimension(3,1.)
        #
        field=MEDCouplingFieldDouble(ON_GAUSS_PT)
        field.setName("MyFieldPG") ; field.setMesh(tmp)
        field.setGaussLocalizationOnType(NORM_SEG2,[-1.,1.],[-0.9,-0.8,0.2,0.4,0.5,0.9],[0.1,0.1,0.1,0.1,0.1,0.5])
        arr=DataArrayDouble(field.getNumberOfTuplesExpected()) ; arr.iota()
        field.setArray(arr)
        field.checkConsistencyLight()
        ####
        fieldOnCell=field.voronoize(1e-12);
        fieldOnCell.checkConsistencyLight()
        self.assertEqual(fieldOnCell.getMesh().getSpaceDimension(),3)
        self.assertEqual(fieldOnCell.getMesh().getMeshDimension(),1)
        assert(fieldOnCell.getArray().isEqual(field.getArray(),1e-12))
        meaRef=field.getMesh().getMeasureField(True).getArray()
        mea=fieldOnCell.getMesh().getMeasureField(True).getArray()
        self.assertEqual(field.getDiscretization().getNbOfGaussLocalization(),1)
        self.assertEqual(field.getDiscretization().getGaussLocalization(0).getNumberOfGaussPt(),6)
        mea.rearrange(6)
        mea2=mea.sumPerTuple()
        self.assertTrue(mea2.isEqual(meaRef,1e-12))
        pass

    def testFieldDoubleConvertToLinear1(self):
        da=DataArrayDouble([0,0, 1,0, 2,0, 3,0, 0.5,0, 1.5,0, 2.5,0, 0,0.5, 0.5,0.5, 1, 0.5, 1.5,0.5, 2,0.5, 3,0.5, 0,1, 1,1, 2,1, 2.5,1, 3,1],18,2)
        da.setInfoOnComponents(["g","h"])
        m=MEDCouplingUMesh("mesh",2)
        m.setCoords(da)
        m.allocateCells()
        m.insertNextCell(NORM_TRI6,[0,1,13,4,9,7])
        m.insertNextCell(NORM_TRI6,[1,2,14,5,10,9])
        m.insertNextCell(NORM_QUAD8,[2,3,17,15,6,12,16,11])
        refPtr=m.getHiddenCppPointer()
        f=MEDCouplingFieldDouble(ON_NODES)
        f.setName("aa")
        f.setMesh(m)
        arr=DataArrayDouble(18*2) ; arr.iota()
        arr.rearrange(2)
        arr.setInfoOnComponents(["bb","ccc"])
        f.setArray(arr)
        f.setTime(0.5,2,3)
        f.checkConsistencyLight()
        #
        f1=f.convertQuadraticCellsToLinear()
        self.assertTrue(f.getMesh().getHiddenCppPointer(),refPtr)
        self.assertTrue(f1.getMesh().getHiddenCppPointer()!=refPtr)
        f1.checkConsistencyLight()
        self.assertEqual(f1.getName(),"aa")
        self.assertEqual(f1.getTypeOfField(),ON_NODES)
        da0=DataArrayDouble([(0,0),(1,0),(2,0),(3,0),(0,1),(1,1),(2,1),(3,1)])
        da0.setInfoOnComponents(["g","h"])
        self.assertTrue(f1.getMesh().getCoords().isEqual(da0,1e-12))
        self.assertTrue(f1.getMesh().getNodalConnectivity().isEqual(DataArrayInt([3,0,1,4,3,1,2,5,4,2,3,7,6])))
        self.assertTrue(f1.getMesh().getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,8,13])))
        da2=DataArrayDouble([(0,1),(2,3),(4,5),(6,7),(26,27),(28,29),(30,31),(34,35)])
        da2.setInfoOnComponents(["bb","ccc"])
        self.assertTrue(f1.getArray().isEqual(da2,1e-12))
        self.assertEqual(f1.getTime(),[0.5,2,3])
        #
        f2=MEDCouplingFieldDouble(ON_CELLS)
        f2.setName("aa")
        f2.setMesh(m)
        arr=DataArrayDouble(3*2) ; arr.iota()
        arr.rearrange(2)
        arr.setInfoOnComponents(["bb","ccc"])
        f2.setArray(arr)
        f2.setTime(0.5,2,3)
        f2.checkConsistencyLight()
        f3=f2.convertQuadraticCellsToLinear()
        self.assertEqual(f2.getMesh().getHiddenCppPointer(),refPtr)
        f3.checkConsistencyLight()
        self.assertTrue(f3.getMesh().getHiddenCppPointer()!=refPtr)
        self.assertTrue(f3.getMesh().getCoords().isEqual(da0,1e-12))
        self.assertTrue(f3.getMesh().getNodalConnectivity().isEqual(DataArrayInt([3,0,1,4,3,1,2,5,4,2,3,7,6])))
        self.assertTrue(f3.getMesh().getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,8,13])))
        self.assertEqual(f2.getArray().getHiddenCppPointer(),f3.getArray().getHiddenCppPointer())
        self.assertEqual(f3.getTime(),[0.5,2,3])
        pass

    def testBuild1DMeshFromCoords1(self):
        da=DataArrayDouble([(3,4),(5,6),(7,8)])
        da.setName("ZeArr")
        da0=da.deepCopy()
        m=MEDCouplingUMesh.Build1DMeshFromCoords(da0)
        m.checkConsistencyLight()
        self.assertEqual(da0.getHiddenCppPointer(),m.getCoords().getHiddenCppPointer())
        self.assertTrue(da.isEqual(da0,1e-12))
        self.assertEqual(m.getName(),da.getName())
        self.assertEqual(m.getMeshDimension(),1)
        self.assertTrue(isinstance(m,MEDCouplingUMesh))
        m1=MEDCoupling1SGTUMesh(m)
        m1.checkConsistencyLight()
        self.assertTrue(m1.getNodalConnectivity().isEqual(DataArrayInt([0,1,1,2])))
        #
        da0.setName("")
        m2=MEDCouplingUMesh.Build1DMeshFromCoords(da0)
        m2.checkConsistencyLight()
        self.assertEqual(da0.getHiddenCppPointer(),m2.getCoords().getHiddenCppPointer())
        self.assertEqual(m2.getName(),"Mesh")
        pass

    def testVoronoi3D_3(self):
        """Non regression test to check MEDCouplingUMesh::clipSingle3DCellByPlane"""
        coo=DataArrayDouble([0.,1.,0.,0.,0.,0.,0.,0.,1.,1.,0.,0.],4,3)
        m=MEDCouplingUMesh("mesh",3)
        m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_TETRA4,[0,2,3,1])
        f=MEDCouplingFieldDouble(ON_GAUSS_PT)
        f.setMesh(m) ; f.setName("field")
        f.setGaussLocalizationOnType(NORM_TETRA4,[0.,1.,0.,0.,0.,0.,0.,0.,1.,1.,0.,0.],[0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.1381966011250105], [0.041667,0.041667,0.041667,0.041667])
        f.setArray(DataArrayDouble([0,1,2,3]))
        f3=f.voronoize(1e-12)
        ref=DataArrayDouble([0.047256836610416179,0.03980327668541684,0.039803276685416833,0.039803276685416833])
        self.assertTrue(f3.getMesh().getMeasureField(False).getArray().isEqual(ref,1e-12))
        self.assertTrue(f3.getArray().isEqual(DataArrayDouble([0,1,2,3]),1e-12))
        pass

    def testVoronoi3D_4(self):
        """Idem testVoronoi3D_3 except that here quadratic cells are considered"""
        coo=DataArrayDouble([0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.5,0.0,0.5,0.5,0.5,0.5,0.0,0.5,0.0,0.0,0.5,0.0,0.5],10,3)
        m=MEDCouplingUMesh("mesh",3)
        m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_TETRA10,[0,1,2,3,4,5,6,7,8,9])
        f=MEDCouplingFieldDouble(ON_GAUSS_PT)
        f.setMesh(m) ; f.setName("field")
        f.setGaussLocalizationOnType(NORM_TETRA10,[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.5,0.0,0.5,0.5,0.5,0.5,0.0,0.5,0.0,0.0,0.5,0.0,0.5],[0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.1381966011250105], [0.041667,0.041667,0.041667,0.041667])
        f.setArray(DataArrayDouble([0,1,2,3]))
        f3=f.voronoize(1e-12)
        ref=DataArrayDouble([0.047256836610416179,0.03980327668541684,0.039803276685416833,0.039803276685416833])
        self.assertTrue(f3.getMesh().getMeasureField(False).getArray().isEqual(ref,1e-12))
        self.assertTrue(f3.getArray().isEqual(DataArrayDouble([0,1,2,3]),1e-12))
        pass

    def testVoronoi3D_5(self):
        """ Cell 0 of Barreau_Elga_V11.rmed and sslv07b.rmed. HEXA8 cut regularly into 8 parts"""
        coo=DataArrayDouble([(0.024,0.024,1.2),(0.024,0.048,1.2),(0.048,0.024,1.2),(0.048,0.048,1.2),(0.024,0.024,1.6),(0.024,0.048,1.6),(0.048,0.024,1.6),(0.048,0.048,1.6)])
        m=MEDCouplingUMesh("",3) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_HEXA8,[0,2,6,4,1,3,7,5])
        f=MEDCouplingFieldDouble(ON_GAUSS_PT) ; f.setMesh(m)
        f.setGaussLocalizationOnType(NORM_HEXA8,[-1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0],[-0.577350269189626, -0.577350269189626, -0.577350269189626, -0.577350269189626, -0.577350269189626, 0.577350269189626, -0.577350269189626, 0.577350269189626, -0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626, 0.577350269189626, -0.577350269189626, -0.577350269189626, 0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626, 0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626, 0.577350269189626],[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        arr=DataArrayDouble(8) ; arr.iota() ; f.setArray(arr)
        f.checkConsistencyLight()
        #
        vol=f.getMesh().getMeasureField(False).getIJ(0,0)
        f2=f.voronoize(1e-12)
        f2.checkConsistencyLight()
        self.assertEqual(f2.getNumberOfTuples(),8)
        volRef=vol/8
        self.assertTrue(f2.getMesh().getMeasureField(False).getArray().isUniform(volRef,1e-12))
        pass

    def testVoronoi3D_6(self):
        """ Cell 0 of brokenshire.med (and pace.med). TETRA10 split into 4 parts"""
        coo=DataArrayDouble([(50.,-50.,200.0),(50.0,-30.,200.0),(30.,-50.,200.0),(50.,-50.,180.0),(50.,-40.,200.0),(40.,-50.,200.0),(50.,-50.,190.0),(40.,-40.,200.0),(50.,-40.,190.0),(40.,-50.,190.0)])
        m=MEDCouplingUMesh("",3) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_TETRA10,[2,0,1,3,5,4,7,9,6,8])
        f=MEDCouplingFieldDouble(ON_GAUSS_PT) ; f.setMesh(m)
        f.setGaussLocalizationOnType(NORM_TETRA10,[0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0.5, 0.5, 0.5, 0.5, 0, 0.5, 0, 0, 0.5, 0, 0.5],[0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.1381966011250105],[0.041666666666666664, 0.041666666666666664, 0.041666666666666664, 0.041666666666666664])
        arr=DataArrayDouble(4) ; arr.iota() ; f.setArray(arr)
        f.checkConsistencyLight()
        f2=f.voronoize(1e-12)
        f2.checkConsistencyLight()
        self.assertEqual(f2.getNumberOfTuples(),4)
        arr=f2.getMesh().getMeasureField(False).getArray()
        self.assertTrue(f2.getMesh().getMeasureField(False).getArray().isEqual(DataArrayDouble([378.0546928833331, 318.42621348333586, 318.4262134833361, 318.4262134833278]),1e-6))
        pass

    def testVoronoi3D_7(self):
        """ sslv07a.rmed. HEXA20 split into 27 parts """
        coo=DataArrayDouble([(-0.5,-0.5,0.0),(-0.25,-0.5,0.0),(0.0,-0.5,0.0),(-0.5,0.0,0.0),(-0.5,-0.25,0.0),(0.0,0.0,0.0),(0.0,-0.25,0.0),(-0.25,0.0,0.0),(-0.5,-0.5,1.0),(-0.25,-0.5,1.0),(0.0,-0.5,1.0),(0.0,-0.25,1.0),(0.0,0.0,1.0),(-0.25,0.0,1.0),(-0.5,0.0,1.0),(-0.5,-0.25,1.0),(-0.5,-0.5,0.5),(0.0,-0.5,0.5),(0.0,0.0,0.5),(-0.5,0.0,0.5)])
        m=MEDCouplingUMesh("",3) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_HEXA20,[0,3,5,2,8,14,12,10,4,7,6,1,15,13,11,9,16,19,18,17])
        f=MEDCouplingFieldDouble(ON_GAUSS_PT) ; f.setMesh(m)
        f.setGaussLocalizationOnType(NORM_HEXA20,
                                     [-1,-1,-1,-1,1,-1,1,1,-1,1,-1,-1,-1,-1,1,-1,1,1,1,1,1,1,-1,1,-1,0,-1,0,1,-1,1,0,-1,0,-1,-1,-1,0,1,0,1,1,1,0,1,0,-1,1,-1,-1,0,-1,1,0,1,1,0,1,-1,0],
                                     [-0.774597,-0.774597,-0.774597,-0.774597,-0.774597,0,-0.774597,-0.774597,0.774597,-0.774597,0,-0.774597,-0.774597,0,0,-0.774597,0,0.774597,-0.774597,0.774597,-0.774597,-0.774597,0.774597,0,-0.774597,0.774597,0.774597,0,-0.774597,-0.774597,0,-0.774597,0,0,-0.774597,0.774597,0,0,-0.774597,0,0,0,0,0,0.774597,0,0.774597,-0.774597,0,0.774597,0,0,0.774597,0.774597,0.774597,-0.774597,-0.774597,0.774597,-0.774597,0,0.774597,-0.774597,0.774597,0.774597,0,-0.774597,0.774597,0,0,0.774597,0,0.774597,0.774597,0.774597,-0.774597,0.774597,0.774597,0,0.774597,0.774597,0.774597],
                                     [0.171468,0.274348,0.171468,0.274348,0.438957,0.274348,0.171468,0.274348,0.171468,0.274348,0.438957,0.274348,0.438957,0.702332,0.438957,0.274348,0.438957,0.274348,0.171468,0.274348,0.171468,0.274348,0.438957,0.274348,0.171468,0.274348,0.171468])
        arr=DataArrayDouble(27) ; arr.iota() ; f.setArray(arr)
        f.checkConsistencyLight()
        f2=f.voronoize(1e-12)
        a=0.007187820185770747 ; b=0.0090870678008658 ; c=0.011488156225861077 ; d=0.014523687548277797
        ref=DataArrayDouble(27) ; ref[::2]=a ; ref[1::2]=b
        ref[[4,10,12,14,16,22]]=c ; ref[13]=d  # 6 cells 4,10,12,14,16,22 are the 6 cells boarding the most inner cell 13
        #
        self.assertTrue(f2.getMesh().getMeasureField(False).getArray().isEqual(ref,1e-7))
        pass

    def testConvertQuadToLin4Gauss_1(self):
        coo=DataArrayDouble([0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.5,0.0,0.5,0.5,0.5,0.5,0.0,0.5,0.0,0.0,0.5,0.0,0.5],10,3)
        m=MEDCouplingUMesh("mesh",3)
        m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_TETRA10,[0,1,2,3,4,5,6,7,8,9])
        f=MEDCouplingFieldDouble(ON_GAUSS_PT)
        f.setMesh(m) ; f.setName("field")
        aaaa=[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.5,0.0,0.5,0.5,0.5,0.5,0.0,0.5,0.0,0.0,0.5,0.0,0.5]
        bbbb=[0.1381966011250105,0.1381966011250105,0.1381966011250105,0.1381966011250105,0.1381966011250105,0.5854101966249685,0.1381966011250105,0.5854101966249685,0.1381966011250105,0.5854101966249685,0.1381966011250105,0.1381966011250105]
        cccc=[0.041667,0.041667,0.041667,0.041667]
        f.setGaussLocalizationOnType(NORM_TETRA10,aaaa,bbbb,cccc)
        f.setArray(DataArrayDouble([0,1,2,3]))
        f.setTime(1.,2,3)
        #
        mcpy=m.deepCopy() ; mcpy.convertQuadraticCellsToLinear() ; mcpy.zipCoords()
        #
        f2=f.convertQuadraticCellsToLinear()
        f2.checkConsistencyLight()
        self.assertTrue(f2.getMesh().isEqual(mcpy,1e-12))
        self.assertTrue(f2.getArray().isEqual(DataArrayDouble([0,1,2,3]),1e-12))
        self.assertEqual(f2.getNbOfGaussLocalization(),1)
        gl=f2.getGaussLocalization(0)
        self.assertEqual(gl.getType(),NORM_TETRA4)
        self.assertTrue(DataArrayDouble(gl.getRefCoords()).isEqual(DataArrayDouble([0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0]),1e-12))
        self.assertTrue(DataArrayDouble(gl.getGaussCoords()).isEqual(DataArrayDouble(bbbb),1e-12))
        self.assertTrue(DataArrayDouble(gl.getWeights()).isEqual(DataArrayDouble(cccc),1e-12))
        self.assertEqual(f2.getName(),"field")
        self.assertEqual(f2.getTime(),[1.,2,3])
        pass

    def testDADCumSum1(self):
        d=DataArrayDouble([3.,2.,4.,5.])
        self.assertTrue(d.cumSum().isEqual(DataArrayDouble([0.,3.,5.,9.,14.]),1e-12))
        d2=DataArrayDouble([])
        self.assertTrue(d2.cumSum().isEqual(DataArrayDouble([0.]),1e-12))
        d.rearrange(2)
        self.assertRaises(InterpKernelException,d.cumSum)
        pass

    def testDAIFromLinkedListOfPairToList1(self):
        d=DataArrayInt([(5,7),(7,3),(3,12),(12,17)])
        zeRes=DataArrayInt([5,7,3,12,17])
        self.assertTrue(d.fromLinkedListOfPairToList().isEqual(zeRes))
        d.rearrange(1)
        self.assertRaises(InterpKernelException,d.fromLinkedListOfPairToList)
        d.rearrange(2)
        self.assertTrue(d.fromLinkedListOfPairToList().isEqual(zeRes))
        d2=DataArrayInt([(5,7)])
        self.assertTrue(d2.fromLinkedListOfPairToList().isEqual(DataArrayInt([5,7])))
        d3=DataArrayInt([(5,7),(7,3),(4,12),(12,17)])
        self.assertRaises(InterpKernelException,d3.fromLinkedListOfPairToList) # not a linked list of pair
        d4=DataArrayInt([(5,7),(7,3),(12,3),(12,17)])
        self.assertRaises(InterpKernelException,d4.fromLinkedListOfPairToList) # not a linked list of pair, but can be repaired !
        d4.sortEachPairToMakeALinkedList()
        self.assertTrue(d4.fromLinkedListOfPairToList().isEqual(zeRes))
        pass

    def testUMeshExplodeIntoEdges1(self):
        m=MEDCouplingCMesh() ; arr=DataArrayDouble(5) ; arr.iota() ; m.setCoords(arr,arr,arr) ; m=m.buildUnstructured()
        self.assertEqual(m.getMeshDimension(),3)
        a0,a1,a2,a3,a4=m.explodeIntoEdges()
        b0,b1,b2,b3,b4=m.explode3DMeshTo1D()
        self.assertTrue(a0.isEqual(b0,1e-12))
        self.assertTrue(a1.isEqual(b1)) ; self.assertTrue(a2.isEqual(b2)) ; self.assertTrue(a3.isEqual(b3)) ; self.assertTrue(a4.isEqual(b4))
        #
        m=MEDCouplingCMesh() ; arr=DataArrayDouble(5) ; arr.iota() ; m.setCoords(arr,arr) ; m=m.buildUnstructured()
        self.assertEqual(m.getMeshDimension(),2)
        a0,a1,a2,a3,a4=m.explodeIntoEdges()
        b0,b1,b2,b3,b4=m.buildDescendingConnectivity()
        self.assertTrue(a0.isEqual(b0,1e-12))
        self.assertTrue(a1.isEqual(b1)) ; self.assertTrue(a2.isEqual(b2)) ; self.assertTrue(a3.isEqual(b3)) ; self.assertTrue(a4.isEqual(b4))
        pass

    def testUMeshComputeEnlargedNeighborsOfNodes(self):
        m=MEDCouplingCMesh() ; arr=DataArrayDouble(4) ; arr.iota() ; m.setCoords(arr,arr) ; m=m.buildUnstructured()
        a,b=m.computeEnlargedNeighborsOfNodes()
        self.assertTrue(a.isEqual(DataArrayInt([1,4,5,0,2,4,5,6,1,3,5,6,7,2,6,7,0,1,5,8,9,0,1,2,4,6,8,9,10,1,2,3,5,7,9,10,11,2,3,6,10,11,4,5,9,12,13,4,5,6,8,10,12,13,14,5,6,7,9,11,13,14,15,6,7,10,14,15,8,9,13,8,9,10,12,14,9,10,11,13,15,10,11,14])))
        self.assertTrue(b.isEqual(DataArrayInt([0,3,8,13,16,21,29,37,42,47,55,63,68,71,76,81,84])))
        pass

    def testDAIfindIdsExt1(self):
        d=DataArrayInt([4,6,-2,3,7,0,10])
        self.assertTrue(d.findIdsGreaterOrEqualTo(3).isEqual(DataArrayInt([0,1,3,4,6])))
        self.assertTrue(d.findIdsGreaterThan(3).isEqual(DataArrayInt([0,1,4,6])))
        self.assertTrue(d.findIdsLowerThan(3).isEqual(DataArrayInt([2,5])))
        self.assertTrue(d.findIdsLowerOrEqualTo(3).isEqual(DataArrayInt([2,3,5])))
        pass

    def testDAFacto1(self):
        """Test focused of new wrapped methods for MEDCouplingFieldInt thanks to code factorization."""
        d=DataArrayDouble(7) ; d.iota()
        m=MEDCouplingUMesh.Build1DMeshFromCoords(d)
        f=MEDCouplingFieldInt(ON_CELLS) ; f.setMesh(m) ; arr=DataArrayInt(6) ; arr.iota() ; f.setArray(arr) ; f.checkConsistencyLight()
        f_0=f[::2] # test is here
        self.assertTrue(f_0.getArray().isEqual(DataArrayInt([0,2,4])))
        self.assertTrue(f_0.getMesh().isEqual(m[[0,2,4]],1e-12))
        #
        f2=MEDCouplingFieldInt(ON_NODES) ; f2.setMesh(m) ; arr=DataArrayInt(7) ; arr.iota() ; f2.setArray(arr) ; f2.checkConsistencyLight()
        f_1=f2[::2] # test is here
        self.assertTrue(f_1.getArray().isEqual(DataArrayInt([0,1,2,3,4,5])))
        m_1=m[[0,2,4]] ; m_1.zipCoords()
        self.assertTrue(f_1.getMesh().isEqual(m_1,1e-12))
        pass

    def testFieldFloatIsOnStage1(self):
        """ My first test with field int."""
        m=MEDCouplingCMesh()
        m.setName("mesh")
        arrX=DataArrayDouble([0,1,2,3])
        m.setCoords(arrX,arrX)
        f=MEDCouplingFieldFloat(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayFloat(8) ; arr.iota() ;f.setArray(arr)
        self.assertRaises(InterpKernelException,f.checkConsistencyLight)
        arr=DataArrayFloat(9) ; arr.iota() ;f.setArray(arr)
        f.checkConsistencyLight()
        f.setTimeUnit("ms")
        self.assertEqual(f.getTimeUnit(),"ms")
        f.setTime(3.2,5,6)
        a,b,c=f.getTime()
        self.assertEqual(b,5)
        self.assertEqual(c,6)
        self.assertEqual(a,3.2,12)
        pass

    def testFieldFloatIsOnStage2(self):
        """ Very important test to check that isEqual of MEDCouplingFieldFloat is OK !"""
        m1=MEDCouplingCMesh() ; m1.setCoords(DataArrayDouble([0,1,2,3]),DataArrayDouble([0,1,2,3,4]))
        m1=m1.buildUnstructured() ; m1.setName("mesh")
        f1=MEDCouplingFieldFloat(ON_CELLS) ; f1.setMesh(m1)
        arr1=DataArrayFloat([(0,1),(2,3),(4,5),(6,7),(8,9),(10,11),(12,13),(14,15),(16,17),(18,19),(20,21),(22,23)]) ; arr1.setInfoOnComponents(["aa","bbb"])
        f1.setArray(arr1) ; f1.setName("f1") ; f1.setTime(2.,3,4)
        #
        m2=MEDCouplingCMesh() ; m2.setCoords(DataArrayDouble([0,1,2,3]),DataArrayDouble([0,1,2,3,4]))
        m2=m2.buildUnstructured() ; m2.setName("mesh")
        f2=MEDCouplingFieldFloat(ON_CELLS) ; f2.setMesh(m2)
        arr2=DataArrayFloat([(0,1),(2,3),(4,5),(6,7),(8,9),(10,11),(12,13),(14,15),(16,17),(18,19),(20,21),(22,23)]) ; arr2.setInfoOnComponents(["aa","bbb"])
        f2.setArray(arr2) ; f2.setName("f1") ; f2.setTime(2.,3,4)
        #
        self.assertTrue(f1.isEqual(f2,1e-12,0.))
        f1.getArray()[:]*=2
        self.assertTrue(not f1.isEqual(f2,1e-12,0.))
        self.assertTrue(not f1.isEqualWithoutConsideringStr(f2,1e-12,0.))
        f1.getArray()[:]/=2
        self.assertTrue(f1.isEqual(f2,1e-12,0.))
        #
        f1.setName("F1")
        self.assertTrue(not f1.isEqual(f2,1e-12,0.))
        f1.setName("f1")
        self.assertTrue(f1.isEqual(f2,1e-12,0.))
        #
        f1.getArray().setInfoOnComponents(["aa","bbbb"])
        self.assertTrue(not f1.isEqual(f2,1e-12,0.))
        self.assertTrue(f1.isEqualWithoutConsideringStr(f2,1e-12,0.))
        f1.getArray().setInfoOnComponents(["aa","bbb"])
        self.assertTrue(f1.isEqual(f2,1e-12,0.))
        #
        f3=f2.deepCopy()
        self.assertTrue(f1.isEqual(f3,1e-12,0.))
        #
        for fd,expected in ((ON_NODES,False),(ON_CELLS,True)):
            f4=MEDCouplingFieldFloat(fd) ; f4.setMesh(m2) ; f4.setTime(2.,3,4)
            arr4=DataArrayFloat([(0,1),(2,3),(4,5),(6,7),(8,9),(10,11),(12,13),(14,15),(16,17),(18,19),(20,21),(22,23)]) ; arr4.setInfoOnComponents(["aa","bbb"])
            f4.setArray(arr4) ; f4.setName("f1")
            self.assertEqual(f1.isEqual(f4,1e-12,0.),expected)
            pass
        pass

    def testLTGTDAD1(self):
        d=DataArrayDouble(10) ; d.iota()
        self.assertTrue(d.findIdsLowerThan(0).empty())
        self.assertTrue(d.findIdsLowerThan(1).isEqual(DataArrayInt([0])))
        d-=5.
        self.assertTrue(d.findIdsStrictlyNegative().isEqual(DataArrayInt([0,1,2,3,4])))
        self.assertTrue(d.findIdsGreaterThan(0.).isEqual(DataArrayInt([6,7,8,9])))
        self.assertTrue(d.convertToFloatArr().isEqual(DataArrayFloat([-5,-4,-3,-2,-1,0,1,2,3,4]),1e-7))
        self.assertTrue(d.convertToFloatArr().convertToDblArr().isEqual(d,1e-12))
        pass

    def testMapII1(self):
        """ Test optimized maps for renumbering. Typical usage local to global in parallel mode"""
        d=DataArrayInt([1003,1007])
        m=d.invertArrayN2O2O2NOptimized()
        d2=DataArrayInt([1003,1003,1007,1003,1007])
        d2.transformWithIndArr(m)
        self.assertTrue(d2.isEqual(DataArrayInt([0,0,1,0,1])))
        pass

    def testDAICheckUniformAndGuess1(self):
        d=DataArrayInt([3,3],1,2)
        self.assertRaises(InterpKernelException,d.checkUniformAndGuess)# non single compo
        d=DataArrayInt([])
        self.assertRaises(InterpKernelException,d.checkUniformAndGuess)# empty
        d=DataArrayInt()
        self.assertRaises(InterpKernelException,d.checkUniformAndGuess)# non allocated
        d=DataArrayInt([3,3,3])
        self.assertEqual(3,d.checkUniformAndGuess())
        d=DataArrayInt([7])
        self.assertEqual(7,d.checkUniformAndGuess())
        d=DataArrayInt([3,4,3])
        self.assertRaises(InterpKernelException,d.checkUniformAndGuess)# non uniform
        pass

    def testUMComputePlaneEquationOf3DFaces1(self):
        """ Consequence of an invalid traduction of matrix inversion transposition."""
        m=MEDCoupling1SGTUMesh("msh",NORM_QUAD4)
        m.setCoords(DataArrayDouble([(0,0,0),(1,0,0),(2,0,0),(0,2,0),(1,2,0),(2,2,0),(0,4,0),(1,4,0),(2,4,0),(0,0,3),(1,0,3),(2,0,3),(0,2,3),(1,2,3),(2,2,3),(0,4,3),(1,4,3),(2,4,3)]))
        m.setNodalConnectivity(DataArrayInt([0,1,4,3,9,12,13,10,0,9,10,1,1,10,13,4,4,13,12,3,3,12,9,0,1,2,5,4,10,13,14,11,1,10,11,2,2,11,14,5,5,14,13,4,3,4,7,6,12,15,16,13,4,13,16,7,7,16,15,6,6,15,12,3,4,5,8,7,13,16,17,14,5,14,17,8,8,17,16,7]))
        m=m.buildUnstructured()
        ref=DataArrayDouble([(0,0,1,0),(0,0,1,-3),(0,1,0,0),(1,0,0,-1),(0,1,0,-2),(1,0,0,0),(0,0,1,0),(0,0,1,-3),(0,1,0,0),(1,0,0,-2),(0,1,0,-2),(0,0,1,0),(0,0,1,-3),(1,0,0,-1),(0,1,0,-4),(1,0,0,0),(0,0,1,0),(0,0,1,-3),(1,0,0,-2),(0,1,0,-4)])
        res=m.computePlaneEquationOf3DFaces()
        self.assertTrue(res.isEqual(ref,1e-12))
        pass

    def testBugInComputationOfEqOfPlane1(self):
        coo=DataArrayDouble([-1.0, 1.0, -0.3872983455657959, -1.0, 1.0, 0.3872983455657959, -1.0, 1.0, 0.693649172782898, 1.0, 1.0, 0.693649172782898, 1.0, 1.0, 0.3872983455657959, 1.0, 1.0, -0.3872983455657959],6,3)
        m=MEDCouplingUMesh("",2)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_POLYGON,[0,1,2,3,4,5])
        self.assertTrue(m.computePlaneEquationOf3DFaces().isEqual(DataArrayDouble([0,1,0,-1],1,4),1e-12))
        pass

    def testSimplifyPolyhedra(self):
        mesh = MEDCouplingUMesh('mesh', 3)
        coo = DataArrayDouble([(-0.01225,-0.0212176,0.02),(-0.00634107,-0.0236652,0.02),(1.50019e-18,-0.0245,0.02),(0.00634107,-0.0236652,0.02),(0.01225,-0.0212176,0.02),(-0.0153864,-0.02665,0),(-0.00714085,-0.02665,0),(1.63184e-18,-0.02665,0),(0.00714085,-0.02665,0),(0.0153864,-0.02665,0),(-0.00714085,-0.02665,0.0101475),(1.63184e-18,-0.02665,0.013145),(0.00714085,-0.02665,0.0101475),(-0.013,-0.0225167,0.02),(-0.0067293,-0.0251141,0.02),(1.59204e-18,-0.026,0.02),(0.0067293,-0.0251141,0.02),(0.013,-0.0225167,0.02),(-0.0161658,-0.028,0),(-0.00750258,-0.028,0),(1.71451e-18,-0.028,0),(0.00750258,-0.028,0),(0.0161658,-0.028,0),(-0.00750258,-0.028,0.0105625),(1.71451e-18,-0.028,0.0136825),(0.00750258,-0.028,0.0105625)])
        mesh.setCoords(coo)
        c = DataArrayInt([31, 13, 14, 15, 16, 17, 4, 3, 2, 1, 0, -1, 18, 5, 6, 7, 8, 9, 22, 21, 20, 19, -1, 19, 23, 18, -1, 23, 14, 13, 18, -1, 20, 24, 23, 19, -1, 24, 15, 14, 23, -1, 21, 25, 24, 20, -1, 25, 16, 15, 24, -1, 22, 25, 21, -1, 22, 17, 16, 25, -1, 9, 4, 17, 22, -1, 8, 12, 9, -1, 12, 3, 4, 9, -1, 7, 11, 12, 8, -1, 11, 2, 3, 12, -1, 6, 10, 11, 7, -1, 10, 1, 2, 11, -1, 5, 10, 6, -1, 5, 0, 1, 10, -1, 18, 13, 0, 5])
        cI = DataArrayInt([0, 108])
        mesh.setConnectivity(c, cI)
        mesh.simplifyPolyhedra(1.0e-8)
        c, cI = mesh.getNodalConnectivity(), mesh.getNodalConnectivityIndex()
        tgt_c = DataArrayInt([31, 23, 18, 19, 20, 21, 22, 25, 24, -1, 12, 9, 8, 7, 6, 5, 10, 11, -1, 13, 14, 15, 16, 17, 4, 3, 2, 1, 0, -1, 18, 5, 6, 7, 8, 9, 22, 21, 20, 19, -1, 23, 14, 13, 18, -1, 24, 15, 14, 23, -1, 25, 16, 15, 24, -1, 22, 17, 16, 25, -1, 9, 4, 17, 22, -1, 12, 3, 4, 9, -1, 11, 2, 3, 12, -1, 10, 1, 2, 11, -1, 5, 0, 1, 10, -1, 18, 13, 0, 5])
        tgt_cI = DataArrayInt([0, 90])
        self.assertEqual(c.getValues(), tgt_c.getValues())
        self.assertEqual(cI.getValues(), tgt_cI.getValues())
        pass

    pass

if __name__ == '__main__':
    unittest.main()
