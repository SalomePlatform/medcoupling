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

from MEDLoader import *
from math import pi,e,sqrt

class MEDLoaderDataForTest:
    def build1DMesh_1(cls):
        coords=[ 0.0, 0.3, 0.75, 1.0, 1.4, 1.3 ]
        conn=[ 0,1, 1,2, 2,3 , 3,4,5]
        mesh=MEDCouplingUMesh.New();
        mesh.setName("1DMesh_1");
        mesh.setMeshDimension(1);
        mesh.allocateCells(4);
        mesh.insertNextCell(NORM_SEG2,2,conn[0:2])
        mesh.insertNextCell(NORM_SEG2,2,conn[2:4])
        mesh.insertNextCell(NORM_SEG2,2,conn[4:6])
        mesh.insertNextCell(NORM_SEG3,3,conn[6:9])
        mesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(coords,6,1);
        myCoords.setInfoOnComponent(0,"tototototototot (m*m*m*m*m*m*m*m)");
        mesh.setCoords(myCoords);
        return mesh;

    def build2DCurveMesh_1(cls):
        coords=[ 0.0,0.0, 0.3,0.3, 0.75,0.75, 1.0,1.0, 1.4,1.4, 1.3,1.3 ]
        conn=[ 0,1, 1,2, 2,3 , 3,4,5]
        mesh=MEDCouplingUMesh.New();
        mesh.setName("2DCurveMesh_1");
        mesh.setMeshDimension(1);
        mesh.allocateCells(4);
        mesh.insertNextCell(NORM_SEG2,2,conn[0:2])
        mesh.insertNextCell(NORM_SEG2,2,conn[2:4])
        mesh.insertNextCell(NORM_SEG2,2,conn[4:6])
        mesh.insertNextCell(NORM_SEG3,3,conn[6:9])
        mesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(coords,6,2);
        mesh.setCoords(myCoords);
        return mesh;

    def build2DMesh_1(cls):
        targetCoords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7,
                      -0.05,0.95, 0.2,1.2, 0.45,0.95]
        targetConn=[1,4,2, 4,5,2, 6,10,8,9,11,7, 0,3,4,1, 6,7,4,3, 7,8,5,4]
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.setMeshDimension(2);
        targetMesh.allocateCells(6);
        targetMesh.setName("2DMesh_1");
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[0:3])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[3:6])
        targetMesh.insertNextCell(NORM_TRI6,6,targetConn[6:12])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[12:16])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[16:20])
        targetMesh.insertNextCell(NORM_POLYGON,4,targetConn[20:24])
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,12,2);
        myCoords.setInfoOnComponent(0,"tototototototot (m)");
        myCoords.setInfoOnComponent(1,"energie (kW)");
        targetMesh.setCoords(myCoords)
        return targetMesh;

    def build2DMesh_2(cls):
        targetCoords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7,
                      -0.05,0.95, 0.2,1.2, 0.45,0.95]
        targetConn=[1,4,2, 4,5,2, 6,10,8,9,11,7, 0,3,4,1, 6,7,4,3, 7,8,5,4]
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.setMeshDimension(2);
        targetMesh.allocateCells(5);
        targetMesh.setName("2DMesh_2");
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[0:3])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[3:6])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[12:16])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[16:20])
        targetMesh.insertNextCell(NORM_TRI6,6,targetConn[6:12])
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,12,2);
        myCoords.setInfoOnComponent(0,"toto (m)");
        myCoords.setInfoOnComponent(1,"energie (kW)");
        targetMesh.setCoords(myCoords);
        return targetMesh;

    def build3DMesh_1(cls):
        coords=[0.,0.,0., 1.,1.,0., 1.,1.25,0., 0.,1.,0., 1.,1.5,0., 2.,0.,0., 2.,1.,0., 1.,2.,0., 0.,2.,0., 3.,1.,0.,
                3.,2.,0., 0.,1.,0., 1.,3.,0., 2.,2.,0., 2.,3.,0.,
                0.,0.,1., 1.,1.,1., 1.,1.25,1., 0.,1.,1., 1.,1.5,1., 2.,0.,1., 2.,1.,1., 1.,2.,1., 0.,2.,1., 3.,1.,1.,
                3.,2.,1., 0.,1.,1., 1.,3.,1., 2.,2.,1., 2.,3.,1.,
                0.,0.,2., 1.,1.,2., 1.,1.25,2., 0.,1.,2., 1.,1.5,2., 2.,0.,2., 2.,1.,2., 1.,2.,2., 0.,2.,2., 3.,1.,2.,
        3.,2.,2., 0.,1.,2., 1.,3.,2., 2.,2.,2., 2.,3.,2.,
                0.,0.,3., 1.,1.,3., 1.,1.25,3., 0.,1.,3., 1.,1.5,3., 2.,0.,3., 2.,1.,3., 1.,2.,3., 0.,2.,3., 3.,1.,3.,
                3.,2.,3., 0.,1.,3., 1.,3.,3., 2.,2.,3., 2.,3.,3.]
        conn=[
            # 0
            0,11,1,3,15,26,16,18,   1,2,4,7,13,6,-1,1,16,21,6,-1,6,21,28,13,-1,13,7,22,28,-1,7,4,19,22,-1,4,2,17,19,-1,2,1,16,17,-1,16,21,28,22,19,17,
            1,6,5,3,16,21,20,18,   13,10,9,6,28,25,24,21,
            11,8,7,4,2,1,-1,11,26,16,1,-1,1,16,17,2,-1,2,17,19,4,-1,4,19,22,7,-1,7,8,23,22,-1,8,11,26,23,-1,26,16,17,19,22,23,
            7,12,14,13,22,27,29,28,
            # 1
            15,26,16,18,30,41,31,33,   16,17,19,22,28,21,-1,16,31,36,21,-1,21,36,43,28,-1,28,22,37,43,-1,22,19,34,37,-1,19,17,32,34,-1,17,16,31,32,-1,31,36,43,37,34,32,
            16,21,20,18,31,36,35,33,   28,25,24,21,43,40,39,36,
            26,23,22,19,17,16,-1,26,41,31,16,-1,16,31,32,17,-1,17,32,34,19,-1,19,34,37,22,-1,22,23,38,37,-1,23,26,41,38,-1,41,31,32,34,37,38,
            22,27,29,28,37,42,44,43,
            # 2
            30,41,31,33,45,56,46,48,  31,32,34,37,43,36,-1,31,46,51,36,-1,36,51,58,43,-1,43,37,52,58,-1,37,34,49,52,-1,34,32,47,49,-1,32,31,46,47,-1,46,51,58,52,49,47,
            31,36,35,33,46,51,50,48,  43,40,39,36,58,55,54,51,
            41,38,37,34,32,31,-1,41,56,46,31,-1,31,46,47,32,-1,32,47,49,34,-1,34,49,52,37,-1,37,38,53,52,-1,38,41,56,53,-1,56,46,47,49,52,53,
            37,42,44,43,52,57,59,58]
        #
        ret=MEDCouplingUMesh.New();
        ret.setName("3DMesh_1");
        ret.setMeshDimension(3);
        ret.allocateCells(18);
        #
        ret.insertNextCell(NORM_HEXA8,8,conn[0:8])
        ret.insertNextCell(NORM_HEXA8,8,conn[51:59])
        ret.insertNextCell(NORM_HEXA8,8,conn[59:67])
        ret.insertNextCell(NORM_HEXA8,8,conn[110:118])
        #
        ret.insertNextCell(NORM_HEXA8,8,conn[118:126])
        ret.insertNextCell(NORM_HEXA8,8,conn[169:177])
        ret.insertNextCell(NORM_HEXA8,8,conn[177:185])
        ret.insertNextCell(NORM_HEXA8,8,conn[228:236])
        #
        ret.insertNextCell(NORM_HEXA8,8,conn[236:244])
        ret.insertNextCell(NORM_HEXA8,8,conn[287:295])
        ret.insertNextCell(NORM_HEXA8,8,conn[295:303])
        ret.insertNextCell(NORM_HEXA8,8,conn[346:354])
        #
        ret.insertNextCell(NORM_POLYHED,43,conn[8:51])
        ret.insertNextCell(NORM_POLYHED,43,conn[67:110])
        ret.insertNextCell(NORM_POLYHED,43,conn[126:169])
        ret.insertNextCell(NORM_POLYHED,43,conn[185:228])
        ret.insertNextCell(NORM_POLYHED,43,conn[244:287])
        ret.insertNextCell(NORM_POLYHED,43,conn[303:346])
        #
        ret.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(coords,60,3);
        myCoords.setInfoOnComponent(0,"titi (m)");
        myCoords.setInfoOnComponent(1,"density power (MW/m^3)");
        myCoords.setInfoOnComponent(2,"t (kW)");
        ret.setCoords(myCoords);
        return ret;
    
    def build3DSurfMesh_1(cls):
        targetCoords=[-0.3,-0.3,-0.3, 0.2,-0.3,-0.3, 0.7,-0.3,-0.3, -0.3,0.2,-0.3, 0.2,0.2,-0.3, 0.7,0.2,-0.3, -0.3,0.7,-0.3, 0.2,0.7,-0.3, 0.7,0.7,-0.3
                      ,-0.05,0.95,-0.3, 0.2,1.2,-0.3, 0.45,0.95,-0.3]
        targetConn=[1,4,2, 4,5,2, 6,10,8,9,11,7, 0,3,4,1, 6,7,4,3, 7,8,5,4]
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.setMeshDimension(2);
        targetMesh.allocateCells(6);
        targetMesh.setName("3DSurfMesh_1");
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[0:3])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[3:6])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[12:16])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[16:20])
        targetMesh.insertNextCell(NORM_TRI6,6,targetConn[6:12])
        targetMesh.insertNextCell(NORM_POLYGON,4,targetConn[20:24])
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,12,3);
        myCoords.setInfoOnComponent(0,"toto (m)");
        myCoords.setInfoOnComponent(2,"ff (km)");#component 1 is not set for test
        targetMesh.setCoords(myCoords);
        return targetMesh;
    
    def build3DMesh_2(cls):
        m3dsurfBase=MEDLoaderDataForTest.build3DSurfMesh_1();
        numbers=[0,1,2,3,5]
        m3dsurf=m3dsurfBase.buildPartOfMySelf(numbers,False);
        m1dBase=MEDLoaderDataForTest.build1DMesh_1();
        numbers2=[0,1,2,3]
        m1d=m1dBase.buildPartOfMySelf(numbers2,False);
        m1d.changeSpaceDimension(3);
        vec=[0.,1.,0.]
        pt=[0.,0.,0.]
        m1d.rotate(pt,vec,-pi/2.);
        ret=m3dsurf.buildExtrudedMesh(m1d,0);
        return ret;

    def buildVecFieldOnCells_1(cls):
        mesh=MEDLoaderDataForTest.build3DSurfMesh_1();
        nbOfCells=mesh.getNumberOfCells();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setName("VectorFieldOnCells");
        f1.setMesh(mesh);
        array=DataArrayDouble.New();
        arr1=[0.,10.,20.,1.,11.,21.,2.,12.,22.,3.,13.,23.,4.,14.,24.,5.,15.,25.]
        array.setValues(arr1,nbOfCells,3);
        array.setInfoOnComponent(0,"power (MW/m^3)");
        array.setInfoOnComponent(1,"density (g/cm^3)");
        array.setInfoOnComponent(2,"temperature (K)");
        f1.setArray(array);
        tmp=array.getPointer();
        f1.setTime(2.,0,1);
        f1.checkCoherency();
        return f1;

    def buildVecFieldOnNodes_1(cls):
        mesh=MEDLoaderDataForTest.build3DSurfMesh_1();
        nbOfNodes=mesh.getNumberOfNodes();
        f1=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME);
        f1.setName("VectorFieldOnNodes");
        f1.setMesh(mesh);
        array=DataArrayDouble.New();
        f1.setArray(array);
        arr1=[70.,80.,90.,71.,81.,91.,72.,82.,92.,73.,83.,93.,74.,84.,94.,75.,85.,95.,
        1000.,10010.,10020.,1001.,10011.,10021.,1002.,10012.,10022.,1003.,10013.,10023.,1004.,10014.,10024.,1005.,10015.,10025.]
        array.setValues(arr1,nbOfNodes,3);
        array.setInfoOnComponent(0,"power (MW/m^3)");
        array.setInfoOnComponent(1,"density (g/cm^3)");
        array.setInfoOnComponent(2,"temperature (K)");
        f1.setTime(2.12,2,3);
        f1.checkCoherency();
        return f1;

    def buildVecFieldOnGauss_1(cls):
        _a=0.446948490915965;
        _b=0.091576213509771;
        _p1=0.11169079483905;
        _p2=0.0549758718227661;
        refCoo1=[ 0.,0., 1.,0., 0.,1. ]
        gsCoo1=[ 2*_b-1, 1-4*_b, 2*_b-1, 2.07*_b-1, 1-4*_b,
                 2*_b-1, 1-4*_a, 2*_a-1, 2*_a-1, 1-4*_a, 2*_a-1, 2*_a-1 ];
        wg1=[ 4*_p2, 4*_p2, 4*_p2, 4*_p1, 4*_p1, 4*_p1 ]
        _refCoo1=refCoo1;
        _gsCoo1=gsCoo1;
        _wg1=wg1;
        m=MEDLoaderDataForTest.build2DMesh_2();
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,ONE_TIME);
        f.setTime(3.14,1,5);
        f.setMesh(m);
        f.setGaussLocalizationOnType(NORM_TRI3,_refCoo1,_gsCoo1,_wg1);
        refCoo2=[-1.0,1.0, -1.0,-1.0, 1.0,-1.0, -1.0,0.0, 0.0,-1.0, 0.0,0.0 ]
        _refCoo2=refCoo2;
        _gsCoo1=_gsCoo1[0:6];
        _gsCoo2=_gsCoo1
        _wg1=_wg1[0:3];
        _wg2=_wg1
        refCoo3=[ 0.,0., 1.,0., 1.,1., 0.,1. ]
        _refCoo3=refCoo3;
        _gsCoo1=_gsCoo1[0:4];
        _wg1=_wg1[0:2];
        f.setGaussLocalizationOnType(NORM_QUAD4,_refCoo3,_gsCoo1,_wg1);
        f.setGaussLocalizationOnType(NORM_TRI6,_refCoo2,_gsCoo2,_wg2);
        array=DataArrayDouble.New();
        array.alloc(19,2);
        ptr=array.getPointer();
        for i in xrange(19*2):
            array.setIJ(0,i,float(i+7));
            pass
        f.setArray(array);
        f.setName("MyFirstFieldOnGaussPoint");
        array.setInfoOnComponent(0,"power (MW/m^3)");
        array.setInfoOnComponent(1,"density");
        f.checkCoherency();
        return f;

    def buildVecFieldOnGaussNE_1(cls):
        m=MEDLoaderDataForTest.build2DMesh_2();
        f=MEDCouplingFieldDouble.New(ON_GAUSS_NE,ONE_TIME);
        f.setTime(3.14,1,5);
        f.setMesh(m);
        array=DataArrayDouble.New();
        array.alloc(20,2);
        for i in xrange(2*20):
            array.setIJ(0,i,float(i+8));
        f.setArray(array);
        array.setInfoOnComponent(0,"power (W)");
        array.setInfoOnComponent(1,"temperature");
        f.setName("MyFieldOnGaussNE");
        f.checkCoherency();
        return f;
    
    build1DMesh_1=classmethod(build1DMesh_1)
    build2DCurveMesh_1=classmethod(build2DCurveMesh_1)
    build2DMesh_1=classmethod(build2DMesh_1)
    build2DMesh_2=classmethod(build2DMesh_2)
    build3DMesh_1=classmethod(build3DMesh_1)
    build3DSurfMesh_1=classmethod(build3DSurfMesh_1)
    build3DMesh_2=classmethod(build3DMesh_2)
    buildVecFieldOnCells_1=classmethod(buildVecFieldOnCells_1)
    buildVecFieldOnNodes_1=classmethod(buildVecFieldOnNodes_1)
    buildVecFieldOnGauss_1=classmethod(buildVecFieldOnGauss_1)
    buildVecFieldOnGaussNE_1=classmethod(buildVecFieldOnGaussNE_1)
    pass
