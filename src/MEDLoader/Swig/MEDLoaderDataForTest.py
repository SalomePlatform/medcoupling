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
# Author : Anthony Geay (EDF R&D)

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
        myCoords.setInfoOnComponent(0,"tototototototot [m*m*m*m*m*m*m*m]");
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
        myCoords.setInfoOnComponent(0,"tototototototot [m]");
        myCoords.setInfoOnComponent(1,"energie [kW]");
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
        myCoords.setInfoOnComponent(0,"toto [m]");
        myCoords.setInfoOnComponent(1,"energie [kW]");
        targetMesh.setCoords(myCoords);
        return targetMesh;

    #this mesh has several cells duplicated ! it is not beautiful but efficient to test file WR.
    def build2DMesh_3(cls):
        targetCoords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7,
                      -0.05,0.95, 0.2,1.2, 0.45,0.95]
        targetConn=[1,4,2, 4,5,2, 6,10,8,9,11,7, 0,3,4,1, 6,7,4,3, 7,8,5,4]
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.setMeshDimension(2);
        targetMesh.allocateCells(13);
        targetMesh.setName("2DMesh_3");
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[0:3])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[3:6])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[0:3])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[3:6])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[0:3])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[3:6])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[12:16])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[16:20])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[12:16])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[16:20])
        targetMesh.insertNextCell(NORM_TRI6,6,targetConn[6:12])
        targetMesh.insertNextCell(NORM_TRI6,6,targetConn[6:12])
        targetMesh.insertNextCell(NORM_TRI6,6,targetConn[6:12])
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,12,2);
        myCoords.setInfoOnComponent(0,"toto [m]");
        myCoords.setInfoOnComponent(1,"energie [kW]");
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
        myCoords.setInfoOnComponent(0,"titi [m]");
        myCoords.setInfoOnComponent(1,"density power [MW/m^3]");
        myCoords.setInfoOnComponent(2,"t [kW]");
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
        myCoords.setInfoOnComponent(0,"toto [m]");
        myCoords.setInfoOnComponent(2,"ff [km]");#component 1 is not set for test
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

    def buildMultiLevelMesh_1(cls):
        coo=[10.,0.,10.,1.25,10.,2.5,10.,3.75,10.,5.,8.75,0.,8.75,1.25,8.75,2.5,8.75,3.75,8.75,5.,7.5,0.,7.5,1.25,7.5,2.5,7.5,3.75,7.5,5.,6.25,0.,6.25,1.25,6.25,2.5,6.25,3.75,6.25,5.,5.,0.,5.,1.25,5.,2.5,5.,3.75,5.,5.,3.75,0.,3.75,1.25,3.75,2.5,3.75,3.75,3.75,5.,2.5,0.,2.5,1.25,2.5,2.5,2.5,3.75,2.5,5.,1.25,0.,1.25,1.25,1.25,2.5,1.25,3.75,1.25,5.,0.,1.25,0.,2.5,0.,3.75,0.,5.,0.,0.,0.,5.,10.,5.,0.,10.,10.,10.,5.,5.,5.,5.,5.,10.,5.,10.,0.625,5.,1.25,5.,1.875,5.,2.5,5.,3.125,5.,3.75,5.,4.375,5.,5.,6.25,5.,7.5,5.,8.75,4.375,10.,3.75,10.,3.125,10.,2.5,10.,1.875,10.,1.25,10.,0.625,10.,0.,8.75,0.,7.5,0.,6.25,4.375,6.25,4.375,7.5,4.375,8.75,3.75,6.25,3.75,7.5,3.75,8.75,3.125,6.25,3.125,7.5,3.125,8.75,2.5,6.25,2.5,7.5,2.5,8.75,1.875,6.25,1.875,7.5,1.875,8.75,1.25,6.25,1.25,7.5,1.25,8.75,0.625,6.25,0.625,7.5,0.625,8.75,5.625,5.,6.25,5.,6.875,5.,7.5,5.,8.125,5.,8.75,5.,9.375,5.,10.,6.25,10.,7.5,10.,8.75,9.375,10.,8.75,10.,8.125,10.,7.5,10.,6.875,10.,6.25,10.,5.625,10.,5.,8.75,5.,7.5,5.,6.25,9.375,6.25,9.375,7.5,9.375,8.75,8.75,6.25,8.75,7.5,8.75,8.75,8.125,6.25,8.125,7.5,8.125,8.75,7.5,6.25,7.5,7.5,7.5,8.75,6.875,6.25,6.875,7.5,6.875,8.75,6.25,6.25,6.25,7.5,6.25,8.75,5.625,6.25,5.625,7.5,5.625,8.75]
        coo2=DataArrayDouble.New()
        coo2.setValues(coo,135,2)
        coo2=coo2.changeNbOfComponents(3,0.)
        coo2.setInfoOnComponent(0,"X [INCONNUE]")
        coo2.setInfoOnComponent(1,"Y [INCONNUE]")
        coo2.setInfoOnComponent(2,"Z [INCONNUE]")
        c2tri=[0,1,6,0,6,5,1,2,6,2,7,6,2,3,8,2,8,7,3,4,8,4,9,8,5,6,11,5,11,10,6,7,11,7,12,11,7,8,13,7,13,12,8,9,13,9,14,13,10,11,16,10,16,15,11,12,16,12,17,16,12,13,18,12,18,17,13,14,18,14,19,18,15,16,21,15,21,20,16,17,21,17,22,21,17,18,23,17,23,22,18,19,23,19,24,23,20,21,26,20,26,25,21,22,26,22,27,26,22,23,28,22,28,27,23,24,28,24,29,28,25,26,31,25,31,30,26,27,31,27,32,31,27,28,33,27,33,32,28,29,33,29,34,33,30,31,36,30,36,35,31,32,36,32,37,36,32,33,38,32,38,37,33,34,38,34,39,38,35,36,40,35,40,44,36,37,40,37,41,40,37,38,42,37,42,41,38,39,42,39,43,42]
        c2quad4=[46,101,114,100,101,102,115,114,102,103,116,115,103,48,104,116,100,114,117,99,114,115,118,117,115,116,119,118,116,104,105,119,99,117,120,98,117,118,121,120,118,119,122,121,119,105,106,122,98,120,123,97,120,121,124,123,121,122,125,124,122,106,107,125,97,123,126,96,123,124,127,126,124,125,128,127,125,107,108,128,96,126,129,95,126,127,130,129,127,128,131,130,128,108,109,131,95,129,132,94,129,130,133,132,130,131,134,133,131,109,110,134,94,132,113,50,132,133,112,113,133,134,111,112,134,110,51,111,49,60,73,59,60,61,74,73,61,62,75,74,62,52,63,75,59,73,76,58,73,74,77,76,74,75,78,77,75,63,64,78,58,76,79,57,76,77,80,79,77,78,81,80,78,64,65,81,57,79,82,56,79,80,83,82,80,81,84,83,81,65,66,84,56,82,85,55,82,83,86,85,83,84,87,86,84,66,67,87,55,85,88,54,85,86,89,88,86,87,90,89,87,67,68,90,54,88,91,53,88,89,92,91,89,90,93,92,90,68,69,93,53,91,72,45,91,92,71,72,92,93,70,71,93,69,47,70]
        m2=MEDCouplingUMesh.New("ma",2)
        m2.setCoords(coo2)
        m2.allocateCells(128)
        nbTri = len(c2tri) // 3
        for i in range(nbTri):
            m2.insertNextCell(NORM_TRI3,3,c2tri[3*i:3*i+3])
            pass
        nbQua = len(c2quad4) // 4
        for i in range(nbQua):
            m2.insertNextCell(NORM_QUAD4,4,c2quad4[4*i:4*i+4])
            pass
        m2.finishInsertingCells()
        m2.setDescription("CREE PAR CODE_ASTER")
        m1=MEDCouplingUMesh.New("ma",1)
        m1.setCoords(coo2)
        c1seg=[0,1,1,2,2,3,3,4,4,9,9,14,14,19,19,24,24,29,29,34,34,39,39,43,43,42,42,41,41,40,40,44,44,35,35,30,30,25,25,20,20,15,15,10,10,5,5,0,43,39,39,34,34,29,29,24,24,19,19,14,14,9,9,4,45,53,53,54,54,55,55,56,56,57,57,58,58,59,59,49,49,60,60,61,61,62,62,52,52,63,63,64,64,65,65,66,66,67,67,68,68,69,69,47,47,70,70,71,71,72,72,45,50,94,94,95,95,96,96,97,97,98,98,99,99,100,100,46,46,101,101,102,102,103,103,48,48,104,104,105,105,106,106,107,107,108,108,109,109,110,110,51,51,111,111,112,112,113,113,50]
        m1.allocateCells(80)
        for i in range(80):
            m1.insertNextCell(NORM_SEG2,2,c1seg[2*i:2*i+2])
            pass
        m1.finishInsertingCells()
        m1.setDescription("CREE PAR CODE_ASTER")
        m0=MEDCouplingUMesh.New("ma",0)
        m0.setCoords(coo2)
        c0pt=[44,0,47,48]
        m0.allocateCells(4)
        for i in range(4):
            m0.insertNextCell(NORM_POINT1,1,[c0pt[i]])
            pass
        m0.finishInsertingCells()
        f2=DataArrayInt.New()
        f2.alloc(128,1)
        f2[:64]=-1
        f2[64:96]=-2
        f2[96:]=-3
        f1=DataArrayInt.New()
        f1.alloc(80,1)
        f1[:4]=-8
        f1[4:12]=-9
        f1[12:16]=-10
        f1[16:24]=-11
        f1[24:28]=-12
        f1[28:32]=-13
        f1[32:40]=-14
        f1[40:44]=-15
        f1[44:52]=-16
        f1[52:56]=-17
        f1[56:64]=-18
        f1[64:68]=-19
        f1[68:76]=-20
        f1[76:]=-21
        f0=DataArrayInt.New()
        f0.setValues([-4,-5,-6,-7],4,1)
        p=DataArrayInt.New()
        p.alloc(135,1)
        p.fillWithZero()
        p1=DataArrayInt.New()
        p1.alloc(13,1)
        p1.iota(1)
        p[[0,4,24,43,44,45,46,47,48,49,50,51,52]]=p1
        n2=DataArrayInt.New()
        n2.alloc(128,1)
        n2.iota(1)
        n1=DataArrayInt.New()
        n1.alloc(80,1)
        n1.iota(133)
        n0=DataArrayInt.New()
        n0.alloc(4,1)
        n0.iota(129)
        fns=['A1A2____________________________', 'A1______________________________', 'A2A4____________________________', 'A2______________________________', 'A3A1____________________________', 'A3C5____________________________', 'A3______________________________', 'A4A3____________________________', 'A4______________________________', 'B1C1____________________________', 'B1______________________________', 'B2B4____________________________', 'B2______________________________', 'B3B1____________________________', 'B3______________________________', 'B4C3____________________________', 'B4______________________________', 'C1C4____________________________', 'C1______________________________', 'C2B2____________________________', 'C2______________________________', 'C3C2____________________________', 'C3______________________________', 'C4B3____________________________', 'C4______________________________', 'C5A4____________________________', 'C5______PMMA____________________', 'FAMILLE_ZERO', 'MESH____APPS____AP1_____________', 'MESH____APPS____AP2_____________', 'MESH____APPS____AP3_____________', 'MESH____APPS____AP4_____________', 'MESH____DALQ1___DALLE___________', 'MESH____DALQ2___DALLE___________', 'MESH____DALT3___DALLE___________']
        fids=[-11, 5, -8, 1, -10, -12, 4, -9, 2, -14, 6, -19, 7, -17, 8, -20, 9, -15, 10, -18, 11, -21, 12, -16, 13, -13, 3, 0, -4, -5, -6, -7, -3, -2, -1]
        grpns=['A1', 'A1A2', 'A2', 'A2A4', 'A3', 'A3A1', 'A3C5', 'A4', 'A4A3', 'AP1', 'AP2', 'AP3', 'AP4', 'APPS', 'B1', 'B1C1', 'B2', 'B2B4', 'B3', 'B3B1', 'B4', 'B4C3', 'C1', 'C1C4', 'C2', 'C2B2', 'C3', 'C3C2', 'C4', 'C4B3', 'C5', 'C5A4', 'DALLE', 'DALQ1', 'DALQ2', 'DALT3', 'MESH', 'PMMA']
        famIdsPerGrp=[[5],[-11],[1],[-8],[4],[-10],[-12],[2],[-9],[-4],[-5],[-6],[-7],[-4,-5,-6,-7],[6],[-14],[7],[-19],[8],[-17],[9],[-20],[10],[-15],[11],[-18],[12],[-21],[13],[-16],[3],[-13],[-3,-2,-1],[-3],[-2],[-1],[-4,-5,-6,-7,-3,-2,-1],[3]]
        return m2,m1,m0,f2,f1,f0,p,n2,n1,n0,fns,fids,grpns,famIdsPerGrp

    def buildMLMeshUnPolyze(cls,tester):
        """Level 0 (meshDim=3) - 2 TETRA4 + 3 PENTA6 + 2 POLYH
        # POLYH #0 becomes 1 TETRA4
        # POLYH #1 becomes HEXA8
        # Level -1 (meshDim=2) - 2 TRI3 + 3 QUAD4 + 4 POLYG
        # POLYG #2 becomes TRI3"""
        meshName="NightmareMesh"
        #
        coords=DataArrayDouble.New(38,3) ; coords.rearrange(1) ; coords.iota(1000.) ; coords.rearrange(3) ; coords.setInfoOnComponents(["X [m]","Y [m]","Z [m]"])
        mesh0=MEDCouplingUMesh(meshName,3)
        type0=[NORM_TETRA4,NORM_TETRA4, NORM_PENTA6,NORM_PENTA6,NORM_PENTA6, NORM_POLYHED,NORM_POLYHED]
        conn0=[[0,1,2,3],[4,5,6,7], [8,9,10,11,12,13],[14,15,16,17,18,19],[20,21,22,23,24,25], [26,27,28,-1,26,29,27,-1,27,29,28,-1,28,29,26],[30,31,32,33,-1,34,37,36,35,-1,30,34,35,31,-1,31,35,36,32,-1,32,36,37,33,-1,33,37,34,30]]
        mesh0.allocateCells(len(type0))
        for typ,nodalConn in zip(type0,conn0):
            mesh0.insertNextCell(typ,nodalConn);
            pass
        mesh0.finishInsertingCells()
        mesh0.setCoords(coords)
        
        meshM1=MEDCouplingUMesh(meshName,2)
        typeM1=[NORM_TRI3,NORM_TRI3, NORM_QUAD4,NORM_QUAD4,NORM_QUAD4, NORM_POLYGON,NORM_POLYGON,NORM_POLYGON,NORM_POLYGON]
        connM1=[[0,1,2],[3,4,5], [6,7,8,9],[10,11,12,13],[14,15,16,17], [18,19,20,21,22],[23,24,25,26,27],[28,29,30],[31,32,33,34,35,36,37]]
        meshM1.allocateCells(len(typeM1))
        for typ,nodalConn in zip(typeM1,connM1):
            meshM1.insertNextCell(typ,nodalConn);
            pass
        meshM1.finishInsertingCells()
        meshM1.setCoords(coords)
        #
        mm=MEDFileUMesh.New()
        mm.setMeshAtLevel(0,mesh0)
        mm.setMeshAtLevel(-1,meshM1)
        grp0_L0=DataArrayInt.New([0,1,5,7]) ; grp0_L0.setName("grp0_L0")
        grp1_L0=DataArrayInt.New([1,2,3,4,6]) ; grp1_L0.setName("grp1_L0")
        tester.assertRaises(InterpKernelException,mm.setGroupsAtLevel,0,[grp0_L0,grp1_L0])# presence of 7 in grp0_L0 (only 7 cells at level 0) -> throw
        grp0_L0=DataArrayInt.New([0,1,5,6]) ; grp0_L0.setName("grp0_L0")
        mm.setGroupsAtLevel(0,[grp0_L0,grp1_L0])
        grp0_LM1=DataArrayInt.New([1,2,3,4,7]) ; grp0_LM1.setName("grp0_LM1")
        grp1_LM1=DataArrayInt.New([2,3,4,5]) ; grp1_LM1.setName("grp1_LM1")
        grp2_LM1=DataArrayInt.New([5,6,7,8]) ; grp2_LM1.setName("grp2_LM1")
        mm.setGroupsAtLevel(-1,[grp0_LM1,grp1_LM1,grp2_LM1])
        grp0_Node=DataArrayInt.New([0,11,15,16]) ; grp0_Node.setName("grp0_Node")
        grp1_Node=DataArrayInt.New([1,2,13,14,16]) ; grp1_Node.setName("grp1_Node")
        mm.setGroupsAtLevel(1,[grp0_Node,grp1_Node])
        #
        tester.assertRaises(InterpKernelException,mm.setRenumFieldArr,0,DataArrayInt.New([0,8,9,4,5,6,7,10]))# to big array
        mm.setRenumFieldArr(0,DataArrayInt.New([0,8,9,4,5,6,7]))
        da=DataArrayInt.New([0,8,9,4,5,6,7,11,12])
        mm.setRenumFieldArr(-1,da)
        mm.setRenumFieldArr(-1,None)
        mm.setRenumFieldArr(-1,da)
        da=DataArrayInt.New(mm.getNumberOfNodes()+1) ; da.iota(8) ; tester.assertRaises(InterpKernelException,mm.setRenumFieldArr,1,da) # to big array more than number of nodes
        da=DataArrayInt.New(mm.getNumberOfNodes()) ; da.iota(8) ; mm.setRenumFieldArr(1,da)
        return mm

    def buildVecFieldOnCells_1(cls):
        mesh=MEDLoaderDataForTest.build3DSurfMesh_1();
        nbOfCells=mesh.getNumberOfCells();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setName("VectorFieldOnCells");
        f1.setMesh(mesh);
        array=DataArrayDouble.New();
        arr1=[0.,10.,20.,1.,11.,21.,2.,12.,22.,3.,13.,23.,4.,14.,24.,5.,15.,25.]
        array.setValues(arr1,nbOfCells,3);
        array.setInfoOnComponent(0,"power [MW/m^3]");
        array.setInfoOnComponent(1,"density [g/cm^3]");
        array.setInfoOnComponent(2,"temperature [K]");
        f1.setArray(array);
        tmp=array.getPointer();
        f1.setTime(2.,0,1);
        f1.checkConsistencyLight();
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
        array.setInfoOnComponent(0,"power [MW/m^3]");
        array.setInfoOnComponent(1,"density [g/cm^3]");
        array.setInfoOnComponent(2,"temperature [K]");
        f1.setTime(2.12,2,3);
        f1.checkConsistencyLight();
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
        for i in range(19 * 2):
            array.setIJ(0,i,float(i+7));
            pass
        f.setArray(array);
        f.setName("MyFirstFieldOnGaussPoint");
        array.setInfoOnComponent(0,"power [MW/m^3]");
        array.setInfoOnComponent(1,"density");
        f.checkConsistencyLight();
        return f;

    def buildVecFieldOnGauss_2(cls):
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
        m=MEDLoaderDataForTest.build2DMesh_3();
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,ONE_TIME);
        f.setTime(3.14,1,5);
        f.setMesh(m);
        di=DataArrayInt.New(); di.setValues([0,2,3],3,1)
        f.setGaussLocalizationOnCells(di,_refCoo1,_gsCoo1,_wg1)
        _wg1[-1]*=2
        f.setGaussLocalizationOnCells([1,5],_refCoo1,_gsCoo1,_wg1);
        _wg1[-1]*=2
        f.setGaussLocalizationOnCells([4],_refCoo1,_gsCoo1,_wg1);
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
        f.setGaussLocalizationOnCells([6,7,8],_refCoo3,_gsCoo1,_wg1);
        _wg1[-1]*=2
        f.setGaussLocalizationOnCells([9],_refCoo3,_gsCoo1,_wg1);
        f.setGaussLocalizationOnType(NORM_TRI6,_refCoo2,_gsCoo2,_wg2);
        array=DataArrayDouble.New();
        array.alloc(53,2);
        ptr=array.getPointer();
        for i in range(53 * 2):
            array.setIJ(0,i,float(i+7));
            pass
        f.setArray(array);
        f.setName("MyFirstFieldOnGaussPoint");
        array.setInfoOnComponent(0,"power [MW/m^3]");
        array.setInfoOnComponent(1,"density");
        f.checkConsistencyLight();
        return f;

    # idem buildVecFieldOnGauss_2 except that different discretizations are sorted inside one type
    def buildVecFieldOnGauss_2_Simpler(cls):
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
        m=MEDLoaderDataForTest.build2DMesh_3();
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,ONE_TIME);
        f.setTime(3.14,1,5);
        f.setMesh(m);
        di=DataArrayInt.New(); di.setValues([0,1,2],3,1)
        f.setGaussLocalizationOnCells(di,_refCoo1,_gsCoo1,_wg1)
        _wg1[-1]*=2
        f.setGaussLocalizationOnCells([3,4],_refCoo1,_gsCoo1,_wg1);
        _wg1[-1]*=2
        f.setGaussLocalizationOnCells([5],_refCoo1,_gsCoo1,_wg1);
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
        f.setGaussLocalizationOnCells([6,7,8],_refCoo3,_gsCoo1,_wg1);
        _wg1[-1]*=2
        f.setGaussLocalizationOnCells([9],_refCoo3,_gsCoo1,_wg1);
        f.setGaussLocalizationOnType(NORM_TRI6,_refCoo2,_gsCoo2,_wg2);
        array=DataArrayDouble.New();
        array.alloc(53,2);
        ptr=array.getPointer();
        for i in range(53 * 2):
            array.setIJ(0,i,float(i+7));
            pass
        f.setArray(array);
        f.setName("MyFirstFieldOnGaussPoint");
        array.setInfoOnComponent(0,"power [MW/m^3]");
        array.setInfoOnComponent(1,"density");
        f.checkConsistencyLight();
        return f;

    def buildVecFieldOnGaussNE_1(cls):
        m=MEDLoaderDataForTest.build2DMesh_2();
        f=MEDCouplingFieldDouble.New(ON_GAUSS_NE,ONE_TIME);
        f.setTime(3.14,1,5);
        f.setMesh(m);
        array=DataArrayDouble.New();
        array.alloc(20,2);
        for i in range(2 * 20):
            array.setIJ(0,i,float(i+8));
        f.setArray(array);
        array.setInfoOnComponent(0,"power [W]");
        array.setInfoOnComponent(1,"temperature");
        f.setName("MyFieldOnGaussNE");
        f.checkConsistencyLight();
        return f;

    def buildACompleteMEDDataStructureWithFieldsOnCells_1(cls):
        coo=DataArrayDouble([0,0,1,0,2,0,0,1,1,1,2,1,0,2,1,2,2,2],9,2)
        m0=MEDCouplingUMesh("mesh",2)
        m0.setCoords(coo)
        m0.allocateCells()
        m0.insertNextCell(NORM_TRI3,[1,4,2])
        m0.insertNextCell(NORM_TRI3,[4,5,2])
        m0.insertNextCell(NORM_QUAD4,[0,3,4,1])
        m0.insertNextCell(NORM_QUAD4,[6,7,4,3])
        m0.insertNextCell(NORM_QUAD4,[7,8,5,4])
        m1=m0.computeSkin()
        mm=MEDFileUMesh()
        #2 levels
        mm.setMeshAtLevel(0,m0) ; mm.setMeshAtLevel(-1,m1)
        #some grps/families on the 2 levels
        grp0=DataArrayInt([0,2,4]); grp0.setName("gr0_0_2_4")
        grp1=DataArrayInt([1,2,3,4]); grp1.setName("gr0_1_2_3_4")
        grp2=DataArrayInt([0,4]); grp2.setName("gr0_0_4")
        mm.setGroupsAtLevel(0,[grp0,grp1,grp2])
        grp3=DataArrayInt([0,1]); grp3.setName("grM1_SegOnTri3")
        grp4=DataArrayInt([2,3,4,5,6,7]); grp4.setName("grM1_SegOnQuad4")
        grp5=DataArrayInt([0,3]); grp5.setName("grM1_bottom")
        mm.setGroupsAtLevel(-1,[grp3,grp4,grp5])
        ms=MEDFileMeshes()
        ms.pushMesh(mm)
        # 3 fields
        fs=MEDFileFields()
        # 1st Field - fNoProfile - no profile on levels 0
        f1Name="fNoProfile"
        timeStepsF1=[(0,-1,0.01),(1,-1,0.02)]
        f1=MEDFileFieldMultiTS()
        for i,(it,order,tim) in enumerate(timeStepsF1):
            f11Tmp=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME)
            f11Tmp.setTime(tim,it,order)
            f11Tmp.setMesh(m0)
            arr=DataArrayDouble(m0.getNumberOfCells(),1) ; arr.iota() ; arr+=1+i ; arr*=0.1
            f11Tmp.setArray(arr)
            f11Tmp.checkConsistencyLight()
            f11Tmp.setName(f1Name)
            f1.appendFieldNoProfileSBT(f11Tmp)
            pass
        fs.pushField(f1)
        # 2nd Field - fNoProfileMultiLevs - no profile on levels 0 and -1
        f2Name="fNoProfileMultiLevs"
        timeStepsF2=[(0,-1,0.),(1,-1,0.1),(2,-1,0.2)]
        f2=MEDFileFieldMultiTS()
        for i,(it,order,tim) in enumerate(timeStepsF2):
            f21Tmp=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME)
            f21Tmp.setTime(tim,it,order)
            f21Tmp.setMesh(m0)
            arr=DataArrayDouble(m0.getNumberOfCells(),1) ; arr.iota() ; arr+=1+i
            f21Tmp.setArray(arr)
            f21Tmp.checkConsistencyLight()
            f21Tmp.setName(f2Name)
            f2.appendFieldNoProfileSBT(f21Tmp)
            f22Tmp=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME)
            f22Tmp.setTime(tim,it,order)
            f22Tmp.setMesh(m1)
            arr=DataArrayDouble(m1.getNumberOfCells(),1) ; arr.iota() ; arr+=100+1+i
            f22Tmp.setArray(arr)
            f22Tmp.checkConsistencyLight()
            f22Tmp.setName(f2Name)
            f2[it,order].setFieldNoProfileSBT(f22Tmp)
            pass
        fs.pushField(f2)
        # 3rd field - fProfileMultiLevs - The most complex one
        f3Name="fProfileMultiLevs"
        timeStepsF3=[(0,-1,0.),(1,-1,10.),(2,-1,20.),(3,-1,30.),]
        f3=MEDFileFieldMultiTS()
        for i,(it,order,tim) in enumerate(timeStepsF3):
            pfl1=DataArrayInt([0,1,3,4]) ; pfl1.setName("pfl1")
            m0Part=m0[pfl1]
            f31Tmp=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME)
            f31Tmp.setTime(tim,it,order)
            f31Tmp.setMesh(m0Part)
            arr=DataArrayDouble(m0Part.getNumberOfCells(),1) ; arr.iota() ; arr+=1000+i+1
            f31Tmp.setArray(arr)
            f31Tmp.checkConsistencyLight()
            f31Tmp.setName(f3Name)
            f3.appendFieldProfile(f31Tmp,mm,0,pfl1)
            pfl2=DataArrayInt([0,3]) ; pfl2.setName("pfl2Bottom")
            m1Part=m1[pfl2]
            f32Tmp=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME)
            f32Tmp.setTime(tim,it,order)
            f32Tmp.setMesh(m1Part)
            arr=DataArrayDouble(m1Part.getNumberOfCells(),1) ; arr.iota() ; arr+=2000+1+i
            f32Tmp.setArray(arr)
            f32Tmp.checkConsistencyLight()
            f32Tmp.setName(f3Name)
            f3[it,order].setFieldProfile(f32Tmp,mm,-1,pfl2)
            pass
        fs.pushField(f3)
        #
        data=MEDFileData() ; data.setMeshes(ms) ; data.setFields(fs)
        return data

    def buildAMEDFileDataWithGroupOnOneFamilyForSauv(self):
        # Coordinates
        coords = [0.,0., 0.,1., 1.,1., 1.,0.]
        # lvl 0 connectivity
        conn2D   = [1,2,3,4]
        # lvl -1 connectivity
        conn1D   = [0,1, 1,2, 2,3, 4,1]
        # lvl 0 mesh
        mesh2D=MEDCouplingUMesh.New()
        mesh2D.setMeshDimension(2)
        mesh2D.allocateCells(1)
        mesh2D.insertNextCell(NORM_QUAD4,4,conn2D)
        mesh2D.finishInsertingCells()
        # lvl -1 mesh
        mesh1D=MEDCouplingUMesh.New()
        mesh1D.setMeshDimension(1)
        mesh1D.allocateCells(4)
        mesh1D.insertNextCell(NORM_SEG2,2,conn1D[0:2])
        mesh1D.insertNextCell(NORM_SEG2,2,conn1D[2:4])
        mesh1D.insertNextCell(NORM_SEG2,2,conn1D[4:6])
        mesh1D.insertNextCell(NORM_SEG2,2,conn1D[6:8])
        mesh1D.finishInsertingCells()
        # assigning coordinates
        meshCoords=DataArrayDouble.New()
        meshCoords.setValues(coords, 4, 2)
        mesh2D.setCoords(meshCoords)
        mesh1D.setCoords(meshCoords)
        # Creating a multi level mesh
        mm = MEDFileUMesh.New()
        mm.setMeshAtLevel(0, mesh2D)
        mm.setMeshAtLevel(-1, mesh1D)
        mm.setName("carre")
        # Creating groups
        # Creating a group with an element on level -1
        grp0_LM1 = DataArrayInt.New([0])
        grp0_LM1.setName("grp0_LM1")
        # Creating a group with all elements on level -1
        grp1_LM1 = DataArrayInt.New([0,1,2,3])
        grp1_LM1.setName("grp1_LM1")
        #
        mm.setGroupsAtLevel(-1,[grp0_LM1,grp1_LM1])
        #
        ms=MEDFileMeshes.New()
        ms.setMeshAtPos(0,mm)
        mfd=MEDFileData.New()
        mfd.setMeshes(ms)
        #
        return mfd
    
    build1DMesh_1=classmethod(build1DMesh_1)
    build2DCurveMesh_1=classmethod(build2DCurveMesh_1)
    build2DMesh_1=classmethod(build2DMesh_1)
    build2DMesh_2=classmethod(build2DMesh_2)
    build2DMesh_3=classmethod(build2DMesh_3)
    build3DMesh_1=classmethod(build3DMesh_1)
    build3DSurfMesh_1=classmethod(build3DSurfMesh_1)
    build3DMesh_2=classmethod(build3DMesh_2)
    buildMLMeshUnPolyze=classmethod(buildMLMeshUnPolyze)
    buildMultiLevelMesh_1=classmethod(buildMultiLevelMesh_1)
    buildVecFieldOnCells_1=classmethod(buildVecFieldOnCells_1)
    buildVecFieldOnNodes_1=classmethod(buildVecFieldOnNodes_1)
    buildVecFieldOnGauss_1=classmethod(buildVecFieldOnGauss_1)
    buildVecFieldOnGauss_2=classmethod(buildVecFieldOnGauss_2)
    buildVecFieldOnGauss_2_Simpler=classmethod(buildVecFieldOnGauss_2_Simpler)
    buildVecFieldOnGaussNE_1=classmethod(buildVecFieldOnGaussNE_1)
    buildACompleteMEDDataStructureWithFieldsOnCells_1=classmethod(buildACompleteMEDDataStructureWithFieldsOnCells_1)
    buildAMEDFileDataWithGroupOnOneFamilyForSauv=classmethod(buildAMEDFileDataWithGroupOnOneFamilyForSauv)
    pass
