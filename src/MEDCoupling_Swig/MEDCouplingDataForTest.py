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

class MEDCouplingDataForTest:
    def build2DTargetMesh_1(cls):
        targetCoords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ];
        targetConn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4];
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
    
    def build2DSourceMesh_1(cls):
        sourceCoords=[-0.3,-0.3, 0.7,-0.3, -0.3,0.7, 0.7,0.7]
        sourceConn=[0,3,1,0,2,3]
        sourceMesh=MEDCouplingUMesh.New("my name of mesh 2D",2);
        sourceMesh.allocateCells(2);
        sourceMesh.insertNextCell(NORM_TRI3,3,sourceConn[0:3]);
        sourceMesh.insertNextCell(NORM_TRI3,3,sourceConn[3:6]);
        sourceMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(sourceCoords,4,2);
        sourceMesh.setCoords(myCoords);
        return sourceMesh;
        
    def build3DTargetMesh_1(cls):
        targetCoords=[ 0., 0., 0., 50., 0., 0. , 200., 0., 0.  , 0., 50., 0., 50., 50., 0. , 200., 50., 0.,   0., 200., 0., 50., 200., 0. , 200., 200., 0. ,
                       0., 0., 50., 50., 0., 50. , 200., 0., 50.  , 0., 50., 50., 50., 50., 50. , 200., 50., 50.,   0., 200., 50., 50., 200., 50. , 200., 200., 50. ,
                       0., 0., 200., 50., 0., 200. , 200., 0., 200.  , 0., 50., 200., 50., 50., 200. , 200., 50., 200.,   0., 200., 200., 50., 200., 200. , 200., 200., 200. ];
        targetConn=[0,1,4,3,9,10,13,12, 1,2,5,4,10,11,14,13, 3,4,7,6,12,13,16,15, 4,5,8,7,13,14,17,16,
                    9,10,13,12,18,19,22,21, 10,11,14,13,19,20,23,22, 12,13,16,15,21,22,25,24, 13,14,17,16,22,23,26,25];
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.setMeshDimension(3);
        targetMesh.allocateCells(12);
        for i in range(8):
            targetMesh.insertNextCell(NORM_HEXA8,8,targetConn[8*i:8*i+8]);
            pass
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,27,3);
        targetMesh.setCoords(myCoords);
        return targetMesh

    def build3DSourceMesh_1(self):
        sourceCoords=[ 0.0, 0.0, 200.0, 0.0, 0.0, 0.0, 0.0, 200.0, 200.0, 0.0, 200.0, 0.0, 200.0, 0.0, 200.0,
                       200.0, 0.0, 0.0, 200.0, 200.0, 200.0, 200.0, 200.0, 0.0, 100.0, 100.0, 100.0]
        sourceConn=[8,1,7,3, 6,0,8,2, 7,4,5,8, 6,8,4,7, 6,8,0,4, 6,8,7,3, 8,1,3,0, 4,1,5,8, 1,7,5,8, 0,3,8,2, 8,1,0,4, 3,6,8,2]
        sourceMesh=MEDCouplingUMesh.New();
        sourceMesh.setMeshDimension(3);
        sourceMesh.allocateCells(12);
        sourceMesh.insertNextCell(NORM_TETRA4,4,sourceConn[0:4])
        sourceMesh.insertNextCell(NORM_TETRA4,4,sourceConn[4:8])
        sourceMesh.insertNextCell(NORM_TETRA4,4,sourceConn[8:12])
        sourceMesh.insertNextCell(NORM_TETRA4,4,sourceConn[12:16])
        sourceMesh.insertNextCell(NORM_TETRA4,4,sourceConn[16:20])
        sourceMesh.insertNextCell(NORM_TETRA4,4,sourceConn[20:24])
        sourceMesh.insertNextCell(NORM_TETRA4,4,sourceConn[24:28])
        sourceMesh.insertNextCell(NORM_TETRA4,4,sourceConn[28:32])
        sourceMesh.insertNextCell(NORM_TETRA4,4,sourceConn[32:36])
        sourceMesh.insertNextCell(NORM_TETRA4,4,sourceConn[36:40])
        sourceMesh.insertNextCell(NORM_TETRA4,4,sourceConn[40:44])
        sourceMesh.insertNextCell(NORM_TETRA4,4,sourceConn[44:48])
        sourceMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(sourceCoords,9,3);
        sourceMesh.setCoords(myCoords);
        return sourceMesh;
        

    def build3DSurfTargetMesh_1(self):
        targetCoords=[-0.3,-0.3,0.5, 0.2,-0.3,1., 0.7,-0.3,1.5, -0.3,0.2,0.5, 0.2,0.2,1., 0.7,0.2,1.5, -0.3,0.7,0.5, 0.2,0.7,1., 0.7,0.7,1.5]
        targetConn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.setMeshDimension(2);
        targetMesh.allocateCells(5);
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[0:4])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[4:7])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[7:10])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[10:14])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[14:18])
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,9,3);
        targetMesh.setCoords(myCoords);
        return targetMesh;

    def build3DExtrudedUMesh_1(self):
        coords=[
            0.,0.,0., 1.,1.,0., 1.,1.25,0., 1.,0.,0., 1.,1.5,0., 2.,0.,0., 2.,1.,0., 1.,2.,0., 0.,2.,0., 3.,1.,0.,
            3.,2.,0., 0.,1.,0., 1.,3.,0., 2.,2.,0., 2.,3.,0.,
            0.,0.,1., 1.,1.,1., 1.,1.25,1., 1.,0.,1., 1.,1.5,1., 2.,0.,1., 2.,1.,1., 1.,2.,1., 0.,2.,1., 3.,1.,1.,
            3.,2.,1., 0.,1.,1., 1.,3.,1., 2.,2.,1., 2.,3.,1.,
            0.,0.,2., 1.,1.,2., 1.,1.25,2., 1.,0.,2., 1.,1.5,2., 2.,0.,2., 2.,1.,2., 1.,2.,2., 0.,2.,2., 3.,1.,2.,
            3.,2.,2., 0.,1.,2., 1.,3.,2., 2.,2.,2., 2.,3.,2.,
            0.,0.,3., 1.,1.,3., 1.,1.25,3., 1.,0.,3., 1.,1.5,3., 2.,0.,3., 2.,1.,3., 1.,2.,3., 0.,2.,3., 3.,1.,3.,
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
        conn2=[7,12,14,13, 11,8,7,4,2,1, 13,10,9,6, 1,6,5,3, 1,2,4,7,13,6, 0,11,1,3]
        #
        ret=MEDCouplingUMesh.New();
        ret.setMeshDimension(3);
        ret.allocateCells(18);
        #
        ret.insertNextCell(NORM_HEXA8,8,conn[0:8]);
        ret.insertNextCell(NORM_POLYHED,43,conn[8:51]);
        ret.insertNextCell(NORM_HEXA8,8,conn[51:59]);
        ret.insertNextCell(NORM_HEXA8,8,conn[59:67]);
        ret.insertNextCell(NORM_POLYHED,43,conn[67:110]);
        ret.insertNextCell(NORM_HEXA8,8,conn[110:118]);
        #
        ret.insertNextCell(NORM_HEXA8,8,conn[118:126]);
        ret.insertNextCell(NORM_POLYHED,43,conn[126:169]);
        ret.insertNextCell(NORM_HEXA8,8,conn[169:177]);
        ret.insertNextCell(NORM_HEXA8,8,conn[177:185]);
        ret.insertNextCell(NORM_POLYHED,43,conn[185:228]);
        ret.insertNextCell(NORM_HEXA8,8,conn[228:236]);
        #
        ret.insertNextCell(NORM_HEXA8,8,conn[236:244]);
        ret.insertNextCell(NORM_POLYHED,43,conn[244:287]);
        ret.insertNextCell(NORM_HEXA8,8,conn[287:295]);
        ret.insertNextCell(NORM_HEXA8,8,conn[295:303]);
        ret.insertNextCell(NORM_POLYHED,43,conn[303:346]);
        ret.insertNextCell(NORM_HEXA8,8,conn[346:354]);
        #
        ret.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(coords,60,3);
        ret.setCoords(myCoords);
        #
        mesh2D=MEDCouplingUMesh.New();
        mesh2D.setMeshDimension(2);
        mesh2D.allocateCells(6);
        mesh2D.insertNextCell(NORM_QUAD4,4,conn2[0:4]);
        mesh2D.insertNextCell(NORM_POLYGON,6,conn2[4:10]);
        mesh2D.insertNextCell(NORM_QUAD4,4,conn2[10:14]);
        mesh2D.insertNextCell(NORM_QUAD4,4,conn2[14:18]);
        mesh2D.insertNextCell(NORM_POLYGON,6,conn2[18:24]);
        mesh2D.insertNextCell(NORM_QUAD4,4,conn2[24:28]);
        mesh2D.finishInsertingCells();
        mesh2D.setCoords(myCoords);
        return ret,mesh2D
    
    def buildCU1DMesh_U(self):
        coords=[ 0.0, 0.3, 0.75, 1.0 ]
        conn=[ 0,1, 1,2, 2,3 ]
        mesh=MEDCouplingUMesh.New();
        mesh.setMeshDimension(1);
        mesh.allocateCells(3);
        mesh.insertNextCell(NORM_SEG2,2,conn[0:2]);
        mesh.insertNextCell(NORM_SEG2,2,conn[2:4]);
        mesh.insertNextCell(NORM_SEG2,2,conn[4:6]);
        mesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(coords,4,1);
        mesh.setCoords(myCoords);
        return mesh;

    def build2DTargetMeshMergeNode_1(self):
        targetCoords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,-0.3, 0.2,-0.3, 0.2,-0.3, 0.2,0.2, 0.2,0.2, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7, 0.2,0.7 ]
        targetConn=[0,9,7,5, 4,6,2, 10,11,8, 9,14,15,7, 17,16,13,6]
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
        myCoords.setValues(targetCoords,18,2);
        targetMesh.setCoords(myCoords);
        return targetMesh;

    def build3DTargetMeshMergeNode_1(self):
        targetCoords=[ 0., 0., 0., 50., 0., 0. , 200., 0., 0.  , 0., 50., 0., 50., 50., 0. , 200., 50., 0.,   0., 200., 0., 50., 200., 0. , 200., 200., 0. ,
                       0., 0., 50., 50., 0., 50. , 200., 0., 50.  , 0., 50., 50., 50., 50., 50. , 200., 50., 50.,   0., 200., 50., 50., 200., 50. , 200., 200., 50. ,
                       0., 0., 200., 50., 0., 200. , 200., 0., 200.  , 0., 50., 200., 50., 50., 200. , 200., 50., 200.,   0., 200., 200., 50., 200., 200. , 200., 200., 200., 50.,0.,0., 50.,0.,0., 50.,0.,0.,  200., 50., 200.]
        targetConn=[0,29,4,3,9,10,13,12, 28,2,5,4,10,11,14,13, 3,4,7,6,12,13,16,15, 4,5,8,7,13,14,17,16,
                    9,10,13,12,18,19,22,21, 10,11,14,13,19,20,23,22, 12,13,16,15,21,22,25,24, 13,14,17,16,22,30,26,25]
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.setMeshDimension(3);
        targetMesh.allocateCells(12);
        for i in range(8):
            targetMesh.insertNextCell(NORM_HEXA8,8,targetConn[8*i:8*(i+1)]);
            pass
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,31,3);
        targetMesh.setCoords(myCoords);
        return targetMesh;

    def build2DTargetMeshMerged_1(self):
        targetCoords=[
            -0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7,
            0.7,-0.3, 1.7,-0.3, 0.7,0.7, 1.7,0.7
            ]
        targetConn=[
            0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4,
            9,12,10,9,11,12
            ]
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.setName("merge");
        targetMesh.setMeshDimension(2);
        targetMesh.allocateCells(10);
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[0:4])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[4:7])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[7:10])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[10:14])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[14:18])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[18:21])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[21:24])
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,13,2);
        targetMesh.setCoords(myCoords);
        return targetMesh;

    def build2DTargetMesh_2(cls):
        targetCoords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        targetConn=[0,3,4, 0,4,1, 1,4,2, 4,5,2, 3,6,4, 6,7,4, 4,7,5, 7,8,5 ]
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.setMeshDimension(2);
        targetMesh.allocateCells(8);
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[0:3])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[3:6])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[6:9])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[9:12])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[12:15])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[15:18])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[18:21])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[21:24])
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,9,2);
        targetMesh.setCoords(myCoords);
        return targetMesh;
    
    def build1DSourceMesh_2(cls):
        ret=MEDCouplingUMesh.New("1DSourceMesh",1);
        ret.allocateCells(4);
        conn=[0,1,2,3,1,2,3,4]
        for i in range(4):
            ret.insertNextCell(NORM_SEG2,2,conn[2*i:2*i+2]);
            pass
        ret.finishInsertingCells();
        myCoords=DataArrayDouble.New([0.3,0.7,0.9,1.0,1.12],5,1);
        ret.setCoords(myCoords);
        return ret

    def build1DTargetMesh_3(cls):
        ret=MEDCouplingUMesh.New("1DMesh_3",1);
        ret.allocateCells(4);
        conn=[0,1,2, 3,4, 6,5,7 ,9,8]
        ret.insertNextCell(NORM_SEG3,3,conn[0:3])
        ret.insertNextCell(NORM_SEG2,2,conn[3:5])
        ret.insertNextCell(NORM_SEG3,3,conn[5:8])
        ret.insertNextCell(NORM_SEG2,2,conn[8:10])
        ret.finishInsertingCells();
        coords=[0.5,1.,0.8,5.,5.21,0.5,1.1,0.7,5.,5.31]
        myCoords=DataArrayDouble.New();
        myCoords.setValues(coords,10,1);
        ret.setCoords(myCoords);
        return ret;

    def build2DCurveTargetMesh_3(cls):
        ret=MEDCouplingUMesh.New("2DCurveMesh_3",1);
        ret.allocateCells(4);
        conn=[0,1,2, 3,4, 6,5,7 ,9,8]
        ret.insertNextCell(NORM_SEG3,3,conn[0:3])
        ret.insertNextCell(NORM_SEG2,2,conn[3:5])
        ret.insertNextCell(NORM_SEG3,3,conn[5:8])
        ret.insertNextCell(NORM_SEG2,2,conn[8:10])
        ret.finishInsertingCells();
        coords=[0.5,0.5,1.,1.,0.8,0.8,5.,5.,5.21,5.21,0.5,0.5,1.1,1.1,0.7,0.7,5.,5.,5.31,5.31]
        myCoords=DataArrayDouble.New();
        myCoords.setValues(coords,10,2);
        ret.setCoords(myCoords);
        return ret;

    def build2DTargetMesh_3(cls):
        ret=MEDCouplingUMesh.New("2DMesh_3",2);
        ret.allocateCells(10);
        conn=[0,1,2, 0,1,3,4, 0,1,3,5,4, 0,1,2,6,7,8, 0,1,3,4,6,9,2,10, 0,2,1, 0,4,3,1, 0,4,5,3,1, 0,2,1,8,7,6, 0,4,3,1,10,2,9,6]
        ret.insertNextCell(NORM_TRI3,3,conn[0:3])
        ret.insertNextCell(NORM_QUAD4,4,conn[3:7])
        ret.insertNextCell(NORM_POLYGON,5,conn[7:12])
        ret.insertNextCell(NORM_TRI6,6,conn[12:18])
        ret.insertNextCell(NORM_QUAD8,8,conn[18:26])
        ret.insertNextCell(NORM_TRI3,3,conn[26:29])
        ret.insertNextCell(NORM_QUAD4,4,conn[29:33])
        ret.insertNextCell(NORM_POLYGON,5,conn[33:38])
        ret.insertNextCell(NORM_TRI6,6,conn[38:44])
        ret.insertNextCell(NORM_QUAD8,8,conn[44:52])
        ret.finishInsertingCells();
        coords=[0.,0.,1.,0.,0.5,1.,1.,1.,0.,1.,0.5,2.,0.5,0.,0.75,0.5,0.25,0.5,1.,0.5,0.,0.5]
        myCoords=DataArrayDouble.New();
        myCoords.setValues(coords,11,2);
        ret.setCoords(myCoords);
        ret.checkConsistencyLight();
        return ret;

    def build2DTargetMesh_4(cls):
        targetCoords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        targetConn=[0,4,5,1, 1,5,3, 5,6,2, 7,8,5,4, 8,9,6,5]
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.setMeshDimension(2);
        targetMesh.allocateCells(5);
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[0:4])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[4:7])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[7:10])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[10:14])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[14:18])
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,10,2);
        targetMesh.setCoords(myCoords);
        return targetMesh;

    def buildMultiFields_1(cls):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m1.setName("m1");
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2.setName("m2");
        vals0=[-0.7,-1.,-2.,-3.,-4.];
        vals1=[0.,1.,2.,3.,4.,0.1,0.2,0.3,0.4];
        vals1_1=[170.,171.,172.,173.,174.,170.1,170.2,170.3,170.4];
        vals2=[5.,6.,7.,8.,9.];
        vals4=[15.,16.,17.,18.,19.];
        d0=DataArrayDouble.New();
        d0.setValues(vals0,5,1);
        d1=DataArrayDouble.New();
        d1.setValues(vals1,9,1);
        d1_1=DataArrayDouble.New();
        d1_1.setValues(vals1_1,9,1);
        d2=DataArrayDouble.New();
        d2.setValues(vals2,5,1);
        d4=DataArrayDouble.New();
        d4.setValues(vals4,5,1);
        d0.setName("d0"); d1.setName("d1"); d1_1.setName("d1_1"); d2.setName("d2"); d4.setName("d4");
        d0.setInfoOnComponent(0,"c1");
        d1.setInfoOnComponent(0,"c6");
        d1_1.setInfoOnComponent(0,"c9");
        d2.setInfoOnComponent(0,"c5");
        d4.setInfoOnComponent(0,"c7");
        f0=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f0.setMesh(m1);
        f0.setArray(d0);
        f0.setTime(0.2,5,6);
        f0.setName("f0");
        f1=MEDCouplingFieldDouble.New(ON_NODES,LINEAR_TIME);
        f1.setMesh(m1);
        f1.setArrays([d1,d1_1]);
        f1.setStartTime(0.7,7,8);
        f1.setEndTime(1.2,9,10);
        f1.setName("f1");
        f2=MEDCouplingFieldDouble.New(ON_CELLS,CONST_ON_TIME_INTERVAL);
        f2.setMesh(m2);
        f2.setArray(d2);
        f2.setTime(1.2,11,12);
        f2.setEndTime(1.5,13,14);
        f2.setName("f2");
        f3=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f3.setMesh(m1);
        f3.setArray(d2);
        f3.setTime(1.7,15,16);
        f3.setName("f3");
        f4=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f4.setMesh(m2);
        f4.setArray(d4);
        f4.setName("f4");
        ret=MEDCouplingMultiFields.New([f0,f1,f2,f3,f4]);
        return ret;

    def buildMultiFields_2(cls):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m1.setName("m1");
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2.setName("m2");
        vals0=[-0.7,-1.,-2.,-3.,-4.];
        vals1=[0.,1.,2.,3.,4.];
        vals1_1=[170.,171.,172.,173.,174.];
        vals2=[5.,6.,7.,8.,9.];
        vals4=[15.,16.,17.,18.,19.];
        d0=DataArrayDouble.New();
        d0.setValues(vals0,5,1);
        d1=DataArrayDouble.New();
        d1.setValues(vals1,5,1);
        d1_1=DataArrayDouble.New();
        d1_1.setValues(vals1_1,5,1);
        d2=DataArrayDouble.New();
        d2.setValues(vals2,5,1);
        d4=DataArrayDouble.New();
        d4.setValues(vals4,5,1);
        d0.setName("d0"); d1.setName("d1"); d1_1.setName("d1_1"); d2.setName("d2"); d4.setName("d4");
        d0.setInfoOnComponent(0,"c1");
        d1.setInfoOnComponent(0,"c6");
        d1_1.setInfoOnComponent(0,"c9");
        d2.setInfoOnComponent(0,"c5");
        d4.setInfoOnComponent(0,"c7");
        f0=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f0.setMesh(m1);
        f0.setArray(d0);
        f0.setTime(0.2,5,6);
        f0.setName("f0");
        f1=MEDCouplingFieldDouble.New(ON_CELLS,LINEAR_TIME);
        f1.setMesh(m1);
        f1.setArrays([d1,d1_1]);
        f1.setStartTime(0.7,7,8);
        f1.setEndTime(1.2,9,10);
        f1.setName("f1");
        f2=MEDCouplingFieldDouble.New(ON_CELLS,CONST_ON_TIME_INTERVAL);
        f2.setMesh(m2);
        f2.setArray(d2);
        f2.setTime(1.2,11,12);
        f2.setEndTime(1.5,13,14);
        f2.setName("f2");
        f3=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f3.setMesh(m1);
        f3.setArray(d2);
        f3.setTime(1.7,15,16);
        f3.setName("f3");
        f4=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f4.setMesh(m2);
        f4.setArray(d4);
        f4.setName("f4");
        return [f0,f1,f2,f3,f4]

    def build1DMultiTypes_1(self):
        mesh=MEDCouplingUMesh.New("Multi1DMesh",1);
        coo=MEDCouplingDataForTest.buildCoordsForMultiTypes_1();
        conn=[0,2, 0,2,1]
        mesh.allocateCells(2);
        mesh.insertNextCell(NORM_SEG2,2,conn[0:2])
        mesh.insertNextCell(NORM_SEG3,3,conn[2:5])
        mesh.finishInsertingCells();
        mesh.setCoords(coo);
        return mesh;

    def build2DMultiTypes_1(self):
        mesh=MEDCouplingUMesh.New("Multi2DMesh",2);
        coo=MEDCouplingDataForTest.buildCoordsForMultiTypes_1();
        conn=[3,4,5, 3,4,5,6,7,8, 0,9,10,11, 0,9,10,11,12,13,14,15]
        mesh.allocateCells(4);
        mesh.insertNextCell(NORM_TRI3,3,conn[0:3])
        mesh.insertNextCell(NORM_TRI6,6,conn[3:9])
        mesh.insertNextCell(NORM_QUAD4,4,conn[9:13])
        mesh.insertNextCell(NORM_QUAD8,8,conn[13:21])
        mesh.finishInsertingCells();
        mesh.setCoords(coo);
        return mesh;

    def build3DMultiTypes_1(self):
        mesh=MEDCouplingUMesh.New("Multi3DMesh",3);
        coo=MEDCouplingDataForTest.buildCoordsForMultiTypes_1();
        conn=[0,16,17,18,
              0,16,17,18,19,20,21,22,23,24,
              0,11,10,9,25,
              0,11,10,9,25,15,14,13,12,26,27,28,29,
              0,30,31,32,33,34,
              0,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
              0,9,10,11,44,45,46,47,
              0,9,10,11,44,45,46,47,12,13,14,15,48,49,50,51,52,53,54,55 ];
        mesh.allocateCells(8);
        mesh.insertNextCell(NORM_TETRA4,4,conn[0:4])
        mesh.insertNextCell(NORM_TETRA10,10,conn[4:14])
        mesh.insertNextCell(NORM_PYRA5,5,conn[14:19])
        mesh.insertNextCell(NORM_PYRA13,13,conn[19:32])
        mesh.insertNextCell(NORM_PENTA6,6,conn[32:38])
        mesh.insertNextCell(NORM_PENTA15,15,conn[38:53])
        mesh.insertNextCell(NORM_HEXA8,8,conn[53:61])
        mesh.insertNextCell(NORM_HEXA20,20,conn[61:81])
        mesh.finishInsertingCells();
        mesh.setCoords(coo);
        return mesh;

    def buildCoordsForMultiTypes_1(self):
        coords=DataArrayDouble.New();
        data=[0.0,0.0,0.0, 0.5,0.5,0.5, 1.0,1.0,1.0, 1.0,1.0,0.0, 2.0,2.5,0.0, 6.0,1.5,0.0, 1.0,2.0,0.0, 4.5,2.5,0.0, 4.0,0.5,0.0, 0.0,4.0,0.0, 4.0,4.0,0.0, 4.0,0.0,0.0, 0.0,2.0,0.0, 2.0,4.0,0.0, 4.0,2.0,0.0, 2.0,0.0,0.0, 0.0,6.0,0.0, 3.0,3.0,0.0, 1.3,3.0,3.0, 0.0,3.0,0.0, 1.5,4.5,0.0, 1.5,1.5,0.0, 0.65,1.5,1.5, 0.65,4.5,1.5, 2.15,3.0,1.5, 2.0,2.0,2.0, 3.0,1.0,1.0, 3.0,3.0,1.0, 1.0,3.0,1.0, 1.0,1.0,1.0, 0.0,3.0,0.0, 2.0,0.0,0.0, 0.0,0.0,6.0, 0.0,3.0,6.0, 3.0,0.0,6.0, 0.0,1.5,0.0, 1.5,1.5,0.0, 1.5,0.0,0.0, 0.0,1.5,6.0, 1.5,1.5,6.0, 1.5,0.0,6.0, 0.0,0.0,3.0, 0.0,3.0,3.0, 3.0,0.0,3.0, 0.0,0.0,4.0, 0.0,4.0,4.0, 4.0,4.0,4.0, 4.0,0.0,4.0, 0.0,2.0,4.0, 2.0,4.0,4.0, 4.0,2.0,4.0, 2.0,0.0,4.0, 0.0,0.0,2.0, 0.0,4.0,2.0, 4.0,4.0,2.0, 4.0,0.0,2.0]
        coords.setValues(data,56,3);
        coords.setInfoOnComponent(0,"X (cm)");
        coords.setInfoOnComponent(1,"Y (cm)");
        coords.setInfoOnComponent(2,"Z (cm)");
        return coords

    def buildHexa8Mesh_1(self):
        mesh=MEDCouplingUMesh.New("Hexa8Only",3);
        coo=DataArrayDouble.New();
        coords=[0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 1.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.5, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.5, 1.0, 0.0, 0.5, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.0, 1.0, 0.5, 0.5, 1.0, 0.5, 1.0, 1.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5, 1.0, 1.0, 0.5, 1.0, 0.0, 1.0, 1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0]
        coo.setValues(coords,27,3);
        conn=[3,12,13,4,0,9,10,1,
              4,13,14,5,1,10,11,2,
              6,15,16,7,3,12,13,4,
              7,16,17,8,4,13,14,5,
              12,21,22,13,9,18,19,10,
              13,22,23,14,10,19,20,11,
              15,24,25,16,12,21,22,13,
              16,25,26,17,13,22,23,14];
        mesh.allocateCells(8);
        for i in range(8):
            mesh.insertNextCell(NORM_HEXA8,8,conn[8*i:8*(i+1)])
            pass
        mesh.finishInsertingCells();
        mesh.setCoords(coo);
        return mesh;

    def buildPointe_1(self):
        mesh=MEDCouplingUMesh.New("Pointe.med",3);
        mesh2=MEDCouplingUMesh.New("Pointe.med",2);
        coords=[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 1.0, 0.0, 2.0, 1.0, -2.0, 0.0, 1.0, 0.0, -2.0, 1.0, 1.0, 1.0, 2.0, -1.0, 1.0, 2.0, -1.0, -1.0, 2.0, 1.0, -1.0, 2.0, 1.0, 1.0, 3.0, -1.0, 1.0, 3.0, -1.0, -1.0, 3.0, 1.0, -1.0, 3.0, 1.0, 1.0, 4.0, -1.0, 1.0, 4.0, -1.0, -1.0, 4.0, 1.0, -1.0, 4.0, 0.0, 0.0, 5.0]
        conn=[0,1,2,5,0,1,3,2,0,1,4,3,0,1,5,4,1,6,3,2,1,7,4,3,1,8,5,4,1,9,2,5,1,6,2,9,1,7,3,6,1,8,4,7,1,9,5,8, 6,7,8,9,1,14,17,16,15,18, 10,11,12,13,6,7,8,9,14,15,16,17,10,11,12,13]
        coo=DataArrayDouble.New();
        coo.setValues(coords,19,3);
        mesh.setCoords(coo);
        mesh2.setCoords(coo);
        mesh.allocateCells(16);
        for i in range(12):
            mesh.insertNextCell(NORM_TETRA4,4,conn[4*i:4*i+4])
            pass
        mesh.insertNextCell(NORM_PYRA5,5,conn[48:53])
        mesh.insertNextCell(NORM_PYRA5,5,conn[53:58])
        mesh.insertNextCell(NORM_HEXA8,8,conn[58:66])
        mesh.insertNextCell(NORM_HEXA8,8,conn[66:74])
        mesh.finishInsertingCells();
        #[1,34,29,23,41,32]
        conn2=[0,5,1,14,18,17,8,7,4,9,5,2, 12,8,9,13,6,7,8,9]
        mesh2.allocateCells(6);
        for i in range(4):
            mesh2.insertNextCell(NORM_TRI3,3,conn2[3*i:3*i+3])
            pass
        mesh2.insertNextCell(NORM_QUAD4,4,conn2[12:16])
        mesh2.insertNextCell(NORM_QUAD4,4,conn2[16:20])
        mesh2.finishInsertingCells();
        return [mesh,mesh2]
    
    # 2D usecase1 for interpolation Gauss Pt-> Gauss Pt. Coming from ASTER : Please, do not touch
    def buildFieldOnGauss_1(self):
        coo=DataArrayDouble([1.0,0.0,1.33333333333333,0.0,1.66666666666667,0.0,0.923879532511287,0.38268343236509006,1.23183937668172,0.510244576486786,1.53979922085214,0.6378057206084831,2.0,0.0,1.8477590650225701,0.7653668647301801,0.9428090415820631,0.9428090415820631,1.1785113019775801,1.1785113019775801,1.4142135623731,1.41421356237309,0.707106781186548,0.707106781186547,0.38268343236509006,0.923879532511287,0.510244576486786,1.23183937668172,0.6378057206084831,1.53979922085214,0.7653668647301801,1.8477590650225701,3.1550283219328204e-17,1.33333333333333,1.16009632455949e-17,1.66666666666667,-2.7620050344068196e-16,2.0,-1.3810025172034098e-16,1.0,-2.0,0.0,-1.53979922085214,0.6378057206084831,-1.66666666666667,0.0,-1.33333333333333,0.0,-0.923879532511287,0.38268343236509006,-1.8477590650225701,0.7653668647301801,-0.9428090415820631,0.9428090415820631,-1.23183937668172,0.510244576486786,-1.83333333333333,0.0,-1.6937791429373599,0.701586292669331,-1.5,0.0,-1.30771370720431,0.26012042935483803,-1.16666666666667,0.0,-1.0778594545965,0.44646400442593803,-1.38578268717091,0.9259503883660041,-1.38581929876693,0.574025148547635,-1.06066017177982,1.06066017177982,-0.8314696123025451,0.5555702330196021,-1.0,0.0,-1.1785113019775801,1.1785113019775801,-0.707106781186548,0.707106781186547,-1.63464213400538,0.325150536693547,-1.9615705608064598,0.390180644032256,-1.47117792060485,0.292635483024192,-0.9807852804032301,0.19509032201612803,-1.524360955888,1.0185454272026,-1.2963624321753402,1.2963624321753402,-1.10862614973673,0.740760310692803,-0.970047881019636,0.6481652718562021,-0.824957911384305,0.824957911384305,-1.4142135623731,1.41421356237309,-1.7981063474059198,0.357665590362902,-1.1442494938037702,0.227605375685483,-1.66293922460509,1.1111404660392,-1.24720441845382,0.833355349529403,-0.7653668647301801,1.8477590650225701,-0.6378057206084831,1.53979922085214,-0.510244576486786,1.23183937668172,-0.701586292669331,1.6937791429373599,-0.574025148547635,1.38581929876693,-0.44646400442593803,1.0778594545965,-0.38268343236509006,0.923879532511287,-0.9259503883660041,1.38578268717091,-0.740760310692803,1.10862614973673,-0.5555702330196021,0.8314696123025451,-0.325150536693547,1.63464213400538,-0.26012042935483803,1.30771370720431,-0.19509032201612803,0.9807852804032301,1.6805133673525298e-18,1.83333333333333,-2.4643915380595496e-16,1.5,-1.4799359654427099e-16,1.16666666666667,-1.1111404660392,1.66293922460509,-0.39018064403225705,1.9615705608064598],73,2)
        coo.setInfoOnComponents(["X [INCONNUE]","Y [INCONNUE]"])
        m=MEDCouplingUMesh("MA1",2)
        m.setDescription("CREE PAR CODE_ASTER") ; m.setTimeUnit("SANS UNITES") ; m.setTime(-1.,-1,-1)
        m.setCoords(coo)
        m.allocateCells()
        conn=[[11,8,13],[11,13,12],[8,9,13],[9,14,13],[9,10,15],[9,15,14],[12,13,19],[13,16,19],[13,14,17],[13,17,16],[14,15,17],[15,18,17],[0,1,4,3],[1,2,5,4],[2,6,7,5],[3,4,8,11],[4,5,9,8],[5,7,10,9],[20,22,21,28,41,51],[21,25,20,29,42,51],[22,23,21,30,43,41],[23,27,21,31,35,43],[23,38,24,32,44,52],[24,27,23,33,31,52],[25,21,50,29,45,53],[21,39,50,34,46,45],[21,27,26,35,47,54],[26,39,21,36,34,54],[27,24,26,33,48,47],[24,40,26,37,49,48],[50,39,56,55,46,62,58,71],[39,26,57,56,36,63,59,62],[26,40,61,57,49,64,60,63],[55,56,17,18,58,65,68,72],[56,57,16,17,59,66,69,65],[57,61,19,16,60,67,70,66]]
        for i in range(0, 12):
            m.insertNextCell(NORM_TRI3,conn[i])
            pass
        for i in range(12, 18):
            m.insertNextCell(NORM_QUAD4,conn[i])
            pass
        for i in range(18, 30):
            m.insertNextCell(NORM_TRI6,conn[i])
            pass
        for i in range(30, 36):
            m.insertNextCell(NORM_QUAD8,conn[i])
            pass
        fff=MEDCouplingFieldDouble.New(ON_GAUSS_PT) ; fff.setName("CH1RB") ; fff.setNature(IntensiveMaximum)
        fff.setMesh(m)
        fff.setGaussLocalizationOnCells(list(range(0, 12)), [0., 0., 1., 0., 0., 1.], [0.3333333333333333, 0.3333333333333333], [0.5])
        fff.setGaussLocalizationOnCells(list(range(12, 18)), [-1., -1., 1., -1., 1., 1., -1., 1.], [-0.577350269189626, -0.577350269189626, 0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626, -0.577350269189626, 0.577350269189626], [1., 1., 1., 1.])
        fff.setGaussLocalizationOnCells(list(range(18, 30)), [0., 0., 1., 0., 0., 1., 0.5, 0., 0.5, 0.5, 0., 0.5], [0.16666666666666666, 0.16666666666666666, 0.6666666666666666, 0.16666666666666666, 0.16666666666666666, 0.6666666666666666], [0.16666666666666666, 0.16666666666666666, 0.16666666666666666])
        fff.setGaussLocalizationOnCells(list(range(30, 36)), [-1., -1., 1., -1., 1., 1., -1., 1., 0., -1., 1., 0., 0., 1., -1., 0.], [-0.774596669241483, -0.774596669241483, 0.774596669241483, -0.774596669241483, 0.774596669241483, 0.774596669241483, -0.774596669241483, 0.774596669241483, 0.0, -0.774596669241483, 0.774596669241483, 0.0, 0.0, 0.774596669241483, -0.774596669241483, 0.0, 0.0, 0.0], [0.30864197530864196, 0.30864197530864196, 0.30864197530864196, 0.30864197530864196, 0.49382716049382713, 0.49382716049382713, 0.49382716049382713, 0.49382716049382713, 0.7901234567901234])
        return MEDCouplingFieldTemplate(fff)

    # 2D usecase2 for interpolation Gauss Pt-> Gauss Pt. Coming from ASTER : Please, do not touch
    def buildFieldOnGauss_2(self):
        coo=DataArrayDouble([1.0,0.0,1.24068268955165,0.15233667925643402,1.25,0.0,1.5,0.0,1.73695576537231,0.21327135095900804,1.75,0.0,2.0,0.0,0.9925461516413221,0.12186934340514702,1.212869657845,0.302402369499585,1.48881922746198,0.182804015107721,1.6980175209829897,0.423363317299419,1.9850923032826397,0.243738686810295,0.9702957262759959,0.241921895599668,1.1669755331215,0.447959936931625,1.4554435894139899,0.362882843399502,1.6337657463701,0.627143911704275,1.94059145255199,0.483843791199335,0.9335804264972021,0.35836794954530005,1.10368449107366,0.586839453482364,1.4003706397458,0.5375519243179501,1.5451582875031202,0.8215752348753091,1.8671608529944002,0.716735899090601,0.882947592858927,0.46947156278589103,1.02394005536124,0.716970545438808,1.32442138928839,0.704207344178836,1.43351607750574,1.00375876361433,1.76589518571785,0.9389431255717821,1.125,0.0,1.11661442059649,0.137103011330791,1.375,0.0,1.4972021976328,0.09157280930228531,1.625,0.0,1.61288749641715,0.198037683033365,1.875,0.0,1.9962695968437298,0.12209707906971401,1.0915826920605,0.272162132549626,1.36475095850682,0.167570347182078,1.47488236134593,0.27335328823822097,1.5767305551984903,0.39312308034946003,1.8610240343274802,0.228505018884652,1.96650981512791,0.364471050984295,1.05027797980935,0.403163943238463,1.3341566236295,0.332642606449543,1.43057542612234,0.45105869925641007,1.5170681930579497,0.5823479180111131,1.8193044867674901,0.45360355424937704,1.9074339014964499,0.601411599008546,0.993316041966293,0.528155508134127,1.28367308643365,0.492755930624788,1.36494190631481,0.6220398639843591,1.43478983839576,0.7628912895270731,1.7504632996822498,0.671939905397438,1.81992254175309,0.8293864853124782,0.921546049825116,0.645273490894927,1.21405294018102,0.6455233988306001,1.27896024653114,0.783747847073923,1.33112207196961,0.93206170907045,1.65552673661049,0.8802591802235451,1.70528032870818,1.0449971294319,0.8191520442889921,0.5735764363510459,1.22872806643349,0.8603646545265691,1.6383040885779798,1.14715287270209,1.24766849802733,0.0763106744185711,0.9981347984218671,0.0610485395348569,1.37243534783007,0.0839417418604282,1.74673589723827,0.106834944186,1.871502747041,0.114466011627857,1.22906863445494,0.227794406865184,0.9832549075639551,0.18223552549214703,1.3519754979004401,0.25057384755170303,1.7206960882369198,0.318912169611258,1.84360295168241,0.341691610297777,1.19214618843528,0.37588224938034104,0.953716950748227,0.300705799504273,1.31136080727881,0.413470474318376,1.6690046638094,0.526235149132478,1.7882192826529297,0.563823374070512,1.13745158859568,0.518366553320299,0.9099612708765431,0.4146932426562391,1.25119674745525,0.570203208652329,1.5924322240339497,0.725713174648418,1.7061773828935198,0.777549829980448,1.06580020544262,0.6531232058949361,0.8526401643540921,0.522498564715949,1.17238022598688,0.7184355264844301,1.12633406089736,0.7886675999826881,1.49212028761966,0.91437248825291,1.59870030816392,0.979684808842404,1.53591008304186,1.07545581815821,1.1229016482246,0.068679606976714,1.6219690474355302,0.0992038767441424,1.10616177100945,0.205014966178666,1.59778922479143,0.29613272892474,1.07293156959176,0.338294024442307,1.5497900449658701,0.488646924194444,1.02370642973611,0.466529897988269,1.47868706517438,0.673876519316388,0.9592201848983541,0.587810885305442,1.3855402670754,0.8490601676634171,0.743144825477394,0.669130606358858,0.9289310318467431,0.836413257948573,1.11471723821609,1.00369590953829,1.30050344458544,1.170978561128,0.656059028990507,0.7547095802227721,0.820073786238134,0.943386975278465,0.984088543485761,1.13206437033416,1.14810330073339,1.32074176538985,0.559192903470747,0.8290375725550421,0.6989911293384331,1.0362969656938,0.8387893552061201,1.24355635883256,0.978587581073807,1.45081575197132,0.453990499739547,0.8910065241883681,0.567488124674433,1.11375815523546,0.6809857496093201,1.3365097862825501,0.794483374544207,1.55926141732964,0.8360379286620679,0.7527719321537151,1.0218241350314199,0.92005458374343,1.20761034140077,1.08733723533314,1.39339654777011,1.25461988692286,0.7380664076143211,0.8490482777506181,0.902081164861948,1.03772567280631,1.06609592210957,1.226403067862,1.2301106793572,1.4150804629177,0.6290920164045901,0.932667269124422,0.7688902422722771,1.13992666226318,0.9086884681399641,1.34718605540194,1.04848669400765,1.5544454485407,0.51073931220699,1.00238233971191,0.624236937141877,1.22513397075901,0.737734562076764,1.4478856018061,0.85123218701165,1.6706372328531902,1.48628965095479,1.33826121271772,1.31211805798101,1.5094191604455398,1.11838580694149,1.6580751451100801,0.907980999479094,1.7820130483767398,0.978260196065517,0.778143295797024,1.17391223527862,0.9337719549564292,1.36956427449172,1.08940061411583,1.56521631370483,1.24502927327524,0.876136580374814,0.891563061442727,1.05136389644978,1.06987567373127,1.22659121252474,1.2481882860198201,1.4018185285997,1.42650089830836,0.7609517862609011,0.991691675364044,0.913142143513081,1.19003001043685,1.06533250076526,1.38836834550966,1.21752285801744,1.5867066805824699,0.6344229537008801,1.07703645055191,0.7613075444410561,1.29244374066229,0.8881921351812321,1.50785103077267,1.01507672592141,1.7232583208830499,0.498436336156558,1.1463250929814,0.5981236033878691,1.37559011157769,0.697810870619181,1.60485513017397,0.7974981378504931,1.8341201487702499,0.42752517915708604,1.17461577598239,0.513030214988503,1.4095389311788602,0.59853525081992,1.6444620863753399,0.6840402866513371,1.87938524157182,0.38477266124137705,1.05715419838415,0.470277697072795,1.29207735358062,0.5557827329042121,1.5270005087771,0.6412877687356291,1.76192366397358,0.34202014332566905,0.9396926207859091,0.782608156852414,0.6225146366376201,0.7009092642998511,0.713250449154182,0.608761429008721,0.7933533402912349,0.507538362960704,0.861629160441526,0.398749068925246,0.917060074385124,-2.0,0.0,-1.75,0.0,-1.5,0.0,-1.25,0.0,-1.9632543668953297,0.38161799075309005,-1.71784757103341,0.333915741908953,-1.4724407751715,0.286213493064817,-1.22703397930958,0.23851124422068104,-1.85436770913357,0.749213186831824,-1.62257174549188,0.655561538477846,-1.39077578185018,0.561909890123868,-1.15897981820848,0.46825824176988995,-1.6773411358908499,1.08927807003005,-1.4676734939044902,0.953118311276297,-1.25800585191814,0.816958552522541,-1.04833820993178,0.680798793768784,-1.4386796006773,1.38931674091799,-1.25884465059264,1.21565214830325,-1.07900970050798,1.0419875556885,-0.8991747504233141,0.868322963073747,-1.0,0.0,-0.981627183447664,0.19080899537654503,-0.9271838545667871,0.374606593415912,-0.838670567945424,0.544639035015027,-0.7193398003386511,0.694658370458997,-1.00375876361433,1.43351607750574,-0.8603646545265691,1.22872806643349,-0.716970545438808,1.02394005536124,-0.5735764363510459,0.8191520442889921,-1.14715287270209,1.6383040885779798,-0.8134732861516011,1.8270909152852002,-0.71178912538265,1.59870455087455,-0.6101049646137,1.3703181864639,-0.50842080384475,1.14193182205325,-0.4067366430758,0.9135454576426011,-0.44990210868773,1.9487401295704703,-0.39366434510176407,1.70514761337416,-0.337426581515798,1.4615550971778501,-0.281188817929831,1.21796258098154,-0.224951054343865,0.974370064785235,-0.06979899340500181,1.9987816540381902,-0.0610741192293767,1.74893394728342,-0.0523492450537515,1.49908624052864,-0.0436243708781263,1.24923853377387,-0.03489949670250091,0.9993908270190961,0.312868930080462,1.97537668119028,0.27376031382040406,1.72845459604149,0.23465169756034704,1.48153251089271,0.19554308130028902,1.23461042574392,0.156434465040231,0.9876883405951381],219,2)
        coo.setInfoOnComponents(["X [INCONNUE]","Y [INCONNUE]"])
        m=MEDCouplingUMesh("MA2",2)
        m.setDescription("CREE PAR CODE_ASTER") ; m.setTimeUnit("SANS UNITES") ; m.setTime(-1.,-1,-1)
        m.setCoords(coo)
        m.allocateCells(0)
        conn=[[198,194,200],[198,200,199],[194,195,200],[195,201,200],[195,196,202],[195,202,201],[196,197,202],[197,203,202],[199,200,205],[199,205,204],[200,201,205],[201,206,205],[201,202,207],[201,207,206],[202,203,207],[203,208,207],[204,205,210],[204,210,209],[205,206,210],[206,211,210],[206,207,212],[206,212,211],[207,208,212],[208,213,212],[209,210,215],[209,215,214],[210,211,215],[211,216,215],[211,212,217],[211,217,216],[212,213,217],[213,218,217],[214,215,157],[214,157,158],[215,216,157],[216,156,157],[216,217,155],[216,155,156],[217,218,155],[218,163,155],[169,170,174,173],[170,171,175,174],[171,172,176,175],[172,189,190,176],[173,174,178,177],[174,175,179,178],[175,176,180,179],[176,190,191,180],[177,178,182,181],[178,179,183,182],[179,180,184,183],[180,191,192,184],[181,182,186,185],[182,183,187,186],[183,184,188,187],[184,192,193,188],[185,186,194,198],[186,187,195,194],[187,188,196,195],[188,193,197,196],[0,2,1,27,62,89],[1,7,0,28,63,89],[2,3,1,29,64,62],[3,9,1,30,36,64],[3,5,4,31,65,90],[4,9,3,32,30,90],[5,6,4,33,66,65],[6,11,4,34,39,66],[7,1,8,28,67,91],[8,12,7,35,68,91],[1,9,8,36,69,67],[9,14,8,37,42,69],[9,4,10,32,70,92],[10,14,9,38,37,92],[4,11,10,39,71,70],[11,16,10,40,45,71],[12,8,13,35,72,93],[13,17,12,41,73,93],[8,14,13,42,74,72],[14,19,13,43,48,74],[14,10,15,38,75,94],[15,19,14,44,43,94],[10,16,15,45,76,75],[16,21,15,46,51,76],[17,13,18,41,77,95],[18,22,17,47,78,95],[13,19,18,48,79,77],[19,24,18,49,54,79],[19,15,20,44,80,96],[20,24,19,50,49,96],[15,21,20,51,81,80],[21,26,20,52,57,81],[22,18,23,47,82,97],[23,59,22,53,83,97],[18,24,23,54,84,82],[24,60,23,55,85,84],[24,20,25,50,86,98],[25,60,24,56,55,98],[20,26,25,57,87,86],[26,61,25,58,88,87],[59,23,100,99,53,135,115,164],[23,60,101,100,85,136,116,135],[60,25,102,101,56,137,117,136],[25,61,131,102,88,138,118,137],[99,100,104,103,115,139,119,165],[100,101,105,104,116,140,120,139],[101,102,106,105,117,141,121,140],[102,131,132,106,118,142,122,141],[103,104,108,107,119,143,123,166],[104,105,109,108,120,144,124,143],[105,106,110,109,121,145,125,144],[106,132,133,110,122,146,126,145],[107,108,112,111,123,147,127,167],[108,109,113,112,124,148,128,147],[109,110,114,113,125,149,129,148],[110,133,134,114,126,150,130,149],[111,112,155,163,127,151,159,168],[112,113,156,155,128,152,160,151],[113,114,157,156,129,153,161,152],[114,134,158,157,130,154,162,153]]
        for i in range(0, 40):
            m.insertNextCell(NORM_TRI3,conn[i])
            pass
        for i in range(40, 60):
            m.insertNextCell(NORM_QUAD4,conn[i])
            pass
        for i in range(60, 100):
            m.insertNextCell(NORM_TRI6,conn[i])
            pass
        for i in range(100, 120):
            m.insertNextCell(NORM_QUAD8,conn[i])
            pass
        fff=MEDCouplingFieldDouble.New(ON_GAUSS_PT) ; fff.setName("CH2RB") ; fff.setNature(IntensiveMaximum)
        fff.setMesh(m)
        fff.setGaussLocalizationOnCells(list(range(0, 40)), [0., 0., 1., 0., 0., 1.], [0.3333333333333333, 0.3333333333333333], [0.5])
        fff.setGaussLocalizationOnCells(list(range(40, 60)), [-1., -1., 1., -1., 1., 1., -1., 1.], [-0.577350269189626, -0.577350269189626, 0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626, -0.577350269189626, 0.577350269189626], [1., 1., 1., 1.])
        fff.setGaussLocalizationOnCells(list(range(60, 100)), [0., 0., 1., 0., 0., 1., 0.5, 0., 0.5, 0.5, 0., 0.5], [0.16666666666666666, 0.16666666666666666, 0.6666666666666666, 0.16666666666666666, 0.16666666666666666, 0.6666666666666666], [0.16666666666666666, 0.16666666666666666, 0.16666666666666666])
        fff.setGaussLocalizationOnCells(list(range(100, 120)), [-1., -1., 1., -1., 1., 1., -1., 1., 0., -1., 1., 0., 0., 1., -1., 0.], [-0.774596669241483, -0.774596669241483, 0.774596669241483, -0.774596669241483, 0.774596669241483, 0.774596669241483, -0.774596669241483, 0.774596669241483, 0.0, -0.774596669241483, 0.774596669241483, 0.0, 0.0, 0.774596669241483, -0.774596669241483, 0.0, 0.0, 0.0], [0.30864197530864196, 0.30864197530864196, 0.30864197530864196, 0.30864197530864196, 0.49382716049382713, 0.49382716049382713, 0.49382716049382713, 0.49382716049382713, 0.7901234567901234])
        return MEDCouplingFieldTemplate(fff)

    # 3D usecase1 for interpolation Gauss Pt-> Gauss Pt. Coming from ASTER : Please, do not touch
    def buildFieldOnGauss_3(self):
        coo=DataArrayDouble([0.,1.,0.,0.,2.,0.,0.,3.,0.,1.,1.,0.,1.,2.,0.,1.,3.,0.,0.,1.,1.,0.,3.,1.,0.5,1.,1.,0.5,3.,1.,1.,1.,1.,1.,3.,1.,0.,0.,0.,0.,1.,0.,1.,0.,0.,1.,1.,0.,0.,0.,1.,0.,0.5,1.,0.,1.,1.,0.5,0.,1.,0.5,0.5,1.,0.5,1.,1.,1.,0.,1.,1.,0.5,1.,1.,1.,1.0],25,3)
        coo.setInfoOnComponents(["X [INCONNUE]","Y [INCONNUE]","Z [INCONNUE]"])
        m=MEDCouplingUMesh("MA1",3)
        m.setDescription("CREE PAR CODE_ASTER") ; m.setTimeUnit("SANS UNITES") ; m.setTime(-1.,-1,-1)
        m.setCoords(coo)
        m.allocateCells(0)
        conn=[[3,10,8,4],[19,22,23,20,14],[0,6,1,3,8,4],[4,8,10,5,9,11],[12,16,17,14,19,20],[14,20,23,15,21,24],[1,2,5,4,6,7,9,8],[12,13,15,14,17,18,21,20]]
        m.insertNextCell(NORM_TETRA4,conn[0])
        m.insertNextCell(NORM_PYRA5,conn[1])
        for i in range(2, 6):
            m.insertNextCell(NORM_PENTA6,conn[i])
            pass
        m.insertNextCell(NORM_HEXA8,conn[6])
        m.insertNextCell(NORM_HEXA8,conn[7])
        fff=MEDCouplingFieldDouble.New(ON_GAUSS_PT) ; fff.setName("CH13") ; fff.setNature(IntensiveMaximum)
        fff.setMesh(m)
        fff.setGaussLocalizationOnCells([0],[0.,1.,0.,0.,0.,0.,0.,0.,1.,1.,0.,0.],[0.25,0.25,0.25],[0.16666666666666666])
        fff.setGaussLocalizationOnCells([1],[1.,0.,0.,0.,-1.,0.,-1.,0.,0.,0.,1.,0.,0.,0.,1.],[0.5,0.,0.1531754163448146,0.,0.5,0.1531754163448146,-0.5,0.,0.1531754163448146,0.,-0.5,0.1531754163448146,0.,0.,0.6372983346207416],[0.1333333333333333,0.1333333333333333,0.1333333333333333,0.1333333333333333,0.1333333333333333])
        fff.setGaussLocalizationOnCells([2,3,4,5],[-1.,1.,0.,-1.,0.,0.,-1.,0.,1.,1.,1.,0.,1.,0.,0.,1.,0.,1.],[-0.577350269189626,0.5,0.5,-0.577350269189626,0.,0.5,-0.577350269189626,0.5,0.,0.577350269189626,0.5,0.5,0.577350269189626,0.,0.5,0.577350269189626,0.5,0.],[0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666])
        fff.setGaussLocalizationOnCells([6,7],[-1.,-1.,-1.,-1.,1.,-1.,1.,1.,-1.,1.,-1.,-1.,-1.,-1.,1.,-1.,1.,1.,1.,1.,1.,1.,-1.,1.],[-0.577350269189626,-0.577350269189626,-0.577350269189626,-0.577350269189626,-0.577350269189626,0.577350269189626,-0.577350269189626,0.577350269189626,-0.577350269189626,-0.577350269189626,0.577350269189626,0.577350269189626,0.577350269189626,-0.577350269189626,-0.577350269189626,0.577350269189626,-0.577350269189626,0.577350269189626,0.577350269189626,0.577350269189626,-0.577350269189626,0.577350269189626,0.577350269189626,0.577350269189626],[1.,1.,1.,1.,1.,1.,1.,1.])
        return MEDCouplingFieldTemplate(fff)

    # 3D usecase2 for interpolation Gauss Pt-> Gauss Pt. Coming from ASTER : Please, do not touch
    def buildFieldOnGauss_4(self):
        coo=DataArrayDouble([0.,2.,0.,0.,1.,0.,0.,0.,0.,1.,2.,0.,1.,1.,0.,1.,0.,0.,0.,2.,1.,0.,0.,1.,0.5, 2.,1.,0.5, 0.,1.,1.,2.,1.,1.,0.,1.,0.,3.,0.,0.,2.,0.,1.,3.,0.,1.,2.,0.,0.,3.,1.,0.,2.5, 1.,0.,2.,1.,0.5, 3.,1.,0.5, 2.5, 1.,0.5, 2.,1.,1.,3.,1.,1.,2.5, 1.,1.,2.,1.0],25,3)
        coo.setInfoOnComponents(["X [INCONNUE]","Y [INCONNUE]","Z [INCONNUE]"])
        m=MEDCouplingUMesh("MA2",3)
        m.setDescription("CREE PAR CODE_ASTER") ; m.setTimeUnit("SANS UNITES") ; m.setTime(-1.,-1,-1)
        m.setCoords(coo)
        m.allocateCells(0)
        conn=[[3,10,8,4],[19,22,23,20,14],[0,6,1,3,8,4],[4,8,10,5,9,11],[12,16,17,14,19,20],[14,20,23,15,21,24],[1,2,5,4,6,7,9,8],[12,13,15,14,17,18,21,20]]
        m.insertNextCell(NORM_TETRA4,conn[0])
        m.insertNextCell(NORM_PYRA5,conn[1])
        for i in range(2, 6):
            m.insertNextCell(NORM_PENTA6,conn[i])
            pass
        m.insertNextCell(NORM_HEXA8,conn[6])
        m.insertNextCell(NORM_HEXA8,conn[7])
        fff=MEDCouplingFieldDouble.New(ON_GAUSS_PT) ; fff.setName("CH23") ; fff.setNature(IntensiveMaximum)
        fff.setMesh(m)
        fff.setGaussLocalizationOnCells([0],[0.,1.,0.,0.,0.,0.,0.,0.,1.,1.,0.,0.],[0.25,0.25,0.25],[0.16666666666666666])
        fff.setGaussLocalizationOnCells([1],[1.,0.,0.,0.,-1.,0.,-1.,0.,0.,0.,1.,0.,0.,0.,1.],[0.5,0.,0.1531754163448146,0.,0.5,0.1531754163448146,-0.5,0.,0.1531754163448146,0.,-0.5,0.1531754163448146,0.,0.,0.6372983346207416],[0.1333333333333333,0.1333333333333333,0.1333333333333333,0.1333333333333333,0.1333333333333333])
        fff.setGaussLocalizationOnCells([2,3,4,5],[-1.,1.,0.,-1.,0.,0.,-1.,0.,1.,1.,1.,0.,1.,0.,0.,1.,0.,1.],[-0.577350269189626,0.5,0.5,-0.577350269189626,0.,0.5,-0.577350269189626,0.5,0.,0.577350269189626,0.5,0.5,0.577350269189626,0.,0.5,0.577350269189626,0.5,0.],[0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666])
        fff.setGaussLocalizationOnCells([6,7],[-1.,-1.,-1.,-1.,1.,-1.,1.,1.,-1.,1.,-1.,-1.,-1.,-1.,1.,-1.,1.,1.,1.,1.,1.,1.,-1.,1.],[-0.577350269189626,-0.577350269189626,-0.577350269189626,-0.577350269189626,-0.577350269189626,0.577350269189626,-0.577350269189626,0.577350269189626,-0.577350269189626,-0.577350269189626,0.577350269189626,0.577350269189626,0.577350269189626,-0.577350269189626,-0.577350269189626,0.577350269189626,-0.577350269189626,0.577350269189626,0.577350269189626,0.577350269189626,-0.577350269189626,0.577350269189626,0.577350269189626,0.577350269189626],[1.,1.,1.,1.,1.,1.,1.,1.])
        return MEDCouplingFieldTemplate(fff)

    def buildCircle(self, center_X, center_Y, radius):
      from cmath import rect
      from math import pi  
  
      c = [rect(radius, i * pi / 4.0) for i in range(8)]
      coords = [c[-1].real,c[-1].imag,  c[3].real,c[3].imag,
                 c[5].real,c[5].imag,  c[1].real,c[1].imag]
      connec = list(range(4))
      baseMesh = MEDCouplingUMesh.New("circle", 2)  
      baseMesh.allocateCells(1)
      meshCoords = DataArrayDouble.New(coords, len(coords) // 2, 2)
      meshCoords += (center_X, center_Y)
      baseMesh.setCoords(meshCoords)
  
      baseMesh.insertNextCell(NORM_QPOLYG, connec)  
      baseMesh.finishInsertingCells()  
      return baseMesh

    def buildCircle2(self, center_X, center_Y, radius):  
      from cmath import rect
      from math import pi  
  
      c = [rect(radius, i * pi / 4.0) for i in range(8)]
      coords = []
      for i in range(8):
          coords.extend([c[i].real,c[i].imag])
      connec = [7,5,3,1,  6,4,2,0]
      baseMesh = MEDCouplingUMesh.New("circle", 2)  
      baseMesh.allocateCells(1)
      meshCoords = DataArrayDouble.New(coords, len(coords) // 2, 2)
      meshCoords += (center_X, center_Y)
      baseMesh.setCoords(meshCoords)
  
      baseMesh.insertNextCell(NORM_QPOLYG, connec)  
      baseMesh.finishInsertingCells()  
      return baseMesh  

    build2DTargetMesh_1=classmethod(build2DTargetMesh_1)
    build2DSourceMesh_1=classmethod(build2DSourceMesh_1)
    build3DTargetMesh_1=classmethod(build3DTargetMesh_1)
    build3DSourceMesh_1=classmethod(build3DSourceMesh_1)
    build3DSurfTargetMesh_1=classmethod(build3DSurfTargetMesh_1)
    build3DExtrudedUMesh_1=classmethod(build3DExtrudedUMesh_1)
    buildCU1DMesh_U=classmethod(buildCU1DMesh_U)
    build2DTargetMeshMergeNode_1=classmethod(build2DTargetMeshMergeNode_1)
    build3DTargetMeshMergeNode_1=classmethod(build3DTargetMeshMergeNode_1)
    build2DTargetMeshMerged_1=classmethod(build2DTargetMeshMerged_1)
    build2DTargetMesh_2=classmethod(build2DTargetMesh_2)
    build1DSourceMesh_2=classmethod(build1DSourceMesh_2)
    build1DTargetMesh_3=classmethod(build1DTargetMesh_3)
    build2DCurveTargetMesh_3=classmethod(build2DCurveTargetMesh_3)
    build2DTargetMesh_3=classmethod(build2DTargetMesh_3)
    build2DTargetMesh_4=classmethod(build2DTargetMesh_4)
    buildMultiFields_1=classmethod(buildMultiFields_1)
    buildMultiFields_2=classmethod(buildMultiFields_2)
    build1DMultiTypes_1=classmethod(build1DMultiTypes_1)
    build2DMultiTypes_1=classmethod(build2DMultiTypes_1)
    build3DMultiTypes_1=classmethod(build3DMultiTypes_1)
    buildCoordsForMultiTypes_1=classmethod(buildCoordsForMultiTypes_1)
    buildHexa8Mesh_1=classmethod(buildHexa8Mesh_1)
    buildPointe_1=classmethod(buildPointe_1)
    buildFieldOnGauss_1=classmethod(buildFieldOnGauss_1)
    buildFieldOnGauss_2=classmethod(buildFieldOnGauss_2)
    buildFieldOnGauss_3=classmethod(buildFieldOnGauss_3)
    buildFieldOnGauss_4=classmethod(buildFieldOnGauss_4)
    buildCircle=classmethod(buildCircle)
    buildCircle2=classmethod(buildCircle2)
    pass



