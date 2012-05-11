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
        for i in xrange(8):
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
        for i in xrange(8):
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
        for i in xrange(4):
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
        ret.checkCoherency();
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
        for i in xrange(8):
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
        for i in xrange(12):
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
        for i in xrange(4):
            mesh2.insertNextCell(NORM_TRI3,3,conn2[3*i:3*i+3])
            pass
        mesh2.insertNextCell(NORM_QUAD4,4,conn2[12:16])
        mesh2.insertNextCell(NORM_QUAD4,4,conn2[16:20])
        mesh2.finishInsertingCells();
        return [mesh,mesh2]

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
    pass
