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

from libMEDCoupling_Swig import *

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
            0.,0.,0., 1.,1.,0., 1.,1.25,0., 0.,1.,0., 1.,1.5,0., 2.,0.,0., 2.,1.,0., 1.,2.,0., 0.,2.,0., 3.,1.,0.,
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

    build2DTargetMesh_1=classmethod(build2DTargetMesh_1)
    build2DSourceMesh_1=classmethod(build2DSourceMesh_1)
    build3DTargetMesh_1=classmethod(build3DTargetMesh_1)
    build3DSurfTargetMesh_1=classmethod(build3DSurfTargetMesh_1)
    build3DExtrudedUMesh_1=classmethod(build3DExtrudedUMesh_1)
    buildCU1DMesh_U=classmethod(buildCU1DMesh_U)
    build2DTargetMeshMergeNode_1=classmethod(build2DTargetMeshMergeNode_1)
    build3DTargetMeshMergeNode_1=classmethod(build3DTargetMeshMergeNode_1)
    build2DTargetMeshMerged_1=classmethod(build2DTargetMeshMerged_1)
    build2DTargetMesh_2=classmethod(build2DTargetMesh_2)
    pass
