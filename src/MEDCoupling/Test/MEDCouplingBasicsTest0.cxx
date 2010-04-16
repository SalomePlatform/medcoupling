//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#include "MEDCouplingBasicsTest.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"

#include "MEDCouplingBasicsTestData1.hxx"

using namespace ParaMEDMEM;

MEDCouplingUMesh *MEDCouplingBasicsTest::build3DSourceMesh_2()
{
  double sourceCoords[84]={100.0, 100.0, 0.0, 100.0, 100.0, 100.0, 100.0, 0.0, 100.0, 100.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 100.0, 100.0, 0.0,
                           0.0, 100.0, 0.0, 0.0, 0.0, 100.0, 100.0, 200.0, 100.0, 0.0, 200.0, 0.0, 100.0, 200.0, 0.0, 0.0, 200.0, 100.0, 200.0,
                           0.0, 100.0, 200.0, 100.0, 0.0, 200.0, 0.0, 0.0, 200.0, 100.0, 100.0, 200.0, 200.0, 0.0, 200.0, 200.0, 200.0, 100.0,
                           0.0, 200.0, 100.00000000833332, 100.00000000833332, 200.0, 0.0, 100.0, 200.0, 0.0, 0.0, 200.0, 100.0, 200.0, 200.0,
                           0.0, 200.0, 200.0, 200.0, 0.0, 200.0, 200.0, 100.0, 200.0, 200.0, 200.0, 149.999999970343, 149.9999999874621, 49.999999881628682};
  
  
  int sourceConn[212]={25, 27, 13, 19, 18, 3, 20, 21, 5, 10, 17, 1, 1, 3, 0, 7, 18, 1, 0, 27, 12, 27, 13, 24, 25, 19, 16, 26, 1, 2, 6, 8, 15, 13, 
                       12, 5, 24, 13, 25, 27, 10, 11, 9, 6, 19, 8, 23, 1, 22, 8, 23, 19, 16, 13, 17, 1, 6, 9, 10, 8, 13, 17, 5, 15, 5, 4, 1, 12, 18,
                       0, 24, 27, 19, 20, 18, 1, 7, 6, 5, 1, 4, 12, 15, 14, 25, 27, 19, 18, 1, 19, 16, 13, 20, 19, 23, 1, 27, 12, 1, 0, 6, 5, 1, 10,
                       4, 5, 1, 7, 12, 27, 1, 13, 5, 15, 4, 12, 19, 16, 26, 22, 13, 5, 17, 1, 1, 3, 7, 2, 13, 5, 1, 12, 18, 1, 3, 0, 8, 23, 2, 9, 3,
                       1, 18, 20, 1, 27, 19, 13, 24, 25, 18, 27, 25, 16, 19, 13, 7, 1, 2, 6, 3, 1, 20, 2, 8, 16, 17, 1, 7, 4, 0, 1, 18, 19, 1, 27,
                       27, 12, 0, 24, 9, 6, 2, 8, 1, 4, 0, 12, 19, 16, 22, 8, 8, 2, 23, 1, 1, 16, 19, 8, 20, 2, 1, 23, 10, 1, 6, 8, 10, 8, 17, 1};
  
  MEDCouplingUMesh *sourceMesh=MEDCouplingUMesh::New();
  sourceMesh->setMeshDimension(3);
  sourceMesh->allocateCells(53);
  for(int i=0;i<53;i++)
    sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,sourceConn+4*i);
  sourceMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(28,3);
  std::copy(sourceCoords,sourceCoords+84,myCoords->getPointer());
  sourceMesh->setCoords(myCoords);
  myCoords->decrRef();
  return sourceMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build3DTargetMesh_2()
{
  double targetCoords[24]={200.0, 200.0, 0.0, 200.0, 200.0, 200.0, 200.0, 0.0, 0.0, 200.0, 0.0, 200.0, 0.0, 200.0, 0.0, 0.0, 200.0, 200.0, 0.0, 0.0, 0.0, 0.0, 0.0, 200.0};
  int targetConn[20]={5, 6, 3, 0, 1, 3, 0, 5, 3, 6, 5, 7, 6, 4, 0, 5, 6, 3, 0, 2};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(3);
  targetMesh->allocateCells(5);
  for(int i=0;i<5;i++)
    targetMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,targetConn+4*i);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(8,3);
  std::copy(targetCoords,targetCoords+24,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build1DTargetMesh_1()
{
  double targetCoords[36]={
    25.,25.,0., 25.,25.,50., 25.,25.,200., 75.,25.,0., 75.,25.,50., 75.,25.,200.,
    25.,125.,0., 25.,125.,50., 25.,125.,200., 125.,125.,0., 125.,125.,50., 125.,125.,200.
  };
  int targetConn[16]={0,1, 1,2, 3,4, 4,5, 6,7, 7,8, 9,10, 10,11};

  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New("my name of mesh 1D",1);
  targetMesh->allocateCells(8);
  for(int i=0;i<8;i++)
    targetMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,targetConn+2*i);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(12,3);
  std::copy(targetCoords,targetCoords+36,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build2DSourceMesh_1()
{
  double sourceCoords[8]={-0.3,-0.3, 0.7,-0.3, -0.3,0.7, 0.7,0.7};
  int sourceConn[6]={0,3,1,0,2,3};
  MEDCouplingUMesh *sourceMesh=MEDCouplingUMesh::New("my name of mesh 2D",2);
  sourceMesh->allocateCells(2);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,sourceConn);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,sourceConn+3);
  sourceMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(4,2);
  std::copy(sourceCoords,sourceCoords+8,myCoords->getPointer());
  sourceMesh->setCoords(myCoords);
  myCoords->decrRef();
  return sourceMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build2DTargetMesh_1()
{
  double targetCoords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  int targetConn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(5);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+4);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+7);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+10);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+14);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(9,2);
  std::copy(targetCoords,targetCoords+18,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build2DTargetMeshPerm_1()
{
  double targetCoords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  int targetConn[18]={0,3,4,1, 1,2,4, 4,5,2, 6,7,4,3, 7,8,5,4};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(5);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+4);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+7);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+10);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+14);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(9,2);
  std::copy(targetCoords,targetCoords+18,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build2DTargetMesh_2()
{
  double targetCoords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  int targetConn[24]={0,3,4, 0,4,1, 1,4,2, 4,5,2, 3,6,4, 6,7,4, 4,7,5, 7,8,5 };
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(8);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+3);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+6);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+9);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+12);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+15);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+18);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+21);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(9,2);
  std::copy(targetCoords,targetCoords+18,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::buildCU1DMesh_U()
{
  double coords[4]={ 0.0, 0.3, 0.75, 1.0 };
  int conn[2*3]={ 0,1, 1,2, 2,3 };
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(1);
  mesh->allocateCells(3);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+2);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+4);
  mesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(4,1);
  std::copy(coords,coords+4,myCoords->getPointer());
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  return mesh;
}
MEDCouplingUMesh *MEDCouplingBasicsTest::buildCU2DMesh_U()
{
  double coords[18]={0.0,0.0, 0.5,0.0, 1.0,0.0, 0.0,0.5, 0.5,0.5, 1.0,0.5, 0.0,1.0, 0.5,1.0, 1.0,1.0 };
  int conn[18]={0,1,4,3, 3,4,7,6, 4,5,8,7, 1,5,4, 1,2,5 };
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+8);
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn+12);
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn+15);
  mesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(9,2);
  std::copy(coords,coords+18,myCoords->getPointer());
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  return mesh;
}
MEDCouplingUMesh *MEDCouplingBasicsTest::buildCU3DMesh_U()
{
  double coords[27*3]=
    {
//   0.0,1.0,0.0 ,0.0,0.3,0.0 ,0.0,0.3,0.3 ,0.3,0.0,0.0 ,0.3,0.3,1.0 ,1.0,0.0,1.0 ,1.0,0.0,0.3 ,0.3,0.0,0.3 ,0.3,1.0,0.3 ,0.0,0.3,1.0 ,0.3,0.0,1.0 ,0.3,0.3,0.3 ,1.0,0.3,1.0 ,1.0,0.0,0.0 ,0.0,0.0,0.0 ,1.0,0.3,0.3 ,0.3,1.0,0.0 ,1.0,1.0,0.3 ,1.0,1.0,1.0 ,0.0,1.0,1.0 ,0.3,0.3,0.0 ,0.0,1.0,0.3 ,0.0,0.0,1.0 ,0.3,1.0,1.0 ,1.0,0.3,0.0 ,0.0,0.0,0.3 ,1.0,1.0,0.0
      0.0,0.0,0.0, 0.3,0.0,0.0, 1.0,0.0,0.0, 0.0,0.3,0.0, 0.3,0.3,0.0, 1.0,0.3,0.0, 0.0,1.0,0.0, 0.3,1.0,0.0, 1.0,1.0,0.0, 0.0,0.0,0.3, 0.3,0.0,0.3, 1.0,0.0,0.3, 0.0,0.3,0.3, 0.3,0.3,0.3, 1.0,0.3,0.3, 0.0,1.0,0.3, 0.3,1.0,0.3, 1.0,1.0,0.3, 0.0,0.0,1.0, 0.3,0.0,1.0, 1.0,0.0,1.0, 0.0,0.3,1.0, 0.3,0.3,1.0, 1.0,0.3,1.0, 0.0,1.0,1.0, 0.3,1.0,1.0, 1.0,1.0,1.0,
    };
  int conn[8*8]=
    {
//       11,15,12,4,8,17,18,23,3,13,6,7,20,24,15,11,14,3,7,25,1,20,11,2,1,20,11,2,0,16,8,21,20,24,15,11,16,26,17,8,25,7,10,22,2,11,4,9,2,11,4,9,21,8,23,19,7,6,5,10,11,15,12,4
      0,3,4,1,9,12,13,10, 1,4,5,2,10,13,14,11, 3,6,7,4,12,15,16,13, 4,7,8,5,13,16,17,14, 9,12,13,10,18,21,22,19, 10,13,14,11,19,22,23,20, 12,15,16,13,21,24,25,22, 13,16,17,14,22,25,26,23
    };
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(3);
  mesh->allocateCells(8);
  for(int i=0;i<8;i++)
    mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+8*i);
  mesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(27,3);
  std::copy(coords,coords+27*3,myCoords->getPointer());
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  return mesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build3DSurfSourceMesh_1()
{
  double sourceCoords[12]={-0.3,-0.3,0.5, 0.7,-0.3,1.5, -0.3,0.7,0.5, 0.7,0.7,1.5};
  int sourceConn[6]={0,3,1,0,2,3};
  MEDCouplingUMesh *sourceMesh=MEDCouplingUMesh::New();
  sourceMesh->setMeshDimension(2);
  sourceMesh->allocateCells(2);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,sourceConn);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,sourceConn+3);
  sourceMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(4,3);
  std::copy(sourceCoords,sourceCoords+12,myCoords->getPointer());
  sourceMesh->setCoords(myCoords);
  myCoords->decrRef();
  return sourceMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build3DSurfSourceMesh_2()
{
  double sourceCoords[12]={-0.3,-0.3,0., 0.7,-0.3,0., -0.3,0.7,0., 0.7,0.7,0.};
  int sourceConn[6]={0,3,1,0,2,3};
  MEDCouplingUMesh *sourceMesh=MEDCouplingUMesh::New();
  sourceMesh->setMeshDimension(2);
  sourceMesh->allocateCells(2);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,sourceConn);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,sourceConn+3);
  sourceMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(4,3);
  std::copy(sourceCoords,sourceCoords+12,myCoords->getPointer());
  sourceMesh->setCoords(myCoords);
  myCoords->decrRef();
  return sourceMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build3DSurfTargetMesh_1()
{
  double targetCoords[27]={-0.3,-0.3,0.5, 0.2,-0.3,1., 0.7,-0.3,1.5, -0.3,0.2,0.5, 0.2,0.2,1., 0.7,0.2,1.5, -0.3,0.7,0.5, 0.2,0.7,1., 0.7,0.7,1.5};
  int targetConn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(5);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+4);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+7);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+10);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+14);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(9,3);
  std::copy(targetCoords,targetCoords+27,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

/*!
 * Idem build3DSurfTargetMesh_1 except that cell id 2 is not correctly numbered.
 */
MEDCouplingUMesh *MEDCouplingBasicsTest::build3DSurfTargetMeshPerm_1()
{
  double targetCoords[27]={-0.3,-0.3,0.5, 0.2,-0.3,1., 0.7,-0.3,1.5, -0.3,0.2,0.5, 0.2,0.2,1., 0.7,0.2,1.5, -0.3,0.7,0.5, 0.2,0.7,1., 0.7,0.7,1.5};
  int targetConn[18]={0,3,4,1, 1,4,2, 4,2,5, 6,7,4,3, 7,8,5,4};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(5);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+4);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+7);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+10);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+14);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(9,3);
  std::copy(targetCoords,targetCoords+27,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build3DSurfTargetMesh_2()
{
  double targetCoords[27]={-0.3,-0.3,0.5, 0.2,-0.3,1., 0.7,-0.3,1.5, -0.3,0.2,0.5, 0.2,0.2,1., 0.7,0.2,1.5, -0.3,0.7,0.5, 0.2,0.7,1., 0.7,0.7,1.5};
  int targetConn[24]={0,3,4, 0,4,1, 1,4,2, 4,5,2, 3,6,4, 6,7,4, 4,7,5, 7,8,5 };
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(8);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+3);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+6);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+9);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+12);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+15);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+18);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+21);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(9,3);
  std::copy(targetCoords,targetCoords+27,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build3DSourceMesh_1()
{
  double sourceCoords[27]={ 0.0, 0.0, 200.0, 0.0, 0.0, 0.0, 0.0, 200.0, 200.0, 0.0, 200.0, 0.0, 200.0, 0.0, 200.0,
                            200.0, 0.0, 0.0, 200.0, 200.0, 200.0, 200.0, 200.0, 0.0, 100.0, 100.0, 100.0 };
  int sourceConn[48]={8,1,7,3, 6,0,8,2, 7,4,5,8, 6,8,4,7, 6,8,0,4, 6,8,7,3, 8,1,3,0, 4,1,5,8, 1,7,5,8, 0,3,8,2, 8,1,0,4, 3,6,8,2};
  MEDCouplingUMesh *sourceMesh=MEDCouplingUMesh::New();
  sourceMesh->setMeshDimension(3);
  sourceMesh->allocateCells(12);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,sourceConn);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,sourceConn+4);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,sourceConn+8);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,sourceConn+12);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,sourceConn+16);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,sourceConn+20);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,sourceConn+24);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,sourceConn+28);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,sourceConn+32);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,sourceConn+36);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,sourceConn+40);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,sourceConn+44);
  sourceMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(9,3);
  std::copy(sourceCoords,sourceCoords+27,myCoords->getPointer());
  sourceMesh->setCoords(myCoords);
  myCoords->decrRef();
  return sourceMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build3DTargetMesh_1()
{
  double targetCoords[81]={ 0., 0., 0., 50., 0., 0. , 200., 0., 0.  , 0., 50., 0., 50., 50., 0. , 200., 50., 0.,   0., 200., 0., 50., 200., 0. , 200., 200., 0. ,
                            0., 0., 50., 50., 0., 50. , 200., 0., 50.  , 0., 50., 50., 50., 50., 50. , 200., 50., 50.,   0., 200., 50., 50., 200., 50. , 200., 200., 50. ,
                            0., 0., 200., 50., 0., 200. , 200., 0., 200.  , 0., 50., 200., 50., 50., 200. , 200., 50., 200.,   0., 200., 200., 50., 200., 200. , 200., 200., 200. };
  int targetConn[64]={0,1,4,3,9,10,13,12, 1,2,5,4,10,11,14,13, 3,4,7,6,12,13,16,15, 4,5,8,7,13,14,17,16,
                      9,10,13,12,18,19,22,21, 10,11,14,13,19,20,23,22, 12,13,16,15,21,22,25,24, 13,14,17,16,22,23,26,25};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(3);
  targetMesh->allocateCells(12);
  for(int i=0;i<8;i++)
    targetMesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,targetConn+8*i);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(27,3);
  std::copy(targetCoords,targetCoords+81,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build2DTargetMeshMergeNode_1()
{
  double targetCoords[36]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,-0.3, 0.2,-0.3, 0.2,-0.3, 0.2,0.2, 0.2,0.2, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7, 0.2,0.7 };
  int targetConn[18]={0,9,7,5, 4,6,2, 10,11,8, 9,14,15,7, 17,16,13,6};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(5);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+4);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+7);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+10);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+14);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(18,2);
  std::copy(targetCoords,targetCoords+36,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build3DTargetMeshMergeNode_1()
{
  double targetCoords[93]={ 0., 0., 0., 50., 0., 0. , 200., 0., 0.  , 0., 50., 0., 50., 50., 0. , 200., 50., 0.,   0., 200., 0., 50., 200., 0. , 200., 200., 0. ,
                            0., 0., 50., 50., 0., 50. , 200., 0., 50.  , 0., 50., 50., 50., 50., 50. , 200., 50., 50.,   0., 200., 50., 50., 200., 50. , 200., 200., 50. ,
                            0., 0., 200., 50., 0., 200. , 200., 0., 200.  , 0., 50., 200., 50., 50., 200. , 200., 50., 200.,   0., 200., 200., 50., 200., 200. , 200., 200., 200., 50.,0.,0., 50.,0.,0., 50.,0.,0.,  200., 50., 200.};
  int targetConn[64]={0,29,4,3,9,10,13,12, 28,2,5,4,10,11,14,13, 3,4,7,6,12,13,16,15, 4,5,8,7,13,14,17,16,
                      9,10,13,12,18,19,22,21, 10,11,14,13,19,20,23,22, 12,13,16,15,21,22,25,24, 13,14,17,16,22,30,26,25};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(3);
  targetMesh->allocateCells(12);
  for(int i=0;i<8;i++)
    targetMesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,targetConn+8*i);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(31,3);
  std::copy(targetCoords,targetCoords+93,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build3DExtrudedUMesh_1(MEDCouplingUMesh *&mesh2D)
{
  double coords[180]={
    0.,0.,0., 1.,1.,0., 1.,1.25,0., 0.,1.,0., 1.,1.5,0., 2.,0.,0., 2.,1.,0., 1.,2.,0., 0.,2.,0., 3.,1.,0.,
    3.,2.,0., 0.,1.,0., 1.,3.,0., 2.,2.,0., 2.,3.,0.,
    0.,0.,1., 1.,1.,1., 1.,1.25,1., 0.,1.,1., 1.,1.5,1., 2.,0.,1., 2.,1.,1., 1.,2.,1., 0.,2.,1., 3.,1.,1.,
    3.,2.,1., 0.,1.,1., 1.,3.,1., 2.,2.,1., 2.,3.,1.,
    0.,0.,2., 1.,1.,2., 1.,1.25,2., 0.,1.,2., 1.,1.5,2., 2.,0.,2., 2.,1.,2., 1.,2.,2., 0.,2.,2., 3.,1.,2.,
    3.,2.,2., 0.,1.,2., 1.,3.,2., 2.,2.,2., 2.,3.,2.,
    0.,0.,3., 1.,1.,3., 1.,1.25,3., 0.,1.,3., 1.,1.5,3., 2.,0.,3., 2.,1.,3., 1.,2.,3., 0.,2.,3., 3.,1.,3.,
    3.,2.,3., 0.,1.,3., 1.,3.,3., 2.,2.,3., 2.,3.,3.};

  int conn[354]={
    // 0
    0,11,1,3,15,26,16,18,   1,2,4,7,13,6,-1,1,16,21,6,-1,6,21,28,13,-1,13,7,22,28,-1,7,4,19,22,-1,4,2,17,19,-1,2,1,16,17,-1,16,21,28,22,19,17,
    1,6,5,3,16,21,20,18,   13,10,9,6,28,25,24,21,
    11,8,7,4,2,1,-1,11,26,16,1,-1,1,16,17,2,-1,2,17,19,4,-1,4,19,22,7,-1,7,8,23,22,-1,8,11,26,23,-1,26,16,17,19,22,23,
    7,12,14,13,22,27,29,28,
    // 1
    15,26,16,18,30,41,31,33,   16,17,19,22,28,21,-1,16,31,36,21,-1,21,36,43,28,-1,28,22,37,43,-1,22,19,34,37,-1,19,17,32,34,-1,17,16,31,32,-1,31,36,43,37,34,32,
    16,21,20,18,31,36,35,33,   28,25,24,21,43,40,39,36,
    26,23,22,19,17,16,-1,26,41,31,16,-1,16,31,32,17,-1,17,32,34,19,-1,19,34,37,22,-1,22,23,38,37,-1,23,26,41,38,-1,41,31,32,34,37,38,
    22,27,29,28,37,42,44,43,
    // 2
    30,41,31,33,45,56,46,48,  31,32,34,37,43,36,-1,31,46,51,36,-1,36,51,58,43,-1,43,37,52,58,-1,37,34,49,52,-1,34,32,47,49,-1,32,31,46,47,-1,46,51,58,52,49,47,
    31,36,35,33,46,51,50,48,  43,40,39,36,58,55,54,51,
    41,38,37,34,32,31,-1,41,56,46,31,-1,31,46,47,32,-1,32,47,49,34,-1,34,49,52,37,-1,37,38,53,52,-1,38,41,56,53,-1,56,46,47,49,52,53,
    37,42,44,43,52,57,59,58
  };
  int conn2[28]={7,12,14,13, 11,8,7,4,2,1, 13,10,9,6, 1,6,5,3, 1,2,4,7,13,6, 0,11,1,3};
  //
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New();
  ret->setMeshDimension(3);
  ret->allocateCells(18);
  //
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+8);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+51);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+59);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+67);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+110);
  //
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+118);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+126);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+169);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+177);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+185);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+228);
  //
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+236);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+244);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+287);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+295);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+303);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+346);
  //
  ret->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(60,3);
  std::copy(coords,coords+180,myCoords->getPointer());
  ret->setCoords(myCoords);
  //
  mesh2D=MEDCouplingUMesh::New();
  mesh2D->setMeshDimension(2);
  mesh2D->allocateCells(6);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn2);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_POLYGON,6,conn2+4);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn2+10);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn2+14);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_POLYGON,6,conn2+18);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn2+24);
  mesh2D->setCoords(myCoords);
  myCoords->decrRef();
  return ret;
}

void MEDCouplingBasicsTest::build3DExtrudedUMesh_2(MEDCouplingUMesh *&meshN, MEDCouplingUMesh *&meshTT)
{
  const double coordsN[270]={
    0, 0, 0, 0.10803000450134277, 0, 0, 0.21606000900268554, 0, 0, 0.28808000564575198, 0, 0, 0.36010002136230468, 0, 0, 0.43212001800537109, 0, 0, 0,
    0.072020001411437995, 0, 0.10803000450134277, 0.072020001411437995, 0, 0.21606000900268554, 0.072020001411437995, 0, 0.28808000564575198, 0.072020001411437995,
    0, 0.36010002136230468, 0.072020001411437995, 0, 0.43212001800537109, 0.072020001411437995, 0, 0, 0.10803000450134277, 0, 0.10803000450134277,
    0.10803000450134277, 0, 0.21606000900268554, 0.10803000450134277, 0, 0.28808000564575198, 0.10803000450134277, 0, 0.36010002136230468, 0.10803000450134277, 0,
    0.43212001800537109, 0.10803000450134277, 0, 0, 0.14404000282287599, 0, 0.10803000450134277, 0.14404000282287599, 0, 0.21606000900268554, 0.14404000282287599, 0,
    0.28808000564575198, 0.14404000282287599, 0, 0.36010002136230468, 0.14404000282287599, 0, 0.43212001800537109, 0.14404000282287599, 0, 0, 0.21606000900268554, 0,
    0.10803000450134277, 0.21606000900268554, 0, 0.21606000900268554, 0.21606000900268554, 0, 0.28808000564575198, 0.21606000900268554, 0, 0.36010002136230468,
    0.21606000900268554, 0, 0.43212001800537109, 0.21606000900268554, 0, 0, 0, 2.1364999389648438, 0.10803000450134277, 0, 2.1364999389648438, 0.21606000900268554,
    0, 2.1364999389648438, 0.28808000564575198, 0, 2.1364999389648438, 0.36010002136230468, 0, 2.1364999389648438, 0.43212001800537109, 0, 2.1364999389648438, 0,
    0.072020001411437995, 2.1364999389648438, 0.10803000450134277, 0.072020001411437995, 2.1364999389648438, 0.21606000900268554, 0.072020001411437995,
    2.1364999389648438, 0.28808000564575198, 0.072020001411437995, 2.1364999389648438, 0.36010002136230468, 0.072020001411437995, 2.1364999389648438,
    0.43212001800537109, 0.072020001411437995, 2.1364999389648438, 0, 0.10803000450134277, 2.1364999389648438, 0.10803000450134277, 0.10803000450134277,
    2.1364999389648438, 0.21606000900268554, 0.10803000450134277, 2.1364999389648438, 0.28808000564575198, 0.10803000450134277, 2.1364999389648438,
    0.36010002136230468, 0.10803000450134277, 2.1364999389648438, 0.43212001800537109, 0.10803000450134277, 2.1364999389648438, 0, 0.14404000282287599,
    2.1364999389648438, 0.10803000450134277, 0.14404000282287599, 2.1364999389648438, 0.21606000900268554, 0.14404000282287599, 2.1364999389648438,
    0.28808000564575198, 0.14404000282287599, 2.1364999389648438, 0.36010002136230468, 0.14404000282287599, 2.1364999389648438, 0.43212001800537109,
    0.14404000282287599, 2.1364999389648438, 0, 0.21606000900268554, 2.1364999389648438, 0.10803000450134277, 0.21606000900268554, 2.1364999389648438,
    0.21606000900268554, 0.21606000900268554, 2.1364999389648438, 0.28808000564575198, 0.21606000900268554, 2.1364999389648438, 0.36010002136230468,
    0.21606000900268554, 2.1364999389648438, 0.43212001800537109, 0.21606000900268554, 2.1364999389648438, 0, 0, 4.2729998779296876, 0.10803000450134277, 0,
    4.2729998779296876, 0.21606000900268554, 0, 4.2729998779296876, 0.28808000564575198, 0, 4.2729998779296876, 0.36010002136230468, 0, 4.2729998779296876,
    0.43212001800537109, 0, 4.2729998779296876, 0, 0.072020001411437995, 4.2729998779296876, 0.10803000450134277, 0.072020001411437995, 4.2729998779296876, 
    0.21606000900268554, 0.072020001411437995, 4.2729998779296876, 0.28808000564575198, 0.072020001411437995, 4.2729998779296876, 0.36010002136230468, 
    0.072020001411437995, 4.2729998779296876, 0.43212001800537109, 0.072020001411437995, 4.2729998779296876, 0, 0.10803000450134277, 4.2729998779296876,
    0.10803000450134277, 0.10803000450134277, 4.2729998779296876, 0.21606000900268554, 0.10803000450134277, 4.2729998779296876, 0.28808000564575198,
    0.10803000450134277, 4.2729998779296876, 0.36010002136230468, 0.10803000450134277, 4.2729998779296876, 0.43212001800537109, 0.10803000450134277, 
    4.2729998779296876, 0, 0.14404000282287599, 4.2729998779296876, 0.10803000450134277, 0.14404000282287599, 4.2729998779296876, 0.21606000900268554,
    0.14404000282287599, 4.2729998779296876, 0.28808000564575198, 0.14404000282287599, 4.2729998779296876, 0.36010002136230468, 0.14404000282287599,
    4.2729998779296876, 0.43212001800537109, 0.14404000282287599, 4.2729998779296876, 0, 0.21606000900268554, 4.2729998779296876, 0.10803000450134277,
    0.21606000900268554, 4.2729998779296876, 0.21606000900268554, 0.21606000900268554, 4.2729998779296876, 0.28808000564575198, 0.21606000900268554,
    4.2729998779296876, 0.36010002136230468, 0.21606000900268554, 4.2729998779296876, 0.43212001800537109, 0.21606000900268554, 4.2729998779296876};
  const int connN[320]={
    0, 1, 7, 6, 30, 31, 37, 36, 1, 2, 8, 7, 31, 32, 38, 37, 2, 3, 9, 8, 32, 33, 39, 38, 3, 4, 10, 9, 33, 34, 40, 39, 4, 5, 11, 10, 34, 35, 41, 40, 6,
    7, 13, 12, 36, 37, 43, 42, 7, 8, 14, 13, 37, 38, 44, 43, 8, 9, 15, 14, 38, 39, 45, 44, 9, 10, 16, 15, 39, 40, 46, 45, 10, 11, 17, 16, 40, 41, 47,
    46, 12, 13, 19, 18, 42, 43, 49, 48, 13, 14, 20, 19, 43, 44, 50, 49, 14, 15, 21, 20, 44, 45, 51, 50, 15, 16, 22, 21, 45, 46, 52, 51, 16, 17, 23,
    22, 46, 47, 53, 52, 18, 19, 25, 24, 48, 49, 55, 54, 19, 20, 26, 25, 49, 50, 56, 55, 20, 21, 27, 26, 50, 51, 57, 56, 21, 22, 28, 27, 51, 52, 58,
    57, 22, 23, 29, 28, 52, 53, 59, 58, 30, 31, 37, 36, 60, 61, 67, 66, 31, 32, 38, 37, 61, 62, 68, 67, 32, 33, 39, 38, 62, 63, 69, 68, 33, 34, 40,
    39, 63, 64, 70, 69, 34, 35, 41, 40, 64, 65, 71, 70, 36, 37, 43, 42, 66, 67, 73, 72, 37, 38, 44, 43, 67, 68, 74, 73, 38, 39, 45, 44, 68, 69, 75,
    74, 39, 40, 46, 45, 69, 70, 76, 75, 40, 41, 47, 46, 70, 71, 77, 76, 42, 43, 49, 48, 72, 73, 79, 78, 43, 44, 50, 49, 73, 74, 80, 79, 44, 45, 51,
    50, 74, 75, 81, 80, 45, 46, 52, 51, 75, 76, 82, 81, 46, 47, 53, 52, 76, 77, 83, 82, 48, 49, 55, 54, 78, 79, 85, 84, 49, 50, 56, 55, 79, 80, 86,
    85, 50, 51, 57, 56, 80, 81, 87, 86, 51, 52, 58, 57, 81, 82, 88, 87, 52, 53, 59, 58, 82, 83, 89, 88};
  meshN=MEDCouplingUMesh::New();
  meshN->setName("meshExtrudedN");
  meshN->setMeshDimension(3);
  meshN->allocateCells(40);
  for(int i=0;i<40;i++)
    meshN->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,connN+8*i);
  meshN->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(90,3);
  std::copy(coordsN,coordsN+270,myCoords->getPointer());
  meshN->setCoords(myCoords);
  myCoords->decrRef();
  //
  meshTT=MEDCouplingUMesh::New();
  meshTT->setName("meshExtrudedTT");
  meshTT->setMeshDimension(3);
  meshTT->allocateCells(200);
  for(int i=0;i<200;i++)
    meshTT->insertNextCell(INTERP_KERNEL::NORM_POLYHED,connITT[i+1]-connITT[i],connTT+connITT[i]);
  meshTT->finishInsertingCells();
  myCoords=DataArrayDouble::New();
  myCoords->alloc(1720,3);
  std::copy(coordsTT,coordsTT+5160,myCoords->getPointer());
  meshTT->setCoords(myCoords);
  myCoords->decrRef();
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build2DTargetMeshMerged_1()
{
  double targetCoords[26]={
    -0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7,
    0.7,-0.3, 1.7,-0.3, 0.7,0.7, 1.7,0.7
  };
  int targetConn[24]={
    0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4,
    9,12,10,9,11,12
  };
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setName("merge");
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(10);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+4);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+7);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+10);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+14);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+18);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+21);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(13,2);
  std::copy(targetCoords,targetCoords+26,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build2DCurveMesh(double dx, double dy)
{
  // 1d mesh:
  //
  //       *
  //      /
  // *---*
  double targetCoords[3*2]=
    {
      0.+dx,0.+dy, 1.+dx,0.+dy, 2.+dx,1.+dy
    };
  int targetConn[2*2]={1,2, 0,1};

  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New("2Dcurve 1D mesh",1);
  targetMesh->allocateCells(2);
  for(int i=0;i<2;i++)
    targetMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,targetConn+2*i);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(3,2);
  std::copy(targetCoords,targetCoords+3*2,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build1DMesh(double dx)
{
  double targetCoords[4]=
    {
      0.+dx, 1.+dx, 3.+dx, 4.+dx
    };
  int targetConn[2*3]={1,2, 0,1, 2,3};

  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New("1D mesh",1);
  targetMesh->allocateCells(3);
  for(int i=0;i<3;i++)
    targetMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,targetConn+2*i);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(4,1);
  std::copy(targetCoords,targetCoords+4,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

double MEDCouplingBasicsTest::sumAll(const std::vector< std::map<int,double> >& matrix)
{
  double ret=0.;
  for(std::vector< std::map<int,double> >::const_iterator iter=matrix.begin();iter!=matrix.end();iter++)
    for(std::map<int,double>::const_iterator iter2=(*iter).begin();iter2!=(*iter).end();iter2++)
      ret+=(*iter2).second;
  return ret;
}
