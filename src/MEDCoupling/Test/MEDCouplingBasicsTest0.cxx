// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
// Author : Anthony Geay (CEA/DEN)

#include "MEDCouplingBasicsTest.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingMappedExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingMultiFields.hxx"

#include "MEDCouplingBasicsTestData1.hxx"

#include "Interpolation2D.txx"
#include "Interpolation2D3D.txx"
#include "Interpolation2D1D.txx"
#include "MEDCouplingNormalizedUnstructuredMesh.txx"
#include "MEDCouplingNormalizedCartesianMesh.txx"

using namespace MEDCoupling;

typedef std::vector<std::map<int,double> > IntersectionMatrix;

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
    0.,0.,0., 1.,1.,0., 1.,1.25,0., 1.,0.,0., 1.,1.5,0., 2.,0.,0., 2.,1.,0., 1.,2.,0., 0.,2.,0., 3.,1.,0.,
    3.,2.,0., 0.,1.,0., 1.,3.,0., 2.,2.,0., 2.,3.,0.,
    0.,0.,1., 1.,1.,1., 1.,1.25,1., 1.,0.,1., 1.,1.5,1., 2.,0.,1., 2.,1.,1., 1.,2.,1., 0.,2.,1., 3.,1.,1.,
    3.,2.,1., 0.,1.,1., 1.,3.,1., 2.,2.,1., 2.,3.,1.,
    0.,0.,2., 1.,1.,2., 1.,1.25,2., 1.,0.,2., 1.,1.5,2., 2.,0.,2., 2.,1.,2., 1.,2.,2., 0.,2.,2., 3.,1.,2.,
    3.,2.,2., 0.,1.,2., 1.,3.,2., 2.,2.,2., 2.,3.,2.,
    0.,0.,3., 1.,1.,3., 1.,1.25,3., 1.,0.,3., 1.,1.5,3., 2.,0.,3., 2.,1.,3., 1.,2.,3., 0.,2.,3., 3.,1.,3.,
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

void MEDCouplingBasicsTest::build3DExtrudedUMesh_2(MEDCouplingUMesh *&meshN, MEDCouplingUMesh *&meshTT, MEDCouplingUMesh *&meshTF)
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
  //
  meshTF=MEDCouplingUMesh::New();
  meshTF->setName("meshExtrudedTF");
  meshTF->setMeshDimension(3);
  meshTF->allocateCells(340);
  for(int i=0;i<320;i++)
    meshTF->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,connTFH8+8*i);
  for(int i=0;i<20;i++)
    meshTF->insertNextCell(INTERP_KERNEL::NORM_POLYHED,connTFPOLH_I[i+1]-connTFPOLH_I[i],connTFPOLH+connTFPOLH_I[i]);
  meshTF->finishInsertingCells();
  myCoords=DataArrayDouble::New();
  myCoords->alloc(567,3);
  std::copy(coordsTF,coordsTF+1701,myCoords->getPointer());
  meshTF->setCoords(myCoords);
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

MEDCouplingUMesh *MEDCouplingBasicsTest::build1DSourceMesh_2()
{
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New("1DSourceMesh",1);
  ret->allocateCells(4);
  int conn[8]={0,1,2,3,1,2,3,4};
  for(int i=0;i<4;i++)
    ret->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+2*i);
  ret->finishInsertingCells();
  double coords[5]={0.3,0.7,0.9,1.0,1.12};
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(5,1);
  std::copy(coords,coords+5,myCoords->getPointer());
  ret->setCoords(myCoords);
  myCoords->decrRef();
  return ret;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build1DTargetMesh_2()
{
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New("1DTargetMesh",1);
  ret->allocateCells(2);
  int conn[4]={1,2,0,1};
  for(int i=0;i<2;i++)
    ret->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+2*i);
  ret->finishInsertingCells();
  double coords[3]={0.5,0.75,1.2};
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(3,1);
  std::copy(coords,coords+3,myCoords->getPointer());
  ret->setCoords(myCoords);
  myCoords->decrRef();
  return ret;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build2DCurveSourceMesh_2()
{
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New("1DSourceMesh",1);
  ret->allocateCells(4);
  int conn[8]={0,1,2,3,1,2,3,4};
  for(int i=0;i<4;i++)
    ret->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+2*i);
  ret->finishInsertingCells();
  double coords[10]={0.3,0.3,0.7,0.7,0.9,0.9,1.0,1.0,1.12,1.12};
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(5,2);
  std::copy(coords,coords+10,myCoords->getPointer());
  ret->setCoords(myCoords);
  myCoords->decrRef();
  return ret;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build2DCurveTargetMesh_2()
{
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New("1DTargetMesh",1);
  ret->allocateCells(2);
  int conn[4]={1,2,0,1};
  for(int i=0;i<2;i++)
    ret->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+2*i);
  ret->finishInsertingCells();
  double coords[6]={0.5,0.5,0.75,0.75,1.2,1.2};
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(3,2);
  std::copy(coords,coords+6,myCoords->getPointer());
  ret->setCoords(myCoords);
  myCoords->decrRef();
  return ret;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build1DTargetMesh_3()
{
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New("1DMesh_3",1);
  ret->allocateCells(4);
  int conn[10]={0,1,2, 3,4, 6,5,7 ,9,8};
  ret->insertNextCell(INTERP_KERNEL::NORM_SEG3,3,conn);
  ret->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+3);
  ret->insertNextCell(INTERP_KERNEL::NORM_SEG3,3,conn+5);
  ret->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+8);
  ret->finishInsertingCells();
  double coords[10]={0.5,1.,0.8,5.,5.21,0.5,1.1,0.7,5.,5.31};
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(10,1);
  std::copy(coords,coords+10,myCoords->getPointer());
  ret->setCoords(myCoords);
  myCoords->decrRef();
  return ret;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build2DCurveTargetMesh_3()
{
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New("2DCurveMesh_3",1);
  ret->allocateCells(4);
  int conn[10]={0,1,2, 3,4, 6,5,7 ,9,8};
  ret->insertNextCell(INTERP_KERNEL::NORM_SEG3,3,conn);
  ret->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+3);
  ret->insertNextCell(INTERP_KERNEL::NORM_SEG3,3,conn+5);
  ret->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+8);
  ret->finishInsertingCells();
  double coords[20]={0.5,0.5,1.,1.,0.8,0.8,5.,5.,5.21,5.21,0.5,0.5,1.1,1.1,0.7,0.7,5.,5.,5.31,5.31};
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(10,2);
  std::copy(coords,coords+20,myCoords->getPointer());
  ret->setCoords(myCoords);
  myCoords->decrRef();
  return ret;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build2DTargetMesh_3()
{
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New("2DMesh_3",2);
  ret->allocateCells(10);
  int conn[52]={
    0,1,2, 0,1,3,4, 0,1,3,5,4, 0,1,2,6,7,8, 0,1,3,4,6,9,2,10,
    0,2,1, 0,4,3,1, 0,4,5,3,1, 0,2,1,8,7,6, 0,4,3,1,10,2,9,6
  };
  ret->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn);
  ret->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+3);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYGON,5,conn+7);
  ret->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,conn+12);
  ret->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8,conn+18);
  ret->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn+26);
  ret->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+29);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYGON,5,conn+33);
  ret->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,conn+38);
  ret->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8,conn+44);
  ret->finishInsertingCells();
  double coords[22]={0.,0.,1.,0.,0.5,1.,1.,1.,0.,1.,0.5,2.,0.5,0.,0.75,0.5,0.25,0.5,1.,0.5,0.,0.5};
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(11,2);
  std::copy(coords,coords+22,myCoords->getPointer());
  ret->setCoords(myCoords);
  myCoords->decrRef();
  ret->checkConsistencyLight();
  return ret;
}

/*!
 * Same as build2DTargetMesh_1 but with more nodes than needed. To check tryToShareSameCoordsPermute method.
 */
MEDCouplingUMesh *MEDCouplingBasicsTest::build2DTargetMesh_4()
{
  double targetCoords[20]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  int targetConn[18]={0,4,5,1, 1,5,3, 5,6,2, 7,8,5,4, 8,9,6,5};
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
  myCoords->alloc(10,2);
  std::copy(targetCoords,targetCoords+20,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build3DTargetMesh_3()
{
  return 0;
}

MEDCouplingMultiFields *MEDCouplingBasicsTest::buildMultiFields_1()
{
  MEDCoupling::MEDCouplingUMesh *m1=build2DTargetMesh_1();
  m1->setName("m1");
  MEDCoupling::MEDCouplingUMesh *m2=build2DTargetMesh_1();
  m2->setName("m2");
  const double vals0[]={-0.7,-1.,-2.,-3.,-4.};
  const double vals1[]={0.,1.,2.,3.,4.,0.1,0.2,0.3,0.4};
  const double vals1_1[]={170.,171.,172.,173.,174.,170.1,170.2,170.3,170.4};
  const double vals2[]={5.,6.,7.,8.,9.};
  const double vals4[]={15.,16.,17.,18.,19.};
  //
  MEDCoupling::DataArrayDouble *d0=MEDCoupling::DataArrayDouble::New(); d0->alloc(5,1); std::copy(vals0,vals0+5,d0->getPointer());
  MEDCoupling::DataArrayDouble *d1=MEDCoupling::DataArrayDouble::New(); d1->alloc(9,1); std::copy(vals1,vals1+9,d1->getPointer());
  MEDCoupling::DataArrayDouble *d1_1=MEDCoupling::DataArrayDouble::New(); d1_1->alloc(9,1); std::copy(vals1_1,vals1_1+9,d1_1->getPointer());
  MEDCoupling::DataArrayDouble *d2=MEDCoupling::DataArrayDouble::New(); d2->alloc(5,1); std::copy(vals2,vals2+5,d2->getPointer());
  MEDCoupling::DataArrayDouble *d4=MEDCoupling::DataArrayDouble::New(); d4->alloc(5,1); std::copy(vals4,vals4+5,d4->getPointer());
  //
  d0->setName("d0"); d1->setName("d1"); d1_1->setName("d1_1"); d2->setName("d2"); d4->setName("d4");
  d0->setInfoOnComponent(0,"c1");
  d1->setInfoOnComponent(0,"c6");
  d1_1->setInfoOnComponent(0,"c9");
  d2->setInfoOnComponent(0,"c5");
  d4->setInfoOnComponent(0,"c7");
  //
  MEDCoupling::MEDCouplingFieldDouble *f0=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::ONE_TIME);
  f0->setMesh(m1);
  f0->setArray(d0);
  f0->setTime(0.2,5,6);
  f0->setName("f0");
  MEDCoupling::MEDCouplingFieldDouble *f1=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES,MEDCoupling::LINEAR_TIME);
  f1->setMesh(m1);
  std::vector<MEDCoupling::DataArrayDouble *> d1s(2); d1s[0]=d1; d1s[1]=d1_1;
  f1->setArrays(d1s);
  f1->setStartTime(0.7,7,8);
  f1->setEndTime(1.2,9,10);
  f1->setName("f1");
  MEDCoupling::MEDCouplingFieldDouble *f2=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::CONST_ON_TIME_INTERVAL);
  f2->setMesh(m2);
  f2->setArray(d2);
  f2->setTime(1.2,11,12);
  f2->setEndTime(1.5,13,14);
  f2->setName("f2");
  MEDCoupling::MEDCouplingFieldDouble *f3=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::ONE_TIME);
  f3->setMesh(m1);
  f3->setArray(d2);
  f3->setTime(1.7,15,16);
  f3->setName("f3");
  MEDCoupling::MEDCouplingFieldDouble *f4=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::NO_TIME);
  f4->setMesh(m2);
  f4->setArray(d4);
  f4->setName("f4");
  //
  std::vector<MEDCoupling::MEDCouplingFieldDouble *> fs(5);
  fs[0]=f0; fs[1]=f1; fs[2]=f2; fs[3]=f3; fs[4]=f4;
  MEDCoupling::MEDCouplingMultiFields *ret=MEDCoupling::MEDCouplingMultiFields::New(fs);
  //
  m1->decrRef();
  m2->decrRef();
  d0->decrRef();
  d1->decrRef();
  d1_1->decrRef();
  d2->decrRef();
  d4->decrRef();
  f0->decrRef();
  f1->decrRef();
  f2->decrRef();
  f3->decrRef();
  f4->decrRef();
  //
  return ret;
}

std::vector<MEDCouplingFieldDouble *> MEDCouplingBasicsTest::buildMultiFields_2()
{
  MEDCoupling::MEDCouplingUMesh *m1=build2DTargetMesh_1();
  m1->setName("m1");
  MEDCoupling::MEDCouplingUMesh *m2=build2DTargetMesh_1();
  m2->setName("m2");
  const double vals0[]={-0.7,-1.,-2.,-3.,-4.};
  const double vals1[]={0.,1.,2.,3.,4.};
  const double vals1_1[]={170.,171.,172.,173.,174.};
  const double vals2[]={5.,6.,7.,8.,9.};
  const double vals4[]={15.,16.,17.,18.,19.};
  //
  MEDCoupling::DataArrayDouble *d0=MEDCoupling::DataArrayDouble::New(); d0->alloc(5,1); std::copy(vals0,vals0+5,d0->getPointer());
  MEDCoupling::DataArrayDouble *d1=MEDCoupling::DataArrayDouble::New(); d1->alloc(5,1); std::copy(vals1,vals1+5,d1->getPointer());
  MEDCoupling::DataArrayDouble *d1_1=MEDCoupling::DataArrayDouble::New(); d1_1->alloc(5,1); std::copy(vals1_1,vals1_1+5,d1_1->getPointer());
  MEDCoupling::DataArrayDouble *d2=MEDCoupling::DataArrayDouble::New(); d2->alloc(5,1); std::copy(vals2,vals2+5,d2->getPointer());
  MEDCoupling::DataArrayDouble *d4=MEDCoupling::DataArrayDouble::New(); d4->alloc(5,1); std::copy(vals4,vals4+5,d4->getPointer());
  //
  d0->setName("d0"); d1->setName("d1"); d1_1->setName("d1_1"); d2->setName("d2"); d4->setName("d4");
  d0->setInfoOnComponent(0,"c1");
  d1->setInfoOnComponent(0,"c6");
  d1_1->setInfoOnComponent(0,"c9");
  d2->setInfoOnComponent(0,"c5");
  d4->setInfoOnComponent(0,"c7");
  //
  MEDCoupling::MEDCouplingFieldDouble *f0=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::ONE_TIME);
  f0->setMesh(m1);
  f0->setArray(d0);
  f0->setTime(0.2,5,6);
  f0->setName("f0");
  MEDCoupling::MEDCouplingFieldDouble *f1=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::LINEAR_TIME);
  f1->setMesh(m1);
  std::vector<MEDCoupling::DataArrayDouble *> d1s(2); d1s[0]=d1; d1s[1]=d1_1;
  f1->setArrays(d1s);
  f1->setStartTime(0.7,7,8);
  f1->setEndTime(1.2,9,10);
  f1->setName("f1");
  MEDCoupling::MEDCouplingFieldDouble *f2=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::CONST_ON_TIME_INTERVAL);
  f2->setMesh(m2);
  f2->setArray(d2);
  f2->setTime(1.2,11,12);
  f2->setEndTime(1.5,13,14);
  f2->setName("f2");
  MEDCoupling::MEDCouplingFieldDouble *f3=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::ONE_TIME);
  f3->setMesh(m1);
  f3->setArray(d2);
  f3->setTime(1.7,15,16);
  f3->setName("f3");
  MEDCoupling::MEDCouplingFieldDouble *f4=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::NO_TIME);
  f4->setMesh(m2);
  f4->setArray(d4);
  f4->setName("f4");
  //
  std::vector<MEDCoupling::MEDCouplingFieldDouble *> fs(5);
  fs[0]=f0; fs[1]=f1; fs[2]=f2; fs[3]=f3; fs[4]=f4;
  m1->decrRef();
  m2->decrRef();
  d0->decrRef();
  d1->decrRef();
  d1_1->decrRef();
  d2->decrRef();
  d4->decrRef();
  //
  return fs;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build1DMultiTypes_1()
{
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New("Multi1DMesh",1);
  DataArrayDouble *coo=buildCoordsForMultiTypes_1();
  const int conn[5]={0,2, 0,2,1};
  mesh->allocateCells(2);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG3,3,conn+2);
  mesh->finishInsertingCells();
  mesh->setCoords(coo);
  coo->decrRef();
  return mesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build2DMultiTypes_1()
{
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New("Multi2DMesh",2);
  DataArrayDouble *coo=buildCoordsForMultiTypes_1();
  const int conn[21]={3,4,5, 3,4,5,6,7,8, 0,9,10,11, 0,9,10,11,12,13,14,15};
  mesh->allocateCells(4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn);
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,conn+3);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+9);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8,conn+13);
  mesh->finishInsertingCells();
  mesh->setCoords(coo);
  coo->decrRef();
  return mesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build3DMultiTypes_1()
{
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New("Multi3DMesh",3);
  DataArrayDouble *coo=buildCoordsForMultiTypes_1();
  const int conn[81]={0,16,17,18,
                      0,16,17,18,19,20,21,22,23,24,
                      0,11,10,9,25,
                      0,11,10,9,25,15,14,13,12,26,27,28,29,
                      0,30,31,32,33,34,
                      0,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
                      0,9,10,11,44,45,46,47,
                      0,9,10,11,44,45,46,47,12,13,14,15,48,49,50,51,52,53,54,55 };
  mesh->allocateCells(8);
  mesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,conn);
  mesh->insertNextCell(INTERP_KERNEL::NORM_TETRA10,10,conn+4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_PYRA5,5,conn+14);
  mesh->insertNextCell(INTERP_KERNEL::NORM_PYRA13,13,conn+19);
  mesh->insertNextCell(INTERP_KERNEL::NORM_PENTA6,6,conn+32);
  mesh->insertNextCell(INTERP_KERNEL::NORM_PENTA15,15,conn+38);
  mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+53);
  mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA20,20,conn+61);
  mesh->finishInsertingCells();
  mesh->setCoords(coo);
  coo->decrRef();
  return mesh;
}

DataArrayDouble *MEDCouplingBasicsTest::buildCoordsForMultiTypes_1()
{
  DataArrayDouble *coords=DataArrayDouble::New();
  coords->alloc(56,3);
  coords->setInfoOnComponent(0,"X (cm)");
  coords->setInfoOnComponent(1,"Y (cm)");
  coords->setInfoOnComponent(2,"Z (cm)");
  const double data[168]={
    0.0, 0.0, 0.0, //#0
    0.5, 0.5, 0.5, //#1
    1.0, 1.0, 1.0, //#2
    1.0, 1.0, 0.0, //#3
    2.0, 2.5, 0.0, //#4
    6.0, 1.5, 0.0, //#5
    1.0, 2.0, 0.0, //#6
    4.5, 2.5, 0.0, //#7
    4.0, 0.5, 0.0, //#8
    0.0, 4.0, 0.0, //#9
    4.0, 4.0, 0.0, //#10
    4.0, 0.0, 0.0, //#11
    0.0, 2.0, 0.0, //#12
    2.0, 4.0, 0.0, //#13
    4.0, 2.0, 0.0, //#14
    2.0, 0.0, 0.0, //#15
    0.0, 6.0, 0.0, //#16
    3.0, 3.0, 0.0, //#17
    1.3, 3.0, 3.0, //#18
    0.0, 3.0, 0.0, //#19
    1.5, 4.5, 0.0, //#20
    1.5, 1.5, 0.0, //#21
    0.65, 1.5, 1.5, //#22
    0.65, 4.5, 1.5, //#23
    2.15, 3.0, 1.5, //#24
    2.0, 2.0, 2.0, //#25
    3.0, 1.0, 1.0, //#26
    3.0, 3.0, 1.0, //#27
    1.0, 3.0, 1.0, //#28
    1.0, 1.0, 1.0, //#29
    0.0, 3.0, 0.0, //#30
    2.0, 0.0, 0.0, //#31
    0.0, 0.0, 6.0, //#32
    0.0, 3.0, 6.0, //#33
    3.0, 0.0, 6.0, //#34
    0.0, 1.5, 0.0, //#35
    1.5, 1.5, 0.0, //#36
    1.5, 0.0, 0.0, //#37
    0.0, 1.5, 6.0, //#38
    1.5, 1.5, 6.0, //#39
    1.5, 0.0, 6.0, //#40
    0.0, 0.0, 3.0, //#41
    0.0, 3.0, 3.0, //#42
    3.0, 0.0, 3.0, //#43
    0.0, 0.0, 4.0, //#44
    0.0, 4.0, 4.0, //#45
    4.0, 4.0, 4.0, //#46
    4.0, 0.0, 4.0, //#47
    0.0, 2.0, 4.0, //#48
    2.0, 4.0, 4.0, //#49
    4.0, 2.0, 4.0, //#50
    2.0, 0.0, 4.0, //#51
    0.0, 0.0, 2.0, //#52
    0.0, 4.0, 2.0, //#53
    4.0, 4.0, 2.0, //#54
    4.0, 0.0, 2.0  //#55
  };
  std::copy(data,data+168,coords->getPointer());
  return coords;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::buildHexa8Mesh_1()
{
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New("Hexa8Only",3);
  DataArrayDouble *coo=DataArrayDouble::New();
  const double coords[81]={0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 1.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.5, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.5, 1.0, 0.0, 0.5, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 0.5, 0.5, 0.0, 1.0, 0.5, 0.5, 1.0, 0.5, 1.0, 1.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5, 1.0, 1.0, 0.5, 1.0, 0.0, 1.0, 1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0};
  coo->alloc(27,3);
  std::copy(coords,coords+81,coo->getPointer());
  const int conn[64]={3,12,13,4,0,9,10,1,
                      4,13,14,5,1,10,11,2,
                      6,15,16,7,3,12,13,4,
                      7,16,17,8,4,13,14,5,
                      12,21,22,13,9,18,19,10,
                      13,22,23,14,10,19,20,11,
                      15,24,25,16,12,21,22,13,
                      16,25,26,17,13,22,23,14};
  mesh->allocateCells(8);
  for(int i=0;i<8;i++)
    mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+8*i);
  mesh->finishInsertingCells();
  mesh->setCoords(coo);
  coo->decrRef();
  return mesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::buildPointe_1(MEDCouplingUMesh *& m1)
{
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New("Pointe.med",3);
  MEDCouplingUMesh *mesh2=MEDCouplingUMesh::New("Pointe.med",2);
  const double coords[57]={0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 1.0, 0.0, 2.0, 1.0, -2.0, 0.0, 1.0, 0.0, -2.0, 1.0, 1.0, 1.0, 2.0, -1.0, 1.0, 2.0, -1.0, -1.0, 2.0, 1.0, -1.0, 2.0, 1.0, 1.0, 3.0, -1.0, 1.0, 3.0, -1.0, -1.0, 3.0, 1.0, -1.0, 3.0, 1.0, 1.0, 4.0, -1.0, 1.0, 4.0, -1.0, -1.0, 4.0, 1.0, -1.0, 4.0, 0.0, 0.0, 5.0};
  const int conn[74]={0,1,2,5,0,1,3,2,0,1,4,3,0,1,5,4,1,6,3,2,1,7,4,3,1,8,5,4,1,9,2,5,1,6,2,9,1,7,3,6,1,8,4,7,1,9,5,8, 6,7,8,9,1,14,17,16,15,18, 10,11,12,13,6,7,8,9,14,15,16,17,10,11,12,13};
  DataArrayDouble *coo=DataArrayDouble::New();
  coo->alloc(19,3);
  std::copy(coords,coords+57,coo->getPointer());
  mesh->setCoords(coo);
  mesh2->setCoords(coo);
  coo->decrRef();
  mesh->allocateCells(16);
  for(int i=0;i<12;i++)
    mesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,conn+4*i);
  mesh->insertNextCell(INTERP_KERNEL::NORM_PYRA5,5,conn+48);
  mesh->insertNextCell(INTERP_KERNEL::NORM_PYRA5,5,conn+53);
  mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+58);
  mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+66);
  mesh->finishInsertingCells();
  //[1,34,29,23,41,32]
  const int conn2[20]={0,5,1,14,18,17,8,7,4,9,5,2, 12,8,9,13,6,7,8,9};
  mesh2->allocateCells(6);
  for(int i=0;i<4;i++)
    mesh2->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn2+3*i);
  mesh2->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn2+12);
  mesh2->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn2+16);
  mesh2->finishInsertingCells();
  m1=mesh2;
  //
  return mesh;
}

double MEDCouplingBasicsTest::sumAll(const std::vector< std::map<int,double> >& matrix)
{
  double ret=0.;
  for(std::vector< std::map<int,double> >::const_iterator iter=matrix.begin();iter!=matrix.end();iter++)
    for(std::map<int,double>::const_iterator iter2=(*iter).begin();iter2!=(*iter).end();iter2++)
      ret+=(*iter2).second;
  return ret;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build2D1DSourceMesh()
{
  double sourceCoords[18]={-17., 3.,  -17., 8.,   -5., 8.,
                            -5., 3.,   -9., 0.,  -13., 3.,
                            -9., 8.,   -7., 0.,   -7., 8.
  };
  int sourceConn[16]={0,1, 1,2, 2,3, 3,0, 3,4, 4,5, 4,6, 7,8};
  MEDCouplingUMesh *sourceMesh=MEDCouplingUMesh::New();
  sourceMesh->setMeshDimension(1);
  sourceMesh->allocateCells(8);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,sourceConn);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,sourceConn+2);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,sourceConn+4);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,sourceConn+6);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,sourceConn+8);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,sourceConn+10);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,sourceConn+12);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,sourceConn+14);
  sourceMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(9,2);
  std::copy(sourceCoords,sourceCoords+18,myCoords->getPointer());
  sourceMesh->setCoords(myCoords);
  myCoords->decrRef();
  return sourceMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build2D1DTargetMesh()
{
  double targetCoords[10]={-17., 0., -17.,6., -9.,6., -9.,0., -5., 3.};
  int targetConn[7]={0,1,2,3, 2,3,4};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(2);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3 ,3,targetConn + 4);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(5,2);
  std::copy(targetCoords,targetCoords+10,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh* MEDCouplingBasicsTest::build2D1DSegSourceMesh(const double shiftX,
                                                                const double inclinationX)
{
  MEDCouplingUMesh *sourceMesh=MEDCouplingUMesh::New();
  sourceMesh->setMeshDimension(1);

  const int nbY = 4;
  const int nbYP1 = nbY + 1;
  sourceMesh->allocateCells(nbY);

  int sourceConn[2];
  for (int iY = 0; iY < nbY; ++iY)
    {
      sourceConn[0] = iY    ;
      sourceConn[1] = iY + 1;
      sourceMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,sourceConn);
    }
  sourceMesh->finishInsertingCells();

  std::vector<double> sourceCoords;
  for (int iY = 0; iY < nbYP1; ++iY)
    {
      sourceCoords.push_back(iY * inclinationX + shiftX);
      sourceCoords.push_back(iY * 4.);
    }
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(nbYP1,2);
  std::copy(sourceCoords.begin(),sourceCoords.end(),myCoords->getPointer());
  sourceMesh->setCoords(myCoords);
  myCoords->decrRef();

  return sourceMesh;
}

MEDCouplingUMesh* MEDCouplingBasicsTest::build2D1DQuadTargetMesh(const double inclinationX)
{
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);

  const int nbX = 5;
  const int nbY = 4;
  const int nbXP1 = nbX + 1;
  const int nbYP1 = nbY + 1;
  targetMesh->allocateCells(nbX * nbY);

  int targetConn[4];
  for (int iX = 0; iX < nbX; ++iX)
    {
      for (int iY = 0; iY < nbY; ++iY)
        {
          targetConn[0] = iY     +  iX      * nbYP1;
          targetConn[1] = iY + 1 +  iX      * nbYP1;
          targetConn[2] = iY + 1 + (iX + 1) * nbYP1;
          targetConn[3] = iY     + (iX + 1) * nbYP1;
          targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn);
        }
    }
  targetMesh->finishInsertingCells();

  std::vector<double> targetCoords;
  for (int iX = 0; iX < nbXP1; ++iX)
    {
      for (int iY = 0; iY < nbYP1; ++iY)
        {
          targetCoords.push_back(iX * 3. + iY * inclinationX);
          targetCoords.push_back(iY * 4.);
        }
    }
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(nbXP1 * nbYP1, 2);
  std::copy(targetCoords.begin(),targetCoords.end(),myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();

  return targetMesh;
}

MEDCouplingUMesh* MEDCouplingBasicsTest::build2D1DTriTargetMesh(const double inclinationX)
{
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);

  const int nbX = 5;
  const int nbY = 4;
  const int nbXP1 = nbX + 1;
  const int nbYP1 = nbY + 1;
  targetMesh->allocateCells(nbX * nbY * 2);

  int targetConn[3];
  for (int iX = 0; iX < nbX; ++iX)
    {
      for (int iY = 0; iY < nbY; ++iY)
        {
          targetConn[0] = iY     +  iX      * nbYP1;
          targetConn[1] = iY + 1 +  iX      * nbYP1;
          targetConn[2] = iY + 1 + (iX + 1) * nbYP1;
          targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn);
          targetConn[0] = iY     +  iX      * nbYP1;
          targetConn[1] = iY + 1 + (iX + 1) * nbYP1;
          targetConn[2] = iY     + (iX + 1) * nbYP1;
          targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn);
        }
    }
  targetMesh->finishInsertingCells();

  std::vector<double> targetCoords;
  for (int iX = 0; iX < nbXP1; ++iX)
    {
      for (int iY = 0; iY < nbYP1; ++iY)
        {
          targetCoords.push_back(iX * 3. + iY * inclinationX);
          targetCoords.push_back(iY * 4.);
        }
    }
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(nbXP1 * nbYP1, 2);
  std::copy(targetCoords.begin(),targetCoords.end(),myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();

  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build3D2DSourceMesh()
{
  double sourceCoords[63]={-12., 6., 10., -12.,10.,  6., -16.,10. , 10.,
                           -20., 0.,  0., -12., 0.,  0., -12., 0. , -4., -20.,0.,-4.,
                           -20., 0., 10., -12., 0., 10., -20.,10. , 10.,
                           -25., 5., -5.,   5., 5., -5.,   5., 5. , 25., -25.,5.,25.,
                           -20., 0., 16., -18., 0., 16., -20., 2.5, 16.,
                           -25., 0., -5.,   5., 0., -5.,   5., 0. , 25., -25.,0.,25.
  };
  int sourceConn[25]={0,1,2, 3,4,5,6, 7,8,9, 10,11,12,13, 14,15,16, 3,4,8,7, 17,18,19,20};
  MEDCouplingUMesh *sourceMesh=MEDCouplingUMesh::New();
  sourceMesh->setMeshDimension(2);
  sourceMesh->allocateCells(7);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3 ,3,sourceConn);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,sourceConn+3);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3 ,3,sourceConn+7);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,sourceConn+10);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3 ,3,sourceConn+14);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,sourceConn+17);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,sourceConn+21);
  sourceMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(21,3);
  std::copy(sourceCoords,sourceCoords+63,myCoords->getPointer());
  sourceMesh->setCoords(myCoords);
  myCoords->decrRef();
  return sourceMesh;
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build3D2DTargetMesh()
{
  double targetCoords[45]={-20., 0., 0., -20.,10., 0., -12.,10., 0.,
                           -12., 0., 0., -20., 0.,10., -20.,10.,10.,
                           -12.,10.,10., -12., 0.,10., -20., 0.,18.,
                           -20.,-5.,10., -20.,-5.,-4., -12.,-5.,-4.,
                           -12.,-5.,10., -20., 0.,-4., -12., 0.,-4.
  };
  int targetConn[20]={4,5,7,8, 0,3,2,1,4,7,6,5, 4,13,14,7,9,10,11,12};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(3);
  targetMesh->allocateCells(3);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,targetConn + 4);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,targetConn + 12);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(15,3);
  std::copy(targetCoords,targetCoords+45,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh* MEDCouplingBasicsTest::build3D2DQuadSourceMesh(const double shiftX,
                                                                 const double inclinationX)
{
  MEDCouplingUMesh *sourceMesh=MEDCouplingUMesh::New();
  sourceMesh->setMeshDimension(2);

  const int nbY = 4;
  const int nbZ = 5;
  const int nbYP1 = nbY + 1;
  const int nbZP1 = nbZ + 1;
  sourceMesh->allocateCells(nbY * nbZ);

  int sourceConn[4];
  for (int iY = 0; iY < nbY; ++iY)
    {
      for (int iZ = 0; iZ < nbZ; ++iZ)
        {
          sourceConn[0] = iZ     +  iY      * nbZP1;
          sourceConn[1] = iZ + 1 +  iY      * nbZP1;
          sourceConn[2] = iZ + 1 + (iY + 1) * nbZP1;
          sourceConn[3] = iZ     + (iY + 1) * nbZP1;
          sourceMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,sourceConn);
        }
    }
  sourceMesh->finishInsertingCells();

  std::vector<double> sourceCoords;
  for (int iY = 0; iY < nbYP1; ++iY)
    {
      for (int iZ = 0; iZ < nbZP1; ++iZ)
        {
            sourceCoords.push_back(iY * inclinationX + shiftX);
            sourceCoords.push_back(iY * 4.);
            sourceCoords.push_back(iZ * 3.);
        }

    }
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(nbYP1 * nbZP1,3);
  std::copy(sourceCoords.begin(),sourceCoords.end(),myCoords->getPointer());
  sourceMesh->setCoords(myCoords);
  myCoords->decrRef();

  return sourceMesh;
}

MEDCouplingUMesh* MEDCouplingBasicsTest::build3D2DTriSourceMesh(const double shiftX,
                                                                const double inclinationX)
{
  MEDCouplingUMesh *sourceMesh=MEDCouplingUMesh::New();
  sourceMesh->setMeshDimension(2);

  const int nbY = 4;
  const int nbZ = 5;
  const int nbYP1 = nbY + 1;
  const int nbZP1 = nbZ + 1;
  sourceMesh->allocateCells(nbY * nbZ * 2);

  int sourceConn[3];
  for (int iY = 0; iY < nbY; ++iY)
    {
      for (int iZ = 0; iZ < nbZ; ++iZ)
        {
        sourceConn[0] = iZ     +  iY      * nbZP1;
        sourceConn[1] = iZ + 1 +  iY      * nbZP1;
        sourceConn[2] = iZ + 1 + (iY + 1) * nbZP1;
        sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,sourceConn);
        sourceConn[0] = iZ     +  iY      * nbZP1;
        sourceConn[1] = iZ     + (iY + 1) * nbZP1;
        sourceConn[2] = iZ + 1 + (iY + 1) * nbZP1;
        sourceMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,sourceConn);
        }
    }
  sourceMesh->finishInsertingCells();

  std::vector<double> sourceCoords;
  for (int iY = 0; iY < nbYP1; ++iY)
    {
      for (int iZ = 0; iZ < nbZP1; ++iZ)
        {
            sourceCoords.push_back(iY * inclinationX + shiftX);
            sourceCoords.push_back(iY * 4.);
            sourceCoords.push_back(iZ * 3.);
        }

    }
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(nbYP1 * nbZP1,3);
  std::copy(sourceCoords.begin(),sourceCoords.end(),myCoords->getPointer());
  sourceMesh->setCoords(myCoords);
  myCoords->decrRef();

  return sourceMesh;
}

MEDCouplingUMesh* MEDCouplingBasicsTest::build3D2DHexaTargetMesh(const double inclinationX)
{
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(3);

  const int nbX = 5;
  const int nbY = 4;
  const int nbZ = 5;
  const int nbXP1 = nbX + 1;
  const int nbYP1 = nbY + 1;
  const int nbZP1 = nbZ + 1;
  targetMesh->allocateCells(nbX * nbY * nbZ);

  int targetConn[8];
  for (int iX = 0; iX < nbX; ++iX)
    {
      for (int iY = 0; iY < nbY; ++iY)
        {
          for (int iZ = 0; iZ < nbZ; ++iZ)
            {
              targetConn[0] = iZ     + ( iY      +  iX      * nbYP1) * nbZP1;
              targetConn[1] = iZ + 1 + ( iY      +  iX      * nbYP1) * nbZP1;
              targetConn[2] = iZ + 1 + ((iY + 1) +  iX      * nbYP1) * nbZP1;
              targetConn[3] = iZ     + ((iY + 1) +  iX      * nbYP1) * nbZP1;
              targetConn[4] = iZ     + ( iY      + (iX + 1) * nbYP1) * nbZP1;
              targetConn[5] = iZ + 1 + ( iY      + (iX + 1) * nbYP1) * nbZP1;
              targetConn[6] = iZ + 1 + ((iY + 1) + (iX + 1) * nbYP1) * nbZP1;
              targetConn[7] = iZ     + ((iY + 1) + (iX + 1) * nbYP1) * nbZP1;
              targetMesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,targetConn);
            }
        }
    }
  targetMesh->finishInsertingCells();

  std::vector<double> targetCoords;
  for (int iX = 0; iX < nbXP1; ++iX)
    {
      for (int iY = 0; iY < nbYP1; ++iY)
        {
          for (int iZ = 0; iZ < nbZP1; ++iZ)
            {
                targetCoords.push_back(iX * 3. + iY * inclinationX);
                targetCoords.push_back(iY * 4.);
                targetCoords.push_back(iZ * 3.);
            }
        }
    }
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(nbXP1 * nbYP1 * nbZP1, 3);
  std::copy(targetCoords.begin(),targetCoords.end(),myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();

  return targetMesh;
}

MEDCouplingUMesh* MEDCouplingBasicsTest::build3D2DTetraTargetMesh(const double inclinationX)
{
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(3);

  const int nbX = 5;
  const int nbY = 4;
  const int nbZ = 5;
  const int nbXP1 = nbX + 1;
  const int nbYP1 = nbY + 1;
  const int nbZP1 = nbZ + 1;
  targetMesh->allocateCells(nbX * nbY * nbZ * 5);

  int targetConn[4];
  for (int iX = 0; iX < nbX; ++iX)
    {
      for (int iY = 0; iY < nbY; ++iY)
        {
          for (int iZ = 0; iZ < nbZ; ++iZ)
            {
              targetConn[0] = iZ     + ( iY      +  iX      * nbYP1) * nbZP1;
              targetConn[1] = iZ + 1 + ( iY      +  iX      * nbYP1) * nbZP1;
              targetConn[2] = iZ + 1 + ( iY      + (iX + 1) * nbYP1) * nbZP1;
              targetConn[3] = iZ + 1 + ((iY + 1) +  iX      * nbYP1) * nbZP1;
              targetMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,targetConn);
              targetConn[0] = iZ     + ( iY      +  iX      * nbYP1) * nbZP1;
              targetConn[1] = iZ     + ( iY      + (iX + 1) * nbYP1) * nbZP1;
              targetConn[2] = iZ + 1 + ( iY      + (iX + 1) * nbYP1) * nbZP1;
              targetConn[3] = iZ     + ((iY + 1) + (iX + 1) * nbYP1) * nbZP1;
              targetMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,targetConn);
              targetConn[0] = iZ     + ( iY      +  iX      * nbYP1) * nbZP1;
              targetConn[1] = iZ     + ((iY + 1) +  iX      * nbYP1) * nbZP1;
              targetConn[2] = iZ     + ((iY + 1) + (iX + 1) * nbYP1) * nbZP1;
              targetConn[3] = iZ + 1 + ((iY + 1) +  iX      * nbYP1) * nbZP1;
              targetMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,targetConn);
              targetConn[0] = iZ + 1 + ( iY      + (iX + 1) * nbYP1) * nbZP1;
              targetConn[1] = iZ + 1 + ((iY + 1) + (iX + 1) * nbYP1) * nbZP1;
              targetConn[2] = iZ     + ((iY + 1) + (iX + 1) * nbYP1) * nbZP1;
              targetConn[3] = iZ + 1 + ((iY + 1) +  iX      * nbYP1) * nbZP1;
              targetMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,targetConn);
              targetConn[0] = iZ     + ( iY      +  iX      * nbYP1) * nbZP1;
              targetConn[1] = iZ + 1 + ((iY + 1) +  iX      * nbYP1) * nbZP1;
              targetConn[2] = iZ + 1 + ( iY      + (iX + 1) * nbYP1) * nbZP1;
              targetConn[3] = iZ     + ((iY + 1) + (iX + 1) * nbYP1) * nbZP1;
              targetMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,targetConn);
            }
        }
    }
  targetMesh->finishInsertingCells();

  std::vector<double> targetCoords;
  for (int iX = 0; iX < nbXP1; ++iX)
    {
      for (int iY = 0; iY < nbYP1; ++iY)
        {
          for (int iZ = 0; iZ < nbZP1; ++iZ)
            {
                targetCoords.push_back(iX * 3. + iY * inclinationX);
                targetCoords.push_back(iY * 4.);
                targetCoords.push_back(iZ * 3.);
            }
        }
    }
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(nbXP1 * nbYP1 * nbZP1, 3);
  std::copy(targetCoords.begin(),targetCoords.end(),myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();

  return targetMesh;
}

int MEDCouplingBasicsTest::countNonZero(const std::vector< std::map<int,double> >& matrix)
{
  int ret=0.;
  for(std::vector< std::map<int,double> >::const_iterator iter=matrix.begin();iter!=matrix.end();iter++)
    for(std::map<int,double>::const_iterator iter2=(*iter).begin();iter2!=(*iter).end();iter2++)
      if (!INTERP_KERNEL::epsilonEqual((*iter2).second, 0.)) ret +=1;
  return ret;
}

void MEDCouplingBasicsTest::test2D1DMeshesIntersection(MEDCouplingUMesh *sourceMesh,
                                                       MEDCouplingUMesh *targetMesh,
                                                       const double correctLength,
                                                       const int correctDuplicateFacesNbr,
                                                       const int correctTotalIntersectFacesNbr)
{
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D1D myInterpolator;
  myInterpolator.setPrecision(1e-12);
  const double prec = 1.0e-5;
  IntersectionMatrix matrix;
  myInterpolator.setIntersectionType(INTERP_KERNEL::Geometric2D);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,matrix,"P0P0");

  std::cout.precision(16);

  const double length = sumAll(matrix);
  LOG(1, "length =  " << length <<"  correctLength = " << correctLength );
  CPPUNIT_ASSERT_DOUBLES_EQUAL(correctLength, length, prec * std::max(correctLength, length));

  INTERP_KERNEL::Interpolation2D3D::DuplicateFacesType duplicateFaces = myInterpolator.retrieveDuplicateFaces();
  int duplicateFacesNbr = duplicateFaces.size();
  LOG(1, "duplicateFacesNbr =  " << duplicateFacesNbr <<"  correctDuplicateFacesNbr = " <<  correctDuplicateFacesNbr);
  CPPUNIT_ASSERT_EQUAL(correctDuplicateFacesNbr, duplicateFacesNbr);

  if (correctTotalIntersectFacesNbr >= 0)
    {
      int totalIntersectFacesNbr = countNonZero(matrix);
      LOG(1, "totalIntersectFacesNbr =  " << totalIntersectFacesNbr <<"  correctTotalIntersectFacesNbr = " << correctTotalIntersectFacesNbr );
      CPPUNIT_ASSERT_EQUAL(correctTotalIntersectFacesNbr, totalIntersectFacesNbr);
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3D2DMeshesIntersection(MEDCouplingUMesh *sourceMesh,
                                                       MEDCouplingUMesh *targetMesh,
                                                       const double correctSurf,
                                                       const int correctDuplicateFacesNbr,
                                                       const int correctTotalIntersectFacesNbr)
{
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D3D myInterpolator;
  myInterpolator.setPrecision(1e-12);
  const double prec = 1.0e-5;
  IntersectionMatrix matrix;
  INTERP_KERNEL::SplittingPolicy sp[] = { INTERP_KERNEL::PLANAR_FACE_5, INTERP_KERNEL::PLANAR_FACE_6, INTERP_KERNEL::GENERAL_24, INTERP_KERNEL::GENERAL_48 };
  for ( size_t i = 0; i < sizeof(sp)/sizeof(sp[0]); ++i )
  {
    myInterpolator.setSplittingPolicy( sp[i] );
    matrix.clear();
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,matrix,"P0P0");

    std::cout.precision(16);

    const double surf = sumAll(matrix);
    LOG(1, "surf =  " << surf <<"  correctSurf = " << correctSurf );
    CPPUNIT_ASSERT_DOUBLES_EQUAL(correctSurf, surf, prec * std::max(correctSurf, surf));

    INTERP_KERNEL::Interpolation2D3D::DuplicateFacesType duplicateFaces = myInterpolator.retrieveDuplicateFaces();
    int duplicateFacesNbr = duplicateFaces.size();
    LOG(1, "duplicateFacesNbr =  " << duplicateFacesNbr <<"  correctDuplicateFacesNbr = " <<  correctDuplicateFacesNbr);
    CPPUNIT_ASSERT_EQUAL(correctDuplicateFacesNbr, duplicateFacesNbr);

    if (correctTotalIntersectFacesNbr >= 0)
      {
        int totalIntersectFacesNbr = countNonZero(matrix);
        LOG(1, "totalIntersectFacesNbr =  " << totalIntersectFacesNbr <<"  correctTotalIntersectFacesNbr = " << correctTotalIntersectFacesNbr );
        CPPUNIT_ASSERT_EQUAL(correctTotalIntersectFacesNbr, totalIntersectFacesNbr);
      }
  }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}
