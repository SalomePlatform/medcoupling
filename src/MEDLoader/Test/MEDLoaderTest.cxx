//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

#include "MEDLoaderTest.hxx"
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"

using namespace ParaMEDMEM;

void MEDLoaderTest::testMesh1DRW()
{
  MEDCouplingUMesh *mesh=build1DMesh_1();
  mesh->checkCoherency();
  MEDLoader::WriteUMesh("file1.med",mesh,true);
  MEDCouplingUMesh *mesh_rw=MEDLoader::ReadUMeshFromFile("file1.med",mesh->getName(),0);
  CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
  mesh_rw->decrRef();
  mesh->decrRef();
}

void MEDLoaderTest::testMesh2DCurveRW()
{
  MEDCouplingUMesh *mesh=build2DCurveMesh_1();
  mesh->checkCoherency();
  MEDLoader::WriteUMesh("file2.med",mesh,true);
  MEDCouplingUMesh *mesh_rw=MEDLoader::ReadUMeshFromFile("file2.med",mesh->getName(),0);
  CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
  mesh_rw->decrRef();
  mesh->decrRef();
}

void MEDLoaderTest::testMesh2DRW()
{
  MEDCouplingUMesh *mesh=build2DMesh_1();
  mesh->checkCoherency();
  MEDLoader::WriteUMesh("file3.med",mesh,true);
  MEDCouplingUMesh *mesh_rw=MEDLoader::ReadUMeshFromFile("file3.med",mesh->getName(),0);
  CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
  mesh_rw->decrRef();
  mesh->decrRef();
}

void MEDLoaderTest::testMesh3DSurfRW()
{
  MEDCouplingUMesh *mesh=build3DSurfMesh_1();
  mesh->checkCoherency();
  MEDLoader::WriteUMesh("file4.med",mesh,true);
  MEDCouplingUMesh *mesh_rw=MEDLoader::ReadUMeshFromFile("file4.med",mesh->getName(),0);
  CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
  mesh_rw->decrRef();
  mesh->decrRef();
}

void MEDLoaderTest::testMesh3DRW()
{
  MEDCouplingUMesh *mesh=build3DMesh_1();
  mesh->checkCoherency();
  MEDLoader::WriteUMesh("file5.med",mesh,true);
  MEDCouplingUMesh *mesh_rw=MEDLoader::ReadUMeshFromFile("file5.med",mesh->getName(),0);
  CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
  mesh_rw->decrRef();
  mesh->decrRef();
}

MEDCouplingUMesh *MEDLoaderTest::build1DMesh_1()
{
  double coords[6]={ 0.0, 0.3, 0.75, 1.0, 1.4, 1.3 };
  int conn[9]={ 0,1, 1,2, 2,3 , 3,4,5};
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setName("1DMesh_1");
  mesh->setMeshDimension(1);
  mesh->allocateCells(4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+2);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG3,3,conn+6);
  mesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(6,1);
  std::copy(coords,coords+6,myCoords->getPointer());
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  return mesh;
}

MEDCouplingUMesh *MEDLoaderTest::build2DCurveMesh_1()
{
  double coords[12]={ 0.0,0.0, 0.3,0.3, 0.75,0.75, 1.0,1.0, 1.4,1.4, 1.3,1.3 };
  int conn[9]={ 0,1, 1,2, 2,3 , 3,4,5};
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setName("2DCurveMesh_1");
  mesh->setMeshDimension(1);
  mesh->allocateCells(4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+2);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG3,3,conn+6);
  mesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(6,2);
  std::copy(coords,coords+12,myCoords->getPointer());
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  return mesh;
}

MEDCouplingUMesh *MEDLoaderTest::build2DMesh_1()
{
  double targetCoords[24]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  int targetConn[24]={1,4,2, 4,5,2, 6,10,8,9,11,7, 0,3,4,1, 6,7,4,3, 7,8,5,4};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(6);
  targetMesh->setName("2DMesh_1");
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+3);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,targetConn+6);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+12);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+16);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+20);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(12,2);
  std::copy(targetCoords,targetCoords+24,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDLoaderTest::build3DSurfMesh_1()
{
  double targetCoords[36]={-0.3,-0.3,-0.3, 0.2,-0.3,-0.3, 0.7,-0.3,-0.3, -0.3,0.2,-0.3, 0.2,0.2,-0.3, 0.7,0.2,-0.3, -0.3,0.7,-0.3, 0.2,0.7,-0.3, 0.7,0.7,-0.3 };
  int targetConn[24]={1,4,2, 4,5,2, 6,10,8,9,11,7, 0,3,4,1, 6,7,4,3, 7,8,5,4};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(6);
  targetMesh->setName("3DSurfMesh_1");
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+3);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,targetConn+6);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+12);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+16);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+20);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(12,3);
  std::copy(targetCoords,targetCoords+36,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDLoaderTest::build3DMesh_1()
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
  //
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New();
  ret->setName("3DMesh_1");
  ret->setMeshDimension(3);
  ret->allocateCells(18);
  //
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+51);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+59);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+110);
  //
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+118);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+169);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+177);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+228);
  //
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+236);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+287);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+295);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+346);
  //
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+8);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+67);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+126);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+185);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+244);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+303);
  //
  ret->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(60,3);
  std::copy(coords,coords+180,myCoords->getPointer());
  ret->setCoords(myCoords);
  myCoords->decrRef();
  return ret;
}
