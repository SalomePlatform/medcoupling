// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#include "MEDCouplingBasicsTest.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingMultiFields.hxx"

void CppExampleFieldDoubleBuildSubPart1()
{
  //! [CppSnippetFieldDoubleBuildSubPart1_1]
  ParaMEDMEM::MEDCouplingUMesh *mesh1=ParaMEDMEM::MEDCouplingBasicsTest::build2DTargetMesh_1();
  ParaMEDMEM::MEDCouplingFieldDouble *f1=ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_CELLS,ParaMEDMEM::ONE_TIME);
  f1->setTime(2.3,5,6);
  f1->setMesh(mesh1);
  ParaMEDMEM::DataArrayDouble *array=ParaMEDMEM::DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),2);
  const double arr1[10]={3.,103.,4.,104.,5.,105.,6.,106.,7.,107.};
  std::copy(arr1,arr1+10,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  //! [CppSnippetFieldDoubleBuildSubPart1_1]
  //! [CppSnippetFieldDoubleBuildSubPart1_2]
  const int part1[3]={2,1,4};
  ParaMEDMEM::MEDCouplingFieldDouble *f2=f1->buildSubPart(part1,part1+3);
  //! [CppSnippetFieldDoubleBuildSubPart1_2]
  f2->zipCoords();
  CPPUNIT_ASSERT_EQUAL(3,f2->getMesh()->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(6,f2->getMesh()->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,f2->getMesh()->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,f2->getMesh()->getMeshDimension());
  ParaMEDMEM::MEDCouplingUMesh *m2C=dynamic_cast<ParaMEDMEM::MEDCouplingUMesh *>(const_cast<ParaMEDMEM::MEDCouplingMesh *>(f2->getMesh()));
  CPPUNIT_ASSERT_EQUAL(13,m2C->getMeshLength());
  const double expected2[12]={0.2, -0.3, 0.7, -0.3, 0.2, 0.2, 0.7, 0.2, 0.2, 0.7, 0.7, 0.7};
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],m2C->getCoords()->getIJ(0,i),1.e-12);
  const double expected3[13]={3,2,3,1,3,0,2,1,4,4,5,3,2};
  CPPUNIT_ASSERT(std::equal(expected3,expected3+13,m2C->getNodalConnectivity()->getConstPointer()));
  const double expected4[4]={0,4,8,13};
  CPPUNIT_ASSERT(std::equal(expected4,expected4+4,m2C->getNodalConnectivityIndex()->getConstPointer()));
  f2->decrRef();
  f1->decrRef();
  //! [CppSnippetFieldDoubleBuildSubPart1_3]
  f1=ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES,ParaMEDMEM::ONE_TIME);
  f1->setTime(2.3,5,6);
  f1->setMesh(mesh1);
  array=ParaMEDMEM::DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfNodes(),2);
  const double arr2[18]={3.,103.,4.,104.,5.,105.,6.,106.,7.,107.,8.,108.,9.,109.,10.,110.,11.,111.};
  std::copy(arr2,arr2+18,array->getPointer());  
  f1->setArray(array);
  array->decrRef();
  //! [CppSnippetFieldDoubleBuildSubPart1_3]
  //! [CppSnippetFieldDoubleBuildSubPart1_4]
  const int part2[2]={1,2};
  f2=f1->buildSubPart(part2,part2+2);
  //! [CppSnippetFieldDoubleBuildSubPart1_4]
  f2->decrRef();
  //idem previous because nodes of cell#4 are not fully present in part3 
  const int part3[2]={1,2};
  ParaMEDMEM::DataArrayInt *arrr=ParaMEDMEM::DataArrayInt::New();
  arrr->alloc(2,1);
  std::copy(part3,part3+2,arrr->getPointer());
  f2=f1->buildSubPart(arrr);
  arrr->decrRef();
  f2->decrRef();
  //
  const int part4[3]={1,2,4};
  f2=f1->buildSubPart(part4,part4+3);
  f2->decrRef();
  //
  f1->decrRef();
  mesh1->decrRef();
  return;
}

void CppSnippetUMeshStdBuild1()
{
  //! [CppSnippetUMeshStdBuild1_1]
  double coords[27]={-0.3,-0.3,0.,   0.2,-0.3,0.,   0.7,-0.3,0.,   -0.3,0.2,0.,   0.2,0.2,0., 
                     0.7,0.2,0.,    -0.3,0.7,0.,    0.2,0.7,0.,     0.7,0.7,0. };
  int nodalConnPerCell[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  //! [CppSnippetUMeshStdBuild1_1]
  //! [CppSnippetUMeshStdBuild1_2]
  ParaMEDMEM::MEDCouplingUMesh *mesh=ParaMEDMEM::MEDCouplingUMesh::New("My2DMesh",2);
  //! [CppSnippetUMeshStdBuild1_2]
  //! [CppSnippetUMeshStdBuild1_3]
  mesh->allocateCells(5);//You can put more than 5 if you want but not less.
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,nodalConnPerCell);
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,nodalConnPerCell+4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,nodalConnPerCell+7);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,nodalConnPerCell+10);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,nodalConnPerCell+14);
  mesh->finishInsertingCells();
  //! [CppSnippetUMeshStdBuild1_3]
  //! [CppSnippetUMeshStdBuild1_4]
  ParaMEDMEM::DataArrayDouble *myCoords=ParaMEDMEM::DataArrayDouble::New();
  myCoords->alloc(9,3);//here myCoords are declared to have 3 components, mesh will deduce that its spaceDim==3. 
  std::copy(coords,coords+27,myCoords->getPointer());
  mesh->setCoords(myCoords);//myCorrds contains 9 tuples, that is to say mesh contains 9 nodes.
  myCoords->decrRef();
  //! [CppSnippetUMeshStdBuild1_4]
  mesh->checkCoherency();
  //! [CppSnippetUMeshStdBuild1_5]
  mesh->decrRef();
  //! [CppSnippetUMeshStdBuild1_5]
}

void CppSnippetUMeshAdvBuild1()
{
  //! [CppSnippetUMeshAdvBuild1_1]
  double coords[27]={-0.3,-0.3,0.,   0.2,-0.3,0.,   0.7,-0.3,0.,   -0.3,0.2,0.,   0.2,0.2,0., 
                     0.7,0.2,0.,    -0.3,0.7,0.,    0.2,0.7,0.,     0.7,0.7,0. };
  int nodalConnPerCell[23]={4,0,3,4,1, 3,1,4,2, 3,4,5,2, 4,6,7,4,3, 4,7,8,5,4};
  int nodalConnPerCellIndex[6]={0,5,9,13,18,23};
  //! [CppSnippetUMeshAdvBuild1_1]
  //! [CppSnippetUMeshAdvBuild1_2]
  ParaMEDMEM::MEDCouplingUMesh *mesh=ParaMEDMEM::MEDCouplingUMesh::New("My2DMesh",2);
  //! [CppSnippetUMeshAdvBuild1_2]
  //! [CppSnippetUMeshAdvBuild1_3]
  ParaMEDMEM::DataArrayInt *nodalConn=ParaMEDMEM::DataArrayInt::New();
  nodalConn->alloc(23,1);
  std::copy(nodalConnPerCell,nodalConnPerCell+23,nodalConn->getPointer());
  ParaMEDMEM::DataArrayInt *nodalConnI=ParaMEDMEM::DataArrayInt::New();
  nodalConnI->alloc(6,1);
  std::copy(nodalConnPerCellIndex,nodalConnPerCellIndex+6,nodalConnI->getPointer());
  mesh->setConnectivity(nodalConn,nodalConnI,true);
  nodalConn->decrRef();// nodalConn DataArrayInt instance is owned by mesh after call to setConnectivity method. No more need here -> decrRef()
  nodalConnI->decrRef();// nodalConnI DataArrayInt instance is owned by mesh after call to setConnectivity method. No more need here -> decrRef()
  //! [CppSnippetUMeshAdvBuild1_3]
  //! [CppSnippetUMeshAdvBuild1_4]
  ParaMEDMEM::DataArrayDouble *myCoords=ParaMEDMEM::DataArrayDouble::New();
  myCoords->alloc(9,3);//here myCoords are declared to have 3 components, mesh will deduce that its spaceDim==3. 
  std::copy(coords,coords+27,myCoords->getPointer());
  mesh->setCoords(myCoords);//myCorrds contains 9 tuples, that is to say mesh contains 9 nodes.
  myCoords->decrRef();
  //! [CppSnippetUMeshAdvBuild1_4]
  mesh->checkCoherency();
  //! [CppSnippetUMeshAdvBuild1_5]
  mesh->decrRef();
  //! [CppSnippetUMeshAdvBuild1_5]
}

int main(int argc, char *argv[])
{
  CppExampleFieldDoubleBuildSubPart1();
  CppSnippetUMeshStdBuild1();
  CppSnippetUMeshAdvBuild1();
  return 0;
}
