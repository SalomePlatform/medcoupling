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
// Author : Anthony Geay (CEA/DEN)

#include "MEDCouplingBasicsTest.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingMultiFields.hxx"


void CppExample_MEDCouplingPointSet_scale()
{
  using namespace ParaMEDMEM;
  //! [CppSnippet_MEDCouplingPointSet_scale_1]
  double coords[4*2]={0.0,0.0, 1.0,0.0, 1.0,1.0, 0.0,1.0}; // 2D coordinates of 4 nodes
  DataArrayDouble *coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 4,2);
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  DataArrayDouble *initCoords = coordsArr->deepCpy();
  //! [CppSnippet_MEDCouplingPointSet_scale_1]
  //! [CppSnippet_MEDCouplingPointSet_scale_2]
  const double center[2] = {0.,0.};
  const double factor = 2.;
  mesh->scale( center, factor );
  //! [CppSnippet_MEDCouplingPointSet_scale_2]
  //! [CppSnippet_MEDCouplingPointSet_scale_3]
  const DataArrayDouble * coordsArr2 = mesh->getCoords();
  CPPUNIT_ASSERT( coordsArr2->isEqualWithoutConsideringStr( *initCoords, 1.0 ));
  CPPUNIT_ASSERT( !coordsArr2->isEqualWithoutConsideringStr( *initCoords, 0.9 ));
  // release data
  mesh->decrRef();
  coordsArr->decrRef();
  initCoords->decrRef();
  //! [CppSnippet_MEDCouplingPointSet_scale_3]
}

void CppExample_MEDCouplingPointSet_translate()
{
  using namespace ParaMEDMEM;
  //! [CppSnippet_MEDCouplingPointSet_translate_1]
  double coords[4*2]={0.0,0.0, 1.0,0.0, 1.0,1.0, 0.0,1.0}; // 2D coordinates of 4 nodes
  DataArrayDouble *coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 4,2);
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  DataArrayDouble *initCoords = coordsArr->deepCpy();
  //! [CppSnippet_MEDCouplingPointSet_translate_1]
  //! [CppSnippet_MEDCouplingPointSet_translate_2]
  double vector[2] = {1.,1.};
  mesh->translate( vector );
  //! [CppSnippet_MEDCouplingPointSet_translate_2]
  //! [CppSnippet_MEDCouplingPointSet_translate_3]
  const DataArrayDouble * coordsArr2 = mesh->getCoords();
  CPPUNIT_ASSERT( coordsArr2->isEqualWithoutConsideringStr( *initCoords, 1.0 ));
  CPPUNIT_ASSERT( !coordsArr2->isEqualWithoutConsideringStr( *initCoords, 0.9 ));
  // release data
  mesh->decrRef();
  coordsArr->decrRef();
  initCoords->decrRef();
  //! [CppSnippet_MEDCouplingPointSet_translate_3]
}

void CppExample_MEDCouplingPointSet_rotate()
{
  using namespace ParaMEDMEM;
  //! [CppSnippet_MEDCouplingPointSet_rotate_1]
  double coords[4*2]={0.0,0.0, 0.1,0.0, 0.1,0.1, 0.0,0.1}; // 2D coordinates of 4 nodes
  DataArrayDouble *coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 4,2);
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingPointSet_rotate_1]
  //! [CppSnippet_MEDCouplingPointSet_rotate_2]
  double center[3] = {0.,0.,0.}; // it suits for 2D as well
  double vector[3] = {0.,0.,1.}; // it is not used in 2D
  mesh->rotate( center, vector, -M_PI/2);
  //! [CppSnippet_MEDCouplingPointSet_rotate_2]
  //! [CppSnippet_MEDCouplingPointSet_rotate_3]
  mesh->changeSpaceDimension(3);
  mesh->rotate( center, vector, +M_PI/2);
  //! [CppSnippet_MEDCouplingPointSet_rotate_3]
  //! [CppSnippet_MEDCouplingPointSet_rotate_4]
  mesh->changeSpaceDimension(2);
  const DataArrayDouble * coordsArr2 = mesh->getCoords();
  CPPUNIT_ASSERT( coordsArr2->isEqualWithoutConsideringStr( *coordsArr, 1e-13 ));
  // release data
  mesh->decrRef();
  coordsArr->decrRef();
  //! [CppSnippet_MEDCouplingPointSet_rotate_4]
}

void CppExample_MEDCouplingPointSet_getBoundingBox()
{
  using namespace ParaMEDMEM;
  //! [CppSnippet_MEDCouplingPointSet_getBoundingBox_1]
  double cc[2*3]={0.0, 0.1, 0.2, // 3D coordinates of 2 nodes
                  2.0, 2.1, 2.2};
  DataArrayDouble *coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(cc, 2,3);
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  coordsArr->decrRef();
  //! [CppSnippet_MEDCouplingPointSet_getBoundingBox_1]
  //! [CppSnippet_MEDCouplingPointSet_getBoundingBox_2]
  double bbox[3][2];
  mesh->getBoundingBox( (double*) bbox );
  mesh->decrRef();

  // check the returned coordinates of extremum points of the bounding box
  for ( int i = 0; i < 2; ++i )   // point id
    for ( int j = 0; j < 3; ++j ) // component
      CPPUNIT_ASSERT_DOUBLES_EQUAL( cc[ i*3 + j ], bbox[j][i], 1e-13);
  //! [CppSnippet_MEDCouplingPointSet_getBoundingBox_2]
}

void CppExample_MEDCouplingPointSet_getNodeIdsNearPoint()
{
  using namespace ParaMEDMEM;
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoint_1]
  // 2D coordinates of 5 nodes
  double coords[5*2]={0.3,-0.301, // #0
                      0.2,-0.3,   // #1
                      0.3,-0.302, // #2
                      1.1,0.0,    // #3
                      0.3,-0.303};// #4
  DataArrayDouble *coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 5,2);
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  coordsArr->decrRef();
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoint_1]
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoint_2]
  double point [2]={0.3, -0.3}; // point close to nodes #0, #2 and #4
  DataArrayInt *ids = mesh->getNodeIdsNearPoint(point, 1e-13);

  // check found ids
  const int expectedIDs[3] = {0,2,4};
  DataArrayInt * okIDs = ids->getIdsEqualList ( expectedIDs, expectedIDs+3 );
  CPPUNIT_ASSERT_EQUAL(3, okIDs->getNumberOfTuples());

  // release data
  mesh->decrRef();
  ids->decrRef();
  okIDs->decrRef();
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoint_2]
}
void CppExample_MEDCouplingPointSet_getNodeIdsNearPoints()
{
  using namespace ParaMEDMEM;
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoints_1]
  // 2D coordinates of 7 nodes
  double coords[7*2]={0.3,-0.301, // #0
                      0.2,-0.3,   // #1
                      0.3,-0.302, // #2
                      1.1,0.0,    // #3
                      1.1,0.0,    // #4
                      1.1,0.002,  // #5
                      0.3,-0.303};// #6
  DataArrayDouble *coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 7,2);
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  coordsArr->decrRef();
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoints_1]
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoints_2]
  const int nbOfPoints = 3;
  double points [nbOfPoints*2]={0.2,-0.301,  // ~ node #1
                                0.0, 0.0,
                                1.1, 0.002}; // ~ nodes #3, #4 and #5
  DataArrayInt *ids, *idsIndex;
  mesh->getNodeIdsNearPoints(points, nbOfPoints, 1e-13,ids,idsIndex);

  // check found ids (i.e. contents of 'ids' array)
  const int expectedIDs[4] = {1, 3, 4, 5};
  DataArrayInt * okIDs = ids->getIdsEqualList ( expectedIDs, expectedIDs+4 );
  CPPUNIT_ASSERT_EQUAL(4, okIDs->getNumberOfTuples());

  // release data
  mesh->decrRef();
  ids->decrRef();
  idsIndex->decrRef();
  okIDs->decrRef();
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoints_2]
}

void CppExample_MEDCouplingPointSet_findCommonNodes()
{
  using namespace ParaMEDMEM;
  //! [CppSnippet_MEDCouplingPointSet_findCommonNodes_1]
  double coords[6*2]={0.3,-0.301, // 0
                      0.2,-0.3,   // 1
                      0.3,-0.302, // 2
                      1.1,0.0,    // 3
                      1.1,0.0,    // 4
                      0.3,-0.303};// 5
  DataArrayDouble *coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 6,2);
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  coordsArr->decrRef();
  //! [CppSnippet_MEDCouplingPointSet_findCommonNodes_1]
  //! [CppSnippet_MEDCouplingPointSet_findCommonNodes_2]
  DataArrayInt *com, *comI;
  mesh->findCommonNodes(1e-13,-1,com,comI);
  CPPUNIT_ASSERT_EQUAL(2, com->getNumberOfTuples());
  mesh->findCommonNodes(0.004,-1,com,comI);
  CPPUNIT_ASSERT_EQUAL(5, com->getNumberOfTuples());
  //! [CppSnippet_MEDCouplingPointSet_findCommonNodes_2]
  mesh->decrRef();
  com->decrRef();
  comI->decrRef();
}

void CppExample_MEDCouplingPointSet_getCoordinatesOfNode()
{
  using namespace ParaMEDMEM;
  //! [CppSnippet_MEDCouplingPointSet_getCoordinatesOfNode_1]
  double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3};
  DataArrayDouble *coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 3,2);
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  coordsArr->decrRef();
  //! [CppSnippet_MEDCouplingPointSet_getCoordinatesOfNode_1]
  //! [CppSnippet_MEDCouplingPointSet_getCoordinatesOfNode_2]
  std::vector<double> coords2;
  mesh->getCoordinatesOfNode(1,coords2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(coords[0],coords2[0],1e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(coords[1],coords2[1],1e-13);
  //! [CppSnippet_MEDCouplingPointSet_getCoordinatesOfNode_2]
  mesh->decrRef();
}

void CppExample_DataArrayInt_buildPermutationArr()
{
  using namespace ParaMEDMEM;
  //! [CppSnippet_DataArrayInt_buildPermutationArr_1]
  DataArrayInt *a=DataArrayInt::New();
  const int vala[5]={4,5,6,7,8};
  a->alloc(5,1);
  std::copy(vala,vala+5,a->getPointer());
  DataArrayInt *b=DataArrayInt::New();
  const int valb[5]={5,4,8,6,7};
  b->alloc(5,1);
  std::copy(valb,valb+5,b->getPointer());
  DataArrayInt *c=a->buildPermutationArr(*b);
  //! [CppSnippet_DataArrayInt_buildPermutationArr_1]
  const int expect1[5]={1,0,4,2,3};
  CPPUNIT_ASSERT_EQUAL(5,c->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,c->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expect1,expect1+5,c->getConstPointer()));
  CPPUNIT_ASSERT(a->isEqualWithoutConsideringStrAndOrder(*b));
  a->decrRef();
  b->decrRef();
  c->decrRef();
}

void CppExample_DataArrayInt_invertArrayO2N2N2O()
{
  using namespace ParaMEDMEM;
  //! [CppSnippet_DataArrayInt_invertArrayO2N2N2O_1]
  const int arr1[6]={2,0,4,1,5,3};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(6,1);
  std::copy(arr1,arr1+6,da->getPointer());
  DataArrayInt *da2=da->invertArrayO2N2N2O(6);
  const int expected1[6]={1,3,0,5,2,4};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],da2->getIJ(i,0));
  //! [CppSnippet_DataArrayInt_invertArrayO2N2N2O_1]
  da->decrRef();
  da2->decrRef();
}

void CppExample_DataArrayInt_invertArrayN2O2O2N()
{
  using namespace ParaMEDMEM;
  //! [CppSnippet_DataArrayInt_invertArrayN2O2O2N_1]
  const int arr1[6]={2,0,4,1,5,3};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(6,1);
  std::copy(arr1,arr1+6,da->getPointer());
  DataArrayInt *da2=da->invertArrayN2O2O2N(6);
  const int expected1[6]={1,3,0,5,2,4};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],da2->getIJ(i,0));
  //! [CppSnippet_DataArrayInt_invertArrayN2O2O2N_1]
  da->decrRef();
  da2->decrRef();
}

void CppExample_DataArrayDouble_getIdsInRange()
{
  using namespace ParaMEDMEM;
  //! [CppSnippet_DataArrayDouble_getIdsInRange_1]
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(10,1);
  da->iota();

  DataArrayInt* da2 = da->getIdsInRange( 2.5, 6 );
  //! [CppSnippet_DataArrayDouble_getIdsInRange_1]
  da->decrRef();
  da2->decrRef();
}

void CppExample_DataArrayDouble_findCommonTuples()
{
  using namespace ParaMEDMEM;
  //! [CppSnippet_DataArrayDouble_findCommonTuples1]
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(6,2);
  const double array2[12]={2.3,2.3, // 0
                           1.2,1.2, // 1
                           1.3,1.3, // 2
                           2.3,2.3, // 3
                           2.301,   // 4
                           2.301,   // 5
                           0.8,0.8};// 6
  std::copy(array2,array2+12,da->getPointer());
  //! [CppSnippet_DataArrayDouble_findCommonTuples1]
  //! [CppSnippet_DataArrayDouble_findCommonTuples2]
  DataArrayInt *c=0,*cI=0;
  da->findCommonTuples(1e-1,-1,c,cI);

  const int expected3[5]={0,3,4,1,2};
  const int expected4[3]={0,3,5};
  CPPUNIT_ASSERT(std::equal(expected3,expected3+5,c->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+3,cI->getConstPointer()));
  c->decrRef();
  cI->decrRef();
  da->decrRef();
  //! [CppSnippet_DataArrayDouble_findCommonTuples2]
}

void CppExample_DataArrayDouble_Meld1()
{
  using namespace ParaMEDMEM;
  //! [CppSnippet_DataArrayDouble_Meld1_1]
  const int sameNbTuples = 7;

  DataArrayDouble *da1=DataArrayDouble::New();
  da1->alloc(sameNbTuples,2);
  da1->fillWithValue(7.);
  da1->setInfoOnComponent(0,"c0da1");
  da1->setInfoOnComponent(1,"c1da1");

  DataArrayDouble *da2=DataArrayDouble::New();
  da2->alloc(sameNbTuples,1);
  da2->iota(0.);
  da2->setInfoOnComponent(0,"c0da2");

  da1->meldWith(da2);
  //! [CppSnippet_DataArrayDouble_Meld1_1]
  //! [CppSnippet_DataArrayDouble_Meld1_2]
  da1->decrRef();
  da2->decrRef();
  //! [CppSnippet_DataArrayDouble_Meld1_2]
}

void CppExample_DataArrayInt_Meld1()
{
  using namespace ParaMEDMEM;
  //! [CppSnippet_DataArrayInt_Meld1_1]
  const int sameNbTuples = 7;

  DataArrayInt *da1=DataArrayInt::New();
  da1->alloc(sameNbTuples,2);
  da1->fillWithValue(7);
  da1->setInfoOnComponent(0,"c0da1");
  da1->setInfoOnComponent(1,"c1da1");

  DataArrayInt *da2=DataArrayInt::New();
  da2->alloc(sameNbTuples,1);
  da2->iota(0);
  da2->setInfoOnComponent(0,"c0da2");

  da1->meldWith(da2);
  //! [CppSnippet_DataArrayInt_Meld1_1]
  //! [CppSnippet_DataArrayInt_Meld1_2]
  da1->decrRef();
  da2->decrRef();
  //! [CppSnippet_DataArrayInt_Meld1_2]
}

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

void CppSnippetCMeshStdBuild1()
{
  //! [CppSnippetCMeshStdBuild1_1]
  double XCoords[9]={-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22};
  double YCoords[7]={0.,0.1,0.37,0.45,0.47,0.49,1.007};
  ParaMEDMEM::DataArrayDouble *arrX=ParaMEDMEM::DataArrayDouble::New();
  arrX->alloc(9,1);
  std::copy(XCoords,XCoords+9,arrX->getPointer());
  arrX->setInfoOnComponent(0,"X [m]");
  ParaMEDMEM::DataArrayDouble *arrY=ParaMEDMEM::DataArrayDouble::New();
  arrY->alloc(7,1);
  std::copy(YCoords,YCoords+7,arrY->getPointer());
  arrY->setInfoOnComponent(0,"Y [m]");
  //! [CppSnippetCMeshStdBuild1_1]
  //! [CppSnippetCMeshStdBuild1_2]
  ParaMEDMEM::MEDCouplingCMesh *mesh=ParaMEDMEM::MEDCouplingCMesh::New("My2D_CMesh");
  mesh->setCoords(arrX,arrY);
  arrX->decrRef();
  arrY->decrRef();
  //! [CppSnippetCMeshStdBuild1_2]
  //! [CppSnippetCMeshStdBuild1_3]
  CPPUNIT_ASSERT_EQUAL(8*6,mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(9*7,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getMeshDimension());
  //! [CppSnippetCMeshStdBuild1_3]
  mesh->decrRef();
  mesh=ParaMEDMEM::MEDCouplingCMesh::New("My2D_CMesh");
  arrX=ParaMEDMEM::DataArrayDouble::New(); arrX->alloc(9,1); std::copy(XCoords,XCoords+9,arrX->getPointer()); arrX->setInfoOnComponent(0,"X [m]");
  arrY=ParaMEDMEM::DataArrayDouble::New(); arrY->alloc(7,1); std::copy(YCoords,YCoords+7,arrY->getPointer()); arrY->setInfoOnComponent(0,"Y [m]");
  //! [CppSnippetCMeshStdBuild1_2bis]
  mesh->setCoordsAt(0,arrX);
  arrX->decrRef();
  mesh->setCoordsAt(1,arrY);
  arrY->decrRef();
  //! [CppSnippetCMeshStdBuild1_2bis]
  CPPUNIT_ASSERT_EQUAL(8*6,mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(9*7,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getMeshDimension());
  //! [CppSnippetCMeshStdBuild1_4]
  mesh->decrRef();
  //! [CppSnippetCMeshStdBuild1_4]
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

void CppSnippetDataArrayBuild1()
{
  //! [CppSnippetDataArrayBuild1_0]
  const int nbOfNodes=12;
  double coords[3*nbOfNodes]={2.,3.,4.,3.,4.,5.,4.,5.,6.,5.,6.,7.,6.,7.,8.,7.,8.,9.,8.,9.,10.,9.,10.,11.,10.,11.,12.,11.,12.,13.,12.,13.,14.,13.,14.,15.};
  //
  ParaMEDMEM::DataArrayDouble *myCoords=0;
  double *tmp=0;
  //! [CppSnippetDataArrayBuild1_0]
  //
  //! [CppSnippetDataArrayBuild1_1]
  myCoords=ParaMEDMEM::DataArrayDouble::New();
  myCoords->useArray(coords,false,ParaMEDMEM::CPP_DEALLOC,nbOfNodes,3);
  //now use myCoords as you need
  //...
  //myCoords is no more useful here : release it
  myCoords->decrRef();
  //! [CppSnippetDataArrayBuild1_1]
  //! [CppSnippetDataArrayBuild1_2]
  myCoords=ParaMEDMEM::DataArrayDouble::New();
  tmp=new double[3*nbOfNodes];
  std::copy(coords,coords+3*nbOfNodes,tmp);
  myCoords->useArray(tmp,true,ParaMEDMEM::CPP_DEALLOC,nbOfNodes,3);
  //now use myCoords as you need
  //...
  //myCoords is no more useful, release it
  myCoords->decrRef();
  //! [CppSnippetDataArrayBuild1_2]
  //! [CppSnippetDataArrayBuild1_3]
  myCoords=ParaMEDMEM::DataArrayDouble::New();
  tmp=(double *)malloc(3*nbOfNodes*sizeof(double));
  std::copy(coords,coords+3*nbOfNodes,tmp);
  myCoords->useArray(tmp,true,ParaMEDMEM::C_DEALLOC,nbOfNodes,3);
  //now use myCoords as you need
  //...
  //myCoords is no more useful here : release it
  myCoords->decrRef();
  //! [CppSnippetDataArrayBuild1_3]
  //! [CppSnippetDataArrayBuild1_4]
  myCoords=ParaMEDMEM::DataArrayDouble::New();
  myCoords->alloc(nbOfNodes,3);
  tmp=myCoords->getPointer();
  std::copy(coords,coords+3*nbOfNodes,tmp);
  myCoords->declareAsNew();//you have modified data pointed by internal pointer notify object
  //now use myCoords as you need
  //...
  //myCoords is no more useful here : release it
  myCoords->decrRef();
  //! [CppSnippetDataArrayBuild1_4]
  myCoords=ParaMEDMEM::DataArrayDouble::New();
  myCoords->alloc(nbOfNodes,3);
  tmp=myCoords->getPointer();
  std::copy(coords,coords+3*nbOfNodes,tmp);
  ParaMEDMEM::DataArrayDouble *myCoordsCpy=0;
  //! [CppSnippetDataArrayBuild1_5]
  myCoordsCpy=myCoords->deepCpy();
  //! [CppSnippetDataArrayBuild1_5]
  //! [CppSnippetDataArrayBuild1_6]
  CPPUNIT_ASSERT(myCoordsCpy->isEqual(*myCoords,1e-12));
  myCoordsCpy->setIJ(0,0,1000.);
  CPPUNIT_ASSERT(!myCoordsCpy->isEqual(*myCoords,1e-12));//myCoordsCpy only has been modified
  //! [CppSnippetDataArrayBuild1_6]
  //! [CppSnippetDataArrayBuild1_7]
  myCoordsCpy->decrRef();
  //! [CppSnippetDataArrayBuild1_7]
  //! [CppSnippetDataArrayBuild1_5bis]
  myCoordsCpy=myCoords->performCpy(true);
  //! [CppSnippetDataArrayBuild1_5bis]
  CPPUNIT_ASSERT(myCoordsCpy->isEqual(*myCoords,1e-12));
  myCoordsCpy->setIJ(0,0,1000.);
  CPPUNIT_ASSERT(!myCoordsCpy->isEqual(*myCoords,1e-12));//myCoordsCpy only has been modified
  myCoordsCpy->decrRef();
  //! [CppSnippetDataArrayBuild1_8]
  myCoordsCpy=myCoords->performCpy(false);
  //! [CppSnippetDataArrayBuild1_8]
  //! [CppSnippetDataArrayBuild1_9]
  CPPUNIT_ASSERT(myCoordsCpy->isEqual(*myCoords,1e-12));
  myCoordsCpy->setIJ(0,0,1000.);
  CPPUNIT_ASSERT(myCoordsCpy->isEqual(*myCoords,1e-12));//myCoords and myCoordsCpy have been modified simultaneously
  //! [CppSnippetDataArrayBuild1_9]
  //! [CppSnippetDataArrayBuild1_10]
  myCoordsCpy->decrRef();
  //! [CppSnippetDataArrayBuild1_10]
  //! [CppSnippetDataArrayBuild1_11]
  myCoordsCpy=ParaMEDMEM::DataArrayDouble::New();
  //! [CppSnippetDataArrayBuild1_11]
  //! [CppSnippetDataArrayBuild1_12]
  myCoordsCpy->cpyFrom(*myCoords);
  //! [CppSnippetDataArrayBuild1_12]
  //! [CppSnippetDataArrayBuild1_13]
  CPPUNIT_ASSERT(myCoordsCpy->isEqual(*myCoords,1e-12));
  myCoordsCpy->setIJ(0,0,2000.);
  CPPUNIT_ASSERT(!myCoordsCpy->isEqual(*myCoords,1e-12));//myCoordsCpy only has been modified
  //! [CppSnippetDataArrayBuild1_13]
  //! [CppSnippetDataArrayBuild1_14]
  myCoordsCpy->decrRef();
  //! [CppSnippetDataArrayBuild1_14]
  myCoords->decrRef();
  //! [CppSnippetDataArrayBuild1_14]
}

void CppSnippetFieldDoubleBuild1()
{
  double XCoords[9]={-0.3,0.07,0.1,0.3,0.45,0.47,0.49,1.,1.22};
  double YCoords[7]={0.07,0.1,0.37,0.45,0.47,0.49,1.007};
  ParaMEDMEM::DataArrayDouble *arrX=ParaMEDMEM::DataArrayDouble::New(); arrX->alloc(9,1); std::copy(XCoords,XCoords+9,arrX->getPointer()); arrX->setInfoOnComponent(0,"X [m]");
  ParaMEDMEM::DataArrayDouble *arrY=ParaMEDMEM::DataArrayDouble::New(); arrY->alloc(7,1); std::copy(YCoords,YCoords+7,arrY->getPointer()); arrY->setInfoOnComponent(0,"Y [m]"); 
  ParaMEDMEM::MEDCouplingCMesh *mesh=ParaMEDMEM::MEDCouplingCMesh::New("My2D_CMesh");
  mesh->setCoords(arrX,arrY); arrX->decrRef(); arrY->decrRef();
  //! [CppSnippetFieldDoubleBuild1_1]
  ParaMEDMEM::MEDCouplingFieldDouble* fieldOnCells=ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_CELLS,ParaMEDMEM::NO_TIME);
  fieldOnCells->setName("MyTensorFieldOnCellNoTime");
  fieldOnCells->setMesh(mesh);
  mesh->decrRef(); // no more need of mesh because mesh has been attached to fieldOnCells
  ParaMEDMEM::DataArrayDouble *array=ParaMEDMEM::DataArrayDouble::New();
  array->alloc(fieldOnCells->getMesh()->getNumberOfCells(),9);//Implicitely fieldOnCells will be a 9 components field.
  array->fillWithValue(7.);
  fieldOnCells->setArray(array);
  array->decrRef();
  // fieldOnCells is now usable
  // ...
  // fieldOnCells is no more useful here : release it
  fieldOnCells->decrRef();
  //! [CppSnippetFieldDoubleBuild1_1]
  arrX=ParaMEDMEM::DataArrayDouble::New(); arrX->alloc(9,1); std::copy(XCoords,XCoords+9,arrX->getPointer()); arrX->setInfoOnComponent(0,"X [m]");
  arrY=ParaMEDMEM::DataArrayDouble::New(); arrY->alloc(7,1); std::copy(YCoords,YCoords+7,arrY->getPointer()); arrY->setInfoOnComponent(0,"Y [m]"); 
  mesh=ParaMEDMEM::MEDCouplingCMesh::New("My2D_CMesh");
  mesh->setCoords(arrX,arrY); arrX->decrRef(); arrY->decrRef();
  //! [CppSnippetFieldDoubleBuild1_2]
  ParaMEDMEM::MEDCouplingFieldDouble *f1=mesh->fillFromAnalytic(ParaMEDMEM::ON_CELLS,1,"x*x+y*y*3+2.*x");//f1 is scalar
  ParaMEDMEM::MEDCouplingFieldDouble *f2=mesh->fillFromAnalytic(ParaMEDMEM::ON_CELLS,1,"cos(x+y/x)");//f2 is scalar too
  ParaMEDMEM::MEDCouplingFieldDouble *f2bis=mesh->fillFromAnalytic(ParaMEDMEM::ON_CELLS,2,"x*x*IVec+3*y*JVec");//f2bis is a vectors field
  ParaMEDMEM::MEDCouplingFieldDouble *f3=(*f1)+(*f2);//f3 scalar
  ParaMEDMEM::MEDCouplingFieldDouble *f4=(*f3)/(*f2);//f4 scalar
  f2bis->applyFunc(1,"sqrt(x*x+y*y)");//f2bis becomes scalar
  ParaMEDMEM::MEDCouplingFieldDouble *f5=(*f2bis)*(*f4);//f5 scalar
  const double pos1[2]={0.48,0.38};
  double res;
  f4->getValueOn(pos1,&res);//f4 is scalar so the returned value is of size 1.
  // ...
  //! [CppSnippetFieldDoubleBuild1_2]
  mesh->decrRef();
  //! [CppSnippetFieldDoubleBuild1_3]
  // f1, f2, f2bis, f3, f4, f5 are no more useful here : release them
  f1->decrRef();
  f2->decrRef();
  f2bis->decrRef();
  f3->decrRef();
  f4->decrRef();
  f5->decrRef();
  //! [CppSnippetFieldDoubleBuild1_3]
}

void CppSnippetFieldDoubleBuild2()
{
  double XCoords[9]={-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22};
  double YCoords[7]={0.,0.1,0.37,0.45,0.47,0.49,1.007};
  ParaMEDMEM::DataArrayDouble *arrX=ParaMEDMEM::DataArrayDouble::New(); arrX->alloc(9,1); std::copy(XCoords,XCoords+9,arrX->getPointer()); arrX->setInfoOnComponent(0,"X [m]");
  ParaMEDMEM::DataArrayDouble *arrY=ParaMEDMEM::DataArrayDouble::New(); arrY->alloc(7,1); std::copy(YCoords,YCoords+7,arrY->getPointer()); arrY->setInfoOnComponent(0,"Y [m]"); 
  ParaMEDMEM::MEDCouplingCMesh *mesh=ParaMEDMEM::MEDCouplingCMesh::New("My2D_CMesh");
  mesh->setCoords(arrX,arrY); arrX->decrRef(); arrY->decrRef();
  //! [CppSnippetFieldDoubleBuild2_1]
  ParaMEDMEM::MEDCouplingFieldDouble* fieldOnNodes=ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES,ParaMEDMEM::NO_TIME);
  fieldOnNodes->setName("MyScalarFieldOnNodeNoTime");
  fieldOnNodes->setMesh(mesh);
  mesh->decrRef(); // no more need of mesh because mesh has been attached to fieldOnNodes
  ParaMEDMEM::DataArrayDouble *array=ParaMEDMEM::DataArrayDouble::New();
  array->alloc(fieldOnNodes->getMesh()->getNumberOfNodes(),1);//Implicitely fieldOnNodes will be a 1 component field.
  array->fillWithValue(8.);
  fieldOnNodes->setArray(array);
  array->decrRef();
  // fieldOnNodes is now usable
  // ...
  // fieldOnNodes is no more useful here : release it
  fieldOnNodes->decrRef();
  //! [CppSnippetFieldDoubleBuild2_1]
}

void CppSnippetFieldDoubleBuild3()
{
  double XCoords[9]={-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22};
  double YCoords[7]={0.,0.1,0.37,0.45,0.47,0.49,1.007};
  ParaMEDMEM::DataArrayDouble *arrX=ParaMEDMEM::DataArrayDouble::New(); arrX->alloc(9,1); std::copy(XCoords,XCoords+9,arrX->getPointer()); arrX->setInfoOnComponent(0,"X [m]");
  ParaMEDMEM::DataArrayDouble *arrY=ParaMEDMEM::DataArrayDouble::New(); arrY->alloc(7,1); std::copy(YCoords,YCoords+7,arrY->getPointer()); arrY->setInfoOnComponent(0,"Y [m]"); 
  ParaMEDMEM::MEDCouplingCMesh *mesh=ParaMEDMEM::MEDCouplingCMesh::New("My2D_CMesh");
  mesh->setCoords(arrX,arrY); arrX->decrRef(); arrY->decrRef();
  //! [CppSnippetFieldDoubleBuild3_1]
  ParaMEDMEM::MEDCouplingFieldDouble* fieldOnCells=ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_CELLS,ParaMEDMEM::ONE_TIME);
  fieldOnCells->setName("MyTensorFieldOnCellNoTime");
  fieldOnCells->setTimeUnit("ms"); // Time unit is ms.
  fieldOnCells->setTime(4.22,2,-1); // Time attached is 4.22 ms, iteration id is 2 and order id (or sub iteration id) is -1
  fieldOnCells->setMesh(mesh);
  mesh->decrRef(); // no more need of mesh because mesh has been attached to fieldOnCells
  ParaMEDMEM::DataArrayDouble *array=ParaMEDMEM::DataArrayDouble::New();
  array->alloc(fieldOnCells->getMesh()->getNumberOfCells(),2);//Implicitely fieldOnCells will be a 2 components field.
  array->fillWithValue(7.);
  fieldOnCells->setArray(array);
  array->decrRef();
  // fieldOnCells is now usable
  // ...
  // fieldOnCells is no more useful here : release it
  fieldOnCells->decrRef();
  //! [CppSnippetFieldDoubleBuild3_1]
}

void CppSnippetFieldDoubleBuild4()
{
  double XCoords[9]={-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22};
  double YCoords[7]={0.,0.1,0.37,0.45,0.47,0.49,1.007};
  ParaMEDMEM::DataArrayDouble *arrX=ParaMEDMEM::DataArrayDouble::New(); arrX->alloc(9,1); std::copy(XCoords,XCoords+9,arrX->getPointer()); arrX->setInfoOnComponent(0,"X [m]");
  ParaMEDMEM::DataArrayDouble *arrY=ParaMEDMEM::DataArrayDouble::New(); arrY->alloc(7,1); std::copy(YCoords,YCoords+7,arrY->getPointer()); arrY->setInfoOnComponent(0,"Y [m]"); 
  ParaMEDMEM::MEDCouplingCMesh *mesh=ParaMEDMEM::MEDCouplingCMesh::New("My2D_CMesh");
  mesh->setCoords(arrX,arrY); arrX->decrRef(); arrY->decrRef();
  //! [CppSnippetFieldDoubleBuild4_1]
  ParaMEDMEM::MEDCouplingFieldDouble* fieldOnNodes=ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES,ParaMEDMEM::CONST_ON_TIME_INTERVAL);
  fieldOnNodes->setName("MyVecFieldOnNodeWithConstTime");
  fieldOnNodes->setTimeUnit("ms"); // Time unit is ms.
  fieldOnNodes->setStartTime(4.22,2,-1);
  fieldOnNodes->setEndTime(6.44,4,-1); // fieldOnNodes is defined in interval [4.22 ms,6.44 ms] 
  fieldOnNodes->setMesh(mesh);
  mesh->decrRef(); // no more need of mesh because mesh has been attached to fieldOnNodes
  ParaMEDMEM::DataArrayDouble *array=ParaMEDMEM::DataArrayDouble::New();
  array->alloc(fieldOnNodes->getMesh()->getNumberOfNodes(),3);//Implicitely fieldOnNodes will be a 3 components field.
  array->fillWithValue(8.);
  fieldOnNodes->setArray(array);
  array->decrRef();
  // fieldOnNodes is now usable
  // ...
  // fieldOnNodes is no more useful here : release it
  fieldOnNodes->decrRef();
  //! [CppSnippetFieldDoubleBuild4_1]
}

int main(int argc, char *argv[])
{
  CppExample_DataArrayInt_buildPermutationArr();
  CppExample_DataArrayInt_invertArrayO2N2N2O();
  CppExample_DataArrayInt_invertArrayN2O2O2N();
  CppExample_DataArrayDouble_getIdsInRange();
  CppExample_DataArrayDouble_findCommonTuples();
  CppExample_DataArrayDouble_Meld1();
  CppExampleFieldDoubleBuildSubPart1();
  CppSnippetUMeshStdBuild1();
  CppSnippetUMeshAdvBuild1();
  CppSnippetDataArrayBuild1();
  CppSnippetCMeshStdBuild1();
  CppSnippetFieldDoubleBuild1();
  CppSnippetFieldDoubleBuild2();
  CppSnippetFieldDoubleBuild3();
  CppSnippetFieldDoubleBuild4();
  return 0;
}
