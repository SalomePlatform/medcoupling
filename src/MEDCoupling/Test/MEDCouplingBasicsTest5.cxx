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

#include "MEDCouplingBasicsTest5.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingMappedExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingGaussLocalization.hxx"
#include "MEDCouplingMultiFields.hxx"
#include "MEDCouplingFieldOverTime.hxx"

#include <cmath>
#include <functional>
#include <iterator>

using namespace MEDCoupling;

void MEDCouplingBasicsTest5::testUMeshTessellate2D1()
{
  double m1Coords[50]={0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1,0.,-1.5,0.5,0.,1.25,0.,0.70710678118654757,0.70710678118654757,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.70710678118654757,0.70710678118654757,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.70710678118654757,-0.70710678118654757,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.70710678118654757,-0.70710678118654757,1.0606601717798214,-1.0606601717798214};
  int m1Conn[56]={0,3,1,13,11,9, 3,4,2,1,14,12,10,11, 5,3,0,15,13,17, 6,4,3,5,16,14,15,18, 5,0,7,17,21,19, 6,5,7,8,18,19,22,20, 0,1,7,9,23,21, 1,2,8,7,10,24,22,23};
  MEDCouplingUMesh *m1=MEDCouplingUMesh::New();
  m1->setMeshDimension(2);
  m1->allocateCells(8);
  m1->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,m1Conn);
  m1->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8,m1Conn+6);
  m1->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,m1Conn+14);
  m1->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8,m1Conn+20);
  m1->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,m1Conn+28);
  m1->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8,m1Conn+34);
  m1->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,m1Conn+42);
  m1->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8,m1Conn+48);
  m1->finishInsertingCells();
  DataArrayDouble *myCoords1=DataArrayDouble::New();
  myCoords1->alloc(25,2);
  std::copy(m1Coords,m1Coords+50,myCoords1->getPointer());
  m1->setCoords(myCoords1);
  myCoords1->decrRef();
  //
  MEDCouplingUMesh *m11=static_cast<MEDCouplingUMesh *>(m1->deepCopy());
  m11->tessellate2D(1.);
  CPPUNIT_ASSERT(m11->getCoords()->isEqual(*m11->getCoords(),1e-12));
  const int expected1[48]={5,0,3,11,1,5,3,4,12,2,1,11,5,5,15,3,0,5,6,16,4,3,15,5,5,5,0,7,19,5,6,5,19,7,8,20,5,0,1,23,7,5,1,2,24,8,7,23};
  const int expected2[9]={0,5,12,17,24,29,36,41,48};
  CPPUNIT_ASSERT_EQUAL(48,(int)m11->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(9,(int)m11->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+48,m11->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+9,m11->getNodalConnectivityIndex()->getConstPointer()));
  m11->decrRef();
  //
  MEDCouplingUMesh *m12=static_cast<MEDCouplingUMesh *>(m1->deepCopy());
  m12->tessellate2D(0.5);
  CPPUNIT_ASSERT_EQUAL(41,m12->getNumberOfNodes());
  const int expected3[60]={5,0,3,25,26,1,5,3,4,27,28,2,1,26,25,5,5,29,30,3,0,5,6,31,32,4,3,30,29,5,5,5,0,7,33,34,5,6,5,34,33,7,8,35,36,5,0,1,37,38,7,5,1,2,39,40,8,7,38,37};
  const int expected4[9]={0,6,15,21,30,36,45,51,60};
  const double expected5[82]={0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1.,0.,-1.5,0.5,0.,1.25,0.,0.7071067811865476,0.7071067811865476,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.7071067811865476,0.7071067811865476,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.7071067811865476,-0.7071067811865476,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.7071067811865476,-0.7071067811865476,1.0606601717798214,-1.0606601717798214,0.479425538604203,0.8775825618903728,0.8414709848078964,0.54030230586814,0.7191383079063044,1.3163738428355591,1.2622064772118446,0.8104534588022099,-0.877582561890373,0.4794255386042027,-0.5403023058681399,0.8414709848078964,-1.3163738428355596,0.7191383079063038,-0.8104534588022098,1.2622064772118446,-0.4794255386042031,-0.8775825618903728,-0.8414709848078965,-0.5403023058681399,-0.7191383079063045,-1.3163738428355591,-1.2622064772118449,-0.8104534588022098,0.8775825618903729,-0.47942553860420295,0.54030230586814,-0.8414709848078964,1.3163738428355594,-0.7191383079063043,0.8104534588022099,-1.2622064772118446};
  for(int i=0;i<82;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected5[i],m12->getCoords()->getIJ(0,i),1e-12);
  CPPUNIT_ASSERT_EQUAL(60,(int)m12->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(9,(int)m12->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+60,m12->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+9,m12->getNodalConnectivityIndex()->getConstPointer()));
  m12->decrRef();
  //
  m1->decrRef();
}

void MEDCouplingBasicsTest5::testUMeshTessellate2DCurve1()
{
  // A quarter of circle:
  double mcoords[6] = {0.4,0.0,   0.0,-0.4,   0.283,-0.283};
  int mconnec[3] = {0,1,2};

  MEDCouplingUMesh *m1 = MEDCouplingUMesh::New();
  m1->setMeshDimension(1);
  m1->allocateCells(1);
  m1->insertNextCell(INTERP_KERNEL::NORM_SEG3, 3, mconnec);

  DataArrayDouble *myCoords = DataArrayDouble::New();
  myCoords->alloc(3,2);
  std::copy(mcoords,mcoords+6,myCoords->getPointer());
  m1->setCoords(myCoords);
  myCoords->decrRef();

  MEDCouplingUMesh *m2 = static_cast<MEDCouplingUMesh *>(m1->deepCopy());
  m2->tessellate2D(0.1);
  CPPUNIT_ASSERT_NO_THROW(m2->checkConsistency(0.0)); // eps param not used
  m1->decrRef();
  m2->decrRef();
}

/*!
 * idem MEDCouplingBasicsTest4::testIntersect2DMeshesTmp3 except that m1 and m2 are permuted on call to MEDCouplingUMesh::Intersect2DMeshes
 */
void MEDCouplingBasicsTest5::testIntersect2DMeshesTmp4()
{
  double m1Coords[50]={0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1,0.,-1.5,0.5,0.,1.25,0.,0.70710678118654757,0.70710678118654757,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.70710678118654757,0.70710678118654757,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.70710678118654757,-0.70710678118654757,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.70710678118654757,-0.70710678118654757,1.0606601717798214,-1.0606601717798214};
  int m1Conn[56]={0,3,1,13,11,9, 3,4,2,1,14,12,10,11, 5,3,0,15,13,17, 6,4,3,5,16,14,15,18, 5,0,7,17,21,19, 6,5,7,8,18,19,22,20, 0,1,7,9,23,21, 1,2,8,7,10,24,22,23};
  MEDCouplingUMesh *m1=MEDCouplingUMesh::New();
  m1->setMeshDimension(2);
  m1->allocateCells(8);
  m1->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,m1Conn);
  m1->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8,m1Conn+6);
  m1->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,m1Conn+14);
  m1->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8,m1Conn+20);
  m1->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,m1Conn+28);
  m1->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8,m1Conn+34);
  m1->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,m1Conn+42);
  m1->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8,m1Conn+48);
  m1->finishInsertingCells();
  DataArrayDouble *myCoords1=DataArrayDouble::New();
  myCoords1->alloc(25,2);
  std::copy(m1Coords,m1Coords+50,myCoords1->getPointer());
  m1->setCoords(myCoords1);
  myCoords1->decrRef();
  //
  double m2Coords[30]={0.,0.,1.1,0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1,-1.1,-1.,0.,-1.,1.1,-1,1.7,-1.};
  int m2Conn[32]={0,3,2,1, 1,2,5,4, 7,6,3,0, 8,9,6,7, 7,0,12,11, 8,7,11,10, 0,1,13,12, 1,4,14,13};
  MEDCouplingUMesh *m2=MEDCouplingUMesh::New();
  m2->setMeshDimension(2);
  m2->allocateCells(8);
  for(int i=0;i<8;i++)
    m2->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,m2Conn+4*i);
  m2->finishInsertingCells();
  DataArrayDouble *myCoords2=DataArrayDouble::New();
  myCoords2->alloc(15,2);
  std::copy(m2Coords,m2Coords+30,myCoords2->getPointer());
  m2->setCoords(myCoords2);
  myCoords2->decrRef();
  //
  DataArrayInt *d1=0,*d2=0;
  MEDCouplingUMesh *m3=MEDCouplingUMesh::Intersect2DMeshes(m2,m1,1e-10,d1,d2);
  m3->unPolyze();
  const int expected1[16]={0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7};
  const int expected2[16]={0,1,1,-1,2,3,3,-1,4,5,5,-1,6,7,7,-1};
 CPPUNIT_ASSERT_EQUAL(16,(int)d1->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(16,(int)d2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(16,(int)m3->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(104,m3->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,m3->getSpaceDimension());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+16,d1->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+16,d2->getConstPointer()));
  const int expected3[136]={6,16,15,18,44,45,46,8,18,2,1,16,47,48,49,50,8,17,1,2,40,51,52,53,54,8,40,5,4,17,55,56,57,58,6,18,15,20,59,60,61,8,20,7,6,18,62,63,64,65,8,41,6,7,21,66,67,68,69,8,21,8,9,41,70,71,72,73,6,20,15,22,74,75,76,8,22,11,7,20,77,78,79,80,8,21,7,11,42,81,82,83,84,8,42,10,8,21,85,86,87,88,6,22,15,16,89,90,91,8,16,1,13,22,92,93,94,95,8,43,13,1,17,96,97,98,99,8,17,4,14,43,100,101,102,103};
  const int expected4[17]={0,7,16,25,34,41,50,59,68,75,84,93,102,109,118,127,136};
  const double expected5[208]={0.,0.,1.1, 0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1.,-1.1,-1.,0.,-1.,1.1,-1.,1.7,-1.,0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1.,0.,-1.5,0.5,0.,1.25,0.,0.7071067811865476,0.7071067811865476,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.7071067811865476,0.7071067811865476,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.7071067811865476,-0.7071067811865476,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.7071067811865476,-0.7071067811865476,1.0606601717798214,-1.0606601717798214,1.1180339887498951,1.,-1.1180339887498951,1.,-1.1180339887498951,-1.,1.1180339887498951,-1.,0.5,0.,0.,0.5,0.7071067811865477,0.7071067811865476,0.55,1.,1.1,0.5,1.05,0.,0.7071067811865477,0.7071067811865475,1.3,0.,1.1,0.5,1.1090169943749475,1.,1.4012585384440737,0.535233134659635,1.4090169943749475,1.,1.7,0.5,1.6,0.,1.4012585384440737,0.535233134659635,0.,0.5,-0.5,0.,-0.7071067811865477,0.7071067811865476,-1.05,0.,-1.1,0.5,-0.55,1.,-0.7071067811865478,0.7071067811865475,-1.1090169943749475,1.,-1.1,0.5,-1.3,0.,-1.4012585384440737,0.5352331346596344,-1.6,0.,-1.7,0.5,-1.4090169943749475,1.,-1.4012585384440737,0.5352331346596344,-0.5,0.,0.,-0.5,-0.7071067811865475,-0.7071067811865477,-0.55,-1.,-1.1,-0.5,-1.05,0.,-0.7071067811865475,-0.7071067811865477,-1.3,0.,-1.1,-0.5,-1.1090169943749475,-1.,-1.4012585384440734,-0.5352331346596354,-1.4090169943749475,-1.,-1.7,-0.5,-1.6,0.,-1.4012585384440732,-0.5352331346596354,0.,-0.5,0.5,0.,0.7071067811865475,-0.7071067811865477,1.05,0.,1.1,-0.5,0.55,-1.,0.7071067811865475,-0.7071067811865477,1.1090169943749475,-1.,1.1,-0.5,1.3,0.,1.4012585384440737,-0.535233134659635,1.6,0.,1.7,-0.5,1.4090169943749475,-1.,1.4012585384440737,-0.535233134659635};
  CPPUNIT_ASSERT_EQUAL(136,(int)m3->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(17,(int)m3->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+136,m3->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+17,m3->getNodalConnectivityIndex()->getConstPointer()));
  for(int i=0;i<208;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected5[i],m3->getCoords()->getIJ(0,i),1e-12);
  d1->decrRef();
  d2->decrRef();
  m3->decrRef();
  //
  m1->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest5::testGetCellIdsCrossingPlane1()
{
  MEDCouplingUMesh *mesh2D=0;
  MEDCouplingUMesh *mesh3D=build3DExtrudedUMesh_1(mesh2D);
  const double vec[3]={-0.07,1.,0.07};
  const double origin[3]={1.524,1.4552,1.74768};
  DataArrayInt *ids1=mesh3D->getCellIdsCrossingPlane(origin,vec,1e-10);
 CPPUNIT_ASSERT_EQUAL(9,(int)ids1->getNumberOfTuples());
  const int expected1[9]={1,3,4,7,9,10,13,15,16};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+9,ids1->getConstPointer()));
  const double vec2[3]={0.,0.,1.};
  DataArrayInt *ids2=mesh3D->getCellIdsCrossingPlane(origin,vec2,1e-10);
  const int expected2[6]={6,7,8,9,10,11};
 CPPUNIT_ASSERT_EQUAL(6,(int)ids2->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected2,expected2+6,ids2->getConstPointer()));
  ids1->decrRef();
  ids2->decrRef();
  mesh3D->decrRef();
  mesh2D->decrRef();
}

void MEDCouplingBasicsTest5::testBuildSlice3D1()
{
  MEDCouplingUMesh *mesh2D=0;
  MEDCouplingUMesh *mesh3D=build3DExtrudedUMesh_1(mesh2D);
  mesh2D->decrRef();
  // First slice in the middle of 3D cells
  const double vec1[3]={-0.07,1.,0.07};
  const double origin1[3]={1.524,1.4552,1.74768};
  DataArrayInt *ids=0;
  MEDCouplingUMesh *slice1=mesh3D->buildSlice3D(origin1,vec1,1e-10,ids);
  const int expected1[9]={1,3,4,7,9,10,13,15,16};
  const int expected2[47]={5,42,41,40,43,44,5,42,46,45,41,5,44,43,40,47,48,5,49,42,44,50,5,49,51,46,42,5,50,44,48,52,5,53,49,50,54,5,53,55,51,49,5,54,50,52,56};
  const int expected3[10]={0,6,11,17,22,27,32,37,42,47};
  const double expected4[171]={1.,1.,0.,1.,1.25,0.,1.,1.5,0.,2.,1.,0.,1.,2.,0.,0.,2.,0.,3.,1.,0.,3.,2.,0.,0.,1.,0.,2.,2.,0.,1.,1.,1.,1.,1.25,1.,1.,1.5,1.,2.,1.,1.,1.,2.,1.,0.,2.,1.,3.,1.,1.,3.,2.,1.,0.,1.,1.,2.,2.,1.,1.,1.,2.,1.,1.25,2.,1.,1.5,2.,2.,1.,2.,1.,2.,2.,0.,2.,2.,3.,1.,2.,3.,2.,2.,0.,1.,2.,2.,2.,2.,1.,1.,3.,1.,1.25,3.,1.,1.5,3.,2.,1.,3.,1.,2.,3.,0.,2.,3.,3.,1.,3.,3.,2.,3.,0.,1.,3.,2.,2.,3.,1.,1.5408576,0.,2.,1.6108576000000001,0.,2.,1.5408576,1.,1.,1.5,0.5836800000000008,1.,1.4708576,1.,3.,1.6808576,0.,3.,1.6108576000000001,1.,0.,1.4708576,0.,0.,1.4008576,1.,2.,1.4708576,2.,1.,1.4008576000000001,2.,3.,1.5408575999999998,2.,0.,1.3308575999999999,2.,2.,1.4008576,3.,1.,1.3308576,3.,3.,1.4708576,3.,0.,1.2608576,3.};
  CPPUNIT_ASSERT_EQUAL(2,slice1->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(3,slice1->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(57,slice1->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(9,(int)slice1->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(9,(int)ids->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(47,(int)slice1->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(10,(int)slice1->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+9,ids->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+47,slice1->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected3,expected3+10,slice1->getNodalConnectivityIndex()->getConstPointer()));
  for(int i=0;i<171;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected4[i],slice1->getCoords()->getIJ(0,i),1e-12);
  ids->decrRef();
  slice1->decrRef();
  // 2nd slice based on already existing nodes of mesh3D.
  const double vec2[3]={0.,3.,1.};
  const double origin2[3]={2.5,1.,3.};
  slice1=mesh3D->buildSlice3D(origin2,vec2,1e-10,ids);
  const int expected5[49]={5,50,10,4,51,5,50,52,7,10,5,51,4,5,53,5,54,50,51,55,56,5,54,57,52,50,5,56,55,51,53,58,5,38,59,56,54,43,5,54,57,46,43,5,38,59,56,58,48};
  const int expected6[10]={0,5,10,15,21,26,32,38,43,49};
  const double expected7[180]={1.,1.,0.,1.,1.25,0.,1.,1.5,0.,2.,1.,0.,1.,2.,0.,0.,2.,0.,3.,1.,0.,3.,2.,0.,0.,1.,0.,1.,3.,0.,2.,2.,0.,2.,3.,0.,1.,1.,1.,1.,1.25,1.,1.,1.5,1.,2.,1.,1.,1.,2.,1.,0.,2.,1.,3.,1.,1.,3.,2.,1.,0.,1.,1.,1.,3.,1.,2.,2.,1.,2.,3.,1.,0.,0.,2.,1.,1.,2.,1.,1.25,2.,1.,0.,2.,1.,1.5,2.,2.,0.,2.,2.,1.,2.,1.,2.,2.,0.,2.,2.,3.,1.,2.,3.,2.,2.,0.,1.,2.,2.,2.,2.,0.,0.,3.,1.,1.,3.,1.,1.25,3.,1.,0.,3.,1.,1.5,3.,2.,0.,3.,2.,1.,3.,1.,2.,3.,0.,2.,3.,3.,1.,3.,3.,2.,3.,0.,1.,3.,2.,2.,3.,2.,1.6666666666666667,1.,1.,1.6666666666666667,1.,3.,1.6666666666666667,1.,0.,1.6666666666666667,1.,2.,1.3333333333333335,2.,1.,1.5,1.5,1.,1.3333333333333333,2.,3.,1.3333333333333335,2.,0.,1.3333333333333335,2.,1.,1.25,2.25};
  CPPUNIT_ASSERT_EQUAL(2,slice1->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(3,slice1->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(60,slice1->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(9,(int)slice1->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(9,(int)ids->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(49,(int)slice1->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(10,(int)slice1->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+9,ids->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected5,expected5+49,slice1->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected6,expected6+10,slice1->getNodalConnectivityIndex()->getConstPointer()));
  for(int i=0;i<180;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected7[i],slice1->getCoords()->getIJ(0,i),1e-12);
  ids->decrRef();
  slice1->decrRef();
  // 3rd slice based on shared face of mesh3D.
  const double vec3[3]={0.,0.,1.};
  const double origin3[3]={2.5,1.,2.};
  slice1=mesh3D->buildSlice3D(origin3,vec3,1e-10,ids);
  const int expected8[12]={6,7,8,9,10,11,12,13,14,15,16,17};
  const int expected9[68]={5,15,26,16,18,5,16,21,28,22,19,17,5,18,20,21,16,5,21,24,25,28,5,26,16,17,19,22,23,5,22,27,29,28,5,15,26,16,18,5,16,21,28,22,19,17,5,18,20,21,16,5,21,24,25,28,5,26,16,17,19,22,23,5,22,27,29,28};
  const int expected10[13]={0,5,12,17,22,29,34,39,46,51,56,63,68};
  const double expected11[135]={0.,0.,1.,1.,1.,1.,1.,1.25, 1.,1.,0.,1.,1.,1.5, 1.,2.,0.,1.,2.,1.,1.,1.,2.,1.,0.,2.,1.,3.,1.,1.,3.,2.,1.,0.,1.,1.,1.,3.,1.,2.,2.,1.,2.,3.,1.,0.,0.,2.,1.,1.,2.,1.,1.25, 2.,1.,0.,2.,1.,1.5, 2.,2.,0.,2.,2.,1.,2.,1.,2.,2.,0.,2.,2.,3.,1.,2.,3.,2.,2.,0.,1.,2.,1.,3.,2.,2.,2.,2.,2.,3.,2.,0.,0.,3.,1.,1.,3.,1.,1.25, 3.,1.,0.,3.,1.,1.5, 3.,2.,0.,3.,2.,1.,3.,1.,2.,3.,0.,2.,3.,3.,1.,3.,3.,2.,3.,0.,1.,3.,1.,3.,3.,2.,2.,3.,2.,3.,3.};
  CPPUNIT_ASSERT_EQUAL(2,slice1->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(3,slice1->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(45,slice1->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(12,(int)slice1->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(12,(int)ids->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(68,(int)slice1->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(13,(int)slice1->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected8,expected8+12,ids->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected9,expected9+68,slice1->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected10,expected10+13,slice1->getNodalConnectivityIndex()->getConstPointer()));
  for(int i=0;i<135;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected11[i],slice1->getCoords()->getIJ(0,i),1e-12);
  ids->decrRef();
  slice1->decrRef();
  //
  mesh3D->decrRef();
}

void MEDCouplingBasicsTest5::testBuildSlice3DSurf1()
{
  MEDCouplingUMesh *mesh2D=0;
  MEDCouplingUMesh *mesh3D=build3DExtrudedUMesh_1(mesh2D);
  mesh2D->decrRef();
  DataArrayInt *a=DataArrayInt::New(),*b=DataArrayInt::New(),*c=DataArrayInt::New(),*d=DataArrayInt::New();
  mesh2D=mesh3D->buildDescendingConnectivity(a,b,c,d);
  a->decrRef(); b->decrRef(); c->decrRef(); d->decrRef();
  mesh3D->decrRef();
  //
  const double vec1[3]={-0.07,1.,0.07};
  const double origin1[3]={1.524,1.4552,1.74768};
  DataArrayInt *ids=0;
  MEDCouplingUMesh *slice1=mesh2D->buildSlice3DSurf(origin1,vec1,1e-10,ids);
  const int expected1[25]={6,8,10,11,13,18,19,21,23,25,26,38,41,43,47,49,52,53,64,67,69,73,75,78,79};
  const int expected2[75]={1,40,41,1,42,41,1,40,43,1,44,43,1,42,44,1,45,41,1,42,46,1,46,45,1,47,40,1,47,48,1,44,48,1,49,42,1,44,50,1,49,50,1,49,51,1,51,46,1,48,52,1,50,52,1,53,49,1,50,54,1,53,54,1,53,55,1,55,51,1,52,56,1,54,56};
  const int expected3[26]={0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75};
  const double expected4[171]={1.,1.,0.,1.,1.25,0.,1.,1.5,0.,2.,1.,0.,1.,2.,0.,0.,2.,0.,3.,1.,0.,3.,2.,0.,0.,1.,0.,2.,2.,0.,1.,1.,1.,1.,1.25,1.,1.,1.5,1.,2.,1.,1.,1.,2.,1.,0.,2.,1.,3.,1.,1.,3.,2.,1.,0.,1.,1.,2.,2.,1.,1.,1.,2.,1.,1.25,2.,1.,1.5,2.,2.,1.,2.,1.,2.,2.,0.,2.,2.,3.,1.,2.,3.,2.,2.,0.,1.,2.,2.,2.,2.,1.,1.,3.,1.,1.25,3.,1.,1.5,3.,2.,1.,3.,1.,2.,3.,0.,2.,3.,3.,1.,3.,3.,2.,3.,0.,1.,3.,2.,2.,3.,1.,1.5408576,0.,2.,1.6108576000000001,0.,2.,1.5408576,1.,1.,1.5,0.5836800000000008,1.,1.4708576,1.,3.,1.6808576,0.,3.,1.6108576000000001,1.,0.,1.4708576,0.,0.,1.4008576,1.,2.,1.4708576,2.,1.,1.4008576000000001,2.,3.,1.5408575999999998,2.,0.,1.3308575999999999,2.,2.,1.4008576,3.,1.,1.3308576,3.,3.,1.4708576,3.,0.,1.2608576,3.};
  CPPUNIT_ASSERT_EQUAL(1,slice1->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(3,slice1->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(57,slice1->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(25,(int)slice1->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(25,(int)ids->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(75,(int)slice1->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(26,(int)slice1->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+25,ids->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+47,slice1->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected3,expected3+26,slice1->getNodalConnectivityIndex()->getConstPointer()));
  for(int i=0;i<171;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected4[i],slice1->getCoords()->getIJ(0,i),1e-12);
  ids->decrRef();
  slice1->decrRef();
  //
  const double vec2[3]={0.,0.,1.};
  const double origin2[3]={2.5,1.,2.};
  slice1=mesh2D->buildSlice3DSurf(origin2,vec2,1e-10,ids);
  const int expected5[68]={32,32,32,32,33,34,35,36,37,38,39,40,41,42,43,43,43,43,43,43,44,44,44,44,45,46,47,47,47,47,48,49,50,51,52,53,53,53,53,53,53,54,54,54,54,55,56,57,59,60,61,62,63,64,65,66,67,68,71,72,74,75,76,77,78,81,82,83};
  const int expected6[204]={1,15,18,1,18,16,1,16,26,1,26,15,1,26,15,1,16,26,1,18,16,1,15,18,1,16,21,1,21,28,1,22,28,1,19,22,1,17,19,1,16,17,1,16,21,1,21,28,1,28,22,1,22,19,1,19,17,1,17,16,1,16,18,1,18,20,1,20,21,1,21,16,1,20,21,1,18,20,1,28,21,1,21,24,1,24,25,1,25,28,1,25,28,1,24,25,1,21,24,1,23,22,1,26,23,1,26,16,1,16,17,1,17,19,1,19,22,1,22,23,1,23,26,1,22,28,1,28,29,1,29,27,1,27,22,1,27,22,1,29,27,1,28,29,1,26,15,1,16,26,1,18,16,1,15,18,1,16,21,1,21,28,1,22,28,1,19,22,1,17,19,1,16,17,1,20,21,1,18,20,1,25,28,1,24,25,1,21,24,1,23,22,1,26,23,1,27,22,1,29,27,1,28,29};
  const int expected7[69]={0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102,105,108,111,114,117,120,123,126,129,132,135,138,141,144,147,150,153,156,159,162,165,168,171,174,177,180,183,186,189,192,195,198,201,204};
  const double expected8[135]={0.,0.,1.,1.,1.,1.,1.,1.25, 1.,1.,0.,1.,1.,1.5, 1.,2.,0.,1.,2.,1.,1.,1.,2.,1.,0.,2.,1.,3.,1.,1.,3.,2.,1.,0.,1.,1.,1.,3.,1.,2.,2.,1.,2.,3.,1.,0.,0.,2.,1.,1.,2.,1.,1.25, 2.,1.,0.,2.,1.,1.5, 2.,2.,0.,2.,2.,1.,2.,1.,2.,2.,0.,2.,2.,3.,1.,2.,3.,2.,2.,0.,1.,2.,1.,3.,2.,2.,2.,2.,2.,3.,2.,0.,0.,3.,1.,1.,3.,1.,1.25, 3.,1.,0.,3.,1.,1.5, 3.,2.,0.,3.,2.,1.,3.,1.,2.,3.,0.,2.,3.,3.,1.,3.,3.,2.,3.,0.,1.,3.,1.,3.,3.,2.,2.,3.,2.,3.,3.};
  CPPUNIT_ASSERT_EQUAL(1,slice1->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(3,slice1->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(45,slice1->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(68,(int)slice1->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(68,(int)ids->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(204,(int)slice1->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(69,(int)slice1->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected5,expected5+68,ids->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected6,expected6+171,slice1->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected7,expected7+69,slice1->getNodalConnectivityIndex()->getConstPointer()));
  for(int i=0;i<135;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected8[i],slice1->getCoords()->getIJ(0,i),1e-12);
  ids->decrRef();
  slice1->decrRef();
  //
  mesh2D->decrRef();
}

void MEDCouplingBasicsTest5::testDataArrayDoubleAdvSetting1()
{
  const double data1[14]={1.,11.,2.,12.,3.,13.,4.,14.,5.,15.,6.,16.,7.,17.};
  const double data2[10]={8.,38.,9.,39.,0.,30.,11.,41.,12.,42.};
  const char *comps[2]={"comp1","comp2"};
  std::vector<std::string> compsCpp(comps,comps+2);
  DataArrayDouble *da=DataArrayDouble::New();
  DataArrayDouble *tmp=0;
  da->setInfoAndChangeNbOfCompo(compsCpp);
  da->setName("da");
  da->alloc(7,2);
  compsCpp.pop_back();
  CPPUNIT_ASSERT_THROW(da->setInfoAndChangeNbOfCompo(compsCpp),INTERP_KERNEL::Exception);
  std::copy(data1,data1+14,da->getPointer());
  //
  std::vector<std::pair<int,int> > p(3);
  p[0].first=0; p[0].second=3; p[1].first=3; p[1].second=5; p[2].first=5; p[2].second=7;
  tmp=dynamic_cast<DataArrayDouble *>(da->selectByTupleRanges(p));
  CPPUNIT_ASSERT(tmp->isEqual(*da,1e-14));
  tmp->decrRef();
  p[0].first=0; p[0].second=2; p[1].first=3; p[1].second=4; p[2].first=5; p[2].second=7;
  tmp=dynamic_cast<DataArrayDouble *>(da->selectByTupleRanges(p));
  const double expected1[10]={1.,11.,2.,12.,4.,14.,6.,16.,7.,17.};
 CPPUNIT_ASSERT_EQUAL(5,(int)tmp->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(2,(int)tmp->getNumberOfComponents());
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],tmp->getIJ(0,i),1e-14);
  tmp->decrRef();
  p[0].first=0; p[0].second=2; p[1].first=0; p[1].second=2; p[2].first=5; p[2].second=6;
  tmp=dynamic_cast<DataArrayDouble *>(da->selectByTupleRanges(p));
  const double expected2[10]={1.,11.,2.,12.,1.,11.,2.,12.,6.,16.};
 CPPUNIT_ASSERT_EQUAL(5,(int)tmp->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(2,(int)tmp->getNumberOfComponents());
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],tmp->getIJ(0,i),1e-14);
  tmp->decrRef();
  p[0].first=0; p[0].second=2; p[1].first=-1; p[1].second=2; p[2].first=5; p[2].second=6;
  CPPUNIT_ASSERT_THROW(da->selectByTupleRanges(p),INTERP_KERNEL::Exception);
  p[0].first=0; p[0].second=2; p[1].first=0; p[1].second=2; p[2].first=5; p[2].second=8;
  CPPUNIT_ASSERT_THROW(da->selectByTupleRanges(p),INTERP_KERNEL::Exception);
  //
  DataArrayDouble *da2=DataArrayDouble::New();
  da2->alloc(5,2);
  std::copy(data2,data2+10,da2->getPointer());
  //
  DataArrayDouble *dac=da->deepCopy();
  dac->setContigPartOfSelectedValuesSlice(1,da2,2,4,1);
  const double expected3[14]={1.,11.,0.,30.,11.,41.,4.,14.,5.,15.,6.,16.,7.,17.};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],dac->getIJ(0,i),1e-14);
  dac->decrRef();
  //
  dac=da->deepCopy();
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValuesSlice(3,da2,0,5,1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValuesSlice(0,da2,4,6,1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValuesSlice(3,da2,5,0,1),INTERP_KERNEL::Exception);
  dac->setContigPartOfSelectedValuesSlice(3,da2,1,5,1);
  const double expected4[14]={1.,11.,2.,12.,3.,13.,9.,39.,0.,30.,11.,41.,12.,42.};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected4[i],dac->getIJ(0,i),1e-14);
  dac->decrRef();
  //
  DataArrayInt *ids=DataArrayInt::New();
  ids->alloc(3,1);
  dac=da->deepCopy();
  ids->setIJ(0,0,2); ids->setIJ(1,0,0); ids->setIJ(2,0,4);
  dac->setContigPartOfSelectedValues(2,da2,ids);
  const double expected5[14]={1.,11.,2.,12.,0.,30.,8.,38.,12.,42.,6.,16.,7.,17.};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected5[i],dac->getIJ(0,i),1e-14);
  dac->decrRef();
  //
  dac=da->deepCopy();
  ids->setIJ(0,0,2); ids->setIJ(1,0,5); ids->setIJ(2,0,4);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues(1,da2,ids),INTERP_KERNEL::Exception);
  ids->setIJ(0,0,2); ids->setIJ(1,0,2); ids->setIJ(2,0,-1);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues(1,da2,ids),INTERP_KERNEL::Exception);
  ids->setIJ(0,0,2); ids->setIJ(1,0,2); ids->setIJ(2,0,1);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues(5,da2,ids),INTERP_KERNEL::Exception);
  dac->decrRef();
  //
  ids->setIJ(0,0,2); ids->setIJ(1,0,2); ids->setIJ(2,0,1);
  dac=da->deepCopy();
  dac->setContigPartOfSelectedValues(4,da2,ids);
  const double expected6[14]={1.,11.,2.,12.,3.,13.,4.,14.,0.,30.,0.,30.,9.,39.};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected6[i],dac->getIJ(0,i),1e-14);
  dac->decrRef();
  ids->decrRef();
  //
  da2->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest5::testDataArrayIntAdvSetting1()
{
  const int data1[14]={1,11,2,12,3,13,4,14,5,15,6,16,7,17};
  const int data2[10]={8,38,9,39,0,30,11,41,12,42};
  const char *comps[2]={"comp1","comp2"};
  std::vector<std::string> compsCpp(comps,comps+2);
  DataArrayInt *da=DataArrayInt::New();
  DataArrayInt *tmp=0;
  da->setInfoAndChangeNbOfCompo(compsCpp);
  da->setName("da");
  da->alloc(7,2);
  compsCpp.pop_back();
  CPPUNIT_ASSERT_THROW(da->setInfoAndChangeNbOfCompo(compsCpp),INTERP_KERNEL::Exception);
  std::copy(data1,data1+14,da->getPointer());
  //
  std::vector<std::pair<int,int> > p(3);
  p[0].first=0; p[0].second=3; p[1].first=3; p[1].second=5; p[2].first=5; p[2].second=7;
  tmp=dynamic_cast<DataArrayInt *>(da->selectByTupleRanges(p));
  CPPUNIT_ASSERT(tmp->isEqual(*da));
  tmp->decrRef();
  p[0].first=0; p[0].second=2; p[1].first=3; p[1].second=4; p[2].first=5; p[2].second=7;
  tmp=dynamic_cast<DataArrayInt *>(da->selectByTupleRanges(p));
  const int expected1[10]={1,11,2,12,4,14,6,16,7,17};
 CPPUNIT_ASSERT_EQUAL(5,(int)tmp->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(2,(int)tmp->getNumberOfComponents());
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],tmp->getIJ(0,i));
  tmp->decrRef();
  p[0].first=0; p[0].second=2; p[1].first=0; p[1].second=2; p[2].first=5; p[2].second=6;
  tmp=dynamic_cast<DataArrayInt *>(da->selectByTupleRanges(p));
  const int expected2[10]={1,11,2,12,1,11,2,12,6,16};
 CPPUNIT_ASSERT_EQUAL(5,(int)tmp->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(2,(int)tmp->getNumberOfComponents());
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],tmp->getIJ(0,i));
  tmp->decrRef();
  p[0].first=0; p[0].second=2; p[1].first=-1; p[1].second=2; p[2].first=5; p[2].second=6;
  CPPUNIT_ASSERT_THROW(da->selectByTupleRanges(p),INTERP_KERNEL::Exception);
  p[0].first=0; p[0].second=2; p[1].first=0; p[1].second=2; p[2].first=5; p[2].second=8;
  CPPUNIT_ASSERT_THROW(da->selectByTupleRanges(p),INTERP_KERNEL::Exception);
  //
  DataArrayInt *da2=DataArrayInt::New();
  da2->alloc(5,2);
  std::copy(data2,data2+10,da2->getPointer());
  //
  DataArrayInt *dac=da->deepCopy();
  dac->setContigPartOfSelectedValuesSlice(1,da2,2,4,1);
  const int expected3[14]={1,11,0,30,11,41,4,14,5,15,6,16,7,17};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_EQUAL(expected3[i],dac->getIJ(0,i));
  dac->decrRef();
  //
  dac=da->deepCopy();
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValuesSlice(3,da2,0,5,1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValuesSlice(0,da2,4,6,1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValuesSlice(3,da2,5,0,1),INTERP_KERNEL::Exception);
  dac->setContigPartOfSelectedValuesSlice(3,da2,1,5,1);
  const int expected4[14]={1,11,2,12,3,13,9,39,0,30,11,41,12,42};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_EQUAL(expected4[i],dac->getIJ(0,i));
  dac->decrRef();
  //
  DataArrayInt *ids=DataArrayInt::New();
  ids->alloc(3,1);
  dac=da->deepCopy();
  ids->setIJ(0,0,2); ids->setIJ(1,0,0); ids->setIJ(2,0,4);
  dac->setContigPartOfSelectedValues(2,da2,ids);
  const int expected5[14]={1,11,2,12,0,30,8,38,12,42,6,16,7,17};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_EQUAL(expected5[i],dac->getIJ(0,i));
  dac->decrRef();
  //
  dac=da->deepCopy();
  ids->setIJ(0,0,2); ids->setIJ(1,0,5); ids->setIJ(2,0,4);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues(1,da2,ids),INTERP_KERNEL::Exception);
  ids->setIJ(0,0,2); ids->setIJ(1,0,2); ids->setIJ(2,0,-1);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues(1,da2,ids),INTERP_KERNEL::Exception);
  ids->setIJ(0,0,2); ids->setIJ(1,0,2); ids->setIJ(2,0,1);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues(5,da2,ids),INTERP_KERNEL::Exception);
  dac->decrRef();
  //
  ids->setIJ(0,0,2); ids->setIJ(1,0,2); ids->setIJ(2,0,1);
  dac=da->deepCopy();
  dac->setContigPartOfSelectedValues(4,da2,ids);
  const int expected6[14]={1,11,2,12,3,13,4,14,0,30,0,30,9,39};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_EQUAL(expected6[i],dac->getIJ(0,i));
  dac->decrRef();
  ids->decrRef();
  //
  da2->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest5::testBuildDescendingConnec2Of3DMesh1()
{
  MEDCouplingUMesh *mesh=build3DSourceMesh_1();
  DataArrayInt *desc=DataArrayInt::New();
  DataArrayInt *descIndx=DataArrayInt::New();
  DataArrayInt *revDesc=DataArrayInt::New();
  DataArrayInt *revDescIndx=DataArrayInt::New();
  //
  MEDCouplingUMesh *mesh2=mesh->buildDescendingConnectivity2(desc,descIndx,revDesc,revDescIndx);
  mesh2->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(2,mesh2->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(30,(int)mesh2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL((std::size_t)31,revDescIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(31,(int)revDescIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)13,descIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(13,(int)descIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)48,desc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(48,(int)desc->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)48,revDesc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(48,(int)revDesc->getNumberOfTuples());
  const int expected1[48]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,-10,15,-5,-13,16,17,-14,18,-4,19,-2,20,21,22,23,24,25,-11,26,-1,-12,-25,-22,27,28,-7,-20,-24,29,-16,-18,30,-8,-28};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+48,desc->getConstPointer()));
  const int expected2[13]={0,4,8,12,16,20,24,28,32,36,40,44,48};
  CPPUNIT_ASSERT(std::equal(expected2,expected2+13,descIndx->getConstPointer()));
  const int expected3[31]={0,2,4,5,7,9,10,12,14,15,17,19,21,23,25,26,28,29,31,32,34,35,37,38,40,42,43,44,46,47,48};
  CPPUNIT_ASSERT(std::equal(expected3,expected3+31,revDescIndx->getConstPointer()));
  const int expected4[48]={0,8,0,6,0,0,5,1,4,1,1,9,1,11,2,2,3,2,7,2,8,3,4,3,5,3,4,10,4,5,11,5,6,10,6,6,9,7,7,10,7,8,8,9,9,11,10,11};
  CPPUNIT_ASSERT(std::equal(expected4,expected4+48,revDesc->getConstPointer()));
  DataArrayInt *conn=mesh2->getNodalConnectivity();
  DataArrayInt *connIndex=mesh2->getNodalConnectivityIndex();
  const int expected5[31]={0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120};
  CPPUNIT_ASSERT(std::equal(expected5,expected5+31,connIndex->getConstPointer()));
  const int expected6[120]={3,8,1,7,3,8,3,1,3,1,3,7,3,7,3,8,3,6,0,8,3,6,2,0,3,0,2,8,3,8,2,6,3,7,4,5,3,7,8,4,3,4,8,5,3,5,8,7,3,6,8,4,3,6,7,8,3,4,7,6,3,8,4,0,3,0,4,6,3,6,3,8,3,7,3,6,3,8,0,1,3,1,0,3,3,3,0,8,3,4,1,5,3,4,8,1,3,1,8,5,3,1,7,5,3,0,2,3,3,3,2,8,3,1,4,0,3,3,2,6};
  CPPUNIT_ASSERT(std::equal(expected6,expected6+120,conn->getConstPointer()));
  //
  desc->decrRef();
  descIndx->decrRef();
  revDesc->decrRef();
  revDescIndx->decrRef();
  mesh2->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest5::testAre2DCellsNotCorrectlyOriented1()
{
  double m1Coords[8]={1.,1.,-1.,-1.,-1.,-1.,1.,-1.};
  int m1Conn[4]={0,3,1,2};
  MEDCouplingUMesh *m1=MEDCouplingUMesh::New();
  m1->setMeshDimension(2);
  m1->allocateCells(1);
  m1->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,m1Conn);
  m1->finishInsertingCells();
  DataArrayDouble *myCoords1=DataArrayDouble::New();
  myCoords1->alloc(4,2);
  std::copy(m1Coords,m1Coords+8,myCoords1->getPointer());
  m1->setCoords(myCoords1);
  myCoords1->decrRef();
  //
  double vec1[3]={0.,0.,1.};
  double *vec2=new double[2];
  for(int i=0;i<18;i++)
    {
      vec2[0]=3.*cos(M_PI/9.*i);
      vec2[1]=3.*sin(M_PI/9.*i);
      MEDCouplingUMesh *m1Cpy=static_cast<MEDCouplingUMesh *>(m1->deepCopy());
      m1Cpy->translate(vec2);
      std::vector<int> res;
      CPPUNIT_ASSERT_THROW(m1Cpy->are2DCellsNotCorrectlyOriented(vec1,false,res),INTERP_KERNEL::Exception);
      res.clear();
      m1Cpy->changeSpaceDimension(3);
      m1Cpy->are2DCellsNotCorrectlyOriented(vec1,false,res);
      CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
      CPPUNIT_ASSERT_EQUAL(0,res[0]);
      m1Cpy->decrRef();
    }
  delete [] vec2;
  //
  m1->decrRef();
}

void MEDCouplingBasicsTest5::testDataArrayAbs1()
{
  DataArrayDouble *d1=DataArrayDouble::New();
  const double val1[12]={2.,-3.,-5.,6.,-7.,-8.,9.,10.,-11.,-12.,-13.,-15.};
  const double expected1[12]={2.,3.,5.,6.,7.,8.,9.,10.,11.,12.,13.,15.};
  d1->alloc(6,2);
  std::copy(val1,val1+12,d1->getPointer());
  MCAuto<DataArrayInt> d2=d1->convertToIntArr();
  //
  d1->abs();
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],d1->getIJ(0,i),1e-14);
  //
  const int expected2[12]={2,3,5,6,7,8,9,10,11,12,13,15};
  d2->abs();
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],d2->getIJ(0,i));
  //
  d1->decrRef();
}

void MEDCouplingBasicsTest5::testGetValueOn3()
{
  const double v[4]={0.,1.,1.5,2.};
  const double v2[5]={0.7,1.25,0.,2.,1.5};
  const double disp[12]={5.,50.,500.,6.,60.,600.,7.,70.,700.,8.,80.,800.};
  MEDCouplingUMesh *m=MEDCouplingUMesh::New("myMesh",1);
  const int nbNodes=4;
  const int nbCells=nbNodes-1;
  m->allocateCells(nbCells);
  DataArrayDouble *coords=DataArrayDouble::New();
  coords->alloc(nbNodes,1);
  std::copy(v,v+nbNodes,coords->getPointer());
  m->setCoords(coords);
  coords->decrRef();
  const int conn[6]={0,1,2,1,2,3};
  m->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
  m->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+2);
  m->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+4);
  m->finishInsertingCells();
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_NODES);
  f->setMesh(m);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(m->getNumberOfNodes(),3);
  std::copy(disp,disp+12,array->getPointer());
  f->setArray(array);
  array->decrRef();
  DataArrayDouble *arr1=f->getValueOnMulti(v2,5);
 CPPUNIT_ASSERT_EQUAL(5,(int)arr1->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(3,(int)arr1->getNumberOfComponents());
  const double expected1[15]={5.7,57.,570.,6.5,65.,650.,5.,50.,500.,8.,80.,800.,7.,70.,700.};
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],arr1->getIJ(0,i),1e-14);
  arr1->decrRef();
  f->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest5::testGetNodeIdsOfCell2()
{
  MEDCouplingCMesh *m1c=MEDCouplingCMesh::New();
  DataArrayDouble *coordsX=DataArrayDouble::New();
  double arrX[5] = { -1., 1., 2., 4., 4.5 };
  coordsX->useArray(arrX,false, CPP_DEALLOC,5,1);
  DataArrayDouble *coordsY=DataArrayDouble::New();
  double arrY[4] = { -2., 2., 4., 8. };
  coordsY->useArray(arrY,false, CPP_DEALLOC,4,1);
  DataArrayDouble *coordsZ=DataArrayDouble::New();
  double arrZ[3] = { -2., 2., 4. };
  coordsZ->useArray(arrZ,false, CPP_DEALLOC,3,1);  
  // test in 1D
  m1c->setCoordsAt(0,coordsX);
  CPPUNIT_ASSERT_EQUAL(4,(int)m1c->getNumberOfCells());
  const int expected1[4][2]={{0,1},{1,2},{2,3},{3,4}};
  for(int i=0;i<4;i++)
    {
      std::vector<int> v;
      m1c->getNodeIdsOfCell(i,v);
      CPPUNIT_ASSERT((int)v.size()==2);
      std::equal(v.begin(),v.end(),expected1[i]);
    }
  // test in 2D
  m1c->setCoordsAt(1,coordsY);
  CPPUNIT_ASSERT_EQUAL(12,(int)m1c->getNumberOfCells());
  const int expected2[12][4]={{0,1,6,5},{1,2,7,6},{2,3,8,7},{3,4,9,8},{4,5,11,10},{5,6,12,11},{6,7,13,12},{7,8,14,13},{8,9,16,15},{9,10,17,16},{10,11,18,17},{11,12,19,18}};
  for(int i=0;i<12;i++)
    {
      std::vector<int> v;
      m1c->getNodeIdsOfCell(i,v);
      CPPUNIT_ASSERT((int)v.size()==4);
      std::equal(v.begin(),v.end(),expected2[i]);
    }
  // test in 3D
  m1c->setCoordsAt(2,coordsZ);
  CPPUNIT_ASSERT_EQUAL(24,(int)m1c->getNumberOfCells());
  const int expected3[24][8]={{0,1,6,5,20,21,26,25},{1,2,7,6,21,22,27,26},{2,3,8,7,22,23,28,27},{3,4,9,8,23,24,29,28},{4,5,11,10,24,25,31,30},{5,6,12,11,25,26,32,31},{6,7,13,12,26,27,33,32},{7,8,14,13,27,28,34,33},{8,9,16,15,28,29,36,35},{9,10,17,16,29,30,37,36},{10,11,18,17,30,31,38,37},{11,12,19,18,31,32,39,38},{20,21,26,25,40,41,46,45},{21,22,27,26,41,42,47,46},{22,23,28,27,42,43,48,47},{23,24,29,28,43,44,49,48},{24,25,31,30,44,45,51,50},{25,26,32,31,45,46,52,51},{26,27,33,32,46,47,53,52},{27,28,34,33,47,48,54,53},{28,29,36,35,48,49,56,55},{29,30,37,36,49,50,57,56},{30,31,38,37,50,51,58,57},{31,32,39,38,51,52,59,58}};
  for(int i=0;i<12;i++)
    {
      std::vector<int> v;
      m1c->getNodeIdsOfCell(i,v);
      CPPUNIT_ASSERT((int)v.size()==8);
      std::equal(v.begin(),v.end(),expected3[i]);
    }
  //
  coordsX->decrRef();
  coordsY->decrRef();
  coordsZ->decrRef();
  m1c->decrRef();
}

void MEDCouplingBasicsTest5::testRenumberNodesInConn1()
{
  double mesh2DCoords[27]={-0.3,-0.3,0., 0.2,-0.3,0., 0.7,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.7,0.2,0., -0.3,0.7,0., 0.2,0.7,0., 0.7,0.7,0. };
  int mesh2DConn[18]={1,4,2, 4,5,2, 0,3,4,1, 6,7,4,3, 7,8,5,4};
  MEDCouplingUMesh *mesh2D=MEDCouplingUMesh::New("mesh",2);
  mesh2D->allocateCells(5);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,mesh2DConn);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,mesh2DConn+3);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,mesh2DConn+6);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,mesh2DConn+10);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,mesh2DConn+14);
  mesh2D->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(9,3);
  std::copy(mesh2DCoords,mesh2DCoords+27,myCoords->getPointer());
  mesh2D->setCoords(myCoords);
  myCoords->decrRef();
  mesh2D->checkConsistencyLight();
  //
  double mesh3DCoords[24]={-0.3,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.2,-0.3,0., -0.3,-0.3,1., -0.3,0.2,1., 0.2,0.2,1., 0.2,-0.3,1. };
  int mesh3DConn[8]={0,1,2,3,4,5,6,7};
  MEDCouplingUMesh *mesh3D=MEDCouplingUMesh::New("mesh",3);
  mesh3D->allocateCells(1);
  mesh3D->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,mesh3DConn);
  mesh3D->finishInsertingCells();
  DataArrayDouble *myCoords3D=DataArrayDouble::New();
  myCoords3D->alloc(8,3);
  std::copy(mesh3DCoords,mesh3DCoords+24,myCoords3D->getPointer());
  mesh3D->setCoords(myCoords3D);
  myCoords3D->decrRef();
  mesh3D->checkConsistencyLight();
  //
  MEDCouplingUMesh *mesh3D_2=dynamic_cast<MEDCouplingUMesh *>(mesh3D->deepCopy());
  MEDCouplingUMesh *mesh2D_2=dynamic_cast<MEDCouplingUMesh *>(mesh2D->deepCopy());
  MEDCouplingUMesh *mesh3D_4=dynamic_cast<MEDCouplingUMesh *>(mesh3D->deepCopy());
  MEDCouplingUMesh *mesh2D_4=dynamic_cast<MEDCouplingUMesh *>(mesh2D->deepCopy());
  DataArrayInt *renumNodes=DataArrayInt::New();
  int oldNbOf3DNodes=mesh3D->getNumberOfNodes();
  renumNodes->alloc(mesh2D->getNumberOfNodes(),1);
  renumNodes->iota(oldNbOf3DNodes);
  DataArrayDouble *coo=DataArrayDouble::Aggregate(mesh3D->getCoords(),mesh2D->getCoords());
  mesh3D->setCoords(coo);
  mesh2D->setCoords(coo);
  coo->decrRef();
  MEDCouplingUMesh *mesh2D_3=dynamic_cast<MEDCouplingUMesh *>(mesh2D->deepCopy());
  mesh2D_3->shiftNodeNumbersInConn(oldNbOf3DNodes);
  mesh2D->renumberNodesInConn(renumNodes->getConstPointer());
  renumNodes->decrRef();
  CPPUNIT_ASSERT(mesh2D_3->isEqual(mesh2D,1e-12));
  mesh2D_3->decrRef();
  //
  DataArrayInt *da1,*da2;
  mesh3D->checkGeoEquivalWith(mesh3D_2,10,1e-12,da1,da2);
  CPPUNIT_ASSERT(da1==0);
 CPPUNIT_ASSERT_EQUAL(8,(int)da2->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da2->getNumberOfComponents());
  const int expected1[8]={8,11,12,9,4,5,6,7};
  for(int i=0;i<8;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],da2->getIJ(i,0));
  da2->decrRef();
  //
  mesh2D->checkGeoEquivalWith(mesh2D_2,10,1e-12,da1,da2);
  CPPUNIT_ASSERT(da1==0);
 CPPUNIT_ASSERT_EQUAL(9,(int)da2->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da2->getNumberOfComponents());
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_EQUAL(8+i,da2->getIJ(i,0));
  da2->decrRef();
  //
  const double vect[3]={1.,0.,0.};
  MEDCouplingUMesh *mesh2D_5=dynamic_cast<MEDCouplingUMesh *>(mesh2D_4->deepCopy());
  mesh2D_5->translate(vect);
  std::vector<MEDCouplingUMesh *> meshes(3);
  meshes[0]=mesh3D_4; meshes[1]=mesh2D_4; meshes[2]=mesh2D_5;
  MEDCouplingUMesh::PutUMeshesOnSameAggregatedCoords(meshes);
  CPPUNIT_ASSERT(mesh3D_4->getCoords()==mesh2D_4->getCoords());
  CPPUNIT_ASSERT(mesh2D_4->getCoords()==mesh2D_5->getCoords());
  mesh3D_4->checkConsistencyLight(); mesh2D_4->checkConsistencyLight(); mesh2D_5->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(26,mesh3D_4->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,mesh3D_4->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(9,(int)mesh3D_4->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(23,(int)mesh2D_4->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(23,(int)mesh2D_5->getNodalConnectivity()->getNumberOfTuples());
  const int expected2[9]={18,0,1,2,3,4,5,6,7};
  const int expected3[23]={3,9,12,10, 3,12,13,10, 4,8,11,12,9, 4,14,15,12,11, 4,15,16,13,12};
  const int expected4[23]={3,18,21,19, 3,21,22,19, 4,17,20,21,18, 4,23,24,21,20, 4,24,25,22,21};
  const double expected5[78]={-0.3,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.2,-0.3,0., -0.3,-0.3,1., -0.3,0.2,1., 0.2,0.2,1., 0.2,-0.3,1., -0.3,-0.3,0., 0.2,-0.3,0., 0.7,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.7,0.2,0., -0.3,0.7,0., 0.2,0.7,0., 0.7,0.7,0., 0.7, -0.3, 0.0, 1.2, -0.3, 0.0, 1.7, -0.3, 0.0, 0.7, 0.2, 0.0, 1.2, 0.2, 0.0, 1.7, 0.2, 0.0, 0.7, 0.7, 0.0, 1.2, 0.7, 0.0, 1.7, 0.7, 0.0};
  CPPUNIT_ASSERT(std::equal(expected2,expected2+9,mesh3D_4->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected3,expected3+23,mesh2D_4->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+23,mesh2D_5->getNodalConnectivity()->getConstPointer()));
  for(int i=0;i<78;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected5[i],mesh3D_4->getCoords()->getIJ(0,i),1e-12);
  //
  MEDCouplingUMesh::MergeNodesOnUMeshesSharingSameCoords(meshes,1e-12);
  mesh3D_4->checkConsistencyLight(); mesh2D_4->checkConsistencyLight(); mesh2D_5->checkConsistencyLight();
  CPPUNIT_ASSERT(mesh3D_4->getCoords()==mesh2D_4->getCoords());
  CPPUNIT_ASSERT(mesh2D_4->getCoords()==mesh2D_5->getCoords());
  CPPUNIT_ASSERT_EQUAL(19,mesh3D_4->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,mesh3D_4->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(9,(int)mesh3D_4->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(23,(int)mesh2D_4->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(23,(int)mesh2D_5->getNodalConnectivity()->getNumberOfTuples());
  const int expected6[9]={18,0,1,2,3,4,5,6,7};
  const int expected7[23]={3,3,2,8, 3,2,9,8, 4,0,1,2,3, 4,10,11,2,1, 4,11,12,9,2};
  const int expected8[23]={3,13,15,14, 3,15,16,14, 4,8,9,15,13, 4,12,17,15,9, 4,17,18,16,15};
  const double expected9[57]={-0.3, -0.3, 0., -0.3, 0.2, 0., 0.2, 0.2, 0., 0.2, -0.3, 0., -0.3, -0.3, 1., -0.3, 0.2, 1., 
                              0.2, 0.2, 1., 0.2, -0.3, 1., 0.7, -0.3, 0., 0.7, 0.2, 0., -0.3, 0.7, 0., 0.2, 0.7, 0., 
                              0.7, 0.7, 0., 1.2, -0.3, 0., 1.7, -0.3, 0., 1.2, 0.2, 0., 1.7, 0.2, 0., 1.2, 0.7, 0., 1.7, 0.7, 0.};
  CPPUNIT_ASSERT(std::equal(expected6,expected6+9,mesh3D_4->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected7,expected7+23,mesh2D_4->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected8,expected8+23,mesh2D_5->getNodalConnectivity()->getConstPointer()));
  for(int i=0;i<57;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected9[i],mesh3D_4->getCoords()->getIJ(0,i),1e-12);
  mesh2D_5->decrRef();
  //
  mesh3D_4->decrRef();
  mesh2D_4->decrRef();
  mesh3D_2->decrRef();
  mesh2D_2->decrRef();
  //
  mesh3D->decrRef();
  mesh2D->decrRef();
}

void MEDCouplingBasicsTest5::testComputeNeighborsOfCells1()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  DataArrayInt *d1=0,*d2=0;
  m->computeNeighborsOfCells(d1,d2);
 CPPUNIT_ASSERT_EQUAL(6,(int)d2->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(10,(int)d1->getNumberOfTuples());
  const int expected1[6]={0,2,4,6,8,10};
  const int expected2[10]={3,1,0,2,4,1,4,0,2,3};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+6,d2->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+10,d1->getConstPointer()));
  d1->decrRef();
  d2->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest5::testCheckButterflyCellsBug1()
{
  double mesh2DCoords[10]={323.85,120.983748908684,317.5,131.982271536747,336.55,120.983748908686,330.2,131.982271536751,323.85,142.98079416481};
  int mesh2DConn[5]={4,1,0,2,3};
  MEDCouplingUMesh *mesh2D=MEDCouplingUMesh::New("mesh",2);
  mesh2D->allocateCells(1);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_POLYGON,5,mesh2DConn);
  mesh2D->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(5,2);
  std::copy(mesh2DCoords,mesh2DCoords+10,myCoords->getPointer());
  mesh2D->setCoords(myCoords);
  myCoords->decrRef();
  mesh2D->checkConsistencyLight();
  //
  std::vector<int> v;
  mesh2D->checkButterflyCells(v);
  CPPUNIT_ASSERT_EQUAL(0,(int)v.size());
  //
  mesh2D->decrRef();
}

void MEDCouplingBasicsTest5::testDataArrayIntRange1()
{
  DataArrayInt *d=DataArrayInt::Range(2,17,7);
  const int expected1[3]={2,9,16};
 CPPUNIT_ASSERT_EQUAL(3,(int)d->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,d->getConstPointer()));
  d->decrRef();
  //
  d=DataArrayInt::Range(2,23,7);
 CPPUNIT_ASSERT_EQUAL(3,(int)d->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,d->getConstPointer()));
  d->decrRef();
  //
  d=DataArrayInt::Range(2,24,7);
  const int expected2[4]={2,9,16,23};
 CPPUNIT_ASSERT_EQUAL(4,(int)d->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected2,expected2+4,d->getConstPointer()));
  d->decrRef();
  //
  d=DataArrayInt::Range(24,2,-7);
  const int expected3[4]={24,17,10,3};
 CPPUNIT_ASSERT_EQUAL(4,(int)d->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+4,d->getConstPointer()));
  d->decrRef();
  //
  d=DataArrayInt::Range(23,2,-7);
  const int expected4[3]={23,16,9};
 CPPUNIT_ASSERT_EQUAL(3,(int)d->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected4,expected4+3,d->getConstPointer()));
  d->decrRef();
  //
  d=DataArrayInt::Range(23,22,-7);
 CPPUNIT_ASSERT_EQUAL(1,(int)d->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)d->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(23,d->getIJ(0,0));
  d->decrRef();
  //
  d=DataArrayInt::Range(22,23,7);
 CPPUNIT_ASSERT_EQUAL(1,(int)d->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)d->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(22,d->getIJ(0,0));
  d->decrRef();
  //
  d=DataArrayInt::Range(22,22,7);
 CPPUNIT_ASSERT_EQUAL(0,(int)d->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)d->getNumberOfComponents());
  d->decrRef();
  //
  d=DataArrayInt::Range(22,22,-7);
 CPPUNIT_ASSERT_EQUAL(0,(int)d->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)d->getNumberOfComponents());
  d->decrRef();
  //
  CPPUNIT_ASSERT_THROW(DataArrayInt::Range(22,23,-7),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(DataArrayInt::Range(23,22,7),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(DataArrayInt::Range(23,22,0),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(DataArrayInt::Range(22,23,0),INTERP_KERNEL::Exception);
}

void MEDCouplingBasicsTest5::testDataArrayDoubleGetMinMaxPerComponent1()
{
  const double values1[12]={1.,2.,3.,-0.9,2.1,3.,1.3,1.7,3.,1.,1.8,3.};
  DataArrayDouble *d1=DataArrayDouble::New();
  double *res=new double[2*3];
  CPPUNIT_ASSERT_THROW(d1->getMinMaxPerComponent(res),INTERP_KERNEL::Exception);
  d1->alloc(4,3);
  std::copy(values1,values1+12,d1->getPointer());
  d1->getMinMaxPerComponent(res);
  const double expected1[6]={-0.9,1.3,1.7,2.1,3.,3.};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],res[i],1e-14);
  delete [] res;
  //
  d1->rearrange(2);
  res=new double[2*2];
  d1->getMinMaxPerComponent(res);
  const double expected2[4]={1.,3.,-0.9,3.};
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],res[i],1e-14);
  delete [] res;
  //
  d1->rearrange(1);
  res=new double[2*1];
  d1->getMinMaxPerComponent(res);
  const double expected3[2]={-0.9,3.};
  for(int i=0;i<2;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],res[i],1e-14);
  delete [] res;
  d1->decrRef();
}

void MEDCouplingBasicsTest5::testDataArrayIntGetHashCode1()
{
  DataArrayInt *d1=DataArrayInt::New(); d1->alloc(3545,1); d1->iota(0);
  DataArrayInt *d2=DataArrayInt::New(); d2->alloc(3545,1); d2->iota(0);
  //
  CPPUNIT_ASSERT_EQUAL(d1->getHashCode(),d2->getHashCode());
  CPPUNIT_ASSERT_EQUAL(232341068,d1->getHashCode());
  d1->setIJ(886,0,6);
  CPPUNIT_ASSERT_EQUAL(232340188,d1->getHashCode());
  //
  d1->decrRef();
  d2->decrRef();
}

void MEDCouplingBasicsTest5::testZipConnectivityPol1()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  const int cells1[3]={2,3,4};
  MEDCouplingPointSet *m2_1=m1->buildPartOfMySelf(cells1,cells1+3,true);
  MEDCouplingUMesh *m2=dynamic_cast<MEDCouplingUMesh *>(m2_1);
  DataArrayInt *arr=0;
  CPPUNIT_ASSERT(m2);
  // no permutation policy 0
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,0,arr));
 CPPUNIT_ASSERT_EQUAL(3,(int)arr->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells1,cells1+3,arr->getConstPointer()));
  arr->decrRef();
  // no permutation policy 1
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,1,arr));
 CPPUNIT_ASSERT_EQUAL(3,(int)arr->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells1,cells1+3,arr->getConstPointer()));
  arr->decrRef();
  // no permutation policy 2
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,2,arr));
 CPPUNIT_ASSERT_EQUAL(3,(int)arr->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells1,cells1+3,arr->getConstPointer()));
  arr->decrRef();
  // some modification into m2
  const int modif1[3]={2,4,5};
  std::copy(modif1,modif1+3,m2->getNodalConnectivity()->getPointer()+1);
  //policy 0 fails because cell0 in m2 has same orientation be not same connectivity
  const int expected1[3]={5,3,4};
  CPPUNIT_ASSERT(!m1->areCellsIncludedIn(m2,0,arr));
 CPPUNIT_ASSERT_EQUAL(3,(int)arr->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,arr->getConstPointer()));
  arr->decrRef();
  //policy 1 succeeds because cell0 in m2 has not exactly the same conn
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,1,arr));
 CPPUNIT_ASSERT_EQUAL(3,(int)arr->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells1,cells1+3,arr->getConstPointer()));
  arr->decrRef();
  //policy 2 succeeds because cell0 in m2 has same nodes in connectivity
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,2,arr));
 CPPUNIT_ASSERT_EQUAL(3,(int)arr->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells1,cells1+3,arr->getConstPointer()));
  arr->decrRef();
  //some new modification into m2
  const int modif2[3]={2,5,4};
  std::copy(modif2,modif2+3,m2->getNodalConnectivity()->getPointer()+1);
  //policy 0 fails because cell0 in m2 has not exactly the same conn
  CPPUNIT_ASSERT(!m1->areCellsIncludedIn(m2,0,arr));
 CPPUNIT_ASSERT_EQUAL(3,(int)arr->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,arr->getConstPointer()));
  arr->decrRef();
  //policy 1 fails too because cell0 in m2 has not same orientation
  CPPUNIT_ASSERT(!m1->areCellsIncludedIn(m2,1,arr));
 CPPUNIT_ASSERT_EQUAL(3,(int)arr->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,arr->getConstPointer()));
  arr->decrRef();
  //policy 2 succeeds because cell0 in m2 has same nodes in connectivity
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,2,arr));
 CPPUNIT_ASSERT_EQUAL(3,(int)arr->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells1,cells1+3,arr->getConstPointer()));
  arr->decrRef();
  m1->decrRef();
  m2->decrRef();
  // Now 1D
  const int cells2[2]={3,2};
  m1=build1DSourceMesh_2();
  m2_1=m1->buildPartOfMySelf(cells2,cells2+2,true);
  m2=dynamic_cast<MEDCouplingUMesh *>(m2_1);
  CPPUNIT_ASSERT(m2);
  arr=0;
  // no permutation policy 0
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,0,arr));
 CPPUNIT_ASSERT_EQUAL(2,(int)arr->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells2,cells2+2,arr->getConstPointer()));
  arr->decrRef();
  // no permutation policy 1
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,1,arr));
 CPPUNIT_ASSERT_EQUAL(2,(int)arr->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells2,cells2+2,arr->getConstPointer()));
  arr->decrRef();
  // no permutation policy 2
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,2,arr));
 CPPUNIT_ASSERT_EQUAL(2,(int)arr->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells2,cells2+2,arr->getConstPointer()));
  arr->decrRef();
  // some modification into m2
  const int modif3[2]={4,3};
  std::copy(modif3,modif3+2,m2->getNodalConnectivity()->getPointer()+1);
  //policy 0 fails because cell0 in m2 has not exactly the same conn
  const int expected2[2]={4,2};
  CPPUNIT_ASSERT(!m1->areCellsIncludedIn(m2,0,arr));
 CPPUNIT_ASSERT_EQUAL(2,(int)arr->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected2,expected2+2,arr->getConstPointer()));
  arr->decrRef();
  //policy 1 fails too because cell0 in m2 has not same orientation
  CPPUNIT_ASSERT(!m1->areCellsIncludedIn(m2,1,arr));
 CPPUNIT_ASSERT_EQUAL(2,(int)arr->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected2,expected2+2,arr->getConstPointer()));
  arr->decrRef();
  //policy 2 succeeds because cell0 in m2 has same nodes in connectivity
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,2,arr));
 CPPUNIT_ASSERT_EQUAL(2,(int)arr->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells2,cells2+2,arr->getConstPointer()));
  arr->decrRef();
  m1->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest5::testConvexEnvelop2D1()
{
  const double coords[662]={7.54758495819e-14,-1.12270326253e-12,8.43143594193,-1.02835845055e-12,4.21571797096,7.30183771609,-4.21571797097,7.30183771609,-8.43143594193,-1.09439981894e-12,-4.21571797097,-7.30183771609,4.21571797097,-7.30183771609,16.8628718839,-1.02835845055e-12,12.6471539129,7.30183771609,8.43143594193,14.6036754322,2.26427548746e-13,14.6036754322,-8.43143594193,14.6036754322,-12.6471539129,7.30183771609,-16.8628718839,-1.39630321727e-12,-12.6471539129,-7.30183771609,-8.43143594193,-14.6036754322,3.7737924791e-14,-14.6036754322,8.43143594193,-14.6036754322,12.6471539129,-7.30183771609,25.2943078258,-1.07553085654e-12,21.0785898548,7.30183771609,16.8628718839,14.6036754322,12.6471539129,21.9055131483,4.21571797096,21.9055131483,-4.21571797097,21.9055131483,-12.6471539129,21.9055131483,-16.8628718839,14.6036754322,-21.0785898548,7.30183771609,-25.2943078258,-1.02835845055e-12,-21.0785898548,-7.30183771609,-16.8628718839,-14.6036754322,-12.6471539129,-21.9055131483,-4.21571797097,-21.9055131483,4.21571797097,-21.9055131483,12.6471539129,-21.9055131483,16.8628718839,-14.6036754322,21.0785898548,-7.30183771609,33.7257437677,-7.45324014622e-13,29.5100257968,7.30183771609,25.2943078258,14.6036754322,21.0785898548,21.9055131483,16.8628718839,29.2073508644,8.43143594193,29.2073508644,-1.20761359331e-12,29.2073508644,-8.43143594193,29.2073508644,-16.8628718839,29.2073508644,-21.0785898548,21.9055131483,-25.2943078258,14.6036754322,-29.5100257968,7.30183771609,-33.7257437677,-7.26455052226e-13,-29.5100257968,-7.30183771609,-25.2943078258,-14.6036754322,-21.0785898548,-21.9055131483,-16.8628718839,-29.2073508644,-8.43143594193,-29.2073508644,4.15117172701e-13,-29.2073508644,8.43143594193,-29.2073508644,16.8628718839,-29.2073508644,21.0785898548,-21.9055131483,25.2943078258,-14.6036754322,29.5100257968,-7.30183771609,42.1571797097,-1.86802727715e-12,37.9414617387,7.30183771609,33.7257437677,14.6036754322,29.5100257968,21.9055131483,25.2943078258,29.2073508644,21.0785898548,36.5091885805,12.6471539129,36.5091885805,4.21571797096,36.5091885805,-4.21571797096,36.5091885805,-12.6471539129,36.5091885805,-21.0785898548,36.5091885805,-25.2943078258,29.2073508644,-29.5100257968,21.9055131483,-33.7257437677,14.6036754322,-37.9414617387,7.30183771609,-42.1571797097,-9.81186044565e-13,-37.9414617387,-7.30183771609,-33.7257437677,-14.6036754322,-29.5100257968,-21.9055131483,-25.2943078258,-29.2073508644,-21.0785898548,-36.5091885805,-12.6471539129,-36.5091885805,-4.21571797097,-36.5091885805,4.21571797097,-36.5091885805,12.6471539129,-36.5091885805,21.0785898548,-36.5091885805,25.2943078258,-29.2073508644,29.5100257968,-21.9055131483,33.7257437677,-14.6036754322,37.9414617387,-7.30183771609,50.5886156516,-6.98151608633e-13,46.3728976806,7.30183771609,42.1571797097,14.6036754322,37.9414617387,21.9055131483,33.7257437677,29.2073508644,29.5100257968,36.5091885805,25.2943078258,43.8110262966,16.8628718839,43.8110262966,8.43143594193,43.8110262966,-1.84915831476e-12,43.8110262966,-8.43143594193,43.8110262966,-16.8628718839,43.8110262966,-25.2943078258,43.8110262966,-29.5100257968,36.5091885805,-33.7257437677,29.2073508644,-37.9414617387,21.9055131483,-42.1571797097,14.6036754322,-46.3728976806,7.30183771609,-50.5886156516,-1.47177906685e-12,-46.3728976806,-7.30183771609,-42.1571797097,-14.6036754322,-37.9414617387,-21.9055131483,-33.7257437677,-29.2073508644,-29.5100257968,-36.5091885805,-25.2943078258,-43.8110262966,-16.8628718839,-43.8110262966,-8.43143594193,-43.8110262966,7.54758495819e-14,-43.8110262966,8.43143594193,-43.8110262966,16.8628718839,-43.8110262966,25.2943078258,-43.8110262966,29.5100257968,-36.5091885805,33.7257437677,-29.2073508644,37.9414617387,-21.9055131483,42.1571797097,-14.6036754322,46.3728976806,-7.30183771609,59.0200515935,-7.9249642061e-13,54.8043336225,7.30183771609,50.5886156516,14.6036754322,46.3728976806,21.9055131483,42.1571797097,29.2073508644,37.9414617387,36.5091885805,33.7257437677,43.8110262966,29.5100257968,51.1128640127,21.0785898548,51.1128640127,12.6471539129,51.1128640127,4.21571797096,51.1128640127,-4.21571797096,51.1128640127,-12.6471539129,51.1128640127,-21.0785898548,51.1128640127,-29.5100257968,51.1128640127,-33.7257437677,43.8110262966,-37.9414617387,36.5091885805,-42.1571797097,29.2073508644,-46.3728976806,21.9055131483,-50.5886156516,14.6036754322,-54.8043336226,7.30183771609,-59.0200515935,-1.31139288649e-12,-54.8043336226,-7.30183771609,-50.5886156516,-14.6036754322,-46.3728976806,-21.9055131483,-42.1571797097,-29.2073508644,-37.9414617387,-36.5091885805,-33.7257437677,-43.8110262966,-29.5100257968,-51.1128640127,-21.0785898548,-51.1128640127,-12.6471539129,-51.1128640127,-4.21571797097,-51.1128640127,4.21571797097,-51.1128640127,12.6471539129,-51.1128640127,21.0785898548,-51.1128640127,29.5100257968,-51.1128640127,33.7257437677,-43.8110262966,37.9414617387,-36.5091885805,42.1571797097,-29.2073508644,46.3728976806,-21.9055131483,50.5886156516,-14.6036754322,54.8043336225,-7.30183771609,67.4514875354,-2.14162723189e-12,63.2357695645,7.30183771609,59.0200515935,14.6036754322,54.8043336226,21.9055131483,50.5886156516,29.2073508644,46.3728976806,36.5091885805,42.1571797097,43.8110262966,37.9414617387,51.1128640127,33.7257437677,58.4147017287,25.2943078258,58.4147017287,16.8628718839,58.4147017287,8.43143594193,58.4147017287,6.79282646237e-13,58.4147017287,-8.43143594193,58.4147017287,-16.8628718839,58.4147017287,-25.2943078258,58.4147017287,-33.7257437677,58.4147017287,-37.9414617387,51.1128640127,-42.1571797097,43.8110262966,-46.3728976806,36.5091885805,-50.5886156516,29.2073508644,-54.8043336226,21.9055131483,-59.0200515935,14.6036754322,-63.2357695645,7.30183771609,-67.4514875354,-1.16044118732e-12,-63.2357695645,-7.30183771609,-59.0200515935,-14.6036754322,-54.8043336226,-21.9055131483,-50.5886156516,-29.2073508644,-46.3728976806,-36.5091885805,-42.1571797097,-43.8110262966,-37.9414617387,-51.1128640127,-33.7257437677,-58.4147017287,-25.2943078258,-58.4147017287,-16.8628718839,-58.4147017287,-8.43143594193,-58.4147017287,-5.66068871864e-14,-58.4147017287,8.43143594193,-58.4147017287,16.8628718839,-58.4147017287,25.2943078258,-58.4147017287,33.7257437677,-58.4147017287,37.9414617387,-51.1128640127,42.1571797097,-43.8110262966,46.3728976806,-36.5091885805,50.5886156516,-29.2073508644,54.8043336226,-21.9055131483,59.0200515935,-14.6036754322,63.2357695645,-7.30183771609,75.8829234774,-2.29257893105e-12,71.6672055064,7.30183771609,67.4514875354,14.6036754322,63.2357695645,21.9055131483,59.0200515935,29.2073508644,54.8043336226,36.5091885805,50.5886156516,43.8110262966,46.3728976806,51.1128640127,42.1571797097,58.4147017287,37.9414617387,65.7165394448,29.5100257968,65.7165394448,21.0785898548,65.7165394448,12.6471539129,65.7165394448,4.21571797097,65.7165394448,-4.21571797096,65.7165394448,-12.6471539129,65.7165394448,-21.0785898548,65.7165394448,-29.5100257968,65.7165394448,-37.9414617387,65.7165394448,-42.1571797097,58.4147017287,-46.3728976806,51.1128640127,-50.5886156516,43.8110262966,-54.8043336226,36.5091885805,-59.0200515935,29.2073508644,-63.2357695645,21.9055131483,-67.4514875354,14.6036754322,-71.6672055064,7.30183771609,-75.8829234774,-1.31139288649e-12,-71.6672055064,-7.30183771609,-67.4514875354,-14.6036754322,-63.2357695645,-21.9055131483,-59.0200515935,-29.2073508644,-54.8043336226,-36.5091885805,-50.5886156516,-43.8110262966,-46.3728976806,-51.1128640127,-42.1571797097,-58.4147017287,-37.9414617387,-65.7165394448,-29.5100257968,-65.7165394448,-21.0785898548,-65.7165394448,-12.6471539129,-65.7165394448,-4.21571797097,-65.7165394448,4.21571797097,-65.7165394448,12.6471539129,-65.7165394448,21.0785898548,-65.7165394448,29.5100257968,-65.7165394448,37.9414617387,-65.7165394448,42.1571797097,-58.4147017287,46.3728976806,-51.1128640127,50.5886156516,-43.8110262966,54.8043336226,-36.5091885805,59.0200515935,-29.2073508644,63.2357695645,-21.9055131483,67.4514875354,-14.6036754322,71.6672055064,-7.30183771609,84.3143594193,-1.49064802924e-12,80.0986414483,7.30183771609,75.8829234774,14.6036754322,71.6672055064,21.9055131483,67.4514875354,29.2073508644,63.2357695645,36.5091885805,59.0200515935,43.8110262966,54.8043336226,51.1128640127,50.5886156516,58.4147017287,46.3728976806,65.7165394448,42.1571797097,73.0183771609,33.7257437677,73.0183771609,25.2943078258,73.0183771609,16.8628718839,73.0183771609,8.43143594193,73.0183771609,2.0755858635e-12,73.0183771609,-8.43143594193,73.0183771609,-16.8628718839,73.0183771609,-25.2943078258,73.0183771609,-33.7257437677,73.0183771609,-42.1571797097,73.0183771609,-46.3728976806,65.7165394448,-50.5886156516,58.4147017287,-54.8043336226,51.1128640127,-59.0200515935,43.8110262966,-63.2357695645,36.5091885805,-67.4514875354,29.2073508644,-71.6672055064,21.9055131483,-75.8829234774,14.6036754322,-80.0986414483,7.30183771609,-84.3143594193,-1.11326878133e-12,-80.0986414483,-7.30183771609,-75.8829234774,-14.6036754322,-71.6672055064,-21.9055131483,-67.4514875354,-29.2073508644,-63.2357695645,-36.5091885805,-59.0200515935,-43.8110262966,-54.8043336226,-51.1128640127,-50.5886156516,-58.4147017287,-46.3728976806,-65.7165394448,-42.1571797097,-73.0183771609,-33.7257437677,-73.0183771609,-25.2943078258,-73.0183771609,-16.8628718839,-73.0183771609,-8.43143594193,-73.0183771609,-5.66068871864e-14,-73.0183771609,8.43143594193,-73.0183771609,16.8628718839,-73.0183771609,25.2943078258,-73.0183771609,33.7257437677,-73.0183771609,42.1571797097,-73.0183771609,46.3728976806,-65.7165394448,50.5886156516,-58.4147017287,54.8043336226,-51.1128640127,59.0200515935,-43.8110262966,63.2357695645,-36.5091885805,67.4514875354,-29.2073508644,71.6672055064,-21.9055131483,75.8829234774,-14.6036754322,80.0986414483,-7.3018377161};
  const int conn[2137]={0,2,3,4,5,6,1,1,8,2,0,6,18,7,2,9,10,3,0,1,8,3,10,11,12,4,0,2,4,3,12,13,14,5,0,5,0,4,14,15,16,6,6,1,0,5,16,17,18,7,20,8,1,18,36,19,8,21,9,2,1,7,20,9,22,23,10,2,8,21,10,23,24,11,3,2,9,11,24,25,26,12,3,10,12,11,26,27,13,4,3,13,12,27,28,29,14,4,14,4,13,29,30,15,5,15,5,14,30,31,32,16,16,6,5,15,32,33,17,17,18,6,16,33,34,35,18,7,1,6,17,35,36,19,38,20,7,36,60,37,20,39,21,8,7,19,38,21,40,22,9,8,20,39,22,41,42,23,9,21,40,23,42,43,24,10,9,22,24,43,44,25,11,10,23,25,44,45,46,26,11,24,26,25,46,47,27,12,11,27,26,47,48,28,13,12,28,27,48,49,50,29,13,29,13,28,50,51,30,14,30,14,29,51,52,31,15,31,15,30,52,53,54,32,32,16,15,31,54,55,33,33,17,16,32,55,56,34,34,35,17,33,56,57,58,35,36,18,17,34,58,59,36,19,7,18,35,59,60,37,62,38,19,60,90,61,38,63,39,20,19,37,62,39,64,40,21,20,38,63,40,65,41,22,21,39,64,41,66,67,42,22,40,65,42,67,68,43,23,22,41,43,68,69,44,24,23,42,44,69,70,45,25,24,43,45,70,71,72,46,25,44,46,45,72,73,47,26,25,47,46,73,74,48,27,26,48,47,74,75,49,28,27,49,48,75,76,77,50,28,50,28,49,77,78,51,29,51,29,50,78,79,52,30,52,30,51,79,80,53,31,53,31,52,80,81,82,54,54,32,31,53,82,83,55,55,33,32,54,83,84,56,56,34,33,55,84,85,57,57,58,34,56,85,86,87,58,59,35,34,57,87,88,59,60,36,35,58,88,89,60,37,19,36,59,89,90,61,92,62,37,90,126,91,62,93,63,38,37,61,92,63,94,64,39,38,62,93,64,95,65,40,39,63,94,65,96,66,41,40,64,95,66,97,98,67,41,65,96,67,98,99,68,42,41,66,68,99,100,69,43,42,67,69,100,101,70,44,43,68,70,101,102,71,45,44,69,71,102,103,104,72,45,70,72,71,104,105,73,46,45,73,72,105,106,74,47,46,74,73,106,107,75,48,47,75,74,107,108,76,49,48,76,75,108,109,110,77,49,77,49,76,110,111,78,50,78,50,77,111,112,79,51,79,51,78,112,113,80,52,80,52,79,113,114,81,53,81,53,80,114,115,116,82,82,54,53,81,116,117,83,83,55,54,82,117,118,84,84,56,55,83,118,119,85,85,57,56,84,119,120,86,86,87,57,85,120,121,122,87,88,58,57,86,122,123,88,89,59,58,87,123,124,89,90,60,59,88,124,125,90,61,37,60,89,125,126,91,128,92,61,126,168,127,92,129,93,62,61,91,128,93,130,94,63,62,92,129,94,131,95,64,63,93,130,95,132,96,65,64,94,131,96,133,97,66,65,95,132,97,134,135,98,66,96,133,98,135,136,99,67,66,97,99,136,137,100,68,67,98,100,137,138,101,69,68,99,101,138,139,102,70,69,100,102,139,140,103,71,70,101,103,140,141,142,104,71,102,104,103,142,143,105,72,71,105,104,143,144,106,73,72,106,105,144,145,107,74,73,107,106,145,146,108,75,74,108,107,146,147,109,76,75,109,108,147,148,149,110,76,110,76,109,149,150,111,77,111,77,110,150,151,112,78,112,78,111,151,152,113,79,113,79,112,152,153,114,80,114,80,113,153,154,115,81,115,81,114,154,155,156,116,116,82,81,115,156,157,117,117,83,82,116,157,158,118,118,84,83,117,158,159,119,119,85,84,118,159,160,120,120,86,85,119,160,161,121,121,122,86,120,161,162,163,122,123,87,86,121,163,164,123,124,88,87,122,164,165,124,125,89,88,123,165,166,125,126,90,89,124,166,167,126,91,61,90,125,167,168,127,170,128,91,168,216,169,128,171,129,92,91,127,170,129,172,130,93,92,128,171,130,173,131,94,93,129,172,131,174,132,95,94,130,173,132,175,133,96,95,131,174,133,176,134,97,96,132,175,134,177,178,135,97,133,176,135,178,179,136,98,97,134,136,179,180,137,99,98,135,137,180,181,138,100,99,136,138,181,182,139,101,100,137,139,182,183,140,102,101,138,140,183,184,141,103,102,139,141,184,185,186,142,103,140,142,141,186,187,143,104,103,143,142,187,188,144,105,104,144,143,188,189,145,106,105,145,144,189,190,146,107,106,146,145,190,191,147,108,107,147,146,191,192,148,109,108,148,147,192,193,194,149,109,149,109,148,194,195,150,110,150,110,149,195,196,151,111,151,111,150,196,197,152,112,152,112,151,197,198,153,113,153,113,152,198,199,154,114,154,114,153,199,200,155,115,155,115,154,200,201,202,156,156,116,115,155,202,203,157,157,117,116,156,203,204,158,158,118,117,157,204,205,159,159,119,118,158,205,206,160,160,120,119,159,206,207,161,161,121,120,160,207,208,162,162,163,121,161,208,209,210,163,164,122,121,162,210,211,164,165,123,122,163,211,212,165,166,124,123,164,212,213,166,167,125,124,165,213,214,167,168,126,125,166,214,215,168,127,91,126,167,215,216,169,218,170,127,216,270,217,170,219,171,128,127,169,218,171,220,172,129,128,170,219,172,221,173,130,129,171,220,173,222,174,131,130,172,221,174,223,175,132,131,173,222,175,224,176,133,132,174,223,176,225,177,134,133,175,224,177,226,227,178,134,176,225,178,227,228,179,135,134,177,179,228,229,180,136,135,178,180,229,230,181,137,136,179,181,230,231,182,138,137,180,182,231,232,183,139,138,181,183,232,233,184,140,139,182,184,233,234,185,141,140,183,185,234,235,236,186,141,184,186,185,236,237,187,142,141,187,186,237,238,188,143,142,188,187,238,239,189,144,143,189,188,239,240,190,145,144,190,189,240,241,191,146,145,191,190,241,242,192,147,146,192,191,242,243,193,148,147,193,192,243,244,245,194,148,194,148,193,245,246,195,149,195,149,194,246,247,196,150,196,150,195,247,248,197,151,197,151,196,248,249,198,152,198,152,197,249,250,199,153,199,153,198,250,251,200,154,200,154,199,251,252,201,155,201,155,200,252,253,254,202,202,156,155,201,254,255,203,203,157,156,202,255,256,204,204,158,157,203,256,257,205,205,159,158,204,257,258,206,206,160,159,205,258,259,207,207,161,160,206,259,260,208,208,162,161,207,260,261,209,209,210,162,208,261,262,263,210,211,163,162,209,263,264,211,212,164,163,210,264,265,212,213,165,164,211,265,266,213,214,166,165,212,266,267,214,215,167,166,213,267,268,215,216,168,167,214,268,269,216,169,127,168,215,269,270,217,272,218,169,270,330,271,218,273,219,170,169,217,272,219,274,220,171,170,218,273,220,275,221,172,171,219,274,221,276,222,173,172,220,275,222,277,223,174,173,221,276,223,278,224,175,174,222,277,224,279,225,176,175,223,278,225,280,226,177,176,224,279,226,281,282,227,177,225,280,227,282,283,228,178,177,226,228,283,284,229,179,178,227,229,284,285,230,180,179,228,230,285,286,231,181,180,229,231,286,287,232,182,181,230,232,287,288,233,183,182,231,233,288,289,234,184,183,232,234,289,290,235,185,184,233,235,290,291,292,236,185,234,236,235,292,293,237,186,185,237,236,293,294,238,187,186,238,237,294,295,239,188,187,239,238,295,296,240,189,188,240,239,296,297,241,190,189,241,240,297,298,242,191,190,242,241,298,299,243,192,191,243,242,299,300,244,193,192,244,243,300,301,302,245,193,245,193,244,302,303,246,194,246,194,245,303,304,247,195,247,195,246,304,305,248,196,248,196,247,305,306,249,197,249,197,248,306,307,250,198,250,198,249,307,308,251,199,251,199,250,308,309,252,200,252,200,251,309,310,253,201,253,201,252,310,311,312,254,254,202,201,253,312,313,255,255,203,202,254,313,314,256,256,204,203,255,314,315,257,257,205,204,256,315,316,258,258,206,205,257,316,317,259,259,207,206,258,317,318,260,260,208,207,259,318,319,261,261,209,208,260,319,320,262,262,263,209,261,320,321,322,263,264,210,209,262,322,323,264,265,211,210,263,323,324,265,266,212,211,264,324,325,266,267,213,212,265,325,326,267,268,214,213,266,326,327,268,269,215,214,267,327,328,269,270,216,215,268,328,329,270,217,169,216,269,329,330,271,272,217,330,273,218,217,271,274,219,218,272,275,220,219,273,276,221,220,274,277,222,221,275,278,223,222,276,279,224,223,277,280,225,224,278,281,226,225,279,281,282,226,280,283,227,226,281,284,228,227,282,285,229,228,283,286,230,229,284,287,231,230,285,288,232,231,286,289,233,232,287,290,234,233,288,291,235,234,289,291,292,235,290,291,293,236,235,292,294,237,236,293,295,238,237,294,296,239,238,295,297,240,239,296,298,241,240,297,299,242,241,298,300,243,242,299,301,244,243,301,300,302,244,244,301,303,245,245,302,304,246,246,303,305,247,247,304,306,248,248,305,307,249,249,306,308,250,250,307,309,251,251,308,310,252,252,309,311,253,311,253,310,312,254,253,311,313,255,254,312,314,256,255,313,315,257,256,314,316,258,257,315,317,259,258,316,318,260,259,317,319,261,260,318,320,262,261,319,321,321,322,262,320,323,263,262,321,324,264,263,322,325,265,264,323,326,266,265,324,327,267,266,325,328,268,267,326,329,269,268,327,330,270,269,328,271,217,270,329};
  const int connI[332]={0,7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126,133,140,147,154,161,168,175,182,189,196,203,210,217,224,231,238,245,252,259,266,273,280,287,294,301,308,315,322,329,336,343,350,357,364,371,378,385,392,399,406,413,420,427,434,441,448,455,462,469,476,483,490,497,504,511,518,525,532,539,546,553,560,567,574,581,588,595,602,609,616,623,630,637,644,651,658,665,672,679,686,693,700,707,714,721,728,735,742,749,756,763,770,777,784,791,798,805,812,819,826,833,840,847,854,861,868,875,882,889,896,903,910,917,924,931,938,945,952,959,966,973,980,987,994,1001,1008,1015,1022,1029,1036,1043,1050,1057,1064,1071,1078,1085,1092,1099,1106,1113,1120,1127,1134,1141,1148,1155,1162,1169,1176,1183,1190,1197,1204,1211,1218,1225,1232,1239,1246,1253,1260,1267,1274,1281,1288,1295,1302,1309,1316,1323,1330,1337,1344,1351,1358,1365,1372,1379,1386,1393,1400,1407,1414,1421,1428,1435,1442,1449,1456,1463,1470,1477,1484,1491,1498,1505,1512,1519,1526,1533,1540,1547,1554,1561,1568,1575,1582,1589,1596,1603,1610,1617,1624,1631,1638,1645,1652,1659,1666,1673,1680,1687,1694,1701,1708,1715,1722,1729,1736,1743,1750,1757,1764,1771,1778,1785,1792,1799,1806,1813,1820,1827,1834,1841,1848,1855,1862,1869,1876,1883,1890,1897,1901,1905,1909,1913,1917,1921,1925,1929,1933,1937,1941,1945,1949,1953,1957,1961,1965,1969,1973,1977,1981,1985,1989,1993,1997,2001,2005,2009,2013,2017,2021,2025,2029,2033,2037,2041,2045,2049,2053,2057,2061,2065,2069,2073,2077,2081,2085,2089,2093,2097,2101,2105,2109,2113,2117,2121,2125,2129,2133,2137};
  //
  MEDCouplingUMesh *m=MEDCouplingUMesh::New("convexhull",2);
  m->allocateCells(331);
  for(int i=0;i<331;i++)
    m->insertNextCell(INTERP_KERNEL::NORM_POLYGON,connI[i+1]-connI[i],conn+connI[i]);
  m->finishInsertingCells();
  DataArrayDouble *coordsDa=DataArrayDouble::New();
  coordsDa->alloc(331,2);
  std::copy(coords,coords+662,coordsDa->getPointer());
  m->setCoords(coordsDa);
  coordsDa->decrRef();
  m->checkConsistencyLight();
  //
  DataArrayInt *da=m->convexEnvelop2D();
  m->checkConsistencyLight();
  CPPUNIT_ASSERT(coordsDa==m->getCoords());
  DataArrayInt *daC=da->buildComplement(331);
  da->decrRef();
  const int expected[58]={271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,302,303,304,305,306,307,308,309,310,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330};
  DataArrayInt *expected2=DataArrayInt::New();
  expected2->alloc(58,1);
  std::copy(expected,expected+58,expected2->getPointer());
  CPPUNIT_ASSERT(expected2->isEqual(*daC));
  //
  expected2->decrRef();
  daC->decrRef();
  //
  MEDCouplingFieldDouble *valsF=m->getMeasureField(ON_CELLS);
  DataArrayDouble *vals=valsF->getArray();
  const double ref[331]={184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,184.69493088478035,-61.564976961404426,-92.34746544254946,-92.34746544259811,-92.34746544253488,-92.3474654425349,-92.34746544180479,-92.34746544253493,-92.3474654419026,-92.34746544190256,-92.34746544253491,-61.564976961404426,-92.34746544254946,-92.34746544259811,-92.34746544253488,-92.3474654425349,-92.34746544180479,-92.34746544253493,-92.3474654419026,-92.34746544190256,-92.34746544253491,-61.564976961404426,-92.34746544254946,-92.34746544259811,-92.34746544253488,-92.3474654425349,-92.34746544180479,-92.34746544253493,-92.3474654419026,-92.34746544190256,-92.34746544253491,61.564976961404426,-92.34746544254946,-92.34746544259811,-92.34746544253488,-92.3474654425349,-92.34746544180479,-92.34746544253493,-92.3474654419026,-92.34746544190256,-92.34746544253491,61.564976961404426,-92.34746544254946,-92.34746544259811,-92.34746544253488,-92.3474654425349,-92.34746544180479,-92.34746544253493,-92.3474654419026,-92.34746544190256,-92.34746544253491,-61.564976961404426,-92.34746544254946,-92.34746544259811,-92.34746544253488,-92.3474654425349,-92.34746544180479,-92.34746544253493,-92.3474654419026,-92.34746544190256,-92.34746544253491};
  DataArrayDouble *ref2=DataArrayDouble::New(); ref2->alloc(331,1); std::copy(ref,ref+331,ref2->getPointer());
  vals->substractEqual(ref2);
  ref2->decrRef();
  vals->abs();
  DataArrayInt *theTest=vals->findIdsInRange(-1.,1e-7);
  CPPUNIT_ASSERT(theTest->isIota(331));
  theTest->decrRef();
  valsF->decrRef();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest5::testDataArraySort1()
{
  DataArrayInt *arr=DataArrayInt::New();
  CPPUNIT_ASSERT_THROW(arr->sort(true),INTERP_KERNEL::Exception);//no allocation
  CPPUNIT_ASSERT_THROW(arr->sort(false),INTERP_KERNEL::Exception);//no allocation
  const int values[6]={2,1,6,5,4,7};
  arr->alloc(3,2);
  CPPUNIT_ASSERT_THROW(arr->sort(true),INTERP_KERNEL::Exception);//no one component
  CPPUNIT_ASSERT_THROW(arr->sort(false),INTERP_KERNEL::Exception);//no one component
  arr->rearrange(1);
  std::copy(values,values+6,arr->getPointer());
  DataArrayInt *arr1=arr->deepCopy();
  DataArrayInt *arr2=arr->deepCopy();
  arr1->sort(true);
  const int expected1[6]={1,2,4,5,6,7};
 CPPUNIT_ASSERT_EQUAL(6,(int)arr1->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr1->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+6,arr1->begin()));
  arr2->sort(false);
  const int expected2[6]={7,6,5,4,2,1};
 CPPUNIT_ASSERT_EQUAL(6,(int)arr2->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)arr2->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected2,expected2+6,arr2->begin()));
  arr1->decrRef();
  arr2->decrRef();
  arr->decrRef();
  //
  DataArrayDouble *ard=DataArrayDouble::New();
  CPPUNIT_ASSERT_THROW(ard->sort(true),INTERP_KERNEL::Exception);//no allocation
  CPPUNIT_ASSERT_THROW(ard->sort(false),INTERP_KERNEL::Exception);//no allocation
  const double valuesD[6]={2.,1.,6.,5.,4.,7.};
  ard->alloc(3,2);
  CPPUNIT_ASSERT_THROW(ard->sort(true),INTERP_KERNEL::Exception);//no one component
  CPPUNIT_ASSERT_THROW(ard->sort(false),INTERP_KERNEL::Exception);//no one component
  ard->rearrange(1);
  std::copy(valuesD,valuesD+6,ard->getPointer());
  DataArrayDouble *ard1=ard->deepCopy();
  DataArrayDouble *ard2=ard->deepCopy();
  ard1->sort(true);
  const double expected3[6]={1.,2.,4.,5.,6.,7.};
  CPPUNIT_ASSERT_EQUAL(6,(int)ard1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)ard1->getNumberOfComponents());
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],ard1->getIJ(i,0),1e-12);
  ard2->sort(false);
  const double expected4[6]={7.,6.,5.,4.,2.,1.};
 CPPUNIT_ASSERT_EQUAL(6,(int)ard2->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)ard2->getNumberOfComponents());
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected4[i],ard2->getIJ(i,0),1e-12);
  ard1->decrRef();
  ard2->decrRef();
  ard->decrRef();
}

void MEDCouplingBasicsTest5::testPartitionBySpreadZone1()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  const int part0[3]={2,3,4};
  const int part1[2]={0,1};
  MEDCouplingUMesh *m1=static_cast<MEDCouplingUMesh *>(m->buildPartOfMySelf(part0,part0+3));
  MEDCouplingUMesh *m2=static_cast<MEDCouplingUMesh *>(m->buildPartOfMySelf(part1,part1+2));
  std::vector<const MEDCouplingUMesh *> v(3); v[0]=m; v[1]=m1; v[2]=m2;
  MEDCouplingUMesh *m4=MEDCouplingUMesh::MergeUMeshes(v);
  const int renum[10]={5,2,9,6,4,7,0,1,3,8};
  m4->renumberCells(renum);
  //
  std::vector<DataArrayInt *> v2=m4->partitionBySpreadZone();
  CPPUNIT_ASSERT_EQUAL(3,(int)v2.size());
  const int expected0[3]={0,1,7};
  const int expected1[5]={2,4,5,6,9};
  const int expected2[2]={3,8};
 CPPUNIT_ASSERT_EQUAL(3,(int)v2[0]->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)v2[0]->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(5,(int)v2[1]->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)v2[1]->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(2,(int)v2[2]->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)v2[2]->getNumberOfComponents());
  //
  CPPUNIT_ASSERT(std::equal(expected0,expected0+3,v2[0]->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected1,expected1+5,v2[1]->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+2,v2[2]->getConstPointer()));
  v2[0]->decrRef();
  v2[1]->decrRef();
  v2[2]->decrRef();
  //
  MEDCouplingUMesh *m5=m4->buildSpreadZonesWithPoly();
  CPPUNIT_ASSERT_EQUAL(3,(int)m5->getNumberOfCells());
  CPPUNIT_ASSERT(m5->getCoords()==m4->getCoords());
  const int expected3[23]={5,15,16,17,14,11,13,12,5,2,1,0,3,6,7,8,5,5,18,21,22,20,19};
  const int expected4[4]={0,8,17,23};
  CPPUNIT_ASSERT_EQUAL(23,(int)m5->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+23,m5->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(4,(int)m5->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected4,expected4+4,m5->getNodalConnectivityIndex()->getConstPointer()));
  //
  m->decrRef();
  m1->decrRef();
  m2->decrRef();
  m4->decrRef();
  m5->decrRef();
}

void MEDCouplingBasicsTest5::testGiveCellsWithType1()
{
  const int expected0[2]={1,2};
  const int expected1[3]={0,3,4};
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  DataArrayInt *da=m->giveCellsWithType(INTERP_KERNEL::NORM_TRI3);
 CPPUNIT_ASSERT_EQUAL(2,(int)da->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected0,expected0+2,da->getConstPointer()));
  da->decrRef();
  //
  da=m->giveCellsWithType(INTERP_KERNEL::NORM_QUAD4);
 CPPUNIT_ASSERT_EQUAL(3,(int)da->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,da->getConstPointer()));
  da->decrRef();
  //
  da=m->giveCellsWithType(INTERP_KERNEL::NORM_TRI6);
 CPPUNIT_ASSERT_EQUAL(0,(int)da->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da->getNumberOfComponents());
  da->decrRef();
  //
  CPPUNIT_ASSERT_THROW(m->giveCellsWithType(INTERP_KERNEL::NORM_SEG2),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(m->giveCellsWithType(INTERP_KERNEL::NORM_HEXA8),INTERP_KERNEL::Exception);
  //
  m->decrRef();
}

void MEDCouplingBasicsTest5::testBuildSlice3D2()
{
  MEDCouplingUMesh *mesh2D=0;
  MEDCouplingUMesh *mesh3D=build3DExtrudedUMesh_1(mesh2D);
  mesh2D->decrRef();
  // First slice in the middle of 3D cells
  const double vec1[3]={-0.07,1.,0.07};
  const double origin1[3]={1.524,1.4552,1.74768};
  DataArrayInt *ids=0;
  MEDCouplingUMesh *slice1=mesh3D->buildSlice3D(origin1,vec1,1e-10,ids);
  //
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f->setTime(4.5,6,7) ; f->setMesh(mesh3D);
  DataArrayDouble *arr=DataArrayDouble::New(); arr->alloc(mesh3D->getNumberOfCells(),2);
  arr->rearrange(1); arr->iota(2.); arr->rearrange(2);
  f->setArray(arr);
  f->checkConsistencyLight();
  const int exp1[9]={1,3,4,7,9,10,13,15,16};
  DataArrayInt *expected1=DataArrayInt::New(); expected1->alloc(9,1); std::copy(exp1,exp1+9,expected1->getPointer());
  CPPUNIT_ASSERT(expected1->isEqual(*ids));
  DataArrayDouble *arr2=arr->selectByTupleIdSafe(expected1->begin(),expected1->end());
  //
  MEDCouplingFieldDouble *f2=f->extractSlice3D(origin1,vec1,1e-10);
  CPPUNIT_ASSERT(f2->getArray()->isEqual(*arr2,1e-12));
  CPPUNIT_ASSERT(slice1->isEqual(f2->getMesh(),1e-12));
  int a,b;
  double c=f2->getTime(a,b);
  CPPUNIT_ASSERT_EQUAL(6,a);
  CPPUNIT_ASSERT_EQUAL(7,b);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.5,c,1e-12);
  //
  ids->decrRef();
  slice1->decrRef();
  arr2->decrRef();
  arr->decrRef();
  f2->decrRef();
  f->decrRef();
  mesh3D->decrRef();
  expected1->decrRef();
}

void MEDCouplingBasicsTest5::testComputeTupleIdsToSelectFromCellIds1()
{
  MEDCouplingUMesh *m=build2DTargetMesh_3();
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_GAUSS_NE,NO_TIME);
  f->setMesh(m);
  DataArrayDouble *arr=DataArrayDouble::New(); arr->alloc(52,2) ; arr->rearrange(1) ; arr->iota(7.); arr->rearrange(2);
  f->setArray(arr);
  //
  const int subPart1[3]={1,5,9};
  MEDCouplingFieldDouble *f2=f->buildSubPart(subPart1,subPart1+3);
  f2->checkConsistencyLight();
  DataArrayInt *cI=m->computeNbOfNodesPerCell();
  cI->computeOffsetsFull();
  const int sel1[3]={1,5,9};
  DataArrayInt *sel=DataArrayInt::New(); sel->useArray(sel1,false,CPP_DEALLOC,3,1);
  DataArrayInt *res=sel->buildExplicitArrByRanges(cI);
  DataArrayDouble *arr2=arr->selectByTupleIdSafe(res->begin(),res->end());
  const double expected1[30]={13.,14.,15.,16.,17.,18.,19.,20.,59.,60.,61.,62.,63.,64.,95.,96.,97.,98.,99.,100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,110.};
  DataArrayDouble *arr3=DataArrayDouble::New(); arr3->useArray(expected1,false,CPP_DEALLOC,15,2);
  CPPUNIT_ASSERT(arr2->isEqual(*arr3,1e-12));
  CPPUNIT_ASSERT(arr2->isEqual(*f2->getArray(),1e-12));
  //
  cI->decrRef();
  arr3->decrRef();
  arr2->decrRef();
  arr->decrRef();
  m->decrRef();
  f->decrRef();
  f2->decrRef();
  sel->decrRef();
  res->decrRef();
}

void MEDCouplingBasicsTest5::testComputeSkin1()
{
  const double input1[5]={2.,3.4,5.6,7.7,8.0};
  const double input2[6]={2.,3.4,5.6,7.7,9.0,14.2};
  DataArrayDouble *arrX=DataArrayDouble::New(); arrX->alloc(5,1); std::copy(input1,input1+5,arrX->getPointer());
  DataArrayDouble *arrY=DataArrayDouble::New(); arrY->alloc(6,1); std::copy(input2,input2+6,arrY->getPointer());
  MEDCouplingCMesh *cmesh=MEDCouplingCMesh::New() ; cmesh->setCoordsAt(0,arrX) ; cmesh->setCoordsAt(1,arrY);
  MEDCouplingUMesh *umesh=cmesh->buildUnstructured();
  cmesh->decrRef(); arrX->decrRef(); arrY->decrRef();
  //
  MEDCouplingUMesh *skin=umesh->computeSkin();
  CPPUNIT_ASSERT_EQUAL(18,(int)skin->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(1,skin->getMeshDimension());
  CPPUNIT_ASSERT(skin->getCoords()==umesh->getCoords());
  const int expected1[19]={0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54};
  const int expected2[54]={1,1,0,1,0,5,1,2,1,1,3,2,1,4,3,1,9,4,1,5,10,1,14,9,1,10,15,1,19,14,1,15,20,1,24,19,1,20,25,1,25,26,1,26,27,1,27,28,1,28,29,1,29,24};
  CPPUNIT_ASSERT_EQUAL((std::size_t)19,skin->getNodalConnectivityIndex()->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+19,skin->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL((std::size_t)54,skin->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(expected2,expected2+54,skin->getNodalConnectivity()->getConstPointer()));
  DataArrayInt *ids=skin->computeFetchedNodeIds();
  const int expected3[18]={0,1,2,3,4,5,9,10,14,15,19,20,24,25,26,27,28,29};
  CPPUNIT_ASSERT_EQUAL((std::size_t)18,ids->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+18,ids->getConstPointer()));
  MEDCouplingUMesh *part=dynamic_cast<MEDCouplingUMesh *>(umesh->buildFacePartOfMySelfNode(ids->begin(),ids->end(),true));
  part->setName(skin->getName().c_str());
  CPPUNIT_ASSERT(part->isEqual(skin,1e-12));
  MEDCouplingUMesh *part2=dynamic_cast<MEDCouplingUMesh *>(part->buildPartOfMySelfSlice(1,18,2,true));
  DataArrayInt *ids2=DataArrayInt::Range(0,18,2);
  part->setPartOfMySelf(ids2->begin(),ids2->end(),*part2);
  ids2->decrRef();
  CPPUNIT_ASSERT(!part->isEqual(skin,1e-12));
  DataArrayInt *trad=part->zipConnectivityTraducer(0);
  CPPUNIT_ASSERT_EQUAL(9,(int)part->getNumberOfCells());
  const int expected4[18]={0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8};
  CPPUNIT_ASSERT(std::equal(expected4,expected4+18,trad->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL((std::size_t)18,trad->getNbOfElems());
  trad->decrRef();
  part->decrRef();
  part2->decrRef();
  //
  ids->decrRef();
  umesh->decrRef();
  skin->decrRef();
}

void MEDCouplingBasicsTest5::testUMeshSetPartOfMySelf2()
{
  // resize with explicit ids list
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  std::set<INTERP_KERNEL::NormalizedCellType> s; s.insert(INTERP_KERNEL::NORM_TRI3); s.insert(INTERP_KERNEL::NORM_QUAD4);
  CPPUNIT_ASSERT(s==m->getAllGeoTypes());
  const int ids1[3]={0,3,4};
  MEDCouplingUMesh *part=static_cast<MEDCouplingUMesh *>(m->buildPartOfMySelf(ids1,ids1+3,true));
  part->simplexize(0)->decrRef();
  const int ids2[3]={1,2,5};
  MEDCouplingUMesh *part2=static_cast<MEDCouplingUMesh *>(part->buildPartOfMySelf(ids2,ids2+3,true));
  m->setPartOfMySelf(ids1,ids1+3,*part2);
  const int expected1[20]={3,0,4,1,3,1,4,2,3,4,5,2,3,6,7,4,3,7,5,4};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+20,m->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL((std::size_t)20,m->getNodalConnectivity()->getNbOfElems());
  const int expected2[6]={0,4,8,12,16,20};
  CPPUNIT_ASSERT(std::equal(expected2,expected2+6,m->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL((std::size_t)6,m->getNodalConnectivityIndex()->getNbOfElems());
  s.clear(); s.insert(INTERP_KERNEL::NORM_TRI3);
  CPPUNIT_ASSERT(s==m->getAllGeoTypes());
  m->decrRef(); part->decrRef(); part2->decrRef();
  // no resize with explicit ids list
  m=build2DTargetMesh_1();
  part=static_cast<MEDCouplingUMesh *>(m->buildPartOfMySelf(ids1,ids1+2,true));
  part->convertAllToPoly();
  m->setPartOfMySelf(ids1+1,ids1+3,*part);
  const int expected3[23]={4,0,3,4,1,3,1,4,2,3,4,5,2,5,0,3,4,1,5,6,7,4,3};
  CPPUNIT_ASSERT(std::equal(expected3,expected3+23,m->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL((std::size_t)23,m->getNodalConnectivity()->getNbOfElems());
  const int expected4[6]={0,5,9,13,18,23};
  CPPUNIT_ASSERT(std::equal(expected4,expected4+6,m->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL((std::size_t)6,m->getNodalConnectivityIndex()->getNbOfElems());
  s.clear(); s.insert(INTERP_KERNEL::NORM_TRI3); s.insert(INTERP_KERNEL::NORM_QUAD4); s.insert(INTERP_KERNEL::NORM_POLYGON);
  CPPUNIT_ASSERT(s==m->getAllGeoTypes());
  m->decrRef(); part->decrRef();
  // resize with range ids
  m=build2DTargetMesh_1();
  part=static_cast<MEDCouplingUMesh *>(m->buildPartOfMySelfSlice(3,5,1,true));
  m->setPartOfMySelfSlice(1,3,1,*part);
  const int expected5[25]={4,0,3,4,1,4,6,7,4,3,4,7,8,5,4,4,6,7,4,3,4,7,8,5,4};
  CPPUNIT_ASSERT(std::equal(expected5,expected5+25,m->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL((std::size_t)25,m->getNodalConnectivity()->getNbOfElems());
  const int expected6[6]={0,5,10,15,20,25};
  CPPUNIT_ASSERT(std::equal(expected6,expected6+6,m->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL((std::size_t)6,m->getNodalConnectivityIndex()->getNbOfElems());
  s.clear(); s.insert(INTERP_KERNEL::NORM_QUAD4);
  CPPUNIT_ASSERT(s==m->getAllGeoTypes());
  m->decrRef(); part->decrRef();
  // no resize with range ids
  m=build2DTargetMesh_1();
  part=static_cast<MEDCouplingUMesh *>(m->buildPartOfMySelfSlice(0,5,3,true));
  part->convertAllToPoly();
  m->setPartOfMySelfSlice(3,5,1,*part);
  const int expected7[23]={4,0,3,4,1,3,1,4,2,3,4,5,2,5,0,3,4,1,5,6,7,4,3};
  CPPUNIT_ASSERT(std::equal(expected7,expected7+23,m->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL((std::size_t)23,m->getNodalConnectivity()->getNbOfElems());
  const int expected8[6]={0,5,9,13,18,23};
  CPPUNIT_ASSERT(std::equal(expected8,expected8+6,m->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL((std::size_t)6,m->getNodalConnectivityIndex()->getNbOfElems());
s.clear(); s.insert(INTERP_KERNEL::NORM_TRI3); s.insert(INTERP_KERNEL::NORM_QUAD4); s.insert(INTERP_KERNEL::NORM_POLYGON);
  CPPUNIT_ASSERT(s==m->getAllGeoTypes());
  m->decrRef(); part->decrRef();
  // no resize with range ids negative direction
  m=build2DTargetMesh_1();
  part=static_cast<MEDCouplingUMesh *>(m->buildPartOfMySelfSlice(3,-1,-3,true));
  part->convertAllToPoly();
  m->setPartOfMySelfSlice(4,2,-1,*part);
  const int expected9[23]={4,0,3,4,1,3,1,4,2,3,4,5,2,5,0,3,4,1,5,6,7,4,3};
  CPPUNIT_ASSERT(std::equal(expected9,expected9+23,m->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL((std::size_t)23,m->getNodalConnectivity()->getNbOfElems());
  const int expected10[6]={0,5,9,13,18,23};
  CPPUNIT_ASSERT(std::equal(expected10,expected10+6,m->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL((std::size_t)6,m->getNodalConnectivityIndex()->getNbOfElems());
  s.clear(); s.insert(INTERP_KERNEL::NORM_TRI3); s.insert(INTERP_KERNEL::NORM_QUAD4); s.insert(INTERP_KERNEL::NORM_POLYGON);
  CPPUNIT_ASSERT(s==m->getAllGeoTypes());
  part->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest5::testUnPolyze3()
{
  const double coord[18]={0.0,0.5,-0.5,-0.5,-0.5,-0.5,0.5,-0.5,-0.5,0.0,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5,0.5};
  const int conn[22]={1,2,5,4,-1,4,3,0,1,-1,2,0,3,5,-1,0,2,1,-1,4,5,3};
  MEDCouplingUMesh *m=MEDCouplingUMesh::New("a mesh",3);
  m->allocateCells(1);
  m->insertNextCell(INTERP_KERNEL::NORM_POLYHED,22,conn);
  m->finishInsertingCells();
  DataArrayDouble *coords=DataArrayDouble::New();
  coords->alloc(6,3);
  std::copy(coord,coord+18,coords->getPointer());
  m->setCoords(coords);
  coords->decrRef();
  m->checkConsistencyLight();
  //
  MEDCouplingFieldDouble *vol=m->getMeasureField(ON_CELLS);
  CPPUNIT_ASSERT_EQUAL(1,(int)vol->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,vol->getArray()->getIJ(0,0),1e-12);
  vol->decrRef();
  //
  m->unPolyze();
  CPPUNIT_ASSERT_EQUAL(1,(int)m->getNumberOfCells());
  std::set<INTERP_KERNEL::NormalizedCellType> s; s.insert(INTERP_KERNEL::NORM_PENTA6);
  CPPUNIT_ASSERT(s==m->getAllGeoTypes());
  //
  const int expected1[2]={0,7};
  const int expected2[7]={16,0,2,1,3,5,4};
  CPPUNIT_ASSERT_EQUAL(2,(int)m->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+2,m->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(7,(int)m->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected2,expected2+7,m->getNodalConnectivity()->getConstPointer()));
  //
  vol=m->getMeasureField(ON_CELLS);
  CPPUNIT_ASSERT_EQUAL(1,(int)vol->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,vol->getArray()->getIJ(0,0),1e-12);
  vol->decrRef();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest5::testKrSpatialDiscretization1()
{
  const double srcPointCoordsX[10]={0.8401877171547095, 0.7830992237586059, 0.9116473579367843, 0.335222755714889, 0.2777747108031878, 0.4773970518621602, 0.3647844727918433, 0.9522297251747128, 0.6357117279599009, 0.1416025553558034};
  const double srcFieldValsOnPoints[10]={2.129892434968836, 2.295320474540621, 1.931948594981134, 2.728013590937196, 2.715603240418478, 2.661778472822935, 2.695696990104364, 1.893710234970982, 2.529628016549284, 2.728432341300668};
  const double targetPointCoordsX[40]={-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,-6.93889390391e-17,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45};
  const double targetFieldValsExpected[40]={2.975379475824351, 2.95613491917003, 2.936890362515361, 2.917645805861018, 2.898401249206574, 2.879156692552137, 2.859912135897732, 2.840667579243201, 2.821423022588731, 2.802178465934342, 2.78293390927989, 2.763689352625457, 2.744444795971001, 2.725209522098197, 2.709077577124666, 2.706677252549218, 2.727467797847971, 2.713338094723676, 2.671342424824244, 2.664877370146978, 2.653840141412181, 2.619607861392791, 2.569777214476479, 2.513263929794591, 2.450732752808528, 2.368313560985155, 2.250909795670307, 2.098194272085416, 1.954257891732065, 1.895040660973802, 1.865256788315972, 1.835475248687992, 1.80569370905998, 1.775912169431971, 1.746130629803976, 1.716349090175918, 1.686567550547855, 1.656786010919941, 1.627004471291988, 1.597222931663817};
  //
  int nbOfInputPoints=10;
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_NODES_KR,ONE_TIME);
  DataArrayDouble *srcArrX=DataArrayDouble::New();
  srcArrX->alloc(nbOfInputPoints,1);
  std::copy(srcPointCoordsX,srcPointCoordsX+nbOfInputPoints,srcArrX->getPointer());
  MEDCouplingCMesh *cmesh=MEDCouplingCMesh::New("aMesh");
  cmesh->setCoordsAt(0,srcArrX);
  MEDCouplingUMesh *umesh=cmesh->buildUnstructured();
  f->setMesh(umesh);
  DataArrayDouble *srcVals=DataArrayDouble::New();
  srcVals->alloc(nbOfInputPoints,1);
  std::copy(srcFieldValsOnPoints,srcFieldValsOnPoints+nbOfInputPoints,srcVals->getPointer());
  f->setArray(srcVals);
  f->checkConsistencyLight();
  //
  double *res0=new double[1];
  f->getValueOn(targetPointCoordsX,res0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(targetFieldValsExpected[0],res0[0],1e-10);
  delete [] res0;
  //
  DataArrayDouble *valuesToTest=f->getValueOnMulti(targetPointCoordsX,40);
 CPPUNIT_ASSERT_EQUAL(40,(int)valuesToTest->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)valuesToTest->getNumberOfComponents());
  for(int i=0;i<40;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(targetFieldValsExpected[i],valuesToTest->getIJ(i,0),1e-10);
  valuesToTest->decrRef();
  //
  cmesh->decrRef();
  umesh->decrRef();
  srcArrX->decrRef();
  srcVals->decrRef();
  f->decrRef();
}

void MEDCouplingBasicsTest5::testDuplicateEachTupleNTimes1()
{
  const double vals0[4]={9.,8.,7.,6.};
  DataArrayDouble *d=DataArrayDouble::New(); d->useArray(vals0,false,CPP_DEALLOC,4,1); d->setInfoOnComponent(0,"mass [kg]"); d->setName("aname");
  DataArrayDouble *d2=d->duplicateEachTupleNTimes(3);
  const double vals1[12]={9.,9.,9.,8.,8.,8.,7.,7.,7.,6.,6.,6.};
  DataArrayDouble *d3=DataArrayDouble::New(); d3->useArray(vals1,false,CPP_DEALLOC,4*3,1); d3->setName("aname"); d3->setInfoOnComponent(0,"mass [kg]");
  CPPUNIT_ASSERT(d2->isEqual(*d2,1e-14)); d3->decrRef();
  d->decrRef();
  d2->decrRef();
  //
  const int vals2[4]={9,8,7,6};
  DataArrayInt *d4=DataArrayInt::New(); d4->useArray(vals2,false,CPP_DEALLOC,4,1); d4->setInfoOnComponent(0,"mass [kg]") ; d4->setName("aname");
  DataArrayInt *d5=d4->duplicateEachTupleNTimes(3);
  const int vals3[12]={9,9,9,8,8,8,7,7,7,6,6,6};
  DataArrayInt *d6=DataArrayInt::New(); d6->useArray(vals3,false,CPP_DEALLOC,4*3,1); d6->setName("aname"); d6->setInfoOnComponent(0,"mass [kg]");
  CPPUNIT_ASSERT(d5->isEqual(*d6)); d6->decrRef();
  d4->decrRef();
  d5->decrRef();
}

void MEDCouplingBasicsTest5::testIntersect2DMeshesTmp5()
{
  // coordinates
  DataArrayDouble *coords=DataArrayDouble::New();
  const double coordsData[376]={41,0,42,0,0,42,0,41,41.5,0,29.698484809834998,29.698484809834994,0,41.5,28.991378028648452,28.991378028648445,-42,0,-41,0,-29.698484809834994,29.698484809834998,-41.5,0,-28.991378028648445,28.991378028648452,0,-42,0,-41,-29.698484809835001,-29.698484809834994,0,-41.5,-28.991378028648455,-28.991378028648445,29.698484809834987,-29.698484809835001,28.991378028648441,-28.991378028648455,43,0,0,43,42.5,0,30.405591591021544,30.40559159102154,0,42.5,-43,0,-30.40559159102154,30.405591591021544,-42.5,0,0,-43,-30.405591591021551,-30.40559159102154,0,-42.5,30.405591591021537,-30.405591591021551,44,0,0,44,43.5,0,31.112698372208094,31.112698372208087,0,43.5,-44,0,-31.112698372208087,31.112698372208094,-43.5,0,0,-44,-31.112698372208097,-31.112698372208087,0,-43.5,31.112698372208083,-31.112698372208097,45,0,0,45,44.5,0,31.81980515339464,31.819805153394636,0,44.5,-45,0,-31.819805153394636,31.81980515339464,-44.5,0,0,-45,-31.819805153394647,-31.819805153394636,0,-44.5,31.819805153394629,-31.819805153394647,47,0,0,47,46,0,33.234018715767739,33.234018715767732,0,46,-47,0,-33.234018715767732,33.234018715767739,-46,0,0,-47,-33.234018715767739,-33.234018715767732,0,-46,33.234018715767725,-33.234018715767739,49,0,0,49,48,0,34.648232278140831,34.648232278140824,0,48,-49,0,-34.648232278140824,34.648232278140831,-48,0,0,-49,-34.648232278140839,-34.648232278140824,0,-48,34.648232278140817,-34.648232278140839,51,0,0,51,50,0,36.062445840513924,36.062445840513924,0,50,-51,0,-36.062445840513924,36.062445840513924,-50,0,0,-51,-36.062445840513931,-36.062445840513924,0,-50,36.062445840513917,-36.062445840513931,53,0,0,53,52,0,37.476659402887023,37.476659402887016,0,52,-53,0,-37.476659402887016,37.476659402887023,-52,0,0,-53,-37.47665940288703,-37.476659402887016,0,-52,37.476659402887009,-37.47665940288703,55,0,0,55,54,0,38.890872965260115,38.890872965260108,0,54,-55,0,-38.890872965260108,38.890872965260115,-54,0,0,-55,-38.890872965260122,-38.890872965260108,0,-54,38.890872965260101,-38.890872965260122,59,0,0,59,57,0,41.719300090006307,41.7193000900063,0,57,-59,0,-41.7193000900063,41.719300090006307,-57,0,0,-59,-41.719300090006314,-41.7193000900063,0,-57,41.719300090006293,-41.719300090006314,63,0,0,63,61,0,44.547727214752499,44.547727214752491,0,61,-63,0,-44.547727214752491,44.547727214752499,-61,0,0,-63,-44.547727214752506,-44.547727214752491,0,-61,44.547727214752484,-44.547727214752506,67,0,0,67,65,0,47.37615433949869,47.376154339498683,0,65,-67,0,-47.376154339498683,47.37615433949869,-65,0,0,-67,-47.376154339498697,-47.376154339498683,0,-65,47.376154339498676,-47.376154339498697,71,0,0,71,69,0,50.204581464244875,50.204581464244868,0,69,-71,0,-50.204581464244868,50.204581464244875,-69,0,0,-71,-50.204581464244889,-50.204581464244868,0,-69,50.20458146424486,-50.204581464244889,75,0,0,75,73,0,53.033008588991066,53.033008588991059,0,73,-75,0,-53.033008588991059,53.033008588991066,-73,0,0,-75,-53.033008588991073,-53.033008588991059,0,-73,53.033008588991052,-53.033008588991073,80,0,0,80,77.5,0,56.568542494923804,56.568542494923797,0,77.5,-80,0,-56.568542494923797,56.568542494923804,-77.5,0,0,-80,-56.568542494923818,-56.568542494923797,0,-77.5,56.56854249492379,-56.568542494923818};
  coords->useArray(coordsData,false,CPP_DEALLOC,188,2);
  coords->setName("");
  DataArrayInt *conn=DataArrayInt::New();
  const int connData[540]={8,0,1,2,3,4,5,6,7,8,3,2,8,9,6,10,11,12,8,9,8,13,14,11,15,16,17,8,14,13,1,0,16,18,4,19,8,1,20,21,2,22,23,24,5,8,2,21,25,8,24,26,27,10,8,8,25,28,13,27,29,30,15,8,13,28,20,1,30,31,22,18,8,20,32,33,21,34,35,36,23,8,21,33,37,25,36,38,39,26,8,25,37,40,28,39,41,42,29,8,28,40,32,20,42,43,34,31,8,32,44,45,33,46,47,48,35,8,33,45,49,37,48,50,51,38,8,37,49,52,40,51,53,54,41,8,40,52,44,32,54,55,46,43,8,44,56,57,45,58,59,60,47,8,45,57,61,49,60,62,63,50,8,49,61,64,52,63,65,66,53,8,52,64,56,44,66,67,58,55,8,56,68,69,57,70,71,72,59,8,57,69,73,61,72,74,75,62,8,61,73,76,64,75,77,78,65,8,64,76,68,56,78,79,70,67,8,68,80,81,69,82,83,84,71,8,69,81,85,73,84,86,87,74,8,73,85,88,76,87,89,90,77,8,76,88,80,68,90,91,82,79,8,80,92,93,81,94,95,96,83,8,81,93,97,85,96,98,99,86,8,85,97,100,88,99,101,102,89,8,88,100,92,80,102,103,94,91,8,92,104,105,93,106,107,108,95,8,93,105,109,97,108,110,111,98,8,97,109,112,100,111,113,114,101,8,100,112,104,92,114,115,106,103,8,104,116,117,105,118,119,120,107,8,105,117,121,109,120,122,123,110,8,109,121,124,112,123,125,126,113,8,112,124,116,104,126,127,118,115,8,116,128,129,117,130,131,132,119,8,117,129,133,121,132,134,135,122,8,121,133,136,124,135,137,138,125,8,124,136,128,116,138,139,130,127,8,128,140,141,129,142,143,144,131,8,129,141,145,133,144,146,147,134,8,133,145,148,136,147,149,150,137,8,136,148,140,128,150,151,142,139,8,140,152,153,141,154,155,156,143,8,141,153,157,145,156,158,159,146,8,145,157,160,148,159,161,162,149,8,148,160,152,140,162,163,154,151,8,152,164,165,153,166,167,168,155,8,153,165,169,157,168,170,171,158,8,157,169,172,160,171,173,174,161,8,160,172,164,152,174,175,166,163,8,164,176,177,165,178,179,180,167,8,165,177,181,169,180,182,183,170,8,169,181,184,172,183,185,186,173,8,172,184,176,164,186,187,178,175};
  conn->useArray(connData,false,CPP_DEALLOC,540,1);
  conn->setName("");
  DataArrayInt *connI=DataArrayInt::New();
  const int connIData[61]={0,9,18,27,36,45,54,63,72,81,90,99,108,117,126,135,144,153,162,171,180,189,198,207,216,225,234,243,252,261,270,279,288,297,306,315,324,333,342,351,360,369,378,387,396,405,414,423,432,441,450,459,468,477,486,495,504,513,522,531,540};
  connI->useArray(connIData,false,CPP_DEALLOC,61,1);
  connI->setName("");
  //
  MEDCouplingUMesh *m1=MEDCouplingUMesh::New("Fix",2);
  m1->setCoords(coords);
  m1->setConnectivity(conn,connI,true);
  coords->decrRef(); conn->decrRef(); connI->decrRef();
  //
  coords=DataArrayDouble::New();
  const double coordsData2[84]={46.5,-2.5,53.5,-2.5,53.5,2.5,46.5,2.5,50,-2.5,53.5,0,50,2.5,46.5,0,60.5,-2.5,60.5,2.5,57,-2.5,60.5,0,57,2.5,53.5,7.5,46.5,7.5,53.5,5,50,7.5,46.5,5,60.5,7.5,60.5,5,57,7.5,-2,47,2,47,2,53,-2,53,0,47,2,50,0,53,-2,50,6,47,6,53,4,47,6,50,4,53,2,59,-2,59,2,56,0,59,-2,56,6,59,6,56,4,59};
  coords->useArray(coordsData2,false,CPP_DEALLOC,42,2);  
  coords->setName("");
  // connectivity
  conn=DataArrayInt::New();
  const int connData2[72]={8,0,1,2,3,4,5,6,7,8,1,8,9,2,10,11,12,5,8,3,2,13,14,6,15,16,17,8,2,9,18,13,12,19,20,15,8,21,22,23,24,25,26,27,28,8,22,29,30,23,31,32,33,26,8,24,23,34,35,27,36,37,38,8,23,30,39,34,33,40,41,36};
  conn->useArray(connData2,false,CPP_DEALLOC,72,1);
  conn->setName("");
  connI=DataArrayInt::New();
  const int connIData2[9]={0,9,18,27,36,45,54,63,72};
  connI->useArray(connIData2,false,CPP_DEALLOC,9,1);
  connI->setName("");
  MEDCouplingUMesh *m2=MEDCouplingUMesh::New("Mobile",2);
  m2->setCoords(coords);
  m2->setConnectivity(conn,connI,true);
  coords->decrRef(); conn->decrRef(); connI->decrRef();
  //
  DataArrayInt *d1=0,*d2=0;
  MEDCouplingUMesh *m3=MEDCouplingUMesh::Intersect2DMeshes(m1,m2,1e-10,d1,d2);
  CPPUNIT_ASSERT_EQUAL(105,(int)m3->getNumberOfCells());
 CPPUNIT_ASSERT_EQUAL(105,(int)d1->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(105,(int)d2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(704,m3->getNumberOfNodes());
  //
  const double areaExpected[105]={-65.18804756198824,-65.18804756198824,-65.18804756198824,-65.18804756198824,-66.75884388878285,-66.75884388878285,-66.7588438887833,-66.75884388878308,-68.32964021557768,-68.32964021557768,-68.32964021557814,-68.32964021557791,-69.9004365423732,-69.9004365423732,-69.90043654237297,-69.90043654237297,-1.194568659706448,-1.0869994447159463,-142.2316939607081,-144.51326206513068,-144.5132620651309,-1.1945686597064424,-143.3186934054243,-5.002264310862817,-10.0261332846393,-3.9727823117092953,-7.290862524642649,-124.504404940456,-3.9727823117093237,-146.82366506060032,-150.79644737231024,-5.002264310862776,-145.79418306144626,-5.00208651738126,-10.054764051268958,-4.001067863263231,-8.027932154428669,-129.99378209314813,-4.001067863263216,-153.07856481622616,-157.0796326794898,-5.0020865173811915,-152.07754616210832,-5.001928880064381,-10.050590216368969,-4.00098721602491,-8.025810856794209,-136.28350081741684,-4.000987216024939,-159.36183077064402,-163.36281798667005,-5.0019288800643285,-158.36088910660442,-1.2991516319851801,-3.702636830195414,-3.7815130030068254,-6.265364371195623,-0.02516260900254963,-0.6553944641345026,-3.975752765070567,-7.368528340442765,-142.57249927881398,-0.02516260900254963,-3.9757527650706095,-165.64508791977525,-169.64600329384803,-1.299151631985167,-3.7026368301953885,-164.6442148316677,-10.00321285677458,-20.08414323176165,-8.001644468035863,-16.042954878437143,-304.0096070742277,-8.00164446803587,-350.1399180412005,-358.1415625092368,-10.003212856774468,-348.13834965246224,-3.794150313030109,-8.65049239704272,-0.02260276689354157,-0.5885167811200915,-370.2185414798688,-0.022602766893559393,-383.2517009710623,-383.2743037379555,-3.7941503130300576,-379.48015342492505,-408.40704496667513,-408.4070449666742,-408.4070449666742,-408.4070449666742,-433.53978619538975,-433.5397861953902,-433.5397861953911,-433.53978619539066,-458.67252742410983,-458.6725274241094,-458.67252742410983,-458.6725274241089,-608.6835766330232,-608.6835766330232,-608.6835766330232,-608.6835766330241};
  const int expected1[105]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16,16,17,18,19,19,20,20,20,20,20,21,21,22,23,23,24,24,24,24,24,25,25,26,27,27,28,28,28,28,28,29,29,30,31,31,32,32,32,32,32,32,32,32,32,33,33,33,34,35,35,35,36,36,36,36,36,37,37,38,39,39,40,40,40,40,40,41,41,42,43,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59};
  const int expected2[105]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,2,-1,-1,-1,0,-1,0,2,4,5,-1,4,-1,-1,0,-1,0,2,4,5,-1,4,-1,-1,0,-1,0,2,4,5,-1,4,-1,-1,0,-1,0,1,2,3,4,5,6,7,-1,4,6,-1,-1,0,1,-1,1,3,6,7,-1,6,-1,-1,1,-1,1,3,6,7,-1,6,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
  MEDCouplingFieldDouble *f3f=m3->getMeasureField(ON_CELLS);
  const double *f3=f3f->getArray()->getConstPointer();
  for(int i=0;i<105;i++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(areaExpected[i],f3[i],1e-10);
      CPPUNIT_ASSERT_EQUAL(expected1[i],d1->getIJ(i,0));
      CPPUNIT_ASSERT_EQUAL(expected2[i],d2->getIJ(i,0));
    }
  //
  f3f->decrRef();
  m3->decrRef();
  d1->decrRef();
  d2->decrRef();
  m2->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest5::testDAIBuildUnique1()
{
  DataArrayInt *d=DataArrayInt::New();
  const int dData[14]={1,2,2,3,3,3,3,4,5,5,7,7,7,19};
  d->useArray(dData,false,CPP_DEALLOC,14,1);
  const int expectedData[7]={1,2,3,4,5,7,19};
  //
  DataArrayInt *e=d->buildUnique();
 CPPUNIT_ASSERT_EQUAL(7,(int)e->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)e->getNumberOfComponents());
  for(int i=0;i<7;i++)
    CPPUNIT_ASSERT_EQUAL(expectedData[i],e->getIJ(i,0));
  //
  e->decrRef();
  d->decrRef();
}

void MEDCouplingBasicsTest5::testDAIPartitionByDifferentValues1()
{
  const int data[9]={1,0,1,2,0,2,2,-3,2};
  const int expected1[4]={-3,0,1,2};
  const int expected2_0[1]={7};
  const int expected2_1[2]={1,4};
  const int expected2_2[2]={0,2};
  const int expected2_3[4]={3,5,6,8};
  DataArrayInt *d=DataArrayInt::New();
  d->useArray(data,false,CPP_DEALLOC,9,1);
  std::vector<int> f;
  static const int nbOfOutputsExpected=4;
  std::vector<DataArrayInt *> e=d->partitionByDifferentValues(f);
  d->decrRef();
  CPPUNIT_ASSERT_EQUAL(nbOfOutputsExpected,(int)e.size());
  CPPUNIT_ASSERT_EQUAL(nbOfOutputsExpected,(int)f.size());
  for(int i=0;i<nbOfOutputsExpected;i++)
    {
      CPPUNIT_ASSERT_EQUAL(expected1[i],f[i]);
    }
  CPPUNIT_ASSERT_EQUAL((std::size_t)1,e[0]->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)2,e[1]->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)2,e[2]->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)4,e[3]->getNbOfElems());
 CPPUNIT_ASSERT_EQUAL(1,(int)e[0]->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(1,(int)e[1]->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(1,(int)e[2]->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(1,(int)e[3]->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected2_0,expected2_0+1,e[0]->begin()));
  CPPUNIT_ASSERT(std::equal(expected2_1,expected2_1+2,e[1]->begin()));
  CPPUNIT_ASSERT(std::equal(expected2_2,expected2_2+2,e[2]->begin()));
  CPPUNIT_ASSERT(std::equal(expected2_3,expected2_3+4,e[3]->begin()));
  e[0]->decrRef(); e[1]->decrRef(); e[2]->decrRef(); e[3]->decrRef();
}

void MEDCouplingBasicsTest5::testDAICheckMonotonic1()
{
  const int data1[6]={-1,0,2,2,4,5};
  const int data2[6]={6,2,0,-8,-9,-56};
  const int data3[6]={-1,0,3,2,4,6};
  const int data4[6]={7,5,2,3,0,-6};
  DataArrayInt *d=DataArrayInt::New();
  d->useArray(data1,false,CPP_DEALLOC,6,1);
  CPPUNIT_ASSERT(d->isMonotonic(true));
  CPPUNIT_ASSERT(!d->isMonotonic(false));
  d->checkMonotonic(true);
  CPPUNIT_ASSERT_THROW(d->checkMonotonic(false),INTERP_KERNEL::Exception);
  d->useArray(data2,false,CPP_DEALLOC,6,1);
  CPPUNIT_ASSERT(d->isMonotonic(false));
  CPPUNIT_ASSERT(!d->isMonotonic(true));
  d->checkMonotonic(false);
  CPPUNIT_ASSERT_THROW(d->checkMonotonic(true),INTERP_KERNEL::Exception);
  d->useArray(data3,false,CPP_DEALLOC,6,1);
  CPPUNIT_ASSERT(!d->isMonotonic(false));
  CPPUNIT_ASSERT(!d->isMonotonic(true));
  CPPUNIT_ASSERT_THROW(d->checkMonotonic(true),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(d->checkMonotonic(false),INTERP_KERNEL::Exception);
  d->useArray(data4,false,CPP_DEALLOC,6,1);
  CPPUNIT_ASSERT(!d->isMonotonic(false));
  CPPUNIT_ASSERT(!d->isMonotonic(true));
  CPPUNIT_ASSERT_THROW(d->checkMonotonic(true),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(d->checkMonotonic(false),INTERP_KERNEL::Exception);
  d->useArray(data4,false,CPP_DEALLOC,0,1);
  CPPUNIT_ASSERT(d->isMonotonic(true));
  CPPUNIT_ASSERT(d->isMonotonic(false));
  d->checkMonotonic(true);
  d->checkMonotonic(false);
  d->useArray(data4,false,CPP_DEALLOC,3,2);//throw because nbComp!=1
  CPPUNIT_ASSERT_THROW(d->isMonotonic(true),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(d->isMonotonic(false),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(d->checkMonotonic(true),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(d->checkMonotonic(false),INTERP_KERNEL::Exception);
  d->decrRef();
}

void MEDCouplingBasicsTest5::testIntersect2DMeshesTmp6()
{
  // coordinates
  DataArrayDouble *coords=DataArrayDouble::New();
  const double coordsData[16]={2.7554552980815448e-15,45,-45,5.5109105961630896e-15,-31.819805153394636,31.81980515339464,2.8779199779962799e-15,47,2.8166876380389124e-15,46,-47,5.7558399559925599e-15,-33.234018715767732,33.234018715767739,-46,5.6333752760778247e-15};
  coords->useArray(coordsData,false,CPP_DEALLOC,8,2);
  // connectivity
  DataArrayInt *conn=DataArrayInt::New();
  const int connData[9]={8,0,3,5,1,4,6,7,2};
  conn->useArray(connData,false,CPP_DEALLOC,9,1);
  DataArrayInt *connI=DataArrayInt::New();
  const int connIData[2]={0,9};
  connI->useArray(connIData,false,CPP_DEALLOC,2,1);
  MEDCouplingUMesh *m1=MEDCouplingUMesh::New("Fixe",2);
  m1->setCoords(coords);
  m1->setConnectivity(conn,connI,true);
  coords->decrRef(); conn->decrRef(); connI->decrRef();
  //
  coords=DataArrayDouble::New();
  const double coordsData2[26]={-7.3800475508445391,41.854329503018846,-3.7041190667754655,42.338274668899189,-3.7041190667754655,45.338274668899189,-7.3800475508445382,44.854329503018839,-5.5473631693521845,42.136406608386956,-3.7041190667754655,43.838274668899189,-5.5420833088100014,45.09630208595901,-7.3800475508445382,43.354329503018839,-3.7041190667754651,52.338274668899189,-7.3800475508445382,51.854329503018839,-3.7041190667754655,48.838274668899189,-5.5420833088100014,52.09630208595901,-7.3800475508445382,48.354329503018839};
  coords->useArray(coordsData2,false,CPP_DEALLOC,13,2);
  // connectivity
  conn=DataArrayInt::New();
  const int connData2[18]={8,0,1,2,3,4,5,6,7,8,3,2,8,9,6,10,11,12};
  conn->useArray(connData2,false,CPP_DEALLOC,18,1);
  connI=DataArrayInt::New();
  const int connIData2[3]={0,9,18};
  connI->useArray(connIData2,false,CPP_DEALLOC,3,1);
  //
  MEDCouplingUMesh *m2=MEDCouplingUMesh::New("Mobile",2);
  m2->setCoords(coords);
  m2->setConnectivity(conn,connI,true);
  coords->decrRef(); conn->decrRef(); connI->decrRef();
  //
  DataArrayInt *d1=0,*d2=0;
  MEDCouplingUMesh *m3=MEDCouplingUMesh::Intersect2DMeshes(m1,m2,1e-10,d1,d2);
  CPPUNIT_ASSERT_EQUAL(4,(int)m3->getNumberOfCells());
 CPPUNIT_ASSERT_EQUAL(4,(int)d1->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(4,(int)d2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(43,m3->getNumberOfNodes());
  bool areMerged=false;
  int newNbOfNodes=-1;
  m3->mergeNodes(1e-12,areMerged,newNbOfNodes)->decrRef();
  CPPUNIT_ASSERT_EQUAL(35,m3->getNumberOfNodes());
  m3->zipCoords();
  CPPUNIT_ASSERT_EQUAL(23,m3->getNumberOfNodes());
  //
  MEDCouplingFieldDouble *f=m3->getMeasureField(true);
  const double *vals=f->getArray()->getConstPointer();
  const double valuesExpected[4]={1.6603638692585716,5.747555728471923,129.68907101754394,7.4162714498559694};
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected[i],vals[i],1e-12);
  f->decrRef();
  //
  m1->decrRef();
  m2->decrRef();
  m3->decrRef();
  d1->decrRef();
  d2->decrRef();
}

void MEDCouplingBasicsTest5::testIntersect2DMeshesTmp7()
{
  double eps = 1.0e-8;
  // coordinates circle - SEE getCircle() on the Python side
  DataArrayDouble *coords1=DataArrayDouble::New();
  const double coordsData1[16]={0.5328427124746189, -0.08284271247461905, -0.03284271247461901, 0.4828427124746191, -0.03284271247461906, -0.082842712474619, 0.5328427124746191, 0.482842712474619};
  coords1->useArray(coordsData1,false,CPP_DEALLOC,8,2);
  // connectivity
  DataArrayInt *conn1=DataArrayInt::New();
  const int connData1[5]={INTERP_KERNEL::NORM_QPOLYG,0,1,2,3};
  conn1->useArray(connData1,false,CPP_DEALLOC,5,1);
  DataArrayInt *connI1=DataArrayInt::New();
  const int connIData1[2]={0,5};
  connI1->useArray(connIData1,false,CPP_DEALLOC,2,1);
  MEDCouplingUMesh *m1=MEDCouplingUMesh::New("circle",2);
  m1->setCoords(coords1);
  m1->setConnectivity(conn1,connI1,true);
  coords1->decrRef(); conn1->decrRef(); connI1->decrRef();

  // square
  DataArrayDouble *coords2=DataArrayDouble::New();
  const double coordsData2[8]={-0.5,-0.5,   -0.5, 0.5, 0.5, 0.5,    0.5,-0.5};
  coords2->useArray(coordsData2,false,CPP_DEALLOC,4,2);
  // connectivity
  DataArrayInt *conn2=DataArrayInt::New();
  const int connData2[5]={INTERP_KERNEL::NORM_POLYGON, 0,1,2,3};
  conn2->useArray(connData2,false,CPP_DEALLOC,5,1);
  DataArrayInt *connI2=DataArrayInt::New();
  const int connIData2[2]={0,5};
  connI2->useArray(connIData2,false,CPP_DEALLOC,2,1);
  MEDCouplingUMesh *m2=MEDCouplingUMesh::New("square",2);
  m2->setCoords(coords2);
  m2->setConnectivity(conn2,connI2,true);
  coords2->decrRef(); conn2->decrRef(); connI2->decrRef();

  DataArrayInt * resToM1 = 0, * resToM2 = 0;
  MEDCouplingUMesh *m_intersec=MEDCouplingUMesh::Intersect2DMeshes(m2, m1, eps, resToM1, resToM2);
  m_intersec->zipCoords();

  const double coo_tgt[34]={-0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, -0.5, -0.03284271247461901, 0.4828427124746191, \
    -0.014575131106459124, 0.5000000000000001, 0.5, -0.11224989991991996, 0.24271243444677046, 0.5, 0.5, 0.19387505004004, \
    -0.04799910280454185, -0.06682678787499614, -0.023843325638122054, 0.4915644577163915, 0.5, -0.30612494995996, 0.0, -0.5,\
    -0.5, 0.0, -0.25728756555322957, 0.5, -0.023843325638122026, 0.49156445771639157, -0.04799910280454181, -0.06682678787499613};
  const int conn_tgt[22]={32, 5, 2, 6, 4, 7, 8, 9, 10, 32, 6, 3, 0, 1, 5, 4, 11, 12, 13, 14, 15, 16};
  const int connI_tgt[3]={0, 9, 22};
  const int res1_tgt[2] = {0, 0};
  const int res2_tgt[2] = {0, -1};

  CPPUNIT_ASSERT(std::equal(conn_tgt,conn_tgt+22,m_intersec->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(connI_tgt,connI_tgt+3,m_intersec->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(res1_tgt,res1_tgt+2,resToM1->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(res2_tgt,res2_tgt+2,resToM2->getConstPointer()));
  for(int i=0;i<34;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(coo_tgt[i],m_intersec->getCoords()->getIJ(0,i),1e-12);
  m1->decrRef(); m2->decrRef(); m_intersec->decrRef();
  resToM1->decrRef(); resToM2->decrRef();
}

void MEDCouplingBasicsTest5::testDAIBuildSubstractionOptimized1()
{
  const int tab1[7]={1,3,5,6,7,9,13};
  const int tab2[3]={3,5,9};
  const int tab3[3]={1,3,5};
  DataArrayInt *da1=DataArrayInt::New(); da1->useArray(tab1,false,CPP_DEALLOC,7,1);
  DataArrayInt *da2=DataArrayInt::New(); da2->useArray(tab2,false,CPP_DEALLOC,3,1);
  DataArrayInt *da3=DataArrayInt::New(); da3->useArray(tab3,false,CPP_DEALLOC,3,1);
  DataArrayInt *da4=DataArrayInt::New(); da4->useArray(tab1,false,CPP_DEALLOC,7,1);
  //
  DataArrayInt *a=0;
  a=da1->buildSubstractionOptimized(da2);
 CPPUNIT_ASSERT_EQUAL(4,(int)a->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)a->getNumberOfComponents());
  const int expected1_0[4]={1,6,7,13};
  CPPUNIT_ASSERT(std::equal(expected1_0,expected1_0+4,a->begin()));
  a->decrRef();
  //
  a=da1->buildSubstractionOptimized(da3);
 CPPUNIT_ASSERT_EQUAL(4,(int)a->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)a->getNumberOfComponents());
  const int expected2_0[4]={6,7,9,13};
  CPPUNIT_ASSERT(std::equal(expected2_0,expected2_0+4,a->begin()));
  a->decrRef();
  //
  a=da1->buildSubstractionOptimized(da4);
 CPPUNIT_ASSERT_EQUAL(0,(int)a->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)a->getNumberOfComponents());
  a->decrRef();
  //
  da1->decrRef();
  da2->decrRef();
  da3->decrRef();
  da4->decrRef();
}

void MEDCouplingBasicsTest5::testDAIIsStrictlyMonotonic1()
{
  const int tab1[7]={1,3,5,6,7,9,13};
  DataArrayInt *da1=DataArrayInt::New(); da1->useArray(tab1,false,CPP_DEALLOC,7,1);
  CPPUNIT_ASSERT(da1->isStrictlyMonotonic(true));
  da1->checkStrictlyMonotonic(true);
  CPPUNIT_ASSERT(da1->isMonotonic(true));
  da1->checkMonotonic(true);
  CPPUNIT_ASSERT(!da1->isStrictlyMonotonic(false));
  CPPUNIT_ASSERT_THROW(da1->checkStrictlyMonotonic(false),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(!da1->isMonotonic(false));
  CPPUNIT_ASSERT_THROW(da1->checkMonotonic(false),INTERP_KERNEL::Exception);
  da1->decrRef();
  //
  int tab2[7]={1,3,5,6,6,9,13};
  da1=DataArrayInt::New(); da1->useArray(tab2,false,CPP_DEALLOC,7,1);
  CPPUNIT_ASSERT(!da1->isStrictlyMonotonic(true));
  CPPUNIT_ASSERT_THROW(da1->checkStrictlyMonotonic(true),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(da1->isMonotonic(true));
  da1->checkMonotonic(true);
  CPPUNIT_ASSERT(!da1->isStrictlyMonotonic(false));
  CPPUNIT_ASSERT_THROW(da1->checkStrictlyMonotonic(false),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(!da1->isMonotonic(false));
  CPPUNIT_ASSERT_THROW(da1->checkMonotonic(false),INTERP_KERNEL::Exception);
  da1->decrRef();
  //
  const int tab3[7]={1,3,5,6,5,9,13};
  da1=DataArrayInt::New(); da1->useArray(tab3,false,CPP_DEALLOC,7,1);
  CPPUNIT_ASSERT(!da1->isStrictlyMonotonic(true));
  CPPUNIT_ASSERT_THROW(da1->checkStrictlyMonotonic(true),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(!da1->isMonotonic(true));
  CPPUNIT_ASSERT_THROW(da1->checkMonotonic(true),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(!da1->isStrictlyMonotonic(false));
  CPPUNIT_ASSERT_THROW(da1->checkStrictlyMonotonic(false),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(!da1->isMonotonic(false));
  CPPUNIT_ASSERT_THROW(da1->checkMonotonic(false),INTERP_KERNEL::Exception);
  da1->decrRef();
  //
  const int tab4[7]={13,9,7,6,5,3,1};
  da1=DataArrayInt::New(); da1->useArray(tab4,false,CPP_DEALLOC,7,1);
  CPPUNIT_ASSERT(!da1->isStrictlyMonotonic(true));
  CPPUNIT_ASSERT_THROW(da1->checkStrictlyMonotonic(true),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(!da1->isMonotonic(true));
  CPPUNIT_ASSERT_THROW(da1->checkMonotonic(true),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(da1->isStrictlyMonotonic(false));
  da1->checkStrictlyMonotonic(false);
  CPPUNIT_ASSERT(da1->isMonotonic(false));
  da1->checkMonotonic(false);
  da1->decrRef();
  //
  const int tab5[7]={13,9,6,6,5,3,1};
  da1=DataArrayInt::New(); da1->useArray(tab5,false,CPP_DEALLOC,7,1);
  CPPUNIT_ASSERT(!da1->isStrictlyMonotonic(true));
  CPPUNIT_ASSERT_THROW(da1->checkStrictlyMonotonic(true),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(!da1->isMonotonic(true));
  CPPUNIT_ASSERT_THROW(da1->checkMonotonic(true),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(!da1->isStrictlyMonotonic(false));
  CPPUNIT_ASSERT_THROW(da1->checkStrictlyMonotonic(false),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(da1->isMonotonic(false));
  da1->checkMonotonic(false);
  da1->decrRef();
  //
  const int tab6[7]={13,9,5,6,5,3,1};
  da1=DataArrayInt::New(); da1->useArray(tab6,false,CPP_DEALLOC,7,1);
  CPPUNIT_ASSERT(!da1->isStrictlyMonotonic(true));
  CPPUNIT_ASSERT_THROW(da1->checkStrictlyMonotonic(true),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(!da1->isMonotonic(true));
  CPPUNIT_ASSERT_THROW(da1->checkMonotonic(true),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(!da1->isStrictlyMonotonic(false));
  CPPUNIT_ASSERT_THROW(da1->checkStrictlyMonotonic(false),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(!da1->isMonotonic(false));
  CPPUNIT_ASSERT_THROW(da1->checkMonotonic(false),INTERP_KERNEL::Exception);
  da1->decrRef();
  //
  da1=DataArrayInt::New(); da1->useArray(tab1,false,CPP_DEALLOC,0,1);
  CPPUNIT_ASSERT(da1->isStrictlyMonotonic(true));
  da1->checkStrictlyMonotonic(true);
  CPPUNIT_ASSERT(da1->isMonotonic(true));
  da1->checkMonotonic(true);
  CPPUNIT_ASSERT(da1->isStrictlyMonotonic(false));
  da1->checkStrictlyMonotonic(false);
  CPPUNIT_ASSERT(da1->isMonotonic(false));
  da1->checkMonotonic(false);
  da1->decrRef();
  //
  da1=DataArrayInt::New(); da1->useArray(tab1,false,CPP_DEALLOC,1,1);
  CPPUNIT_ASSERT(da1->isStrictlyMonotonic(true));
  da1->checkStrictlyMonotonic(true);
  CPPUNIT_ASSERT(da1->isMonotonic(true));
  da1->checkMonotonic(true);
  CPPUNIT_ASSERT(da1->isStrictlyMonotonic(false));
  da1->checkStrictlyMonotonic(false);
  CPPUNIT_ASSERT(da1->isMonotonic(false));
  da1->checkMonotonic(false);
  da1->decrRef();
}

void MEDCouplingBasicsTest5::testSimplexize3()
{
  const int conn[24]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  MEDCouplingUMesh *m=MEDCouplingUMesh::New("toto",3);
  m->allocateCells(0);
  m->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,conn+0);
  m->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+4);
  m->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+12);
  m->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,conn+20);
  const double coords[72]={0.,0.,0.,0.,1.,0.,1.,0.,0.,0.,0.,1.,2.,0.,0.,2.,1.,0.,3.,1.,0.,3.,0.,0.,2.,0.,1.,2.,1.,1.,3.,1.,1.,3.,0.,1.,4.,0.,0.,4.,1.,0.,5.,1.,0.,5.,0.,0.,4.,0.,1.,4.,1.,1.,5.,1.,1.,5.,0.,1.,6.,0.,0.,6.,1.,0.,7.,0.,0.,6.,0.,1.};
  DataArrayDouble *c=DataArrayDouble::New();
  c->useArray(coords,false,CPP_DEALLOC,24,3);
  m->setCoords(c);
  c->decrRef();
  m->checkConsistency();
  //
  MEDCouplingUMesh *m1=static_cast<MEDCouplingUMesh *>(m->deepCopy());
  DataArrayInt *d1=m1->simplexize(INTERP_KERNEL::PLANAR_FACE_5);
  m1->checkConsistency();
  MEDCouplingFieldDouble *f1=m1->getMeasureField(ON_CELLS);
  const double vol1Expected[12]={1./6, 1./6, 1./6,1./6, 1./6, 1./3,1./6, 1./6, 1./6, 1./6, 1./3, 1./6};
  CPPUNIT_ASSERT_EQUAL(1,(int)f1->getArray()->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(12,(int)f1->getArray()->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vol1Expected[i],f1->getIJ(i,0),1e-12);
  const int connExpected1[60]={14,0,1,2,3,14,4,9,5,6,14,4,8,9,11,14,4,7,11,6,14,9,11,10,6,14,4,9,6,11,14,12,17,13,14,14,12,16,17,19,14,12,15,19,14,14,17,19,18,14,14,12,17,14,19,14,20,21,22,23};
  const int connIExpected1[13]={0,5,10,15,20,25,30,35,40,45,50,55,60};
  const int n2o1[12]={0,1,1,1,1,1,2,2,2,2,2,3};
  CPPUNIT_ASSERT_EQUAL(1,(int)m1->getNodalConnectivity()->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(60,(int)m1->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)m1->getNodalConnectivityIndex()->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(13,(int)m1->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(connExpected1,connExpected1+60,m1->getNodalConnectivity()->begin()));
  CPPUNIT_ASSERT(std::equal(connIExpected1,connIExpected1+13,m1->getNodalConnectivityIndex()->begin()));
 CPPUNIT_ASSERT_EQUAL(1,(int)d1->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(12,(int)d1->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(n2o1,n2o1+12,d1->begin()));
  f1->decrRef();
  m1->decrRef();
  d1->decrRef();
  //
  MEDCouplingUMesh *m2=static_cast<MEDCouplingUMesh *>(m->deepCopy());
  DataArrayInt *d2=m2->simplexize(INTERP_KERNEL::PLANAR_FACE_6);
  m2->checkConsistency();
  MEDCouplingFieldDouble *f2=m2->getMeasureField(ON_CELLS);
  const double vol2Expected[14]={1./6, 1./6, 1./6,1./6, 1./6, 1./6,1./6,1./6, 1./6, 1./6, 1./6, 1./6,1./6,1./6};
  CPPUNIT_ASSERT_EQUAL(1,(int)f2->getArray()->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(14,(int)f2->getArray()->getNumberOfTuples());
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vol2Expected[i],f2->getIJ(i,0),1e-12);
  const int connExpected2[70]={14,0,1,2,3,14,4,9,5,10,14,4,5,6,10,14,4,8,9,10,14,4,11,8,10,14,4,6,7,10,14,4,7,11,10,14,12,17,13,18,14,12,13,14,18,14,12,16,17,18,14,12,19,16,18,14,12,14,15,18,14,12,15,19,18,14,20,21,22,23};
  const int connIExpected2[15]={0,5,10,15,20,25,30,35,40,45,50,55,60,65,70};
  const int n2o2[14]={0,1,1,1,1,1,1,2,2,2,2,2,2,3};
  CPPUNIT_ASSERT_EQUAL(1,(int)m2->getNodalConnectivity()->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(70,(int)m2->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)m2->getNodalConnectivityIndex()->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(15,(int)m2->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(connExpected2,connExpected2+70,m2->getNodalConnectivity()->begin()));
  CPPUNIT_ASSERT(std::equal(connIExpected2,connIExpected2+15,m2->getNodalConnectivityIndex()->begin()));
 CPPUNIT_ASSERT_EQUAL(1,(int)d2->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(14,(int)d2->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(n2o2,n2o2+14,d2->begin()));
  f2->decrRef();
  m2->decrRef();
  d2->decrRef();
  //
  m->decrRef();
}
