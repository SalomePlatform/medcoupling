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

#include "MEDCouplingBasicsTest5.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingGaussLocalization.hxx"
#include "MEDCouplingMultiFields.hxx"
#include "MEDCouplingFieldOverTime.hxx"

#include <cmath>
#include <functional>
#include <iterator>

using namespace ParaMEDMEM;

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
  MEDCouplingUMesh *m11=static_cast<MEDCouplingUMesh *>(m1->deepCpy());
  m11->tessellate2D(1.);
  CPPUNIT_ASSERT(m11->getCoords()->isEqual(*m11->getCoords(),1e-12));
  const int expected1[48]={5,0,3,11,1,5,3,4,12,2,1,11,5,5,15,3,0,5,6,16,4,3,15,5,5,5,0,7,19,5,6,5,19,7,8,20,5,0,1,23,7,5,1,2,24,8,7,23};
  const int expected2[9]={0,5,12,17,24,29,36,41,48};
  CPPUNIT_ASSERT_EQUAL(48,m11->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(9,m11->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+48,m11->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+9,m11->getNodalConnectivityIndex()->getConstPointer()));
  m11->decrRef();
  //
  MEDCouplingUMesh *m12=static_cast<MEDCouplingUMesh *>(m1->deepCpy());
  m12->tessellate2D(0.5);
  CPPUNIT_ASSERT_EQUAL(41,m12->getNumberOfNodes());
  const int expected3[60]={5,0,3,25,26,1,5,3,4,27,28,2,1,26,25,5,5,29,30,3,0,5,6,31,32,4,3,30,29,5,5,5,0,7,33,34,5,6,5,34,33,7,8,35,36,5,0,1,37,38,7,5,1,2,39,40,8,7,38,37};
  const int expected4[9]={0,6,15,21,30,36,45,51,60};
  const double expected5[82]={0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1.,0.,-1.5,0.5,0.,1.25,0.,0.7071067811865476,0.7071067811865476,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.7071067811865476,0.7071067811865476,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.7071067811865476,-0.7071067811865476,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.7071067811865476,-0.7071067811865476,1.0606601717798214,-1.0606601717798214,0.479425538604203,0.8775825618903728,0.8414709848078964,0.54030230586814,0.7191383079063044,1.3163738428355591,1.2622064772118446,0.8104534588022099,-0.877582561890373,0.4794255386042027,-0.5403023058681399,0.8414709848078964,-1.3163738428355596,0.7191383079063038,-0.8104534588022098,1.2622064772118446,-0.4794255386042031,-0.8775825618903728,-0.8414709848078965,-0.5403023058681399,-0.7191383079063045,-1.3163738428355591,-1.2622064772118449,-0.8104534588022098,0.8775825618903729,-0.47942553860420295,0.54030230586814,-0.8414709848078964,1.3163738428355594,-0.7191383079063043,0.8104534588022099,-1.2622064772118446};
  for(int i=0;i<82;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected5[i],m12->getCoords()->getIJ(0,i),1e-12);
  CPPUNIT_ASSERT_EQUAL(60,m12->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(9,m12->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+60,m12->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+9,m12->getNodalConnectivityIndex()->getConstPointer()));
  m12->decrRef();
  //
  m1->decrRef();
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
  const int expected1[12]={0,0,1,2,2,3,4,4,5,6,6,7};
  const int expected2[12]={0,1,1,2,3,3,4,5,5,6,7,7};
  CPPUNIT_ASSERT_EQUAL(12,d1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(12,d2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(12,m3->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(88,m3->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,m3->getSpaceDimension());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+12,d1->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+12,d2->getConstPointer()));
  const int expected3[100]={6,16,15,18,44,45,46,8,18,2,1,16,47,48,49,50,8,17,1,2,40,51,52,53,54,6,18,15,20,55,56,57,8,20,7,6,18,58,59,60,61,8,41,6,7,21,62,63,64,65,6,20,15,22,66,67,68,8,22,11,7,20,69,70,71,72,8,21,7,11,42,73,74,75,76,6,22,15,16,77,78,79,8,16,1,13,22,80,81,82,83,8,43,13,1,17,84,85,86,87};
  const int expected4[13]={0,7,16,25,32,41,50,57,66,75,82,91,100};
  const double expected5[176]={0.,0.,1.1,0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1.,-1.1,-1.,0.,-1.,1.1,-1.,1.7,-1.,0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1.,0.,-1.5,0.5,0.,1.25,0.,0.7071067811865476,0.7071067811865476,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.7071067811865476,0.7071067811865476,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.7071067811865476,-0.7071067811865476,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.7071067811865476,-0.7071067811865476,1.0606601717798214,-1.0606601717798214,1.1180339887498951,1.,-1.1180339887498951,1.,-1.1180339887498951,-1.,1.1180339887498951,-1.,0.5,0.,0.,0.5,0.7071067811865477,0.7071067811865476,0.55,1.,1.1,0.5,1.05,0.,0.7071067811865477,0.7071067811865477,1.3,0.,1.1,0.5,1.1090169943749475,1.,1.4012585384440737,0.535233134659635,0.,0.5,-0.5,0.,-0.7071067811865477,0.7071067811865476,-1.05,0.,-1.1,0.5,-0.55,1.,-0.7071067811865477,0.7071067811865477,-1.1090169943749475,1.,-1.1,0.5,-1.3,0.,-1.4012585384440737,0.5352331346596344,-0.5,0.,0.,-0.5,-0.7071067811865475,-0.7071067811865477,-0.55,-1.,-1.1,-0.5,-1.05,0.,-0.7071067811865479,-0.7071067811865476,-1.3,0.,-1.1,-0.5,-1.1090169943749475,-1.,-1.4012585384440734,-0.5352331346596354,0.,-0.5,0.5,0.,0.7071067811865475,-0.7071067811865477,1.05,0.,1.1,-0.5,0.55,-1.,0.7071067811865477,-0.7071067811865476,1.1090169943749475,-1.,1.1,-0.5,1.3,0.,1.4012585384440737,-0.535233134659635};
  CPPUNIT_ASSERT_EQUAL(100,m3->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(13,m3->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+100,m3->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+13,m3->getNodalConnectivityIndex()->getConstPointer()));
  for(int i=0;i<176;i++)
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
  CPPUNIT_ASSERT_EQUAL(9,ids1->getNumberOfTuples());
  const int expected1[9]={1,3,4,7,9,10,13,15,16};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+9,ids1->getConstPointer()));
  const double vec2[3]={0.,0.,1.};
  DataArrayInt *ids2=mesh3D->getCellIdsCrossingPlane(origin,vec2,1e-10);
  const int expected2[6]={6,7,8,9,10,11};
  CPPUNIT_ASSERT_EQUAL(6,ids2->getNumberOfTuples());
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
  CPPUNIT_ASSERT_EQUAL(9,slice1->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(9,ids->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(47,slice1->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(10,slice1->getNodalConnectivityIndex()->getNumberOfTuples());
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
  CPPUNIT_ASSERT_EQUAL(9,slice1->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(9,ids->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(49,slice1->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(10,slice1->getNodalConnectivityIndex()->getNumberOfTuples());
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
  CPPUNIT_ASSERT_EQUAL(12,slice1->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(12,ids->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(68,slice1->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(13,slice1->getNodalConnectivityIndex()->getNumberOfTuples());
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
  CPPUNIT_ASSERT_EQUAL(25,slice1->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(25,ids->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(75,slice1->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(26,slice1->getNodalConnectivityIndex()->getNumberOfTuples());
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
  CPPUNIT_ASSERT_EQUAL(68,slice1->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(68,ids->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(204,slice1->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(69,slice1->getNodalConnectivityIndex()->getNumberOfTuples());
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
  tmp=da->selectByTupleRanges(p);
  CPPUNIT_ASSERT(tmp->isEqual(*da,1e-14));
  tmp->decrRef();
  p[0].first=0; p[0].second=2; p[1].first=3; p[1].second=4; p[2].first=5; p[2].second=7;
  tmp=da->selectByTupleRanges(p);
  const double expected1[10]={1.,11.,2.,12.,4.,14.,6.,16.,7.,17.};
  CPPUNIT_ASSERT_EQUAL(5,tmp->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,tmp->getNumberOfComponents());
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],tmp->getIJ(0,i),1e-14);
  tmp->decrRef();
  p[0].first=0; p[0].second=2; p[1].first=0; p[1].second=2; p[2].first=5; p[2].second=6;
  tmp=da->selectByTupleRanges(p);
  const double expected2[10]={1.,11.,2.,12.,1.,11.,2.,12.,6.,16.};
  CPPUNIT_ASSERT_EQUAL(5,tmp->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,tmp->getNumberOfComponents());
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
  DataArrayDouble *dac=da->deepCpy();
  dac->setContigPartOfSelectedValues2(1,da2,2,4,1);
  const double expected3[14]={1.,11.,0.,30.,11.,41.,4.,14.,5.,15.,6.,16.,7.,17.};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],dac->getIJ(0,i),1e-14);
  dac->decrRef();
  //
  dac=da->deepCpy();
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues2(3,da2,0,5,1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues2(0,da2,4,6,1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues2(3,da2,5,0,1),INTERP_KERNEL::Exception);
  dac->setContigPartOfSelectedValues2(3,da2,1,5,1);
  const double expected4[14]={1.,11.,2.,12.,3.,13.,9.,39.,0.,30.,11.,41.,12.,42.};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected4[i],dac->getIJ(0,i),1e-14);
  dac->decrRef();
  //
  DataArrayInt *ids=DataArrayInt::New();
  ids->alloc(3,1);
  dac=da->deepCpy();
  ids->setIJ(0,0,2); ids->setIJ(1,0,0); ids->setIJ(2,0,4);
  dac->setContigPartOfSelectedValues(2,da2,ids);
  const double expected5[14]={1.,11.,2.,12.,0.,30.,8.,38.,12.,42.,6.,16.,7.,17.};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected5[i],dac->getIJ(0,i),1e-14);
  dac->decrRef();
  //
  dac=da->deepCpy();
  ids->setIJ(0,0,2); ids->setIJ(1,0,5); ids->setIJ(2,0,4);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues(1,da2,ids),INTERP_KERNEL::Exception);
  ids->setIJ(0,0,2); ids->setIJ(1,0,2); ids->setIJ(2,0,-1);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues(1,da2,ids),INTERP_KERNEL::Exception);
  ids->setIJ(0,0,2); ids->setIJ(1,0,2); ids->setIJ(2,0,1);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues(5,da2,ids),INTERP_KERNEL::Exception);
  dac->decrRef();
  //
  ids->setIJ(0,0,2); ids->setIJ(1,0,2); ids->setIJ(2,0,1);
  dac=da->deepCpy();
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
  tmp=da->selectByTupleRanges(p);
  CPPUNIT_ASSERT(tmp->isEqual(*da));
  tmp->decrRef();
  p[0].first=0; p[0].second=2; p[1].first=3; p[1].second=4; p[2].first=5; p[2].second=7;
  tmp=da->selectByTupleRanges(p);
  const int expected1[10]={1,11,2,12,4,14,6,16,7,17};
  CPPUNIT_ASSERT_EQUAL(5,tmp->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,tmp->getNumberOfComponents());
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],tmp->getIJ(0,i));
  tmp->decrRef();
  p[0].first=0; p[0].second=2; p[1].first=0; p[1].second=2; p[2].first=5; p[2].second=6;
  tmp=da->selectByTupleRanges(p);
  const int expected2[10]={1,11,2,12,1,11,2,12,6,16};
  CPPUNIT_ASSERT_EQUAL(5,tmp->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,tmp->getNumberOfComponents());
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
  DataArrayInt *dac=da->deepCpy();
  dac->setContigPartOfSelectedValues2(1,da2,2,4,1);
  const int expected3[14]={1,11,0,30,11,41,4,14,5,15,6,16,7,17};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_EQUAL(expected3[i],dac->getIJ(0,i));
  dac->decrRef();
  //
  dac=da->deepCpy();
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues2(3,da2,0,5,1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues2(0,da2,4,6,1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues2(3,da2,5,0,1),INTERP_KERNEL::Exception);
  dac->setContigPartOfSelectedValues2(3,da2,1,5,1);
  const int expected4[14]={1,11,2,12,3,13,9,39,0,30,11,41,12,42};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_EQUAL(expected4[i],dac->getIJ(0,i));
  dac->decrRef();
  //
  DataArrayInt *ids=DataArrayInt::New();
  ids->alloc(3,1);
  dac=da->deepCpy();
  ids->setIJ(0,0,2); ids->setIJ(1,0,0); ids->setIJ(2,0,4);
  dac->setContigPartOfSelectedValues(2,da2,ids);
  const int expected5[14]={1,11,2,12,0,30,8,38,12,42,6,16,7,17};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_EQUAL(expected5[i],dac->getIJ(0,i));
  dac->decrRef();
  //
  dac=da->deepCpy();
  ids->setIJ(0,0,2); ids->setIJ(1,0,5); ids->setIJ(2,0,4);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues(1,da2,ids),INTERP_KERNEL::Exception);
  ids->setIJ(0,0,2); ids->setIJ(1,0,2); ids->setIJ(2,0,-1);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues(1,da2,ids),INTERP_KERNEL::Exception);
  ids->setIJ(0,0,2); ids->setIJ(1,0,2); ids->setIJ(2,0,1);
  CPPUNIT_ASSERT_THROW(dac->setContigPartOfSelectedValues(5,da2,ids),INTERP_KERNEL::Exception);
  dac->decrRef();
  //
  ids->setIJ(0,0,2); ids->setIJ(1,0,2); ids->setIJ(2,0,1);
  dac=da->deepCpy();
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
  mesh2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(2,mesh2->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(30,mesh2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(31,revDescIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(31,revDescIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(13,descIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(13,descIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(48,desc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(48,desc->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(48,revDesc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(48,revDesc->getNumberOfTuples());
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
      MEDCouplingUMesh *m1Cpy=static_cast<MEDCouplingUMesh *>(m1->deepCpy());
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
  DataArrayInt *d2=d1->convertToIntArr();
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
  d2->decrRef();
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
  CPPUNIT_ASSERT_EQUAL(5,arr1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,arr1->getNumberOfComponents());
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
  CPPUNIT_ASSERT_EQUAL(4,m1c->getNumberOfCells());
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
  CPPUNIT_ASSERT_EQUAL(12,m1c->getNumberOfCells());
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
  CPPUNIT_ASSERT_EQUAL(24,m1c->getNumberOfCells());
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
  mesh2D->checkCoherency();
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
  mesh3D->checkCoherency();
  //
  MEDCouplingUMesh *mesh3D_2=dynamic_cast<MEDCouplingUMesh *>(mesh3D->deepCpy());
  MEDCouplingUMesh *mesh2D_2=dynamic_cast<MEDCouplingUMesh *>(mesh2D->deepCpy());
  MEDCouplingUMesh *mesh3D_4=dynamic_cast<MEDCouplingUMesh *>(mesh3D->deepCpy());
  MEDCouplingUMesh *mesh2D_4=dynamic_cast<MEDCouplingUMesh *>(mesh2D->deepCpy());
  DataArrayInt *renumNodes=DataArrayInt::New();
  int oldNbOf3DNodes=mesh3D->getNumberOfNodes();
  renumNodes->alloc(mesh2D->getNumberOfNodes(),1);
  renumNodes->iota(oldNbOf3DNodes);
  DataArrayDouble *coo=DataArrayDouble::Aggregate(mesh3D->getCoords(),mesh2D->getCoords());
  mesh3D->setCoords(coo);
  mesh2D->setCoords(coo);
  coo->decrRef();
  MEDCouplingUMesh *mesh2D_3=dynamic_cast<MEDCouplingUMesh *>(mesh2D->deepCpy());
  mesh2D_3->shiftNodeNumbersInConn(oldNbOf3DNodes);
  mesh2D->renumberNodesInConn(renumNodes->getConstPointer());
  renumNodes->decrRef();
  CPPUNIT_ASSERT(mesh2D_3->isEqual(mesh2D,1e-12));
  mesh2D_3->decrRef();
  //
  DataArrayInt *da1,*da2;
  mesh3D->checkGeoEquivalWith(mesh3D_2,10,1e-12,da1,da2);
  CPPUNIT_ASSERT(da1==0);
  CPPUNIT_ASSERT_EQUAL(8,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da2->getNumberOfComponents());
  const int expected1[8]={8,11,12,9,4,5,6,7};
  for(int i=0;i<8;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],da2->getIJ(i,0));
  da2->decrRef();
  //
  mesh2D->checkGeoEquivalWith(mesh2D_2,10,1e-12,da1,da2);
  CPPUNIT_ASSERT(da1==0);
  CPPUNIT_ASSERT_EQUAL(9,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da2->getNumberOfComponents());
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_EQUAL(8+i,da2->getIJ(i,0));
  da2->decrRef();
  //
  const double vect[3]={1.,0.,0.};
  MEDCouplingUMesh *mesh2D_5=dynamic_cast<MEDCouplingUMesh *>(mesh2D_4->deepCpy());
  mesh2D_5->translate(vect);
  std::vector<MEDCouplingUMesh *> meshes(3);
  meshes[0]=mesh3D_4; meshes[1]=mesh2D_4; meshes[2]=mesh2D_5;
  MEDCouplingUMesh::PutUMeshesOnSameAggregatedCoords(meshes);
  CPPUNIT_ASSERT(mesh3D_4->getCoords()==mesh2D_4->getCoords());
  CPPUNIT_ASSERT(mesh2D_4->getCoords()==mesh2D_5->getCoords());
  mesh3D_4->checkCoherency(); mesh2D_4->checkCoherency(); mesh2D_5->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(26,mesh3D_4->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,mesh3D_4->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(9,mesh3D_4->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(23,mesh2D_4->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(23,mesh2D_5->getNodalConnectivity()->getNumberOfTuples());
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
  mesh3D_4->checkCoherency(); mesh2D_4->checkCoherency(); mesh2D_5->checkCoherency();
  CPPUNIT_ASSERT(mesh3D_4->getCoords()==mesh2D_4->getCoords());
  CPPUNIT_ASSERT(mesh2D_4->getCoords()==mesh2D_5->getCoords());
  CPPUNIT_ASSERT_EQUAL(19,mesh3D_4->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,mesh3D_4->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(9,mesh3D_4->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(23,mesh2D_4->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(23,mesh2D_5->getNodalConnectivity()->getNumberOfTuples());
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
  CPPUNIT_ASSERT_EQUAL(6,d2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(10,d1->getNumberOfTuples());
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
  mesh2D->checkCoherency();
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
  CPPUNIT_ASSERT_EQUAL(3,d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,d->getConstPointer()));
  d->decrRef();
  //
  d=DataArrayInt::Range(2,23,7);
  CPPUNIT_ASSERT_EQUAL(3,d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,d->getConstPointer()));
  d->decrRef();
  //
  d=DataArrayInt::Range(2,24,7);
  const int expected2[4]={2,9,16,23};
  CPPUNIT_ASSERT_EQUAL(4,d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected2,expected2+4,d->getConstPointer()));
  d->decrRef();
  //
  d=DataArrayInt::Range(24,2,-7);
  const int expected3[4]={24,17,10,3};
  CPPUNIT_ASSERT_EQUAL(4,d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+4,d->getConstPointer()));
  d->decrRef();
  //
  d=DataArrayInt::Range(23,2,-7);
  const int expected4[3]={23,16,9};
  CPPUNIT_ASSERT_EQUAL(3,d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected4,expected4+3,d->getConstPointer()));
  d->decrRef();
  //
  d=DataArrayInt::Range(23,22,-7);
  CPPUNIT_ASSERT_EQUAL(1,d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,d->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(23,d->getIJ(0,0));
  d->decrRef();
  //
  d=DataArrayInt::Range(22,23,7);
  CPPUNIT_ASSERT_EQUAL(1,d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,d->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(22,d->getIJ(0,0));
  d->decrRef();
  //
  d=DataArrayInt::Range(22,22,7);
  CPPUNIT_ASSERT_EQUAL(0,d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,d->getNumberOfComponents());
  d->decrRef();
  //
  d=DataArrayInt::Range(22,22,-7);
  CPPUNIT_ASSERT_EQUAL(0,d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,d->getNumberOfComponents());
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
  CPPUNIT_ASSERT_EQUAL(3,arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells1,cells1+3,arr->getConstPointer()));
  arr->decrRef();
  // no permutation policy 1
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,1,arr));
  CPPUNIT_ASSERT_EQUAL(3,arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells1,cells1+3,arr->getConstPointer()));
  arr->decrRef();
  // no permutation policy 2
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,2,arr));
  CPPUNIT_ASSERT_EQUAL(3,arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells1,cells1+3,arr->getConstPointer()));
  arr->decrRef();
  // some modification into m2
  const int modif1[3]={2,4,5};
  std::copy(modif1,modif1+3,m2->getNodalConnectivity()->getPointer()+1);
  //policy 0 fails because cell0 in m2 has same orientation be not same connectivity
  const int expected1[3]={5,3,4};
  CPPUNIT_ASSERT(!m1->areCellsIncludedIn(m2,0,arr));
  CPPUNIT_ASSERT_EQUAL(3,arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,arr->getConstPointer()));
  arr->decrRef();
  //policy 1 succeeds because cell0 in m2 has not exactly the same conn
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,1,arr));
  CPPUNIT_ASSERT_EQUAL(3,arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells1,cells1+3,arr->getConstPointer()));
  arr->decrRef();
  //policy 2 succeeds because cell0 in m2 has same nodes in connectivity
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,2,arr));
  CPPUNIT_ASSERT_EQUAL(3,arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells1,cells1+3,arr->getConstPointer()));
  arr->decrRef();
  //some new modification into m2
  const int modif2[3]={2,5,4};
  std::copy(modif2,modif2+3,m2->getNodalConnectivity()->getPointer()+1);
  //policy 0 fails because cell0 in m2 has not exactly the same conn
  CPPUNIT_ASSERT(!m1->areCellsIncludedIn(m2,0,arr));
  CPPUNIT_ASSERT_EQUAL(3,arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,arr->getConstPointer()));
  arr->decrRef();
  //policy 1 fails too because cell0 in m2 has not same orientation
  CPPUNIT_ASSERT(!m1->areCellsIncludedIn(m2,1,arr));
  CPPUNIT_ASSERT_EQUAL(3,arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,arr->getConstPointer()));
  arr->decrRef();
  //policy 2 succeeds because cell0 in m2 has same nodes in connectivity
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,2,arr));
  CPPUNIT_ASSERT_EQUAL(3,arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr->getNumberOfComponents());
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
  CPPUNIT_ASSERT_EQUAL(2,arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells2,cells2+2,arr->getConstPointer()));
  arr->decrRef();
  // no permutation policy 1
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,1,arr));
  CPPUNIT_ASSERT_EQUAL(2,arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells2,cells2+2,arr->getConstPointer()));
  arr->decrRef();
  // no permutation policy 2
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,2,arr));
  CPPUNIT_ASSERT_EQUAL(2,arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(cells2,cells2+2,arr->getConstPointer()));
  arr->decrRef();
  // some modification into m2
  const int modif3[2]={4,3};
  std::copy(modif3,modif3+2,m2->getNodalConnectivity()->getPointer()+1);
  //policy 0 fails because cell0 in m2 has not exactly the same conn
  const int expected2[2]={4,2};
  CPPUNIT_ASSERT(!m1->areCellsIncludedIn(m2,0,arr));
  CPPUNIT_ASSERT_EQUAL(2,arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected2,expected2+2,arr->getConstPointer()));
  arr->decrRef();
  //policy 1 fails too because cell0 in m2 has not same orientation
  CPPUNIT_ASSERT(!m1->areCellsIncludedIn(m2,1,arr));
  CPPUNIT_ASSERT_EQUAL(2,arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected2,expected2+2,arr->getConstPointer()));
  arr->decrRef();
  //policy 2 succeeds because cell0 in m2 has same nodes in connectivity
  CPPUNIT_ASSERT(m1->areCellsIncludedIn(m2,2,arr));
  CPPUNIT_ASSERT_EQUAL(2,arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr->getNumberOfComponents());
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
  m->checkCoherency();
  //
  DataArrayInt *da=m->convexEnvelop2D();
  m->checkCoherency();
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
  DataArrayInt *theTest=vals->getIdsInRange(-1.,1e-7);
  CPPUNIT_ASSERT(theTest->isIdentity());
  CPPUNIT_ASSERT_EQUAL(331,theTest->getNumberOfTuples());
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
  DataArrayInt *arr1=arr->deepCpy();
  DataArrayInt *arr2=arr->deepCpy();
  arr1->sort(true);
  const int expected1[6]={1,2,4,5,6,7};
  CPPUNIT_ASSERT_EQUAL(6,arr1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr1->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+6,arr1->begin()));
  arr2->sort(false);
  const int expected2[6]={7,6,5,4,2,1};
  CPPUNIT_ASSERT_EQUAL(6,arr2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,arr2->getNumberOfComponents());
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
  DataArrayDouble *ard1=ard->deepCpy();
  DataArrayDouble *ard2=ard->deepCpy();
  ard1->sort(true);
  const double expected3[6]={1.,2.,4.,5.,6.,7.};
  CPPUNIT_ASSERT_EQUAL(6,ard1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,ard1->getNumberOfComponents());
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],ard1->getIJ(i,0),1e-12);
  ard2->sort(false);
  const double expected4[6]={7.,6.,5.,4.,2.,1.};
  CPPUNIT_ASSERT_EQUAL(6,ard2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,ard2->getNumberOfComponents());
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
  CPPUNIT_ASSERT_EQUAL(3,v2[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,v2[0]->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,v2[1]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,v2[1]->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(2,v2[2]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,v2[2]->getNumberOfComponents());
  //
  CPPUNIT_ASSERT(std::equal(expected0,expected0+3,v2[0]->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected1,expected1+5,v2[1]->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+2,v2[2]->getConstPointer()));
  v2[0]->decrRef();
  v2[1]->decrRef();
  v2[2]->decrRef();
  //
  MEDCouplingUMesh *m5=m4->buildSpreadZonesWithPoly();
  CPPUNIT_ASSERT_EQUAL(3,m5->getNumberOfCells());
  CPPUNIT_ASSERT(m5->getCoords()==m4->getCoords());
  const int expected3[23]={5,15,16,17,14,11,13,12,5,2,1,0,3,6,7,8,5,5,18,21,22,20,19};
  const int expected4[4]={0,8,17,23};
  CPPUNIT_ASSERT_EQUAL(23,m5->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+23,m5->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(4,m5->getNodalConnectivityIndex()->getNumberOfTuples());
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
  CPPUNIT_ASSERT_EQUAL(2,da->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected0,expected0+2,da->getConstPointer()));
  da->decrRef();
  //
  da=m->giveCellsWithType(INTERP_KERNEL::NORM_QUAD4);
  CPPUNIT_ASSERT_EQUAL(3,da->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,da->getConstPointer()));
  da->decrRef();
  //
  da=m->giveCellsWithType(INTERP_KERNEL::NORM_TRI6);
  CPPUNIT_ASSERT_EQUAL(0,da->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da->getNumberOfComponents());
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
  f->checkCoherency();
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
  f2->checkCoherency();
  DataArrayInt *cI=m->computeNbOfNodesPerCell();
  cI->computeOffsets2();
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
  CPPUNIT_ASSERT_EQUAL(18,skin->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(1,skin->getMeshDimension());
  CPPUNIT_ASSERT(skin->getCoords()==umesh->getCoords());
  const int expected1[19]={0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54};
  const int expected2[54]={1,1,0,1,0,5,1,2,1,1,3,2,1,4,3,1,9,4,1,5,10,1,14,9,1,10,15,1,19,14,1,15,20,1,24,19,1,20,25,1,25,26,1,26,27,1,27,28,1,28,29,1,29,24};
  CPPUNIT_ASSERT_EQUAL(19,skin->getNodalConnectivityIndex()->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+19,skin->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(54,skin->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(expected2,expected2+54,skin->getNodalConnectivity()->getConstPointer()));
  DataArrayInt *ids=skin->computeFetchedNodeIds();
  const int expected3[18]={0,1,2,3,4,5,9,10,14,15,19,20,24,25,26,27,28,29};
  CPPUNIT_ASSERT_EQUAL(18,ids->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+18,ids->getConstPointer()));
  MEDCouplingUMesh *part=dynamic_cast<MEDCouplingUMesh *>(umesh->buildFacePartOfMySelfNode(ids->begin(),ids->end(),true));
  part->setName(skin->getName());
  CPPUNIT_ASSERT(part->isEqual(skin,1e-12));
  MEDCouplingUMesh *part2=dynamic_cast<MEDCouplingUMesh *>(part->buildPartOfMySelf2(1,18,2,true));
  DataArrayInt *ids2=DataArrayInt::Range(0,18,2);
  part->setPartOfMySelf(ids2->begin(),ids2->end(),*part2);
  ids2->decrRef();
  CPPUNIT_ASSERT(!part->isEqual(skin,1e-12));
  DataArrayInt *trad=part->zipConnectivityTraducer(0);
  CPPUNIT_ASSERT_EQUAL(9,part->getNumberOfCells());
  const int expected4[18]={0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8};
  CPPUNIT_ASSERT(std::equal(expected4,expected4+18,trad->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(18,trad->getNbOfElems());
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
  CPPUNIT_ASSERT(s==m->getAllTypes());
  const int ids1[3]={0,3,4};
  MEDCouplingUMesh *part=static_cast<MEDCouplingUMesh *>(m->buildPartOfMySelf(ids1,ids1+3,true));
  part->simplexize(0)->decrRef();
  const int ids2[3]={1,2,5};
  MEDCouplingUMesh *part2=static_cast<MEDCouplingUMesh *>(part->buildPartOfMySelf(ids2,ids2+3,true));
  m->setPartOfMySelf(ids1,ids1+3,*part2);
  const int expected1[20]={3,0,4,1,3,1,4,2,3,4,5,2,3,6,7,4,3,7,5,4};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+20,m->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(20,m->getNodalConnectivity()->getNbOfElems());
  const int expected2[6]={0,4,8,12,16,20};
  CPPUNIT_ASSERT(std::equal(expected2,expected2+6,m->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(6,m->getNodalConnectivityIndex()->getNbOfElems());
  s.clear(); s.insert(INTERP_KERNEL::NORM_TRI3);
  CPPUNIT_ASSERT(s==m->getAllTypes());
  m->decrRef(); part->decrRef(); part2->decrRef();
  // no resize with explicit ids list
  m=build2DTargetMesh_1();
  part=static_cast<MEDCouplingUMesh *>(m->buildPartOfMySelf(ids1,ids1+2,true));
  part->convertAllToPoly();
  m->setPartOfMySelf(ids1+1,ids1+3,*part);
  const int expected3[23]={4,0,3,4,1,3,1,4,2,3,4,5,2,5,0,3,4,1,5,6,7,4,3};
  CPPUNIT_ASSERT(std::equal(expected3,expected3+23,m->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(23,m->getNodalConnectivity()->getNbOfElems());
  const int expected4[6]={0,5,9,13,18,23};
  CPPUNIT_ASSERT(std::equal(expected4,expected4+6,m->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(6,m->getNodalConnectivityIndex()->getNbOfElems());
  s.clear(); s.insert(INTERP_KERNEL::NORM_TRI3); s.insert(INTERP_KERNEL::NORM_QUAD4); s.insert(INTERP_KERNEL::NORM_POLYGON);
  CPPUNIT_ASSERT(s==m->getAllTypes());
  m->decrRef(); part->decrRef();
  // resize with range ids
  m=build2DTargetMesh_1();
  part=static_cast<MEDCouplingUMesh *>(m->buildPartOfMySelf2(3,5,1,true));
  m->setPartOfMySelf2(1,3,1,*part);
  const int expected5[25]={4,0,3,4,1,4,6,7,4,3,4,7,8,5,4,4,6,7,4,3,4,7,8,5,4};
  CPPUNIT_ASSERT(std::equal(expected5,expected5+25,m->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(25,m->getNodalConnectivity()->getNbOfElems());
  const int expected6[6]={0,5,10,15,20,25};
  CPPUNIT_ASSERT(std::equal(expected6,expected6+6,m->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(6,m->getNodalConnectivityIndex()->getNbOfElems());
  s.clear(); s.insert(INTERP_KERNEL::NORM_QUAD4);
  CPPUNIT_ASSERT(s==m->getAllTypes());
  m->decrRef(); part->decrRef();
  // no resize with range ids
  m=build2DTargetMesh_1();
  part=static_cast<MEDCouplingUMesh *>(m->buildPartOfMySelf2(0,5,3,true));
  part->convertAllToPoly();
  m->setPartOfMySelf2(3,5,1,*part);
  const int expected7[23]={4,0,3,4,1,3,1,4,2,3,4,5,2,5,0,3,4,1,5,6,7,4,3};
  CPPUNIT_ASSERT(std::equal(expected7,expected7+23,m->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(23,m->getNodalConnectivity()->getNbOfElems());
  const int expected8[6]={0,5,9,13,18,23};
  CPPUNIT_ASSERT(std::equal(expected8,expected8+6,m->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(6,m->getNodalConnectivityIndex()->getNbOfElems());
s.clear(); s.insert(INTERP_KERNEL::NORM_TRI3); s.insert(INTERP_KERNEL::NORM_QUAD4); s.insert(INTERP_KERNEL::NORM_POLYGON);
  CPPUNIT_ASSERT(s==m->getAllTypes());
  m->decrRef(); part->decrRef();
  // no resize with range ids negative direction
  m=build2DTargetMesh_1();
  part=static_cast<MEDCouplingUMesh *>(m->buildPartOfMySelf2(3,-1,-3,true));
  part->convertAllToPoly();
  m->setPartOfMySelf2(4,2,-1,*part);
  const int expected9[23]={4,0,3,4,1,3,1,4,2,3,4,5,2,5,0,3,4,1,5,6,7,4,3};
  CPPUNIT_ASSERT(std::equal(expected9,expected9+23,m->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(23,m->getNodalConnectivity()->getNbOfElems());
  const int expected10[6]={0,5,9,13,18,23};
  CPPUNIT_ASSERT(std::equal(expected10,expected10+6,m->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(6,m->getNodalConnectivityIndex()->getNbOfElems());
  s.clear(); s.insert(INTERP_KERNEL::NORM_TRI3); s.insert(INTERP_KERNEL::NORM_QUAD4); s.insert(INTERP_KERNEL::NORM_POLYGON);
  CPPUNIT_ASSERT(s==m->getAllTypes());
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
  m->checkCoherency();
  //
  MEDCouplingFieldDouble *vol=m->getMeasureField(ON_CELLS);
  CPPUNIT_ASSERT_EQUAL(1,vol->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,vol->getArray()->getIJ(0,0),1e-12);
  vol->decrRef();
  //
  m->unPolyze();
  CPPUNIT_ASSERT_EQUAL(1,m->getNumberOfCells());
  std::set<INTERP_KERNEL::NormalizedCellType> s; s.insert(INTERP_KERNEL::NORM_PENTA6);
  CPPUNIT_ASSERT(s==m->getAllTypes());
  //
  const int expected1[2]={0,7};
  const int expected2[7]={16,0,2,1,3,5,4};
  CPPUNIT_ASSERT_EQUAL(2,m->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+2,m->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(7,m->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected2,expected2+7,m->getNodalConnectivity()->getConstPointer()));
  //
  vol=m->getMeasureField(ON_CELLS);
  CPPUNIT_ASSERT_EQUAL(1,vol->getArray()->getNumberOfTuples());
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
  f->checkCoherency();
  //
  double *res0=new double[1];
  f->getValueOn(targetPointCoordsX,res0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(targetFieldValsExpected[0],res0[0],1e-10);
  delete [] res0;
  //
  cmesh->decrRef();
  umesh->decrRef();
  srcArrX->decrRef();
  srcVals->decrRef();
  f->decrRef();
}
