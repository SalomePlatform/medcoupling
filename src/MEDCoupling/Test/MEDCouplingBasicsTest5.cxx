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
