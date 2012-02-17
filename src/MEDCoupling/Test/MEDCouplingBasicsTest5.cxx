// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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
