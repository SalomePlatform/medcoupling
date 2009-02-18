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
#include "MEDCouplingFieldDouble.hxx"
#include "MemArray.hxx"
#include "Interpolation2D.txx"
#include "Interpolation3DSurf.txx"

#include "MEDCouplingNormalizedUnstructuredMesh.txx"

#include <cmath>

using namespace std;
using namespace ParaMEDMEM;

void MEDCouplingBasicsTest::testMesh()
{
  const int nbOfCells=6;
  const int nbOfNodes=12;
  
  double coords[3*nbOfNodes]={ 
    0.024155, 0.04183768725682622, -0.305, 0.04831000000000001, -1.015761910347357e-17, -0.305, 0.09662000000000001, -1.832979297858306e-18, 
    -0.305, 0.120775, 0.04183768725682623, -0.305, 0.09662000000000001, 0.08367537451365245, -0.305, 0.04831000000000001, 
    0.08367537451365246, -0.305, 0.024155, 0.04183768725682622, -0.2863, 0.04831000000000001, -1.015761910347357e-17, -0.2863, 
    0.09662000000000001, -1.832979297858306e-18, -0.2863, 0.120775, 0.04183768725682623, -0.2863, 0.09662000000000001, 0.08367537451365245, 
    -0.2863, 0.04831000000000001, 0.08367537451365246, -0.2863, };
  
  int tab4[4*nbOfCells]={ 
    1, 2, 8, 7, 2, 3, 9, 8, 3, 4, 10, 9, 4, 5, 11, 10, 5, 0, 6, 11, 
    0, 1, 7, 6, };
  
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(8);
  const int *curConn=tab4;
  for(int i=0;i<nbOfCells;i++,curConn+=4)
    mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,curConn);
  mesh->finishInsertingCells();
  CPPUNIT_ASSERT_EQUAL(30,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(nbOfCells,mesh->getNumberOfCells());
  //test 0 - no copy no ownership
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->useArray(coords,false,CPP_DEALLOC,nbOfNodes,3);
  mesh->setCoords(myCoords);
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  CPPUNIT_ASSERT_EQUAL(nbOfCells,mesh->getNumberOfCells());
  mesh->checkCoherency();
  //test 1 - no copy ownership C++
  myCoords=DataArrayDouble::New();
  double *tmp=new double[3*nbOfNodes];
  copy(coords,coords+3*nbOfNodes,tmp);
  myCoords->useArray(tmp,true,CPP_DEALLOC,nbOfNodes,3);
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  CPPUNIT_ASSERT_EQUAL(nbOfCells,mesh->getNumberOfCells());
  mesh->checkCoherency();
  //test 2 - no copy ownership C
  myCoords=DataArrayDouble::New();
  tmp=(double *)malloc(3*nbOfNodes*sizeof(double));
  copy(coords,coords+3*nbOfNodes,tmp);
  myCoords->useArray(tmp,true,C_DEALLOC,nbOfNodes,3);
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  CPPUNIT_ASSERT_EQUAL(nbOfNodes,mesh->getNumberOfNodes());
  mesh->checkCoherency();
  //test 3 - copy.
  myCoords=DataArrayDouble::New();
  myCoords->alloc(nbOfNodes,3);
  tmp=myCoords->getPointer();
  copy(coords,coords+3*nbOfNodes,tmp);
  // test 3 bis deepcopy
  DataArrayDouble *myCoords2=DataArrayDouble::New();
  *myCoords2=*myCoords;
  myCoords2->decrRef();
  //
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  CPPUNIT_ASSERT_EQUAL(nbOfNodes,mesh->getNumberOfNodes());
  mesh->checkCoherency();
  //test 4 - Field on cells
  MEDCouplingFieldDouble *fieldOnCells=MEDCouplingFieldDouble::New(ON_CELLS);
  fieldOnCells->setMesh(mesh);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(nbOfCells,9);
  fieldOnCells->setArray(array);
  tmp=array->getPointer();
  array->decrRef();
  fill(tmp,tmp+9*nbOfCells,7.);
  fieldOnCells->declareAsNew();
  fieldOnCells->checkCoherency();
  fieldOnCells->decrRef();
  //clean-up
  mesh->decrRef();
}

void MEDCouplingBasicsTest::test2DInterpP0P0_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  vector<map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[3]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Convex, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<3;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
      CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[0][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[2][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[4][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[4][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test2DInterpP0P1_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  vector<map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[3]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Convex, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<3;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");
      CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[0][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333329,res[2][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[3][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[4][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[4][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[5][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333329,res[6][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[7][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[8][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[8][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.25,sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test2DInterpP1P0_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  vector<map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[2]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<2;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
      CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[3][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333333,res[1][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333333,res[2][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666667,res[3][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[2][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[3][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[4][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DSurfInterpP0P0_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  vector<map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[3]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Convex, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<3;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
      CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[2][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[3][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[4][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[4][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.*sqrt(2.),sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DSurfInterpP0P1_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  vector<map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[2]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<2;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");
      CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[0][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333329*sqrt(2.),res[2][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666*sqrt(2.),res[3][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666*sqrt(2.),res[4][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666*sqrt(2.),res[4][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[5][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333329*sqrt(2.),res[6][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666*sqrt(2.),res[7][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[8][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[8][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.25*sqrt(2.),sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DSurfInterpP1P0_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  vector<map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[2]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<2;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
      CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[3][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333333*sqrt(2.),res[1][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.083333333333333333*sqrt(2.),res[2][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666667*sqrt(2.),res[3][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[2][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[3][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[4][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.*sqrt(2.),sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DInterpP0P0_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  //clean up
  sourceMesh->decrRef();
}

MEDCouplingUMesh *MEDCouplingBasicsTest::build2DSourceMesh_1()
{
  double sourceCoords[8]={-0.3,-0.3, 0.7,-0.3, -0.3,0.7, 0.7,0.7};
  int sourceConn[6]={0,3,1,0,2,3};
  MEDCouplingUMesh *sourceMesh=MEDCouplingUMesh::New();
  sourceMesh->setMeshDimension(2);
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
  std::copy(sourceCoords,sourceCoords+12,myCoords->getPointer());
  sourceMesh->setCoords(myCoords);
  myCoords->decrRef();
  return sourceMesh;
}

double MEDCouplingBasicsTest::sumAll(const std::vector< std::map<int,double> >& matrix)
{
  double ret=0.;
  for(std::vector< std::map<int,double> >::const_iterator iter=matrix.begin();iter!=matrix.end();iter++)
    for(std::map<int,double>::const_iterator iter2=(*iter).begin();iter2!=(*iter).end();iter2++)
      ret+=(*iter2).second;
  return ret;
}
