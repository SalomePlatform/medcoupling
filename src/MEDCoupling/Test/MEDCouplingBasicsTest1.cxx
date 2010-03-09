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
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"

#include <cmath>
#include <functional>

using namespace std;
using namespace ParaMEDMEM;

void MEDCouplingBasicsTest::testArray()
{
  int tmp1[6]={7,6,5,4,3,2};
  const int tmp2[3]={8,9,10};
  {
    MemArray<int> mem;
    mem.useArray(tmp1,false,CPP_DEALLOC,6);
    CPPUNIT_ASSERT(tmp1==mem.getConstPointer());
    CPPUNIT_ASSERT_THROW(mem.getPointer(),INTERP_KERNEL::Exception);
    CPPUNIT_ASSERT_THROW(mem[2]=7,INTERP_KERNEL::Exception);
    CPPUNIT_ASSERT_THROW(mem.writeOnPlace(0,12,tmp2,3),INTERP_KERNEL::Exception);
    mem.writeOnPlace(4,12,tmp2,3);
  }
  {
    int *tmp3=new int[6];
    std::copy(tmp1,tmp1+6,tmp3);
    MemArray<int> mem2;
    mem2.useArray(tmp3,true,CPP_DEALLOC,6);
    CPPUNIT_ASSERT(tmp3==mem2.getConstPointer());
    CPPUNIT_ASSERT(tmp3==mem2.getPointer());
    CPPUNIT_ASSERT_EQUAL(5,mem2[2]);
    mem2[2]=7;
    CPPUNIT_ASSERT_EQUAL(7,mem2[2]);
    mem2.writeOnPlace(0,12,tmp2,3);
    CPPUNIT_ASSERT_EQUAL(9,mem2[2]);
    CPPUNIT_ASSERT_EQUAL(12,mem2[0]);
    mem2.writeOnPlace(4,12,tmp2,3);
  }
}

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
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  // test clone not recursively
  MEDCouplingUMesh *mesh2=mesh->clone(false);
  CPPUNIT_ASSERT(mesh2!=mesh);
  mesh2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(nbOfCells,mesh2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(nbOfNodes,mesh2->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,mesh2->getSpaceDimension());
  CPPUNIT_ASSERT(mesh!=mesh2);
  CPPUNIT_ASSERT(mesh->getCoords()==mesh2->getCoords());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.2863,mesh2->getCoords()->getIJ(11,2),1e-14);
  CPPUNIT_ASSERT(mesh->getNodalConnectivity()==mesh2->getNodalConnectivity());
  CPPUNIT_ASSERT_EQUAL(3,mesh2->getNodalConnectivity()->getIJ(7,0));
  CPPUNIT_ASSERT(mesh->getNodalConnectivityIndex()==mesh2->getNodalConnectivityIndex());
  CPPUNIT_ASSERT_EQUAL(15,mesh2->getNodalConnectivityIndex()->getIJ(3,0));
  mesh2->decrRef();
  // test clone not recursively
  MEDCouplingUMesh *mesh3=mesh->clone(true);
  CPPUNIT_ASSERT(mesh3!=mesh);
  mesh3->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(nbOfCells,mesh3->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(nbOfNodes,mesh3->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,mesh3->getSpaceDimension());
  CPPUNIT_ASSERT(mesh!=mesh3);
  CPPUNIT_ASSERT(mesh->getCoords()!=mesh3->getCoords());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.2863,mesh3->getCoords()->getIJ(11,2),1e-14);
  CPPUNIT_ASSERT(mesh->getNodalConnectivity()!=mesh3->getNodalConnectivity());
  CPPUNIT_ASSERT_EQUAL(3,mesh3->getNodalConnectivity()->getIJ(7,0));
  CPPUNIT_ASSERT(mesh->getNodalConnectivityIndex()!=mesh3->getNodalConnectivityIndex());
  CPPUNIT_ASSERT_EQUAL(15,mesh3->getNodalConnectivityIndex()->getIJ(3,0));
  mesh3->decrRef();
  //test 4 - Field on cells
  MEDCouplingFieldDouble *fieldOnCells=MEDCouplingFieldDouble::New(ON_CELLS);
  fieldOnCells->setMesh(mesh);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(nbOfCells,9);
  fieldOnCells->setArray(array);
  tmp=array->getPointer();
  array->decrRef();
  fill(tmp,tmp+9*nbOfCells,7.);
  //content of field changed -> declare it.
  fieldOnCells->declareAsNew();
  fieldOnCells->checkCoherency();
  // testing clone of fields - no recursive
  MEDCouplingFieldDouble *fieldOnCells2=fieldOnCells->clone(false);
  CPPUNIT_ASSERT(fieldOnCells2!=fieldOnCells);
  fieldOnCells2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(nbOfCells,fieldOnCells2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(9,fieldOnCells2->getNumberOfComponents());
  CPPUNIT_ASSERT(fieldOnCells2->getArray()==fieldOnCells->getArray());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,fieldOnCells2->getArray()->getIJ(3,7),1e-14);
  CPPUNIT_ASSERT(fieldOnCells2->getMesh()==fieldOnCells->getMesh());
  // testing clone of fields - recursive
  MEDCouplingFieldDouble *fieldOnCells3=fieldOnCells->clone(true);
  CPPUNIT_ASSERT(fieldOnCells3!=fieldOnCells);
  fieldOnCells3->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(nbOfCells,fieldOnCells3->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(9,fieldOnCells3->getNumberOfComponents());
  CPPUNIT_ASSERT(fieldOnCells3->getArray()!=fieldOnCells->getArray());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,fieldOnCells3->getArray()->getIJ(3,7),1e-14);
  CPPUNIT_ASSERT(fieldOnCells3->getMesh()==fieldOnCells->getMesh());
  fieldOnCells2->decrRef();
  fieldOnCells3->decrRef();
  //
  fieldOnCells->decrRef();
  //clean-up
  mesh->decrRef();
}

void MEDCouplingBasicsTest::testMeshPointsCloud()
{
  double targetCoords[27]={-0.3,-0.3,0.5, 0.2,-0.3,1., 0.7,-0.3,1.5, -0.3,0.2,0.5, 0.2,0.2,1., 0.7,0.2,1.5, -0.3,0.7,0.5, 0.2,0.7,1., 0.7,0.7,1.5};
  int *targetConn=0;
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(0);
  targetMesh->allocateCells(8);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT0,0,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT0,0,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT0,0,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT0,0,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT0,0,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT0,0,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT0,0,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT0,0,targetConn);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(9,3);
  std::copy(targetCoords,targetCoords+27,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  //
  targetMesh->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(3,targetMesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(8,targetMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(9,targetMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(0,targetMesh->getMeshDimension());
  //
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::testMeshM1D()
{
  MEDCouplingUMesh *meshM1D=MEDCouplingUMesh::New();
  CPPUNIT_ASSERT_THROW(meshM1D->getMeshDimension(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(meshM1D->getNumberOfNodes(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(meshM1D->getNumberOfCells(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(meshM1D->setMeshDimension(-2),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(meshM1D->setMeshDimension(-10),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(meshM1D->checkCoherency(),INTERP_KERNEL::Exception);
  meshM1D->setMeshDimension(-1);
  meshM1D->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(-1,meshM1D->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(1,meshM1D->getNumberOfCells());
  CPPUNIT_ASSERT_THROW(meshM1D->getNumberOfNodes(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(meshM1D->getSpaceDimension(),INTERP_KERNEL::Exception);
  MEDCouplingUMesh *cpy=meshM1D->clone(true);
  CPPUNIT_ASSERT(cpy->isEqual(meshM1D,1e-12));
  cpy->decrRef();
  MEDCouplingFieldDouble *fieldOnCells=MEDCouplingFieldDouble::New(ON_CELLS);
  fieldOnCells->setMesh(meshM1D);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(1,6);
  fieldOnCells->setArray(array);
  double *tmp=array->getPointer();
  array->decrRef();
  fill(tmp,tmp+6,7.);
  fieldOnCells->checkCoherency();
  //
  fieldOnCells->decrRef();
  meshM1D->decrRef();
}

void MEDCouplingBasicsTest::testDeepCopy()
{
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(5,3);
  fill(array->getPointer(),array->getPointer()+5*3,7.);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,array->getIJ(3,2),1e-14);
  double *tmp1=array->getPointer();
  DataArrayDouble *array2=array->deepCopy();
  double *tmp2=array2->getPointer();
  CPPUNIT_ASSERT(tmp1!=tmp2);
  array->decrRef();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,array2->getIJ(3,2),1e-14);
  array2->decrRef();
  //
  DataArrayInt *array3=DataArrayInt::New();
  array3->alloc(5,3);
  fill(array3->getPointer(),array3->getPointer()+5*3,17);
  CPPUNIT_ASSERT_EQUAL(17,array3->getIJ(3,2));
  int *tmp3=array3->getPointer();
  DataArrayInt *array4=array3->deepCopy();
  int *tmp4=array4->getPointer();
  CPPUNIT_ASSERT(tmp3!=tmp4);
  array3->decrRef();
  CPPUNIT_ASSERT_EQUAL(17,array4->getIJ(3,2));
  array4->decrRef();
}

void MEDCouplingBasicsTest::testRevNodal()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  DataArrayInt *revNodal=DataArrayInt::New();
  DataArrayInt *revNodalIndx=DataArrayInt::New();
  //
  mesh->getReverseNodalConnectivity(revNodal,revNodalIndx);
  const int revNodalExpected[18]={0,0,1,1,2,0,3,0,1,2,3,4,2,4,3,3,4,4};
  const int revNodalIndexExpected[10]={0,1,3,5,7,12,14,15,17,18};
  CPPUNIT_ASSERT_EQUAL(18,revNodal->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(10,revNodalIndx->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(revNodalExpected,revNodalExpected+18,revNodal->getPointer()));
  CPPUNIT_ASSERT(std::equal(revNodalIndexExpected,revNodalIndexExpected+10,revNodalIndx->getPointer()));
  //
  revNodal->decrRef();
  revNodalIndx->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest::testConvertToPolyTypes()
{
  ////// 2D
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  //
  const int elts[2]={1,3};
  std::vector<int> eltsV(elts,elts+2);
  mesh->convertToPolyTypes(eltsV);
  mesh->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(5,mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(23,mesh->getNodalConnectivity()->getNumberOfTuples());
  const int *pt=mesh->getNodalConnectivity()->getConstPointer();
  const int expected1[23]={4, 0, 3, 4, 1, 5, 1, 4, 2, 3, 4, 5, 2, 5, 6, 7, 4, 3, 4, 7, 8, 5, 4};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+23,pt));
  //
  mesh->decrRef();
  ////// 3D
  mesh=build3DTargetMesh_1();
  mesh->convertToPolyTypes(eltsV);
  mesh->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(8,mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(114,mesh->getNodalConnectivity()->getNumberOfTuples());
  mesh->convertToPolyTypes(eltsV);
  mesh->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(8,mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(114,mesh->getNodalConnectivity()->getNumberOfTuples());
  //
  mesh->decrRef();
}

void MEDCouplingBasicsTest::testDescConn2D()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  DataArrayInt *desc=DataArrayInt::New();
  DataArrayInt *descIndx=DataArrayInt::New();
  DataArrayInt *revDesc=DataArrayInt::New();
  DataArrayInt *revDescIndx=DataArrayInt::New();
  //
  MEDCouplingUMesh *mesh2=mesh->buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  mesh2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(1,mesh2->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(13,mesh2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(14,revDescIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(14,revDescIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(6,descIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(6,descIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(18,desc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(18,desc->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(18,revDesc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(18,revDesc->getNumberOfTuples());
  const int expected1[18]={0,1,2,3, 2,4,5, 6,7,4, 8,9,1,10, 11,12,6,9};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+18,desc->getConstPointer()));
  const int expected2[6]={0,4,7,10,14,18};
  CPPUNIT_ASSERT(std::equal(expected2,expected2+6,descIndx->getConstPointer()));
  const int expected3[14]={0,1,3,5,6,8,9,11,12,13,15,16,17,18};
  CPPUNIT_ASSERT(std::equal(expected3,expected3+14,revDescIndx->getConstPointer()));
  const int expected4[18]={0, 0,3, 0,1, 0, 1,2, 1, 2,4, 2, 3, 3,4, 3, 4, 4};
  CPPUNIT_ASSERT(std::equal(expected4,expected4+18,revDesc->getConstPointer()));
  DataArrayInt *conn=mesh2->getNodalConnectivity();
  DataArrayInt *connIndex=mesh2->getNodalConnectivityIndex();
  const int expected5[14]={0,3,6,9,12,15,18,21,24,27,30,33,36,39};
  CPPUNIT_ASSERT(std::equal(expected5,expected5+14,connIndex->getConstPointer()));
  const int expected6[39]={1, 0, 3, 1, 3, 4, 1, 4, 1, 1, 1, 0, 1, 4, 2, 1, 2, 1, 1, 4, 5, 1, 5, 2, 1, 6, 7, 1, 7, 4, 1, 3, 6, 1, 7, 8, 1, 8, 5};
  CPPUNIT_ASSERT(std::equal(expected6,expected6+39,conn->getConstPointer()));
  //
  desc->decrRef();
  descIndx->decrRef();
  revDesc->decrRef();
  revDescIndx->decrRef();
  mesh2->decrRef();
  //
  const int elts[2]={1,3};
  std::vector<int> eltsV(elts,elts+2);
  mesh->convertToPolyTypes(eltsV);
  mesh->checkCoherency();
  //
  desc=DataArrayInt::New();
  descIndx=DataArrayInt::New();
  revDesc=DataArrayInt::New();
  revDescIndx=DataArrayInt::New();
  //
  mesh2=mesh->buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  mesh2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(1,mesh2->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(13,mesh2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(14,revDescIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(14,revDescIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(6,descIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(6,descIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(18,desc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(18,desc->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(18,revDesc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(18,revDesc->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+18,desc->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+6,descIndx->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected3,expected3+14,revDescIndx->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+18,revDesc->getConstPointer()));
  conn=mesh2->getNodalConnectivity();
  connIndex=mesh2->getNodalConnectivityIndex();
  CPPUNIT_ASSERT(std::equal(expected5,expected5+14,connIndex->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected6,expected6+39,conn->getConstPointer()));
  //
  desc->decrRef();
  descIndx->decrRef();
  revDesc->decrRef();
  revDescIndx->decrRef();
  mesh2->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest::testDescConn3D()
{
  MEDCouplingUMesh *mesh=build3DTargetMesh_1();
  DataArrayInt *desc=DataArrayInt::New();
  DataArrayInt *descIndx=DataArrayInt::New();
  DataArrayInt *revDesc=DataArrayInt::New();
  DataArrayInt *revDescIndx=DataArrayInt::New();
  //
  MEDCouplingUMesh *mesh2=mesh->buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  mesh2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(2,mesh2->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(36,mesh2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(37,revDescIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(37,revDescIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(9,descIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(9,descIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(48,desc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(48,desc->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(48,revDesc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(48,revDesc->getNumberOfTuples());
  const int expected1[9]={0, 6, 12, 18, 24, 30, 36, 42, 48};
  const int expected2[48]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 3, 11, 12, 4, 13, 14, 15, 16, 17, 10, 18, 19, 13, 1, 20, 21, 22, 23, 24, 7, 25, 26, 27, 28, 22, 12, 29, 23, 30, 31, 32, 17, 33, 28, 34, 35, 30};
  const int expected3[37]={0, 1, 3, 4, 6, 8, 9, 10, 12, 13, 14, 16, 17, 19, 21, 22, 23, 24, 26, 27, 28, 29, 30, 32, 34, 35, 36, 37, 38, 40, 41, 43, 44, 45, 46, 47, 48};
  const int expected4[48]={0, 0, 4, 0, 0, 1, 0, 2, 0, 1, 1, 5, 1, 1, 1, 3, 2, 2, 6, 2, 3, 2, 2, 3, 3, 7, 3, 3, 4, 4, 4, 5, 4, 6, 4, 5, 5, 5, 5, 7, 6, 6, 7, 6, 6, 7, 7, 7};
  const int expected5[37]={0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180};
  const int expected6[180]={4, 0, 1, 4, 3, 4, 9, 12, 13, 10, 4, 0, 9, 10, 1, 4, 1, 10, 13, 4, 4, 4, 13, 12, 3, 4, 3, 12, 9, 0, 4, 1, 2, 5, 4, 4, 10, 13, 14, 11, 4, 1, 10, 11, 2, 4, 2, 11, 14,
                            5, 4, 5, 14, 13, 4, 4, 3, 4, 7, 6, 4, 12, 15, 16, 13, 4, 4, 13, 16, 7, 4, 7, 16, 15, 6, 4, 6, 15, 12, 3, 4, 4, 5, 8, 7, 4, 13, 16, 17, 14, 4, 5, 14, 17, 8, 4, 8,
                            17, 16, 7, 4, 18, 21, 22, 19, 4, 9, 18, 19, 10, 4, 10, 19, 22, 13, 4, 13, 22, 21, 12, 4, 12, 21, 18, 9, 4, 19, 22, 23, 20, 4, 10, 19, 20, 11, 4, 11, 20, 23, 14, 4,
                            14, 23, 22, 13, 4, 21, 24, 25, 22, 4, 13, 22, 25, 16, 4, 16, 25, 24, 15, 4, 15, 24, 21, 12, 4, 22, 25, 26, 23, 4, 14, 23, 26, 17, 4, 17, 26, 25, 16};
  const int expected7[180]={4, 0, 1, 4, 3, 4, 9, 12, 13, 10, 4, 0, 9, 10, 1, 4, 1, 10, 13, 4, 4, 4, 13, 12, 3, 4, 3, 12, 9, 0, 5, 1, 2, 5, 4, 5, 10, 13, 14, 11, 5, 1, 10, 11, 2, 5, 2, 11, 14,
                            5, 5, 5, 14, 13, 4, 4, 3, 4, 7, 6, 4, 12, 15, 16, 13, 4, 4, 13, 16, 7, 4, 7, 16, 15, 6, 4, 6, 15, 12, 3, 5, 4, 5, 8, 7, 5, 13, 16, 17, 14, 5, 5, 14, 17, 8, 5, 8,
                            17, 16, 7, 4, 18, 21, 22, 19, 4, 9, 18, 19, 10, 4, 10, 19, 22, 13, 4, 13, 22, 21, 12, 4, 12, 21, 18, 9, 4, 19, 22, 23, 20, 4, 10, 19, 20, 11, 4, 11, 20, 23, 14, 4,
                            14, 23, 22, 13, 4, 21, 24, 25, 22, 4, 13, 22, 25, 16, 4, 16, 25, 24, 15, 4, 15, 24, 21, 12, 4, 22, 25, 26, 23, 4, 14, 23, 26, 17, 4, 17, 26, 25, 16};

  CPPUNIT_ASSERT(std::equal(expected1,expected1+9,descIndx->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+48,desc->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected3,expected3+37,revDescIndx->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+48,revDesc->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected5,expected5+37,mesh2->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected6,expected6+180,mesh2->getNodalConnectivity()->getConstPointer()));
  //
  desc->decrRef();
  descIndx->decrRef();
  revDesc->decrRef();
  revDescIndx->decrRef();
  mesh2->decrRef();
  //
  const int elts[2]={1,3};
  std::vector<int> eltsV(elts,elts+2);
  mesh->convertToPolyTypes(eltsV);
  mesh->checkCoherency();
  desc=DataArrayInt::New();
  descIndx=DataArrayInt::New();
  revDesc=DataArrayInt::New();
  revDescIndx=DataArrayInt::New();
  mesh2=mesh->buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  mesh2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(2,mesh2->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(36,mesh2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(37,revDescIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(37,revDescIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(9,descIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(9,descIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(48,desc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(48,desc->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(48,revDesc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(48,revDesc->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+9,descIndx->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+48,desc->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected3,expected3+37,revDescIndx->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+48,revDesc->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected5,expected5+37,mesh2->getNodalConnectivityIndex()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected7,expected7+180,mesh2->getNodalConnectivity()->getConstPointer()));
  //
  desc->decrRef();
  descIndx->decrRef();
  revDesc->decrRef();
  revDescIndx->decrRef();
  mesh2->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest::testFindBoundaryNodes()
{
  MEDCouplingUMesh *mesh=build3DTargetMesh_1();
  std::vector<int> boundaryNodes;
  mesh->findBoundaryNodes(boundaryNodes);
  CPPUNIT_ASSERT_EQUAL(26,(int)boundaryNodes.size());
  const int expected1[26]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+26,boundaryNodes.begin()));
  mesh->decrRef();
}

void MEDCouplingBasicsTest::testBoundaryMesh()
{
  MEDCouplingUMesh *mesh=build3DTargetMesh_1();
  MEDCouplingPointSet *mesh2=mesh->buildBoundaryMesh(false);
  CPPUNIT_ASSERT_EQUAL(24,mesh2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(26,mesh2->getNumberOfNodes());
  mesh2->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest::testBuildPartOfMySelf()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  mesh->setName("Toto");
  const int tab1[2]={0,4};
  const int tab2[3]={0,2,3};
  //
  MEDCouplingPointSet *subMeshSimple=mesh->buildPartOfMySelf(tab1,tab1+2,true);
  MEDCouplingUMesh *subMesh=dynamic_cast<MEDCouplingUMesh *>(subMeshSimple);
  CPPUNIT_ASSERT(subMesh);
  std::string name(subMesh->getName());
  CPPUNIT_ASSERT_EQUAL(2,(int)mesh->getAllTypes().size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,*mesh->getAllTypes().begin());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*(++(mesh->getAllTypes().begin())));
  CPPUNIT_ASSERT_EQUAL(1,(int)subMesh->getAllTypes().size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*subMesh->getAllTypes().begin());
  CPPUNIT_ASSERT(name=="PartOf_Toto");
  CPPUNIT_ASSERT(mesh->getCoords()==subMesh->getCoords());
  CPPUNIT_ASSERT_EQUAL(2,subMesh->getNumberOfCells());
  const int subConn[10]={4,0,3,4,1,4,7,8,5,4};
  const int subConnIndex[3]={0,5,10};
  CPPUNIT_ASSERT_EQUAL(10,subMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(3,subMesh->getNodalConnectivityIndex()->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(subConn,subConn+10,subMesh->getNodalConnectivity()->getPointer()));
  CPPUNIT_ASSERT(std::equal(subConnIndex,subConnIndex+3,subMesh->getNodalConnectivityIndex()->getPointer()));
  subMesh->decrRef();
  //
  subMeshSimple=mesh->buildPartOfMySelf(tab2,tab2+3,true);
  subMesh=dynamic_cast<MEDCouplingUMesh *>(subMeshSimple);
  CPPUNIT_ASSERT(subMesh);
  name=subMesh->getName();
  CPPUNIT_ASSERT_EQUAL(2,(int)subMesh->getAllTypes().size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,*subMesh->getAllTypes().begin());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*(++(subMesh->getAllTypes().begin())));
  CPPUNIT_ASSERT(name=="PartOf_Toto");
  CPPUNIT_ASSERT(mesh->getCoords()==subMesh->getCoords());
  CPPUNIT_ASSERT_EQUAL(3,subMesh->getNumberOfCells());
  const int subConn2[14]={4,0,3,4,1,3,4,5,2,4,6,7,4,3};
  const int subConnIndex2[4]={0,5,9,14};
  CPPUNIT_ASSERT_EQUAL(14,subMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(4,subMesh->getNodalConnectivityIndex()->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(subConn2,subConn2+14,subMesh->getNodalConnectivity()->getPointer()));
  CPPUNIT_ASSERT(std::equal(subConnIndex2,subConnIndex2+4,subMesh->getNodalConnectivityIndex()->getPointer()));
  subMesh->decrRef();
  //
  mesh->decrRef();
}

void MEDCouplingBasicsTest::testBuildPartOfMySelfNode()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  const int tab1[2]={5,7};
  MEDCouplingPointSet *subMeshSimple=mesh->buildPartOfMySelfNode(tab1,tab1+2,true);
  MEDCouplingUMesh *subMesh=dynamic_cast<MEDCouplingUMesh *>(subMeshSimple);
  CPPUNIT_ASSERT(subMesh);
  CPPUNIT_ASSERT_EQUAL(1,(int)subMesh->getAllTypes().size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*subMesh->getAllTypes().begin());
  CPPUNIT_ASSERT_EQUAL(1,subMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(5,subMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(2,subMesh->getNodalConnectivityIndex()->getNbOfElems());
  const int subConn[5]={4,7,8,5,4};
  const int subConnIndex[3]={0,5};
  CPPUNIT_ASSERT(std::equal(subConn,subConn+5,subMesh->getNodalConnectivity()->getPointer()));
  CPPUNIT_ASSERT(std::equal(subConnIndex,subConnIndex+2,subMesh->getNodalConnectivityIndex()->getPointer()));
  CPPUNIT_ASSERT(subMesh->getCoords()==mesh->getCoords());
  subMeshSimple->decrRef();
  //
  subMeshSimple=mesh->buildPartOfMySelfNode(tab1,tab1+2,false);
  subMesh=dynamic_cast<MEDCouplingUMesh *>(subMeshSimple);
  CPPUNIT_ASSERT(subMesh);
  CPPUNIT_ASSERT_EQUAL(2,(int)subMesh->getAllTypes().size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,*subMesh->getAllTypes().begin());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*(++subMesh->getAllTypes().begin()));
  CPPUNIT_ASSERT_EQUAL(3,subMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(14,subMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(4,subMesh->getNodalConnectivityIndex()->getNbOfElems());
  const int subConn2[14]={3,4,5,2,4,6,7,4,3,4,7,8,5,4};
  const int subConnIndex2[4]={0,4,9,14};
  CPPUNIT_ASSERT(std::equal(subConn2,subConn2+14,subMesh->getNodalConnectivity()->getPointer()));
  CPPUNIT_ASSERT(std::equal(subConnIndex2,subConnIndex2+4,subMesh->getNodalConnectivityIndex()->getPointer()));
  CPPUNIT_ASSERT(subMesh->getCoords()==mesh->getCoords());
  subMeshSimple->decrRef();
  //testing the case where length of tab2 is greater than max number of node per cell.
  const int tab2[7]={0,3,2,1,4,5,6};
  subMeshSimple=mesh->buildPartOfMySelfNode(tab2,tab2+7,true);
  subMesh=dynamic_cast<MEDCouplingUMesh *>(subMeshSimple);
  CPPUNIT_ASSERT(subMesh);
  CPPUNIT_ASSERT_EQUAL(2,(int)subMesh->getAllTypes().size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,*subMesh->getAllTypes().begin());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*(++subMesh->getAllTypes().begin()));
  CPPUNIT_ASSERT_EQUAL(3,subMesh->getNumberOfCells());
  subMeshSimple->decrRef();
  //
  mesh->decrRef();
}

void MEDCouplingBasicsTest::testZipCoords()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  CPPUNIT_ASSERT_EQUAL(2,(int)mesh->getAllTypes().size());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(9,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(5,mesh->getNumberOfCells());
  std::vector<int> oldConn(mesh->getNodalConnectivity()->getNbOfElems());
  std::vector<int> oldConnIndex(mesh->getNumberOfCells()+1);
  std::copy(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+oldConn.size(),oldConn.begin());
  std::copy(mesh->getNodalConnectivityIndex()->getPointer(),mesh->getNodalConnectivityIndex()->getPointer()+mesh->getNumberOfCells()+1,oldConnIndex.begin());
  DataArrayDouble *oldCoords=mesh->getCoords();
  oldCoords->incrRef();
  mesh->zipCoords();
  CPPUNIT_ASSERT_EQUAL(2,(int)mesh->getAllTypes().size());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(9,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(5,mesh->getNumberOfCells());
  CPPUNIT_ASSERT(mesh->getCoords()!=oldCoords);
  CPPUNIT_ASSERT(std::equal(mesh->getCoords()->getPointer(),mesh->getCoords()->getPointer()+2*9,oldCoords->getPointer()));
  CPPUNIT_ASSERT(std::equal(oldConn.begin(),oldConn.end(),mesh->getNodalConnectivity()->getPointer()));
  CPPUNIT_ASSERT(std::equal(oldConnIndex.begin(),oldConnIndex.end(),mesh->getNodalConnectivityIndex()->getPointer()));
  oldCoords->decrRef();
  //
  const int tab1[2]={0,4};
  MEDCouplingPointSet *subMeshPtSet=mesh->buildPartOfMySelf(tab1,tab1+2,true);
  MEDCouplingUMesh *subMesh=dynamic_cast<MEDCouplingUMesh *>(subMeshPtSet);
  CPPUNIT_ASSERT(subMesh);
  DataArrayInt *traducer=subMesh->zipCoordsTraducer();
  const int expectedTraducer[9]={0,1,-1,2,3,4,-1,5,6};
  CPPUNIT_ASSERT(std::equal(expectedTraducer,expectedTraducer+9,traducer->getPointer()));
  traducer->decrRef();
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*subMesh->getAllTypes().begin());
  CPPUNIT_ASSERT_EQUAL(2,subMesh->getNumberOfCells());
  const int subConn[10]={4,0,2,3,1,4,5,6,4,3};
  const int subConnIndex[3]={0,5,10};
  CPPUNIT_ASSERT_EQUAL(7,subMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(10,subMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(3,subMesh->getNodalConnectivityIndex()->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(subConn,subConn+10,subMesh->getNodalConnectivity()->getPointer()));
  CPPUNIT_ASSERT(std::equal(subConnIndex,subConnIndex+3,subMesh->getNodalConnectivityIndex()->getPointer()));
  subMesh->decrRef();
  //
  subMeshPtSet=mesh->buildPartOfMySelf(tab1,tab1+2,false);
  subMesh=dynamic_cast<MEDCouplingUMesh *>(subMeshPtSet);
  CPPUNIT_ASSERT(subMesh);
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*subMesh->getAllTypes().begin());
  CPPUNIT_ASSERT_EQUAL(2,subMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(7,subMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(10,subMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(3,subMesh->getNodalConnectivityIndex()->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(subConn,subConn+10,subMesh->getNodalConnectivity()->getPointer()));
  CPPUNIT_ASSERT(std::equal(subConnIndex,subConnIndex+3,subMesh->getNodalConnectivityIndex()->getPointer()));
  subMesh->decrRef();
  //
  mesh->decrRef();
}

void MEDCouplingBasicsTest::testEqualMesh()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  MEDCouplingUMesh *mesh2=build2DTargetMesh_1();
  //
  CPPUNIT_ASSERT(mesh1->isEqual(mesh1,1e-12));
  //
  CPPUNIT_ASSERT(mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(mesh2->isEqual(mesh1,1e-12));
  double *pt=mesh2->getCoords()->getPointer();
  double tmp=pt[1];
  pt[1]=5.999;
  CPPUNIT_ASSERT(!mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(!mesh2->isEqual(mesh1,1e-12));
  pt[1]=tmp;
  CPPUNIT_ASSERT(mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(mesh2->isEqual(mesh1,1e-12));
  //
  int *pt2=mesh1->getNodalConnectivity()->getPointer();
  pt2[5]++;
  CPPUNIT_ASSERT(!mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(!mesh2->isEqual(mesh1,1e-12));
  pt2[5]--;
  CPPUNIT_ASSERT(mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(mesh2->isEqual(mesh1,1e-12));
  //
  pt2=mesh1->getNodalConnectivityIndex()->getPointer();
  pt2[1]++;
  CPPUNIT_ASSERT(!mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(!mesh2->isEqual(mesh1,1e-12));
  pt2[1]--;
  CPPUNIT_ASSERT(mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(mesh2->isEqual(mesh1,1e-12));
  //
  std::string tmp3=mesh1->getName();
  mesh1->setName("lllll");
  CPPUNIT_ASSERT(!mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(!mesh2->isEqual(mesh1,1e-12));
  mesh1->setName(tmp3.c_str());
  CPPUNIT_ASSERT(mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(mesh2->isEqual(mesh1,1e-12));
  //
  tmp3=mesh2->getCoords()->getInfoOnComponent(1);
  mesh2->getCoords()->setInfoOnComponent(1,"kkkkkk");
  CPPUNIT_ASSERT(!mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(!mesh2->isEqual(mesh1,1e-12));
  mesh2->getCoords()->setInfoOnComponent(1,tmp3.c_str());
  CPPUNIT_ASSERT(mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(mesh2->isEqual(mesh1,1e-12));
  //
  mesh1->decrRef();
  mesh2->decrRef();
}

void MEDCouplingBasicsTest::testEqualFieldDouble()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  MEDCouplingUMesh *mesh2=build2DTargetMesh_1();
  //
  MEDCouplingFieldDouble *fieldOnCells1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  fieldOnCells1->setMesh(mesh1);
  MEDCouplingFieldDouble *fieldOnCells2=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  fieldOnCells2->setMesh(mesh2);
  //
  CPPUNIT_ASSERT(fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  fieldOnCells2->decrRef();
  //
  MEDCouplingFieldDouble *fieldOnNodes1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  CPPUNIT_ASSERT(!fieldOnCells1->isEqual(fieldOnNodes1,1e-12,1e-15));
  CPPUNIT_ASSERT(!fieldOnNodes1->isEqual(fieldOnCells1,1e-12,1e-15));
  fieldOnNodes1->decrRef();
  //
  fieldOnCells2=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  CPPUNIT_ASSERT(!fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(!fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  fieldOnCells1->decrRef();
  fieldOnCells1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  CPPUNIT_ASSERT(fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  fieldOnCells1->setTime(4.,6,7);
  CPPUNIT_ASSERT(!fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(!fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  fieldOnCells2->setTime(4.,6,7);
  CPPUNIT_ASSERT(fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  fieldOnCells1->setName("Power");
  CPPUNIT_ASSERT(!fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(!fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  fieldOnCells2->setName("Power");
  CPPUNIT_ASSERT(fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  //
  fieldOnCells1->setMesh(mesh1);
  CPPUNIT_ASSERT(!fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(!fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  fieldOnCells2->setMesh(mesh1);
  CPPUNIT_ASSERT(fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  DataArrayDouble *arr=DataArrayDouble::New();
  arr->setName("popo");
  arr->alloc(mesh1->getNumberOfCells(),3);
  double *pt=arr->getPointer();
  std::fill(pt,pt+mesh1->getNumberOfCells()*3,6.);
  fieldOnCells1->setArray(arr);
  CPPUNIT_ASSERT(!fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(!fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  fieldOnCells2->setArray(arr);
  arr->decrRef();
  CPPUNIT_ASSERT(fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  //
  DataArrayDouble *arr2=arr->deepCopy();
  fieldOnCells2->setArray(arr2);
  arr2->decrRef();
  CPPUNIT_ASSERT(fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  pt[4]=6.1;
  CPPUNIT_ASSERT(!fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(!fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  pt[4]=6.;
  CPPUNIT_ASSERT(fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  arr2->setName("popo2");
  CPPUNIT_ASSERT(!fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(!fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  //
  arr2->setName("popo");
  CPPUNIT_ASSERT(fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  //
  arr2->setInfoOnComponent(2,"jjj");
  CPPUNIT_ASSERT(!fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(!fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  arr->setInfoOnComponent(2,"jjj");
  CPPUNIT_ASSERT(fieldOnCells1->isEqual(fieldOnCells2,1e-12,1e-15));
  CPPUNIT_ASSERT(fieldOnCells2->isEqual(fieldOnCells1,1e-12,1e-15));
  //
  fieldOnCells1->decrRef();
  fieldOnCells2->decrRef();
  //
  mesh1->decrRef();
  mesh2->decrRef();
}

void MEDCouplingBasicsTest::testNatureChecking()
{
  MEDCouplingFieldDouble *field=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  field->setNature(Integral);
  field->setNature(ConservativeVolumic);
  field->setNature(IntegralGlobConstraint);
  field->decrRef();
  field=MEDCouplingFieldDouble::New(ON_NODES,NO_TIME);
  field->setNature(ConservativeVolumic);
  CPPUNIT_ASSERT_THROW(field->setNature(Integral),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(field->setNature(IntegralGlobConstraint),INTERP_KERNEL::Exception);
  field->decrRef();
}

void MEDCouplingBasicsTest::testBuildSubMeshData()
{
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //check buildSubMesh on field on cells
  MEDCouplingFieldDouble *fieldCells=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  fieldCells->setMesh(targetMesh);
  const int elts[3]={1,2,4};
  DataArrayInt *di;
  MEDCouplingMesh *ret1=fieldCells->buildSubMeshData(elts,elts+3,di);
  CPPUNIT_ASSERT_EQUAL(3,ret1->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(6,ret1->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,di->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,di->getNumberOfComponents());
  const int *toCheck=di->getConstPointer();
  CPPUNIT_ASSERT(std::equal(elts,elts+3,toCheck));
  MEDCouplingUMesh *ret1DC=dynamic_cast<MEDCouplingUMesh *>(ret1);
  CPPUNIT_ASSERT(ret1DC);
  ret1->decrRef();
  di->decrRef();
  fieldCells->decrRef();
  //check buildSubMesh on field on nodes
  MEDCouplingFieldDouble *fieldNodes=MEDCouplingFieldDouble::New(ON_NODES,NO_TIME);
  fieldNodes->setMesh(targetMesh);
  MEDCouplingMesh *ret2=fieldNodes->buildSubMeshData(elts,elts+3,di);
  MEDCouplingUMesh *ret2DC=dynamic_cast<MEDCouplingUMesh *>(ret2);
  CPPUNIT_ASSERT(ret2DC);
  CPPUNIT_ASSERT_EQUAL(3,ret2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(6,ret2->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(6,di->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,di->getNumberOfComponents());
  toCheck=di->getConstPointer();
  const int expected[6]={1,2,4,5,7,8};
  CPPUNIT_ASSERT(std::equal(expected,expected+6,toCheck));
  ret2->decrRef();
  di->decrRef();
  fieldNodes->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::testExtrudedMesh1()
{
  MEDCouplingUMesh *mesh2D=0;
  MEDCouplingUMesh *mesh3D=build3DExtrudedUMesh_1(mesh2D);
  MEDCouplingExtrudedMesh *ext=MEDCouplingExtrudedMesh::New(mesh3D,mesh2D,1);
  CPPUNIT_ASSERT_EQUAL(18,ext->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(60,ext->getNumberOfNodes());
  DataArrayInt *ids3D=ext->getMesh3DIds();
  const int ids3DExpected[18]={5,4,3,2,1,0, 11,10,9,8,7,6, 17,16,15,14,13,12};
  CPPUNIT_ASSERT_EQUAL(18,ids3D->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,ids3D->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(ids3DExpected,ids3DExpected+18,ids3D->getConstPointer()));
  MEDCouplingUMesh *mesh1D=ext->getMesh1D();
  CPPUNIT_ASSERT_EQUAL(4,mesh1D->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,mesh1D->getNumberOfCells());
  const double mesh1DExpected[12]={0.66666666666666663, 1.4583333333333333, 0, 0.66666666666666663, 1.4583333333333333, 1, 0.66666666666666663, 1.4583333333333333, 2, 0.66666666666666663, 1.4583333333333333, 3};
  DataArrayDouble *mesh1DCoords=mesh1D->getCoords();
  CPPUNIT_ASSERT_EQUAL(4,mesh1DCoords->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,mesh1DCoords->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(mesh1DExpected,mesh1DExpected+12,mesh1DCoords->getConstPointer()));
  DataArrayInt *conn1D=mesh1D->getNodalConnectivity();
  CPPUNIT_ASSERT_EQUAL(9,conn1D->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,conn1D->getNumberOfComponents());
  const int conn1DExpected[9]={1,0,1,1,1,2,1,2,3};
  CPPUNIT_ASSERT(std::equal(conn1DExpected,conn1DExpected+9,conn1D->getConstPointer()));
  ext->decrRef();
  mesh3D->decrRef();
  mesh2D->decrRef();
}

void MEDCouplingBasicsTest::testFindCommonNodes()
{
  DataArrayInt *comm,*commI;
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  targetMesh->findCommonNodes(comm,commI,1e-10);
  CPPUNIT_ASSERT_EQUAL(1,commI->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,comm->getNumberOfTuples());
  int newNbOfNodes;
  DataArrayInt *o2n=targetMesh->buildNewNumberingFromCommNodesFrmt(comm,commI,newNbOfNodes);
  CPPUNIT_ASSERT_EQUAL(27,newNbOfNodes);
  CPPUNIT_ASSERT_EQUAL(27,o2n->getNumberOfTuples());
  const int o2nExp1[27]=
    {
      0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
      21,22,23,24,25,26
    };
  CPPUNIT_ASSERT(std::equal(o2nExp1,o2nExp1+27,o2n->getConstPointer()));
  o2n->decrRef();
  comm->decrRef();
  commI->decrRef();
  targetMesh->decrRef();
  //
  targetMesh=build3DTargetMeshMergeNode_1();
  CPPUNIT_ASSERT_EQUAL(31,targetMesh->getNumberOfNodes());
  targetMesh->findCommonNodes(comm,commI,1e-10);
  CPPUNIT_ASSERT_EQUAL(3,commI->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(6,comm->getNumberOfTuples());
  const int commExpected[6]={1,27,28,29,23,30};
  const int commIExpected[3]={0,4,6};
  CPPUNIT_ASSERT(std::equal(commExpected,commExpected+6,comm->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(commIExpected,commIExpected+3,commI->getConstPointer()));
  o2n=targetMesh->buildNewNumberingFromCommNodesFrmt(comm,commI,newNbOfNodes);
  CPPUNIT_ASSERT_EQUAL(31,o2n->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(27,newNbOfNodes);
  const int o2nExp2[31]=
    {
      0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
      21,22,23,24,25,26,1,1,1,23
    };
  CPPUNIT_ASSERT(std::equal(o2nExp2,o2nExp2+31,o2n->getConstPointer()));
  o2n->decrRef();
  comm->decrRef();
  commI->decrRef();
  targetMesh->decrRef();
  //
  targetMesh=build3DTargetMesh_1();
  bool areNodesMerged;
  unsigned int time=targetMesh->getTimeOfThis();
  o2n=targetMesh->mergeNodes(1e-10,areNodesMerged);
  targetMesh->updateTime();
  CPPUNIT_ASSERT(time==targetMesh->getTimeOfThis());
  CPPUNIT_ASSERT(!areNodesMerged);
  targetMesh->decrRef();
  o2n->decrRef();
  //
  targetMesh=build3DTargetMeshMergeNode_1();
  time=targetMesh->getTimeOfThis();
  o2n=targetMesh->mergeNodes(1e-10,areNodesMerged);
  targetMesh->updateTime();
  CPPUNIT_ASSERT(time!=targetMesh->getTimeOfThis());
  CPPUNIT_ASSERT(areNodesMerged);
  int connExp[72]={18,0,1,4,3,9,10,13,12, 18,1,2,5,4,10,11,14,13, 18,3,4,7,6,12,13,16,15,
                   18,4,5,8,7,13,14,17,16,
                   18,9,10,13,12,18,19,22,21, 18,10,11,14,13,19,20,23,22, 18,12,13,16,15,21,22,25,24,
                   18,13,14,17,16,22,23,26,25};
  CPPUNIT_ASSERT_EQUAL(72,targetMesh->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(connExp,connExp+72,targetMesh->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(27,targetMesh->getCoords()->getNumberOfTuples());
  double coordsExp[81]={ 0., 0., 0., 50., 0., 0. , 200., 0., 0.  , 0., 50., 0., 50., 50., 0. ,
                         200., 50., 0.,   0., 200., 0., 50., 200., 0. , 200., 200., 0. ,
                         0., 0., 50., 50., 0., 50. , 200., 0., 50.  , 0., 50., 50., 50.,
                         50., 50. , 200., 50., 50.,   0., 200., 50., 50., 200., 50. ,
                         200., 200., 50. , 0., 0., 200., 50., 0., 200. , 200., 0., 200.  
                         , 0., 50., 200., 50., 50., 200. , 200., 50., 200., 
                         0., 200., 200., 50., 200., 200. , 200., 200., 200. };
  CPPUNIT_ASSERT(std::equal(coordsExp,coordsExp+81,targetMesh->getCoords()->getConstPointer()));
  targetMesh->decrRef();
  o2n->decrRef();
  //2D
  targetMesh=build2DTargetMeshMergeNode_1();
  CPPUNIT_ASSERT_EQUAL(18,targetMesh->getNumberOfNodes());
  time=targetMesh->getTimeOfThis();
  o2n=targetMesh->mergeNodes(1e-10,areNodesMerged);
  CPPUNIT_ASSERT(time!=targetMesh->getTimeOfThis());
  CPPUNIT_ASSERT(areNodesMerged);
  CPPUNIT_ASSERT_EQUAL(9,targetMesh->getNumberOfNodes());
  int connExp2[23]={4,0,4,3,1, 3,1,3,2, 3,3,5,2, 4,4,6,7,3, 4,7,8,5,3};
  CPPUNIT_ASSERT_EQUAL(23,targetMesh->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(connExp2,connExp2+23,targetMesh->getNodalConnectivity()->getConstPointer()));
  double coordsExp2[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, 0.2,0.2, -0.3,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7};
  CPPUNIT_ASSERT_EQUAL(9,targetMesh->getCoords()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(coordsExp2,coordsExp2+18,targetMesh->getCoords()->getConstPointer()));
  targetMesh->decrRef();
  o2n->decrRef();
}

void MEDCouplingBasicsTest::testCheckButterflyCells()
{
  std::vector<int> cells;
  MEDCouplingUMesh *sourceMesh=build2DTargetMesh_1();
  sourceMesh->checkButterflyCells(cells);
  CPPUNIT_ASSERT(cells.empty());
  int *pt=sourceMesh->getNodalConnectivity()->getPointer();
  std::swap(pt[15],pt[16]);
  sourceMesh->checkButterflyCells(cells);
  CPPUNIT_ASSERT_EQUAL(1,(int)cells.size());
  CPPUNIT_ASSERT_EQUAL(3,cells[0]);
  cells.clear();
  std::swap(pt[15],pt[16]);
  sourceMesh->checkButterflyCells(cells);
  CPPUNIT_ASSERT(cells.empty());
  sourceMesh->decrRef();
  // 3D surf
  sourceMesh=build3DSurfTargetMesh_1();
  sourceMesh->checkButterflyCells(cells);
  CPPUNIT_ASSERT(cells.empty());
  pt=sourceMesh->getNodalConnectivity()->getPointer();
  std::swap(pt[15],pt[16]);
  sourceMesh->checkButterflyCells(cells);
  CPPUNIT_ASSERT_EQUAL(1,(int)cells.size());
  CPPUNIT_ASSERT_EQUAL(3,cells[0]);
  cells.clear();
  std::swap(pt[15],pt[16]);
  sourceMesh->checkButterflyCells(cells);
  CPPUNIT_ASSERT(cells.empty());
  sourceMesh->decrRef();
}

void MEDCouplingBasicsTest::testMergeMesh1()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  MEDCouplingUMesh *m2=build2DSourceMesh_1();
  const double vec[2]={1.,0.};
  m2->translate(vec);
  MEDCouplingMesh *m3=m1->mergeMyselfWith(m2);
  MEDCouplingUMesh *m3C=dynamic_cast<MEDCouplingUMesh *>(m3);
  CPPUNIT_ASSERT(m3C);
  m3->checkCoherency();
  MEDCouplingUMesh *m4=build2DTargetMeshMerged_1();
  CPPUNIT_ASSERT(m3->isEqual(m4,1.e-12));
  m4->decrRef();
  bool isMerged;
  DataArrayInt *da=m3C->mergeNodes(1.e-12,isMerged);
  CPPUNIT_ASSERT_EQUAL(11,m3C->getNumberOfNodes());
  CPPUNIT_ASSERT(isMerged);
  da->decrRef();
  m3->decrRef();
  m1->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest::testMergeField1()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  MEDCouplingUMesh *m2=build2DSourceMesh_1();
  const double vec[2]={1.,0.};
  m2->translate(vec);
  MEDCouplingFieldDouble *f1=m1->getMeasureField(true);
  MEDCouplingFieldDouble *f2=m2->getMeasureField(true);
  MEDCouplingFieldDouble *f3=MEDCouplingFieldDouble::mergeFields(f1,f2);
  f3->checkCoherency();
  MEDCouplingUMesh *m4=build2DTargetMeshMerged_1();
  CPPUNIT_ASSERT(f3->getMesh()->isEqual(m4,1.e-12));
  std::string name=f3->getName();
  CPPUNIT_ASSERT(name=="MeasureOfMesh_");
  CPPUNIT_ASSERT(f3->getTypeOfField()==ON_CELLS);
  CPPUNIT_ASSERT(f3->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_EQUAL(1,f3->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(7,f3->getNumberOfTuples());
  double values[7]={0.25,0.125,0.125,0.25,0.25,0.5,0.5};
  const double *tmp=f3->getArray()->getConstPointer();
  std::transform(tmp,tmp+7,values,values,std::minus<double>());
  std::transform(values,values+7,values,std::ptr_fun<double,double>(fabs));
  double max=*std::max_element(values,values+7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  m4->decrRef();
  f3->decrRef();
  f1->decrRef();
  f2->decrRef();
  m1->decrRef();
  m2->decrRef();
}

bool func1(const double *pt, double *res);
bool func2(const double *pt, double *res);
bool func3(const double *pt, double *res);
bool func4(const double *pt, double *res);

bool func1(const double *pt, double *res)
{
  res[0]=pt[0]+pt[1];
  return true;
}

bool func2(const double *pt, double *res)
{
  res[0]=pt[0]+pt[1];
  res[1]=2.*(pt[0]+pt[1]);
  return true;
}

bool func3(const double *pt, double *res)
{
  if(fabs(pt[0]-0.2)<1e-12)
    return false;
  res[0]=1./(pt[0]-0.2);
  return true;
}

void MEDCouplingBasicsTest::testFillFromAnalytic()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=m->fillFromAnalytic(ON_CELLS,1,func1);
  f1->checkCoherency();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_CELLS);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfTuples());
  double values1[5]={-0.1,0.23333333333333336,0.56666666666666665,0.4,0.9};
  const double *tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+5,values1,values1,std::minus<double>());
  std::transform(values1,values1+5,values1,std::ptr_fun<double,double>(fabs));
  double max=*std::max_element(values1,values1+5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  //
  f1=m->fillFromAnalytic(ON_NODES,1,func1);
  f1->checkCoherency();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  double values2[9]={-0.6,-0.1,0.4,-0.1,0.4,0.9,0.4,0.9,1.4};
  tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+9,values2,values2,std::minus<double>());
  std::transform(values2,values2+9,values2,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values2,values2+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  //
  f1=m->fillFromAnalytic(ON_NODES,2,func2);
  f1->checkCoherency();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_EQUAL(2,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  double values3[18]={-0.6,-1.2,-0.1,-0.2,0.4,0.8,-0.1,-0.2,0.4,0.8,0.9,1.8,0.4,0.8,0.9,1.8,1.4,2.8};
  tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+18,values3,values3,std::minus<double>());
  std::transform(values3,values3+18,values3,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values3,values3+18);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  double values4[2];
  f1->accumulate(values4);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.6,values4[0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.2,values4[1],1.e-12);
  f1->measureAccumulate(true,values4);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,values4[0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,values4[1],1.e-12);
  f1->decrRef();
  //
  CPPUNIT_ASSERT_THROW(f1=m->fillFromAnalytic(ON_NODES,1,func3),INTERP_KERNEL::Exception);
  //
  m->decrRef();
}

void MEDCouplingBasicsTest::testFillFromAnalytic2()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=m->fillFromAnalytic(ON_CELLS,1,"y+x");
  f1->checkCoherency();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_CELLS);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfTuples());
  double values1[5]={-0.1,0.23333333333333336,0.56666666666666665,0.4,0.9};
  const double *tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+5,values1,values1,std::minus<double>());
  std::transform(values1,values1+5,values1,std::ptr_fun<double,double>(fabs));
  double max=*std::max_element(values1,values1+5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  //
  f1=m->fillFromAnalytic(ON_NODES,1,"y+2*x");
  f1->checkCoherency();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  double values2[9]={-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1};
  tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+9,values2,values2,std::minus<double>());
  std::transform(values2,values2+9,values2,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values2,values2+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  f1=m->fillFromAnalytic(ON_NODES,1,"2.*x+y");
  f1->checkCoherency();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  tmp=f1->getArray()->getConstPointer();
  double values2Bis[9]={-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1};
  std::transform(tmp,tmp+9,values2Bis,values2Bis,std::minus<double>());
  std::transform(values2,values2+9,values2Bis,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values2Bis,values2Bis+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  //
  f1=m->fillFromAnalytic(ON_NODES,2,"(x+y)*IVec+2*(x+y)*JVec");
  f1->checkCoherency();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_EQUAL(2,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  double values3[18]={-0.6,-1.2,-0.1,-0.2,0.4,0.8,-0.1,-0.2,0.4,0.8,0.9,1.8,0.4,0.8,0.9,1.8,1.4,2.8};
  tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+18,values3,values3,std::minus<double>());
  std::transform(values3,values3+18,values3,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values3,values3+18);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  double values4[2];
  f1->accumulate(values4);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.6,values4[0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.2,values4[1],1.e-12);
  f1->measureAccumulate(true,values4);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,values4[0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,values4[1],1.e-12);
  f1->decrRef();
  //
  CPPUNIT_ASSERT_THROW(f1=m->fillFromAnalytic(ON_NODES,1,"1./(x-0.2)"),INTERP_KERNEL::Exception);
  //
  m->decrRef();
}

void MEDCouplingBasicsTest::testApplyFunc()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=m->fillFromAnalytic(ON_NODES,2,func2);
  f1->checkCoherency();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_EQUAL(2,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  f1->applyFunc(1,func1);
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  double values1[9]={-1.8,-0.3,1.2,-0.3,1.2,2.7,1.2,2.7,4.2};
  const double *tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+9,values1,values1,std::minus<double>());
  std::transform(values1,values1+9,values1,std::ptr_fun<double,double>(fabs));
  double max=*std::max_element(values1,values1+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest::testApplyFunc2()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=m->fillFromAnalytic(ON_NODES,2,func2);
  f1->checkCoherency();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_EQUAL(2,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  //
  MEDCouplingFieldDouble *f2=f1->clone(true);
  f2->applyFunc("abs(u)^2.4+2*u");
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_EQUAL(2,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  double values2[18]={-0.9065304805418678, -0.85105859001709905, -0.19601892829446504, -0.37898777756476987,
                      0.91090317490482353, 2.1853504664669781, -0.19601892829446504, -0.37898777756476987,
                      0.91090317490482353, 2.1853504664669781, 2.5765725275664879, 7.6987743736515295,
                      0.91090317490482353, 2.1853504664669781, 2.5765725275664879, 7.6987743736515295,
                      5.0423700574830965, 17.435300118916864};
  const double *tmp=f2->getArray()->getConstPointer();
  std::transform(tmp,tmp+18,values2,values2,std::minus<double>());
  std::transform(values2,values2+18,values2,std::ptr_fun<double,double>(fabs));
  double max=*std::max_element(values2,values2+18);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f2->decrRef();
  //
  f1->applyFunc(1,"x+y");
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  double values1[9]={-1.8,-0.3,1.2,-0.3,1.2,2.7,1.2,2.7,4.2};
  tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+9,values1,values1,std::minus<double>());
  std::transform(values1,values1+9,values1,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values1,values1+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest::testOperationsOnFields()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=m->fillFromAnalytic(ON_NODES,1,func1);
  MEDCouplingFieldDouble *f2=m->fillFromAnalytic(ON_NODES,1,func1);
  f1->checkCoherency();
  f2->checkCoherency();
  MEDCouplingFieldDouble *f3=(*f1)+(*f2);
  f3->checkCoherency();
  CPPUNIT_ASSERT(f3->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f3->getTimeDiscretization()==NO_TIME);
  double values1[9]={-1.2,-0.2,0.8,-0.2,0.8,1.8,0.8,1.8,2.8};
  const double *tmp=f3->getArray()->getConstPointer();
  std::transform(tmp,tmp+9,values1,values1,std::minus<double>());
  std::transform(values1,values1+9,values1,std::ptr_fun<double,double>(fabs));
  double max=*std::max_element(values1,values1+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f3->decrRef();
  //
  f3=(*f1)*(*f2);
  f3->checkCoherency();
  CPPUNIT_ASSERT(f3->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f3->getTimeDiscretization()==NO_TIME);
  double values2[9]={0.36,0.01,0.16,0.01,0.16,0.81,0.16,0.81,1.96};
  tmp=f3->getArray()->getConstPointer();
  std::transform(tmp,tmp+9,values2,values2,std::minus<double>());
  std::transform(values2,values2+9,values2,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values2,values2+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f3->decrRef();
  //
  f3=(*f1)+(*f2);
  MEDCouplingFieldDouble *f4=(*f1)-(*f3);
  f4->checkCoherency();
  CPPUNIT_ASSERT(f4->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f4->getTimeDiscretization()==NO_TIME);
  double values3[9]={0.6,0.1,-0.4,0.1,-0.4,-0.9,-0.4,-0.9,-1.4};
  tmp=f4->getArray()->getConstPointer();
  std::transform(tmp,tmp+9,values3,values3,std::minus<double>());
  std::transform(values3,values3+9,values3,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values3,values3+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f3->decrRef();
  f4->decrRef();
  //
  f3=(*f1)+(*f2);
  f4=(*f3)/(*f2);
  f4->checkCoherency();
  CPPUNIT_ASSERT(f4->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f4->getTimeDiscretization()==NO_TIME);
  tmp=f4->getArray()->getConstPointer();
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,tmp[i],1.e-12);
  f3->decrRef();
  f4->decrRef();
  //
  f4=f2->buildNewTimeReprFromThis(ONE_TIME,false);
  f4->checkCoherency();
  CPPUNIT_ASSERT(f4->getArray()==f2->getArray());
  CPPUNIT_ASSERT(f4->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f4->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_THROW(f3=(*f1)+(*f4),INTERP_KERNEL::Exception);
  MEDCouplingFieldDouble *f5=f4->buildNewTimeReprFromThis(NO_TIME,false);
  CPPUNIT_ASSERT(f4->getArray()==f5->getArray());
  CPPUNIT_ASSERT(f5->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f5->getTimeDiscretization()==NO_TIME);
  f3=(*f1)+(*f5);
  tmp=f3->getArray()->getConstPointer();
  double values4[9]={-1.2,-0.2,0.8,-0.2,0.8,1.8,0.8,1.8,2.8};
  std::transform(tmp,tmp+9,values4,values4,std::minus<double>());
  std::transform(values4,values4+9,values4,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values4,values4+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f5->decrRef();
  f4->decrRef();
  f3->decrRef();
  //
  f4=f2->buildNewTimeReprFromThis(ONE_TIME,true);
  f4->checkCoherency();
  CPPUNIT_ASSERT(f4->getArray()!=f2->getArray());
  CPPUNIT_ASSERT(f4->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f4->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_THROW(f3=(*f1)+(*f4),INTERP_KERNEL::Exception);
  f5=f4->buildNewTimeReprFromThis(NO_TIME,true);
  CPPUNIT_ASSERT(f4->getArray()!=f5->getArray());
  CPPUNIT_ASSERT(f2->getArray()!=f5->getArray());
  CPPUNIT_ASSERT(f5->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f5->getTimeDiscretization()==NO_TIME);
  f3=(*f1)+(*f5);
  tmp=f3->getArray()->getConstPointer();
  double values5[9]={-1.2,-0.2,0.8,-0.2,0.8,1.8,0.8,1.8,2.8};
  std::transform(tmp,tmp+9,values5,values5,std::minus<double>());
  std::transform(values5,values5+9,values5,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values5,values5+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f5->decrRef();
  f4->decrRef();
  f3->decrRef();
  //
  f1->decrRef();
  f2->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest::testOperationsOnFields2()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  MEDCouplingFieldDouble *f1=m->fillFromAnalytic(ON_NODES,1,"x+y+z");
  MEDCouplingFieldDouble *f2=m->fillFromAnalytic(ON_NODES,1,"a*a+b+c*c");
  MEDCouplingFieldDouble *f3=(*f1)/(*f2);
  f3->checkCoherency();
  CPPUNIT_ASSERT(f3->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f3->getTimeDiscretization()==NO_TIME);
  const double expected1[9]={-2.4999999999999991, 1.2162162162162162, 0.77868852459016391,
                             0.7407407407407407, 1.129032258064516, 0.81632653061224492,
                             0.86538461538461531, 1.0919540229885056, 0.84302325581395343};
  CPPUNIT_ASSERT_EQUAL(1,f3->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f3->getNumberOfTuples());
  const double *val=f3->getArray()->getConstPointer();
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],val[i],1.e-12);
  f3->decrRef();
  f1->decrRef();
  f2->decrRef();
  m->decrRef();
}

bool func4(const double *pt, double *res)
{
  res[0]=pt[0]+pt[1]+pt[2];
  return true;
}

void MEDCouplingBasicsTest::testMergeNodesOnField()
{
  double *tmp;
  MEDCouplingUMesh *targetMesh=build3DTargetMeshMergeNode_1();
  MEDCouplingFieldDouble *f1=targetMesh->fillFromAnalytic(ON_NODES,1,func4);
  f1->mergeNodes(1e-10);
  f1->decrRef();
  targetMesh->decrRef();
  //
  targetMesh=build3DTargetMeshMergeNode_1();
  f1=targetMesh->fillFromAnalytic(ON_NODES,1,func4);
  tmp=f1->getArray()->getPointer();
  tmp[0]=1000.;
  f1->mergeNodes(1e-10);
  f1->decrRef();
  targetMesh->decrRef();
  //
  targetMesh=build3DTargetMeshMergeNode_1();
  f1=targetMesh->fillFromAnalytic(ON_NODES,1,func4);
  tmp=f1->getArray()->getPointer();
  tmp[1]=1000.;
  CPPUNIT_ASSERT_THROW(f1->mergeNodes(1e-10),INTERP_KERNEL::Exception);
  f1->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::testCheckConsecutiveCellTypes()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  CPPUNIT_ASSERT(sourceMesh->checkConsecutiveCellTypes());
  CPPUNIT_ASSERT(!targetMesh->checkConsecutiveCellTypes());
  targetMesh->decrRef();
  sourceMesh->decrRef();
}

void MEDCouplingBasicsTest::testBuildOrthogonalField()
{
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  MEDCouplingFieldDouble *field=targetMesh->buildOrthogonalField();
  double expected[3]={0.70710678118654746,0.,-0.70710678118654746};
  CPPUNIT_ASSERT_EQUAL(5,field->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,field->getNumberOfComponents());
  const double *vals=field->getArray()->getConstPointer();
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[i%3],vals[i],1e-12);
  field->decrRef();
  targetMesh->decrRef();
  // testing 
  double targetCoords[12]={0.,0.,0.,0.5,0.,0.5,1.,0.,1.,0.,1.,0.};
  int targetConn[4]={0,1,2,3};
  targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(1);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(4,3);
  std::copy(targetCoords,targetCoords+12,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  field=targetMesh->buildOrthogonalField();
  CPPUNIT_ASSERT_EQUAL(1,field->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,field->getNumberOfComponents());
  vals=field->getArray()->getConstPointer();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.70710678118654746,vals[0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,vals[1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.70710678118654746,vals[2],1e-12);
  field->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::testGetCellsContainingPoint()
{
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  double pos[12]={0.,0.,0.4,0.4,0.,0.4,0.1,0.1,0.25,0.,0.65,0.};
  std::vector<int> t1,t2;
  //2D basic
  targetMesh->getCellsContainingPoints(pos,6,1e-12,t1,t2);
  CPPUNIT_ASSERT_EQUAL(6,(int)t1.size());
  CPPUNIT_ASSERT_EQUAL(7,(int)t2.size());
  const int expectedValues1[6]={0,4,3,0,1,2};
  const int expectedValues2[7]={0,1,2,3,4,5,6};
  CPPUNIT_ASSERT(std::equal(t1.begin(),t1.end(),expectedValues1));
  CPPUNIT_ASSERT(std::equal(t2.begin(),t2.end(),expectedValues2));
  //2D with no help of bounding box.
  double center[2]={0.2,0.2};
  MEDCouplingPointSet::rotate2DAlg(center,0.78539816339744830962,6,pos);
  targetMesh->rotate(center,0,0.78539816339744830962);
  t1.clear(); t2.clear();
  targetMesh->getCellsContainingPoints(pos,6,1e-12,t1,t2);
  CPPUNIT_ASSERT_EQUAL(6,(int)t1.size());
  CPPUNIT_ASSERT_EQUAL(7,(int)t2.size());
  CPPUNIT_ASSERT(std::equal(t1.begin(),t1.end(),expectedValues1));
  CPPUNIT_ASSERT(std::equal(t2.begin(),t2.end(),expectedValues2));
  //2D outside
  const double pos1bis[2]={-0.3303300858899107,-0.11819805153394641};
  CPPUNIT_ASSERT_EQUAL(-1,targetMesh->getCellContainingPoint(pos1bis,1e-12));
  targetMesh->decrRef();
  //test limits 2D
  targetMesh=build2DTargetMesh_1();
  const double pos2[2]={0.2,-0.05};
  t1.clear();
  targetMesh->getCellsContainingPoint(pos2,1e-12,t1);
  CPPUNIT_ASSERT_EQUAL(2,(int)t1.size());
  const int expectedValues3[2]={0,1};
  CPPUNIT_ASSERT(std::equal(t1.begin(),t1.end(),expectedValues3));
  const double pos3[2]={0.2,0.2};
  t1.clear();
  targetMesh->getCellsContainingPoint(pos3,1e-12,t1);
  CPPUNIT_ASSERT_EQUAL(5,(int)t1.size());
  const int expectedValues4[5]={0,1,2,3,4};
  CPPUNIT_ASSERT(std::equal(t1.begin(),t1.end(),expectedValues4));
  CPPUNIT_ASSERT_EQUAL(0,targetMesh->getCellContainingPoint(pos3,1e-12));
  targetMesh->decrRef();
  //3D
  targetMesh=build3DTargetMesh_1();
  const double pos4[3]={25.,25.,25.};
  CPPUNIT_ASSERT_EQUAL(0,targetMesh->getCellContainingPoint(pos4,1e-12));
  const double pos5[3]={50.,50.,50.};
  t1.clear();
  targetMesh->getCellsContainingPoint(pos5,1e-12,t1);
  CPPUNIT_ASSERT_EQUAL(8,(int)t1.size());
  const int expectedValues5[8]={0,1,2,3,4,5,6,7};
  CPPUNIT_ASSERT(std::equal(t1.begin(),t1.end(),expectedValues5));
  const double pos6[3]={0., 50., 0.};
  t1.clear();
  targetMesh->getCellsContainingPoint(pos6,1e-12,t1);
  CPPUNIT_ASSERT_EQUAL(2,(int)t1.size());
  const int expectedValues6[2]={0,2};
  CPPUNIT_ASSERT(std::equal(t1.begin(),t1.end(),expectedValues6));
  //3D outside
  const double pos7[3]={-1.0,-1.0,0.};
  CPPUNIT_ASSERT_EQUAL(-1,targetMesh->getCellContainingPoint(pos7,1e-12));
  //3D outside 2
  const double center2[3]={0.,0.,0.};
  const double vec2[3]={0.,-1.,0.};
  targetMesh->rotate(center2,vec2,0.78539816339744830962);
  const double pos8[3]={-25,25.,12.};
  CPPUNIT_ASSERT_EQUAL(-1,targetMesh->getCellContainingPoint(pos8,1e-12));
  //
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::testGetValueOn1()
{
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  MEDCouplingFieldDouble *fieldOnCells=MEDCouplingFieldDouble::New(ON_CELLS);
  int nbOfCells=targetMesh->getNumberOfCells();
  fieldOnCells->setMesh(targetMesh);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(nbOfCells,2);
  fieldOnCells->setArray(array);
  double *tmp=array->getPointer();
  for(int i=0;i<nbOfCells;i++)
    { tmp[2*i]=7.+(double)i; tmp[2*i+1]=17.+(double)i; }
  array->decrRef();
  //
  const double pos1[2]={0.25,0.};
  double res[2];
  fieldOnCells->getValueOn(pos1,res);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.,res[0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(18.,res[1],1e-12);
  //
  fieldOnCells->decrRef();
  targetMesh->decrRef();
  //
  targetMesh=build2DSourceMesh_1();
  MEDCouplingFieldDouble *fieldOnNodes=MEDCouplingFieldDouble::New(ON_NODES);
  int nbOfNodes=targetMesh->getNumberOfNodes();
  fieldOnNodes->setMesh(targetMesh);
  array=DataArrayDouble::New();
  array->alloc(nbOfNodes,2);
  fieldOnNodes->setArray(array);
  tmp=array->getPointer();
  for(int i=0;i<nbOfNodes;i++)
    { tmp[2*i]=17.+(double)i; tmp[2*i+1]=27.+(double)i; }
  array->decrRef();
  //
  const double pos2[2]={-0.13333333333333333,-0.13333333333333333};
  fieldOnNodes->getValueOn(pos2,res);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(17.5,res[0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(27.5,res[1],1e-12);
  const double pos3[2]={0.033333333333333326,0.36666666666666664};
  fieldOnNodes->getValueOn(pos3,res);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(18.666666666666667,res[0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(28.666666666666667,res[1],1e-12);
  //
  fieldOnNodes->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::testCMesh0()
{
  MEDCouplingCMesh* mesh=MEDCouplingCMesh::New();
  DataArrayDouble* coordsX=DataArrayDouble::New();
  double arrX[4] = { -1., 1., 2., 4. };
  coordsX->useArray(arrX,false, CPP_DEALLOC,4,1);
  DataArrayDouble* coordsY=DataArrayDouble::New();
  double arrY[4] = { -2., 2., 4., 8. };
  coordsY->useArray(arrY,false, CPP_DEALLOC,4,1);
  DataArrayDouble* coordsZ=DataArrayDouble::New();
  double arrZ[4] = { -3., 3., 6., 12. };
  coordsZ->useArray(arrZ,false, CPP_DEALLOC,4,1);
  mesh->setCoords(coordsX,coordsY,coordsZ);
  coordsX->decrRef();
  coordsY->decrRef();
  coordsZ->decrRef();
  //
  MEDCouplingFieldDouble *fieldOnNodes=mesh->fillFromAnalytic(ON_NODES,1,"x+y/2.+z/3.");
  CPPUNIT_ASSERT_EQUAL(1,fieldOnNodes->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(64,fieldOnNodes->getNumberOfTuples());
  const double expected1[64]={-3., -1., 0., 2., -1., 1., 2., 4., 0., 2., 3., 5., 2., 4., 5., 7., -1., 1., 2.,
                              4., 1., 3., 4., 6., 2., 4., 5., 7., 4., 6., 7., 9., 0., 2., 3., 5., 2., 4., 5.,
                              7., 3., 5., 6., 8., 5., 7., 8., 10., 2., 4., 5.,
                              7., 4., 6., 7., 9., 5., 7., 8., 10., 7., 9., 10., 12.};
  const double *val=fieldOnNodes->getArray()->getConstPointer();
  for(int i=0;i<64;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],val[i],1e-12);
  double res;
  fieldOnNodes->getValueOnPos(1,3,2,&res);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,res,1e-12);
  fieldOnNodes->decrRef();
  //
  MEDCouplingFieldDouble *fieldOnCells=mesh->fillFromAnalytic(ON_CELLS,1,"x+y/2.+z/3.");
  CPPUNIT_ASSERT_EQUAL(1,fieldOnCells->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(27,fieldOnCells->getNumberOfTuples());
  val=fieldOnCells->getArray()->getConstPointer();
  const double expected2[27]={0, 1.5, 3, 1.5, 3, 4.5, 3, 4.5, 6, 1.5, 3, 4.5, 3, 4.5,
                              6, 4.5, 6, 7.5, 3, 4.5, 6, 4.5, 6, 7.5, 6, 7.5, 9};
  for(int i=0;i<27;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],val[i],1e-12);
  fieldOnCells->getValueOnPos(1,2,1,&res);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(6.,res,1e-12);
  fieldOnCells->decrRef();
  //
  mesh->decrRef();
}

void MEDCouplingBasicsTest::testScale()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  const double pos[2]={0.2,0.2};
  mesh->scale(pos,0.5);
  const double expected1[18]={-0.05,-0.05, 0.2,-0.05, 0.45,-0.05, -0.05,0.2, 0.2,0.2, 0.45,0.2,
                              -0.05,0.45, 0.2,0.45, 0.45,0.45};
  const double *val=mesh->getCoords()->getConstPointer();
  for(int i=0;i<18;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],val[i],1e-12);
  //
  mesh->decrRef();
}
