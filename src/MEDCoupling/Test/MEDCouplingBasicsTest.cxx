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
#include "MEDCouplingMemArray.hxx"
#include "Interpolation2D.txx"
#include "Interpolation3DSurf.txx"
#include "Interpolation3D.txx"
#include "InterpolationCC.txx"
#include <InterpolationCU2D.txx>

#include "MEDCouplingNormalizedUnstructuredMesh.txx"
#include "MEDCouplingNormalizedCartesianMesh.txx"

#include <cmath>

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

void MEDCouplingBasicsTest::test2DInterpP0P0PL_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  vector<map<int,double> > res;
  //
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  //
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test2DInterpP0P0PL_2()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  std::vector<int> cellsIds(targetMesh->getNumberOfCells());
  for(int i=0;i<targetMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  targetMesh->convertToPolyTypes(cellsIds);
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  vector<map<int,double> > res;
  //
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  //
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test2DInterpP0P0PL_3()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  std::vector<int> cellsIds(sourceMesh->getNumberOfCells());
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  sourceMesh->convertToPolyTypes(cellsIds);
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  vector<map<int,double> > res;
  //
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  //
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test2DInterpP0P0PL_4()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  std::vector<int> cellsIds(sourceMesh->getNumberOfCells());
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  sourceMesh->convertToPolyTypes(cellsIds);
  cellsIds.resize(targetMesh->getNumberOfCells());
  for(int i=0;i<targetMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  targetMesh->convertToPolyTypes(cellsIds);
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  vector<map<int,double> > res;
  //
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  //
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

void MEDCouplingBasicsTest::test2DInterpP0P1PL_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");
  CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[8][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[8][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.,sumAll(res),1e-12);
  res.clear();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test2DInterpP0P1PL_2()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  std::vector<int> cellsIds(sourceMesh->getNumberOfCells());
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  sourceMesh->convertToPolyTypes(cellsIds);
  //
  cellsIds.resize(targetMesh->getNumberOfCells());
  for(int i=0;i<targetMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  targetMesh->convertToPolyTypes(cellsIds);
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");
  CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[8][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[8][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.,sumAll(res),1e-12);
  res.clear();
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

void MEDCouplingBasicsTest::test2DInterpP1P0PL_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[0][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.333333333333333333,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[1][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666666,res[1][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666666,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[2][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.333333333333333333,res[2][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[3][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[4][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5,res[4][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  res.clear();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test2DInterpP1P0PL_2()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  std::vector<int >cellsIds(targetMesh->getNumberOfCells());
  for(int i=0;i<targetMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  targetMesh->convertToPolyTypes(cellsIds);
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[0][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.333333333333333333,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[1][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666666,res[1][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666666,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[2][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.333333333333333333,res[2][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[3][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[4][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5,res[4][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  res.clear();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test2DInterpP1P1_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_2();
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
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P1");
      CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08333333333333334,res[0][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.05416666666666665,res[1][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666666,res[1][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08333333333333334,res[2][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.05416666666666665,res[3][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666668,res[3][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1416666666666666,res[4][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02499999999999999,res[4][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02499999999999999,res[4][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.09999999999999999,res[4][3],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666666,res[5][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.09583333333333333,res[5][3],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08333333333333333,res[6][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666667,res[7][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.09583333333333331,res[7][3],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.04166666666666668,res[8][3],1.e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test2DInterpP1P1PL_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P1");
  CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,res[0][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,res[2][1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][2],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.,res[4][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.,res[4][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][2],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][2],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,res[8][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(25.,sumAll(res),1e-12);
  res.clear();
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

void MEDCouplingBasicsTest::test3DSurfInterpP0P0PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  res.clear();
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

void MEDCouplingBasicsTest::test3DSurfInterpP0P1PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");
  CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[8][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[8][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.,sumAll(res),1e-12);
  res.clear();
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

void MEDCouplingBasicsTest::test3DSurfInterpP1P0PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[0][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.333333333333333333,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[1][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666666,res[1][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.166666666666666666,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[2][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.333333333333333333,res[2][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[3][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res[4][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5,res[4][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,sumAll(res),1e-12);
  res.clear();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DSurfInterpP1P1_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_2();
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
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P1");
      CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08333333333333334*sqrt(2.),res[0][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.05416666666666665*sqrt(2.),res[1][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666666*sqrt(2.),res[1][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08333333333333334*sqrt(2.),res[2][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.05416666666666665*sqrt(2.),res[3][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666668*sqrt(2.),res[3][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1416666666666666*sqrt(2.),res[4][0],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02499999999999999*sqrt(2.),res[4][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02499999999999999*sqrt(2.),res[4][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.09999999999999999*sqrt(2.),res[4][3],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666666*sqrt(2.),res[5][1],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.09583333333333333*sqrt(2.),res[5][3],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08333333333333333*sqrt(2.),res[6][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.02916666666666667*sqrt(2.),res[7][2],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.09583333333333331*sqrt(2.),res[7][3],1.e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.04166666666666668*sqrt(2.),res[8][3],1.e-12);
      res.clear();
    }
  //
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DSurfInterpP1P1PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P1");
  CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,res[0][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,res[2][1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][2],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.,res[4][0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.,res[4][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][2],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][2],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,res[8][3],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(25.,sumAll(res),1e-12);
  res.clear();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DSurfInterpP0P0_2()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMeshPerm_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::Triangulation);
  {
    myInterpolator.setOrientation(2);
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
  {
    myInterpolator.setOrientation(0);
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
    CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[1][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.125*sqrt(2.),res[2][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[3][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[4][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[4][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.75*sqrt(2.),sumAll(res),1e-12);
    res.clear();
  }
  {
    myInterpolator.setOrientation(1);
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
    CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[1][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[3][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[4][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[4][1],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.875*sqrt(2.),sumAll(res),1e-12);
    res.clear();
  }
  {
    myInterpolator.setOrientation(-1);
    myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
    CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[2][0],1e-12);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),sumAll(res),1e-12);
    res.clear();
  }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

/*!
 * Test of precision option implemented by Fabien that represents distance of "barycenter" to the other cell.
 */
void MEDCouplingBasicsTest::test3DSurfInterpP0P0_3()
{
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  vector<map<int,double> > res;
  double vecTrans[3]={0.,0.,1.e-10};
  double vec[3]={0.,-1.,0.};
  double pt[3]={-0.3,-0.3,5.e-11};
  const int N=32;
  const double deltaA=M_PI/N;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::Triangulation);
  myInterpolator.setMaxDistance3DSurfIntersect(1e-9);
  for(int i=0;i<N;i++)
    {
      res.clear();
      MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_2();
      sourceMesh->rotate(pt,vec,i*deltaA);
      MEDCouplingUMesh *targetMesh=build3DSurfSourceMesh_2();
      targetMesh->translate(vecTrans);
      targetMesh->rotate(pt,vec,i*deltaA);
      MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
      MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
      CPPUNIT_ASSERT_EQUAL(2,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
      sourceMesh->decrRef();
      targetMesh->decrRef();
    }
  //
  myInterpolator.setMaxDistance3DSurfIntersect(1e-11);
  for(int i=0;i<N;i++)
    {
      res.clear();
      MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_2();
      sourceMesh->rotate(pt,vec,i*deltaA);
      MEDCouplingUMesh *targetMesh=build3DSurfSourceMesh_2();
      targetMesh->translate(vecTrans);
      targetMesh->rotate(pt,vec,i*deltaA);
      MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
      MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
      CPPUNIT_ASSERT_EQUAL(2,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,sumAll(res),1e-12);
      sourceMesh->decrRef();
      targetMesh->decrRef();
    }
  //
  res.clear();
  myInterpolator.setMaxDistance3DSurfIntersect(-1.);//unactivate fabien lookup
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_2();
  MEDCouplingUMesh *targetMesh=build3DSurfSourceMesh_2();
  targetMesh->translate(vecTrans);
  myInterpolator.setBoundingBoxAdjustment(1e-11);
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper0(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper0(targetMesh);
  myInterpolator.interpolateMeshes(sourceWrapper0,targetWrapper0,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(2,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,sumAll(res),1e-12);
  sourceMesh->decrRef();
  targetMesh->decrRef();
  //
  res.clear();
  sourceMesh=build3DSurfSourceMesh_2();
  targetMesh=build3DSurfSourceMesh_2();
  targetMesh->translate(vecTrans);
  myInterpolator.setBoundingBoxAdjustment(1e-9);
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper1(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper1(targetMesh);
  myInterpolator.interpolateMeshes(sourceWrapper1,targetWrapper1,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(2,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
  sourceMesh->decrRef();
  targetMesh->decrRef();
  //keeping the same bbox adj == 1.e-11 but trying rotation
  res.clear();
  sourceMesh=build3DSurfSourceMesh_2();
  sourceMesh->rotate(pt,vec,M_PI/4.);
  targetMesh=build3DSurfSourceMesh_2();
  targetMesh->translate(vecTrans);
  targetMesh->rotate(pt,vec,M_PI/4.);
  myInterpolator.setBoundingBoxAdjustment(1e-11);
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper2(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper2(targetMesh);
  myInterpolator.interpolateMeshes(sourceWrapper2,targetWrapper2,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(2,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DInterpP0P0_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.e6,sumAll(res),1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[0][0],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(41666.66666666667,res[0][6],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[0][7],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[0][8],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[0][10],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(41666.66666666667,res[1][2],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[1][7],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[1][8],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[2][0],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[2][5],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333333,res[2][6],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[2][9],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[2][11],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(395833.3333333333,res[3][0],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333333,res[3][2],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333331,res[3][3],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[3][5],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(395833.3333333333,res[3][8],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[4][1],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[4][4],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333333,res[4][6],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333333,res[4][9],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[4][10],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333333,res[5][2],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(20833.33333333331,res[5][3],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[5][4],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(395833.3333333333,res[5][7],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(395833.3333333333,res[5][10],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[6][1],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(250000,res[6][6],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(541666.6666666667,res[6][9],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[6][11],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(83333.33333333331,res[7][0],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(479166.6666666667,res[7][1],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(333333.3333333333,res[7][2],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(624999.9999999997,res[7][3],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(479166.6666666667,res[7][4],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(479166.6666666667,res[7][5],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(83333.33333333333,res[7][6],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(83333.33333333331,res[7][7],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(83333.33333333333,res[7][8],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(83333.33333333333,res[7][9],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(83333.33333333331,res[7][10],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(479166.6666666667,res[7][11],1e-7);
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DInterpP0P0PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][10],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][10],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][10],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][9],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][5],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][11],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(21.,sumAll(res),1e-12);
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DInterpP0P0PL_2()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  std::vector<int> cellsIds(targetMesh->getNumberOfCells());
  for(int i=0;i<targetMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  targetMesh->convertToPolyTypes(cellsIds);
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][10],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][10],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][10],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][9],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][5],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][11],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(21.,sumAll(res),1e-12);
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DInterpP0P0PL_3()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  std::vector<int> cellsIds(sourceMesh->getNumberOfCells());
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  sourceMesh->convertToPolyTypes(cellsIds);
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][10],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][10],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][10],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][9],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][5],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][11],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(21.,sumAll(res),1e-12);
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DInterpP0P0PL_4()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  std::vector<int> cellsIds(sourceMesh->getNumberOfCells());
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    cellsIds[i]=i;
  sourceMesh->convertToPolyTypes(cellsIds);
  cellsIds.resize(targetMesh->getNumberOfCells());
  for(int j=0;j<targetMesh->getNumberOfCells();j++)
    cellsIds[j]=j;
  targetMesh->convertToPolyTypes(cellsIds);
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][10],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][10],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][10],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][9],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][5],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][11],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(21.,sumAll(res),1e-12);
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DInterpP0P1_1()
{
  MEDCouplingUMesh *sourceMesh=build3DTargetMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSourceMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");
  CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(244444.4444444445,res[0][4],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333333,res[0][5],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(291666.6666666666,res[0][6],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[0][7],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(125000,res[1][0],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(140277.7777777778,res[1][1],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(119444.4444444444,res[1][2],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[1][3],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(119444.4444444444,res[1][4],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[1][5],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(26388.88888888889,res[1][6],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(348611.1111111111,res[2][6],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888888,res[2][7],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(244444.4444444444,res[3][2],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333334,res[3][3],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(291666.6666666666,res[3][6],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[3][7],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(536111.111111111,res[4][5],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(297222.2222222221,res[4][7],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(223611.1111111111,res[5][1],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(125000,res[5][3],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(125000,res[5][5],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(26388.88888888892,res[5][7],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[6][7],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(536111.1111111109,res[7][3],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(297222.2222222221,res[7][7],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(11111.1111111111,res[8][1],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(11111.11111111111,res[8][2],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666666,res[8][3],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(11111.11111111111,res[8][4],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[8][5],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[8][6],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1466666.666666668,res[8][7],1e-7);
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DInterpP0P1PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DTargetMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSourceMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P1");
  CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][5],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[8][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(9.,sumAll(res),1e-12);
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DInterpP1P0_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
  CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(125000,res[0][1],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(140277.7777777778,res[1][1],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(223611.1111111111,res[1][5],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(11111.1111111111,res[1][8],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(119444.4444444444,res[2][1],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(244444.4444444445,res[2][3],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(11111.11111111111,res[2][8],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[3][1],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333333,res[3][3],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(125000,res[3][5],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(536111.1111111109,res[3][7],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666667,res[3][8],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(244444.4444444445,res[4][0],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(119444.4444444445,res[4][1],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(11111.11111111111,res[4][8],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(145833.3333333333,res[5][0],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[5][1],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(536111.1111111109,res[5][4],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(125000,res[5][5],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666666,res[5][8],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(291666.6666666666,res[6][0],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(26388.88888888889,res[6][1],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(348611.1111111112,res[6][2],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(291666.6666666667,res[6][3],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(166666.6666666666,res[6][8],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[7][0],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[7][2],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(151388.8888888889,res[7][3],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(297222.2222222221,res[7][4],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(26388.88888888892,res[7][5],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[7][6],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(297222.2222222222,res[7][7],1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1466666.666666668,res[7][8],1e-7);
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DInterpP1P0PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
  CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.75,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.25,res[0][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[1][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][5],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[1][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[2][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[2][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[3][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[3][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[4][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[4][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[5][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,res[5][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[6][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[6][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[6][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[6][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.25,res[7][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.75,res[7][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(21.,sumAll(res),1e-12);
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DInterpP1P1_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_2();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_2();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P1");
  CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
  double res3D[8][28]= {{124999.999883775978, 245370.370390364464, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 203703.703634892299, 187500.000094145857, 0.0, 0.0, 4629.6296266718, 0.0, 215277.777751402784, 209722.222322299582, 0.0, 0.0, 0.0, 0.0, 104166.666590829205, 121296.296368812196, 0.0, 250000.000003472145},
                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 120370.370368827047, 0.0, 0.0, 38888.888897777797, 0.0, 0.0, 45370.3703701697596, 0.0, 0.0, 45370.3703701697596, 83333.3333263888926, 0.0},
                        {0.0, 0.0, 0.0, 97222.2222222221753, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 97222.2222222221608, 0.0, 97222.2222222222044, 41666.6666666666642, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {0.0, 277777.777787084982, 199074.074074073927, 0.0, 0.0, 0.0, 4629.62962962962774, 0.0, 321759.259254934732, 83333.3333333333139, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4629.62962667180363, 0.0, 0.0, 251388.88888319055, 194444.444454861077, 0.0, 79629.6296194135939, 250000.000003472145, 0.0, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 85185.1851851851534, 4629.62962962962774, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 118518.518518518511, 0.0, 41666.6666666666642, 83333.3333333333285, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {0.0, 324074.07407629228, 0.0, 0.0, 0.0, 247685.185185184964, 6481.48148148147993, 0.0, 173611.11111196311, 0.0, 164814.814814814832, 0.0, 4629.62962962962865, 208333.33333418527, 0.0, 83333.3333333333285, 203703.703697273799, 249999.999999999767, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {125000.000000000015, 423611.111111110775, 134259.259259259241, 194444.444444444351, 164814.814814814745, 164351.851851851825, 203703.703703703592, 249999.999999999825, 0.0, 0.0, 0.0, 0.0, 6481.48148148147902, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 118518.518518518453, 0.0, 4629.62962962962956, 83333.3333333333139, 85185.1851851851825, 41666.6666666666642, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  int i=0;
  double sum = 0;
  //cout.precision(18);
  for(std::vector<std::map<int,double> >::const_iterator iter1=res.begin();iter1!=res.end();iter1++,i++)
    {
      //cout<< "res3D[" <<i<< "][]={";
      for(int j=0;j<28;j++)
        {
          std::map<int,double>::const_iterator iter2=(*iter1).find(j);
          if(iter2!=(*iter1).end())
            {
              //cout<< iter2->second<< ", ";
              sum += iter2->second;
              CPPUNIT_ASSERT_DOUBLES_EQUAL(res3D[i][j],(*iter2).second,1.e-5);
            }
          else
            {
              //cout << "0.0, ";
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res3D[i][j],1e-14);
            }
        }
      //cout << "}" << endl;
    }
  //cout << "Sum = " << sum << endl;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8000000,sum,1.e-5);
  //clean-up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DInterpP1P1PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_2();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_2();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P1");
  CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(20.,res[0][24],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,res[1][26],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][21],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(24.,res[3][23],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][14],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(24.,res[5][17],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(24.,res[6][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][11],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(97.,sumAll(res),1e-12);
  //clean-up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DInterpP0P0Empty()
{
  MEDCouplingUMesh *sourceMesh=MEDCouplingUMesh::New();
  sourceMesh->setMeshDimension(2);
  sourceMesh->allocateCells(0);
  sourceMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(0,0);
  sourceMesh->setCoords(myCoords);
  myCoords->decrRef();
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(0);
  targetMesh->finishInsertingCells();
  myCoords=DataArrayDouble::New();
  myCoords->alloc(0,2);
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::testInterpolationCC()
{
  double arr1[3] = { 0/2., 1/2., 2/2. };
  double arr2[4] = { 0/3, 1/3., 2/3., 3/3. };
  MEDCouplingCMesh* mesh[2];
  for ( int i = 0; i < 2; ++i )
    {
      const double* arr = i ? arr1 : arr2;
      const int nb_coord = i ? 3 : 4;
      DataArrayDouble* coords = DataArrayDouble::New();
      coords->useArray( arr, /*ownership=*/false, CPP_DEALLOC, nb_coord, 1 );

      mesh[i] = MEDCouplingCMesh::New();
      mesh[i]->setCoords( coords, coords, coords );
      coords->decrRef();
    }
  MEDCouplingNormalizedCartesianMesh<3> targetWrapper(mesh[1]);
  MEDCouplingNormalizedCartesianMesh<3> sourceWrapper(mesh[0]);
  CPPUNIT_ASSERT_EQUAL( 27,int( sourceWrapper.getNumberOfElements()));
  CPPUNIT_ASSERT_EQUAL( 3, int( sourceWrapper.nbCellsAlongAxis(0)));
  CPPUNIT_ASSERT_EQUAL( 3, int( sourceWrapper.nbCellsAlongAxis(1)));
  CPPUNIT_ASSERT_EQUAL( 3, int( sourceWrapper.nbCellsAlongAxis(2)));
  CPPUNIT_ASSERT_THROW( sourceWrapper.nbCellsAlongAxis(3), INTERP_KERNEL::Exception);

  INTERP_KERNEL::InterpolationCC myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");

  CPPUNIT_ASSERT_EQUAL(8,int( res.size()));
  CPPUNIT_ASSERT_EQUAL(8,int( res[0].size()));
  const double precis = 1e-7;
  set<double> vals;
  double sum = 0;
  for ( int i = 0; i < (int)res.size(); ++i )
    for ( map<int,double>::iterator s_v = res[i].begin(); s_v != res[i].end(); ++s_v)
      {
        sum += s_v->second;
        double vvv;
#ifdef WNT
        double vv = s_v->second / precis;
        if(vv>=0.0)
          {
            vvv = floor(vv+0.5);
          }
        else
          {
            vvv = ceil(vv-0.5);
          }
#else
        vvv = round( s_v->second / precis );
#endif
        vals.insert( precis * vvv );
      }
  //cout << "tgt: " << i << " src: " << s_v->first << " - w: " << s_v->second << endl;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, sum, precis );

  set<double>::iterator v = vals.begin();
  CPPUNIT_ASSERT_EQUAL( 4, int( vals.size()) );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00462963, *v++, precis );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.00925926, *v++, precis );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.01851850, *v++, precis );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.03703700, *v++, precis );

  mesh[0]->decrRef();
  mesh[1]->decrRef();
}

void MEDCouplingBasicsTest::testInterpolationCU2D()
{
  MEDCouplingCMesh* meshC = MEDCouplingCMesh::New();
  DataArrayDouble* coords = DataArrayDouble::New();
  double arr[4] = { 0/3, 1/3., 2/3., 3/3. };
  coords->useArray( arr, /*ownership=*/false, CPP_DEALLOC, 4, 1 );
  meshC->setCoords( coords, coords );
  coords->decrRef();

  MEDCouplingUMesh * meshU = buildCU2DMesh_U();

  MEDCouplingNormalizedCartesianMesh<2>      sourceWrapper(meshC);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(meshU);
  INTERP_KERNEL::InterpolationCU2D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");

//   cout.precision(18);
//   for ( int i = 0; i < (int)res.size(); ++i )
//     for ( map<int,double>::iterator s_v = res[i].begin(); s_v != res[i].end(); ++s_v)
//     {
//       cout << "CPPUNIT_ASSERT_DOUBLES_EQUAL( "<<s_v->second<<" ,res["<<i<<"]["<<s_v->first<<"],precis);"<<endl;
//     }

  const double precis = 1e-7;
  CPPUNIT_ASSERT_EQUAL(5,int( res.size()));
  double sum = sumAll(res);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 1, sum, precis );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0555556 ,res[0][1],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0555556 ,res[0][3],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0277778 ,res[0][4],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0555556 ,res[1][3],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0277778 ,res[1][4],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1111111 ,res[1][6],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0555556 ,res[1][7],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0277778 ,res[2][4],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0555556 ,res[2][5],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0555556 ,res[2][7],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1111111 ,res[2][8],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0416667 ,res[3][1],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0138889 ,res[3][2],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0277778 ,res[3][4],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0416667 ,res[3][5],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0138889 ,res[4][1],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0972222 ,res[4][2],precis);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0138889 ,res[4][5],precis);


  meshC->decrRef();
  meshU->decrRef();
}

void MEDCouplingBasicsTest::test2DInterpP0IntegralUniform()
{
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  vector<map<int,double> > res;
  CPPUNIT_ASSERT_EQUAL(5,myInterpolator.fromIntegralUniform(targetWrapper,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
  res.clear();
  CPPUNIT_ASSERT_EQUAL(1,myInterpolator.toIntegralUniform(targetWrapper,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[3][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
  res.clear();
  targetMesh->decrRef();
  //
  targetMesh=build2DTargetMeshPerm_1();
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper2(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator2;
  CPPUNIT_ASSERT(myInterpolator2.getMeasureAbsStatus());
  CPPUNIT_ASSERT_EQUAL(5,myInterpolator2.fromIntegralUniform(targetWrapper2,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
  res.clear();
  myInterpolator2.setMeasureAbsStatus(false);
  CPPUNIT_ASSERT(!myInterpolator2.getMeasureAbsStatus());
  CPPUNIT_ASSERT_EQUAL(5,myInterpolator2.fromIntegralUniform(targetWrapper2,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.125,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,res[0][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.75,sumAll(res),1e-12);
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DSurfInterpP0IntegralUniform()
{
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  vector<map<int,double> > res;
  CPPUNIT_ASSERT_EQUAL(5,myInterpolator.fromIntegralUniform(targetWrapper,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[0][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[0][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.*sqrt(2.),sumAll(res),1e-12);
  res.clear();
  CPPUNIT_ASSERT_EQUAL(1,myInterpolator.toIntegralUniform(targetWrapper,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[3][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25*sqrt(2.),res[4][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.*sqrt(2.),sumAll(res),1e-12);
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DInterpP0IntegralUniform()
{
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  INTERP_KERNEL::Interpolation3D myInterpolator;
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  vector<map<int,double> > res;
  CPPUNIT_ASSERT_EQUAL(8,myInterpolator.fromIntegralUniform(targetWrapper,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(125000.,res[0][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(375000.,res[0][1],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(375000.,res[0][2],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1125000.,res[0][3],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(375000.,res[0][4],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1125000.,res[0][5],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1125000.,res[0][6],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3375000.,res[0][7],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8000000.,sumAll(res),1e-6);
  res.clear();
  CPPUNIT_ASSERT_EQUAL(1,myInterpolator.toIntegralUniform(targetWrapper,res,"P0"));
  CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(125000.,res[0][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(375000.,res[1][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(375000.,res[2][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1125000.,res[3][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(375000.,res[4][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1125000.,res[5][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1125000.,res[6][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3375000.,res[7][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8000000.,sumAll(res),1e-6);
  res.clear();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test2DInterpP1IntegralUniform()
{
  MEDCouplingUMesh *targetMesh=build2DSourceMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  vector<map<int,double> > res;
  CPPUNIT_ASSERT_EQUAL(4,myInterpolator.fromIntegralUniform(targetWrapper,res,"P1"));
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.33333333333333331,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[0][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[0][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.33333333333333331,res[0][3],1e-12);
  res.clear();
  CPPUNIT_ASSERT_EQUAL(1,myInterpolator.toIntegralUniform(targetWrapper,res,"P1"));
  CPPUNIT_ASSERT_EQUAL(4,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.33333333333333331,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[1][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[2][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.33333333333333331,res[3][0],1e-12);
  res.clear();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DInterpP1IntegralUniform()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(sourceMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  vector<map<int,double> > res;
  CPPUNIT_ASSERT_EQUAL(9,myInterpolator.fromIntegralUniform(targetWrapper,res,"P1"));
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[0][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[0][1],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(500000.,res[0][2],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[0][3],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[0][4],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(500000.,res[0][5],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[0][6],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[0][7],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2000000.,res[0][8],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8000000.,sumAll(res),1e-6);
  res.clear();
  CPPUNIT_ASSERT_EQUAL(1,myInterpolator.toIntegralUniform(targetWrapper,res,"P1"));
  CPPUNIT_ASSERT_EQUAL(9,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[0][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[1][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(500000.,res[2][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[3][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[4][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(500000.,res[5][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[6][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(833333.333333333,res[7][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2000000.,res[8][0],1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8000000.,sumAll(res),1e-6);
  sourceMesh->decrRef();
}

void MEDCouplingBasicsTest::test2DInterpP1P0Bary_1()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation2D myInterpolator;
  myInterpolator.setP1P0BaryMethod(true);
  vector<map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[2]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<2;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
      CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666669,res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343,res[0][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343,res[0][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[0][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625,res[1][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343,res[1][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343,res[2][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625,res[2][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[2][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625,res[3][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,res[3][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625,res[3][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664,res[4][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343,res[4][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343,res[4][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666,res[4][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DSurfInterpP1P0Bary_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSurfSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,2> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,2> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3DSurf myInterpolator;
  myInterpolator.setP1P0BaryMethod(true);
  vector<map<int,double> > res;
  INTERP_KERNEL::IntersectionType types[2]={INTERP_KERNEL::Triangulation, INTERP_KERNEL::Geometric2D};
  for(int i=0;i<2;i++)
    {
      myInterpolator.setPrecision(1e-12);
      myInterpolator.setIntersectionType(types[i]);
      myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
      CPPUNIT_ASSERT_EQUAL(5,(int)res.size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666669*sqrt(2.),res[0][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343*sqrt(2.),res[0][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343*sqrt(2.),res[0][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[0][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[1][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625*sqrt(2.),res[1][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343*sqrt(2.),res[1][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343*sqrt(2.),res[2][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625*sqrt(2.),res[2][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[2][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625*sqrt(2.),res[3][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125*sqrt(2.),res[3][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0625*sqrt(2.),res[3][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.041666666666666664*sqrt(2.),res[4][0],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343*sqrt(2.),res[4][1],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.020833333333333343*sqrt(2.),res[4][2],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.16666666666666666*sqrt(2.),res[4][3],1e-12);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.*sqrt(2.),sumAll(res),1e-12);
      res.clear();
    }
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

#include <iomanip>
void MEDCouplingBasicsTest::test3DInterpP1P0Bary_1()
{
  MEDCouplingUMesh *sourceMesh=build3DSourceMesh_2();
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_2();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  myInterpolator.setP1P0BaryMethod(true);
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P1P0");
  CPPUNIT_ASSERT_EQUAL(5,(int)res.size());

  double res3D[5][28]={{104166.66658918398, 885416.666685817763, 135416.666666666541, 36458.3333333335031, 31249.9999999999018, 145833.333333333256, 41666.6666666667516, 124999.999999999971, 177083.333326388849, 0.0, 31249.9999999999636, 0.0, 41666.666620792399, 159722.22229009436, 0.0, 0.0, 41666.6666631944681, 125000, 43499.2283723790752, 164351.851924000395, 36458.3333372396883, 0.0, 0.0, 125000.000001736029, 34722.2221800900952, 13599.5370788455439, 0.0, 167438.27159690368},
                       {0.0, 41666.6664479170649, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 125000.000161457952, 0.0, 0.0, 0.0, 0.0, 111111.11112005508, 0.0, 0.0, 291666.666656249959, 41666.6666666666933, 6944.4444415638809, 270833.333520485845, 0.0, 0.0, 124999.999989583303, 41666.6665798612958, 20833.3333186342825, 145833.333354303701, 83333.3333263888198, 27777.7777501651799},
                       {0.0, 93750.0000000000728, 125000.000000000058, 0.0, 0.0, 72916.666666666526, 291666.666666666628, 41666.6666666667152, 197916.66666666657, 166666.666666666802, 218750.000000000116, 41666.6666666665697, 0.0, 0.0, 0.0, 0.0, 0.0, 41666.6666666666861, 0.0, 0.0, 0.0, 0.0, 0.0, 41666.6666666666642, 0.0, 0.0, 0.0, 0.0},
                       {72916.6666484848247, 82465.2777799315081, 0.0, 0.0, 217447.916666666686, 197916.666666666802, 0.0, 41666.6666666666715, 0.0, 0.0, 0.0, 0.0, 290364.583310396119, 125000.000018181803, 41666.6666666666351, 166666.666666666599, 0.0, 41666.6666666665551, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 27777.7777734705051, 0.0, 0.0, 27777.7778028684952},
                       {72916.6666461071727, 172309.027782170655, 70312.5000000000437, 253906.250000000029, 0.0, 0.0, 0.0, 41666.666666666657, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 258246.527775988478, 71180.5555571812583, 253906.250006944378, 41666.6666666666861, 0.0, 41666.6666649305407, 20833.3333186342534, 6944.44445267237552, 0.0, 27777.7777953707919}};

  double sum = 0;
  int i=0;
  for(std::vector<std::map<int,double> >::const_iterator iter1=res.begin();iter1!=res.end();iter1++,i++)
    {
      for(int j=0;j<28;j++)
        {
          std::map<int,double>::const_iterator iter2=(*iter1).find(j);
          if(iter2!=(*iter1).end())
            {
              sum += iter2->second;
              CPPUNIT_ASSERT_DOUBLES_EQUAL(res3D[i][j],(*iter2).second,1.e-5);
            }
          else
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,res3D[i][j],1e-14);
            }
        }
    }
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8000000,sum,1.e-5);
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::test3DTo1DInterpP0P0PL_1()
{
  MEDCouplingUMesh *sourceMesh=build3DTargetMesh_1();
  MEDCouplingUMesh *targetMesh=build1DTargetMesh_1();
  //
  MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(sourceMesh);
  MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(targetMesh);
  INTERP_KERNEL::Interpolation3D myInterpolator;
  vector<map<int,double> > res;
  myInterpolator.setPrecision(1e-12);
  myInterpolator.setIntersectionType(INTERP_KERNEL::PointLocator);
  myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,res,"P0P0");
  CPPUNIT_ASSERT_EQUAL(8,(int)res.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[0][0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[1][4],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[2][1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[3][5],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[4][2],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[5][6],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[6][3],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,res[7][7],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.,sumAll(res),1e-12);
  //
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

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

double MEDCouplingBasicsTest::sumAll(const std::vector< std::map<int,double> >& matrix)
{
  double ret=0.;
  for(std::vector< std::map<int,double> >::const_iterator iter=matrix.begin();iter!=matrix.end();iter++)
    for(std::map<int,double>::const_iterator iter2=(*iter).begin();iter2!=(*iter).end();iter2++)
      ret+=(*iter2).second;
  return ret;
}
