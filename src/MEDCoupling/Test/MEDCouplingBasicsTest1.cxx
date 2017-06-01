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

#include "MEDCouplingBasicsTest1.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingMappedExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"

#include <sstream>
#include <cmath>
#include <algorithm>
#include <functional>

#ifdef WIN32
#include "MEDCouplingMemArray.txx"
#endif

using namespace MEDCoupling;

void MEDCouplingBasicsTest1::testArray()
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

void MEDCouplingBasicsTest1::testArray2()
{
  DataArrayDouble *arr=DataArrayDouble::New();
  arr->alloc(3,4);
  double *tmp=arr->getPointer();
  const double arrRef[12]={12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.};
  std::copy(arrRef,arrRef+12,tmp);
  arr->setInfoOnComponent(0,"ggg");
  arr->setInfoOnComponent(1,"hhhh");
  arr->setInfoOnComponent(2,"jj");
  arr->setInfoOnComponent(3,"kkkkkk");
  MCAuto<DataArrayInt> arr2(arr->convertToIntArr());
  MCAuto<DataArrayDouble> arr3(arr2->convertToDblArr());
  CPPUNIT_ASSERT(arr->isEqual(*arr3,1e-14));
  arr->decrRef();
}

void MEDCouplingBasicsTest1::testArray3()
{
  DataArrayInt *arr1=DataArrayInt::New();
  arr1->alloc(7,2);
  int *tmp=arr1->getPointer();
  const int arr1Ref[14]={0,10,1,11,2,12,3,13,4,14,5,15,6,16};
  std::copy(arr1Ref,arr1Ref+14,tmp);
  CPPUNIT_ASSERT_EQUAL(7,(int)arr1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)arr1->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(arr1Ref,arr1Ref+14,arr1->getConstPointer()));
  DataArrayInt *arr2=arr1->subArray(3);
  CPPUNIT_ASSERT_EQUAL(4,(int)arr2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)arr2->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(arr1Ref+6,arr1Ref+14,arr2->getConstPointer()));
  arr2->decrRef();
  DataArrayInt *arr3=arr1->subArray(2,5);
  CPPUNIT_ASSERT_EQUAL(3,(int)arr3->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)arr3->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(arr1Ref+4,arr1Ref+10,arr3->getConstPointer()));
  arr1->decrRef();
  arr3->decrRef();
  //
  DataArrayDouble *arr4=DataArrayDouble::New();
  arr4->alloc(7,2);
  double *tmp2=arr4->getPointer();
  const double arr4Ref[14]={0.8,10.8,1.9,11.9,2.1,12.1,3.2,13.2,4.3,14.3,5.4,15.4,6.5,16.5};
  std::copy(arr4Ref,arr4Ref+14,tmp2);
  CPPUNIT_ASSERT_EQUAL(7,(int)arr4->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)arr4->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(arr4Ref,arr4Ref+14,arr4->getConstPointer()));
  DataArrayDouble *arr5=arr4->subArray(3);
  CPPUNIT_ASSERT_EQUAL(4,(int)arr5->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)arr5->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(arr4Ref+6,arr4Ref+14,arr5->getConstPointer()));
  arr5->decrRef();
  DataArrayDouble *arr6=arr4->subArray(2,5);
  CPPUNIT_ASSERT_EQUAL(3,(int)arr6->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)arr6->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(arr4Ref+4,arr4Ref+10,arr6->getConstPointer()));
  arr4->decrRef();
  arr6->decrRef();
}

void MEDCouplingBasicsTest1::testMesh()
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
  CPPUNIT_ASSERT_EQUAL(MEDCouplingMesh::GetNumberOfNodesOfGeometricType(INTERP_KERNEL::NORM_TRI3),3);
  CPPUNIT_ASSERT(MEDCouplingMesh::IsStaticGeometricType(INTERP_KERNEL::NORM_TRI3));
  CPPUNIT_ASSERT(MEDCouplingMesh::IsLinearGeometricType(INTERP_KERNEL::NORM_TRI3));
  CPPUNIT_ASSERT_EQUAL(MEDCouplingMesh::GetDimensionOfGeometricType(INTERP_KERNEL::NORM_TRI3),2);
  CPPUNIT_ASSERT_EQUAL(std::string(MEDCouplingMesh::GetReprOfGeometricType(INTERP_KERNEL::NORM_TRI3)),std::string("NORM_TRI3"));
  CPPUNIT_ASSERT_THROW(MEDCouplingMesh::GetNumberOfNodesOfGeometricType(INTERP_KERNEL::NORM_POLYGON),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(!MEDCouplingMesh::IsStaticGeometricType(INTERP_KERNEL::NORM_POLYGON));
  CPPUNIT_ASSERT(MEDCouplingMesh::IsLinearGeometricType(INTERP_KERNEL::NORM_POLYGON));
  CPPUNIT_ASSERT_EQUAL(MEDCouplingMesh::GetDimensionOfGeometricType(INTERP_KERNEL::NORM_POLYGON),2);
  CPPUNIT_ASSERT_EQUAL(std::string(MEDCouplingMesh::GetReprOfGeometricType(INTERP_KERNEL::NORM_POLYGON)),std::string("NORM_POLYGON"));
  CPPUNIT_ASSERT_EQUAL(MEDCouplingMesh::GetNumberOfNodesOfGeometricType(INTERP_KERNEL::NORM_TRI6),6);
  CPPUNIT_ASSERT(MEDCouplingMesh::IsStaticGeometricType(INTERP_KERNEL::NORM_TRI6));
  CPPUNIT_ASSERT(!MEDCouplingMesh::IsLinearGeometricType(INTERP_KERNEL::NORM_TRI6));
  CPPUNIT_ASSERT_EQUAL(MEDCouplingMesh::GetDimensionOfGeometricType(INTERP_KERNEL::NORM_TRI6),2);
  CPPUNIT_ASSERT_EQUAL(std::string(MEDCouplingMesh::GetReprOfGeometricType(INTERP_KERNEL::NORM_TRI6)),std::string("NORM_TRI6"));
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(8);
  const int *curConn=tab4;
  for(int i=0;i<nbOfCells;i++,curConn+=4)
    mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,curConn);
  mesh->finishInsertingCells();
  CPPUNIT_ASSERT_EQUAL((std::size_t)30,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(nbOfCells,(int)mesh->getNumberOfCells());
  //test 0 - no copy no ownership
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->useArray(coords,false,CPP_DEALLOC,nbOfNodes,3);
  mesh->setCoords(myCoords);
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  CPPUNIT_ASSERT_EQUAL(nbOfCells,(int)mesh->getNumberOfCells());
  mesh->checkConsistencyLight();
  //test 1 - no copy ownership C++
  myCoords=DataArrayDouble::New();
  double *tmp=new double[3*nbOfNodes];
  std::copy(coords,coords+3*nbOfNodes,tmp);
  myCoords->useArray(tmp,true,CPP_DEALLOC,nbOfNodes,3);
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  CPPUNIT_ASSERT_EQUAL(nbOfCells,(int)mesh->getNumberOfCells());
  mesh->checkConsistencyLight();
  //test 2 - no copy ownership C
  myCoords=DataArrayDouble::New();
  tmp=(double *)malloc(3*nbOfNodes*sizeof(double));
  std::copy(coords,coords+3*nbOfNodes,tmp);
  myCoords->useArray(tmp,true,C_DEALLOC,nbOfNodes,3);
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  CPPUNIT_ASSERT_EQUAL(nbOfNodes,mesh->getNumberOfNodes());
  mesh->checkConsistencyLight();
  //test 3 - copy.
  myCoords=DataArrayDouble::New();
  myCoords->alloc(nbOfNodes,3);
  tmp=myCoords->getPointer();
  std::copy(coords,coords+3*nbOfNodes,tmp);
  // test 3 bis deepcopy
  DataArrayDouble *myCoords2=DataArrayDouble::New();
  *myCoords2=*myCoords;
  myCoords2->decrRef();
  //
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  CPPUNIT_ASSERT_EQUAL(nbOfNodes,mesh->getNumberOfNodes());
  mesh->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  // test clone not recursively
  MEDCouplingUMesh *mesh2=mesh->clone(false);
  CPPUNIT_ASSERT(mesh2!=mesh);
  mesh2->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(nbOfCells,(int)mesh2->getNumberOfCells());
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
  mesh3->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(nbOfCells,(int)mesh3->getNumberOfCells());
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
  std::fill(tmp,tmp+9*nbOfCells,7.);
  //content of field changed -> declare it.
  fieldOnCells->declareAsNew();
  fieldOnCells->checkConsistencyLight();
  // testing clone of fields - no recursive
  MEDCouplingFieldDouble *fieldOnCells2=fieldOnCells->clone(false);
  CPPUNIT_ASSERT(fieldOnCells2!=fieldOnCells);
  fieldOnCells2->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(nbOfCells,(int)fieldOnCells2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(9,(int)fieldOnCells2->getNumberOfComponents());
  CPPUNIT_ASSERT(fieldOnCells2->getArray()==fieldOnCells->getArray());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,fieldOnCells2->getArray()->getIJ(3,7),1e-14);
  CPPUNIT_ASSERT(fieldOnCells2->getMesh()==fieldOnCells->getMesh());
  // testing clone of fields - recursive
  MEDCouplingFieldDouble *fieldOnCells3=fieldOnCells->clone(true);
  CPPUNIT_ASSERT(fieldOnCells3!=fieldOnCells);
  fieldOnCells3->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(nbOfCells,(int)fieldOnCells3->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(9,(int)fieldOnCells3->getNumberOfComponents());
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

void MEDCouplingBasicsTest1::testMeshPointsCloud()
{
  double targetCoords[27]={-0.3,-0.3,0.5, 0.2,-0.3,1., 0.7,-0.3,1.5, -0.3,0.2,0.5, 0.2,0.2,1., 0.7,0.2,1.5, -0.3,0.7,0.5, 0.2,0.7,1., 0.7,0.7,1.5};
  const int targetConn[]={0,1,2,3,4,5,7,6};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(0);
  targetMesh->allocateCells(8);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT1,1,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT1,1,targetConn+1);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT1,1,targetConn+2);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT1,1,targetConn+3);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT1,1,targetConn+4);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT1,1,targetConn+5);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT1,1,targetConn+6);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT1,1,targetConn+7);
  targetMesh->finishInsertingCells();
  CPPUNIT_ASSERT_THROW(targetMesh->checkConsistencyLight(),INTERP_KERNEL::Exception);
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(9,3);
  std::copy(targetCoords,targetCoords+27,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  //
  targetMesh->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(3,targetMesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(8,(int)targetMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(9,targetMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(0,targetMesh->getMeshDimension());
  //
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest1::testMeshM1D()
{
  MEDCouplingUMesh *meshM1D=MEDCouplingUMesh::New();
  CPPUNIT_ASSERT_THROW(meshM1D->getMeshDimension(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(meshM1D->getNumberOfNodes(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(meshM1D->getNumberOfCells(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(meshM1D->setMeshDimension(-2),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(meshM1D->setMeshDimension(-10),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(meshM1D->checkConsistencyLight(),INTERP_KERNEL::Exception);
  meshM1D->setMeshDimension(-1);
  meshM1D->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(-1,meshM1D->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(1,(int)meshM1D->getNumberOfCells());
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
  std::fill(tmp,tmp+6,7.);
  fieldOnCells->checkConsistencyLight();
  //
  fieldOnCells->decrRef();
  meshM1D->decrRef();
}

void MEDCouplingBasicsTest1::testDeepCopy()
{
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(5,3);
  std::fill(array->getPointer(),array->getPointer()+5*3,7.);
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
  std::fill(array3->getPointer(),array3->getPointer()+5*3,17);
  CPPUNIT_ASSERT_EQUAL(17,array3->getIJ(3,2));
  int *tmp3=array3->getPointer();
  DataArrayInt *array4=array3->deepCopy();
  int *tmp4=array4->getPointer();
  CPPUNIT_ASSERT(tmp3!=tmp4);
  array3->decrRef();
  CPPUNIT_ASSERT_EQUAL(17,array4->getIJ(3,2));
  array4->decrRef();
}

void MEDCouplingBasicsTest1::testRevNodal()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  DataArrayInt *revNodal=DataArrayInt::New();
  DataArrayInt *revNodalIndx=DataArrayInt::New();
  //
  mesh->getReverseNodalConnectivity(revNodal,revNodalIndx);
  const int revNodalExpected[18]={0,0,1,1,2,0,3,0,1,2,3,4,2,4,3,3,4,4};
  const int revNodalIndexExpected[10]={0,1,3,5,7,12,14,15,17,18};
  CPPUNIT_ASSERT_EQUAL((std::size_t)18,revNodal->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)10,revNodalIndx->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(revNodalExpected,revNodalExpected+18,revNodal->getPointer()));
  CPPUNIT_ASSERT(std::equal(revNodalIndexExpected,revNodalIndexExpected+10,revNodalIndx->getPointer()));
  //
  revNodal->decrRef();
  revNodalIndx->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest1::testConvertToPolyTypes()
{
  ////// 2D
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  //
  const int elts[2]={1,3};
  std::vector<int> eltsV(elts,elts+2);
  mesh->convertToPolyTypes(&eltsV[0],&eltsV[0]+eltsV.size());
  mesh->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(5,(int)mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(23,(int)mesh->getNodalConnectivity()->getNumberOfTuples());
  const int *pt=mesh->getNodalConnectivity()->getConstPointer();
  const int expected1[23]={4, 0, 3, 4, 1, 5, 1, 4, 2, 3, 4, 5, 2, 5, 6, 7, 4, 3, 4, 7, 8, 5, 4};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+23,pt));
  //
  mesh->decrRef();
  ////// 3D
  mesh=build3DTargetMesh_1();
  mesh->convertToPolyTypes(&eltsV[0],&eltsV[0]+eltsV.size());
  mesh->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(8,(int)mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(114,(int)mesh->getNodalConnectivity()->getNumberOfTuples());
  mesh->convertToPolyTypes(&eltsV[0],&eltsV[0]+eltsV.size());
  mesh->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(8,(int)mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(114,(int)mesh->getNodalConnectivity()->getNumberOfTuples());
  //
  mesh->decrRef();
}

void MEDCouplingBasicsTest1::testDescConn2D()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  DataArrayInt *desc=DataArrayInt::New();
  DataArrayInt *descIndx=DataArrayInt::New();
  DataArrayInt *revDesc=DataArrayInt::New();
  DataArrayInt *revDescIndx=DataArrayInt::New();
  //
  MEDCouplingUMesh *mesh2=mesh->buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  mesh2->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(1,mesh2->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(13,(int)mesh2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL((std::size_t)14,revDescIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(14,(int)revDescIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)6,descIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(6,(int)descIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)18,desc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(18,(int)desc->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)18,revDesc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(18,(int)revDesc->getNumberOfTuples());
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
  mesh->convertToPolyTypes(&eltsV[0],&eltsV[0]+eltsV.size());
  mesh->checkConsistencyLight();
  //
  desc=DataArrayInt::New();
  descIndx=DataArrayInt::New();
  revDesc=DataArrayInt::New();
  revDescIndx=DataArrayInt::New();
  //
  mesh2=mesh->buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  mesh2->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(1,mesh2->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(13,(int)mesh2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL((std::size_t)14,revDescIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(14,(int)revDescIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)6,descIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(6,(int)descIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)18,desc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(18,(int)desc->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)18,revDesc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(18,(int)revDesc->getNumberOfTuples());
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

void MEDCouplingBasicsTest1::testDescConn3D()
{
  MEDCouplingUMesh *mesh=build3DTargetMesh_1();
  DataArrayInt *desc=DataArrayInt::New();
  DataArrayInt *descIndx=DataArrayInt::New();
  DataArrayInt *revDesc=DataArrayInt::New();
  DataArrayInt *revDescIndx=DataArrayInt::New();
  //
  MEDCouplingUMesh *mesh2=mesh->buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  mesh2->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(2,mesh2->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(36,(int)mesh2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL((std::size_t)37,revDescIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(37,(int)revDescIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)9,descIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(9,(int)descIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)48,desc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(48,(int)desc->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)48,revDesc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(48,(int)revDesc->getNumberOfTuples());
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
  mesh->convertToPolyTypes(&eltsV[0],&eltsV[0]+eltsV.size());
  mesh->checkConsistencyLight();
  desc=DataArrayInt::New();
  descIndx=DataArrayInt::New();
  revDesc=DataArrayInt::New();
  revDescIndx=DataArrayInt::New();
  mesh2=mesh->buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  mesh2->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(2,mesh2->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(36,(int)mesh2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL((std::size_t)37,revDescIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(37,(int)revDescIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)9,descIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(9,(int)descIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)48,desc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(48,(int)desc->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)48,revDesc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(48,(int)revDesc->getNumberOfTuples());
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

void MEDCouplingBasicsTest1::testFindBoundaryNodes()
{
  MEDCouplingUMesh *mesh=build3DTargetMesh_1();
  DataArrayInt *boundaryNodes=mesh->findBoundaryNodes();
  CPPUNIT_ASSERT_EQUAL(26,(int)boundaryNodes->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)boundaryNodes->getNumberOfComponents());
  const int expected1[26]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+26,boundaryNodes->begin()));
  boundaryNodes->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest1::testBoundaryMesh()
{
  MEDCouplingUMesh *mesh=build3DTargetMesh_1();
  MEDCouplingPointSet *mesh2=mesh->buildBoundaryMesh(false);
  CPPUNIT_ASSERT_EQUAL(24,(int)mesh2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(26,mesh2->getNumberOfNodes());
  mesh2->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest1::testBuildPartOfMySelf()
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
  CPPUNIT_ASSERT_EQUAL(2,(int)mesh->getAllGeoTypes().size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,*mesh->getAllGeoTypes().begin());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*(++(mesh->getAllGeoTypes().begin())));
  CPPUNIT_ASSERT_EQUAL(1,(int)subMesh->getAllGeoTypes().size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*subMesh->getAllGeoTypes().begin());
  CPPUNIT_ASSERT(name=="Toto");
  CPPUNIT_ASSERT(mesh->getCoords()==subMesh->getCoords());
  CPPUNIT_ASSERT_EQUAL(2,(int)subMesh->getNumberOfCells());
  const int subConn[10]={4,0,3,4,1,4,7,8,5,4};
  const int subConnIndex[3]={0,5,10};
  CPPUNIT_ASSERT_EQUAL((std::size_t)10,subMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)3,subMesh->getNodalConnectivityIndex()->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(subConn,subConn+10,subMesh->getNodalConnectivity()->getPointer()));
  CPPUNIT_ASSERT(std::equal(subConnIndex,subConnIndex+3,subMesh->getNodalConnectivityIndex()->getPointer()));
  subMesh->decrRef();
  //
  subMeshSimple=mesh->buildPartOfMySelf(tab2,tab2+3,true);
  subMesh=dynamic_cast<MEDCouplingUMesh *>(subMeshSimple);
  CPPUNIT_ASSERT(subMesh);
  name=subMesh->getName();
  CPPUNIT_ASSERT_EQUAL(2,(int)subMesh->getAllGeoTypes().size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,*subMesh->getAllGeoTypes().begin());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*(++(subMesh->getAllGeoTypes().begin())));
  CPPUNIT_ASSERT(name=="Toto");
  CPPUNIT_ASSERT(mesh->getCoords()==subMesh->getCoords());
  CPPUNIT_ASSERT_EQUAL(3,(int)subMesh->getNumberOfCells());
  const int subConn2[14]={4,0,3,4,1,3,4,5,2,4,6,7,4,3};
  const int subConnIndex2[4]={0,5,9,14};
  CPPUNIT_ASSERT_EQUAL((std::size_t)14,subMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)4,subMesh->getNodalConnectivityIndex()->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(subConn2,subConn2+14,subMesh->getNodalConnectivity()->getPointer()));
  CPPUNIT_ASSERT(std::equal(subConnIndex2,subConnIndex2+4,subMesh->getNodalConnectivityIndex()->getPointer()));
  const int tab3[3]={0,1,2};
  MEDCouplingPointSet *subMeshSimple2=subMeshSimple->buildPartOfMySelf(tab3,tab3+3,true);
  subMesh->decrRef();
  name=subMeshSimple2->getName();
  CPPUNIT_ASSERT(name=="Toto");
  subMeshSimple2->decrRef();
  //
  mesh->decrRef();
}

void MEDCouplingBasicsTest1::testBuildPartOfMySelfNode()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  const int tab1[4]={5,7,8,4};
  MEDCouplingPointSet *subMeshSimple=mesh->buildPartOfMySelfNode(tab1,tab1+4,true);
  MEDCouplingUMesh *subMesh=dynamic_cast<MEDCouplingUMesh *>(subMeshSimple);
  CPPUNIT_ASSERT(subMesh);
  CPPUNIT_ASSERT_EQUAL(1,(int)subMesh->getAllGeoTypes().size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*subMesh->getAllGeoTypes().begin());
  CPPUNIT_ASSERT_EQUAL(1,(int)subMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL((std::size_t)5,subMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)2,subMesh->getNodalConnectivityIndex()->getNbOfElems());
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
  CPPUNIT_ASSERT_EQUAL(2,(int)subMesh->getAllGeoTypes().size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,*subMesh->getAllGeoTypes().begin());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*(++subMesh->getAllGeoTypes().begin()));
  CPPUNIT_ASSERT_EQUAL(3,(int)subMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL((std::size_t)14,subMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)4,subMesh->getNodalConnectivityIndex()->getNbOfElems());
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
  CPPUNIT_ASSERT_EQUAL(2,(int)subMesh->getAllGeoTypes().size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,*subMesh->getAllGeoTypes().begin());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*(++subMesh->getAllGeoTypes().begin()));
  CPPUNIT_ASSERT_EQUAL(3,(int)subMesh->getNumberOfCells());
  subMeshSimple->decrRef();
  //
  mesh->decrRef();
}

void MEDCouplingBasicsTest1::testZipCoords()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  CPPUNIT_ASSERT_EQUAL(2,(int)mesh->getAllGeoTypes().size());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(9,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(5,(int)mesh->getNumberOfCells());
  std::vector<int> oldConn(mesh->getNodalConnectivity()->getNbOfElems());
  std::vector<int> oldConnIndex(mesh->getNumberOfCells()+1);
  std::copy(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+oldConn.size(),oldConn.begin());
  std::copy(mesh->getNodalConnectivityIndex()->getPointer(),mesh->getNodalConnectivityIndex()->getPointer()+mesh->getNumberOfCells()+1,oldConnIndex.begin());
  DataArrayDouble *oldCoords=mesh->getCoords();
  oldCoords->incrRef();
  mesh->zipCoords();
  CPPUNIT_ASSERT_EQUAL(2,(int)mesh->getAllGeoTypes().size());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(9,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(5,(int)mesh->getNumberOfCells());
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
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*subMesh->getAllGeoTypes().begin());
  CPPUNIT_ASSERT_EQUAL(2,(int)subMesh->getNumberOfCells());
  const int subConn[10]={4,0,2,3,1,4,5,6,4,3};
  const int subConnIndex[3]={0,5,10};
  CPPUNIT_ASSERT_EQUAL(7,subMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL((std::size_t)10,subMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)3,subMesh->getNodalConnectivityIndex()->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(subConn,subConn+10,subMesh->getNodalConnectivity()->getPointer()));
  CPPUNIT_ASSERT(std::equal(subConnIndex,subConnIndex+3,subMesh->getNodalConnectivityIndex()->getPointer()));
  subMesh->decrRef();
  //
  subMeshPtSet=mesh->buildPartOfMySelf(tab1,tab1+2,false);
  subMesh=dynamic_cast<MEDCouplingUMesh *>(subMeshPtSet);
  CPPUNIT_ASSERT(subMesh);
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,*subMesh->getAllGeoTypes().begin());
  CPPUNIT_ASSERT_EQUAL(2,(int)subMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(7,subMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL((std::size_t)10,subMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)3,subMesh->getNodalConnectivityIndex()->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(subConn,subConn+10,subMesh->getNodalConnectivity()->getPointer()));
  CPPUNIT_ASSERT(std::equal(subConnIndex,subConnIndex+3,subMesh->getNodalConnectivityIndex()->getPointer()));
  subMesh->decrRef();
  //
  mesh->decrRef();
}

void MEDCouplingBasicsTest1::testZipConnectivity()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  MEDCouplingUMesh *m2=build2DTargetMesh_1();
  int cells1[3]={2,3,4};
  MEDCouplingPointSet *m3_1=m2->buildPartOfMySelf(cells1,cells1+3,true);
  MEDCouplingUMesh *m3=dynamic_cast<MEDCouplingUMesh *>(m3_1);
  CPPUNIT_ASSERT(m3);
  m2->decrRef();
  MEDCouplingUMesh *m4=build2DSourceMesh_1();
  MEDCouplingUMesh *m5=MEDCouplingUMesh::MergeUMeshes(m1,m3);
  m1->decrRef();
  m3->decrRef();
  MEDCouplingUMesh *m6=MEDCouplingUMesh::MergeUMeshes(m5,m4);
  m4->decrRef();
  m5->decrRef();
  //
  bool areNodesMerged;
  int newNbOfNodes;
  CPPUNIT_ASSERT_EQUAL(10,(int)m6->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(22,m6->getNumberOfNodes());
  DataArrayInt *arr=m6->mergeNodes(1e-13,areNodesMerged,newNbOfNodes);
  arr->decrRef();
  CPPUNIT_ASSERT(areNodesMerged);
  CPPUNIT_ASSERT_EQUAL(10,(int)m6->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(9,m6->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(9,newNbOfNodes);
  //
  arr=m6->zipConnectivityTraducer(0);
  CPPUNIT_ASSERT_EQUAL(7,(int)m6->getNumberOfCells());
  arr->decrRef();
  MEDCouplingUMesh *m7=m6->clone(true);
  arr=m6->zipConnectivityTraducer(0);
  CPPUNIT_ASSERT(m7->isEqual(m6,1e-12));
  CPPUNIT_ASSERT_EQUAL(7,(int)m6->getNumberOfCells());
  arr->decrRef();
  //
  m7->decrRef();
  m6->decrRef();
}

void MEDCouplingBasicsTest1::testEqualMesh()
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

void MEDCouplingBasicsTest1::testEqualFieldDouble()
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

void MEDCouplingBasicsTest1::testNatureChecking()
{
  MEDCouplingFieldDouble *field=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  field->setNature(ExtensiveMaximum);
  field->setNature(IntensiveMaximum);
  field->setNature(ExtensiveConservation);
  field->decrRef();
  field=MEDCouplingFieldDouble::New(ON_NODES,NO_TIME);
  field->setNature(IntensiveMaximum);
  CPPUNIT_ASSERT_THROW(field->setNature(ExtensiveMaximum),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(field->setNature(ExtensiveConservation),INTERP_KERNEL::Exception);
  field->decrRef();
}

void MEDCouplingBasicsTest1::testBuildSubMeshData()
{
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  //check buildSubMesh on field on cells
  MEDCouplingFieldDouble *fieldCells=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  fieldCells->setMesh(targetMesh);
  const int elts[3]={1,2,4};
  DataArrayInt *di;
  MEDCouplingMesh *ret1=fieldCells->buildSubMeshData(elts,elts+3,di);
  CPPUNIT_ASSERT_EQUAL(3,(int)ret1->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(9,ret1->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,(int)di->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)di->getNumberOfComponents());
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
  CPPUNIT_ASSERT_EQUAL(3,(int)ret2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(6,ret2->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(6,(int)di->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)di->getNumberOfComponents());
  toCheck=di->getConstPointer();
  const int expected[6]={1,2,4,5,7,8};
  CPPUNIT_ASSERT(std::equal(expected,expected+6,toCheck));
  ret2->decrRef();
  di->decrRef();
  fieldNodes->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest1::testExtrudedMesh1()
{
  MEDCouplingUMesh *mesh2D=0;
  MEDCouplingUMesh *mesh3D=build3DExtrudedUMesh_1(mesh2D);
  MEDCouplingMappedExtrudedMesh *ext=MEDCouplingMappedExtrudedMesh::New(mesh3D,mesh2D,1);
  CPPUNIT_ASSERT_EQUAL(18,(int)ext->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(60,ext->getNumberOfNodes());
  DataArrayInt *ids3D=ext->getMesh3DIds();
  const int ids3DExpected[18]={5,4,3,2,1,0, 11,10,9,8,7,6, 17,16,15,14,13,12};
  CPPUNIT_ASSERT_EQUAL(18,(int)ids3D->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)ids3D->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(ids3DExpected,ids3DExpected+18,ids3D->getConstPointer()));
  MEDCouplingUMesh *mesh1D=ext->getMesh1D();
  CPPUNIT_ASSERT_EQUAL(4,mesh1D->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,(int)mesh1D->getNumberOfCells());
  const double mesh1DExpected[12]={0.66666666666666663, 1.4583333333333333, 0, 0.66666666666666663, 1.4583333333333333, 1, 0.66666666666666663, 1.4583333333333333, 2, 0.66666666666666663, 1.4583333333333333, 3};
  DataArrayDouble *mesh1DCoords=mesh1D->getCoords();
  CPPUNIT_ASSERT_EQUAL(4,(int)mesh1DCoords->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,(int)mesh1DCoords->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(mesh1DExpected,mesh1DExpected+12,mesh1DCoords->getConstPointer()));
  DataArrayInt *conn1D=mesh1D->getNodalConnectivity();
  CPPUNIT_ASSERT_EQUAL(9,(int)conn1D->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)conn1D->getNumberOfComponents());
  const int conn1DExpected[9]={1,0,1,1,1,2,1,2,3};
  CPPUNIT_ASSERT(std::equal(conn1DExpected,conn1DExpected+9,conn1D->getConstPointer()));
  ext->decrRef();
  mesh3D->decrRef();
  mesh2D->decrRef();
}

void MEDCouplingBasicsTest1::testExtrudedMesh2()
{
  MEDCouplingUMesh *mN,*mTT,*mTF;
  build3DExtrudedUMesh_2(mN,mTT,mTF);
  //
  bool b=false;
  int newNbOfNodes;
  DataArrayInt *da=mTT->mergeNodes(1e-12,b,newNbOfNodes);
  CPPUNIT_ASSERT(b);
  da->decrRef();
  std::vector<int> n;
  double pt[3]={300.,300.,0.};
  double v[3]={0.,0.,2.};
  mTT->findNodesOnPlane(pt,v,1e-12,n);
  CPPUNIT_ASSERT_EQUAL(43,(int)n.size());
  MEDCouplingUMesh *mTT3dSurf=(MEDCouplingUMesh *)mTT->buildFacePartOfMySelfNode(&n[0],&n[0]+n.size(),true);
  MEDCouplingMappedExtrudedMesh *meTT=MEDCouplingMappedExtrudedMesh::New(mTT,mTT3dSurf,0);
  CPPUNIT_ASSERT_EQUAL(200,(int)meTT->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(10,(int)meTT->getMesh2D()->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(20,(int)meTT->getMesh1D()->getNumberOfCells());
  mTT3dSurf->decrRef();
  //
  b=false;
  da=mN->mergeNodes(1e-12,b,newNbOfNodes);
  da->decrRef();
  CPPUNIT_ASSERT(!b);
  n.clear();
  mN->findNodesOnPlane(pt,v,1e-12,n);
  CPPUNIT_ASSERT_EQUAL(30,(int)n.size());
  MEDCouplingUMesh *mN3dSurf=(MEDCouplingUMesh *)mN->buildFacePartOfMySelfNode(&n[0],&n[0]+n.size(),true);
  MEDCouplingMappedExtrudedMesh *meN=MEDCouplingMappedExtrudedMesh::New(mN,mN3dSurf,0);
  CPPUNIT_ASSERT_EQUAL(40,(int)meN->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(20,(int)meN->getMesh2D()->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(2,(int)meN->getMesh1D()->getNumberOfCells());
  mN3dSurf->decrRef();
  //
  b=false;
  da=mTF->mergeNodes(1e-12,b,newNbOfNodes);
  da->decrRef();
  CPPUNIT_ASSERT(!b);
  n.clear();
  mTF->findNodesOnPlane(pt,v,1e-12,n);
  CPPUNIT_ASSERT_EQUAL(27,(int)n.size());
  MEDCouplingUMesh *mTF3dSurf=(MEDCouplingUMesh *)mTF->buildFacePartOfMySelfNode(&n[0],&n[0]+n.size(),true);
  MEDCouplingMappedExtrudedMesh *meTF=MEDCouplingMappedExtrudedMesh::New(mTF,mTF3dSurf,0);
  CPPUNIT_ASSERT_EQUAL(340,(int)meTF->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(17,(int)meTF->getMesh2D()->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(20,(int)meTF->getMesh1D()->getNumberOfCells());
  mTF3dSurf->decrRef();
  //
  meTT->decrRef();
  meN->decrRef();
  meTF->decrRef();
  //
  mN->decrRef();
  mTT->decrRef();
  mTF->decrRef();
}

/*!
 * This test check MEDCouplingUMesh::buildExtrudedMesh method.
 */
void MEDCouplingBasicsTest1::testExtrudedMesh3()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  m1->changeSpaceDimension(3);
  MEDCouplingUMesh *m2=buildCU1DMesh_U();
  m2->changeSpaceDimension(3);
  double center[3]={0.,0.,0.};
  double vector[3]={0,1,0};
  m2->rotate(center,vector,-M_PI/2.);
  MEDCouplingUMesh *m3=m1->buildExtrudedMesh(m2,0);
  //
  MEDCouplingMappedExtrudedMesh *m4=MEDCouplingMappedExtrudedMesh::New(m3,m1,0);
  CPPUNIT_ASSERT_EQUAL(15,(int)m4->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(5,(int)m4->getMesh2D()->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(3,(int)m4->getMesh1D()->getNumberOfCells());
  const int *m3DIds=m4->getMesh3DIds()->getConstPointer();
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_EQUAL(i,m3DIds[i]);
  m4->decrRef();
  //some random in cells to check that extrusion alg find it correctly
  const int expected1[15]={1,3,2,0,6,5,7,10,11,8,12,9,14,13,4};
  m3->renumberCells(expected1,false);
  m4=MEDCouplingMappedExtrudedMesh::New(m3,m1,0);
  CPPUNIT_ASSERT_EQUAL(15,(int)m4->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(5,(int)m4->getMesh2D()->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(3,(int)m4->getMesh1D()->getNumberOfCells());
  m3DIds=m4->getMesh3DIds()->getConstPointer();
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],m3DIds[i]);
  m4->decrRef();
  m3->decrRef();
  //play with polygons and polyedrons
  std::vector<int> cells(2); cells[0]=2; cells[1]=3;
  m1->convertToPolyTypes(&cells[0],&cells[0]+cells.size());
  m3=m1->buildExtrudedMesh(m2,0);
  CPPUNIT_ASSERT_EQUAL((int)INTERP_KERNEL::NORM_HEXA8,(int)m3->getTypeOfCell(0));
  CPPUNIT_ASSERT_EQUAL((int)INTERP_KERNEL::NORM_PENTA6,(int)m3->getTypeOfCell(1));
  CPPUNIT_ASSERT_EQUAL((int)INTERP_KERNEL::NORM_POLYHED,(int)m3->getTypeOfCell(2));
  CPPUNIT_ASSERT_EQUAL((int)INTERP_KERNEL::NORM_POLYHED,(int)m3->getTypeOfCell(3));
  CPPUNIT_ASSERT_EQUAL((int)INTERP_KERNEL::NORM_HEXA8,(int)m3->getTypeOfCell(4));
  m3->renumberCells(expected1,false);
  m4=MEDCouplingMappedExtrudedMesh::New(m3,m1,0);
  CPPUNIT_ASSERT_EQUAL(15,(int)m4->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(5,(int)m4->getMesh2D()->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(3,(int)m4->getMesh1D()->getNumberOfCells());
  m3DIds=m4->getMesh3DIds()->getConstPointer();
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],m3DIds[i]);
  m4->decrRef();
  m3->decrRef();
  //
  m2->decrRef();
  m1->decrRef();
}

/*!
 * This test check MEDCouplingUMesh::buildExtrudedMesh method, but also, MEDCouplingMappedExtrudedMesh following methods :
 * getCellContainingPoint getMeasureField getNodeIdsOfCell getCoordinateOfNode getTypeOfCell build3DUnstructuredMesh.
 */
void MEDCouplingBasicsTest1::testExtrudedMesh4()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  std::vector<int> cells(2); cells[0]=2; cells[1]=4;
  m1->convertToPolyTypes(&cells[0],&cells[0]+cells.size());
  m1->changeSpaceDimension(3);
  MEDCouplingUMesh *m2=buildCU1DMesh_U();
  m2->changeSpaceDimension(3);
  double center[3]={0.,0.,0.};
  double vector[3]={0.,1.,0.};
  m2->rotate(center,vector,-M_PI/2.);
  m1->zipCoords();
  MEDCouplingUMesh *m3=m1->buildExtrudedMesh(m2,0);
  const int expected1[15]= {1,3,2,0,6,5,7,10,11,8,12,9,14,13,4};
  const int rexpected1[15]={3, 0, 2, 1, 14, 5, 4, 6, 9, 11, 7, 8, 10, 13, 12};
  m3->renumberCells(expected1,false);
  MEDCouplingMappedExtrudedMesh *m4=MEDCouplingMappedExtrudedMesh::New(m3,m1,0);
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_HEXA8,m4->getTypeOfCell(0));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_HEXA8,m4->getTypeOfCell(1));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_POLYHED,m4->getTypeOfCell(2));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_PENTA6,m4->getTypeOfCell(7));
  MEDCouplingFieldDouble *f=m4->getMeasureField(true);
  DataArrayDouble *arr=f->getArray();
  CPPUNIT_ASSERT_EQUAL(15,(int)arr->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
  const double *arrPtr=arr->getConstPointer();
  const double expected2[15]={0.075,0.0375,0.0375,0.075,0.075,   0.1125,0.05625,0.05625,0.1125,0.1125,   0.0625,0.03125,0.03125,0.0625,0.0625};
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[rexpected1[i]],arrPtr[i],1e-16);
  f->decrRef();
  MEDCouplingUMesh *m5=m4->build3DUnstructuredMesh();
  m5->zipCoords();
  CPPUNIT_ASSERT(m5->isEqual(m3,1e-12));
  f=m5->getMeasureField(true);
  arr=f->getArray();
  arrPtr=arr->getConstPointer();
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[rexpected1[i]],arrPtr[i],1e-16);
  f->decrRef();
  m5->decrRef();
  //
  m4->decrRef();
  m3->decrRef();
  m2->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest1::testFindCommonNodes()
{
  DataArrayInt *comm,*commI;
  MEDCouplingUMesh *targetMesh=build3DTargetMesh_1();
  targetMesh->findCommonNodes(1e-10,-1,comm,commI);
  CPPUNIT_ASSERT_EQUAL(1,(int)commI->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,(int)comm->getNumberOfTuples());
  int newNbOfNodes;
  DataArrayInt *o2n=targetMesh->buildNewNumberingFromCommonNodesFormat(comm,commI,newNbOfNodes);
  CPPUNIT_ASSERT_EQUAL(27,newNbOfNodes);
  CPPUNIT_ASSERT_EQUAL(27,(int)o2n->getNumberOfTuples());
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
  targetMesh->findCommonNodes(1e-10,-1,comm,commI);
  CPPUNIT_ASSERT_EQUAL(3,(int)commI->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(6,(int)comm->getNumberOfTuples());
  const int commExpected[6]={1,27,28,29,23,30};
  const int commIExpected[3]={0,4,6};
  CPPUNIT_ASSERT(std::equal(commExpected,commExpected+6,comm->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(commIExpected,commIExpected+3,commI->getConstPointer()));
  o2n=targetMesh->buildNewNumberingFromCommonNodesFormat(comm,commI,newNbOfNodes);
  CPPUNIT_ASSERT_EQUAL(31,(int)o2n->getNumberOfTuples());
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
  o2n=targetMesh->mergeNodes(1e-10,areNodesMerged,newNbOfNodes);
  targetMesh->updateTime();
  CPPUNIT_ASSERT(time==targetMesh->getTimeOfThis());
  CPPUNIT_ASSERT(!areNodesMerged);
  targetMesh->decrRef();
  o2n->decrRef();
  //
  targetMesh=build3DTargetMeshMergeNode_1();
  time=targetMesh->getTimeOfThis();
  o2n=targetMesh->mergeNodes(1e-10,areNodesMerged,newNbOfNodes);
  targetMesh->updateTime();
  CPPUNIT_ASSERT(time!=targetMesh->getTimeOfThis());
  CPPUNIT_ASSERT(areNodesMerged);
  int connExp[72]={18,0,1,4,3,9,10,13,12, 18,1,2,5,4,10,11,14,13, 18,3,4,7,6,12,13,16,15,
                   18,4,5,8,7,13,14,17,16,
                   18,9,10,13,12,18,19,22,21, 18,10,11,14,13,19,20,23,22, 18,12,13,16,15,21,22,25,24,
                   18,13,14,17,16,22,23,26,25};
  CPPUNIT_ASSERT_EQUAL(72,(int)targetMesh->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(connExp,connExp+72,targetMesh->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT_EQUAL(27,(int)targetMesh->getCoords()->getNumberOfTuples());
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
  o2n=targetMesh->mergeNodes(1e-10,areNodesMerged,newNbOfNodes);
  CPPUNIT_ASSERT(time!=targetMesh->getTimeOfThis());
  CPPUNIT_ASSERT(areNodesMerged);
  CPPUNIT_ASSERT_EQUAL(9,targetMesh->getNumberOfNodes());
  int connExp2[23]={4,0,4,3,1, 3,1,3,2, 3,3,5,2, 4,4,6,7,3, 4,7,8,5,3};
  CPPUNIT_ASSERT_EQUAL(23,(int)targetMesh->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(connExp2,connExp2+23,targetMesh->getNodalConnectivity()->getConstPointer()));
  double coordsExp2[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, 0.2,0.2, -0.3,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7};
  CPPUNIT_ASSERT_EQUAL(9,(int)targetMesh->getCoords()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(coordsExp2,coordsExp2+18,targetMesh->getCoords()->getConstPointer()));
  targetMesh->decrRef();
  o2n->decrRef();
}

void MEDCouplingBasicsTest1::testCheckButterflyCells()
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

void MEDCouplingBasicsTest1::testMergeMesh1()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  MEDCouplingUMesh *m2=build2DSourceMesh_1();
  const double vec[2]={1.,0.};
  m2->translate(vec);
  MEDCouplingMesh *m3=m1->mergeMyselfWith(m2);
  MEDCouplingUMesh *m3C=dynamic_cast<MEDCouplingUMesh *>(m3);
  CPPUNIT_ASSERT(m3C);
  m3->checkConsistencyLight();
  MEDCouplingUMesh *m4=build2DTargetMeshMerged_1();
  CPPUNIT_ASSERT(m3->isEqual(m4,1.e-12));
  m4->decrRef();
  bool isMerged;
  int newNbOfNodes;
  DataArrayInt *da=m3C->mergeNodes(1.e-12,isMerged,newNbOfNodes);
  CPPUNIT_ASSERT_EQUAL(11,m3C->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(11,newNbOfNodes);
  CPPUNIT_ASSERT(isMerged);
  da->decrRef();
  m3->decrRef();
  m1->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest1::testMergeMeshOnSameCoords1()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  MEDCouplingUMesh *m2=build2DTargetMesh_1();
  std::vector<int> cells(5);
  for(int i=0;i<5;i++)
    cells[i]=i;
  m2->convertToPolyTypes(&cells[0],&cells[0]+cells.size());
  m1->tryToShareSameCoords(*m2,1e-12);
  MEDCouplingUMesh *m3=build2DTargetMesh_1();
  m3->tryToShareSameCoords(*m2,1e-12);
  std::vector<const MEDCouplingUMesh *> meshes;
  meshes.push_back(m1); meshes.push_back(m2); meshes.push_back(m3);
  MEDCouplingUMesh *m4=MEDCouplingUMesh::MergeUMeshesOnSameCoords(meshes);
  m4->checkConsistencyLight();
  CPPUNIT_ASSERT(m4->getCoords()==m1->getCoords());
  CPPUNIT_ASSERT_EQUAL(15,(int)m4->getNumberOfCells());
  const int cells1[5]={0,1,2,3,4};
  MEDCouplingPointSet *m1_1=m4->buildPartOfMySelf(cells1,cells1+5,true);
  m1_1->setName(m1->getName().c_str());
  CPPUNIT_ASSERT(m1->isEqual(m1_1,1e-12));
  const int cells2[5]={5,6,7,8,9};
  MEDCouplingPointSet *m2_1=m4->buildPartOfMySelf(cells2,cells2+5,true);
  m2_1->setName(m2->getName().c_str());
  CPPUNIT_ASSERT(m2->isEqual(m2_1,1e-12));
  const int cells3[5]={10,11,12,13,14};
  MEDCouplingPointSet *m3_1=m4->buildPartOfMySelf(cells3,cells3+5,true);
  m3_1->setName(m3->getName().c_str());
  CPPUNIT_ASSERT(m3->isEqual(m3_1,1e-12));
  m1_1->decrRef(); m2_1->decrRef(); m3_1->decrRef();
  //
  m4->decrRef();
  m1->decrRef();
  m2->decrRef();
  m3->decrRef();
}

void MEDCouplingBasicsTest1::testMergeField1()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  MEDCouplingUMesh *m2=build2DSourceMesh_1();
  const double vec[2]={1.,0.};
  m2->translate(vec);
  MEDCouplingFieldDouble *f1=m1->getMeasureField(true);
  MEDCouplingFieldDouble *f2=m2->getMeasureField(true);
  MEDCouplingFieldDouble *f3=MEDCouplingFieldDouble::MergeFields(f1,f2);
  f3->checkConsistencyLight();
  MEDCouplingUMesh *m4=build2DTargetMeshMerged_1();
  CPPUNIT_ASSERT(f3->getMesh()->isEqual(m4,1.e-12));
  std::string name=f3->getName();
  CPPUNIT_ASSERT(name=="MeasureOfMesh_");
  CPPUNIT_ASSERT(f3->getTypeOfField()==ON_CELLS);
  CPPUNIT_ASSERT(f3->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_EQUAL(1,(int)f3->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(7,(int)f3->getNumberOfTuples());
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

void MEDCouplingBasicsTest1::testFillFromAnalytic()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  m->setTime(3.4,5,6); m->setTimeUnit("us");
  int a,b;
  MEDCouplingFieldDouble *f1=m->fillFromAnalytic(ON_CELLS,1,func1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.4,f1->getTime(a,b),1.e-14);
  CPPUNIT_ASSERT_EQUAL(5,a); CPPUNIT_ASSERT_EQUAL(6,b);
  CPPUNIT_ASSERT_EQUAL(std::string(f1->getTimeUnit()),std::string("us"));
  f1->checkConsistencyLight();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_CELLS);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_EQUAL(1,(int)f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,(int)f1->getNumberOfTuples());
  double values1[5]={-0.1,0.23333333333333336,0.56666666666666665,0.4,0.9};
  const double *tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+5,values1,values1,std::minus<double>());
  std::transform(values1,values1+5,values1,std::ptr_fun<double,double>(fabs));
  double max=*std::max_element(values1,values1+5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  //
  f1=m->fillFromAnalytic(ON_NODES,1,func1);
  f1->checkConsistencyLight();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_EQUAL(1,(int)f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,(int)f1->getNumberOfTuples());
  double values2[9]={-0.6,-0.1,0.4,-0.1,0.4,0.9,0.4,0.9,1.4};
  tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+9,values2,values2,std::minus<double>());
  std::transform(values2,values2+9,values2,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values2,values2+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  //
  f1=m->fillFromAnalytic(ON_NODES,2,func2);
  f1->checkConsistencyLight();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_EQUAL(2,(int)f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,(int)f1->getNumberOfTuples());
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
  f1->integral(true,values4);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,values4[0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,values4[1],1.e-12);
  f1->decrRef();
  //
  CPPUNIT_ASSERT_THROW(f1=m->fillFromAnalytic(ON_NODES,1,func3),INTERP_KERNEL::Exception);
  //
  m->decrRef();
}

void MEDCouplingBasicsTest1::testFillFromAnalytic2()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=m->fillFromAnalytic(ON_CELLS,1,"y+x");
  f1->checkConsistencyLight();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_CELLS);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_EQUAL(1,(int)f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,(int)f1->getNumberOfTuples());
  double values1[5]={-0.1,0.23333333333333336,0.56666666666666665,0.4,0.9};
  const double *tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+5,values1,values1,std::minus<double>());
  std::transform(values1,values1+5,values1,std::ptr_fun<double,double>(fabs));
  double max=*std::max_element(values1,values1+5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  //
  f1=m->fillFromAnalytic(ON_NODES,1,"y+2*x");
  f1->checkConsistencyLight();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_EQUAL(1,(int)f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,(int)f1->getNumberOfTuples());
  double values2[9]={-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1};
  tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+9,values2,values2,std::minus<double>());
  std::transform(values2,values2+9,values2,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values2,values2+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  f1=m->fillFromAnalytic(ON_NODES,1,"2.*x+y");
  f1->checkConsistencyLight();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_EQUAL(1,(int)f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,(int)f1->getNumberOfTuples());
  tmp=f1->getArray()->getConstPointer();
  double values2Bis[9]={-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1};
  std::transform(tmp,tmp+9,values2Bis,values2Bis,std::minus<double>());
  std::transform(values2,values2+9,values2Bis,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values2Bis,values2Bis+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  //
  f1=m->fillFromAnalytic(ON_NODES,2,"(x+y)*IVec+2*(x+y)*JVec");
  f1->checkConsistencyLight();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_EQUAL(2,(int)f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,(int)f1->getNumberOfTuples());
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
  f1->integral(true,values4);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,values4[0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,values4[1],1.e-12);
  f1->decrRef();
  //
  CPPUNIT_ASSERT_THROW(f1=m->fillFromAnalytic(ON_NODES,1,"1./(x-0.2)"),INTERP_KERNEL::Exception);
  //
  m->decrRef();
}

void MEDCouplingBasicsTest1::testApplyFunc()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=m->fillFromAnalytic(ON_NODES,2,func2);
  f1->checkConsistencyLight();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_EQUAL(2,(int)f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,(int)f1->getNumberOfTuples());
  f1->applyFunc(1,func1);
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_EQUAL(1,(int)f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,(int)f1->getNumberOfTuples());
  double values1[9]={-1.8,-0.3,1.2,-0.3,1.2,2.7,1.2,2.7,4.2};
  const double *tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+9,values1,values1,std::minus<double>());
  std::transform(values1,values1+9,values1,std::ptr_fun<double,double>(fabs));
  double max=*std::max_element(values1,values1+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest1::testApplyFunc2()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=m->fillFromAnalytic(ON_NODES,2,func2);
  f1->checkConsistencyLight();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_EQUAL(2,(int)f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,(int)f1->getNumberOfTuples());
  //
  MEDCouplingFieldDouble *f2=f1->clone(true);
  CPPUNIT_ASSERT_THROW(f2->applyFunc(1,"a+b+c+d"),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(f2->applyFunc(1,"a/0"),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(f2->applyFunc("a/0"),INTERP_KERNEL::Exception);
  f2->applyFunc("abs(u)^2.4+2*u");
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_EQUAL(2,(int)f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,(int)f1->getNumberOfTuples());
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
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_EQUAL(1,(int)f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,(int)f1->getNumberOfTuples());
  double values1[9]={-1.8,-0.3,1.2,-0.3,1.2,2.7,1.2,2.7,4.2};
  tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+9,values1,values1,std::minus<double>());
  std::transform(values1,values1+9,values1,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values1,values1+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest1::testOperationsOnFields()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=m->fillFromAnalytic(ON_NODES,1,func1);
  MEDCouplingFieldDouble *f2=m->fillFromAnalytic(ON_NODES,1,func1);
  f1->checkConsistencyLight();
  f2->checkConsistencyLight();
  MEDCouplingFieldDouble *f3=(*f1)+(*f2);
  f3->checkConsistencyLight();
  CPPUNIT_ASSERT(f3->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f3->getTimeDiscretization()==ONE_TIME);
  double values1[9]={-1.2,-0.2,0.8,-0.2,0.8,1.8,0.8,1.8,2.8};
  const double *tmp=f3->getArray()->getConstPointer();
  std::transform(tmp,tmp+9,values1,values1,std::minus<double>());
  std::transform(values1,values1+9,values1,std::ptr_fun<double,double>(fabs));
  double max=*std::max_element(values1,values1+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f3->decrRef();
  //
  f3=(*f1)*(*f2);
  f3->checkConsistencyLight();
  CPPUNIT_ASSERT(f3->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f3->getTimeDiscretization()==ONE_TIME);
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
  f4->checkConsistencyLight();
  CPPUNIT_ASSERT(f4->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f4->getTimeDiscretization()==ONE_TIME);
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
  f4->checkConsistencyLight();
  CPPUNIT_ASSERT(f4->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f4->getTimeDiscretization()==ONE_TIME);
  tmp=f4->getArray()->getConstPointer();
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,tmp[i],1.e-12);
  f3->decrRef();
  f4->decrRef();
  //
  f4=f2->buildNewTimeReprFromThis(NO_TIME,false);
  f4->checkConsistencyLight();
  CPPUNIT_ASSERT(f4->getArray()==f2->getArray());
  CPPUNIT_ASSERT(f4->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f4->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_THROW(f3=(*f1)+(*f4),INTERP_KERNEL::Exception);
  MEDCouplingFieldDouble *f5=f4->buildNewTimeReprFromThis(ONE_TIME,false);
  CPPUNIT_ASSERT(f4->getArray()==f5->getArray());
  CPPUNIT_ASSERT(f5->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f5->getTimeDiscretization()==ONE_TIME);
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
  f4=f2->buildNewTimeReprFromThis(NO_TIME,true);
  f4->checkConsistencyLight();
  CPPUNIT_ASSERT(f4->getArray()!=f2->getArray());
  CPPUNIT_ASSERT(f4->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f4->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_THROW(f3=(*f1)+(*f4),INTERP_KERNEL::Exception);
  f5=f4->buildNewTimeReprFromThis(ONE_TIME,true);
  CPPUNIT_ASSERT(f4->getArray()!=f5->getArray());
  CPPUNIT_ASSERT(f2->getArray()!=f5->getArray());
  CPPUNIT_ASSERT(f5->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f5->getTimeDiscretization()==ONE_TIME);
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

void MEDCouplingBasicsTest1::testOperationsOnFields2()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  m->setTime(3.4,5,6); m->setTimeUnit("us");
  int a,b;
  MEDCouplingFieldDouble *f1=m->fillFromAnalytic(ON_NODES,1,"x+y+z");
  MEDCouplingFieldDouble *f2=m->fillFromAnalytic(ON_NODES,1,"a*a+b+c*c");
  MEDCouplingFieldDouble *f3=(*f1)/(*f2);
  f3->checkConsistencyLight();
  CPPUNIT_ASSERT(f3->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f3->getTimeDiscretization()==ONE_TIME);
  const double expected1[9]={-2.4999999999999991, 1.2162162162162162, 0.77868852459016391,
                             0.7407407407407407, 1.129032258064516, 0.81632653061224492,
                             0.86538461538461531, 1.0919540229885056, 0.84302325581395343};
  CPPUNIT_ASSERT_EQUAL(1,(int)f3->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,(int)f3->getNumberOfTuples());
  const double *val=f3->getArray()->getConstPointer();
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],val[i],1.e-12);
  f3->decrRef();
  f1->decrRef();
  f2->decrRef();
  //
  f1=m->buildOrthogonalField();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.4,f1->getTime(a,b),1.e-14);
  CPPUNIT_ASSERT_EQUAL(5,a); CPPUNIT_ASSERT_EQUAL(6,b);
  CPPUNIT_ASSERT_EQUAL(std::string(f1->getTimeUnit()),std::string("us"));
  f2=m->fillFromAnalytic(ON_CELLS,1,"x");
  f3=(*f1)*(*f2);
  const double expected2[15]={-0.035355339059327376,0.,0.035355339059327376, 0.2592724864350674,0.,-0.2592724864350674, 0.37712361663282529,0.,-0.37712361663282529, -0.035355339059327376,0.,0.035355339059327376, 0.31819805153394637,0.,-0.31819805153394637};
  val=f3->getArray()->getConstPointer();
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],val[i],1.e-12);
  f3->decrRef();
  //
  f3=(*f2)*(*f1);
  val=f3->getArray()->getConstPointer();
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],val[i],1.e-12);
  f3->decrRef();
  //
  f1->decrRef();
  f2->decrRef();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest1::testOperationsOnFields3()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  MEDCouplingFieldDouble *f1=m->fillFromAnalytic(ON_NODES,1,"x+y+z");
  MEDCouplingFieldDouble *f2=m->fillFromAnalytic(ON_NODES,1,"a*a+b+c*c");
  (*f1)/=(*f2);
  f1->checkConsistencyLight();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==ONE_TIME);
  const double expected1[9]={-2.4999999999999991, 1.2162162162162162, 0.77868852459016391,
                             0.7407407407407407, 1.129032258064516, 0.81632653061224492,
                             0.86538461538461531, 1.0919540229885056, 0.84302325581395343};
  CPPUNIT_ASSERT_EQUAL(1,(int)f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,(int)f1->getNumberOfTuples());
  const double *val=f1->getArray()->getConstPointer();
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],val[i],1.e-12);
  f1->decrRef();
  f2->decrRef();
  //
  f1=m->buildOrthogonalField();
  f2=m->fillFromAnalytic(ON_CELLS,1,"x");
  (*f1)*=(*f2);
  const double expected2[15]={-0.035355339059327376,0.,0.035355339059327376, 0.2592724864350674,0.,-0.2592724864350674, 0.37712361663282529,0.,-0.37712361663282529, -0.035355339059327376,0.,0.035355339059327376, 0.31819805153394637,0.,-0.31819805153394637};
  val=f1->getArray()->getConstPointer();
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],val[i],1.e-12);
  f1->decrRef();
  //
  f1=m->buildOrthogonalField();
  CPPUNIT_ASSERT_THROW((*f2)*=(*f1),INTERP_KERNEL::Exception);
  f1->decrRef();
  f2->decrRef();
  //
  m->decrRef();
}

/*!
 * Check of LINEAR_TIME and CONST_ON_TIME_INTERVAL policies
 */
void MEDCouplingBasicsTest1::testOperationsOnFields4()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  int nbOfCells=m->getNumberOfCells();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,CONST_ON_TIME_INTERVAL);
  f1->setMesh(m);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(nbOfCells,3);
  f1->setArray(array);
  CPPUNIT_ASSERT_THROW(f1->setEndArray(array),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(f1->getEndArray(),INTERP_KERNEL::Exception);
  array->decrRef();
  double *tmp=array->getPointer();
  const double arr1[15]={0.,10.,20.,1.,11.,21.,2.,12.,22.,3.,13.,23.,4.,14.,24.};
  const double arr2[15]={5.,15.,25.,6.,16.,26.,7.,17.,27.,8.,18.,28.,9.,19.,29.};
  std::copy(arr1,arr1+15,tmp);
  f1->setStartTime(2.,0,0);
  f1->setEndTime(3.,0,0);
  f1->checkConsistencyLight();
  double res[3];
  const double pos[2]={0.3,-0.2};
  f1->getValueOn(pos,res);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[3],res[0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[4],res[1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[5],res[2],1.e-12);
  std::fill(res,res+3,0.);
  f1->getValueOn(pos,2.2,res);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[3],res[0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[4],res[1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[5],res[2],1.e-12);
  std::fill(res,res+3,0.);
  CPPUNIT_ASSERT_THROW(f1->getValueOn(pos,3.2,res),INTERP_KERNEL::Exception);
  MEDCouplingFieldDouble *f2=MEDCouplingFieldDouble::New(ON_CELLS,LINEAR_TIME);
  f2->setMesh(m);
  f2->setArray(f1->getArray());
  f2->setStartTime(2.,3,0);
  f2->setEndTime(4.,13,0);
  CPPUNIT_ASSERT_THROW(f2->checkConsistencyLight(),INTERP_KERNEL::Exception);
  DataArrayDouble *array2=DataArrayDouble::New();
  array2->alloc(nbOfCells,3);
  tmp=array2->getPointer();
  std::copy(arr2,arr2+15,tmp);
  f2->setEndArray(array2);
  array2->decrRef();
  f2->checkConsistencyLight();
  //
  std::fill(res,res+3,0.);
  f2->getValueOn(pos,3.21,res);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.025,res[0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(14.025,res[1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(24.025,res[2],1.e-12);
  MEDCouplingFieldDouble *f3=f2->clone(true);
  CPPUNIT_ASSERT(f2->isEqual(f3,1e-12,1e-12));
  f3->getEndArray()->getPointer()[0]=5.001;
  CPPUNIT_ASSERT(!f2->isEqual(f3,1e-12,1e-12));
  CPPUNIT_ASSERT(f2->isEqual(f3,1e-12,1e-2));
  f3->setStartTime(2.1,3,0);
  CPPUNIT_ASSERT(!f2->isEqual(f3,1e-12,1e-2));
  f3->setStartTime(2.,3,0);
  CPPUNIT_ASSERT(f2->isEqual(f3,1e-12,1e-2));
  f3->setStartTime(2.,4,0);
  CPPUNIT_ASSERT(!f2->isEqual(f3,1e-12,1e-2));
  f3->setStartTime(2.,3,1);
  CPPUNIT_ASSERT(!f2->isEqual(f3,1e-12,1e-2));
  f3->setStartTime(2.,3,0);
  CPPUNIT_ASSERT(f2->isEqual(f3,1e-12,1e-2));
  f3->setEndTime(4.1,13,0);
  CPPUNIT_ASSERT(!f2->isEqual(f3,1e-12,1e-2));
  f3->setEndTime(4.,13,0);
  CPPUNIT_ASSERT(f2->isEqual(f3,1e-12,1e-2));
  f3->setEndTime(4.,14,0);
  CPPUNIT_ASSERT(!f2->isEqual(f3,1e-12,1e-2));
  f3->setEndTime(4.,13,1);
  CPPUNIT_ASSERT(!f2->isEqual(f3,1e-12,1e-2));
  f3->setEndTime(4.,13,0);
  CPPUNIT_ASSERT(f2->isEqual(f3,1e-12,1e-2));
  f3->decrRef();
  MEDCouplingFieldDouble *f4=(*f2)+(*f2);
  std::fill(res,res+3,0.);
  f4->getValueOn(pos,3.21,res);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.05,res[0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(28.05,res[1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(48.05,res[2],1.e-12);
  (*f4)+=*f2;
  std::fill(res,res+3,0.);
  f4->getValueOn(pos,3.21,res);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.075,res[0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42.075,res[1],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(72.075,res[2],1.e-12);
  f4->decrRef();
  //
  f2->decrRef();
  f1->decrRef();
  m->decrRef();
}

bool func4(const double *pt, double *res)
{
  res[0]=pt[0]+pt[1]+pt[2];
  return true;
}

void MEDCouplingBasicsTest1::testMergeNodesOnField()
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

void MEDCouplingBasicsTest1::testCheckConsecutiveCellTypes()
{
  MEDCouplingUMesh *sourceMesh=build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  CPPUNIT_ASSERT(sourceMesh->checkConsecutiveCellTypes());
  const INTERP_KERNEL::NormalizedCellType order1[]={INTERP_KERNEL::NORM_TRI3, INTERP_KERNEL::NORM_QUAD4};
  const INTERP_KERNEL::NormalizedCellType order2[]={INTERP_KERNEL::NORM_QUAD4, INTERP_KERNEL::NORM_TRI3};
  CPPUNIT_ASSERT(!targetMesh->checkConsecutiveCellTypes());
  CPPUNIT_ASSERT(!targetMesh->checkConsecutiveCellTypesAndOrder(order1,order1+2));
  CPPUNIT_ASSERT(!targetMesh->checkConsecutiveCellTypesAndOrder(order2,order2+2));
  DataArrayInt *da=targetMesh->getRenumArrForConsecutiveCellTypesSpec(order1,order1+2);
  CPPUNIT_ASSERT_EQUAL(5,(int)da->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)da->getNumberOfComponents());
  const int expected1[5]={2,0,1,3,4};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+5,da->getConstPointer()));
  da->decrRef();
  da=targetMesh->getRenumArrForConsecutiveCellTypesSpec(order2,order2+2);
  CPPUNIT_ASSERT_EQUAL(5,(int)da->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)da->getNumberOfComponents());
  const int expected2[5]={0,3,4,1,2};
  CPPUNIT_ASSERT(std::equal(expected2,expected2+5,da->getConstPointer()));
  da->decrRef();
  const int renumber1[5]={4,0,1,2,3};
  targetMesh->renumberCells(renumber1,false);
  CPPUNIT_ASSERT(targetMesh->checkConsecutiveCellTypes());
  CPPUNIT_ASSERT(targetMesh->checkConsecutiveCellTypesAndOrder(order1,order1+2));
  CPPUNIT_ASSERT(!targetMesh->checkConsecutiveCellTypesAndOrder(order2,order2+2));
  targetMesh->decrRef();
  sourceMesh->decrRef();
}

void MEDCouplingBasicsTest1::testRearrange2ConsecutiveCellTypes()
{
  MEDCouplingUMesh *m1_1=build2DSourceMesh_1();
  MEDCouplingUMesh *m2_1=build2DTargetMesh_1();
  DataArrayInt *arr1=m1_1->rearrange2ConsecutiveCellTypes();
  MEDCouplingUMesh *m1_2=build2DSourceMesh_1();
  CPPUNIT_ASSERT(m1_2->isEqual(m1_1,1e-12));
  const int expected1[2]={0,1};
  CPPUNIT_ASSERT_EQUAL(2,(int)arr1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)arr1->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+2,arr1->getConstPointer()));
  arr1->decrRef();
  const int expected2[5]={0,3,4,1,2};
  arr1=m2_1->rearrange2ConsecutiveCellTypes();
  CPPUNIT_ASSERT_EQUAL(5,(int)arr1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)arr1->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected2,expected2+5,arr1->getConstPointer()));
  MEDCouplingUMesh *m2_2=build2DTargetMesh_1();
  CPPUNIT_ASSERT_EQUAL(5,(int)arr1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)arr1->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected2,expected2+5,arr1->getConstPointer()));
  CPPUNIT_ASSERT(!m2_2->isEqual(m2_1,1e-12));
  m2_2->renumberCells(expected2,false);
  CPPUNIT_ASSERT(m2_2->isEqual(m2_1,1e-12));
  arr1->decrRef();
  m1_1->decrRef();
  m1_2->decrRef();
  m2_1->decrRef();
  m2_2->decrRef();
}

void MEDCouplingBasicsTest1::testSplitByType()
{
  MEDCouplingUMesh *m1=build3DSurfTargetMesh_1();
  std::vector<MEDCouplingUMesh *> v=m1->splitByType();
  CPPUNIT_ASSERT_EQUAL(3,(int)v.size());
  std::vector<const MEDCouplingUMesh *> v2(v.begin(),v.end());
  MEDCouplingUMesh *m2=MEDCouplingUMesh::MergeUMeshesOnSameCoords(v2);
  m2->setName(m1->getName().c_str());
  CPPUNIT_ASSERT(m1->isEqual(m2,1.e-12));
  for(std::vector<MEDCouplingUMesh *>::const_iterator iter=v.begin();iter!=v.end();iter++)
    (*iter)->decrRef();
  m2->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest1::testFuseUMeshesOnSameCoords()
{
  std::vector<const MEDCouplingUMesh *> meshes;
  MEDCouplingUMesh *m2=build2DTargetMesh_1();
  int cells1[3]={2,3,4};
  MEDCouplingPointSet *m3_1=m2->buildPartOfMySelf(cells1,cells1+3,true);
  MEDCouplingUMesh *m3=dynamic_cast<MEDCouplingUMesh *>(m3_1);
  CPPUNIT_ASSERT(m3);
  meshes.push_back(m3);
  int cells2[3]={1,2,4};
  MEDCouplingPointSet *m4_1=m2->buildPartOfMySelf(cells2,cells2+3,true);
  MEDCouplingUMesh *m4=dynamic_cast<MEDCouplingUMesh *>(m4_1);
  CPPUNIT_ASSERT(m4);
  meshes.push_back(m4);
  int cells3[2]={1,2};
  MEDCouplingPointSet *m5_1=m2->buildPartOfMySelf(cells3,cells3+2,true);
  MEDCouplingUMesh *m5=dynamic_cast<MEDCouplingUMesh *>(m5_1);
  CPPUNIT_ASSERT(m5);
  meshes.push_back(m5);
  m2->decrRef();
  //
  std::vector<DataArrayInt *> corr;
  MEDCouplingUMesh *m7=MEDCouplingUMesh::FuseUMeshesOnSameCoords(meshes,0,corr);
  CPPUNIT_ASSERT_EQUAL(4,(int)m7->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(3,(int)corr.size());
  const int expectedVals1[3]={3,3,2};
  const int expectedVals2[3][3]={{0,1,2},{3,0,2},{3,0,111111}};
  for(int i=0;i<3;i++)
    {
      DataArrayInt *arr=corr[i];
      CPPUNIT_ASSERT_EQUAL(1,(int)arr->getNumberOfComponents());
      int nbOfVals=expectedVals1[i];
      CPPUNIT_ASSERT_EQUAL(nbOfVals,(int)arr->getNumberOfTuples());
      const int *vals=arr->getConstPointer();
      for(int j=0;j<nbOfVals;j++)
        CPPUNIT_ASSERT_EQUAL(expectedVals2[i][j],vals[j]);
    }
  std::vector< std::vector<int> > fidsOfGroups;
  std::vector<const DataArrayInt *> corr2(corr.begin(),corr.end());
  DataArrayInt *arr2=DataArrayInt::MakePartition(corr2,m7->getNumberOfCells(),fidsOfGroups);
  const int fidExp[4]={5,1,3,4};
  const int fidsGrp[3][3]={{1,3,5},{3,4,5},{4,5,23344}};
  CPPUNIT_ASSERT_EQUAL(3,(int)fidsOfGroups.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)arr2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(4,(int)arr2->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(fidExp,fidExp+4,arr2->getConstPointer()));
  for(int i=0;i<3;i++)
    {
      int nbOfVals=expectedVals1[i];
      CPPUNIT_ASSERT_EQUAL(nbOfVals,(int)fidsOfGroups[i].size());
      CPPUNIT_ASSERT(std::equal(fidsOfGroups[i].begin(),fidsOfGroups[i].end(),fidsGrp[i]));
    }
  for(std::vector<DataArrayInt *>::iterator iter=corr.begin();iter!=corr.end();iter++)
    (*iter)->decrRef();
  arr2->decrRef();
  m7->decrRef();
  //
  m3->decrRef();
  m4->decrRef();
  m5->decrRef();
}

void MEDCouplingBasicsTest1::testFuseUMeshesOnSameCoords2()
{
  MEDCouplingUMesh *m2;
  MEDCouplingUMesh *m1=build3DExtrudedUMesh_1(m2);
  m2->decrRef();
  const int part1[5]={2,3,6,4,10};
  MEDCouplingUMesh *m3=(MEDCouplingUMesh *)m1->buildPartOfMySelf(part1,part1+5,true);
  const int part2[4]={5,6,4,7};
  MEDCouplingUMesh *m4=(MEDCouplingUMesh *)m1->buildPartOfMySelf(part2,part2+4,true);
  std::vector<const MEDCouplingUMesh *> meshes;
  meshes.push_back(m1);
  meshes.push_back(m3);
  meshes.push_back(m3);
  meshes.push_back(m4);
  std::vector<DataArrayInt *> corr;
  MEDCouplingUMesh *m5=MEDCouplingUMesh::FuseUMeshesOnSameCoords(meshes,0,corr);
  CPPUNIT_ASSERT_EQUAL(18,(int)m5->getNumberOfCells());
  std::vector<DataArrayInt *>::iterator it=corr.begin();
  const int exp1[4]={18,5,5,4};
  const int exp2[4][18]={
    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17},
    {2,3,6,4,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {2,3,6,4,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {5,6,4,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}
  };
  int i=0;
  for(;it!=corr.end();it++,i++)
    {
      int sz=(*it)->getNumberOfTuples();
      CPPUNIT_ASSERT_EQUAL(exp1[i],sz);
      CPPUNIT_ASSERT(std::equal(exp2[i],exp2[i]+sz,(*it)->getConstPointer()));
    }
  for(it=corr.begin();it!=corr.end();it++)
    (*it)->decrRef();
  m5->decrRef();
  m4->decrRef();
  m3->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest1::testBuildOrthogonalField()
{
  MEDCouplingUMesh *targetMesh=build3DSurfTargetMesh_1();
  MEDCouplingFieldDouble *field=targetMesh->buildOrthogonalField();
  double expected[3]={0.70710678118654746,0.,-0.70710678118654746};
  CPPUNIT_ASSERT_EQUAL(5,(int)field->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,(int)field->getNumberOfComponents());
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
  CPPUNIT_ASSERT_EQUAL(1,(int)field->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,(int)field->getNumberOfComponents());
  vals=field->getArray()->getConstPointer();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.70710678118654746,vals[0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,vals[1],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.70710678118654746,vals[2],1e-12);
  field->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest1::testGetCellsContainingPoint()
{
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_1();
  double pos[12]={0.,0.,0.4,0.4,0.,0.4,0.1,0.1,0.25,0.,0.65,0.};
  MCAuto<DataArrayInt> t1,t2;
  //2D basic
  targetMesh->getCellsContainingPoints(pos,6,1e-12,t1,t2);
  CPPUNIT_ASSERT_EQUAL(6,(int)t1->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(7,(int)t2->getNbOfElems());
  const int expectedValues1[6]={0,4,3,0,1,2};
  const int expectedValues2[7]={0,1,2,3,4,5,6};
  CPPUNIT_ASSERT(std::equal(t1->begin(),t1->end(),expectedValues1));
  CPPUNIT_ASSERT(std::equal(t2->begin(),t2->end(),expectedValues2));
  //2D with no help of bounding box.
  double center[2]={0.2,0.2};
  DataArrayDouble::Rotate2DAlg(center,0.78539816339744830962,6,pos,pos);
  targetMesh->rotate(center,0,0.78539816339744830962);
  targetMesh->getCellsContainingPoints(pos,6,1e-12,t1,t2);
  CPPUNIT_ASSERT_EQUAL(6,(int)t1->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(7,(int)t2->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(t1->begin(),t1->end(),expectedValues1));
  CPPUNIT_ASSERT(std::equal(t2->begin(),t2->end(),expectedValues2));
  //2D outside
  const double pos1bis[2]={-0.3303300858899107,-0.11819805153394641};
  CPPUNIT_ASSERT_EQUAL(-1,targetMesh->getCellContainingPoint(pos1bis,1e-12));
  targetMesh->decrRef();
  //test limits 2D
  targetMesh=build2DTargetMesh_1();
  const double pos2[2]={0.2,-0.05};
  std::vector<int> t11;
  t11.clear();
  targetMesh->getCellsContainingPoint(pos2,1e-12,t11);
  CPPUNIT_ASSERT_EQUAL(2,(int)t11.size());
  const int expectedValues3[2]={0,1};
  CPPUNIT_ASSERT(std::equal(t11.begin(),t11.end(),expectedValues3));
  const double pos3[2]={0.2,0.2};
  t11.clear();
  targetMesh->getCellsContainingPoint(pos3,1e-12,t11);
  CPPUNIT_ASSERT_EQUAL(5,(int)t11.size());
  const int expectedValues4[5]={0,1,2,3,4};
  CPPUNIT_ASSERT(std::equal(t11.begin(),t11.end(),expectedValues4));
  CPPUNIT_ASSERT_EQUAL(0,targetMesh->getCellContainingPoint(pos3,1e-12));
  targetMesh->decrRef();
  //3D
  targetMesh=build3DTargetMesh_1();
  const double pos4[3]={25.,25.,25.};
  CPPUNIT_ASSERT_EQUAL(0,targetMesh->getCellContainingPoint(pos4,1e-12));
  const double pos5[3]={50.,50.,50.};
  t11.clear();
  targetMesh->getCellsContainingPoint(pos5,1e-12,t11);
  CPPUNIT_ASSERT_EQUAL(8,(int)t11.size());
  const int expectedValues5[8]={0,1,2,3,4,5,6,7};
  CPPUNIT_ASSERT(std::equal(t11.begin(),t11.end(),expectedValues5));
  const double pos6[3]={0., 50., 0.};
  t11.clear();
  targetMesh->getCellsContainingPoint(pos6,1e-12,t11);
  CPPUNIT_ASSERT_EQUAL(2,(int)t11.size());
  const int expectedValues6[2]={0,2};
  CPPUNIT_ASSERT(std::equal(t11.begin(),t11.end(),expectedValues6));
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

void MEDCouplingBasicsTest1::testGetValueOn1()
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

void MEDCouplingBasicsTest1::testCMesh0()
{
  MEDCouplingCMesh* mesh=MEDCouplingCMesh::New();
  MEDCouplingCMesh* meshEmpty=mesh->clone(true);
  CPPUNIT_ASSERT(meshEmpty->isEqual(mesh,1e-12));
  
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
  CPPUNIT_ASSERT_EQUAL(1,(int)fieldOnNodes->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(64,(int)fieldOnNodes->getNumberOfTuples());
  const double expected1[64]={-3., -1., 0., 2., -1., 1., 2., 4., 0., 2., 3., 5., 2., 4., 5., 7., -1., 1., 2.,
                              4., 1., 3., 4., 6., 2., 4., 5., 7., 4., 6., 7., 9., 0., 2., 3., 5., 2., 4., 5.,
                              7., 3., 5., 6., 8., 5., 7., 8., 10., 2., 4., 5.,
                              7., 4., 6., 7., 9., 5., 7., 8., 10., 7., 9., 10., 12.};
  const double *val=fieldOnNodes->getArray()->getConstPointer();
  for(int i=0;i<64;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],val[i],1e-12);
  double res[1];  //size fieldOnNodes->getNumberOfComponents()
  fieldOnNodes->getValueOnPos(1,3,2,&res[0]);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,res[0],1e-12);
  fieldOnNodes->decrRef();
  //
  MEDCouplingFieldDouble *fieldOnCells=mesh->fillFromAnalytic(ON_CELLS,1,"x+y/2.+z/3.");
  CPPUNIT_ASSERT_EQUAL(1,(int)fieldOnCells->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(27,(int)fieldOnCells->getNumberOfTuples());
  val=fieldOnCells->getArray()->getConstPointer();
  const double expected2[27]={0, 1.5, 3, 1.5, 3, 4.5, 3, 4.5, 6, 1.5, 3, 4.5, 3, 4.5,
                              6, 4.5, 6, 7.5, 3, 4.5, 6, 4.5, 6, 7.5, 6, 7.5, 9};
  for(int i=0;i<27;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],val[i],1e-12);
  fieldOnCells->getValueOnPos(1,2,1,&res[0]);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(6.,res[0],1e-12);
  fieldOnCells->decrRef();
  //
  MEDCouplingMesh* meshDeepCopy=mesh->deepCopy();
  MEDCouplingCMesh* meshClone=mesh->clone(false);
  
  CPPUNIT_ASSERT_THROW(meshEmpty->copyTinyStringsFrom(0),INTERP_KERNEL::Exception);
  meshEmpty->copyTinyStringsFrom(mesh);
  //no data in meshEmpty, expected false
  CPPUNIT_ASSERT(!meshEmpty->isEqual(mesh,1e-12));
  
  CPPUNIT_ASSERT(meshDeepCopy->isEqual(mesh,1e-12));
  meshDeepCopy->copyTinyStringsFrom(mesh);
  CPPUNIT_ASSERT(meshDeepCopy->isEqual(mesh,1e-12));
  CPPUNIT_ASSERT(meshClone->isEqual(mesh,1e-12));
  
  CPPUNIT_ASSERT_EQUAL(CARTESIAN,mesh->getType());
  CPPUNIT_ASSERT_EQUAL(CARTESIAN,meshEmpty->getType());
  CPPUNIT_ASSERT_EQUAL(CARTESIAN,meshDeepCopy->getType());
  CPPUNIT_ASSERT_EQUAL(CARTESIAN,meshClone->getType());
  
  mesh->decrRef();
  meshEmpty->decrRef();
  meshDeepCopy->decrRef();
  meshClone->decrRef();
}

void MEDCouplingBasicsTest1::testCMesh1()
{
  MEDCouplingCMesh *mesh1,*mesh2,*mesh3;
  mesh1=MEDCouplingCMesh::New();
  DataArrayDouble* coordsX1=DataArrayDouble::New();
  double arrX1[4] = { -1., 1., 2., 4. };
  coordsX1->useArray(arrX1,false, CPP_DEALLOC,4,1);
  DataArrayDouble* coordsY1=DataArrayDouble::New();
  double arrY1[4] = { -2., 2., 4., 8. };
  coordsY1->useArray(arrY1,false, CPP_DEALLOC,4,1);
  DataArrayDouble* coordsZ1=DataArrayDouble::New();
  double arrZ1[4] = { -3., 3., 6., 12. };
  coordsZ1->useArray(arrZ1,false, CPP_DEALLOC,4,1);
  mesh1->setCoords(coordsX1,coordsY1,coordsZ1);
  
  mesh2=MEDCouplingCMesh::New();
  DataArrayDouble* coordsX2=DataArrayDouble::New();
  double arrX2[4] = { -1., 1., 2., 4. };
  coordsX2->useArray(arrX2,false, CPP_DEALLOC,4,1);
  DataArrayDouble* coordsY2=DataArrayDouble::New();
  double arrY2[4] = { -2., 2., 4., 8. };
  coordsY2->useArray(arrY2,false, CPP_DEALLOC,4,1);
  DataArrayDouble* coordsZ2=DataArrayDouble::New();
  double arrZ2[4] = { -3., 3., 6., 12.+1e-6 };   //here is not equal
  coordsZ2->useArray(arrZ2,false, CPP_DEALLOC,4,1);
  mesh2->setCoords(coordsX2,coordsY2,coordsZ2);
  
  mesh3=MEDCouplingCMesh::New();
  DataArrayDouble* coordsX3=DataArrayDouble::New();
  double arrX3[1] = { -1.};
  coordsX3->useArray(arrX3,false, CPP_DEALLOC,1,1);
  DataArrayDouble* coordsY3=DataArrayDouble::New();
  double arrY3[1] = { -2.};
  coordsY3->useArray(arrY3,false, CPP_DEALLOC,1,1);
  DataArrayDouble* coordsZ3=DataArrayDouble::New();
  double arrZ3[1] = { -3.};
  coordsZ3->useArray(arrZ3,false, CPP_DEALLOC,1,1);
  mesh3->setCoords(coordsX3,coordsY3,coordsZ3);
  
  CPPUNIT_ASSERT_EQUAL(3,mesh1->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,mesh1->getMeshDimension());
  
  CPPUNIT_ASSERT(!mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(!mesh2->isEqual(mesh1,1e-12));
  CPPUNIT_ASSERT(!mesh2->isEqualWithoutConsideringStr(mesh1,1e-12));
  CPPUNIT_ASSERT(mesh1->isEqual(mesh2,1e-5));
  CPPUNIT_ASSERT(!mesh1->isEqual(mesh2,1e-7));
  
  CPPUNIT_ASSERT_THROW(mesh3->checkConsistency(1e-12),INTERP_KERNEL::Exception);
  mesh1->checkConsistency(1e-12);
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_HEXA8,mesh1->getTypeOfCell(1));
  
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_HEXA8,*((mesh1->getAllGeoTypes()).begin()));
  CPPUNIT_ASSERT_EQUAL(27,(int)mesh1->getNumberOfCellsWithType(INTERP_KERNEL::NORM_HEXA8));
  CPPUNIT_ASSERT_THROW(mesh1->getNumberOfCellsWithType(INTERP_KERNEL::NORM_QUAD4),INTERP_KERNEL::Exception);
  
  std::vector<double> coo;
  mesh1->getCoordinatesOfNode(0, coo);
  CPPUNIT_ASSERT_EQUAL(3,(int) coo.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.,coo[0],1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-2.,coo[1],1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.,coo[2],1e-14);
  coo.clear();
  mesh1->getCoordinatesOfNode(63, coo);
  CPPUNIT_ASSERT_EQUAL(3,(int) coo.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.,coo[0],1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.,coo[1],1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(12.,coo[2],1e-14);
  
  std::string repr;
  repr=mesh1->simpleRepr();
  repr=mesh1->advancedRepr();
  CPPUNIT_ASSERT(!(repr.find("Cartesian")==std::string::npos));
  CPPUNIT_ASSERT(!(repr.find("Number of components : 1")==std::string::npos));
  CPPUNIT_ASSERT(!(repr.find("Number of tuples : 4")==std::string::npos));
  CPPUNIT_ASSERT(!(repr.find("Z Array :")==std::string::npos));
  coordsX1->decrRef();
  coordsY1->decrRef();
  coordsZ1->decrRef();
  coordsX2->decrRef();
  coordsY2->decrRef();
  coordsZ2->decrRef();
  coordsX3->decrRef();
  coordsY3->decrRef();
  coordsZ3->decrRef();
  mesh1->decrRef();
  mesh2->decrRef();
  mesh3->decrRef();
}
  
void MEDCouplingBasicsTest1::testCMesh2()
{
  MEDCouplingCMesh *mesh1;
  mesh1=MEDCouplingCMesh::New();
  DataArrayDouble* coordsX1=DataArrayDouble::New();
  double arrX1[4] = { -1., 1., 2., 4. };
  coordsX1->useArray(arrX1,false, CPP_DEALLOC,4,1);
  DataArrayDouble* coordsY1=DataArrayDouble::New();
  double arrY1[4] = { -2., 2., 4., 8. };
  coordsY1->useArray(arrY1,false, CPP_DEALLOC,4,1);
  DataArrayDouble* coordsZ1=DataArrayDouble::New();
  double arrZ1[4] = { -3., 3., 6., 12. };
  coordsZ1->useArray(arrZ1,false, CPP_DEALLOC,4,1);
  mesh1->setCoords(coordsX1,coordsY1,coordsZ1);
  
  std::vector<int> dis=mesh1->getDistributionOfTypes();
  CPPUNIT_ASSERT_EQUAL(3,(int) dis.size());
  CPPUNIT_ASSERT_EQUAL((int) INTERP_KERNEL::NORM_HEXA8,dis[0]);
  CPPUNIT_ASSERT_EQUAL(27,dis[1]);
  CPPUNIT_ASSERT_EQUAL(-1,dis[2]);
  
  std::vector<const DataArrayInt *> idsPerType;
  CPPUNIT_ASSERT(!(mesh1->checkTypeConsistencyAndContig(dis, idsPerType)));
  dis[0]=(int) INTERP_KERNEL::NORM_QUAD4;
  CPPUNIT_ASSERT_THROW(mesh1->checkTypeConsistencyAndContig(dis, idsPerType),INTERP_KERNEL::Exception);
  
  dis[0]=(int) INTERP_KERNEL::NORM_HEXA8;
  dis[2]=0;
  DataArrayInt *ids=DataArrayInt::New();
  ids->alloc(10,1);
  ids->fillWithValue(23);
  idsPerType.push_back(ids);
  DataArrayInt* check=mesh1->checkTypeConsistencyAndContig(dis, idsPerType);
  CPPUNIT_ASSERT(check);
  CPPUNIT_ASSERT(check->isEqual(*ids));
  
  std::vector<int> code;
  std::vector<DataArrayInt *> idsInPflPerType;
  std::vector<DataArrayInt *> pfls;
  mesh1->splitProfilePerType(ids,code,idsInPflPerType,pfls);
  CPPUNIT_ASSERT_EQUAL(3,(int)code.size());
  CPPUNIT_ASSERT_EQUAL((int) INTERP_KERNEL::NORM_HEXA8,code[0]);
  CPPUNIT_ASSERT_EQUAL(10,code[1]);
  CPPUNIT_ASSERT_EQUAL(0,code[2]);
  CPPUNIT_ASSERT_EQUAL(1,(int)idsInPflPerType.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)pfls.size());
  DataArrayInt *exp=DataArrayInt::New(); exp->alloc(10,1); exp->iota(0);
  CPPUNIT_ASSERT(idsInPflPerType[0]->isEqual(*exp));
  exp->decrRef();
  CPPUNIT_ASSERT(pfls[0]->isEqual(*ids));
  idsInPflPerType[0]->decrRef();
  pfls[0]->decrRef();

  ids->decrRef();
  check->decrRef();
  int cells1[4]={0,1,25,26};
  MEDCouplingUMesh *partMesh1=
    dynamic_cast<MEDCouplingUMesh *>(mesh1->buildPart(cells1,cells1+4));
  CPPUNIT_ASSERT(partMesh1);
  CPPUNIT_ASSERT_EQUAL(4,(int)partMesh1->getNumberOfCellsWithType(INTERP_KERNEL::NORM_HEXA8));
  CPPUNIT_ASSERT_EQUAL(64,mesh1->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(64,partMesh1->getNumberOfNodes());
  
  int cells2[2]={25,26};
  DataArrayInt* arr1;
  MEDCouplingCMesh *partMesh2=
    dynamic_cast<MEDCouplingCMesh *>(mesh1->buildPartAndReduceNodes(cells2,cells2+2,arr1));
  CPPUNIT_ASSERT(partMesh2);
  CPPUNIT_ASSERT_EQUAL(2,(int)partMesh2->getNumberOfCellsWithType(INTERP_KERNEL::NORM_HEXA8));
  CPPUNIT_ASSERT_EQUAL(12,partMesh2->getNumberOfNodes());
  
  int cells3[2]={2,3};
  DataArrayInt* arr2;
  MEDCouplingUMesh *partMesh3=
    dynamic_cast<MEDCouplingUMesh *>(partMesh1->buildPartAndReduceNodes(cells3,cells3+2,arr2));
  CPPUNIT_ASSERT(partMesh3);
  CPPUNIT_ASSERT_EQUAL(2,(int)partMesh3->getNumberOfCellsWithType(INTERP_KERNEL::NORM_HEXA8));
  CPPUNIT_ASSERT_EQUAL(12,partMesh3->getNumberOfNodes());
  
  CPPUNIT_ASSERT_THROW(mesh1->simplexize(0),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(mesh1->getMeasureFieldOnNode(true),INTERP_KERNEL::Exception);

  double bbox1[6];
  double bbox2[6];
  mesh1->getBoundingBox(bbox1);
  partMesh1->getBoundingBox(bbox2);
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(bbox1[i],bbox2[i],1e-12);
  partMesh3->getBoundingBox(bbox1);
  partMesh2->getBoundingBox(bbox2);
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(bbox1[i],bbox2[i],1e-12);
  
  CPPUNIT_ASSERT_THROW(mesh1->buildOrthogonalField(),INTERP_KERNEL::Exception);
  MEDCouplingCMesh *mesh2d=MEDCouplingCMesh::New();
  mesh2d->setCoords(coordsX1,coordsY1);
  MEDCouplingFieldDouble *f1=mesh2d->buildOrthogonalField();
  
  std::vector<double> tinyInfoD;
  std::vector<int> tinyInfo;
  std::vector<std::string> littleStrings;
  mesh2d->getTinySerializationInformation(tinyInfoD, tinyInfo, littleStrings);
  CPPUNIT_ASSERT_EQUAL(5,(int)tinyInfo.size());
  CPPUNIT_ASSERT_EQUAL(4,(int)tinyInfo[0]);   //x
  CPPUNIT_ASSERT_EQUAL(4,(int)tinyInfo[1]);   //y
  CPPUNIT_ASSERT_EQUAL(-1,(int)tinyInfo[2]);  //z
  CPPUNIT_ASSERT_EQUAL(-1,(int)tinyInfo[3]);  //it
  CPPUNIT_ASSERT_EQUAL(-1,(int)tinyInfo[4]);   //order
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,tinyInfoD[0],1e-14); //time
  DataArrayInt* d1=DataArrayInt::New();
  DataArrayDouble* d2=DataArrayDouble::New();
  mesh2d->resizeForUnserialization(tinyInfo, d1, d2, littleStrings);
  CPPUNIT_ASSERT_EQUAL(0,(int)d1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(8,(int)d2->getNumberOfTuples());
 
  partMesh1->decrRef();
  partMesh2->decrRef();
  partMesh3->decrRef();
  mesh2d->decrRef();
  arr1->decrRef();
  arr2->decrRef();
  f1->decrRef();
  d1->decrRef();
  d2->decrRef();
  coordsX1->decrRef();
  coordsY1->decrRef();
  coordsZ1->decrRef();
  mesh1->decrRef();
}

void MEDCouplingBasicsTest1::testScale()
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

void MEDCouplingBasicsTest1::testTryToShareSameCoords()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  MEDCouplingUMesh *m2=build2DTargetMesh_1();
  CPPUNIT_ASSERT(m1->getCoords()!=m2->getCoords());
  m1->tryToShareSameCoords(*m2,1e-12);
  CPPUNIT_ASSERT(m1->getCoords()==m2->getCoords());
  m1->tryToShareSameCoords(*m2,1e-12);
  CPPUNIT_ASSERT(m1->getCoords()==m2->getCoords());
  m2->tryToShareSameCoords(*m1,1e-12);
  CPPUNIT_ASSERT(m1->getCoords()==m2->getCoords());
  m1->decrRef();
  m2->decrRef();
  //
  m1=build2DTargetMesh_1();
  m2=build2DTargetMesh_2();
  CPPUNIT_ASSERT(m1->getCoords()!=m2->getCoords());
  m1->tryToShareSameCoords(*m2,1e-12);
  CPPUNIT_ASSERT(m1->getCoords()==m2->getCoords());
  m1->tryToShareSameCoords(*m2,1e-12);
  CPPUNIT_ASSERT(m1->getCoords()==m2->getCoords());
  m2->tryToShareSameCoords(*m1,1e-12);
  CPPUNIT_ASSERT(m1->getCoords()==m2->getCoords());
  m1->decrRef();
  m2->decrRef();
  //
  m1=build2DTargetMesh_1();
  m2=build2DSourceMesh_1();
  CPPUNIT_ASSERT(m1->getCoords()!=m2->getCoords());
  CPPUNIT_ASSERT_THROW(m1->tryToShareSameCoords(*m2,1e-12),INTERP_KERNEL::Exception);
  m1->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest1::testFindNodeOnPlane()
{
  MEDCouplingUMesh *mesh=build3DTargetMesh_1();
  std::vector<int> n;
  double pt[3]={300.,300.,0.};
  double v[3]={0.,0.,2.};
  mesh->findNodesOnPlane(pt,v,1e-12,n);
  CPPUNIT_ASSERT_EQUAL(9,(int)n.size());
  MEDCouplingUMesh *m3dSurf=(MEDCouplingUMesh *)mesh->buildFacePartOfMySelfNode(&n[0],&n[0]+n.size(),true);
  MEDCouplingMappedExtrudedMesh *me=MEDCouplingMappedExtrudedMesh::New(mesh,m3dSurf,0);
  const DataArrayInt *da=me->getMesh3DIds();
  CPPUNIT_ASSERT_EQUAL(8,(int)me->getNumberOfCells());
  const int expected[8]={0,1,2,3,4,5,6,7};
  const int *val=da->getConstPointer();
  for(int i=0;i<8;i++)
    CPPUNIT_ASSERT_EQUAL(expected[i],val[i]);
  me->decrRef();
  m3dSurf->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest1::testRenumberCells()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  MEDCouplingUMesh *m2=build3DSurfTargetMesh_1();
  CPPUNIT_ASSERT(m->isEqual(m2,0));
  const int arr[5]={12,3,25,2,26};
  m->renumberCells(arr,true);
  CPPUNIT_ASSERT(!m->isEqual(m2,0));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,m->getTypeOfCell(0));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,m->getTypeOfCell(1));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,m->getTypeOfCell(2));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,m->getTypeOfCell(3));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_QUAD4,m->getTypeOfCell(4));
  const int arr2[5]={5,-1,-5,4,8};
  m->renumberCells(arr2,true);
  CPPUNIT_ASSERT(m->isEqual(m2,0));
  m->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest1::testChangeSpaceDimension()
{
  MEDCouplingUMesh *m1=build3DSurfTargetMesh_1();
  MEDCouplingUMesh *m2=build2DTargetMesh_1();
  //
  CPPUNIT_ASSERT_EQUAL(3,m1->getSpaceDimension());
  m1->changeSpaceDimension(2);
  CPPUNIT_ASSERT_EQUAL(2,m1->getSpaceDimension());
  m1->setName(m2->getName().c_str());
  CPPUNIT_ASSERT(m1->isEqual(m2,1e-12));
  m1->changeSpaceDimension(3);
  CPPUNIT_ASSERT_EQUAL(3,m1->getSpaceDimension());
  const double expected[27]={-0.3,-0.3,0., 0.2,-0.3,0., 0.7,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.7,0.2,0., -0.3,0.7,0., 0.2,0.7,0., 0.7,0.7,0.};
  const double *val=m1->getCoords()->getConstPointer();
  for(int i=0;i<27;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[i],val[i],1e-14);
  //
  m1->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest1::testSetConnectivity()
{
  MEDCouplingUMesh *m1 = build1DTargetMesh_1();

  DataArrayInt * conn = DataArrayInt::New();
  DataArrayInt * connI = DataArrayInt::New();
  m1->setConnectivity(conn, connI, true); // was SEG-Faulting with empty arrays
  conn->decrRef();
  connI->decrRef();
  m1->decrRef();
}
