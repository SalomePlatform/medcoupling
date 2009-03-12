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

void MEDCouplingBasicsTest::testBuildPartOfMySelf()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  mesh->setName("Toto");
  const int tab1[2]={0,4};
  const int tab2[3]={0,2,3};
  //
  MEDCouplingUMesh *subMesh=mesh->buildPartOfMySelf(tab1,tab1+2,true);
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
  subMesh=mesh->buildPartOfMySelf(tab2,tab2+3,true);
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
  MEDCouplingUMesh *subMesh=mesh->buildPartOfMySelf(tab1,tab1+2,true);
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
  subMesh=mesh->buildPartOfMySelf(tab1,tab1+2,false);
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
  
  //
  mesh1->decrRef();
  mesh2->decrRef();
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
