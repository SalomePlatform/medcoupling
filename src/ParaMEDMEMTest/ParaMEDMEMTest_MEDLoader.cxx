// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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

#include "ParaMEDMEMTest.hxx"
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include <cppunit/TestAssert.h>

#include <algorithm>
#include <numeric>
#include <iostream>
#include <iterator>

using namespace std;
using namespace INTERP_KERNEL;
using namespace ParaMEDMEM;

void ParaMEDMEMTest::testMEDLoaderRead1()
{
  string fileName=getResourceFile("pointe.med");
  vector<string> meshNames=MEDLoader::GetMeshNames(fileName.c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)meshNames.size());
  MEDCouplingUMesh *mesh=MEDLoader::ReadUMeshFromFile(fileName.c_str(),meshNames[0].c_str(),0);
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(16,mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,(int)mesh->getAllGeoTypes().size());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,mesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,mesh->getTypeOfCell(12));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,mesh->getTypeOfCell(13));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,mesh->getTypeOfCell(14));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,mesh->getTypeOfCell(15));
  CPPUNIT_ASSERT_EQUAL((std::size_t)90,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(701,std::accumulate(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+90,0));
  CPPUNIT_ASSERT_EQUAL(705,std::accumulate(mesh->getNodalConnectivityIndex()->getPointer(),mesh->getNodalConnectivityIndex()->getPointer()+17,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(46.,std::accumulate(mesh->getCoords()->getPointer(),mesh->getCoords()->getPointer()+57,0),1e-12);
  mesh->decrRef();
  //
  vector<string> families=MEDLoader::GetMeshFamiliesNames(fileName.c_str(),meshNames[0].c_str());
  CPPUNIT_ASSERT_EQUAL(8,(int)families.size());
  CPPUNIT_ASSERT(families[2]=="FAMILLE_ELEMENT_3");
  //
  vector<string> families2;
  families2.push_back(families[2]);
  mesh=MEDLoader::ReadUMeshFromFamilies(fileName.c_str(),meshNames[0].c_str(),0,families2);
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,(int)mesh->getAllGeoTypes().size());
  CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,mesh->getTypeOfCell(0));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,mesh->getTypeOfCell(1));
  CPPUNIT_ASSERT_EQUAL((std::size_t)11,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(132,std::accumulate(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+11,0));
  CPPUNIT_ASSERT_EQUAL(16,std::accumulate(mesh->getNodalConnectivityIndex()->getPointer(),mesh->getNodalConnectivityIndex()->getPointer()+3,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(46.,std::accumulate(mesh->getCoords()->getPointer(),mesh->getCoords()->getPointer()+57,0),1e-12);
  mesh->decrRef();
  //
  vector<string> groups=MEDLoader::GetMeshGroupsNames(fileName.c_str(),meshNames[0].c_str());
  CPPUNIT_ASSERT_EQUAL(5,(int)groups.size());
  CPPUNIT_ASSERT(groups[0]=="groupe1");
  CPPUNIT_ASSERT(groups[1]=="groupe2");
  CPPUNIT_ASSERT(groups[2]=="groupe3");
  CPPUNIT_ASSERT(groups[3]=="groupe4");
  CPPUNIT_ASSERT(groups[4]=="groupe5");
  vector<string> groups2;
  groups2.push_back(groups[0]);
  mesh=MEDLoader::ReadUMeshFromGroups(fileName.c_str(),meshNames[0].c_str(),0,groups2);
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(7,mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,(int)mesh->getAllGeoTypes().size());
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,mesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,mesh->getTypeOfCell(6));
  CPPUNIT_ASSERT_EQUAL((std::size_t)36,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(254,std::accumulate(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+36,0));
  CPPUNIT_ASSERT_EQUAL(141,std::accumulate(mesh->getNodalConnectivityIndex()->getPointer(),mesh->getNodalConnectivityIndex()->getPointer()+8,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(46.,std::accumulate(mesh->getCoords()->getPointer(),mesh->getCoords()->getPointer()+57,0),1e-12);
  mesh->decrRef();
  //
  std::vector<std::string> fieldsName=MEDLoader::GetCellFieldNamesOnMesh(fileName.c_str(),meshNames[0].c_str());
  CPPUNIT_ASSERT_EQUAL(2,(int)fieldsName.size());
  CPPUNIT_ASSERT(fieldsName[0]=="fieldcelldoublescalar");
  CPPUNIT_ASSERT(fieldsName[1]=="fieldcelldoublevector");
  std::vector<std::pair<int,int> > its0=MEDLoader::GetCellFieldIterations(fileName.c_str(),meshNames[0].c_str(),fieldsName[0].c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)its0.size());
  CPPUNIT_ASSERT_EQUAL(-1,its0[0].first);
  CPPUNIT_ASSERT_EQUAL(-1,its0[0].second);
  std::vector<std::pair<int,int> > its1=MEDLoader::GetCellFieldIterations(fileName.c_str(),meshNames[0].c_str(),fieldsName[1].c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)its1.size());
  CPPUNIT_ASSERT_EQUAL(-1,its1[0].first);
  CPPUNIT_ASSERT_EQUAL(-1,its1[0].second);
  //
  MEDCouplingFieldDouble *field0=MEDLoader::ReadFieldCell(fileName.c_str(),meshNames[0].c_str(),0,fieldsName[0].c_str(),its0[0].first,its0[0].second);
  field0->checkCoherency();
  CPPUNIT_ASSERT(field0->getName()==fieldsName[0]);
  CPPUNIT_ASSERT_EQUAL(1,field0->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(16,field0->getNumberOfTuples());
  const double expectedValues[16]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,2.,3.,3.,2.};
  double diffValue[16];
  std::transform(field0->getArray()->getPointer(),field0->getArray()->getPointer()+16,expectedValues,diffValue,std::minus<double>());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::max_element(diffValue,diffValue+16),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::min_element(diffValue,diffValue+16),1e-12);
  const MEDCouplingUMesh *constMesh=dynamic_cast<const MEDCouplingUMesh *>(field0->getMesh());
  CPPUNIT_ASSERT(constMesh);
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(16,constMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,constMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,(int)constMesh->getAllGeoTypes().size());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,constMesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(12));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(13));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(14));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(15));
  CPPUNIT_ASSERT_EQUAL((std::size_t)90,constMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(701,std::accumulate(constMesh->getNodalConnectivity()->getConstPointer(),constMesh->getNodalConnectivity()->getConstPointer()+90,0));
  CPPUNIT_ASSERT_EQUAL(705,std::accumulate(constMesh->getNodalConnectivityIndex()->getConstPointer(),constMesh->getNodalConnectivityIndex()->getConstPointer()+17,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(46.,std::accumulate(constMesh->getCoords()->getConstPointer(),constMesh->getCoords()->getConstPointer()+57,0),1e-12);
  field0->decrRef();
  //
  MEDCouplingFieldDouble *field1=MEDLoader::ReadFieldCell(fileName.c_str(),meshNames[0].c_str(),0,fieldsName[1].c_str(),its1[0].first,its1[0].second);
  field1->checkCoherency();
  CPPUNIT_ASSERT(field1->getName()==fieldsName[1]);
  CPPUNIT_ASSERT_EQUAL(3,field1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(16,field1->getNumberOfTuples());
  const double expectedValues2[48]={1.,0.,1.,1.,0.,1.,1.,0.,1.,2.,1.,0.,2.,1.,0.,2.,1.,0.,3.,0.,1.,3.,0.,1.,3.,0.,1.,4.,1.,0.,4.,1.,0.,4.,1.,0.,5.,0.,0.,6.,1.,1.,6.,0.,0.,5.,1.,1.};
  double diffValue2[48];
  std::transform(field1->getArray()->getPointer(),field1->getArray()->getPointer()+48,expectedValues2,diffValue2,std::minus<double>());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::max_element(diffValue2,diffValue2+48),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::min_element(diffValue2,diffValue2+48),1e-12);
  constMesh=dynamic_cast<const MEDCouplingUMesh *>(field1->getMesh());
  CPPUNIT_ASSERT(constMesh);
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(16,constMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,constMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,(int)constMesh->getAllGeoTypes().size());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,constMesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(12));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(13));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(14));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(15));
  CPPUNIT_ASSERT_EQUAL((std::size_t)90,constMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(701,std::accumulate(constMesh->getNodalConnectivity()->getConstPointer(),constMesh->getNodalConnectivity()->getConstPointer()+90,0));
  CPPUNIT_ASSERT_EQUAL(705,std::accumulate(constMesh->getNodalConnectivityIndex()->getConstPointer(),constMesh->getNodalConnectivityIndex()->getConstPointer()+17,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(46.,std::accumulate(constMesh->getCoords()->getConstPointer(),constMesh->getCoords()->getConstPointer()+57,0),1e-12);
  field1->decrRef();
  //fields on nodes
  std::vector<std::string> fieldsNameNode=MEDLoader::GetNodeFieldNamesOnMesh(fileName.c_str(),meshNames[0].c_str());
  CPPUNIT_ASSERT_EQUAL(2,(int)fieldsNameNode.size());
  CPPUNIT_ASSERT(fieldsNameNode[0]=="fieldnodedouble");
  CPPUNIT_ASSERT(fieldsNameNode[1]=="fieldnodeint");
  std::vector<std::pair<int,int> > its0Node=MEDLoader::GetNodeFieldIterations(fileName.c_str(),meshNames[0].c_str(),fieldsNameNode[0].c_str());
  CPPUNIT_ASSERT_EQUAL(3,(int)its0Node.size());
  CPPUNIT_ASSERT_EQUAL(-1,its0Node[0].first);
  CPPUNIT_ASSERT_EQUAL(-1,its0Node[0].second);
  CPPUNIT_ASSERT_EQUAL(1,its0Node[1].first);
  CPPUNIT_ASSERT_EQUAL(-1,its0Node[1].second);
  CPPUNIT_ASSERT_EQUAL(2,its0Node[2].first);
  CPPUNIT_ASSERT_EQUAL(-1,its0Node[2].second);
  MEDCouplingFieldDouble *field0Nodes=MEDLoader::ReadFieldNode(fileName.c_str(),meshNames[0].c_str(),0,fieldsNameNode[0].c_str(),its0Node[0].first,its0Node[0].second);
  field0Nodes->checkCoherency();
  CPPUNIT_ASSERT(field0Nodes->getName()==fieldsNameNode[0]);
  CPPUNIT_ASSERT_EQUAL(1,field0Nodes->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(19,field0Nodes->getNumberOfTuples());
  const double expectedValues3[19]={1.,1.,1.,2.,2.,2.,3.,3.,3.,4.,4.,4.,5.,5.,5.,6.,6.,6.,7.};
  double diffValue3[19];
  std::transform(field0Nodes->getArray()->getPointer(),field0Nodes->getArray()->getPointer()+19,expectedValues3,diffValue3,std::minus<double>());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::max_element(diffValue3,diffValue3+19),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::min_element(diffValue3,diffValue3+19),1e-12);
  constMesh=dynamic_cast<const MEDCouplingUMesh *>(field0Nodes->getMesh());
  CPPUNIT_ASSERT(constMesh);
  field0Nodes->decrRef();
  //
  field0Nodes=MEDLoader::ReadFieldNode(fileName.c_str(),meshNames[0].c_str(),0,fieldsNameNode[0].c_str(),its0Node[2].first,its0Node[2].second);
  field0Nodes->checkCoherency();
  CPPUNIT_ASSERT(field0Nodes->getName()==fieldsNameNode[0]);
  CPPUNIT_ASSERT_EQUAL(1,field0Nodes->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(19,field0Nodes->getNumberOfTuples());
  const double expectedValues4[19]={1.,2.,2.,2.,3.,3.,3.,4.,4.,4.,5.,5.,5.,6.,6.,6.,7.,7.,7.};
  std::transform(field0Nodes->getArray()->getPointer(),field0Nodes->getArray()->getPointer()+19,expectedValues4,diffValue3,std::minus<double>());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::max_element(diffValue3,diffValue3+19),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::min_element(diffValue3,diffValue3+19),1e-12);
  constMesh=dynamic_cast<const MEDCouplingUMesh *>(field0Nodes->getMesh());
  CPPUNIT_ASSERT(constMesh);
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(16,constMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,constMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,(int)constMesh->getAllGeoTypes().size());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,constMesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(12));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(13));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(14));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(15));
  CPPUNIT_ASSERT_EQUAL((std::size_t)90,constMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(701,std::accumulate(constMesh->getNodalConnectivity()->getConstPointer(),constMesh->getNodalConnectivity()->getConstPointer()+90,0));
  CPPUNIT_ASSERT_EQUAL(705,std::accumulate(constMesh->getNodalConnectivityIndex()->getConstPointer(),constMesh->getNodalConnectivityIndex()->getConstPointer()+17,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(46.,std::accumulate(constMesh->getCoords()->getConstPointer(),constMesh->getCoords()->getConstPointer()+57,0),1e-12);
  field0Nodes->decrRef();
  //
  field0Nodes=MEDLoader::ReadFieldNode(fileName.c_str(),meshNames[0].c_str(),0,fieldsNameNode[0].c_str(),its0Node[0].first,its0Node[0].second);
  field0Nodes->checkCoherency();
  CPPUNIT_ASSERT(field0Nodes->getName()==fieldsNameNode[0]);
  CPPUNIT_ASSERT_EQUAL(1,field0Nodes->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(19,field0Nodes->getNumberOfTuples());
  const double expectedValues5[19]={1.,1.,1.,2.,2.,2.,3.,3.,3.,4.,4.,4.,5.,5.,5.,6.,6.,6.,7.};
  std::transform(field0Nodes->getArray()->getPointer(),field0Nodes->getArray()->getPointer()+19,expectedValues5,diffValue3,std::minus<double>());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::max_element(diffValue3,diffValue3+19),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::min_element(diffValue3,diffValue3+19),1e-12);
  constMesh=dynamic_cast<const MEDCouplingUMesh *>(field0Nodes->getMesh());
  CPPUNIT_ASSERT(constMesh);
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(16,constMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,constMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,(int)constMesh->getAllGeoTypes().size());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,constMesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(12));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(13));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(14));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(15));
  CPPUNIT_ASSERT_EQUAL((std::size_t)90,constMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(701,std::accumulate(constMesh->getNodalConnectivity()->getConstPointer(),constMesh->getNodalConnectivity()->getConstPointer()+90,0));
  CPPUNIT_ASSERT_EQUAL(705,std::accumulate(constMesh->getNodalConnectivityIndex()->getConstPointer(),constMesh->getNodalConnectivityIndex()->getConstPointer()+17,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(46.,std::accumulate(constMesh->getCoords()->getConstPointer(),constMesh->getCoords()->getConstPointer()+57,0),1e-12);
  field0Nodes->decrRef();
}

void ParaMEDMEMTest::testMEDLoaderPolygonRead()
{
  string fileName=getResourceFile("polygones.med");
  vector<string> meshNames=MEDLoader::GetMeshNames(fileName.c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)meshNames.size());
  CPPUNIT_ASSERT(meshNames[0]=="Bord");
  MEDCouplingUMesh *mesh=MEDLoader::ReadUMeshFromFile(fileName.c_str(),meshNames[0].c_str(),0);
  mesh->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(538,mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(579,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,(int)mesh->getAllGeoTypes().size());
  for(int i=0;i<514;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(i));
  for(int i=514;i<538;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_POLYGON,mesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,std::accumulate(mesh->getCoords()->getPointer(),mesh->getCoords()->getPointer()+1737,0),1e-12);
  const double expectedVals1[12]={1.4851585216522212,-0.5,0.,1.4851585216522212,-0.4,0.,1.4851585216522212,-0.3,0., 1.5741585216522211, -0.5, 0. };
  double diffValue1[12];
  std::transform(mesh->getCoords()->getPointer(),mesh->getCoords()->getPointer()+12,expectedVals1,diffValue1,std::minus<double>());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::max_element(diffValue1,diffValue1+12),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::min_element(diffValue1,diffValue1+12),1e-12);
  CPPUNIT_ASSERT_EQUAL((std::size_t)2768,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(651050,std::accumulate(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+2768,0));
  CPPUNIT_ASSERT_EQUAL(725943,std::accumulate(mesh->getNodalConnectivityIndex()->getPointer(),mesh->getNodalConnectivityIndex()->getPointer()+539,0));
  mesh->decrRef();
  //
  std::vector<std::string> fieldsName=MEDLoader::GetCellFieldNamesOnMesh(fileName.c_str(),meshNames[0].c_str());
  CPPUNIT_ASSERT_EQUAL(3,(int)fieldsName.size());
  CPPUNIT_ASSERT(fieldsName[0]=="bord_:_distorsion");
  CPPUNIT_ASSERT(fieldsName[1]=="bord_:_familles");
  CPPUNIT_ASSERT(fieldsName[2]=="bord_:_non-ortho");
  std::vector<std::pair<int,int> > its0=MEDLoader::GetCellFieldIterations(fileName.c_str(),meshNames[0].c_str(),fieldsName[0].c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)its0.size());
  MEDCouplingFieldDouble *field=MEDLoader::ReadFieldCell(fileName.c_str(),meshNames[0].c_str(),0,fieldsName[0].c_str(),its0[0].first,its0[0].second);
  field->checkCoherency();
  CPPUNIT_ASSERT(field->getName()==fieldsName[0]);
  CPPUNIT_ASSERT_EQUAL(1,field->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(538,field->getNumberOfTuples());
  const MEDCouplingUMesh *constMesh=dynamic_cast<const MEDCouplingUMesh *>(field->getMesh());
  CPPUNIT_ASSERT(constMesh);
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,constMesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(538,constMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(579,constMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,(int)constMesh->getAllGeoTypes().size());
  for(int i=0;i<514;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,constMesh->getTypeOfCell(i));
  for(int i=514;i<538;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_POLYGON,constMesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,std::accumulate(constMesh->getCoords()->getConstPointer(),constMesh->getCoords()->getConstPointer()+1737,0),1e-12);
  std::transform(constMesh->getCoords()->getConstPointer(),constMesh->getCoords()->getConstPointer()+12,expectedVals1,diffValue1,std::minus<double>());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::max_element(diffValue1,diffValue1+12),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::min_element(diffValue1,diffValue1+12),1e-12);
  CPPUNIT_ASSERT_EQUAL((std::size_t)2768,constMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(651050,std::accumulate(constMesh->getNodalConnectivity()->getConstPointer(),constMesh->getNodalConnectivity()->getConstPointer()+2768,0));
  CPPUNIT_ASSERT_EQUAL(725943,std::accumulate(constMesh->getNodalConnectivityIndex()->getConstPointer(),constMesh->getNodalConnectivityIndex()->getConstPointer()+539,0));
  const double *values=field->getArray()->getPointer();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.87214203182918,std::accumulate(values,values+538,0.),1e-12);
  field->decrRef();
}

void ParaMEDMEMTest::testMEDLoaderPolyhedronRead()
{
  string fileName=getResourceFile("poly3D.med");
  vector<string> meshNames=MEDLoader::GetMeshNames(fileName.c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)meshNames.size());
  CPPUNIT_ASSERT(meshNames[0]=="poly3D");
  MEDCouplingUMesh *mesh=MEDLoader::ReadUMeshFromFile(fileName.c_str(),meshNames[0].c_str(),0);
  mesh->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(3,mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,(int)mesh->getAllGeoTypes().size());
  CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,mesh->getTypeOfCell(0));
  CPPUNIT_ASSERT_EQUAL(NORM_POLYHED,mesh->getTypeOfCell(1));
  CPPUNIT_ASSERT_EQUAL(NORM_POLYHED,mesh->getTypeOfCell(2));
  CPPUNIT_ASSERT_EQUAL((std::size_t)98,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(725,std::accumulate(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+98,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(110.,std::accumulate(mesh->getCoords()->getPointer(),mesh->getCoords()->getPointer()+57,0),1e-12);
  CPPUNIT_ASSERT_EQUAL(155,std::accumulate(mesh->getNodalConnectivityIndex()->getPointer(),mesh->getNodalConnectivityIndex()->getPointer()+4,0));
  mesh->decrRef();
  //
  mesh=MEDLoader::ReadUMeshFromFile(fileName.c_str(),meshNames[0].c_str(),-1);
  mesh->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(17,mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,(int)mesh->getAllGeoTypes().size());
  CPPUNIT_ASSERT_EQUAL(NORM_POLYGON,mesh->getTypeOfCell(0));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(1));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(2));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(3));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(4));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(5));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(6));
  CPPUNIT_ASSERT_EQUAL(NORM_TRI3,mesh->getTypeOfCell(7));
  CPPUNIT_ASSERT_EQUAL(NORM_POLYGON,mesh->getTypeOfCell(8));
  CPPUNIT_ASSERT_EQUAL(NORM_POLYGON,mesh->getTypeOfCell(9));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(10));
  CPPUNIT_ASSERT_EQUAL(NORM_TRI3,mesh->getTypeOfCell(11));
  CPPUNIT_ASSERT_EQUAL(NORM_TRI3,mesh->getTypeOfCell(12));
  CPPUNIT_ASSERT_EQUAL(NORM_TRI3,mesh->getTypeOfCell(13));
  CPPUNIT_ASSERT_EQUAL(NORM_TRI3,mesh->getTypeOfCell(14));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(15));
  CPPUNIT_ASSERT_EQUAL(NORM_TRI3,mesh->getTypeOfCell(16));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(110.,std::accumulate(mesh->getCoords()->getPointer(),mesh->getCoords()->getPointer()+57,0),1e-12);
  CPPUNIT_ASSERT_EQUAL((std::size_t)83,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(619,std::accumulate(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+83,0));
  mesh->decrRef();
  //
  vector<string> families=MEDLoader::GetMeshFamiliesNames(fileName.c_str(),meshNames[0].c_str());
  CPPUNIT_ASSERT_EQUAL(4,(int)families.size());
  CPPUNIT_ASSERT(families[0]=="FAMILLE_FACE_POLYGONS3");
  CPPUNIT_ASSERT(families[1]=="FAMILLE_FACE_QUAD41");
  CPPUNIT_ASSERT(families[2]=="FAMILLE_FACE_TRIA32");
  CPPUNIT_ASSERT(families[3]=="FAMILLE_ZERO");
  vector<string> families2;
  families2.push_back(families[0]);
  mesh=MEDLoader::ReadUMeshFromFamilies(fileName.c_str(),meshNames[0].c_str(),-1,families2);
  mesh->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(3,mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(1,(int)mesh->getAllGeoTypes().size());
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_POLYGON,mesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_EQUAL((std::size_t)19,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(117,std::accumulate(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+19,0));
  mesh->decrRef();
  //
  mesh=MEDLoader::ReadUMeshFromFamilies(fileName.c_str(),meshNames[0].c_str(),0,families2);
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(0,mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(0,(int)mesh->getAllGeoTypes().size());
  mesh->decrRef();
}
