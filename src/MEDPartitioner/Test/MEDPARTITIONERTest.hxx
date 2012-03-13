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

#ifndef __MEDPARTITIONERTEST_HXX__
#define __MEDPARTITIONERTEST_HXX__

#include <cppunit/extensions/HelperMacros.h>

#include <set>
#include <string>
#include <iostream>

#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"

class MEDPARTITIONERTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( MEDPARTITIONERTest );
  CPPUNIT_TEST( testSmallSize );
  CPPUNIT_TEST( testMedianSize );
  //CPPUNIT_TEST( deleteTestMeshes );
  CPPUNIT_TEST_SUITE_END();

public:

  //global use
  int _ni;  //nb of hexa9
  int _nj;
  int _nk;
  int _ntot;
  std::string _fileName; //initial test mesh file med CUBE3D
  std::string _fileNameWithFaces; //initial test mesh file med CUBE3D plus a set of faces
  std::string _fileName2; //initial test mesh file med CARRE3D
  std::string _fileNameHugeXml;
  int _nbTargetHuge;
  std::string _meshName; //initial test mesh file med
  int _verbose;
  
  //for utils
  void setSize(int ni, int nj, int nk);
  void setSmallSize();
  void setMedianSize();
  void setbigSize();
  ParaMEDMEM::MEDCouplingUMesh * buildCUBE3DMesh();
  ParaMEDMEM::MEDCouplingUMesh * buildFACE3DMesh();
  ParaMEDMEM::MEDCouplingUMesh * buildCARRE3DMesh();
  ParaMEDMEM::MEDCouplingFieldDouble * buildVecFieldOnCells(std::string myfileName);
  ParaMEDMEM::MEDCouplingFieldDouble * buildVecFieldOnNodes();
  void createTestMeshWithoutField();
  void createTestMeshWithVecFieldOnCells();
  void createTestMeshWithVecFieldOnNodes();
  void verifyTestMeshWithVecFieldOnNodes();
  void verifyMedpartitionerOnSmallSizeForMesh();
  void verifyMedpartitionerOnSmallSizeForFieldOnCells();
  void verifyMedpartitionerOnSmallSizeForFieldOnGaussNe();
  void createTestMeshes();
  void createHugeTestMesh(int ni, int nj, int nk, int nbx, int nby, int nbz, int nbTarget);
  void launchMedpartitionerOnTestMeshes();
  void launchMedpartitionerOnHugeTestMeshes();
  void deleteTestMeshes();
  
  //for CPPUNIT_TEST
  void setUp();
  void tearDown();
  void testSmallSize();
  void testMedianSize();
};

#endif
