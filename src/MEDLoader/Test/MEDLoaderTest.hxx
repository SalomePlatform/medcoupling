//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

#ifndef __MEDLOADERTEST_HXX__
#define __MEDLOADERTEST_HXX__

#include <cppunit/extensions/HelperMacros.h>

namespace ParaMEDMEM
{
  class MEDCouplingUMesh;

  class MEDLoaderTest : public CppUnit::TestFixture
  {
    CPPUNIT_TEST_SUITE(MEDLoaderTest);
    CPPUNIT_TEST( testMesh1DRW );
    CPPUNIT_TEST( testMesh2DCurveRW );
    CPPUNIT_TEST( testMesh2DRW );
    CPPUNIT_TEST( testMesh3DSurfRW );
    CPPUNIT_TEST( testMesh3DRW );
    CPPUNIT_TEST_SUITE_END();
  public:
    void testMesh1DRW();
    void testMesh2DCurveRW();
    void testMesh2DRW();
    void testMesh3DSurfRW();
    void testMesh3DRW();
  private:
    MEDCouplingUMesh *build1DMesh_1();
    MEDCouplingUMesh *build2DCurveMesh_1();
    MEDCouplingUMesh *build2DMesh_1();
    MEDCouplingUMesh *build3DSurfMesh_1();
    MEDCouplingUMesh *build3DMesh_1();
  };
}

#endif
