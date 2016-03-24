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

#ifndef __MEDCOUPLINGREMAPPERTEST_HXX__
#define __MEDCOUPLINGREMAPPERTEST_HXX__

#include <cppunit/extensions/HelperMacros.h>

#include <map>
#include <vector>

namespace MEDCoupling
{
  class MEDCouplingUMesh;

  class MEDCouplingRemapperTest : public CppUnit::TestFixture
  {
    CPPUNIT_TEST_SUITE(MEDCouplingRemapperTest);
    CPPUNIT_TEST( test2DInterpP0P0_1 );
    CPPUNIT_TEST( test2DInterpP0P0R_1 );
    CPPUNIT_TEST( test1DInterp_1 );
    CPPUNIT_TEST( test2DInterpMultiMethods );
    CPPUNIT_TEST( testMultiDimCombi );
    CPPUNIT_TEST( testNatureOfField );
    CPPUNIT_TEST( testExtruded );
    CPPUNIT_TEST( testExtruded2 );
    CPPUNIT_TEST( testPrepareEx1 );
    CPPUNIT_TEST( testPartialTransfer1 );
    CPPUNIT_TEST( testBugNonRegression1 );
    CPPUNIT_TEST_SUITE_END();
  public:
    void test2DInterpP0P0_1();
    void test2DInterpP0P0R_1();
    void test1DInterp_1();
    void test2DInterpMultiMethods();
    void testMultiDimCombi();
    void testNatureOfField();
    void testExtruded();
    void testExtruded2();
    void testPrepareEx1();
    void testPartialTransfer1();
    //
    void testBugNonRegression1();
  private:
    static MEDCouplingUMesh *build1DTargetMesh_2();
    static MEDCouplingUMesh *build2DTargetMesh_3();
    static MEDCouplingUMesh *build3DExtrudedUMesh_1(MEDCouplingUMesh *&mesh2D);
  };
}

#endif
