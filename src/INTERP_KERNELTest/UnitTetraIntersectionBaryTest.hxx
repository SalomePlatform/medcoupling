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

// File      : UnitTetraIntersectionBaryTests.hxx
// Created   : Thu Nov  6 17:11:27 2008
// Author    : Edward AGAPOV (eap)
//
#ifndef __UNITTETRAINTERSECTIONBARYTEST_HXX__
#define __UNITTETRAINTERSECTIONBARYTEST_HXX__

#include <cppunit/extensions/HelperMacros.h>

#include "InterpKernelTestExport.hxx"

namespace INTERP_TEST
{
  /**
   * \brief Test suite testing UnitTetraIntersectionBary class.
   *
   */
  class INTERPKERNELTEST_EXPORT UnitTetraIntersectionBaryTest : public CppUnit::TestFixture
  {
    CPPUNIT_TEST_SUITE( UnitTetraIntersectionBaryTest );
    CPPUNIT_TEST( test_UnitTetraIntersectionBary_13 );
    CPPUNIT_TEST( test_UnitTetraIntersectionBary_12 );
    CPPUNIT_TEST( test_UnitTetraIntersectionBary_1 );
    CPPUNIT_TEST( test_UnitTetraIntersectionBary_2 );
    CPPUNIT_TEST( test_UnitTetraIntersectionBary_3 );
    CPPUNIT_TEST( test_UnitTetraIntersectionBary_4 );
    CPPUNIT_TEST( test_UnitTetraIntersectionBary_5 );
    CPPUNIT_TEST( test_UnitTetraIntersectionBary_6 );
    CPPUNIT_TEST( test_UnitTetraIntersectionBary_7 );
    CPPUNIT_TEST( test_UnitTetraIntersectionBary_8 );
    CPPUNIT_TEST( test_UnitTetraIntersectionBary_9 );
    CPPUNIT_TEST( test_UnitTetraIntersectionBary_10 );
    CPPUNIT_TEST( test_UnitTetraIntersectionBary_11 );
    CPPUNIT_TEST( test_TetraAffineTransform_reverseApply );
    CPPUNIT_TEST( test_barycentric_coords );
    CPPUNIT_TEST( test_cuboid_mapped_coords_3D );
    CPPUNIT_TEST( test_quad_mapped_coords_2D );
    CPPUNIT_TEST_SUITE_END();
  public:
    void test_UnitTetraIntersectionBary_1();
    void test_UnitTetraIntersectionBary_2();
    void test_UnitTetraIntersectionBary_3();
    void test_UnitTetraIntersectionBary_4();
    void test_UnitTetraIntersectionBary_5();
    void test_UnitTetraIntersectionBary_6();
    void test_UnitTetraIntersectionBary_7();
    void test_UnitTetraIntersectionBary_8();
    void test_UnitTetraIntersectionBary_9();
    void test_UnitTetraIntersectionBary_10();
    void test_UnitTetraIntersectionBary_11();
    void test_UnitTetraIntersectionBary_12();
    void test_UnitTetraIntersectionBary_13();
    void test_TetraAffineTransform_reverseApply();
    void test_barycentric_coords();
    void test_cuboid_mapped_coords_3D();
    void test_quad_mapped_coords_2D();
  };
}

#endif
