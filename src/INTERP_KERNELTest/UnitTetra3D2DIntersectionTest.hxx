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

#ifndef __UNITTETRA3D2DINTERSECTIONTEST_HXX__
#define __UNITTETRA3D2DINTERSECTIONTEST_HXX__

#include <cppunit/extensions/HelperMacros.h>

#include "InterpKernelTestExport.hxx"

namespace INTERP_TEST
{
  /**
   * \brief Test suite testing UnitTetra3D2DIntersection class.
   *
   */
  class INTERPKERNELTEST_EXPORT UnitTetra3D2DIntersectionTest : public CppUnit::TestFixture
  {
    CPPUNIT_TEST_SUITE( UnitTetra3D2DIntersectionTest );
    CPPUNIT_TEST( test_UnitTetra3D2DIntersection_1 );
    CPPUNIT_TEST( test_UnitTetra3D2DIntersection_2 );
    CPPUNIT_TEST( test_UnitTetra3D2DIntersection_3 );
    CPPUNIT_TEST_SUITE_END();
  public:
    void test_UnitTetra3D2DIntersection_1();
    void test_UnitTetra3D2DIntersection_2();
    void test_UnitTetra3D2DIntersection_3();
  };
}

#endif
