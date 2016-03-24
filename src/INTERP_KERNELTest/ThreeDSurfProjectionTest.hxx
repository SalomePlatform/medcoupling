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

#ifndef __THREEDSURFPROJECTIONTEST_HXX__
#define __THREEDSURFPROJECTIONTEST_HXX__

#include <cppunit/extensions/HelperMacros.h>

#include "InterpKernelTestExport.hxx"

namespace INTERP_TEST 
{
  /**
   * \brief Class dedicated of the test of the preprocessing of 3D surf cells before performing invoking 2D algorithms.
   */
  class INTERPKERNELTEST_EXPORT ThreeDSurfProjectionTest : public CppUnit::TestFixture
  {
    CPPUNIT_TEST_SUITE( ThreeDSurfProjectionTest );
    CPPUNIT_TEST ( test1 );
    CPPUNIT_TEST ( test2 );
    CPPUNIT_TEST_SUITE_END();
  public:
    void test1();
    void test2();
  };
}

#endif
