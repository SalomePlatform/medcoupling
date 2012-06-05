// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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

#ifndef __TU_INTERPOLATIONOPTIONS_HXX__
#define __TU_INTERPOLATIONOPTIONS_HXX__

#include <cppunit/extensions/HelperMacros.h>

#include "InterpKernelTestExport.hxx"
#include "InterpolationOptions.hxx"
#include "MEDMEM_Field.hxx"

namespace INTERP_TEST
{

  /**
   * \brief Test suite testing some of the low level methods of TransformedTriangle.
   *
   */
  class INTERPKERNELTEST_EXPORT InterpolationOptionsTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( InterpolationOptionsTest );
    CPPUNIT_TEST( test_InterpolationOptions );
    CPPUNIT_TEST_SUITE_END();

   
  public:
    void setUp();

    void tearDown();

    // tests
    void test_InterpolationOptions();

  private:
   
  };




}



#endif
