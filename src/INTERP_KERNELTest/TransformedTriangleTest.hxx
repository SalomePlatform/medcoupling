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

#ifndef __TU_TRANSFORMED_TRIANGLE_HXX__
#define __TU_TRANSFORMED_TRIANGLE_HXX__

#include <cppunit/extensions/HelperMacros.h>

#include "InterpKernelTestExport.hxx"
#include "TransformedTriangle.hxx"

#define ERR_TOL 1.0e-8

namespace INTERP_TEST
{

  /**
   * \brief Test suite testing some of the low level methods of TransformedTriangle.
   *
   */
  class INTERPKERNELTEST_EXPORT TransformedTriangleTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( TransformedTriangleTest );
    CPPUNIT_TEST( test_constructor );
    CPPUNIT_TEST( test_calcUnstableC );
    CPPUNIT_TEST( test_calcUnstableT );
    //removed because the test fails to enter the desired code branch
 //   CPPUNIT_TEST( test_calcStableC_Consistency );
    CPPUNIT_TEST_SUITE_END();

    typedef INTERP_KERNEL::TransformedTriangle::TriSegment TriSegment;
    typedef INTERP_KERNEL::TransformedTriangle::DoubleProduct DoubleProduct;

  public:
    void setUp();

    void tearDown();

    // tests
    void test_constructor();

    void test_calcUnstableC(); 

    void test_calcUnstableT();

    void test_calcStableC_Consistency();

    double p1[3], q1[3], r1[3];
    double hp1, hq1, hr1;
    double Hp1, Hq1, Hr1;

    double p2[3], q2[3], r2[3];
    double hp2, hq2, hr2;
    double Hp2, Hq2, Hr2;

    double stable_c2[24];
  
  private:
    INTERP_KERNEL::TransformedTriangle* tri1;
    INTERP_KERNEL::TransformedTriangle* tri2;

  };




}



#endif
