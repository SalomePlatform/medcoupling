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

#ifndef __TU_TRANSFORMED_TRIANGLE_INTERSECT_HXX__
#define __TU_TRANSFORMED_TRIANGLE_INTERSECT_HXX__

#include <cppunit/extensions/HelperMacros.h>

#include "InterpKernelTestExport.hxx"
#include "TransformedTriangle.hxx"

namespace INTERP_TEST
{

  class INTERPKERNELTEST_EXPORT TransformedTriangleIntersectTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( TransformedTriangleIntersectTest );

    CPPUNIT_TEST( testTriangle1 );
    CPPUNIT_TEST( testTriangle2 );
    CPPUNIT_TEST( testTriangle3 );
    CPPUNIT_TEST( testTriangle4 );
    CPPUNIT_TEST( testTriangle5 );
    CPPUNIT_TEST( testTriangle6 );
    CPPUNIT_TEST( testTriangle7 );
    CPPUNIT_TEST( testTriangle8 );
    CPPUNIT_TEST( testTriangle9 );
    CPPUNIT_TEST( testTriangle10 );
    CPPUNIT_TEST( testTriangle11 );
    CPPUNIT_TEST( testTriangle12 );
    CPPUNIT_TEST( testTriangle13 );

    CPPUNIT_TEST_SUITE_END();

    typedef INTERP_KERNEL::TransformedTriangle::TriSegment TriSegment;
    typedef INTERP_KERNEL::TransformedTriangle::DoubleProduct DoubleProduct;

  public:

    void testTriangle1();
  
    void testTriangle2();

    void testTriangle3();

    void testTriangle4();

    void testTriangle5();

    void testTriangle6();

    void testTriangle7();

    void testTriangle8();

    void testTriangle9();
  
    void testTriangle10();
  
    void testTriangle11();
  
    void testTriangle12();

    void testTriangle13();

  private:
 
  };

}






#endif
