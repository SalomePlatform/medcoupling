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

#ifndef __TU_BB_TREE_HXX__
#define __TU_BB_TREE_HXX__

#include <cppunit/extensions/HelperMacros.h>

#include "InterpKernelTestExport.hxx"
#include "BBTree.txx"

namespace INTERP_TEST
{

  /**
   * \brief Test suite testing some of the low level methods of TransformedTriangle.
   *
   */
  class INTERPKERNELTEST_EXPORT BBTreeTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( BBTreeTest );
    CPPUNIT_TEST( test_BBTree );
    CPPUNIT_TEST( test_DirectedBB_1D );
    CPPUNIT_TEST( test_DirectedBB_2D );
    CPPUNIT_TEST( test_DirectedBB_3D );
    CPPUNIT_TEST_SUITE_END();

   
  public:
    void setUp();

    void tearDown();

    // tests
    void test_BBTree();
    void test_DirectedBB_1D();
    void test_DirectedBB_2D();
    void test_DirectedBB_3D();

  };




}



#endif
