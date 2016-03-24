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

#ifndef _EXPREVALINTERPTEST_HXX_
#define _EXPREVALINTERPTEST_HXX_

#include <cppunit/extensions/HelperMacros.h>

#include "InterpKernelTestExport.hxx"

namespace INTERP_TEST
{
  class INTERPKERNELTEST_EXPORT ExprEvalInterpTest : public CppUnit::TestFixture
  {
    CPPUNIT_TEST_SUITE( ExprEvalInterpTest );
    CPPUNIT_TEST( testBuildStringFromFortran );
    CPPUNIT_TEST( testDeleteWhiteSpaces );
    CPPUNIT_TEST( testInterpreter0 );
    CPPUNIT_TEST( testInterpreter1 );
    CPPUNIT_TEST( testInterpreter2 );
    CPPUNIT_TEST( testInterpreterUnit0 );
    CPPUNIT_TEST( testInterpreterUnit1 );
    CPPUNIT_TEST( testInterpreter3 );
    CPPUNIT_TEST( testInterpreter4 );
    CPPUNIT_TEST( testInterpreter5 );
    CPPUNIT_TEST( testInterpreter6 );
    CPPUNIT_TEST_SUITE_END();
  public:
    void setUp() { }
    void tearDown() { }
    void cleanUp() { }
    void testBuildStringFromFortran();
    void testDeleteWhiteSpaces();
    void testInterpreter0();
    void testInterpreter1();
    void testInterpreter2();
    void testInterpreter3();
    void testInterpreter4();
    void testInterpreter5();
    void testInterpreter6();
    void testInterpreterUnit0();
    void testInterpreterUnit1();
  };
}

#endif
