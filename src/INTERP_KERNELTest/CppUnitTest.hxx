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

#ifndef __TU_TEST_CPPUNIT_HXX__
#define __TU_TEST_CPPUNIT_HXX__

#include <cppunit/extensions/HelperMacros.h>

#include "InterpKernelTestExport.hxx"

/**
 * \brief Class tested by TestBogusClass : not very useful
 */
class INTERPKERNELTEST_EXPORT BogusClass {
  friend class TestBogusClass;

public:
  BogusClass(double _x) : x(_x) {;} 
 
  double getX() { return x; }

private: 
  double x;
};
  
/**
 * \brief Class used to figure out CppUnit : not very useful
 *
 */
class INTERPKERNELTEST_EXPORT TestBogusClass : public CppUnit::TestFixture
{

  CPPUNIT_TEST_SUITE( TestBogusClass );
  CPPUNIT_TEST( test1 );
  CPPUNIT_TEST( test2 );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
    obj = new BogusClass(3.14);
  }

  void tearDown() {
    delete obj;
  }

  void test1() {
    // test something
    CPPUNIT_ASSERT(obj->x == 3.14);
    CPPUNIT_ASSERT(obj->getX() == obj->x);
  }

  void test2() {
    // test something else
    obj->x += 2.6;
    CPPUNIT_ASSERT(obj->getX() > 3.14);
  }

private:
  BogusClass* obj;

};








#endif
