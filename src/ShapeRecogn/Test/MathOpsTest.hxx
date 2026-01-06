// Copyright (C) 2024-2026  CEA, EDF
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

#pragma once

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

namespace MEDCoupling
{
class MathOpsTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(MathOpsTest);
    CPPUNIT_TEST(testLstsq);
    CPPUNIT_TEST(testLstsq2);
    CPPUNIT_TEST(testLstsqBig);
    CPPUNIT_TEST(testComputeCov);
    CPPUNIT_TEST(testComputePCAFirstAxis);
    CPPUNIT_TEST(testComputeAngles);
    CPPUNIT_TEST(testComputeBaseFromNormal);
    CPPUNIT_TEST_SUITE_END();

   public:
    void testLstsq();
    void testLstsq2();
    void testLstsqBig();
    void testComputeCov();
    void testComputePCAFirstAxis();
    void testComputeAngles();
    void testComputeBaseFromNormal();
};
}  // namespace MEDCoupling
