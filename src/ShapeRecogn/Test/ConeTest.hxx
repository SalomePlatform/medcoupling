// Copyright (C) 2007-2024  CEA, EDF
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

#ifndef __CONETEST_HXX__
#define __CONETEST_HXX__

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

namespace MEDCoupling
{
    class ShapeRecognMeshBuilder;
    class Areas;

    class ConeTest : public CppUnit::TestFixture
    {
        CPPUNIT_TEST_SUITE(ConeTest);
        CPPUNIT_TEST(testNumberOfAreas);
        CPPUNIT_TEST(testComputePlaneProperties);
        CPPUNIT_TEST(testComputeCylinderProperties);
        CPPUNIT_TEST(testComputeConeProperties);
        CPPUNIT_TEST(testFirstArea);
        CPPUNIT_TEST(testSecondArea);
        CPPUNIT_TEST(testThirdArea);
        CPPUNIT_TEST_SUITE_END();

    public:
        void setUp() override;
        void tearDown() override;

        void testComputePlaneProperties();
        void testComputeCylinderProperties();
        void testComputeConeProperties();

        void testNumberOfAreas();
        void testFirstArea();
        void testSecondArea();
        void testThirdArea();

    private:
        ShapeRecognMeshBuilder *srMesh = 0;
        const Areas *areas;
    };
};

#endif // __CONETEST_HXX__
