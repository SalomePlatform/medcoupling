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

#include "SphereTest.hxx"

#include "ShapeRecognMeshBuilder.hxx"
#include "ShapeRecognMesh.hxx"
#include "Areas.hxx"
#include "MathOps.hxx"
#include "TestInterpKernelUtils.hxx"  // getResourceFile()

#include "ShapeRecognTest.hxx"

using namespace MEDCoupling;

void
SphereTest::setUp()
{
    std::string file = INTERP_TEST::getResourceFile("ShapeRecognSphere.med", 3);
    srMesh = BuildShapeRecognMeshBuilderFromFile(file);
    srMesh->recognize();
    areas = srMesh->getAreas();
}

void
SphereTest::tearDown()
{
    areas = 0;
}

void
SphereTest::testArea()
{
    CPPUNIT_ASSERT_EQUAL(1, (int)areas->getNumberOfAreas());
    // 8 double nodes so 147 - 6 nodes
    CPPUNIT_ASSERT_EQUAL(141, (int)areas->getNumberOfNodes(0));
    CPPUNIT_ASSERT_EQUAL((int)PrimitiveType::Sphere, (int)areas->getPrimitiveType(0));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, areas->getRadius(0), 1E-2);
    std::array<double, 3> centerRef = {5.3, -6.7, -9.02};
    std::array<double, 3> center = areas->getCenter(0);
    for (size_t j = 0; j < 3; ++j) CPPUNIT_ASSERT_DOUBLES_EQUAL(centerRef[j], center[j], 1E-2);
}
