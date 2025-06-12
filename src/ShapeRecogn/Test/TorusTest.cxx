// Copyright (C) 2024  CEA, EDF
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

#include "TorusTest.hxx"

#include "ShapeRecognMeshBuilder.hxx"
#include "Areas.hxx"
#include "MathOps.hxx"
#include "TestInterpKernelUtils.hxx"  // getResourceFile()
#include "ShapeRecognMesh.hxx"

#include "ShapeRecognTest.hxx"

using namespace MEDCoupling;

void
TorusTest::setUp()
{
    std::string file = INTERP_TEST::getResourceFile("ShapeRecognTorus.med", 3);
    srMesh = BuildShapeRecognMeshBuilderFromFile(file);
    srMesh->recognize();
    areas = srMesh->getAreas();
}

void
TorusTest::tearDown()
{
    areas = 0;
}

void
TorusTest::testArea()
{
    CPPUNIT_ASSERT_EQUAL(275, (int)srMesh->getNodes()->getNbNodes());
    CPPUNIT_ASSERT_EQUAL(1, (int)areas->getNumberOfAreas());
    CPPUNIT_ASSERT_EQUAL((int)PrimitiveType::Torus, (int)areas->getPrimitiveType(0));
    // Some nodes are unknown
    CPPUNIT_ASSERT_EQUAL(272, (int)areas->getNumberOfNodes(0));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.843297, areas->getMinorRadius(0), 1E-2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.156428, areas->getRadius(0), 1E-1);
    std::array<double, 3> centerRef = {7.687022, -3.726887, -9.02};
    std::array<double, 3> center = areas->getCenter(0);
    for (size_t j = 0; j < 3; ++j) CPPUNIT_ASSERT_DOUBLES_EQUAL(centerRef[j], center[j], 1E-2);
}
