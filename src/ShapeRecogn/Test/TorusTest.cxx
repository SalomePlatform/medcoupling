#include "TorusTest.hxx"

#include "ShapeRecognMeshBuilder.hxx"
#include "Areas.hxx"
#include "MathOps.hxx"
#include "TestInterpKernelUtils.hxx" // getResourceFile()

using namespace MEDCoupling;

void TorusTest::setUp()
{
    std::string file = INTERP_TEST::getResourceFile("ShapeRecognTorus.med", 3);
    srMesh = new ShapeRecognMeshBuilder(file);
    srMesh->recognize();
    areas = srMesh->getAreas();
}

void TorusTest::tearDown()
{
    if (srMesh != 0)
        delete srMesh;
    areas = 0;
}

void TorusTest::testArea()
{
    CPPUNIT_ASSERT_EQUAL(275, (int)srMesh->getNodes()->getNbNodes());
    CPPUNIT_ASSERT_EQUAL(1, (int)areas->getNumberOfAreas());
    CPPUNIT_ASSERT_EQUAL(PrimitiveType::Torus, areas->getPrimitiveType(0));
    // Some nodes are unknown
    CPPUNIT_ASSERT_EQUAL(272, (int)areas->getNumberOfNodes(0));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.843297, areas->getMinorRadius(0), 1E-2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.156428, areas->getRadius(0), 1E-1);
    std::array<double, 3> centerRef = {7.687022, -3.726887, -9.02};
    std::array<double, 3> center = areas->getCenter(0);
    for (size_t j = 0; j < 3; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(centerRef[j], center[j], 1E-2);
}
