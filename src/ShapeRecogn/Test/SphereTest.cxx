#include "SphereTest.hxx"

#include "ShapeRecognMeshBuilder.hxx"
#include "Areas.hxx"
#include "MathOps.hxx"
#include "TestInterpKernelUtils.hxx" // getResourceFile()

using namespace MEDCoupling;

void SphereTest::setUp()
{
    std::string file = INTERP_TEST::getResourceFile("ShapeRecognSphere.med", 3);
    srMesh = new ShapeRecognMeshBuilder(file);
    srMesh->recognize();
    areas = srMesh->getAreas();
}

void SphereTest::tearDown()
{
    if (srMesh != 0)
        delete srMesh;
    areas = 0;
}

void SphereTest::testArea()
{
    CPPUNIT_ASSERT_EQUAL(1, (int)areas->getNumberOfAreas());
    // 8 double nodes so 147 - 6 nodes
    CPPUNIT_ASSERT_EQUAL(141, (int)areas->getNumberOfNodes(0));
    CPPUNIT_ASSERT_EQUAL(PrimitiveType::Sphere, areas->getPrimitiveType(0));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, areas->getRadius(0), 1E-2);
    std::array<double, 3> centerRef = {5.3, -6.7, -9.02};
    std::array<double, 3> center = areas->getCenter(0);
    for (size_t j = 0; j < 3; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(centerRef[j], center[j], 1E-2);
}
