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

#include "MathOpsTest.hxx"
#include "MathOps.hxx"

using namespace MEDCoupling;

void MathOpsTest::testLstsq()
{
    std::vector<double> a = {
        0.69473263, 0.83318004, 0.60822673, 0.59243878, 0.82872553,
        0.84048546, 0.95698819, 0.02171218, 0.27683381, 0.20628928,
        0.80895323, 0.4207767, 0.37468575, 0.86258204, 0.42571846};
    std::vector<double> b = {
        0.91167508, 0.95878824, 0.7234827, 0.51753917, 0.18603306};
    std::vector<double> x = MathOps::lstsq(a, b);
    std::array<double, 3> xRef = {0.35719095, 0.5134345, 0.26039343};
    for (size_t i = 0; i < 3; ++i)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(xRef[i], x[i], 1E-6);
}

void MathOpsTest::testLstsq2()
{
    std::vector<double> a = {
        0.4564562, 0.3517006,
        0.28928215, 0.72309086,
        0.05944836, 0.56024464};
    std::vector<double> b = {0.98902712, 0.46791812};
    std::vector<double> x = MathOps::lstsq(a, b);
    std::array<double, 3> xRef = {2.10752524, 0.2636243, -0.82807416};
    for (size_t i = 0; i < 3; ++i)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(xRef[i], x[i], 1E-6);
}

void MathOpsTest::testLstsqBig()
{
    std::vector<double> a(5000000 * 3, 1.0);
    std::vector<double> b(5000000, 1.0);
    std::vector<double> x = MathOps::lstsq(a, b);
    std::array<double, 3> xRef = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    for (size_t i = 2; i < 3; ++i)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(xRef[i], x[i], 1E-6);
}

void MathOpsTest::testComputeCov()
{
    std::vector<double> coordinates{
        1, 2, -3, -8, 6, -3, 9, 2, -1, -10, -11, -120};
    std::array<double, 9>
        covMatrix = MathOps::computeCov(coordinates);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(76.666666, covMatrix[0], 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(26.666666, covMatrix[1], 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(319.333333, covMatrix[2], 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(26.666666, covMatrix[3], 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(54.916666, covMatrix[4], 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(420.75, covMatrix[5], 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(319.333333, covMatrix[6], 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(420.75, covMatrix[7], 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3462.25, covMatrix[8], 1E-6);
}

void MathOpsTest::testComputePCAFirstAxis()
{
    std::vector<double> coordinates{
        1, 2, -3, -8, 6, -3, 9, 2, -1, -10, -11, -120};
    std::array<double, 3>
        axis = MathOps::computePCAFirstAxis(coordinates);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.09198798, axis[0], 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.11994164, axis[1], 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.9885101, axis[2], 1E-6);
}

void MathOpsTest::testComputeAngles()
{
    std::vector<double> directions{
        1, 2, -3, -8, 6, -3, 9, 2, -1, -10, -11, -120};
    std::array<double, 3> axis{1, 2, 3};
    std::vector<double> angles = MathOps::computeAngles(
        directions, axis);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.86054803, angles[0], 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.69914333, angles[1], 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.27845478, angles[2], 1E-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.61880376, angles[3], 1E-6);
}

void MathOpsTest::testComputeBaseFromNormal()
{
    std::array<double, 3> normal = {1.0, 2.0, 3.0};
    std::array<double, 6> base = MathOps::computeBaseFromNormal(normal);
    std::array<double, 6> baseRef = {
        -0.53452248, 0.77454192, -0.33818712,
        -0.80178373, -0.33818712, 0.49271932};
    for (size_t i = 0; i < 6; ++i)
    {
        std::ostringstream message;
        message << "Mismatch at index " << i;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message.str().c_str(), baseRef[i], base[i], 1E-6);
    }
}
