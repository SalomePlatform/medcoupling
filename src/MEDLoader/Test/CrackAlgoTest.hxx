// Copyright (C) 2007-2026  CEA, EDF
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
// Author : Aymeric SONOLET (CEA/DES)

#ifndef __CRACKALGOTEST_HXX__
#define __CRACKALGOTEST_HXX__

#include <cppunit/extensions/HelperMacros.h>
#include <string>
#include <vector>
#include <utility>

#include <MCIdType.hxx>

namespace MEDCoupling
{

class MEDFileUMesh;
class MEDCouplingUMesh;

class CrackAlgoTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(CrackAlgoTest);
    CPPUNIT_TEST(test3DHalfCrossCut);
    CPPUNIT_TEST(test3DFullCrossCut);
    CPPUNIT_TEST(test3DHalfCut);
    CPPUNIT_TEST(test3DFullCut);
    CPPUNIT_TEST(test3DFullFullCut);
    CPPUNIT_TEST(test3DAngleCut);
    CPPUNIT_TEST(test2DGrp);
    CPPUNIT_TEST(test2DMeshNonConnexCut);
    CPPUNIT_TEST(test1DMesh);
    CPPUNIT_TEST(testInnerCrossCut);
    CPPUNIT_TEST_SUITE_END();

   public:
    void test3DHalfCrossCut();
    void test3DFullCrossCut();
    void test3DAngleCut();
    void test3DFullCut();
    void test3DFullFullCut();
    void test3DHalfCut();
    void test2DGrp();
    void test2DMeshNonConnexCut();
    void test1DMesh();
    void testInnerCrossCut();

   private:
    using Connections = std::vector<std::pair<mcIdType, mcIdType>>;
    static std::pair<bool, bool> TestCrack(
        MEDFileUMesh *mm_init, const std::string &grp_name, const std::string &test_name
    );

    static bool CheckM0Mesh(
        const MEDFileUMesh *mm, const Connections &c2cBrokenConnection, const Connections &c2cPreservedConnection
    );

    static bool CheckM1Mesh(const MEDCouplingUMesh *f2dup_before, const MEDCouplingUMesh *f2dup_after);

    static Connections GetC2CBroken(const MEDFileUMesh *mm, const std::string &grp_name);

    static std::vector<std::pair<mcIdType, mcIdType>> GetC2CPreserved(
        const MEDFileUMesh *mm, const std::string &grp_name
    );

    static MEDFileUMesh *make2x2Voxel();

    static MEDFileUMesh *make4x4Voxel();

    static MEDFileUMesh *make2DMesh();

    static MEDFileUMesh *make2DMesh2();

    static MEDFileUMesh *make1DMesh();
};
}  // namespace MEDCoupling
#endif  // __CRACKALGOTEST_HXX__
