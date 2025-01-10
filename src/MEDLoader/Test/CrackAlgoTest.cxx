// Copyright (C) 2007-2025  CEA, EDF
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

#include "CrackAlgoTest.hxx"

#include <string>
#include <algorithm>
#include <map>
#include <set>
#include <vector>
#include <utility>

#include <InterpKernelException.hxx>
#include <MEDFileMesh.hxx>
#include <MEDCouplingMemArray.hxx>
#include <MEDCouplingCMesh.hxx>
#include <MEDCouplingUMesh.hxx>
#include <MCAuto.hxx>
#include <MCType.hxx>

using std::vector;
using std::string;
using std::pair;
using std::copy;

using MEDCoupling::CrackAlgoTest;
using MEDCoupling::MEDFileUMesh;
using MEDCoupling::MEDCouplingCMesh;
using MEDCoupling::MEDCouplingUMesh;
using MEDCoupling::DataArrayIdType;
using MEDCoupling::MCAuto;

using MFU = MCAuto<MEDFileUMesh>;
using MCC = MCAuto<MEDCouplingCMesh>;
using MCU = MCAuto<MEDCouplingUMesh>;
using DAI = MCAuto<DataArrayIdType>;


/*
 *            TOP z=2              MIDDLE z=1           BOTTOM z=0
 *   2 +--------+--------+   +--------+--------+   +--------+--------+
 *     |        |        |   |        |        |   |        o        |
 *     |        |        |   |        |        |   |        o        |
 *     |        |        |   |        |        |   |       +-+       |
 *   1 +--------+--------+   +--------+--------+   +ooooooo|o|ooooooo+
 *     |        |        |   |        |        |   |       +-+       |
 *     |        |        |   |        |        |   |        o        |
 *     |        |        |   |        |        |   |        o        |
 *   0 +--------+--------+   +--------+--------+   +--------+--------+
 *
 *     0        1        2
 */
void CrackAlgoTest::test3DHalfCrossCut(
) {
    MFU mf = make2x2Voxel();

    const MCU m1 = mf->getMeshAtLevel(-1);
    double box[6] {1.0, 1.0, 1.0, 1.0, 0.1, 0.9};
    MCAuto<DataArrayIdType> g1 = m1->getCellsInBoundingBox(box, 0.01);
    g1->setName("M1");
    mf->setGroupsAtLevel(-1, vector<const DataArrayIdType*>{g1});

    CPPUNIT_ASSERT(g1->getNumberOfTuples() == 4);
    const auto res = TestCrack(mf, "M1", "halfCrossCut");
    CPPUNIT_ASSERT(res.first);
    CPPUNIT_ASSERT(res.second);
}

/*
 *            TOP z=2              MIDDLE z=1           BOTTOM z=0
 *   2 +--------+--------+   ooooooooooooooooooo   +--------+--------+
 *     |        o        |   ooooooooooooooooooo   |        o        |
 *     |        o        |   ooooooooooooooooooo   |        o        |
 *     |       +-+       |   oooooooo+-+oooooooo   |       +-+       |
 *   1 +ooooooo|o|ooooooo+   oooooooo|o|oooooooo   +ooooooo|o|ooooooo+
 *     |       +-+       |   oooooooo+-+oooooooo   |       +-+       |
 *     |        o        |   ooooooooooooooooooo   |        o        |
 *     |        o        |   ooooooooooooooooooo   |        o        |
 *   0 +--------+--------+   ooooooooooooooooooo   +--------+--------+
 *
 *     0        1        2
 */
void CrackAlgoTest::test3DFullFullCut(
) {
    MFU mf = make2x2Voxel();

    const MCU m1 = mf->getMeshAtLevel(-1);
    double box[6] {1.0, 1.0, 1.0, 1.0, 0.1, 1.9};
    DAI g1 = m1->getCellsInBoundingBox(box, 0.01);
    g1->setName("M1");
    mf->setGroupsAtLevel(-1, vector<const DataArrayIdType*>{g1});

    const auto res = TestCrack(mf, "M1", "fullFullCut");
    CPPUNIT_ASSERT(res.first);
    CPPUNIT_ASSERT(res.second);

    CPPUNIT_ASSERT(g1->getNumberOfTuples() == 12);
}

/*
 *            TOP z=2              MIDDLE z=1           BOTTOM z=0
 *   2 +--------+--------+   +--------+--------+   +--------+--------+
 *     |        o        |   |        o        |   |        o        |
 *     |        o        |   |        o        |   |        o        |
 *     |       +-+       |   |        o        |   |       +-+       |
 *   1 +ooooooo|o|ooooooo+   +ooooooooooooooooo+   +ooooooo|o|ooooooo+
 *     |       +-+       |   |        o        |   |       +-+       |
 *     |        o        |   |        o        |   |        o        |
 *     |        o        |   |        o        |   |        o        |
 *   0 +--------+--------+   +--------+--------+   +--------+--------+
 *
 *     0        1        2
 */
void CrackAlgoTest::test3DFullCrossCut(
) {
    MFU mf = make2x2Voxel();

    const MCU m1 = mf->getMeshAtLevel(-1);
    double box[6] {1.0, 1.0, 1.0, 1.0, 0.1, 0.9};
    MCAuto<DataArrayIdType> g1 = m1->getCellsInBoundingBox(box, 0.01);
    double box2[6] {1.0, 1.0, 1.0, 1.0, 1.1, 1.9};
    MCAuto<DataArrayIdType> g2 = m1->getCellsInBoundingBox(box, 0.01);
    g1->aggregate(g2);
    g1->setName("M1");
    mf->setGroupsAtLevel(-1, vector<const DataArrayIdType*>{g1});

    CPPUNIT_ASSERT(g1->getNumberOfTuples() == 8);

    const auto res = TestCrack(mf, "M1", "fullCrossCut");
    CPPUNIT_ASSERT(res.first);
    CPPUNIT_ASSERT(res.second);
}

/*
 *            TOP z=2              MIDDLE z=1           BOTTOM z=0
 *   2 +--------+--------+   +--------+--------+   +--------+--------+
 *     |        |        |   |        |        |   |        |        |
 *     |        |        |   |        |        |   |        |        |
 *     |        |        |   |        |        |   |  +-+   |   +-+  |
 *   1 +--------+--------+   +--------+--------+   +oo|o|ooooooo|o|oo+
 *     |        |        |   |        |        |   |  +-+   |   +-+  |
 *     |        |        |   |        |        |   |        |        |
 *     |        |        |   |        |        |   |        |        |
 *   0 +--------+--------+   +--------+--------+   +--------+--------+
 *
 *     0        1        2
 */
void CrackAlgoTest::test3DHalfCut(
) {
    MFU mf = make2x2Voxel();

    const MCU m1 = mf->getMeshAtLevel(-1);
    double box1[6] {0.5, 0.5, 1.0, 1.0, 0.1, 0.9};
    DAI g1 = m1->getCellsInBoundingBox(box1, 0.01);
    double box2[6] {1.5, 1.5, 1.0, 1.0, 0.1, 0.9};
    DAI g2 = m1->getCellsInBoundingBox(box2, 0.01);

    g1->aggregate(g2);

    g1->setName("M1");
    mf->setGroupsAtLevel(-1, vector<const DataArrayIdType*>{g1});

    CPPUNIT_ASSERT(g1->getNumberOfTuples() == 2);
    const auto res = TestCrack(mf, "M1", "halfCut");
    CPPUNIT_ASSERT(res.first);
    CPPUNIT_ASSERT(res.second);
}

/*
 *            TOP z=2              MIDDLE z=1           BOTTOM z=0
 *   2 +--------+--------+   +--------+--------+   +--------+--------+
 *     |        |        |   |        |        |   |        |        |
 *     |        |        |   |        |        |   |        |        |
 *     |  +-+   |   +-+  |   |        |        |   |  +-+   |   +-+  |
 *   1 +oo|o|ooooooo|o|oo+   +ooooooooooooooooo+   +oo|o|ooooooo|o|oo+
 *     |  +-+   |   +-+  |   |        |        |   |  +-+   |   +-+  |
 *     |        |        |   |        |        |   |        |        |
 *     |        |        |   |        |        |   |        |        |
 *   0 +--------+--------+   +--------+--------+   +--------+--------+
 *
 *     0        1        2
 */
void CrackAlgoTest::test3DFullCut(
) {
    MFU mf = make2x2Voxel();

    const MCU m1 = mf->getMeshAtLevel(-1);
    double box1[6] {0.5, 0.5, 1.0, 1.0, 0.1, 0.9};
    DAI g1 = m1->getCellsInBoundingBox(box1, 0.01);
    double box2[6] {1.5, 1.5, 1.0, 1.0, 0.1, 0.9};
    DAI g2 = m1->getCellsInBoundingBox(box2, 0.01);
    double box3[6] {0.5, 0.5, 1.0, 1.0, 1.1, 1.9};
    DAI g3 = m1->getCellsInBoundingBox(box1, 0.01);
    double box4[6] {1.5, 1.5, 1.0, 1.0, 1.1, 1.9};
    DAI g4 = m1->getCellsInBoundingBox(box2, 0.01);

    g1->aggregate(g2);
    g1->aggregate(g3);
    g1->aggregate(g4);

    g1->setName("M1");
    mf->setGroupsAtLevel(-1, vector<const DataArrayIdType*>{g1});

    CPPUNIT_ASSERT(g1->getNumberOfTuples() == 4);
    const auto res = TestCrack(mf, "M1", "fullCut");
    CPPUNIT_ASSERT(res.first);
    CPPUNIT_ASSERT(res.second);
}

/*
 *            TOP z=2              MIDDLE z=1           BOTTOM z=0
 *   2 +--------+--------+   +--------+--------+   +--------+--------+
 *     |        |        |   |        |        |   |        |        |
 *     |        |        |   |        |        |   |        |        |
 *     | +---+  |        |   |        |        |   |        |        |
 *   1 oo|ooo|oo+--------+   +--------+--------+   +--------+--------+
 *     | +---+  |        |   |        |        |   |        |        |
 *     |        |        |   |        |        |   |        |        |
 *     |        |        |   |        |        |   |        |        |
 *   0 +--------+--------+   +--------+--------+   +--------+--------+
 *
 *     0        1        2
 */
void CrackAlgoTest::test3DAngleCut(
) {
    MFU mf = make2x2Voxel();

    const MCU m1 = mf->getMeshAtLevel(-1);
    double box[6] {0.1, 0.9, 1.0, 1.0, 1.1, 1.9};
    MCAuto<DataArrayIdType> g1 = m1->getCellsInBoundingBox(box, 0.01);
    g1->setName("M1");
    mf->setGroupsAtLevel(-1, vector<const DataArrayIdType*>{g1});

    CPPUNIT_ASSERT(g1->getNumberOfTuples() == 1);
    const auto res = TestCrack(mf, "M1", "angleCut");
    CPPUNIT_ASSERT(res.first);
    CPPUNIT_ASSERT(res.second);
}

void CrackAlgoTest::testInnerCrossCut(
) {
    MFU mf = make4x4Voxel();

    const MCU m1 = mf->getMeshAtLevel(-1);
    double box[6] {2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
    MCAuto<DataArrayIdType> g1 = m1->getCellsInBoundingBox(box, 0.01);
    g1->setName("M1");
    mf->setGroupsAtLevel(-1, vector<const DataArrayIdType*>{g1});

    CPPUNIT_ASSERT(g1->getNumberOfTuples() == 12);
    const auto res = TestCrack(mf, "M1", "innerCrossCut");
    CPPUNIT_ASSERT(res.first);
    CPPUNIT_ASSERT(res.second);
}

void CrackAlgoTest::test2DGrp(
) {
    MCAuto<MEDFileUMesh> mf = make2DMesh();
    const auto res = TestCrack(mf, "Grp", "2dGrp");
    CPPUNIT_ASSERT(res.first);
    CPPUNIT_ASSERT(res.second);
}

void CrackAlgoTest::test1DMesh(
) {
    MFU mf = make1DMesh();

    const MCU m1 = mf->getMeshAtLevel(-1);
    double box[2] {1.9, 2.1};
    MCAuto<DataArrayIdType> g1 = m1->getCellsInBoundingBox(box, 0.01);
    g1->setName("M1");
    mf->setGroupsAtLevel(-1, vector<const DataArrayIdType*>{g1});

    CPPUNIT_ASSERT(g1->getNumberOfTuples() == 1);
    const auto res = TestCrack(mf, "M1", "1DMesh");
    CPPUNIT_ASSERT(res.first);
    CPPUNIT_ASSERT(res.second);
}

void CrackAlgoTest::test2DMeshNonConnexCut(
) {
    MCAuto<MEDFileUMesh> mf = make2DMesh2();
    const auto res = TestCrack(mf, "Grp", "2dGrpNonConnex");
    CPPUNIT_ASSERT(res.first);
    CPPUNIT_ASSERT(res.second);
}

MEDFileUMesh * CrackAlgoTest::make2x2Voxel(
) {
    double coords[3]  {0., 1., 2.};
    MCAuto<DataArrayDouble> x_coords = DataArrayDouble::New();
    x_coords->alloc(3);
    copy(coords, coords + 3, x_coords->getPointer());

    MCC m_cmesh = MEDCouplingCMesh::New();
    m_cmesh->setCoords(x_coords, x_coords, x_coords);
    MCU m0 = m_cmesh->buildUnstructured();
    m0->setName("UMesh");

    DAI desc(DataArrayIdType::New()),
        descIdx(DataArrayIdType::New()),
        revDesc(DataArrayIdType::New()),
        revDescIdx(DataArrayIdType::New());
    MCU m1 = m0->buildDescendingConnectivity(desc, descIdx, revDesc, revDescIdx);

    MFU mf = MEDFileUMesh::New();
    mf->setCoords(m0->getCoords());
    mf->setMeshAtLevel(0, m0);
    mf->setMeshAtLevel(-1, m1);

    return mf.retn();
}

MEDFileUMesh * CrackAlgoTest::make1DMesh(
) {
    double coords[5]  {0., 1., 2., 3., 4.};
    MCAuto<DataArrayDouble> x_coords = DataArrayDouble::New();
    x_coords->alloc(5);
    copy(coords, coords + 5, x_coords->getPointer());

    MCC m_cmesh = MEDCouplingCMesh::New();
    m_cmesh->setCoords(x_coords);
    MCU m0 = m_cmesh->buildUnstructured();
    m0->setName("UMesh");

    DAI desc(DataArrayIdType::New()),
        descIdx(DataArrayIdType::New()),
        revDesc(DataArrayIdType::New()),
        revDescIdx(DataArrayIdType::New());
    MCU m1 = m0->buildDescendingConnectivity(desc, descIdx, revDesc, revDescIdx);

    MFU mf = MEDFileUMesh::New();
    mf->setCoords(m0->getCoords());
    mf->setMeshAtLevel(0, m0);
    mf->setMeshAtLevel(-1, m1);

    return mf.retn();
}

MEDFileUMesh * CrackAlgoTest::make4x4Voxel(
) {
    double coords[5]  {0., 1., 2., 3., 4.};
    MCAuto<DataArrayDouble> x_coords = DataArrayDouble::New();
    x_coords->alloc(5);
    copy(coords, coords + 5, x_coords->getPointer());

    MCAuto<MEDCouplingCMesh> m_cmesh = MEDCouplingCMesh::New();
    m_cmesh->setCoords(x_coords, x_coords, x_coords);
    MCAuto<MEDCouplingUMesh> m0 = m_cmesh->buildUnstructured();
    m0->setName("UMesh");

    DAI desc(DataArrayIdType::New()), descIdx(DataArrayIdType::New()), revDesc(DataArrayIdType::New()), revDescIdx(DataArrayIdType::New());
    MCU m1 = m0->buildDescendingConnectivity(desc, descIdx, revDesc, revDescIdx);

    MEDFileUMesh * mf = MEDFileUMesh::New();
    mf->setCoords(m0->getCoords());
    mf->setMeshAtLevel(0, m0);
    mf->setMeshAtLevel(-1, m1);

    return mf;
}

MEDFileUMesh * CrackAlgoTest::make2DMesh(
) {
    using DAD = MCAuto<DataArrayDouble>;

    double c0[6] {0.0, 1.1, 2.3, 3.6, 5.0, 6.5};
    DAD coords0 = DataArrayDouble::New();
    coords0->alloc(6);
    copy(c0, c0 + 6, coords0->rwBegin());

    double c1[5] {0.0, 1.1, 2.3, 3.6, 5.0};
    DAD coords1 = DataArrayDouble::New();
    coords1->alloc(5);
    copy(c1, c1+ 5, coords1->rwBegin());

    MCC m = MEDCouplingCMesh::New();
    m->setCoordsAt(0, coords0);
    m->setCoordsAt(1, coords1);

    MCU m0 = m->buildUnstructured();
    m0->setName("duplicate");

    DAI desc(DataArrayIdType::New()),
        descIdx(DataArrayIdType::New()),
        revDesc(DataArrayIdType::New()),
        revDescIdx(DataArrayIdType::New());
    MCU m2 = m0->buildDescendingConnectivity(desc, descIdx, revDesc, revDescIdx);
    m2->setName(m0->getName());

    mcIdType ids_[17] {8,11,14,20,21,22,23,24,25,26,31,32,33,34,35,36,37};
    MCU m3 = m2->buildPartOfMySelf(ids_, ids_+17);

    mcIdType grp_[3] {4, 6, 8};
    DAI grp = DataArrayIdType::New();
    grp->alloc(3);
    copy(grp_, grp_+3, grp->rwBegin());
    grp->setName("Grp");

    mcIdType grp2_[2] {9, 16};
    DAI grp2 = DataArrayIdType::New();
    grp2->alloc(2);
    copy(grp2_, grp2_+2, grp2->rwBegin());
    grp2->setName("Grp2");

    MEDFileUMesh * mm = MEDFileUMesh::New();
    mm->setMeshAtLevel(0, m0);
    mm->setMeshAtLevel(-1, m3);
    mm->setGroupsAtLevel(-1, {grp, grp2});

    mcIdType grpNode_[3] {4, 21, 23};
    DAI grpNode = DataArrayIdType::New();
    grpNode->alloc(3);
    copy(grpNode_, grpNode_+3, grpNode->rwBegin());
    grpNode->setName("GrpNode");
    mm->setGroupsAtLevel(1, {grpNode});

    return mm;
}

MEDFileUMesh * CrackAlgoTest::make2DMesh2(
) {
    using DAD = MCAuto<DataArrayDouble>;

    double c0[5] {0.0, 1.1, 2.3, 3.6, 5.0};
    DAD coords0 = DataArrayDouble::New();
    coords0->alloc(5);
    copy(c0, c0 + 5, coords0->rwBegin());

    double c1[3] {0.0, 1.0, 2.0};
    DAD coords1 = DataArrayDouble::New();
    coords1->alloc(3);
    copy(c1, c1 + 3, coords1->rwBegin());

    MCC m = MEDCouplingCMesh::New();
    m->setCoordsAt(0, coords0);
    m->setCoordsAt(1, coords1);

    MCU m0 = m->buildUnstructured();
    m0->setName("simple");

    DAI desc(DataArrayIdType::New()),
        descIdx(DataArrayIdType::New()),
        revDesc(DataArrayIdType::New()),
        revDescIdx(DataArrayIdType::New());
    MCU m2 = m0->buildDescendingConnectivity(desc, descIdx, revDesc, revDescIdx);
    m2->setName(m0->getName());

    mcIdType grp_[2] {3, 19};
    DAI grp = DataArrayIdType::New();
    grp->alloc(2);
    copy(grp_, grp_+2, grp->rwBegin());
    grp->setName("Grp");

    MEDFileUMesh * mm = MEDFileUMesh::New();
    mm->setMeshAtLevel(0, m0);
    mm->setMeshAtLevel(-1, m2);
    mm->setGroupsAtLevel(-1, {grp});

    return mm;
}


pair<bool, bool> CrackAlgoTest::TestCrack(
    MEDFileUMesh * mm_init,
    const string & grp_name,
    const string & test_name
) {
    const auto c2cBroken = GetC2CBroken(mm_init, grp_name);
    const auto c2cPreserved = GetC2CPreserved(mm_init, grp_name);

    const MCU f2dup = mm_init->getGroup(-1, grp_name);

    mm_init->write(test_name + "_in.med", 2);
    const auto cellOld2NewNode = mm_init->crackAlong(grp_name);
    mm_init->write(test_name + "_out.med", 2);

    const MCU f2dup_b = mm_init->getGroup(-1, grp_name);

    const auto res = pair<bool, bool> {
        CheckM0Mesh(mm_init, c2cBroken, c2cPreserved),
        CheckM1Mesh(f2dup, f2dup_b)
        };

    mm_init->openCrack(cellOld2NewNode, 0.9);
    mm_init->write(test_name + "_cracked.med", 2);

    return res;
}

bool CrackAlgoTest::CheckM1Mesh(
    const MEDCouplingUMesh * f2dup_before,
    const MEDCouplingUMesh * f2dup_after
) {
    const mcIdType nbFaces_0 = f2dup_before->getNumberOfCells();
    const mcIdType nbFaces_1 = f2dup_after->getNumberOfCells();

    MCU f2dup_before_copy = f2dup_after->deepCopy();
    f2dup_before_copy->zipCoords();

    MCU f2dup_after_copy = f2dup_after->deepCopy();
    bool nodesAreMerged = false;
    mcIdType newNbOfNodes = 0;
    DAI o2n = f2dup_after_copy->mergeNodes(1e-12, nodesAreMerged, newNbOfNodes);
    f2dup_after_copy->zipCoords();

    DataArrayIdType * cellCor;
    DataArrayIdType * nodeCor;
    f2dup_after_copy->checkGeoEquivalWith(f2dup_before_copy, 12, 1e-12, cellCor, nodeCor);
    DAI cellCorrI(cellCor);
    DAI nodeCorrI(nodeCor);

    return (2 * nbFaces_0 == nbFaces_1);
}

vector<pair<mcIdType, mcIdType>>
CrackAlgoTest::GetC2CBroken(
    const MEDFileUMesh * mm,
    const string & grp_name
) {
    vector<pair<mcIdType, mcIdType>> res;
    MCU m0 = mm->getMeshAtLevel(0);
    DAI desc(DataArrayIdType::New()),
        descIdx(DataArrayIdType::New()),
        revDesc(DataArrayIdType::New()),
        revDescIdx(DataArrayIdType::New());
    MCU mf = m0->buildDescendingConnectivity(desc, descIdx, revDesc, revDescIdx);
    MCU f2dupMesh = mm->getGroup(-1, grp_name);
    DataArrayIdType * f2dup_p;
    mf->areCellsIncludedIn(f2dupMesh, 2, f2dup_p);
    DAI f2dup(f2dup_p);

    for (const auto & f : *f2dup) {
        const auto & cellsId_start = revDescIdx->begin()[f];
        const auto & cellsId_end = revDescIdx->begin()[f+1];
        if (cellsId_end - cellsId_start != 2)
            throw INTERP_KERNEL::Exception(
                "crackAlong: A face of group M1 does not have two neighbors.");
        const auto & cell0 = revDesc->begin()[cellsId_start];
        const auto & cell1 = revDesc->begin()[cellsId_start + 1];
        res.push_back(pair<mcIdType, mcIdType>(cell0, cell1));
    }
    return res;
}

vector<pair<mcIdType, mcIdType>>
CrackAlgoTest::GetC2CPreserved(
    const MEDFileUMesh * mm,
    const string & grp_name
) {
    vector<pair<mcIdType, mcIdType>> res;
    MCU m0 = mm->getMeshAtLevel(0);
    DAI desc(DataArrayIdType::New()),
        descIdx(DataArrayIdType::New()),
        revDesc(DataArrayIdType::New()),
        revDescIdx(DataArrayIdType::New());
    MCU mf = m0->buildDescendingConnectivity(desc, descIdx, revDesc, revDescIdx);
    MCU f2dupMesh = mm->getGroup(-1, grp_name);

    DataArrayIdType * f2dup_p;
    mf->areCellsIncludedIn(f2dupMesh, 2, f2dup_p);
    DAI f2dup(f2dup_p);

    mcIdType nbOfFaces = mf->getNumberOfCells();
    for (mcIdType f = 0; f < nbOfFaces; f++) {
        if (std::find(f2dup->begin(), f2dup->end(), f) != f2dup->end())
            continue;

        const auto & cellsId_start = revDescIdx->begin()[f];
        const auto & cellsId_end = revDescIdx->begin()[f+1];
        if (cellsId_end - cellsId_start == 1)
            continue;
        if (cellsId_end - cellsId_start > 2)
            throw INTERP_KERNEL::Exception(
                "crackAlong: A face of the mf has more than two neighbors.");

        const auto & cell0 = revDesc->begin()[cellsId_start];
        const auto & cell1 = revDesc->begin()[cellsId_start + 1];
        res.push_back(pair<mcIdType, mcIdType>(cell0, cell1));
    }
    return res;
}

bool
CrackAlgoTest::CheckM0Mesh(
    const MEDFileUMesh * mm,
    const Connections& c2cBrokenConnection,
    const Connections& c2cPreservedConnection
) {
    MCU m0 = mm->getMeshAtLevel(0);
    DAI desc(DataArrayIdType::New()), descIdx(DataArrayIdType::New()), revDesc(DataArrayIdType::New()), revDescIdx(DataArrayIdType::New());
    MCU mf = m0->buildDescendingConnectivity(desc, descIdx, revDesc, revDescIdx);

    std::map<mcIdType, std::set<mcIdType>> c2c;

    const mcIdType * descIdx_ptr = descIdx->begin();
    const mcIdType * desc_ptr = desc->begin();
    const mcIdType * revDescIdx_ptr = revDescIdx->begin();
    const mcIdType * revDesc_ptr = revDesc->begin();
    for (mcIdType cell = 0; cell < m0->getNumberOfCells(); cell++) {
        auto & neighbors = c2c[cell];
        for (
            mcIdType descId = descIdx_ptr[cell];
            descId < descIdx_ptr[cell+1];
            descId++
        ) {
            const auto & face = desc_ptr[descId];
            for (
                mcIdType revDescId = revDescIdx_ptr[face];
                revDescId < revDescIdx_ptr[face+1];
                revDescId++
            ) {
                const auto & otherCell = revDesc_ptr[revDescId];
                if (otherCell != cell)
                    neighbors.insert(otherCell);
            }
        }
    }

    bool res = true;
    for (const auto & conn : c2cBrokenConnection) {
        const auto & cell0 = conn.first;
        const auto & cell1 = conn.second;
        // NOTE: all cells are in c2c so cell0 and cell1 are in c2c
        if (c2c.at(cell0).find(cell1) != c2c.at(cell0).end())
            res = false;
        if (c2c.at(cell1).find(cell0) != c2c.at(cell1).end())
            res = false;
    }

    for (const auto & conn : c2cPreservedConnection) {
        const auto & cell0 = conn.first;
        const auto & cell1 = conn.second;
        // NOTE: all cells are in c2c so cell0 and cell1 are in c2c
        if (c2c.at(cell0).find(cell1) == c2c.at(cell0).end())
            res = false;
        if (c2c.at(cell1).find(cell0) == c2c.at(cell1).end())
            res = false;
    }
    return res;
}
