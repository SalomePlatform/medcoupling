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

#include "ParaMEDMEMTest.hxx"
#include "MEDLoader.hxx"

#include "ParaMEDFileMesh.hxx"
#include "MEDFileMesh.hxx"
#include "MEDFileField1TS.hxx"
#include "TestInterpKernelUtils.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include <cppunit/TestAssert.h>

#include <algorithm>
#include <numeric>
#include <iostream>
#include <iterator>

using namespace MEDCoupling;

/*
 * Generate a 2D mesh that is supposed to match the part that will be loaded by each proc in testParallelLoad1
 */
MEDCouplingUMesh *
genLocMesh2D(int rk)
{
    int nxTot = 4, nyTot = 2;
    int nx = 2, ny = 2;
    MCAuto<MEDCouplingCMesh> msh = MEDCouplingCMesh::New("mesh");
    MCAuto<DataArrayDouble> dax = DataArrayDouble::New();
    dax->alloc(nx + 1, 1);
    MCAuto<DataArrayDouble> day = DataArrayDouble::New();
    day->alloc(ny + 1, 1);
    dax->iota();
    day->iota();
    if (rk == 0)
    {
        std::transform(dax->begin(), dax->end(), dax->rwBegin(), [nxTot](const int &c) { return c / (float)nxTot; });
        std::transform(day->begin(), day->end(), day->rwBegin(), [nyTot](const int &c) { return c / (float)nyTot; });
    }
    else
    {
        std::transform(
            dax->begin(), dax->end(), dax->rwBegin(), [nxTot](const int &c) { return c / (float)nxTot + 0.5; }
        );
        std::transform(day->begin(), day->end(), day->rwBegin(), [nyTot](const int &c) { return c / (float)nyTot; });
    }
    msh->setCoords(dax, day);
    MCAuto<MEDCouplingUMesh> ret = msh->buildUnstructured();
    return ret.retn();
}

/*
 * Generate a 2D mesh that is supposed to match the part that will be loaded by proc0 in testParallelLoad2
 */
MEDCouplingUMesh *
genLocMeshMultipleTypes1()
{
    MCAuto<MEDCouplingUMesh> ret = MEDCouplingUMesh::New("mesh", 2);
    double coords[10] = {0., 1., 0., 2., 1., 2., 0., 3., 1., 3.};
    DataArrayDouble *myCoords = DataArrayDouble::New();
    myCoords->alloc(5, 2);
    std::copy(coords, coords + 10, myCoords->getPointer());
    ret->setCoords(myCoords);
    myCoords->decrRef();
    mcIdType conn[7] = {0, 2, 1, 1, 2, 4, 3};
    ret->allocateCells(2);
    ret->insertNextCell(INTERP_KERNEL::NORM_TRI3, 3, conn);
    ret->insertNextCell(INTERP_KERNEL::NORM_QUAD4, 4, conn + 3);
    ret->finishInsertingCells();
    return ret.retn();
}

/*
 * Generate a 2D mesh that is supposed to match the part that will be loaded by proc1 in testParallelLoad2
 */
MEDCouplingUMesh *
genLocMeshMultipleTypes2()
{
    MCAuto<MEDCouplingUMesh> ret = MEDCouplingUMesh::New("mesh", 2);
    double coords[10] = {0., 0., 1., 0., 0., 1., 1., 1., 1., 2.};
    DataArrayDouble *myCoords = DataArrayDouble::New();
    myCoords->alloc(5, 2);
    std::copy(coords, coords + 10, myCoords->getPointer());
    ret->setCoords(myCoords);
    myCoords->decrRef();
    mcIdType conn[7] = {2, 3, 4, 0, 1, 3, 2};
    ret->allocateCells(2);
    ret->insertNextCell(INTERP_KERNEL::NORM_TRI3, 3, conn);
    ret->insertNextCell(INTERP_KERNEL::NORM_QUAD4, 4, conn + 3);
    ret->finishInsertingCells();
    return ret.retn();
}

/*
 * Generate a 2D mesh that is supposed to match the part that will be loaded by proc2 in testParallelLoad2
 */
MEDCouplingUMesh *
genLocMeshMultipleTypes3()
{
    MCAuto<MEDCouplingUMesh> ret = MEDCouplingUMesh::New("mesh", 2);
    double coords[16] = {1., 0., 2., 0., 1., 1., 2., 1., 1., 2., 2., 2., 1., 3., 2., 3.};
    DataArrayDouble *myCoords = DataArrayDouble::New();
    myCoords->alloc(8, 2);
    std::copy(coords, coords + 16, myCoords->getPointer());
    ret->setCoords(myCoords);
    myCoords->decrRef();
    mcIdType conn[14] = {0, 1, 3, 0, 3, 2, 2, 3, 5, 4, 4, 5, 7, 6};
    ret->allocateCells(4);
    ret->insertNextCell(INTERP_KERNEL::NORM_TRI3, 3, conn);
    ret->insertNextCell(INTERP_KERNEL::NORM_TRI3, 3, conn + 3);
    ret->insertNextCell(INTERP_KERNEL::NORM_QUAD4, 4, conn + 6);
    ret->insertNextCell(INTERP_KERNEL::NORM_QUAD4, 4, conn + 10);
    ret->finishInsertingCells();
    return ret.retn();
}

/*
 * Generate a 2D mesh that is supposed to match the part that will be loaded by proc0 in testParallelLoad6
 */
MEDCouplingUMesh *
genPartialLocMeshMultipleTypes1()
{
    MCAuto<MEDCouplingUMesh> ret = MEDCouplingUMesh::New("mesh", 2);
    double coords[8] = {0., 2., 1., 2., 0., 3., 1., 3.};
    DataArrayDouble *myCoords = DataArrayDouble::New();
    myCoords->alloc(4, 2);
    std::copy(coords, coords + 8, myCoords->getPointer());
    ret->setCoords(myCoords);
    myCoords->decrRef();
    mcIdType conn[4] = {0, 1, 3, 2};
    ret->allocateCells(1);
    ret->insertNextCell(INTERP_KERNEL::NORM_QUAD4, 4, conn);
    ret->finishInsertingCells();
    return ret.retn();
}

/*
 * Generate a 2D mesh that is supposed to match the part that will be loaded by proc1 in testParallelLoad6
 */
MEDCouplingUMesh *
genPartialLocMeshMultipleTypes2()
{
    MCAuto<MEDCouplingUMesh> ret = MEDCouplingUMesh::New("mesh", 2);
    double coords[6] = {0., 1., 1., 1., 1., 2.};
    DataArrayDouble *myCoords = DataArrayDouble::New();
    myCoords->alloc(3, 2);
    std::copy(coords, coords + 6, myCoords->getPointer());
    ret->setCoords(myCoords);
    myCoords->decrRef();
    mcIdType conn[3] = {0, 1, 2};
    ret->allocateCells(1);
    ret->insertNextCell(INTERP_KERNEL::NORM_TRI3, 3, conn);
    ret->finishInsertingCells();
    return ret.retn();
}

/*
 * Generate a 2D mesh that is supposed to match the part that will be loaded by proc2 in testParallelLoad6
 */
MEDCouplingUMesh *
genPartialLocMeshMultipleTypes3()
{
    MCAuto<MEDCouplingUMesh> ret = MEDCouplingUMesh::New("mesh", 2);
    double coords[12] = {1., 0., 2., 0., 1., 1., 2., 1., 1., 2., 2., 2.};
    DataArrayDouble *myCoords = DataArrayDouble::New();
    myCoords->alloc(6, 2);
    std::copy(coords, coords + 12, myCoords->getPointer());
    ret->setCoords(myCoords);
    myCoords->decrRef();
    mcIdType conn[7] = {0, 1, 3, 2, 3, 5, 4};
    ret->allocateCells(2);
    ret->insertNextCell(INTERP_KERNEL::NORM_TRI3, 3, conn);
    ret->insertNextCell(INTERP_KERNEL::NORM_QUAD4, 4, conn + 3);
    ret->finishInsertingCells();
    return ret.retn();
}

/*
 * Generate a 2D field that is supposed to match the local field loaded by each proc in testParallelLoad4
 */
MEDCouplingFieldDouble *
genLocFieldCells(int rank)
{
    MCAuto<MEDCouplingUMesh> mesh = genLocMesh2D(rank);
    MCAuto<MEDCouplingFieldDouble> f1 = MEDCouplingFieldDouble::New(ON_CELLS, ONE_TIME);
    f1->setName("field");
    f1->setMesh(mesh);

    MCAuto<DataArrayDouble> array(DataArrayDouble::New());
    array->alloc(4, 2);
    std::vector<double> values;
    if (rank == 0)
        values = {0., 10., 20., 30., 80., 90., 100., 110.};
    else
        values = {40., 50., 60., 70., 120., 130., 140., 150.};
    std::copy(values.data(), values.data() + 8, array->getPointer());
    array->setInfoOnComponent(0, "");
    f1->setArray(array);
    return f1.retn();
}

/*
 * Generate a 2D field that is supposed to match the local field loaded by each proc in testParallelLoad5
 */
MEDCouplingFieldDouble *
genLocFieldNodes(int rank)
{
    MCAuto<MEDCouplingUMesh> mesh = genLocMesh2D(rank);
    MCAuto<MEDCouplingFieldDouble> f1 = MEDCouplingFieldDouble::New(ON_NODES, ONE_TIME);
    f1->setName("field");
    f1->setMesh(mesh);

    MCAuto<DataArrayDouble> array(DataArrayDouble::New());
    array->alloc(9, 2);
    std::vector<double> values;
    if (rank == 0)
        values = {0., 10., 20., 30., 40., 50., 100., 110., 120., 130., 140., 150., 200., 210., 220., 230., 240., 250.};
    else
        values = {40., 50., 60., 70., 80., 90., 140., 150., 160., 170., 180., 190., 240., 250., 260., 270., 280., 290.};
    std::copy(values.data(), values.data() + 18, array->getPointer());
    array->setInfoOnComponent(0, "");
    f1->setArray(array);
    return f1.retn();
}

/*!
 * Test case to load a simple 2D cartesian mesh in parallel on 2 procs
 */
void
ParaMEDMEMTest::testParallelLoad1()
{
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //
    if (size != 2)
        return;

    std::map<INTERP_KERNEL::NormalizedCellType, std::vector<mcIdType>> distrib;
    if (rank == 0)
        distrib = {{INTERP_KERNEL::NORM_QUAD4, {0, 1, 4, 5} /*c++ type of indexing: index starts from zero!*/}};
    else
        distrib = {{INTERP_KERNEL::NORM_QUAD4, {2, 3, 6, 7}}};

    std::string filename = INTERP_TEST::getResourceFile("SimpleTest2D.med");
    MCAuto<MEDFileUMesh> mu = ParaMEDFileUMesh::ParaNew(distrib, MPI_COMM_WORLD, MPI_INFO_NULL, filename, "mesh");
    MCAuto<MEDCouplingUMesh> mesh = mu->getMeshAtLevel(0);
    MCAuto<MEDCouplingUMesh> meshRef = genLocMesh2D(rank);
    CPPUNIT_ASSERT(mesh->isEqual(meshRef, 1e-12));
    MPI_Barrier(MPI_COMM_WORLD);
}

/*!
 * Test case to load a 2D mesh made of squares and triangles in parallel on 3 procs.
 * Each proc is going to load a part of the mesh.
 */
void
ParaMEDMEMTest::testParallelLoad2()
{
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //
    if (size != 3)
        return;

    std::map<INTERP_KERNEL::NormalizedCellType, std::vector<mcIdType>> distrib;
    // independant numerotation for each geometric type!
    if (rank == 0)
        distrib = {{INTERP_KERNEL::NORM_TRI3, {3}}, {INTERP_KERNEL::NORM_QUAD4, {2}}};
    else if (rank == 1)
        distrib = {{INTERP_KERNEL::NORM_TRI3, {2}}, {INTERP_KERNEL::NORM_QUAD4, {0}}};
    else
        distrib = {{INTERP_KERNEL::NORM_TRI3, {0, 1}}, {INTERP_KERNEL::NORM_QUAD4, {1, 3}}};

    std::string filename = INTERP_TEST::getResourceFile("Test2DMultiGeoType.med");
    // partial loading
    MCAuto<MEDFileUMesh> mu = ParaMEDFileUMesh::ParaNew(distrib, MPI_COMM_WORLD, MPI_INFO_NULL, filename, "mesh");
    MCAuto<MEDCouplingUMesh> mesh = mu->getMeshAtLevel(0);

    MEDCouplingUMesh *meshRef;
    if (rank == 0)
        meshRef = genLocMeshMultipleTypes1();
    else if (rank == 1)
        meshRef = genLocMeshMultipleTypes2();
    else
        meshRef = genLocMeshMultipleTypes3();
    // checking that all 3 procs have correctly loaded their part
    int equal = (int)mesh->isEqual(meshRef, 1e-12);
    int allEqual = -1;
    MPI_Allreduce(&equal, &allEqual, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CPPUNIT_ASSERT(allEqual == 3);
    meshRef->decrRef();

    MPI_Barrier(MPI_COMM_WORLD);
}

/*!
 * Test case to load a 3D box meshed with tetras in parallel on 2 procs
 */
void
ParaMEDMEMTest::testParallelLoad3()
{
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //
    if (size != 2)
        return;

    std::map<INTERP_KERNEL::NormalizedCellType, std::vector<mcIdType>> distrib;
    if (rank == 0)
    {
        std::vector<mcIdType> distribCells = {0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,
                                              14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  48,  49,  50,  51,
                                              52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  65,
                                              66,  67,  68,  69,  70,  71,  96,  97,  98,  99,  100, 101, 102, 103,
                                              104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117,
                                              118, 119, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155,
                                              156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167};
        distrib = {{INTERP_KERNEL::NORM_TETRA4, distribCells}};
    }
    else
    {
        std::vector<mcIdType> distribCells = {24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,
                                              38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  72,  73,  74,  75,
                                              76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,
                                              90,  91,  92,  93,  94,  95,  120, 121, 122, 123, 124, 125, 126, 127,
                                              128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141,
                                              142, 143, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179,
                                              180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191};
        distrib = {{INTERP_KERNEL::NORM_TETRA4, distribCells}};
    }

    std::string filename = INTERP_TEST::getResourceFile("SimpleTest3D.med");
    MCAuto<MEDFileUMesh> mu = ParaMEDFileUMesh::ParaNew(distrib, MPI_COMM_WORLD, MPI_INFO_NULL, filename, "mesh");
    MCAuto<MEDCouplingUMesh> mesh = mu->getMeshAtLevel(0);
    CPPUNIT_ASSERT_EQUAL(96, (int)mesh->getNumberOfCells());

    // checking nodal connectivity
    double nodalConnec[480] = {
        14, 1,  7,  18, 24, 14, 7,  6,  18, 24, 14, 6,  0,  18, 24, 14, 0,  1,  18, 24, 14, 1,  0,  19, 24, 14, 0,
        2,  19, 24, 14, 2,  3,  19, 24, 14, 3,  1,  19, 24, 14, 1,  3,  20, 24, 14, 3,  9,  20, 24, 14, 9,  7,  20,
        24, 14, 7,  1,  20, 24, 14, 0,  6,  21, 24, 14, 6,  8,  21, 24, 14, 8,  2,  21, 24, 14, 2,  0,  21, 24, 14,
        7,  9,  22, 24, 14, 9,  8,  22, 24, 14, 8,  6,  22, 24, 14, 6,  7,  22, 24, 14, 2,  8,  23, 24, 14, 8,  9,
        23, 24, 14, 9,  3,  23, 24, 14, 3,  2,  23, 24, 14, 3,  9,  25, 31, 14, 9,  8,  25, 31, 14, 8,  2,  25, 31,
        14, 2,  3,  25, 31, 14, 3,  2,  26, 31, 14, 2,  4,  26, 31, 14, 4,  5,  26, 31, 14, 5,  3,  26, 31, 14, 3,
        5,  27, 31, 14, 5,  11, 27, 31, 14, 11, 9,  27, 31, 14, 9,  3,  27, 31, 14, 2,  8,  28, 31, 14, 8,  10, 28,
        31, 14, 10, 4,  28, 31, 14, 4,  2,  28, 31, 14, 9,  11, 29, 31, 14, 11, 10, 29, 31, 14, 10, 8,  29, 31, 14,
        8,  9,  29, 31, 14, 4,  10, 30, 31, 14, 10, 11, 30, 31, 14, 11, 5,  30, 31, 14, 5,  4,  30, 31, 14, 7,  13,
        32, 38, 14, 13, 12, 32, 38, 14, 12, 6,  32, 38, 14, 6,  7,  32, 38, 14, 7,  6,  33, 38, 14, 6,  8,  33, 38,
        14, 8,  9,  33, 38, 14, 9,  7,  33, 38, 14, 7,  9,  34, 38, 14, 9,  15, 34, 38, 14, 15, 13, 34, 38, 14, 13,
        7,  34, 38, 14, 6,  12, 35, 38, 14, 12, 14, 35, 38, 14, 14, 8,  35, 38, 14, 8,  6,  35, 38, 14, 13, 15, 36,
        38, 14, 15, 14, 36, 38, 14, 14, 12, 36, 38, 14, 12, 13, 36, 38, 14, 8,  14, 37, 38, 14, 14, 15, 37, 38, 14,
        15, 9,  37, 38, 14, 9,  8,  37, 38, 14, 9,  15, 39, 45, 14, 15, 14, 39, 45, 14, 14, 8,  39, 45, 14, 8,  9,
        39, 45, 14, 9,  8,  40, 45, 14, 8,  10, 40, 45, 14, 10, 11, 40, 45, 14, 11, 9,  40, 45, 14, 9,  11, 41, 45,
        14, 11, 17, 41, 45, 14, 17, 15, 41, 45, 14, 15, 9,  41, 45, 14, 8,  14, 42, 45, 14, 14, 16, 42, 45, 14, 16,
        10, 42, 45, 14, 10, 8,  42, 45, 14, 15, 17, 43, 45, 14, 17, 16, 43, 45, 14, 16, 14, 43, 45, 14, 14, 15, 43,
        45, 14, 10, 16, 44, 45, 14, 16, 17, 44, 45, 14, 17, 11, 44, 45, 14, 11, 10, 44, 45
    };
    const mcIdType *nc = mesh->getNodalConnectivity()->getConstPointer();
    CPPUNIT_ASSERT_EQUAL(480, (int)mesh->getNodalConnectivity()->getNumberOfTuples());
    CPPUNIT_ASSERT(std::equal(nodalConnec, nodalConnec + 480, nc));

    double nodalConnecInd[97] = {0,   5,   10,  15,  20,  25,  30,  35,  40,  45,  50,  55,  60,  65,  70,  75,  80,
                                 85,  90,  95,  100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165,
                                 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250,
                                 255, 260, 265, 270, 275, 280, 285, 290, 295, 300, 305, 310, 315, 320, 325, 330, 335,
                                 340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 390, 395, 400, 405, 410, 415, 420,
                                 425, 430, 435, 440, 445, 450, 455, 460, 465, 470, 475, 480};
    const mcIdType *ncIndx = mesh->getNodalConnectivityIndex()->getConstPointer();
    CPPUNIT_ASSERT_EQUAL(97, (int)mesh->getNodalConnectivityIndex()->getNumberOfTuples());
    CPPUNIT_ASSERT(std::equal(nodalConnecInd, nodalConnecInd + 97, ncIndx));

    // checking coords
    std::vector<double> coords(138);
    if (rank == 0)
        coords = {0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0.0, 4.0, 0.0, 2.0, 4.0, 0.0, 0.0, 0.0,
                  2.0, 2.0, 0.0, 2.0, 0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 4.0, 2.0, 2.0, 4.0, 2.0, 0.0, 0.0, 4.0, 2.0,
                  0.0, 4.0, 0.0, 2.0, 4.0, 2.0, 2.0, 4.0, 0.0, 4.0, 4.0, 2.0, 4.0, 4.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0,
                  2.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 3.0,
                  0.0, 2.0, 3.0, 1.0, 0.0, 3.0, 1.0, 1.0, 3.0, 2.0, 1.0, 4.0, 1.0, 1.0, 3.0, 1.0, 1.0, 0.0, 3.0, 1.0,
                  1.0, 2.0, 2.0, 1.0, 3.0, 0.0, 1.0, 3.0, 1.0, 1.0, 4.0, 1.0, 2.0, 3.0, 1.0, 1.0, 3.0, 1.0, 2.0, 3.0,
                  1.0, 3.0, 2.0, 2.0, 3.0, 3.0, 0.0, 3.0, 3.0, 1.0, 3.0, 4.0, 1.0, 4.0, 3.0, 1.0, 3.0, 3.0};
    else
        coords = {2.0, 0.0, 0.0, 4.0, 0.0, 0.0, 2.0, 2.0, 0.0, 4.0, 2.0, 0.0, 2.0, 4.0, 0.0, 4.0, 4.0, 0.0, 2.0, 0.0,
                  2.0, 4.0, 0.0, 2.0, 2.0, 2.0, 2.0, 4.0, 2.0, 2.0, 2.0, 4.0, 2.0, 4.0, 4.0, 2.0, 2.0, 0.0, 4.0, 4.0,
                  0.0, 4.0, 2.0, 2.0, 4.0, 4.0, 2.0, 4.0, 2.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.0, 0.0, 1.0, 3.0, 1.0, 0.0,
                  4.0, 1.0, 1.0, 2.0, 1.0, 1.0, 3.0, 1.0, 2.0, 3.0, 2.0, 1.0, 3.0, 1.0, 1.0, 3.0, 2.0, 1.0, 3.0, 3.0,
                  0.0, 4.0, 3.0, 1.0, 2.0, 3.0, 1.0, 3.0, 3.0, 2.0, 3.0, 4.0, 1.0, 3.0, 3.0, 1.0, 3.0, 0.0, 3.0, 3.0,
                  1.0, 2.0, 4.0, 1.0, 3.0, 2.0, 1.0, 3.0, 3.0, 1.0, 4.0, 3.0, 2.0, 3.0, 3.0, 1.0, 3.0, 3.0, 2.0, 3.0,
                  3.0, 3.0, 2.0, 4.0, 3.0, 3.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, 3.0, 4.0, 3.0, 3.0, 3.0, 3.0};
    const double *coo = mesh->getCoords()->getConstPointer();
    CPPUNIT_ASSERT_EQUAL(46, (int)mesh->getCoords()->getNumberOfTuples());
    CPPUNIT_ASSERT(std::equal(coords.data(), coords.data() + 138, coo));

    MPI_Barrier(MPI_COMM_WORLD);
}

/*!
 * Test case to load a field located on cells in parallel on 2 procs.
 */
void
ParaMEDMEMTest::testParallelLoad4()
{
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //
    if (size != 2)
        return;

    std::vector<mcIdType> distrib;
    if (rank == 0)
        distrib = {0, 1, 4, 5};  // c++ type of indexing: index starts from zero!
    else
        distrib = {2, 3, 6, 7};

    std::string filename = INTERP_TEST::getResourceFile("SimpleTest2D.med");
    MCAuto<MEDFileField1TS> f1TS = ParaMEDFileField1TS::ParaNew(
        MPI_COMM_WORLD, MPI_INFO_NULL, filename, "fieldOnCells", "mesh", distrib, ON_CELLS, INTERP_KERNEL::NORM_QUAD4
    );
    MCAuto<MEDCouplingFieldDouble> fieldRef = genLocFieldCells(rank);
    CPPUNIT_ASSERT(f1TS->getUndergroundDataArray()->isEqual(*fieldRef->getArray(), 1e-12));
    MPI_Barrier(MPI_COMM_WORLD);
}

/*!
 * Test case to load a field located on nodes in parallel on 2 procs.
 */
void
ParaMEDMEMTest::testParallelLoad5()
{
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //
    if (size != 2)
        return;

    std::vector<mcIdType> distrib;
    if (rank == 0)
        distrib = {0, 1, 2, 5, 6, 7, 10, 11, 12};  // c++ type of indexing: index starts from zero!
    else
        distrib = {2, 3, 4, 7, 8, 9, 12, 13, 14};

    std::string filename = INTERP_TEST::getResourceFile("SimpleTest2D.med");
    // for fields on nodes, geometrical type is not needed
    MCAuto<MEDFileField1TS> f1TS = ParaMEDFileField1TS::ParaNew(
        MPI_COMM_WORLD, MPI_INFO_NULL, filename, "fieldOnNodes", "mesh", distrib, ON_NODES, INTERP_KERNEL::NORM_ERROR
    );
    MCAuto<MEDCouplingFieldDouble> fieldRef = genLocFieldNodes(rank);
    CPPUNIT_ASSERT(f1TS->getUndergroundDataArray()->isEqual(*fieldRef->getArray(), 1e-12));
    MPI_Barrier(MPI_COMM_WORLD);
}

/*!
 * Test case to load a 2D mesh with multiple geometric types in parallel on 3 procs.
 * Some procs may have empty partition to load for some types
 */
void
ParaMEDMEMTest::testParallelLoad6()
{
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //
    if (size != 3)
        return;

    std::map<INTERP_KERNEL::NormalizedCellType, std::vector<mcIdType>> distrib;
    // independant numerotation for each geometric type!
    if (rank == 0)
        distrib = {{INTERP_KERNEL::NORM_TRI3, {}}, {INTERP_KERNEL::NORM_QUAD4, {2}}};
    else if (rank == 1)
        distrib = {{INTERP_KERNEL::NORM_TRI3, {2}}, {INTERP_KERNEL::NORM_QUAD4, {}}};
    else
        distrib = {{INTERP_KERNEL::NORM_TRI3, {0}}, {INTERP_KERNEL::NORM_QUAD4, {1}}};

    std::string filename = INTERP_TEST::getResourceFile("Test2DMultiGeoType.med");
    MCAuto<MEDFileUMesh> mu = ParaMEDFileUMesh::ParaNew(distrib, MPI_COMM_WORLD, MPI_INFO_NULL, filename, "mesh");
    MCAuto<MEDCouplingUMesh> mesh = mu->getMeshAtLevel(0);

    MEDCouplingUMesh *meshRef;
    if (rank == 0)
        meshRef = genPartialLocMeshMultipleTypes1();
    else if (rank == 1)
        meshRef = genPartialLocMeshMultipleTypes2();
    else
        meshRef = genPartialLocMeshMultipleTypes3();
    // checking that all 3 procs have correctly loaded their part
    int equal = (int)mesh->isEqual(meshRef, 1e-12);
    int allEqual = -1;
    MPI_Allreduce(&equal, &allEqual, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    CPPUNIT_ASSERT(allEqual == 3);
    meshRef->decrRef();

    MPI_Barrier(MPI_COMM_WORLD);
}
