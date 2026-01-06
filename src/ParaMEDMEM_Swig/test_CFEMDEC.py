#!/usr/bin/env python3
#  -*- coding: utf-8 -*-
# Copyright (C) 2025-2026  CEA, EDF
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

"""
See EDF31187 : Test P1P1 on 7 procs. 3 source procs and 4 target procs.

This tests checks that DEC is OK for both sides (src->target and target->src)
Matrix computation is required for both sides
"""

import medcoupling as mc
from mpi4py import MPI

globalComm = MPI.COMM_WORLD

size = globalComm.size
rank = globalComm.rank

procs_source = [0, 1, 2]
procs_target = [3, 4, 5, 6]

idec = mc.CFEMDEC(procs_source, procs_target)


def createSourceFieldExNihilo(src_mesh):
    src_field = mc.MEDCouplingFieldDouble(mc.ON_NODES_FE)
    src_field.setMesh(src_mesh)
    #
    coords = src_mesh.getCoords()
    x = coords[:, 0]
    y = coords[:, 1]
    value = x**2 + y**2
    #
    src_field.setArray(value)
    src_field.setName("Field")
    src_field.setNature(mc.IntensiveMaximum)
    value.setInfoOnComponents(["ABC"])
    return src_field


def createFieldForParaView(zeField):
    retField = mc.MEDCouplingFieldDouble(mc.ON_NODES)
    retField.setName(zeField.getName())
    retField.setMesh(zeField.getMesh())
    retField.setArray(zeField.getArray())
    return retField


def computeInSequentialReferenceField(src_field, trg_Mesh):
    trgFt = mc.MEDCouplingFieldTemplate(mc.ON_NODES_FE)
    trgFt.setMesh(trg_Mesh)
    rem = mc.MEDCouplingRemapper()
    rem.setIntersectionType(mc.PointLocator)
    srcFt = mc.MEDCouplingFieldTemplate(src_field)
    rem.prepareEx(srcFt, trgFt)
    trg_field = rem.transferField(src_field, 1e300)
    return trg_field


def addGhostCells(fieldGlob, meshLoc):
    fetched_node_ids = meshLoc.computeFetchedNodeIds()
    ghostCells = fieldGlob.getMesh().getCellIdsLyingOnNodes(fetched_node_ids, False)
    initialNbOfNodes = fieldGlob.getMesh().getNumberOfNodes()
    ghostCells_mesh = fieldGlob.getMesh()[ghostCells]
    o2nIDS = ghostCells_mesh.zipCoordsTraducer()
    globalNodeIds = o2nIDS.invertArrayO2N2N2O(ghostCells_mesh.getNumberOfNodes())
    arrWithGhost = fieldGlob.getArray()[globalNodeIds]
    arrWithGhost.copyStringInfoFrom(fieldGlob.getArray())
    tabO2N = mc.DataArrayInt(initialNbOfNodes)
    tabO2N.iota()
    tabO2N[1::2] += 1000000000
    # Force a renumbering with high value to be sure that memory consumtion is OK
    globalNodeIds = tabO2N[globalNodeIds]
    return globalNodeIds, ghostCells_mesh, arrWithGhost


def buildLocalFieldFromGlobal(arrWithGhost, ghostCells_mesh):
    field_on_local = mc.MEDCouplingFieldDouble(mc.ON_NODES)
    field_on_local.setMesh(ghostCells_mesh)
    field_on_local.setName("globalnodeids_field")

    field_on_local.setArray(arrWithGhost)
    return field_on_local


if rank in procs_source:
    src_mesh = mc.MEDFileMesh.New("test_CFEMDEC.med", "src")[0]
    src_field = createSourceFieldExNihilo(src_mesh)
    #
    src_cellsIds_submesh_part0 = [
        24,
        25,
        4,
        30,
        9,
        35,
        14,
        40,
        19,
        45,
        0,
        26,
        5,
        31,
        10,
        36,
        15,
        41,
        20,
        46,
    ]
    src_cellsIds_submesh_part1 = [
        1,
        27,
        6,
        32,
        11,
        37,
        16,
        42,
        21,
        47,
        2,
        28,
        7,
        33,
        12,
        38,
        17,
        43,
        22,
        48,
    ]
    src_cellsIds_submesh_part2 = [3, 29, 8, 34, 13, 39, 18, 44, 23, 49]
    src_cell_groups = [
        src_cellsIds_submesh_part0,
        src_cellsIds_submesh_part1,
        src_cellsIds_submesh_part2,
    ]
    #
    src_submeshes = src_mesh[src_cell_groups[rank]]
    src_globalNodeIds, src_ghostCells_mesh, src_arrWithGhost = addGhostCells(
        src_field, src_submeshes
    )
    src_field_on_local = buildLocalFieldFromGlobal(
        src_arrWithGhost, src_ghostCells_mesh
    )

    # All data are prepared let s go for the test

    # first basic test. Scalar src -> trg
    idec.attachLocalMesh(src_ghostCells_mesh, src_globalNodeIds)
    idec.sendToTarget(src_field_on_local)

    # second test vector field src -> trg
    arr2 = mc.DataArrayDouble(src_field_on_local.getArray().getNumberOfTuples(), 2)
    arr2[:, 0] = src_field_on_local.getArray()
    arr2[:, 1] = arr2[:, 0] * 2
    arr2.setInfoOnComponents(["DEF", "GHIJK"])
    src_field_on_local.setArray(arr2)
    idec.sendToTarget(src_field_on_local)

    # third test scalar field trg -> src
    expected_values_on_whole_src_mesh = mc.DataArrayDouble(
        [
            20000.0,
            19999.999999999996,
            20000.0,
            20000.0,
            14148.571428571424,
            10765.714285714195,
            10765.714285714412,
            14148.571428571424,
            14148.571428571424,
            10765.714285714193,
            10765.71428571441,
            14148.571428571422,
            14148.571428571422,
            10765.714285714192,
            10765.714285714412,
            14148.571428571424,
            14148.571428571424,
            10765.714285714192,
            10765.714285714414,
            14148.571428571424,
            7868.681089602829,
            4867.991833215152,
            4344.494057912346,
            8028.916998049083,
            4878.38401782005,
            1337.8705692098215,
            1564.290580656142,
            4916.716392647952,
            4853.406143945703,
            1264.246476298977,
            1678.2209184671142,
            4902.165551383907,
            8052.863022697784,
            4671.9689984247125,
            4908.030552914453,
            8008.660190314067,
        ]
    )
    zeResu3 = idec.receiveFromTarget()
    assert zeResu3.getArray().getInfoOnComponents() == ["LMNOP"]
    a, n2o = src_mesh.getCoords().areIncludedInMe(zeResu3.getMesh().getCoords(), 1e-12)
    assert a
    assert expected_values_on_whole_src_mesh[n2o].isEqualWithoutConsideringStr(
        zeResu3.getArray(), 1e-11
    )

if rank in procs_target:
    # computed with computeInSequentialReferenceField
    expected_values_on_whole_target_mesh = mc.DataArrayDouble(
        [
            20000.0,
            20000.0,
            19999.999999999996,
            20000.0,
            15428.57142857149,
            12228.571428571518,
            10400.0,
            10400.0,
            12228.571428571246,
            15428.571428571244,
            15428.57142857149,
            12228.571428571518,
            10400.0,
            10400.0,
            12228.571428571246,
            15428.571428571246,
            15428.571428571488,
            12228.571428571517,
            10399.999999999998,
            10400.0,
            12228.57142857125,
            15428.571428571246,
            15428.571428571488,
            12228.571428571518,
            10399.999999999998,
            10400.0,
            12228.571428571246,
            15428.571428571246,
            9036.51525644838,
            6589.302297134688,
            5995.698676861471,
            6841.132031926167,
            9603.880571402078,
            9599.487442051013,
            6282.661530571854,
            9376.084921754715,
            9870.81512141087,
            5836.508790702477,
            9206.0374584893,
            9186.849454375477,
            6939.199661303852,
            6691.314512608502,
            7709.042737936002,
            6813.120840231214,
            7142.388603259622,
            6969.223286406583,
            6506.737808081094,
            8069.986820492897,
            2757.416083531578,
            2854.5200072667694,
            4931.1974413357775,
            5411.149737477748,
            2649.9939831637894,
            3388.4160359651282,
            4826.53771440446,
            4457.152560636136,
            2712.2232402537034,
            3335.925702172402,
            3679.3996445719995,
            2582.1283228800535,
            800.0000000000002,
            1275.8376494997938,
            1569.187581760969,
            800.0000000000001,
            1221.9603216011838,
        ]
    )
    #
    trg_cellsIds_submesh_part0 = [
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        19,
        27,
        44,
        46,
        47,
        48,
        49,
        50,
        57,
        62,
        68,
        69,
        72,
    ]
    trg_cellsIds_submesh_part1 = [
        0,
        11,
        31,
        12,
        29,
        13,
        14,
        40,
        30,
        32,
        28,
        33,
        34,
        45,
        51,
        52,
        53,
        54,
        55,
        56,
        58,
        67,
    ]
    trg_cellsIds_submesh_part2 = [17, 18, 26, 71, 61, 39, 42, 59, 65, 43, 73]
    trg_cellsIds_submesh_part3 = [
        15,
        16,
        20,
        21,
        22,
        23,
        24,
        25,
        35,
        36,
        37,
        38,
        41,
        60,
        63,
        64,
        66,
        70,
        74,
        75,
        76,
        77,
        78,
        79,
        80,
        81,
        82,
        83,
        84,
        85,
        86,
        87,
        88,
        89,
        90,
        91,
        92,
        93,
        94,
        95,
        96,
        97,
        98,
        99,
    ]
    trg_cell_groups = [
        trg_cellsIds_submesh_part0,
        trg_cellsIds_submesh_part1,
        trg_cellsIds_submesh_part2,
        trg_cellsIds_submesh_part3,
    ]
    trg_Mesh = mc.MEDFileMesh.New("test_CFEMDEC.med", "trg")[0]
    expected_field_target = mc.MEDCouplingFieldDouble(mc.ON_NODES)
    expected_field_target.setArray(expected_values_on_whole_target_mesh)
    expected_field_target.setMesh(trg_Mesh)
    #
    proc_trg_mesh = trg_Mesh[trg_cell_groups[procs_target.index(rank)]]
    #
    trg_globalNodeIds, trg_ghostCells_mesh, trg_arrWithGhost = addGhostCells(
        expected_field_target, proc_trg_mesh
    )
    expected_trg_field_on_local = buildLocalFieldFromGlobal(
        trg_arrWithGhost, trg_ghostCells_mesh
    )
    # All data are prepared let s go for the test

    # first basic test. Scalar src -> trg
    idec.attachLocalMesh(trg_ghostCells_mesh, trg_globalNodeIds)
    zeResu = idec.receiveFromSource()
    assert expected_trg_field_on_local.getArray().isEqualWithoutConsideringStr(
        zeResu.getArray(), 1e-11
    )
    assert zeResu.getArray().getInfoOnComponents() == ["ABC"]
    # createFieldForParaView(zeResu).writeVTK(f"trg_field_array{rank}.vtu")

    # second test vector field src -> trg
    zeResu2 = idec.receiveFromSource()
    assert zeResu2.getArray().getNumberOfComponents() == 2
    assert zeResu2.getArray().getInfoOnComponents() == ["DEF", "GHIJK"]
    assert expected_trg_field_on_local.getArray().isEqualWithoutConsideringStr(
        zeResu2.getArray()[:, 0], 1e-11
    )
    assert expected_trg_field_on_local.getArray().isEqualWithoutConsideringStr(
        zeResu2.getArray()[:, 1] / 2, 1e-11
    )

    # third test scalar field trg -> src
    zeResu.getArray().setInfoOnComponents(["LMNOP"])
    idec.sendToSource(zeResu)

# mpirun -np 7 xterm -e "gdb -x cmd.gdb --args python3 test_CFEMDEC.py"
