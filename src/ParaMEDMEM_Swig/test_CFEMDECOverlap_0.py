#!/usr/bin/env python3
#  -*- coding: utf-8 -*-
# Copyright (C) 2026  CEA, EDF
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

# See EDF35712


def MyAssert(clue):
    if not clue:
        raise RuntimeError("My assertion failed !")


import medcoupling as mc

from mpi4py import MPI


def MyAssert(v):
    if not v:
        raise RuntimeError("Assertion failed !")


globalComm = MPI.COMM_WORLD

size = globalComm.size
rank = globalComm.rank

if size != 3:
    raise RuntimeError("Expected to be lanched with 3 procs !")

fieldName = "Toto"


def generateSrcMeshExNihilo(nbNodesFine: int) -> mc.MEDCouplingUMesh:
    srcArr = mc.DataArrayDouble(nbNodesFine)
    srcArr.iota()
    srcArr /= nbNodesFine - 1
    mSrc = mc.MEDCouplingCMesh()
    mSrc.setCoords(srcArr, srcArr)
    mSrc = mSrc.buildUnstructured()
    #
    srcIds2_0 = mSrc.computeIsoBarycenterOfNodesPerCell()[:, 0].findIdsLowerThan(1 / 3)
    srcIds2_1 = mSrc.computeIsoBarycenterOfNodesPerCell()[:, 1].findIdsLowerThan(1 / 3)
    srcIds2 = mc.DataArrayInt.BuildUnion([srcIds2_0, srcIds2_1])
    srcIds0_1 = srcIds2.buildComplement(mSrc.getNumberOfCells())
    srcIds0 = srcIds0_1[
        mSrc.computeIsoBarycenterOfNodesPerCell()[srcIds0_1, 0].findIdsLowerThan(2 / 3)
    ]
    srcIds1 = srcIds0_1[
        mSrc.computeIsoBarycenterOfNodesPerCell()[srcIds0_1, 0].findIdsGreaterOrEqualTo(
            2 / 3
        )
    ]
    #
    checkSrcIds = mc.DataArrayInt.Aggregate([srcIds0, srcIds1, srcIds2])
    checkSrcIds.sort()
    MyAssert(checkSrcIds.isIota(mSrc.getNumberOfCells()))
    return mSrc, [srcIds0, srcIds1, srcIds2]


def generateTrgMeshExNihilo(nbNodesCoarse: int) -> mc.MEDCouplingUMesh:
    trgArr = mc.DataArrayDouble(nbNodesCoarse)
    trgArr.iota()
    trgArr /= nbNodesCoarse - 1
    mTrg = mc.MEDCouplingCMesh()
    mTrg.setCoords(trgArr, trgArr)
    mTrg = mTrg.buildUnstructured()
    #
    trgIds0_0 = mTrg.computeIsoBarycenterOfNodesPerCell()[:, 0].findIdsGreaterThan(0.8)
    trgIds0 = trgIds0_0[
        mTrg[trgIds0_0].computeIsoBarycenterOfNodesPerCell()[:, 1].findIdsLowerThan(0.2)
    ]
    trgIds_1_2 = trgIds0.buildComplement(mTrg.getNumberOfCells())
    bary_1_2 = mTrg[trgIds_1_2].computeIsoBarycenterOfNodesPerCell()
    thres = bary_1_2[:, 1] - (1.0 - bary_1_2[:, 0])
    trgIds1 = trgIds_1_2[thres.findIdsGreaterOrEqualTo(0.0)]
    trgIds2 = trgIds_1_2[thres.findIdsLowerThan(0.0)]
    #
    checkTrgIds = mc.DataArrayInt.Aggregate([trgIds0, trgIds1, trgIds2])
    checkTrgIds.sort()
    MyAssert(checkTrgIds.isIota(mTrg.getNumberOfCells()))
    return mTrg, [trgIds0, trgIds1, trgIds2]


def createSourceFieldExNihilo(
    srcMesh: mc.MEDCouplingUMesh,
) -> mc.MEDCouplingFieldDouble:
    srcField = mc.MEDCouplingFieldDouble(mc.ON_NODES_FE)
    srcField.setMesh(srcMesh)
    srcField.setName(fieldName)
    #
    coords = srcMesh.getCoords()
    x = coords[:, 0]
    y = coords[:, 1]
    value = x**2 + y**2
    #
    srcField.setArray(value)
    srcField.setNature(mc.IntensiveMaximum)
    value.setInfoOnComponents(["ABC"])
    return srcField


def createFieldForParaView(
    zeField: mc.MEDCouplingFieldDouble,
) -> mc.MEDCouplingFieldDouble:
    retField = mc.MEDCouplingFieldDouble(mc.ON_NODES)
    retField.setName(zeField.getName())
    retField.setMesh(zeField.getMesh())
    retField.setArray(zeField.getArray())
    return retField


def addGhostCellsMesh(meshGlob: mc.MEDCouplingUMesh, meshLoc: mc.MEDCouplingUMesh):
    """
    meshLoc and meshGlob are supposed to lie on same set of nodes.
    return:
    - mesh without orphan nodes meshLoc + layer of ghost cells (taken from meshGlob)
    - global node ids for OverlapCFEMDEC. Voluntary this array is not directly numbering in meshGlob
    - global node ids in meshGlob node numbering
    """
    fetched_node_ids = meshLoc.computeFetchedNodeIds()
    ghostCells = meshGlob.getCellIdsLyingOnNodes(fetched_node_ids, False)
    ghostCells_mesh = meshGlob[ghostCells]
    o2nIDS = ghostCells_mesh.zipCoordsTraducer()
    globalNodeIdsNonModified = o2nIDS.invertArrayO2N2N2O(
        ghostCells_mesh.getNumberOfNodes()
    )
    #
    initialNbOfNodes = meshGlob.getNumberOfNodes()
    tabO2N = mc.DataArrayInt(initialNbOfNodes)
    tabO2N.iota()
    tabO2N[1::2] += 1000000000
    # Force a renumbering with high value to be sure that memory consumtion is OK
    globalNodeIds = tabO2N[globalNodeIdsNonModified]
    return ghostCells_mesh, globalNodeIdsNonModified, globalNodeIds


def addGhostCells(
    srcFieldGlob: mc.MEDCouplingFieldDouble, srcMesh: mc.MEDCouplingUMesh
):
    """
    return :
    - field reduced to srcMesh + layer of ghost added without orphan nodes
    - global node ids for OverlapCFEMDEC. Voluntary this array is not directly numbering in meshGlob
    """
    ghostCells_mesh, globalNodeIdsNonModified, globalNodeIds = addGhostCellsMesh(
        srcFieldGlob.getMesh(), srcMesh
    )
    arrWithGhost = srcFieldGlob.getArray()[globalNodeIdsNonModified]
    arrWithGhost.copyStringInfoFrom(srcFieldGlob.getArray())
    #
    srcFieldLoc = mc.MEDCouplingFieldDouble(mc.ON_NODES)
    srcFieldLoc.setName(srcFieldGlob.getName())
    srcFieldLoc.setMesh(ghostCells_mesh)
    srcFieldLoc.setArray(arrWithGhost)
    return srcFieldLoc, globalNodeIds


def computeInSequentialReferenceField(src_field, trg_Mesh) -> mc.MEDCouplingFieldDouble:
    trgFt = mc.MEDCouplingFieldTemplate(mc.ON_NODES_FE)
    trgFt.setMesh(trg_Mesh)
    rem = mc.MEDCouplingRemapper()
    rem.setIntersectionType(mc.PointLocator)
    srcFt = mc.MEDCouplingFieldTemplate(src_field)
    rem.prepareEx(srcFt, trgFt)
    trg_field = rem.transferField(src_field, 1e300)
    return trg_field


def createFieldForParaView(zeField) -> mc.MEDCouplingFieldDouble:
    retField = mc.MEDCouplingFieldDouble(mc.ON_NODES)
    retField.setName(zeField.getName())
    retField.setMesh(zeField.getMesh())
    retField.setArray(zeField.getArray())
    return retField


mSrcGlobal, srcIds = generateSrcMeshExNihilo(nbNodesFine=101)
mTrgGlobal, trgIds = generateTrgMeshExNihilo(nbNodesCoarse=16)

srcMesh = mSrcGlobal[srcIds[rank]]
srcField = createSourceFieldExNihilo(srcMesh)
srcFieldLoc, srcGlobNodeIds = addGhostCells(srcField, srcMesh)

trgMesh = mTrgGlobal[trgIds[rank]]
trgMeshWithGhost, trgGlobalNodeIds, trgGlobalNodeIds4Dec = addGhostCellsMesh(
    mTrgGlobal, trgMesh
)

# compute reference field
trgFieldRef = computeInSequentialReferenceField(
    createSourceFieldExNihilo(mSrcGlobal), mTrgGlobal
)
# createFieldForParaView( trgFieldRef ).writeVTK("res.vtu")

dec = mc.OverlapCFEMDEC(list(range(size)))
dec.attachSourceField(srcFieldLoc, srcGlobNodeIds)
dec.attachTargetMesh(trgMeshWithGhost, trgGlobalNodeIds4Dec)
dec.synchronize()
resTrgField = dec.computeTargetField()

MyAssert(resTrgField.getName() == fieldName)
MyAssert(
    resTrgField.getMesh().getHiddenCppPointer()
    == trgMeshWithGhost.getHiddenCppPointer()
)
MyAssert(
    resTrgField.getArray().isEqual(trgFieldRef.getArray()[trgGlobalNodeIds], 1e-12)
)  # <- big test is here

srcFieldLoc.getArray()[:] *= 2
resTrgField = dec.computeTargetField()
MyAssert(resTrgField.getName() == fieldName)
MyAssert(
    resTrgField.getMesh().getHiddenCppPointer()
    == trgMeshWithGhost.getHiddenCppPointer()
)
MyAssert(
    resTrgField.getArray().isEqual(2 * trgFieldRef.getArray()[trgGlobalNodeIds], 1e-12)
)
