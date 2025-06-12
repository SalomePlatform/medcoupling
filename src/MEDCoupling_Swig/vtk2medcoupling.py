#! /usr/bin/env python3
#  -*- coding: utf-8 -*-
# Copyright (C) 2020-2025  CEA, EDF
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
# Author : Anthony GEAY (EDF R&D)

import vtk
from vtk.util import numpy_support
import medcoupling as mc
import numpy as np
import logging

logger = None


def getLogger():
    global logger
    if logger is None:
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger()
    return logger


def setLogLevel(lev):
    getLogger().setLevel(lev)


def mesh_convertor(fileName):
    # vtk.vtkDataSetReader()
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fileName)
    reader.Update()
    ug = reader.GetOutput()
    return mesh_convertor_mem(ug), ug


def mesh_convertor_mem(ug, dimRequested=-1):
    """
    Main entry
    """
    from distutils.version import LooseVersion

    vtkVersion = vtk.vtkVersion().GetVTKVersion()
    getLogger().debug(f"VTK version : {vtkVersion}")
    if LooseVersion(vtkVersion) <= LooseVersion("9.3.0"):
        return mesh_convertor_mem_93(ug, dimRequested)
    else:
        return mesh_convertor_mem_94(ug, dimRequested)


def get_coarse_data_from_VTK_unstructuredGrid93(ug):
    facesLoc = numpy_support.vtk_to_numpy(ug.GetFaceLocations())
    faces = numpy_support.vtk_to_numpy(ug.GetFaces())
    return facesLoc, faces


def generate_polyhedrons_mesh_from_VTK_UnstructuredGrid_93(ug):
    """
    Returns from vtkUnstructuredGrid instance MEDCouplingMesh containing only polyhedra cells
    """
    facesLoc, faces = get_coarse_data_from_VTK_unstructuredGrid93(ug)
    facesLocMed = mc.DataArrayInt(facesLoc)
    facesMed = mc.DataArrayInt(faces)
    cellIdsPolyh = facesLocMed.findIdsNotEqual(-1)
    facesLocCpy = facesLocMed[:]
    #
    polyhStart = facesLocCpy[cellIdsPolyh]
    polyhStart.pushBackSilent(len(facesMed))
    a, b = mc.DataArrayInt64.FromVTKInternalReprOfPolyedra(facesMed, polyhStart)
    m = mc.MEDCoupling1DGTUMesh("mesh", mc.NORM_POLYHED)
    m.setNodalConnectivity(a, b)
    return m.buildUnstructured()


def generate_polyhedrons_mesh_from_VTK_UnstructuredGrid_94(ug):
    facesLoc = (
        ug.GetPolyhedronFaceLocations()
    )  # vtkCellArray.SetData(offsets,connectivity)
    facesLoc_o = mc.DataArrayInt(numpy_support.vtk_to_numpy(facesLoc.GetOffsetsArray()))
    facesLoc_c = mc.DataArrayInt(
        numpy_support.vtk_to_numpy(facesLoc.GetConnectivityArray())
    )

    faces = ug.GetPolyhedronFaces()
    faces_o = mc.DataArrayInt(numpy_support.vtk_to_numpy(faces.GetOffsetsArray()))
    faces_c = mc.DataArrayInt(numpy_support.vtk_to_numpy(faces.GetConnectivityArray()))

    nbFacesPerCell = facesLoc_o.deltaShiftIndex()
    polyhedronsCellIds = nbFacesPerCell.findIdsNotEqual(0)

    _, facesLoc_o_eff = mc.DataArrayInt.ExtractFromIndexedArrays(
        polyhedronsCellIds, facesLoc_c, facesLoc_o
    )

    a, b = mc.DataArrayInt64.FromVTK94InternalReprOfPolyedra(
        faces_c, faces_o, facesLoc_o_eff
    )

    m = mc.MEDCoupling1DGTUMesh("mesh", mc.NORM_POLYHED)
    m.setNodalConnectivity(a, b)
    return m.buildUnstructured()


def from_93_to_94_coarse_VTK_data(ug):
    facesLoc, faces, _ = get_coarse_data_from_VTK_unstructuredGrid93(ug)
    facesLocMed = mc.DataArrayInt(facesLoc)
    facesMed = mc.DataArrayInt(faces)
    cellIdsPolyh = facesLocMed.findIdsNotEqual(-1)
    nbFacesPerPolyh = facesMed[facesLocMed[cellIdsPolyh]]
    facesLoc_o = mc.DataArrayInt(len(facesLocMed))
    facesLoc_o[:] = 0
    facesLoc_o[cellIdsPolyh] = nbFacesPerPolyh
    facesLoc_o.computeOffsetsFull()
    faceLocMedNoMinusOne = facesLocMed[cellIdsPolyh]
    faceLocMedNoMinusOne.pushBackSilent(len(facesMed))
    faces_c, faces_o = mc.DataArrayInt64.FromVTK93To94FacesInternaReprOfPolyedra(
        facesMed, faceLocMedNoMinusOne
    )
    nbFacesInPolyhedra = facesLoc_o[-1]
    facesLoc_c = mc.DataArrayInt(nbFacesInPolyhedra)
    facesLoc_c.iota()
    return facesLoc_o, facesLoc_c, faces_o, faces_c


def get_coarse_data_from_VTK_unstructuredGrid_94(ug):
    facesLoc = ug.GetPolyhedronFaceLocations()
    facesLoc_o = numpy_support.vtk_to_numpy(facesLoc.GetOffsetsArray())
    facesLoc_c = numpy_support.vtk_to_numpy(facesLoc.GetConnectivityArray())

    faces = ug.GetPolyhedronFaces()
    faces_o = numpy_support.vtk_to_numpy(faces.GetOffsetsArray())
    faces_c = numpy_support.vtk_to_numpy(faces.GetConnectivityArray())
    return facesLoc_o, facesLoc_c, faces_o, faces_c


def mesh_convertor_mem_gen(ug, dimRequested, callback_for_polyedrons):
    coords = mc.DataArrayDouble(
        np.array(numpy_support.vtk_to_numpy(ug.GetPoints().GetData()), dtype=np.float64)
    )
    cells = ug.GetCells()
    offsets = numpy_support.vtk_to_numpy(cells.GetOffsetsArray())
    ci = mc.DataArrayInt(np.array(offsets, dtype=np.int64))[:]
    conn = mc.DataArrayInt(
        np.array(
            numpy_support.vtk_to_numpy(cells.GetConnectivityArray()), dtype=np.int64
        )
    )
    types = numpy_support.vtk_to_numpy(ug.GetCellTypesArray())
    ct = mc.DataArrayInt(
        np.array(types, dtype="int{}".format(mc.MEDCouplingSizeOfIDs()))
    )[:]
    vtk2med = mc.DataArrayInt(mc.vtk2med_cell_types())
    ctdim = ct[:]
    ct.transformWithIndArr(vtk2med)

    activeType = vtk2med.findIdsNotEqual(-1)
    vtk2dim = vtk2med[:]
    vtk2dim[activeType] = mc.DataArrayInt(
        [
            mc.MEDCouplingUMesh.GetDimensionOfGeometricType(int(elt))
            for elt in vtk2med[activeType]
        ]
    )
    ctdim.transformWithIndArr(vtk2dim)
    dims = ctdim.getDifferentValues().getValues()
    getLogger().debug(f"Available dimensions : {dims}")
    zeDim = dimRequested
    if zeDim < 0:
        zeDim = max(dims)
    else:
        if zeDim not in dims:
            raise RuntimeError(
                f"{zeDim} not in detected dimension of vtkUnstructuredGrid instance : {dims}"
            )
    getLogger().debug(f"Filtering on dimension : {zeDim}")
    selectedCells = ctdim.findIdsEqual(zeDim)
    conn1, ci1 = mc.DataArrayInt.ExtractFromIndexedArrays(selectedCells, conn, ci)
    ct1 = ct[selectedCells]

    classicalCells = ct1.findIdsNotEqual(mc.NORM_POLYHED)
    isPostWithPoly = ct1.presenceOfValue(mc.NORM_POLYHED)
    resu = []
    resu2 = []
    if not classicalCells.empty():
        getLogger().debug(f"Presence of {len(classicalCells)} classical cells")
        conn2, ci2 = mc.DataArrayInt.ExtractFromIndexedArrays(
            classicalCells, conn1, ci1
        )
        ct2 = ct1[classicalCells]
        m = mc.MEDCouplingUMesh("mesh", zeDim)
        # insert geo type in connectivity
        a = ci2.deltaShiftIndex()
        a2 = a + 1
        a2.computeOffsetsFull()
        b = mc.DataArrayInt(len(a))
        b[:] = 1
        c = mc.DataArrayInt.Meld([b, a])
        c.rearrange(1)
        c.computeOffsetsFull()
        d = mc.DataArrayInt(len(classicalCells))
        d.iota()
        d *= 2
        d += 1
        d2 = d.buildExplicitArrByRanges(c)
        d3 = mc.DataArrayInt(len(conn2) + len(classicalCells))
        d3[:] = -1
        idsForGeoType = d2.buildComplement(len(conn2) + len(classicalCells))
        d3[idsForGeoType] = ct2
        d3[d2] = conn2
        #
        m.setConnectivity(d3, a2, True)
        m.setCoords(coords)
        resu2.append(classicalCells)
        resu.append(m)

    if isPostWithPoly:
        getLogger().debug(f"Presence of polyhedrons cells")
        resu2.append(ct1.findIdsEqual(mc.NORM_POLYHED))
        m = callback_for_polyedrons(ug)
        m.setCoords(coords)
        resu.append(m)

    resu2 = mc.DataArrayInt.Aggregate(resu2)
    resu = mc.MEDCouplingUMesh.MergeUMeshesOnSameCoords(resu)
    ren = resu2.invertArrayO2N2N2O(len(resu2))
    ret = resu[ren]
    return ret


def mesh_convertor_mem_93(ug, dimRequested=-1):
    getLogger().debug(f"Reader 9.3")
    return mesh_convertor_mem_gen(
        ug, dimRequested, generate_polyhedrons_mesh_from_VTK_UnstructuredGrid_93
    )


def mesh_convertor_mem_94(ug, dimRequested=-1):
    getLogger().debug(f"Reader 9.4")
    return mesh_convertor_mem_gen(
        ug, dimRequested, generate_polyhedrons_mesh_from_VTK_UnstructuredGrid_94
    )


def mesh_convertor_mem_93_legacy(ug):
    def patchForPolyedra(polyhedCellIds, ug, mesh):
        """
        Method in charge to change the connectivity of polyedra contained in mesh using ug vtkUnstructuredGrid.

        :param in polyhedCellIds: mc.DataArrayInt of cells ids in mesh to be patched
        :param in ug: vtkUnstructuredGrid containing polyhedra
        :param in-out mesh: mc.MEDCouplingUMesh. 3D Mesh whose polyedra cells connectivity will be modified
        """
        facesLoc = mc.DataArrayInt(numpy_support.vtk_to_numpy(ug.GetFaceLocations()))
        faces = mc.DataArrayInt(numpy_support.vtk_to_numpy(ug.GetFaces()))
        facesLoc = facesLoc[polyhedCellIds]
        facesLoc = mc.DataArrayInt.Aggregate([facesLoc, mc.DataArrayInt([len(faces)])])
        connForPoly, facesLoc = mc.DataArrayInt.FromVTKInternalReprOfPolyedra(
            faces, facesLoc
        )
        meshPoly = mc.MEDCoupling1DGTUMesh(mesh.getName(), mc.NORM_POLYHED)
        meshPoly.setCoords(mesh.getCoords())
        meshPoly.setNodalConnectivity(connForPoly, facesLoc)
        mesh[polyhedCellIds] = meshPoly.buildUnstructured()
        pass

    from distutils.version import LooseVersion

    #
    pts = numpy_support.vtk_to_numpy(ug.GetPoints().GetData())
    #
    cla = numpy_support.vtk_to_numpy(ug.GetCellLocationsArray())
    ctvtk = numpy_support.vtk_to_numpy(ug.GetCellTypesArray())
    conn = numpy_support.vtk_to_numpy(ug.GetCells().GetData())
    #
    ct = mc.DataArrayInt(
        np.array(ctvtk, dtype="int{}".format(mc.MEDCouplingSizeOfIDs()))
    )[:]
    c = mc.DataArrayInt(conn)[:]
    ci = mc.DataArrayInt(cla)[:]
    # for pv580
    if LooseVersion(vtk.vtkVersion().GetVTKVersion()) >= LooseVersion("8.90.0"):
        ci = ci.deltaShiftIndex() + 1
        ci.computeOffsetsFull()
    #
    vtk2med = mc.DataArrayInt(mc.vtk2med_cell_types())
    #
    ct.transformWithIndArr(vtk2med)
    c[ci] = ct
    ci = mc.DataArrayInt.Aggregate([ci, mc.DataArrayInt([len(c)])])
    #
    gtv = ct.getDifferentValues().getValues()
    dimArray = mc.DataArrayInt(len(gtv))
    dimArray[:] = -1
    for i, gt in enumerate(gtv):
        dimArray[i] = mc.MEDCouplingUMesh.GetDimensionOfGeometricType(gt)
    dim = dimArray.getMaxValueInArray()
    if not dimArray.isUniform(dim):
        raise RuntimeError(
            "The input vtkUnstructuredGrid instance is not a single level mesh ! need to split it !"
        )
    #
    m = mc.MEDCouplingUMesh("mesh", dim)
    m.setCoords(mc.DataArrayDouble(np.array(pts, dtype=np.float64)))
    m.setConnectivity(c, ci, True)
    m.checkConsistencyLight()
    #
    if m.getMeshDimension() == 3:
        polyhedCellIds = ct.findIdsEqual(mc.NORM_POLYHED)
        if not polyhedCellIds.empty():
            patchForPolyedra(polyhedCellIds, ug, m)
    #
    return m
