// Copyright (C) 2026  CEA, EDF
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

%pythoncode %{
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    import numpy as np
    import pyvista as pv

MC_DIM = {
    NORM_POINT1: 0,
    NORM_QUAD4: 2,
    NORM_TRI3: 2,
    NORM_HEXA8: 3,
    NORM_SEG2: 1,
    NORM_SEG3: 1,
    NORM_TETRA4: 3,
    NORM_POLYGON: 2,
    NORM_QPOLYG: 3,  # only UnstructuredGrid
    NORM_POLYHED: 3,
}



def _mesh_to_pv(mesh):
    import numpy as np
    def _mc_ph_to_vtk_fast(connectivity: np.ndarray, offsets: np.ndarray) -> np.ndarray:
        if len(connectivity) == 0:
            return np.array([], dtype=int)
        cell_length = offsets[1:] - offsets[:-1]
        face_delims = connectivity == -1  # all face seperator
        face_delims[offsets[:-1]] = True  # all cell seperator

        face_offsets = np.r_[np.flatnonzero(face_delims), len(connectivity)]

        # Compute number of face per cell
        cell_delim_in_face_offsets = face_offsets.searchsorted(offsets)
        num_faces = cell_delim_in_face_offsets[1:] - cell_delim_in_face_offsets[:-1]

        # Replace -1 and cell_type with face length
        face_length = face_offsets[1:] - face_offsets[:-1] - 1
        connectivity[face_delims] = face_length

        inds = np.empty((2 * len(cell_length),), dtype=int)
        inds[::2] = offsets[:-1]
        inds[1::2] = offsets[:-1]

        vals = np.empty((2 * len(cell_length),), dtype=int)
        vals[::2] = cell_length + 1
        vals[1::2] = num_faces

        return np.insert(connectivity, obj=inds, values=vals)


    def _to_unstructured(
        mesh: MEDCouplingUMesh, coords: np.ndarray
    ) -> "pv.UnstructuredGrid":
        import pyvista as pv
        MC_TO_PV_CELLTYPE = {
            NORM_POINT1: pv.CellType.VERTEX,
            NORM_QUAD4: pv.CellType.QUAD,
            NORM_TRI3: pv.CellType.TRIANGLE,
            NORM_HEXA8: pv.CellType.HEXAHEDRON,
            NORM_SEG2: pv.CellType.LINE,
            NORM_SEG3: pv.CellType.QUADRATIC_EDGE,
            NORM_TETRA4: pv.CellType.TETRA,
            NORM_POLYGON: pv.CellType.POLYGON,
            NORM_QPOLYG: pv.CellType.QUADRATIC_POLYGON,
            NORM_POLYHED: pv.CellType.POLYHEDRON,
        }

        offsets = mesh.getNodalConnectivityIndex().toNumPyArray()
        cell_length = offsets[1:] - offsets[:-1] - 1

        cell_types_idx = offsets[:-1]
        connectivity = np.array(mesh.getNodalConnectivity().toNumPyArray())
        mc_cell_types = np.array(connectivity[cell_types_idx])

        # Split off polyhedrons
        ind_l = np.searchsorted(mc_cell_types, NORM_POLYHED, side="left")
        ind_r = np.searchsorted(mc_cell_types, NORM_POLYHED, side="right")
        cell_types = np.array(mc_cell_types)

        # Non polyhedron case
        types_idx = dict()
        for mc_type in MC_TO_PV_CELLTYPE:
            types_idx[mc_type] = cell_types == mc_type
        for mc_type, pv_type in MC_TO_PV_CELLTYPE.items():
            cell_types[types_idx[mc_type]] = pv_type

        connectivity_npl = connectivity[: offsets[ind_l]]
        connectivity_npl[cell_types_idx[:ind_l]] = cell_length[:ind_l]
        connectivity_npr = connectivity[offsets[ind_r] :]
        connectivity_npr[cell_types_idx[ind_r:]] = cell_length[ind_r:]

        # Polyhedron case
        connectivity_ph = _mc_ph_to_vtk_fast(
            connectivity[offsets[ind_l] : offsets[ind_r]],
            offsets=offsets[ind_l : ind_r + 1] - offsets[ind_l],
        )

        # Combination
        connectivity = np.r_[connectivity_npl, connectivity_npr, connectivity_ph]

        pv_mesh = pv.UnstructuredGrid(connectivity, cell_types, coords)
        return pv_mesh


    def _to_polydata(mesh: MEDCouplingUMesh, coords: np.ndarray) -> "pv.PolyData":
        import pyvista as pv
        offsets = mesh.getNodalConnectivityIndex().toNumPyArray()
        cell_length = offsets[1:] - offsets[:-1] - 1
        cell_types_idx = offsets[:-1]

        connectivity = np.array(mesh.getNodalConnectivity().toNumPyArray())
        any_type = connectivity[0]
        connectivity[cell_types_idx] = cell_length

        if MC_DIM[any_type] == 0:
            return pv.PolyData(coords, verts=connectivity)
        if MC_DIM[any_type] == 1:
            return pv.PolyData(coords, lines=connectivity)
        if MC_DIM[any_type] == 2:
            return pv.PolyData(coords, faces=connectivity)
        else:
            raise Exception(f"{any_type} if not in known types: {MC_DIM=}")

    # Prepare mesh : make copy and sort cell types
    mesh = mesh.deepCopyConnectivityOnly()
    permut = mesh.sortCellsInMEDFileFrmt()

    coords = np.array(mesh.getCoords().toNumPyArray())
    coords = np.c_[coords, np.zeros((coords.shape[0], 3 - coords.shape[1]))]

    offsets = mesh.getNodalConnectivityIndex().toNumPyArray()
    connectivity = np.array(mesh.getNodalConnectivity().toNumPyArray())
    last_cell_type = connectivity[offsets[-2]]
    max_type = MC_DIM[last_cell_type]

    if max_type == 3:
        pv_mesh = _to_unstructured(mesh, coords)
    elif max_type < 3:
        pv_mesh = _to_polydata(mesh, coords)
    else:
        raise Exception(f"The coords shape is not valid: {coords.shape=}")
    return permut, pv_mesh

def _field_to_pv(field) -> "pv.PolyData | pv.UnstructuredGrid":
    import numpy as np
    mesh = field.getMesh()
    permut, pv_mesh = _mesh_to_pv(mesh)
    fname = mesh.getName()
    farr = field.getArray().toNumPyArray()
    if fname is None or fname == "":
        fname = "Field"
    if field.getTypeOfField() == ON_CELLS:
        permut = permut.toNumPyArray()
        p = np.empty_like(permut)
        p[permut] = np.arange(permut.size)
        pv_mesh.cell_data[fname] = farr[p]
    elif field.getTypeOfField() == ON_NODES:
        pv_mesh.point_data[fname] = farr
    else:
        raise NotImplementedError(f"{field.getTypeOfField()=} is not supported.")
    return pv_mesh

def _mesh_to_pv2(mesh) -> "pv.PolyData | pv.UnstructuredGrid":
    _, pv_mesh = _mesh_to_pv(mesh)
    return pv_mesh

MEDCouplingUMesh.toPyvista = _mesh_to_pv2
del _mesh_to_pv2
MEDCouplingFieldDouble.toPyvista = _field_to_pv
del _field_to_pv
%}
