# -*- coding: utf-8 -*-
# Copyright (C) 2007-2025  CEA, EDF
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
import unittest
import medcoupling as mc

EPS = 1.0e-8
COORDS = [
    4.0, -4.4408920985e-16,
    2.82842712475,  2.82842712475,
    2.44929359829e-16,  4.0,
    -2.82842712475, 2.82842712475,
    -4.0, 1.12948887387e-15,
    -2.82842712475, -2.82842712475,
    -6.43249059871e-16, -4.0,
    2.82842712475,  -2.82842712475,
    0.0, 0.0,
    0.798636, 0.601815,
    1.597271, 1.203630,
    2.395907, 1.805445,
    3.194542, 2.407260,
]

# quadratic 2D: [v0..v(num_verts-1), m01, m12, ...]
QPOLYG_CURVED   = [0, 2, 4, 6, 1, 3, 5, 7]
QUAD8_CURVED    = [0, 2, 4, 6, 1, 3, 5, 7]
TRI6_CURVED     = [0, 2, 4, 1, 3, 8]
QPOLYG_STRAIGHT = [8, 10, 12, 9, 11, 10]

def coords_has_point(da_coords, pt, tol=1e-9):
    n = da_coords.getNumberOfTuples()
    x, y = pt
    for i in range(n):
        if abs(da_coords.getIJ(i,0)-x) <= tol and abs(da_coords.getIJ(i,1)-y) <= tol:
            return True
    return False

def build_single_cell_mesh(elem_type, conn, xy):
    mesh = mc.MEDCouplingUMesh.New("mesh", 2)
    mesh.allocateCells(1)
    da = mc.DataArrayDouble.New(xy, len(xy)//2, 2)
    mesh.setCoords(da)
    mesh.insertNextCell(elem_type, conn)
    return mesh

class ColinearizeKeepConformTest(unittest.TestCase):
    def run_pipeline(self, elem_type, conn, xy):
        mesh = build_single_cell_mesh(elem_type, conn, xy)
        removed = mesh.colinearizeKeepingConform2D(EPS)
        if removed.getNumberOfTuples():
            mesh.mergeNodes(EPS)
            mesh.conformize2D(EPS)
            mesh.zipCoords()
        return mesh

    def assert_original_vertices_preserved(self, elem_type, conn, xy):
        conn_len = len(conn)
        num_verts = (conn_len - 1)//2 if (conn_len % 2)==1 else conn_len//2  
        verts = conn[:num_verts]

        mesh_out = self.run_pipeline(elem_type, conn, xy)
        out_coords = mesh_out.getCoords()

        for v in verts:
            px, py = xy[2*v], xy[2*v+1]
            self.assertTrue(coords_has_point(out_coords, (px, py), 1e-8),
                            msg=f"Original vertex {v} was removed")

    def test_qpolyg_curved_vertices_preserved(self):
        self.assert_original_vertices_preserved(mc.NORM_QPOLYG, QPOLYG_CURVED, COORDS)

    def test_quad8_curved_vertices_preserved(self):
        self.assert_original_vertices_preserved(mc.NORM_QUAD8, QUAD8_CURVED, COORDS)

    def test_tri6_curved_vertices_preserved(self):
        self.assert_original_vertices_preserved(mc.NORM_TRI6,  TRI6_CURVED, COORDS)

    def test_qpolyg_straight_pipeline_ok(self):
        mesh = self.run_pipeline(mc.NORM_QPOLYG, QPOLYG_STRAIGHT, COORDS)
        self.assertEqual(2, mesh.getSpaceDimension())
        self.assertGreaterEqual(mesh.getNumberOfCells(), 1)
        self.assertGreaterEqual(mesh.getNumberOfNodes(), 3)

if __name__ == "__main__":
    unittest.main()
