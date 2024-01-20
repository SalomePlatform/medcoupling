#!/usr/bin/env python
#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2020-2024  CEA, EDF
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

__doc__ = """Here a mesh is initially loaded into two parts in the 2 procs along Y axis ( low part and upper part ).
Then a redistribution of cells is performed along X axis ( left part and right part )
"""

import medcoupling as mc

from mpi4py import MPI

def GenerateGlobalMesh():
    mesh = mc.MEDCouplingCMesh()
    arr = mc.DataArrayDouble(11) ; arr.iota()
    mesh.setCoords(arr,arr)
    mesh = mesh.buildUnstructured()
    return mesh

def MeshPart0():
    mesh = GenerateGlobalMesh()
    p0 = mesh[list(range(50))]
    nodeIds = p0.computeFetchedNodeIds()
    p0.zipCoords()
    return mc.ParaUMesh(p0,mc.DataArrayInt(list(range(50))),nodeIds)
    
def MeshPart1():
    mesh = GenerateGlobalMesh()
    cellIds = mc.DataArrayInt(list(range(50,100,1)))
    p1 = mesh[cellIds]
    nodeIds = p1.computeFetchedNodeIds()
    p1.zipCoords()
    return mc.ParaUMesh(p1,cellIds,nodeIds)

workPerProc = {0 : MeshPart0 , 1 : MeshPart1}
distribPerProc = { 0 : [0, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37, 38, 44, 45, 46, 47, 48, 49, 55, 56, 57, 58, 59, 60, 66, 67, 68, 69, 70, 71, 77, 78, 79, 80, 81, 82, 88, 89, 90, 91, 92, 93, 99, 100, 101, 102, 103, 104, 110, 111, 112, 113, 114, 115],
1 : [6, 7, 8, 9, 10, 17, 18, 19, 20, 21, 28, 29, 30, 31, 32, 39, 40, 41, 42, 43, 50, 51, 52, 53, 54, 61, 62, 63, 64, 65, 72, 73, 74, 75, 76, 83, 84, 85, 86, 87, 94, 95, 96, 97, 98, 105, 106, 107, 108, 109, 116, 117, 118, 119, 120]}

cellsExpected = {0 : [0, 1, 2, 3, 4, 10, 11, 12, 13, 14, 20, 21, 22, 23, 24, 30, 31, 32, 33, 34, 40, 41, 42, 43, 44, 50, 51, 52, 53, 54, 60, 61, 62, 63, 64, 70, 71, 72, 73, 74, 80, 81, 82, 83, 84, 90, 91, 92, 93, 94],
1 : [6, 7, 8, 9, 16, 17, 18, 19, 26, 27, 28, 29, 36, 37, 38, 39, 46, 47, 48, 49, 56, 57, 58, 59, 66, 67, 68, 69, 76, 77, 78, 79, 86, 87, 88, 89, 96, 97, 98, 99]}

if MPI.COMM_WORLD.size != 2 :
    raise RuntimeError("Expected to be lanched with 2 procs !")

def test():
    pmesh = workPerProc[MPI.COMM_WORLD.rank]()
    pnodeids = mc.DataArrayInt(distribPerProc[MPI.COMM_WORLD.rank])
    cells = pmesh.getCellIdsLyingOnNodes(pnodeids,True)
    assert(cells.isEqual( mc.DataArrayInt(cellsExpected[MPI.COMM_WORLD.rank]) ))
    pmesh_red = pmesh.redistributeCells(cells)
    expected_mesh = GenerateGlobalMesh()[cells] ; expected_mesh.zipCoords()
    assert( expected_mesh.isEqual(pmesh_red.getMesh(),1e-12) )

if __name__ == "__main__":
    test()
