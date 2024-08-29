#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2023-2024  CEA, EDF
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

def MEDFileUMeshReduceToCells(self, level, keepCells, removeOrphanNodes=True):
    """
    Method returning a new MEDFileUMesh, restriction of self to level and keepCell cells at this level.
    This method also 

    :param level: Specifies the top level of the returned MEDFileUMesh expected
    :param keepCells: A DataArrayInt specifying cell ids at level level of self
    :param removeOrphanNodes: Specifies if orphan nodes should be removed at the end
    
    see also MEDFileUMesh.extractPart
    """
    import MEDLoader as ml
    subLevs = [l for l in self.getNonEmptyLevels() if l<=level]
    subMeshes = [self[lev] for lev in subLevs]
    allFamilyFields = [self.getFamilyFieldAtLevel(lev) for lev in subLevs]
    allRefMesh = subMeshes[0]
    refMesh = allRefMesh[keepCells]

    mmOut = ml.MEDFileUMesh()
    # level 0
    mmOut[0] = refMesh
    mmOut.setFamilyFieldArr(0,allFamilyFields[0][keepCells])

    # subLevels
    for curLev,meshLev,famFieldLev in zip(subLevs[1:],subMeshes[1:],allFamilyFields[1:]):
        allMeshLev,d,di, rd,rdi = allRefMesh.explodeMeshTo( curLev-level )
        a,b = allMeshLev.areCellsIncludedIn(meshLev,2)
        if not a:
            raise RuntimeError("Error in mesh {}")
        dlev,dlevi = ml.DataArrayInt.ExtractFromIndexedArrays( keepCells, d,di )
        dlev2 = dlev.buildUniqueNotSorted()
        cellsToKeepLev = ml.DataArrayInt.BuildIntersection([dlev2,b])
        cellsToKeepLev = b.indicesOfSubPart(cellsToKeepLev)
        cellsToKeepLev.sort()
        mmOut[curLev] = meshLev[cellsToKeepLev]
        mmOut.setFamilyFieldArr(curLev,famFieldLev[cellsToKeepLev])

    allFamNodes = mmOut.getFamilyFieldAtLevel(1)
    if allFamNodes:
        mmOut.setFamilyFieldArr(1,allFamNodes[:])

    if removeOrphanNodes:
        mmOut.zipCoords()

    mmOut.copyFamGrpMapsFrom(self)
    return mmOut
