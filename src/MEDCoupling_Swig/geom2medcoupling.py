#! /usr/bin/env python3
#  -*- coding: utf-8 -*-
# Copyright (C) 2021  CEA/DEN, EDF R&D
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

import medcoupling as mc

def __to_geomshape_3D(mcmesh):
    """
    Precondition mcmesh is a MEDCouplingUMesh containing exactly one linear 3D cell.
    """
    import salome
    salome.standalone()
    salome.salome_init()
    import GEOM
    from salome.geom import geomBuilder
    geompy = geomBuilder.New()

    mcmesh2 = mcmesh.deepCopyConnectivityOnly()
    vertices = [ geompy.MakeVertex(*list(elt)) for elt in mcmesh2.getCoords()]

    mcfaces = mcmesh2.buildDescendingConnectivity()[0]
    shell_1 = geompy.MakeShell(
        [
            geompy.MakeFaceWires(
                [
                    geompy.MakePolyline([ vertices[nodeidx] for nodeidx in elt.getAllConn()[1:] ], True)
                ]
                , 1
            )
        for elt in mcfaces ]
        )
    return geompy.MakeSolid([shell_1])

def to_geomshape(mcmesh):
    """
    Method converting a unique 3D linear cell in a MEDCoupling mesh to GEOM Shape solid

    :param mcmesh: Mesh with single 3D linear cell.
    :type mcmesh: mc.MEDCouplingUMesh
    :return: GEOM shape
    :rtype: GEOM_Object
    """
    if not isinstance(mcmesh,mc.MEDCouplingUMesh):
        raise RuntimeError("Input mesh is expected to be of type mc.MEDCouplingUMesh !")
    if mcmesh.getNumberOfCells() != 1:
        raise RuntimeError("Input mesh is expected to contain exactly one cell !")
    if not mc.MEDCouplingUMesh.IsLinearGeometricType( mc.MEDCoupling1SGTUMesh(mcmesh).getCellModelEnum() ) :
        raise RuntimeError("The unique cell in the mesh is expected to be linear !")
    dico = { 3 : __to_geomshape_3D }
    mdim = mcmesh.getMeshDimension()
    if mdim not in dico:
        raise RuntimeError( "Input mesh is expected to have mesh dimension in {}".format(list(dico.keys())) )
    if mcmesh.getSpaceDimension() != 3:
        mcmesh.changeSpaceDimension(3,0.0)
    return (dico[mdim])(mcmesh)
