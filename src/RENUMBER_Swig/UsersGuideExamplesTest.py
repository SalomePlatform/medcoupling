#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2019  CEA/DEN, EDF R&D
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

import sys

if sys.platform == "win32":
    from MEDCouplingCompat import *
else:
    from medcoupling import *

from math import pi, sqrt

if sys.platform == "win32":
    import MEDCouplingCompat as MEDCoupling
else:
    import MEDCoupling

from MEDCouplingDataForTest import MEDCouplingDataForTest
m=MEDCouplingDataForTest.build2DTargetMesh_1();
#! [UG_Optimization_0]
from MEDRenumber import RenumberingFactory
ren=RenumberingFactory("BOOST")
a,b=m.computeNeighborsOfCells()
n2o,_=ren.renumber(a,b)
mrenum=m[n2o]
#! [UG_Optimization_0]

#! [UG_Optimization_1]
import MEDPartitioner
# prepare a MEDPartitioner
a,b=m.computeNeighborsOfCells()
sk=MEDCouplingSkyLineArray(b,a)
g=MEDPartitioner.MEDPartitioner.Graph(sk)
# compute partitioning into 4 parts
g.partGraph(4)
# get the 1st of parts of m
procIdOnCells=g.getPartition().getValuesArray()
p0=procIdOnCells.findIdsEqual(0)
part0=m[p0]
#! [UG_Optimization_1]
#! [UG_Optimization_2]
boundary_nodes_part0=part0.findBoundaryNodes()
boundary_cells_part0=p0[part0.getCellIdsLyingOnNodes(boundary_nodes_part0,False)]
# starting from knowledge of neighborhood it s possible to know the neighbors of boundary_cells_part0
neighbors_boundary_cells_part0=MEDCouplingUMesh.ExtractFromIndexedArrays(boundary_cells_part0,a,b)[0]
neighbors_boundary_cells_part0.sort()
neighbors_boundary_cells_part0=neighbors_boundary_cells_part0.buildUnique()
#
layer_of_part0=neighbors_boundary_cells_part0.buildSubstraction(p0)
#
whole_part_with_layer=DataArrayInt.Aggregate([p0,layer_of_part0])
whole_part_with_layer.sort()
part0_with_layer=m[whole_part_with_layer]
#! [UG_Optimization_2]

