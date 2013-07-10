#  -*- coding: iso-8859-1 -*-
# This illustrates the use of getValueOn in the case of hexahedron
# meshes (for which a temporary limitation implies the usage of the
# method simplexize that split the hexahedron into simplex,
# i.e. triangles and tetrahedrons).
#
# (gboulant, nov. 2012)
from MEDLoader import *
from MEDCoupling import *

import os
filename = "timeseries.med"
filepath = os.path.join(os.path.abspath(os.path.dirname(__file__)),filename)

rmedfilename = filepath

# Load the meshe data
meshname = "Grid_80x80"
fieldname = "Pulse"
dimrestriction = 0 # no restriction
msource = MEDLoader.ReadUMeshFromFile(rmedfilename,meshname,dimrestriction)

# WARN: In the current version of MEDCoupling, the getValueOn works
# only with simplex cells (triangles, tetrahedron). This is not a
# technical problem, but a question of specification of the
# interpolation to be performed in the case of other cells, in
# particular in the case of hexahedron meshes.
#
# A temporary solution (with good numerical results) is to split
# hexahedrons into simplex (before the association of the mesh to the
# field) using the method simplexize.
policy = 0
msource.simplexize(policy)

# Load the field data at iteration 0
iteration = 0
order = -1
fieldOnNodes = MEDLoader.ReadField(ON_NODES,rmedfilename,
                                   meshname,dimrestriction,
                                   fieldname,iteration,order)


fieldOnNodes.setMesh(msource)

# Get the value of field at coordinates x,y
x=0.5
y=0.5
fieldValue = fieldOnNodes.getValueOn([x,y])
print fieldValue
