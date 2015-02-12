#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2012-2015  CEA/DEN, EDF R&D
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

# This illustrates how to make a mesh partition using the value of a
# field defined on this mesh (for example to extract the cells where
# the field takes a value greater than a threshold L.
# (Anthony Geay, nov. 2012)

# WRN: this use case does not require a med input file because the
# data (mesh and field) to work with are created from scratch. 

from MEDCoupling import *

# =======================================================
# Creation of the input data (mesh and field) 
# =======================================================
#
# We prepare the input field from scratch instead of load it from a
# file, but there is no difference
from MEDCouplingDataForTest import MEDCouplingDataForTest
m3D=MEDCouplingDataForTest.build3DTargetMesh_1()
m3D.setName("m3D")
a=DataArrayDouble.New([1.,1.,1.,-1.,-1.,-1.,-1.,-1.])
field=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
field.setMesh(m3D)
field.setArray(a)
field.checkCoherency()
field.setName("f")

# Save the field (and associated mesh) 
from MEDLoader import MEDLoader
MEDLoader.WriteField("partition_input.med",field,True)

# =======================================================
# Determine the border skin mesh
# =======================================================
#
# We have to determine the 2D mesh that delimits the volume where the
# field is greater than a threshold L from the volume where the field
# is lower than this threshold (in this example L=0).
#
# WRN: This works in SALOME V660 only
#
# _T1A
L=0.
arr = field.getArray()
ids = arr.getIdsInRange(L,1e300)
m3DSub = field.getMesh()[ids]
skin = m3DSub.computeSkin()
MEDLoader.WriteUMesh("partition_skin.med",skin,True);
# _T1B

# =======================================================
# Compare to the result in SALOME V650
# =======================================================
# SALOME V650 requires a more complicated syntax.
m2D,desc,descI,revDesc,revDescI=m3DSub.buildDescendingConnectivity()
numberOf3DVolSharing=revDescI.deltaShiftIndex()
ids2D=numberOf3DVolSharing.getIdsEqual(1)
skin_V650=m2D[ids2D]
# We can check if the two skins are identical
print "Are two meshes equal between V660 and V650 ?",skin.isEqual(skin_V650,1e-12)
