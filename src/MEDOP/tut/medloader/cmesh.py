#!/usr/bin/env python
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

from MEDLoader import MEDLoader

import os
filename = "madnex_field.med"
filepath = os.path.join(os.path.abspath(os.path.dirname(__file__)),filename)

meshName="my_field_RG"
fieldName="my_field"
itNumber=0
itOrder=0

# Load as an unstructured mesh
meshDimRelToMax = 0 # 0 = no restriction
umesh = MEDLoader.ReadUMeshFromFile(filepath,meshName,meshDimRelToMax)
print "umesh is structured: %s"%umesh.isStructured()

# Load as a structured mesh explicitly
# _T2A
from MEDLoader import MEDFileCMesh
medfile = MEDFileCMesh.New(filepath,meshName)
cmesh = medfile.getMesh()
# Note that the getMesh method is a short way to the method:
#cmesh = medfile.getGenMeshAtLevel(0,False)
print "cmesh is structured: %s"%cmesh.isStructured()
# _T2B

# Load and let MEDLoader decide what is nature of the mesh
# _T1A
from MEDLoader import MEDFileMesh
medfile = MEDFileMesh.New(filepath,meshName)
print medfile.advancedRepr()
meshDimRelToMax = 0 # 0 = no restriction
mesh = medfile.getGenMeshAtLevel(meshDimRelToMax)
print "mesh is structured: %s"%mesh.isStructured()
# _T1B


# Write the mesh to another file
# _T3A
outputfilepath="output.med"
mode=0
medfile.write(outputfilepath,mode)
# _T3B

# test to reload the mesh
medfile = MEDFileCMesh.New(outputfilepath,meshName)
cmesh = medfile.getMesh()
print "cmesh is structured: %s"%cmesh.isStructured()

# Q: Is it possible to know if a mesh is structured or unstructured
# without loading the mesh.
