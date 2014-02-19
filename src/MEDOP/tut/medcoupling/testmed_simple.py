#!/usr/bin/env python
#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2011-2014  CEA/DEN, EDF R&D
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

# This simple use case illustrates the basic usage of MEDCoupling and
# MEDLoader to create a cartesian mesh, define a field on this mesh,
# and save all the stuff in a med file.
# (gboulant - 27/06/2011)

import MEDCoupling as MC
import MEDLoader as ML

#
# ===============================================================
# Creating a 512x512 cartesian mesh
# ===============================================================
#
# The size is the number of discrete values in a direction, and then
# corresponds to the number of cells in that direction.
size=8
#size=512

# The mesh is created using MEDCoupling. The code below creates a
# cartesian mesh as a sizexsize grid

# >>>
# WARNING: remember the problem of tics and spaces. The data values
# are considered as values defined on cells. With size values in a
# direction, we have to create size+1 mesh nodes in that direction.
# <<<

cmesh=MC.MEDCouplingCMesh.New();
cmesh.setName("512x512 cartesian mesh")

sizeX = size
nbNodesX = sizeX+1
stepX = 0.1
arrX = [float(i * stepX) for i in range(nbNodesX)]
print "Size of arrX = %d"%len(arrX)

coordsX=MC.DataArrayDouble.New()
coordsX.setValues(arrX,nbNodesX,1)

sizeY = size
nbNodesY = sizeY+1
stepY = 0.1
arrY=[float(i * stepY) for i in range(nbNodesY)]
coordsY=MC.DataArrayDouble.New()
coordsY.setValues(arrY,sizeY,1)

cmesh.setCoords(coordsX,coordsY)
print cmesh.getSpaceDimension()
#print cmesh

# WARN: In the current state of development of MEDLoader, only
# unstructured meshes are supported for writting function in med
# files. We just have to convert the cartesian mesh in an unstructured
# mesh before creating the field.
umesh=cmesh.buildUnstructured();
umesh.setName("512x512 unstructured mesh")

# This can be used to save the mesh only (can be visualize using
# SMESH).
meshFileName = "umesh.med"
ML.MEDLoader.WriteUMesh(meshFileName,umesh,True);

# Alternatively, you can use a MEDFileMesh to write the mesh in a
# file.
medFileCMesh = ML.MEDFileCMesh.New()
medFileCMesh.setMesh(cmesh)
medFileCMesh.setName(cmesh.getName())
meshFileName = "cmesh.med"
mode = 2
medFileCMesh.write(meshFileName,mode)

#
# ===============================================================
# Creating a scalar field on the 512x512 mesh
# ===============================================================
#
# For the simple test, we create a field that varies in space as
# field(x,y)=x+y where x and y are coordinates on the mesh

# --- Field on cells

# Create the field
field = MC.MEDCouplingFieldDouble.New(MC.ON_CELLS);
field.setName("AnalyticField_onCells");
field.setMesh(umesh);

nbComponents=1 # Only one single component for a scalar field
fillFunction="x+y"
field.fillFromAnalytic(nbComponents,fillFunction);

# The MEDLoader can be used to save all the stuff in a med file. You
# just have to specify the field and the MEDLoader will save the
# underlying mesh.
createFromScratch=True
ML.MEDLoader.WriteField("fieldtest.med",field,createFromScratch)

# --- Field on nodes

field = MC.MEDCouplingFieldDouble.New(MC.ON_NODES);
field.setName("AnalyticField_onNodes");
field.setMesh(umesh);
field.fillFromAnalytic(nbComponents,fillFunction);
createFromScratch=False
ML.MEDLoader.WriteField("fieldtest.med",field,createFromScratch)


#
# ===============================================================
# Creating a scalar field, working with numpy
# ===============================================================
#

# We start by creating a numpy matrix
import numpy
rows=[]
for irow in range(sizeY):
    row = numpy.arange(irow*sizeY,irow*sizeY+sizeX,dtype='float64')
    rows.append(row)

marray = numpy.vstack(rows)

# Then, we can reshape the matrix in a 1D vector that concatenate all
# the rows
data=marray.reshape(1,sizeX*sizeY)[0]
# Finally, we can create a simple list as required by the MEDCoupling
# DataArrayDouble. Note also the usage of float type because
# MEDCoupling works only with real numbers
listdata=list(data)

# Create the field using the list obtained from the numpy array
fieldWithNumpy = MC.MEDCouplingFieldDouble.New(MC.ON_CELLS);
fieldWithNumpy.setName("Numpy Field");
fieldWithNumpy.setMesh(umesh);

nbCells=sizeX*sizeY
dataArray=MC.DataArrayDouble.New();
dataArray.setValues(listdata,nbCells,nbComponents)
fieldWithNumpy.setArray(dataArray);

createFromScratch=False
ML.MEDLoader.WriteField("fieldtest.med",fieldWithNumpy,createFromScratch)


