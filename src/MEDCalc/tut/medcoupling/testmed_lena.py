#!/usr/bin/env python
#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2011-2015  CEA/DEN, EDF R&D
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

# This use case illustrates the usage of PIL (Python Imaging Library)
# combined with the MEDCoupling and MEDLoader modules to save an image
# as a field in a med file.
# (gboulant - 27/06/2011)

import MEDCoupling as MC
import MEDLoader as ML

#
# ===============================================================
# We first get data from the test image to render as a field
# ===============================================================
#
import scipy, numpy
# The data field array may be created from the lena image
#image = scipy.lena()
# We could either read a real image using the PIL python package.
from scipy.misc import pilutil
image = pilutil.imread("images/avatar.png",True)


#from PIL import Image
#im=Image.open("images/irm.png")
#im=Image.open("images/lena.png")
#image=pilutil.fromimage(im,True)
#image=numpy.asarray(im)
#print image

dim  = len(image.shape)
print "Image space dimension = %d"%dim
sizeX = image.shape[1]
sizeY = image.shape[0]

# The sizes defined the number of pixel in a direction, then the
# number of cells to create in the mesh in that direction.

# We must reshape the matrix of pixel in a 1D vector that concatenates
# all the rows, and then convert this vector in a simple list of
# double as required by the MEDCoupling field specification.
import numpy
imageDataNArray       = image.reshape(1,sizeX*sizeY)[0]
print imageDataNArray

imageDataNArrayDouble = numpy.array(imageDataNArray, dtype='float64')
imageDataArrayDouble  = list(imageDataNArrayDouble)

#
# ===============================================================
# Creating a cartesian mesh with a grid of the size of the image
# ===============================================================
#

# >>>
# WARNING: remember the problem of tics and spaces. The data values
# are considered as values defined on cells. With size values in a
# direction, we have to create size+1 mesh nodes in that direction.
# <<<

# The mesh is created using MEDCoupling
cmesh=MC.MEDCouplingCMesh.New();
cmesh.setName("imagemesh")

# We use an arbitrary step between cells (the value does not matter)
stepX = 0.1
nbNodesX = sizeX+1
arrX = [float(i * stepX) for i in range(nbNodesX)]
coordsX=MC.DataArrayDouble.New()
coordsX.setValues(arrX,nbNodesX,1)

stepY = 0.1
nbNodesY = sizeY+1
arrY=[float(i * stepY) for i in range(nbNodesY)]
coordsY=MC.DataArrayDouble.New()
coordsY.setValues(arrY,nbNodesY,1)

cmesh.setCoords(coordsX,coordsY)
print "Imagem mesh dimension: %d"%cmesh.getSpaceDimension()

# WARN: In the current state of development of MEDLoader, only
# unstructured meshes are supported for writting function in med
# files. We just have to convert the cartesian mesh in an unstructured
# mesh before creating the field.
umesh=cmesh.buildUnstructured();
umesh.setName("imagemesh")

#
# ===============================================================
# Creating a scalar field on the mesh using image data
# ===============================================================
#

# Create the field using MEDCoupling
field = MC.MEDCouplingFieldDouble.New(MC.ON_CELLS,MC.ONE_TIME);
field.setName("imagefield");
field.setMesh(umesh);
# OPTIONAL: We set an arbitrary time step for test purpose
field.setIteration(3);
field.setOrder(0)

dataArray=MC.DataArrayDouble.New();
nbCells = sizeX*sizeY
nbComponents=1 # For a scalar field

# This example shows haw to initialize all cell with the same
# value. Just create an array of size nbCells
# dataArray.setValues(nbCells*[3.4],nbCells,nbComponents)

dataArray.setValues(imageDataArrayDouble,nbCells,nbComponents)
field.setArray(dataArray);

# The MEDLoader can be used to save all the stuff in a med file. You
# just have to specify the field and the MEDLoader will save the
# underlying mesh.
createFromScratch=True
ML.MEDLoader.WriteField("fieldimage.med",field,createFromScratch)
