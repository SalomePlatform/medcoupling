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

# This script illustrates how to create a matrix of pixels from data
# read in an image file. At the end, the cells of the matrix
# corresponds to the cells of a cartesian mesh that could hold a field
# whose value is the value of the pixel.
# (gboulant - 13/11/2011)

from PIL import Image
from PIL import ImageOps
import numpy

def image2matrix():
    # Load the image
    #img=Image.open("images/avatar.png")
    img=Image.open("images/tests.pgm")
    
    # Get a grayscale version
    imgbw=ImageOps.grayscale(img)
    
    # Save the image (optionnal)
    imgbw.save(fp="testsbw.pgm")
    
    # Get the data
    imgdata=imgbw.getdata()
    width,height=imgbw.size
    print list(imgdata)
    print width,height

    # Convert the data in a matrix using numpy
    tab=numpy.array(imgdata,dtype='float64')
    print list(tab)
    print tab
    nbRows=height
    nbCols=width
    matrix=numpy.reshape(tab,(nbRows,nbCols))
    # Note that in the reshape function, the height (sizeY) of the image
    # is specified first, because it corresponds to the number of rows.
    print matrix
    print list(matrix)

import MEDCoupling as MC
import MEDLoader as ML
def createMesh(meshname, sizeX, sizeY):
    """
    Creating a cartesian mesh with a grid of the size of the image.
    sizeX and sizeY should be respectively the width and heigth of the
    image.
    """
    # >>>
    # WARNING: remember the problem of tics and spaces. The data values
    # are considered as values defined on cells. With size values in a
    # direction, we have to create size+1 mesh nodes in that direction.
    # <<<
    
    # The mesh is created using MEDCoupling
    cmesh=MC.MEDCouplingCMesh.New();
    cmesh.setName(meshname)
    
    # We use an arbitrary step between cells (the value does not matter)
    stepX = 0.1
    nbNodesX = sizeX+1
    arrX = [float(i * stepX) for i in range(nbNodesX)]
    coordsX=MC.DataArrayDouble.New()
    coordsX.setValues(arrX,nbNodesX,1)

    # For the Y dimension, we have to reverse the coordinates (the
    # first pixel is at y=height and not at y=0).
    stepY = 0.1
    nbNodesY = sizeY+1
    lengthY = sizeY*stepY
    arrY=[float(lengthY - i * stepY) for i in range(nbNodesY)]
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
    
    return umesh

def createField(fieldname, mesh, image):
    """
    Creating a scalar field on the mesh using image data
    """    
    # Create the field using MEDCoupling
    field = MC.MEDCouplingFieldDouble.New(MC.ON_CELLS,MC.ONE_TIME);
    field.setName(fieldname);
    field.setMesh(mesh);
    # OPTIONAL: We set an arbitrary time step for test purpose
    field.setIteration(0);
    field.setOrder(0)

    imagedata=list(image.getdata())
    width,height=image.size
    nbCells = width*height
    dataArray=MC.DataArrayDouble.New();
    nbComponents=1 # For a scalar field
    
    dataArray.setValues(imagedata,nbCells,nbComponents)
    field.setArray(dataArray);
    
    return field

def image2med():
    img=Image.open("images/avatar.png")
    #img=Image.open("images/irm.png")
    imgbw=ImageOps.grayscale(img)
    # We keep only the grayscale. Maybe, it could be usefull to get
    # the RGB scales each on one component of the field.
    
    width,height=imgbw.size
    mesh=createMesh("mesh",width,height)
    field=createField("field",mesh,imgbw)
    
    createFromScratch=True
    ML.MEDLoader.WriteField("image.med",field,createFromScratch)


# ===================================================================
    
if __name__ == "__main__":
    #image2matrix()
    image2med()
