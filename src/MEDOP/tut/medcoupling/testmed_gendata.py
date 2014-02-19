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

# This script illustrates the basic usage of MEDCoupling and MEDLoader
# to generate test data files for various cases of med operation. It
# illustrates also the usage of numpy to specify the values of the
# fields when defined on a cartesian mesh (grid).
# (gboulant - 11/07/2011)

import MEDCoupling as MC
import MEDLoader as ML

import numpy

#
# ===============================================================
# Helper functions to create meshes
# ===============================================================
#

def createGridMesh(meshName, nbCellsX, nbCellsY):
    """
    The mesh is created using MEDCoupling. The code below creates a
    cartesian mesh as a grid with nbCellsX segments in the X direction
    and nbCellsY in the Y direction (nb. cells = nbCellsX * nbCellsY)
    """
    print "Creating grid mesh of size %sx%s"%(nbCellsX, nbCellsY)
    cmesh=MC.MEDCouplingCMesh.New();

    # Create X coordinates
    nbNodesX = nbCellsX+1
    stepX = 0.1
    arrX = [float(i * stepX) for i in range(nbNodesX)]    
    coordsX=MC.DataArrayDouble.New()
    coordsX.setValues(arrX,nbNodesX,1)

    # Create Y coordinates
    nbNodesY = nbCellsY+1
    stepY = 0.1
    arrY=[float(i * stepY) for i in range(nbNodesY)]
    coordsY=MC.DataArrayDouble.New()
    coordsY.setValues(arrY,nbNodesY,1)

    # Create the grid
    cmesh.setCoords(coordsX,coordsY)
    cmesh.setName(meshName)

    return cmesh

def unstructuredMesh(cartesianMesh):
    """
    Convert the cartesian mesh in unstructured mesh for the need of
    write function of MEDLoader
    """
    print "Creating unstructured mesh from %s"%(cartesianMesh.getName())
    umesh=cartesianMesh.buildUnstructured();
    umesh.setName(cartesianMesh.getName())
    return umesh

#
# ===============================================================
# Creating a cartesian mesh
# ===============================================================
#
# The size is the number of discrete values in a direction, and then
# corresponds to the number of cells in that direction.
size=80
#size=512


# >>>
# WARNING: remember the problem of tics and spaces. The parameter
# "size" is considered to be a number of cells (intervals). The number
# of nodes in that direction is size+1.
# <<<

nbCellsX = size
nbNodesX = nbCellsX+1

nbCellsY = size # The size could be different than the X size
nbNodesY = nbCellsY+1

meshName = "Grid_%sx%s"%(nbCellsX, nbCellsY)
cmesh = createGridMesh(meshName, nbCellsX, nbCellsY)
umesh = unstructuredMesh(cmesh)
medFileName="gendata.med"
ML.MEDLoader.WriteUMesh(medFileName,umesh,True);

#
# ===============================================================
# Creating a scalar field, working with numpy
# ===============================================================
#

def createField(fieldName,gridMesh,
                numpy2Darray,typeOfField=MC.ON_CELLS,
                iteration=0):
    """
    The number of values for the fields is deduced from the sizes of
    the numpy array. If typeOfField is ON_CELLS, the size is considered
    as the number of cells, otherwise it's considered as the number of
    nodes. In any case, it must be consistent with the dimensions of
    the numpy 2D array.
    """
    print "Creating field %s with iteration=%s"%(fieldName,iteration)

    # The sizes are deduced from the numpy array. Note that if
    # typeOfField is ON_CELLS, then the size should correspond to the
    # number of cells, while if typeOfField is ON_NODES, then the size
    # should correspond to the number of nodes
    [sizeX,sizeY] = numpy2Darray.shape

    # We first have to reshape the 2D numpy array in a 1D vector that
    # concatenate all the rows
    data=numpy2Darray.reshape(1,sizeX*sizeY)[0]
    # Then, we can create a simple list as required by the MEDCoupling
    # DataArrayDouble. Note also the usage of float type because
    # MEDCoupling works only with real numbers
    listdata=list(data)
    
    # Create the field using the list obtained from the numpy array
    field = MC.MEDCouplingFieldDouble.New(typeOfField,MC.ONE_TIME);
    field.setName(fieldName);
    field.setMesh(gridMesh);
    field.setIteration(iteration)
    field.setTimeValue(float(iteration))
    
    nbComponents=1 # Only one single component for a scalar field
    nbCells=sizeX*sizeY
    dataArray=MC.DataArrayDouble.New();
    dataArray.setValues(listdata,nbCells,nbComponents)
    field.setArray(dataArray);

    return field

def writeField(fieldName, numpy2Darray,
               typeOfField=MC.ON_CELLS,
               iteration=0):

    field = createField(fieldName, umesh, numpy2Darray,
                        typeOfField, iteration)
    createFromScratch=False
    ML.MEDLoader.WriteField(medFileName,field,createFromScratch)
    

def createTestNumpy2DArray(sizeX, sizeY):
    """
    This illustrates how to create a numpy 2D array for input of the
    createField function.
    """
    rows=[]
    for irow in range(sizeY):
        row = numpy.arange(start = irow*sizeY,
                           stop  = irow*sizeY+sizeX,
                           step  = 1,
                           dtype='float64')
        rows.append(row)

    numpy2Darray = numpy.vstack(rows)
    return numpy2Darray
    
def createTestFieldOnCells():
    # Test field on cells
    numpy2Darray = createTestNumpy2DArray(sizeX=nbCellsX, sizeY=nbCellsY)
    writeField("FieldOnCells", numpy2Darray,
               typeOfField=MC.ON_CELLS)

def createTestFieldOnNodes():
    # Test field on nodes
    numpy2Darray = createTestNumpy2DArray(sizeX=nbNodesX, sizeY=nbNodesY)
    writeField("FieldOnNodes", numpy2Darray,
               typeOfField=MC.ON_NODES)
    

#
# =================================================
# Creating a time series
# =================================================
#

# -------------------------------------------------
# Simple demo of the principles
# -------------------------------------------------

# In these functions, (x,y) are the indexes of the element in the
# numpy array. Note that theses indexes maps the indexes of the
# cartesian mesh.

# A function can be a simple python function ...
def f1(x,y):
    z = 10*x
    print "x=%s\ny=%s\nz=%s"%(x,y,z)
    return z

# ... but also a more sophisticated callable object, for example to
# defines some parameters
class Function(object):
    def __init__(self, sizeX, sizeY, param):
        self.sizeX = sizeX
        self.sizeY = sizeY
        self.param = param

    def function(self, x,y):
        z = self.param*x
        print "x=%s\ny=%s\nz=%s"%(x,y,z)
        return z

    def __call__(self, x,y):
        return self.function(x,y)

fOnNodes=Function(sizeX=nbNodesX, sizeY=nbNodesY, param=10)
fOnCells=Function(sizeX=nbCellsX, sizeY=nbCellsY, param=3)

def createFunctionField_01():
    sizeX=nbNodesX
    sizeY=nbNodesY
    typeOfField=MC.ON_NODES
    f=fOnNodes
    numpy2Darray = numpy.fromfunction(f,(sizeX,sizeY),dtype='float64')
    writeField("FieldOnNodesUsingFunc", numpy2Darray,typeOfField)

    f=fOnCells
    sizeX=nbCellsX
    sizeY=nbCellsY
    typeOfField=MC.ON_CELLS
    numpy2Darray = numpy.fromfunction(f,(sizeX,sizeY),dtype='float64')
    writeField("FieldOnCellsUsingFunc", numpy2Darray,typeOfField)


# -------------------------------------------------
# Using the pyfunctions package to generate data
# -------------------------------------------------

def createNumpy2DArrayWithFunc(sizeX, sizeY, function):
    """
    @function : a callable than can be used as a function of X.
    Typically function should be an instance of Function object
    defined in pyfunctions.functions.
    """

    # X coordinates should range between 0 and 1 to use the normalized
    # functions. We have to generate sizeX points:
    step=1./sizeX
    arrX=[float(i * step) for i in range(sizeX)]

    values = function(arrX)

    # Then on can create the base row for the numpy 2D array
    rowX = numpy.array(values)
    # and replicate this row along the Y axis
    rows=[]
    for irow in range(sizeY):
        rows.append(rowX)

    numpy2Darray = numpy.vstack(rows)
    return numpy2Darray

from pyfunctions.functions import FuncStiffPulse
def createNumpy2DArrayWithFuncStiff(sizeX, sizeY):
    f=FuncStiffPulse(xlimit=0.3,stiffness=30,nbPeriods=10)
    return createNumpy2DArrayWithFunc(sizeX, sizeY, f)
    
def createFunctionField_02():
    sizeX=nbCellsX
    sizeY=nbCellsY
    typeOfField=MC.ON_CELLS
    numpy2Darray = createNumpy2DArrayWithFuncStiff(sizeX,sizeY)
    writeField("FieldOnCellsUsingFunc02", numpy2Darray,typeOfField)

    sizeX=nbNodesX
    sizeY=nbNodesY
    typeOfField=MC.ON_NODES
    numpy2Darray = createNumpy2DArrayWithFuncStiff(sizeX,sizeY)
    writeField("FieldOnNodesUsingFunc02", numpy2Darray,typeOfField)

#
# =================================================
# Functions to create custom fields for MEDOP tests
# =================================================
#
def createTimeSeries():
    """
    Create a single med file with a single mesh and a field defined on
    several time steps (time series).
    """
    meshName = "Grid_%sx%s"%(nbCellsX, nbCellsY)
    cmesh = createGridMesh(meshName, nbCellsX, nbCellsY)
    umesh = unstructuredMesh(cmesh)
    medFileName="timeseries.med"
    ML.MEDLoader.WriteUMesh(medFileName,umesh,True);

    sizeX=nbNodesX
    sizeY=nbNodesY
    typeOfField=MC.ON_NODES

    nbIterations=10
    pulseStiffNess = 20
    pulseNbPeriods = 10
    for iteration in range(nbIterations):
        xlimit = float(iteration)/float(nbIterations)
        f=FuncStiffPulse(xlimit,stiffness=pulseStiffNess,nbPeriods=pulseNbPeriods)
        numpy2Darray = createNumpy2DArrayWithFunc(sizeX,sizeY,f)
        field = createField("Pulse",umesh,numpy2Darray,typeOfField,iteration)
        ML.MEDLoader.WriteField(medFileName,field,False)

from pyfunctions.functions import FuncStiffExp
def createParametrics():
    """
    Create 2 med files containing each a mesh (identical) and a field
    defined on this mesh in each file.
    """
    meshName = "Grid_%sx%s_01"%(nbCellsX, nbCellsY)
    cmesh = createGridMesh(meshName, nbCellsX, nbCellsY)
    umesh = unstructuredMesh(cmesh)
    
    sizeX=nbNodesX
    sizeY=nbNodesY
    typeOfField=MC.ON_NODES

    medFileName="parametric_01.med"
    ML.MEDLoader.WriteUMesh(medFileName,umesh,True);
    f=FuncStiffExp(xlimit=0.3,stiffness=30)
    numpy2Darray = createNumpy2DArrayWithFunc(sizeX,sizeY,f)
    fieldName = "StiffExp_01"
    field = createField(fieldName,umesh, numpy2Darray,typeOfField)
    ML.MEDLoader.WriteField(medFileName,field,False)

    medFileName="parametric_02.med"
    umesh.setName("Grid_%sx%s_02"%(nbCellsX, nbCellsY))
    ML.MEDLoader.WriteUMesh(medFileName,umesh,True);
    f=FuncStiffExp(xlimit=0.4,stiffness=30)
    numpy2Darray = createNumpy2DArrayWithFunc(sizeX,sizeY,f)
    fieldName = "StiffExp_02"
    field = createField(fieldName,umesh, numpy2Darray,typeOfField)
    ML.MEDLoader.WriteField(medFileName,field,False)

def createParametrics_demo():
    """
    Create 2 med files containing each a mesh (identical) and a field
    defined on this mesh in each file.
    """
    meshName = "mesh1"
    cmesh = createGridMesh(meshName, nbCellsX, nbCellsY)
    umesh = unstructuredMesh(cmesh)
    
    sizeX=nbNodesX
    sizeY=nbNodesY
    typeOfField=MC.ON_NODES

    listIteration = [0,1,2,3,4]

    medFileName="parametric_01.med"
    ML.MEDLoader.WriteUMesh(medFileName,umesh,True);
    fieldName = "field1"
    for iteration in listIteration:
        #f=FuncStiffPulse(xlimit=0.3+0.1*iteration,stiffness=10,nbPeriods=5)
        f=FuncStiffExp(xlimit=0.3+0.1*iteration,stiffness=10)
        numpy2Darray = createNumpy2DArrayWithFunc(sizeX,sizeY,f)
        field = createField(fieldName,umesh, numpy2Darray,typeOfField,iteration)
        ML.MEDLoader.WriteField(medFileName,field,False)
    
    medFileName="parametric_02.med"
    umesh.setName("mesh2")
    ML.MEDLoader.WriteUMesh(medFileName,umesh,True);
    fieldName = "field2"
    for iteration in listIteration:
        #f=FuncStiffPulse(xlimit=0.3+0.1*iteration,stiffness=10,nbPeriods=6)
        f=FuncStiffExp(xlimit=0.3+0.1*iteration,stiffness=15)
        numpy2Darray = createNumpy2DArrayWithFunc(sizeX,sizeY,f)
        field = createField(fieldName,umesh, numpy2Darray,typeOfField,iteration)
        ML.MEDLoader.WriteField(medFileName,field,False)



#
# =================================================
# Main runner
# =================================================
#
if __name__ == "__main__":
    #createTestFieldOnCells()
    #createTestFieldOnNodes()
    #createFunctionField_01()
    #createFunctionField_02()
    #createTimeSeries()
    createParametrics_demo()
