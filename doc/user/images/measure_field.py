# -*- coding: utf-8 -*-

###
### This script is intended to be launched in a new SALOME study
###

import os
import salome

salome.salome_init()
theStudy = salome.myStudy

import sys
if sys.platform == "win32":
    from MEDCouplingCompat import *
else:
    from MEDCoupling import *

from MEDLoader import WriteMesh, WriteFieldUsingAlreadyWrittenMesh

medfile="mesh1.med"
XCoords=[-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22] # 9 values along X
YCoords=[0.,0.1,0.37,0.45,0.47,0.49,1.007] # 7 values along Y
arrX=DataArrayDouble(XCoords)
arrX.setInfoOnComponent(0,"X [m]")
arrY=DataArrayDouble(YCoords)
arrY.setInfoOnComponent(0,"Y [m]")
mesh=MEDCouplingCMesh("CMesh")
mesh.setCoords(arrX,arrY)
WriteMesh(medfile,mesh,True)

f=mesh.getMeasureField(True)
WriteFieldUsingAlreadyWrittenMesh(medfile,f)

import iparameters
ipar = iparameters.IParameters(salome.myStudy.GetCommonParameters("Interface Applicative", 1), True)

#Set up visual properties:
ipar.setProperty("AP_ACTIVE_VIEW", "ParaView_0_0")
ipar.setProperty("AP_WORKSTACK_INFO", "0000000100000000000000020100000001000003b5000000040000000100000001000000080000001800500061007200610056006900650077005f0030005f00300000000102")
ipar.setProperty("AP_ACTIVE_MODULE", "ParaViS")
ipar.setProperty("AP_SAVEPOINT_NAME", "GUI state: 1")
#Set up lists:
# fill list AP_VIEWERS_LIST
ipar.append("AP_VIEWERS_LIST", "ParaView_1")
# fill list ParaView_1
ipar.append("ParaView_1", "ParaView scene:33 - viewer:1")
ipar.append("ParaView_1", "empty")
# fill list AP_MODULES_LIST
ipar.append("AP_MODULES_LIST", "ParaViS")

###
### PARAVIS component
###

import pvsimple
pvsimple.ShowParaviewView()
#### import the simple module from the paraview
from pvsimple import *
#### disable automatic camera reset on 'Show'
pvsimple._DisableFirstRenderCameraReset()

# create a new 'MED Reader'
mesh1med = MEDReader(FileName=medfile)

# Properties modified on mesh1med
mesh1med.AllArrays = ['TS0/CMesh/ComSup0/MeasureOfMesh_CMesh@@][@@P0']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [935, 531]

# show data in view
mesh1medDisplay = Show(mesh1med, renderView1)

# trace defaults for the display properties.
mesh1medDisplay.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(mesh1medDisplay, ('CELLS', 'MeasureOfMesh_CMesh'))

# rescale color and/or opacity maps used to include current data range
mesh1medDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
mesh1medDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'MeasureOfMesh_CMesh'
measureOfMesh_CMeshLUT = GetColorTransferFunction('MeasureOfMesh_CMesh')

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [1.7629894333642653, 1.5037968367059755, -3.5223589503586297]
renderView1.CameraFocalPoint = [1.7629894333642653, 1.5037968367059755, 0.0]
renderView1.CameraParallelScale = 1.6150499279094872


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
  iparameters.getSession().restoreVisualState(1)
