# -*- coding: utf-8 -*-

###
### This script is intended to be launched in a new SALOME study
###

import os
import salome

salome.salome_init()

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
ipar.append("ParaView_1", "ParaView scene:2 - viewer:1")
ipar.append("ParaView_1", "empty")
# fill list AP_MODULES_LIST
ipar.append("AP_MODULES_LIST", "ParaViS")

import sys
if sys.platform == "win32":
    from MEDCouplingCompat import *
else:
    from MEDCoupling import *

from MEDLoader import WriteMesh, WriteFieldUsingAlreadyWrittenMesh

medfile1="mesh1.med"
medfile2="mesh2.med"

m4=MEDCouplingCMesh("box")
coo=DataArrayDouble(range(7))
m4.setCoords(coo[:5],coo[:5],coo)
m4=m4.buildUnstructured()
valsArr1=m4.computeCellCenterOfMass()
valsArr1.applyFunc(1,"sqrt(X*X+Y*Y+Z*Z)")
field4 = MEDCouplingFieldDouble(ON_CELLS)
field4.setArray(valsArr1)
field4.setMesh(m4)
field4.setName("field4")

WriteMesh(medfile1,m4,True)
WriteFieldUsingAlreadyWrittenMesh(medfile1,field4)

origin=[0,0,2]
normvec=[-4,-4,6]
slice4=field4.extractSlice3D(origin,normvec,1e-10)
slice4.getMesh().translate([6,0,0])

WriteMesh(medfile2,slice4.getMesh(),True)
WriteFieldUsingAlreadyWrittenMesh(medfile2,slice4)



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
mesh1med = MEDReader(FileName=medfile1)

# Properties modified on mesh1med
mesh1med.AllArrays = ['TS0/box/ComSup0/box@@][@@P0', 'TS0/box/ComSup0/field4@@][@@P0']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [935, 531]

# show data in view
mesh1medDisplay = Show(mesh1med, renderView1)

# trace defaults for the display properties.
mesh1medDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(mesh1medDisplay, ('CELLS', 'field4', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
mesh1medDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
mesh1medDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'field4'
field4LUT = GetColorTransferFunction('field4')

# create a new 'MED Reader'
mesh2med = MEDReader(FileName=medfile2)

# Properties modified on mesh2med
mesh2med.AllArrays = ['TS0/Slice3D/ComSup0/Slice3D@@][@@P0', 'TS0/Slice3D/ComSup0/field4@@][@@P0']

# show data in view
mesh2medDisplay = Show(mesh2med, renderView1)

# trace defaults for the display properties.
mesh2medDisplay.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(mesh2medDisplay, ('CELLS', 'field4', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
mesh2medDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
mesh2medDisplay.SetScalarBarVisibility(renderView1, True)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [13.435468781360525, -42.4334252478128, 12.397177265244004]
renderView1.CameraFocalPoint = [11.468912675367923, -26.954985015354552, 9.183145534397903]
renderView1.CameraViewUp = [-0.0347897943004739, 0.19894659037092524, 0.9793926303542998]
renderView1.CameraParallelScale = 4.123105625617661


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
  iparameters.getSession().restoreVisualState(1)
