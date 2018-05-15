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
medfile = os.path.join( os.getenv("MEDCOUPLING_ROOT_DIR"),"share","resources","med", "pointe.med")
pointemed = MEDReader(FileName=medfile)

# Properties modified on pointemed
pointemed.AllArrays = ['TS0/maa1/ComSup0/fieldnodeint@@][@@P1']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [935, 561]

# show data in view
pointemedDisplay = Show(pointemed, renderView1)

# trace defaults for the display properties.
pointemedDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(pointemedDisplay, ('POINTS', 'fieldnodeint', 'comp1'))

# rescale color and/or opacity maps used to include current data range
pointemedDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
pointemedDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'fieldnodeint'
fieldnodeintLUT = GetColorTransferFunction('fieldnodeint')

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [9.22339478264647, 20.564112667735554, 4.873867168433736]
renderView1.CameraFocalPoint = [3.2873912814202906, 7.329436317878308, 3.34609088265869]
renderView1.CameraViewUp = [0.08023261286174106, -0.14974521444076697, 0.9854639001939477]
renderView1.CameraParallelScale = 3.774917217635375


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
  iparameters.getSession().restoreVisualState(1)
