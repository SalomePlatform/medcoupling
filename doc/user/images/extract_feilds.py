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

medfile = os.path.join( os.getenv("MEDCOUPLING_ROOT_DIR"),"share","resources","med", "pointe.med")

from MEDLoader import ReadField, WriteField, WriteMesh
f=ReadField(medfile,"fieldnodeint") # field on 19 nodes
f4 = f[(range(19/2))]
f4.getMesh().translate( [5,0,0] )

import tempfile
medfile2=tempfile.NamedTemporaryFile().name + ".med"

WriteMesh(medfile2,f4.getMesh(), True )
WriteField(medfile2,f4,False)

# ###
### PARAVIS component
###

import pvsimple
pvsimple.ShowParaviewView()
#### import the simple module from the paraview
from pvsimple import *
#### disable automatic camera reset on 'Show'
pvsimple._DisableFirstRenderCameraReset()

# create a new 'MED Reader'
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

# set scalar coloring
ColorBy(pointemedDisplay, ('POINTS', 'fieldnodeint', 'comp1'))

# rescale color and/or opacity maps used to include current data range
pointemedDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
pointemedDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'fieldnodeint'
fieldnodeintLUT = GetColorTransferFunction('fieldnodeint')


# create a new 'MED Reader'
partmed = MEDReader(FileName='/data/eap/S8/MEDCOUPLING_BUILD/doc/part.med')

# Properties modified on partmed
partmed.AllArrays = ['TS0/maa1/ComSup0/fieldnodeint@@][@@P1', 'TS0/maa1/ComSup0/maa1@@][@@P0']

# show data in view
partmedDisplay = Show(partmed, renderView1)

# trace defaults for the display properties.
partmedDisplay.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(partmed, renderView1)

# destroy partmed
Delete(partmed)
del partmed

# create a new 'MED Reader'
partmed = MEDReader(FileName='/data/eap/S8/MEDCOUPLING_BUILD/doc/part.med')

# Properties modified on partmed
partmed.AllArrays = ['TS0/maa1/ComSup0/fieldnodeint@@][@@P1', 'TS0/maa1/ComSup0/maa1@@][@@P0']

# show data in view
partmedDisplay = Show(partmed, renderView1)

# trace defaults for the display properties.
partmedDisplay.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# destroy partmed
Delete(partmed)
del partmed

# create a new 'MED Reader'
partmed = MEDReader(FileName=medfile2)

# Properties modified on partmed
partmed.AllArrays = ['TS0/maa1/ComSup0/fieldnodeint@@][@@P1']

# show data in view
partmedDisplay = Show(partmed, renderView1)

# trace defaults for the display properties.
partmedDisplay.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(partmedDisplay, ('POINTS', 'fieldnodeint', 'comp1'))

# rescale color and/or opacity maps used to include current data range
partmedDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
partmedDisplay.SetScalarBarVisibility(renderView1, True)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [10.652899780874208, -21.421142047941988, 4.509321505454196]
renderView1.CameraFocalPoint = [-1.941979639934524, 14.064247480193199, 0.8656523817352476]
renderView1.CameraViewUp = [-0.033512904575591126, 0.09030855881511615, 0.9953498125943683]
renderView1.CameraParallelScale = 3.774917217635375


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
  iparameters.getSession().restoreVisualState(1)
