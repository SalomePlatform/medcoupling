# -*- coding: utf-8 -*-

###
### This script is intended to be launched in a new SALOME study
###

import sys
if sys.platform == "win32":
    from MEDCouplingCompat import *
else:
    from MEDCoupling import *

coordsArr=DataArrayDouble(range(3))
mesh1=MEDCouplingCMesh("mesh")
mesh1.setCoords(coordsArr,coordsArr,coordsArr[:1])
mesh1 = mesh1.buildUnstructured()

import math
r = 3.
nb = 10
a = 5.
coords = []
for i in range( nb ):
  x = r * math.cos( math.radians( i*a )) - r
  z = r * math.sin( math.radians( i*a ))
  coords.extend([ x, 0, z ])

m2=MEDCouplingCurveLinearMesh("myCurveLinearMesh")
m2.setNodeGridStructure([1,nb])
coo=DataArrayDouble(coords,nb,3)
m2.setCoords(coo)
m2 = m2.buildUnstructured()

m3=mesh1.buildExtrudedMesh(m2,1)

skin=m3.computeSkin()
skin.setName("skin")
ortho = skin.buildOrthogonalField()
ortho.setName("ortho_field")

from MEDLoader import WriteField, WriteMesh
medfile="mesh1.med"
WriteMesh(medfile,ortho.getMesh(),True)
WriteField(medfile,ortho,False)


import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

import iparameters
ipar = iparameters.IParameters(salome.myStudy.GetCommonParameters("Interface Applicative", 1), True)

#Set up visual properties:
ipar.setProperty("AP_ACTIVE_VIEW", "ParaView_0_0")
ipar.setProperty("AP_WORKSTACK_INFO", "00000001000000000000000201000000010000038f000000040000000100000001000000080000001800500061007200610056006900650077005f0030005f00300000000102")
ipar.setProperty("AP_ACTIVE_MODULE", "ParaViS")
ipar.setProperty("AP_SAVEPOINT_NAME", "GUI state: 1")
#Set up lists:
# fill list AP_VIEWERS_LIST
ipar.append("AP_VIEWERS_LIST", "ParaView_1")
# fill list ParaView_1
ipar.append("ParaView_1", "ParaView scene:7 - viewer:1")
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

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [897, 531]

# show data in view
mesh1medDisplay = Show(mesh1med, renderView1)

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# change representation type
mesh1medDisplay.SetRepresentationType('Surface With Edges')

# create a new 'Glyph'
glyph1 = Glyph(Input=mesh1med,GlyphType='Arrow')

# Properties modified on glyph1
glyph1.Vectors = ['CELLS', 'ortho_field']

# show data in view
glyph1Display = Show(glyph1, renderView1)

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on glyph1
glyph1.GlyphMode = 'All Points'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on glyph1
glyph1.ScaleFactor = 0.7071067811865466

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(mesh1med)

# Properties modified on renderView1
renderView1.Background = [1.0, 1.0, 1.0]

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [-9.40023601342297, -11.82197473872875, 6.269507101837803]
renderView1.CameraFocalPoint = [-3.714211157684638, -4.502747059961789, 3.6997577635620735]
renderView1.CameraViewUp = [0.6941544280720486, -0.3114076799016199, 0.6489798817268972]
renderView1.CameraParallelScale = 2.4893170029349183


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
  iparameters.getSession().restoreVisualState(1)
