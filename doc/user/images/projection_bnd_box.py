# -*- coding: utf-8 -*-

###
### This script is intended to be launched in a new SALOME study
###

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

import iparameters
ipar = iparameters.IParameters(salome.myStudy.GetCommonParameters("Interface Applicative", 1), True)

#Set up visual properties:
ipar.setProperty("AP_ACTIVE_VIEW", "OCCViewer_0_0")
ipar.setProperty("AP_WORKSTACK_INFO", "000000010000000000000002010000000100000426000000040000000100000001000000080000001a004f00430043005600690065007700650072005f0030005f00300000000102")
ipar.setProperty("AP_ACTIVE_MODULE", "Geometry")
ipar.setProperty("AP_SAVEPOINT_NAME", "GUI state: 1")
#Set up lists:
# fill list AP_VIEWERS_LIST
ipar.append("AP_VIEWERS_LIST", "OCCViewer_1")
# fill list OCCViewer_1
ipar.append("OCCViewer_1", "OCC scene:1 - viewer:1")
ipar.append("OCCViewer_1", "0|-1|scale=3.425984215998e+2*projX=4.501122113466e-1*projY=-8.801114765472e-1*projZ=1.510059139458e-1*twist=6.274342135079e+0*atX=2.565570424306e-1*atY=3.358900466217e-1*atZ=5.981017875427e-2*eyeX=3.015682642359e-1*eyeY=2.478788976555e-1*eyeZ=7.491077037387e-2*scaleX=1.000000000000e+0*scaleY=1.000000000000e+0*scaleZ=1.000000000000e+0*isVisible=0*size=1.30*gtIsVisible=0*gtDrawNameX=1*gtDrawNameY=1*gtDrawNameZ=1*gtNameX=X*gtNameY=Y*gtNameZ=Z*gtNameColorRX=255*gtNameColorGX=0*gtNameColorBX=0*gtNameColorRY=0*gtNameColorGY=255*gtNameColorBY=0*gtNameColorRZ=0*gtNameColorGZ=0*gtNameColorBZ=255*gtDrawValuesX=1*gtDrawValuesY=1*gtDrawValuesZ=1*gtNbValuesX=3*gtNbValuesY=3*gtNbValuesZ=3*gtOffsetX=2*gtOffsetY=2*gtOffsetZ=2*gtColorRX=255*gtColorGX=0*gtColorBX=0*gtColorRY=0*gtColorGY=255*gtColorBY=0*gtColorRZ=0*gtColorGZ=0*gtColorBZ=255*gtDrawTickmarksX=1*gtDrawTickmarksY=1*gtDrawTickmarksZ=1*gtTickmarkLengthX=5*gtTickmarkLengthY=5*gtTickmarkLengthZ=5*lightSource=lightType~1;lightX~0;lightY~0;lightZ~-1;lightColorR~1;lightColorG~1;lightColorB~1;lightHeadlight~1;*background=bt$1;fn$;tm$0;ts$false;c1$#ffffff;c2$#698fff;gt$1;gr$")
# fill list AP_MODULES_LIST
ipar.append("AP_MODULES_LIST", "Geometry")

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
sk = geompy.Sketcher2D()
sk.addPoint(0.000000, 0.000000)
sk.addSegmentAbsolute(1.000000, 0.000000)
sk.addSegmentAbsolute(0.000000, 1.000000)
sk.close()
geomObj_1 = geompy.MakeMarker(0, 0, 0, 1, 0, 0, 0, 1, 0)
Sketch_1 = sk.wire(geomObj_1)
Triangle_src = geompy.MakeFaceWires([Sketch_1], 1)

Triangle_tgt = geompy.MakeTranslation(Triangle_src, 0, 0, 0.1)

sk = geompy.Sketcher2D()
sk.addPoint(-0.150000, -0.150000)
sk.addSegmentAbsolute(1.150000, -0.150000)
sk.addSegmentAbsolute(1.150000, 1.150000)
sk.addSegmentAbsolute(-0.150000, 1.150000)
sk.close()
geomObj_4 = geompy.MakeMarker(0, 0, 0, 1, 0, 0, 0, 1, 0)
Sketch_2 = sk.wire(geomObj_4)

Bnd_box = geompy.MakePrismVecH(Sketch_2, OZ, 0.3)
geompy.TranslateDXDYDZ(Bnd_box, 0, 0, -0.15)

Bnd_box.SetColor(SALOMEDS.Color(1,0,0))
Triangle_src.SetColor(SALOMEDS.Color(0,0.333333,1))
Triangle_tgt.SetColor(SALOMEDS.Color(0,1,0.498039))

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Sketch_1, 'Sketch_1' )
geompy.addToStudy( Triangle_src, 'Triangle_src' )
geompy.addToStudy( Triangle_tgt, 'Triangle_tgt' )
geompy.addToStudy( Sketch_2, 'Sketch_2' )
geompy.addToStudy( Bnd_box, 'Bnd_box' )

### Store presentation parameters of displayed objects
import iparameters
ipar = iparameters.IParameters(theStudy.GetModuleParameters("Interface Applicative", "GEOM", 1))

#Set up entries:
# set up entry GEOM_6 (Triangle_src) parameters
objId = geompy.getObjectID(Triangle_src)
ipar.setParameter(objId, "OCCViewer_0_Visibility", "On")
ipar.setParameter(objId, "OCCViewer_0_DisplayMode", "2")
ipar.setParameter(objId, "OCCViewer_0_Color", "0:0.333333:1")
ipar.setParameter(objId, "OCCViewer_0_Transparency", "0")
ipar.setParameter(objId, "OCCViewer_0_TopLevelFlag", "false")
ipar.setParameter(objId, "OCCViewer_0_Isos", "0:0")
ipar.setParameter(objId, "OCCViewer_0_VectorMode", "false")
ipar.setParameter(objId, "OCCViewer_0_VerticesMode", "false")
ipar.setParameter(objId, "OCCViewer_0_NameMode", "false")
ipar.setParameter(objId, "OCCViewer_0_DeflectionCoeff", "0.001")
ipar.setParameter(objId, "OCCViewer_0_MarkerType", "7:3")
ipar.setParameter(objId, "OCCViewer_0_Material", "Physical=0:FrontShininess=0.13:BackShininess=0.13:Transparency=0:Ambient=1:AmbientColor=#333333:FrontAmbientCoefficient=0.3:BackAmbientCoefficient=0.25:Diffuse=1:DiffuseColor=#000000:FrontDiffuseCoefficient=0.5:BackDiffuseCoefficient=0.4:Specular=1:SpecularColor=#ffffff:FrontSpecularCoefficient=0.3:BackSpecularCoefficient=0.3:Emissive=0:EmissiveColor=#000000:FrontEmissiveCoefficient=0:BackEmissiveCoefficient=0")
ipar.setParameter(objId, "OCCViewer_0_EdgeWidth", "1")
ipar.setParameter(objId, "OCCViewer_0_IsosWidth", "1")
# set up entry GEOM_7 (Triangle_tgt) parameters
objId = geompy.getObjectID(Triangle_tgt)
ipar.setParameter(objId, "OCCViewer_0_Visibility", "On")
ipar.setParameter(objId, "OCCViewer_0_DisplayMode", "2")
ipar.setParameter(objId, "OCCViewer_0_Color", "0:1:0.498039")
ipar.setParameter(objId, "OCCViewer_0_Transparency", "0")
ipar.setParameter(objId, "OCCViewer_0_TopLevelFlag", "false")
ipar.setParameter(objId, "OCCViewer_0_Isos", "0:0")
ipar.setParameter(objId, "OCCViewer_0_VectorMode", "false")
ipar.setParameter(objId, "OCCViewer_0_VerticesMode", "false")
ipar.setParameter(objId, "OCCViewer_0_NameMode", "false")
ipar.setParameter(objId, "OCCViewer_0_DeflectionCoeff", "0.001")
ipar.setParameter(objId, "OCCViewer_0_MarkerType", "7:3")
ipar.setParameter(objId, "OCCViewer_0_Material", "Physical=0:FrontShininess=0.13:BackShininess=0.13:Transparency=0:Ambient=1:AmbientColor=#333333:FrontAmbientCoefficient=0.3:BackAmbientCoefficient=0.25:Diffuse=1:DiffuseColor=#000000:FrontDiffuseCoefficient=0.5:BackDiffuseCoefficient=0.4:Specular=1:SpecularColor=#ffffff:FrontSpecularCoefficient=0.3:BackSpecularCoefficient=0.3:Emissive=0:EmissiveColor=#000000:FrontEmissiveCoefficient=0:BackEmissiveCoefficient=0")
ipar.setParameter(objId, "OCCViewer_0_EdgeWidth", "1")
ipar.setParameter(objId, "OCCViewer_0_IsosWidth", "1")
# set up entry GEOM_9 (Bnd_box) parameters
objId = geompy.getObjectID(Bnd_box)
ipar.setParameter(objId, "OCCViewer_0_Visibility", "On")
ipar.setParameter(objId, "OCCViewer_0_DisplayMode", "0")
ipar.setParameter(objId, "OCCViewer_0_Color", "1:0:0")
ipar.setParameter(objId, "OCCViewer_0_Transparency", "0")
ipar.setParameter(objId, "OCCViewer_0_TopLevelFlag", "false")
ipar.setParameter(objId, "OCCViewer_0_Isos", "0:0")
ipar.setParameter(objId, "OCCViewer_0_VectorMode", "false")
ipar.setParameter(objId, "OCCViewer_0_VerticesMode", "false")
ipar.setParameter(objId, "OCCViewer_0_NameMode", "false")
ipar.setParameter(objId, "OCCViewer_0_DeflectionCoeff", "0.001")
ipar.setParameter(objId, "OCCViewer_0_MarkerType", "7:3")
ipar.setParameter(objId, "OCCViewer_0_Material", "Physical=0:FrontShininess=0.13:BackShininess=0.13:Transparency=0:Ambient=1:AmbientColor=#333333:FrontAmbientCoefficient=0.3:BackAmbientCoefficient=0.25:Diffuse=1:DiffuseColor=#000000:FrontDiffuseCoefficient=0.5:BackDiffuseCoefficient=0.4:Specular=1:SpecularColor=#ffffff:FrontSpecularCoefficient=0.3:BackSpecularCoefficient=0.3:Emissive=0:EmissiveColor=#000000:FrontEmissiveCoefficient=0:BackEmissiveCoefficient=0")
ipar.setParameter(objId, "OCCViewer_0_EdgeWidth", "1")
ipar.setParameter(objId, "OCCViewer_0_IsosWidth", "1")


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
  iparameters.getSession().restoreVisualState(1)
