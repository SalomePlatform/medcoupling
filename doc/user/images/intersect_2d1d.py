# -*- coding: utf-8 -*-

###
### This script is intended to be launched in a new SALOME study
###

import sys
import salome
import math

salome.salome_init()
theStudy = salome.myStudy

import iparameters
ipar = iparameters.IParameters(salome.myStudy.GetCommonParameters("Interface Applicative", 1), True)

#Set up visual properties:
ipar.setProperty("AP_ACTIVE_VIEW", "VTKViewer_0_0")
ipar.setProperty("AP_WORKSTACK_INFO", "0000000100000000000000020100000001000003a0000000040000000100000002000000080000001a00560054004b005600690065007700650072005f0030005f00300000000202")
ipar.setProperty("AP_ACTIVE_MODULE", "Mesh")
ipar.setProperty("AP_SAVEPOINT_NAME", "GUI state: 2")
#Set up lists:
# fill list AP_VIEWERS_LIST
ipar.append("AP_VIEWERS_LIST", "VTKViewer_1")
# fill list VTKViewer_1
ipar.append("VTKViewer_1", "VTK scene:1 - viewer:1")
ipar.append("VTKViewer_1", """<?xml version="1.0"?>
<ViewState>
    <Position X="4.41007" Y="1.25267" Z="-97.8303"/>
    <FocalPoint X="2.00346" Y="1.42923" Z="0.314897"/>
    <ViewUp X="0.000135221" Y="-0.999998" Z="0.00180225"/>
    <ViewScale Parallel="5.42135" X="1" Y="1" Z="1"/>
    <DisplayCubeAxis Show="0"/>
    <GraduatedAxis Axis="X">
        <Title isVisible="1" Text="X" Font="0" Bold="0" Italic="0" Shadow="0">
            <Color R="1" G="0" B="0"/>
        </Title>
        <Labels isVisible="1" Number="3" Offset="2" Font="0" Bold="0" Italic="0" Shadow="0">
            <Color R="1" G="0" B="0"/>
        </Labels>
        <TickMarks isVisible="1" Length="5"/>
    </GraduatedAxis>
    <GraduatedAxis Axis="Y">
        <Title isVisible="1" Text="Y" Font="0" Bold="0" Italic="0" Shadow="0">
            <Color R="0" G="1" B="0"/>
        </Title>
        <Labels isVisible="1" Number="3" Offset="2" Font="0" Bold="0" Italic="0" Shadow="0">
            <Color R="0" G="1" B="0"/>
        </Labels>
        <TickMarks isVisible="1" Length="5"/>
    </GraduatedAxis>
    <GraduatedAxis Axis="Z">
        <Title isVisible="1" Text="Z" Font="0" Bold="0" Italic="0" Shadow="0">
            <Color R="0" G="0" B="1"/>
        </Title>
        <Labels isVisible="1" Number="3" Offset="2" Font="0" Bold="0" Italic="0" Shadow="0">
            <Color R="0" G="0" B="1"/>
        </Labels>
        <TickMarks isVisible="1" Length="5"/>
    </GraduatedAxis>
    <Trihedron isShown="0" Size="100"/>
    <Background Value="bt=1;fn=;tm=0;ts=false;c1=#ffffff;c2=#000000;gt=-1;gr="/>
</ViewState>
""")
# fill list AP_MODULES_LIST
ipar.append("AP_MODULES_LIST", "Mesh")

if sys.platform == "win32":
    from MEDCouplingCompat import *
else:
    from MEDCoupling import *

coordsArr=DataArrayDouble(range(6))
mesh2d=MEDCouplingCMesh("mesh2d")
mesh2d.setCoords(coordsArr,coordsArr[:2])
mesh2d=mesh2d.buildUnstructured()

mesh1d=MEDCouplingCMesh("mesh1d")
mesh1d.setCoords(coordsArr,coordsArr[:1])
mesh1d=mesh1d.buildUnstructured()
mesh1d.rotate( [2.3,0], math.radians( 25 ))
mesh1d.translate( [0.2,0.4] )

m2d,m1d,a2d,a1d=MEDCouplingUMesh.Intersect2DMeshWith1DLine( mesh2d, mesh1d, 1e-12 )
m1d.setName("m1d")

from MEDLoader import WriteMesh
WriteMesh("mesh1.med",mesh2d,True)
WriteMesh("mesh1.med",mesh1d,False)
WriteMesh("mesh1.med",m2d,False)
WriteMesh("mesh1.med",m1d,False)

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)
([m_1d, m_2d, mesh_1d, mesh_2d], status) = smesh.CreateMeshesFromMED(r'mesh1.med')

m_1d.TranslateObject( m_1d, [0,2,0], False )
m_2d.TranslateObject( m_2d, [0,2,0], False )

m_1d.MakeGroup   ( "nodes", SMESH.NODE, SMESH.FT_NodeConnectivityNumber,">",0)
mesh_1d.MakeGroup( "all nodes", SMESH.NODE, SMESH.FT_RangeOfIds,"=","1-20")


### Store presentation parameters of displayed objects
import iparameters
ipar = iparameters.IParameters(theStudy.GetModuleParameters("Interface Applicative", "SMESH", 1))

#Set up entries:
# set up entry SMESH_4 (merge) parameters
ipar.setParameter("SMESH_4", "VTKViewer_0_Visibility", "On")
ipar.setParameter("SMESH_4", "VTKViewer_0_Representation", "2")
ipar.setParameter("SMESH_4", "VTKViewer_0_IsShrunk", "1")
ipar.setParameter("SMESH_4", "VTKViewer_0_Entities", "e:0:f:1:v:0:0d:0:b:0")
ipar.setParameter("SMESH_4", "VTKViewer_0_Colors", "surface:0:0.666667:1:backsurface:100:volume:1:0:0.666667:-100:edge:0:0.666667:1:node:1:0:0:outline:0:0.27451:0:elem0d:0:1:0:ball:0:0.333333:1:orientation:1:1:1")
ipar.setParameter("SMESH_4", "VTKViewer_0_Sizes", "line:1:outline:1:elem0d:5:ball:10:1:shrink:0.75:orientation:0.1:0")
ipar.setParameter("SMESH_4", "VTKViewer_0_PointMarker", "std:1:9")
ipar.setParameter("SMESH_4", "VTKViewer_0_Opacity", "1")
ipar.setParameter("SMESH_4", "VTKViewer_0_ClippingPlane", "Off")
# set up entry SMESH_3 (m1d) parameters
ipar.setParameter("SMESH_3", "VTKViewer_0_Visibility", "On")
ipar.setParameter("SMESH_3", "VTKViewer_0_Representation", "1")
ipar.setParameter("SMESH_3", "VTKViewer_0_IsShrunk", "0")
ipar.setParameter("SMESH_3", "VTKViewer_0_Entities", "e:1:f:0:v:0:0d:0:b:0")
ipar.setParameter("SMESH_3", "VTKViewer_0_Colors", "surface:0:0.666667:1:backsurface:100:volume:1:0:0.666667:-100:edge:0:0.666667:1:node:1:0:0:outline:0:0.27451:0:elem0d:0:1:0:ball:0:0.333333:1:orientation:1:1:1")
ipar.setParameter("SMESH_3", "VTKViewer_0_Sizes", "line:1:outline:1:elem0d:5:ball:10:1:shrink:0.75:orientation:0.1:0")
ipar.setParameter("SMESH_3", "VTKViewer_0_PointMarker", "std:1:9")
ipar.setParameter("SMESH_3", "VTKViewer_0_Opacity", "1")
ipar.setParameter("SMESH_3", "VTKViewer_0_ClippingPlane", "Off")
# set up entry SMESH_5 (mesh1d) parameters
ipar.setParameter("SMESH_5", "VTKViewer_0_Visibility", "On")
ipar.setParameter("SMESH_5", "VTKViewer_0_Representation", "1")
ipar.setParameter("SMESH_5", "VTKViewer_0_IsShrunk", "0")
ipar.setParameter("SMESH_5", "VTKViewer_0_Entities", "e:1:f:0:v:0:0d:0:b:0")
ipar.setParameter("SMESH_5", "VTKViewer_0_Colors", "surface:0:0.666667:1:backsurface:100:volume:1:0:0.666667:-100:edge:0:0.666667:1:node:1:0:0:outline:0:0.27451:0:elem0d:0:1:0:ball:0:0.333333:1:orientation:1:1:1")
ipar.setParameter("SMESH_5", "VTKViewer_0_Sizes", "line:1:outline:1:elem0d:5:ball:10:1:shrink:0.75:orientation:0.1:0")
ipar.setParameter("SMESH_5", "VTKViewer_0_PointMarker", "std:1:9")
ipar.setParameter("SMESH_5", "VTKViewer_0_Opacity", "1")
ipar.setParameter("SMESH_5", "VTKViewer_0_ClippingPlane", "Off")
# set up entry SMESH_6 (mesh2d) parameters
ipar.setParameter("SMESH_6", "VTKViewer_0_Visibility", "On")
ipar.setParameter("SMESH_6", "VTKViewer_0_Representation", "2")
ipar.setParameter("SMESH_6", "VTKViewer_0_IsShrunk", "1")
ipar.setParameter("SMESH_6", "VTKViewer_0_Entities", "e:0:f:1:v:0:0d:0:b:0")
ipar.setParameter("SMESH_6", "VTKViewer_0_Colors", "surface:0:0.666667:1:backsurface:100:volume:1:0:0.666667:-100:edge:0:0.666667:1:node:1:0:0:outline:0:0.27451:0:elem0d:0:1:0:ball:0:0.333333:1:orientation:1:1:1")
ipar.setParameter("SMESH_6", "VTKViewer_0_Sizes", "line:1:outline:1:elem0d:5:ball:10:1:shrink:0.75:orientation:0.1:0")
ipar.setParameter("SMESH_6", "VTKViewer_0_PointMarker", "std:1:9")
ipar.setParameter("SMESH_6", "VTKViewer_0_Opacity", "1")
ipar.setParameter("SMESH_6", "VTKViewer_0_ClippingPlane", "Off")
# set up entry SMESH_3:11:1 (nodes) parameters
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_Visibility", "On")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_Representation", "0")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_IsShrunk", "0")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_Entities", "e:1:f:1:v:1:0d:1:b:1")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_Colors", "surface:0:0.666667:1:backsurface:100:volume:1:0:0.666667:-100:edge:0:0.666667:1:node:1:0.666667:0:outline:0:0.27451:0:elem0d:0:1:0:ball:0:0.333333:1:orientation:1:1:1")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_Sizes", "line:1:outline:1:elem0d:5:ball:10:1:shrink:0.75:orientation:0.1:0")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_PointMarker", "std:1:9")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_Opacity", "1")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_ClippingPlane", "Off")
# set up entry SMESH_5:11:1 (all nodes) parameters
ipar.setParameter("SMESH_5:11:1", "VTKViewer_0_Visibility", "On")
ipar.setParameter("SMESH_5:11:1", "VTKViewer_0_Representation", "0")
ipar.setParameter("SMESH_5:11:1", "VTKViewer_0_IsShrunk", "0")
ipar.setParameter("SMESH_5:11:1", "VTKViewer_0_Entities", "e:1:f:1:v:1:0d:1:b:1")
ipar.setParameter("SMESH_5:11:1", "VTKViewer_0_Colors", "surface:0:0.666667:1:backsurface:100:volume:1:0:0.666667:-100:edge:0:0.666667:1:node:1:0.666667:0:outline:0:0.27451:0:elem0d:0:1:0:ball:0:0.333333:1:orientation:1:1:1")
ipar.setParameter("SMESH_5:11:1", "VTKViewer_0_Sizes", "line:1:outline:1:elem0d:5:ball:10:1:shrink:0.75:orientation:0.1:0")
ipar.setParameter("SMESH_5:11:1", "VTKViewer_0_PointMarker", "std:1:9")
ipar.setParameter("SMESH_5:11:1", "VTKViewer_0_Opacity", "1")
ipar.setParameter("SMESH_5:11:1", "VTKViewer_0_ClippingPlane", "Off")


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
  iparameters.getSession().restoreVisualState(1)
