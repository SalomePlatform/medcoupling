# -*- coding: utf-8 -*-

###
### This script is intended to be launched in a new SALOME study
###

import os
import salome

salome.salome_init()
theStudy = salome.myStudy

import iparameters
ipar = iparameters.IParameters(salome.myStudy.GetCommonParameters("Interface Applicative", 1), True)

#Set up visual properties:
ipar.setProperty("AP_ACTIVE_VIEW", "VTKViewer_0_0")
ipar.setProperty("AP_WORKSTACK_INFO", "0000000100000000000000020100000001000003a0000000040000000100000001000000080000001a00560054004b005600690065007700650072005f0030005f00300000000102")
ipar.setProperty("AP_ACTIVE_MODULE", "Mesh")
ipar.setProperty("AP_SAVEPOINT_NAME", "GUI state: 2")
#Set up lists:
# fill list AP_VIEWERS_LIST
ipar.append("AP_VIEWERS_LIST", "VTKViewer_1")
# fill list VTKViewer_1
ipar.append("VTKViewer_1", "VTK scene:2 - viewer:1")
ipar.append("VTKViewer_1", """<?xml version="1.0"?>
<ViewState>
    <Position X="1.81423" Y="1.44665" Z="20.0926"/>
    <FocalPoint X="1.81423" Y="1.33011" Z="-0.0048134"/>
    <ViewUp X="0" Y="0.999983" Z="-0.0057984"/>
    <ViewScale Parallel="6.01999" X="1" Y="1" Z="1"/>
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

import sys
if sys.platform == "win32":
    from MEDCouplingCompat import *
else:
    from MEDCoupling import *
from MEDLoader import WriteMesh

arr=DataArrayDouble(range(2))
mc=MEDCouplingCMesh("m1")
mc.setCoords(arr,arr)
m1=mc.buildUnstructured()
m2=m1.deepCopy()
m2.translate([1.5,0])
m1.convertLinearCellsToQuadratic(0)
m2.convertLinearCellsToQuadratic(1)
m2.setName("m2")
WriteMesh("mesh1.med",m1,True)
WriteMesh("mesh1.med",m2,False)

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(salome.myStudy)
([mesh_1,mesh_2], status) = smesh.CreateMeshesFromMED(r'mesh1.med')
g1 = mesh_1.MakeGroupByIds("all nodes", SMESH.NODE, range(1,10))
g2 = mesh_2.MakeGroupByIds("all nodes", SMESH.NODE, range(1,10))

### Store presentation parameters of displayed objects
import iparameters
ipar = iparameters.IParameters(theStudy.GetModuleParameters("Interface Applicative", "SMESH", 1))

#Set up entries:
# set up entry SMESH_3 (m1) parameters
ipar.setParameter("SMESH_3", "VTKViewer_0_Visibility", "On")
ipar.setParameter("SMESH_3", "VTKViewer_0_Representation", "1")
ipar.setParameter("SMESH_3", "VTKViewer_0_IsShrunk", "0")
ipar.setParameter("SMESH_3", "VTKViewer_0_Entities", "e:0:f:1:v:0:0d:0:b:0")
ipar.setParameter("SMESH_3", "VTKViewer_0_Colors", "surface:0:0.666667:1:backsurface:100:volume:1:0:0.666667:-100:edge:0:0.666667:1:node:1:0:0:outline:0:0.27451:0:elem0d:0:1:0:ball:0:0.333333:1:orientation:1:1:1")
ipar.setParameter("SMESH_3", "VTKViewer_0_Sizes", "line:1:outline:1:elem0d:5:ball:10:1:shrink:0.75:orientation:0.1:0")
ipar.setParameter("SMESH_3", "VTKViewer_0_PointMarker", "std:1:9")
ipar.setParameter("SMESH_3", "VTKViewer_0_Opacity", "1")
ipar.setParameter("SMESH_3", "VTKViewer_0_ClippingPlane", "Off")
# set up entry SMESH_3:11:1 (all nodes) parameters
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_Visibility", "On")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_Representation", "0")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_IsShrunk", "0")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_Entities", "e:1:f:1:v:1:0d:1:b:1")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_Colors", "surface:0:0.666667:1:backsurface:100:volume:1:0:0.666667:-100:edge:0:0.666667:1:node:1:0.666667:0:outline:0:0.27451:0:elem0d:0:1:0:ball:0:0.333333:1:orientation:1:1:1")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_Sizes", "line:1:outline:1:elem0d:5:ball:10:1:shrink:0.75:orientation:0.1:0")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_PointMarker", "std:1:9")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_Opacity", "1")
ipar.setParameter("SMESH_3:11:1", "VTKViewer_0_ClippingPlane", "Off")
# set up entry SMESH_4 (m2) parameters
ipar.setParameter("SMESH_4", "VTKViewer_0_Visibility", "On")
ipar.setParameter("SMESH_4", "VTKViewer_0_Representation", "1")
ipar.setParameter("SMESH_4", "VTKViewer_0_IsShrunk", "0")
ipar.setParameter("SMESH_4", "VTKViewer_0_Entities", "e:1:f:1:v:1:0d:1:b:1")
ipar.setParameter("SMESH_4", "VTKViewer_0_Colors", "surface:0:0.666667:1:backsurface:100:volume:1:0:0.666667:-100:edge:0:0.666667:1:node:1:0:0:outline:0:0.27451:0:elem0d:0:1:0:ball:0:0.333333:1:orientation:1:1:1")
ipar.setParameter("SMESH_4", "VTKViewer_0_Sizes", "line:1:outline:1:elem0d:5:ball:10:1:shrink:0.75:orientation:0.1:0")
ipar.setParameter("SMESH_4", "VTKViewer_0_PointMarker", "std:1:9")
ipar.setParameter("SMESH_4", "VTKViewer_0_Opacity", "1")
ipar.setParameter("SMESH_4", "VTKViewer_0_ClippingPlane", "Off")
# set up entry SMESH_4:11:1 (all nodes) parameters
ipar.setParameter("SMESH_4:11:1", "VTKViewer_0_Visibility", "On")
ipar.setParameter("SMESH_4:11:1", "VTKViewer_0_Representation", "0")
ipar.setParameter("SMESH_4:11:1", "VTKViewer_0_IsShrunk", "0")
ipar.setParameter("SMESH_4:11:1", "VTKViewer_0_Entities", "e:1:f:1:v:1:0d:1:b:1")
ipar.setParameter("SMESH_4:11:1", "VTKViewer_0_Colors", "surface:0:0.666667:1:backsurface:100:volume:1:0:0.666667:-100:edge:0:0.666667:1:node:1:0.666667:0:outline:0:0.27451:0:elem0d:0:1:0:ball:0:0.333333:1:orientation:1:1:1")
ipar.setParameter("SMESH_4:11:1", "VTKViewer_0_Sizes", "line:1:outline:1:elem0d:5:ball:10:1:shrink:0.75:orientation:0.1:0")
ipar.setParameter("SMESH_4:11:1", "VTKViewer_0_PointMarker", "std:1:9")
ipar.setParameter("SMESH_4:11:1", "VTKViewer_0_Opacity", "1")
ipar.setParameter("SMESH_4:11:1", "VTKViewer_0_ClippingPlane", "Off")


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
  iparameters.getSession().restoreVisualState(1)
