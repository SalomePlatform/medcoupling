# -*- coding: utf-8 -*-

###
### This script is intended to be launched in a new SALOME study
###

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)
sys.path.insert( 0, r'/misc/dn27/users_Linux/eap/salome/tmp')

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
ipar.append("VTKViewer_1", "VTK scene:1 - viewer:1")
ipar.append("VTKViewer_1", """<?xml version="1.0"?>
<ViewState>
    <Position X="6.64164" Y="2.49199" Z="130.355"/>
    <FocalPoint X="6.64164" Y="2.49199" Z="0"/>
    <ViewUp X="0" Y="1" Z="0"/>
    <ViewScale Parallel="28.0904" X="1" Y="1" Z="1"/>
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

from MEDLoader import ReadMeshFromFile, WriteMesh

coordsArr=DataArrayDouble(range(8))
m=MEDCouplingCMesh("mesh")
m.setCoords(coordsArr,coordsArr)
m = m.buildUnstructured()
#m.rotate([3.5,3.5],math.pi/4.)

#WriteMesh("mesh1.med",m,True)

from MEDCoupling import MEDCouplingSkyLineArray
import MEDPartitioner
a,b=m.computeNeighborsOfCells()
sk=MEDCouplingSkyLineArray(b,a)
g=MEDPartitioner.MEDPartitioner.Graph(sk)
g.partGraph(4)
procIdOnCells=g.getPartition().getValuesArray()
p0=procIdOnCells.findIdsEqual(0)
part0=m[p0]

WriteMesh("part.med",part0,True)

boundary_nodes_part0=part0.findBoundaryNodes()
boundary_cells_part0=p0[part0.getCellIdsLyingOnNodes(boundary_nodes_part0,False)]
# starting from knowledge of neighborhood it s possible to know the neighbors of boundary_cells_part0
neighbors_boundary_cells_part0=MEDCouplingUMesh.ExtractFromIndexedArrays(boundary_cells_part0,a,b)[0]
neighbors_boundary_cells_part0.sort()
neighbors_boundary_cells_part0=neighbors_boundary_cells_part0.buildUnique()
#
layer_of_part0=neighbors_boundary_cells_part0.buildSubstraction(p0)
#
whole_part_with_layer=DataArrayInt.Aggregate([p0,layer_of_part0])
whole_part_with_layer.sort()
part0_with_layer=m[whole_part_with_layer]
part0_with_layer.setName("part0_with_layer")

WriteMesh("part.med",part0_with_layer,False)

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)
#([Mesh_1], status) = smesh.CreateMeshesFromMED(r'mesh1.med')
([Mesh_2,Mesh_2_WL], status) = smesh.CreateMeshesFromMED(r'part.med')
#Mesh_2.TranslateObject( Mesh_2,       [ 9, 0, 0 ], 0 )
Mesh_2_WL.TranslateObject( Mesh_2_WL, [ 5, 0, 0 ], 0 )

### Store presentation parameters of displayed objects
import iparameters
ipar = iparameters.IParameters(theStudy.GetModuleParameters("Interface Applicative", "SMESH", 1))

#Set up entries:
# set up entry SMESH_3 (mesh) parameters
ipar.setParameter("SMESH_3", "VTKViewer_0_Visibility", "On")
ipar.setParameter("SMESH_3", "VTKViewer_0_Representation", "2")
ipar.setParameter("SMESH_3", "VTKViewer_0_IsShrunk", "0")
ipar.setParameter("SMESH_3", "VTKViewer_0_Entities", "e:0:f:1:v:0:0d:0:b:0")
ipar.setParameter("SMESH_3", "VTKViewer_0_Colors", "surface:0:0.666667:1:backsurface:100:volume:1:0:0.666667:-100:edge:0:0.666667:1:node:1:0:0:outline:0:0.27451:0:elem0d:0:1:0:ball:0:0.333333:1:orientation:1:1:1")
ipar.setParameter("SMESH_3", "VTKViewer_0_Sizes", "line:1:outline:1:elem0d:5:ball:10:1:shrink:0.75:orientation:0.1:0")
ipar.setParameter("SMESH_3", "VTKViewer_0_PointMarker", "std:1:9")
ipar.setParameter("SMESH_3", "VTKViewer_0_Opacity", "1")
ipar.setParameter("SMESH_3", "VTKViewer_0_ClippingPlane", "Off")
# set up entry SMESH_4 (mesh) parameters
ipar.setParameter("SMESH_4", "VTKViewer_0_Visibility", "On")
ipar.setParameter("SMESH_4", "VTKViewer_0_Representation", "2")
ipar.setParameter("SMESH_4", "VTKViewer_0_IsShrunk", "0")
ipar.setParameter("SMESH_4", "VTKViewer_0_Entities", "e:0:f:1:v:0:0d:0:b:0")
ipar.setParameter("SMESH_4", "VTKViewer_0_Colors", "surface:0:0.666667:1:backsurface:100:volume:1:0:0.666667:-100:edge:0:0.666667:1:node:1:0:0:outline:0:0.27451:0:elem0d:0:1:0:ball:0:0.333333:1:orientation:1:1:1")
ipar.setParameter("SMESH_4", "VTKViewer_0_Sizes", "line:1:outline:1:elem0d:5:ball:10:1:shrink:0.75:orientation:0.1:0")
ipar.setParameter("SMESH_4", "VTKViewer_0_PointMarker", "std:1:9")
ipar.setParameter("SMESH_4", "VTKViewer_0_Opacity", "1")
ipar.setParameter("SMESH_4", "VTKViewer_0_ClippingPlane", "Off")


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
  iparameters.getSession().restoreVisualState(1)

import libSALOME_Swig
gui = libSALOME_Swig.SALOMEGUI_Swig()
gui.AddIObject( salome.ObjectToID( Mesh_2.GetMesh() ))
gui.AddIObject( salome.ObjectToID( Mesh_2_WL.GetMesh() ))
