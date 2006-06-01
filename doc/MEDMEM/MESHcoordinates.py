# Copyright (C) 2005  OPEN CASCADE, CEA, EDF R&D, LEG
#           PRINCIPIA R&D, EADS CCR, Lip6, BV, CEDRAT
# 
from libMEDMEM_Swig import *

MedFile = "pointe.med"
meshName = "maa1"

myMesh = MESH(MED_DRIVER,MedFile,meshName)

name = myMesh.getName()

print "Mesh name : ",name
spaceDimension = myMesh.getSpaceDimension()
numberOfNodes = myMesh.getNumberOfNodes()
print "Space Dimension : ",spaceDimension
print "Number of Nodes : ",numberOfNodes

print "Show Nodes Coordinates :"
print "Name :"
coordinatesNames = myMesh.getCoordinatesNames()
for i in range(spaceDimension):
    coordinateName = coordinatesNames[i]
    print " - ",coordinateName

print "Unit :"
coordinatesUnits = myMesh.getCoordinatesUnits()
for i in range(spaceDimension):
    coordinateUnit = coordinatesUnits[i]
    print " - ",coordinateUnit

coordinates = myMesh.getCoordinates(MED_FULL_INTERLACE)
for i in range(numberOfNodes):
    print "Node ",(i+1)," : ",coordinates[i*spaceDimension:(i+1)*spaceDimension]
