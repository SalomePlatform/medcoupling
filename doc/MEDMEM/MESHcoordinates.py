#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2011  CEA/DEN, EDF R&D, OPEN CASCADE
#
# Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
# CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
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
