#!/usr/bin/env python
# Copyright (C) 2012-2015  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
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

# This script illustrates the basic features of MEDLoader
# (gboulant, 17 nov 2012)
import os
filename = "timeseries.med"
filepath = os.path.join(os.path.abspath(os.path.dirname(__file__)),filename)

# _T1A
from MEDLoader import MEDLoader
meshNames = MEDLoader.GetMeshNames(filepath)
# _T1B
meshName=meshNames[0]
# _T2A
fieldNames = MEDLoader.GetAllFieldNamesOnMesh(filepath,meshName)
# _T2B
fieldName=fieldNames[0]
# _T3A
listOfTypes = MEDLoader.GetTypesOfField(filepath,meshName,fieldName)
# _T3B
typeOfDiscretization=listOfTypes[0]
# _T4A
fieldIterations = MEDLoader.GetFieldIterations(typeOfDiscretization,
                                               filepath,
                                               meshName,
                                               fieldName)
# _T4B

iteration = fieldIterations[0]
iterationNumber = iteration[0]
iterationOrder  = iteration[1]

dimrestriction = 0
# _T5A
mesh = MEDLoader.ReadUMeshFromFile(filepath, meshName, dimrestriction)
# _T5B
# _T6A
field = MEDLoader.ReadField(typeOfDiscretization,
                            filepath, meshName, dimrestriction,
                            fieldName, iterationNumber, iterationOrder)
# _T6B
