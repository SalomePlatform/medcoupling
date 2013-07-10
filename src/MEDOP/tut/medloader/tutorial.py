#!/usr/bin/env python
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
