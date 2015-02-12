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

# _T1A
import collections
def tree():
    return collections.defaultdict(tree)

fieldTree = tree()
meshDict = {}
# _T1B

import os
filename = "timeseries.med"
filepath = os.path.join(os.path.abspath(os.path.dirname(__file__)),filename)

# _T2A
from MEDLoader import MEDLoader
meshNames = MEDLoader.GetMeshNames(filepath)

meshDimRelToMax = 0 # 0 = no restriction

for meshName in meshNames:
    mesh = MEDLoader.ReadUMeshFromFile(filepath,meshName,meshDimRelToMax)
    meshDict[meshName] = mesh

    fieldNames = MEDLoader.GetAllFieldNamesOnMesh(filepath,meshName)
    for fieldName in fieldNames:
        listOfTypes = MEDLoader.GetTypesOfField(filepath,meshName,fieldName)
        for typeOfDiscretization in listOfTypes:
            fieldIterations = MEDLoader.GetFieldIterations(typeOfDiscretization,
                                                           filepath,
                                                           meshName,
                                                           fieldName)
            for fieldIteration in fieldIterations:
                itNumber = fieldIteration[0]
                itOrder  = fieldIteration[1]

                field = MEDLoader.ReadField(typeOfDiscretization,
                                            filepath,
                                            meshName,
                                            meshDimRelToMax,
                                            fieldName,
                                            itNumber,
                                            itOrder)

                fieldTree\
                           [meshName]\
                           [fieldName]\
                           [typeOfDiscretization]\
                           [itNumber][itOrder] = field
# _T2B

# Q: use a list of structures whose an attribute could be a
# MEDCoupling field? Or a tree that you cross using attribute and
# whose leaves are the MEDCoupling fields?
# R: I think that the default structure should be a simple list that
# store objects whith properties that corresponds to the metadata (and
# if loaded the MEDCouplingField or Mesh). Then for specific request,
# a BTree could be create to organize the search (for example if we
# request all the fields for a given iteration step, then we should
# use the iteration step as a first classifaction switch of the tree

print fieldTree.keys()

# _T3A
for meshName in fieldTree.keys():
    print "%s"%meshName
    for fieldName in fieldTree[meshName].keys():
        print "  %s"%fieldName
        for fieldType in fieldTree[meshName][fieldName].keys():
            print "    %s"%fieldType
            for itNumber in fieldTree[meshName][fieldName][fieldType].keys():
                for itOrder in fieldTree[meshName][fieldName][fieldType][itNumber].keys():
                    print "      (%s,%s)"%(itNumber,itOrder)
                    print fieldTree[meshName][fieldName][fieldType][itNumber][itOrder]
# _T3B
