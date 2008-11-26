#  Copyright (C) 2007-2008  CEA/DEN, EDF R&D, OPEN CASCADE
#
#  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
#  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#
from libMEDMEM_Swig import *

MedFile = "pointe.med"
#MedFile = "carre_quad4_3.med"
#MedFile = "polyedres.med"
#MedFile = "polygones.med"
meshName = "maa1"
#meshName = "CARRE_EN_QUAD4"
#meshName = "Erreur orientation"
#meshName = "Bord"

myMesh = MESH(MED_DRIVER,MedFile,meshName)
myMesh.read()

nameMesh = myMesh.getName()

print "Mesh name : ",nameMesh

numberOfTypes = myMesh.getNumberOfTypes(MED_CELL)
print "Show Connectivity (Nodal) : "

# This example use access with a specified medGeometryElement through
# CELLMODEL class

for i in range(numberOfTypes):
    cellType = myMesh.getCellType(MED_CELL,i)
    nameType = cellType.getName()
    type = cellType.getType()
    numberOfElements = myMesh.getNumberOfElements(MED_CELL,type)
    numberOfNodesPerCell = cellType.getNumberOfNodes()
    connectivity = myMesh.getConnectivity(MED_FULL_INTERLACE,
                                          MED_NODAL,MED_CELL,type)
    print "For Type ",nameType," : "
    for j in range(numberOfElements):
        print "Element ",(j+1)," : ",connectivity[j*numberOfNodesPerCell:
                                                  (j+1)*numberOfNodesPerCell]

print "Show Reverse Nodal Connectivity :"

# This example use global access with index array

numberOfNodes = myMesh.getNumberOfNodes()

reverseNodalConnectivity = myMesh.getReverseConnectivity(MED_NODAL)
reverseNodalConnectivityIndex = myMesh.getReverseConnectivityIndex(MED_NODAL)

for i in range(numberOfNodes):
    indexBegin = reverseNodalConnectivityIndex[i]
    indexEnd = reverseNodalConnectivityIndex[i+1]

    # Index value begin at 1 so (index-1) is in fact used here

    print "Node ",(i+1)," : ",reverseNodalConnectivity[(indexBegin-1):
                                                       (indexEnd-1)]

print "Show Connectivity (Descending) :"

# This example use global access with index array

numberOfElements = myMesh.getNumberOfElements(MED_CELL,MED_ALL_ELEMENTS)
descendingConnectivity = myMesh.getConnectivity(MED_FULL_INTERLACE,
                                                MED_DESCENDING,MED_CELL,
                                                MED_ALL_ELEMENTS)
descendingConnectivityIndex = myMesh.getConnectivityIndex(MED_DESCENDING,
                                                          MED_CELL)

for i in range(numberOfElements):
    indexBegin = descendingConnectivityIndex[i]
    indexEnd = descendingConnectivityIndex[i+1]

    # Index value begin at 1 so (index-1) is in fact used here

    print "Element ",(i+1)," : ",descendingConnectivity[(indexBegin-1):
                                                        (indexEnd-1)]

print "Show Reverse Descending Connectivity :"

# This example use global access with index array

meshDimension = myMesh.getMeshDimension()

if (meshDimension == 1):
    print "ERROR : Mesh Dimension = 1"
    print "Then the Reverse Descending Connectivity could not be seen"
else:
    if (meshDimension == 2):
        constituent = "Edge"
        constituentEntity = MED_EDGE

    if (meshDimension == 3):
        constituent = "Face"
        constituentEntity = MED_FACE

    numberOfConstituents = myMesh.getNumberOfElements(constituentEntity,
                                                      MED_ALL_ELEMENTS)
    reverseDescendingConnectivity = myMesh.getReverseConnectivity(
        MED_DESCENDING)
    reverseDescendingConnectivityIndex = myMesh.getReverseConnectivityIndex(
        MED_DESCENDING)

    for i in range(numberOfConstituents):
        indexBegin = reverseDescendingConnectivityIndex[i]
        indexEnd = reverseDescendingConnectivityIndex[i+1]

        # Index value begin at 1 so (index-1) is in fact used here

        print constituent," : ",(i+1)," : ",reverseDescendingConnectivity[
            (indexBegin-1):(indexEnd-1)]

    print "Show ",constituent," Connectivity (Nodal) :"

    constituentConnectivity = myMesh.getConnectivity(MED_FULL_INTERLACE,
                                                     MED_NODAL,
                                                     constituentEntity,
                                                     MED_ALL_ELEMENTS)
    constituentConnectivityIndex = myMesh.getConnectivityIndex(MED_NODAL,
                                                               constituentEntity)

    for i in range(numberOfConstituents):
        indexBegin = constituentConnectivityIndex[i]
        indexEnd = constituentConnectivityIndex[i+1]

        # Index value begin at 1 so (index-1) is in fact used here

        print constituent," : ",(i+1)," : ",constituentConnectivity[
            (indexBegin-1):(indexEnd-1)]
        pass
    pass

nbPolygons = myMesh.getNumberOfPolygons()
if nbPolygons > 0 :
    print ""
    print "     Show Connectivity (Nodal) of POLYGONS:"
    print ""
    connectivity = myMesh.getPolygonsConnectivity(MED_NODAL,MED_CELL)
    index = myMesh.getPolygonsConnectivityIndex(MED_NODAL,MED_CELL)
    for j in range(nbPolygons):
        print "       Polygon",(j+1)," ",connectivity[ index[j]-1 : index[j+1]-1 ]
        pass
    pass

nbPolyhedrons = myMesh.getNumberOfPolyhedron()
if nbPolyhedrons > 0 :
    print ""
    print "     Show Connectivity (Nodal) of POLYHEDRONS:"
    print ""
    connectivity = myMesh.getPolyhedronConnectivity(MED_NODAL)
    fIndex = myMesh.getPolyhedronFacesIndex()
    index = myMesh.getPolyhedronIndex(MED_NODAL)
    for j in range(nbPolyhedrons):
        print     "       Polyhedra",(j+1)
        iF1, iF2 = index[ j ]-1, index[ j+1 ]-1
        for f in range( iF2 - iF1 ):
            iN1, iN2 = fIndex[ iF1+f ]-1, fIndex[ iF1+f+1 ]-1
            print "         Face",f+1," ",connectivity[ iN1 : iN2 ]
            pass
        pass
    pass
