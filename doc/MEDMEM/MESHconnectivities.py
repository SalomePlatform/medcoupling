from libMEDMEM_Swig import *

MedFile = "pointe.med"
#MedFile = "carre_quad4_3.med"
meshName = "maa1"
#meshName = "CARRE_EN_QUAD4"

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
