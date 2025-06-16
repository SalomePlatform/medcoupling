#  -*- coding: utf-8 -*-
# Copyright (C) 2025  CEA, EDF
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

# EDF32699
# See https://ericca.uqtr.ca/fr11.7/man_u/u3/u3.01.00.pdf for documentation of format


def _getGeoTypeDict():
    import MEDLoader as ml

    dico = {
        "TRIA3": ml.NORM_TRI3,
        "TRIA6": ml.NORM_TRI6,
        "TRIA7": ml.NORM_TRI7,
        "QUAD4": ml.NORM_QUAD4,
        "QUAD8": ml.NORM_QUAD8,
        "QUAD9": ml.NORM_QUAD9,
        "SEG2": ml.NORM_SEG2,
        "SEG3": ml.NORM_SEG3,
        "SEG4": ml.NORM_SEG4,
        "POI1": ml.NORM_POINT1,
        "HEXA8": ml.NORM_HEXA8,
        "HEXA20": ml.NORM_HEXA20,
        "HEXA27": ml.NORM_HEXA27,
        "PENTA6": ml.NORM_PENTA6,
        "PENTA15": ml.NORM_PENTA15,
        "PENTA18": ml.NORM_PENTA18,
        "TETRA4": ml.NORM_TETRA4,
        "TETRA10": ml.NORM_TETRA10,
        "PYRAM5": ml.NORM_PYRA5,
        "PYRAM13": ml.NORM_PYRA13,
    }
    return dico


def _getReorderArray(mailGt: str):
    """
    This method gives components permutation array in order to goes from mail to med format convention
    :input mailGt: str of geometric type in mail format
    :return: list of components to reorder.
    """
    import MEDLoader as ml

    typesToReorder = {
        "PENTA15": [0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11],
        "PENTA18": [0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16, 17],
        "HEXA20": [
            0,
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            16,
            17,
            18,
            19,
            12,
            13,
            14,
            15,
        ],
        "HEXA27": [
            0,
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            16,
            17,
            18,
            19,
            12,
            13,
            14,
            15,
            20,
            21,
            22,
            23,
            24,
            25,
            26,
        ],
    }
    if mailGt not in typesToReorder:
        return list(
            range(
                ml.MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(
                    _getGeoTypeDict()[mailGt]
                )
            )
        )
    else:
        return typesToReorder[mailGt]


def _buildUmesh(coords, elements, title, groupesNo, groupesMa, mapGeoCellName):
    """
    Creates a MEDCouplingUMesh object from coordinate list and element connectivity.

    coords : list of [x,y(,z)] coordinates
    elements : dict of element_type → list of lists of node indices
    """
    import MEDLoader as ml

    # Convert node coordinates to medCoupling array
    coordsMC = ml.DataArrayDouble(coords)

    # Get the dimension of each geometric type
    allDims = set(
        [
            ml.MEDCouplingUMesh.GetDimensionOfGeometricType(_getGeoTypeDict()[elt])
            for elt in elements
        ]
    )

    # Get max dimension
    maxDim = max(
        [
            ml.MEDCouplingUMesh.GetDimensionOfGeometricType(_getGeoTypeDict()[elt])
            for elt in elements
        ]
    )
    # Create an offset dictionary for each geometry
    dicoOffset = {elt: 0 for elt in elements}

    mm = ml.MEDFileUMesh()

    dicoCellDim = {}  # Dictionary for geometric types and their dimensions
    mapGeoIdMm = {}  # Dictionary for geometric types and their mm IDs

    for dim in reversed(list(allDims)):
        # Filter the geo types of the current dimension
        gts = [
            (elt, _getGeoTypeDict()[elt])
            for elt in elements
            if ml.MEDCouplingUMesh.GetDimensionOfGeometricType(_getGeoTypeDict()[elt])
            == dim
        ]
        ms = []
        offset = 0
        for mailGt, mcGt in sorted(gts, key=lambda x: _getGeoTypeDict()[x[0]]):
            connQuad4 = ml.DataArrayInt(elements[mailGt])
            connQuad4 = connQuad4[:, _getReorderArray(mailGt)]
            connQuad4.rearrange(1)
            mQuad4 = ml.MEDCoupling1SGTUMesh("", mcGt)
            mQuad4.setCoords(coordsMC)
            mQuad4.setNodalConnectivity(connQuad4)
            ms.append(mQuad4.buildUnstructured())
            dicoOffset[mailGt] = offset
            i = offset
            offset += mQuad4.getNumberOfCells()
            dicoCellDim[mailGt] = dim
            IdCellMm = []
            while i < offset:
                IdCellMm.append(i)
                i = i + 1

            mapGeoIdMm[mailGt] = IdCellMm

        mm[dim - maxDim] = ml.MEDCouplingUMesh.MergeUMeshesOnSameCoords(ms)

    mapIdNameCell = {}
    for key in mapGeoIdMm:
        # Using get() to provide a default value if the key is missing
        list1 = mapGeoIdMm.get(key, [])
        list2 = mapGeoCellName.get(key, [])

        # Merge the lists with zip, but we manage the size difference if necessary
        mapIdNameCell[key] = list(zip(list2, list1))  # Pairs of corresponding elements

    mm.setName(title)

    dimPerCell = {}
    idMedPerCell = {}
    for gtMail in mapIdNameCell:
        for cellName, medId in mapIdNameCell[gtMail]:
            dimPerCell[cellName] = ml.MEDCouplingUMesh.GetDimensionOfGeometricType(
                _getGeoTypeDict()[gtMail]
            )
            idMedPerCell[cellName] = medId
    for cellGrp in groupesMa:
        elems = groupesMa[cellGrp]
        dimPerCellVect = [dimPerCell[elem] for elem in elems]
        idMedVect = [idMedPerCell[elem] for elem in elems]
        dimsOfCellGrp = set(dimPerCellVect)
        for dimOfCellGrp in dimsOfCellGrp:
            idsOfCellsWithDim = ml.DataArrayInt(dimPerCellVect).findIdsEqual(
                dimOfCellGrp
            )
            arr = ml.DataArrayInt(idMedVect)[idsOfCellsWithDim]
            arr.setName(cellGrp)
            mm.addGroup(dimOfCellGrp - maxDim, arr)

    for nodesGrp in groupesNo:
        elems = groupesNo[nodesGrp]
        arr = ml.DataArrayInt(elems)
        arr.setName(nodesGrp)
        arr.sort()
        mm.addGroup(1, arr)

    return mm


def _parseMeshFile(fichier):
    # Open and read all lines from the mesh file
    with open(fichier, "r") as f:
        lines = iter(f.readlines())  # Create an iterator over the lines of the file

    # Initialize variables and flags
    coords, nodesMap, elements = [], {}, {}
    # coords : List of node coordinates
    # nodesMap : ex: 'N1' → 0
    # element : Dict of element_type → list of connectivity lists
    groupesNo, groupesMa = {}, {}  # dict of nodes and mesh group
    mapGeoCellName = {}  # Map of element type → list of cells IDs
    title, dim = "", 0

    # Parsing state variables
    mode = None
    currentElemType = ""
    currentElemList = []
    currentGroup = []
    groupName = ""

    # Loop through each line in the file
    lineNum = 0
    for line in lines:
        lineNum += 1
        line = line.strip()
        if not line or line.startswith("#"):  # Skip empty lines and comments
            continue

        # Handle mode switching based on section headers
        if _is_mode_switch(line):
            mode, dim, currentElemType, currentElemList = _update_mode(
                line, elements, mapGeoCellName, dim
            )
            if mode in ["groupNo", "groupMa"]:
                groupName, currentGroup = line, []  # Reset group data
            continue
        elif line.startswith("FINSF"):
            _handleGroupEnd(mode, groupName, currentGroup, groupesNo, groupesMa)
            mode, groupName, currentGroup = None, "", []  # Reset mode and group data
            continue

        # Dispatch parsing based on current mode
        if mode == "title":
            title = line.strip()  # Store mesh title
        elif mode == "coords":
            _parseCoordsLine(line, dim, coords, nodesMap)  # Parse node coordinates
        elif mode == "elements":
            lineNum = _parseElementsLine(
                line, lines, currentElemType, elements, currentElemList, lineNum
            )  # Parse element connectivity
        elif mode == "groupNo" or mode == "groupMa":
            groupName, currentGroup = _parseGroupLine(
                line, groupName, currentGroup, lineNum
            )  # Parse group definitions

    # Transforming connectivities into integer indices
    for etype in elements:
        elements[etype] = [[nodesMap[n] for n in conn] for conn in elements[etype]]
    # Transform node groups into indices
    for group in groupesNo:
        groupesNo[group] = [nodesMap[n] for n in groupesNo[group]]
    #     # Transformer les groupes de mailles en indices
    # for groupName in groupesMa:
    #     groupesMa[groupName] = [int(m[1:]) for m in groupesMa[groupName] if m.startswith('M')]

    return coords, elements, title, groupesNo, groupesMa, mapGeoCellName


def _is_mode_switch(line):
    #  Check if the current line indicates a mode change (section header).
    return (
        line.startswith("TITRE")
        or line.startswith("COOR_")
        or line in _getGeoTypeDict()
        or line.startswith("GROUP_NO")
        or line.startswith("GROUP_MA")
    )


def _update_mode(line, elements, mapGeoCellName, dim):
    #  Update parsing mode and initialize relevant data structures.
    if line.startswith("TITRE"):
        return "title", dim, "", []
    elif line.startswith("COOR_2D"):
        dim = 2
        return "coords", dim, "", []
    elif line.startswith("COOR_3D"):
        dim = 3
        return "coords", dim, "", []
    elif line in _getGeoTypeDict():
        elements.setdefault(line, [])
        mapGeoCellName[line] = []
        return "elements", dim, line, mapGeoCellName[line]
    elif line.startswith("GROUP_NO"):
        return "groupNo", dim, "", []
    elif line.startswith("GROUP_MA"):
        return "groupMa", dim, "", []
    return None, dim, "", []


def _parseCoordsLine(line, dim, coords, nodesMap):
    # Parse a single coordinate line and update coords and node map.
    # if len(parts) < (1 + dim): genrate an error? ??
    parts = line.split()
    nid = parts[0]
    coord = list(map(float, parts[1 : 1 + dim]))
    nodesMap[nid] = len(coords)
    coords.append(coord)


def _parseElementsLine(line, lines, elementType, elements, cell_names, lineNum):
    # Parse a line from an element block, handling multi-line node connectivity.
    import MEDLoader as ml

    parts = line.split()
    if not parts:
        return
    name = parts[0]
    conn = parts[1:]
    nodes_needed = ml.MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(
        _getGeoTypeDict()[elementType]
    )
    while len(conn) < nodes_needed:
        conn += next(lines).strip().split()
        lineNum += 1
    elements[elementType].append(conn)
    cell_names.append(name)
    return lineNum


def _parseGroupLine(line, groupName, groupData, lineNum):
    #  Parse a line inside a GROUP_NO or GROUP_MA section.
    import re

    if groupName == "GROUP_NO" or groupName == "GROUP_MA":
        groupName = line
        return groupName, []
    elif groupName.startswith("GROUP_NO") or groupName.startswith("GROUP_MA"):
        pat = "[\s]*NOM[\s]*\=[\s]*"
        if not re.search(pat, groupName):
            raise ValueError(
                f" Group name is missing or badly formatted on line {lineNum - 1}"
            )
        parts = re.split(pat, groupName)
        groupName = parts[1].strip()
    groupData.extend(line.split())
    return groupName, groupData


def _handleGroupEnd(mode, groupName, groupData, groupesNo, groupesMa):
    # Finalize group and store it in the correct dictionary
    if not groupName:
        return
    if mode == "groupNo":
        groupesNo[groupName] = groupData[:]
    elif mode == "groupMa":
        groupesMa[groupName] = groupData[:]


def LoadMailFileInMEDFileUMeshInstance(inputMailFilePath: str):
    """
    :return: MEDFileUMesh instance representing MED file structure in memory
    """
    # Parse mesh file into Python data structures
    coords, elements, title, groupesNo, groupesMa, mapGeoCellName = _parseMeshFile(
        inputMailFilePath
    )

    # Build MEDCouplingUMesh from data structures
    return _buildUmesh(coords, elements, title, groupesNo, groupesMa, mapGeoCellName)


def ConvertFromMailToMEDFile(inputMailFilePath: str, outputMedFilePath: str):
    mm = LoadMailFileInMEDFileUMeshInstance(inputMailFilePath)
    mm.write(outputMedFilePath, 2)
    return outputMedFilePath
