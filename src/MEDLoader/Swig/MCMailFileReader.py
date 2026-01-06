#  -*- coding: utf-8 -*-
# Copyright (C) 2025-2026  CEA, EDF
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
        "TETRA4": [0, 2, 1, 3],
        "TETRA10": [0, 2, 1, 3, 6, 5, 4, 7, 9, 8],
        "PYRAM5": [0, 3, 2, 1, 4],
        "PYRAM13": [0, 3, 2, 1, 4, 8, 7, 6, 5, 9, 12, 11, 10],
        "PENTA6": [0, 2, 1, 3, 5, 4],
        "PENTA15": [0, 2, 1, 3, 5, 4, 8, 7, 6, 14, 13, 12, 9, 11, 10],
        "PENTA18": [0, 2, 1, 3, 5, 4, 8, 7, 6, 14, 13, 12, 9, 11, 10, 17, 16, 15],
        "HEXA8": [0, 3, 2, 1, 4, 7, 6, 5],
        "HEXA20": [
            0,
            3,
            2,
            1,
            4,
            7,
            6,
            5,
            11,
            10,
            9,
            8,
            19,
            18,
            17,
            16,
            12,
            15,
            14,
            13,
        ],
        "HEXA27": [
            0,
            3,
            2,
            1,
            4,
            7,
            6,
            5,
            11,
            10,
            9,
            8,
            19,
            18,
            17,
            16,
            12,
            15,
            14,
            13,
            20,
            24,
            23,
            22,
            21,
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
    import collections

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
    renumDim = {}  # Dictionary containing if necessary renumbering of cells o2n

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
            connGt = ml.DataArrayInt(elements[mailGt])
            connGt = connGt[:, _getReorderArray(mailGt)]
            connGt.rearrange(1)
            mGt = ml.MEDCoupling1SGTUMesh("", mcGt)
            mGt.setCoords(coordsMC)
            mGt.setNodalConnectivity(connGt)
            ms.append(mGt.buildUnstructured())
            dicoOffset[mailGt] = offset
            i = offset
            offset += mGt.getNumberOfCells()
            dicoCellDim[mailGt] = dim
            IdCellMm = []
            while i < offset:
                IdCellMm.append(i)
                i = i + 1

            mapGeoIdMm[mailGt] = IdCellMm

        # see EDF33583 for zzzz366a.mail
        um = ml.MEDCouplingUMesh.MergeUMeshesOnSameCoords(ms)

        o2n = um.sortCellsInMEDFileFrmt()

        if o2n:
            n2o = o2n.invertArrayO2N2N2O(um.getNumberOfCells())
            renumDim[dim] = n2o

        mm[dim - maxDim] = um

    mapIdNameCell = {}
    for key in mapGeoIdMm:
        # Using get() to provide a default value if the key is missing
        list1 = mapGeoIdMm.get(key, [])
        list2 = mapGeoCellName.get(key, [])

        # Merge the lists with zip, but we manage the size difference if necessary
        mapIdNameCell[key] = list(zip(list2, list1))  # Pairs of corresponding elements

    meshName = title[:64].strip()  # 64 == MED_NAME_SIZE

    mm.setName(meshName)

    dimPerCell = {}
    idMedPerCell = {}
    for gtMail in mapIdNameCell:
        for cellName, medId in mapIdNameCell[gtMail]:
            dimPerCell[cellName] = ml.MEDCouplingUMesh.GetDimensionOfGeometricType(
                _getGeoTypeDict()[gtMail]
            )
            idMedPerCell[cellName] = medId
    #
    grps = collections.defaultdict(list)
    #
    for cellGrp in groupesMa:
        elems = groupesMa[cellGrp]
        dimPerCellVect = [dimPerCell[elem] for elem in elems]
        idMedVect = [idMedPerCell[elem] for elem in elems]
        dimsOfCellGrp = set(dimPerCellVect)
        for dimOfCellGrp in dimsOfCellGrp:
            idsOfCellsWithDim = ml.DataArrayInt(dimPerCellVect).findIdsEqual(
                dimOfCellGrp
            )
            arr0 = ml.DataArrayInt(idMedVect)[idsOfCellsWithDim]
            arr = ml.DataArrayInt(sorted(arr0.getValues()))

            #
            if dimOfCellGrp in renumDim:
                if renumDim[dimOfCellGrp]:
                    arr.transformWithIndArr(renumDim[dimOfCellGrp])
                    arr.sort()
            #
            arr.setName(cellGrp)
            grps[dimOfCellGrp - maxDim].append(arr)

    for dimRelative in grps:
        arrs = grps[dimRelative]
        mm.setGroupsAtLevel(dimRelative, arrs)

    for nodesGrp in groupesNo:
        elems = groupesNo[nodesGrp]
        arr = ml.DataArrayInt(elems)
        arr.setName(nodesGrp)
        arr.sort()
        mm.addGroup(1, arr)

    return mm


def _zip_line(line):
    # Restituer la ligne sans champs et les champs à part
    spline = line.split()
    fields = dict(i.split("=") for i in spline if ("=" in i and len(i.split("=")) == 2))
    zipped = " ".join((i for i in spline if "=" not in i)).strip()
    return zipped, fields


def _zip_block(block):
    strip = {" =": "=", "= ": "="}
    # Supprimer l'ensemble des espaces avant et après le =
    # Afin d'identifier les champs
    while any(key in block for key in strip.keys()):
        for k, sk in strip.items():
            block = block.replace(k, sk)
    lines = (line for line in block.split("\n") if len(line) > 0)

    fields = {}
    zipped_lines = []
    for item in lines:
        zipped, flds = _zip_line(item)
        if len(zipped) > 0:
            zipped_lines.append(zipped)
        fields.update(flds)

    block_name = zipped_lines[0]
    block_body = " ".join(zipped_lines[1:]).split()
    return block_name, block_body, fields


def _remove_comments_and_split(fstream):
    comment = "%"
    txt = "\n".join(line.partition(comment)[0].strip() for line in fstream)
    return txt.split("FINSF")


def _parseCoordsBlock(bname, bbody):
    """Parse coordinates in a block"""
    space_dim = int(bname.strip("COOR_").strip("D"))
    n = space_dim + 1

    coords = []
    nodesMap = {}

    for i in range(0, len(bbody), n):
        nodesMap[bbody[i]] = len(coords)
        coor_str = [s.replace("D", "E").replace("d", "e") for s in bbody[i + 1 : i + n]]
        coords.append(list(map(float, coor_str)))

    return coords, nodesMap


def _getGrpName(bfields):
    """Get name for group"""
    try:
        grp_name = bfields["NOM"]
    except:
        val = bfields[list(bfields.keys())[0]]
        raise ValueError(f" Group name '{val}' is missing or badly formatted.")
    return grp_name


def _parseElementsBlock(bname, bbody):
    """Parse all elements in a block"""
    import MEDLoader as ml

    etype = bname
    nb_nodes = ml.MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(
        _getGeoTypeDict()[etype]
    )
    n = nb_nodes + 1

    elements = []
    cellsName = []

    for i in range(0, len(bbody), n):
        cellsName.append(bbody[i])
        elements.append(bbody[i + 1 : i + n])

    return elements, cellsName


def _parseMeshFile(filename):
    # Open and read all lines from the mesh file
    with open(filename, "r") as f:
        mesh_blocks = _remove_comments_and_split(f)

    # Initialize variables and flags
    coords, nodesMap, elements = [], {}, {}
    # coords : List of node coordinates
    # nodesMap : ex: 'N1' → 0
    # element : Dict of element_type → list of connectivity lists
    groupesNo, groupesMa = {}, {}  # dict of nodes and mesh group
    mapGeoCellName = {}  # Map of element type → list of cells IDs
    title = "mesh"

    for block in mesh_blocks:
        bname, bbody, bfields = _zip_block(block)

        if bname in ("COOR_3D", "COOR_2D"):
            coords, nodesMap = _parseCoordsBlock(bname, bbody)

        elif bname in ("GROUP_MA",):
            if len(bfields) > 0:
                grp_name = _getGrpName(bfields)
                groupesMa[grp_name] = bbody
            else:
                grp_name = bbody[0]
                groupesMa[grp_name] = bbody[1:]

        elif bname in ("GROUP_NO",):
            if len(bfields) > 0:
                grp_name = _getGrpName(bfields)
                groupesNo[grp_name] = bbody
            else:
                grp_name = bbody[0]
                groupesNo[grp_name] = bbody[1:]

        elif bname in _getGeoTypeDict():
            etype = bname
            new_elem, new_name = _parseElementsBlock(etype, bbody)

            elements.setdefault(etype, [])
            elements[etype].extend(new_elem)
            mapGeoCellName.setdefault(etype, [])
            mapGeoCellName[etype].extend(new_name)
        elif bname in ("TITRE"):
            if len(bbody) > 0:
                title = ""
                if len(bfields) > 0:
                    key = list(bfields.keys())[0]
                    title = key + " = " + bfields[key]
                for word in bbody:
                    title += " " + word
        else:
            pass

    # Transforming connectivities into integer indices
    for etype in elements:
        elements[etype] = [[nodesMap[n] for n in conn] for conn in elements[etype]]
    # Transform node groups into indices
    for group in groupesNo:
        groupesNo[group] = sorted([nodesMap[n] for n in groupesNo[group]])

    return coords, elements, title, groupesNo, groupesMa, mapGeoCellName


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
