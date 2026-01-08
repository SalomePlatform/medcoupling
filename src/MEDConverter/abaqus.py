#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (C) 2023-2026  CEA, EDF
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

import time
import os.path as osp
import numpy as np
import sys

from .logger import logger
from .MEDConverterMesh import MEDConverterMesh
from .cells import CellsTypeConverter
from .connectivity import ConnectivityRenumberer


class AbaqusNode:
    def __init__(self, node_id=-1, node_coordinates=[]):
        self.id = node_id
        self.coordinates = node_coordinates

    def __repr__(self):
        return "<Node> Id: {0}, Coordinates: {1}".format(self.id, self.coordinates)

    def __str__(self):
        return "<Node> Id: {0}, Coordinates: {1}".format(self.id, self.coordinates)

    def setId(self, node_id):
        self.id = node_id

    def getId(self):
        return self.id

    def setCoordinates(self, node_coordinates):
        self.coordinates = node_coordinates

    def getCoordinates(self):
        if len(self.coordinates) != 3:
            raise RuntimeError("Coordinates have to have 3 items")
        return self.coordinates


class AbaqusElement:
    def __init__(self, elem_type=None, elem_id=None, elem_nodes=None, multilevel=False):
        self.id = elem_id
        if elem_nodes is not None:
            self.nodes = elem_nodes
        else:
            self.nodes = []
        self.type = elem_type
        self.multilevel = multilevel

    def __repr__(self):
        return "<Element> Id: {0}, Type: {1}, Nodes: {2}".format(
            self.id, self.type, self.nodes
        )

    def __str__(self):
        return "<Element> Id: {0}, Type: {1}, Nodes: {2}".format(
            self.id, self.type, self.nodes
        )

    def setId(self, elem_id):
        self.id = elem_id

    def getId(self):
        return self.id

    def setNodes(self, elem_nodes):
        self.nodes = elem_nodes

    def getNodes(self):
        return self.nodes

    def setType(self, elem_type):
        self.type = elem_type

    def getType(self):
        return self.type


class AbaqusGroup:
    def __init__(self, name=None, instance=None, multilevel=False, group=None):
        self.name = name
        self.instance = instance
        self.multilevel = multilevel
        if group is not None:
            self.group = group
        else:
            self.group = []

    def __repr__(self):
        return "<Group> Name: {0}, Instance: {1}, Group: {2}".format(
            self.name, self.instance, self.group
        )

    def __str__(self):
        return "<Group> Name: {0}, Instance: {1}, Group: {2}".format(
            self.name, self.instance, self.group
        )

    def setName(self, name):
        self.name = name

    def getName(self):
        return self.name

    def setInstance(self, instance):
        self.instance = instance

    def getInstance(self):
        return self.instance

    def addGroup(self, group):
        self.group += group

    def setGroup(self, group):
        self.group = group

    def getGroup(self):
        return self.group


class AbaqusSurface:
    def __init__(self, name, type, elem):
        # Create connectivity for faces - the abaqus's element used is just
        # for the med conversion after and has no sense
        self.faces = {
            "TRI3": [3, ["B21", 2, [1, 2]], ["B21", 2, [2, 3]], ["B21", 2, [3, 1]]],
            "TRI6": [
                3,
                ["B22", 3, [1, 2, 4]],
                ["B22", 3, [2, 3, 5]],
                ["B22", 3, [3, 1, 6]],
            ],
            "QUAD4": [
                4,
                ["B21", 2, [1, 2]],
                ["B21", 2, [2, 3]],
                ["B21", 2, [3, 4]],
                ["B21", 2, [4, 1]],
            ],
            "QUAD8": [
                4,
                ["B22", 2, [1, 2, 5]],
                ["B22", 2, [2, 3, 6]],
                ["B22", 2, [3, 4, 7]],
                ["B22", 2, [4, 1, 8]],
            ],
            "QUAD9": [
                4,
                ["B22", 2, [1, 2, 5]],
                ["B22", 2, [2, 3, 6]],
                ["B22", 2, [3, 4, 7]],
                ["B22", 2, [4, 1, 8]],
            ],
            "TETRA4": [
                4,
                ["CPE3", 3, [1, 2, 3]],
                ["CPE3", 3, [1, 4, 2]],
                ["CPE3", 3, [2, 4, 3]],
                ["CPE3", 3, [3, 4, 1]],
            ],
            "TETRA10": [
                4,
                ["CPE6", 6, [1, 2, 3, 5, 6, 7]],
                ["CPE6", 6, [1, 4, 2, 8, 9, 5]],
                ["CPE6", 6, [2, 4, 3, 9, 10, 6]],
                ["CPE6", 6, [3, 4, 1, 10, 8, 7]],
            ],
            "HEXA8": [
                6,
                ["CPE4", 4, [1, 2, 3, 4]],
                ["CPE4", 4, [5, 8, 7, 6]],
                ["CPE4", 4, [1, 5, 6, 2]],
                ["CPE4", 4, [2, 6, 7, 3]],
                ["CPE4", 4, [3, 7, 8, 4]],
                ["CPE4", 4, [4, 8, 5, 1]],
            ],
            "HEXA20": [
                6,
                ["CPE8", 8, [1, 2, 3, 4, 9, 10, 11, 12]],
                ["CPE8", 8, [5, 8, 7, 6, 16, 15, 14, 13]],
                ["CPE8", 8, [1, 5, 6, 2, 17, 13, 18, 9]],
                ["CPE8", 8, [2, 6, 7, 3, 18, 14, 19, 10]],
                ["CPE8", 8, [3, 7, 8, 4, 19, 15, 20, 11]],
                ["CPE8", 8, [4, 8, 5, 1, 20, 16, 17, 12]],
            ],
            "HEXA27": [
                6,
                ["CPE9", 9, [1, 2, 3, 4, 9, 10, 11, 12, 22]],
                ["CPE9", 9, [5, 8, 7, 6, 16, 15, 14, 13, 23]],
                ["CPE9", 9, [1, 5, 6, 2, 17, 13, 18, 9, 24]],
                ["CPE9", 9, [2, 6, 7, 3, 18, 14, 19, 10, 25]],
                ["CPE9", 9, [3, 7, 8, 4, 19, 15, 20, 11, 26]],
                ["CPE9", 9, [4, 8, 5, 1, 20, 16, 17, 12, 27]],
            ],
            "PYRA5": [
                5,
                ["CPE4", 4, [1, 2, 3, 4]],
                ["CPE3", 3, [1, 5, 2]],
                ["CPE3", 3, [2, 5, 3]],
                ["CPE3", 3, [3, 5, 4]],
                ["CPE3", 3, [4, 5, 1]],
            ],
            "PENTA6": [
                5,
                ["CPE3", 3, [1, 2, 3]],
                ["CPE3", 3, [4, 6, 5]],
                ["CPE4", 4, [1, 4, 5, 2]],
                ["CPE4", 4, [2, 5, 6, 3]],
                ["CPE4", 4, [3, 6, 4, 1]],
            ],
            "PENTA15": [
                5,
                ["CPE6", 6, [1, 2, 3, 7, 8, 9]],
                ["CPE6", 6, [4, 6, 5, 12, 11, 10]],
                ["CPE8", 8, [1, 4, 5, 2, 13, 10, 14, 7]],
                ["CPE8", 8, [2, 5, 6, 3, 14, 11, 15, 8]],
                ["CPE8", 8, [3, 6, 4, 1, 15, 12, 13, 9]],
            ],
            "PENTA18": [
                5,
                ["CPE6", 6, [1, 2, 3, 7, 8, 9]],
                ["CPE6", 6, [4, 6, 5, 12, 11, 10]],
                ["CPE9", 9, [1, 4, 5, 2, 13, 10, 14, 7, 16]],
                ["CPE9", 9, [2, 5, 6, 3, 14, 11, 15, 8, 17]],
                ["CPE9", 9, [3, 6, 4, 1, 15, 12, 13, 9, 18]],
            ],
        }

        self.conv = CellsTypeConverter._abaqus_to_med
        self.name = name
        self.elem = elem
        self.type = type

    def getName(self):
        return self.name

    def getType(self):
        return self.type

    def getSurface(self):
        return self.elem

    def createElement(self, cell, faceId, newID):
        """Create AbaqusElement as a surface of a given cell"""

        med_type = self.conv[cell.getType()]
        face = self.faces[med_type][faceId]
        cell_nodes = cell.getNodes()
        index_nodes = [cell_nodes[node - 1] for node in face[2]]

        return AbaqusElement(face[0], newID, index_nodes, False)


class AbaqusPart:
    def __init__(self):
        self.name = " "
        self.Nodes = []
        self.Elements = []
        self.Elset = []
        self.Nset = []
        self.ElsetName = {}
        self.NsetName = {}
        self.Surfaces = []
        self.Parts = []

    def setName(self, name):
        self.name = name

    def getName(self):
        return self.name


class AbaqusInstance:
    def __init__(self):
        self.name = " "
        self.PartName = " "
        self.Nodes = []
        self.Elements = []
        self.Elset = []
        self.Nset = []
        self.ElsetName = {}
        self.NsetName = {}
        self.Surfaces = []
        self.translation = None
        self.rotation = None
        self.Parts = []

    def setName(self, name):
        self.name = name

    def getName(self):
        return self.name

    def setPartName(self, part_name):
        self.PartName = part_name

    def getPartName(self):
        return self.PartName

    def _addPart(self, Part):
        self.Nodes += Part.Nodes
        self.Elements += Part.Elements
        self.Elset += Part.Elset
        self.Nset += Part.Nset

        def check_name(dic1, dic2):
            for name in dic1.keys():
                if name in dic2.keys():
                    raise RuntimeError("Group %s already exists" % name)

        check_name(Part.ElsetName, self.ElsetName)
        self.ElsetName.update(Part.ElsetName)
        check_name(Part.NsetName, self.NsetName)
        self.NsetName.update(Part.NsetName)
        for part in Part.Parts:
            self._addPart(part)

    def setPart(self, Part):
        self.setPartName(Part.getName())
        self._addPart(Part)


class AbaqusAssembly:
    def __init__(self):
        self.name = " "
        self.Instance = []
        self.Nodes = []
        self.Elements = []
        self.Elset = []
        self.Nset = []
        self.ElsetName = {}
        self.NsetName = {}
        self.Parts = []
        self.Surfaces = []

    def setName(self, name):
        self.name = name

    def getName(self):
        return self.name

    def addInstance(self, Instance):
        self.Instance.append(Instance)

    def getInstance(self):
        return self.Instance

    def getPart(self, name):
        for part in self.Parts:
            if name == part.getName():
                return part

        raise RuntimeError("Part nod found: " + name)


class AbaqusNumbering:
    def __init__(self):
        self.name = " "
        self.corresponding_nodes = None
        self.corresponding_elems = None

    def setName(self, name):
        self.name = name

    def getName(self):
        return self.name


class AbaqusMesh:
    def __init__(self):
        self.name = " "
        self.Nodes = []
        self.Elements = []
        self.Elset = []
        self.Nset = []
        self.ElsetName = {}
        self.NsetName = {}
        self.nodesOffset = 0
        self.elemsOffset = 0
        self.surfOffset = 0
        self.Numbering = []

    def setName(self, name):
        self.name = name

    def getName(self):
        return self.name

    def translation(self, point, translation):
        return np.array(point) + np.array(translation)

    def matrix_rotation(self, axis, angle_radian):
        u = np.array(axis)
        u = u / np.linalg.norm(u)

        c = np.cos(angle_radian)
        s = np.sin(angle_radian)

        row1 = [
            u[0] * u[0] * (1 - c) + c,
            u[0] * u[1] * (1 - c) - u[2] * s,
            u[0] * u[2] * (1 - c) + u[1] * s,
        ]
        row2 = [
            u[0] * u[1] * (1 - c) + u[2] * s,
            u[1] * u[1] * (1 - c) + c,
            u[1] * u[2] * (1 - c) - u[0] * s,
        ]
        row3 = [
            u[0] * u[2] * (1 - c) - u[1] * s,
            u[1] * u[2] * (1 - c) + u[0] * s,
            u[2] * u[2] * (1 - c) + c,
        ]

        mrot = np.array([row1, row2, row3])

        return mrot

    def rotation(self, point, center, matrix_rotation):
        return matrix_rotation @ (np.array(point) - np.array(center)) + np.array(center)

    def geometric_transfo(
        self, point, translation=None, center=None, matrix_rotation=None
    ):
        if translation is not None:
            transla = self.translation(point, translation)
        else:
            transla = point

        if matrix_rotation is not None:
            rota = self.rotation(transla, center, matrix_rotation)
        else:
            rota = transla

        return tuple(rota)

    def addNodes(self, Nodes, translation=None, rotation_param=None):
        # object for translation + rotation (optimization for large mesh)
        if translation is not None:
            translation = np.array(translation)

        if rotation_param is not None:
            assert translation is not None
            center = np.array(rotation_param[0:3])
            axis = np.array(rotation_param[3:6])
            angle = rotation_param[6]
            angle_radian = np.radians(angle)
            mrot = self.matrix_rotation(axis, angle_radian)
        else:
            center, mrot = None, None

        corresponding_nodes = {}
        for idx, node in enumerate(Nodes):
            if node.getId() in corresponding_nodes:
                raise RuntimeError(
                    "Two nodes with identical id: {0}".format(node.getId())
                )
            else:
                corresponding_nodes[node.getId()] = self.nodesOffset + idx

            # add Node
            node_id = corresponding_nodes[node.getId()]
            coor = node.getCoordinates()
            new_coor = self.geometric_transfo(coor, translation, center, mrot)
            self.Nodes.append(AbaqusNode(node_id, new_coor))

        self.nodesOffset += len(Nodes)

        return corresponding_nodes

    def addElements(self, Elements, corresponding_nodes):
        corresponding_elems = {}
        for elem in Elements:
            self.elemsOffset += 1
            if elem.getId() in corresponding_elems:
                raise RuntimeError(
                    "Two elements with identical id: {0}".format(elem.getId())
                )
            else:
                corresponding_elems[elem.getId()] = self.elemsOffset

            if elem.multilevel:
                list_nodes = []
                for node in elem.getNodes():
                    try:
                        node_local_id = int(node)
                        list_nodes.append(corresponding_nodes[node_local_id])
                    except:
                        entries = [x.strip() for x in node.strip().split(".")]
                        assert len(entries) == 2
                        subinstance = entries[0]
                        node_local_id = int(entries[1])

                        l_find = False
                        for nume in self.Numbering:
                            if subinstance == nume.getName():
                                global_id = self.getGlobalId(
                                    nume.corresponding_nodes, node_local_id
                                )

                                if global_id is None:
                                    raise RuntimeError(
                                        "Create Element: node %s is \
                                        not in the mesh"
                                        % node
                                    )
                                else:
                                    list_nodes.append(global_id)

                                l_find = True
                                break

                        if not l_find:
                            raise RuntimeError(
                                "Create Element: node %s is \
                                    not in the mesh"
                                % {node}
                            )
            else:
                nodes_elem = map(int, elem.getNodes())
                try:
                    list_nodes = tuple(corresponding_nodes[k] for k in nodes_elem)
                except KeyError:
                    msg = (
                        "Element %s of type %s can not be converted (nodes not finded)"
                        % (
                            elem.getId(),
                            elem.getType(),
                        )
                    )
                    raise RuntimeError(msg)
            elem_id = corresponding_elems[elem.getId()]
            self.Elements.append(AbaqusElement(elem.getType(), elem_id, list_nodes))

        return corresponding_elems

    def addSurface(self, Surfaces, Elset, corresponding_elems):
        for surfs in Surfaces:
            elemSurf = []
            if surfs.getType() == "ELEMENT":
                for surf in surfs.getSurface():
                    l_global_grp = False
                    try:
                        listElem = [int(surf[0])]
                    except:
                        name_grp = surf[0]
                        l_find = False
                        for group in self.Elset:
                            if name_grp == group.getName():
                                listElem = group.getGroup()
                                l_find = True
                                l_global_grp = True
                                break
                        if not l_find:
                            for group in Elset:
                                if name_grp == group.getName():
                                    listElem = group.getGroup()
                                    l_find = True
                                    l_global_grp = False
                                    break
                        if not l_find:
                            raise RuntimeError("Group not find")

                    if len(surf) == 1 or surf[1] in ("SPOS", "SNEG"):
                        l_create_elem = False
                    else:
                        l_create_elem = True

                    for cell_loc_id in listElem:
                        if l_global_grp:
                            cell_id = cell_loc_id
                        else:
                            cell_id = corresponding_elems[cell_loc_id]

                        if l_create_elem:
                            self.surfOffset += 1
                            self.elemsOffset += 1
                            global_id = self.elemsOffset
                            cell = self.Elements[cell_id - 1]
                            surf_id = int(surf[1][1:])

                            self.Elements.append(
                                surfs.createElement(cell, surf_id, self.surfOffset)
                            )

                            if self.surfOffset in corresponding_elems:
                                raise RuntimeError(
                                    "Two elements with identical id: {0}".format(
                                        self.surfOffset
                                    )
                                )
                            else:
                                corresponding_elems[self.surfOffset] = self.elemsOffset
                        else:
                            global_id = cell_id

                        elemSurf.append(global_id)

                if surfs.getName() in self.ElsetName:
                    raise RuntimeError(
                        "Two surfaces with identical name: {0}".format(surfs.getName())
                    )
                self.Elset.append(AbaqusGroup(surfs.getName(), "xxx", False, elemSurf))
            else:
                logger.debug("Ignore SURFACE keyword")

    def fuseCommonGroup(self, Groups, GroupsName, Group):
        name = Group.getName()
        if name in GroupsName:
            Groups[GroupsName[name]].addGroup(Group.getGroup())
        else:
            Groups.append(Group)
            GroupsName[name] = len(Groups) - 1

    def getGlobalId(self, corresponding, local_id):
        if local_id in corresponding:
            global_id = corresponding[local_id]
        else:
            global_id = None

        return global_id

    def addGroups(self, typeGrp, Groups, corresponding):
        for Group in Groups:
            # add group
            name = Group.getName()
            instance = Group.getInstance()
            list_clean = []
            if Group.multilevel:
                for k in Group.getGroup():
                    entries = [x.strip() for x in k.strip().split(".")]
                    assert len(entries) == 2
                    subinstance = entries[0]
                    local_id = int(entries[1])
                    l_find = False
                    for nume in self.Numbering:
                        if subinstance == nume.getName():
                            if typeGrp == "NSET":
                                global_id = self.getGlobalId(
                                    nume.corresponding_nodes, local_id
                                )
                            elif typeGrp == "ELSET":
                                global_id = self.getGlobalId(
                                    nume.corresponding_elems, local_id
                                )
                            else:
                                raise RuntimeError("Unknown type of group")

                            if global_id is None:
                                raise RuntimeError(
                                    "Create Group %s: element %s is \
                                    not in the mesh"
                                    % {name, k}
                                )
                            else:
                                list_clean.append(global_id)

                            l_find = True
                            break

                    if not l_find:
                        raise RuntimeError(
                            "Create Group %s: element %s is \
                                    not in the mesh"
                            % {name, k}
                        )

                list_item = tuple(int(k) for k in list_clean)
            else:
                items = map(int, Group.getGroup())
                for k in items:
                    if k in corresponding:
                        list_clean.append(k)
                    else:
                        raise RuntimeError(
                            "Create Group %s: element %s is \
                                    not in the mesh"
                            % {name, k}
                        )

                list_item = tuple(corresponding[k] for k in list_clean)

            if len(list_item) == 0:
                raise RuntimeError("No items in group: " + name)

            if typeGrp == "NSET":
                self.fuseCommonGroup(
                    self.Nset,
                    self.NsetName,
                    AbaqusGroup(name, instance, False, list_item),
                )
            elif typeGrp == "ELSET":
                self.fuseCommonGroup(
                    self.Elset,
                    self.ElsetName,
                    AbaqusGroup(name, instance, False, list_item),
                )
            else:
                raise RuntimeError("Unknown type of group")

    def addFromEntities(self, Entities, translation=None, rotation_param=None):
        tic = time.perf_counter()
        corresponding_nodes = self.addNodes(Entities.Nodes, translation, rotation_param)
        toc = time.perf_counter()
        logger.debug(
            "-> Number of nodes : %d (in %0.4f seconds)"
            % (len(Entities.Nodes), toc - tic)
        )

        tic = time.perf_counter()
        corresponding_elems = self.addElements(Entities.Elements, corresponding_nodes)
        self.addSurface(Entities.Surfaces, Entities.Elset, corresponding_elems)
        toc = time.perf_counter()
        logger.debug(
            "-> Number of elements : %d (in %0.4f seconds)"
            % (len(Entities.Elements), toc - tic)
        )

        tic = time.perf_counter()
        self.addGroups("NSET", Entities.Nset, corresponding_nodes)
        toc = time.perf_counter()
        logger.debug(
            "-> Number of groups of nodes : %d (in %0.4f seconds)"
            % (len(Entities.Nset), toc - tic)
        )

        tic = time.perf_counter()
        self.addGroups("ELSET", Entities.Elset, corresponding_elems)
        toc = time.perf_counter()
        logger.debug(
            "-> Number of groups of elements : %d (in %0.4f seconds)"
            % (len(Entities.Elset), toc - tic)
        )

        localNumbering = AbaqusNumbering()
        localNumbering.setName(Entities.getName())
        localNumbering.corresponding_nodes = corresponding_nodes
        localNumbering.corresponding_elems = corresponding_elems
        self.Numbering.append(localNumbering)

    def addGroupsInRightPlace(self, Assembly):
        new_Nset = []
        for group in Assembly.Nset:
            instance_name = group.getInstance()
            if instance_name != "":
                for Instance in Assembly.Instance:
                    if Instance.getName() == instance_name:
                        Instance.Nset.append(group)
            else:
                new_Nset.append(group)

        Assembly.Nset = new_Nset

        new_Elset = []
        for group in Assembly.Elset:
            instance_name = group.getInstance()
            if instance_name != "":
                for Instance in Assembly.Instance:
                    if Instance.getName() == instance_name:
                        Instance.Elset.append(group)
            else:
                new_Elset.append(group)

        Assembly.Elset = new_Elset

    def assemble(self, Assembly):
        if len(Assembly.Parts) > 0:
            if len(Assembly.Instance) == 0:
                # create an instance with all parts
                for part in Assembly.Parts:
                    Instance = AbaqusInstance()

                    Instance.setName("Instance_" + part.getName())
                    Instance.setPart(part)

                    Assembly.addInstance(Instance)

        self.surfOffset = self._estimateNbElem(Assembly)
        self.addGroupsInRightPlace(Assembly)

        # loop on instance of Assembly
        for Instance in Assembly.Instance:
            logger.debug("Processing Instance: " + Instance.getName())
            self.addFromEntities(Instance, Instance.translation, Instance.rotation)

        # add others objects in assembly
        logger.debug("Processing rest of the mesh: ")
        self.addFromEntities(Assembly)

    def _estimateNbElem(self, Assembly):
        nb_elem = 0

        # loop on instance of Assembly
        for Instance in Assembly.Instance:
            nb_elem += len(Instance.Elements)

        nb_elem += len(Assembly.Elements)

        return nb_elem


class MEDConverterAbaqus(MEDConverterMesh):
    @staticmethod
    def convert_abaqus_to_med(filename_abaqus, verbose=False):
        tic = time.perf_counter()
        c = MEDConverterAbaqus()
        c.verbose = verbose
        c.read_abaqus_mesh(filename_abaqus)
        c.create_UMesh()
        toc = time.perf_counter()
        logger.debug("Mesh converted (in %0.4f seconds)" % (toc - tic))
        return c.umesh

    def __init__(self):
        super(MEDConverterAbaqus, self).__init__()
        self.abaqusmesh = None
        self.nbAssembly = 0

    def read_abaqus_mesh(self, filename):
        logger.debug("Reading ABAQUS mesh file : %s" % filename)
        self._reset_structures()

        Assembly = AbaqusAssembly()

        # Lecture du fichier .inp où les blocs sont separés par des *Instance et *End Instance
        with open(filename, "r", encoding=self._get_file_encoding(filename)) as file:
            self.file = []
            self.filename = filename
            self.mesh_name = self._read_meshname(filename)
            # a priori, this is a 3D mesh
            self.space_dim = 3

            logger.debug("Mesh name : %s" % self.mesh_name)
            logger.debug("Space Dimension : 3")
            logger.debug("Beginning to parse mesh file")
            tic = time.perf_counter()

            recur_level_default = sys.getrecursionlimit()
            sys.setrecursionlimit(10000)

            for line in file:
                self.line = line
                self._read_data(file, Assembly)

            sys.setrecursionlimit(recur_level_default)

            # close open file
            for fich in self.file:
                fich.close()

            toc = time.perf_counter()
            logger.debug("Ending to parse mesh file in %0.4f seconds" % (toc - tic))

        # create Abaqus mesh
        logger.debug(" ")
        logger.debug("Creating ABAQUS mesh:")
        tic = time.perf_counter()

        mesh = AbaqusMesh()
        mesh.setName(self.mesh_name)
        mesh.assemble(Assembly)

        del Assembly

        toc = time.perf_counter()

        logger.debug("End creating ABAQUS mesh in %0.4f seconds" % (toc - tic))

        logger.debug("Statistics of the mesh : " + mesh.getName())
        logger.debug("-> Number of nodes : %d" % (len(mesh.Nodes)))
        logger.debug("-> Number of elements : %d" % (len(mesh.Elements)))
        logger.debug("-> Number of groups of nodes : %d" % (len(mesh.Nset)))
        logger.debug("-> Number of groups of elements : %d" % (len(mesh.Elset)))

        # nodes of the mesh (collection of double)
        for node in mesh.Nodes:
            self.add_node(node.getId(), node.getCoordinates())

        # Les elements, triés par dimension
        e_conv = CellsTypeConverter("ABAQUS")
        c_renum = ConnectivityRenumberer("ABAQUS")

        groups_by_elem = {}

        for elem in mesh.Elements:
            idx_element_abaqus = elem.getId()
            element_abaqus_type = elem.getType()
            elements_nodes_abaqus = tuple(map(int, elem.getNodes()))

            element_medcoupling_type = e_conv.external_to_medcoupling(
                element_abaqus_type
            )
            element_nodes_med = c_renum.external_to_medcoupling(
                element_medcoupling_type, elements_nodes_abaqus
            )

            self.add_cell(
                idx_element_abaqus, element_medcoupling_type, element_nodes_med
            )

            if element_abaqus_type not in groups_by_elem:
                groups_by_elem[element_abaqus_type] = []
            groups_by_elem[element_abaqus_type].append(idx_element_abaqus)

        # Les groups
        # Nodes' group
        for group in mesh.Nset:
            group_name = group.getName()
            group_nodes_abaqus = map(int, group.getGroup())
            self.add_group_nodes(group_name, group_nodes_abaqus)

        # Element's group
        for group in mesh.Elset:
            group_name = group.getName()
            group_element_abaqus = map(int, group.getGroup())
            self.add_group_cells(group_name, group_element_abaqus)

        # Add a group by cell type:
        for group_name, group_element_abaqus in groups_by_elem.items():
            if group_name not in mesh.ElsetName:
                self.add_group_cells(group_name, group_element_abaqus)

    def _read_meshname(self, filename):
        return osp.splitext(osp.basename(filename))[0]

    def _read_data(self, file, Entities):
        self.line = self.line.strip()
        # print("LINE: ", self.line )
        if self.line.upper().startswith("*NODE"):
            self._read_nodes(file, Entities.Nodes, Entities.Nset, Entities.NsetName)
            # print("Nodes")
            # print(Entities.Nodes)
            self._read_data(file, Entities)
        elif self.line.upper().startswith("*ELEMENT"):
            self._read_cells(
                file, Entities.Elements, Entities.Elset, Entities.ElsetName
            )
            # print("Cells")
            # print(Entities.Elements)
            self._read_data(file, Entities)
        elif self.line.upper().startswith("*NSET"):
            self._read_group(file, "NSET", Entities.Nset, Entities.NsetName)
            # print("Nset")
            # print(Entities.Nset)
            self._read_data(file, Entities)
        elif self.line.upper().startswith("*ELSET"):
            self._read_group(file, "ELSET", Entities.Elset, Entities.ElsetName)
            # print("Elset")
            # print(Entities.Elset)
            self._read_data(file, Entities)
        elif self.line.upper().startswith("*INCLUDE"):
            self._read_include_file(file, Entities)
        elif self.line.upper().startswith("*INSTANCE"):
            self._read_instance(file, Entities)
        elif self.line.upper().startswith("*PART"):
            self._read_part(file, Entities.Parts)
        elif self.line.upper().startswith("*ASSEMBLY"):
            self._read_assembly(file, Entities)
            self.nbAssembly += 1

            if self.nbAssembly > 1:
                raise RuntimeError("Only one Assembly allowed")
        elif self.line.upper().startswith("*SURFACE"):
            self._read_surfaces(file, Entities.Surfaces)
            self._read_data(file, Entities)
        elif self.line.upper().startswith("*NGEN"):
            raise RuntimeError("Keyword not supported: NGEN")
        elif self.line.upper().startswith("*NFILL"):
            raise RuntimeError("Keyword not supported: NFILL")
        elif self.line.upper().startswith("*NMAP"):
            raise RuntimeError("Keyword not supported: NMAP")
        elif self.line.upper().startswith("*NCOPY"):
            raise RuntimeError("Keyword not supported: NCOPY")

    def _read_nodes(self, file, Nodes, Nset, NsetName):
        # get informations about nodes
        params_map = self._get_param_map(self.line)

        # this is not a list of node
        if "*NODE" not in params_map:
            self.line = file.readline()
            return

        logger.debug("-> Reading Nodes")

        # create directly a group from the list of nodes
        if "NSET" in params_map:
            create_nset = True
            list_nodes = []
        else:
            create_nset = False

        # the coordinates are in an external file
        if "INPUT" in params_map:
            l_extern_file = True
            filename_node = osp.dirname(self.filename) + "/" + params_map["INPUT"]
            file_to_read = open(filename_node, "r")
            self.file.append(file_to_read)
        else:
            l_extern_file = False
            file_to_read = file

        # loop on nodes
        while True:
            self.line = file_to_read.readline()
            l_process_line = True
            if l_extern_file:
                if self.line.startswith("*"):
                    if self.line.upper().startswith("*NODE"):
                        self._read_nodes(file_to_read, Nodes, Nset, NsetName)
                    else:
                        l_process_line = False

                if self.line == "":
                    break
                elif self.line in ["\n", "\r\n"]:
                    break
            elif self.breakLoop(self.line):
                break

            if not self.line.startswith("**"):
                if l_process_line:
                    entries = self._read_continuous_line(file_to_read, ",")
                    # read id and coordinatines
                    nid, x = int(entries[0]), entries[1:]
                    # fill with zero if not enougth coordinates
                    if len(x) < 3:
                        for i in range(0, 3 - len(x)):
                            x.append("0.0")
                    assert len(x) == 3

                    Nodes.append(AbaqusNode(nid, [float(xx) for xx in x]))

                    # add node in the group
                    if create_nset:
                        list_nodes.append(nid)

                    if self.line.lstrip().startswith("*"):
                        break

        # add group in Nset
        if create_nset:
            name = params_map["NSET"]
            if name in NsetName:
                Nset[NsetName[name]].addGroup(list_nodes)
            else:
                Nset.append(AbaqusGroup(name, "", False, list_nodes))
                NsetName[name] = len(Nset) - 1

        if l_extern_file:
            file_to_read.close()
            self.line = file.readline()

    # Read a list of element
    def _read_cells(self, file, Elements, Elset, ElsetName):
        # get informations about elements
        params_map = self._get_param_map(self.line)

        # this is not a list of element
        if "*ELEMENT" not in params_map:
            self.line = file.readline()
            return

        logger.debug("-> Reading Elements : " + params_map["TYPE"])

        # create directly a group from the list of elements
        if "ELSET" in params_map:
            create_elset = True
            list_elem = []
        else:
            create_elset = False

        # the elements are in an external file
        if "INPUT" in params_map:
            l_extern_file = True
            filename_elem = osp.dirname(self.filename) + "/" + params_map["INPUT"]
            file_to_read = open(filename_elem, "r")
            self.file.append(file_to_read)
        else:
            l_extern_file = False
            file_to_read = file

        # get type of element to create
        etype = params_map["TYPE"].upper()

        # loop on list of elements
        while True:
            self.line = file_to_read.readline()

            l_process_line = True
            if l_extern_file:
                if self.line.startswith("*"):
                    if self.line.upper().startswith("ELEMENT"):
                        self._read_cells(file_to_read, Elements, Elset, ElsetName)
                    else:
                        l_process_line = False

                if self.line == "":
                    break
                elif self.line in ["\n", "\r\n"]:
                    break
            elif self.breakLoop(self.line):
                break

            multilevel = False
            if not self.line.startswith("**"):
                if l_process_line:
                    entries = self._read_continuous_line(file_to_read, ",")
                    # get id and list of nodes
                    eid, nodes = int(entries[0]), entries[1:]
                    try:
                        index_nodes = [int(n) for n in nodes]
                    except:
                        index_nodes = nodes
                        multilevel = True
                    # add element
                    Elements.append(AbaqusElement(etype, eid, index_nodes, multilevel))
                    # add element in the group
                    if create_elset:
                        list_elem.append(eid)

                    if self.line.lstrip().startswith("*"):
                        break

        # add group in Elset
        if create_elset:
            name = params_map["ELSET"]

            if name in ElsetName:
                Elset[ElsetName[name]].addGroup(list_elem)
            else:
                Elset.append(AbaqusGroup(name, "", False, list_elem))
                ElsetName[name] = len(Elset) - 1

        if l_extern_file:
            file_to_read.close()
            self.line = file.readline()

    # Read a list of element
    def _read_surfaces(self, file, Surfaces):
        # get informations about elements
        params_map = self._get_param_map(self.line)

        # this is not a list of element
        if "TYPE" not in params_map:
            logger.debug("TYPE is not present for SURFACE: Ignore keyword")
            self.line = file.readline()
            return

        logger.debug(
            "-> Reading Surface : %s (%s)" % (params_map["NAME"], params_map["TYPE"])
        )

        elem = []
        # loop on list of elements
        while True:
            self.line = file.readline()

            if self.line.lstrip().startswith("*"):
                break

            elem.append(self._read_continuous_line(file, ","))

        Surfaces.append(
            AbaqusSurface(params_map["NAME"], params_map["TYPE"].upper(), elem)
        )

    def _read_group(self, file, typyeGroup, Group, GroupName):
        # find type of element
        params_map = self._get_param_map(self.line)

        # this is not a group
        if typyeGroup not in params_map:
            self.line = file.readline()
            return

        # to generate groups
        if "GENERATE" in params_map.keys():
            generate = True
        else:
            generate = False

        # to generate groups
        if "INSTANCE" in params_map.keys():
            instance = params_map["INSTANCE"]
        else:
            instance = ""

        logger.debug(
            "-> Reading Group: " + params_map[typyeGroup] + " (" + typyeGroup + ")"
        )

        list_item = []
        multilevel = False
        while True:
            self.line = file.readline()
            if self.breakLoop(self.line):
                break

            if not self.line.startswith("**"):
                entries = [x.strip() for x in self.line.strip().rstrip(",").split(",")]

                try:
                    int(entries[0])
                    l_list_grp = False
                except:
                    l_list_grp = True
                if l_list_grp:
                    l_find = False
                    for grp_name in entries:
                        try:
                            index_group = GroupName[grp_name]
                            list_item += Group[index_group].getGroup()
                            l_find = True
                        except:
                            pass

                    if not l_find:
                        multilevel = True
                        list_item += entries
                else:
                    if generate:
                        # default value is 1
                        if len(entries) == 2:
                            entries.append("1")
                        # generate elements in group
                        # first element, last_element, step

                        list_item += [
                            int(n)
                            for n in range(
                                int(entries[0]), int(entries[1]) + 1, int(entries[2])
                            )
                        ]
                    else:
                        # read directely list of elements
                        list_item += [int(n) for n in entries]

        if len(list_item) == 0:
            raise RuntimeError("No items for this group: " + params_map[typyeGroup])
        # add group
        name = params_map[typyeGroup]
        if name in GroupName:
            if Group[GroupName[name]].getInstance() == instance:
                Group[GroupName[name]].addGroup(list_item)
            else:
                Group.append(AbaqusGroup(name, instance, multilevel, list_item))
                GroupName[name] = [GroupName[name], len(Group) - 1]
        else:
            Group.append(AbaqusGroup(name, instance, multilevel, list_item))
            GroupName[name] = len(Group) - 1

    # Read an included file
    def _read_include_file(self, file, Entities):
        # get informations about file
        params_map = self._get_param_map(self.line)

        # this is not an included file
        if "INPUT" not in params_map:
            return

        # open external file
        filename_elem = osp.dirname(self.filename) + "/" + params_map["INPUT"]
        file_to_read = open(filename_elem, "r")
        self.file.append(file_to_read)

        logger.debug("-> Reading included file: " + filename_elem)

        # read external file
        for self.line in file_to_read:
            self._read_data(file_to_read, Entities)

        file_to_read.close()

    def _read_part(self, file, Parts):
        # get informations about part
        params_map = self._get_param_map(self.line)

        # this is not an included file
        if "*PART" not in params_map:
            return

        Part = AbaqusPart()

        Part.setName(params_map["NAME"])

        logger.debug("-> Reading Part: " + Part.getName())

        l_finish = False
        for line in file:
            self.line = line
            self._read_data(file, Part)
            if self.line.strip().upper().startswith("*END PART"):
                l_finish = True
                break

        if not l_finish:
            raise RuntimeError("Not Find: End Part")

        Parts.append(Part)

    def _read_assembly(self, file, Assembly):
        # get informations about part
        params_map = self._get_param_map(self.line)

        # this is not an included file
        if "*ASSEMBLY" not in params_map:
            return

        Assembly.setName(params_map["NAME"])

        l_finish = False
        for line in file:
            self.line = line
            self._read_data(file, Assembly)
            if self.line.strip().upper().startswith("*END ASSEMBLY"):
                l_finish = True
                break

        if not l_finish:
            raise RuntimeError("Not Find: End Assembly")

    def _read_instance(self, file, Assembly):
        # get informations about part
        params_map = self._get_param_map(self.line)

        # this is not an included file
        if "*INSTANCE" not in params_map:
            return

        Instance = AbaqusInstance()

        Instance.setName(params_map["NAME"])
        Instance.setPart(Assembly.getPart(params_map["PART"]))

        logger.debug("-> Reading Instance: " + Instance.getName())

        l_finish = False
        l_first_line = True
        for line in file:
            self.line = line
            if l_first_line:
                # read translation
                if not self.line.strip().startswith("*"):
                    Instance.translation = [
                        float(x.strip())
                        for x in self.line.strip().rstrip(",").split(",")
                    ]
                    assert len(Instance.translation) == 3
                    self.line = file.readline()
                    if not self.line.strip().startswith("*"):
                        Instance.rotation = [
                            float(x.strip())
                            for x in self.line.strip().rstrip(",").split(",")
                        ]
                        assert len(Instance.rotation) == 7
                        self.line = file.readline()

                l_first_line = False

            self._read_data(file, Instance)

            if self.line.strip().upper().startswith("*END INSTANCE"):
                l_finish = True
                break

        if not l_finish:
            raise RuntimeError("Not Find: End Instance")

        Assembly.addInstance(Instance)

    def _get_param_map(self, word, required_keys=None):
        """
        get the optional arguments on a line
        Example
        -------
        >>> word = 'elset,instance=dummy2,generate'
        >>> params = get_param_map(word, required_keys=['instance'])
        params = {
            'elset' : None,
            'instance' : 'dummy2,
            'generate' : None,
        }
        """
        if required_keys is None:
            required_keys = []
        words = word.split(",")
        param_map = {}
        for wordi in words:
            if "=" not in wordi:
                key = wordi.strip().upper()
                value = None
            else:
                sword = wordi.split("=")
                assert len(sword) == 2, sword
                key = (sword[0].strip()).upper()
                value = sword[1].strip()
            param_map[key] = value

        msg = ""
        for key in required_keys:
            if key not in param_map:
                msg += "%r not found in %r\n" % (key, word)
        if msg:
            raise RuntimeError(msg)
        return param_map

    # read a string which are in more that one line. If terminates by separator
    def _read_continuous_line(self, file, separator):
        # read the line
        entries = [
            x.strip() for x in self.line.strip().rstrip(separator).split(separator)
        ]

        # more than one line to read
        if self.line.rstrip().endswith(separator):
            self.line = file.readline()

            if not self.line.lstrip().startswith("*"):
                entries += self._read_continuous_line(file, separator)

        return entries

    def breakLoop(self, line):
        if line.startswith("*") and not line.startswith("**"):
            return True
        elif line == "":
            return True
        elif line in ["\n", "\r\n"]:
            return True

        return False

    def splitAndCleanLine(self, my_string, separator):
        return [x.strip() for x in my_string.split(separator)]
