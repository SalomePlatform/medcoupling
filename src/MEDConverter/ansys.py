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
from collections import OrderedDict

from .logger import logger
from .MEDConverterMesh import MEDConverterMesh
from .cells import CellsTypeConverter
from .connectivity import ConnectivityRenumberer

# doc: https://www.mm.bme.hu/~gyebro/files/ans_help_v182/ans_cmd/Hlp_C_CM.html
# doc: http://oss.jishulink.com/caenet/forums/upload/2013/11/25/389/21437609302438.pdf


class AnsysCell:
    def __init__(
        self,
        elem_type=None,
        elem_id=None,
        elem_nodes=None,
        sec_id=None,
        elem_rep=None,
        elem_const=None,
        elem_tension=None,
    ):
        self.id = elem_id
        if elem_nodes != None:
            self.nodes = elem_nodes
        else:
            self.nodes = []
        self.type = elem_type
        self.sec = sec_id
        self.rep = elem_rep
        self.const = elem_const
        self.tension = elem_tension

    def __repr__(self):
        return "<Cell> Id: {0}, Type: {1}, Nodes: {2}".format(
            self.id, self.type, self.nodes
        )

    def __str__(self):
        return "<Cell> Id: {0}, Type: {1}, Nodes: {2}".format(
            self.id, self.type, self.nodes
        )


class AnsysGroup:
    def __init__(self, name=None, typeg=None, group=None):
        self.name = name
        self.type = typeg
        if group != None:
            self.elems = group
        else:
            self.elems = []

    def __repr__(self):
        return "<Group> Name: {0}, Instance: {1}, Group: {2}".format(
            self.name, self.type, self.elems
        )

    def __str__(self):
        return "<Group> Name: {0}, Instance: {1}, Group: {2}".format(
            self.name, self.type, self.elems
        )


class Section:
    def __init__(self, sec_type, sec_subtype, sec_data=None, option=0, courbure=0.0):
        self.type = sec_type
        self.subtype = sec_subtype

    def setData(self, sec_data):
        self.data = sec_data


class Repere:
    def __init__(self, rep_type):
        self.type = rep_type
        self.orig = []
        self.angle = []

    def getRep(self, coord_nodes):
        if self.type == "CS":
            x1 = coord_nodes[self.angle[0]][0] - coord_nodes[self.orig[0]][0]
            x2 = coord_nodes[self.angle[0]][1] - coord_nodes[self.orig[0]][1]
            x3 = coord_nodes[self.angle[0]][2] - coord_nodes[self.orig[0]][2]
            y1 = coord_nodes[self.angle[1]][0] - coord_nodes[self.orig[0]][0]
            y2 = coord_nodes[self.angle[1]][1] - coord_nodes[self.orig[0]][1]
            y3 = coord_nodes[self.angle[1]][2] - coord_nodes[self.orig[0]][2]
            return (x1, x2, x3, y1, y2, y3)

        elif self.type == "LOCAL" or "CLOCAL":
            return (self.angle[0], self.angle[2], self.angle[1])


class MEDConverterAnsys(MEDConverterMesh):
    @staticmethod
    def convert_ansys_to_med(filename_ansys, verbose=False):
        tic = time.perf_counter()
        c = MEDConverterAnsys()
        c.verbose = verbose
        c.read_ansys_mesh(filename_ansys)
        c.create_UMesh()
        toc = time.perf_counter()
        logger.debug("Mesh converted (in %0.4f seconds)" % (toc - tic))
        return c.umesh

    def __init__(self):
        super(MEDConverterAnsys, self).__init__()
        self.ansysmesh = None

    def read_ansys_mesh(self, filename):
        logger.debug("Read ANSYS mesh.")

        self._reset_structures()
        Cells, Groups, title = [], [], None
        (nb_total_cells, last_idx_sec) = (0, 0)
        time_nodes, time_cell, time_groups = (0.0, 0.0, 0.0)
        Esel = []
        ElemEntities = {}
        ElemAnsys = {}
        ElemOpt = {}
        Sect = {}
        nodes = {}
        GROUPSMODELE = {}
        CMELEM = {}

        tic = time.perf_counter()
        # Lecture du fichier .cdb où les blocs sont separés par des BEGIN_* et END_*
        with open(filename, "r", encoding=self._get_file_encoding(filename)) as file:
            for line in file:
                strip_line = line.strip().upper()

                if strip_line.startswith("NBLOCK"):
                    tic0 = time.perf_counter()
                    nline = line.split(",")
                    dim = int(nline[1])
                    if dim == 6:
                        self.space_dim = 3
                    else:
                        self.space_dim = dim
                    nodes = self.__read_nodes(file, nodes)
                    toc0 = time.perf_counter()
                    time_nodes += toc0 - tic0
                elif strip_line.startswith("EBLOCK"):
                    tic0 = time.perf_counter()
                    nb_total_cells += int(line.split(",")[4])
                    self.__read_cells(file, Cells)
                    toc0 = time.perf_counter()
                    time_cell += toc0 - tic0
                elif strip_line.startswith("CMBLOCK"):
                    tic0 = time.perf_counter()
                    self.__read_groups(file, line, Groups)
                    toc0 = time.perf_counter()
                    time_groups += toc0 - tic0
                elif strip_line.startswith("ET,"):
                    sspline = strip_line.split(",")
                    ElemAnsys[int(sspline[1])] = int(sspline[2])
                elif strip_line.startswith("ETBLOCK"):
                    ElemAnsys = self.__read_elem(file, ElemAnsys)
                elif strip_line.startswith("ESEL,"):
                    sspline = strip_line.split(",")
                    assert sspline[1] in ("S", "ALL", "A")
                    if sspline[1] in ("S", "ALL"):
                        Esel = []

                    if sspline[1] in ("S", "A") and sspline[2] == "TYPE":
                        if len(sspline) == 6:
                            for i in range(int(sspline[4]), int(sspline[5]) + 1):
                                Esel.append(i)
                        else:
                            Esel.append(int(sspline[4]))
                elif strip_line.startswith("CM,"):
                    sspline = strip_line.split(",")
                    cname = sspline[1]
                    entity = sspline[2]
                    if entity == "ELEM":
                        CMELEM[cname] = Esel
                elif strip_line.startswith("KEYOP"):
                    sspline = strip_line.split(",")
                    ElemOpt[int(sspline[1])] = [int(sspline[2]), int(sspline[3])]
                elif strip_line.startswith("SECTYPE"):
                    sspline = strip_line.split(",")
                    Sect[int(sspline[1])] = Section(sspline[2], sspline[3].strip())
                    last_idx_sec = int(sspline[1])
                elif strip_line.startswith("SECDATA"):
                    sspline = strip_line.split(",")
                    if sspline[-1] == "":
                        sspline.pop(-1)
                    data = [float(i) for i in sspline[1:]]
                    Sect[last_idx_sec].setData(data)
                elif strip_line.startswith("SECBLOCK"):
                    sspline = strip_line.split(",")
                    for i in range(int(sspline[1])):
                        nextLine = next(file)
                        next_strip = nextLine.strip()
                        snext = next_strip.split(",")
                        data = [float(i) for i in snext[0:-1]]
                        Sect[last_idx_sec].setData(data)
                        line = nextLine
                elif strip_line.startswith("SECCONTROL"):
                    sspline = strip_line.split(",")
                    if len(sspline) > 2:
                        tmp = float(sspline[2].strip())
                        Sect[last_idx_sec].option = int(tmp)
                elif strip_line.startswith("ESYS"):
                    sspline = strip_line.split(",")
                elif "/TITLE" in line:
                    # Gestion du titre
                    title = line.split(",")[1].strip().replace("\n", "")

        # assert nb_total_nodes == len(self.nodes)
        assert nb_total_cells == len(Cells)
        # Réupération du nom du maillage
        self.mesh_name = title or osp.splitext(osp.split(filename)[-1])[0]

        toc = time.perf_counter()
        logger.debug(
            " File name : %s (parsed in %0.4f seconds)" % (filename, toc - tic)
        )
        logger.debug(
            " -> nodes: %d (parsed in %0.4f seconds)" % (len(self.nodes), time_nodes)
        )
        logger.debug(
            " -> cells: %d (parsed in %0.4f seconds)" % (nb_total_cells, time_cell)
        )
        logger.debug(
            " -> groups: %d (parsed in %0.4f seconds)"
            % (len(Groups) + len(CMELEM), time_groups)
        )

        logger.debug(" Mesh name : %s" % self.mesh_name)
        logger.debug(" Space Dimension : %d" % self.space_dim)
        logger.debug(" Number of nodes : %d" % (len(self.nodes)))

        # Les elements
        tic = time.perf_counter()
        e_conv = CellsTypeConverter("ANSYS")
        c_renum = ConnectivityRenumberer("ANSYS")
        for cell in Cells:
            if cell.type not in ElemEntities:
                ElemEntities[cell.type] = []
            ElemEntities[cell.type].append(cell.id)
            element_ansys_type = ElemAnsys[cell.type]
            # some trick for few cells (remove last node)
            element_ansys_test = str(element_ansys_type) + "_" + str(len(cell.nodes))
            if element_ansys_test in (
                "16_3",
                "18_3",
                "188_3",
                "189_4",
                "288_3",
                "289_4",
            ):
                nb_nodes = len(cell.nodes) - 1
                # logger.debug("Présence de noeuds orphelins")
            else:
                nb_nodes = len(cell.nodes)

            if "189" in element_ansys_test:
                nb_nodes = nb_nodes - 1
                logger.debug(
                    "Présence de BEAM189 : Passage d'une maille support SEG3 à SEG2"
                )
                del cell.nodes[2]

            elements_nodes_ansys = list(OrderedDict.fromkeys(cell.nodes[:nb_nodes]))

            if element_ansys_type == 200:
                element_ansys_type = "_".join(
                    map(
                        str,
                        (
                            element_ansys_type,
                            len(elements_nodes_ansys),
                            ElemOpt[cell.type][1],
                        ),
                    )
                )
            else:
                element_ansys_type = "_".join(
                    map(str, (element_ansys_type, len(elements_nodes_ansys)))
                )

            element_medcoupling_type = e_conv.external_to_medcoupling(
                element_ansys_type
            )
            element_nodes_med = c_renum.external_to_medcoupling(
                element_medcoupling_type, elements_nodes_ansys
            )

            # add group
            element_group = element_ansys_type.split("_")[0]
            if element_group in GROUPSMODELE:
                GROUPSMODELE[element_group].append(cell.id)
            else:
                GROUPSMODELE[element_group] = [cell.id]

            self.add_cell(cell.id, element_medcoupling_type, element_nodes_med)

        toc = time.perf_counter()
        logger.debug(" Load %d cells (in %0.4f seconds)" % (len(Cells), toc - tic))

        # Les groups
        tic = time.perf_counter()
        for group in Groups:
            values = []
            for elem in group.elems:
                if elem > 0:
                    values.append(elem)
                else:
                    values += range(values[-1] + 1, (-elem) + 1)

            if group.type == "NODE":
                self.add_group_nodes(group.name.strip(), values)
            elif group.type == "ELEM":
                self.add_group_cells(group.name.strip(), values)
            else:
                raise RuntimeError("Unknown group's type")

        # CMELEM
        for name, elem in CMELEM.items():
            values = []
            for ent in elem:
                values += ElemEntities[ent]

            self.add_group_cells(name.strip(), values)

        for group in GROUPSMODELE:
            rname = "Grp_FE_" + group.replace("-", "_").strip()
            self.add_group_cells(rname, GROUPSMODELE[group])

        toc = time.perf_counter()
        logger.debug(
            " Load %d groups (in %0.4f seconds)"
            % (len(Groups) + len(CMELEM), toc - tic)
        )

    def getCoor(self, line, firstStr, longFloat):
        # le premier decimal commence a la colonne firstStrg
        rline = line.rstrip()[firstStr:]
        elems = [
            float(rline[i : i + longFloat]) for i in range(0, len(rline), longFloat)
        ]

        nbElem = len(elems)

        if nbElem >= 3:
            return elems[0:3]
        else:
            return elems + [0.0] * (3 - nbElem)

    def __read_elem(self, file, elem):
        while True:
            line = file.readline()
            strip_line = line.strip()
            if strip_line.startswith("("):
                [nbElem, LongInt] = self.elem_format(strip_line)
            elif strip_line.startswith("-1"):
                break
            else:
                rline = line.rstrip()
                enum = [
                    int(rline[i : i + LongInt]) for i in range(0, len(rline), LongInt)
                ]
                assert len(enum) <= nbElem
                elem_id = enum[0]
                elem_type = enum[1]
                if elem_id in elem:
                    assert elem[elem_id] == elem_type
                else:
                    elem[elem_id] = elem_type
        return elem

    def __read_nodes(self, file, nodes):
        while True:
            line = file.readline()
            strip_line = line.strip()
            if strip_line.startswith("("):
                [firstStr, LongFloat] = self.node_format(strip_line)
            elif strip_line.startswith("N,") or strip_line.startswith("-1"):
                break
            else:
                spline = strip_line.split()
                self.add_node(int(spline[0]), self.getCoor(line, firstStr, LongFloat))
                nodes[int(spline[0])] = self.getCoor(line, firstStr, LongFloat)
        return nodes

    def __read_cells(self, file, Cells):
        l_new_cell = True
        while True:
            line = file.readline()
            rline = line.rstrip()
            strip_line = rline.lstrip()
            if strip_line.startswith("("):
                nbElem, LongInt = self.cell_format(strip_line)
            elif strip_line.startswith("-1"):
                break
            else:
                enum = [
                    int(rline[i : i + LongInt]) for i in range(0, len(rline), LongInt)
                ]
                assert len(enum) <= nbElem

                if l_new_cell:
                    cnodes = enum[11:]
                    nb_nodes = enum[8]
                    cid = enum[10]
                    ctype = enum[1]
                    sec = enum[3]
                    rep = enum[4]
                    const = enum[2]
                    if len(cnodes) < nb_nodes:
                        l_new_cell = False
                else:
                    cnodes += enum
                    if len(cnodes) == nb_nodes:
                        l_new_cell = True

                if l_new_cell:
                    assert len(cnodes) == nb_nodes
                    Cells.append(AnsysCell(ctype, cid, cnodes, sec, rep, const))

    def __read_groups(self, file, line, Groups):
        spline = line.split(",")
        gname = spline[1]
        gtype = spline[2]
        nb_elem = int(spline[3].split()[0])
        elems = []
        while True:
            line = file.readline()
            rline = line.rstrip()
            strip_line = rline.lstrip()
            if strip_line.startswith("("):
                nbElem, LongInt = self.cell_format(strip_line)
            else:
                elems += [
                    int(rline[i : i + LongInt]) for i in range(0, len(rline), LongInt)
                ]

                if len(elems) == nb_elem:
                    Groups.append(AnsysGroup(gname, gtype, elems))
                    break
                elif len(elems) > nb_elem:
                    raise RuntimeError("Wrong reading of groups")

    def decode_format(self, line):
        return line.strip().lstrip("(").rstrip(")").split(",")

    def elem_format(self, line):
        format = self.decode_format(line)
        s0 = format[0].split("i")
        s1 = format[1].split("a")
        nbElem = int(s0[0]) + int(s1[0])
        long = int(s0[1])
        assert long == int(s0[1])

        return [nbElem, long]

    def node_format(self, line):
        format = self.decode_format(line)
        s0 = format[0].split("i")
        firstStr = int(s0[0]) * int(s0[1])
        long = int(format[1].split("e")[1].split(".")[0])

        return [firstStr, long]

    def cell_format(self, line):
        format = self.decode_format(line)
        s0 = format[0].split("i")
        nbElem = int(s0[0])
        long = int(s0[1])

        return [nbElem, long]
