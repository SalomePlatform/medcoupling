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


dicoMod = {
    "5_2": "MMA-3D",
    "5_8": "TMA-3D",
    "11": "MBA-BARRE",
    "13_3_3": "MMA-AXIS",
    "13_3": "MMA-C_PLAN",
    "13_2_3": "TMA-PLAN",
    "13_2": "TMA-PLAN",
    "14": "MDI-DIS_T",
    "14_0": "MDI-DIS_T",
    "14_2": "MDD-2D_DIS_T",
    "16": "MPO-POU_D_T",
    "16_2": "NAN-undefined",
    "16_3": "NAN-undefined",
    "18": "MPO-POU_D_T",
    "18_2": "NAN-undefined",
    "18_3": "NAN-undefined",
    "21": "MDI-DIS_TR",
    "21_0": "MDI-DIS_TR",
    "21_2": "MDI-DIS_T",
    "21_3": "MDD-2D_DIS_TR",
    "21_4": "MDD-2D_DIS_T",
    "25": "MMA-AXIS",
    "25_3": "MMA-AXIS",
    "25_4": "MMA-D_PLAN",
    "27": "NAN-undefined",
    "29_3": "MMA-AXIS",
    "29": "MMA-D_PLAN_HHM",
    "30": "AMA-3D",
    "31": "TMA-AXIS",
    "33": "TMA-AXIS",
    "34": "TMA-AXIS",
    "35": "TMA-PLAN",
    "39": "NAN-undefined",
    "40": "NAN-undefined",
    "40_0": "NAN-undefined",
    "40_1": "NAN-undefined",
    "40_2": "NAN-undefined",
    "40_3": "NAN-undefined",
    "40_4": "NAN-undefined",
    "40_5": "NAN-undefined",
    "40_6": "NAN-undefined",
    "40_7": "NAN-undefined",
    "40_8": "NAN-undefined",
    "42": "MMA-C_PLAN",
    "42_0_3": "MMA-AXIS",
    "42_2_3": "MMA-AXIS",
    "42_0": "MMA-C_PLAN",
    "42_2": "MMA-D_PLAN",
    "45": "MMA-3D",
    "45_0": "MMA-3D",
    "45_1": "MMA-3D",
    "47": "TMA-PLAN",
    "47_1": "TMA-PLAN",
    "55": "TMA-PLAN",
    "59": "MPO-POU_D_T",
    "59_1": "MCA-CABLE",
    "59_0": "MDD-2D_POU_D_TR",
    "59_2": "MDD-2D_POU_D_TR",
    "61": "MCO-COQUE_AXIS",
    "63": "MCO-DKT",
    "65": "MMA-3D",
    "68": "TMA-AXIS",
    "70": "TMA-3D",
    "71": "NAN-undefined",
    "75": "TMA-PLAN",
    "77": "TMA-PLAN",
    "78": "TMA-PLAN",
    "82": "MMA-C_PLAN",
    "83": "MMA-C_PLAN",
    "87": "TMA-3D",
    "90": "TMA-3D",
    "92": "MMA-3D",
    "95": "MMA-3D",
    "95_0": "MMA-3D",
    "95_1": "MMA-3D",
    "96": "NAN-undefined",
    "98": "MMA-3D",
    "98_2": "MMA-3D",
    "98_8": "TMA-3D",
    "111": "NAN-undefined",
    "111_3": "TMA-3D",
    "116": "TMA-AXIS",
    "116_0": "TMA-AXIS",
    "116_1": "TMA-AXIS",
    "116_2": "NAN-undefined",
    "116_3": "NAN-undefined",
    "120": "NAN-undefined",
    "121": "NAN-undefined",
    "122": "NAN-undefined",
    "123": "NAN-undefined",
    "129": "NAN-undefined",
    "131": "TCO-COQUE",
    "132": "TCO-COQUE",
    "136": "MMA-D_PLAN_HM",
    "138": "NAN-undefined",
    "43": "MCO-Q4G",
    "43_0": "MCO-Q4G",
    "43_1": "MCO-Q4G",
    "143": "MCO-Q4G",
    "143_0": "MCO-Q4G",
    "143_1": "MCO-Q4G",
    "151": "TMA-AXIS",
    "152": "TMA-PLAN",
    "153": "MMA-AXIS",
    "154": "NAN-undefined",
    "156": "MMA-AXIS",
    "157": "TCO-COQUE",
    "160": "MBA-BARRE",
    "161": "MPO-POU_D_T",
    "162": "NAN-undefined",
    "163": "NAN-undefined",
    "164": "MMA-3D",
    "164_0": "NAN-undefined",
    "164_1": "NAN-undefined",
    "164_2": "NAN-undefined",
    "165": "MDI-DIS_T",
    "165_0": "MDI-DIS_T",
    "165_1": "MDI-DIS_TR",
    "166": "MDI-DIS_T",
    "167": "MCA-CABLE",
    "168": "MMA-3D",
    "169": "NAN-undefined",
    "169_3": "NAN-undefined",
    "170": "NAN-undefined",
    "171": "NAN-undefined",
    "171_0": "NAN-undefined",
    "171_1": "NAN-undefined",
    "171_2": "NAN-undefined",
    "171_7": "NAN-undefined",
    "171_8": "NAN-undefined",
    "171_9": "NAN-undefined",
    "171_10": "NAN-undefined",
    "172": "NAN-undefined",
    "172_0": "NAN-undefined",
    "172_1": "NAN-undefined",
    "172_2": "NAN-undefined",
    "172_7": "NAN-undefined",
    "172_8": "NAN-undefined",
    "172_9": "NAN-undefined",
    "172_10": "NAN-undefined",
    "173": "NAN-undefined",
    "173_0": "NAN-undefined",
    "173_1": "NAN-undefined",
    "173_2": "NAN-undefined",
    "173_8": "NAN-undefined",
    "173_9": "NAN-undefined",
    "173_10": "NAN-undefined",
    "174": "NAN-undefined",
    "174_0": "NAN-undefined",
    "174_1": "NAN-undefined",
    "174_2": "NAN-undefined",
    "174_8": "NAN-undefined",
    "174_9": "NAN-undefined",
    "174_10": "NAN-undefined",
    "175": "NAN-undefined",
    "175_0": "NAN-undefined",
    "175_1": "NAN-undefined",
    "175_2": "NAN-undefined",
    "175_8": "NAN-undefined",
    "175_9": "NAN-undefined",
    "175_10": "NAN-undefined",
    "176": "NAN-undefined",
    "177": "NAN-undefined",
    "180_0": "MBA-BARRE",
    "180_1": "MCA-CABLE",
    "181": "MCO-Q4G",
    "181_0": "MCO-Q4G",
    "181_1": "MCO-Q4G",
    "182": "MMA-C_PLAN",
    "182_0": "MMA-AXIS",
    "182_2": "MMA-D_PLAN",
    "182_3": "MMA-AXIS",
    "182_5": "MMA-C_PLAN",
    "183": "MMA-C_PLAN",
    "183_0": "MMA-C_PLAN",
    "183_2": "MMA-C_PLAN",
    "183_3": "MMA-C_PLAN",
    "183_5": "MMA-C_PLAN",
    "184": "MPO-POU_D_T",
    "185": "MMA-3D",
    "186": "MMA-3D_SI",
    "186_0": "MMA-3D_SI",
    "186_1": "MMA-3D",
    "187": "MMA-3D",
    "4": "MPO-POU_D_T",
    "4_3": "MPO-POU_D_T",
    "44": "MPO-POU_D_T",
    "44_3": "MPO-POU_D_T",
    "188": "MPO-POU_D_T",
    "189": "MPO-POU_D_T",
    "189_3": "MPO-POU_D_T",
    "190": "MCO-DKT",
    "192": "M-PLAN_JOINT",
    "195": "M-3D_JOINT",
    "200": "NAN-undefined",
    "202": "MMA-C_PLAN",
    "202_0": "MMA-C_PLAN",
    "202_2": "MMA-D_PLAN",
    "202_3": "MMA-C_PLAN",
    "204": "MMA-3D",
    "212_3": "NAN-undefined",
    "212": "NAN-undefined",
    "213_3": "M-D_PLAN_HM",
    "213": "M-D_PLAN_HHM",
    "215": "NAN-undefined",
    "216": "M-3D_HM",
    "217": "NAN-undefined",
    "218": "M-D_PLAN_THH",
    "218_0": "NAN-undefined",
    "218_1": "M-D_PLAN_THH",
    "220": "AMA-3D",
    "221": "AMA-3D",
    "223": "NAN-undefined",
    "223_11": "NAN-undefined",
    "223_100001": "NAN-undefined",
    "223_100010": "NAN-undefined",
    "223_100011": "NAN-undefined",
    "226": "NAN-undefined",
    "226_11": "NAN-undefined",
    "226_100001": "NAN-undefined",
    "226_100010": "NAN-undefined",
    "226_100011": "NAN-undefined",
    "227": "NAN-undefined",
    "227_11": "NAN-undefined",
    "227_100001": "NAN-undefined",
    "227_100010": "NAN-undefined",
    "227_100011": "NAN-undefined",
    "230": "NAN-undefined",
    "231": "NAN-undefined",
    "232": "NAN-undefined",
    "233": "NAN-undefined",
    "236": "NAN-undefined",
    "237": "NAN-undefined",
    "238": "NAN-undefined",
    "239": "NAN-undefined",
    "240": "NAN-undefined",
    "251": "NAN-undefined",
    "252": "NAN-undefined",
    "278": "NAN-undefined",
    "279": "NAN-undefined",
    "281": "MCO-COQUE_3D",
    "281_0": "MCO-COQUE_3D",
    "281_1": "MCO-COQUE_3D",
    "285": "NAN-undefined",
    "288": "MPO-POU_D_T",
    "289_3": "MPO-TUYAU_3M",
    "290": "MPO-TUYAU_3M",
}

dicoOpt = {
    "4": None,
    "5": 1,
    "11": None,
    "13": 1,
    "14": 3,
    "16": None,
    "18": None,
    "21": 3,
    "25": None,
    "27": None,
    "29": None,
    "30": None,
    "31": None,
    "33": None,
    "34": None,
    "35": None,
    "39": 4,
    "40": 3,
    "42": 3,
    "44": None,
    "45": 2,
    "47": 1,
    "55": None,
    "59": 1,
    "61": None,
    "63": None,
    "65": None,
    "68": None,
    "70": None,
    "71": None,
    "75": None,
    "77": None,
    "78": None,
    "82": None,
    "83": None,
    "87": None,
    "90": None,
    "92": None,
    "95": 11,
    "96": None,
    "98": 1,
    "111": 1,
    "116": 1,
    "120": None,
    "121": None,
    "122": None,
    "123": None,
    "129": None,
    "131": None,
    "132": None,
    "136": None,
    "138": None,
    "143": 1,
    "151": None,
    "152": None,
    "153": None,
    "154": None,
    "156": None,
    "157": None,
    "160": None,
    "161": None,
    "162": None,
    "163": None,
    "164": 1,
    "165": 1,
    "166": None,
    "167": None,
    "168": None,
    "169": None,
    "170": None,
    "171": 1,
    "172": 1,
    "173": 1,
    "174": 1,
    "175": 1,
    "176": None,
    "177": None,
    "180": None,
    "43": 1,
    "181": 1,
    "182": 3,
    "183": 3,
    "184": None,
    "185": None,
    "186": None,
    "187": None,
    "188": None,
    "189": None,
    "190": None,
    "192": None,
    "195": None,
    "200": None,
    "202": 3,
    "204": None,
    "212": None,
    "213": None,
    "215": None,
    "216": None,
    "217": None,
    "218": 1,
    "220": None,
    "221": None,
    "223": 1,
    "226": 1,
    "227": 1,
    "230": None,
    "231": None,
    "232": None,
    "233": None,
    "236": None,
    "237": None,
    "238": None,
    "239": None,
    "240": None,
    "251": None,
    "252": None,
    "278": None,
    "279": None,
    "281": 1,
    "285": None,
    "288": None,
    "289": None,
    "290": None,
}


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

            # add cell in medcoupling format
            if (
                cell.type in ElemOpt
                and ElemOpt[cell.type][0] == dicoOpt[str(element_ansys_type)]
            ):
                element_group_type = (
                    str(element_ansys_type) + "_" + str(ElemOpt[cell.type][1])
                )
            else:
                element_group_type = str(element_ansys_type)
            if nb_nodes == 3:
                element_group_type = element_group_type + "_" + str(nb_nodes)
            if element_group_type == "180":
                element_group_type = (
                    element_group_type + "_" + str(Sect[cell.sec].option)
                )
            try:
                element_group = dicoMod[element_group_type]
            except (ValueError, TypeError):
                msg = "Erreur: Option de l'élément non traitée"
                raise RuntimeError(msg)

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
                self.add_group_nodes(group.name, values)
            elif group.type == "ELEM":
                self.add_group_cells(group.name, values)
            else:
                raise RuntimeError("Unknown group's type")

        # CMELEM
        for name, elem in CMELEM.items():
            values = []
            for ent in elem:
                values += ElemEntities[ent]

            self.add_group_cells(name, values)

        for group in GROUPSMODELE:
            rname = group.replace("-", "_")
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
