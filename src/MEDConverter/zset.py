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

from .logger import logger
from .MEDConverterMesh import MEDConverterMesh
from .cells import CellsTypeConverter, GroupCellsTypeConverter
from .connectivity import ConnectivityRenumberer

ZSET_MAX_LINE_SIZE = 21
ZSET_NODES_SHIFT = 1  # La numérotation ZSET des noeuds démarre à 1
ZSET_CELLS_SHIFT = 1  # La numérotation ZSET des élements démarre à 1


class MEDConverterZset(MEDConverterMesh):
    @staticmethod
    def convert_zset_to_med(filename_zset, verbose=False):

        tic = time.perf_counter()
        c = MEDConverterZset()
        c.verbose = verbose
        c.read_zset_mesh(filename_zset)
        c.create_UMesh()
        toc = time.perf_counter()
        logger.debug("Mesh converted (in %0.4f seconds)" % (toc - tic))
        return c.umesh

    def __init__(self):
        super().__init__()
        self.zsetmesh = None

    def read_zset_mesh(self, filename):
        logger.debug("Read ZSET mesh.")

        self._reset_structures()

        NODES, ELEMENTS, GROUPS = [], [], []

        flag = {"NODES": 0, "ELEMENTS": 0, "GROUPS": 0}

        tic = time.perf_counter()
        with open(filename, "r", encoding=self._get_file_encoding(filename)) as f:
            self.mesh_name = osp.splitext(osp.split(filename)[-1])[0]

            for line in f:
                if flag["NODES"] == 1:
                    NODES.append(line)
                elif flag["ELEMENTS"] == 1:
                    ELEMENTS.append(line)
                elif flag["GROUPS"] == 1:
                    GROUPS.append(line)

                if "**node" in line:
                    flag["NODES"] = 1

                elif "**element" in line:
                    flag["ELEMENTS"] = 1
                    flag["NODES"] = 0

                elif "***group" in line:
                    flag["ELEMENTS"] = 0
                    flag["GROUPS"] = 1

                elif "***return" in line:
                    flag["GROUPS"] = 0

        nb_nodes, self.space_dim = map(int, NODES[0].split())
        nb_elements = int(ELEMENTS[0].split()[0])
        toc = time.perf_counter()

        logger.debug(
            " File name : %s (parsed in %0.4f seconds)" % (filename, toc - tic)
        )
        logger.debug(" Mesh name : %s" % self.mesh_name)
        logger.debug(" Space Dimension : %d" % self.space_dim)

        # Les noeuds
        tic = time.perf_counter()
        for line in NODES[1:-1]:
            spline = line.split()
            idx_zset = int(spline[0])
            coords = tuple(map(float, spline[-self.space_dim :]))
            self.add_node(idx_zset, coords)
        toc = time.perf_counter()
        logger.debug(" Load %d nodes (in %0.4f seconds)" % (len(NODES) - 2, toc - tic))

        # Les elements
        e_conv = CellsTypeConverter("ZSET")
        g_conv = GroupCellsTypeConverter("ZSET")
        c_renum = ConnectivityRenumberer("ZSET")

        tic = time.perf_counter()
        for line in ELEMENTS[1:-1]:
            spline = line.split()
            idx_element_zset = int(spline[0])
            element_zset_type = spline[1]
            elements_nodes_zset = tuple(map(int, spline[2:]))

            element_medcoupling_type = e_conv.external_to_medcoupling(element_zset_type)
            element_nodes_med = c_renum.external_to_medcoupling(
                element_medcoupling_type, elements_nodes_zset
            )

            self.add_cell(idx_element_zset, element_medcoupling_type, element_nodes_med)

        toc = time.perf_counter()
        logger.debug(
            " Load %d cells (in %0.4f seconds)" % (len(ELEMENTS) - 2, toc - tic)
        )

        # Les groupes
        tic = time.perf_counter()
        groups = {}
        for line in GROUPS[:-1]:
            if "**" in line:
                type_grp, name_grp = line.strip("**").split()
                if type_grp not in groups:
                    groups[type_grp] = {}
                if name_grp not in groups[type_grp]:
                    groups[type_grp][name_grp] = []
            else:
                groups[type_grp][name_grp].append(line.split())

        toc = time.perf_counter()
        logger.debug(" Parse groups (in %0.4f seconds)" % (toc - tic))

        tic = time.perf_counter()
        nodes_groups = groups.get("nset", {})
        for group_name, items in sorted(nodes_groups.items()):
            values = (int(i) for line in items for i in line)
            self.add_group_nodes(group_name, values)
        toc = time.perf_counter()
        logger.debug(" Add %d nset (in %0.4f seconds)" % (len(nodes_groups), toc - tic))

        tic = time.perf_counter()
        cells_groups = groups.get("elset", {})
        for group_name, items in sorted(cells_groups.items()):
            values = (int(i) for line in items for i in line)
            self.add_group_cells(group_name, values)
        toc = time.perf_counter()
        logger.debug(
            " Add %d elset (in %0.4f seconds)" % (len(cells_groups), toc - tic)
        )

        tic = time.perf_counter()
        faset = groups.get("faset", {})
        for faset_name, items in sorted(faset.items()):
            self._add_bset(faset_name, items, g_conv, c_renum)
        toc = time.perf_counter()
        logger.debug(" Add %d faset (in %0.4f seconds)" % (len(faset), toc - tic))

        tic = time.perf_counter()
        liset = groups.get("liset", {})
        for liset_name, items in sorted(liset.items()):
            self._add_bset(liset_name, items, g_conv, c_renum)
        toc = time.perf_counter()
        logger.debug(" Add %d liset (in %0.4f seconds)" % (len(liset), toc - tic))

    def _add_bset(self, bset_name, bset_items, g_conv, c_renum):
        max_idx_elements = max(
            i for dim in self.corresponding_cells.values() for i in dim
        )
        values = []
        for i, spline in enumerate((j for j in bset_items if bool(j))):
            idx_element_zset = max_idx_elements + i + 1
            element_zset_type = spline[0]
            element_nodes_zset = tuple(map(int, spline[1:]))

            element_medcoupling_type = g_conv.external_to_medcoupling(element_zset_type)
            element_nodes_med = c_renum.external_to_medcoupling(
                element_medcoupling_type, element_nodes_zset
            )

            self.add_cell(idx_element_zset, element_medcoupling_type, element_nodes_med)
            values.append(idx_element_zset)
        self.add_group_cells(bset_name, values)
