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
from .cells import CellsTypeConverter
from .connectivity import ConnectivityRenumberer

SYSTUS_NODES_SHIFT = 1  # La numérotation SYSTUS des noeuds démarre à 1
SYSTUS_CELLS_SHIFT = 1  # La numérotation SYSTUS des élements démarre à 1


class MEDConverterSystus(MEDConverterMesh):
    @staticmethod
    def convert_systus_to_med(filename_systus, skip_types, verbose=False):

        tic = time.perf_counter()
        c = MEDConverterSystus()
        c.verbose = verbose
        c.read_systus_mesh(filename_systus, skip_types)
        c.create_UMesh()
        toc = time.perf_counter()
        logger.debug("Mesh converted (in %0.4f seconds)" % (toc - tic))
        return c.umesh

    def __init__(self):
        super().__init__()
        self.systusmesh = None

    def read_systus_mesh(self, filename, skip_types=[]):
        logger.debug("Read SYSTUS mesh.")

        self._reset_structures()
        NODES, ELEMENTS, GROUPS = [], [], []

        flag = {"NODES": 0, "ELEMENTS": 0, "GROUPS": 0}

        tic = time.perf_counter()
        # Lecture du fichier .ASC où les blocs sont separés par des BEGIN_* et END_*
        with open(filename, "r", encoding=self._get_file_encoding(filename)) as f:
            type_systus, num_systus = next(f).split()[:2]
            if not (type_systus in ("1VSD",) and num_systus in ("0",)):
                msg = "Cannot convert SYSTUS mesh tagged '{}' '{}'. Supported SYSTUS mesh must be tagged '1VSD' '0'.".format(
                    type_systus, num_systus
                )
                raise RuntimeError(msg)

            # Lecture du nom du maillage si disponible
            line_1 = next(f).strip()
            self.mesh_name = line_1 or osp.splitext(osp.split(filename)[-1])[0]

            for line in f:

                if flag["NODES"] == 1:
                    NODES.append(line)
                elif flag["ELEMENTS"] == 1:
                    ELEMENTS.append(line)
                elif flag["GROUPS"] == 1:
                    GROUPS.append(line)

                if "BEGIN_NODES" in line:
                    flag["NODES"] = 1
                    self.space_dim = int(line.split()[2])
                elif "END_NODES" in line:
                    flag["NODES"] = 0

                elif "BEGIN_ELEMENTS" in line:
                    flag["ELEMENTS"] = 1
                elif "END_ELEMENTS" in line:
                    flag["ELEMENTS"] = 0

                elif "BEGIN_GROUPS" in line:
                    flag["GROUPS"] = 1
                elif "END_GROUPS" in line:
                    flag["GROUPS"] = 0

        toc = time.perf_counter()
        logger.debug(
            " File name : %s (parsed in %0.4f seconds)" % (filename, toc - tic)
        )
        logger.debug(" Mesh name : %s" % self.mesh_name)
        logger.debug(" Space Dimension : %d" % self.space_dim)

        # Les noeuds
        tic = time.perf_counter()
        for line in NODES[:-1]:
            spline = line.split()
            idx_systus = int(spline[0])
            coords = tuple(map(float, spline[-self.space_dim :]))
            self.add_node(idx_systus, coords)
        toc = time.perf_counter()
        logger.debug(" Load %d nodes (in %0.4f seconds)" % (len(NODES) - 1, toc - tic))

        # Les elements
        e_conv = CellsTypeConverter("SYSTUS")
        c_renum = ConnectivityRenumberer("SYSTUS")
        skipped = 0

        tic = time.perf_counter()
        for line in ELEMENTS[:-1]:
            spline = line.split()
            idx_element_systus = int(spline[0])
            element_systus_type = "%04d" % int(spline[1])

            # Skip some type of cells:
            if element_systus_type in skip_types:
                skipped += 1
                continue

            elements_nodes_systus = tuple(map(int, spline[5:]))

            element_medcoupling_type = e_conv.external_to_medcoupling(
                element_systus_type
            )
            element_nodes_med = c_renum.external_to_medcoupling(
                element_medcoupling_type, elements_nodes_systus
            )

            self.add_cell(
                idx_element_systus, element_medcoupling_type, element_nodes_med
            )

        toc = time.perf_counter()

        logger.debug(
            " Load %d cells (in %0.4f seconds)"
            % (len(ELEMENTS) - 1 - skipped, toc - tic)
        )
        if skipped:
            logger.debug(" Skip %d cells" % skipped)

        tic = time.perf_counter()
        # Les groups
        for line in GROUPS[:-1]:
            spline = line.split()
            values = map(int, line.split('"')[-1].split())
            group_name = spline[1]
            group_tag_systus = spline[2]

            if group_tag_systus == "1":
                self.add_group_nodes(group_name, values)
            else:
                self.add_group_cells(group_name, values)
        toc = time.perf_counter()
        logger.debug(
            " Load %d groups (in %0.4f seconds)" % (len(GROUPS) - 1, toc - tic)
        )

    def write_systus_mesh(self, filename):
        tic = time.perf_counter()
        with open(filename, "w") as f:
            f.write(self.systusmesh)
        toc = time.perf_counter()
        logger.debug(
            "Write SYSTUS mesh file : %s (in %0.4f seconds)" % (filename, toc - tic)
        )

    def create_systus_mesh(self):
        self.systusmesh = None
        logger.debug("Create SYSTUS mesh.")

        if not (bool(self.cells_continuous) or bool(self.groups_e_continuous)):
            tic = time.perf_counter()
            self._make_continuous()
            toc = time.perf_counter()
            logger.debug("Make continuous (in %0.4f seconds)" % (toc - tic))

        logger.debug(" Mesh name : %s" % self.mesh_name)
        logger.debug(" Space Dimension : %d" % self.space_dim)

        # Noeuds
        tic = time.perf_counter()
        nb_nodes = len(self.nodes)
        nodes_lines = (
            "%d 0 0 0 0 0 " % (i + SYSTUS_NODES_SHIFT) + " ".join(map(repr, node))
            for i, node in enumerate(self.nodes)
        )
        toc = time.perf_counter()
        logger.debug(" Add %d nodes (in %0.4f seconds)" % (nb_nodes, toc - tic))

        # Elements et groupes
        elements_lines = []
        groups_lines = []
        groups_e_ids = {}
        groups_n_ids = {}

        c_renum = ConnectivityRenumberer("SYSTUS")
        e_conv = CellsTypeConverter("SYSTUS")

        tic = time.perf_counter()
        for j, (medcoupling_type, element_nodes_med) in self.cells_continuous.items():
            systus_type = e_conv.medcoupling_to_external(medcoupling_type)
            element_nodes_med = [i + SYSTUS_NODES_SHIFT for i in element_nodes_med]
            element_nodes_asc = c_renum.medcoupling_to_external(
                medcoupling_type, element_nodes_med
            )
            elements_lines.append(
                "%d %s 0 0 0 " % (j + SYSTUS_CELLS_SHIFT, systus_type)
                + " ".join(map(str, element_nodes_asc))
            )

        nb_elements = len(self.cells_continuous)
        toc = time.perf_counter()
        logger.debug(" Add %d cells (in %0.4f seconds)" % (nb_elements, toc - tic))

        tic = time.perf_counter()
        for group, values in self.groups_e_continuous.items():
            groups_e_ids[group] = [i + SYSTUS_CELLS_SHIFT for i in values]

        for group, values in self.groups_n.items():
            groups_n_ids[group] = [i + SYSTUS_NODES_SHIFT for i in values]

        id_groups = 1  # La numérotation des groupes systus est incrementale et commune à tout type de groupe
        for name in sorted(groups_e_ids.keys()):
            group_e = groups_e_ids[name]
            group_line = '%d %s 2 0 "COLLECTOR_ID %d"  ""  "" %s' % (
                id_groups,
                name,
                id_groups,
                " ".join(map(str, group_e)),
            )
            id_groups += 1
            groups_lines.append(group_line)

        for name in sorted(groups_n_ids.keys()):
            group_n = groups_n_ids[name]
            group_line = '%d %s 1 0 "COLLECTOR_ID %d"  ""  "" %s' % (
                id_groups,
                name,
                id_groups,
                " ".join(map(str, group_n)),
            )
            id_groups += 1
            groups_lines.append(group_line)

        nb_groups = id_groups - 1
        toc = time.perf_counter()
        logger.debug(" Add %d groups (in %0.4f seconds)" % (nb_groups, toc - tic))

        tic = time.perf_counter()
        # Entete du fichier
        txt_header = """1VSD 0 {0} {0}
{1}
 100000 4 {2} {3} 0 3 6 0 0
BEGIN_INFORMATIONS
{1}
 4 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 3 3 9 0 0 0 0 9 0 0 0 0 6 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0
END_INFORMATIONS
""".format(
            *[time.strftime("%y%m%d %H%M%S"), self.mesh_name, nb_nodes, nb_elements]
        )

        txt_nodes = "BEGIN_NODES %d %d\n%s\nEND_NODES\n" % (
            nb_nodes,
            self.space_dim,
            "\n".join(nodes_lines),
        )
        txt_elements = "BEGIN_ELEMENTS %d\n%s\nEND_ELEMENTS\n" % (
            nb_elements,
            "\n".join(elements_lines),
        )
        txt_groups = (
            "BEGIN_GROUPS %d\n%s\nEND_GROUPS\n" % (nb_groups, "\n".join(groups_lines))
            if groups_lines
            else ""
        )
        self.systusmesh = "".join((txt_header, txt_nodes, txt_elements, txt_groups))
        toc = time.perf_counter()
        logger.debug(" Assembly file (in %0.4f seconds)" % (toc - tic))
