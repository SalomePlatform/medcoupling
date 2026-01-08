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
import re
import os.path as osp

from .utilities import chunks
from .logger import logger
from .MEDConverterMesh import MEDConverterMesh
from .cells import CellsTypeConverter
from .connectivity import ConnectivityRenumberer

ASTER_MAX_LINE_SIZE = 8
ASTER_NODES_SHIFT = 1  # La numérotation ASTER des noeuds démarre à 1
ASTER_CELLS_SHIFT = 1  # La numérotation ASTER des élements démarre à 1


class MEDConverterAster(MEDConverterMesh):
    @staticmethod
    def convert_aster_to_med(filename_aster, verbose=False):

        tic = time.perf_counter()
        c = MEDConverterAster()
        c.verbose = verbose
        c.read_aster_mesh(filename_aster)
        c.create_UMesh()
        toc = time.perf_counter()
        logger.debug("Mesh converted (in %0.4f seconds)" % (toc - tic))
        return c.umesh

    def __init__(self):
        super(MEDConverterAster, self).__init__()
        self.astermesh = None

    def read_aster_mesh(self, filename):
        logger.debug("Read ASTER mesh.")

        self._reset_structures()

        def zip_line(line):
            # Restituer la ligne sans champs et les champs à part
            spline = line.split()
            fields = dict(
                i.split("=") for i in spline if ("=" in i and len(i.split("=")) == 2)
            )
            zipped = " ".join((i for i in spline if "=" not in i)).strip()
            return zipped, fields

        def zip_block(block):
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
                zipped, flds = zip_line(item)
                if len(zipped) > 0:
                    zipped_lines.append(zipped)
                fields.update(flds)

            block_name = zipped_lines[0]
            block_body = " ".join(zipped_lines[1:]).split()
            return block_name, block_body, fields

        def remove_comments_and_split(fstream):
            comment = "%"
            txt = "\n".join(line.partition(comment)[0].strip() for line in fstream)
            return txt.split("FINSF")

        tic = time.perf_counter()

        with open(filename, "r", encoding=self._get_file_encoding(filename)) as f:
            mesh_blocks = remove_comments_and_split(f)
        toc = time.perf_counter()
        logger.debug(
            " File name : %s (splitted in %0.4f seconds)" % (filename, toc - tic)
        )

        tic = time.perf_counter()
        NODES, ELEMENTS, GROUPS_N, GROUPS_M = [], {}, {}, {}

        self.mesh_name = osp.splitext(osp.split(filename)[-1])[0]
        for block in mesh_blocks:
            bname, bbody, bfields = zip_block(block)

            if bname in ("COOR_3D", "COOR_2D"):
                self.space_dim = int(bname.strip("COOR_").strip("D"))
                NODES = list(chunks(bbody, self.space_dim + 1))

            elif bname in ("GROUP_MA",):
                if "NOM" in bfields:
                    grp_name = bfields["NOM"]
                    GROUPS_M[grp_name] = bbody
                else:
                    for k, v in bfields.items():
                        raise ValueError(f"Group use an unknown keyword: {k}={v}.")
                    grp_name = bbody[0]
                    GROUPS_M[grp_name] = bbody[1:]

            elif bname in ("GROUP_NO",):
                if "NOM" in bfields:
                    grp_name = bfields["NOM"]
                    GROUPS_N[grp_name] = bbody
                else:
                    for k, v in bfields.items():
                        raise ValueError(f"Group use an unknown keyword: {k}={v}.")
                    grp_name = bbody[0]
                    GROUPS_N[grp_name] = bbody[1:]

            elif bname in CellsTypeConverter._aster_to_med.keys():
                etype = bname
                nb_nodes = int(re.findall(r"\d+", etype)[0])
                ELEMENTS.setdefault(etype, [])
                ELEMENTS[etype].extend(list(chunks(bbody, 1 + nb_nodes)))

            elif bname in ("TITRE",):
                parts = []
                for k, v in bfields.items():
                    parts.append(f"{k} = {v}")
                parts.extend(bbody)
                new_name = " ".join(parts).strip()
                if len(new_name) > 0:
                    self.mesh_name = new_name
            else:
                pass

        toc = time.perf_counter()
        logger.debug(
            " Mesh name : %s (parsed in %0.4f seconds)" % (self.mesh_name, toc - tic)
        )
        logger.debug(" Space Dimension : %d" % self.space_dim)

        # Les noeuds
        strip_fortran_notation = lambda s: s.replace("d", "e").replace("D", "E")
        tic = time.perf_counter()
        for spline in NODES:
            idx_aster = spline[0]
            coords = tuple(
                float(strip_fortran_notation(c)) for c in spline[-self.space_dim :]
            )
            self.add_node(idx_aster, coords)
        toc = time.perf_counter()
        logger.debug(" Load %d nodes (in %0.4f seconds)" % (len(NODES), toc - tic))

        # Les elements
        e_conv = CellsTypeConverter("ASTER")
        c_renum = ConnectivityRenumberer("ASTER")

        tic = time.perf_counter()
        nb_elements = 0
        for element_aster_type, cells in ELEMENTS.items():
            for spline in cells:
                idx_element_aster = spline[0]
                elements_nodes_aster = spline[1:]

                element_medcoupling_type = e_conv.external_to_medcoupling(
                    element_aster_type
                )
                element_nodes_med = c_renum.external_to_medcoupling(
                    element_medcoupling_type, elements_nodes_aster
                )

                self.add_cell(
                    idx_element_aster, element_medcoupling_type, element_nodes_med
                )
                nb_elements += 1

        toc = time.perf_counter()
        logger.debug(" Load %d cells (in %0.4f seconds)" % (nb_elements, toc - tic))

        # Les groups
        tic = time.perf_counter()
        nb_groups = 0
        for group_name, values in GROUPS_N.items():
            self.add_group_nodes(group_name, values)
            nb_groups += 1
        toc = time.perf_counter()
        logger.debug(
            " Load %d groups of nodes (in %0.4f seconds)" % (nb_groups, toc - tic)
        )

        tic = time.perf_counter()
        nb_groups = 0
        for group_name, values in GROUPS_M.items():
            self.add_group_cells(group_name, values)
            nb_groups += 1
        toc = time.perf_counter()
        logger.debug(
            " Load %d groups of cells (in %0.4f seconds)" % (nb_groups, toc - tic)
        )
