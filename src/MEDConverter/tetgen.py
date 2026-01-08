#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import os.path as osp

from .logger import logger
from .MEDConverterMesh import MEDConverterMesh
from .cells import CellsTypeConverter
from .connectivity import ConnectivityRenumberer


class MEDConverterTetgen(MEDConverterMesh):
    @staticmethod
    def convert_tetgen_to_med(filename_tetgen, verbose=False):

        tic = time.perf_counter()
        c = MEDConverterTetgen()
        c.verbose = verbose
        c.read_tetgen_mesh(filename_tetgen)
        c.create_UMesh()
        toc = time.perf_counter()
        logger.debug("Mesh converted (in %0.4f seconds)" % (toc - tic))
        return c.umesh

    def __init__(self):
        super(MEDConverterTetgen, self).__init__()

    def read_tetgen_mesh(self, filename):
        logger.debug("Read TETGEN mesh.")

        self._reset_structures()

        GROUPS_M = {}

        with open(filename, "r", encoding=self._get_file_encoding(filename)) as file:
            self.filename = filename
            self.mesh_name = osp.splitext(osp.split(filename)[-1])[0]
            self.space_dim = 3

            logger.debug("Mesh name : %s" % self.mesh_name)
            logger.debug("Space Dimension : 3")
            logger.debug("Beginning to parse mesh file")

            # read nodes
            tic = time.perf_counter()
            nb_nodes = int(file.readline())
            for i_node in range(nb_nodes):
                idx_tetgen = i_node + 1
                spline = file.readline().split()
                coords = tuple(map(float, spline))
                self.add_node(idx_tetgen, coords)
            toc = time.perf_counter()
            logger.debug(" Load %d nodes (in %0.4f seconds)" % (nb_nodes, toc - tic))

            # Les elements
            e_conv = CellsTypeConverter("TETGEN")
            c_renum = ConnectivityRenumberer("TETGEN")

            # read tetra
            tic = time.perf_counter()
            idx_tetgen = 0
            nb_tetra = int(file.readline())
            for i_tet in range(nb_tetra):
                idx_tetgen += 1
                spline = file.readline().split()

                nodes_tet = tuple(map(int, spline[1:]))
                tet_medcoupling_type = e_conv.external_to_medcoupling("TETRA4")
                nodes_med = c_renum.external_to_medcoupling(
                    tet_medcoupling_type, nodes_tet
                )
                self.add_cell(idx_tetgen, tet_medcoupling_type, nodes_med)

                group_id = "C" + spline[0]
                if group_id in GROUPS_M:
                    GROUPS_M[group_id].append(idx_tetgen)
                else:
                    GROUPS_M[group_id] = [idx_tetgen]

            toc = time.perf_counter()
            logger.debug(
                " Load %d tetrahedra (in %0.4f seconds)" % (nb_tetra, toc - tic)
            )

            # read triangles
            tic = time.perf_counter()
            nb_tri = int(file.readline())
            for i_tri in range(nb_tri):
                idx_tetgen += 1
                spline = file.readline().split()

                nodes_tet = tuple(map(int, spline[1:]))
                tet_medcoupling_type = e_conv.external_to_medcoupling("TRIA3")
                nodes_med = c_renum.external_to_medcoupling(
                    tet_medcoupling_type, nodes_tet
                )
                self.add_cell(idx_tetgen, tet_medcoupling_type, nodes_med)

                group_id = "F" + spline[0]
                if group_id in GROUPS_M:
                    GROUPS_M[group_id].append(idx_tetgen)
                else:
                    GROUPS_M[group_id] = [idx_tetgen]

            toc = time.perf_counter()
            logger.debug(" Load %d triangles (in %0.4f seconds)" % (nb_tri, toc - tic))

            tic = time.perf_counter()
            nb_groups = 0
            for group_name, values in GROUPS_M.items():
                self.add_group_cells(group_name, values)
                nb_groups += 1

            toc = time.perf_counter()
            logger.debug(" Load %d groups (in %0.4f seconds)" % (nb_groups, toc - tic))
