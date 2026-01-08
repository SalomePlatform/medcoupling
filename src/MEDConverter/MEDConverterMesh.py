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
import logging
from collections import OrderedDict

from .logger import logger


class MEDConverterMesh:
    @classmethod
    def from_mesh(cls, other):
        obj = cls()

        attr_to_import = (
            "mesh_name",
            "space_dim",
            "nodes",
            "cells",
            "groups_e",
            "groups_n",
            "verbose",
            "medmesh",
            "_corresponding_nodes",
            "_corresponding_cells",
            "groups_e_continuous",
            "cells_continuous",
        )

        for attr in attr_to_import:
            setattr(obj, attr, getattr(other, attr))

        return obj

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, verbose):
        self._verbose = verbose
        if verbose:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)

    @property
    def mesh_name(self):
        return self._mesh_name

    @mesh_name.setter
    def mesh_name(self, name):
        not_allowed_symbols_in_name = ("/",)
        if any(symbol in name for symbol in not_allowed_symbols_in_name):
            msg = f"The following symbols are not allowed in MED mesh name: {not_allowed_symbols_in_name}. It is replaced by '-'."
            logger.info(msg)
            for symbol in not_allowed_symbols_in_name:
                name = name.replace(symbol, "-")

        MED_NAME_SIZE = 64
        if len(name) > MED_NAME_SIZE:
            msg = f"Mesh name '{name}' is too long {len(name)}>{MED_NAME_SIZE}. Only {MED_NAME_SIZE} first characters are used."
            logger.info(msg)

        self._mesh_name = name[: MED_NAME_SIZE - 1].strip().strip('"').strip("'")

    @property
    def dimensions(self):
        return sorted(self.cells.keys())[::-1]

    @property
    def max_dim_cells(self):
        return max(self.cells.keys())

    @property
    def levels(self):
        mdim = int(self.max_dim_cells[0])
        return {"%dD" % i: i - mdim for i in range(mdim, -1, -1)}

    @property
    def corresponding_cells(self):
        return self._corresponding_cells

    @property
    def corresponding_cells_reversed(self):
        return {
            dim: {item: key for key, item in values.items()}
            for dim, values in self._corresponding_cells.items()
        }

    @property
    def corresponding_nodes(self):
        return self._corresponding_nodes

    @property
    def corresponding_nodes_reversed(self):
        return {item: key for key, item in self._corresponding_nodes.items()}

    def _reset_structures(self):
        self._mesh_name = ""
        self.space_dim = 0
        self.nodes = []
        self.cells = OrderedDict()  # Par niveau
        self.groups_e = OrderedDict()  # Par niveau
        self.groups_n = OrderedDict()

        self.groups_e_continuous = OrderedDict()  # Numérotation globale
        self.cells_continuous = OrderedDict()  # Numérotation globale

        self.umesh = None

        self._corresponding_nodes = {}
        self._corresponding_cells = {}

    def __init__(self):
        self._reset_structures()
        self._verbose = False

    def _get_file_encoding(self, filename):

        encodings = "utf8 latin_1 cp437".split()

        for enc in encodings:
            try:
                with open(filename, mode="r", encoding=enc) as f:
                    f.read()
                return enc
            except UnicodeDecodeError:
                continue

        msg = "File encoding is not among : %s" % (", ".join(encodings))
        raise RuntimeError(msg)

    def _check_group_name(self, name):
        MED_LNAME_SIZE = 80
        if len(name) > MED_LNAME_SIZE:
            msg = "Group name '%s' is too long %d>%d" % (
                name,
                len(name),
                MED_LNAME_SIZE,
            )
            raise RuntimeError(msg)

    def add_node(self, idx, coords):
        self._corresponding_nodes[idx] = len(self.nodes)
        self.nodes.append(coords)

    def add_cell(self, idx, medcoupling_cell_type, cell_nodes):
        import MEDLoader as ml

        cell_dim = ml.MEDCouplingUMesh.GetDimensionOfGeometricType(
            medcoupling_cell_type
        )
        cell_nodes_med = tuple(self._corresponding_nodes[k] for k in cell_nodes)

        key = "%dD" % cell_dim
        if not key in self.cells:
            self.cells[key] = []
        if not key in self._corresponding_cells:
            self._corresponding_cells[key] = {}

        self.cells[key].append((medcoupling_cell_type, cell_nodes_med))
        assert idx not in self._corresponding_cells[key]
        self._corresponding_cells[key][idx] = len(self._corresponding_cells[key])

    def add_group_nodes(self, group_name, group_nodes):
        self._check_group_name(group_name)
        self.groups_n[group_name] = tuple(
            self._corresponding_nodes[k] for k in group_nodes
        )

    def add_group_cells(self, group_name, group_cells):
        self._check_group_name(group_name)
        for cell in group_cells:
            for dim in self.cells.keys():
                if cell in self._corresponding_cells[dim]:
                    if not dim in self.groups_e:
                        self.groups_e[dim] = {}
                    if not group_name in self.groups_e[dim]:
                        self.groups_e[dim][group_name] = []
                    self.groups_e[dim][group_name].append(
                        self._corresponding_cells[dim][cell]
                    )

    def _make_continuous(self):

        self.groups_e_continuous = OrderedDict()
        self.cells_continuous = OrderedDict()

        cells_shift = 0  # Variable pour la creation d'une numérotation globale
        for dim in sorted(self.cells.keys())[::-1]:
            for j, (medcoupling_cell_type, element_nodes_med) in enumerate(
                self.cells[dim]
            ):
                self.cells_continuous[cells_shift + j] = (
                    medcoupling_cell_type,
                    element_nodes_med,
                )

            if dim in self.groups_e:
                for group, values in self.groups_e[dim].items():
                    if group in self.groups_e_continuous:
                        for v in values:
                            self.groups_e_continuous[group].append(cells_shift + v)
                    else:
                        self.groups_e_continuous[group] = [
                            cells_shift + v for v in values
                        ]

            cells_shift += j + 1

    def create_UMesh(self):
        import MEDLoader as ml

        logger.debug("Create MED mesh.")

        self.umesh = ml.MEDFileUMesh()
        coords = ml.DataArrayDouble(self.nodes)

        # Les clés de elements correspondent aux dimensions dans le maillage
        for dim in self.dimensions:
            level = self.levels[dim]
            logger.debug(" Level : %d" % level)

            tic = time.perf_counter()
            mesh_at_current_level = ml.MEDCouplingUMesh(self.mesh_name, int(dim[0]))
            mesh_at_current_level.setCoords(coords)
            number_of_elements_at_level = len(self.cells[dim])
            mesh_at_current_level.allocateCells(number_of_elements_at_level)
            toc = time.perf_counter()
            logger.debug("  Set nodes (in %0.4f seconds)" % (toc - tic))
            tic = time.perf_counter()

            # Elements par niveau, avec renumerotation au passage
            for medcoupling_type, element_nodes_med in self.cells[dim]:
                number_of_nodes_current_element = (
                    ml.MEDCouplingUMesh.GetNumberOfNodesOfGeometricType(
                        medcoupling_type
                    )
                )
                mesh_at_current_level.insertNextCell(
                    medcoupling_type, number_of_nodes_current_element, element_nodes_med
                )

            mesh_at_current_level.finishInsertingCells()
            o2n = mesh_at_current_level.sortCellsInMEDFileFrmt()
            mesh_at_current_level.checkConsistencyLight()
            self.umesh.setMeshAtLevel(level, mesh_at_current_level)
            toc = time.perf_counter()
            logger.debug(
                "  Add %d elements (in %0.4f seconds)"
                % (number_of_elements_at_level, toc - tic)
            )
            tic = time.perf_counter()

            # Groupes d'elements par niveau
            try:
                groups_e_at_level = []
                for group_name, group_elements in self.groups_e[dim].items():
                    group_medcoupling = ml.DataArrayInt(group_elements)
                    group_medcoupling.transformWithIndArr(o2n)
                    group_medcoupling.setName(group_name.strip('"').strip("'"))
                    group_medcoupling.sort()
                    groups_e_at_level.append(group_medcoupling)
                self.umesh.setGroupsAtLevel(level, groups_e_at_level)
            except KeyError:
                # On peut ne pas avoir de groupes de mailles d'une certaine dimension
                pass
            toc = time.perf_counter()
            logger.debug(
                "  Add %d groups of elements (in %0.4f seconds)"
                % (len(groups_e_at_level), toc - tic)
            )

        tic = time.perf_counter()

        # Groupes de noeuds
        groups_n_at_level = []
        for group_name, group_nodes in self.groups_n.items():
            group_medcoupling = ml.DataArrayInt(group_nodes)
            group_medcoupling.setName(group_name.strip('"').strip("'"))
            group_medcoupling.sort()
            groups_n_at_level.append(group_medcoupling)
        self.umesh.setGroupsAtLevel(
            1, groups_n_at_level
        )  # Groupes de noeuds au niveau 1
        self.umesh.setName(self.mesh_name)
        toc = time.perf_counter()
        logger.debug(" Level : 1")
        logger.debug(
            "  Add %d groups of nodes (in %0.4f seconds)"
            % (len(groups_n_at_level), toc - tic)
        )

        tic = time.perf_counter()
        self.umesh.rearrangeFamilies()
        toc = time.perf_counter()
        logger.debug(" Sort families (in %0.4f seconds)" % (toc - tic))

        return self.umesh
