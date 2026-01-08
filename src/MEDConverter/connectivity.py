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

from collections import OrderedDict

class ConnectivityRenumberer:

    # Les items d'une liste indiquent à quelle position MED se trouve le noeud SYSTUS correspondant à l'index dans la liste.
    # e.g. pour le QUAD8 :
    # Le noeud 0 SYSTUS correspond au noeud 0 MED
    # Le noeud 1 SYSTUS correspond au noeud 4 MED
    # Le noeud 2 SYSTUS correspond au noeud 1 MED
    # Le noeud 3 SYSTUS correspond au noeud 5 MED
    # Le noeud 4 SYSTUS correspond au noeud 2 MED
    # Le noeud 5 SYSTUS correspond au noeud 6 MED
    # Le noeud 6 SYSTUS correspond au noeud 3 MED
    # Le noeud 7 SYSTUS correspond au noeud 7 MED

    # fmt: off
    _systus = {
        "POINT1": [0],
        "SEG2": range(2),
        "TRI3": range(3),
        "QUAD4": range(4),
        "TETRA4": [0, 2, 1, 3],
        "HEXA8": [0, 3, 2, 1, 4, 7, 6, 5],
        "PYRA5": [0, 3, 2, 1, 4],
        "PENTA6": [0, 2, 1, 3, 5, 4],
        "SEG3": [0, 2, 1],
        "TRI6": [0, 3, 1, 4, 2, 5],
        "QUAD8": [0, 4, 1, 5, 2, 6, 3, 7],
        "TETRA10": [0, 6, 2, 5, 1, 4, 7, 9, 8, 3],
        "HEXA20": [0, 11, 3, 10, 2, 9, 1, 8, 16, 19, 18, 17, 4, 15, 7, 14, 6, 13, 5, 12],
        "PYRA13": [0, 8, 3, 7, 2, 6, 1, 5, 9, 12, 11, 10, 4],
        "PENTA15": [0, 8, 2, 7, 1, 6, 12, 14, 13, 3, 11, 5, 10, 4, 9],
        "SEG4": [0, 2, 3, 1],
    }

    _abaqus = {
        "POINT1": [0],
        "SEG2": range(2),
        "SEG3": [0, 2, 1],
        "TRI3": range(3),
        "TRI6": range(6),
        "QUAD4": range(4),
        "QUAD8": range(8),
        "QUAD9": range(9),
        "TETRA4": [0, 2, 1, 3],
        "TETRA10": [0, 2, 1, 3, 6, 5, 4, 7, 9, 8],
        "HEXA8": [0, 3, 2, 1, 4, 7, 6, 5],
        "HEXA20": [0, 3, 2, 1, 4, 7, 6, 5, 11, 10, 9, 8, 15, 14, 13, 12, 16, 19, 18, 17],
        "HEXA27": [0, 3, 2, 1, 4, 7, 6, 5, 11, 10, 9, 8, 15, 14, 13, 12, 16, 19, 18, 17, 26, 20, 25, 24, 23, 22, 21],
        "PENTA6": [0, 2, 1, 3, 5, 4],
        "PENTA15": [0, 2, 1, 3, 5, 4, 8, 7, 6, 11, 10, 9, 12, 14, 13],
        "PENTA18": [0, 2, 1, 3, 5, 4, 8, 7, 6, 11, 10, 9, 12, 14, 13, 17, 16, 15],
        "PYRA5": [0, 3, 2, 1, 4],
        "PYRA13": [0, 3, 2, 1, 4, 8, 7, 6, 5, 9, 12, 11, 10],
    }

    _aster = {
        "POINT1": [0],
        "SEG2": range(2),
        "TRI3": range(3),
        "QUAD4": range(4),
        "TETRA4": [0, 2, 1, 3],
        "HEXA8": [0, 3, 2, 1, 4, 7, 6, 5],
        "PYRA5": [0, 3, 2, 1, 4],
        "PENTA6": [0, 2, 1, 3, 5, 4],
        "SEG3": range(3),
        "TRI6": range(6),
        "QUAD8": range(8),
        "TETRA10": [0, 2, 1, 3, 6, 5, 4, 7, 9, 8],
        "HEXA20": [0, 3, 2, 1, 4, 7, 6, 5, 11, 10, 9, 8, 16, 19, 18, 17, 15, 14, 13, 12],
        "PYRA13": [0, 3, 2, 1, 4, 8, 7, 6, 5, 9, 12, 11, 10],
        "PENTA15": [0, 2, 1, 3, 5, 4, 8, 7, 6, 12, 14, 13, 11, 10, 9],
        "SEG4": range(4),
        "TRI7": range(7),
        "QUAD9": range(9),
        "PENTA18": [0, 2, 1, 3, 5, 4, 8, 7, 6, 12, 14, 13, 11, 10, 9, 17, 16, 15],
        "HEXA27": [0, 3, 2, 1, 4, 7, 6, 5, 11, 10, 9, 8, 16, 19, 18, 17, 15, 14, 13, 12, 20, 24, 23, 22, 21, 25, 26],
    }

    _tetgen = {"TRI3": range(3), "TETRA4": range(4)}

    _zset = {
        "SEG2": range(2),
        "SEG3": [0, 2, 1],
        "TRI3": range(3),
        "TRI6": [2, 5, 0, 3, 1, 4],
        "QUAD4": [0, 1, 2, 3],
        "QUAD8": [0, 4, 1, 5, 2, 6, 3, 7],
        "TETRA4": [1, 3, 0, 2],
        "TETRA10": [1, 3, 0, 8, 7, 4, 5, 9, 6, 2],
        "HEXA8": [0, 3, 2, 1, 4, 7, 6, 5],
        "HEXA20": [0, 11, 3, 10, 2, 9, 1, 8, 16, 19, 18, 17, 4, 15, 7, 14, 6, 13, 5, 12],
        "PENTA6": [1, 0, 2, 4, 3, 5],
        "PENTA15": [1, 6, 0, 8, 2, 7, 13, 12, 14, 4, 9, 3, 11, 5, 10],
        "PYRA5": [0, 3, 2, 1, 4],
        "PYRA13": [0, 3, 2, 1, 4, 8, 7, 6, 5, 9, 12, 11, 10],
    }

    _ansys = {
        "POINT1": [0],
        "SEG2": range(2),
        "TRI3": range(3),
        "QUAD4": range(4),
        "HEXA8": range(8),
        "PENTA6": range(6),
        "TETRA4": range(4),
        "PYRA5": range(5),
        "SEG3": range(3),
        "TRI6": range(6),
        "TETRA10": range(10),
        "QUAD8": range(8),
        "PYRA13": range(13),
        "PENTA15": range(15),
        "HEXA20": range(20),
    }

    _med_types = "POINT1 SEG2 TRI3 QUAD4 TETRA4 HEXA8 PYRA5 PENTA6 SEG3 TRI6 QUAD8 TETRA10 HEXA20 PYRA13 PENTA15 SEG4 TRI7 QUAD9 PENTA18 HEXA27".split()
    # fmt: on

    def __init__(self, code):
        import MEDLoader as ml

        self._connectivity_med_to_external = {}
        self._connectivity_external_to_med = {}

        try:
            connectivity = getattr(self, "_{}".format(code.lower()))
            assert set(connectivity.keys()) <= set(self._med_types)

            for elem, nodes in connectivity.items():
                elem_mc = getattr(ml, "NORM_%s" % elem)
                self._connectivity_med_to_external[elem_mc] = OrderedDict()
                self._connectivity_external_to_med[elem_mc] = OrderedDict()
                tmp = {}

                for i, val in enumerate(nodes):
                    tmp[val] = i
                    self._connectivity_med_to_external[elem_mc][i] = val
                for i in sorted(tmp):
                    self._connectivity_external_to_med[elem_mc][i] = tmp[i]

        except AttributeError as exc:
            raise RuntimeError("Unknown connectivity {}".format(code)) from exc

    def external_to_medcoupling(self, medcoupling_type, nodes):
        try:
            return tuple(
                nodes[self._connectivity_external_to_med[medcoupling_type][i]]
                for i in self._connectivity_external_to_med[medcoupling_type]
            )
        except KeyError as exc:
            raise RuntimeError(
                "Unsupported element type %s" % medcoupling_type
            ) from exc

    def medcoupling_to_external(self, medcoupling_type, nodes):
        try:
            return tuple(
                nodes[self._connectivity_med_to_external[medcoupling_type][i]]
                for i in self._connectivity_med_to_external[medcoupling_type]
            )
        except KeyError as exc:
            raise RuntimeError(
                "Unsupported element type %s" % medcoupling_type
            ) from exc
