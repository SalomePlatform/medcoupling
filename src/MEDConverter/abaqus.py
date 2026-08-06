#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (C) 2023-2026 CEA, EDF
# SPDX-License-Identifier: LGPL-2.1-or-later

"""Abaqus input mesh reader for MEDConverter.

The reader intentionally converts mesh topology, sets and discrete surfaces.
Analysis/model cards that do not alter the mesh are skipped. Unsupported
keywords are ignored with an explicit warning.
"""

from collections import OrderedDict
from dataclasses import dataclass, field
import math
import os
import time

from .logger import logger
from .MEDConverterMesh import MEDConverterMesh
from .cells import CellsTypeConverter
from .connectivity import ConnectivityRenumberer


@dataclass
class AbaqusNode:
    id: int
    coordinates: tuple

    def getId(self):
        return self.id

    def getCoordinates(self):
        if len(self.coordinates) != 3:
            raise RuntimeError("Coordinates have to have 3 items")
        return self.coordinates


@dataclass
class AbaqusElement:
    type: str
    id: int
    nodes: tuple

    def getId(self):
        return self.id

    def getType(self):
        return self.type

    def getNodes(self):
        return self.nodes


@dataclass
class AbaqusGroup:
    name: str
    instance: str = ""
    group: list = field(default_factory=list)

    def getName(self):
        return self.name

    def getInstance(self):
        return self.instance

    def getGroup(self):
        return self.group

    def addGroup(self, values):
        self.group.extend(values)


@dataclass
class AbaqusSurface:
    name: str
    type: str
    entries: list
    instance: str = ""


@dataclass
class AbaqusEntity:
    name: str = ""
    nodes: OrderedDict = field(default_factory=OrderedDict)
    elements: OrderedDict = field(default_factory=OrderedDict)
    nsets: OrderedDict = field(default_factory=OrderedDict)
    elsets: OrderedDict = field(default_factory=OrderedDict)
    surfaces: OrderedDict = field(default_factory=OrderedDict)
    pending_nsets: list = field(default_factory=list)


@dataclass
class AbaqusInstance:
    name: str
    part_name: str
    translation: tuple = None
    rotation: tuple = None


# Face indices follow Abaqus local numbering. Values are internal topology
# aliases resolved centrally by CellsTypeConverter("ABAQUS").
FACES = {
    "TRI3": [
        ("SEG2", (1, 2)),
        ("SEG2", (2, 3)),
        ("SEG2", (3, 1)),
    ],
    "TRI6": [
        ("SEG3", (1, 2, 4)),
        ("SEG3", (2, 3, 5)),
        ("SEG3", (3, 1, 6)),
    ],
    "QUAD4": [
        ("SEG2", (1, 2)),
        ("SEG2", (2, 3)),
        ("SEG2", (3, 4)),
        ("SEG2", (4, 1)),
    ],
    "QUAD8": [
        ("SEG3", (1, 2, 5)),
        ("SEG3", (2, 3, 6)),
        ("SEG3", (3, 4, 7)),
        ("SEG3", (4, 1, 8)),
    ],
    "QUAD9": [
        ("SEG3", (1, 2, 5)),
        ("SEG3", (2, 3, 6)),
        ("SEG3", (3, 4, 7)),
        ("SEG3", (4, 1, 8)),
    ],
    "TETRA4": [
        ("TRI3", (1, 2, 3)),
        ("TRI3", (1, 4, 2)),
        ("TRI3", (2, 4, 3)),
        ("TRI3", (3, 4, 1)),
    ],
    "TETRA10": [
        ("TRI6", (1, 2, 3, 5, 6, 7)),
        ("TRI6", (1, 4, 2, 8, 9, 5)),
        ("TRI6", (2, 4, 3, 9, 10, 6)),
        ("TRI6", (3, 4, 1, 10, 8, 7)),
    ],
    "HEXA8": [
        ("QUAD4", (1, 2, 3, 4)),
        ("QUAD4", (5, 8, 7, 6)),
        ("QUAD4", (1, 5, 6, 2)),
        ("QUAD4", (2, 6, 7, 3)),
        ("QUAD4", (3, 7, 8, 4)),
        ("QUAD4", (4, 8, 5, 1)),
    ],
    "HEXA20": [
        ("QUAD8", (1, 2, 3, 4, 9, 10, 11, 12)),
        ("QUAD8", (5, 8, 7, 6, 16, 15, 14, 13)),
        ("QUAD8", (1, 5, 6, 2, 17, 13, 18, 9)),
        ("QUAD8", (2, 6, 7, 3, 18, 14, 19, 10)),
        ("QUAD8", (3, 7, 8, 4, 19, 15, 20, 11)),
        ("QUAD8", (4, 8, 5, 1, 20, 16, 17, 12)),
    ],
    "HEXA27": [
        ("QUAD9", (1, 2, 3, 4, 9, 10, 11, 12, 22)),
        ("QUAD9", (5, 8, 7, 6, 16, 15, 14, 13, 23)),
        ("QUAD9", (1, 5, 6, 2, 17, 13, 18, 9, 24)),
        ("QUAD9", (2, 6, 7, 3, 18, 14, 19, 10, 25)),
        ("QUAD9", (3, 7, 8, 4, 19, 15, 20, 11, 26)),
        ("QUAD9", (4, 8, 5, 1, 20, 16, 17, 12, 27)),
    ],
    "PYRA5": [
        ("QUAD4", (1, 2, 3, 4)),
        ("TRI3", (1, 5, 2)),
        ("TRI3", (2, 5, 3)),
        ("TRI3", (3, 5, 4)),
        ("TRI3", (4, 5, 1)),
    ],
    "PYRA13": [
        ("QUAD8", (1, 2, 3, 4, 6, 7, 8, 9)),
        ("TRI6", (1, 5, 2, 10, 11, 6)),
        ("TRI6", (2, 5, 3, 11, 12, 7)),
        ("TRI6", (3, 5, 4, 12, 13, 8)),
        ("TRI6", (4, 5, 1, 13, 10, 9)),
    ],
    "PENTA6": [
        ("TRI3", (1, 2, 3)),
        ("TRI3", (4, 6, 5)),
        ("QUAD4", (1, 4, 5, 2)),
        ("QUAD4", (2, 5, 6, 3)),
        ("QUAD4", (3, 6, 4, 1)),
    ],
    "PENTA15": [
        ("TRI6", (1, 2, 3, 7, 8, 9)),
        ("TRI6", (4, 6, 5, 12, 11, 10)),
        ("QUAD8", (1, 4, 5, 2, 13, 10, 14, 7)),
        ("QUAD8", (2, 5, 6, 3, 14, 11, 15, 8)),
        ("QUAD8", (3, 6, 4, 1, 15, 12, 13, 9)),
    ],
    "PENTA18": [
        ("TRI6", (1, 2, 3, 7, 8, 9)),
        ("TRI6", (4, 6, 5, 12, 11, 10)),
        ("QUAD9", (1, 4, 5, 2, 13, 10, 14, 7, 16)),
        ("QUAD9", (2, 5, 6, 3, 14, 11, 15, 8, 17)),
        ("QUAD9", (3, 6, 4, 1, 15, 12, 13, 9, 18)),
    ],
}


class MEDConverterAbaqus(MEDConverterMesh):
    _mesh_generation_keywords = {
        "*NGEN",
        "*NFILL",
        "*NMAP",
        "*NCOPY",
        "*ELGEN",
        "*ELCOPY",
    }

    def __init__(self):
        super().__init__()
        self.parts = OrderedDict()
        self.instances = []
        # Mesh data can be declared directly between *INSTANCE and
        # *END INSTANCE, independently of the referenced *PART.
        self.instance_entities = OrderedDict()
        self.root = AbaqusEntity("assembly")
        self._include_stack = []
        # Containers remain case-insensitive internally, while these maps
        # preserve the spelling used in the Abaqus file for MED group names.
        self._nset_names = {}
        self._elset_names = {}
        self._type_converter = CellsTypeConverter("ABAQUS")
        self._connectivity_converter = ConnectivityRenumberer("ABAQUS")
        self._element_info_cache = {}
        self._surface_topology_cache = {}
        self._surface_med_type_cache = {}

    @staticmethod
    def convert_abaqus_to_med(filename_abaqus, verbose=False):
        tic = time.perf_counter()
        converter = MEDConverterAbaqus()
        converter.verbose = verbose
        converter.read_abaqus_mesh(filename_abaqus)
        converter.create_UMesh()
        logger.debug("Mesh converted (in %0.4f seconds)" % (time.perf_counter() - tic))
        return converter.umesh

    @staticmethod
    def _keyword(line):
        stripped = line.lstrip()
        return stripped.startswith("*") and not stripped.startswith("**")

    @staticmethod
    def _params(line):
        result = {}
        for index, token in enumerate(line.strip().split(",")):
            token = token.strip()
            if "=" in token:
                key, value = token.split("=", 1)
                result[key.strip().upper()] = value.strip()
            else:
                result[token.upper()] = None
        return result

    @staticmethod
    def _data_tokens(line):
        return [
            token.strip()
            for token in line.strip().rstrip(",").split(",")
            if token.strip()
        ]

    def _load_lines(self, filename):
        path = os.path.realpath(filename)
        if path in self._include_stack:
            raise RuntimeError("Cyclic Abaqus include detected: %s" % path)
        self._include_stack.append(path)
        try:
            with open(path, "r", encoding=self._get_file_encoding(path)) as stream:
                source = stream.readlines()
            result = []
            i = 0
            while i < len(source):
                line = source[i]
                upper = line.lstrip().upper()
                if upper.startswith("*INCLUDE"):
                    params = self._params(line)
                    input_name = params.get("INPUT") or params.get("FILE")
                    if not input_name:
                        raise RuntimeError("*INCLUDE without INPUT in %s" % path)
                    include = (
                        input_name
                        if os.path.isabs(input_name)
                        else os.path.join(os.path.dirname(path), input_name)
                    )
                    result.extend(self._load_lines(include))
                else:
                    result.append((path, i + 1, line.rstrip("\n\r")))

                    # Abaqus also permits external data files directly on
                    # *NODE and *ELEMENT cards.  These files contain only data
                    # records, not an *INCLUDE card, and must therefore be
                    # inserted immediately after their owning keyword.
                    if upper.startswith("*NODE") or upper.startswith("*ELEMENT"):
                        params = self._params(line)
                        input_name = params.get("INPUT") or params.get("FILE")
                        if input_name:
                            include = (
                                input_name
                                if os.path.isabs(input_name)
                                else os.path.join(os.path.dirname(path), input_name)
                            )
                            result.extend(self._load_external_data(include))
                i += 1
            return result
        finally:
            self._include_stack.pop()

    def _load_external_data(self, filename):
        """Load a data-only file referenced by INPUT on *NODE/*ELEMENT."""
        path = os.path.realpath(filename)
        if path in self._include_stack:
            raise RuntimeError("Cyclic Abaqus include detected: %s" % path)
        self._include_stack.append(path)
        try:
            with open(path, "r", encoding=self._get_file_encoding(path)) as stream:
                return [
                    (path, lineno, line.rstrip("\n\r"))
                    for lineno, line in enumerate(stream, 1)
                ]
        finally:
            self._include_stack.pop()

    def read_abaqus_mesh(self, filename):
        logger.debug("Reading ABAQUS mesh file : %s" % filename)
        self._reset_structures()
        self.filename = filename
        self.mesh_name = os.path.splitext(os.path.basename(filename))[0]
        self.space_dim = 3
        lines = self._load_lines(filename)
        self._parse(lines)
        self._build_mesh()

    def _data_block(self, lines, index):
        rows = []
        i = index + 1
        while i < len(lines) and not self._keyword(lines[i][2]):
            text = lines[i][2].strip()
            if text and not text.startswith("**"):
                rows.append((lines[i][0], lines[i][1], text))
            i += 1
        return rows, i

    def _skip_data_block(self, lines, index):
        """Return the next keyword index without storing ignored data rows."""
        i = index + 1
        keyword = self._keyword
        while i < len(lines) and not keyword(lines[i][2]):
            i += 1
        return i

    @staticmethod
    def _warn_unsupported_keyword(path, lineno, keyword):
        logger.warning(
            "Unsupported Abaqus keyword ignored at %s:%d: %s",
            path,
            lineno,
            keyword,
        )

    @staticmethod
    def _logical_rows(rows):
        """Join Abaqus physical lines terminated by a comma.

        Abaqus uses a trailing comma to continue node and element records on
        the next physical line.  Returning the first source location keeps
        diagnostics attached to the start of the logical record.
        """
        logical = []
        pending = []
        first_path = None
        first_lineno = None
        for path, lineno, row in rows:
            if not pending:
                first_path, first_lineno = path, lineno
            pending.extend(MEDConverterAbaqus._data_tokens(row))
            if not row.rstrip().endswith(","):
                logical.append((first_path, first_lineno, pending))
                pending = []
        if pending:
            logical.append((first_path, first_lineno, pending))
        return logical

    def _element_info(self, element_type):
        element_type = element_type.upper()
        info = self._element_info_cache.get(element_type)
        if info is None:
            med_type = self._type_converter.external_to_medcoupling(element_type)
            permutation = tuple(
                self._connectivity_converter._connectivity_external_to_med[med_type]
            )
            info = med_type, permutation
            self._element_info_cache[element_type] = info
        return info

    def _element_rows(self, rows, element_type):
        """Build Abaqus element records using the topology node count.

        A trailing comma is not sufficient to identify a continuation: some
        Abaqus writers also put one on complete records.  The connectivity
        permutation is the authoritative expected size.
        """
        try:
            _, permutation = self._element_info(element_type)
            expected_nodes = len(permutation)
        except (KeyError, TypeError, AttributeError):
            expected_nodes = None

        logical = []
        pending = []
        first_path = None
        first_lineno = None
        for path, lineno, row in rows:
            tokens = MEDConverterAbaqus._data_tokens(row)
            if not pending:
                first_path, first_lineno = path, lineno
            pending.extend(tokens)
            if expected_nodes is None:
                complete = not row.rstrip().endswith(",")
            else:
                complete = len(pending) >= expected_nodes + 1
            if complete:
                if expected_nodes is not None and len(pending) != expected_nodes + 1:
                    raise RuntimeError(
                        "Invalid %s connectivity at %s:%d: expected %d nodes, "
                        "got %d"
                        % (
                            element_type,
                            first_path,
                            first_lineno,
                            expected_nodes,
                            len(pending) - 1,
                        )
                    )
                logical.append((first_path, first_lineno, pending))
                pending = []
        if pending:
            raise RuntimeError(
                "Incomplete %s connectivity at %s:%d"
                % (
                    element_type,
                    first_path,
                    first_lineno,
                )
            )
        return logical

    def _remember_set_name(self, entity, is_node, name):
        table = self._nset_names if is_node else self._elset_names
        table[(id(entity), name.upper())] = name

    def _display_set_name(self, entity, is_node, key):
        table = self._nset_names if is_node else self._elset_names
        return table.get((id(entity), key.upper()), key)

    def _parse(self, lines):
        current = self.root
        in_assembly = False
        i = 0
        while i < len(lines):
            path, lineno, raw = lines[i]
            text = raw.strip()
            if not text or text.startswith("**") or not self._keyword(raw):
                i += 1
                continue
            params = self._params(text)
            keyword = next(iter(params))
            if keyword in self._mesh_generation_keywords:
                self._warn_unsupported_keyword(path, lineno, keyword)
                _, i = self._skip_data_block(lines, i)
                continue
            if keyword == "*PART":
                name = params.get("NAME")
                if not name:
                    raise RuntimeError("*PART without NAME at %s:%d" % (path, lineno))
                current = AbaqusEntity(name)
                self.parts[name.upper()] = current
                i += 1
                continue
            if keyword == "*END PART":
                current = self.root
                i += 1
                continue
            if keyword == "*ASSEMBLY":
                in_assembly = True
                current = self.root
                i += 1
                continue
            if keyword == "*END ASSEMBLY":
                in_assembly = False
                current = self.root
                i += 1
                continue
            if keyword == "*INSTANCE":
                rows, next_i = self._data_block(lines, i)
                name, part = params.get("NAME"), params.get("PART")
                if not name or not part:
                    raise RuntimeError(
                        "*INSTANCE requires NAME and PART at %s:%d" % (path, lineno)
                    )
                numbers = [
                    [float(x) for x in self._data_tokens(row[2])] for row in rows
                ]
                translation = (
                    tuple(numbers[0]) if numbers and len(numbers[0]) == 3 else None
                )
                rotation = None
                for values in numbers:
                    if len(values) == 7:
                        rotation = tuple(values)
                    elif len(values) not in (3,):
                        raise RuntimeError(
                            "Invalid instance transform at %s:%d" % (path, lineno)
                        )
                self.instances.append(AbaqusInstance(name, part, translation, rotation))
                current = AbaqusEntity(name)
                self.instance_entities[name.upper()] = current
                i = next_i
                continue
            if keyword == "*END INSTANCE":
                current = self.root
                i += 1
                continue
            if keyword == "*NODE":
                rows, i = self._data_block(lines, i)
                ids = []
                for row_path, rowno, values in self._logical_rows(rows):
                    if len(values) < 2:
                        raise RuntimeError("Invalid node at %s:%d" % (row_path, rowno))
                    node_id = int(values[0])
                    coords = [float(x) for x in values[1:]]
                    coords += [0.0] * (3 - len(coords))
                    if len(coords) != 3:
                        raise RuntimeError("Invalid node at %s:%d" % (row_path, rowno))
                    current.nodes[node_id] = AbaqusNode(node_id, tuple(coords))
                    ids.append(node_id)
                if params.get("NSET"):
                    self._remember_set_name(current, True, params["NSET"])
                    self._merge_set(current.nsets, params["NSET"], ids)
                continue
            if keyword == "*ELEMENT":
                rows, i = self._data_block(lines, i)
                etype = params.get("TYPE")
                if not etype:
                    raise RuntimeError(
                        "*ELEMENT without TYPE at %s:%d" % (path, lineno)
                    )
                ids = []
                for row_path, rowno, values in self._element_rows(rows, etype):
                    if len(values) < 2:
                        raise RuntimeError(
                            "Invalid element at %s:%d" % (row_path, rowno)
                        )
                    element_id = int(values[0])
                    nodes = tuple(
                        value.upper() if "." in value else int(value)
                        for value in values[1:]
                    )
                    current.elements[element_id] = AbaqusElement(
                        etype.upper(), element_id, nodes
                    )
                    ids.append(element_id)
                if params.get("ELSET"):
                    self._remember_set_name(current, False, params["ELSET"])
                    self._merge_set(current.elsets, params["ELSET"], ids)
                continue
            if keyword in ("*NSET", "*ELSET"):
                rows, i = self._data_block(lines, i)
                is_node = keyword == "*NSET"
                set_name = params.get("NSET" if is_node else "ELSET")
                if not set_name:
                    raise RuntimeError(
                        "%s without name at %s:%d" % (keyword, path, lineno)
                    )
                target = current.nsets if is_node else current.elsets
                self._remember_set_name(current, is_node, set_name)
                if is_node and (params.get("ELSET") or params.get("SURFACE")):
                    current.pending_nsets.append(
                        (
                            set_name,
                            "ELSET" if params.get("ELSET") else "SURFACE",
                            params.get("ELSET") or params.get("SURFACE"),
                        )
                    )
                    continue
                values = []
                for _, _, row in rows:
                    values.extend(self._data_tokens(row))
                expanded = []
                if "GENERATE" in params:
                    for row_path, rowno, row in rows:
                        values = self._data_tokens(row)
                        if len(values) not in (2, 3):
                            raise RuntimeError(
                                "Invalid %s GENERATE range at %s:%d: "
                                "expected first,last[,step], got %r"
                                % (keyword, row_path, rowno, values)
                            )
                        first = int(values[0])
                        last = int(values[1])
                        step = int(values[2]) if len(values) == 3 else 1
                        if step == 0:
                            raise RuntimeError(
                                "Invalid %s GENERATE range at %s:%d: "
                                "step must not be zero" % (keyword, row_path, rowno)
                            )
                        expanded.extend(range(first, last + 1, step))
                else:
                    values = []
                    for _, _, row in rows:
                        values.extend(self._data_tokens(row))
                    for value in values:
                        try:
                            expanded.append(int(value))
                        except ValueError:
                            source = target.get(value.upper())
                            if source is None:
                                expanded.append(value)
                            else:
                                expanded.extend(source)
                # At assembly level, INSTANCE= scopes unqualified labels to
                # that instance. Keep qualified labels until global numbering
                # is available in _build_mesh().
                instance_name = params.get("INSTANCE")
                if instance_name:
                    expanded = [
                        value if "." in str(value) else "%s.%s" % (instance_name, value)
                        for value in expanded
                    ]
                if "UNSORTED" not in params and all(
                    isinstance(x, int) for x in expanded
                ):
                    expanded = sorted(set(expanded))
                self._merge_set(target, set_name, expanded)
                continue
            if keyword == "*SURFACE":
                rows, i = self._data_block(lines, i)
                name = params.get("NAME")
                surface_type = (params.get("TYPE") or "ELEMENT").upper()
                if surface_type not in ("NODE", "ELEMENT"):
                    logger.warning(
                        "Unsupported Abaqus analytical surface ignored at %s:%d: "
                        "NAME=%s, TYPE=%s",
                        path,
                        lineno,
                        name or "<unnamed>",
                        surface_type,
                    )
                    continue
                entries = [self._data_tokens(row) for _, _, row in rows]
                current.surfaces[name.upper()] = AbaqusSurface(
                    name, surface_type, entries, params.get("INSTANCE") or ""
                )
                continue
            _, i = self._data_block(lines, i)

    @staticmethod
    def _merge_set(container, name, values):
        key = name.upper()
        container.setdefault(key, []).extend(values)

    @staticmethod
    def _rotation_matrix(rotation):
        if rotation is None:
            return None, None
        center = rotation[:3]
        axis = [rotation[3 + i] - center[i] for i in range(3)]
        norm = math.sqrt(sum(x * x for x in axis))
        if norm == 0:
            raise RuntimeError("Abaqus instance rotation axis must not be zero")
        x, y, z = [v / norm for v in axis]
        a = math.radians(rotation[6])
        c, s, q = math.cos(a), math.sin(a), 1 - math.cos(a)
        return center, (
            (c + x * x * q, x * y * q - z * s, x * z * q + y * s),
            (y * x * q + z * s, c + y * y * q, y * z * q - x * s),
            (z * x * q - y * s, z * y * q + x * s, c + z * z * q),
        )

    @classmethod
    def _transform(cls, point, translation, rotation):
        p = list(point)
        if translation:
            p = [p[i] + translation[i] for i in range(3)]
        center, matrix = cls._rotation_matrix(rotation)
        if matrix:
            v = [p[i] - center[i] for i in range(3)]
            p = [
                center[i] + sum(matrix[i][j] * v[j] for j in range(3)) for i in range(3)
            ]
        return tuple(p)

    @staticmethod
    def _transform_precomputed(point, translation, center, matrix):
        x, y, z = point
        if translation is not None:
            x += translation[0]
            y += translation[1]
            z += translation[2]
        if matrix is None:
            return x, y, z
        cx, cy, cz = center
        vx, vy, vz = x - cx, y - cy, z - cz
        return (
            cx + matrix[0][0] * vx + matrix[0][1] * vy + matrix[0][2] * vz,
            cy + matrix[1][0] * vx + matrix[1][1] * vy + matrix[1][2] * vz,
            cz + matrix[2][0] * vx + matrix[2][1] * vy + matrix[2][2] * vz,
        )

    @staticmethod
    def _group_values(groups, name):
        """Return a group by name without changing its original case."""
        wanted = name.upper()
        for group_name, values in groups.items():
            if group_name.upper() == wanted:
                return values
        return []

    def _build_mesh(self):
        type_converter = self._type_converter
        node_groups, cell_groups = OrderedDict(), OrderedDict()
        element_records = {}
        qualified_nodes = {}
        qualified_cells = {}
        next_node, next_cell = 1, 1

        entities = []
        if self.parts:
            for instance in self.instances or [
                AbaqusInstance("Instance_" + p.name, p.name)
                for p in self.parts.values()
            ]:
                part = self.parts.get(instance.part_name.upper())
                if part is None:
                    raise RuntimeError("Part not found: %s" % instance.part_name)
                inline = self.instance_entities.get(instance.name.upper())
                # Abaqus permits nodes, elements and sets directly inside an
                # instance. Use that entity when it carries mesh content;
                # otherwise instantiate the referenced part as before.
                if inline is not None and (
                    inline.nodes
                    or inline.elements
                    or inline.nsets
                    or inline.elsets
                    or inline.surfaces
                ):
                    entity = inline
                else:
                    entity = part
                entities.append(
                    (instance.name, entity, instance.translation, instance.rotation)
                )
        entities.append(("", self.root, None, None))

        for instance_name, entity, translation, rotation in entities:
            node_map, cell_map = {}, {}
            cells_by_type = OrderedDict()
            rotation_center, rotation_matrix = self._rotation_matrix(rotation)
            for local_id, node in entity.nodes.items():
                node_map[local_id] = next_node
                if instance_name:
                    qualified_nodes["%s.%s" % (instance_name.upper(), local_id)] = (
                        next_node
                    )
                coordinates = self._transform_precomputed(
                    node.coordinates,
                    translation,
                    rotation_center,
                    rotation_matrix,
                )
                self.add_node(next_node, coordinates)
                next_node += 1
            for local_id, element in entity.elements.items():
                external_nodes = []
                for node_label in element.nodes:
                    if isinstance(node_label, int):
                        global_node = node_map.get(node_label)
                    else:
                        global_node = qualified_nodes.get(node_label)
                    if global_node is None:
                        raise RuntimeError(
                            "Element %s references missing node %s"
                            % (local_id, node_label)
                        )
                    external_nodes.append(global_node)
                med_type, _ = self._element_info(element.type)
                # Keep ConnectivityRenumberer as the authoritative conversion.
                # Its internal table must not be applied directly as a forward
                # permutation: doing so changes TETRA, PYRA and HEXA ordering.
                med_nodes = self._connectivity_converter.external_to_medcoupling(
                    med_type, tuple(external_nodes)
                )
                cell_map[local_id] = next_cell
                if instance_name:
                    qualified_cells["%s.%s" % (instance_name.upper(), local_id)] = (
                        next_cell
                    )
                self.add_cell(next_cell, med_type, med_nodes)
                element_records[next_cell] = (element.type, external_nodes)
                cells_by_type.setdefault(element.type, []).append(next_cell)
                next_cell += 1
            prefix = ""
            for name, values in entity.nsets.items():
                resolved = []
                for value in values:
                    label = str(value)
                    if "." in label:
                        global_id = qualified_nodes.get(label.upper())
                    else:
                        try:
                            global_id = node_map.get(int(label))
                        except ValueError:
                            global_id = None
                    if global_id is not None:
                        resolved.append(global_id)
                    else:
                        logger.warning(
                            "Node reference ignored in NSET %s: %s",
                            self._display_set_name(entity, True, name),
                            label,
                        )
                display_name = self._display_set_name(entity, True, name)
                node_groups.setdefault(prefix + display_name, []).extend(resolved)
            for name, values in entity.elsets.items():
                resolved = []
                for value in values:
                    label = str(value)
                    if "." in label:
                        global_id = qualified_cells.get(label.upper())
                    else:
                        try:
                            global_id = cell_map.get(int(label))
                        except ValueError:
                            global_id = None
                    if global_id is not None:
                        resolved.append(global_id)
                    else:
                        logger.warning(
                            "Element reference ignored in ELSET %s: %s",
                            self._display_set_name(entity, False, name),
                            label,
                        )
                display_name = self._display_set_name(entity, False, name)
                cell_groups.setdefault(prefix + display_name, []).extend(resolved)

            # Always create one uniform Grp_FE_<type> group for every Abaqus
            # element type, even when an explicit ELSET covers the same cells.
            for element_type, type_cells in cells_by_type.items():
                cell_groups.setdefault(prefix + "Grp_FE_" + element_type, []).extend(
                    type_cells
                )

            self._resolve_surfaces(
                entity,
                prefix,
                node_map,
                cell_map,
                element_records,
                node_groups,
                cell_groups,
                type_converter,
                next_cell,
            )
            next_cell = max(element_records, default=next_cell - 1) + 1
            for set_name, source_type, source_name in entity.pending_nsets:
                if source_type == "ELSET":
                    ids = self._group_values(cell_groups, prefix + source_name)
                    nodes = []
                    for cid in ids:
                        nodes.extend(element_records[cid][1])
                    node_groups[prefix + set_name.upper()] = list(dict.fromkeys(nodes))
                elif source_type == "SURFACE":
                    ids = self._group_values(cell_groups, prefix + source_name)
                    nodes = []
                    for cid in ids:
                        nodes.extend(element_records[cid][1])
                    node_groups[prefix + set_name.upper()] = list(dict.fromkeys(nodes))

        for name, values in node_groups.items():
            values = list(dict.fromkeys(values))
            if values:
                self.add_group_nodes(name, values)
        for name, values in cell_groups.items():
            values = list(dict.fromkeys(values))
            if values:
                self.add_group_cells(name, values)

    def _resolve_surfaces(
        self,
        entity,
        prefix,
        node_map,
        cell_map,
        records,
        node_groups,
        cell_groups,
        type_converter,
        next_cell,
    ):
        for _, surface in entity.surfaces.items():
            name = prefix + surface.name
            if surface.type == "NODE":
                nodes = []
                for entry in surface.entries:
                    label = entry[0]
                    try:
                        nodes.append(node_map[int(label)])
                    except ValueError:
                        nodes.extend(self._group_values(node_groups, prefix + label))
                node_groups[name] = list(dict.fromkeys(nodes))
                continue
            if surface.type != "ELEMENT":
                logger.warning(
                    "Unsupported Abaqus surface type ignored: %s (surface %s)",
                    surface.type,
                    surface.name,
                )
                continue
            generated = []
            for entry in surface.entries:
                if not entry:
                    continue
                label = entry[0]
                face = entry[1].upper() if len(entry) > 1 else None
                try:
                    parents = [cell_map[int(label)]]
                except ValueError:
                    parents = self._group_values(cell_groups, prefix + label)
                for parent in parents:
                    if face in (None, "SPOS", "SNEG"):
                        generated.append(parent)
                        continue
                    # Abaqus uses S<n> for solid faces and E<n> for
                    # shell/2D-element edges. Both map to the corresponding
                    # entry of FACES for the converted MED topology.
                    if (
                        len(face) < 2
                        or face[0] not in ("S", "E")
                        or not face[1:].isdigit()
                    ):
                        logger.warning(
                            "Unsupported Abaqus surface selector ignored: "
                            "%s (element type %s)",
                            face,
                            records[parent][0],
                        )
                        continue
                    etype, nodes = records[parent]
                    available = self._surface_topology_cache.get(etype)
                    if available is None:
                        med_name = CellsTypeConverter._abaqus_to_med.get(etype)
                        available = FACES.get(med_name, ())
                        self._surface_topology_cache[etype] = available
                    face_id = int(face[1:])
                    if not available or not (1 <= face_id <= len(available)):
                        logger.warning(
                            "Unsupported Abaqus surface selector ignored: " "%s for %s",
                            face,
                            etype,
                        )
                        continue
                    surface_type, local = available[face_id - 1]
                    face_nodes = tuple(nodes[i - 1] for i in local)
                    med_type = self._surface_med_type_cache.get(surface_type)
                    if med_type is None:
                        med_type = self._type_converter.external_to_medcoupling(
                            surface_type
                        )
                        self._surface_med_type_cache[surface_type] = med_type
                    self.add_cell(next_cell, med_type, face_nodes)
                    records[next_cell] = (surface_type, face_nodes)
                    generated.append(next_cell)
                    next_cell += 1
            cell_groups[name] = generated
