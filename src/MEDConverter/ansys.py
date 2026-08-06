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

"""ANSYS CDB mesh reader for MEDConverter.

Scope
-----
This reader intentionally keeps only data required to build a MED mesh:

* nodes and their global coordinates;
* finite-element type entities and topology-relevant KEYOPT values;
* element connectivities;
* node and element components/groups.

Material identifiers, real constants, sections, element coordinate systems,
section orientations, loads and analysis properties are deliberately ignored.

Supported CDB representations
-----------------------------
* classic blocked NBLOCK and EBLOCK;
* EBLOCK with the COMPACT label, when a current TYPE context is available;
* ET and ETBLOCK element-type declarations;
* unblocked N, EN, E and EMORE mesh commands;
* CMBLOCK components;
* NSEL, ESEL, CM, CMSEL, CMGRP, CMEDIT, CMDELE and ALLSEL.
"""

from collections import OrderedDict
from dataclasses import dataclass, field
import os.path as osp
import re
import time

from .logger import logger
from .MEDConverterMesh import MEDConverterMesh
from .cells import CellsTypeConverter
from .connectivity import ConnectivityRenumberer


@dataclass
class AnsysElementType:
    entity_id: int
    ansys_type: int
    # Only topology-relevant options are retained. At present this is used by
    # MESH200 KEYOPT(1). No physical/property KEYOPT is preserved.
    topology_options: dict = field(default_factory=dict)


@dataclass
class AnsysCell:
    id: int
    type_entity: int
    nodes: list


@dataclass
class AnsysComponent:
    name: str
    entity: str
    ids: list = field(default_factory=list)
    children: list = field(default_factory=list)


class MEDConverterAnsys(MEDConverterMesh):
    """Convert ANSYS CDB mesh topology and groups to MED."""

    # Number of geometrical nodes for variants carrying one trailing
    # orientation node. The orientation node itself is intentionally ignored.
    _SUPPORT_NODE_COUNT = {
        (16, 3): 2,
        (18, 3): 2,
        (188, 3): 2,
        (288, 3): 2,
        (189, 4): 3,
        (289, 4): 3,
    }

    @staticmethod
    def convert_ansys_to_med(filename_ansys, verbose=False):
        tic = time.perf_counter()
        converter = MEDConverterAnsys()
        converter.verbose = verbose
        converter.read_ansys_mesh(filename_ansys)
        converter.create_UMesh()
        logger.debug(
            "Mesh converted (in %0.4f seconds)",
            time.perf_counter() - tic,
        )
        return converter.umesh

    def __init__(self):
        super().__init__()
        self.ansysmesh = None
        self._warned = set()

    def _warn_once(self, key, message, *args):
        if key not in self._warned:
            logger.warning(message, *args)
            self._warned.add(key)

    @staticmethod
    def _split(line):
        return [field.strip() for field in line.strip().split(",")]

    @staticmethod
    def _to_int(value, default=None):
        value = value.strip()
        if not value:
            return default
        return int(float(value))

    @staticmethod
    def _fortran_float(value):
        value = value.strip().replace("D", "E").replace("d", "e")
        return float(value) if value else 0.0

    @staticmethod
    def _fixed_ints(line, width):
        raw = line.rstrip("\r\n")
        values = []
        for start in range(0, len(raw), width):
            token = raw[start : start + width].strip()
            if token:
                values.append(int(token))
        return values

    @staticmethod
    def _parse_int_format(line):
        match = re.search(r"(\d+)\s*[iI]\s*(\d+)", line)
        if not match:
            raise RuntimeError("Unsupported ANSYS integer format: %s" % line.strip())
        return int(match.group(1)), int(match.group(2))

    @staticmethod
    def _parse_node_format(line):
        integer = re.search(r"(\d+)\s*[iI]\s*(\d+)", line)
        real = re.search(r"(?:\d+\s*)?[eE]\s*(\d+)\.", line)
        if not integer or not real:
            raise RuntimeError("Unsupported ANSYS node format: %s" % line.strip())
        first_real = int(integer.group(1)) * int(integer.group(2))
        return first_real, int(real.group(1))

    @staticmethod
    def _node_coordinates(line, first_real, real_width):
        raw = line.rstrip("\r\n")
        return tuple(
            MEDConverterAnsys._fortran_float(
                raw[
                    first_real
                    + index * real_width : first_real
                    + (index + 1) * real_width
                ]
            )
            for index in range(3)
        )

    @staticmethod
    def _clean_group_name(name):
        return name.strip()

    def _read_nblock(self, stream, nodes):
        first_real = real_width = None
        while True:
            line = stream.readline()
            if not line:
                raise RuntimeError("Unexpected end of file inside ANSYS NBLOCK")
            stripped = line.lstrip()
            upper = stripped.upper()
            if stripped.startswith("("):
                first_real, real_width = self._parse_node_format(stripped)
                continue
            if stripped.startswith("-1") or upper.startswith("N,"):
                return
            if not stripped.strip() or stripped.startswith("!"):
                continue
            if first_real is None:
                raise RuntimeError("ANSYS NBLOCK data precedes its format line")
            node_id = int(stripped.split(None, 1)[0])
            if node_id in nodes:
                raise RuntimeError("Duplicate ANSYS node %s" % node_id)
            nodes[node_id] = self._node_coordinates(line, first_real, real_width)

    def _append_cell(self, cells, cell_ids, element_id, type_entity, nodes):
        if element_id in cell_ids:
            raise RuntimeError("Duplicate ANSYS element %s" % element_id)
        if type_entity is None:
            raise RuntimeError(
                "ANSYS element %s has no element TYPE context" % element_id
            )
        if not nodes:
            raise RuntimeError(
                "ANSYS element %s has an empty connectivity" % element_id
            )
        cells.append(AnsysCell(element_id, type_entity, list(nodes)))
        cell_ids.add(element_id)

    def _read_classic_eblock(self, stream, cells, cell_ids):
        fields_per_line = integer_width = None
        pending = None
        while True:
            line = stream.readline()
            if not line:
                raise RuntimeError("Unexpected end of file inside ANSYS EBLOCK")
            stripped = line.lstrip()
            if stripped.startswith("("):
                fields_per_line, integer_width = self._parse_int_format(stripped)
                continue
            if stripped.startswith("-1"):
                if pending is not None:
                    raise RuntimeError(
                        "Incomplete ANSYS element %s at end of EBLOCK" % pending["id"]
                    )
                return
            if not stripped.strip() or stripped.startswith("!"):
                continue
            if integer_width is None:
                raise RuntimeError("ANSYS EBLOCK data precedes its format line")
            values = self._fixed_ints(line, integer_width)
            if fields_per_line is not None and len(values) > fields_per_line:
                raise RuntimeError("Too many fields in ANSYS EBLOCK record")

            if pending is None:
                if len(values) < 11:
                    raise RuntimeError("Invalid classic ANSYS EBLOCK record")
                pending = {
                    "type": values[1],
                    "count": values[8],
                    "id": values[10],
                    "nodes": list(values[11:]),
                }
            else:
                pending["nodes"].extend(values)

            if len(pending["nodes"]) > pending["count"]:
                raise RuntimeError(
                    "ANSYS element %s has too many node entries" % pending["id"]
                )
            if len(pending["nodes"]) == pending["count"]:
                self._append_cell(
                    cells,
                    cell_ids,
                    pending["id"],
                    pending["type"],
                    pending["nodes"],
                )
                pending = None

    def _read_compact_eblock(self, stream, cells, cell_ids, current_type_entity):
        """Read EBLOCK,COMPACT records.

        A compact record contains the element identifier followed by node
        identifiers. Because it does not carry the classic attribute header,
        a current TYPE command must identify the element-type entity.
        """
        integer_width = None
        while True:
            line = stream.readline()
            if not line:
                raise RuntimeError("Unexpected end of file inside compact ANSYS EBLOCK")
            stripped = line.lstrip()
            if stripped.startswith("("):
                _, integer_width = self._parse_int_format(stripped)
                continue
            if stripped.startswith("-1"):
                return
            if not stripped.strip() or stripped.startswith("!"):
                continue
            if integer_width is None:
                # Accept whitespace-separated compact records as well.
                values = [int(value) for value in stripped.split()]
            else:
                values = self._fixed_ints(line, integer_width)
            if len(values) < 2:
                raise RuntimeError("Invalid compact ANSYS EBLOCK record")
            self._append_cell(
                cells, cell_ids, values[0], current_type_entity, values[1:]
            )

    def _read_etblock(self, stream, element_types):
        """Read modern ETBLOCK records including embedded KEYOPT values."""
        integer_width = None
        character_fields = 0
        character_width = None
        while True:
            line = stream.readline()
            if not line:
                raise RuntimeError("Unexpected end of file inside ETBLOCK")
            stripped = line.lstrip()
            if stripped.startswith("("):
                integer_match = re.search(r"(\d+)\s*[iI]\s*(\d+)", stripped)
                char_match = re.search(r"(\d+)\s*[aA]\s*(\d+)", stripped)
                if not integer_match:
                    raise RuntimeError(
                        "Unsupported ANSYS ETBLOCK format: %s" % stripped.strip()
                    )
                integer_width = int(integer_match.group(2))
                if char_match:
                    character_fields = int(char_match.group(1))
                    character_width = int(char_match.group(2))
                continue
            if stripped.startswith("-1"):
                return
            if not stripped.strip() or stripped.startswith("!"):
                continue
            if integer_width is None:
                raise RuntimeError("ANSYS ETBLOCK data precedes its format line")

            raw = line.rstrip("\r\n")
            entity_id = int(raw[0:integer_width])
            ansys_type = int(raw[integer_width : 2 * integer_width])
            descriptor = AnsysElementType(entity_id, ansys_type)

            # In the common (2i9,19a9) representation, each A-field contains
            # a numeric KEYOPT value. Only KEYOPT(1) is retained, because it
            # determines the MESH200 topology. Other properties are ignored.
            if character_fields and character_width:
                offset = 2 * integer_width
                options = []
                for index in range(character_fields):
                    token = raw[
                        offset
                        + index * character_width : offset
                        + (index + 1) * character_width
                    ].strip()
                    if token:
                        try:
                            options.append(int(float(token)))
                        except ValueError:
                            options.append(0)
                    else:
                        options.append(0)
                if options:
                    descriptor.topology_options[1] = options[0]
            element_types[entity_id] = descriptor

    def _read_cmblock(self, stream, command):
        fields = self._split(command)
        if len(fields) < 4:
            raise RuntimeError("Invalid ANSYS CMBLOCK command")
        name = fields[1]
        entity = fields[2].upper()

        # The declared CMBLOCK count is the number of encoded integer
        # fields, not the number of identifiers obtained after expanding
        # negative range terminators. Example: 42302, -42316 contains two
        # encoded fields but represents identifiers 42302 through 42316.
        expected_encoded = int(fields[3].split()[0])
        width = None
        encoded_values = []

        while len(encoded_values) < expected_encoded:
            line = stream.readline()
            if not line:
                raise RuntimeError("Unexpected end of file inside CMBLOCK")
            stripped = line.lstrip()
            if stripped.startswith("("):
                _, width = self._parse_int_format(stripped)
                continue
            if not stripped.strip() or stripped.startswith("!"):
                continue
            if width is None:
                raise RuntimeError("ANSYS CMBLOCK data precedes its format line")

            encoded_values.extend(self._fixed_ints(line, width))
            if len(encoded_values) > expected_encoded:
                raise RuntimeError(
                    "Too many encoded entries in ANSYS component %s" % name
                )

        values = []
        for encoded in encoded_values:
            if encoded >= 0:
                values.append(encoded)
                continue

            if not values:
                raise RuntimeError(
                    "Invalid compressed ANSYS component %s: "
                    "negative range terminator appears first" % name
                )

            previous = values[-1]
            range_end = -encoded
            if range_end < previous:
                raise RuntimeError(
                    "Invalid compressed ANSYS component %s: "
                    "decreasing range %s to %s" % (name, previous, range_end)
                )
            values.extend(range(previous + 1, range_end + 1))

        return AnsysComponent(name, entity, values)

    @staticmethod
    def _selection_range(fields):
        first = MEDConverterAnsys._to_int(fields[4], None) if len(fields) > 4 else None
        last = MEDConverterAnsys._to_int(fields[5], first) if len(fields) > 5 else first
        step = MEDConverterAnsys._to_int(fields[6], 1) if len(fields) > 6 else 1
        if step is None or step <= 0:
            raise RuntimeError("Selection increment must be positive")
        return first, last, step

    def _apply_id_selection(self, fields, universe, selected):
        operation = fields[1].upper() if len(fields) > 1 and fields[1] else "S"
        if operation == "ALL":
            return set(universe)
        if operation == "NONE":
            return set()
        if operation == "INVE":
            return set(universe) - selected
        if operation == "STAT":
            return selected
        if operation not in {"S", "R", "A", "U"}:
            self._warn_once(
                ("SELECT_OPERATION", operation),
                "Unsupported ANSYS selection operation ignored: %s",
                operation,
            )
            return selected
        first, last, step = self._selection_range(fields)
        if first is None:
            return selected
        matches = set(range(first, last + 1, step)) & set(universe)
        if operation == "S":
            return matches
        if operation == "R":
            return selected & matches
        if operation == "A":
            return selected | matches
        return selected - matches

    def _apply_nsel(self, fields, nodes, selected):
        item = fields[2].upper() if len(fields) > 2 and fields[2] else "NODE"
        if item in {"NODE", ""}:
            return self._apply_id_selection(fields, nodes.keys(), selected)
        if item == "LOC":
            operation = fields[1].upper() if len(fields) > 1 else "S"
            axis = fields[3].upper() if len(fields) > 3 else ""
            axis_index = {"X": 0, "Y": 1, "Z": 2}.get(axis)
            if axis_index is None:
                self._warn_once(
                    ("NSEL_AXIS", axis),
                    "Unsupported ANSYS NSEL LOC axis ignored: %s",
                    axis,
                )
                return selected
            first, last, _ = self._selection_range(fields)
            if first is None:
                return selected
            lower, upper = sorted((float(first), float(last)))
            tolerance = max(1.0, abs(lower), abs(upper)) * 1.0e-12
            matches = {
                node_id
                for node_id, coordinates in nodes.items()
                if lower - tolerance <= coordinates[axis_index] <= upper + tolerance
            }
            if operation == "S":
                return matches
            if operation == "R":
                return selected & matches
            if operation == "A":
                return selected | matches
            if operation == "U":
                return selected - matches
            if operation == "ALL":
                return set(nodes)
            if operation == "NONE":
                return set()
            if operation == "INVE":
                return set(nodes) - selected
            return selected
        self._warn_once(
            ("NSEL_ITEM", item),
            "Unsupported ANSYS NSEL item ignored: %s",
            item,
        )
        return selected

    def _apply_esel(self, fields, cells, selected):
        item = fields[2].upper() if len(fields) > 2 and fields[2] else "ELEM"
        if item in {"ELEM", ""}:
            return self._apply_id_selection(
                fields, (cell.id for cell in cells), selected
            )
        if item == "TYPE":
            operation = fields[1].upper() if len(fields) > 1 else "S"
            first, last, step = self._selection_range(fields)
            if first is None:
                return selected
            wanted = set(range(first, last + 1, step))
            matches = {cell.id for cell in cells if cell.type_entity in wanted}
            if operation == "S":
                return matches
            if operation == "R":
                return selected & matches
            if operation == "A":
                return selected | matches
            if operation == "U":
                return selected - matches
            if operation == "ALL":
                return {cell.id for cell in cells}
            if operation == "NONE":
                return set()
            if operation == "INVE":
                return {cell.id for cell in cells} - selected
            return selected
        self._warn_once(
            ("ESEL_ITEM", item),
            "Unsupported ANSYS ESEL item ignored because properties are not read: %s",
            item,
        )
        return selected

    @staticmethod
    def _component_ids(component, components, entity, visited=None):
        if component is None:
            return set()
        if visited is None:
            visited = set()
        key = component.name.upper()
        if key in visited:
            raise RuntimeError(
                "Cyclic ANSYS component assembly detected at %s" % component.name
            )
        visited = set(visited)
        visited.add(key)
        if component.entity == entity:
            return set(component.ids)
        if component.entity != "GROUP":
            return set()
        result = set()
        for child_name in component.children:
            child = components.get(child_name.upper())
            result.update(
                MEDConverterAnsys._component_ids(child, components, entity, visited)
            )
        return result

    def _apply_cmsel(self, fields, components, selected_nodes, selected_elements):
        operation = fields[1].upper() if len(fields) > 1 and fields[1] else "S"
        name = fields[2].upper() if len(fields) > 2 else ""
        component = components.get(name)
        if component is None:
            self._warn_once(
                ("CMSEL_UNKNOWN", name),
                "Unknown ANSYS component in CMSEL ignored: %s",
                name,
            )
            return selected_nodes, selected_elements

        node_ids = self._component_ids(component, components, "NODE")
        element_ids = self._component_ids(component, components, "ELEM")

        def apply(current, values):
            if operation == "S":
                return set(values)
            if operation == "A":
                return current | values
            if operation == "R":
                return current & values
            if operation == "U":
                return current - values
            return current

        return apply(selected_nodes, node_ids), apply(selected_elements, element_ids)

    def _support_nodes(self, ansys_type, raw_nodes):
        count = self._SUPPORT_NODE_COUNT.get(
            (ansys_type, len(raw_nodes)), len(raw_nodes)
        )
        # Repeated positions encode degenerate solid/shell topologies.
        return list(dict.fromkeys(raw_nodes[:count]))

    @staticmethod
    def _topology_key(descriptor, support_nodes):
        if descriptor.ansys_type == 200:
            option = descriptor.topology_options.get(1)
            if option is None:
                raise RuntimeError(
                    "MESH200 type entity %s has no topology KEYOPT(1)"
                    % descriptor.entity_id
                )
            return "200_%s_%s" % (len(support_nodes), option)
        return "%s_%s" % (descriptor.ansys_type, len(support_nodes))

    def read_ansys_mesh(self, filename):
        self._reset_structures()
        self.space_dim = 3
        self.mesh_name = osp.splitext(osp.basename(filename))[0]

        nodes = OrderedDict()
        cells = []
        cell_ids = set()
        element_types = {}
        components = OrderedDict()
        selected_nodes = set()
        selected_elements = set()
        current_type_entity = None
        pending_unblocked = None
        declared_cells = 0
        classic_cells_read = 0

        with open(filename, "r", encoding=self._get_file_encoding(filename)) as stream:
            for line in stream:
                stripped = line.strip()
                upper = stripped.upper()
                if not stripped or stripped.startswith("!"):
                    continue

                # E starts an unblocked element and EMORE extends it. Any
                # other command closes the pending element before that command
                # is interpreted, so selections/components can see it.
                if pending_unblocked is not None and not upper.startswith("EMORE,"):
                    self._append_cell(
                        cells,
                        cell_ids,
                        pending_unblocked[0],
                        pending_unblocked[1],
                        pending_unblocked[2],
                    )
                    pending_unblocked = None

                if upper.startswith("/TITLE"):
                    fields = self._split(line)
                    if len(fields) > 1 and fields[1]:
                        self.mesh_name = fields[1]
                elif upper.startswith("NBLOCK"):
                    self._read_nblock(stream, nodes)
                elif upper.startswith("EBLOCK"):
                    fields = self._split(line)
                    if "COMPACT" in upper:
                        self._read_compact_eblock(
                            stream, cells, cell_ids, current_type_entity
                        )
                    else:
                        if len(fields) <= 4 or not fields[4]:
                            raise RuntimeError(
                                "Invalid classic ANSYS EBLOCK header: %s" % stripped
                            )
                        declared_cells += int(fields[4])
                        before = len(cells)
                        self._read_classic_eblock(stream, cells, cell_ids)
                        classic_cells_read += len(cells) - before
                elif upper.startswith("ETBLOCK"):
                    self._read_etblock(stream, element_types)
                elif upper.startswith("ET,"):
                    fields = self._split(line)
                    entity_id = int(fields[1])
                    element_types[entity_id] = AnsysElementType(
                        entity_id, int(fields[2])
                    )
                elif upper.startswith("KEYOP"):
                    fields = self._split(line)
                    if len(fields) >= 4:
                        entity_id = int(fields[1])
                        option_id = int(fields[2])
                        value = int(float(fields[3]))
                        # Only MESH200 KEYOPT(1) influences topology.
                        if option_id == 1:
                            descriptor = element_types.get(entity_id)
                            if descriptor is None:
                                descriptor = AnsysElementType(entity_id, 200)
                                element_types[entity_id] = descriptor
                            descriptor.topology_options[1] = value
                elif upper.startswith("TYPE,"):
                    fields = self._split(line)
                    current_type_entity = int(fields[1])
                elif upper.startswith("N,"):
                    fields = self._split(line)

                    # Classic NBLOCK exports may be followed by one or more
                    # control records such as:
                    #   N,UNBL,LOC,-1
                    #   N,R5.1,LOC,-1
                    # They are not unblocked node definitions. An actual
                    # unblocked N command has a numeric node identifier in
                    # its second field.
                    try:
                        node_id = int(fields[1])
                    except (IndexError, ValueError):
                        continue

                    coords = tuple(
                        (
                            self._fortran_float(fields[index])
                            if index < len(fields)
                            else 0.0
                        )
                        for index in range(2, 5)
                    )
                    if node_id in nodes:
                        raise RuntimeError("Duplicate ANSYS node %s" % node_id)
                    nodes[node_id] = coords
                elif upper.startswith("EN,"):
                    fields = self._split(line)

                    # Un EBLOCK classique peut être suivi de :
                    #
                    # EN,UNBL,ATTR,-1
                    #
                    # Ce n'est pas une définition d'élément.
                    try:
                        element_id = int(fields[1])
                    except (IndexError, ValueError):
                        continue

                    connectivity = []

                    for value in fields[2:]:
                        if not value:
                            continue

                    try:
                        connectivity.append(int(value))
                    except ValueError as exc:
                        raise RuntimeError(
                            "Invalid ANSYS EN connectivity "
                            "for element %s: %s"
                            % (
                                element_id,
                                value,
                            )
                        ) from exc

                    self._append_cell(
                        cells,
                        cell_ids,
                        element_id,
                        current_type_entity,
                        connectivity,
                    )
                elif upper.startswith("E,"):
                    fields = self._split(line)
                    node_values = [int(value) for value in fields[1:] if value]
                    if not node_values:
                        raise RuntimeError("Invalid ANSYS E command")
                    next_id = max(cell_ids, default=0) + 1
                    pending_unblocked = [next_id, current_type_entity, node_values]
                elif upper.startswith("EMORE,"):
                    if pending_unblocked is None:
                        raise RuntimeError("ANSYS EMORE without preceding E")
                    fields = self._split(line)
                    pending_unblocked[2].extend(
                        int(value) for value in fields[1:] if value
                    )
                elif upper.startswith("CMBLOCK"):
                    component = self._read_cmblock(stream, line)
                    components[component.name.upper()] = component
                elif upper.startswith("NSEL,"):
                    selected_nodes = self._apply_nsel(
                        self._split(line), nodes, selected_nodes
                    )
                elif upper.startswith("ESEL,"):
                    selected_elements = self._apply_esel(
                        self._split(line), cells, selected_elements
                    )
                elif upper.startswith("ALLSEL"):
                    selected_nodes = set(nodes)
                    selected_elements = {cell.id for cell in cells}
                elif upper.startswith("CM,"):
                    fields = self._split(line)
                    if len(fields) >= 3:
                        name, entity = fields[1], fields[2].upper()
                        if entity == "NODE":
                            ids = sorted(selected_nodes)
                        elif entity == "ELEM":
                            ids = sorted(selected_elements)
                        else:
                            self._warn_once(
                                ("CM_ENTITY", entity),
                                "Unsupported ANSYS CM entity ignored: %s",
                                entity,
                            )
                            continue
                        components[name.upper()] = AnsysComponent(name, entity, ids)
                elif upper.startswith("CMSEL,"):
                    selected_nodes, selected_elements = self._apply_cmsel(
                        self._split(line),
                        components,
                        selected_nodes,
                        selected_elements,
                    )
                elif upper.startswith("CMGRP,"):
                    fields = self._split(line)
                    if len(fields) >= 3:
                        name = fields[1]
                        children = [value.upper() for value in fields[2:] if value]
                        components[name.upper()] = AnsysComponent(
                            name, "GROUP", children=children
                        )
                elif upper.startswith("CMEDIT,"):
                    fields = self._split(line)
                    if len(fields) >= 4:
                        name = fields[1].upper()
                        operation = fields[2].upper()
                        children = [value.upper() for value in fields[3:] if value]
                        component = components.get(name)
                        if component is None or component.entity != "GROUP":
                            self._warn_once(
                                ("CMEDIT_UNKNOWN", name),
                                "Unknown ANSYS component assembly ignored: %s",
                                name,
                            )
                        elif operation in {"ADD", "A"}:
                            component.children.extend(
                                child
                                for child in children
                                if child not in component.children
                            )
                        elif operation in {"DELE", "DELETE", "D"}:
                            component.children = [
                                child
                                for child in component.children
                                if child not in children
                            ]
                        else:
                            self._warn_once(
                                ("CMEDIT_OPERATION", operation),
                                "Unsupported ANSYS CMEDIT operation ignored: %s",
                                operation,
                            )
                elif upper.startswith("CMDELE,"):
                    fields = self._split(line)
                    if len(fields) > 1:
                        components.pop(fields[1].upper(), None)

        if pending_unblocked is not None:
            self._append_cell(
                cells,
                cell_ids,
                pending_unblocked[0],
                pending_unblocked[1],
                pending_unblocked[2],
            )
        if declared_cells and declared_cells != classic_cells_read:
            raise RuntimeError(
                "Classic ANSYS EBLOCK declared %d elements, but %d classic "
                "elements were read" % (declared_cells, classic_cells_read)
            )

        for node_id, coordinates in nodes.items():
            self.add_node(node_id, coordinates)

        type_converter = CellsTypeConverter("ANSYS")
        connectivity_converter = ConnectivityRenumberer("ANSYS")
        groups_by_type = OrderedDict()

        for cell in cells:
            descriptor = element_types.get(cell.type_entity)
            if descriptor is None:
                raise RuntimeError(
                    "ANSYS element %s references undefined type entity %s"
                    % (cell.id, cell.type_entity)
                )
            support_nodes = self._support_nodes(descriptor.ansys_type, cell.nodes)
            missing = [node for node in support_nodes if node not in nodes]
            if missing:
                raise RuntimeError(
                    "ANSYS element %s references missing nodes %s" % (cell.id, missing)
                )
            topology_key = self._topology_key(descriptor, support_nodes)
            try:
                med_type = type_converter.external_to_medcoupling(topology_key)
            except (KeyError, RuntimeError) as exc:
                raise RuntimeError(
                    "Unsupported ANSYS element topology %s for element %s"
                    % (topology_key, cell.id)
                ) from exc
            med_nodes = connectivity_converter.external_to_medcoupling(
                med_type, support_nodes
            )
            self.add_cell(cell.id, med_type, med_nodes)
            groups_by_type.setdefault(str(descriptor.ansys_type), []).append(cell.id)

        # Flatten component assemblies by entity type for MED groups.
        materialized = list(components.values())
        for component in materialized:
            if component.entity == "GROUP":
                node_ids = sorted(self._component_ids(component, components, "NODE"))
                element_ids = sorted(self._component_ids(component, components, "ELEM"))
                if node_ids:
                    self.add_group_nodes(component.name, node_ids)
                if element_ids:
                    self.add_group_cells(component.name, element_ids)
                continue

            ids = list(dict.fromkeys(component.ids))
            if component.entity == "NODE":
                missing = [value for value in ids if value not in nodes]
                if missing:
                    raise RuntimeError(
                        "ANSYS node component %s references missing nodes %s"
                        % (component.name, missing[:10])
                    )
                if ids:
                    self.add_group_nodes(component.name, ids)
            elif component.entity == "ELEM":
                missing = [value for value in ids if value not in cell_ids]
                if missing:
                    raise RuntimeError(
                        "ANSYS element component %s references missing elements %s"
                        % (component.name, missing[:10])
                    )
                if ids:
                    self.add_group_cells(component.name, ids)

        for ansys_type, values in groups_by_type.items():
            self.add_group_cells("Grp_FE_" + ansys_type, values)
