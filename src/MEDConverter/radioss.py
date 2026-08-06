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

"""Native Radioss Starter mesh reader for MEDConverter.

Supported element blocks:
    /BRICK, /BRIC20, /TETRA4, /TETRA10, /PENTA6,
    /SHELL, /SH3N, /QUAD, /TRIA, /BEAM, /TRUSS, /SPRING.

Includes are resolved recursively relative to their parent file. Cyclic
includes are prevented with a visited-file set. Group dependency resolution is
strictly bounded and fails explicitly on cycles or unsupported selectors.
"""

from collections import OrderedDict
from pathlib import Path
import re
import time

from .logger import logger
from .MEDConverterMesh import MEDConverterMesh
from .cells import CellsTypeConverter
from .connectivity import ConnectivityRenumberer


class MEDConverterRadioss(MEDConverterMesh):
    # keyword: (topological node count, integer node fields in the card)
    _element_specs = OrderedDict(
        (
            ("BRICK", (8, 8)),
            ("BRIC20", (20, 20)),
            ("TETRA4", (4, 4)),
            ("TETRA10", (10, 10)),
            ("PENTA6", (6, 6)),
            ("SHEL16", (16, 16)),
            ("SHELL", (4, 4)),
            ("SH3N", (3, 3)),
            ("QUAD", (4, 4)),
            ("TRIA", (3, 3)),
            # BEAM and SPRING may contain optional orientation/auxiliary
            # nodes. Only the two topological end nodes are required for the
            # MED SEG2 connectivity; any additional fields are ignored.
            ("BEAM", (2, 2)),
            ("TRUSS", (2, 2)),
            ("SPRING", (2, 2)),
        )
    )

    _element_header = re.compile(
        r"^/(%s)/(\d+)(?:/\d+)?\s*$" % "|".join(_element_specs),
        re.I,
    )
    _part_header = re.compile(
        r"^/PART/(\d+)(?:/\d+)?\s*$",
        re.I,
    )
    _property_header = re.compile(
        r"^/PROP/([A-Z0-9_]+)/(\d+)(?:/\d+)?\s*$",
        re.I,
    )

    # Canonical property family and zero-based index of the finite-element
    # formulation flag in the first numerical line after the property title.
    #
    # Examples:
    #   /PROP/SHELL or /PROP/TYPE1  -> SHELL, Ishell at index 0
    #   /PROP/SOLID or /PROP/TYPE14 -> SOLID, Isolid at index 0
    #
    # Only property cards whose formulation position is explicitly known are
    # listed here. Unknown cards retain their raw keyword and do not append a
    # guessed formulation value.
    _property_formulation_specs = {
        "SHELL": ("SHELL", 0),
        "TYPE1": ("SHELL", 0),
        "SOLID": ("SOLID", 0),
        "TYPE14": ("SOLID", 0),
    }
    _group_header = re.compile(
        r"^/(GRNOD|GRPART|GRBEAM|GRBRIC|GRQUAD|GRSH3N|"
        r"GRSHEL|GRSPRI|GRTRIA|GRTRUS)/"
        r"([A-Z0-9_]+)/(\d+)(?:/\d+)?\s*$",
        re.I,
    )
    _surface_part_ext_header = re.compile(
        r"^/SURF/PART/EXT/(\d+)(?:/\d+)?\s*$",
        re.I,
    )
    _surface_group_header = re.compile(
        r"^/SURF/(GRSHEL|GRSH3N|SURF|PART|MAT|PROP)/(\d+)" r"(?:/\d+)?\s*$",
        re.I,
    )
    _surface_seg_header = re.compile(
        r"^/SURF/SEG/(\d+)(?:/\d+)?\s*$",
        re.I,
    )
    _admas_header = re.compile(
        r"^/ADMAS/(\d+)/(\d+)(?:/\d+)?\s*$",
        re.I,
    )

    # Surface-cell specifications in native Radioss connectivity. Each entry
    # is (MED type name, full face-node indices, corner-node indices). Corner
    # nodes identify coincident faces; full indices retain quadratic nodes.
    _surface_face_specs = {
        "BRICK": (
            ("QUAD4", (0, 3, 2, 1), (0, 1, 2, 3)),
            ("QUAD4", (4, 5, 6, 7), (4, 5, 6, 7)),
            ("QUAD4", (0, 1, 5, 4), (0, 1, 5, 4)),
            ("QUAD4", (1, 2, 6, 5), (1, 2, 6, 5)),
            ("QUAD4", (2, 3, 7, 6), (2, 3, 7, 6)),
            ("QUAD4", (3, 0, 4, 7), (3, 0, 4, 7)),
        ),
        "BRIC20": (
            ("QUAD8", (0, 3, 2, 1, 11, 10, 9, 8), (0, 1, 2, 3)),
            ("QUAD8", (4, 5, 6, 7, 12, 13, 14, 15), (4, 5, 6, 7)),
            ("QUAD8", (0, 1, 5, 4, 8, 17, 12, 16), (0, 1, 5, 4)),
            ("QUAD8", (1, 2, 6, 5, 9, 18, 13, 17), (1, 2, 6, 5)),
            ("QUAD8", (2, 3, 7, 6, 10, 19, 14, 18), (2, 3, 7, 6)),
            ("QUAD8", (3, 0, 4, 7, 11, 16, 15, 19), (3, 0, 4, 7)),
        ),
        "TETRA4": (
            ("TRI3", (0, 2, 1), (0, 1, 2)),
            ("TRI3", (0, 1, 3), (0, 1, 3)),
            ("TRI3", (1, 2, 3), (1, 2, 3)),
            ("TRI3", (2, 0, 3), (2, 0, 3)),
        ),
        "TETRA10": (
            ("TRI6", (0, 2, 1, 6, 5, 4), (0, 1, 2)),
            ("TRI6", (0, 1, 3, 4, 8, 7), (0, 1, 3)),
            ("TRI6", (1, 2, 3, 5, 9, 8), (1, 2, 3)),
            ("TRI6", (2, 0, 3, 6, 7, 9), (2, 0, 3)),
        ),
        "SHEL16": (
            ("QUAD4", (0, 3, 2, 1), (0, 1, 2, 3)),
            ("QUAD4", (4, 5, 6, 7), (4, 5, 6, 7)),
            ("QUAD4", (0, 1, 5, 4), (0, 1, 5, 4)),
            ("QUAD4", (1, 2, 6, 5), (1, 2, 6, 5)),
            ("QUAD4", (2, 3, 7, 6), (2, 3, 7, 6)),
            ("QUAD4", (3, 0, 4, 7), (3, 0, 4, 7)),
        ),
        "PENTA6": (
            ("TRI3", (0, 2, 1), (0, 1, 2)),
            ("TRI3", (3, 4, 5), (3, 4, 5)),
            ("QUAD4", (0, 1, 4, 3), (0, 1, 4, 3)),
            ("QUAD4", (1, 2, 5, 4), (1, 2, 5, 4)),
            ("QUAD4", (2, 0, 3, 5), (2, 0, 3, 5)),
        ),
    }
    _boundary_card_header = re.compile(
        r"^/(BCS|CLOAD|GRAV|IMPDISP|INIVEL/TRA|INTER/TYPE10|"
        r"INTER/TYPE24|RBODY|MONVOL/[A-Z0-9_]+)/(\d+)"
        r"(?:/\d+)?\s*$",
        re.I,
    )

    _include = re.compile(
        r"^#include\s+(.+?)\s*$",
        re.I,
    )

    _group_family = {
        "GRBEAM": "BEAM",
        "GRBRIC": "SOLID",
        "GRQUAD": "QUAD",
        "GRSH3N": "SH3N",
        "GRSHEL": "SHELL",
        "GRSPRI": "SPRING",
        "GRTRIA": "TRIA",
        "GRTRUS": "TRUSS",
    }

    _direct_selector = {
        "GRBEAM": "BEAM",
        "GRBRIC": "BRIC",
        "GRQUAD": "QUAD",
        "GRSH3N": "SH3N",
        "GRSHEL": "SHEL",
        "GRSPRI": "SPRI",
        "GRTRIA": "TRIA",
        "GRTRUS": "TRUS",
        "GRPART": "PART",
        "GRNOD": "NODE",
    }

    def __init__(self):
        super(MEDConverterRadioss, self).__init__()

    @staticmethod
    def convert_radioss_to_med(filename_radioss, verbose=False):
        tic = time.perf_counter()
        converter = MEDConverterRadioss()
        converter.verbose = verbose
        converter.read_radioss_mesh(filename_radioss)
        converter.create_UMesh()
        logger.debug("Mesh converted (in %0.4f seconds)" % (time.perf_counter() - tic))
        return converter.umesh

    @staticmethod
    def _read(path):
        for encoding in ("utf8", "latin_1", "cp437"):
            try:
                return Path(path).read_text(encoding=encoding).splitlines()
            except UnicodeDecodeError:
                continue
        raise RuntimeError("Unable to determine encoding of %s" % path)

    @staticmethod
    def _data(lines, start):
        for line in lines[start + 1 :]:
            value = line.strip()
            if value.startswith("/"):
                break
            if value.lower().startswith("#enddata"):
                break
            if value and not value.startswith("#"):
                yield line

    @classmethod
    def _resolve_include(cls, parent, name):
        parent_directory = Path(parent).resolve().parent
        name = name.strip().strip('"').strip("'")

        candidates = [parent_directory / name]
        stem = Path(name).stem
        candidates.extend(
            parent_directory / (stem + extension)
            for extension in (".txt", ".inc", ".rad", "")
        )

        for candidate in candidates:
            if candidate.is_file():
                return candidate.resolve()

        raise FileNotFoundError(
            "Radioss include not found: %s (from %s)" % (name, parent_directory)
        )

    @classmethod
    def _deck(cls, filename_radioss):
        output = []
        visited = set()
        active_stack = []

        def visit(path):
            path = Path(path).resolve()

            if path in active_stack:
                cycle = active_stack[active_stack.index(path) :] + [path]
                raise RuntimeError(
                    "Cyclic Radioss include dependency: %s"
                    % " -> ".join(map(str, cycle))
                )

            if path in visited:
                return

            active_stack.append(path)
            visited.add(path)
            lines = cls._read(path)
            output.extend(lines)

            for line in lines:
                match = cls._include.match(line.strip())
                if match:
                    include_path = cls._resolve_include(
                        path,
                        match.group(1),
                    )
                    visit(include_path)

            active_stack.pop()

        visit(filename_radioss)
        return output

    @staticmethod
    def _node(line):
        line = line.ljust(70)
        try:
            node_id = int(line[:10])
            coordinates = tuple(
                float(line[index : index + 20]) for index in (10, 30, 50)
            )
        except ValueError as exc:
            raise RuntimeError("Invalid Radioss /NODE line: %s" % line) from exc

        return node_id, coordinates

    @staticmethod
    def _title_and_ints(lines, index):
        """Read a Radioss group title followed by integer identifiers.

        Radioss decks encountered in practice use both forms::

            /GRNOD/NODE/21
            #title
            BCS_1
            365 366 767

        and::

            /GRNOD/NODE/21
            BCS_1
            365 366 767

        The first non-comment, non-integer line is therefore accepted as the
        implicit title. Once a title or identifiers have been read, any other
        non-integer data line is reported instead of being silently ignored.
        """
        title = ""
        values = []
        title_next = False

        for line_number in range(index + 1, len(lines)):
            item = lines[line_number].strip()

            if item.startswith("/"):
                break
            if item.lower().startswith("#enddata"):
                break
            if not item:
                continue

            if item.lower() == "#title":
                title_next = True
                continue

            if item.startswith("#"):
                continue

            if title_next:
                title = item
                title_next = False
                continue

            fields = item.split()
            try:
                identifiers = [int(field) for field in fields]
            except ValueError as exc:
                if not title and not values:
                    # Implicit Radioss title, for example ``BCS_1``.
                    title = item
                    continue

                raise RuntimeError(
                    "Invalid non-integer line in Radioss group at line %d: %s"
                    % (line_number + 1, item)
                ) from exc

            values.extend(identifiers)

        if title_next:
            raise RuntimeError(
                "Missing title after #title in Radioss group at line %d" % (index + 1)
            )

        return title, values

    @staticmethod
    def _generated(values, increment=False):
        result = []
        stride = 3 if increment else 2

        if len(values) % stride:
            raise RuntimeError("Invalid generated group definition: %s" % values)

        for index in range(0, len(values), stride):
            first = values[index]
            last = values[index + 1]
            step = values[index + 2] if increment else 1

            if step == 0:
                raise RuntimeError("A generated group cannot use a zero increment")

            stop = last + (1 if step > 0 else -1)
            result.extend(range(first, stop, step))

        return result

    @staticmethod
    def _element_records(keyword, rows, input_count):
        if keyword in ("BRIC20", "TETRA10", "SHEL16"):
            tokens = [token for row in rows for token in row.split()]
            width = 1 + input_count

            if len(tokens) % width:
                raise RuntimeError("Invalid multiline Radioss /%s block" % keyword)

            return [
                tokens[index : index + width] for index in range(0, len(tokens), width)
            ]

        # For variable-width records such as /BEAM and /SPRING, only the
        # required mesh fields are tested here. Optional orientation nodes,
        # auxiliary nodes, vectors, skew IDs, thicknesses, and other trailing
        # values remain in the row but are deliberately ignored when the
        # topological connectivity is extracted.
        return [row.split() for row in rows if len(row.split()) >= 1 + input_count]

    @staticmethod
    def _apply_signed_identifiers(values):
        """Apply Radioss add/remove semantics in linear time.

        Positive identifiers are inserted once while preserving their first
        insertion order. A negative identifier removes its positive
        counterpart. ``dict`` preserves insertion order and provides average
        O(1) membership, insertion and deletion, unlike the former list-based
        implementation whose duplicate checks were O(n^2).
        """
        selected = {}

        for value in values:
            if value < 0:
                selected.pop(-value, None)
            else:
                selected.setdefault(value, None)

        return list(selected)

    @staticmethod
    def _card_content(lines, index):
        """Return non-empty, non-comment content until the next card."""
        content = []
        for line in lines[index + 1 :]:
            item = line.strip()
            if item.startswith("/") or item.lower().startswith("#enddata"):
                break
            if item and not item.startswith("#"):
                content.append(item)
        return content

    @staticmethod
    def _integer_at(fields, position):
        if len(fields) <= position:
            return None
        try:
            return int(fields[position])
        except ValueError:
            return None

    @classmethod
    def _boundary_references_from_card(cls, header, lines, index):
        """Extract GRNOD/SURF references from supported BC/load cards."""
        match = cls._boundary_card_header.match(header)
        if not match:
            return []

        keyword = match.group(1).upper()
        card_id = int(match.group(2))
        content = cls._card_content(lines, index)
        references = []

        def add(entity, entity_id):
            if entity_id not in (None, 0):
                references.append((keyword, card_id, entity, entity_id))

        if keyword == "BCS":
            for item in content:
                fields = item.split()
                group_id = cls._integer_at(fields, 3)
                if group_id is not None:
                    add("GRNOD", group_id)
                    break

        elif keyword in ("CLOAD", "GRAV", "IMPDISP"):
            for item in content:
                fields = item.split()
                group_id = cls._integer_at(fields, 4)
                if group_id is not None:
                    add("GRNOD", group_id)
                    break

        elif keyword == "INIVEL/TRA":
            for item in content:
                fields = item.split()
                group_id = cls._integer_at(fields, 3)
                if group_id is not None:
                    add("GRNOD", group_id)
                    break

        elif keyword == "INTER/TYPE10":
            for item in content:
                fields = item.split()
                group_id = cls._integer_at(fields, 0)
                surface_id = cls._integer_at(fields, 1)
                if group_id is not None and surface_id is not None:
                    add("GRNOD", group_id)
                    add("SURF", surface_id)
                    break

        elif keyword == "INTER/TYPE24":
            # First numerical record: surf_id1, surf_id2, ...
            for item in content:
                fields = item.split()
                surface_1 = cls._integer_at(fields, 0)
                surface_2 = cls._integer_at(fields, 1)
                if surface_1 is not None and surface_2 is not None:
                    add("SURF", surface_1)
                    add("SURF", surface_2)
                    break

            # Optional Grnod_id is on the physical line following its comment.
            for line_number in range(index + 1, len(lines) - 1):
                item = lines[line_number].strip()
                if item.startswith("/"):
                    break
                if "grnod_id" in item.lower():
                    candidate = lines[line_number + 1][:10].strip()
                    if candidate:
                        try:
                            add("GRNOD", int(candidate))
                        except ValueError:
                            pass
                    break

        elif keyword == "RBODY":
            for item in content:
                fields = item.split()
                group_id = cls._integer_at(fields, 5)
                surface_id = cls._integer_at(fields, 8)
                if group_id is not None:
                    add("GRNOD", group_id)
                    add("SURF", surface_id)
                    break

        elif keyword.startswith("MONVOL/"):
            for item in content:
                fields = item.split()
                surface_id = cls._integer_at(fields, 0)
                if surface_id is not None:
                    add("SURF", surface_id)
                    break

        return references

    @staticmethod
    def _normalize_element_topology(keyword, identifiers):
        """Return the actual MED-compatible topology of a Radioss element."""
        if keyword == "BRICK":
            unique_nodes = list(dict.fromkeys(i for i in identifiers if i != 0))
            derived = {4: "TETRA4", 5: "PYRA5", 6: "PENTA6", 8: "BRICK"}
            actual_keyword = derived.get(len(unique_nodes))
            if actual_keyword is None:
                raise RuntimeError(
                    "Unsupported degenerated /BRICK with %d distinct nodes: %s"
                    % (len(unique_nodes), identifiers)
                )
            return actual_keyword, tuple(unique_nodes)

        if keyword == "SHEL16":
            # MEDCoupling has no 16-node thick-shell cell. Preserve the exact
            # external geometry with the eight corner nodes as HEXA8. The eight
            # additional interpolation nodes are not representable in MED.
            corners = tuple(i for i in identifiers[:8] if i != 0)
            if len(corners) != 8:
                raise RuntimeError("Invalid /SHEL16 corner connectivity")
            return "SHEL16", corners

        topology = tuple(i for i in identifiers if i != 0)
        return keyword, topology

    def read_radioss_mesh(self, filename_radioss):
        source = Path(filename_radioss).resolve()
        self.filename = str(source)
        self.mesh_name = source.stem.replace("_0000", "")
        self.space_dim = 3

        lines = self._deck(source)
        nodes = OrderedDict()
        elements = []
        part_names = {}
        part_property_ids = {}
        part_material_ids = {}
        property_signatures = {}
        admas_definitions = []
        raw_groups = []
        raw_surfaces = {}
        boundary_references = []

        for index, line in enumerate(lines):
            value = line.strip()

            boundary_references.extend(
                self._boundary_references_from_card(
                    value,
                    lines,
                    index,
                )
            )

            if value.upper() == "/NODE":
                for row in self._data(lines, index):
                    node_id, coordinates = self._node(row)
                    if node_id in nodes:
                        raise RuntimeError("Duplicate node ID %d" % node_id)
                    nodes[node_id] = coordinates
                continue

            match = self._element_header.match(value)
            if match:
                keyword = match.group(1).upper()
                part_id = int(match.group(2))
                topology_count, input_count = self._element_specs[keyword]
                rows = list(self._data(lines, index))
                records = self._element_records(
                    keyword,
                    rows,
                    input_count,
                )

                for fields in records:
                    try:
                        cell_id = int(fields[0])
                        identifiers = tuple(map(int, fields[1 : 1 + input_count]))
                    except ValueError as exc:
                        raise RuntimeError(
                            "Invalid Radioss /%s mesh fields: %s" % (keyword, fields)
                        ) from exc

                    actual_keyword, topology = self._normalize_element_topology(
                        keyword,
                        identifiers[:topology_count],
                    )
                    elements.append((cell_id, actual_keyword, part_id, topology))
                continue

            match = self._property_header.match(value)
            if match:
                raw_property_type = match.group(1).upper()
                property_id = int(match.group(2))
                property_data = list(self._data(lines, index))

                family, formulation_index = self._property_formulation_specs.get(
                    raw_property_type,
                    (raw_property_type, None),
                )

                formulation = None
                if formulation_index is not None:
                    # property_data[0] is the title. The next non-comment data
                    # line starts with Ishell for SHELL/TYPE1 and Isolid for
                    # SOLID/TYPE14.
                    if len(property_data) < 2:
                        raise RuntimeError(
                            "Property /PROP/%s/%d has no formulation data"
                            % (raw_property_type, property_id)
                        )

                    fields = property_data[1].split()
                    if len(fields) <= formulation_index:
                        raise RuntimeError(
                            "Cannot read formulation from /PROP/%s/%d"
                            % (raw_property_type, property_id)
                        )

                    try:
                        formulation = int(fields[formulation_index])
                    except ValueError as exc:
                        raise RuntimeError(
                            "Invalid formulation in /PROP/%s/%d: %s"
                            % (
                                raw_property_type,
                                property_id,
                                property_data[1],
                            )
                        ) from exc

                signature = (family, formulation)
                previous_signature = property_signatures.get(property_id)
                if previous_signature is not None and previous_signature != signature:
                    raise RuntimeError(
                        "Property ID %d is declared with two signatures: "
                        "%s and %s"
                        % (
                            property_id,
                            previous_signature,
                            signature,
                        )
                    )

                property_signatures[property_id] = signature
                continue

            match = self._part_header.match(value)
            if match:
                part_id = int(match.group(1))
                data = list(self._data(lines, index))

                if data:
                    part_names[part_id] = data[0].strip()

                if len(data) >= 2:
                    fields = data[1].split()
                    if fields:
                        try:
                            part_property_ids[part_id] = int(fields[0])
                            if len(fields) >= 2:
                                part_material_ids[part_id] = int(fields[1])
                        except ValueError as exc:
                            raise RuntimeError(
                                "Invalid prop_ID in /PART/%d: %s" % (part_id, data[1])
                            ) from exc
                continue

            match = self._admas_header.match(value)
            if match:
                admas_type = int(match.group(1))
                admas_id = int(match.group(2))
                data = list(self._data(lines, index))
                title = data[0].strip() if data else "ADMAS_%d" % admas_id
                node_ids = []
                group_ids = []
                surface_ids = []
                if admas_type in (0, 1) and len(data) >= 2:
                    fields = data[1].split()
                    if len(fields) >= 2:
                        group_ids.append(int(fields[1]))
                elif admas_type == 2 and len(data) >= 2:
                    fields = data[1].split()
                    if len(fields) >= 2:
                        surface_ids.append(int(fields[1]))
                elif admas_type == 5:
                    for row in data[1:]:
                        fields = row.split()
                        if len(fields) >= 2:
                            node_ids.append(int(fields[1]))
                admas_definitions.append(
                    (admas_id, title, node_ids, group_ids, surface_ids)
                )
                continue

            match = self._surface_seg_header.match(value)
            if match:
                surface_id = int(match.group(1))
                data = list(self._data(lines, index))
                title = data[0].strip() if data else "SURF_%d" % surface_id
                segments = []
                for row in data[1:]:
                    fields = list(map(int, row.split()))
                    if len(fields) >= 4:
                        nodes_on_segment = fields[1:5]
                        nodes_on_segment = [n for n in nodes_on_segment if n != 0]
                        if len(nodes_on_segment) in (3, 4):
                            segments.append(tuple(nodes_on_segment))
                raw_surfaces[surface_id] = ("SEG", title, segments)
                continue

            match = self._surface_group_header.match(value)
            if match:
                selector = match.group(1).upper()
                surface_id = int(match.group(2))
                title, item_ids = self._title_and_ints(lines, index)
                raw_surfaces[surface_id] = (
                    selector,
                    title or "SURF_%d" % surface_id,
                    item_ids,
                )
                continue

            match = self._surface_part_ext_header.match(value)
            if match:
                surface_id = int(match.group(1))
                title, part_ids = self._title_and_ints(lines, index)
                raw_surfaces[surface_id] = (
                    "PART/EXT",
                    title or "SURF_%d" % surface_id,
                    part_ids,
                )
                continue

            match = self._group_header.match(value)
            if match:
                family = match.group(1).upper()
                selector = match.group(2).upper()
                group_id = int(match.group(3))
                title, values = self._title_and_ints(
                    lines,
                    index,
                )
                raw_groups.append(
                    (
                        family,
                        selector,
                        group_id,
                        title or "%s_%d" % (family, group_id),
                        values,
                    )
                )

        if not nodes:
            raise RuntimeError("The Radioss deck contains no /NODE block")
        if not elements:
            raise RuntimeError("The Radioss deck contains no supported element block")

        for node_id, coordinates in nodes.items():
            self.add_node(node_id, coordinates)

        type_converter = CellsTypeConverter("RADIOSS")
        connectivity_converter = ConnectivityRenumberer("RADIOSS")

        element_ids = set()
        part_cells = OrderedDict()
        family_cells = OrderedDict()
        element_nodes = {}
        element_keywords = {}

        for cell_id, keyword, part_id, connectivity in elements:
            if cell_id in element_ids:
                raise RuntimeError(
                    "Element IDs must be globally unique; duplicate %d" % cell_id
                )
            element_ids.add(cell_id)

            missing = [node_id for node_id in connectivity if node_id not in nodes]
            if missing:
                raise RuntimeError(
                    "Element %d references missing nodes %s" % (cell_id, missing)
                )

            med_type = type_converter.external_to_medcoupling(keyword)
            connectivity_med = connectivity_converter.external_to_medcoupling(
                med_type,
                connectivity,
            )
            self.add_cell(
                cell_id,
                med_type,
                connectivity_med,
            )

            part_cells.setdefault(part_id, []).append(cell_id)
            family = (
                "SOLID"
                if keyword
                in (
                    "BRICK",
                    "BRIC20",
                    "TETRA4",
                    "TETRA10",
                    "PENTA6",
                )
                else keyword
            )
            family_cells.setdefault(family, []).append(cell_id)
            element_nodes[cell_id] = connectivity
            element_keywords[cell_id] = keyword

        used_names = set()

        def unique(name):
            base = name
            suffix = 2
            while name in used_names:
                name = "%s_%d" % (base, suffix)
                suffix += 1
            used_names.add(name)
            return name

        for part_id, identifiers in part_cells.items():
            name = part_names.get(
                part_id,
                "PART_%d" % part_id,
            )
            self.add_group_cells(
                unique(name),
                identifiers,
            )

        # Aggregate by the finite-element property signature rather than by
        # prop_ID. Thus two /PROP/SHELL cards with Ishell=24 contribute to the
        # same MED group Grp_FE_SHELL_24.
        property_signature_cells = OrderedDict()
        missing_property_ids = set()

        for part_id, identifiers in part_cells.items():
            property_id = part_property_ids.get(part_id)
            if property_id is None:
                continue

            signature = property_signatures.get(property_id)
            if signature is None:
                missing_property_ids.add(property_id)
                continue

            property_signature_cells.setdefault(
                signature,
                [],
            ).extend(identifiers)

        for signature, identifiers in property_signature_cells.items():
            family, formulation = signature
            identifiers = list(dict.fromkeys(identifiers))

            if formulation is None:
                group_name = "Grp_FE_%s" % family
            else:
                group_name = "Grp_FE_%s_%d" % (
                    family,
                    formulation,
                )

            self.add_group_cells(
                unique(group_name),
                identifiers,
            )

        if missing_property_ids:
            logger.info(
                "No /PROP card found for prop_ID(s): %s. "
                "No Grp_FE group was created for these properties."
                % ", ".join(map(str, sorted(missing_property_ids)))
            )

        surface_nodes, surface_cells = self._create_part_ext_surface_cells(
            raw_surfaces,
            raw_groups,
            part_cells,
            part_property_ids,
            part_material_ids,
            element_nodes,
            element_keywords,
            element_ids,
            type_converter,
            unique,
        )

        resolved_groups = self._resolve_groups(
            raw_groups=raw_groups,
            nodes=nodes,
            element_ids=element_ids,
            part_cells=part_cells,
            family_cells=family_cells,
            element_nodes=element_nodes,
            surface_nodes=surface_nodes,
            unique_name=unique,
        )

        self._create_admas_point_cells(
            admas_definitions,
            resolved_groups,
            surface_nodes,
            surface_cells,
            nodes,
            element_ids,
            type_converter,
            unique,
        )

        self._validate_boundary_references(
            boundary_references,
            resolved_groups,
            surface_cells,
        )

        logger.debug("Radioss nodes: %d" % len(nodes))
        logger.debug("Radioss cells: %d" % len(elements))
        logger.debug("Radioss explicit groups: %d" % len(raw_groups))

    def _resolve_groups(
        self,
        raw_groups,
        nodes,
        element_ids,
        part_cells,
        family_cells,
        element_nodes,
        surface_nodes,
        unique_name,
    ):
        resolved = {}
        pending = list(raw_groups)
        maximum_passes = len(pending) + 1
        pass_number = 0

        while pending:
            pass_number += 1

            if pass_number > maximum_passes:
                unresolved = [
                    "%s/%s/%d" % (family, selector, group_id)
                    for family, selector, group_id, _, _ in pending
                ]
                raise RuntimeError(
                    "Maximum number of Radioss group-resolution passes "
                    "exceeded. Possible cyclic dependency: %s" % ", ".join(unresolved)
                )

            previous_pending_size = len(pending)
            rest = []

            logger.debug(
                "Radioss group-resolution pass %d: "
                "%d unresolved group(s)" % (pass_number, previous_pending_size)
            )

            for family, selector, group_id, name, values in pending:
                key = (family, group_id)
                direct_selector = self._direct_selector.get(family)
                result = None

                if selector in (direct_selector, "NODENS"):
                    result = list(values)

                elif selector == "GENE":
                    result = self._generated(
                        values,
                        increment=False,
                    )

                elif selector == "GEN_INCR":
                    result = self._generated(
                        values,
                        increment=True,
                    )

                elif selector == "SURF" and family == "GRNOD":
                    result = []
                    all_surfaces_known = all(
                        abs(surface_id) in surface_nodes for surface_id in values
                    )
                    if all_surfaces_known:
                        for signed_surface_id in values:
                            nodes_on_surface = surface_nodes[abs(signed_surface_id)]
                            if signed_surface_id > 0:
                                result.extend(nodes_on_surface)
                            else:
                                result.extend(-node_id for node_id in nodes_on_surface)
                    else:
                        result = None

                elif selector == "PART":
                    result = self._resolve_part_selector(
                        family,
                        values,
                        part_cells,
                        family_cells,
                        element_nodes,
                    )

                elif selector == family:
                    dependencies = [(family, abs(identifier)) for identifier in values]

                    if all(dependency in resolved for dependency in dependencies):
                        result = []
                        for signed_identifier in values:
                            dependency = (
                                family,
                                abs(signed_identifier),
                            )
                            dependency_values = resolved[dependency]

                            if signed_identifier > 0:
                                result.extend(dependency_values)
                            else:
                                # Keep the signed identifiers. The common
                                # linear-time reducer below performs removal
                                # without repeatedly rebuilding large lists.
                                result.extend(
                                    -identifier for identifier in dependency_values
                                )

                if result is None:
                    rest.append(
                        (
                            family,
                            selector,
                            group_id,
                            name,
                            values,
                        )
                    )
                    continue

                selected = self._apply_signed_identifiers(result)

                if family == "GRNOD":
                    selected = [
                        identifier for identifier in selected if identifier in nodes
                    ]
                    resolved[key] = selected
                    if selected:
                        self.add_group_nodes(
                            unique_name(name),
                            selected,
                        )

                elif family == "GRPART":
                    selected_cells = []
                    for part_id in selected:
                        selected_cells.extend(part_cells.get(part_id, []))
                    selected_cells = list(dict.fromkeys(selected_cells))
                    resolved[key] = selected_cells
                    if selected_cells:
                        self.add_group_cells(
                            unique_name(name),
                            selected_cells,
                        )

                else:
                    selected = [
                        identifier
                        for identifier in selected
                        if identifier in element_ids
                    ]
                    resolved[key] = selected
                    if selected:
                        self.add_group_cells(
                            unique_name(name),
                            selected,
                        )

            pending = rest

            if pending and len(pending) >= previous_pending_size:
                unresolved = [
                    "%s/%s/%d" % (family, selector, group_id)
                    for family, selector, group_id, _, _ in pending
                ]
                raise RuntimeError(
                    "Unable to resolve Radioss groups. Possible cyclic "
                    "dependency or unsupported selector: %s" % ", ".join(unresolved)
                )

        return resolved

    def _create_admas_point_cells(
        self,
        definitions,
        resolved_groups,
        surface_nodes,
        surface_cells,
        nodes,
        element_ids,
        type_converter,
        unique_name,
    ):
        """Materialize nodal /ADMAS definitions as MED POINT1 cells."""
        occupied_cell_ids = set(element_ids)
        for cell_ids in surface_cells.values():
            occupied_cell_ids.update(cell_ids)
        next_cell_id = max(occupied_cell_ids) + 1 if occupied_cell_ids else 1
        point_by_node = {}

        for admas_id, title, direct_nodes, group_ids, surface_ids in definitions:
            selected_nodes = list(direct_nodes)
            for group_id in group_ids:
                selected_nodes.extend(resolved_groups.get(("GRNOD", group_id), []))
            for surface_id in surface_ids:
                selected_nodes.extend(surface_nodes.get(surface_id, []))
            selected_nodes = [n for n in dict.fromkeys(selected_nodes) if n in nodes]

            cell_ids = []
            for node_id in selected_nodes:
                cell_id = point_by_node.get(node_id)
                if cell_id is None:
                    while next_cell_id in occupied_cell_ids:
                        next_cell_id += 1
                    cell_id = next_cell_id
                    next_cell_id += 1
                    med_type = type_converter.external_to_medcoupling("ADMAS")
                    self.add_cell(cell_id, med_type, (node_id,))
                    point_by_node[node_id] = cell_id
                    occupied_cell_ids.add(cell_id)
                cell_ids.append(cell_id)
            if cell_ids:
                self.add_group_cells(unique_name(title), cell_ids)

    @staticmethod
    def _validate_boundary_references(
        references,
        resolved_groups,
        surface_cells,
    ):
        """Ensure every BC/load reference resolves to a non-empty MED group."""
        errors = []

        for keyword, card_id, entity, entity_id in references:
            if entity == "GRNOD":
                values = resolved_groups.get(("GRNOD", entity_id))
            elif entity == "SURF":
                values = surface_cells.get(entity_id)
            else:
                continue

            if values is None:
                errors.append(
                    "/%s/%d references missing %s %d"
                    % (keyword, card_id, entity, entity_id)
                )
            elif not values:
                errors.append(
                    "/%s/%d references empty %s %d"
                    % (keyword, card_id, entity, entity_id)
                )

        if errors:
            raise RuntimeError(
                "Invalid Radioss boundary-condition references:\n  - "
                + "\n  - ".join(errors)
            )

        if references:
            logger.debug(
                "Validated %d Radioss boundary-condition reference(s)" % len(references)
            )

    def _create_part_ext_surface_cells(
        self,
        raw_surfaces,
        raw_groups,
        part_cells,
        part_property_ids,
        part_material_ids,
        element_nodes,
        element_keywords,
        element_ids,
        type_converter,
        unique_name,
    ):
        """Create 2D MED cells and groups for /SURF/PART/EXT surfaces.

        A face occurring once in the selected volume set is external. Faces
        shared by two selected volumes are internal and omitted. Coincident
        faces reused by several Radioss surfaces are represented by one MED
        cell and may belong to several MED groups.
        """
        surface_nodes = {}
        surface_cells = {}
        face_cell_cache = {}
        next_cell_id = max(element_ids) + 1 if element_ids else 1

        for surface_id, (selector, surface_name, items) in raw_surfaces.items():
            if selector == "SEG":
                cell_ids = []
                nodes_on_surface = {}
                for segment in items:
                    med_name = "TRI3" if len(segment) == 3 else "QUAD4"
                    cache_key = (med_name, tuple(sorted(segment)))
                    cell_id = face_cell_cache.get(cache_key)
                    if cell_id is None:
                        med_type = type_converter.external_to_medcoupling(med_name)
                        cell_id = next_cell_id
                        next_cell_id += 1
                        self.add_cell(cell_id, med_type, segment)
                        face_cell_cache[cache_key] = cell_id
                    cell_ids.append(cell_id)
                    for node_id in segment:
                        nodes_on_surface.setdefault(node_id, None)
                if cell_ids:
                    self.add_group_cells(unique_name(surface_name), cell_ids)
                surface_cells[surface_id] = list(dict.fromkeys(cell_ids))
                surface_nodes[surface_id] = list(nodes_on_surface)
                continue

            if selector != "PART/EXT":
                continue
            signed_part_ids = items
            selected_cells = {}
            for signed_part_id in signed_part_ids:
                cells = part_cells.get(abs(signed_part_id), [])
                if signed_part_id > 0:
                    for cell_id in cells:
                        selected_cells.setdefault(cell_id, None)
                else:
                    for cell_id in cells:
                        selected_cells.pop(cell_id, None)

            occurrences = {}
            descriptions = {}
            surface_cell_ids = []
            ordered_nodes = {}

            for source_cell_id in selected_cells:
                keyword = element_keywords[source_cell_id]
                connectivity = element_nodes[source_cell_id]

                if keyword in ("SHELL", "SH3N", "QUAD", "TRIA"):
                    # A shell part is already a 2D surface mesh. Reuse its
                    # native cell instead of trying to extract a lower-
                    # dimensional boundary or duplicating the cell.
                    surface_cell_ids.append(source_cell_id)
                    for node_id in connectivity:
                        ordered_nodes.setdefault(node_id, None)
                    continue

                specifications = self._surface_face_specs.get(keyword)
                if specifications is None:
                    raise RuntimeError(
                        "Cannot create /SURF/PART/EXT cells from element type %s"
                        % keyword
                    )

                for med_name, full_indices, corner_indices in specifications:
                    corners = tuple(connectivity[i] for i in corner_indices)
                    face_key = tuple(sorted(corners))
                    occurrences[face_key] = occurrences.get(face_key, 0) + 1
                    descriptions.setdefault(
                        face_key,
                        (
                            med_name,
                            tuple(connectivity[i] for i in full_indices),
                        ),
                    )

            # Only volume elements require generation of new 2D boundary
            # cells. Faces occurring twice are internal to the selected volume
            # set and are excluded from /SURF/PART/EXT.
            for face_key, count in occurrences.items():
                if count != 1:
                    continue

                med_name, face_connectivity = descriptions[face_key]
                cache_key = (med_name, tuple(sorted(face_connectivity)))
                surface_cell_id = face_cell_cache.get(cache_key)

                if surface_cell_id is None:
                    med_type = type_converter.external_to_medcoupling(med_name)
                    surface_cell_id = next_cell_id
                    next_cell_id += 1
                    self.add_cell(
                        surface_cell_id,
                        med_type,
                        face_connectivity,
                    )
                    face_cell_cache[cache_key] = surface_cell_id

                surface_cell_ids.append(surface_cell_id)
                for node_id in face_connectivity:
                    ordered_nodes.setdefault(node_id, None)

            # A cell may be introduced through several signed part selectors;
            # keep one occurrence while preserving the Radioss order.
            surface_cell_ids = list(dict.fromkeys(surface_cell_ids))

            if surface_cell_ids:
                self.add_group_cells(
                    unique_name(surface_name),
                    surface_cell_ids,
                )
            surface_nodes[surface_id] = list(ordered_nodes)
            surface_cells[surface_id] = list(surface_cell_ids)

        # Resolve surfaces assembled from shell element groups and from
        # previously resolved surfaces. These cards do not create cells; they
        # create additional MED groups over existing 2D cells.
        direct_element_groups = {}
        for family, selector, group_id, _, values in raw_groups:
            if (family, selector) in (("GRSHEL", "SHEL"), ("GRSH3N", "SH3N")):
                direct_element_groups[(family, group_id)] = [
                    cell_id
                    for cell_id in self._apply_signed_identifiers(values)
                    if cell_id in element_ids
                ]

        pending = {
            surface_id: definition
            for surface_id, definition in raw_surfaces.items()
            if definition[0] not in ("PART/EXT", "SEG")
        }

        while pending:
            previous_size = len(pending)
            rest = {}

            for surface_id, (selector, surface_name, item_ids) in pending.items():
                cell_ids = []
                resolvable = True

                if selector in ("GRSHEL", "GRSH3N"):
                    for signed_group_id in item_ids:
                        key = (selector, abs(signed_group_id))
                        group_cells = direct_element_groups.get(key)
                        if group_cells is None:
                            resolvable = False
                            break
                        if signed_group_id > 0:
                            cell_ids.extend(group_cells)
                        else:
                            removals = set(group_cells)
                            cell_ids = [c for c in cell_ids if c not in removals]

                elif selector in ("PART", "MAT", "PROP"):
                    selected_part_ids = []
                    for signed_id in item_ids:
                        identifier = abs(signed_id)
                        if selector == "PART":
                            matches = [identifier]
                        elif selector == "MAT":
                            matches = [
                                p
                                for p, m in part_material_ids.items()
                                if m == identifier
                            ]
                        else:
                            matches = [
                                p
                                for p, prop in part_property_ids.items()
                                if prop == identifier
                            ]
                        if signed_id > 0:
                            selected_part_ids.extend(matches)
                        else:
                            removals = set(matches)
                            selected_part_ids = [
                                p for p in selected_part_ids if p not in removals
                            ]
                    for part_id in dict.fromkeys(selected_part_ids):
                        for cell_id in part_cells.get(part_id, []):
                            if element_keywords[cell_id] in (
                                "SHELL",
                                "SH3N",
                                "QUAD",
                                "TRIA",
                            ):
                                cell_ids.append(cell_id)

                elif selector == "SURF":
                    for signed_surface_id in item_ids:
                        group_cells = surface_cells.get(abs(signed_surface_id))
                        if group_cells is None:
                            resolvable = False
                            break
                        if signed_surface_id > 0:
                            cell_ids.extend(group_cells)
                        else:
                            removals = set(group_cells)
                            cell_ids = [c for c in cell_ids if c not in removals]
                else:
                    resolvable = False

                if not resolvable:
                    rest[surface_id] = (selector, surface_name, item_ids)
                    continue

                cell_ids = list(dict.fromkeys(cell_ids))
                nodes = {}
                for cell_id in cell_ids:
                    for node_id in element_nodes[cell_id]:
                        nodes.setdefault(node_id, None)

                if cell_ids:
                    self.add_group_cells(unique_name(surface_name), cell_ids)
                surface_cells[surface_id] = cell_ids
                surface_nodes[surface_id] = list(nodes)

            if rest and len(rest) >= previous_size:
                unresolved = [
                    "/SURF/%s/%d" % (selector, surface_id)
                    for surface_id, (selector, _, _) in rest.items()
                ]
                raise RuntimeError(
                    "Unable to resolve Radioss surfaces: %s" % ", ".join(unresolved)
                )
            pending = rest

        return surface_nodes, surface_cells

    def _resolve_part_selector(
        self,
        family,
        values,
        part_cells,
        family_cells,
        element_nodes,
    ):
        if family == "GRPART":
            return list(values)

        if family == "GRNOD":
            result = []
            for signed_part_id in values:
                part_nodes = []
                for cell_id in part_cells.get(
                    abs(signed_part_id),
                    [],
                ):
                    part_nodes.extend(element_nodes[cell_id])

                if signed_part_id > 0:
                    result.extend(part_nodes)
                else:
                    result.extend(-node_id for node_id in part_nodes)
            return result

        target_family = self._group_family[family]
        target_cells = set(family_cells.get(target_family, []))
        result = []

        for signed_part_id in values:
            selected_cells = [
                cell_id
                for cell_id in part_cells.get(
                    abs(signed_part_id),
                    [],
                )
                if cell_id in target_cells
            ]

            if signed_part_id > 0:
                result.extend(selected_cells)
            else:
                result.extend(-cell_id for cell_id in selected_cells)

        return result
