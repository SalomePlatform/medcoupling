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


class CellsTypeConverter:

    _systus_to_med = {
        "0001": "POINT1",
        "1002": "SEG2",
        "2003": "TRI3",
        "2004": "QUAD4",
        "3004": "TETRA4",
        "3008": "HEXA8",
        "3005": "PYRA5",
        "3006": "PENTA6",
        "1003": "SEG3",
        "2006": "TRI6",
        "2008": "QUAD8",
        "3010": "TETRA10",
        "3020": "HEXA20",
        "3013": "PYRA13",
        "3015": "PENTA15",
        "1004": "SEG4",
    }

    _abaqus_to_med = OrderedDict(
        (
            ("Node", "POINT1"),
            # Mass element - 0D
            ("MASS", "POINT1"),
            # Frame Element
            ("FRAME2D", "SEG2"),
            ("FRAME3D", "SEG2"),
            # Joint Element
            ("JOINTC", "SEG2"),
            # Spring element
            ("SPRING1", "POINT1"),
            ("SPRING2", "SEG2"),
            ("SPRINGA", "SEG2"),
            # Truss Element
            ("T3D2", "SEG2"),
            ("T3D3", "SEG3"),
            ("T3D3H", "SEG3"),
            # Elbow Element
            ("ELBOW31", "SEG2"),
            ("ELBOW32", "SEG3"),
            ("ELBOW31B", "SEG2"),
            ("ELBOW31C", "SEG2"),
            # Diffusive heat transfer elements
            ("DC1D2", "SEG2"),
            ("DC1D3", "SEG3"),
            # Forced convection heat transfer element
            ("DCC1D2", "SEG2"),
            ("DCC1D2D", "SEG2"),
            # Coupled thermal-electrical elements
            ("DC1D2E", "SEG2"),
            ("DC1D3E", "SEG3"),
            # Axis symmetric thermal element
            ("DCCAX4", "QUAD4"),
            # Acoustic element
            ("AC1D2", "SEG2"),
            ("AC1D3", "SEG3"),
            # Rigid element
            ("R2D2", "SEG2"),
            ("R2D3", "SEG3"),
            ("RB2D2", "SEG2"),
            ("RB3D2", "SEG2"),
            ("RB2D3", "SEG3"),
            ("RAX2", "SEG2"),
            ("RAX3", "SEG3"),
            ("R3D3", "TRI3"),
            ("R3D4", "QUAD4"),
            ("R3D6", "TRI6"),
            ("R3D8", "QUAD8"),
            # Beam element
            ("B21", "SEG2"),
            ("B22", "SEG3"),
            ("B23", "SEG2"),
            ("B21H", "SEG2"),
            ("B22H", "SEG3"),
            ("B23H", "SEG2"),
            ("B31", "SEG2"),
            ("B32", "SEG3"),
            ("B33", "SEG2"),
            ("B31H", "SEG2"),
            ("B32H", "SEG3"),
            ("B33H", "SEG2"),
            # PIPE element
            ("PIPE21", "SEG2"),
            ("PIPE22", "SEG3"),
            ("PIPE21H", "SEG2"),
            ("PIPE22H", "SEG3"),
            ("PIPE31", "SEG2"),
            ("PIPE32", "SEG3"),
            ("PIPE31H", "SEG2"),
            ("PIPE32H", "SEG3"),
            # Membrane Element
            ("M3D3", "TRI3"),
            ("M3D4", "QUAD4"),
            ("M3D6", "TRI6"),
            ("M3D8", "QUAD8"),
            ("M3D9", "QUAD9"),
            ("MAX1", "SEG2"),
            ("MAX2", "SEG3"),
            ("MGAX1", "SEG2"),
            ("MGAX2", "SEG3"),
            # Eulerian Element
            ("EC3D8R", "HEXA8"),
            # Reduced Structural Element
            ("S3R", "TRI3"),
            ("S3RS", "TRI3"),
            ("S4R", "QUAD4"),
            ("S4RS", "QUAD4"),
            ("S4R5", "QUAD4"),
            ("S4RSW", "QUAD4"),
            ("S8R", "QUAD8"),
            ("S8RS", "QUAD8"),
            ("S9R5", "QUAD9"),
            # Structural Element
            ("S3", "TRI3"),
            ("S4", "QUAD4"),
            # Reduced Continuum Element
            ("CPE4R", "QUAD4"),
            ("CPE4RH", "QUAD4"),
            ("CPE8R", "QUAD8"),
            ("CPE8RH", "QUAD8"),
            ("CPS4R", "QUAD4"),
            ("CPS8R", "QUAD8"),
            ("C3D8R", "HEXA8"),
            ("C3D8RH", "HEXA8"),
            ("C3D20R", "HEXA20"),
            ("C3D20RH", "HEXA20"),
            ("C3D27R", "HEXA27"),
            ("C3D27RH", "HEXA27"),
            # Plane strain element
            ("CPE3", "TRI3"),
            ("CPE6", "TRI6"),
            ("CPE4", "QUAD4"),
            ("CPE8", "QUAD8"),
            ("CPE9", "QUAD9"),
            ("CPE3H", "TRI3"),
            ("CPE6H", "TRI6"),
            ("CPE4H", "QUAD4"),
            ("CPE8H", "QUAD8"),
            ("CPE9H", "QUAD9"),
            ("CPE6M", "TRI6"),
            ("CPE6MH", "TRI6"),
            ("CPE4I", "QUAD4"),
            ("CPE4IH", "QUAD4"),
            # Plane stress element
            ("CPS3", "TRI3"),
            ("CPS6", "TRI6"),
            ("CPS6M", "TRI6"),
            ("CPS4", "QUAD4"),
            ("CPS4I", "QUAD4"),
            ("CPS8", "QUAD8"),
            # Axi element
            ("CAX3", "TRI3"),
            ("CAX6", "TRI6"),
            ("CAX4", "QUAD4"),
            ("CAX8", "QUAD8"),
            ("CAX9", "QUAD9"),
            ("CAX3H", "TRI3"),
            ("CAX6H", "TRI6"),
            ("CAX4H", "QUAD4"),
            ("CAX8H", "QUAD8"),
            ("CAX9H", "QUAD9"),
            ("CAX6M", "TRI6"),
            ("CAX6MH", "TRI6"),
            ("CAX4I", "QUAD4"),
            ("CAX4IH", "QUAD4"),
            ("CAX4R", "QUAD4"),
            ("CAX4RH", "QUAD4"),
            ("CAX8R", "QUAD8"),
            ("CAX8RH", "QUAD8"),
            # Continuum Element - Hybrid element
            ("C3D4H", "TETRA4"),
            ("C3D10H", "TETRA10"),
            ("C3D10M", "TETRA10"),
            ("C3D10MH", "TETRA10"),
            ("C3D5H", "PYRA5"),
            ("C3D13H", "PYRA13"),
            ("C3D6H", "PENTA6"),
            ("C3D15H", "PENTA15"),
            ("C3D15VH", "PENTA18"),
            ("C3D8I", "HEXA8"),
            ("C3D8IH", "HEXA8"),
            ("C3D8H", "HEXA8"),
            ("C3D20H", "HEXA20"),
            ("C3D27H", "HEXA27"),
            # Continuum Element, must be declared last
            ("C3D4", "TETRA4"),
            ("C3D10", "TETRA10"),
            ("C3D5", "PYRA5"),
            ("C3D13", "PYRA13"),
            ("C3D6", "PENTA6"),
            ("C3D15", "PENTA15"),
            ("C3D15V", "PENTA18"),
            ("C3D8", "HEXA8"),
            ("C3D20", "HEXA20"),
            ("C3D27", "HEXA27"),
        )
    )

    _aster_to_med = {
        "POI1": "POINT1",
        "SEG2": "SEG2",
        "TRIA3": "TRI3",
        "QUAD4": "QUAD4",
        "TETRA4": "TETRA4",
        "HEXA8": "HEXA8",
        "PYRAM5": "PYRA5",
        "PENTA6": "PENTA6",
        "SEG3": "SEG3",
        "TRIA6": "TRI6",
        "QUAD8": "QUAD8",
        "TETRA10": "TETRA10",
        "HEXA20": "HEXA20",
        "PYRAM13": "PYRA13",
        "PENTA15": "PENTA15",
        "SEG4": "SEG4",
        "TRIA7": "TRI7",
        "QUAD9": "QUAD9",
        "PENTA18": "PENTA18",
        "HEXA27": "HEXA27",
    }

    _zset_to_med = OrderedDict(
        (
            ("l2d2", "SEG2"),
            ("l2d3", "SEG3"),
            ("c2d3r", "TRI3"),
            ("c2d3", "TRI3"),
            ("c2d6r", "TRI6"),
            ("c2d6", "TRI6"),
            ("c2d4r", "QUAD4"),
            ("c2d4", "QUAD4"),
            ("c2d8r", "QUAD8"),
            ("c2d8", "QUAD8"),
            ("c3d4r", "TETRA4"),
            ("c3d4", "TETRA4"),
            ("c3d10r", "TETRA10"),
            ("c3d10_4", "TETRA10"),
            ("c3d10", "TETRA10"),
            ("c3d8r", "HEXA8"),
            ("c3d8", "HEXA8"),
            ("c3d20r", "HEXA20"),
            ("c3d20", "HEXA20"),
            ("c3d6r", "PENTA6"),
            ("c3d6", "PENTA6"),
            ("c3d15_9", "PENTA15"),
            ("c3d15r", "PENTA15"),
            ("c3d15", "PENTA15"),
            ("c3d5_6", "PYRA5"),
            ("c3d5_27", "PYRA5"),
            ("c3d5", "PYRA5"),
            ("c3d13_27", "PYRA13"),
            ("c3d13r", "PYRA13"),
            ("c3d13", "PYRA13"),
        )
    )

    _tetgen_to_med = {"TRIA3": "TRI3", "TETRA4": "TETRA4"}

    _ansys_to_med = OrderedDict(
        (
            # Mesh200
            ("200_2_2", "SEG2"),
            ("200_3_3", "SEG3"),
            ("200_3_4", "TRI3"),
            ("200_6_5", "TRI6"),
            ("200_4_6", "QUAD4"),
            ("200_8_7", "QUAD8"),
            ("200_4_8", "TETRA4"),
            ("200_10_9", "TETRA10"),
            ("200_8_10", "HEXA8"),
            ("200_27_11", "HEXA27"),
            # Mass element - 0D
            ("21_1", "POINT1"),
            ("71_1", "POINT1"),
            ("166_1", "POINT1"),
            # SOLID5
            ("5_8", "HEXA8"),
            ("5_6", "PENTA6"),
            # LINK11
            ("11_2", "SEG2"),
            # PLANE13
            ("13_4", "QUAD4"),
            ("13_3", "TRI3"),
            # COMBIN14
            ("14_2", "SEG2"),
            # PIPE16
            ("16_2", "SEG2"),
            ("16_3", "SEG3"),
            # PIPE18
            ("18_2", "SEG2"),
            ("18_3", "SEG3"),
            # PLANE25
            ("25_4", "QUAD4"),
            ("25_3", "TRI3"),
            # MATRIX27
            ("27_2", "SEG2"),
            # FLUID29
            ("29_4", "QUAD4"),
            ("29_3", "TRI3"),
            # FLUID30
            ("30_8", "HEXA8"),
            ("30_4", "TETRA4"),
            ("30_6", "PENTA6"),
            ("30_5", "PYRA5"),
            # LINK31
            ("31_2", "SEG2"),
            # LINK33
            ("33_2", "SEG2"),
            # LINK34
            ("34_2", "SEG2"),
            # PLANE35
            ("35_6", "TRI6"),
            # COMBIN39
            ("39_2", "SEG2"),
            # COMBIN40
            ("40_2", "SEG2"),
            # PLANE42
            ("42_4", "QUAD4"),
            # SOLID45
            ("45_8", "HEXA8"),
            ("45_6", "PENTA6"),
            # INFIN47
            ("47_4", "QUAD4"),
            ("47_3", "TRI3"),
            # PLANE55
            ("55_4", "QUAD4"),
            ("55_3", "TRI3"),
            # PIPE59
            ("59_2", "SEG2"),
            # SHELL61
            ("61_2", "SEG2"),
            # SHELL63
            ("63_4", "QUAD4"),
            ("63_3", "TRI3"),
            # SOLID65
            ("65_8", "HEXA8"),
            ("65_6", "PENTA6"),
            # LINK68
            ("68_2", "SEG2"),
            # SOLID70
            ("70_8", "HEXA8"),
            ("70_4", "TETRA4"),
            ("70_6", "PENTA6"),
            ("70_5", "PYRA5"),
            # PLANE75
            ("75_4", "QUAD4"),
            ("75_3", "TRI3"),
            # PLANE77
            ("77_4", "QUAD8"),
            ("77_3", "TRI6"),
            # PLANE78
            ("78_4", "QUAD8"),
            ("78_3", "TRI6"),
            # PLANE82
            ("82_8", "QUAD8"),
            ("82_6", "TRI6"),
            # PLANE83
            ("83_4", "QUAD8"),
            ("83_3", "TRI6"),
            # SOLID87
            ("87_10", "TETRA10"),
            # SOLID90
            ("90_20", "HEXA20"),
            ("90_10", "TETRA10"),
            ("90_13", "PYRA13"),
            ("90_15", "PENTA15"),
            # SOLID92
            ("92_10", "TETRA10"),
            # SOLID95
            ("95_20", "HEXA20"),
            ("95_10", "TETRA10"),
            ("95_13", "PYRA13"),
            ("95_15", "PENTA15"),
            # SOLID96
            ("96_8", "HEXA8"),
            ("96_4", "TETRA4"),
            ("96_6", "PENTA6"),
            ("96_5", "PYRA5"),
            # SOLID98
            ("98_10", "TETRA10"),
            # INFIN111
            ("111_20", "HEXA20"),
            ("111_8", "HEXA8"),
            ("111_15", "PENTA15"),
            ("111_6", "PENTA6"),
            # FLUID116
            ("116_2", "SEG2"),
            # SOLID120
            ("120_20", "HEXA20"),
            ("120_10", "TETRA10"),
            ("120_13", "PYRA13"),
            ("120_15", "PENTA15"),
            # PLANE121
            ("121_8", "QUAD8"),
            ("121_6", "TRI6"),
            # SOLID122
            ("122_20", "HEXA20"),
            ("122_10", "TETRA10"),
            ("122_13", "PYRA13"),
            ("122_15", "PENTA15"),
            # SOLID123
            ("123_10", "TETRA10"),
            # FLUID129
            ("129_2", "SEG2"),
            # SHELL131
            ("131_4", "QUAD4"),
            ("131_3", "TRI3"),
            # SHELL132
            ("132_8", "QUAD8"),
            ("132_6", "TRI6"),
            # FLUID136
            ("136_4", "QUAD4"),
            ("136_3", "TRI3"),
            ("136_8", "QUAD8"),
            ("136_6", "TRI6"),
            # FLUID138
            ("138_2", "SEG2"),
            # SHELL43 : n'existe plus mais ce comporte comme un SHELL181
            ("43_4", "QUAD4"),
            ("43_3", "TRI3"),
            # SHELL143 : n'existe plus mais ce comporte comme un SHELL181
            ("143_4", "QUAD4"),
            ("143_3", "TRI3"),
            # SURF151
            ("151_2", "SEG2"),
            ("151_3", "SEG3"),
            # SURF152
            ("152_3", "TRI3"),
            ("152_4", "QUAD4"),
            ("152_6", "TRI6"),
            ("152_8", "QUAD8"),
            # SURF153
            ("153_2", "SEG2"),
            ("153_3", "SEG3"),
            # SURF156
            ("156_2", "SEG2"),
            # SURF154
            ("154_3", "TRI3"),
            ("154_6", "TRI6"),
            ("154_4", "QUAD4"),
            ("154_8", "QUAD8"),
            # SHELL157
            ("157_4", "QUAD4"),
            ("157_3", "TRI3"),
            # LINK160
            ("160_2", "SEG2"),
            # BEAM161
            ("161_2", "SEG2"),
            # SHELL162
            ("162_4", "QUAD4"),
            ("162_3", "TRI3"),
            # SHELL163
            ("163_4", "QUAD4"),
            ("163_3", "TRI3"),
            # SOLID164
            ("164_8", "HEXA8"),
            ("164_6", "PENTA6"),
            ("164_4", "TETRA4"),
            ("164_5", "PYRA5"),
            # COMBIN165
            ("165_2", "SEG2"),
            # LINK167
            ("167_2", "SEG2"),
            # SOLID168
            ("168_10", "TETRA10"),
            # TARGE169
            ("169_1", "POINT1"),
            ("169_2", "SEG2"),
            ("169_3", "SEG3"),
            # TARGE170
            ("170_1", "POINT1"),
            ("170_2", "SEG2"),
            ("170_3", "TRI3"),
            ("170_4", "QUAD4"),
            ("170_6", "TRI6"),
            ("170_8", "QUAD8"),
            # CONTACT171
            ("171_2", "SEG2"),
            # CONTACT172
            ("172_2", "SEG2"),
            ("172_3", "SEG3"),
            # CONTACT173
            ("173_4", "QUAD4"),
            ("173_3", "TRI3"),
            # CONTACT174
            ("174_4", "QUAD4"),
            ("174_3", "TRI3"),
            ("174_8", "QUAD8"),
            ("174_6", "TRI6"),
            # CONTACT175
            ("175_1", "POINT1"),
            # CONTACT176
            ("176_2", "SEG2"),
            ("176_3", "SEG3"),
            # CONTACT177
            ("177_2", "SEG2"),
            ("177_3", "SEG3"),
            # LINK180
            ("180_2", "SEG2"),
            # SHELL181
            ("181_4", "QUAD4"),
            ("181_3", "TRI3"),
            # PLANE182
            ("182_4", "QUAD4"),
            ("182_3", "TRI3"),
            # PLANE183
            ("183_8", "QUAD8"),
            ("183_6", "TRI6"),
            # MPC184
            ("184_2", "SEG2"),
            # SOLID185
            ("185_8", "HEXA8"),
            ("185_6", "PENTA6"),
            ("185_4", "TETRA4"),
            ("185_5", "PYRA5"),
            # SOLID186
            ("186_20", "HEXA20"),
            ("186_10", "TETRA10"),
            ("186_13", "PYRA13"),
            ("186_15", "PENTA15"),
            # SOLID187
            ("187_10", "TETRA10"),
            # BEAM4
            ("4_2", "SEG2"),
            ("4_3", "SEG3"),
            # BEAM44
            ("44_2", "SEG2"),
            ("44_3", "SEG3"),
            # BEAM188
            ("188_2", "SEG2"),
            # BEAM189
            ("189_2", "SEG2"),
            ("189_3", "SEG3"),
            # SOLSH190
            ("190_8", "HEXA8"),
            ("190_6", "PENTA6"),
            # INTER192
            ("192_4", "QUAD4"),
            # INTER195
            ("195_8", "HEXA8"),
            # INTER202
            ("202_4", "QUAD4"),
            # INTER204
            ("204_5", "HEXA8"),
            # CPT212
            ("212_4", "QUAD4"),
            ("212_3", "TRI3"),
            # CPT213
            ("213_8", "QUAD8"),
            ("213_3", "TRI6"),
            # CPT215
            ("215_8", "HEXA8"),
            ("215_6", "PENTA6"),
            ("215_4", "TETRA4"),
            # CPT216
            ("216_20", "HEXA20"),
            ("216_10", "TETRA10"),
            ("216_13", "PYRA13"),
            ("216_15", "PENTA15"),
            # CPT217
            ("217_10", "TETRA10"),
            # FLUID218
            ("218_4", "QUAD4"),
            ("218_3", "TRI3"),
            # FLUID220
            ("220_20", "HEXA20"),
            ("220_10", "TETRA10"),
            ("220_13", "PYRA13"),
            ("220_15", "PENTA15"),
            # FLUID221
            ("221_10", "TETRA10"),
            # PLANE223
            ("223_8", "QUAD8"),
            ("223_6", "TRI6"),
            # SOLID226
            ("226_20", "HEXA20"),
            ("226_10", "TETRA10"),
            ("226_13", "PYRA13"),
            ("226_15", "PENTA15"),
            # SOLID227
            ("227_10", "TETRA10"),
            # PLANE230
            ("230_8", "QUAD8"),
            ("230_6", "TRI6"),
            # SOLID231
            ("231_20", "HEXA20"),
            ("231_10", "TETRA10"),
            ("231_13", "PYRA13"),
            ("231_15", "PENTA15"),
            # SOLID232
            ("232_10", "TETRA10"),
            # PLANE233
            ("233_8", "QUAD8"),
            ("233_6", "TRI6"),
            # SOLID236
            ("236_20", "HEXA20"),
            ("236_10", "TETRA10"),
            ("236_13", "PYRA13"),
            ("236_15", "PENTA15"),
            # SOLID237
            ("237_10", "TETRA10"),
            # PLANE238
            ("238_8", "QUAD8"),
            ("238_6", "TRI6"),
            # SOLID239
            ("239_20", "HEXA20"),
            ("239_10", "TETRA10"),
            ("239_13", "PYRA13"),
            ("239_15", "PENTA15"),
            # SOLID240
            ("240_10", "TETRA10"),
            # SURF251
            ("251_2", "SEG2"),
            # SURF252
            ("252_4", "QUAD4"),
            ("252_3", "TRI3"),
            # SOLID278
            ("278_8", "HEXA8"),
            ("278_6", "PENTA6"),
            ("278_4", "TETRA4"),
            ("278_5", "PYRA5"),
            # SOLID279
            ("279_20", "HEXA20"),
            ("279_10", "TETRA10"),
            ("279_13", "PYRA13"),
            ("279_15", "PENTA15"),
            # SHELL281
            ("281_8", "QUAD8"),
            ("281_6", "TRI6"),
            # SOLID285
            ("285_4", "TETRA4"),
            # PIPE288
            ("288_2", "SEG2"),
            # PIPE289
            ("289_3", "SEG3"),
            # ELBOW290
            ("290_3", "SEG3"),
        )
    )

    _med_types = "POINT1 SEG2 TRI3 QUAD4 TETRA4 HEXA8 PYRA5 PENTA6 SEG3 TRI6 QUAD8 TETRA10 HEXA20 PYRA13 PENTA15 SEG4 TRI7 QUAD9 PENTA18 HEXA27".split()

    def __init__(self, code):
        import MEDLoader as ml

        self.code = code.lower()

        try:
            data = getattr(self, "_{}_to_med".format(self.code))
        except AttributeError as exc:
            raise RuntimeError("Unknown format '{}'".format(code)) from exc

        assert set(data.values()) <= set(self._med_types)
        mdata = {i: getattr(ml, "NORM_%s" % k) for i, k in data.items()}

        self._external_to_medcoupling = {i: k for i, k in mdata.items()}
        self._medcoupling_to_external = {k: i for i, k in mdata.items()}

        if "systus" in self.code:
            self._f_e2m = self._systus_to_mc
            self._f_m2e = self._to_ext
        else:
            self._f_e2m = self._to_mc
            self._f_m2e = self._to_ext

    def external_to_medcoupling(self, external_type):
        return self._f_e2m(external_type)

    def medcoupling_to_external(self, medcoupling_type):
        return self._f_m2e(medcoupling_type)

    # Specific functions
    def _systus_to_mc(self, systus_type):
        dim, stype, nb_nodes = systus_type[0], systus_type[1], systus_type[-2:]

        if not stype in ("0", "1", "2", "3"):
            raise RuntimeError(
                "Cannot convert {} type '{}'".format(*(self.code.title(), systus_type))
            )

        return self._to_mc("0".join((dim, nb_nodes)))

    # Generic functions
    def _to_mc(self, external_type):
        try:
            return self._external_to_medcoupling[external_type]
        except KeyError as exc:
            raise RuntimeError(
                "{} type '{}' unknown.".format(*(self.code.title(), external_type))
            ) from exc

    def _to_ext(self, medcoupling_type):
        try:
            return self._medcoupling_to_external[medcoupling_type]
        except KeyError as exc:
            raise RuntimeError(
                "MedCoupling type '{}' unknown.".format(medcoupling_type)
            ) from exc


class GroupCellsTypeConverter(CellsTypeConverter):

    _systus_to_med = {}
    _abaqus_to_med = {}
    _aster_to_med = {}
    _ansys_to_med = {}

    _zset_to_med = OrderedDict(
        (
            ("line", "SEG2"),
            ("quad", "SEG3"),
            ("t3", "TRI3"),
            ("t6", "TRI6"),
            ("q4", "QUAD4"),
            ("q8", "QUAD8"),
        )
    )

    def __init__(self, code):
        super().__init__(code)
