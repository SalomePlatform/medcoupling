// Copyright (C) 2026  CEA, EDF
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#include "MEDLoaderGeoTypesConverter.hxx"

namespace MEDCoupling
{
med_geometry_type typmai[MED_N_CELL_FIXED_GEO] = {MED_POINT1,  MED_SEG2,    MED_SEG3,     MED_SEG4,      MED_TRIA3,
                                                  MED_QUAD4,   MED_TRIA6,   MED_TRIA7,    MED_QUAD8,     MED_QUAD9,
                                                  MED_TETRA4,  MED_PYRA5,   MED_PENTA6,   MED_HEXA8,     MED_OCTA12,
                                                  MED_TETRA10, MED_PYRA13,  MED_PENTA15,  MED_PENTA18,   MED_HEXA20,
                                                  MED_HEXA27,  MED_POLYGON, MED_POLYGON2, MED_POLYHEDRON};

INTERP_KERNEL::NormalizedCellType typmai2[MED_N_CELL_FIXED_GEO] = {
    INTERP_KERNEL::NORM_POINT1, INTERP_KERNEL::NORM_SEG2,    INTERP_KERNEL::NORM_SEG3,    INTERP_KERNEL::NORM_SEG4,
    INTERP_KERNEL::NORM_TRI3,   INTERP_KERNEL::NORM_QUAD4,   INTERP_KERNEL::NORM_TRI6,    INTERP_KERNEL::NORM_TRI7,
    INTERP_KERNEL::NORM_QUAD8,  INTERP_KERNEL::NORM_QUAD9,   INTERP_KERNEL::NORM_TETRA4,  INTERP_KERNEL::NORM_PYRA5,
    INTERP_KERNEL::NORM_PENTA6, INTERP_KERNEL::NORM_HEXA8,   INTERP_KERNEL::NORM_HEXGP12, INTERP_KERNEL::NORM_TETRA10,
    INTERP_KERNEL::NORM_PYRA13, INTERP_KERNEL::NORM_PENTA15, INTERP_KERNEL::NORM_PENTA18, INTERP_KERNEL::NORM_HEXA20,
    INTERP_KERNEL::NORM_HEXA27, INTERP_KERNEL::NORM_POLYGON, INTERP_KERNEL::NORM_QPOLYG,  INTERP_KERNEL::NORM_POLYHED
};

med_geometry_type typmai3[INTERP_KERNEL::NORM_MAXTYPE] = {
    MED_POINT1,      // 0
    MED_SEG2,        // 1
    MED_SEG3,        // 2
    MED_TRIA3,       // 3
    MED_QUAD4,       // 4
    MED_POLYGON,     // 5
    MED_TRIA6,       // 6
    MED_TRIA7,       // 7
    MED_QUAD8,       // 8
    MED_QUAD9,       // 9
    MED_SEG4,        // 10
    MED_NONE,        // 11
    MED_NONE,        // 12
    MED_NONE,        // 13
    MED_TETRA4,      // 14
    MED_PYRA5,       // 15
    MED_PENTA6,      // 16
    MED_NONE,        // 17
    MED_HEXA8,       // 18
    MED_NONE,        // 19
    MED_TETRA10,     // 20
    MED_NONE,        // 21
    MED_OCTA12,      // 22
    MED_PYRA13,      // 23
    MED_NONE,        // 24
    MED_PENTA15,     // 25
    MED_NONE,        // 26
    MED_HEXA27,      // 27
    MED_PENTA18,     // 28
    MED_NONE,        // 29
    MED_HEXA20,      // 30
    MED_POLYHEDRON,  // 31
    MED_POLYGON2     // 32
};

}  // namespace MEDCoupling
