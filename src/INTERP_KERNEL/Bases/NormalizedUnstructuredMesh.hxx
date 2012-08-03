// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#ifndef __NORMALIZEDUNSTRUCTUREDMESH_HXX__
#define __NORMALIZEDUNSTRUCTUREDMESH_HXX__

namespace INTERP_KERNEL
{
  typedef enum
    {
      ALL_C_MODE       ,
      ALL_FORTRAN_MODE
    } NumberingPolicy;


  typedef enum
    {
      NORM_POINT1  =  0,
      NORM_SEG2    =  1,
      NORM_SEG3    =  2,
      NORM_SEG4    =  10,
      NORM_POLYL   =  33,
      NORM_TRI3    =  3,
      NORM_QUAD4   =  4,
      NORM_POLYGON =  5,
      NORM_TRI6    =  6,
      NORM_TRI7    =  7,
      NORM_QUAD8   =  8,
      NORM_QUAD9   =  9,
      NORM_QPOLYG  =  32,
      //
      NORM_TETRA4  = 14,
      NORM_PYRA5   = 15,
      NORM_PENTA6  = 16,
      NORM_HEXA8   = 18,
      NORM_TETRA10 = 20,
      NORM_HEXGP12 = 22,
      NORM_PYRA13  = 23,
      NORM_PENTA15 = 25,
      NORM_HEXA20  = 30,
      NORM_HEXA27  = 27,
      NORM_POLYHED = 31,
      NORM_ERROR   = 40,
      NORM_MAXTYPE = 33
    } NormalizedCellType;

  class GenericMesh
  {};
}

#endif
