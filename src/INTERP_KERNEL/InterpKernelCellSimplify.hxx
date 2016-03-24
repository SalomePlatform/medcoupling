// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (CEA/DEN)

#ifndef __INTERPKERNELCELLSIMPLIFY_HXX__
#define __INTERPKERNELCELLSIMPLIFY_HXX__

#include "INTERPKERNELDefines.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpKernelException.hxx"

namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT CellSimplify
  {
  public:
    static INTERP_KERNEL::NormalizedCellType simplifyDegeneratedCell(INTERP_KERNEL::NormalizedCellType type, const int *conn, int lgth, int *retConn, int& retLgth);
    static int *getFullPolyh3DCell(INTERP_KERNEL::NormalizedCellType type, const int *conn, int lgth,
                                   int& retNbOfFaces, int& retLgth);
    static INTERP_KERNEL::NormalizedCellType tryToUnPoly2D(bool isQuad, const int *conn, int lgth, int *retConn, int& retLgth);
    static INTERP_KERNEL::NormalizedCellType tryToUnPoly3D(const int *conn, int nbOfFaces, int lgth, int *retConn, int& retLgth);
    static INTERP_KERNEL::NormalizedCellType tryToUnPolyHex8(const int *conn, int nbOfFaces, int lgth, int *retConn, int& retLgth);
    static INTERP_KERNEL::NormalizedCellType tryToUnPolyHexp12(const int *conn, int nbOfFaces, int lgth, int *retConn, int& retLgth);
    static INTERP_KERNEL::NormalizedCellType tryToUnPolyPenta6(const int *conn, int nbOfFaces, int lgth, int *retConn, int& retLgth);
    static INTERP_KERNEL::NormalizedCellType tryToUnPolyPyra5(const int *conn, int nbOfFaces, int lgth, int *retConn, int& retLgth);
    static INTERP_KERNEL::NormalizedCellType tryToUnPolyTetra4(const int *conn, int nbOfFaces, int lgth, int *retConn, int& retLgth);
    static bool tryToArrangeOppositeFace(const int *conn, int lgth, int lgthBaseFace, const int *baseFace, const int *oppFaceId, int nbOfFaces, int *retConnOfOppFace);
    static bool isWellOriented(const int *baseFace, int *retConn, const int *sideFace, int lgthBaseFace);
    static bool orientOppositeFace(const int *baseFace, int *retConn, const int *sideFace, int lgthBaseFace);
  };
}

#endif
