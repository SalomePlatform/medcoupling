// Copyright (C) 2007-2023  CEA, EDF
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
#include "MCIdType.hxx"

namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT CellSimplify
  {
  public:
    static INTERP_KERNEL::NormalizedCellType simplifyDegeneratedCell(INTERP_KERNEL::NormalizedCellType type, const mcIdType *conn, mcIdType lgth, mcIdType *retConn, mcIdType& retLgth);
    static mcIdType *getFullPolyh3DCell(INTERP_KERNEL::NormalizedCellType type, const mcIdType *conn, mcIdType lgth,
                                   mcIdType& retNbOfFaces, mcIdType& retLgth);
    static INTERP_KERNEL::NormalizedCellType tryToUnPoly2D(bool isQuad, const mcIdType *conn, mcIdType lgth, mcIdType *retConn, mcIdType& retLgth);
    static INTERP_KERNEL::NormalizedCellType tryToUnPoly3D(const mcIdType *conn, mcIdType nbOfFaces, mcIdType lgth, mcIdType *retConn, mcIdType& retLgth);
    static INTERP_KERNEL::NormalizedCellType tryToUnPolyHex8(const mcIdType *conn, mcIdType nbOfFaces, mcIdType lgth, mcIdType *retConn, mcIdType& retLgth);
    static INTERP_KERNEL::NormalizedCellType tryToUnPolyHexp12(const mcIdType *conn, mcIdType nbOfFaces, mcIdType lgth, mcIdType *retConn, mcIdType& retLgth);
    static INTERP_KERNEL::NormalizedCellType tryToUnPolyPenta6(const mcIdType *conn, mcIdType nbOfFaces, mcIdType lgth, mcIdType *retConn, mcIdType& retLgth);
    static INTERP_KERNEL::NormalizedCellType tryToUnPolyPyra5(const mcIdType *conn, mcIdType nbOfFaces, mcIdType lgth, mcIdType *retConn, mcIdType& retLgth);
    static INTERP_KERNEL::NormalizedCellType tryToUnPolyTetra4(const mcIdType *conn, mcIdType nbOfFaces, mcIdType lgth, mcIdType *retConn, mcIdType& retLgth);
    static bool tryToArrangeOppositeFace(const mcIdType *conn, mcIdType lgth, mcIdType lgthBaseFace, const mcIdType *baseFace, const mcIdType *oppFaceId, mcIdType nbOfFaces, mcIdType *retConnOfOppFace);
    static bool isWellOriented(const mcIdType *baseFace, mcIdType *retConn, const mcIdType *sideFace, mcIdType lgthBaseFace);
    static bool orientOppositeFace(const mcIdType *baseFace, mcIdType *retConn, const mcIdType *sideFace, mcIdType lgthBaseFace);
    static bool isFlatCell(const mcIdType* conn, mcIdType pos, mcIdType lgth, NormalizedCellType type);
  };
}

#endif
