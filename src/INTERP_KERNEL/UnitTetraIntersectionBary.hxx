//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
// File      : UnitTetraIntersectionBary.hxx
// Created   : Tue Dec  9 16:06:33 2008
// Author    : Edward AGAPOV (eap)
#ifndef __UNITTETRAINTERSECTIONBARY_HXX__
#define __UNITTETRAINTERSECTIONBARY_HXX__

#include "TransformedTriangle.hxx"
#include "INTERPKERNELDefines.hxx"

#include <vector>
#include <list>

namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT UnitTetraIntersectionBary : protected TransformedTriangle
  {
  public:
    UnitTetraIntersectionBary();

    void init();
    /*!
     * \brief Stores a part of triangle common with the unit tetrahedron
     *  \param triangle - triangle side of other cell, whose calculateIntersectionVolume()
     *                    must have already been called
     */
    void addSide(const TransformedTriangle& triangle);

    /*!
     * \brief Computes and return coordinates of barycentre
     */
    bool getBary(double* baryCenter);

    /*!
     * \brief Returns volume of intersection
     *  \retval double - 
     */
    inline double getVolume() const { return _int_volume; }

    virtual ~UnitTetraIntersectionBary();

  private:

    int addSideFaces();

    void setTriangleOnSide(int i);

    void clearPolygons(bool andFaces=false);

    /// volume of intersection
    double  _int_volume;

    /// faces of intersection polyhedron
    std::list< std::vector< double* > > _faces;
  };

}

#endif
