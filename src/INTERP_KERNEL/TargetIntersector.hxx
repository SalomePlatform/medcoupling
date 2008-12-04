//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D, OPEN CASCADE
//
//  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
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
#ifndef __TARGETINTERSECTOR__HXX__
#define __TARGETINTERSECTOR__HXX__

#include <INTERPKERNEL_defines.hxx>

namespace INTERP_KERNEL
{
  /**
   * \brief Abstract base class of Intersector classes. 
   * These classes represent a target element and calculate its intersection
   * with the source elements.
   */
  template<class ConnType>
  class INTERPKERNEL_EXPORT TargetIntersector
  {
  public:

    /// Virtual destructor
    virtual ~TargetIntersector() {}
    
    /**
     * Calculate the volume of the intersection between target cell 
     * and the given source cell.
     *
     * @param srcCell     global number of the source cell
     */
    virtual double intersectSourceCell(ConnType srcCell) = 0;
  };
};

#endif
