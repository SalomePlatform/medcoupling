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

#ifndef __DEC_HXX__
#define __DEC_HXX__

#include "MEDCouplingFieldDouble.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "DECOptions.hxx"

namespace MEDCoupling
{
  class CommInterface;

  /*!
   * DEC stands for Data Exchange Channel. See the page \ref para-dec for more on this.
   *
   * This class is purely abstract. See the derivations:
   * - \ref DisjointDEC-det "DisjointDEC"
   * - \ref NonCoincidentDEC "NonCoincidentDEC"
   * - \ref OverlapDEC "OverlapDEC"
   */
  class DEC : public DECOptions
  {
  public:
    DEC();
    void copyFrom(const DEC& other);
    virtual void synchronize() = 0;
    virtual void sendRecvData(bool way=true) = 0;
    virtual ~DEC();
  protected:
    const CommInterface* _comm_interface;
  };
}

#endif
