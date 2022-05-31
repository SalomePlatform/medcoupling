// Copyright (C) 2007-2022  CEA/DEN, EDF R&D
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

#ifndef __EXPLICITMAPPING_HXX__
#define __EXPLICITMAPPING_HXX__

#include "MCIdType.hxx"

#include <vector>
#include <map>
#include <set>

namespace MEDCoupling
{
  /*!
   * Internal class, not part of the public API.
   *
   * Used by the ExplicitCoincidentDEC.
   */
  class ExplicitMapping
  {
  public:
    ExplicitMapping();
    ~ExplicitMapping();
    
    void pushBackElem(std::pair<int,mcIdType> idistant);
    void  setDistantElem(mcIdType ilocal, std::pair<int,mcIdType> idistant);
    int nbDistantDomains();
    std::pair <int,mcIdType> getDistantNumbering(mcIdType ielem) const;
    
    int getDistantDomain(int i);
    int getNbDistantElems(int i);
    mcIdType* serialize(int idproc);
    void unserialize(int nbprocs, int* sizes,int nbtarget, int* targetrank, mcIdType* commbuffer);
    
    int* getBufferIndex() const { return _buffer_index; }
    int* getCounts() const { return _send_counts; }
  private:
    std::vector <std::pair<int,mcIdType> > _mapping;
    std::set<int> _distant_domains;
    int* _numbers;
    int* _domains;
    mcIdType* _comm_buffer;
    int* _buffer_index;
    int* _send_counts;

    void computeNumbers();
  };
}

#endif
