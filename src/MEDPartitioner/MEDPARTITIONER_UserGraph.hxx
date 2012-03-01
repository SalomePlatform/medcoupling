//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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
#ifndef MEDPARTITIONER_USERGRAPH_HXX_
#define MEDPARTITIONER_USERGRAPH_HXX_

#include "MEDPARTITIONER.hxx"
#include "MEDPARTITIONER_Graph.hxx"
#include <string>
namespace MEDPARTITIONER
{
  class MEDSKYLINEARRAY;
  class ParaDomainSelector;
  class MEDPARTITIONER_EXPORT UserGraph : public Graph
  {
  public:
    UserGraph(MEDPARTITIONER::MEDSKYLINEARRAY*, const int*, int);
    virtual ~UserGraph();
    void partGraph(int, const std::string& options=std::string(""), ParaDomainSelector* sel=0);
  };
}

#endif /*MEDPARTITIONER_USERGRAPH_HXX_*/