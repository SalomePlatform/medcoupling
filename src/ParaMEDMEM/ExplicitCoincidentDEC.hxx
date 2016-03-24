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

#ifndef __EXPLICITCOINCIDENTDEC_HXX__
#define __EXPLICITCOINCIDENTDEC_HXX__

#include "DisjointDEC.hxx"
#include "ExplicitMapping.hxx"
#include "ExplicitTopology.hxx"

#include <map>

namespace MEDCoupling
{
  class BlockTopology;

  class ExplicitCoincidentDEC : public DisjointDEC
  {
  public:
    ExplicitCoincidentDEC();
    virtual ~ExplicitCoincidentDEC();
    void synchronize();
    void broadcastTopology(BlockTopology*&, int tag);
    void broadcastTopology(const ExplicitTopology* toposend, ExplicitTopology* toporecv, int tag);
    void transferMappingToSource();
    void prepareSourceDE();
    void prepareTargetDE();
    void recvData();
    void sendData();
  private:  
    ExplicitTopology* _toposource;
    ExplicitTopology* _topotarget;
    ProcessorGroup* _targetgroup;
    ProcessorGroup* _sourcegroup;
    int* _sendcounts;
    int* _recvcounts;
    int* _senddispls;
    int* _recvdispls;
    double* _recvbuffer;
    double* _sendbuffer;
    std::map<int,std::pair<int,int> > _distant_elems;
    ExplicitMapping _explicit_mapping;
  }; 
}

#endif
