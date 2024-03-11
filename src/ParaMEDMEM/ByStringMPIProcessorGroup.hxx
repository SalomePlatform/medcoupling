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

#ifndef __BYSTRINGMPIPROCESSORGROUP_HXX__
#define __BYSTRINGMPIPROCESSORGROUP_HXX__

#include "MPIProcessorGroup.hxx"

namespace MEDCoupling
{
  class CommInterface;

  class ByStringMPIProcessorGroup : public MPIProcessorGroup
  {
  public:
    ByStringMPIProcessorGroup(const CommInterface& interface);
    ByStringMPIProcessorGroup(const CommInterface& interface, const std::string& simCodeTag, const MPI_Comm& world_comm=MPI_COMM_WORLD);
    ByStringMPIProcessorGroup(const ByStringMPIProcessorGroup& other);
    virtual ~ByStringMPIProcessorGroup();
    virtual ByStringMPIProcessorGroup *deepCopy() const;
    
  };
}

#endif
