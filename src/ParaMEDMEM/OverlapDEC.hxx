// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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

#ifndef __OVERLAPDEC_HXX__
#define __OVERLAPDEC_HXX__

#include "DEC.hxx"
#include "InterpolationOptions.hxx"

#include <mpi.h>

namespace ParaMEDMEM
{
  class OverlapInterpolationMatrix;
  class ProcessorGroup;
  class ParaFIELD;

  class OverlapDEC : public DEC, public INTERP_KERNEL::InterpolationOptions
  {
  public:
    OverlapDEC(const std::set<int>& procIds,const MPI_Comm& world_comm=MPI_COMM_WORLD);
    virtual ~OverlapDEC();
    void sendRecvData(bool way=true);
    void sendData();
    void recvData();
    void synchronize();
    void attachSourceLocalField(ParaFIELD *field, bool ownPt=false);
    void attachTargetLocalField(ParaFIELD *field, bool ownPt=false);
    ProcessorGroup *getGrp() { return _group; }
    bool isInGroup() const;
  private:
    bool _own_group;
    OverlapInterpolationMatrix* _interpolation_matrix;
    ProcessorGroup *_group;

    ParaFIELD *_source_field;
    bool _own_source_field;
    ParaFIELD *_target_field;
    bool _own_target_field;
    MPI_Comm _comm;
  };
}

#endif
