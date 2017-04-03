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

#ifndef __OVERLAPDEC_HXX__
#define __OVERLAPDEC_HXX__

#include "DEC.hxx"
#include "InterpolationOptions.hxx"

#include <mpi.h>
#include <string>

namespace ICoCo {
  class MEDField;
}

namespace MEDCoupling
{
  class OverlapInterpolationMatrix;
  class OverlapElementLocator;
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
    void attachSourceLocalField(MEDCouplingFieldDouble *field);
    void attachTargetLocalField(MEDCouplingFieldDouble *field);
    void attachSourceLocalField(ICoCo::MEDField *field);
    void attachTargetLocalField(ICoCo::MEDField *field);
    ProcessorGroup *getGroup() { return _group; }
    bool isInGroup() const;

    void setDefaultValue(double val) {_default_field_value = val;}
    //! 0 means initial algo from Antho, 1 or 2 means Adrien's algo (2 should be better). Make your choice :-))
    void setWorkSharingAlgo(int method)  { _load_balancing_algo = method; }

    void debugPrintWorkSharing(std::ostream & ostr) const;
  private:
    int _load_balancing_algo;

    bool _own_group;
    OverlapInterpolationMatrix* _interpolation_matrix;
    OverlapElementLocator* _locator;
    ProcessorGroup *_group;

    double _default_field_value;

    ParaFIELD *_source_field;
    bool _own_source_field;
    ParaFIELD *_target_field;
    bool _own_target_field;
    MPI_Comm _comm;
  };
}

#endif
