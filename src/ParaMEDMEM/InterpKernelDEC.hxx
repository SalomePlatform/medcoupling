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

#ifndef __INTERPKERNELDEC_HXX__
#define __INTERPKERNELDEC_HXX__

#include "DisjointDEC.hxx"
#include "MxN_Mapping.hxx"
#include "InterpolationOptions.hxx"

namespace MEDCoupling
{
  class InterpolationMatrix;

  class InterpKernelDEC : public DisjointDEC, public INTERP_KERNEL::InterpolationOptions
  {
  public:  
    InterpKernelDEC();
    InterpKernelDEC(ProcessorGroup& source_group, ProcessorGroup& target_group);
    InterpKernelDEC(const std::set<int>& src_ids, const std::set<int>& trg_ids,
                    const MPI_Comm& world_comm=MPI_COMM_WORLD);
    virtual ~InterpKernelDEC();
    void synchronize();
    void recvData();
    void recvData(double time);
    void sendData();
    void sendData(double time , double deltatime);
    void prepareSourceDE() { }
    void prepareTargetDE() { }
  private :
    //Number of distant points to be located locally 
    int _nb_distant_points;
    //coordinates of distant points 
    const double* _distant_coords;
    //local element number containing the distant points  
    const int* _distant_locations; 
    InterpolationMatrix* _interpolation_matrix;
  };
}

#endif
