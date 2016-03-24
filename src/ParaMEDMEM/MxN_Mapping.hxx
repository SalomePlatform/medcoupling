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

#ifndef __MXN_MAPPING_HXX__
#define __MXN_MAPPING_HXX__

#include "MEDCouplingFieldDouble.hxx"
#include "MPIAccessDEC.hxx"
#include "DECOptions.hxx"

#include <vector>

namespace MEDCoupling
{

  class ProcessorGroup;

  /*!
   * Internal class, not part of the public API.
   *
   * Used by InterpolationMatrix. This class manages the mapping between a given processor and part
   * of the mesh (cell ids).
   */
  class MxN_Mapping : public DECOptions
  {
  public:
    MxN_Mapping(const ProcessorGroup& source_group, const ProcessorGroup& target_group, const DECOptions& dec_options);
    virtual ~MxN_Mapping();
    void addElementFromSource(int distant_proc, int distant_elem);
    void prepareSendRecv();
    void sendRecv(MEDCouplingFieldDouble& field);
    void sendRecv(double* sendfield, MEDCouplingFieldDouble& field) const ;
    void reverseSendRecv(double* recvfield, MEDCouplingFieldDouble& field) const ;
 
    //
    const std::vector<std::pair<int,int> >& getSendingIds() const { return _sending_ids; }
    const std::vector<int>& getSendProcsOffsets() const { return _send_proc_offsets; }
    void initialize();

    MPIAccessDEC* getAccessDEC(){ return _access_DEC; }
  private :
    ProcessorGroup* _union_group;
    MPIAccessDEC * _access_DEC;
    int _nb_comps;
    std::vector<std::pair<int,int> > _sending_ids;
    std::vector<int> _recv_ids;
    std::vector<int> _send_proc_offsets;
    std::vector<int> _recv_proc_offsets;
  };

  std::ostream & operator<< (std::ostream &,const AllToAllMethod &);

}

#endif
