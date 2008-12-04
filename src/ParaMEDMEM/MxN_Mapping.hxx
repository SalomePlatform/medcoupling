//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#ifndef MXN_MAPPING_HXX_
#define MXN_MAPPING_HXX_

#include <vector>

#include "MEDMEM_Field.hxx"
#include "MPI_AccessDEC.hxx"
#include "DECOptions.hxx"

namespace ParaMEDMEM
{

class ProcessorGroup;

	class MxN_Mapping : public DECOptions
{
public:
	MxN_Mapping();
  MxN_Mapping(const ProcessorGroup& local_group, const ProcessorGroup& distant_group, const DECOptions& dec_options);
	virtual ~MxN_Mapping();
  void addElementFromSource(int distant_proc, int distant_elem);
  void prepareSendRecv();
  void sendRecv(MEDMEM::FIELD<double>& field);
  void sendRecv(double* field, MEDMEM::FIELD<double>& field) const ;
	void reverseSendRecv(double* field, MEDMEM::FIELD<double>& field) const ;
 
	MPI_AccessDEC* getAccessDEC(){return _accessDEC;}
private :
//  ProcessorGroup& _local_group;
//  ProcessorGroup& _distant_group;
  ProcessorGroup* _union_group;
  MPI_AccessDEC * _accessDEC;
  int _nb_comps;
  std::vector<pair<int,int> > _sending_ids;
  std::vector<int> _recv_ids;
  std::vector<int> _send_proc_offsets;
  std::vector<int> _recv_proc_offsets;
};

ostream & operator<< (ostream &,const AllToAllMethod &);

}

#endif /*MXN_MAPPING_HXX_*/
