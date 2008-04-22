#ifndef MXN_MAPPING_HXX_
#define MXN_MAPPING_HXX_

#include <vector>

#include "MEDMEM_Field.hxx"
#include "MPI_AccessDEC.hxx"


namespace ParaMEDMEM
{
typedef enum{Native,PointToPoint} AllToAllMethod;

class ProcessorGroup;

class MxN_Mapping
{
public:
	MxN_Mapping();
  MxN_Mapping(const ProcessorGroup& local_group, const ProcessorGroup& distant_group);
	virtual ~MxN_Mapping();
  void addElementFromSource(int distant_proc, int distant_elem);
  void prepareSendRecv();
  void sendRecv(MEDMEM::FIELD<double>& field);
  void sendRecv(double* field, MEDMEM::FIELD<double>& field) const ;
	void reverseSendRecv(double* field, MEDMEM::FIELD<double>& field) const ;
  void setAllToAllMethod(const AllToAllMethod& method){_allToAllMethod = method ;};
	MPI_AccessDEC* getAccessDEC(){return _accessDEC;}
private :
//  ProcessorGroup& _local_group;
//  ProcessorGroup& _distant_group;
  ProcessorGroup* _union_group;
  MPI_AccessDEC * _accessDEC;
  AllToAllMethod _allToAllMethod ;
  int _nb_comps;
  std::vector<pair<int,int> > _sending_ids;
  std::vector<int> _recv_ids;
  std::vector<int> _send_proc_offsets;
  std::vector<int> _recv_proc_offsets;
};

ostream & operator<< (ostream &,const AllToAllMethod &);

}

#endif /*MXN_MAPPING_HXX_*/
