#ifndef ExplicitCOINCIDENTDEC_HXX_
#define ExplicitCOINCIDENTDEC_HXX_

#include "DEC.hxx"
#include "ExplicitMapping.hxx"
#include "ExplicitTopology.hxx"
#include <map>

namespace ParaMEDMEM
{
  class DEC;
  class BlockTopology;
  //  class ExplicitMapping;
  class ExplicitCoincidentDEC: public DEC
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
  private :
    
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

#endif /*ExplicitCOINCIDENTDEC_HXX_*/
	
