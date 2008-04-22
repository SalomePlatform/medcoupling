#ifndef PROCESSORGROUP_HXX_
#define PROCESSORGROUP_HXX_
#include <set>
#include "CommInterface.hxx"


namespace ParaMEDMEM
{
class CommInterface;

class ProcessorGroup
{
public:
  
  ProcessorGroup(const CommInterface& interface):_comm_interface(interface){}
  ProcessorGroup(const CommInterface& interface, std::set<int> proc_ids):
    _comm_interface(interface),_proc_ids(proc_ids){}
  ProcessorGroup (const ProcessorGroup& proc_group, std::set<int> proc_ids):
    _comm_interface(proc_group.getCommInterface()){}
  ProcessorGroup (const CommInterface& interface, int start, int end);
  virtual ~ProcessorGroup(){}
  virtual ProcessorGroup* fuse (const ProcessorGroup&) const=0;
  virtual void intersect (ProcessorGroup&)=0;
  bool contains(int rank) const {return _proc_ids.find(rank)!=_proc_ids.end();};
  virtual bool containsMyRank() const=0;
  int size() const  {return _proc_ids.size();}
  const CommInterface& getCommInterface()const {return _comm_interface;};
  virtual int myRank() const =0;
  virtual int translateRank(const ProcessorGroup*, int) const =0;
  virtual ProcessorGroup* createComplementProcGroup() const =0;
  virtual ProcessorGroup* createProcGroup() const=0;
  virtual const std::set<int>& getProcIDs()const  {return _proc_ids;} 
protected:
  const CommInterface _comm_interface;
  std::set<int> _proc_ids;
};
  
}

#endif /*PROCESSORGROUP_HXX_*/
