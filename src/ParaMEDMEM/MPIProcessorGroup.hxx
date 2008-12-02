#ifndef MPIPROCESSORGROUP_HXX_
#define MPIPROCESSORGROUP_HXX_

#include <set>
#include <mpi.h>

using namespace std;
namespace ParaMEDMEM
{
class ProcessorGroup;
class CommInterface;

class MPIProcessorGroup:public ProcessorGroup
{
public:
  MPIProcessorGroup(const CommInterface& interface);
  MPIProcessorGroup(const CommInterface& interface, set<int> proc_ids, const MPI_Comm& world_comm=MPI_COMM_WORLD);
  MPIProcessorGroup (const ProcessorGroup& proc_group, set<int> proc_ids);
  MPIProcessorGroup(const CommInterface& interface,int pstart, int pend);
  virtual ~MPIProcessorGroup();
  virtual ProcessorGroup* fuse (const ProcessorGroup&) const;
  void intersect (ProcessorGroup&){};
  int myRank() const {int rank; MPI_Comm_rank(_comm,&rank); return rank;}
  bool containsMyRank() const { int rank; MPI_Group_rank(_group, &rank); return (rank!=MPI_UNDEFINED);}
  int translateRank(const ProcessorGroup* group, int rank) const;
  const MPI_Comm* getComm() const {return &_comm;}
  ProcessorGroup* createComplementProcGroup() const;
  ProcessorGroup* createProcGroup() const;
  MPI_Comm getWorldComm() { return _world_comm;}
  
private:
  MPI_Group _group;
  const MPI_Comm _world_comm;
  MPI_Comm _comm;
};

}

#endif /*MPIPROCESSORGROUP_HXX_*/
