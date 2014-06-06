// Data exchange channel for trio
// TrioDEC.h
// version 0.0 06/06/2014

#ifndef _TrioDEC_included_
#define _TrioDEC_included_

#include "InterpKernelDEC.hxx"

namespace ICoCo
{
  class MEDField;
  class TrioField;
}

namespace ParaMEDMEM
{
  class TrioDEC : public InterpKernelDEC
  {
  public:  
    TrioDEC();
    TrioDEC(ProcessorGroup& source_group, ProcessorGroup& target_group);
    TrioDEC(const std::set<int>& src_ids, const std::set<int>& trg_ids, const MPI_Comm& world_comm=MPI_COMM_WORLD);
    void attachLocalField(ICoCo::TrioField *field);
    virtual ~TrioDEC();
  private:
    void releaseInternalPointer();
  private :
    ICoCo::MEDField *_traduced_field;
  };
}

#endif
