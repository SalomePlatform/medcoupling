// Data exchange channel for trio
// TrioDEC.cxx
// version 0.0 06/06/2014

#include "TrioDEC.hxx"

#include "ICoCoMEDField.hxx"
#include "ICoCoTrioField.hxx"

using namespace ParaMEDMEM;
using namespace ICoCo;
 
TrioDEC::TrioDEC():_traduced_field(0)
{
}

TrioDEC::TrioDEC(ProcessorGroup& source_group, ProcessorGroup& target_group):InterpKernelDEC(source_group,target_group),_traduced_field(0)
{
}

TrioDEC::TrioDEC(const std::set<int>& src_ids, const std::set<int>& trg_ids, const MPI_Comm& world_comm):InterpKernelDEC(src_ids,trg_ids,world_comm),_traduced_field(0)
{
}

void TrioDEC::attachLocalField(ICoCo::TrioField *field)
{
  if(!field)
    throw INTERP_KERNEL::Exception("TrioDEC::attachLocalField : The input trio Field is NULL !");
  releaseInternalPointer();
  _traduced_field=field->build_medfield();
  DisjointDEC::attachLocalField(_traduced_field);
}

void TrioDEC::releaseInternalPointer()
{
  if(_traduced_field)
    delete _traduced_field;
  _traduced_field=0;
}

TrioDEC::~TrioDEC()
{
  releaseInternalPointer();
}
