#include "ProcessorGroup.hxx"
#include "MEDMEM_Exception.hxx"
namespace ParaMEDMEM
{

  ProcessorGroup::ProcessorGroup (const CommInterface& interface, int start, int end):_comm_interface(interface)
  {
    if (start>end) throw MEDMEM::MEDEXCEPTION("wrong call to Processor group constructor");
    for (int i=start; i<=end;i++)
      _proc_ids.insert(i);
  }
  
}
