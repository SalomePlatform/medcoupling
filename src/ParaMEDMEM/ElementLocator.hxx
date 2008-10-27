#ifndef ELEMENTLOCATOR_HXX_
#define ELEMENTLOCATOR_HXX_

#include <vector>
#include <set>
#include "InterpolationOptions.hxx"

namespace MEDMEM
{
  class MESH;
}
namespace ParaMEDMEM
{
class ParaMESH;
class ProcessorGroup;
class ParaSUPPORT;
class InterpolationMatrix;


	class ElementLocator: public INTERP_KERNEL::InterpolationOptions
{
public:
	ElementLocator(const ParaMESH& mesh, const ProcessorGroup& distant_group);
  ElementLocator(const ParaSUPPORT& support, const ProcessorGroup& distant_group);
	virtual ~ElementLocator();
  void exchangeMesh(int idistantrank, MEDMEM::MESH*& distant_mesh, int*& distant_ids);
 
private:
  MEDMEM::MESH* _local_mesh;
  std::vector<MEDMEM::MESH*> _distant_meshes;
  double* _domain_bounding_boxes;
  const ProcessorGroup& _distant_group;
  const ProcessorGroup& _local_group;
  ProcessorGroup* _union_group;
  std::vector<int> _distant_proc_ids;
  
  void _computeBoundingBoxes();
  bool _intersectsBoundingBox(int irank);
  bool _intersectsBoundingBox(double* bb1, double* bb2, int dim);
  void _exchangeMesh(MEDMEM::MESH* local_mesh, MEDMEM::MESH*& distant_mesh, int iproc_distant, const int* distant_ids_send, int*& distant_ids_recv);
  MEDMEM::MESH* _meshFromElems(std::set<int>& elems);
};

}

#endif /*ELEMENTLOCATOR_HXX_*/
