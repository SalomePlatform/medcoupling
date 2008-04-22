#ifndef STRUCTUREDPARASUPPORT_HXX_
#define STRUCTUREDPARASUPPORT_HXX_

#include "ParaSUPPORT.hxx"
#include "MEDMEM_define.hxx"

using namespace MED_EN;
namespace MEDMEM
{
	class SUPPORT;
}

namespace ParaMEDMEM
{
class BlockTopology;
class ParaGRID;
class ParaMESH;

class StructuredParaSUPPORT:public ParaSUPPORT
{
public:
	
	StructuredParaSUPPORT(const ParaGRID* const grid, const MED_EN::medEntityMesh entity);
	StructuredParaSUPPORT(const ParaMESH* const mesh, const MED_EN::medEntityMesh entity);
	
	virtual ~StructuredParaSUPPORT();
	const Topology* getTopology() const {return _block_topology;}
	const MEDMEM::SUPPORT* getSupport() {return _support;}
	const ParaMESH* getParaMesh()const {return _mesh;}
	
private:
	const BlockTopology* const  _block_topology;
	const ParaGRID* const _grid;
  //	const ParaMESH* const _mesh;
	const MED_EN::medEntityMesh _entity;
  //	const MEDMEM::SUPPORT* _support;
	
};

}
#endif /*STRUCTUREDPARASUPPORT_HXX_*/
