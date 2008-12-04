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
