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
#ifndef UNSTRUCTUREDPARASUPPORT_HXX_
#define UNSTRUCTUREDPARASUPPORT_HXX_

#include "ParaSUPPORT.hxx"
#include "MEDMEM_define.hxx"
#include "Topology.hxx"
#include "ExplicitTopology.hxx"

using namespace MED_EN;
namespace MEDMEM
{
	class SUPPORT;
}

namespace ParaMEDMEM
{
  class Topology;
  class ExplicitTopology;
  class ParaMESH;
  

  class UnstructuredParaSUPPORT:public ParaSUPPORT
  {
  public:
    
    UnstructuredParaSUPPORT(const ParaMESH* const mesh, const MEDMEM::SUPPORT* support );
    UnstructuredParaSUPPORT(const MEDMEM::SUPPORT* support , const ProcessorGroup& group);
    
    UnstructuredParaSUPPORT(const ParaMESH* const mesh, const MED_EN::medEntityMesh entity);
    virtual ~UnstructuredParaSUPPORT(); 
    const Topology* getTopology() const
      {return _explicit_topology;}

  private:
    const Topology*  _explicit_topology;
    const MED_EN::medEntityMesh _entity;

	
  };
  
}
#endif /*STRUCTUREDPARASUPPORT_HXX_*/
