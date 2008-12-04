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
#ifndef PARASUPPORT_HXX_
#define PARASUPPORT_HXX_

#include "MEDMEM_Exception.hxx"
namespace MEDMEM
{
  class SUPPORT;
}
namespace ParaMEDMEM
{
  class Topology;
  class ParaMESH;
  class ProcessorGroup;
  class ParaSUPPORT
  {
  public:
    ParaSUPPORT();
    ParaSUPPORT(const ParaMESH* mesh, const MEDMEM::SUPPORT* support):
      _support(support), _mesh(mesh), _has_support_ownership(false), _has_mesh_ownership(false){};
		ParaSUPPORT(const ParaMESH* mesh):_has_mesh_ownership(false){};
		ParaSUPPORT(const MEDMEM::SUPPORT&, const ProcessorGroup& proc_group);
		virtual ~ParaSUPPORT();
    virtual const Topology* getTopology() const =0;
    virtual const MEDMEM::SUPPORT* getSupport() const {return _support;}
    virtual const ParaMESH* getMesh() const {return _mesh;}
		virtual const int* getGlobalNumbering() const;

  protected :
		
    const MEDMEM::SUPPORT* _support;
    const ParaMESH* _mesh;
		bool _has_support_ownership;
    bool _has_mesh_ownership;
  };
  
}


#endif /*PARASUPPORT_HXX_*/
