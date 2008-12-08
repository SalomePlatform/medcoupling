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
#include "ParaMESH.hxx"
#include "ParaSUPPORT.hxx"
#include "ParaMESH.hxx"
#include "MEDMEM_Support.hxx"

namespace ParaMEDMEM
{

  ParaSUPPORT::ParaSUPPORT()
  {
  }
  
  ParaSUPPORT::ParaSUPPORT(const MEDMEM::SUPPORT& support, const ProcessorGroup& proc_group):
  _support(&support), 
  _has_support_ownership(false),
  _has_mesh_ownership(true)
   {
    _mesh = new ParaMESH(*(support.getMesh()),  proc_group, "mesh from support");
  } 

  ParaSUPPORT::~ParaSUPPORT()
  {
		if (_has_support_ownership)
			{
				delete _support;
				_support=0;
			}
    if (_has_mesh_ownership)
			{
				delete _mesh;
				_mesh=0;
			}
  }

	const int* ParaSUPPORT::getGlobalNumbering() const
	{
		if (! _support->isOnAllElements())
			throw MEDMEM::MEDEXCEPTION("GlobalNumbering can only be retrieved on supports on all elements");
		return _mesh->getGlobalNumbering(_support->getEntity());
	}

}

