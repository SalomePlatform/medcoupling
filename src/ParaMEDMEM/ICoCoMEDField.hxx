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
#ifndef _ICOCOMEDFIELD_HXX_
#define _ICOCOMEDFIELD_HXX_

#include <vector>
#include "ICoCoField.hxx"
#include "MEDMEM_Field.hxx"

namespace ParaMEDMEM
{
	class ParaMESH;
	class ParaFIELD;
	class ParaSUPPORT;
	class ComponentTopology;
	class ProcessorGroup;
}
namespace MEDMEM
{
	class CONNECTZONE;
	class MESHING;
	class SUPPORT;
}
namespace ICoCo
{
	class TrioField;
	
  class MEDField:public ICoCo::Field
    {
    public:
      
      MEDField(){};
      MEDField(ParaMEDMEM::ParaMESH* mesh, ParaMEDMEM::ParaFIELD* field);
			MEDField(TrioField& , const ParaMEDMEM::ProcessorGroup& group);
			virtual ~MEDField();
      ParaMEDMEM::ParaFIELD* getField() const  {return _field;}
      ParaMEDMEM::ParaMESH* getMesh()const {return _mesh;}

    private:
      ParaMEDMEM::ParaMESH* _mesh;
      ParaMEDMEM::ParaFIELD* _field;
      MEDMEM::FIELD<double>* _local_field;
			bool _has_field_ownership;
			MEDMEM::MESHING* _local_mesh;
			MEDMEM::SUPPORT* _local_support;
			ParaMEDMEM::ParaSUPPORT* _support;
			ParaMEDMEM::ComponentTopology* _comp_topology;
			
      std::vector<MEDMEM::CONNECTZONE*> _connectzones;
    };
};

#endif
