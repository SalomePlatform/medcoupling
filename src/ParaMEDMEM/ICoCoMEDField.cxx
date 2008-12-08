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
#include "ICoCoMEDField.hxx"
#include "ICoCoTrioField.hxx"
#include "ProcessorGroup.hxx"
#include "UnstructuredParaSUPPORT.hxx"
#include "ParaFIELD.hxx"
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Meshing.hxx"
#include "MEDMEM_Support.hxx"
#include "MEDMEM_ArrayInterface.hxx"
#include "MEDMEM_Array.hxx"
#include "MEDMEM_Field.hxx"

namespace ICoCo
{

	/*! Constructor directly attaching a ParaMESH and a ParaFIELD
		the object does not take the control the objects pointed by 
		\a mesh and \a field.
	*/
		
	MEDField::MEDField(ParaMEDMEM::ParaMESH* mesh, ParaMEDMEM::ParaFIELD* field): 
		_mesh(mesh),
		_field(field),
		_has_field_ownership(false),
		_local_mesh(0),
		_local_support(0),
		_support(0),
		_comp_topology(0)
	{
	}
	
	MEDField::MEDField(TrioField& triofield, const ParaMEDMEM::ProcessorGroup& group):
	_has_field_ownership(true)
	{
		_local_mesh = new MEDMEM::MESHING();
		string system="CARTESIAN";
		_local_mesh->setCoordinates(triofield._space_dim,
																triofield._nbnodes,
																triofield._coords,
																system,
																MED_EN::MED_FULL_INTERLACE);
		_local_mesh->setNumberOfTypes(1, MED_EN::MED_CELL);
		MED_EN::medGeometryElement elemtype;
		switch (triofield._mesh_dim)
			{
			case 2:
				switch (triofield._nodes_per_elem)
					{
					case 3:
						elemtype=MED_TRIA3;
						break;
					case 4 : 
						elemtype=MED_QUAD4;
						break;
					default:
						throw MEDMEM::MEDEXCEPTION("incompatible Trio field - wrong nb of nodes per elem");
					}
				break;
			case 3:
				switch (triofield._nodes_per_elem)
					{
					case 4:
						elemtype=MED_TETRA4;
						break;
					case 8 : 
						elemtype=MED_HEXA8;
						break;
					default:
						throw MEDMEM::MEDEXCEPTION("incompatible Trio field - wrong nb of nodes per elem");
					}
				break;
			default:
				throw MEDMEM::MEDEXCEPTION("incompatible Trio field - wrong mesh dimension");
			}
		_local_mesh->setTypes (&elemtype, MED_EN::MED_CELL);
		_local_mesh->setNumberOfElements(&triofield._nb_elems, MED_EN::MED_CELL);

		//creating a connectivity table that complies to MED (1 indexing)
		//and passing it to _local_mesh
		int* conn=new int[triofield._nb_elems*triofield._nodes_per_elem];
		for (int i=0; i<triofield._nb_elems*triofield._nodes_per_elem;i++)
		  conn[i]=(triofield._connectivity)[i]+1;
		_local_mesh->setConnectivity(conn, MED_EN::MED_CELL, elemtype);
		delete[] conn;


		_local_mesh->setMeshDimension(triofield._mesh_dim);
		
		_mesh=new ParaMEDMEM::ParaMESH(*_local_mesh, group, "support for trio field");
		_local_support = new MEDMEM::SUPPORT(_local_mesh,"support on all cells", MED_EN::MED_CELL);
		_support=new ParaMEDMEM::UnstructuredParaSUPPORT(_local_support,group);
		_comp_topology=new ParaMEDMEM::ComponentTopology(triofield._nb_field_components);
		
    //field on the sending end
		if (triofield._field!=0)
			{
        _field =  new ParaMEDMEM::ParaFIELD(_support, *_comp_topology );
        _field->getField()->setName(triofield.getName());
				_field->getField()->setTime(triofield._time1);
				_field->getField()->setIterationNumber(triofield._itnumber);
				_field->getField()->setOrderNumber(0);
				for (int i =0; i<triofield._nb_elems; i++)
					for (int j=0; j<triofield._nb_field_components; j++)
						{
							_field->getField()->setValueIJ(i+1,j+1,triofield._field[i*triofield._nb_field_components+j]);
						}
			}
    //field on the receiving end
		else
			{
        #if 0
        _local_field=new MEDMEM::FIELD<double>();
        _local_field->setNumberOfComponents(triofield._nb_field_components);
        _local_field->setSupport(_local_support);
        
        MEDMEM::MEDMEM_ArrayInterface<double,FullInterlace,MEDMEM::NoGauss>::Array* array=
        new MEDMEM::MEDMEM_ArrayInterface<double,FullInterlace,MEDMEM::NoGauss>::Array
        (triofield._field, _local_field->getNumberOfComponents(),_local_field->getNumberOfValues(),true,true);
        
        _local_field->setArray(array);
        #endif
        _field =  new ParaMEDMEM::ParaFIELD(_support, *_comp_topology );
        _field->getField()->setName(triofield.getName());
        _field->getField()->setTime(triofield._time1);
        _field->getField()->setIterationNumber(triofield._itnumber);
        _field->getField()->setOrderNumber(0);
				// the trio field points to the pointer inside the MED field
				triofield._field=const_cast<double*> (_field->getField()->getValue());
				for (int i=0; i<triofield._nb_field_components*triofield._nb_elems;i++)
					triofield._field[i]=0.0;
			}
	
	}
	MEDField::~MEDField()
	{
		delete _local_mesh;
		delete _local_support;
		   _local_support=0;
		delete _comp_topology;
		delete _support;
		_support=0;
		if (_has_field_ownership)
			{
				delete _field;
				delete _mesh;
			}
	}
};
