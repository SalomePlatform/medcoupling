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
#include "ICoCoMEDField.hxx"
#include "ICoCoTrioField.hxx"
#include "ProcessorGroup.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "NormalizedUnstructuredMesh.hxx"

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
    _support(0),
    _comp_topology(0)
  {
  }
  
  MEDField::MEDField(TrioField& triofield, const ParaMEDMEM::ProcessorGroup& group):
    _has_field_ownership(true),_support(0)
  {
    _local_mesh = ParaMEDMEM::MEDCouplingUMesh::New();
    _local_mesh->setMeshDimension(triofield._space_dim);
    ParaMEDMEM::DataArrayDouble *myCoords=ParaMEDMEM::DataArrayDouble::New();
    myCoords->alloc(triofield._nbnodes,triofield._space_dim);
    _local_mesh->setCoords(myCoords);
    myCoords->decrRef();
    double *ptr=myCoords->getPointer();
    std::copy(triofield._coords,triofield._coords+triofield._space_dim*triofield._nbnodes,ptr);
    _local_mesh->allocateCells(triofield._nb_elems);
    INTERP_KERNEL::NormalizedCellType elemtype;
    switch (triofield._mesh_dim)
      {
      case 2:
        switch (triofield._nodes_per_elem)
          {
          case 3:
            elemtype=INTERP_KERNEL::NORM_TRI3;
            break;
          case 4 : 
            elemtype=INTERP_KERNEL::NORM_QUAD4;
            break;
          default:
            throw INTERP_KERNEL::Exception("incompatible Trio field - wrong nb of nodes per elem");
          }
        break;
      case 3:
        switch (triofield._nodes_per_elem)
          {
          case 4:
            elemtype=INTERP_KERNEL::NORM_TETRA4;
            break;
          case 8 : 
            elemtype=INTERP_KERNEL::NORM_HEXA8;
            break;
          default:
            throw INTERP_KERNEL::Exception("incompatible Trio field - wrong nb of nodes per elem");
          }
        break;
      default:
        throw INTERP_KERNEL::Exception("incompatible Trio field - wrong mesh dimension");
      }
    //creating a connectivity table that complies to MED (1 indexing)
    //and passing it to _local_mesh
    int* conn=new int[triofield._nodes_per_elem];
    for (int i=0; i<triofield._nb_elems;i++)
    {
      for(int j=0;j<triofield._nodes_per_elem;j++)
        {
          conn[j]=(triofield._connectivity)[i*triofield._nodes_per_elem+j];
	}
      _local_mesh->insertNextCell(elemtype,triofield._nodes_per_elem,conn);
    }
    delete[] conn;
    
    _local_mesh->setMeshDimension(triofield._mesh_dim);
    _local_mesh->finishInsertingCells();
    _mesh=new ParaMEDMEM::ParaMESH(_local_mesh, group, "support for trio field");
    //_support=new ParaMEDMEM::UnstructuredParaSUPPORT(_local_support,group);
    _comp_topology=new ParaMEDMEM::ComponentTopology(triofield._nb_field_components);
    
    //field on the sending end
    if (triofield._field!=0)
      {
        _field =  new ParaMEDMEM::ParaFIELD(ParaMEDMEM::ON_CELLS,ParaMEDMEM::ONE_TIME,_mesh, *_comp_topology );
        ParaMEDMEM::DataArrayDouble *fieldArr=_field->getField()->getArray();
        _field->getField()->setName(triofield.getName().c_str());
        _field->getField()->setTime(triofield._time1,0,triofield._itnumber);
        for (int i =0; i<triofield._nb_elems; i++)
          for (int j=0; j<triofield._nb_field_components; j++)
            {
              fieldArr->setIJ(i,j,triofield._field[i*triofield._nb_field_components+j]);
            }
      }
    //field on the receiving end
    else
      {
     //   _field =  new ParaMEDMEM::ParaFIELD(ParaMEDMEM::ON_CELLS,_support, *_comp_topology );
        _field =  new ParaMEDMEM::ParaFIELD(ParaMEDMEM::ON_CELLS,ParaMEDMEM::ONE_TIME,_mesh, *_comp_topology );
        _field->getField()->setName(triofield.getName().c_str());
        _field->getField()->setTime(triofield._time1,0,triofield._itnumber);
        // the trio field points to the pointer inside the MED field
        triofield._field=const_cast<double*> (_field->getField()->getArray()->getPointer());
        for (int i=0; i<triofield._nb_field_components*triofield._nb_elems;i++)
          triofield._field[i]=0.0;
      }
    _local_mesh->decrRef();
  }

  MEDField::~MEDField()
  {
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
