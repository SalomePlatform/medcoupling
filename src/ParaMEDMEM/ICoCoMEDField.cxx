// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#include "ICoCoMEDField.hxx"
#include "ICoCoTrioField.hxx"
#include "ProcessorGroup.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "NormalizedUnstructuredMesh.hxx"

namespace ICoCo
{

  /*! Constructor directly attaching a MEDCouplingUMesh and a MEDCouplingFieldDouble
    the object does not take the control the objects pointed by 
    \a mesh and \a field.
  */
    
  MEDField::MEDField(ParaMEDMEM::MEDCouplingUMesh* mesh, ParaMEDMEM::MEDCouplingFieldDouble* field): 
    _mesh(mesh),
    _field(field)
  {
    if(_mesh)
      _mesh->incrRef();
    if(_field)
      _field->incrRef();
  }
  
  MEDField::MEDField(TrioField& triofield)
  {
    _mesh = ParaMEDMEM::MEDCouplingUMesh::New();
    _mesh->setMeshDimension(triofield._space_dim);
    ParaMEDMEM::DataArrayDouble *myCoords=ParaMEDMEM::DataArrayDouble::New();
    myCoords->alloc(triofield._nbnodes,triofield._space_dim);
    _mesh->setCoords(myCoords);
    myCoords->decrRef();
    double *ptr=myCoords->getPointer();
    std::copy(triofield._coords,triofield._coords+triofield._space_dim*triofield._nbnodes,ptr);
    _mesh->allocateCells(triofield._nb_elems);
    INTERP_KERNEL::NormalizedCellType elemtype;
    switch (triofield._mesh_dim)
      {
      case 1:
        switch (triofield._nodes_per_elem)
          {
          case 2:
            elemtype=INTERP_KERNEL::NORM_SEG2;
            break;
          default:
            throw INTERP_KERNEL::Exception("incompatible Trio field - wrong nb of nodes per elem");
          }
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
    //and passing it to _mesh
    int* conn=new int[triofield._nodes_per_elem];
    _mesh->setMeshDimension(triofield._mesh_dim);
    for (int i=0; i<triofield._nb_elems;i++)
      {
        for(int j=0;j<triofield._nodes_per_elem;j++)
          {
            conn[j]=(triofield._connectivity)[i*triofield._nodes_per_elem+j];
          }
        _mesh->insertNextCell(elemtype,triofield._nodes_per_elem,conn);
      }
    delete[] conn;
    
    _mesh->finishInsertingCells();
    
    //field on the sending end
    int nb_case=triofield.nb_values();
    if (triofield._type==0)
      {
        _field =  ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_CELLS,ParaMEDMEM::ONE_TIME);
      }
    else
      {
        _field =  ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES,ParaMEDMEM::ONE_TIME ); 
      }
    _field->setMesh(_mesh);
    _field->setNature(ParaMEDMEM::ConservativeVolumic);
    ParaMEDMEM::DataArrayDouble *fieldArr=ParaMEDMEM::DataArrayDouble::New();
    fieldArr->alloc(_field->getNumberOfTuples(),triofield._nb_field_components);
    _field->setName(triofield.getName().c_str());
    std::string meshName("SupportOf_"); meshName+=_field->getName();
    _mesh->setName(meshName.c_str());
    _field->setTime(triofield._time1,0,triofield._itnumber);
    if (triofield._field!=0)
      {
        for (int i =0; i<nb_case; i++)
          for (int j=0; j<triofield._nb_field_components; j++)
            {
              fieldArr->setIJ(i,j,triofield._field[i*triofield._nb_field_components+j]);
            }
      }
    //field on the receiving end
    else
      {
        // the trio field points to the pointer inside the MED field
        triofield._field=const_cast<double*> (fieldArr->getPointer());
        for (int i=0; i<triofield._nb_field_components*nb_case;i++)
          triofield._field[i]=0.0;
      }
    _field->setArray(fieldArr);
    fieldArr->decrRef();
  }

  MEDField::~MEDField()
  {
    if(_field)
      _field->decrRef();
    if(_mesh)
      _mesh->decrRef();
  }
}
