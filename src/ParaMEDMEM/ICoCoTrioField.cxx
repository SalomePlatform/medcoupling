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

// ICoCo file common to several codes
// ICoCoTrioField.cxx
// version 1.2 10/05/2010

#include <ICoCoTrioField.hxx>

#include "ICoCoMEDField.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include <string.h>
#include <iostream>
#include <iomanip>

using namespace ICoCo;
using namespace std;

// Default constructor
TrioField::TrioField() :
  _type(0),
  _mesh_dim(0),
  _space_dim(0),
  _nbnodes(0),
  _nodes_per_elem(0),
  _nb_elems(0),
  _itnumber(0),
  _connectivity(0),
  _coords(0),
  _time1(0.),
  _time2(0.),
  _nb_field_components(0),
  _field(0),
  _has_field_ownership(false) { }

// Copy constructor
TrioField::TrioField(const TrioField& OtherField):_connectivity(0),_coords(0),_field(0) {
  (*this)=OtherField;
}

// Destructor
TrioField::~TrioField() {
  clear();
}

// After the call to clear(), all pointers are null and field ownership is false.
// Arrays are deleted if necessary
void TrioField::clear() {
  if (_connectivity)
    delete[] _connectivity;
  if (_coords)
    delete[] _coords;
  if (_field && _has_field_ownership)
    delete[] _field;
  _connectivity=0;
  _coords=0;
  _field=0;
  _has_field_ownership=false;
}

// Returns the number of value locations
// The size of field is nb_values()*_nb_field_components
int TrioField::nb_values() const {
  if (_type==0)
    return _nb_elems;
  else if (_type==1)
    return _nbnodes;
  throw 0;
  return -1;
}

// Save field to a .field file (loadable by visit!)
void TrioField::save(ostream& os) const{

  os << setprecision(12);
  os << getName() << endl;
  os << _type << endl;
  os << _mesh_dim << endl;
  os << _space_dim << endl;
  os << _nbnodes << endl;
  os << _nodes_per_elem << endl;
  os << _nb_elems << endl;
    
  os<< _itnumber<<endl;
  for (int i=0;i<_nb_elems;i++) {
    for (int j=0;j<_nodes_per_elem;j++)
      os << " " << _connectivity[i*_nodes_per_elem+j];
    os<<endl;
  }
  
  for (int i=0;i<_nbnodes;i++) {
    for (int j=0;j<_space_dim;j++)
      os << " " << _coords[i*_space_dim+j] ;
    os << endl;
  }
  
  os << _time1 << endl;
  os << _time2 << endl;
  os << _nb_field_components << endl;

  if (_field) {
    os << 1 << endl;
    for (int i=0;i<nb_values();i++) {
      for (int j=0;j<_nb_field_components;j++)
        os << " " << _field[i*_nb_field_components+j];
      os << endl;
    }
  }
  else
    os << 0 << endl;

  os << _has_field_ownership << endl;

}

// Restore field from a .field file
void TrioField::restore(istream& in) {

  string name;
  in >> name; 
  setName(name);
  in >> _type;
  in >> _mesh_dim;
  in >> _space_dim;
  in >> _nbnodes;
  in >> _nodes_per_elem;
  in >> _nb_elems;
    
  in >> _itnumber;
  if (_connectivity)
    delete [] _connectivity;
  _connectivity=new int[_nodes_per_elem*_nb_elems];
  for (int i=0;i<_nb_elems;i++) {
    for (int j=0;j<_nodes_per_elem;j++)
      in >> _connectivity[i*_nodes_per_elem+j];
  }
  if (_coords)
    delete [] _coords;
  _coords=new double[_nbnodes*_space_dim];
  for (int i=0;i<_nbnodes;i++) {
    for (int j=0;j<_space_dim;j++)
      in >> _coords[i*_space_dim+j];
  }
  
  in >> _time1;
  in >> _time2;
  in >> _nb_field_components;
  int test;
  in >> test;
  if (test) {
    if (_field)
      delete [] _field;
    _field=new double[_nb_field_components*nb_values()];
    for (int i=0;i<nb_values();i++) {
      for (int j=0;j<_nb_field_components;j++)
        in>> _field[i*_nb_field_components+j];
    }
  }
  else
    _field=0;

  in >> _has_field_ownership;
}


// After the call to set_standalone(), field ownership is true and field is allocated
// to the size _nb_field_components*nb_values().
// The values of the field have been copied if necessary.
void TrioField::set_standalone() {
  if (!_field) {
    _field=new double[_nb_field_components*nb_values()];
    _has_field_ownership=true;
    
  }
  else if (!_has_field_ownership) {
    double *tmp_field=new double[_nb_field_components*nb_values()];
    memcpy(tmp_field,_field,_nb_field_components*nb_values()*sizeof(double));
    _field=tmp_field;
    _has_field_ownership=true;
  }
}

// Used to simulate a 0D geometry (Cathare/Trio for example).
void TrioField::dummy_geom() {
  _type=0;
  _mesh_dim=2;
  _space_dim=2;
  _nbnodes=3;
  _nodes_per_elem=3;
  _nb_elems=1;
  _itnumber=0;
  if (_connectivity)
    delete[] _connectivity;
  _connectivity=new int[3];
  _connectivity[0]=0;
  _connectivity[1]=1;
  _connectivity[2]=2;
  if (_coords)
    delete[] _coords;
  _coords=new double[6];
  _coords[0]=0;
  _coords[1]=0;
  _coords[2]=1;
  _coords[3]=0;
  _coords[4]=0;
  _coords[5]=1;
  _time1=0;
  _time2=1;
  _nb_field_components=1;
  if (_field && _has_field_ownership)
    delete[] _field;
  _has_field_ownership=false;
  _field=0;
}

// Overloading operator = for TrioField
// This becomes an exact copy of NewField.
// If NewField._has_field_ownership is false, they point to the same values.
// Otherwise the values are copied.
TrioField& TrioField::operator=(const TrioField& NewField){

  clear();

  _type=NewField._type;
  _mesh_dim=NewField._mesh_dim;
  _space_dim=NewField._space_dim;
  _nbnodes=NewField._nbnodes;
  _nodes_per_elem=NewField._nodes_per_elem;
  _nb_elems=NewField._nb_elems;
  _itnumber=NewField._itnumber;
  _time1=NewField._time1;
  _time2=NewField._time2;
  _nb_field_components=NewField._nb_field_components;

  if (!NewField._connectivity)
    _connectivity=0;
  else {
    _connectivity=new int[_nodes_per_elem*_nb_elems];
    memcpy( _connectivity,NewField._connectivity,_nodes_per_elem*_nb_elems*sizeof(int));
  }

  if (!NewField._coords)
    _coords=0;
  else {
    _coords=new double[_nbnodes*_space_dim];
    memcpy( _coords,NewField._coords,_nbnodes*_space_dim*sizeof(double));
  }

  //Copie des valeurs du champ
  _has_field_ownership=NewField._has_field_ownership;
  if (_has_field_ownership) {
    _field=new double[nb_values()*_nb_field_components];
    memcpy(_field,NewField._field,nb_values()*_nb_field_components*sizeof(double));
  }
  else
    _field=NewField._field;

  return(*this);

}

/*!
 * This method is non const only due to this->_field that can be modified (to point to the same zone than returned object).
 * So \b warning, to access to \a this->_field only when the returned object is alive.
 */
MEDField *TrioField::build_medfield()
{
  ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingUMesh> mesh(ParaMEDMEM::MEDCouplingUMesh::New("",_mesh_dim));
  ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> coo(ParaMEDMEM::DataArrayDouble::New()); coo->alloc(_nbnodes,_space_dim);
  mesh->setCoords(coo);
  double *ptr(coo->getPointer());
  std::copy(_coords,_coords+_space_dim*_nbnodes,ptr);
  mesh->allocateCells(_nb_elems);
  INTERP_KERNEL::NormalizedCellType elemtype;
  switch(_mesh_dim)
  {
    case 1:
      {
        switch (_nodes_per_elem)
        {
          case 2:
            elemtype=INTERP_KERNEL::NORM_SEG2;
            break;
          default:
            throw INTERP_KERNEL::Exception("incompatible Trio field - wrong nb of nodes per elem");
        }
        break;
      }
    case 2:
      {
        switch (_nodes_per_elem)
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
      }
    case 3:
      {
        switch (_nodes_per_elem)
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
  }
  //creating a connectivity table that complies to MED (1 indexing)
  //and passing it to _mesh
  ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> field;
  int *conn(new int[_nodes_per_elem]);
  for (int i=0; i<_nb_elems;i++)
    {
      for(int j=0;j<_nodes_per_elem;j++)
        {
          conn[j]=_connectivity[i*_nodes_per_elem+j];
        }
      mesh->insertNextCell(elemtype,_nodes_per_elem,conn);
    }
  delete [] conn;
  mesh->finishInsertingCells();
  //
  //field on the sending end
  int nb_case=nb_values();
  if (_type==0)
    {
      field =  ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_CELLS,ParaMEDMEM::ONE_TIME);
    }
  else
    {
      field =  ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES,ParaMEDMEM::ONE_TIME );
    }
  field->setMesh(mesh);
  field->setNature(ParaMEDMEM::ConservativeVolumic);
  ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> fieldArr(ParaMEDMEM::DataArrayDouble::New());
  fieldArr->alloc(field->getNumberOfTuplesExpected(),_nb_field_components);
  field->setName(getName().c_str());
  std::string meshName("SupportOf_"); meshName+=getName();
  mesh->setName(meshName.c_str());
  field->setTime(_time1,0,_itnumber);
  if (_field!=0)
    {
      for (int i =0; i<nb_case; i++)
        for (int j=0; j<_nb_field_components; j++)
          {
            fieldArr->setIJ(i,j,_field[i*_nb_field_components+j]);
          }
    }
  //field on the receiving end
  else
    {
      // the trio field points to the pointer inside the MED field
      _field=fieldArr->getPointer();
      for (int i=0; i<_nb_field_components*nb_case;i++)
        _field[i]=0.0;
    }
  field->setArray(fieldArr);
  return new MEDField(field);
}



