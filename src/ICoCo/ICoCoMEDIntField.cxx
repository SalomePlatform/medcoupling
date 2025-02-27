// Copyright (C) 2007-2025  CEA, EDF
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

#include "ICoCoMEDIntField.h"
#include "MEDCouplingFieldInt32.hxx"

namespace ICoCo
{

  MEDIntField::MEDIntField() : _field(0) {}

  /*! Constructor directly attaching a MEDCouplingFieldInt
    the object does not take the control the objects pointed by
    \a field.
   */
  MEDIntField::MEDIntField(MEDCoupling::MEDCouplingFieldInt32 *field):_field(field)
  {
    if(_field)
      {
        _field->incrRef();
        setName(_field->getName());
      }
    else
      setName("");
  }

  MEDIntField::MEDIntField(const MEDIntField& field):_field(field.getMCField())
  {
    if(_field)
      _field->incrRef();
    setName(field.getName());
  }

  MEDIntField::~MEDIntField()
  {
    if(_field)
      _field->decrRef();
  }


  MEDIntField& MEDIntField::operator=(const MEDIntField& field)
  {
    if (_field)
      _field->decrRef();

    _field=field.getMCField();
    if(_field)
      _field->incrRef();
    setName(field.getName());
    return *this;
  }

  MEDCoupling::MEDCouplingFieldInt32 *MEDIntField::getMCField() const
  {
    return _field;
  }

  void MEDIntField::setMCField(MEDCoupling::MEDCouplingFieldInt32 * f)
  {
    if(_field)
      _field->decrRef();
    _field = f;
    if(f != nullptr)
      {
        _field->incrRef();
        setName(_field->getName());
      }
    else
      setName("");
  }

}
