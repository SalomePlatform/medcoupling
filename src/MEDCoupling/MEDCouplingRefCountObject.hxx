// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#ifndef __PARAMEDMEM_MEDCOUPLINGREFCOUNTOBJECT_HXX__
#define __PARAMEDMEM_MEDCOUPLINGREFCOUNTOBJECT_HXX__

#include "MEDCoupling.hxx"

namespace ParaMEDMEM
{
  typedef enum
    {
      C_DEALLOC = 2,
      CPP_DEALLOC = 3
    } DeallocType;

  typedef enum
    {
      ON_CELLS = 0,
      ON_NODES = 1,
      ON_GAUSS_PT = 2,
      ON_GAUSS_NE = 3
    } TypeOfField;

  typedef enum
    {
      NO_TIME = 4,
      ONE_TIME = 5,
      LINEAR_TIME = 6,
      CONST_ON_TIME_INTERVAL = 7
    } TypeOfTimeDiscretization;

  typedef bool (*FunctionToEvaluate)(const double *pos, double *res);

  class MEDCOUPLING_EXPORT RefCountObject
  {
  protected:
    RefCountObject();
    RefCountObject(const RefCountObject& other);
  public:
    bool decrRef() const;
    void incrRef() const;
  protected:
    virtual ~RefCountObject();
  private:
    mutable int _cnt;
  };
}

#endif
