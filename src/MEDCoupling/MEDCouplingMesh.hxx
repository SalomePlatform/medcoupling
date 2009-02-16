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
#ifndef __PARAMEDMEM_MEDCOUPLINGMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGMESH_HXX__

#include "RefCountObject.hxx"
#include "InterpKernelException.hxx"

namespace ParaMEDMEM
{
  class MEDCouplingMesh : public RefCountObject
  {
  public:
    void setName(const char *name) { _name=name; }
    const char *getName() const { return _name.c_str(); }
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception) = 0;
    virtual bool isStructured() const = 0;
    virtual int getNumberOfCells() const = 0;
    virtual int getNumberOfNodes() const = 0;
    virtual int getSpaceDimension() const = 0;
    virtual int getMeshDimension() const = 0;
  protected:
    virtual ~MEDCouplingMesh() { }
  private:
    std::string _name;
  };
}

#endif
