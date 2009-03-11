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
#ifndef __PARAMEDMEM_MEDCOUPLINGFIELD_HXX__
#define __PARAMEDMEM_MEDCOUPLINGFIELD_HXX__

#include "RefCountObject.hxx"
#include "InterpKernelException.hxx"

#include <string>

namespace ParaMEDMEM
{
  class MEDCouplingMesh;

  class MEDCouplingField : public RefCountObject
  {
  public:
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception) = 0;
    void setMesh(ParaMEDMEM::MEDCouplingMesh *mesh);
    void setTime(double val) { _time=val; }
    double getTime() const { return _time; }
    void setDtIt(int dt, int it) { _dt=dt; _it=it; }
    void getDtIt(int& dt, int& it) { dt=_dt; it=_it; }
    ParaMEDMEM::MEDCouplingMesh *getMesh() const { return _mesh; }
    void setName(const char *name) { _name=name; }
    void setDescription(const char *desc) { _desc=desc; }
    const char *getName() const { return _name.c_str(); }
    TypeOfField getEntity() const { return _type; }
  protected:
    void updateTime();
  protected:
    MEDCouplingField(TypeOfField type):_time(0.),_dt(-1),_it(-1),_mesh(0),_type(type) { }
    MEDCouplingField(const MEDCouplingField& other);
    virtual ~MEDCouplingField();
  protected:
    std::string _name;
    std::string _desc;
    double _time;
    int _dt;
    int _it;
    MEDCouplingMesh *_mesh;
    const TypeOfField _type;
  };
}

#endif
