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

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "InterpKernelException.hxx"

namespace ParaMEDMEM
{
  typedef enum
    {
      UNSTRUCTURED = 5,
      UNSTRUCTURED_DESC = 6,
      CARTESIAN = 7,
      EXTRUDED = 8
    } MEDCouplingMeshType;

  class DataArrayDouble;
  class MEDCouplingFieldDouble;

  class MEDCOUPLING_EXPORT MEDCouplingMesh : public RefCountObject, public TimeLabel
  {
  public:
    void setName(const char *name) { _name=name; }
    const char *getName() const { return _name.c_str(); }
    virtual MEDCouplingMeshType getType() const = 0;
    virtual bool isEqual(const MEDCouplingMesh *other, double prec) const { return _name==other->_name; }
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception) = 0;
    virtual bool isStructured() const = 0;
    virtual int getNumberOfCells() const = 0;
    virtual int getNumberOfNodes() const = 0;
    virtual int getSpaceDimension() const = 0;
    virtual int getMeshDimension() const = 0;
    virtual DataArrayDouble *getCoordinatesAndOwner() const = 0;
    virtual DataArrayDouble *getBarycenterAndOwner() const = 0;
    // tools
    virtual void getBoundingBox(double *bbox) const = 0;
    virtual MEDCouplingFieldDouble *getMeasureField(bool isAbs) const = 0;
    virtual MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const = 0;
    virtual int getElementContainingPoint(const double *pos, double eps) const = 0;
    virtual MEDCouplingFieldDouble *fillFromAnalytic(TypeOfField t, int nbOfComp, FunctionToEvaluate func) const;
    virtual MEDCouplingFieldDouble *fillFromAnalytic(TypeOfField t, int nbOfComp, const char *func) const;
    virtual MEDCouplingFieldDouble *buildOrthogonalField() const = 0;
    virtual void rotate(const double *center, const double *vector, double angle) = 0;
    virtual void translate(const double *vector) = 0;
    virtual MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const = 0;
    virtual bool areCompatible(const MEDCouplingMesh *other) const;
    static MEDCouplingMesh *mergeMeshes(const MEDCouplingMesh *mesh1, const MEDCouplingMesh *mesh2);
  protected:
    MEDCouplingMesh() { }
    MEDCouplingMesh(const MEDCouplingMesh& other):_name(other._name) { }
    virtual ~MEDCouplingMesh() { }
  private:
    std::string _name;
  };
}

#endif
