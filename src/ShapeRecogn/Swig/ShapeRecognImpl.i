// Copyright (C) 2024  CEA, EDF
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

%{
#include <vector>
#include <string>
#include "ShapeRecognMesh.hxx"
#include "ShapeRecognMeshBuilder.hxx"
#include "PrimitiveType.hxx"
#include "Areas.hxx"
 
 #include <type_traits>		
%}

%feature("unref") ShapeRecognMesh "$this->decrRef();"

%newobject Areas::getNodeIds;

%newobject ShapeRecognMeshBuilder::recognize;

%newobject ShapeRecognMesh::getNodeK1;
%newobject ShapeRecognMesh::getNodeK2;
%newobject ShapeRecognMesh::getNodePrimitiveType;
%newobject ShapeRecognMesh::getNodeNormal;
%newobject ShapeRecognMesh::getAreaId;
%newobject ShapeRecognMesh::getAreaPrimitiveType;
%newobject ShapeRecognMesh::getAreaNormal;
%newobject ShapeRecognMesh::getMinorRadius;
%newobject ShapeRecognMesh::getRadius;
%newobject ShapeRecognMesh::getAngle;
%newobject ShapeRecognMesh::getCenter;
%newobject ShapeRecognMesh::getAxis;
%newobject ShapeRecognMesh::getApex;
%newobject ShapeRecognMeshBuilder::buildNodeWeakDirections;
%newobject ShapeRecognMeshBuilder::buildNodeMainDirections;
%newobject ShapeRecognMeshBuilder::buildAreaAxisPoint;
%newobject ShapeRecognMeshBuilder::buildAreaAffinePoint;

%rename (ConvertStringToPrimitive) ConvertStringToPrimitiveSwig;
%rename (ConvertPrimitiveToString) ConvertPrimitiveToStringSwig;

class Areas
{
public:
    double getMinorRadius(mcIdType areaId) const;
    double getRadius(mcIdType areaId) const;
    double getAngle(mcIdType areaId) const;
    //
    bool isEmpty(mcIdType areaId) const;
    size_t getNumberOfAreas() const;
    size_t getNumberOfNodes(mcIdType areaId) const;
    int getPrimitiveType(mcIdType areaId) const;
    std::string getPrimitiveTypeName(mcIdType areaId) const;
    %extend
    {
      DataArrayIdType *getNodeIds(mcIdType areaId) const
      {
        const std::vector<mcIdType> &res = self->getNodeIds(areaId);
        MCAuto< DataArrayIdType > ret( DataArrayIdType::New() );
        ret->alloc(res.size(),1);
        std::copy(res.begin(),res.end(),ret->getPointer());
        return ret.retn();
      }
      
      std::vector<double> getAxis(mcIdType areaId) const
      {
        std::array<double, 3> tmp(self->getAxis(areaId));
        return {tmp.cbegin(),tmp.cend()};
      }

      std::vector<double> getAxisPoint(mcIdType areaId) const
      {
        std::array<double, 3> tmp(self->getAxisPoint(areaId));
        return {tmp.cbegin(),tmp.cend()};
      }

      std::vector<double> getNormal(mcIdType areaId) const
      {
        std::array<double, 3> tmp(self->getNormal(areaId));
        return {tmp.cbegin(),tmp.cend()};
      }
      
      std::vector<double> getAffinePoint(mcIdType areaId) const
      {
        std::array<double, 3> tmp(self->getAffinePoint(areaId));
        return {tmp.cbegin(),tmp.cend()};
      }
      
      std::vector<double> getCenter(mcIdType areaId) const
      {
        std::array<double, 3> tmp(self->getCenter(areaId));
        return {tmp.cbegin(),tmp.cend()};
      }

      std::vector<double> getApex(mcIdType areaId) const
      {
        std::array<double, 3> tmp(self->getApex(areaId));
        return {tmp.cbegin(),tmp.cend()};
      }
    }
private:
    Areas();
    ~Areas();
};

using namespace MEDCoupling;

std::vector<std::string> AllManagedPrimitivesStr();

%inline
{
  std::string ConvertPrimitiveToStringSwig(int type)
  {
    return ConvertPrimitiveToString(static_cast<PrimitiveType>(type));
  }
  
  int ConvertStringToPrimitiveSwig(const std::string& type)
  {
    return static_cast<std::underlying_type_t<PrimitiveType>>( ConvertStringToPrimitive(type) );
  }
}

class ShapeRecognMesh : public RefCountObject
{
public:
    ~ShapeRecognMesh();
    %extend
    {
      // Node properties
      MEDCouplingFieldDouble *getNodeK1() const
      {
          MEDCouplingFieldDouble *ret = const_cast<MEDCouplingFieldDouble *>( self->getNodeK1() );
          ret->incrRef();
          return ret;
      }

      MEDCouplingFieldDouble *getNodeK2() const
      {
          MEDCouplingFieldDouble *ret = const_cast<MEDCouplingFieldDouble *>( self->getNodeK2() );
          ret->incrRef();
          return ret;
      }

      MEDCouplingFieldInt32 *getNodePrimitiveType() const
      {
          MEDCouplingFieldInt32 *ret = const_cast<MEDCouplingFieldInt32 *>( self->getNodePrimitiveType() );
          ret->incrRef();
          return ret;
      }
      
      MEDCouplingFieldDouble *getNodeNormal() const
      {
          MEDCouplingFieldDouble *ret = const_cast<MEDCouplingFieldDouble *>( self->getNodeNormal() );
          ret->incrRef();
          return ret;
      }

      // Area properties

      MEDCouplingFieldInt32 *getAreaId() const
      {
          MEDCouplingFieldInt32 *ret = const_cast<MEDCouplingFieldInt32 *>( self->getAreaId() );
          ret->incrRef();
          return ret;
      }
      
      MEDCouplingFieldInt32 *getAreaPrimitiveType() const
      {
          MEDCouplingFieldInt32 *ret = const_cast<MEDCouplingFieldInt32 *>( self->getAreaPrimitiveType() );
          ret->incrRef();
          return ret;
      }
      
      MEDCouplingFieldDouble *getAreaNormal() const
      {
          MEDCouplingFieldDouble *ret = const_cast<MEDCouplingFieldDouble *>( self->getAreaNormal() );
          ret->incrRef();
          return ret;
      }
      
      MEDCouplingFieldDouble *getMinorRadius() const
      {
          MEDCouplingFieldDouble *ret = const_cast<MEDCouplingFieldDouble *>( self->getMinorRadius() );
          ret->incrRef();
          return ret;
      }
      
      MEDCouplingFieldDouble *getRadius() const
      {
          MEDCouplingFieldDouble *ret = const_cast<MEDCouplingFieldDouble *>( self->getRadius() );
          ret->incrRef();
          return ret;
      }

      MEDCouplingFieldDouble *getAngle() const
      {
          MEDCouplingFieldDouble *ret = const_cast<MEDCouplingFieldDouble *>( self->getAngle() );
          ret->incrRef();
          return ret;
      }
      
      MEDCouplingFieldDouble *getCenter() const
      {
          MEDCouplingFieldDouble *ret = const_cast<MEDCouplingFieldDouble *>( self->getCenter() );
          ret->incrRef();
          return ret;
      }
      
      MEDCouplingFieldDouble *getAxis() const
      {
          MEDCouplingFieldDouble *ret = const_cast<MEDCouplingFieldDouble *>( self->getAxis() );
          ret->incrRef();
          return ret;
      }

      MEDCouplingFieldDouble *getApex() const
      {
          MEDCouplingFieldDouble *ret = const_cast<MEDCouplingFieldDouble *>( self->getApex() );
          ret->incrRef();
          return ret;
      }
    }
private:
    ShapeRecognMesh();
};

class ShapeRecognMeshBuilder
{
public:
    ShapeRecognMeshBuilder(MEDCoupling::MEDCouplingUMesh *mesh);
    ~ShapeRecognMeshBuilder();
    const Areas *getAreas() const;
    %extend
    {
      MEDCouplingFieldDouble *buildNodeWeakDirections() const
      {
          MCAuto<MEDCouplingFieldDouble> ret = self->buildNodeWeakDirections();
          ret->incrRef();
          return ret;
      }

      MEDCouplingFieldDouble *buildNodeMainDirections() const
      {
          MCAuto<MEDCouplingFieldDouble> ret = self->buildNodeMainDirections();
          ret->incrRef();
          return ret;
      }

      MEDCouplingFieldDouble *buildAreaAxisPoint() const
      {
          MCAuto<MEDCouplingFieldDouble> ret = self->buildAreaAxisPoint();
          ret->incrRef();
          return ret;
      }

      MEDCouplingFieldDouble *buildAreaAffinePoint() const
      {
          MCAuto<MEDCouplingFieldDouble> ret = self->buildAreaAffinePoint();
          ret->incrRef();
          return ret;
      }
    
      ShapeRecognMesh *recognize()
      {
          MCAuto<ShapeRecognMesh> ret = self->recognize();
          return ret.retn();
      }
    }
};
