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
#include "MEDCouplingPointSet.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingUMeshDesc.hxx"
#include "MEDCouplingMemArray.hxx"

using namespace ParaMEDMEM;

MEDCouplingPointSet::MEDCouplingPointSet():_coords(0)
{
}

MEDCouplingPointSet::MEDCouplingPointSet(const MEDCouplingPointSet& other, bool deepCpy):MEDCouplingMesh(other),_coords(0)
{
  if(other._coords)
    _coords=other._coords->performCpy(deepCpy);
}

MEDCouplingPointSet::~MEDCouplingPointSet()
{
  if(_coords)
    _coords->decrRef();
}

int MEDCouplingPointSet::getNumberOfNodes() const
{
  if(_coords)
    return _coords->getNumberOfTuples();
  else
    throw INTERP_KERNEL::Exception("Unable to get number of nodes because no coordinates specified !");
}

int MEDCouplingPointSet::getSpaceDimension() const
{
  if(_coords)
    return _coords->getNumberOfComponents();
  else
    throw INTERP_KERNEL::Exception("Unable to get space dimension because no coordinates specified !");
}

void MEDCouplingPointSet::updateTime()
{
  if(_coords)
    {
      updateTimeWith(*_coords);
    }
}

bool MEDCouplingPointSet::isStructured() const
{
  return false;
}

void MEDCouplingPointSet::setCoords(DataArrayDouble *coords)
{
  if( coords != _coords )
    {
      if (_coords)
        _coords->decrRef();
      _coords=coords;
      if(_coords)
        _coords->incrRef();
      declareAsNew();
    }
}

bool MEDCouplingPointSet::areCoordsEqual(const MEDCouplingPointSet& other, double prec) const
{
  if(_coords==0 && other._coords==0)
    return true;
  if(_coords==0 || other._coords==0)
    return false;
  if(_coords==other._coords)
    return true;
  return _coords->isEqual(*other._coords,prec);
}

void MEDCouplingPointSet::getBoundingBox(double *bbox) const
{
  int dim=getSpaceDimension();
  for (int idim=0; idim<dim; idim++)
    {
      bbox[idim*2]=std::numeric_limits<double>::max();
      bbox[idim*2+1]=-std::numeric_limits<double>::max();
    } 
  const double *coords=_coords->getConstPointer();
  int nbnodes=getNumberOfNodes();
  for (int i=0; i<nbnodes; i++)
    {
      for (int idim=0; idim<dim;idim++)
        {
          if ( bbox[idim*2] > coords[i*dim+idim] )
            {
              bbox[idim*2] = coords[i*dim+idim] ;
            }
          if ( bbox[idim*2+1] < coords[i*dim+idim] )
            {
              bbox[idim*2+1] = coords[i*dim+idim] ;
            }
        }
    }
}

MEDCouplingPointSet *MEDCouplingPointSet::buildInstanceFromMeshType(MEDCouplingMeshType type)
{
  switch(type)
    {
    case UNSTRUCTURED:
      return MEDCouplingUMesh::New();
    case UNSTRUCTURED_DESC:
      return MEDCouplingUMeshDesc::New();
    default:
      throw INTERP_KERNEL::Exception("Invalid type of mesh specified");
    }
}

void MEDCouplingPointSet::getTinySerializationInformation(std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const
{
  if(_coords)
    {
      int spaceDim=getSpaceDimension();
      littleStrings.resize(spaceDim+1);
      littleStrings[0]=getName();
      for(int i=0;i<spaceDim;i++)
        littleStrings[i+1]=getCoords()->getInfoOnComponent(i);
      tinyInfo.clear();
      tinyInfo.push_back(getType());
      tinyInfo.push_back(spaceDim);
      tinyInfo.push_back(getNumberOfNodes());
    }
  else
    {
      littleStrings.resize(1);
      littleStrings[0]=getName();
      tinyInfo.clear();
      tinyInfo.push_back(getType());
      tinyInfo.push_back(-1);
      tinyInfo.push_back(-1);
    }
}

void MEDCouplingPointSet::serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const
{
  if(_coords)
    {
      a2=getCoords();
      a2->incrRef();
    }
  else
    a2=0;
}

void MEDCouplingPointSet::resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings)
{
  if(tinyInfo[2]>=0 && tinyInfo[1]>=1)
    {
      a2->alloc(tinyInfo[2],tinyInfo[1]);
      littleStrings.resize(tinyInfo[1]+1);
    }
  else
    {
      littleStrings.resize(1);
    }
}

void MEDCouplingPointSet::unserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, const std::vector<std::string>& littleStrings)
{
  if(tinyInfo[2]>=0 && tinyInfo[1]>=1)
    {
      setCoords(a2);
      setName(littleStrings[0].c_str());
      for(int i=0;i<tinyInfo[1];i++)
        getCoords()->setInfoOnComponent(i,littleStrings[i+1].c_str());
    }
  else
    setName(littleStrings[0].c_str());
}

// =============================================
// Intersect Bounding Box given 2 Bounding Boxes
// =============================================
bool MEDCouplingPointSet::intersectsBoundingBox(const double* bb1, const double* bb2, int dim, double eps)
{
  double bbtemp[2*dim];
  double deltamax=0.0;

  for (int i=0; i< dim; i++)
    {
      double delta = bb1[2*i+1]-bb1[2*i];
      if ( delta > deltamax )
        {
          deltamax = delta ;
        }
    }
  for (int i=0; i<dim; i++)
    {
      bbtemp[i*2]=bb1[i*2]-deltamax*eps;
      bbtemp[i*2+1]=bb1[i*2+1]+deltamax*eps;
    }
  
  for (int idim=0; idim < dim; idim++)
    {
      bool intersects = (bbtemp[idim*2]<bb2[idim*2+1])
        && (bb2[idim*2]<bbtemp[idim*2+1]) ;
      if (!intersects) return false; 
    }
  return true;
}

