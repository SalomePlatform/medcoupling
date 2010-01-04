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

#include <cmath>
#include <limits>
#include <numeric>

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

void MEDCouplingPointSet::zipCoords()
{
  checkFullyDefined();
  DataArrayInt *traducer=zipCoordsTraducer();
  traducer->decrRef();
}

void MEDCouplingPointSet::rotate(const double *center, const double *vector, double angle)
{
  int spaceDim=getSpaceDimension();
  if(spaceDim==3)
    rotate3D(center,vector,angle);
  else if(spaceDim==2)
    rotate2D(center,angle);
  else
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::rotate : invalid space dim for rotation must be 2 or 3");
  _coords->declareAsNew();
  updateTime();
}

void MEDCouplingPointSet::translate(const double *vector)
{
  double *coords=_coords->getPointer();
  int nbNodes=getNumberOfNodes();
  int dim=getSpaceDimension();
  for(int i=0; i<nbNodes; i++)
    for(int idim=0; idim<dim;idim++)
      coords[i*dim+idim]+=vector[idim];
  _coords->declareAsNew();
  updateTime();
}

void MEDCouplingPointSet::tryToShareSameCoords(MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception)
{
  if(_coords==other._coords)
    return ;
  if(!_coords)
    throw INTERP_KERNEL::Exception("Current instance has no coords whereas other has !");
  if(!other._coords)
    throw INTERP_KERNEL::Exception("Other instance has no coords whereas current has !");
  if(!_coords->isEqual(*other._coords,epsilon))
    throw INTERP_KERNEL::Exception("Coords are not the same !");
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
  double* bbtemp = new double[2*dim];
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
      if (!intersects)
        {
          delete [] bbtemp;
          return false; 
        }
    }
  delete [] bbtemp;
  return true;
}

/*!
 * 'This' is expected to be of spaceDim==3. Idem for 'center' and 'vect'
 */
void MEDCouplingPointSet::rotate3D(const double *center, const double *vect, double angle)
{
  double sina=sin(angle);
  double cosa=cos(angle);
  double vectorNorm[3];
  double matrix[9];
  double matrixTmp[9];
  double norm=sqrt(vect[0]*vect[0]+vect[1]*vect[1]+vect[2]*vect[2]);
  std::transform(vect,vect+3,vectorNorm,std::bind2nd(std::multiplies<double>(),1/norm));
  //rotation matrix computation
  matrix[0]=cosa; matrix[1]=0.; matrix[2]=0.; matrix[3]=0.; matrix[4]=cosa; matrix[5]=0.; matrix[6]=0.; matrix[7]=0.; matrix[8]=cosa;
  matrixTmp[0]=vectorNorm[0]*vectorNorm[0]; matrixTmp[1]=vectorNorm[0]*vectorNorm[1]; matrixTmp[2]=vectorNorm[0]*vectorNorm[2];
  matrixTmp[3]=vectorNorm[1]*vectorNorm[0]; matrixTmp[4]=vectorNorm[1]*vectorNorm[1]; matrixTmp[5]=vectorNorm[1]*vectorNorm[2];
  matrixTmp[6]=vectorNorm[2]*vectorNorm[0]; matrixTmp[7]=vectorNorm[2]*vectorNorm[1]; matrixTmp[8]=vectorNorm[2]*vectorNorm[2];
  std::transform(matrixTmp,matrixTmp+9,matrixTmp,std::bind2nd(std::multiplies<double>(),1-cosa));
  std::transform(matrix,matrix+9,matrixTmp,matrix,std::plus<double>());
  matrixTmp[0]=0.; matrixTmp[1]=-vectorNorm[2]; matrixTmp[2]=vectorNorm[1];
  matrixTmp[3]=vectorNorm[2]; matrixTmp[4]=0.; matrixTmp[5]=-vectorNorm[0];
  matrixTmp[6]=-vectorNorm[1]; matrixTmp[7]=vectorNorm[0]; matrixTmp[8]=0.;
  std::transform(matrixTmp,matrixTmp+9,matrixTmp,std::bind2nd(std::multiplies<double>(),sina));
  std::transform(matrix,matrix+9,matrixTmp,matrix,std::plus<double>());
  //rotation matrix computed.
  double *coords=_coords->getPointer();
  int nbNodes=getNumberOfNodes();
  double tmp[3];
  for(int i=0; i<nbNodes; i++)
    {
      std::transform(coords+i*3,coords+(i+1)*3,center,tmp,std::minus<double>());
      coords[i*3]=matrix[0]*tmp[0]+matrix[1]*tmp[1]+matrix[2]*tmp[2]+center[0];
      coords[i*3+1]=matrix[3]*tmp[0]+matrix[4]*tmp[1]+matrix[5]*tmp[2]+center[1];
      coords[i*3+2]=matrix[6]*tmp[0]+matrix[7]*tmp[1]+matrix[8]*tmp[2]+center[2];
    }
}

/*!
 * 'This' is expected to be of spaceDim==2. Idem for 'center' and 'vect'
 */
void MEDCouplingPointSet::rotate2D(const double *center, double angle)
{
  double cosa=cos(angle);
  double sina=sin(angle);
  double matrix[4];
  matrix[0]=cosa; matrix[1]=-sina; matrix[2]=sina; matrix[3]=cosa;
  double *coords=_coords->getPointer();
  int nbNodes=getNumberOfNodes();
  double tmp[2];
  for(int i=0; i<nbNodes; i++)
    {
      std::transform(coords+i*2,coords+(i+1)*2,center,tmp,std::minus<double>());
      coords[i*2]=matrix[0]*tmp[0]+matrix[1]*tmp[1]+center[0];
      coords[i*2+1]=matrix[2]*tmp[0]+matrix[3]*tmp[1]+center[1];
    }
}
