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
#include "MEDCouplingPointSet.txx"
#include "PlanarIntersector.txx"
#include "InterpKernelGeo2DQuadraticPolygon.hxx"
#include "InterpKernelGeo2DNode.hxx"

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

DataArrayDouble *MEDCouplingPointSet::getCoordinatesAndOwner() const
{
  if(_coords)
    _coords->incrRef();
  return _coords;
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

/*!
 * @param comm out parameter (not inout)
 * @param commIndex out parameter (not inout)
 */
void MEDCouplingPointSet::findCommonNodes(DataArrayInt *&comm, DataArrayInt *&commIndex, double prec) const
{
  comm=DataArrayInt::New();
  commIndex=DataArrayInt::New();
  //
  int nbNodesOld=getNumberOfNodes();
  int spaceDim=getSpaceDimension();
  std::vector<double> bbox(2*nbNodesOld*spaceDim);
  const double *coordsPtr=_coords->getConstPointer();
  for(int i=0;i<nbNodesOld;i++)
    {
      for(int j=0;j<spaceDim;j++)
        {
          bbox[2*spaceDim*i+2*j]=coordsPtr[spaceDim*i+j];
          bbox[2*spaceDim*i+2*j+1]=coordsPtr[spaceDim*i+j];
        }
    }
  //
  std::vector<int> c,cI(1);
  switch(spaceDim)
    {
    case 3:
      findCommonNodesAlg<3>(bbox,nbNodesOld,prec,c,cI);
      break;
    case 2:
      findCommonNodesAlg<2>(bbox,nbNodesOld,prec,c,cI);
      break;
    case 1:
      findCommonNodesAlg<1>(bbox,nbNodesOld,prec,c,cI);
      break;
    default:
      throw INTERP_KERNEL::Exception("Unexpected spacedim of coords. Must be 1,2 or 3.");
    }
  commIndex->alloc(cI.size(),1);
  std::copy(cI.begin(),cI.end(),commIndex->getPointer());
  comm->alloc(cI.back(),1);
  std::copy(c.begin(),c.end(),comm->getPointer());
}

/*!
 * @param comm in param in the same format than one returned by findCommonNodes method.
 * @param commI in param in the same format than one returned by findCommonNodes method.
 * @return the old to new correspondance array.
 */
DataArrayInt *MEDCouplingPointSet::buildNewNumberingFromCommNodesFrmt(const DataArrayInt *comm, const DataArrayInt *commIndex,
                                                                      int& newNbOfNodes) const
{
  DataArrayInt *ret=DataArrayInt::New();
  int nbNodesOld=getNumberOfNodes();
  ret->alloc(nbNodesOld,1);
  std::fill(ret->getPointer(),ret->getPointer()+nbNodesOld,-1);
  int *retPtr=ret->getPointer();
  std::vector<int> commRemain(comm->getConstPointer(),comm->getConstPointer()+comm->getNumberOfTuples());
  std::vector<int> commIRemain(commIndex->getConstPointer(),commIndex->getConstPointer()+commIndex->getNumberOfTuples());
  int newNb=0;
  for(int iNode=0;iNode<nbNodesOld;iNode++)
    {
      if(retPtr[iNode]!=-1)
        continue;
      if(commRemain.empty())
        {
          retPtr[iNode]=newNb++;
          continue;
        }
      if(commRemain[0]!=iNode)
        retPtr[iNode]=newNb;
      else
        {
          for(std::vector<int>::const_iterator iNode2=commRemain.begin();
              iNode2!=commRemain.begin()+commIRemain[1];iNode2++)
            retPtr[*iNode2]=newNb;
          int delta=commIRemain[1];
          commRemain.erase(commRemain.begin(),commRemain.begin()+commIRemain[1]);
          commIRemain.erase(commIRemain.begin());
          std::transform(commIRemain.begin(),commIRemain.end(),commIRemain.begin(),std::bind2nd(std::minus<int>(),delta));
        }
      newNb++;
    }
  newNbOfNodes=newNb;
  return ret;
}

void MEDCouplingPointSet::renumberNodes(const int *newNodeNumbers, int newNbOfNodes)
{
  DataArrayDouble *newCoords=DataArrayDouble::New();
  int spaceDim=getSpaceDimension();
  newCoords->alloc(newNbOfNodes,spaceDim);
  newCoords->copyStringInfoFrom(*_coords);
  int oldNbOfNodes=getNumberOfNodes();
  double *ptToFill=newCoords->getPointer();
  const double *oldCoordsPtr=_coords->getConstPointer();
  for(int i=0;i<oldNbOfNodes;i++)
    std::copy(oldCoordsPtr+i*spaceDim,oldCoordsPtr+(i+1)*spaceDim,ptToFill+newNodeNumbers[i]*spaceDim);
  setCoords(newCoords);
  newCoords->decrRef();
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

DataArrayDouble *MEDCouplingPointSet::mergeNodesArray(const MEDCouplingPointSet *m1, const MEDCouplingPointSet *m2)
{
  int spaceDim=m1->getSpaceDimension();
  if(spaceDim!=m2->getSpaceDimension())
    throw INTERP_KERNEL::Exception("Mismatch in SpaceDim during call of mergeNodesArray !");
  return DataArrayDouble::aggregate(m1->getCoords(),m2->getCoords());
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
  double *coords=_coords->getPointer();
  int nbNodes=getNumberOfNodes();
  rotate3DAlg(center,vect,angle,nbNodes,coords);
}

void MEDCouplingPointSet::rotate3DAlg(const double *center, const double *vect, double angle, int nbNodes, double *coords)
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
  double *coords=_coords->getPointer();
  int nbNodes=getNumberOfNodes();
  rotate2DAlg(center,angle,nbNodes,coords);
}

void MEDCouplingPointSet::rotate2DAlg(const double *center, double angle, int nbNodes, double *coords)
{
  double cosa=cos(angle);
  double sina=sin(angle);
  double matrix[4];
  matrix[0]=cosa; matrix[1]=-sina; matrix[2]=sina; matrix[3]=cosa;
  double tmp[2];
  for(int i=0; i<nbNodes; i++)
    {
      std::transform(coords+i*2,coords+(i+1)*2,center,tmp,std::minus<double>());
      coords[i*2]=matrix[0]*tmp[0]+matrix[1]*tmp[1]+center[0];
      coords[i*2+1]=matrix[2]*tmp[0]+matrix[3]*tmp[1]+center[1];
    }
}

class DummyClsMCPS
{
public:
  static const int MY_SPACEDIM=3;
  static const int MY_MESHDIM=2;
  typedef int MyConnType;
  static const INTERP_KERNEL::NumberingPolicy My_numPol=INTERP_KERNEL::ALL_C_MODE;
};

void MEDCouplingPointSet::project2DCellOnXY(const int *startConn, const int *endConn, std::vector<double>& res) const
{
  const double *coords=_coords->getConstPointer();
  int spaceDim=getSpaceDimension();
  for(const int *it=startConn;it!=endConn;it++)
    res.insert(res.end(),coords+spaceDim*(*it),coords+spaceDim*(*it+1));
  if(spaceDim==2)
    return ;
  if(spaceDim==3)
    {
      std::vector<double> cpy(res);
      int nbNodes=endConn-startConn;
      INTERP_KERNEL::PlanarIntersector<DummyClsMCPS,int>::projection(&res[0],&cpy[0],nbNodes,nbNodes,1.e-12,0.,0.,true);
      res.resize(2*nbNodes);
      for(int i=0;i<nbNodes;i++)
        {
          res[2*i]=cpy[3*i];
          res[2*i+1]=cpy[3*i+1];
        }
      return ;
    }
  throw INTERP_KERNEL::Exception("Invalid spacedim for project2DCellOnXY !");
}

bool MEDCouplingPointSet::isButterfly2DCell(const std::vector<double>& res, bool isQuad)
{
  int nbOfNodes=res.size()/2;
  std::vector<INTERP_KERNEL::Node *> nodes(nbOfNodes);
  for(int i=0;i<nbOfNodes;i++)
    {
      INTERP_KERNEL::Node *tmp=new INTERP_KERNEL::Node(res[2*i],res[2*i+1]);
      nodes[i]=tmp;
    }
  INTERP_KERNEL::QuadraticPolygon *pol=0;
  if(isQuad)
    pol=INTERP_KERNEL::QuadraticPolygon::buildArcCirclePolygon(nodes);
  else
    pol=INTERP_KERNEL::QuadraticPolygon::buildLinearPolygon(nodes);
  bool ret=pol->isButterfly();
  delete pol;
  return ret;
}
