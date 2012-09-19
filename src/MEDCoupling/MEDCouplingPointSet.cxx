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
// Author : Anthony Geay (CEA/DEN)

#include "MEDCouplingPointSet.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingUMeshDesc.hxx"
#include "MEDCouplingMemArray.hxx"
#include "PlanarIntersector.txx"
#include "InterpKernelGeo2DQuadraticPolygon.hxx"
#include "InterpKernelGeo2DNode.hxx"
#include "DirectedBoundingBox.hxx"
#include "InterpKernelAutoPtr.hxx"

#include <cmath>
#include <limits>
#include <numeric>

using namespace ParaMEDMEM;

MEDCouplingPointSet::MEDCouplingPointSet():_coords(0)
{
}

MEDCouplingPointSet::MEDCouplingPointSet(const MEDCouplingPointSet& other, bool deepCopy):MEDCouplingMesh(other),_coords(0)
{
  if(other._coords)
    _coords=other._coords->performCpy(deepCopy);
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

void MEDCouplingPointSet::updateTime() const
{
  if(_coords)
    {
      updateTimeWith(*_coords);
    }
}

void MEDCouplingPointSet::setCoords(const DataArrayDouble *coords)
{
  if( coords != _coords )
    {
      if (_coords)
        _coords->decrRef();
      _coords=const_cast<DataArrayDouble *>(coords);
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

/*!
 * This method copyies all tiny strings from other (name and components name).
 * @throw if other and this have not same mesh type.
 */
void MEDCouplingPointSet::copyTinyStringsFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingPointSet *otherC=dynamic_cast<const MEDCouplingPointSet *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::copyTinyStringsFrom : meshes have not same type !");
  MEDCouplingMesh::copyTinyStringsFrom(other);
  if(_coords && otherC->_coords)
    _coords->copyStringInfoFrom(*otherC->_coords);
}

bool MEDCouplingPointSet::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::isEqualIfNotWhy : null mesh instance in input !");
  const MEDCouplingPointSet *otherC=dynamic_cast<const MEDCouplingPointSet *>(other);
  if(!otherC)
    {
      reason="mesh given in input is not castable in MEDCouplingPointSet !";
      return false;
    }
  if(!MEDCouplingMesh::isEqualIfNotWhy(other,prec,reason))
    return false;
  if(!areCoordsEqualIfNotWhy(*otherC,prec,reason))
    return false;
  return true;
}

bool MEDCouplingPointSet::isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const
{
  const MEDCouplingPointSet *otherC=dynamic_cast<const MEDCouplingPointSet *>(other);
  if(!otherC)
    return false;
  if(!areCoordsEqualWithoutConsideringStr(*otherC,prec))
    return false;
  return true;
}

bool MEDCouplingPointSet::areCoordsEqualIfNotWhy(const MEDCouplingPointSet& other, double prec, std::string& reason) const
{
  if(_coords==0 && other._coords==0)
    return true;
  if(_coords==0 || other._coords==0)
    {
      reason="Only one PointSet between the two this and other has coordinate defined !";
      return false;
    }
  if(_coords==other._coords)
    return true;
  bool ret=_coords->isEqualIfNotWhy(*other._coords,prec,reason);
  if(!ret)
    reason.insert(0,"Coordinates DataArray do not match : ");
  return ret;
}

bool MEDCouplingPointSet::areCoordsEqual(const MEDCouplingPointSet& other, double prec) const
{
  std::string tmp;
  return areCoordsEqualIfNotWhy(other,prec,tmp);
}

bool MEDCouplingPointSet::areCoordsEqualWithoutConsideringStr(const MEDCouplingPointSet& other, double prec) const
{
  if(_coords==0 && other._coords==0)
    return true;
  if(_coords==0 || other._coords==0)
    return false;
  if(_coords==other._coords)
    return true;
  return _coords->isEqualWithoutConsideringStr(*other._coords,prec);
}

/*!
 * Returns coordinates of node with id 'nodeId' and append it in 'coo'.
 */
void MEDCouplingPointSet::getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const throw(INTERP_KERNEL::Exception)
{
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::getCoordinatesOfNode : no coordinates array set !");
  int nbNodes=getNumberOfNodes();
  if(nodeId>=0 && nodeId<nbNodes)
    {
      const double *cooPtr=_coords->getConstPointer();
      int spaceDim=getSpaceDimension();
      coo.insert(coo.end(),cooPtr+spaceDim*nodeId,cooPtr+spaceDim*(nodeId+1));
    }
  else
    {
      std::ostringstream oss; oss << "MEDCouplingPointSet::getCoordinatesOfNode : request of nodeId \"" << nodeId << "\" but it should be in [0,"<< nbNodes << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * This method is typically the base method used for implementation of mergeNodes. This method computes this permutation array using as input,
 * This method is const ! So this method simply computes the array, no permutation of nodes is done.
 * a precision 'precision' and a 'limitNodeId' that is the node id so that every nodes which id is strictly lower than 'limitNodeId' will not be merged.
 * To desactivate this advanced feature put -1 to this argument.
 * @param areNodesMerged output parameter that states if some nodes have been "merged" in returned array
 * @param newNbOfNodes output parameter too this is the maximal id in returned array to avoid to recompute it.
 */
DataArrayInt *MEDCouplingPointSet::buildPermArrayForMergeNode(double precision, int limitNodeId, bool& areNodesMerged, int& newNbOfNodes) const
{
  DataArrayInt *comm,*commI;
  findCommonNodes(precision,limitNodeId,comm,commI);
  int oldNbOfNodes=getNumberOfNodes();
  DataArrayInt *ret=buildNewNumberingFromCommonNodesFormat(comm,commI,newNbOfNodes);
  areNodesMerged=(oldNbOfNodes!=newNbOfNodes);
  comm->decrRef();
  commI->decrRef();
  return ret;
}

/*!
 * This methods searches for each node if there are any nodes in _coords that are less far than 'prec' from n1. if any, these nodes are stored in out params
 * comm and commIndex.
 * @param limitNodeId is the limit node id. All nodes which id is strictly lower than 'limitNodeId' will not be merged each other.
 * @param comm out parameter (not inout)
 * @param commIndex out parameter (not inout)
 */
void MEDCouplingPointSet::findCommonNodes(double prec, int limitNodeId, DataArrayInt *&comm, DataArrayInt *&commIndex) const
{
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::findCommonNodes : no coords specified !");
  _coords->findCommonTuples(prec,limitNodeId,comm,commIndex);
}

std::vector<int> MEDCouplingPointSet::getNodeIdsNearPoint(const double *pos, double eps) const throw(INTERP_KERNEL::Exception)
{
  std::vector<int> c,cI;
  getNodeIdsNearPoints(pos,1,eps,c,cI);
  return c;
}

/*!
 * Given a point given by its position 'pos' this method finds the set of node ids that are a a distance lower than eps.
 * Position 'pos' is expected to be of size getSpaceDimension()*nbOfNodes. If not the behabiour is not warranted.
 * This method throws an exception if no coordiantes are set.
 */
void MEDCouplingPointSet::getNodeIdsNearPoints(const double *pos, int nbOfNodes, double eps, std::vector<int>& c, std::vector<int>& cI) const throw(INTERP_KERNEL::Exception)
{
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::getNodeIdsNearPoint : no coordiantes set !");
  int spaceDim=getSpaceDimension();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> points=DataArrayDouble::New();
  points->useArray(pos,false,CPP_DEALLOC,nbOfNodes,spaceDim);
  _coords->computeTupleIdsNearTuples(points,eps,c,cI);
}

/*!
 * @param comm in param in the same format than one returned by findCommonNodes method.
 * @param commI in param in the same format than one returned by findCommonNodes method.
 * @return the old to new correspondance array.
 */
DataArrayInt *MEDCouplingPointSet::buildNewNumberingFromCommonNodesFormat(const DataArrayInt *comm, const DataArrayInt *commIndex,
                                                                          int& newNbOfNodes) const
{
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::buildNewNumberingFromCommonNodesFormat : no coords specified !");
  return DataArrayInt::BuildOld2NewArrayFromSurjectiveFormat2(getNumberOfNodes(),comm,commIndex,newNbOfNodes);
}

/*
 * This method renumber 'this' using 'newNodeNumbers' array of size this->getNumberOfNodes.
 * newNbOfNodes specifies the *std::max_element(newNodeNumbers,newNodeNumbers+this->getNumberOfNodes())
 * This value is asked because often known by the caller of this method.
 * @param newNodeNumbers array specifying the new numbering in old2New convention..
 * @param newNbOfNodes the new number of nodes.
 */
void MEDCouplingPointSet::renumberNodes(const int *newNodeNumbers, int newNbOfNodes)
{
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::renumberNodes : no coords specified !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> newCoords=_coords->renumberAndReduce(newNodeNumbers,newNbOfNodes);
  setCoords(newCoords);
}

/*
 * This method renumber 'this' using 'newNodeNumbers' array of size this->getNumberOfNodes.
 * newNbOfNodes specifies the *std::max_element(newNodeNumbers,newNodeNumbers+this->getNumberOfNodes())
 * This value is asked because often known by the caller of this method.
 * Contrary to ParaMEDMEM::MEDCouplingPointSet::renumberNodes method for merged nodes the barycenter of them is computed here.
 *
 * @param newNodeNumbers array specifying the new numbering.
 * @param newNbOfNodes the new number of nodes.
 */
void MEDCouplingPointSet::renumberNodes2(const int *newNodeNumbers, int newNbOfNodes)
{
  DataArrayDouble *newCoords=DataArrayDouble::New();
  std::vector<int> div(newNbOfNodes);
  int spaceDim=getSpaceDimension();
  newCoords->alloc(newNbOfNodes,spaceDim);
  newCoords->copyStringInfoFrom(*_coords);
  newCoords->fillWithZero();
  int oldNbOfNodes=getNumberOfNodes();
  double *ptToFill=newCoords->getPointer();
  const double *oldCoordsPtr=_coords->getConstPointer();
  for(int i=0;i<oldNbOfNodes;i++)
    {
      std::transform(oldCoordsPtr+i*spaceDim,oldCoordsPtr+(i+1)*spaceDim,ptToFill+newNodeNumbers[i]*spaceDim,
                     ptToFill+newNodeNumbers[i]*spaceDim,std::plus<double>());
      div[newNodeNumbers[i]]++;
    }
  for(int i=0;i<newNbOfNodes;i++)
    ptToFill=std::transform(ptToFill,ptToFill+spaceDim,ptToFill,std::bind2nd(std::multiplies<double>(),1./(double)div[i]));
  setCoords(newCoords);
  newCoords->decrRef();
}

/*!
 * This method fills bbox params like that : bbox[0]=XMin, bbox[1]=XMax, bbox[2]=YMin...
 * The returned bounding box is arranged along trihedron.
 * @param bbox out array of size 2*this->getSpaceDimension().
 */
void MEDCouplingPointSet::getBoundingBox(double *bbox) const throw(INTERP_KERNEL::Exception)
{
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::getBoundingBox : Coordinates not set !");
  _coords->getMinMaxPerComponent(bbox);
}

/*!
 * This method removes useless nodes in coords.
 */
void MEDCouplingPointSet::zipCoords()
{
  checkFullyDefined();
  DataArrayInt *traducer=zipCoordsTraducer();
  traducer->decrRef();
}

struct MEDCouplingCompAbs
{
  bool operator()(double x, double y) { return std::abs(x)<std::abs(y);}
};

/*!
 * This method expects that _coords attribute is set.
 * @return the carateristic dimension of point set. This caracteristic dimension is the max of difference 
 * @exception If _coords attribute not set.
 */
double MEDCouplingPointSet::getCaracteristicDimension() const
{
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::getCaracteristicDimension : Coordinates not set !");
  const double *coords=_coords->getConstPointer();
  int nbOfValues=_coords->getNbOfElems();
  return std::abs(*std::max_element(coords,coords+nbOfValues,MEDCouplingCompAbs()));
}

/*!
 * This method recenter coordinates of nodes in \b this in order to be centered at the origin to benefit about the advantages of the precision to be around the box
 * around origin of 'radius' 1.
 *
 * \param [in] eps absolute epsilon. under that value of delta between max and min no scale is performed.
 *
 * \warning this method is non const and alterates coordinates in \b this without modifying.
 */
void MEDCouplingPointSet::recenterForMaxPrecision(double eps) throw(INTERP_KERNEL::Exception)
{
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::recenterForMaxPrecision : Coordinates not set !");
  _coords->recenterForMaxPrecision(eps);
  updateTime();
}

/*!
 * Non const method that operates a rotation of 'this'.
 * If spaceDim==2 'vector' parameter is ignored (and could be 0) and the rotation is done around 'center' with angle specified by 'angle'.
 * If spaceDim==3 the rotation axe is defined by ('center','vector') and the angle is 'angle'.
 * @param center an array of size getSpaceDimension().
 * @param vector in array of size getSpaceDimension().
 * @param angle angle of rotation in radian.
 */
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

/*!
 * Non const method that operates a translation of 'this'.
 * @param vector in array of size getSpaceDimension().
 */
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

/*!
 * Non const method that operates a scale on 'this' with 'point' as reference point of scale and with factor 'factor'.
 * @param point in array of size getSpaceDimension().
 * @param factor factor of the scaling
 */
void MEDCouplingPointSet::scale(const double *point, double factor)
{
  double *coords=_coords->getPointer();
  int nbNodes=getNumberOfNodes();
  int dim=getSpaceDimension();
  double *tmp=new double[dim];
  for(int i=0;i<nbNodes;i++)
    {
      std::transform(coords+i*dim,coords+(i+1)*dim,point,coords+i*dim,std::minus<double>());
      std::transform(coords+i*dim,coords+(i+1)*dim,coords+i*dim,std::bind2nd(std::multiplies<double>(),factor));
      std::transform(coords+i*dim,coords+(i+1)*dim,point,coords+i*dim,std::plus<double>());
    }
  delete [] tmp;
  _coords->declareAsNew();
  updateTime();
}

/*!
 * This method is only available for already defined coordinates.
 * If not an INTERP_KERNEL::Exception is thrown. The 'newSpaceDim' input must be greater or equal to 1.
 * This method simply convert this to newSpaceDim space :
 * - by putting a 0. for each \f$ i^{th} \f$ components of each coord of nodes so that i>=getSpaceDim(), if 'newSpaceDim'>getSpaceDimsion()
 * - by ignoring each \f$ i^{th} \f$ components of each coord of nodes so that i >= 'newSpaceDim', 'newSpaceDim'<getSpaceDimension()
 * If newSpaceDim==getSpaceDim() nothing is done by this method.
 */
void MEDCouplingPointSet::changeSpaceDimension(int newSpaceDim, double dftValue) throw(INTERP_KERNEL::Exception)
{
  if(getCoords()==0)
    throw INTERP_KERNEL::Exception("changeSpaceDimension must be called on an MEDCouplingPointSet instance with coordinates set !");
  if(newSpaceDim<1)
    throw INTERP_KERNEL::Exception("changeSpaceDimension must be called a newSpaceDim >=1 !");
  int oldSpaceDim=getSpaceDimension();
  if(newSpaceDim==oldSpaceDim)
    return ;
  DataArrayDouble *newCoords=getCoords()->changeNbOfComponents(newSpaceDim,dftValue);
  setCoords(newCoords);
  newCoords->decrRef();
  updateTime();
}

/*!
 * This method try to substitute this->_coords with other._coords if arrays match.
 * This method potentially modifies 'this' if it succeeds, otherway an exception is thrown.
 */
void MEDCouplingPointSet::tryToShareSameCoords(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception)
{
  if(_coords==other._coords)
    return ;
  if(!_coords)
    throw INTERP_KERNEL::Exception("Current instance has no coords whereas other has !");
  if(!other._coords)
    throw INTERP_KERNEL::Exception("Other instance has no coords whereas current has !");
  if(!_coords->isEqualWithoutConsideringStr(*other._coords,epsilon))
    throw INTERP_KERNEL::Exception("Coords are not the same !");
  setCoords(other._coords);
}

/*!
 * This method duplicates the nodes whose ids are in [\b nodeIdsToDuplicateBg, \b nodeIdsToDuplicateEnd) and put the result of their duplication at the end
 * of existing node ids.
 * 
 * \param [in] nodeIdsToDuplicateBg begin of node ids (included) to be duplicated in connectivity only
 * \param [in] nodeIdsToDuplicateEnd end of node ids (excluded) to be duplicated in connectivity only
 */
void MEDCouplingPointSet::duplicateNodesInCoords(const int *nodeIdsToDuplicateBg, const int *nodeIdsToDuplicateEnd) throw(INTERP_KERNEL::Exception)
{
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::duplicateNodesInCoords : no coords set !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> newCoords=_coords->selectByTupleIdSafe(nodeIdsToDuplicateBg,nodeIdsToDuplicateEnd);
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> newCoords2=DataArrayDouble::Aggregate(_coords,newCoords);
  setCoords(newCoords2);
}

/*!
 * This method is expecting to be called for meshes so that getSpaceDimension() returns 3.
 * This method returns in 'nodes' output all the nodes that are at a distance lower than epsilon from plane
 * defined by the point 'pt' and the vector 'vec'.
 * @param pt points to an array of size 3 and represents a point that owns to plane.
 * @param vec points to an array of size 3 and represents the normal vector of the plane. The norm of the vector is not compulsory equal to 1. But norm must be greater than 10*abs(eps)
 * @param eps is the maximal distance around the plane where node in this->_coords will be picked.
 * @param nodes is the output of the method. The vector is not compulsory empty before call. The nodes that fulfills the condition will be added at the end of the nodes.
 */
void MEDCouplingPointSet::findNodesOnPlane(const double *pt, const double *vec, double eps, std::vector<int>& nodes) const throw(INTERP_KERNEL::Exception)
{
  if(getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::findNodesOnPlane : Invalid spacedim to be applied on this ! Must be equal to 3 !");
  int nbOfNodes=getNumberOfNodes();
  double a=vec[0],b=vec[1],c=vec[2],d=-pt[0]*vec[0]-pt[1]*vec[1]-pt[2]*vec[2];
  double deno=sqrt(a*a+b*b+c*c);
  const double *work=_coords->getConstPointer();
  for(int i=0;i<nbOfNodes;i++)
    {
      if(std::abs(a*work[0]+b*work[1]+c*work[2]+d)/deno<eps)
        nodes.push_back(i);
      work+=3;
    }
}

/*!
 * This method is expecting to be called for meshes so that getSpaceDimension() returns 2 or 3.
 * This method returns in 'nodes' output all the nodes that are at a distance lower than epsilon from a line
 * defined by the point 'pt' and the vector 'vec'. 'pt' and 'vec' are expected to have a dimension equal to space dimension of 'this'
 * @param pt points to an array of size this->getSpaceDimension and represents a point that owns to plane.
 * @param vec points to an array of size this->getSpaceDimension and represents the direction vector of the line. The norm of the vector is not compulsory equal to 1.
 *            But norm must be greater than 10*abs(eps)
 * @param eps is the maximal distance around the plane where node in this->_coords will be picked.
 * @param nodes is the output of the method. The vector is not compulsory empty before call. The nodes that fulfills the condition will be added at the end of the nodes.
 */
void MEDCouplingPointSet::findNodesOnLine(const double *pt, const double *vec, double eps, std::vector<int>& nodes) const throw(INTERP_KERNEL::Exception)
{
  int spaceDim=getSpaceDimension();
  if(spaceDim!=2 && spaceDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::findNodesOnLine : Invalid spacedim to be applied on this ! Must be equal to 2 or 3 !");
  int nbOfNodes=getNumberOfNodes();
  double den=0.;
  for(int i=0;i<spaceDim;i++)
    den+=vec[i]*vec[i];
  double deno=sqrt(den);
  if(deno<10.*eps)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::findNodesOnLine : Invalid given direction vector ! Norm is too small !");
  INTERP_KERNEL::AutoPtr<double> vecn=new double[spaceDim];
  for(int i=0;i<spaceDim;i++)
    vecn[i]=vec[i]/deno;
  const double *work=_coords->getConstPointer();
  if(spaceDim==2)
    {
      for(int i=0;i<nbOfNodes;i++)
        {
          if(std::abs(vecn[0]*(work[1]-pt[1])-vecn[1]*(work[0]-pt[0]))<eps)
            nodes.push_back(i);
          work+=2;
        }
    }
  else
    {
      for(int i=0;i<nbOfNodes;i++)
        {
          double a=vecn[0]*(work[1]-pt[1])-vecn[1]*(work[0]-pt[0]);
          double b=vecn[1]*(work[2]-pt[2])-vecn[2]*(work[1]-pt[1]);
          double c=vecn[2]*(work[0]-pt[0])-vecn[0]*(work[2]-pt[2]);
          if(std::sqrt(a*a+b*b+c*c)<eps)
            nodes.push_back(i);
          work+=3;
        }
    }
}

/*!
 * merge _coords arrays of m1 and m2 and returns the union. The returned instance is newly created with ref count == 1.
 */
DataArrayDouble *MEDCouplingPointSet::MergeNodesArray(const MEDCouplingPointSet *m1, const MEDCouplingPointSet *m2) throw(INTERP_KERNEL::Exception)
{
  int spaceDim=m1->getSpaceDimension();
  if(spaceDim!=m2->getSpaceDimension())
    throw INTERP_KERNEL::Exception("Mismatch in SpaceDim during call of MergeNodesArray !");
  return DataArrayDouble::Aggregate(m1->getCoords(),m2->getCoords());
}

DataArrayDouble *MEDCouplingPointSet::MergeNodesArray(const std::vector<const MEDCouplingPointSet *>& ms) throw(INTERP_KERNEL::Exception)
{
  if(ms.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::MergeNodesArray : input array must be NON EMPTY !");
  std::vector<const MEDCouplingPointSet *>::const_iterator it=ms.begin();
  std::vector<const DataArrayDouble *> coo(ms.size());
  int spaceDim=(*it)->getSpaceDimension();
  coo[0]=(*it++)->getCoords();
  for(int i=1;it!=ms.end();it++,i++)
    {
      const DataArrayDouble *tmp=(*it)->getCoords();
      if(tmp)
        {
          if((*it)->getSpaceDimension()==spaceDim)
            coo[i]=tmp;
          else
            throw INTERP_KERNEL::Exception("Mismatch in SpaceDim during call of MergeNodesArray !");
        }
      else
        throw INTERP_KERNEL::Exception("Empty coords detected during call of MergeNodesArray !");
    }
  return DataArrayDouble::Aggregate(coo);
}

/*!
 * Factory to build new instance of instanciable subclasses of MEDCouplingPointSet.
 * This method is used during unserialization process.
 */
MEDCouplingPointSet *MEDCouplingPointSet::BuildInstanceFromMeshType(MEDCouplingMeshType type)
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

/*!
 * First step of serialization process. Used by ParaMEDMEM and MEDCouplingCorba to transfert data between process.
 */
void MEDCouplingPointSet::getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const
{
  int it,order;
  double time=getTime(it,order);
  if(_coords)
    {
      int spaceDim=getSpaceDimension();
      littleStrings.resize(spaceDim+4);
      littleStrings[0]=getName();
      littleStrings[1]=getDescription();
      littleStrings[2]=_coords->getName();
      littleStrings[3]=getTimeUnit();
      for(int i=0;i<spaceDim;i++)
        littleStrings[i+4]=getCoords()->getInfoOnComponent(i);
      tinyInfo.clear();
      tinyInfo.push_back(getType());
      tinyInfo.push_back(spaceDim);
      tinyInfo.push_back(getNumberOfNodes());
      tinyInfo.push_back(it);
      tinyInfo.push_back(order);
      tinyInfoD.push_back(time);
    }
  else
    {
      littleStrings.resize(3);
      littleStrings[0]=getName();
      littleStrings[1]=getDescription();
      littleStrings[2]=getTimeUnit();
      tinyInfo.clear();
      tinyInfo.push_back(getType());
      tinyInfo.push_back(-1);
      tinyInfo.push_back(-1);
      tinyInfo.push_back(it);
      tinyInfo.push_back(order);
      tinyInfoD.push_back(time);
    }
}

/*!
 * Third and final step of serialization process.
 */
void MEDCouplingPointSet::serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const
{
  if(_coords)
    {
      a2=const_cast<DataArrayDouble *>(getCoords());
      a2->incrRef();
    }
  else
    a2=0;
}

/*!
 * Second step of serialization process.
 * @param tinyInfo must be equal to the result given by getTinySerializationInformation method.
 */
void MEDCouplingPointSet::resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const
{
  if(tinyInfo[2]>=0 && tinyInfo[1]>=1)
    {
      a2->alloc(tinyInfo[2],tinyInfo[1]);
      littleStrings.resize(tinyInfo[1]+4);
    }
  else
    {
      littleStrings.resize(3);
    }
}

/*!
 * Second and final unserialization process.
 * @param tinyInfo must be equal to the result given by getTinySerializationInformation method.
 */
void MEDCouplingPointSet::unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2, const std::vector<std::string>& littleStrings)
{
  if(tinyInfo[2]>=0 && tinyInfo[1]>=1)
    {
      setCoords(a2);
      setName(littleStrings[0].c_str());
      setDescription(littleStrings[1].c_str());
      a2->setName(littleStrings[2].c_str());
      setTimeUnit(littleStrings[3].c_str());
      for(int i=0;i<tinyInfo[1];i++)
        getCoords()->setInfoOnComponent(i,littleStrings[i+4].c_str());
      setTime(tinyInfoD[0],tinyInfo[3],tinyInfo[4]);
    }
  else
    {
      setName(littleStrings[0].c_str());
      setDescription(littleStrings[1].c_str());
      setTimeUnit(littleStrings[2].c_str());
      setTime(tinyInfoD[0],tinyInfo[3],tinyInfo[4]);
    }
}

/*!
 * Intersect Bounding Box given 2 Bounding Boxes.
 */
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
 * Intersect 2 given Bounding Boxes.
 */
bool MEDCouplingPointSet::intersectsBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bb1, const double* bb2, int dim, double eps)
{
  double* bbtemp = new double[2*dim];
  double deltamax=0.0;

  for (int i=0; i< dim; i++)
    {
      double delta = bb2[2*i+1]-bb2[2*i];
      if ( delta > deltamax )
        {
          deltamax = delta ;
        }
    }
  for (int i=0; i<dim; i++)
    {
      bbtemp[i*2]=bb2[i*2]-deltamax*eps;
      bbtemp[i*2+1]=bb2[i*2+1]+deltamax*eps;
    }
  
  bool intersects = !bb1.isDisjointWith( bbtemp );
  delete [] bbtemp;
  return intersects;
}

/*!
 * 'This' is expected to be of spaceDim==3. Idem for 'center' and 'vect'
 */
void MEDCouplingPointSet::rotate3D(const double *center, const double *vect, double angle)
{
  double *coords=_coords->getPointer();
  int nbNodes=getNumberOfNodes();
  Rotate3DAlg(center,vect,angle,nbNodes,coords);
}

/*!
 * Low static method that operates 3D rotation of 'nbNodes' 3D nodes whose coordinates are arranged in 'coords'
 * around an axe ('center','vect') and with angle 'angle'.
 */
void MEDCouplingPointSet::Rotate3DAlg(const double *center, const double *vect, double angle, int nbNodes, double *coords)
{
  if(!center || !vect)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::Rotate3DAlg : null vector in input !");
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
 * This method implements pure virtual method MEDCouplingMesh::buildPart.
 * This method build a part of 'this' by simply keeping cells whose ids are in ['start','end').
 * The coords are kept unchanged contrary to pure virtual method MEDCouplingMesh::buildPartAndReduceNodes.
 * The returned mesh has to be managed by the caller.
 */
MEDCouplingMesh *MEDCouplingPointSet::buildPart(const int *start, const int *end) const
{
  return buildPartOfMySelf(start,end,true);
}

/*!
 * This method implements pure virtual method MEDCouplingMesh::buildPartAndReduceNodes.
 * This method build a part of 'this' by simply keeping cells whose ids are in ['start','end') \b and potentially reduces the nodes set
 * behind returned mesh. This cause an overhead but it is lesser in memory.
 * This method returns an array too. This array allows to the caller to know the mapping between nodeids in 'this' and nodeids in 
 * returned mesh. This is quite usefull for MEDCouplingFieldDouble on nodes for example...
 * 'arr' is in old2New format of size ret->getNumberOfCells like MEDCouplingUMesh::zipCoordsTraducer is.
 * The returned mesh has to be managed by the caller.
 */
MEDCouplingMesh *MEDCouplingPointSet::buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const
{
  MEDCouplingPointSet *ret=buildPartOfMySelf(start,end,true);
  arr=ret->zipCoordsTraducer();
  return ret;
}


/*!
 * 'This' is expected to be of spaceDim==2. Idem for 'center' and 'vect'
 */
void MEDCouplingPointSet::rotate2D(const double *center, double angle)
{
  double *coords=_coords->getPointer();
  int nbNodes=getNumberOfNodes();
  Rotate2DAlg(center,angle,nbNodes,coords);
}

/*!
 * Low static method that operates 3D rotation of 'nbNodes' 3D nodes whose coordinates are arranged in 'coords'
 * around the center point 'center' and with angle 'angle'.
 */
void MEDCouplingPointSet::Rotate2DAlg(const double *center, double angle, int nbNodes, double *coords)
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

/// @cond INTERNAL

class DummyClsMCPS
{
public:
  static const int MY_SPACEDIM=3;
  static const int MY_MESHDIM=2;
  typedef int MyConnType;
  static const INTERP_KERNEL::NumberingPolicy My_numPol=INTERP_KERNEL::ALL_C_MODE;
};

/// @endcond

/*!
 * res should be an empty vector before calling this method.
 * This method returns all the node coordinates included in _coords which ids are in [startConn;endConn) and put it into 'res' vector.
 * If spaceDim==3 a projection will be done for each nodes on the middle plane containing these all nodes in [startConn;endConn).
 * And after each projected nodes are moved to Oxy plane in order to consider these nodes as 2D nodes.
 */
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
      int nbNodes=(int)std::distance(startConn,endConn);
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

/*!
 * low level method that checks that the 2D cell is not a butterfly cell.
 */
bool MEDCouplingPointSet::isButterfly2DCell(const std::vector<double>& res, bool isQuad, double eps)
{
  std::size_t nbOfNodes=res.size()/2;
  std::vector<INTERP_KERNEL::Node *> nodes(nbOfNodes);
  for(std::size_t i=0;i<nbOfNodes;i++)
    {
      INTERP_KERNEL::Node *tmp=new INTERP_KERNEL::Node(res[2*i],res[2*i+1]);
      nodes[i]=tmp;
    }
  INTERP_KERNEL::QUADRATIC_PLANAR::_precision=eps;
  INTERP_KERNEL::QUADRATIC_PLANAR::_arc_detection_precision=eps;
  INTERP_KERNEL::QuadraticPolygon *pol=0;
  if(isQuad)
    pol=INTERP_KERNEL::QuadraticPolygon::BuildArcCirclePolygon(nodes);
  else
    pol=INTERP_KERNEL::QuadraticPolygon::BuildLinearPolygon(nodes);
  bool ret=pol->isButterflyAbs();
  delete pol;
  return ret;
}
