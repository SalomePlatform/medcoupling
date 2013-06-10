// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
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

std::size_t MEDCouplingPointSet::getHeapMemorySize() const
{
  std::size_t ret=0;
  if(_coords)
    ret+=_coords->getHeapMemorySize();
  return MEDCouplingMesh::getHeapMemorySize()+ret;
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

/*!
 * Returns a pointer to the array of point coordinates held by \a this.
 *  \return DataArrayDouble * - the pointer to the array of point coordinates. The
 *          caller is to delete this array using decrRef() as it is no more needed.
 */
DataArrayDouble *MEDCouplingPointSet::getCoordinatesAndOwner() const
{
  if(_coords)
    _coords->incrRef();
  return _coords;
}

/*!
 * Copies string attributes from an \a other mesh. The copied strings are
 * - mesh name
 * - mesh description
 * - time units
 * - textual data of the coordinates array (name and components info)
 *
 *  \param [in] other - the mesh to copy string attributes from.
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

/*!
 * Checks equality of point coordinates with coordinates of an \a other mesh.
 *        None textual data is considered.
 *  \param [in] other - the mesh to compare coordinates with \a this one.
 *  \param [in] prec - precision value to compare coordinates.
 *  \return bool - \a true if coordinates of points are equal, \a false else.
 */
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

/*!
 * Checks equality of point coordinates with \a other point coordinates.
 *        Textual data (name and components info) \b is compared as well.
 *  \param [in] other - the point coordinates to compare with \a this one.
 *  \param [in] prec - precision value to compare coordinates.
 *  \return bool - \a true if coordinates of points are equal, \a false else.
 */
bool MEDCouplingPointSet::areCoordsEqual(const MEDCouplingPointSet& other, double prec) const
{
  std::string tmp;
  return areCoordsEqualIfNotWhy(other,prec,tmp);
}

/*!
 * Checks equality of point coordinates with \a other point coordinates.
 *        None textual data is considered.
 *  \param [in] other - the point coordinates to compare with \a this one.
 *  \param [in] prec - precision value to compare coordinates.
 *  \return bool - \a true if coordinates of points are equal, \a false else.
 */
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
 * Returns coordinates of \a nodeId-th node.
 *  \param [in] nodeId - the ID of the node of interest.
 *  \param [in, out] coo - the array filled with coordinates of the \a nodeId-th
 *         node. This array is not cleared before filling in, the coordinates are
 *         appended to its end.
 *  \throw If the coordinates array is not set.
 *  \throw If \a nodeId is not a valid index for the coordinates array.
 *
 *  \ref cpp_mcpointset_getcoordinatesofnode "Here is a C++ example".<br>
 *  \ref  py_mcpointset_getcoordinatesofnode "Here is a Python example".
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
 * Finds nodes equal within \a precision and returns an array describing the 
 * permutation to remove duplicated nodes.
 *  \param [in] precision - minimal absolute distance between two nodes at which they are
 *              considered not coincident.
 *  \param [in] limitNodeId - limit node id. If all nodes within a group of coincident
 *              nodes have id strictly lower than \a limitTupleId then they are not
 *              returned. Put -1 to this parameter to have all nodes returned.
 *  \param [out] areNodesMerged - is set to \a true if any coincident nodes found.
 *  \param [out] newNbOfNodes - returns number of unique nodes.
 *  \return DataArrayInt * - the permutation array in "Old to New" mode. For more 
 *          info on "Old to New" mode see \ref MEDCouplingArrayRenumbering. The caller
 *          is to delete this array using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 */
DataArrayInt *MEDCouplingPointSet::buildPermArrayForMergeNode(double precision, int limitNodeId, bool& areNodesMerged, int& newNbOfNodes) const
{
  DataArrayInt *comm,*commI;
  findCommonNodes(precision,limitNodeId,comm,commI);
  int oldNbOfNodes=getNumberOfNodes();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=buildNewNumberingFromCommonNodesFormat(comm,commI,newNbOfNodes);
  areNodesMerged=(oldNbOfNodes!=newNbOfNodes);
  comm->decrRef();
  commI->decrRef();
  return ret.retn();
}

/*!
 * Finds nodes coincident within \a prec tolerance.
 * Ids of coincident nodes are stored in output arrays.
 * A pair of arrays (\a comm, \a commIndex) is called "Surjective Format 2".
 *  \param [in] prec - minimal absolute distance (using infinite norm) between two nodes at which they are
 *              considered not coincident.
 *  \param [in] limitNodeId - limit node id. If all nodes within a group of coincident
 *              nodes have id strictly lower than \a limitTupleId then they are not
 *              returned. Put -1 to this parameter to have all nodes treated.
 *  \param [out] comm - the array holding ids of coincident nodes.
 *               \a comm->getNumberOfComponents() == 1. 
 *               \a comm->getNumberOfTuples() == \a commIndex->back(). The caller
 *               is to delete this array using decrRef() as it is no more needed.
 *  \param [out] commIndex - the array dividing all ids stored in \a comm into
 *               groups of (ids of) coincident nodes. Its every value is a tuple
 *               index where a next group of nodes begins. For example the second
 *               group of nodes in \a comm is described by following range of indices:
 *               [ \a commIndex[1], \a commIndex[2] ). \a commIndex->getNumberOfTuples()-1
 *               gives the number of groups of coincident nodes. The caller
 *               is to delete this array using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *
 *  \ref cpp_mcpointset_findcommonnodes "Here is a C++ example".<br>
 *  \ref  py_mcpointset_findcommonnodes "Here is a Python example".
 */
void MEDCouplingPointSet::findCommonNodes(double prec, int limitNodeId, DataArrayInt *&comm, DataArrayInt *&commIndex) const
{
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::findCommonNodes : no coords specified !");
  _coords->findCommonTuples(prec,limitNodeId,comm,commIndex);
}

/*!
 * Finds nodes located at distances lower that \a eps from a given point.
 *  \param [in] pos - pointer to coordinates of the point.  This array is expected to
 *         be of length \a this->getSpaceDimension() at least, else the
 *         behavior is not warranted.
 *  \param [in] eps - the lowest distance between a point and a node (using infinite norm) at which the node is
 *          not returned by this method.
 *  \return DataArrayInt * - a new instance of DataArrayInt holding ids of nodes
 *          close to the point. The caller is to delete this
 *          array using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *
 *  \ref cpp_mcpointset_getnodeidsnearpoint "Here is a C++ example".<br>
 *  \ref  py_mcpointset_getnodeidsnearpoint "Here is a Python example".
 */
DataArrayInt *MEDCouplingPointSet::getNodeIdsNearPoint(const double *pos, double eps) const throw(INTERP_KERNEL::Exception)
{
  DataArrayInt *c=0,*cI=0;
  getNodeIdsNearPoints(pos,1,eps,c,cI);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cITmp(cI);
  return c;
}

/*!
 * Finds nodes located at distances lower that \a eps from given points.
 *  \param [in] pos - pointer to coordinates of the points. This array is expected to
 *         be of length \a nbOfPoints * \a this->getSpaceDimension() at least, else the
 *         behavior is not warranted.
 *  \param [in] nbOfPoints - number of points whose coordinates are given by \a pos
 *         parameter. 
 *  \param [in] eps - the lowest distance between (using infinite norm) a point and a node at which the node is
 *         not returned by this method.
 *  \param [out] c - array returning ids of nodes located closer than \a eps to the
 *         given points. The caller
 *         is to delete this array using decrRef() as it is no more needed.
 *  \param [out] cI - for each i-th given point, the array specifies tuples of \a c
 *         holding ids of nodes close to the i-th point. <br>The i-th value of \a cI is an 
 *         index of tuple of \a c holding id of a first (if any) node close to the
 *         i-th given point. Difference between the i-th and (i+1)-th value of \a cI
 *         (i.e. \a cI[ i+1 ] - \a cI[ i ]) defines number of nodes close to the i-th
 *         point (that can be zero!). For example, the group of nodes close to the
 *         second point is described by following range of indices [ \a cI[1], \a cI[2] ).
 *         The caller is to delete this array using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *
 *  \ref cpp_mcpointset_getnodeidsnearpoints "Here is a C++ example".<br>
 *  \ref  py_mcpointset_getnodeidsnearpoints "Here is a Python example".
 */
void MEDCouplingPointSet::getNodeIdsNearPoints(const double *pos, int nbOfPoints, double eps, DataArrayInt *& c, DataArrayInt *& cI) const throw(INTERP_KERNEL::Exception)
{
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::getNodeIdsNearPoint : no coordiantes set !");
  int spaceDim=getSpaceDimension();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> points=DataArrayDouble::New();
  points->useArray(pos,false,CPP_DEALLOC,nbOfPoints,spaceDim);
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
  return DataArrayInt::BuildOld2NewArrayFromSurjectiveFormat2(getNumberOfNodes(),comm->begin(),commIndex->begin(),commIndex->end(),newNbOfNodes);
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
 * Computes the minimum box bounding all nodes. The edges of the box are parallel to
 * the Cartesian coordinate axes. The bounding box is described by coordinates of its
 * two extremum points with minimal and maximal coordinates.
 *  \param [out] bbox - array filled with coordinates of extremum points in "no
 *         interlace" mode, i.e. xMin, xMax, yMin, yMax, zMin, zMax (if in 3D). This
 *         array, of length 2 * \a this->getSpaceDimension() at least, is to be
 *         pre-allocated by the caller.
 *  \throw If the coordinates array is not set.
 *
 *  \ref cpp_mcpointset_getBoundingBox "Here is a C++ example".<br>
 *  \ref  py_mcpointset_getBoundingBox "Here is a Python example".
 */
void MEDCouplingPointSet::getBoundingBox(double *bbox) const throw(INTERP_KERNEL::Exception)
{
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::getBoundingBox : Coordinates not set !");
  _coords->getMinMaxPerComponent(bbox);
}

/*!
 * Removes "free" nodes, i.e. nodes not used to define any element.
 *  \throw If the coordinates array is not set.
 *  \throw If the elements are not defined.
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
 * Returns the carateristic dimension of \a this point set, that is a maximal
 * absolute values of node coordinates.
 *  \throw If the coordinates array is not set.
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
 * \warning this method is non const and alterates coordinates in \b this without modifying.
 * \param [in] eps absolute epsilon. under that value of delta between max and min no scale is performed.
 *
 */
void MEDCouplingPointSet::recenterForMaxPrecision(double eps) throw(INTERP_KERNEL::Exception)
{
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::recenterForMaxPrecision : Coordinates not set !");
  _coords->recenterForMaxPrecision(eps);
  updateTime();
}

/*!
 * Rotates \a this set of nodes by \a angle around either an axis (in 3D) or a point
 * (in 2D). 
 *  \param [in] center - coordinates either of an origin of rotation axis (in 3D) or
 *         of center of rotation (in 2D). This array is to be of size \a
 *         this->getSpaceDimension() at least.
 *  \param [in] vector - 3 components of a vector defining direction of the rotation
 *         axis in 3D. In 2D this parameter is not used.
 *  \param [in] angle - the rotation angle in radians.
 *  \throw If the coordinates array is not set.
 *  \throw If \a this->getSpaceDimension() != 2 && \a this->getSpaceDimension() != 3.
 *  \throw If \a center == NULL
 *  \throw If \a vector == NULL && \a this->getSpaceDimension() == 3.
 *  \throw If Magnitude of \a vector is zero.
 *
 *  \ref cpp_mcpointset_rotate "Here is a C++ example".<br>
 *  \ref  py_mcpointset_rotate "Here is a Python example".
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
 * Translates \a this set of nodes. 
 *  \param [in] vector - components of a translation vector. This array is to be of
 *         size \a this->getSpaceDimension() at least. 
 *  \throw If the coordinates array is not set.
 *  \throw If \a vector == NULL.
 *
 *  \ref cpp_mcpointset_translate "Here is a C++ example".<br>
 *  \ref  py_mcpointset_translate "Here is a Python example".
 */
void MEDCouplingPointSet::translate(const double *vector)
{
  if(!vector)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::translate : NULL input vector !");
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::translate : no coordinates set !");
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
 * Applies scaling transformation to \a this set of nodes. 
 *  \param [in] point - coordinates of a scaling center. This array is to be of
 *         size \a this->getSpaceDimension() at least.
 *  \param [in] factor - a scale factor.
 *  \throw If the coordinates array is not set.
 *  \throw If \a point == NULL.
 *
 *  \ref cpp_mcpointset_scale "Here is a C++ example".<br>
 *  \ref  py_mcpointset_scale "Here is a Python example".
 */
void MEDCouplingPointSet::scale(const double *point, double factor)
{
  if(!point)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::scale : NULL input point !");
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::scale : no coordinates set !");
  double *coords=_coords->getPointer();
  int nbNodes=getNumberOfNodes();
  int dim=getSpaceDimension();
  for(int i=0;i<nbNodes;i++)
    {
      std::transform(coords+i*dim,coords+(i+1)*dim,point,coords+i*dim,std::minus<double>());
      std::transform(coords+i*dim,coords+(i+1)*dim,coords+i*dim,std::bind2nd(std::multiplies<double>(),factor));
      std::transform(coords+i*dim,coords+(i+1)*dim,point,coords+i*dim,std::plus<double>());
    }
  _coords->declareAsNew();
  updateTime();
}

/*!
 * Converts \a this set of points to an other dimension by changing number of
 * components of point coordinates. If the dimension increases, added components
 * are filled with \a dftValue. If the dimension decreases, last components are lost.
 * If the new dimension is same as \a this->getSpaceDimension(), nothing is done.
 *  \param [in] newSpaceDim - the new space dimension.
 *  \param [in] dftValue - the value to assign to added components of point coordinates
 *         (if the dimension increases).
 *  \throw If the coordinates array is not set.
 *  \throw If \a newSpaceDim < 1.
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
 * Substitutes \a this->_coords with \a other._coords provided that coordinates of
 * the two point sets match with a specified precision, else an exception is thrown.
 *  \param [in] other - the other point set whose coordinates array will be used by
 *         \a this point set in case of their equality.
 *  \param [in] epsilon - the precision used to compare coordinates.
 *  \throw If the coordinates array of \a this is not set.
 *  \throw If the coordinates array of \a other is not set.
 *  \throw If the coordinates of \a this and \a other do not match.
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
 * Finds nodes located at distance lower that \a eps from a specified plane.
 *  \param [in] pt - 3 components of a point defining location of the plane.
 *  \param [in] vec - 3 components of a normal vector to the plane. Vector magnitude
 *         must be greater than 10*\a eps.
 *  \param [in] eps - maximal distance of a node from the plane at which the node is
 *         considered to lie on the plane.
 *  \param [in,out] nodes - a vector returning ids of found nodes. This vector is not
 *         cleared before filling in.
 *  \throw If the coordinates array is not set.
 *  \throw If \a pt == NULL.
 *  \throw If \a vec == NULL.
 *  \throw If the magnitude of \a vec is zero.
 *  \throw If \a this->getSpaceDimension() != 3.
 */
void MEDCouplingPointSet::findNodesOnPlane(const double *pt, const double *vec, double eps, std::vector<int>& nodes) const throw(INTERP_KERNEL::Exception)
{
  if(getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::findNodesOnPlane : Invalid spacedim to be applied on this ! Must be equal to 3 !");
  if(!pt)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::findNodesOnPlane : NULL point pointer specified !");
  if(!vec)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::findNodesOnPlane : NULL vector pointer specified !");
  int nbOfNodes=getNumberOfNodes();
  double a=vec[0],b=vec[1],c=vec[2],d=-pt[0]*vec[0]-pt[1]*vec[1]-pt[2]*vec[2];
  double deno=sqrt(a*a+b*b+c*c);
  if(deno<std::numeric_limits<double>::min())
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::findNodesOnPlane : vector pointer specified has norm equal to 0. !");
  const double *work=_coords->getConstPointer();
  for(int i=0;i<nbOfNodes;i++)
    {
      if(std::abs(a*work[0]+b*work[1]+c*work[2]+d)/deno<eps)
        nodes.push_back(i);
      work+=3;
    }
}

/*!
 * Finds nodes located at distance lower that \a eps from a specified line in 2D and 3D.
 *  \param [in] pt - components of coordinates of an initial point of the line. This
 *         array is to be of size \a this->getSpaceDimension() at least.
 *  \param [in] vec - components of a vector defining the line direction. This array
 *         is to be of size \a this->getSpaceDimension() at least. Vector magnitude 
 *         must be greater than 10*\a eps.
 *  \param [in] eps - maximal distance of a node from the line at which the node is
 *         considered to lie on the line.
 *  \param [in,out] nodes - a vector returning ids of found nodes. This vector is not
 *         cleared before filling in.
 *  \throw If the coordinates array is not set.
 *  \throw If \a pt == NULL.
 *  \throw If \a vec == NULL.
 *  \throw If the magnitude of \a vec is zero.
 *  \throw If ( \a this->getSpaceDimension() != 3 && \a this->getSpaceDimension() != 2 ).
 */
void MEDCouplingPointSet::findNodesOnLine(const double *pt, const double *vec, double eps, std::vector<int>& nodes) const throw(INTERP_KERNEL::Exception)
{
  int spaceDim=getSpaceDimension();
  if(spaceDim!=2 && spaceDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::findNodesOnLine : Invalid spacedim to be applied on this ! Must be equal to 2 or 3 !");
  if(!pt)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::findNodesOnLine : NULL point pointer specified !");
  if(!vec)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::findNodesOnLine : NULL vector pointer specified !");
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
 * Returns a new array of node coordinates by concatenating node coordinates of two
 * given point sets, so that (1) the number of nodes in the result array is a sum of the
 * number of nodes of given point sets and (2) the number of component in the result array
 * is same as that of each of given point sets. Info on components is copied from the first
 * of the given point set. Space dimension of the given point sets must be the same. 
 *  \param [in] m1 - a point set whose coordinates will be included in the result array.
 *  \param [in] m2 - another point set whose coordinates will be included in the
 *         result array. 
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If both \a m1 and \a m2 are NULL.
 *  \throw If \a m1->getSpaceDimension() != \a m2->getSpaceDimension().
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

void MEDCouplingPointSet::checkCoherency() const throw(INTERP_KERNEL::Exception)
{
  if(!_coords)
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::checkCoherency : no coordinates set !");
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
  if(norm<std::numeric_limits<double>::min())
    throw INTERP_KERNEL::Exception("MEDCouplingPointSet::Rotate3DAlg : magnitude of input vector is too close of 0. !");
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
 * Creates a new MEDCouplingMesh containing a part of cells of \a this mesh. The new
 * mesh shares a coordinates array with \a this one. The cells to include to the
 * result mesh are specified by an array of cell ids.
 *  \param [in] start - an array of cell ids to include to the result mesh.
 *  \param [in] end - specifies the end of the array \a start, so that
 *              the last value of \a start is \a end[ -1 ].
 *  \return MEDCouplingMesh * - a new instance of MEDCouplingMesh. The caller is to
 *         delete this mesh using decrRef() as it is no more needed. 
 */
MEDCouplingMesh *MEDCouplingPointSet::buildPart(const int *start, const int *end) const
{
  return buildPartOfMySelf(start,end,true);
}

/*!
 * Creates a new MEDCouplingMesh containing a part of cells of \a this mesh. The
 * cells to include to the result mesh are specified by an array of cell ids. 
 * <br> This method additionally returns a renumbering map in "Old to New" mode
 * which allows the caller to know the mapping between nodes in \a this and the result mesh.
 *  \param [in] start - an array of cell ids to include to the result mesh.
 *  \param [in] end - specifies the end of the array \a start, so that
 *              the last value of \a start is \a end[ -1 ].
 *  \param [out] arr - a new DataArrayInt that is the "Old to New" renumbering
 *         map. The caller is to delete this array using decrRef() as it is no more needed.
 *  \return MEDCouplingMesh * - a new instance of MEDCouplingMesh. The caller is to
 *         delete this mesh using decrRef() as it is no more needed. 
 */
MEDCouplingMesh *MEDCouplingPointSet::buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingPointSet> ret=buildPartOfMySelf(start,end,true);
  arr=ret->zipCoordsTraducer();
  return ret.retn();
}

/*!
 * This method specialized the MEDCouplingMesh::buildPartRange
 *
 * \sa MEDCouplingUMesh::buildPartOfMySelf2
 */
MEDCouplingMesh *MEDCouplingPointSet::buildPartRange(int beginCellIds, int endCellIds, int stepCellIds) const throw(INTERP_KERNEL::Exception)
{
  return buildPartOfMySelf2(beginCellIds,endCellIds,stepCellIds,true);
}

/*!
 * This method specialized the MEDCouplingMesh::buildPartRangeAndReduceNodes
 *
 * \param [out] beginOut valid only if \a arr not NULL !
 * \param [out] endOut valid only if \a arr not NULL !
 * \param [out] stepOut valid only if \a arr not NULL !
 * \param [out] arr correspondance old to new in node ids.
 * 
 * \sa MEDCouplingUMesh::buildPartOfMySelf2
 */
MEDCouplingMesh *MEDCouplingPointSet::buildPartRangeAndReduceNodes(int beginCellIds, int endCellIds, int stepCellIds, int& beginOut, int& endOut, int& stepOut, DataArrayInt*& arr) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingPointSet> ret=buildPartOfMySelf2(beginCellIds,endCellIds,stepCellIds,true);
  arr=ret->zipCoordsTraducer();
  return ret.retn();
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
