// Copyright (C) 2007-2024  CEA, EDF
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
// Author : Anthony Geay (EDF R&D)

#include "MEDCouplingUMesh.txx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingSkyLineArray.hxx"
#include "CellModel.hxx"
#include "VolSurfUser.txx"
#include "InterpolationUtils.hxx"
#include "PointLocatorAlgos.txx"
#include "BBTree.txx"
#include "BBTreeDst.txx"
#include "SplitterTetra.hxx"
#include "DiameterCalculator.hxx"
#include "DirectedBoundingBox.hxx"
#include "InterpKernelMatrixTools.hxx"
#include "InterpKernelMeshQuality.hxx"
#include "InterpKernelCellSimplify.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelAutoPtr.hxx"
#include "InterpKernelGeo2DNode.hxx"
#include "InterpKernelGeo2DEdgeLin.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelGeo2DQuadraticPolygon.hxx"
#include "OrientationInverter.hxx"
#include "MEDCouplingUMesh_internal.hxx"

#include <sstream>
#include <fstream>
#include <numeric>
#include <memory>
#include <cstring>
#include <limits>
#include <list>

using namespace MEDCoupling;

double MEDCouplingUMesh::EPS_FOR_POLYH_ORIENTATION=1.e-14;

/// @cond INTERNAL

const INTERP_KERNEL::NormalizedCellType MEDCouplingUMesh::MEDMEM_ORDER[N_MEDMEM_ORDER] = { INTERP_KERNEL::NORM_POINT1, INTERP_KERNEL::NORM_SEG2, INTERP_KERNEL::NORM_SEG3, INTERP_KERNEL::NORM_SEG4, INTERP_KERNEL::NORM_POLYL, INTERP_KERNEL::NORM_TRI3, INTERP_KERNEL::NORM_QUAD4, INTERP_KERNEL::NORM_TRI6, INTERP_KERNEL::NORM_TRI7, INTERP_KERNEL::NORM_QUAD8, INTERP_KERNEL::NORM_QUAD9, INTERP_KERNEL::NORM_POLYGON, INTERP_KERNEL::NORM_QPOLYG, INTERP_KERNEL::NORM_TETRA4, INTERP_KERNEL::NORM_PYRA5, INTERP_KERNEL::NORM_PENTA6, INTERP_KERNEL::NORM_HEXA8, INTERP_KERNEL::NORM_HEXGP12, INTERP_KERNEL::NORM_TETRA10, INTERP_KERNEL::NORM_PYRA13, INTERP_KERNEL::NORM_PENTA15, INTERP_KERNEL::NORM_PENTA18, INTERP_KERNEL::NORM_HEXA20, INTERP_KERNEL::NORM_HEXA27, INTERP_KERNEL::NORM_POLYHED };
/// @endcond

MEDCouplingUMesh *MEDCouplingUMesh::New()
{
  return new MEDCouplingUMesh;
}

MEDCouplingUMesh *MEDCouplingUMesh::New(const std::string& meshName, int meshDim)
{
  MEDCouplingUMesh *ret=new MEDCouplingUMesh;
  ret->setName(meshName);
  ret->setMeshDimension(meshDim);
  return ret;
}

/*!
 * Returns a new MEDCouplingUMesh which is a full copy of \a this one. No data is shared
 * between \a this and the new mesh.
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingMesh. The caller is to
 *          delete this mesh using decrRef() as it is no more needed.
 */
MEDCouplingUMesh *MEDCouplingUMesh::deepCopy() const
{
  return clone(true);
}


/*!
 * Returns a new MEDCouplingUMesh which is a copy of \a this one.
 *  \param [in] recDeepCpy - if \a true, the copy is deep, else all data arrays of \a
 * this mesh are shared by the new mesh.
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingMesh. The caller is to
 *          delete this mesh using decrRef() as it is no more needed.
 */
MEDCouplingUMesh *MEDCouplingUMesh::clone(bool recDeepCpy) const
{
  return new MEDCouplingUMesh(*this,recDeepCpy);
}

/*!
 * This method behaves mostly like MEDCouplingUMesh::deepCopy method, except that only nodal connectivity arrays are deeply copied.
 * The coordinates are shared between \a this and the returned instance.
 *
 * \return MEDCouplingUMesh * - A new object instance holding the copy of \a this (deep for connectivity, shallow for coordiantes)
 * \sa MEDCouplingUMesh::deepCopy
 */
MEDCouplingUMesh *MEDCouplingUMesh::deepCopyConnectivityOnly() const
{
  checkConnectivityFullyDefined();
  MCAuto<MEDCouplingUMesh> ret=clone(false);
  MCAuto<DataArrayIdType> c(getNodalConnectivity()->deepCopy()),ci(getNodalConnectivityIndex()->deepCopy());
  ret->setConnectivity(c,ci);
  return ret.retn();
}

void MEDCouplingUMesh::shallowCopyConnectivityFrom(const MEDCouplingPointSet *other)
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::shallowCopyConnectivityFrom : input pointer is null !");
  const MEDCouplingUMesh *otherC=dynamic_cast<const MEDCouplingUMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::shallowCopyConnectivityFrom : input pointer is not an MEDCouplingUMesh instance !");
  MEDCouplingUMesh *otherC2=const_cast<MEDCouplingUMesh *>(otherC);//sorry :(
  setConnectivity(otherC2->getNodalConnectivity(),otherC2->getNodalConnectivityIndex(),true);
}

std::size_t MEDCouplingUMesh::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(MEDCouplingPointSet::getHeapMemorySizeWithoutChildren());
  return ret;
}

std::vector<const BigMemoryObject *> MEDCouplingUMesh::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(MEDCouplingPointSet::getDirectChildrenWithNull());
  ret.push_back(_nodal_connec);
  ret.push_back(_nodal_connec_index);
  return ret;
}

void MEDCouplingUMesh::updateTime() const
{
  MEDCouplingPointSet::updateTime();
  if(_nodal_connec)
    {
      updateTimeWith(*_nodal_connec);
    }
  if(_nodal_connec_index)
    {
      updateTimeWith(*_nodal_connec_index);
    }
}

MEDCouplingUMesh::MEDCouplingUMesh():_mesh_dim(-2),_nodal_connec(0),_nodal_connec_index(0)
{
}

/*!
 * Checks if \a this mesh is well defined. If no exception is thrown by this method,
 * then \a this mesh is most probably is writable, exchangeable and available for most
 * of algorithms. When a mesh is constructed from scratch, it is a good habit to call
 * this method to check that all is in order with \a this mesh.
 *  \throw If the mesh dimension is not set.
 *  \throw If the coordinates array is not set (if mesh dimension != -1 ).
 *  \throw If \a this mesh contains elements of dimension different from the mesh dimension.
 *  \throw If the connectivity data array has more than one component.
 *  \throw If the connectivity data array has a named component.
 *  \throw If the connectivity index data array has more than one component.
 *  \throw If the connectivity index data array has a named component.
 */
void MEDCouplingUMesh::checkConsistencyLight() const
{
  if(_mesh_dim<-1)
    throw INTERP_KERNEL::Exception("No mesh dimension specified !");
  if(_mesh_dim!=-1)
    MEDCouplingPointSet::checkConsistencyLight();
  for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iter=_types.begin();iter!=_types.end();iter++)
    {
      if(ToIdType(INTERP_KERNEL::CellModel::GetCellModel(*iter).getDimension())!=_mesh_dim)
        {
          std::ostringstream message;
          message << "Mesh invalid because dimension is " << _mesh_dim << " and there is presence of cell(s) with type " << (*iter);
          throw INTERP_KERNEL::Exception(message.str().c_str());
        }
    }
  if(_nodal_connec)
    {
      if(_nodal_connec->getNumberOfComponents()!=1)
        throw INTERP_KERNEL::Exception("Nodal connectivity array is expected to be with number of components set to one !");
      if(_nodal_connec->getInfoOnComponent(0)!="")
        throw INTERP_KERNEL::Exception("Nodal connectivity array is expected to have no info on its single component !");
    }
  else
    if(_mesh_dim!=-1)
      throw INTERP_KERNEL::Exception("Nodal connectivity array is not defined !");
  if(_nodal_connec_index)
    {
      if(_nodal_connec_index->getNumberOfComponents()!=1)
        throw INTERP_KERNEL::Exception("Nodal connectivity index array is expected to be with number of components set to one !");
      if(_nodal_connec_index->getInfoOnComponent(0)!="")
        throw INTERP_KERNEL::Exception("Nodal connectivity index array is expected to have no info on its single component !");
    }
  else
    if(_mesh_dim!=-1)
      throw INTERP_KERNEL::Exception("Nodal connectivity index array is not defined !");
}

/*!
 * Checks if \a this mesh is well defined. If no exception is thrown by this method,
 * then \a this mesh is informatically clean, most probably is writable, exchangeable and available for all
 * algorithms. <br> In addition to the checks performed by checkConsistencyLight(), this
 * method thoroughly checks the nodal connectivity. For more geometrical checking
 * checkGeomConsistency method is better than this.
 * 
 * \sa MEDCouplingUMesh::checkGeomConsistency
 * 
 *  \param [in] eps - a not used parameter.
 *  \throw If the mesh dimension is not set.
 *  \throw If the coordinates array is not set (if mesh dimension != -1 ).
 *  \throw If \a this mesh contains elements of dimension different from the mesh dimension.
 *  \throw If the connectivity data array has more than one component.
 *  \throw If the connectivity data array has a named component.
 *  \throw If the connectivity index data array has more than one component.
 *  \throw If the connectivity index data array has a named component.
 *  \throw If number of nodes defining an element does not correspond to the type of element.
 *  \throw If the nodal connectivity includes an invalid node id.
 */
void MEDCouplingUMesh::checkConsistency(double eps) const
{
  checkConsistencyLight();
  if(_mesh_dim==-1)
    return ;
  int meshDim=getMeshDimension();
  mcIdType nbOfNodes=getNumberOfNodes();
  mcIdType nbOfCells=getNumberOfCells();
  const mcIdType *ptr=_nodal_connec->getConstPointer();
  const mcIdType *ptrI=_nodal_connec_index->getConstPointer();
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)ptr[ptrI[i]]);
      if(ToIdType(cm.getDimension())!=meshDim)
        {
          std::ostringstream oss;
          oss << "MEDCouplingUMesh::checkConsistency : cell << #" << i<< " with type Type " << cm.getRepr() << " in 'this' whereas meshdim == " << meshDim << " !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      mcIdType nbOfNodesInCell=ptrI[i+1]-ptrI[i]-1;
      if(!cm.isDynamic())
        if(nbOfNodesInCell!=ToIdType(cm.getNumberOfNodes()))
          {
            std::ostringstream oss;
            oss << "MEDCouplingUMesh::checkConsistency : cell #" << i << " with static Type '" << cm.getRepr() << "' has " <<  cm.getNumberOfNodes();
            oss << " nodes whereas in connectivity there is " << nbOfNodesInCell << " nodes ! Looks very bad !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      if(cm.isQuadratic() && cm.isDynamic() && meshDim == 2)
        if (nbOfNodesInCell % 2 || nbOfNodesInCell < 4)
          {
            std::ostringstream oss;
            oss << "MEDCouplingUMesh::checkConsistency : cell #" << i << " with quadratic type '" << cm.getRepr() << "' has " <<  nbOfNodesInCell;
            oss << " nodes. This should be even, and greater or equal than 4!! Looks very bad!";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      for(const mcIdType *w=ptr+ptrI[i]+1;w!=ptr+ptrI[i+1];w++)
        {
          mcIdType nodeId=*w;
          if(nodeId>=0)
            {
              if(nodeId>=nbOfNodes)
                {
                  std::ostringstream oss; oss << "Cell #" << i << " is built with node #" << nodeId << " whereas there are only " << nbOfNodes << " nodes in the mesh !";
                  throw INTERP_KERNEL::Exception(oss.str());
                }
            }
          else if(nodeId<-1)
            {
              std::ostringstream oss; oss << "Cell #" << i << " is built with node #" << nodeId << " in connectivity ! sounds bad !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
          else
            {
              if((INTERP_KERNEL::NormalizedCellType)(ptr[ptrI[i]])!=INTERP_KERNEL::NORM_POLYHED)
                {
                  std::ostringstream oss; oss << "Cell #" << i << " is built with node #-1 in connectivity ! sounds bad !";
                  throw INTERP_KERNEL::Exception(oss.str());
                }
            }
        }
    }
}

/*!
 * This method adds some geometrical checks in addition to the informatical check of checkConsistency method.
 * This method in particular checks that a same node is not repeated several times in a cell.
 * 
 *  \throw If there is a presence a multiple same node ID in nodal connectivity of cell.
 */
void MEDCouplingUMesh::checkGeomConsistency(double eps) const
{
  this->checkConsistency(eps);
  auto nbOfCells(getNumberOfCells());
  const mcIdType *ptr(_nodal_connec->begin()),*ptrI(_nodal_connec_index->begin());
  for(auto icell = 0 ; icell < nbOfCells ; ++icell)
  {
    std::set<mcIdType> s(ptr+ptrI[icell]+1,ptr+ptrI[icell+1]);
    if(ToIdType(s.size())==ptrI[icell+1]-ptrI[icell]-1)
      continue;
    std::ostringstream oss; oss << "MEDCouplingUMesh::checkGeomConsistency : for cell #" << icell << " presence of multiple same nodeID !";
    throw INTERP_KERNEL::Exception(oss.str());
  }
}


/*!
 * Sets dimension of \a this mesh. The mesh dimension in general depends on types of
 * elements contained in the mesh. For more info on the mesh dimension see
 * \ref MEDCouplingUMeshPage.
 *  \param [in] meshDim - a new mesh dimension.
 *  \throw If \a meshDim is invalid. A valid range is <em> -1 <= meshDim <= 3</em>.
 */
void MEDCouplingUMesh::setMeshDimension(int meshDim)
{
  if(meshDim<-1 || meshDim>3)
    throw INTERP_KERNEL::Exception("Invalid meshDim specified ! Must be greater or equal to -1 and lower or equal to 3 !");
  _mesh_dim=meshDim;
  declareAsNew();
}

/*!
 * Allocates memory to store an estimation of the given number of cells.
 * The closer the estimation to the number of cells effectively inserted, the less need the library requires
 * to reallocate memory. If the number of cells to be inserted is not known simply assign 0 to this parameter.
 * If a nodal connectivity previously existed before the call of this method, it will be reset.
 *
 *  \param [in] nbOfCells - estimation of the number of cell \a this mesh will contain.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref medcouplingcppexamplesUmeshStdBuild1 "Here is a C++ example".<br>
 *  \ref medcouplingpyexamplesUmeshStdBuild1 "Here is a Python example".
 *  \endif
 */
void MEDCouplingUMesh::allocateCells(mcIdType nbOfCells)
{
  if(nbOfCells<0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::allocateCells : the input number of cells should be >= 0 !");
  if(_nodal_connec_index)
    {
      _nodal_connec_index->decrRef();
    }
  if(_nodal_connec)
    {
      _nodal_connec->decrRef();
    }
  _nodal_connec_index=DataArrayIdType::New();
  _nodal_connec_index->reserve(nbOfCells+1);
  _nodal_connec_index->pushBackSilent(0);
  _nodal_connec=DataArrayIdType::New();
  _nodal_connec->reserve(2*nbOfCells);
  _types.clear();
  declareAsNew();
}

/*!
 * Appends a cell to the connectivity array. For deeper understanding what is
 * happening see \ref MEDCouplingUMeshNodalConnectivity.
 *  \param [in] type - type of cell to add.
 *  \param [in] size - number of nodes constituting this cell.
 *  \param [in] nodalConnOfCell - the connectivity of the cell to add.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref medcouplingcppexamplesUmeshStdBuild1 "Here is a C++ example".<br>
 *  \ref medcouplingpyexamplesUmeshStdBuild1 "Here is a Python example".
 *  \endif
 */
void MEDCouplingUMesh::insertNextCell(INTERP_KERNEL::NormalizedCellType type, mcIdType size, const mcIdType *nodalConnOfCell)
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  if(_nodal_connec_index==0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::insertNextCell : nodal connectivity not set ! invoke allocateCells before calling insertNextCell !");
  if(ToIdType(cm.getDimension())==_mesh_dim)
    {
      if(!cm.isDynamic())
        if(size!=ToIdType(cm.getNumberOfNodes()))
          {
            std::ostringstream oss; oss << "MEDCouplingUMesh::insertNextCell : Trying to push a " << cm.getRepr() << " cell with a size of " << size;
            oss << " ! Expecting " << cm.getNumberOfNodes() << " !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      mcIdType idx=_nodal_connec_index->back();
      mcIdType val=idx+size+1;
      _nodal_connec_index->pushBackSilent(val);
      _nodal_connec->writeOnPlace(idx,type,nodalConnOfCell,size);
      _types.insert(type);
    }
  else
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::insertNextCell : cell type " << cm.getRepr() << " has a dimension " << cm.getDimension();
      oss << " whereas Mesh Dimension of current UMesh instance is set to " << _mesh_dim << " ! Please invoke \"setMeshDimension\" method before or invoke ";
      oss << "\"MEDCouplingUMesh::New\" static method with 2 parameters name and meshDimension !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
}

/*!
 * Compacts data arrays to release unused memory. This method is to be called after
 * finishing cell insertion using \a this->insertNextCell().
 *
 *  \if ENABLE_EXAMPLES
 *  \ref medcouplingcppexamplesUmeshStdBuild1 "Here is a C++ example".<br>
 *  \ref medcouplingpyexamplesUmeshStdBuild1 "Here is a Python example".
 *  \endif
 */
void MEDCouplingUMesh::finishInsertingCells()
{
  _nodal_connec->pack();
  _nodal_connec_index->pack();
  _nodal_connec->declareAsNew();
  _nodal_connec_index->declareAsNew();
  updateTime();
}

/*!
 * Entry point for iteration over cells of this. Warning the returned cell iterator should be deallocated.
 * Useful for python users.
 */
MEDCouplingUMeshCellIterator *MEDCouplingUMesh::cellIterator()
{
  return new MEDCouplingUMeshCellIterator(this);
}

/*!
 * Entry point for iteration over cells groups geo types per geotypes. Warning the returned cell iterator should be deallocated.
 * If \a this is not so that the cells are grouped by geo types, this method will throw an exception.
 * In this case MEDCouplingUMesh::sortCellsInMEDFileFrmt or MEDCouplingUMesh::rearrange2ConsecutiveCellTypes methods for example can be called before invoking this method.
 * Useful for python users.
 */
MEDCouplingUMeshCellByTypeEntry *MEDCouplingUMesh::cellsByType()
{
  if(!checkConsecutiveCellTypes())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::cellsByType : this mesh is not sorted by type !");
  return new MEDCouplingUMeshCellByTypeEntry(this);
}

/*!
 * Returns a set of all cell types available in \a this mesh.
 * \return std::set<INTERP_KERNEL::NormalizedCellType> - the set of cell types.
 * \warning this method does not throw any exception even if \a this is not defined.
 * \sa MEDCouplingUMesh::getAllGeoTypesSorted
 */
std::set<INTERP_KERNEL::NormalizedCellType> MEDCouplingUMesh::getAllGeoTypes() const
{
  return _types;
}

/*!
 * This method returns the sorted list of geometric types in \a this.
 * Sorted means in the same order than the cells in \a this. A single entry in return vector means the maximal chunk of consecutive cells in \a this
 * having the same geometric type. So a same geometric type can appear more than once if the cells are not sorted per geometric type.
 *
 * \throw if connectivity in \a this is not correctly defined.
 *
 * \sa MEDCouplingMesh::getAllGeoTypes
 */
std::vector<INTERP_KERNEL::NormalizedCellType> MEDCouplingUMesh::getAllGeoTypesSorted() const
{
  std::vector<INTERP_KERNEL::NormalizedCellType> ret;
  checkConnectivityFullyDefined();
  mcIdType nbOfCells=getNumberOfCells();
  if(nbOfCells==0)
    return ret;
  if(getNodalConnectivityArrayLen()<1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getAllGeoTypesSorted : the connectivity in this seems invalid !");
  const mcIdType *c(_nodal_connec->begin()),*ci(_nodal_connec_index->begin());
  ret.push_back((INTERP_KERNEL::NormalizedCellType)c[*ci++]);
  for(mcIdType i=1;i<nbOfCells;i++,ci++)
    if(ret.back()!=((INTERP_KERNEL::NormalizedCellType)c[*ci]))
      ret.push_back((INTERP_KERNEL::NormalizedCellType)c[*ci]);
  return ret;
}

/*!
 * This method is a method that compares \a this and \a other.
 * This method compares \b all attributes, even names and component names.
 */
bool MEDCouplingUMesh::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::isEqualIfNotWhy : input other pointer is null !");
  std::ostringstream oss; oss.precision(15);
  const MEDCouplingUMesh *otherC=dynamic_cast<const MEDCouplingUMesh *>(other);
  if(!otherC)
    {
      reason="mesh given in input is not castable in MEDCouplingUMesh !";
      return false;
    }
  if(!MEDCouplingPointSet::isEqualIfNotWhy(other,prec,reason))
    return false;
  if(_mesh_dim!=otherC->_mesh_dim)
    {
      oss << "umesh dimension mismatch : this mesh dimension=" << _mesh_dim << " other mesh dimension=" <<  otherC->_mesh_dim;
      reason=oss.str();
      return false;
    }
  if(_types!=otherC->_types)
    {
      oss << "umesh geometric type mismatch :\nThis geometric types are :";
      for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iter=_types.begin();iter!=_types.end();iter++)
        { const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(*iter); oss << cm.getRepr() << ", "; }
      oss << "\nOther geometric types are :";
      for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iter=otherC->_types.begin();iter!=otherC->_types.end();iter++)
        { const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(*iter); oss << cm.getRepr() << ", "; }
      reason=oss.str();
      return false;
    }
  if(_nodal_connec!=0 || otherC->_nodal_connec!=0)
    if(_nodal_connec==0 || otherC->_nodal_connec==0)
      {
        reason="Only one UMesh between the two this and other has its nodal connectivity DataArrayInt defined !";
        return false;
      }
  if(_nodal_connec!=otherC->_nodal_connec)
    if(!_nodal_connec->isEqualIfNotWhy(*otherC->_nodal_connec,reason))
      {
        reason.insert(0,"Nodal connectivity DataArrayInt differ : ");
        return false;
      }
  if(_nodal_connec_index!=0 || otherC->_nodal_connec_index!=0)
    if(_nodal_connec_index==0 || otherC->_nodal_connec_index==0)
      {
        reason="Only one UMesh between the two this and other has its nodal connectivity index DataArrayInt defined !";
        return false;
      }
  if(_nodal_connec_index!=otherC->_nodal_connec_index)
    if(!_nodal_connec_index->isEqualIfNotWhy(*otherC->_nodal_connec_index,reason))
      {
        reason.insert(0,"Nodal connectivity index DataArrayInt differ : ");
        return false;
      }
  return true;
}

/*!
 * Checks if data arrays of this mesh (node coordinates, nodal
 * connectivity of cells, etc) of two meshes are same. Textual data like name etc. are
 * not considered.
 *  \param [in] other - the mesh to compare with.
 *  \param [in] prec - precision value used to compare node coordinates.
 *  \return bool - \a true if the two meshes are same.
 */
bool MEDCouplingUMesh::isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const
{
  const MEDCouplingUMesh *otherC=dynamic_cast<const MEDCouplingUMesh *>(other);
  if(!otherC)
    return false;
  if(!MEDCouplingPointSet::isEqualWithoutConsideringStr(other,prec))
    return false;
  if(_mesh_dim!=otherC->_mesh_dim)
    return false;
  if(_types!=otherC->_types)
    return false;
  if(_nodal_connec!=0 || otherC->_nodal_connec!=0)
    if(_nodal_connec==0 || otherC->_nodal_connec==0)
      return false;
  if(_nodal_connec!=otherC->_nodal_connec)
    if(!_nodal_connec->isEqualWithoutConsideringStr(*otherC->_nodal_connec))
      return false;
  if(_nodal_connec_index!=0 || otherC->_nodal_connec_index!=0)
    if(_nodal_connec_index==0 || otherC->_nodal_connec_index==0)
      return false;
  if(_nodal_connec_index!=otherC->_nodal_connec_index)
    if(!_nodal_connec_index->isEqualWithoutConsideringStr(*otherC->_nodal_connec_index))
      return false;
  return true;
}

/*!
 * Checks if \a this and \a other meshes are geometrically equivalent with high
 * probability, else an exception is thrown. The meshes are considered equivalent if
 * (1) meshes contain the same number of nodes and the same number of elements of the
 * same types (2) three cells of the two meshes (first, last and middle) are based
 * on coincident nodes (with a specified precision).
 *  \param [in] other - the mesh to compare with.
 *  \param [in] prec - the precision used to compare nodes of the two meshes.
 *  \throw If the two meshes do not match.
 */
void MEDCouplingUMesh::checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const
{
  MEDCouplingPointSet::checkFastEquivalWith(other,prec);
  const MEDCouplingUMesh *otherC=dynamic_cast<const MEDCouplingUMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::checkFastEquivalWith : Two meshes are not not unstructured !");
}

/*!
 * Returns the reverse nodal connectivity. The reverse nodal connectivity enumerates
 * cells each node belongs to.
 * \warning For speed reasons, this method does not check if node ids in the nodal
 *          connectivity correspond to the size of node coordinates array.
 * \param [in,out] revNodal - an array holding ids of cells sharing each node.
 * \param [in,out] revNodalIndx - an array, of length \a this->getNumberOfNodes() + 1,
 *        dividing cell ids in \a revNodal into groups each referring to one
 *        node. Its every element (except the last one) is an index pointing to the
 *         first id of a group of cells. For example cells sharing the node #1 are
 *        described by following range of indices:
 *        [ \a revNodalIndx[1], \a revNodalIndx[2] ) and the cell ids are
 *        \a revNodal[ \a revNodalIndx[1] ], \a revNodal[ \a revNodalIndx[1] + 1], ...
 *        Number of cells sharing the *i*-th node is
 *        \a revNodalIndx[ *i*+1 ] - \a revNodalIndx[ *i* ].
 * \throw If the coordinates array is not set.
 * \throw If the nodal connectivity of cells is not defined.
 *
 * \if ENABLE_EXAMPLES
 * \ref cpp_mcumesh_getReverseNodalConnectivity "Here is a C++ example".<br>
 * \ref  py_mcumesh_getReverseNodalConnectivity "Here is a Python example".
 * \endif
 */
void MEDCouplingUMesh::getReverseNodalConnectivity(DataArrayIdType *revNodal, DataArrayIdType *revNodalIndx) const
{
  checkFullyDefined();
  mcIdType nbOfNodes(getNumberOfNodes());
  mcIdType *revNodalIndxPtr=(mcIdType *)malloc((nbOfNodes+1)*sizeof(mcIdType));
  revNodalIndx->useArray(revNodalIndxPtr,true,DeallocType::C_DEALLOC,nbOfNodes+1,1);
  std::fill(revNodalIndxPtr,revNodalIndxPtr+nbOfNodes+1,0);
  const mcIdType *conn(_nodal_connec->begin()),*connIndex(_nodal_connec_index->begin());
  mcIdType nbOfCells(getNumberOfCells()),nbOfEltsInRevNodal(0);
  for(mcIdType eltId=0;eltId<nbOfCells;eltId++)
    {
      const mcIdType *strtNdlConnOfCurCell(conn+connIndex[eltId]+1),*endNdlConnOfCurCell(conn+connIndex[eltId+1]);
      for(const mcIdType *iter=strtNdlConnOfCurCell;iter!=endNdlConnOfCurCell;iter++)
        if(*iter>=0)//for polyhedrons
          {
            nbOfEltsInRevNodal++;
            revNodalIndxPtr[(*iter)+1]++;
          }
    }
  std::transform(revNodalIndxPtr+1,revNodalIndxPtr+nbOfNodes+1,revNodalIndxPtr,revNodalIndxPtr+1,std::plus<mcIdType>());
  mcIdType *revNodalPtr=(mcIdType *)malloc(nbOfEltsInRevNodal*sizeof(mcIdType));
  revNodal->useArray(revNodalPtr,true,DeallocType::C_DEALLOC,nbOfEltsInRevNodal,1);
  std::fill(revNodalPtr,revNodalPtr+nbOfEltsInRevNodal,-1);
  for(mcIdType eltId=0;eltId<nbOfCells;eltId++)
    {
      const mcIdType *strtNdlConnOfCurCell=conn+connIndex[eltId]+1;
      const mcIdType *endNdlConnOfCurCell=conn+connIndex[eltId+1];
      for(const mcIdType *iter=strtNdlConnOfCurCell;iter!=endNdlConnOfCurCell;iter++)
        if(*iter>=0)//for polyhedrons
          *std::find_if(revNodalPtr+revNodalIndxPtr[*iter],revNodalPtr+revNodalIndxPtr[*iter+1],std::bind(std::equal_to<mcIdType>(),std::placeholders::_1,-1))=eltId;
    }
}

/*!
 * Creates a new MEDCouplingUMesh containing cells, of dimension one less than \a
 * this->getMeshDimension(), that bound cells of \a this mesh. In addition arrays
 * describing correspondence between cells of \a this and the result meshes are
 * returned. The arrays \a desc and \a descIndx (\ref numbering-indirect) describe the descending connectivity,
 * i.e. enumerate cells of the result mesh bounding each cell of \a this mesh. The
 * arrays \a revDesc and \a revDescIndx (\ref numbering-indirect) describe the reverse descending connectivity,
 * i.e. enumerate cells of  \a this mesh bounded by each cell of the result mesh.
 * \warning For speed reasons, this method does not check if node ids in the nodal
 *          connectivity correspond to the size of node coordinates array.
 * \warning Cells of the result mesh are \b not sorted by geometric type, hence,
 *          to write this mesh to the MED file, its cells must be sorted using
 *          sortCellsInMEDFileFrmt().
 *  \param [in,out] desc - the array containing cell ids of the result mesh bounding
 *         each cell of \a this mesh.
 *  \param [in,out] descIndx - the array, of length \a this->getNumberOfCells() + 1,
 *        dividing cell ids in \a desc into groups each referring to one
 *        cell of \a this mesh. Its every element (except the last one) is an index
 *        pointing to the first id of a group of cells. For example cells of the
 *        result mesh bounding the cell #1 of \a this mesh are described by following
 *        range of indices:
 *        [ \a descIndx[1], \a descIndx[2] ) and the cell ids are
 *        \a desc[ \a descIndx[1] ], \a desc[ \a descIndx[1] + 1], ...
 *        Number of cells of the result mesh sharing the *i*-th cell of \a this mesh is
 *        \a descIndx[ *i*+1 ] - \a descIndx[ *i* ].
 *  \param [in,out] revDesc - the array containing cell ids of \a this mesh bounded
 *         by each cell of the result mesh.
 *  \param [in,out] revDescIndx - the array, of length one more than number of cells
 *        in the result mesh,
 *        dividing cell ids in \a revDesc into groups each referring to one
 *        cell of the result mesh the same way as \a descIndx divides \a desc.
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. The caller is to
 *        delete this mesh using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is node defined.
 *  \throw If \a desc == NULL || \a descIndx == NULL || \a revDesc == NULL || \a
 *         revDescIndx == NULL.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_buildDescendingConnectivity "Here is a C++ example".<br>
 *  \ref  py_mcumesh_buildDescendingConnectivity "Here is a Python example".
 *  \endif
 * \sa buildDescendingConnectivity2()
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildDescendingConnectivity(DataArrayIdType *desc, DataArrayIdType *descIndx, DataArrayIdType *revDesc, DataArrayIdType *revDescIndx) const
{
  return buildDescendingConnectivityGen<MinusOneSonsGenerator>(desc,descIndx,revDesc,revDescIndx,MEDCouplingFastNbrer);
}

/*!
 * This method is a generalization of buildDescendingConnectivity method. This method return mesh containing subcells of this at level specified by \a targetDeltaLevel
 * and return descending and reverse descending correspondances to this.
 *
 * \param [in] targetDeltaLevel target delta level compared to \a this mesh dimension. This parameter is expected to be lower than zero.
 * 
 * \throw If targetDeltaLevel is greater or equal to zero
 * \throw If targetDeltaLevel is lower than -meshDim
 * \sa MEDCouplingUMesh::buildDescendingConnectivity, MEDCouplingUMesh::explode3DMeshTo1D
 */
MCAuto<MEDCouplingUMesh> MEDCouplingUMesh::explodeMeshTo(int targetDeltaLevel, MCAuto<DataArrayIdType>& desc, MCAuto<DataArrayIdType>& descIndx, MCAuto<DataArrayIdType>& revDesc, MCAuto<DataArrayIdType>& revDescIndx) const
{
  this->checkConsistencyLight();
  if( targetDeltaLevel >= 0 )
    THROW_IK_EXCEPTION("Input parameter targetDeltaLevel is expected to be lower than zero !");
  if( targetDeltaLevel == -1 )
  {
    desc = DataArrayIdType::New(); descIndx = DataArrayIdType::New(); revDesc = DataArrayIdType::New(); revDescIndx = DataArrayIdType::New();
    MCAuto<MEDCouplingUMesh> ret( this->buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx) );
    return ret;
  }
  if( targetDeltaLevel == -2 && this->getMeshDimension() == 3 )
  {
    desc = DataArrayIdType::New(); descIndx = DataArrayIdType::New(); revDesc = DataArrayIdType::New(); revDescIndx = DataArrayIdType::New();
    MCAuto<MEDCouplingUMesh> ret( this->explode3DMeshTo1D(desc,descIndx,revDesc,revDescIndx) );
    return ret;
  }
  if( targetDeltaLevel == -this->getMeshDimension() )
  {
    MCAuto<MEDCouplingUMesh> ret = MEDCouplingUMesh::Build0DMeshFromCoords( const_cast<DataArrayDouble *>( this->getCoords() ) );
    MEDCouplingUMesh::DeleteCellTypeInIndexedArray(getNodalConnectivity(),getNodalConnectivityIndex(),desc,descIndx);
    revDesc = DataArrayIdType::New(); revDescIndx = DataArrayIdType::New();
    this->getReverseNodalConnectivity(revDesc,revDescIndx);
    return ret;
  }
  THROW_IK_EXCEPTION("Not valid input targetDeltaLevel regarding mesh dimension");
}

/*!
 * \a this has to have a mesh dimension equal to 3. If it is not the case an INTERP_KERNEL::Exception will be thrown.
 * This behaves exactly as MEDCouplingUMesh::buildDescendingConnectivity does except that this method compute directly the transition from mesh dimension 3 to sub edges (dimension 1)
 * in one shot. That is to say that this method is equivalent to 2 successive calls to MEDCouplingUMesh::buildDescendingConnectivity.
 * This method returns 4 arrays and a mesh as MEDCouplingUMesh::buildDescendingConnectivity does.
 * \sa MEDCouplingUMesh::buildDescendingConnectivity
 */
MEDCouplingUMesh *MEDCouplingUMesh::explode3DMeshTo1D(DataArrayIdType *desc, DataArrayIdType *descIndx, DataArrayIdType *revDesc, DataArrayIdType *revDescIndx) const
{
  checkFullyDefined();
  if(getMeshDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::explode3DMeshTo1D : This has to have a mesh dimension to 3 !");
  return buildDescendingConnectivityGen<MinusTwoSonsGenerator>(desc,descIndx,revDesc,revDescIndx,MEDCouplingFastNbrer);
}

/*!
 * This method computes the micro edges constituting each cell in \a this. Micro edge is an edge for non quadratic cells. Micro edge is an half edge for quadratic cells.
 * This method works for both meshes with mesh dimension equal to 2 or 3. Dynamical cells are not supported (polygons, polyhedrons...)
 * 
 * \sa explode3DMeshTo1D, buildDescendingConnectiviy
 */
MEDCouplingUMesh *MEDCouplingUMesh::explodeMeshIntoMicroEdges(DataArrayIdType *desc, DataArrayIdType *descIndx, DataArrayIdType *revDesc, DataArrayIdType *revDescIndx) const
{
   checkFullyDefined();
   switch(getMeshDimension())
     {
     case 2:
       return buildDescendingConnectivityGen<MicroEdgesGenerator2D>(desc,descIndx,revDesc,revDescIndx,MEDCouplingFastNbrer);
     case 3:
       return buildDescendingConnectivityGen<MicroEdgesGenerator2D>(desc,descIndx,revDesc,revDescIndx,MEDCouplingFastNbrer);
     default:
       throw INTERP_KERNEL::Exception("MEDCouplingUMesh::explodeMeshIntoMicroEdges : Only 2D and 3D supported !");
     }
}

/*!
 * Creates a new MEDCouplingUMesh containing cells, of dimension one less than \a
 * this->getMeshDimension(), that bound cells of \a this mesh. In
 * addition arrays describing correspondence between cells of \a this and the result
 * meshes are returned. The arrays \a desc and \a descIndx (\ref numbering-indirect) describe the descending
 * connectivity, i.e. enumerate cells of the result mesh bounding each cell of \a this
 *  mesh. This method differs from buildDescendingConnectivity() in that apart
 * from cell ids, \a desc returns mutual orientation of cells in \a this and the
 * result meshes. So a positive id means that order of nodes in corresponding cells
 * of two meshes is same, and a negative id means a reverse order of nodes. Since a
 * cell with id #0 can't be negative, the array \a desc returns ids in FORTRAN mode,
 * i.e. cell ids are one-based.
 * Arrays \a revDesc and \a revDescIndx (\ref numbering-indirect) describe the reverse descending connectivity,
 * i.e. enumerate cells of  \a this mesh bounded by each cell of the result mesh.
 * \warning For speed reasons, this method does not check if node ids in the nodal
 *          connectivity correspond to the size of node coordinates array.
 * \warning Cells of the result mesh are \b not sorted by geometric type, hence,
 *          to write this mesh to the MED file, its cells must be sorted using
 *          sortCellsInMEDFileFrmt().
 *  \param [in,out] desc - the array containing cell ids of the result mesh bounding
 *         each cell of \a this mesh.
 *  \param [in,out] descIndx - the array, of length \a this->getNumberOfCells() + 1,
 *        dividing cell ids in \a desc into groups each referring to one
 *        cell of \a this mesh. Its every element (except the last one) is an index
 *        pointing to the first id of a group of cells. For example cells of the
 *        result mesh bounding the cell #1 of \a this mesh are described by following
 *        range of indices:
 *        [ \a descIndx[1], \a descIndx[2] ) and the cell ids are
 *        \a desc[ \a descIndx[1] ], \a desc[ \a descIndx[1] + 1], ...
 *        Number of cells of the result mesh sharing the *i*-th cell of \a this mesh is
 *        \a descIndx[ *i*+1 ] - \a descIndx[ *i* ].
 *  \param [in,out] revDesc - the array containing cell ids of \a this mesh bounded
 *         by each cell of the result mesh.
 *  \param [in,out] revDescIndx - the array, of length one more than number of cells
 *        in the result mesh,
 *        dividing cell ids in \a revDesc into groups each referring to one
 *        cell of the result mesh the same way as \a descIndx divides \a desc.
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. This result mesh
 *        shares the node coordinates array with \a this mesh. The caller is to
 *        delete this mesh using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is node defined.
 *  \throw If \a desc == NULL || \a descIndx == NULL || \a revDesc == NULL || \a
 *         revDescIndx == NULL.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_buildDescendingConnectivity2 "Here is a C++ example".<br>
 *  \ref  py_mcumesh_buildDescendingConnectivity2 "Here is a Python example".
 *  \endif
 * \sa buildDescendingConnectivity()
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildDescendingConnectivity2(DataArrayIdType *desc, DataArrayIdType *descIndx, DataArrayIdType *revDesc, DataArrayIdType *revDescIndx) const
{
  return buildDescendingConnectivityGen<MinusOneSonsGenerator>(desc,descIndx,revDesc,revDescIndx,MEDCouplingOrientationSensitiveNbrer);
}

/*!
 * \b WARNING this method do the assumption that connectivity lies on the coordinates set.
 * For speed reasons no check of this will be done. This method calls
 * MEDCouplingUMesh::buildDescendingConnectivity to compute the result.
 * This method lists for every cell in \b this its neighbor \b cells. To compute the result
 * only connectivities are considered.
 * The neighbor cells of a given cell having id 'cellId' are neighbors[neighborsIndx[cellId]:neighborsIndx[cellId+1]].
 * The format of return is hence \ref numbering-indirect.
 *
 * \param [out] neighbors is an array storing all the neighbors of all cells in \b this. This array is newly
 * allocated and should be dealt by the caller. \b neighborsIndx 2nd output
 * parameter allows to select the right part in this array (\ref numbering-indirect). The number of tuples
 * is equal to the last values in \b neighborsIndx.
 * \param [out] neighborsIndx is an array of size this->getNumberOfCells()+1 newly allocated and should be
 * dealt by the caller. This arrays allow to use the first output parameter \b neighbors (\ref numbering-indirect).
 */
void MEDCouplingUMesh::computeNeighborsOfCells(DataArrayIdType *&neighbors, DataArrayIdType *&neighborsIndx) const
{
  MCAuto<DataArrayIdType> desc=DataArrayIdType::New();
  MCAuto<DataArrayIdType> descIndx=DataArrayIdType::New();
  MCAuto<DataArrayIdType> revDesc=DataArrayIdType::New();
  MCAuto<DataArrayIdType> revDescIndx=DataArrayIdType::New();
  MCAuto<MEDCouplingUMesh> meshDM1=buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  meshDM1=0;
  ComputeNeighborsOfCellsAdv(desc,descIndx,revDesc,revDescIndx,neighbors,neighborsIndx);
}

/**
 * Given a set of identifiers indexed by the node IDs of the mesh (and given in the (\ref numbering-indirect format) ,
 * re-arrange the data to produce a set indexed by cell IDs. The mapping between a node ID and a cell ID is done using the connectivity
 * of the mesh (e.g. a triangular element will receive the information from its three vertices).
 * Doublons are eliminated. If present in the inital dataset, the ID of the cell itself is also remooved.
 *
 * \param [in] nodeNeigh a set of identifiers (mcIdType) stored by node index (\ref numbering-indirect format)
 * \param [in] nodeNeighI a set of identifiers (mcIdType) stored by node index (\ref numbering-indirect format)
 * \param [out] cellNeigh This array is newly allocated and should be dealt by the caller. It contains the initial identifiers
 * provided in the input parameters but stored now by cell index (See 2nd output parameter and \ref numbering-indirect).
 * \param [out] cellNeighI is an array of size this->getNumberOfCells()+1 newly allocated and should be
 * dealt by the caller. This arrays allow to use the first output parameter \b neighbors (\ref numbering-indirect).
 *
 * \raise if the number of tuples in nodeNeighI is not equal to the number of nodes in the mesh.
 */
void MEDCouplingUMesh::computeCellNeighborhoodFromNodesOne(const DataArrayIdType *nodeNeigh, const DataArrayIdType *nodeNeighI,
                                                           MCAuto<DataArrayIdType>& cellNeigh, MCAuto<DataArrayIdType>& cellNeighIndex) const
{
  if(!nodeNeigh || !nodeNeighI)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::computeCellNeighborhoodFromNodesOne : null pointer !");
  checkConsistencyLight();
  nodeNeigh->checkAllocated(); nodeNeighI->checkAllocated();
  nodeNeigh->checkNbOfComps(1,"MEDCouplingUMesh::computeCellNeighborhoodFromNodesOne : node neigh");
  nodeNeighI->checkNbOfComps(1,"MEDCouplingUMesh::computeCellNeighborhoodFromNodesOne : node neigh index");
  nodeNeighI->checkNbOfTuples(1+getNumberOfNodes(),"MEDCouplingUMesh::computeCellNeighborhoodFromNodesOne : invalid length");
  mcIdType nbCells=getNumberOfCells();
  const mcIdType *c(_nodal_connec->begin()),*ci(_nodal_connec_index->begin()),*ne(nodeNeigh->begin()),*nei(nodeNeighI->begin());
  cellNeigh=DataArrayIdType::New(); cellNeigh->alloc(0,1); cellNeighIndex=DataArrayIdType::New(); cellNeighIndex->alloc(1,1); cellNeighIndex->setIJ(0,0,0);
  for(mcIdType i=0;i<nbCells;i++)
    {
      std::set<mcIdType> s;
      for(const mcIdType *it=c+ci[i]+1;it!=c+ci[i+1];it++)
        if(*it>=0)  // avoid -1 in polygons or polyedrons
          s.insert(ne+nei[*it],ne+nei[*it+1]);
      s.erase(i);
      cellNeigh->insertAtTheEnd(s.begin(),s.end());
      cellNeighIndex->pushBackSilent(cellNeigh->getNumberOfTuples());
    }
}

/*!
 * This method is called by MEDCouplingUMesh::computeNeighborsOfCells. This methods performs the algorithm
 * of MEDCouplingUMesh::computeNeighborsOfCells.
 * This method is useful for users that want to reduce along a criterion the set of neighbours cell. This is
 * typically the case to extract a set a neighbours,
 * excluding a set of meshdim-1 cells in input descending connectivity.
 * Typically \b desc, \b descIndx, \b revDesc and \b revDescIndx (\ref numbering-indirect) input params are
 * the result of MEDCouplingUMesh::buildDescendingConnectivity.
 * This method lists cell by cell in \b this which are its neighbors. To compute the result only connectivities
 * are considered.
 * The neighbor cells of cell having id 'cellId' are neighbors[neighborsIndx[cellId]:neighborsIndx[cellId+1]].
 *
 * \param [in] desc descending connectivity array.
 * \param [in] descIndx descending connectivity index array used to walk through \b desc (\ref numbering-indirect).
 * \param [in] revDesc reverse descending connectivity array.
 * \param [in] revDescIndx reverse descending connectivity index array used to walk through \b revDesc (\ref numbering-indirect).
 * \param [out] neighbors is an array storing all the neighbors of all cells in \b this. This array is newly allocated and should be dealt by the caller. \b neighborsIndx 2nd output
 *                        parameter allows to select the right part in this array. The number of tuples is equal to the last values in \b neighborsIndx.
 * \param [out] neighborsIndx is an array of size this->getNumberOfCells()+1 newly allocated and should be dealt by the caller. This arrays allow to use the first output parameter \b neighbors.
 */
void MEDCouplingUMesh::ComputeNeighborsOfCellsAdv(const DataArrayIdType *desc, const DataArrayIdType *descIndx, const DataArrayIdType *revDesc, const DataArrayIdType *revDescIndx,
                                                  DataArrayIdType *&neighbors, DataArrayIdType *&neighborsIndx)
{
  if(!desc || !descIndx || !revDesc || !revDescIndx)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ComputeNeighborsOfCellsAdv some input array is empty !");
  const mcIdType *descPtr=desc->begin();
  const mcIdType *descIPtr=descIndx->begin();
  const mcIdType *revDescPtr=revDesc->begin();
  const mcIdType *revDescIPtr=revDescIndx->begin();
  //
  mcIdType nbCells=descIndx->getNumberOfTuples()-1;
  MCAuto<DataArrayIdType> out0=DataArrayIdType::New();
  MCAuto<DataArrayIdType> out1=DataArrayIdType::New(); out1->alloc(nbCells+1,1);
  mcIdType *out1Ptr=out1->getPointer();
  *out1Ptr++=0;
  out0->reserve(desc->getNumberOfTuples());
  for(mcIdType i=0;i<nbCells;i++,descIPtr++,out1Ptr++)
    {
      for(const mcIdType *w1=descPtr+descIPtr[0];w1!=descPtr+descIPtr[1];w1++)
        {
          std::set<mcIdType> s(revDescPtr+revDescIPtr[*w1],revDescPtr+revDescIPtr[(*w1)+1]);
          s.erase(i);
          out0->insertAtTheEnd(s.begin(),s.end());
        }
      *out1Ptr=out0->getNumberOfTuples();
    }
  neighbors=out0.retn();
  neighborsIndx=out1.retn();
}

/*!
 * Explodes \a this into edges whatever its dimension.
 */
MCAuto<MEDCouplingUMesh> MEDCouplingUMesh::explodeIntoEdges(MCAuto<DataArrayIdType>& desc, MCAuto<DataArrayIdType>& descIndex, MCAuto<DataArrayIdType>& revDesc, MCAuto<DataArrayIdType>& revDescIndx) const
{
  checkFullyDefined();
  int mdim(getMeshDimension());
  desc=DataArrayIdType::New(); descIndex=DataArrayIdType::New(); revDesc=DataArrayIdType::New(); revDescIndx=DataArrayIdType::New();
  MCAuto<MEDCouplingUMesh> mesh1D;
  switch(mdim)
  {
    case 3:
      {
        mesh1D=explode3DMeshTo1D(desc,descIndex,revDesc,revDescIndx);
        break;
      }
    case 2:
      {
        mesh1D=buildDescendingConnectivity(desc,descIndex,revDesc,revDescIndx);
        break;
      }
    default:
      {
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::computeNeighborsOfNodes : Mesh dimension supported are [3,2] !");
      }
  }
  return mesh1D;
}

/*!
 * \b WARNING this method do the assumption that connectivity lies on the coordinates set.
 * For speed reasons no check of this will be done. This method calls
 * MEDCouplingUMesh::buildDescendingConnectivity to compute the result.
 * This method lists for every node in \b this its neighbor \b nodes. To compute the result
 * only connectivities are considered.
 * The neighbor nodes of node having id 'nodeId' are neighbors[neighborsIndx[cellId]:neighborsIndx[cellId+1]].
 *
 * \param [out] neighbors is an array storing all the neighbors of all nodes in \b this. This array
 * is newly allocated and should be dealt by the caller. \b neighborsIndx 2nd output
 * parameter allows to select the right part in this array (\ref numbering-indirect).
 * The number of tuples is equal to the last values in \b neighborsIndx.
 * \param [out] neighborsIdx is an array of size this->getNumberOfCells()+1 newly allocated and should
 * be dealt by the caller. This arrays allow to use the first output parameter \b neighbors.
 *
 * \sa MEDCouplingUMesh::computeEnlargedNeighborsOfNodes
 */
void MEDCouplingUMesh::computeNeighborsOfNodes(DataArrayIdType *&neighbors, DataArrayIdType *&neighborsIdx) const
{
  checkFullyDefined();
  mcIdType mdim(getMeshDimension()),nbNodes(getNumberOfNodes());
  MCAuto<DataArrayIdType> desc(DataArrayIdType::New()),descIndx(DataArrayIdType::New()),revDesc(DataArrayIdType::New()),revDescIndx(DataArrayIdType::New());
  MCConstAuto<MEDCouplingUMesh> mesh1D;
  switch(mdim)
  {
    case 3:
      {
        mesh1D=explode3DMeshTo1D(desc,descIndx,revDesc,revDescIndx);
        break;
      }
    case 2:
      {
        mesh1D=buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
        break;
      }
    case 1:
      {
        mesh1D.takeRef(this);
        break;
      }
    default:
      {
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::computeNeighborsOfNodes : Mesh dimension supported are [3,2,1] !");
      }
  }
  desc=DataArrayIdType::New(); descIndx=DataArrayIdType::New(); revDesc=0; revDescIndx=0;
  mesh1D->getReverseNodalConnectivity(desc,descIndx);
  MCAuto<DataArrayIdType> ret0(DataArrayIdType::New());
  ret0->alloc(desc->getNumberOfTuples(),1);
  mcIdType *r0Pt(ret0->getPointer());
  const mcIdType *c1DPtr(mesh1D->getNodalConnectivity()->begin()),*rn(desc->begin()),*rni(descIndx->begin());
  for(mcIdType i=0;i<nbNodes;i++,rni++)
    {
      for(const mcIdType *oneDCellIt=rn+rni[0];oneDCellIt!=rn+rni[1];oneDCellIt++)
        *r0Pt++=c1DPtr[3*(*oneDCellIt)+1]==i?c1DPtr[3*(*oneDCellIt)+2]:c1DPtr[3*(*oneDCellIt)+1];
    }
  neighbors=ret0.retn();
  neighborsIdx=descIndx.retn();
}

/*!
 * Computes enlarged neighbors for each nodes in \a this. The behavior of this method is close to MEDCouplingUMesh::computeNeighborsOfNodes except that the neighborhood of each node is wider here.
 * A node j is considered to be in the neighborhood of i if and only if there is a cell in \a this containing in its nodal connectivity both i and j.
 * This method is useful to find ghost cells of a part of a mesh with a code based on fields on nodes.
 *
 * \sa MEDCouplingUMesh::computeNeighborsOfNodes
 */
void MEDCouplingUMesh::computeEnlargedNeighborsOfNodes(MCAuto<DataArrayIdType> &neighbors, MCAuto<DataArrayIdType>& neighborsIdx) const
{
  checkFullyDefined();
  mcIdType nbOfNodes(getNumberOfNodes());
  const mcIdType *conn(_nodal_connec->begin()),*connIndex(_nodal_connec_index->begin());
  mcIdType nbOfCells=getNumberOfCells();
  std::vector< std::set<mcIdType> > st0(nbOfNodes);
  for(mcIdType eltId=0;eltId<nbOfCells;eltId++)
    {
      const mcIdType *strtNdlConnOfCurCell(conn+connIndex[eltId]+1),*endNdlConnOfCurCell(conn+connIndex[eltId+1]);
      std::set<mcIdType> s(strtNdlConnOfCurCell,endNdlConnOfCurCell); s.erase(-1); //for polyhedrons
      for(std::set<mcIdType>::const_iterator iter2=s.begin();iter2!=s.end();iter2++)
        st0[*iter2].insert(s.begin(),s.end());
    }
  neighborsIdx=DataArrayIdType::New(); neighborsIdx->alloc(nbOfNodes+1,1); neighborsIdx->setIJ(0,0,0);
  {
    mcIdType *neighIdx(neighborsIdx->getPointer());
    for(std::vector< std::set<mcIdType> >::const_iterator it=st0.begin();it!=st0.end();it++,neighIdx++)
      {
        if ((*it).empty())
          neighIdx[1]=neighIdx[0];
        else
          neighIdx[1]=neighIdx[0]+ToIdType((*it).size())-1;
      }
  }
  neighbors=DataArrayIdType::New(); neighbors->alloc(neighborsIdx->back(),1);
  {
    const mcIdType *neighIdx(neighborsIdx->begin());
    mcIdType *neigh(neighbors->getPointer()),nodeId(0);
    for(std::vector< std::set<mcIdType> >::const_iterator it=st0.begin();it!=st0.end();it++,neighIdx++,nodeId++)
      {
        std::set<mcIdType> s(*it); s.erase(nodeId);
        std::copy(s.begin(),s.end(),neigh+*neighIdx);
      }
  }
}

/*!
 * Converts specified cells to either polygons (if \a this is a 2D mesh) or
 * polyhedrons (if \a this is a 3D mesh). The cells to convert are specified by an
 * array of cell ids. Pay attention that after conversion all algorithms work slower
 * with \a this mesh than before conversion. <br> If an exception is thrown during the
 * conversion due presence of invalid ids in the array of cells to convert, as a
 * result \a this mesh contains some already converted elements. In this case the 2D
 * mesh remains valid but 3D mesh becomes \b inconsistent!
 *  \warning This method can significantly modify the order of geometric types in \a this,
 *          hence, to write this mesh to the MED file, its cells must be sorted using
 *          sortCellsInMEDFileFrmt().
 *  \param [in] cellIdsToConvertBg - the array holding ids of cells to convert.
 *  \param [in] cellIdsToConvertEnd - a pointer to the last-plus-one-th element of \a
 *         cellIdsToConvertBg.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is node defined.
 *  \throw If dimension of \a this mesh is not either 2 or 3.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_convertToPolyTypes "Here is a C++ example".<br>
 *  \ref  py_mcumesh_convertToPolyTypes "Here is a Python example".
 *  \endif
 */
void MEDCouplingUMesh::convertToPolyTypes(const mcIdType *cellIdsToConvertBg, const mcIdType *cellIdsToConvertEnd)
{
  checkFullyDefined();
  int dim=getMeshDimension();
  if(dim<2 || dim>3)
    throw INTERP_KERNEL::Exception("Invalid mesh dimension : must be 2 or 3 !");
  mcIdType nbOfCells=getNumberOfCells();
  if(dim==2)
    {
      const mcIdType *connIndex=_nodal_connec_index->begin();
      mcIdType *conn=_nodal_connec->getPointer();
      for(const mcIdType *iter=cellIdsToConvertBg;iter!=cellIdsToConvertEnd;iter++)
        {
          if(*iter>=0 && *iter<nbOfCells)
            {
              const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[connIndex[*iter]]);
              if(!cm.isQuadratic())
                conn[connIndex[*iter]]=INTERP_KERNEL::NORM_POLYGON;
              else
                conn[connIndex[*iter]]=INTERP_KERNEL::NORM_QPOLYG;
            }
          else
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::convertToPolyTypes : On rank #" << std::distance(cellIdsToConvertBg,iter) << " value is " << *iter << " which is not";
              oss << " in range [0," << nbOfCells << ") !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
        }
    }
  else
    {
      mcIdType *connIndex(_nodal_connec_index->getPointer());
      const mcIdType *connOld(_nodal_connec->getConstPointer());
      MCAuto<DataArrayIdType> connNew(DataArrayIdType::New()),connNewI(DataArrayIdType::New()); connNew->alloc(0,1); connNewI->alloc(1,1); connNewI->setIJ(0,0,0);
      std::vector<bool> toBeDone(nbOfCells,false);
      for(const mcIdType *iter=cellIdsToConvertBg;iter!=cellIdsToConvertEnd;iter++)
        {
          if(*iter>=0 && *iter<nbOfCells)
            toBeDone[*iter]=true;
          else
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::convertToPolyTypes : On rank #" << std::distance(cellIdsToConvertBg,iter) << " value is " << *iter << " which is not";
              oss << " in range [0," << nbOfCells << ") !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
        }
      for(mcIdType cellId=0;cellId<nbOfCells;cellId++)
        {
          mcIdType pos(connIndex[cellId]),posP1(connIndex[cellId+1]);
          mcIdType lgthOld(posP1-pos-1);
          if(toBeDone[cellId])
            {
              const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)connOld[pos]);
              unsigned nbOfFaces(cm.getNumberOfSons2(connOld+pos+1,lgthOld));
              mcIdType *tmp(new mcIdType[nbOfFaces*lgthOld+1]);
              mcIdType *work=tmp; *work++=INTERP_KERNEL::NORM_POLYHED;
              for(unsigned j=0;j<nbOfFaces;j++)
                {
                  INTERP_KERNEL::NormalizedCellType type;
                  unsigned offset=cm.fillSonCellNodalConnectivity2(j,connOld+pos+1,lgthOld,work,type);
                  work+=offset;
                  *work++=-1;
                }
              std::size_t newLgth(std::distance(tmp,work)-1);//-1 for last -1
              connNew->pushBackValsSilent(tmp,tmp+newLgth);
              connNewI->pushBackSilent(connNewI->back()+ToIdType(newLgth));
              delete [] tmp;
            }
          else
            {
              connNew->pushBackValsSilent(connOld+pos,connOld+posP1);
              connNewI->pushBackSilent(connNewI->back()+posP1-pos);
            }
        }
      setConnectivity(connNew,connNewI,false);//false because computeTypes called just behind.
    }
  computeTypes();
}

/*!
 * Converts all cells to either polygons (if \a this is a 2D mesh) or
 * polyhedrons (if \a this is a 3D mesh).
 *  \warning As this method is purely for user-friendliness and no optimization is
 *          done to avoid construction of a useless vector, this method can be costly
 *          in memory.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is node defined.
 *  \throw If dimension of \a this mesh is not either 2 or 3.
 */
void MEDCouplingUMesh::convertAllToPoly()
{
  mcIdType nbOfCells=getNumberOfCells();
  std::vector<mcIdType> cellIds(nbOfCells);
  for(mcIdType i=0;i<nbOfCells;i++)
    cellIds[i]=i;
  convertToPolyTypes(&cellIds[0],&cellIds[0]+ToIdType(cellIds.size()));
}

/*!
 * Fixes nodal connectivity of invalid cells of type NORM_POLYHED. This method
 * expects that all NORM_POLYHED cells have connectivity similar to that of prismatic
 * volumes like NORM_HEXA8, NORM_PENTA6 etc., i.e. the first half of nodes describes a
 * base facet of the volume and the second half of nodes describes an opposite facet
 * having the same number of nodes as the base one. This method converts such
 * connectivity to a valid polyhedral format where connectivity of each facet is
 * explicitly described and connectivity of facets are separated by -1. If \a this mesh
 * contains a NORM_POLYHED cell with a valid connectivity, or an invalid connectivity is
 * not as expected, an exception is thrown and the mesh remains unchanged. Care of
 * a correct orientation of the first facet of a polyhedron, else orientation of a
 * corrected cell is reverse.<br>
 * This method is useful to build an extruded unstructured mesh with polyhedrons as
 * it releases the user from boring description of polyhedra connectivity in the valid
 * format.
 *  \throw If \a this->getMeshDimension() != 3.
 *  \throw If \a this->getSpaceDimension() != 3.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If the coordinates array is not set.
 *  \throw If \a this mesh contains polyhedrons with the valid connectivity.
 *  \throw If \a this mesh contains polyhedrons with odd number of nodes.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_arePolyhedronsNotCorrectlyOriented "Here is a C++ example".<br>
 *  \ref  py_mcumesh_arePolyhedronsNotCorrectlyOriented "Here is a Python example".
 *  \endif
 */
void MEDCouplingUMesh::convertExtrudedPolyhedra()
{
  checkFullyDefined();
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::convertExtrudedPolyhedra works on umeshes with meshdim equal to 3 and spaceDim equal to 3 too!");
  mcIdType nbOfCells=getNumberOfCells();
  MCAuto<DataArrayIdType> newCi=DataArrayIdType::New();
  newCi->alloc(nbOfCells+1,1);
  mcIdType *newci=newCi->getPointer();
  const mcIdType *ci=_nodal_connec_index->getConstPointer();
  const mcIdType *c=_nodal_connec->getConstPointer();
  newci[0]=0;
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)c[ci[i]];
      if(type==INTERP_KERNEL::NORM_POLYHED)
        {
          if(std::count(c+ci[i]+1,c+ci[i+1],-1)!=0)
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::convertExtrudedPolyhedra : cell # " << i << " is a polhedron BUT it has NOT exactly 1 face !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
          std::size_t n2=std::distance(c+ci[i]+1,c+ci[i+1]);
          if(n2%2!=0)
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::convertExtrudedPolyhedra : cell # " << i << " is a polhedron with 1 face but there is a mismatch of number of nodes in face should be even !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
          mcIdType n1=ToIdType(n2/2);
          newci[i+1]=7*n1+2+newci[i];//6*n1 (nodal length) + n1+2 (number of faces) - 1 (number of '-1' separator is equal to number of faces -1) + 1 (for cell type)
        }
      else
        newci[i+1]=(ci[i+1]-ci[i])+newci[i];
    }
  MCAuto<DataArrayIdType> newC=DataArrayIdType::New();
  newC->alloc(newci[nbOfCells],1);
  mcIdType *newc=newC->getPointer();
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)c[ci[i]];
      if(type==INTERP_KERNEL::NORM_POLYHED)
        {
          std::size_t n1=std::distance(c+ci[i]+1,c+ci[i+1])/2;
          newc=std::copy(c+ci[i],c+ci[i]+n1+1,newc);
          *newc++=-1;
          for(std::size_t j=0;j<n1;j++)
            {
              newc[j]=c[ci[i]+1+n1+(n1-j)%n1];
              newc[n1+5*j]=-1;
              newc[n1+5*j+1]=c[ci[i]+1+j];
              newc[n1+5*j+2]=c[ci[i]+1+j+n1];
              newc[n1+5*j+3]=c[ci[i]+1+(j+1)%n1+n1];
              newc[n1+5*j+4]=c[ci[i]+1+(j+1)%n1];
            }
          newc+=n1*6;
        }
      else
        newc=std::copy(c+ci[i],c+ci[i+1],newc);
    }
  _nodal_connec_index->decrRef(); _nodal_connec_index=newCi.retn();
  _nodal_connec->decrRef(); _nodal_connec=newC.retn();
}


/*!
 * Converts all polygons (if \a this is a 2D mesh) or polyhedrons (if \a this is a 3D
 * mesh) to cells of classical types. This method is opposite to convertToPolyTypes().
 * \warning Cells of the result mesh are \b not sorted by geometric type, hence,
 *          to write this mesh to the MED file, its cells must be sorted using
 *          sortCellsInMEDFileFrmt().
 * \warning Cells (and most notably polyhedrons) must be correctly oriented for this to work
 *          properly. See orientCorrectlyPolyhedrons() and arePolyhedronsNotCorrectlyOriented().
 * \return \c true if at least one cell has been converted, \c false else. In the
 *         last case the nodal connectivity remains unchanged.
 * \throw If the coordinates array is not set.
 * \throw If the nodal connectivity of cells is not defined.
 * \throw If \a this->getMeshDimension() < 0.
 */
bool MEDCouplingUMesh::unPolyze()
{
  checkFullyDefined();
  int mdim=getMeshDimension();
  if(mdim<0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::unPolyze works on umeshes with meshdim equals to 0, 1 2 or 3 !");
  if(mdim<=1)
    return false;
  mcIdType nbOfCells=getNumberOfCells();
  if(nbOfCells<1)
    return false;
  mcIdType initMeshLgth=getNodalConnectivityArrayLen();
  mcIdType *conn=_nodal_connec->getPointer();
  mcIdType *index=_nodal_connec_index->getPointer();
  mcIdType posOfCurCell=0;
  mcIdType newPos=0;
  mcIdType lgthOfCurCell;
  bool ret=false;
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      lgthOfCurCell=index[i+1]-posOfCurCell;
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)conn[posOfCurCell];
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      INTERP_KERNEL::NormalizedCellType newType=INTERP_KERNEL::NORM_ERROR;
      mcIdType newLgth=0;
      if(cm.isDynamic())
        {
          switch(cm.getDimension())
          {
            case 2:
              {
                INTERP_KERNEL::AutoPtr<mcIdType> tmp=new mcIdType[lgthOfCurCell-1];
                std::copy(conn+posOfCurCell+1,conn+posOfCurCell+lgthOfCurCell,(mcIdType *)tmp);
                newType=INTERP_KERNEL::CellSimplify::tryToUnPoly2D(cm.isQuadratic(),tmp,lgthOfCurCell-1,conn+newPos+1,newLgth);
                break;
              }
            case 3:
              {
                mcIdType nbOfFaces,lgthOfPolyhConn;
                INTERP_KERNEL::AutoPtr<mcIdType> zipFullReprOfPolyh=INTERP_KERNEL::CellSimplify::getFullPolyh3DCell(type,conn+posOfCurCell+1,lgthOfCurCell-1,nbOfFaces,lgthOfPolyhConn);
                newType=INTERP_KERNEL::CellSimplify::tryToUnPoly3D(zipFullReprOfPolyh,nbOfFaces,lgthOfPolyhConn,conn+newPos+1,newLgth);
                break;
              }
         /*   case 1:  // Not supported yet
              {
                newType=(lgthOfCurCell==3)?INTERP_KERNEL::NORM_SEG2:INTERP_KERNEL::NORM_POLYL;
                break;
              }
         */
          }
          ret=ret || (newType!=type);
          conn[newPos]=newType;
          newPos+=newLgth+1;
          posOfCurCell=index[i+1];
          index[i+1]=newPos;
        }
      else
        {
          std::copy(conn+posOfCurCell,conn+posOfCurCell+lgthOfCurCell,conn+newPos);
          newPos+=lgthOfCurCell;
          posOfCurCell+=lgthOfCurCell;
          index[i+1]=newPos;
        }
    }
  if(newPos!=initMeshLgth)
    _nodal_connec->reAlloc(newPos);
  if(ret)
    computeTypes();
  return ret;
}

/*!
 * This method expects that spaceDimension is equal to 3 and meshDimension equal to 3.
 * This method performs operation only on polyhedrons in \b this. If no polyhedrons exists in \b this, \b this remains unchanged.
 * This method allows to merge if any coplanar 3DSurf cells that may appear in some polyhedrons cells.
 *
 * \b WARNING: this method will not modify edges connectivity! Take a look at colinearizeEdges for that.
 *
 * \param [in] eps is a relative precision that allows to establish if some 3D plane are coplanar or not. This epsilon is used to recenter around origin to have maximal
 *             precision.
 */
void MEDCouplingUMesh::simplifyPolyhedra(double eps)
{
  checkFullyDefined();
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::simplifyPolyhedra : works on meshdimension 3 and spaceDimension 3 !");
  MCAuto<DataArrayDouble> coords=getCoords()->deepCopy();
  coords->recenterForMaxPrecision(eps);
  //
  mcIdType nbOfCells=getNumberOfCells();
  const mcIdType *conn=_nodal_connec->getConstPointer();
  const mcIdType *index=_nodal_connec_index->getConstPointer();
  MCAuto<DataArrayIdType> connINew=DataArrayIdType::New();
  connINew->alloc(nbOfCells+1,1);
  mcIdType *connINewPtr=connINew->getPointer(); *connINewPtr++=0;
  MCAuto<DataArrayIdType> connNew=DataArrayIdType::New(); connNew->alloc(0,1);
  MCAuto<DataArrayIdType> E_Fi(DataArrayIdType::New()), E_F(DataArrayIdType::New()), F_Ei(DataArrayIdType::New()), F_E(DataArrayIdType::New());
  MCAuto<MEDCouplingUMesh> m_faces(buildDescendingConnectivity(E_F, E_Fi, F_E, F_Ei));
  bool changed=false;
  for(mcIdType i=0;i<nbOfCells;i++,connINewPtr++)
    {
      if(conn[index[i]]==ToIdType(INTERP_KERNEL::NORM_POLYHED))
        {
          SimplifyPolyhedronCell(eps,coords, i,connNew, m_faces, E_Fi, E_F, F_Ei, F_E);
          changed=true;
        }
      else
        connNew->insertAtTheEnd(conn+index[i],conn+index[i+1]);
      *connINewPtr=connNew->getNumberOfTuples();
    }
  if(changed)
    setConnectivity(connNew,connINew,false);
}

/*!
 * This method expects that spaceDimension is equal to 3 and meshDimension equal to 3.
 * This method performs operation only on polyhedrons in \b this. If no polyhedrons exists in \b this, \b this remains unchanged.
 * This method allows to simplify edges of polyhedron cells so that consecutive colinear segments (with intermediate points
 * not used by any other cell) are merged together.
 *
 * \param [in] eps is a relative precision that allows to establish if two consecutive 3D segments are colinear or not.
 *
 * \sa simplifyPolyhedra
 */
void MEDCouplingUMesh::colinearizeEdges(double eps)
{
  //
  // Thanks to Antoine Gerschenfeld (CEA) for contributing this method!
  //
  using DAI = MCAuto<DataArrayIdType>;
  checkFullyDefined();
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::colinearizeEdges() : works with meshdim=3 and spaceDim=3!");
  double seps = sqrt(1-eps);
  // Computing connectivities and correspondances : elements -> segments -> points
  DAI E_Fi(DataArrayIdType::New()), E_F(DataArrayIdType::New()), F_Ei(DataArrayIdType::New()), F_E(DataArrayIdType::New()),
         F_Si(DataArrayIdType::New()), F_S(DataArrayIdType::New()), S_Fi(DataArrayIdType::New()), S_F(DataArrayIdType::New()),
         S_Pi(DataArrayIdType::New()), S_P(DataArrayIdType::New()), P_Si(DataArrayIdType::New()), P_S(DataArrayIdType::New());
  MCAuto<MEDCouplingUMesh> m_f(buildDescendingConnectivity(E_F, E_Fi, F_E, F_Ei)),
         m_s(m_f->buildDescendingConnectivity(F_S, F_Si, S_F, S_Fi)),
         m_p(m_s->buildDescendingConnectivity(S_P, S_Pi, P_S, P_Si)); // E: elem, F: faces, S: segments (edges), P: points (vertices)
  const mcIdType *S_Pp(S_P->begin()), *S_Pip(S_Pi->begin()), *P_Sp(P_S->begin()), *P_Sip(P_Si->begin());
  std::set<mcIdType> pt_rem;
  const mcIdType *m_pi = m_p->getNodalConnectivityIndex()->begin(),
                 *m_pc = m_p->getNodalConnectivity()->begin();
  double (*coord)[3] = (double (*)[3]) getCoords()->begin();
  // Find all points only connected to exaclty 2 segments - they are the candidates for elimination
  // Note that in 3D this can only happen for polyhedrons (when this happens at all)
  DAI dsi = P_Si->deltaShiftIndex();
  DAI cand = dsi->findIdsEqual(2);
  for (const mcIdType& i: *cand)  // i is a point to be potentially eliminated, shared by 2 segs only
    {
      double n2[2] = {0., 0.}, scal = 0.; // n2 is a squared norm, scal is a scalar product
      mcIdType p[2][2];                   // p[j][k] is the ID (in the coord array) of the k-th point of the j-th segment
      for (mcIdType j = 0; j < 2; j++)
        for (mcIdType k = 0; k < 2; k++)
          {
            mcIdType off1 = P_Sip[i] + j;   // offset to get ID of the j-th seg (around the i-th point) in the point->seg correspondance
            mcIdType pt_id = P_Sp[off1] + k; // ID of the k-th point of the j-th seg in the point->seg correspondance
            mcIdType pt_id2 = S_Pp[S_Pip[pt_id]]; // ID of the point in the point mesh
            p[j][k] = m_pc[m_pi[pt_id2] + 1];  // Absolute ID, as read from the connectvity (+1 to skip type: NORM_POINT1)
            // Just for fun, as initially written by Antoine :-)
            // p[j][k] = m_pc[m_pi[S_P->getIJ(S_Pi->getIJ(P_S->getIJ(P_Si->getIJ(i, 0) + j, 0), 0) + k, 0)] + 1];
          }
      // Geometric test on scalar product
      for (int d = 0; d < 3; d++) // dimension
        {
          for (int j = 0; j < 2; j++)
            n2[j] += std::pow(coord[p[j][1]][d] - coord[p[j][0]][d], 2);
          scal += (coord[p[1][1]][d] - coord[p[1][0]][d]) * (coord[p[0][1]][d] - coord[p[0][0]][d]);
        }
      if (scal * scal > seps * n2[0] * n2[1]) // seps is a sqrt for homogeneity
        pt_rem.insert(m_pc[m_pi[i] + 1]);  // point should be removed
    }
  // Clean connectivity by filtering points to be removed:
  DataArrayIdType *old_index = getNodalConnectivityIndex(), *old_conn = getNodalConnectivity();
  DAI new_index(DataArrayIdType::New()), new_conn(DataArrayIdType::New());
  const mcIdType *old_index_p(old_index->begin()), *old_conn_p(old_conn->begin());
  for (mcIdType i = 0; i < getNumberOfCells(); i++)
    {
      new_index->pushBackSilent(new_conn->getNbOfElems());
      for (mcIdType j = old_index_p[i]; j < old_index_p[i + 1]; j++)
        {
          // Keep point if it is not to be removed, or if is first in connectivity (TODO this last check could be removed?)
          if (std::find(pt_rem.begin(), pt_rem.end(), old_conn_p[j]) == pt_rem.end() || j == old_index_p[i])
            new_conn->pushBackSilent(old_conn_p[j]);
        }
    }
  new_index->pushBackSilent(new_conn->getNbOfElems());
  setConnectivity(new_conn, new_index);
}

/*!
 * This method returns all node ids used in the connectivity of \b this. The data array returned has to be dealt by the caller.
 * The returned node ids are sorted ascendingly. This method is close to MEDCouplingUMesh::getNodeIdsInUse except
 * the format of the returned DataArrayIdType instance.
 *
 * \return a newly allocated DataArrayIdType sorted ascendingly of fetched node ids.
 * \sa MEDCouplingUMesh::getNodeIdsInUse, areAllNodesFetched
 */
DataArrayIdType *MEDCouplingUMesh::computeFetchedNodeIds() const
{
  checkConnectivityFullyDefined();
  const mcIdType *maxEltPt(std::max_element(_nodal_connec->begin(),_nodal_connec->end()));
  mcIdType maxElt(maxEltPt==_nodal_connec->end()?0:std::abs(*maxEltPt)+1);
  std::vector<bool> retS(maxElt,false);
  computeNodeIdsAlg(retS);
  return DataArrayIdType::BuildListOfSwitchedOn(retS);
}

/*!
 * \param [in,out] nodeIdsInUse an array of size typically equal to nbOfNodes.
 * \sa MEDCouplingUMesh::getNodeIdsInUse, areAllNodesFetched
 */
void MEDCouplingUMesh::computeNodeIdsAlg(std::vector<bool>& nodeIdsInUse) const
{
  mcIdType nbOfNodes=ToIdType(nodeIdsInUse.size()),
           nbOfCells=getNumberOfCells();
  const mcIdType *connIndex(_nodal_connec_index->getConstPointer()),*conn(_nodal_connec->getConstPointer());
  for(mcIdType i=0;i<nbOfCells;i++)
    for(mcIdType j=connIndex[i]+1;j<connIndex[i+1];j++)
      if(conn[j]>=0)
        {
          if(conn[j]<nbOfNodes)
            nodeIdsInUse[conn[j]]=true;
          else
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::computeNodeIdsAlg : In cell #" << i  << " presence of node id " <<  conn[j] << " not in [0," << nbOfNodes << ") !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
        }
}

/// @cond INTERNAL

struct MEDCouplingAccVisit
{
  MEDCouplingAccVisit():_new_nb_of_nodes(0) { }
  mcIdType operator()(mcIdType val) { if(val!=-1) return _new_nb_of_nodes++; else return -1; }
  mcIdType _new_nb_of_nodes;
};

/// @endcond

/*!
 * Finds nodes not used in any cell and returns an array giving a new id to every node
 * by excluding the unused nodes, for which the array holds -1. The result array is
 * a mapping in "Old to New" mode.
 *  \param [out] nbrOfNodesInUse - number of node ids present in the nodal connectivity.
 *  \return DataArrayIdType * - a new instance of DataArrayIdType. Its length is \a
 *          this->getNumberOfNodes(). It holds for each node of \a this mesh either -1
 *          if the node is unused or a new id else. The caller is to delete this
 *          array using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If the nodal connectivity includes an invalid id.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_getNodeIdsInUse "Here is a C++ example".<br>
 *  \ref  py_mcumesh_getNodeIdsInUse "Here is a Python example".
 *  \endif
 * \sa computeFetchedNodeIds, computeNodeIdsAlg()
 */
DataArrayIdType *MEDCouplingUMesh::getNodeIdsInUse(mcIdType& nbrOfNodesInUse) const
{
  nbrOfNodesInUse=-1;
  mcIdType nbOfNodes(getNumberOfNodes());
  MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
  ret->alloc(nbOfNodes,1);
  mcIdType *traducer=ret->getPointer();
  std::fill(traducer,traducer+nbOfNodes,-1);
  mcIdType nbOfCells=getNumberOfCells();
  const mcIdType *connIndex=_nodal_connec_index->getConstPointer();
  const mcIdType *conn=_nodal_connec->getConstPointer();
  for(mcIdType i=0;i<nbOfCells;i++)
    for(mcIdType j=connIndex[i]+1;j<connIndex[i+1];j++)
      if(conn[j]>=0)
        {
          if(conn[j]<nbOfNodes)
            traducer[conn[j]]=1;
          else
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::getNodeIdsInUse : In cell #" << i  << " presence of node id " <<  conn[j] << " not in [0," << nbOfNodes << ") !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
        }
  nbrOfNodesInUse=ToIdType(std::count(traducer,traducer+nbOfNodes,1));
  std::transform(traducer,traducer+nbOfNodes,traducer,MEDCouplingAccVisit());
  return ret.retn();
}

/*!
 * This method returns a newly allocated array containing this->getNumberOfCells() tuples and 1 component.
 * For each cell in \b this the number of nodes constituting cell is computed.
 * For each polyhedron cell, the sum of the number of nodes of each face constituting polyhedron cell is returned.
 * So for pohyhedrons some nodes can be counted several times in the returned result.
 *
 * \return a newly allocated array
 * \sa MEDCouplingUMesh::computeEffectiveNbOfNodesPerCell
 */
DataArrayIdType *MEDCouplingUMesh::computeNbOfNodesPerCell() const
{
  checkConnectivityFullyDefined();
  mcIdType nbOfCells=getNumberOfCells();
  MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
  ret->alloc(nbOfCells,1);
  mcIdType *retPtr=ret->getPointer();
  const mcIdType *conn=getNodalConnectivity()->getConstPointer();
  const mcIdType *connI=getNodalConnectivityIndex()->getConstPointer();
  for(mcIdType i=0;i<nbOfCells;i++,retPtr++)
    {
      if(conn[connI[i]]!=ToIdType(INTERP_KERNEL::NORM_POLYHED))
        *retPtr=connI[i+1]-connI[i]-1;
      else
        *retPtr=connI[i+1]-connI[i]-1-ToIdType(std::count(conn+connI[i]+1,conn+connI[i+1],-1));
    }
  return ret.retn();
}

/*!
 * This method computes effective number of nodes per cell. That is to say nodes appearing several times in nodal connectivity of a cell,
 * will be counted only once here whereas it will be counted several times in MEDCouplingUMesh::computeNbOfNodesPerCell method.
 *
 * \return DataArrayIdType * - new object to be deallocated by the caller.
 * \sa MEDCouplingUMesh::computeNbOfNodesPerCell
 */
DataArrayIdType *MEDCouplingUMesh::computeEffectiveNbOfNodesPerCell() const
{
  checkConnectivityFullyDefined();
  mcIdType nbOfCells=getNumberOfCells();
  MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
  ret->alloc(nbOfCells,1);
  mcIdType *retPtr=ret->getPointer();
  const mcIdType *conn=getNodalConnectivity()->getConstPointer();
  const mcIdType *connI=getNodalConnectivityIndex()->getConstPointer();
  for(mcIdType i=0;i<nbOfCells;i++,retPtr++)
    {
      std::set<mcIdType> s(conn+connI[i]+1,conn+connI[i+1]);
      if(conn[connI[i]]!=ToIdType(INTERP_KERNEL::NORM_POLYHED))
        *retPtr=ToIdType(s.size());
      else
        {
          s.erase(-1);
          *retPtr=ToIdType(s.size());
        }
    }
  return ret.retn();
}

/*!
 * This method returns a newly allocated array containing this->getNumberOfCells() tuples and 1 component.
 * For each cell in \b this the number of faces constituting (entity of dimension this->getMeshDimension()-1) cell is computed.
 *
 * \return a newly allocated array
 */
DataArrayIdType *MEDCouplingUMesh::computeNbOfFacesPerCell() const
{
  checkConnectivityFullyDefined();
  mcIdType nbOfCells=getNumberOfCells();
  MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
  ret->alloc(nbOfCells,1);
  mcIdType *retPtr=ret->getPointer();
  const mcIdType *conn=getNodalConnectivity()->getConstPointer();
  const mcIdType *connI=getNodalConnectivityIndex()->getConstPointer();
  for(mcIdType i=0;i<nbOfCells;i++,retPtr++,connI++)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[*connI]);
      *retPtr=cm.getNumberOfSons2(conn+connI[0]+1,connI[1]-connI[0]-1);
    }
  return ret.retn();
}

/*!
 * Removes unused nodes (the node coordinates array is shorten) and returns an array
 * mapping between new and old node ids in "Old to New" mode. -1 values in the returned
 * array mean that the corresponding old node is no more used.
 *  \return DataArrayIdType * - a new instance of DataArrayIdType of length \a
 *           this->getNumberOfNodes() before call of this method. The caller is to
 *           delete this array using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If the nodal connectivity includes an invalid id.
 *  \sa areAllNodesFetched
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_zipCoordsTraducer "Here is a C++ example".<br>
 *  \ref  py_mcumesh_zipCoordsTraducer "Here is a Python example".
 *  \endif
 */
DataArrayIdType *MEDCouplingUMesh::zipCoordsTraducer()
{
  return MEDCouplingPointSet::zipCoordsTraducer();
}

/*!
 * This method stands if 'cell1' and 'cell2' are equals regarding 'compType' policy.
 * The semantic of 'compType' is specified in MEDCouplingPointSet::zipConnectivityTraducer method.
 */
int MEDCouplingUMesh::AreCellsEqual(const mcIdType *conn, const mcIdType *connI, mcIdType cell1, mcIdType cell2, int compType)
{
  switch(compType)
  {
    case 0:
      return AreCellsEqualPolicy0(conn,connI,cell1,cell2);
    case 1:
      return AreCellsEqualPolicy1(conn,connI,cell1,cell2);
    case 2:
      return AreCellsEqualPolicy2(conn,connI,cell1,cell2);
    case 3:
      return AreCellsEqualPolicy2NoType(conn,connI,cell1,cell2);
    case 7:
      return AreCellsEqualPolicy7(conn,connI,cell1,cell2);
  }
  throw INTERP_KERNEL::Exception("Unknown comparison asked ! Must be in 0,1,2,3 or 7.");
}

/*!
 * This method is the last step of the MEDCouplingPointSet::zipConnectivityTraducer with policy 0.
 */
int MEDCouplingUMesh::AreCellsEqualPolicy0(const mcIdType *conn, const mcIdType *connI, mcIdType cell1, mcIdType cell2)
{
  if(connI[cell1+1]-connI[cell1]==connI[cell2+1]-connI[cell2])
    return std::equal(conn+connI[cell1]+1,conn+connI[cell1+1],conn+connI[cell2]+1)?1:0;
  return 0;
}

/*!
 * This method is the last step of the MEDCouplingPointSet::zipConnectivityTraducer with policy 1.
 */
int MEDCouplingUMesh::AreCellsEqualPolicy1(const mcIdType *conn, const mcIdType *connI, mcIdType cell1, mcIdType cell2)
{
  mcIdType sz=connI[cell1+1]-connI[cell1];
  if(sz==connI[cell2+1]-connI[cell2])
    {
      if(conn[connI[cell1]]==conn[connI[cell2]])
        {
          const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[connI[cell1]]);
          unsigned dim=cm.getDimension();
          if(dim!=3)
            {
              if(dim!=1)
                {
                  mcIdType sz1=2*(sz-1);
                  INTERP_KERNEL::AutoPtr<mcIdType> tmp=new mcIdType[sz1];
                  mcIdType *work=std::copy(conn+connI[cell1]+1,conn+connI[cell1+1],(mcIdType *)tmp);
                  std::copy(conn+connI[cell1]+1,conn+connI[cell1+1],work);
                  work=std::search((mcIdType *)tmp,(mcIdType *)tmp+sz1,conn+connI[cell2]+1,conn+connI[cell2+1]);
                  return work!=tmp+sz1?1:0;
                }
              else
                return std::equal(conn+connI[cell1]+1,conn+connI[cell1+1],conn+connI[cell2]+1)?1:0;//case of SEG2 and SEG3
            }
          else
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::AreCellsEqualPolicy1 : not implemented yet for meshdim == 3 !");
        }
    }
  return 0;
}

/*!
 * This method is the last step of the MEDCouplingPointSet::zipConnectivityTraducer with policy 2.
 */
int MEDCouplingUMesh::AreCellsEqualPolicy2(const mcIdType *conn, const mcIdType *connI, mcIdType cell1, mcIdType cell2)
{
  if(connI[cell1+1]-connI[cell1]==connI[cell2+1]-connI[cell2])
    {
      if(conn[connI[cell1]]==conn[connI[cell2]])
        {
          std::set<mcIdType> s1(conn+connI[cell1]+1,conn+connI[cell1+1]);
          std::set<mcIdType> s2(conn+connI[cell2]+1,conn+connI[cell2+1]);
          return s1==s2?1:0;
        }
    }
  return 0;
}

/*!
 * This method is less restrictive than AreCellsEqualPolicy2. Here the geometric type is absolutely not taken into account !
 */
int MEDCouplingUMesh::AreCellsEqualPolicy2NoType(const mcIdType *conn, const mcIdType *connI, mcIdType cell1, mcIdType cell2)
{
  if(connI[cell1+1]-connI[cell1]==connI[cell2+1]-connI[cell2])
    {
      std::set<mcIdType> s1(conn+connI[cell1]+1,conn+connI[cell1+1]);
      std::set<mcIdType> s2(conn+connI[cell2]+1,conn+connI[cell2+1]);
      return s1==s2?1:0;
    }
  return 0;
}

/*!
 * This method is the last step of the MEDCouplingPointSet::zipConnectivityTraducer with policy 7.
 */
int MEDCouplingUMesh::AreCellsEqualPolicy7(const mcIdType *conn, const mcIdType *connI, mcIdType cell1, mcIdType cell2)
{
  mcIdType sz=connI[cell1+1]-connI[cell1];
  if(sz==connI[cell2+1]-connI[cell2])
    {
      if(conn[connI[cell1]]==conn[connI[cell2]])
        {
          const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[connI[cell1]]);
          unsigned dim=cm.getDimension();
          if(dim!=3)
            {
              if(dim!=1)
                {
                  mcIdType sz1=2*(sz-1);
                  INTERP_KERNEL::AutoPtr<mcIdType> tmp=new mcIdType[sz1];
                  mcIdType *work=std::copy(conn+connI[cell1]+1,conn+connI[cell1+1],(mcIdType *)tmp);
                  std::copy(conn+connI[cell1]+1,conn+connI[cell1+1],work);
                  work=std::search((mcIdType *)tmp,(mcIdType *)tmp+sz1,conn+connI[cell2]+1,conn+connI[cell2+1]);
                  if(work!=tmp+sz1)
                    return 1;
                  else
                    {
                      std::reverse_iterator<mcIdType *> it1((mcIdType *)tmp+sz1);
                      std::reverse_iterator<mcIdType *> it2((mcIdType *)tmp);
                      if(std::search(it1,it2,conn+connI[cell2]+1,conn+connI[cell2+1])!=it2)
                        return 2;
                      else
                        return 0;
                    }
                }
              else
                {//case of SEG2 and SEG3
                  if(std::equal(conn+connI[cell1]+1,conn+connI[cell1+1],conn+connI[cell2]+1))
                    return 1;
                  if(!cm.isQuadratic())
                    {
                      std::reverse_iterator<const mcIdType *> it1(conn+connI[cell1+1]);
                      std::reverse_iterator<const mcIdType *> it2(conn+connI[cell1]+1);
                      if(std::equal(it1,it2,conn+connI[cell2]+1))
                        return 2;
                      return 0;
                    }
                  else
                    {
                      if(conn[connI[cell1]+1]==conn[connI[cell2]+2] && conn[connI[cell1]+2]==conn[connI[cell2]+1] && conn[connI[cell1]+3]==conn[connI[cell2]+3])
                        return 2;
                      return 0;
                    }
                }
            }
          else
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::AreCellsEqualPolicy7 : not implemented yet for meshdim == 3 !");
        }
    }
  return 0;
}


/*!
 * This method find cells that are equal (regarding \a compType) in \a this. The comparison is specified by \a compType (see zipConnectivityTraducer).
 * This method keeps the coordinates of \a this. The comparison starts at rank \a startCellId cell id (included).
 * If \a startCellId is equal to 0 algorithm starts at cell #0 and for each cell candidates being searched have cell id higher than current cellId.
 * If \a startCellId is greater than 0 algorithm starts at cell #startCellId but for each cell all candidates are considered.
 * This method is time consuming.
 *
 * \param [in] compType input specifying the technique used to compare cells each other.
 *   - 0 : exactly. A cell is detected to be the same if and only if the connectivity is exactly the same without permutation and types same too. This is the strongest policy.
 *   - 1 : permutation same orientation. cell1 and cell2 are considered equal if the connectivity of cell2 can be deduced by those of cell1 by direct permutation (with exactly the same orientation)
 * and their type equal. For 1D mesh the policy 1 is equivalent to 0.
 *   - 2 : nodal. cell1 and cell2 are equal if and only if cell1 and cell2 have same type and have the same nodes constituting connectivity. This is the laziest policy. This policy
 * can be used for users not sensitive to orientation of cell
 * \param [in] startCellId specifies the cellId starting from which the equality computation will be carried out. By default it is 0, which it means that all cells in \a this will be scanned.
 * \param [out] commonCellsArr common cells ids (\ref numbering-indirect)
 * \param [out] commonCellsIArr common cells ids (\ref numbering-indirect)
 *
 */
void MEDCouplingUMesh::findCommonCells(int compType, mcIdType startCellId, DataArrayIdType *& commonCellsArr, DataArrayIdType *& commonCellsIArr) const
{
  MCAuto<DataArrayIdType> revNodal=DataArrayIdType::New(),revNodalI=DataArrayIdType::New();
  getReverseNodalConnectivity(revNodal,revNodalI);
  FindCommonCellsAlg(compType,startCellId,_nodal_connec,_nodal_connec_index,revNodal,revNodalI,commonCellsArr,commonCellsIArr);
}

void MEDCouplingUMesh::FindCommonCellsAlg(int compType, mcIdType startCellId, const DataArrayIdType *nodal, const DataArrayIdType *nodalI, const DataArrayIdType *revNodal, const DataArrayIdType *revNodalI,
                                          DataArrayIdType *& commonCellsArr, DataArrayIdType *& commonCellsIArr)
{
  MCAuto<DataArrayIdType> commonCells=DataArrayIdType::New(),commonCellsI=DataArrayIdType::New(); commonCells->alloc(0,1);
  mcIdType nbOfCells=nodalI->getNumberOfTuples()-1;
  commonCellsI->reserve(1); commonCellsI->pushBackSilent(0);
  const mcIdType *revNodalPtr=revNodal->getConstPointer(),*revNodalIPtr=revNodalI->getConstPointer();
  const mcIdType *connPtr=nodal->getConstPointer(),*connIPtr=nodalI->getConstPointer();
  std::vector<bool> isFetched(nbOfCells,false);
  if(startCellId==0)
    {
      for(mcIdType i=startCellId;i<nbOfCells;i++)
        {
          if(!isFetched[i])
            {
              const mcIdType *connOfNode=std::find_if(connPtr+connIPtr[i]+1,connPtr+connIPtr[i+1],std::bind(std::not_equal_to<mcIdType>(),std::placeholders::_1,-1));
              std::vector<mcIdType> v,v2;
              if(connOfNode!=connPtr+connIPtr[i+1])
                {
                  const mcIdType *locRevNodal=std::find(revNodalPtr+revNodalIPtr[*connOfNode],revNodalPtr+revNodalIPtr[*connOfNode+1],i);
                  v2.insert(v2.end(),locRevNodal,revNodalPtr+revNodalIPtr[*connOfNode+1]);
                  connOfNode++;
                }
              for(;connOfNode!=connPtr+connIPtr[i+1] && v2.size()>1;connOfNode++)
                if(*connOfNode>=0)
                  {
                    v=v2;
                    const mcIdType *locRevNodal=std::find(revNodalPtr+revNodalIPtr[*connOfNode],revNodalPtr+revNodalIPtr[*connOfNode+1],i);
                    std::vector<mcIdType>::iterator it=std::set_intersection(v.begin(),v.end(),locRevNodal,revNodalPtr+revNodalIPtr[*connOfNode+1],v2.begin());
                    v2.resize(std::distance(v2.begin(),it));
                  }
              if(v2.size()>1)
                {
                  if(AreCellsEqualInPool(v2,compType,connPtr,connIPtr,commonCells))
                    {
                      mcIdType pos=commonCellsI->back();
                      commonCellsI->pushBackSilent(commonCells->getNumberOfTuples());
                      for(const mcIdType *it=commonCells->begin()+pos;it!=commonCells->end();it++)
                        isFetched[*it]=true;
                    }
                }
            }
        }
    }
  else
    {
      for(mcIdType i=startCellId;i<nbOfCells;i++)
        {
          if(!isFetched[i])
            {
              const mcIdType *connOfNode=std::find_if(connPtr+connIPtr[i]+1,connPtr+connIPtr[i+1],std::bind(std::not_equal_to<mcIdType>(),std::placeholders::_1,-1));
              // v2 contains the result of successive intersections using rev nodal on on each node of cell #i
              std::vector<mcIdType> v,v2;
              if(connOfNode!=connPtr+connIPtr[i+1])
                {
                  v2.insert(v2.end(),revNodalPtr+revNodalIPtr[*connOfNode],revNodalPtr+revNodalIPtr[*connOfNode+1]);
                  connOfNode++;
                }
              for(;connOfNode!=connPtr+connIPtr[i+1] && v2.size()>1;connOfNode++)
                if(*connOfNode>=0)
                  {
                    v=v2;
                    std::vector<mcIdType>::iterator it=std::set_intersection(v.begin(),v.end(),revNodalPtr+revNodalIPtr[*connOfNode],revNodalPtr+revNodalIPtr[*connOfNode+1],v2.begin());
                    v2.resize(std::distance(v2.begin(),it));
                  }
              // v2 contains now candidates. Problem candidates are sorted using id rank.
              if(v2.size()>1)
                {
                  if(v2[0]!=i)
                  {
                    auto it(std::find(v2.begin(),v2.end(),i));
                    std::swap(*v2.begin(),*it);
                  }
                  if(AreCellsEqualInPool(v2,compType,connPtr,connIPtr,commonCells))
                    {
                      mcIdType newPos(commonCells->getNumberOfTuples());
                      mcIdType pos(commonCellsI->back());
                      std::sort(commonCells->getPointerSilent()+pos,commonCells->getPointerSilent()+newPos);
                      commonCellsI->pushBackSilent(newPos);
                      for(const mcIdType *it=commonCells->begin()+pos;it!=commonCells->end();it++)
                        isFetched[*it]=true;
                    }
                }
            }
        }
    }
  commonCellsArr=commonCells.retn();
  commonCellsIArr=commonCellsI.retn();
}

/*!
 * Checks if \a this mesh includes all cells of an \a other mesh, and returns an array
 * giving for each cell of the \a other an id of a cell in \a this mesh. A value larger
 * than \a this->getNumberOfCells() in the returned array means that there is no
 * corresponding cell in \a this mesh.
 * It is expected that \a this and \a other meshes share the same node coordinates
 * array, if it is not so an exception is thrown.
 *  \param [in] other - the mesh to compare with.
 *  \param [in] compType - specifies a cell comparison technique. For meaning of its
 *         valid values [0,1,2], see zipConnectivityTraducer().
 *  \param [out] arr - a new instance of DataArrayIdType returning correspondence
 *         between cells of the two meshes. It contains \a other->getNumberOfCells()
 *         values. The caller is to delete this array using
 *         decrRef() as it is no more needed.
 *  \return bool - \c true if all cells of \a other mesh are present in the \a this
 *         mesh.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_areCellsIncludedIn "Here is a C++ example".<br>
 *  \ref  py_mcumesh_areCellsIncludedIn "Here is a Python example".
 *  \endif
 *  \sa checkDeepEquivalOnSameNodesWith()
 *  \sa checkGeoEquivalWith()
 */
bool MEDCouplingUMesh::areCellsIncludedIn(const MEDCouplingUMesh *other, int compType, DataArrayIdType *& arr) const
{
  MCAuto<MEDCouplingUMesh> mesh=MergeUMeshesOnSameCoords(this,other);
  mcIdType nbOfCells=getNumberOfCells();
  static const int possibleCompType[]={0,1,2};
  if(std::find(possibleCompType,possibleCompType+sizeof(possibleCompType)/sizeof(int),compType)==possibleCompType+sizeof(possibleCompType)/sizeof(int))
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::areCellsIncludedIn : only following policies are possible : ";
      std::copy(possibleCompType,possibleCompType+sizeof(possibleCompType)/sizeof(int),std::ostream_iterator<int>(oss," "));
      oss << " !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  //
  if(other->getNumberOfCells()==0)
  {
    MCAuto<DataArrayIdType> dftRet(DataArrayIdType::New()); dftRet->alloc(0,1); arr=dftRet.retn(); arr->setName(other->getName());
    return true;
  }
  DataArrayIdType *commonCells(nullptr),*commonCellsI(nullptr);
  mesh->findCommonCells(compType,nbOfCells,commonCells,commonCellsI);
  MCAuto<DataArrayIdType> commonCellsTmp(commonCells),commonCellsITmp(commonCellsI);
  mcIdType newNbOfCells=-1;
  MCAuto<DataArrayIdType> o2n = DataArrayIdType::ConvertIndexArrayToO2N(ToIdType(mesh->getNumberOfCells()),commonCells->begin(),commonCellsI->begin(),commonCellsI->end(),newNbOfCells);
  MCAuto<DataArrayIdType> p0(o2n->selectByTupleIdSafeSlice(0,nbOfCells,1));
  mcIdType maxPart(p0->getMaxValueInArray());
  bool ret(maxPart==newNbOfCells-1);
  MCAuto<DataArrayIdType> p1(p0->invertArrayO2N2N2O(newNbOfCells));
  // fill p1 array in case of presence of cells in other not in this
  mcIdType *pt(p1->getPointer());
  for(mcIdType i = maxPart ; i < newNbOfCells-1 ; ++i )
    pt[i+1] = i+1;
  //
  MCAuto<DataArrayIdType> p2(o2n->subArray(nbOfCells));
  p2->transformWithIndArr(p1->begin(),p1->end()); p2->setName(other->getName());
  arr = p2.retn();
  return ret;
}

/*!
 * This method makes the assumption that \a this and \a other share the same coords. If not an exception will be thrown !
 * This method tries to determine if \b other is fully included in \b this.
 * The main difference is that this method is not expected to throw exception.
 * This method has two outputs :
 *
 * \param other other mesh
 * \param arr is an output parameter that returns a \b newly created instance. This array is of size 'other->getNumberOfCells()'.
 * \return If \a other is fully included in 'this 'true is returned. If not false is returned.
 */
bool MEDCouplingUMesh::areCellsIncludedInPolicy7(const MEDCouplingUMesh *other, DataArrayIdType *& arr) const
{
  MCAuto<MEDCouplingUMesh> mesh=MergeUMeshesOnSameCoords(this,other);
  DataArrayIdType *commonCells=0,*commonCellsI=0;
  mcIdType thisNbCells=getNumberOfCells();
  mesh->findCommonCells(7,thisNbCells,commonCells,commonCellsI);
  MCAuto<DataArrayIdType> commonCellsTmp(commonCells),commonCellsITmp(commonCellsI);
  const mcIdType *commonCellsPtr=commonCells->getConstPointer(),*commonCellsIPtr=commonCellsI->getConstPointer();
  mcIdType otherNbCells=other->getNumberOfCells();
  MCAuto<DataArrayIdType> arr2=DataArrayIdType::New();
  arr2->alloc(otherNbCells,1);
  arr2->fillWithZero();
  mcIdType *arr2Ptr=arr2->getPointer();
  mcIdType nbOfCommon=commonCellsI->getNumberOfTuples()-1;
  for(mcIdType i=0;i<nbOfCommon;i++)
    {
      mcIdType start=commonCellsPtr[commonCellsIPtr[i]];
      if(start<thisNbCells)
        {
          for(mcIdType j=commonCellsIPtr[i]+1;j!=commonCellsIPtr[i+1];j++)
            {
              mcIdType sig=commonCellsPtr[j]>0?1:-1;
              mcIdType val=std::abs(commonCellsPtr[j])-1;
              if(val>=thisNbCells)
                arr2Ptr[val-thisNbCells]=sig*(start+1);
            }
        }
    }
  arr2->setName(other->getName());
  if(arr2->presenceOfValue(0))
    return false;
  arr=arr2.retn();
  return true;
}

MEDCouplingUMesh *MEDCouplingUMesh::mergeMyselfWithOnSameCoords(const MEDCouplingPointSet *other) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::mergeMyselfWithOnSameCoords : input other is null !");
  const MEDCouplingUMesh *otherC=dynamic_cast<const MEDCouplingUMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::mergeMyselfWithOnSameCoords : the input other mesh is not of type unstructured !");
  std::vector<const MEDCouplingUMesh *> ms(2);
  ms[0]=this;
  ms[1]=otherC;
  return MergeUMeshesOnSameCoords(ms);
}

/*!
 * Build a sub part of \b this lying or not on the same coordinates than \b this (regarding value of \b keepCoords).
 * By default coordinates are kept. This method is close to MEDCouplingUMesh::buildPartOfMySelf except that here input
 * cellIds is not given explicitly but by a range python like.
 *
 * \param start starting ID
 * \param end end ID (excluded)
 * \param step step size
 * \param keepCoords that specifies if you want or not to keep coords as this or zip it (see MEDCoupling::MEDCouplingUMesh::zipCoords). If true zipCoords is \b NOT called, if false, zipCoords is called.
 * \return a newly allocated
 *
 * \warning This method modifies can generate an unstructured mesh whose cells are not sorted by geometric type order.
 * In view of the MED file writing, a renumbering of cells of returned unstructured mesh (using MEDCouplingUMesh::sortCellsInMEDFileFrmt) should be necessary.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildPartOfMySelfSlice(mcIdType start, mcIdType end, mcIdType step, bool keepCoords) const
{
  if(getMeshDimension()!=-1)
    return static_cast<MEDCouplingUMesh *>(MEDCouplingPointSet::buildPartOfMySelfSlice(start,end,step,keepCoords));
  else
    {
      mcIdType newNbOfCells=DataArray::GetNumberOfItemGivenBESRelative(start,end,step,"MEDCouplingUMesh::buildPartOfMySelfSlice for -1 dimension mesh ");
      if(newNbOfCells!=1)
        throw INTERP_KERNEL::Exception("-1D mesh has only one cell !");
      if(start!=0)
        throw INTERP_KERNEL::Exception("-1D mesh has only one cell : 0 !");
      incrRef();
      return const_cast<MEDCouplingUMesh *>(this);
    }
}

/*!
 * Creates a new MEDCouplingUMesh containing specified cells of \a this mesh.
 * The result mesh shares or not the node coordinates array with \a this mesh depending
 * on \a keepCoords parameter.
 *  \warning Cells of the result mesh can be \b not sorted by geometric type, hence,
 *           to write this mesh to the MED file, its cells must be sorted using
 *           sortCellsInMEDFileFrmt().
 *  \param [in] begin - an array of cell ids to include to the new mesh.
 *  \param [in] end - a pointer to last-plus-one-th element of \a begin.
 *  \param [in] keepCoords - if \c true, the result mesh shares the node coordinates
 *         array of \a this mesh, else "free" nodes are removed from the result mesh
 *         by calling zipCoords().
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. The caller is
 *         to delete this mesh using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If any cell id in the array \a begin is not valid.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_buildPartOfMySelf "Here is a C++ example".<br>
 *  \ref  py_mcumesh_buildPartOfMySelf "Here is a Python example".
 *  \endif
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildPartOfMySelf(const mcIdType *begin, const mcIdType *end, bool keepCoords) const
{
  if(getMeshDimension()!=-1)
    return static_cast<MEDCouplingUMesh *>(MEDCouplingPointSet::buildPartOfMySelf(begin,end,keepCoords));
  else
    {
      if(end-begin!=1)
        throw INTERP_KERNEL::Exception("-1D mesh has only one cell !");
      if(begin[0]!=0)
        throw INTERP_KERNEL::Exception("-1D mesh has only one cell : 0 !");
      incrRef();
      return const_cast<MEDCouplingUMesh *>(this);
    }
}

/*!
 * This method operates only on nodal connectivity on \b this. Coordinates of \b this is completely ignored here.
 *
 * This method allows to partially modify some cells in \b this (whose list is specified by [ \b cellIdsBg, \b cellIdsEnd ) ) with cells coming in \b otherOnSameCoordsThanThis.
 * Size of [ \b cellIdsBg, \b cellIdsEnd ) ) must be equal to the number of cells of otherOnSameCoordsThanThis.
 * The number of cells of \b this will remain the same with this method.
 *
 * \param [in] cellIdsBg begin of cell ids (included) of cells in this to assign
 * \param [in] cellIdsEnd end of cell ids (excluded) of cells in this to assign
 * \param [in] otherOnSameCoordsThanThis an another mesh with same meshdimension than \b this with exactly the same number of cells than cell ids list in [\b cellIdsBg, \b cellIdsEnd ).
 *             Coordinate pointer of \b this and those of \b otherOnSameCoordsThanThis must be the same
 */
void MEDCouplingUMesh::setPartOfMySelf(const mcIdType *cellIdsBg, const mcIdType *cellIdsEnd, const MEDCouplingUMesh& otherOnSameCoordsThanThis)
{
  checkConnectivityFullyDefined();
  otherOnSameCoordsThanThis.checkConnectivityFullyDefined();
  if(getCoords()!=otherOnSameCoordsThanThis.getCoords())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::setPartOfMySelf : coordinates pointer are not the same ! Invoke setCoords or call tryToShareSameCoords method !");
  if(getMeshDimension()!=otherOnSameCoordsThanThis.getMeshDimension())
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::setPartOfMySelf : Mismatch of meshdimensions ! this is equal to " << getMeshDimension();
      oss << ", whereas other mesh dimension is set equal to " << otherOnSameCoordsThanThis.getMeshDimension() << " !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  mcIdType nbOfCellsToModify( ToIdType((std::distance(cellIdsBg,cellIdsEnd))));
  if(nbOfCellsToModify!=otherOnSameCoordsThanThis.getNumberOfCells())
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::setPartOfMySelf : cells ids length (" <<  nbOfCellsToModify << ") do not match the number of cells of other mesh (" << otherOnSameCoordsThanThis.getNumberOfCells() << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  mcIdType nbOfCells(getNumberOfCells());
  bool easyAssign(true);
  const mcIdType *connI(_nodal_connec_index->begin());
  const mcIdType *connIOther=otherOnSameCoordsThanThis._nodal_connec_index->begin();
  for(const mcIdType *it=cellIdsBg;it!=cellIdsEnd && easyAssign;it++,connIOther++)
    {
      if(*it>=0 && *it<nbOfCells)
        {
          easyAssign=(connIOther[1]-connIOther[0])==(connI[*it+1]-connI[*it]);
        }
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::setPartOfMySelf : On pos #" << std::distance(cellIdsBg,it) << " id is equal to " << *it << " which is not in [0," << nbOfCells << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  if(easyAssign)
    {
      DataArrayIdType::SetPartOfIndexedArraysSameIdx(cellIdsBg,cellIdsEnd,_nodal_connec,_nodal_connec_index,otherOnSameCoordsThanThis._nodal_connec,otherOnSameCoordsThanThis._nodal_connec_index);
      computeTypes();
    }
  else
    {
      DataArrayIdType *arrOut=0,*arrIOut=0;
      DataArrayIdType::SetPartOfIndexedArrays(cellIdsBg,cellIdsEnd,_nodal_connec,_nodal_connec_index,otherOnSameCoordsThanThis._nodal_connec,otherOnSameCoordsThanThis._nodal_connec_index,
                                               arrOut,arrIOut);
      MCAuto<DataArrayIdType> arrOutAuto(arrOut),arrIOutAuto(arrIOut);
      setConnectivity(arrOut,arrIOut,true);
    }
}

void MEDCouplingUMesh::setPartOfMySelfSlice(mcIdType start, mcIdType end, mcIdType step, const MEDCouplingUMesh& otherOnSameCoordsThanThis)
{
  checkConnectivityFullyDefined();
  otherOnSameCoordsThanThis.checkConnectivityFullyDefined();
  if(getCoords()!=otherOnSameCoordsThanThis.getCoords())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::setPartOfMySelfSlice : coordinates pointer are not the same ! Invoke setCoords or call tryToShareSameCoords method !");
  if(getMeshDimension()!=otherOnSameCoordsThanThis.getMeshDimension())
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::setPartOfMySelfSlice : Mismatch of meshdimensions ! this is equal to " << getMeshDimension();
      oss << ", whereas other mesh dimension is set equal to " << otherOnSameCoordsThanThis.getMeshDimension() << " !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  mcIdType nbOfCellsToModify=DataArray::GetNumberOfItemGivenBESRelative(start,end,step,"MEDCouplingUMesh::setPartOfMySelfSlice : ");
  if(nbOfCellsToModify!=otherOnSameCoordsThanThis.getNumberOfCells())
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::setPartOfMySelfSlice : cells ids length (" <<  nbOfCellsToModify << ") do not match the number of cells of other mesh (" << otherOnSameCoordsThanThis.getNumberOfCells() << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  mcIdType nbOfCells=getNumberOfCells();
  bool easyAssign=true;
  const mcIdType *connI=_nodal_connec_index->getConstPointer();
  const mcIdType *connIOther=otherOnSameCoordsThanThis._nodal_connec_index->getConstPointer();
  mcIdType it=start;
  for(mcIdType i=0;i<nbOfCellsToModify && easyAssign;i++,it+=step,connIOther++)
    {
      if(it>=0 && it<nbOfCells)
        {
          easyAssign=(connIOther[1]-connIOther[0])==(connI[it+1]-connI[it]);
        }
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::setPartOfMySelfSlice : On pos #" << i << " id is equal to " << it << " which is not in [0," << nbOfCells << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  if(easyAssign)
    {
      DataArrayIdType::SetPartOfIndexedArraysSameIdxSlice(start,end,step,_nodal_connec,_nodal_connec_index,otherOnSameCoordsThanThis._nodal_connec,otherOnSameCoordsThanThis._nodal_connec_index);
      computeTypes();
    }
  else
    {
      DataArrayIdType *arrOut=0,*arrIOut=0;
      DataArrayIdType::SetPartOfIndexedArraysSlice(start,end,step,_nodal_connec,_nodal_connec_index,otherOnSameCoordsThanThis._nodal_connec,otherOnSameCoordsThanThis._nodal_connec_index,
                                                arrOut,arrIOut);
      MCAuto<DataArrayIdType> arrOutAuto(arrOut),arrIOutAuto(arrIOut);
      setConnectivity(arrOut,arrIOut,true);
    }
}


/*!
 * Creates a new MEDCouplingUMesh containing cells, of dimension one less than \a
 * this->getMeshDimension(), that bound some cells of \a this mesh.
 * The cells of lower dimension to include to the result mesh are selected basing on
 * specified node ids and the value of \a fullyIn parameter. If \a fullyIn ==\c true, a
 * cell is copied if its all nodes are in the array \a begin of node ids. If \a fullyIn
 * ==\c false, a cell is copied if any its node is in the array of node ids. The
 * created mesh shares the node coordinates array with \a this mesh.
 *  \param [in] begin - the array of node ids.
 *  \param [in] end - a pointer to the (last+1)-th element of \a begin.
 *  \param [in] fullyIn - if \c true, then cells whose all nodes are in the
 *         array \a begin are added, else cells whose any node is in the
 *         array \a begin are added.
 *  \return MEDCouplingUMesh * - new instance of MEDCouplingUMesh. The caller is
 *         to delete this mesh using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If any node id in \a begin is not valid.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_buildFacePartOfMySelfNode "Here is a C++ example".<br>
 *  \ref  py_mcumesh_buildFacePartOfMySelfNode "Here is a Python example".
 *  \endif
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildFacePartOfMySelfNode(const mcIdType *begin, const mcIdType *end, bool fullyIn) const
{
  MCAuto<DataArrayIdType> desc,descIndx,revDesc,revDescIndx;
  desc=DataArrayIdType::New(); descIndx=DataArrayIdType::New(); revDesc=DataArrayIdType::New(); revDescIndx=DataArrayIdType::New();
  MCAuto<MEDCouplingUMesh> subMesh=buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  desc=0; descIndx=0; revDesc=0; revDescIndx=0;
  return static_cast<MEDCouplingUMesh*>(subMesh->buildPartOfMySelfNode(begin,end,fullyIn));
}

/*!
 * Creates a new MEDCouplingUMesh containing cells, of dimension one less than \a
 * this->getMeshDimension(), which bound only one cell of \a this mesh.
 *  \param [in] keepCoords - if \c true, the result mesh shares the node coordinates
 *         array of \a this mesh, else "free" nodes are removed from the result mesh
 *         by calling zipCoords().
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. The caller is
 *         to delete this mesh using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_buildBoundaryMesh "Here is a C++ example".<br>
 *  \ref  py_mcumesh_buildBoundaryMesh "Here is a Python example".
 *  \endif
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildBoundaryMesh(bool keepCoords) const
{
  DataArrayIdType *desc=DataArrayIdType::New();
  DataArrayIdType *descIndx=DataArrayIdType::New();
  DataArrayIdType *revDesc=DataArrayIdType::New();
  DataArrayIdType *revDescIndx=DataArrayIdType::New();
  //
  MCAuto<MEDCouplingUMesh> meshDM1=buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  revDesc->decrRef();
  desc->decrRef();
  descIndx->decrRef();
  mcIdType nbOfCells=meshDM1->getNumberOfCells();
  const mcIdType *revDescIndxC=revDescIndx->getConstPointer();
  std::vector<mcIdType> boundaryCells;
  for(mcIdType i=0;i<nbOfCells;i++)
    if(revDescIndxC[i+1]-revDescIndxC[i]==1)
      boundaryCells.push_back(i);
  revDescIndx->decrRef();
  MEDCouplingUMesh *ret=meshDM1->buildPartOfMySelf(&boundaryCells[0],&boundaryCells[0]+boundaryCells.size(),keepCoords);
  return ret;
}

/*!
 * This method returns a newly created DataArrayIdType instance containing ids of cells located in boundary.
 * A cell is detected to be on boundary if it contains one or more than one face having only one father.
 * This method makes the assumption that \a this is fully defined (coords,connectivity). If not an exception will be thrown.
 */
DataArrayIdType *MEDCouplingUMesh::findCellIdsOnBoundary() const
{
  checkFullyDefined();
  MCAuto<DataArrayIdType> ret2(DataArrayIdType::New());

  if (getNumberOfCells() == 0)
    {
      ret2->alloc(0,1);
      return ret2.retn();
    }

  MCAuto<DataArrayIdType> desc(DataArrayIdType::New()), descIndx(DataArrayIdType::New()), revDesc(DataArrayIdType::New()), revDescIndx(DataArrayIdType::New());
  //
  buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx)->decrRef();
  desc=(DataArrayIdType*)0; descIndx=(DataArrayIdType*)0;
  //
  MCAuto<DataArrayIdType> tmp=revDescIndx->deltaShiftIndex();
  MCAuto<DataArrayIdType> faceIds=tmp->findIdsEqual(1); tmp=(DataArrayIdType*)0;
  const mcIdType *revDescPtr=revDesc->getConstPointer();
  const mcIdType *revDescIndxPtr=revDescIndx->getConstPointer();
  mcIdType nbOfCells=getNumberOfCells();
  std::vector<bool> ret1(nbOfCells,false);
  mcIdType sz=0;
  for(const mcIdType *pt=faceIds->begin();pt!=faceIds->end();pt++)
    if(!ret1[revDescPtr[revDescIndxPtr[*pt]]])
      { ret1[revDescPtr[revDescIndxPtr[*pt]]]=true; sz++; }
  //
  ret2->alloc(sz,1);
  mcIdType *ret2Ptr=ret2->getPointer();
  sz=0;
  for(std::vector<bool>::const_iterator it=ret1.begin();it!=ret1.end();it++,sz++)
    if(*it)
      *ret2Ptr++=sz;
  ret2->setName("BoundaryCells");
  return ret2.retn();
}

/*!
 * This method finds in \b this the cell ids that lie on mesh \b otherDimM1OnSameCoords.
 * \b this and \b otherDimM1OnSameCoords have to lie on the same coordinate array pointer. The coherency of that coords array with connectivity
 * of \b this and \b otherDimM1OnSameCoords is not important here because this method works only on connectivity.
 * this->getMeshDimension() - 1 must be equal to otherDimM1OnSameCoords.getMeshDimension()
 *
 * s0 is the cell ids set in \b this lying on at least one node in the fetched nodes in \b otherDimM1OnSameCoords.
 * This method also returns the cells ids set s1 which contains the cell ids in \b this for which one of the dim-1 constituent
 * equals a cell in \b otherDimM1OnSameCoords.
 *
 * \throw if \b otherDimM1OnSameCoords is not part of constituent of \b this, or if coordinate pointer of \b this and \b otherDimM1OnSameCoords
 *        are not same, or if this->getMeshDimension()-1!=otherDimM1OnSameCoords.getMeshDimension()
 *
 * \param [in] otherDimM1OnSameCoords other mesh
 * \param [out] cellIdsRk0 a newly allocated array containing the cell ids of s0 (which are cell ids of \b this) in the above algorithm.
 * \param [out] cellIdsRk1 a newly allocated array containing the cell ids of s1 \b indexed into the \b cellIdsRk0 subset. To get the absolute ids of s1, simply invoke
 *              cellIdsRk1->transformWithIndArr(cellIdsRk0->begin(),cellIdsRk0->end());
 */
void MEDCouplingUMesh::findCellIdsLyingOn(const MEDCouplingUMesh& otherDimM1OnSameCoords, DataArrayIdType *&cellIdsRk0, DataArrayIdType *&cellIdsRk1) const
{
  if(getCoords()!=otherDimM1OnSameCoords.getCoords())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findCellIdsLyingOn : coordinates pointer are not the same ! Use tryToShareSameCoords method !");
  checkConnectivityFullyDefined();
  otherDimM1OnSameCoords.checkConnectivityFullyDefined();
  if(getMeshDimension()-1!=otherDimM1OnSameCoords.getMeshDimension())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findCellIdsLyingOn : invalid mesh dimension of input mesh regarding meshdimesion of this !");
  MCAuto<DataArrayIdType> fetchedNodeIds1=otherDimM1OnSameCoords.computeFetchedNodeIds();
  MCAuto<DataArrayIdType> s0arr=getCellIdsLyingOnNodes(fetchedNodeIds1->begin(),fetchedNodeIds1->end(),false);
  MCAuto<MEDCouplingUMesh> thisPart=static_cast<MEDCouplingUMesh *>(buildPartOfMySelf(s0arr->begin(),s0arr->end(),true));
  MCAuto<DataArrayIdType> descThisPart=DataArrayIdType::New(),descIThisPart=DataArrayIdType::New(),revDescThisPart=DataArrayIdType::New(),revDescIThisPart=DataArrayIdType::New();
  MCAuto<MEDCouplingUMesh> thisPartConsti=thisPart->buildDescendingConnectivity(descThisPart,descIThisPart,revDescThisPart,revDescIThisPart);
  const mcIdType *revDescThisPartPtr=revDescThisPart->getConstPointer(),*revDescIThisPartPtr=revDescIThisPart->getConstPointer();
  DataArrayIdType *idsOtherInConsti=0;
  bool b=thisPartConsti->areCellsIncludedIn(&otherDimM1OnSameCoords,2,idsOtherInConsti);
  MCAuto<DataArrayIdType> idsOtherInConstiAuto(idsOtherInConsti);
  if(!b)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findCellIdsLyingOn : the given mdim-1 mesh in other is not a constituent of this !");
  std::set<mcIdType> s1;
  for(const mcIdType *idOther=idsOtherInConsti->begin();idOther!=idsOtherInConsti->end();idOther++)
    s1.insert(revDescThisPartPtr+revDescIThisPartPtr[*idOther],revDescThisPartPtr+revDescIThisPartPtr[*idOther+1]);
  MCAuto<DataArrayIdType> s1arr_renum1=DataArrayIdType::New(); s1arr_renum1->alloc(s1.size(),1); std::copy(s1.begin(),s1.end(),s1arr_renum1->getPointer());
  s1arr_renum1->sort();
  cellIdsRk0=s0arr.retn();
  //cellIdsRk1=s_renum1.retn();
  cellIdsRk1=s1arr_renum1.retn();
}

/*!
 * This method computes the skin of \b this. That is to say the consituting meshdim-1 mesh is built and only the boundary subpart is
 * returned. This subpart of meshdim-1 mesh is built using meshdim-1 cells in it shared only one cell in \b this.
 *
 * \return a newly allocated mesh lying on the same coordinates than \b this. The caller has to deal with returned mesh.
 */
MEDCouplingUMesh *MEDCouplingUMesh::computeSkin() const
{
  MCAuto<DataArrayIdType> desc=DataArrayIdType::New();
  MCAuto<DataArrayIdType> descIndx=DataArrayIdType::New();
  MCAuto<DataArrayIdType> revDesc=DataArrayIdType::New();
  MCAuto<DataArrayIdType> revDescIndx=DataArrayIdType::New();
  //
  MCAuto<MEDCouplingUMesh> meshDM1=buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  revDesc=0; desc=0; descIndx=0;
  MCAuto<DataArrayIdType> revDescIndx2=revDescIndx->deltaShiftIndex();
  MCAuto<DataArrayIdType> part=revDescIndx2->findIdsEqual(1);
  return static_cast<MEDCouplingUMesh *>(meshDM1->buildPartOfMySelf(part->begin(),part->end(),true));
}

/*!
 * Finds nodes lying on the boundary of \a this mesh.
 *  \return DataArrayIdType * - a new instance of DataArrayIdType holding ids of found
 *          nodes. The caller is to delete this array using decrRef() as it is no
 *          more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is node defined.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_findBoundaryNodes "Here is a C++ example".<br>
 *  \ref  py_mcumesh_findBoundaryNodes "Here is a Python example".
 *  \endif
 */
DataArrayIdType *MEDCouplingUMesh::findBoundaryNodes() const
{
  MCAuto<MEDCouplingUMesh> skin=computeSkin();
  return skin->computeFetchedNodeIds();
}

MEDCouplingUMesh *MEDCouplingUMesh::buildUnstructured() const
{
  incrRef();
  return const_cast<MEDCouplingUMesh *>(this);
}

/*!
 * This method expects that \b this and \b otherDimM1OnSameCoords share the same coordinates array.
 * otherDimM1OnSameCoords->getMeshDimension() is expected to be equal to this->getMeshDimension()-1.
 * This method searches for nodes needed to be duplicated. Those are the nodes fetched by \b otherDimM1OnSameCoords which are not part of the boundary of \b otherDimM1OnSameCoords.
 * If a node is in the boundary of \b this \b and in the boundary of \b otherDimM1OnSameCoords, it will be duplicated.
 * When the set of node ids \b nodeIdsToDuplicate is computed, cell ids in \b this is searched so that their connectivity includes at least 1 node in \b nodeIdsToDuplicate.
 *
 * \param [in] crackingMesh a mesh lying on the same coords than \b this and with a mesh dimension equal to those of \b this minus 1. WARNING this input
 *             parameter is altered during the call.
 * \return node ids which need to be duplicated following the algorithm explained above.
 *
 */
DataArrayIdType* MEDCouplingUMesh::findNodesToDuplicate(const MEDCouplingUMesh& crackingMesh) const
{
  // DEBUG NOTE: in case of issue with the algorithm in this method, see Python script in resources/dev
  // which mimicks the C++
  using DAInt = MCAuto<DataArrayIdType>;
  using MCUMesh = MCAuto<MEDCouplingUMesh>;

  checkFullyDefined();
  crackingMesh.checkFullyDefined();
  if(getCoords()!=crackingMesh.getCoords())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findNodesToDuplicate : meshes do not share the same coords array !");
  if(crackingMesh.getMeshDimension()!=getMeshDimension()-1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findNodesToDuplicate : the mesh given in other parameter must have this->getMeshDimension()-1 !");

  // Clean the M1 group (cracking mesh): the M1 cells which are part of M0 boundary are irrelevant (we can't create a crack on the boundary of M0!)
  MCUMesh m0skin = computeSkin();
  DataArrayIdType *idsToKeepP;
  m0skin->areCellsIncludedIn(&crackingMesh,2, idsToKeepP);
  DAInt idsToKeep(idsToKeepP);
  DAInt ids2 = idsToKeep->findIdsNotInRange(0, m0skin->getNumberOfCells());  // discard cells on the skin of M0
  MCUMesh otherDimM1OnSameCoords =static_cast<MEDCouplingUMesh *>(crackingMesh.buildPartOfMySelf(ids2->begin(), ids2->end(), true));

  if (!otherDimM1OnSameCoords->getNumberOfCells())
    return MCAuto<DataArrayIdType>(DataArrayIdType::New()).retn();

  // Checking star-shaped M1 group:
  DAInt dt0=DataArrayIdType::New(),dit0=DataArrayIdType::New(),rdt0=DataArrayIdType::New(),rdit0=DataArrayIdType::New();
  MCUMesh meshM2 = otherDimM1OnSameCoords->buildDescendingConnectivity(dt0, dit0, rdt0, rdit0); // 2D: a mesh of points, 3D: a mesh of segs
  DAInt dsi = rdit0->deltaShiftIndex();
  DAInt idsTmp0 = dsi->findIdsNotInRange(-1, 3);  // for 2D: if a point is connected to more than 2 segs. For 3D: if a seg is connected to more than two faces.
  if(idsTmp0->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("MEDFileUMesh::buildInnerBoundaryAlongM1Group: group is too complex: some points (or edges) have more than two connected segments (or faces)!");
  dt0=0; dit0=0; rdt0=0; rdit0=0; idsTmp0=0;

  // Get extreme nodes from the group (they won't be duplicated except if they also lie on bound of M0 -- see below),
  // ie nodes belonging to the boundary "cells" (might be points) of M1
  DAInt xtremIdsM2 = dsi->findIdsEqual(1); dsi = 0;
  MCUMesh meshM2Part = static_cast<MEDCouplingUMesh *>(meshM2->buildPartOfMySelf(xtremIdsM2->begin(), xtremIdsM2->end(),true));
  DAInt xtrem = meshM2Part->computeFetchedNodeIds();
  // Remove from the list points on the boundary of the M0 mesh (those need duplication!).
  //    In 3D, some points on the boundary of M0 will NOT be duplicated (where as in 2D, points on the boundary of M0 are always duplicated)
  //    Think of a partial (plane) crack in a cube: the points at the tip of the crack and not located inside the volume of the cube are not duplicated
  //    although they are technically on the skin of the cube.
  DAInt fNodes = m0skin->computeFetchedNodeIds();
  DAInt notDup = 0;
  if (getMeshDimension() == 3)
    {
      DAInt dnu1=DataArrayIdType::New(), dnu2=DataArrayIdType::New(), dnu3=DataArrayIdType::New(), dnu4=DataArrayIdType::New();
      MCUMesh m0skinDesc = m0skin->buildDescendingConnectivity(dnu1, dnu2, dnu3, dnu4); // all segments of the skin of the 3D (M0) mesh
      dnu1=0;dnu2=0;dnu3=0;dnu4=0;
      DataArrayIdType * corresp=0;
      meshM2->areCellsIncludedIn(m0skinDesc,2,corresp);
      // validIds is the list of segments which are on both the skin of *this*, and in the segments of the M1 group
      // In the cube example above, this is a U-shaped polyline.
      DAInt validIds = corresp->findIdsInRange(0, meshM2->getNumberOfCells());
      corresp->decrRef();
      if (validIds->getNumberOfTuples())
        {
          // Build the set of segments which are: in the desc mesh of the skin of the 3D mesh (M0) **and** in the desc mesh of the M1 group:
          // (the U-shaped polyline described above)
          MCUMesh m1IntersecSkin = static_cast<MEDCouplingUMesh *>(m0skinDesc->buildPartOfMySelf(validIds->begin(), validIds->end(), true));
          // Its boundary nodes should no be duplicated (this is for example the tip of the crack inside the cube described above)
          DAInt notDuplSkin = m1IntersecSkin->findBoundaryNodes();
          DAInt fNodes1 = fNodes->buildSubstraction(notDuplSkin);

          // Specific logic to handle singular points :
          //   - a point on this U-shape line used in a cell which has no face in common with M1 is deemed singular.
          //   - indeed, if duplicated, such a point would lead to the duplication of a cell which has no face touching M1 ! The
          //   algorithm would be duplicating too much ...
          // This is a costly algorithm so only go into it if a simple (non sufficient) criteria is met: a node connected to more than 3 segs in meshM2:
          dnu1=DataArrayIdType::New(), dnu2=DataArrayIdType::New(), dnu3=DataArrayIdType::New(), rdit0=DataArrayIdType::New();
          MCUMesh meshM2Desc = meshM2->buildDescendingConnectivity(dnu1, dnu2, dnu3, rdit0);  // a mesh made of node cells
          dnu1=0;dnu2=0;dnu3=0;
          dsi = rdit0->deltaShiftIndex();  rdit0=0;
          DAInt singPoints = dsi->findIdsNotInRange(-1,4) ;    dsi=0;// points connected to (strictly) more than 3 segments
          if (singPoints->getNumberOfTuples())
            {
              DAInt boundNodes = m1IntersecSkin->computeFetchedNodeIds();
              // If a point on this U-shape line is connected to cells which do not share any face with M1, then it
              // should not be duplicated
              //    1. Extract N D cells touching U-shape line:
              DAInt cellsAroundBN = getCellIdsLyingOnNodes(boundNodes->begin(), boundNodes->end(), false);  // false= take cell in, even if not all nodes are in dupl
              MCUMesh mAroundBN = static_cast<MEDCouplingUMesh *>(this->buildPartOfMySelf(cellsAroundBN->begin(), cellsAroundBN->end(), true));
              DAInt descBN=DataArrayIdType::New(), descIBN=DataArrayIdType::New(), revDescBN=DataArrayIdType::New(), revDescIBN=DataArrayIdType::New();
              MCUMesh mAroundBNDesc = mAroundBN->buildDescendingConnectivity(descBN,descIBN,revDescBN,revDescIBN);
              //    2. Identify cells in sub-mesh mAroundBN which have a face in common with M1
              DataArrayIdType *idsOfM1BNt;
              mAroundBNDesc->areCellsIncludedIn(otherDimM1OnSameCoords,2, idsOfM1BNt);
              DAInt idsOfM1BN(idsOfM1BNt);
              mcIdType nCells=mAroundBN->getNumberOfCells(), nCellsDesc=mAroundBNDesc->getNumberOfCells();
              DAInt idsTouch=DataArrayIdType::New(); idsTouch->alloc(0,1);
              const mcIdType *revDescIBNP=revDescIBN->begin(), *revDescBNP=revDescBN->begin();
              for(const auto& v: *idsOfM1BN)
                {
                  if (v >= nCellsDesc)    // Keep valid match only
                    continue;
                  mcIdType idx0 = revDescIBNP[v];
                  // Keep the two cells on either side of the face v of M1:
                  mcIdType c1=revDescBNP[idx0], c2=revDescBNP[idx0+1];
                  idsTouch->pushBackSilent(c1);  idsTouch->pushBackSilent(c2);
                }
              //    3. Build complement
              DAInt idsTouchCompl = idsTouch->buildComplement(nCells);
              MCUMesh mAroundBNStrict = static_cast<MEDCouplingUMesh *>(mAroundBN->buildPartOfMySelf(idsTouchCompl->begin(), idsTouchCompl->end(), true));
              DAInt nod3 = mAroundBNStrict->computeFetchedNodeIds();
              DAInt inters = boundNodes->buildIntersection(nod3);
              fNodes1 = fNodes1->buildSubstraction(inters);  // reminder: fNodes1 represent nodes that need dupl.
            }
          notDup = xtrem->buildSubstraction(fNodes1);
        }
      else  // if (validIds-> ...)
        notDup = xtrem->buildSubstraction(fNodes);
    }
  else  // if (3D ...)
    notDup = xtrem->buildSubstraction(fNodes);

  DAInt m1Nodes = otherDimM1OnSameCoords->computeFetchedNodeIds();
  DAInt dupl = m1Nodes->buildSubstraction(notDup);
  return dupl.retn();
}


/*!
 * This method expects that \b this and \b otherDimM1OnSameCoords share the same coordinates array.
 * otherDimM1OnSameCoords->getMeshDimension() is expected to be equal to this->getMeshDimension()-1.
 * This method is part of the MEDFileUMesh::buildInnerBoundaryAlongM1Group() algorithm.
 * Given a set of nodes to duplicate, this method identifies which cells should have their connectivity modified
 * to produce the inner boundary. It is typically called after findNodesToDuplicate().
 *
 * \param [in] otherDimM1OnSameCoords a mesh lying on the same coords than \b this and with a mesh dimension equal to those of \b this minus 1.
 * \param [in] nodeIdsToDuplicateBg node ids needed to be duplicated, as returned by findNodesToDuplicate.
 * \param [in] nodeIdsToDuplicateEnd node ids needed to be duplicated, as returned by findNodesToDuplicate.
 * \param [out] cellIdsNeededToBeRenum cell ids in \b this in which the renumber of nodes should be performed.
 * \param [out] cellIdsNotModified cell ids in \b this that lies on \b otherDimM1OnSameCoords mesh whose connectivity do \b not need to be modified as it is the case for \b cellIdsNeededToBeRenum.
 *
 */
void MEDCouplingUMesh::findCellsToRenumber(const MEDCouplingUMesh& otherDimM1OnSameCoords, const mcIdType *nodeIdsToDuplicateBg, const mcIdType *nodeIdsToDuplicateEnd,
                                           DataArrayIdType *& cellIdsNeededToBeRenum, DataArrayIdType *& cellIdsNotModified) const
{
  using DAInt = MCAuto<DataArrayIdType>;
  using MCUMesh = MCAuto<MEDCouplingUMesh>;

  checkFullyDefined();
  otherDimM1OnSameCoords.checkFullyDefined();
  if(getCoords()!=otherDimM1OnSameCoords.getCoords())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findCellsToRenumber: meshes do not share the same coords array !");
  if(otherDimM1OnSameCoords.getMeshDimension()!=getMeshDimension()-1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findCellsToRenumber: the mesh given in other parameter must have this->getMeshDimension()-1 !");

  // Degenerated case - no nodes to duplicate
  if (nodeIdsToDuplicateBg == nodeIdsToDuplicateEnd)
    {
      cellIdsNeededToBeRenum = DataArrayIdType::New(); cellIdsNeededToBeRenum->alloc(0,1);
      cellIdsNotModified = DataArrayIdType::New(); cellIdsNotModified->alloc(0,1);
      return;
    }

  // Compute cell IDs of the mesh with cells that touch the M1 group with a least one node:
  DAInt cellsAroundGroupLarge = getCellIdsLyingOnNodes(nodeIdsToDuplicateBg, nodeIdsToDuplicateEnd, false);  // false= take cell in, even if not all nodes are in dupl
  MCUMesh mAroundGrpLarge=static_cast<MEDCouplingUMesh *>(buildPartOfMySelf(cellsAroundGroupLarge->begin(),cellsAroundGroupLarge->end(),true));
  mcIdType nCellsLarge=cellsAroundGroupLarge->getNumberOfTuples();
  DAInt descL=DataArrayIdType::New(),descIL=DataArrayIdType::New(),revDescL=DataArrayIdType::New(),revDescIL=DataArrayIdType::New();
  MCUMesh mArGrpLargeDesc=mAroundGrpLarge->buildDescendingConnectivity(descL,descIL,revDescL,revDescIL);
  const mcIdType *descILP=descIL->begin(), *descLP=descL->begin();
  DataArrayIdType *idsOfM1t;
  mArGrpLargeDesc->areCellsIncludedIn(&otherDimM1OnSameCoords,2, idsOfM1t);
  DAInt idsOfM1Large(idsOfM1t);
  mcIdType nL = mArGrpLargeDesc->getNumberOfCells();

  // Computation of the neighbor information of the mesh WITH the crack (some neighbor links are removed):
  //     In the neighbor information remove the connection between high dimension cells and its low level constituents which are part
  //     of the frontier given in parameter (i.e. the cells of low dimension from the group delimiting the crack):
  DAInt descLTrunc = descL->deepCopy(), descILTrunc = descIL->deepCopy();
  DataArrayIdType::RemoveIdsFromIndexedArrays(idsOfM1Large->begin(), idsOfM1Large->end(),descLTrunc,descILTrunc);
  DataArrayIdType *neight=0, *neighIt=0;
  MEDCouplingUMesh::ComputeNeighborsOfCellsAdv(descLTrunc,descILTrunc,revDescL,revDescIL, neight, neighIt);
  DAInt neighL(neight), neighIL(neighIt);

  DAInt hitCellsLarge = DataArrayIdType::New(); hitCellsLarge->alloc(nCellsLarge,1);
  hitCellsLarge->fillWithValue(0);  // 0 : not hit, +1: one side of the crack, -1: other side of the crack,
  mcIdType* hitCellsLargeP = hitCellsLarge->rwBegin();

  // Now loop on the faces of the M1 group and fill spread zones on either side of the crack:
  const mcIdType *revDescILP=revDescIL->begin(), *revDescLP=revDescL->begin();
  for(const auto& v: *idsOfM1Large)
    {
      if (v >= nL) continue;   // Keep valid match only - see doc of areCellsIncludedIn()
      mcIdType idx0 = revDescILP[v];
      // Retrieve the two cells on either side of the face v of M1:
      mcIdType c1=revDescLP[idx0], c2=revDescLP[idx0+1];
      std::map<mcIdType, mcIdType> toOther = {{c1, c2}, {c2, c1}};
      // Handle the spread zones on the two sides of the crack:
      for (const auto c: {c1, c2})
        {
          if (hitCellsLargeP[c]) continue;
          // Identify connex zone around this cell - if we find a value already assigned there, use it.
          mcIdType dnu;
          DAInt spreadZone = MEDCouplingUMesh::ComputeSpreadZoneGraduallyFromSeed(&c, &c+1, neighL,neighIL, -1, dnu);
          std::set<mcIdType> sv;
          for (const mcIdType& s: *spreadZone)
            if (hitCellsLargeP[s]) sv.insert(hitCellsLargeP[s]);
          if (sv.size() > 1)
            // Strange: we find in the same spread zone a +1 and -1 !
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findCellsToRenumber: internal error #0 - conflicting values - should not happen!");
          // If a valid value was found, use it:
          mcIdType val = sv.size()==1 ? *sv.begin() : 0;
          // Hopefully this does not conflict with an potential value on the other side:
          mcIdType other = toOther[c];
          if (hitCellsLargeP[other])
            {
              if(val && hitCellsLargeP[other] != -val)
                throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findCellsToRenumber: internal error #1 - conflicting values - should not happen!");;
              // We do not yet have a value, but other side has one. Use it!
              if(!val) val = -hitCellsLargeP[other];
            }
          // Cover first initialisation:
          if (!val) val = 1;
          // And finally, fill the current spread zone:
          for(const mcIdType& s: *spreadZone) hitCellsLargeP[s] = val;
        }
    }

  DAInt cellsRet1 = hitCellsLarge->findIdsEqual(1);
  DAInt cellsRet2 = hitCellsLarge->findIdsEqual(-1);

  if (cellsRet1->getNumberOfTuples() + cellsRet2->getNumberOfTuples() != cellsAroundGroupLarge->getNumberOfTuples())
    {
      DAInt nonHitCells = hitCellsLarge->findIdsEqual(0); // variable kept for debug ...
      throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findCellsToRenumber: Some cells not hit - Internal error should not happen");
    }
  cellsRet1->transformWithIndArr(cellsAroundGroupLarge->begin(),cellsAroundGroupLarge->end());
  cellsRet2->transformWithIndArr(cellsAroundGroupLarge->begin(),cellsAroundGroupLarge->end());
  //
  cellIdsNeededToBeRenum=cellsRet1.retn();
  cellIdsNotModified=cellsRet2.retn();
}

/*!
 * This method operates a modification of the connectivity and coords in \b this.
 * Every time that a node id in [ \b nodeIdsToDuplicateBg, \b nodeIdsToDuplicateEnd ) will append in nodal connectivity of \b this
 * its ids will be modified to id this->getNumberOfNodes()+std::distance(nodeIdsToDuplicateBg,std::find(nodeIdsToDuplicateBg,nodeIdsToDuplicateEnd,id)).
 * More explicitly the renumber array in nodes is not explicitly given in old2new to avoid to build a big array of renumbering whereas typically few node ids needs to be
 * renumbered. The node id nodeIdsToDuplicateBg[0] will have id this->getNumberOfNodes()+0, node id nodeIdsToDuplicateBg[1] will have id this->getNumberOfNodes()+1,
 * node id nodeIdsToDuplicateBg[2] will have id this->getNumberOfNodes()+2...
 *
 * As a consequence nodal connectivity array length will remain unchanged by this method, and nodal connectivity index array will remain unchanged by this method.
 *
 * \param [in] nodeIdsToDuplicateBg begin of node ids (included) to be duplicated in connectivity only
 * \param [in] nodeIdsToDuplicateEnd end of node ids (excluded) to be duplicated in connectivity only
 */
void MEDCouplingUMesh::duplicateNodes(const mcIdType *nodeIdsToDuplicateBg, const mcIdType *nodeIdsToDuplicateEnd)
{
  mcIdType nbOfNodes=getNumberOfNodes();
  duplicateNodesInCoords(nodeIdsToDuplicateBg,nodeIdsToDuplicateEnd);
  duplicateNodesInConn(nodeIdsToDuplicateBg,nodeIdsToDuplicateEnd,nbOfNodes);
}

/*!
 * This method renumbers only nodal connectivity in \a this. The renumbering is only an offset applied. So this method is a specialization of
 * \a renumberNodesInConn. \b WARNING, this method does not check that the resulting node ids in the nodal connectivity is in a valid range !
 *
 * \param [in] offset - specifies the offset to be applied on each element of connectivity.
 *
 * \sa renumberNodesInConn
 */
void MEDCouplingUMesh::renumberNodesWithOffsetInConn(mcIdType offset)
{
  checkConnectivityFullyDefined();
  mcIdType *conn(getNodalConnectivity()->getPointer());
  const mcIdType *connIndex(getNodalConnectivityIndex()->getConstPointer());
  mcIdType nbOfCells=getNumberOfCells();
  for(mcIdType i=0;i<nbOfCells;i++)
    for(mcIdType iconn=connIndex[i]+1;iconn!=connIndex[i+1];iconn++)
      {
        mcIdType& node=conn[iconn];
        if(node>=0)//avoid polyhedron separator
          {
            node+=offset;
          }
      }
  _nodal_connec->declareAsNew();
  updateTime();
}

/*!
 *  Same than renumberNodesInConn(const mcIdType *) except that here the format of old-to-new traducer is using map instead
 *  of array. This method is dedicated for renumbering from a big set of nodes the a tiny set of nodes which is the case during extraction
 *  of a big mesh.
 */
void MEDCouplingUMesh::renumberNodesInConn(const INTERP_KERNEL::HashMap<mcIdType,mcIdType>& newNodeNumbersO2N)
{
  this->renumberNodesInConnT< INTERP_KERNEL::HashMap<mcIdType,mcIdType> >(newNodeNumbersO2N);
}

/*!
 *  Same than renumberNodesInConn(const mcIdType *) except that here the format of old-to-new traducer is using map instead
 *  of array. This method is dedicated for renumbering from a big set of nodes the a tiny set of nodes which is the case during extraction
 *  of a big mesh.
 */
void MEDCouplingUMesh::renumberNodesInConn(const std::map<mcIdType,mcIdType>& newNodeNumbersO2N)
{
  this->renumberNodesInConnT< std::map<mcIdType,mcIdType> >(newNodeNumbersO2N);
}

/*!
 * Changes ids of nodes within the nodal connectivity arrays according to a permutation
 * array in "Old to New" mode. The node coordinates array is \b not changed by this method.
 * This method is a generalization of shiftNodeNumbersInConn().
 *  \warning This method performs no check of validity of new ids. **Use it with care !**
 *  \param [in] newNodeNumbersO2N - a permutation array, of length \a
 *         this->getNumberOfNodes(), in "Old to New" mode.
 *         See \ref numbering for more info on renumbering modes.
 *  \throw If the nodal connectivity of cells is not defined.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_renumberNodesInConn "Here is a C++ example".<br>
 *  \ref  py_mcumesh_renumberNodesInConn "Here is a Python example".
 *  \endif
 */
void MEDCouplingUMesh::renumberNodesInConn(const mcIdType *newNodeNumbersO2N)
{
  checkConnectivityFullyDefined();
  mcIdType *conn=getNodalConnectivity()->getPointer();
  const mcIdType *connIndex=getNodalConnectivityIndex()->getConstPointer();
  mcIdType nbOfCells=getNumberOfCells();
  for(mcIdType i=0;i<nbOfCells;i++)
    for(mcIdType iconn=connIndex[i]+1;iconn!=connIndex[i+1];iconn++)
      {
        mcIdType& node=conn[iconn];
        if(node>=0)//avoid polyhedron separator
          {
            node=newNodeNumbersO2N[node];
          }
      }
  _nodal_connec->declareAsNew();
  updateTime();
}

/*!
 * This method renumbers nodes \b in \b connectivity \b only \b without \b any \b reference \b to \b coords.
 * This method performs no check on the fact that new coordinate ids are valid. \b Use \b it \b with \b care !
 * This method is an specialization of \ref MEDCoupling::MEDCouplingUMesh::renumberNodesInConn "renumberNodesInConn method".
 *
 * \param [in] delta specifies the shift size applied to nodeId in nodal connectivity in \b this.
 */
void MEDCouplingUMesh::shiftNodeNumbersInConn(mcIdType delta)
{
  checkConnectivityFullyDefined();
  mcIdType *conn=getNodalConnectivity()->getPointer();
  const mcIdType *connIndex=getNodalConnectivityIndex()->getConstPointer();
  mcIdType nbOfCells=getNumberOfCells();
  for(mcIdType i=0;i<nbOfCells;i++)
    for(mcIdType iconn=connIndex[i]+1;iconn!=connIndex[i+1];iconn++)
      {
        mcIdType& node=conn[iconn];
        if(node>=0)//avoid polyhedron separator
          {
            node+=delta;
          }
      }
  _nodal_connec->declareAsNew();
  updateTime();
}

/*!
 * This method operates a modification of the connectivity in \b this.
 * Coordinates are \b NOT considered here and will remain unchanged by this method. this->_coords can ever been null for the needs of this method.
 * Every time that a node id in [ \b nodeIdsToDuplicateBg, \b nodeIdsToDuplicateEnd ) will append in nodal connectivity of \b this
 * its ids will be modified to id offset+std::distance(nodeIdsToDuplicateBg,std::find(nodeIdsToDuplicateBg,nodeIdsToDuplicateEnd,id)).
 * More explicitly the renumber array in nodes is not explicitly given in old2new to avoid to build a big array of renumbering whereas typically few node ids needs to be
 * renumbered. The node id nodeIdsToDuplicateBg[0] will have id offset+0, node id nodeIdsToDuplicateBg[1] will have id offset+1,
 * node id nodeIdsToDuplicateBg[2] will have id offset+2...
 *
 * As a consequence nodal connectivity array length will remain unchanged by this method, and nodal connectivity index array will remain unchanged by this method.
 * As an another consequense after the call of this method \b this can be transiently non cohrent.
 *
 * \param [in] nodeIdsToDuplicateBg begin of node ids (included) to be duplicated in connectivity only
 * \param [in] nodeIdsToDuplicateEnd end of node ids (excluded) to be duplicated in connectivity only
 * \param [in] offset the offset applied to all node ids in connectivity that are in [ \a nodeIdsToDuplicateBg, \a nodeIdsToDuplicateEnd ).
 */
void MEDCouplingUMesh::duplicateNodesInConn(const mcIdType *nodeIdsToDuplicateBg, const mcIdType *nodeIdsToDuplicateEnd, mcIdType offset)
{
  checkConnectivityFullyDefined();
  std::map<mcIdType,mcIdType> m;
  mcIdType val=offset;
  for(const mcIdType *work=nodeIdsToDuplicateBg;work!=nodeIdsToDuplicateEnd;work++,val++)
    m[*work]=val;
  mcIdType *conn=getNodalConnectivity()->getPointer();
  const mcIdType *connIndex=getNodalConnectivityIndex()->getConstPointer();
  mcIdType nbOfCells=getNumberOfCells();
  for(mcIdType i=0;i<nbOfCells;i++)
    for(mcIdType iconn=connIndex[i]+1;iconn!=connIndex[i+1];iconn++)
      {
        mcIdType& node=conn[iconn];
        if(node>=0)//avoid polyhedron separator
          {
            std::map<mcIdType,mcIdType>::iterator it=m.find(node);
            if(it!=m.end())
              node=(*it).second;
          }
      }
  updateTime();
}

/*!
 * This method renumbers cells of \a this using the array specified by [old2NewBg;old2NewBg+getNumberOfCells())
 *
 * Contrary to MEDCouplingPointSet::renumberNodes, this method makes a permutation without any fuse of cell.
 * After the call of this method the number of cells remains the same as before.
 *
 * If 'check' equals true the method will check that any elements in [ \a old2NewBg; \a old2NewEnd ) is unique ; if not
 * an INTERP_KERNEL::Exception will be thrown. When 'check' equals true [ \a old2NewBg ; \a old2NewEnd ) is not expected to
 * be strictly in [0;this->getNumberOfCells()).
 *
 * If 'check' equals false the method will not check the content of [ \a old2NewBg ; \a old2NewEnd ).
 * To avoid any throw of SIGSEGV when 'check' equals false, the elements in [ \a old2NewBg ; \a old2NewEnd ) should be unique and
 * should be contained in[0;this->getNumberOfCells()).
 *
 * \param [in] old2NewBg is expected to be a dynamically allocated pointer of size at least equal to this->getNumberOfCells()
 * \param check whether to check content of old2NewBg
 */
void MEDCouplingUMesh::renumberCells(const mcIdType *old2NewBg, bool check)
{
  checkConnectivityFullyDefined();
  mcIdType nbCells=getNumberOfCells();
  const mcIdType *array=old2NewBg;
  if(check)
    array=DataArrayIdType::CheckAndPreparePermutation(old2NewBg,old2NewBg+nbCells);
  //
  const mcIdType *conn=_nodal_connec->getConstPointer();
  const mcIdType *connI=_nodal_connec_index->getConstPointer();
  MCAuto<DataArrayIdType> o2n=DataArrayIdType::New(); o2n->useArray(array,false,DeallocType::C_DEALLOC,nbCells,1);
  MCAuto<DataArrayIdType> n2o=o2n->invertArrayO2N2N2O(nbCells);
  const mcIdType *n2oPtr=n2o->begin();
  MCAuto<DataArrayIdType> newConn=DataArrayIdType::New();
  newConn->alloc(_nodal_connec->getNumberOfTuples(),_nodal_connec->getNumberOfComponents());
  newConn->copyStringInfoFrom(*_nodal_connec);
  MCAuto<DataArrayIdType> newConnI=DataArrayIdType::New();
  newConnI->alloc(_nodal_connec_index->getNumberOfTuples(),_nodal_connec_index->getNumberOfComponents());
  newConnI->copyStringInfoFrom(*_nodal_connec_index);
  //
  mcIdType *newC=newConn->getPointer();
  mcIdType *newCI=newConnI->getPointer();
  mcIdType loc=0;
  newCI[0]=loc;
  for(mcIdType i=0;i<nbCells;i++)
    {
      mcIdType pos=n2oPtr[i];
      mcIdType nbOfElts=connI[pos+1]-connI[pos];
      newC=std::copy(conn+connI[pos],conn+connI[pos+1],newC);
      loc+=nbOfElts;
      newCI[i+1]=loc;
    }
  //
  setConnectivity(newConn,newConnI);
  if(check)
    free(const_cast<mcIdType *>(array));
}

/*!
 * Finds cells whose bounding boxes intersect a given bounding box.
 *  \param [in] bbox - an array defining the bounding box via coordinates of its
 *         extremum points in "no interlace" mode, i.e. xMin, xMax, yMin, yMax, zMin,
 *         zMax (if in 3D).
 *  \param [in] eps - a factor used to increase size of the bounding box of cell
 *         before comparing it with \a bbox. This factor is multiplied by the maximal
 *         extent of the bounding box of cell to produce an addition to this bounding box.
 *  \return DataArrayIdType * - a new instance of DataArrayIdType holding ids for found
 *         cells. The caller is to delete this array using decrRef() as it is no more
 *         needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_getCellsInBoundingBox "Here is a C++ example".<br>
 *  \ref  py_mcumesh_getCellsInBoundingBox "Here is a Python example".
 *  \endif
 */
DataArrayIdType *MEDCouplingUMesh::getCellsInBoundingBox(const double *bbox, double eps) const
{
  MCAuto<DataArrayIdType> elems=DataArrayIdType::New(); elems->alloc(0,1);
  if(getMeshDimension()==-1)
    {
      elems->pushBackSilent(0);
      return elems.retn();
    }
  int dim=getSpaceDimension();
  INTERP_KERNEL::AutoPtr<double> elem_bb=new double[2*dim];
  const mcIdType* conn      = getNodalConnectivity()->getConstPointer();
  const mcIdType* conn_index= getNodalConnectivityIndex()->getConstPointer();
  const double* coords = getCoords()->getConstPointer();
  mcIdType nbOfCells=getNumberOfCells();
  for ( mcIdType ielem=0; ielem<nbOfCells;ielem++ )
    {
      for (int i=0; i<dim; i++)
        {
          elem_bb[i*2]=std::numeric_limits<double>::max();
          elem_bb[i*2+1]=-std::numeric_limits<double>::max();
        }

      for (mcIdType inode=conn_index[ielem]+1; inode<conn_index[ielem+1]; inode++)//+1 due to offset of cell type.
        {
          mcIdType node= conn[inode];
          if(node>=0)//avoid polyhedron separator
            {
              for (int idim=0; idim<dim; idim++)
                {
                  if ( coords[node*dim+idim] < elem_bb[idim*2] )
                    {
                      elem_bb[idim*2] = coords[node*dim+idim] ;
                    }
                  if ( coords[node*dim+idim] > elem_bb[idim*2+1] )
                    {
                      elem_bb[idim*2+1] = coords[node*dim+idim] ;
                    }
                }
            }
        }
      if (intersectsBoundingBox(elem_bb, bbox, dim, eps))
        elems->pushBackSilent(ielem);
    }
  return elems.retn();
}

/*!
 * Given a boundary box 'bbox' returns elements 'elems' contained in this 'bbox' or touching 'bbox' (within 'eps' distance).
 * Warning 'elems' is incremented during the call so if elems is not empty before call returned elements will be
 * added in 'elems' parameter.
 */
DataArrayIdType *MEDCouplingUMesh::getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps)
{
  MCAuto<DataArrayIdType> elems=DataArrayIdType::New(); elems->alloc(0,1);
  if(getMeshDimension()==-1)
    {
      elems->pushBackSilent(0);
      return elems.retn();
    }
  int dim=getSpaceDimension();
  INTERP_KERNEL::AutoPtr<double> elem_bb=new double[2*dim];
  const mcIdType* conn      = getNodalConnectivity()->getConstPointer();
  const mcIdType* conn_index= getNodalConnectivityIndex()->getConstPointer();
  const double* coords = getCoords()->getConstPointer();
  mcIdType nbOfCells=getNumberOfCells();
  for ( mcIdType ielem=0; ielem<nbOfCells;ielem++ )
    {
      for (int i=0; i<dim; i++)
        {
          elem_bb[i*2]=std::numeric_limits<double>::max();
          elem_bb[i*2+1]=-std::numeric_limits<double>::max();
        }

      for (mcIdType inode=conn_index[ielem]+1; inode<conn_index[ielem+1]; inode++)//+1 due to offset of cell type.
        {
          mcIdType node= conn[inode];
          if(node>=0)//avoid polyhedron separator
            {
              for (int idim=0; idim<dim; idim++)
                {
                  if ( coords[node*dim+idim] < elem_bb[idim*2] )
                    {
                      elem_bb[idim*2] = coords[node*dim+idim] ;
                    }
                  if ( coords[node*dim+idim] > elem_bb[idim*2+1] )
                    {
                      elem_bb[idim*2+1] = coords[node*dim+idim] ;
                    }
                }
            }
        }
      if(intersectsBoundingBox(bbox, elem_bb, dim, eps))
        elems->pushBackSilent(ielem);
    }
  return elems.retn();
}

/*!
 * Returns a type of a cell by its id.
 *  \param [in] cellId - the id of the cell of interest.
 *  \return INTERP_KERNEL::NormalizedCellType - enumeration item describing the cell type.
 *  \throw If \a cellId is invalid. Valid range is [0, \a this->getNumberOfCells() ).
 */
INTERP_KERNEL::NormalizedCellType MEDCouplingUMesh::getTypeOfCell(mcIdType cellId) const
{
  const mcIdType *ptI(_nodal_connec_index->begin()),*pt(_nodal_connec->begin());
  if(cellId<_nodal_connec_index->getNbOfElems()-1)
    return (INTERP_KERNEL::NormalizedCellType) pt[ptI[cellId]];
  else
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::getTypeOfCell : Requesting type of cell #" << cellId << " but it should be in [0," << _nodal_connec_index->getNbOfElems()-1 << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
}

/*!
 * This method returns a newly allocated array containing cell ids (ascendingly sorted) whose geometric type are equal to type.
 * This method does not throw exception if geometric type \a type is not in \a this.
 * This method throws an INTERP_KERNEL::Exception if meshdimension of \b this is not equal to those of \b type.
 * The coordinates array is not considered here.
 *
 * \param [in] type the geometric type
 * \return cell ids in this having geometric type \a type.
 */
DataArrayIdType *MEDCouplingUMesh::giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const
{

  MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
  ret->alloc(0,1);
  checkConnectivityFullyDefined();
  mcIdType nbCells=getNumberOfCells();
  int mdim=getMeshDimension();
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  if(mdim!=ToIdType(cm.getDimension()))
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::giveCellsWithType : Mismatch between mesh dimension and dimension of the cell !");
  const mcIdType *ptI=_nodal_connec_index->getConstPointer();
  const mcIdType *pt=_nodal_connec->getConstPointer();
  for(mcIdType i=0;i<nbCells;i++)
    {
      if((INTERP_KERNEL::NormalizedCellType)pt[ptI[i]]==type)
        ret->pushBackSilent(i);
    }
  return ret.retn();
}

/*!
 * Returns nb of cells having the geometric type \a type. No throw if no cells in \a this has the geometric type \a type.
 */
mcIdType MEDCouplingUMesh::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const
{
  const mcIdType *ptI(_nodal_connec_index->begin()),*pt(_nodal_connec->begin());
  mcIdType nbOfCells(getNumberOfCells()),ret(0);
  for(mcIdType i=0;i<nbOfCells;i++)
    if((INTERP_KERNEL::NormalizedCellType) pt[ptI[i]]==type)
      ret++;
  return ret;
}

/*!
 * Returns the nodal connectivity of a given cell.
 * The separator of faces within polyhedron connectivity (-1) is not returned, thus
 * all returned node ids can be used in getCoordinatesOfNode().
 *  \param [in] cellId - an id of the cell of interest.
 *  \param [in,out] conn - a vector where the node ids are appended. It is not
 *         cleared before the appending.
 *  \throw If \a cellId is invalid. Valid range is [0, \a this->getNumberOfCells() ).
 */
void MEDCouplingUMesh::getNodeIdsOfCell(mcIdType cellId, std::vector<mcIdType>& conn) const
{
  const mcIdType *ptI(_nodal_connec_index->begin()),*pt(_nodal_connec->begin());
  for(const mcIdType *w=pt+ptI[cellId]+1;w!=pt+ptI[cellId+1];w++)
    if(*w>=0)
      conn.push_back(*w);
}

std::string MEDCouplingUMesh::simpleRepr() const
{
  static const char msg0[]="No coordinates specified !";
  std::ostringstream ret;
  ret << "Unstructured mesh with name : \"" << getName() << "\"\n";
  ret << "Description of mesh : \"" << getDescription() << "\"\n";
  int tmpp1,tmpp2;
  double tt=getTime(tmpp1,tmpp2);
  ret << "Time attached to the mesh [unit] : " << tt << " [" << getTimeUnit() << "]\n";
  ret << "Iteration : " << tmpp1  << " Order : " << tmpp2 << "\n";
  if(_mesh_dim>=-1)
    { ret << "Mesh dimension : " << _mesh_dim << "\nSpace dimension : "; }
  else
    { ret << " Mesh dimension has not been set or is invalid !"; }
  if(_coords!=0)
    {
      const int spaceDim=getSpaceDimension();
      ret << spaceDim << "\nInfo attached on space dimension : ";
      for(int i=0;i<spaceDim;i++)
        ret << "\"" << _coords->getInfoOnComponent(i) << "\" ";
      ret << "\n";
    }
  else
    ret << msg0 << "\n";
  ret << "Number of nodes : ";
  if(_coords!=0)
    ret << getNumberOfNodes() << "\n";
  else
    ret << msg0 << "\n";
  ret << "Number of cells : ";
  if(_nodal_connec!=0 && _nodal_connec_index!=0)
    ret << getNumberOfCells() << "\n";
  else
    ret << "No connectivity specified !" << "\n";
  ret << "Cell types present : ";
  for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iter=_types.begin();iter!=_types.end();iter++)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(*iter);
      ret << cm.getRepr() << " ";
    }
  ret << "\n";
  return ret.str();
}

std::string MEDCouplingUMesh::advancedRepr() const
{
  std::ostringstream ret;
  ret << simpleRepr();
  ret << "\nCoordinates array : \n___________________\n\n";
  if(_coords)
    _coords->reprWithoutNameStream(ret);
  else
    ret << "No array set !\n";
  ret << "\n\nConnectivity arrays : \n_____________________\n\n";
  reprConnectivityOfThisLL(ret);
  return ret.str();
}

/*!
 * This method returns a C++ code that is a dump of \a this.
 * This method will throw if this is not fully defined.
 */
std::string MEDCouplingUMesh::cppRepr() const
{
  static const char coordsName[]="coords";
  static const char connName[]="conn";
  static const char connIName[]="connI";
  checkFullyDefined();
  std::ostringstream ret; ret << "// coordinates" << std::endl;
  _coords->reprCppStream(coordsName,ret); ret << std::endl << "// connectivity" << std::endl;
  _nodal_connec->reprCppStream(connName,ret); ret << std::endl;
  _nodal_connec_index->reprCppStream(connIName,ret); ret << std::endl;
  ret << "MEDCouplingUMesh *mesh=MEDCouplingUMesh::New(\"" << getName() << "\"," << getMeshDimension() << ");" << std::endl;
  ret << "mesh->setCoords(" << coordsName << ");" << std::endl;
  ret << "mesh->setConnectivity(" << connName << "," << connIName << ",true);" << std::endl;
  ret << coordsName << "->decrRef(); " << connName << "->decrRef(); " << connIName << "->decrRef();" << std::endl;
  return ret.str();
}

std::string MEDCouplingUMesh::reprConnectivityOfThis() const
{
  std::ostringstream ret;
  reprConnectivityOfThisLL(ret);
  return ret.str();
}

/*!
 * This method builds a newly allocated instance (with the same name than \a this) that the caller has the responsibility to deal with.
 * This method returns an instance with all arrays allocated (connectivity, connectivity index, coordinates)
 * but with length of these arrays set to 0. It allows to define an "empty" mesh (with nor cells nor nodes but compliant with
 * some algos).
 *
 * This method expects that \a this has a mesh dimension set and higher or equal to 0. If not an exception will be thrown.
 * This method analyzes the 3 arrays of \a this. For each the following behaviour is done : if the array is null a newly one is created
 * with number of tuples set to 0, if not the array is taken as this in the returned instance.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildSetInstanceFromThis(std::size_t spaceDim) const
{
  int mdim=getMeshDimension();
  if(mdim<0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSetInstanceFromThis : invalid mesh dimension ! Should be >= 0 !");
  MCAuto<MEDCouplingUMesh> ret=MEDCouplingUMesh::New(getName(),mdim);
  MCAuto<DataArrayIdType> tmp1,tmp2;
  bool needToCpyCT=true;
  if(!_nodal_connec)
    {
      tmp1=DataArrayIdType::New(); tmp1->alloc(0,1);
      needToCpyCT=false;
    }
  else
    {
      tmp1=_nodal_connec;
      tmp1->incrRef();
    }
  if(!_nodal_connec_index)
    {
      tmp2=DataArrayIdType::New(); tmp2->alloc(1,1); tmp2->setIJ(0,0,0);
      needToCpyCT=false;
    }
  else
    {
      tmp2=_nodal_connec_index;
      tmp2->incrRef();
    }
  ret->setConnectivity(tmp1,tmp2,false);
  if(needToCpyCT)
    ret->_types=_types;
  if(!_coords)
    {
      MCAuto<DataArrayDouble> coords=DataArrayDouble::New(); coords->alloc(0,spaceDim);
      ret->setCoords(coords);
    }
  else
    ret->setCoords(_coords);
  return ret.retn();
}

mcIdType MEDCouplingUMesh::getNumberOfNodesInCell(mcIdType cellId) const
{
  const mcIdType *ptI=_nodal_connec_index->getConstPointer();
  const mcIdType *pt=_nodal_connec->getConstPointer();
  if(pt[ptI[cellId]]!=INTERP_KERNEL::NORM_POLYHED)
    return ptI[cellId+1]-ptI[cellId]-1;
  else
    return ToIdType(std::count_if(pt+ptI[cellId]+1,pt+ptI[cellId+1],std::bind(std::not_equal_to<mcIdType>(),std::placeholders::_1,-1)));
}

/*!
 * Returns types of cells of the specified part of \a this mesh.
 * This method avoids computing sub-mesh explicitly to get its types.
 *  \param [in] begin - an array of cell ids of interest.
 *  \param [in] end - the end of \a begin, i.e. a pointer to its (last+1)-th element.
 *  \return std::set<INTERP_KERNEL::NormalizedCellType> - a set of enumeration items
 *         describing the cell types.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \sa getAllGeoTypes()
 */
std::set<INTERP_KERNEL::NormalizedCellType> MEDCouplingUMesh::getTypesOfPart(const mcIdType *begin, const mcIdType *end) const
{
  checkFullyDefined();
  std::set<INTERP_KERNEL::NormalizedCellType> ret;
  const mcIdType *conn=_nodal_connec->getConstPointer();
  const mcIdType *connIndex=_nodal_connec_index->getConstPointer();
  for(const mcIdType *w=begin;w!=end;w++)
    ret.insert((INTERP_KERNEL::NormalizedCellType)conn[connIndex[*w]]);
  return ret;
}

/*!
 * Defines the nodal connectivity using given connectivity arrays in \ref numbering-indirect format.
 * Optionally updates
 * a set of types of cells constituting \a this mesh.
 * This method is for advanced users having prepared their connectivity before. For
 * more info on using this method see \ref MEDCouplingUMeshAdvBuild.
 *  \param [in] conn - the nodal connectivity array.
 *  \param [in] connIndex - the nodal connectivity index array.
 *  \param [in] isComputingTypes - if \c true, the set of types constituting \a this
 *         mesh is updated.
 */
void MEDCouplingUMesh::setConnectivity(DataArrayIdType *conn, DataArrayIdType *connIndex, bool isComputingTypes)
{
  DataArrayIdType::SetArrayIn(conn,_nodal_connec);
  DataArrayIdType::SetArrayIn(connIndex,_nodal_connec_index);
  if(isComputingTypes)
    computeTypes();
  declareAsNew();
}

/*!
 * Copy constructor. If 'deepCopy' is false \a this is a shallow copy of other.
 * If 'deeCpy' is true all arrays (coordinates and connectivities) are deeply copied.
 */
MEDCouplingUMesh::MEDCouplingUMesh(const MEDCouplingUMesh& other, bool deepCpy):MEDCouplingPointSet(other,deepCpy),_mesh_dim(other._mesh_dim),
    _nodal_connec(0),_nodal_connec_index(0),
    _types(other._types)
{
  if(other._nodal_connec)
    _nodal_connec=other._nodal_connec->performCopyOrIncrRef(deepCpy);
  if(other._nodal_connec_index)
    _nodal_connec_index=other._nodal_connec_index->performCopyOrIncrRef(deepCpy);
}

MEDCouplingUMesh::~MEDCouplingUMesh()
{
  if(_nodal_connec)
    _nodal_connec->decrRef();
  if(_nodal_connec_index)
    _nodal_connec_index->decrRef();
}

/*!
 * Recomputes a set of cell types of \a this mesh. For more info see
 * \ref MEDCouplingUMeshNodalConnectivity.
 */
void MEDCouplingUMesh::computeTypes()
{
  ComputeAllTypesInternal(_types,_nodal_connec,_nodal_connec_index);
}


/*!
 * Returns a number of cells constituting \a this mesh.
 *  \return mcIdType - the number of cells in \a this mesh.
 *  \throw If the nodal connectivity of cells is not defined.
 */
mcIdType MEDCouplingUMesh::getNumberOfCells() const
{
  if(_nodal_connec_index)
    return _nodal_connec_index->getNumberOfTuples()-1;
  else
    if(_mesh_dim==-1)
      return 1;
    else
      throw INTERP_KERNEL::Exception("Unable to get number of cells because no connectivity specified !");
}

/*!
 * Returns a dimension of \a this mesh, i.e. a dimension of cells constituting \a this
 * mesh. For more info see \ref meshes.
 *  \return int - the dimension of \a this mesh.
 *  \throw If the mesh dimension is not defined using setMeshDimension().
 */
int MEDCouplingUMesh::getMeshDimension() const
{
  if(_mesh_dim<-1)
    throw INTERP_KERNEL::Exception("No mesh dimension specified !");
  return _mesh_dim;
}

/*!
 * Returns a length of the nodal connectivity array.
 * This method is for test reason. Normally the integer returned is not useable by
 * user.  For more info see \ref MEDCouplingUMeshNodalConnectivity.
 *  \return mcIdType - the length of the nodal connectivity array.
 */
mcIdType MEDCouplingUMesh::getNodalConnectivityArrayLen() const
{
  return _nodal_connec->getNbOfElems();
}

/*!
 * First step of serialization process. Used by ParaMEDMEM and MEDCouplingCorba to transfert data between process.
 */
void MEDCouplingUMesh::getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<mcIdType>& tinyInfo, std::vector<std::string>& littleStrings) const
{
  MEDCouplingPointSet::getTinySerializationInformation(tinyInfoD,tinyInfo,littleStrings);
  tinyInfo.push_back(ToIdType(getMeshDimension()));
  tinyInfo.push_back(getNumberOfCells());
  if(_nodal_connec)
    tinyInfo.push_back(getNodalConnectivityArrayLen());
  else
    tinyInfo.push_back(-1);
}

/*!
 * First step of unserialization process.
 */
bool MEDCouplingUMesh::isEmptyMesh(const std::vector<mcIdType>& tinyInfo) const
{
  return tinyInfo[6]<=0;
}

/*!
 * Second step of serialization process.
 * \param tinyInfo must be equal to the result given by getTinySerializationInformation method.
 * \param a1 DataArrayDouble
 * \param a2 DataArrayDouble
 * \param littleStrings string vector
 */
void MEDCouplingUMesh::resizeForUnserialization(const std::vector<mcIdType>& tinyInfo, DataArrayIdType *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const
{
  MEDCouplingPointSet::resizeForUnserialization(tinyInfo,a1,a2,littleStrings);
  if(tinyInfo[5]!=-1)
    a1->alloc(tinyInfo[7]+tinyInfo[6]+1,1);
}

/*!
 * Third and final step of serialization process.
 */
void MEDCouplingUMesh::serialize(DataArrayIdType *&a1, DataArrayDouble *&a2) const
{
  MEDCouplingPointSet::serialize(a1,a2);
  if(getMeshDimension()>-1)
    {
      a1=DataArrayIdType::New();
      a1->alloc(getNodalConnectivityArrayLen()+getNumberOfCells()+1,1);
      mcIdType *ptA1=a1->getPointer();
      const mcIdType *conn=getNodalConnectivity()->getConstPointer();
      const mcIdType *index=getNodalConnectivityIndex()->getConstPointer();
      ptA1=std::copy(index,index+getNumberOfCells()+1,ptA1);
      std::copy(conn,conn+getNodalConnectivityArrayLen(),ptA1);
    }
  else
    a1=0;
}

/*!
 * Second and final unserialization process.
 * \param tinyInfo must be equal to the result given by getTinySerializationInformation method.
 */
void MEDCouplingUMesh::unserialization(const std::vector<double>& tinyInfoD, const std::vector<mcIdType>& tinyInfo, const DataArrayIdType *a1, DataArrayDouble *a2, const std::vector<std::string>& littleStrings)
{
  MEDCouplingPointSet::unserialization(tinyInfoD,tinyInfo,a1,a2,littleStrings);
  setMeshDimension(FromIdType<int>(tinyInfo[5]));
  if(tinyInfo[7]!=-1)
    {
      // Connectivity
      const mcIdType *recvBuffer=a1->getConstPointer();
      MCAuto<DataArrayIdType> myConnecIndex=DataArrayIdType::New();
      myConnecIndex->alloc(tinyInfo[6]+1,1);
      std::copy(recvBuffer,recvBuffer+tinyInfo[6]+1,myConnecIndex->getPointer());
      MCAuto<DataArrayIdType> myConnec=DataArrayIdType::New();
      myConnec->alloc(tinyInfo[7],1);
      std::copy(recvBuffer+tinyInfo[6]+1,recvBuffer+tinyInfo[6]+1+tinyInfo[7],myConnec->getPointer());
      setConnectivity(myConnec, myConnecIndex);
    }
}



/*!
 * Returns a new MEDCouplingFieldDouble containing volumes of cells constituting \a this
 * mesh.<br>
 * For 1D cells, the returned field contains lengths.<br>
 * For 2D cells, the returned field contains areas.<br>
 * For 3D cells, the returned field contains volumes.
 *  \param [in] isAbs - if \c true, the computed cell volume does not reflect cell
 *         orientation, i.e. the volume is always positive.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble on cells
 *         and one time . The caller is to delete this field using decrRef() as it is no
 *         more needed.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::getMeasureField(bool isAbs) const
{
  std::string name="MeasureOfMesh_";
  name+=getName();
  mcIdType nbelem=getNumberOfCells();
  MCAuto<MEDCouplingFieldDouble> field=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  field->setName(name);
  MCAuto<DataArrayDouble> array=DataArrayDouble::New();
  array->alloc(nbelem,1);
  double *area_vol=array->getPointer();
  field->setArray(array) ; array=0;
  field->setMesh(const_cast<MEDCouplingUMesh *>(this));
  field->synchronizeTimeWithMesh();
  if(getMeshDimension()!=-1)
    {
      mcIdType ipt;
      INTERP_KERNEL::NormalizedCellType type;
      int dim_space=getSpaceDimension();
      const double *coords=getCoords()->getConstPointer();
      const mcIdType *connec=getNodalConnectivity()->getConstPointer();
      const mcIdType *connec_index=getNodalConnectivityIndex()->getConstPointer();
      for(mcIdType iel=0;iel<nbelem;iel++)
        {
          ipt=connec_index[iel];
          type=(INTERP_KERNEL::NormalizedCellType)connec[ipt];
          area_vol[iel]=INTERP_KERNEL::computeVolSurfOfCell2<mcIdType,INTERP_KERNEL::ALL_C_MODE>(type,connec+ipt+1,connec_index[iel+1]-ipt-1,coords,dim_space);
        }
      if(isAbs)
        std::transform(area_vol,area_vol+nbelem,area_vol,[](double c){return fabs(c);});
    }
  else
    {
      area_vol[0]=std::numeric_limits<double>::max();
    }
  return field.retn();
}

/*!
 * Returns a new DataArrayDouble containing volumes of specified cells of \a this
 * mesh.<br>
 * For 1D cells, the returned array contains lengths.<br>
 * For 2D cells, the returned array contains areas.<br>
 * For 3D cells, the returned array contains volumes.
 * This method avoids building explicitly a part of \a this mesh to perform the work.
 *  \param [in] isAbs - if \c true, the computed cell volume does not reflect cell
 *         orientation, i.e. the volume is always positive.
 *  \param [in] begin - an array of cell ids of interest.
 *  \param [in] end - the end of \a begin, i.e. a pointer to its (last+1)-th element.
 *  \return DataArrayDouble * - a new instance of DataArrayDouble. The caller is to
 *          delete this array using decrRef() as it is no more needed.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_getPartMeasureField "Here is a C++ example".<br>
 *  \ref  py_mcumesh_getPartMeasureField "Here is a Python example".
 *  \endif
 *  \sa getMeasureField()
 */
DataArrayDouble *MEDCouplingUMesh::getPartMeasureField(bool isAbs, const mcIdType *begin, const mcIdType *end) const
{
  std::string name="PartMeasureOfMesh_";
  name+=getName();
  std::size_t nbelem=std::distance(begin,end);
  MCAuto<DataArrayDouble> array=DataArrayDouble::New();
  array->setName(name);
  array->alloc(nbelem,1);
  double *area_vol=array->getPointer();
  if(getMeshDimension()!=-1)
    {
      mcIdType ipt;
      INTERP_KERNEL::NormalizedCellType type;
      int dim_space=getSpaceDimension();
      const double *coords=getCoords()->getConstPointer();
      const mcIdType *connec=getNodalConnectivity()->getConstPointer();
      const mcIdType *connec_index=getNodalConnectivityIndex()->getConstPointer();
      for(const mcIdType *iel=begin;iel!=end;iel++)
        {
          ipt=connec_index[*iel];
          type=(INTERP_KERNEL::NormalizedCellType)connec[ipt];
          *area_vol++=INTERP_KERNEL::computeVolSurfOfCell2<mcIdType,INTERP_KERNEL::ALL_C_MODE>(type,connec+ipt+1,connec_index[*iel+1]-ipt-1,coords,dim_space);
        }
      if(isAbs)
        std::transform(array->getPointer(),area_vol,array->getPointer(),[](double c){return fabs(c);});
    }
  else
    {
      area_vol[0]=std::numeric_limits<double>::max();
    }
  return array.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble containing volumes of cells of a dual mesh of
 * \a this one. The returned field contains the dual cell volume for each corresponding
 * node in \a this mesh. In other words, the field returns the getMeasureField() of
 *  the dual mesh in P1 sens of \a this.<br>
 * For 1D cells, the returned field contains lengths.<br>
 * For 2D cells, the returned field contains areas.<br>
 * For 3D cells, the returned field contains volumes.
 * This method is useful to check "P1*" conservative interpolators.
 *  \param [in] isAbs - if \c true, the computed cell volume does not reflect cell
 *         orientation, i.e. the volume is always positive.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble on
 *          nodes and one time. The caller is to delete this array using decrRef() as
 *          it is no more needed.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::getMeasureFieldOnNode(bool isAbs) const
{
  MCAuto<MEDCouplingFieldDouble> tmp=getMeasureField(isAbs);
  std::string name="MeasureOnNodeOfMesh_";
  name+=getName();
  mcIdType nbNodes=getNumberOfNodes();
  MCAuto<DataArrayDouble> nnpc;
  {
    MCAuto<DataArrayIdType> tmp2(computeNbOfNodesPerCell());
    nnpc=tmp2->convertToDblArr();
  }
  std::for_each(nnpc->rwBegin(),nnpc->rwEnd(),[](double& v) { v=1./v; });
  const double *nnpcPtr(nnpc->begin());
  MCAuto<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(ON_NODES);
  MCAuto<DataArrayDouble> array=DataArrayDouble::New();
  array->alloc(nbNodes,1);
  double *valsToFill=array->getPointer();
  std::fill(valsToFill,valsToFill+nbNodes,0.);
  const double *values=tmp->getArray()->getConstPointer();
  MCAuto<DataArrayIdType> da=DataArrayIdType::New();
  MCAuto<DataArrayIdType> daInd=DataArrayIdType::New();
  getReverseNodalConnectivity(da,daInd);
  const mcIdType *daPtr=da->getConstPointer();
  const mcIdType *daIPtr=daInd->getConstPointer();
  for(mcIdType i=0;i<nbNodes;i++)
    for(const mcIdType *cell=daPtr+daIPtr[i];cell!=daPtr+daIPtr[i+1];cell++)
      valsToFill[i]+=nnpcPtr[*cell]*values[*cell];
  ret->setMesh(this);
  ret->setArray(array);
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble holding normal vectors to cells of \a this
 * mesh. The returned normal vectors to each cell have a norm2 equal to 1.
 * The computed vectors have <em> this->getMeshDimension()+1 </em> components
 * and are normalized.
 * <br> \a this can be either
 * - a  2D mesh in 2D or 3D space or
 * - an 1D mesh in 2D space.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble on
 *          cells and one time. The caller is to delete this field using decrRef() as
 *          it is no more needed.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If the coordinates array is not set.
 *  \throw If the mesh dimension is not set.
 *  \throw If the mesh and space dimension is not as specified above.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::buildOrthogonalField() const
{
  if((getMeshDimension()!=2) && (getMeshDimension()!=1 || getSpaceDimension()!=2))
    throw INTERP_KERNEL::Exception("Expected a umesh with ( meshDim == 2 spaceDim == 2 or 3 ) or ( meshDim == 1 spaceDim == 2 ) !");
  MCAuto<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  MCAuto<DataArrayDouble> array=DataArrayDouble::New();
  mcIdType nbOfCells=getNumberOfCells();
  int nbComp=getMeshDimension()+1;
  array->alloc(nbOfCells,nbComp);
  double *vals=array->getPointer();
  const mcIdType *connI=_nodal_connec_index->getConstPointer();
  const mcIdType *conn=_nodal_connec->getConstPointer();
  const double *coords=_coords->getConstPointer();
  if(getMeshDimension()==2)
    {
      if(getSpaceDimension()==3)
        {
          MCAuto<DataArrayDouble> loc=computeCellCenterOfMass();
          const double *locPtr=loc->getConstPointer();
          for(mcIdType i=0;i<nbOfCells;i++,vals+=3)
            {
              mcIdType offset=connI[i];
              INTERP_KERNEL::crossprod<3>(locPtr+3*i,coords+3*conn[offset+1],coords+3*conn[offset+2],vals);
              double n=INTERP_KERNEL::norm<3>(vals);
              std::transform(vals,vals+3,vals,std::bind(std::multiplies<double>(),std::placeholders::_1,1./n));
            }
        }
      else
        {
          MCAuto<MEDCouplingFieldDouble> isAbs=getMeasureField(false);
          const double *isAbsPtr=isAbs->getArray()->begin();
          for(mcIdType i=0;i<nbOfCells;i++,isAbsPtr++)
            { vals[3*i]=0.; vals[3*i+1]=0.; vals[3*i+2]=*isAbsPtr>0.?1.:-1.; }
        }
    }
  else//meshdimension==1
    {
      double tmp[2];
      for(mcIdType i=0;i<nbOfCells;i++)
        {
          mcIdType offset=connI[i];
          std::transform(coords+2*conn[offset+2],coords+2*conn[offset+2]+2,coords+2*conn[offset+1],tmp,std::minus<double>());
          double n=INTERP_KERNEL::norm<2>(tmp);
          std::transform(tmp,tmp+2,tmp,std::bind(std::multiplies<double>(),std::placeholders::_1,1./n));
          *vals++=-tmp[1];
          *vals++=tmp[0];
        }
    }
  ret->setArray(array);
  ret->setMesh(this);
  ret->synchronizeTimeWithSupport();
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble holding normal vectors to specified cells of
 * \a this mesh. The computed vectors have <em> this->getMeshDimension()+1 </em> components
 * and are normalized.
 * <br> \a this can be either
 * - a  2D mesh in 2D or 3D space or
 * - an 1D mesh in 2D space.
 *
 * This method avoids building explicitly a part of \a this mesh to perform the work.
 *  \param [in] begin - an array of cell ids of interest.
 *  \param [in] end - the end of \a begin, i.e. a pointer to its (last+1)-th element.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble on
 *          cells and one time. The caller is to delete this field using decrRef() as
 *          it is no more needed.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If the coordinates array is not set.
 *  \throw If the mesh dimension is not set.
 *  \throw If the mesh and space dimension is not as specified above.
 *  \sa buildOrthogonalField()
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_buildPartOrthogonalField "Here is a C++ example".<br>
 *  \ref  py_mcumesh_buildPartOrthogonalField "Here is a Python example".
 *  \endif
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::buildPartOrthogonalField(const mcIdType *begin, const mcIdType *end) const
{
  if((getMeshDimension()!=2) && (getMeshDimension()!=1 || getSpaceDimension()!=2))
    throw INTERP_KERNEL::Exception("Expected a umesh with ( meshDim == 2 spaceDim == 2 or 3 ) or ( meshDim == 1 spaceDim == 2 ) !");
  MCAuto<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  MCAuto<DataArrayDouble> array=DataArrayDouble::New();
  std::size_t nbelems=std::distance(begin,end);
  int nbComp=getMeshDimension()+1;
  array->alloc(nbelems,nbComp);
  double *vals=array->getPointer();
  const mcIdType *connI=_nodal_connec_index->getConstPointer();
  const mcIdType *conn=_nodal_connec->getConstPointer();
  const double *coords=_coords->getConstPointer();
  if(getMeshDimension()==2)
    {
      if(getSpaceDimension()==3)
        {
          MCAuto<DataArrayDouble> loc=getPartBarycenterAndOwner(begin,end);
          const double *locPtr=loc->getConstPointer();
          for(const mcIdType *i=begin;i!=end;i++,vals+=3,locPtr+=3)
            {
              mcIdType offset=connI[*i];
              INTERP_KERNEL::crossprod<3>(locPtr,coords+3*conn[offset+1],coords+3*conn[offset+2],vals);
              double n=INTERP_KERNEL::norm<3>(vals);
              std::transform(vals,vals+3,vals,std::bind(std::multiplies<double>(),std::placeholders::_1,1./n));
            }
        }
      else
        {
          for(std::size_t i=0;i<nbelems;i++)
            { vals[3*i]=0.; vals[3*i+1]=0.; vals[3*i+2]=1.; }
        }
    }
  else//meshdimension==1
    {
      double tmp[2];
      for(const mcIdType *i=begin;i!=end;i++)
        {
          mcIdType offset=connI[*i];
          std::transform(coords+2*conn[offset+2],coords+2*conn[offset+2]+2,coords+2*conn[offset+1],tmp,std::minus<double>());
          double n=INTERP_KERNEL::norm<2>(tmp);
          std::transform(tmp,tmp+2,tmp,std::bind(std::multiplies<double>(),std::placeholders::_1,1./n));
          *vals++=-tmp[1];
          *vals++=tmp[0];
        }
    }
  ret->setArray(array);
  ret->setMesh(this);
  ret->synchronizeTimeWithSupport();
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble holding a direction vector for each SEG2 in \a
 * this 1D mesh. The computed vectors have <em> this->getSpaceDimension() </em> components
 * and are \b not normalized.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble on
 *          cells and one time. The caller is to delete this field using decrRef() as
 *          it is no more needed.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If the coordinates array is not set.
 *  \throw If \a this->getMeshDimension() != 1.
 *  \throw If \a this mesh includes cells of type other than SEG2.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::buildDirectionVectorField() const
{
  if(getMeshDimension()!=1)
    throw INTERP_KERNEL::Exception("Expected a umesh with meshDim == 1 for buildDirectionVectorField !");
  if(_types.size()!=1 || *(_types.begin())!=INTERP_KERNEL::NORM_SEG2)
    throw INTERP_KERNEL::Exception("Expected a umesh with only NORM_SEG2 type of elements for buildDirectionVectorField !");
  MCAuto<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  MCAuto<DataArrayDouble> array=DataArrayDouble::New();
  mcIdType nbOfCells=getNumberOfCells();
  int spaceDim=getSpaceDimension();
  array->alloc(nbOfCells,spaceDim);
  double *pt=array->getPointer();
  const double *coo=getCoords()->getConstPointer();
  std::vector<mcIdType> conn;
  conn.reserve(2);
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      conn.resize(0);
      getNodeIdsOfCell(i,conn);
      pt=std::transform(coo+conn[1]*spaceDim,coo+(conn[1]+1)*spaceDim,coo+conn[0]*spaceDim,pt,std::minus<double>());
    }
  ret->setArray(array);
  ret->setMesh(this);
  ret->synchronizeTimeWithSupport();
  return ret.retn();
}

/*!
 * Creates a 2D mesh by cutting \a this 3D mesh with a plane. In addition to the mesh,
 * returns a new DataArrayIdType, of length equal to the number of 2D cells in the result
 * mesh, holding, for each cell in the result mesh, an id of a 3D cell it comes
 * from. If a result face is shared by two 3D cells, then the face in included twice in
 * the result mesh.
 *  \param [in] origin - 3 components of a point defining location of the plane.
 *  \param [in] vec - 3 components of a vector normal to the plane. Vector magnitude
 *         must be greater than 1e-6.
 *  \param [in] eps - half-thickness of the plane.
 *  \param [out] cellIds - a new instance of DataArrayIdType holding ids of 3D cells
 *         producing correspondent 2D cells. The caller is to delete this array
 *         using decrRef() as it is no more needed.
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. This mesh does
 *         not share the node coordinates array with \a this mesh. The caller is to
 *         delete this mesh using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If \a this->getMeshDimension() != 3 or \a this->getSpaceDimension() != 3.
 *  \throw If magnitude of \a vec is less than 1e-6.
 *  \throw If the plane does not intersect any 3D cell of \a this mesh.
 *  \throw If \a this includes quadratic cells.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildSlice3D(const double *origin, const double *vec, double eps, DataArrayIdType *&cellIds) const
{
  checkFullyDefined();
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSlice3D works on umeshes with meshdim equal to 3 and spaceDim equal to 3 too!");
  MCAuto<DataArrayIdType> candidates=getCellIdsCrossingPlane(origin,vec,eps);
  if(candidates->empty())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSlice3D : No 3D cells in this intercepts the specified plane considering bounding boxes !");
  std::vector<mcIdType> nodes;
  DataArrayIdType *cellIds1D=0;
  MCAuto<MEDCouplingUMesh> subMesh=static_cast<MEDCouplingUMesh*>(buildPartOfMySelf(candidates->begin(),candidates->end(),false));
  subMesh->findNodesOnPlane(origin,vec,eps,nodes);
  MCAuto<DataArrayIdType> desc1=DataArrayIdType::New(),desc2=DataArrayIdType::New();
  MCAuto<DataArrayIdType> descIndx1=DataArrayIdType::New(),descIndx2=DataArrayIdType::New();
  MCAuto<DataArrayIdType> revDesc1=DataArrayIdType::New(),revDesc2=DataArrayIdType::New();
  MCAuto<DataArrayIdType> revDescIndx1=DataArrayIdType::New(),revDescIndx2=DataArrayIdType::New();
  MCAuto<MEDCouplingUMesh> mDesc2=subMesh->buildDescendingConnectivity(desc2,descIndx2,revDesc2,revDescIndx2);//meshDim==2 spaceDim==3
  revDesc2=0; revDescIndx2=0;
  MCAuto<MEDCouplingUMesh> mDesc1=mDesc2->buildDescendingConnectivity(desc1,descIndx1,revDesc1,revDescIndx1);//meshDim==1 spaceDim==3
  revDesc1=0; revDescIndx1=0;
  //Marking all 1D cells that contained at least one node located on the plane
  //the intersection between those cells and the plane, which consist of the nodes previously tagged, thus don't need to be computed afterwards
  //(if said intersection is computed in MEDCouplingUMesh::split3DCurveWithPlane, then we might create additional nodes
  //due to accuracy errors when the needed nodes already exist)
  mDesc1->fillCellIdsToKeepFromNodeIds(&nodes[0],&nodes[0]+nodes.size(),false,cellIds1D);
  MCAuto<DataArrayIdType> cellIds1DTmp(cellIds1D);
  //
  std::vector<mcIdType> cut3DCurve(mDesc1->getNumberOfCells(),-2);
  for(const mcIdType *it=cellIds1D->begin();it!=cellIds1D->end();it++)
    cut3DCurve[*it]=-1;
  mDesc1->split3DCurveWithPlane(origin,vec,eps,cut3DCurve);
  std::vector< std::pair<mcIdType,mcIdType> > cut3DSurf(mDesc2->getNumberOfCells());
  AssemblyForSplitFrom3DCurve(cut3DCurve,nodes,mDesc2->getNodalConnectivity()->getConstPointer(),mDesc2->getNodalConnectivityIndex()->getConstPointer(),
                              mDesc1->getNodalConnectivity()->getConstPointer(),mDesc1->getNodalConnectivityIndex()->getConstPointer(),
                              desc1->getConstPointer(),descIndx1->getConstPointer(),cut3DSurf);
  MCAuto<DataArrayIdType> conn(DataArrayIdType::New()),connI(DataArrayIdType::New()),cellIds2(DataArrayIdType::New());
  connI->pushBackSilent(0); conn->alloc(0,1); cellIds2->alloc(0,1);
  subMesh->assemblyForSplitFrom3DSurf(cut3DSurf,desc2->getConstPointer(),descIndx2->getConstPointer(),conn,connI,cellIds2);
  if(cellIds2->empty())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSlice3D : No 3D cells in this intercepts the specified plane !");
  MCAuto<MEDCouplingUMesh> ret=MEDCouplingUMesh::New("Slice3D",2);
  ret->setCoords(mDesc1->getCoords());
  ret->setConnectivity(conn,connI,true);
  cellIds=candidates->selectByTupleId(cellIds2->begin(),cellIds2->end());
  return ret.retn();
}

/*!
 * Creates an 1D mesh by cutting \a this 2D mesh in 3D space with a plane. In
addition to the mesh, returns a new DataArrayIdType, of length equal to the number of 1D cells in the result mesh, holding, for each cell in the result mesh, an id of a 2D cell it comes
from. If a result segment is shared by two 2D cells, then the segment in included twice in
the result mesh.
 *  \param [in] origin - 3 components of a point defining location of the plane.
 *  \param [in] vec - 3 components of a vector normal to the plane. Vector magnitude
 *         must be greater than 1e-6.
 *  \param [in] eps - half-thickness of the plane.
 *  \param [out] cellIds - a new instance of DataArrayIdType holding ids of faces
 *         producing correspondent segments. The caller is to delete this array
 *         using decrRef() as it is no more needed.
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. This is an 1D
 *         mesh in 3D space. This mesh does not share the node coordinates array with
 *         \a this mesh. The caller is to delete this mesh using decrRef() as it is
 *         no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If \a this->getMeshDimension() != 2 or \a this->getSpaceDimension() != 3.
 *  \throw If magnitude of \a vec is less than 1e-6.
 *  \throw If the plane does not intersect any 2D cell of \a this mesh.
 *  \throw If \a this includes quadratic cells.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildSlice3DSurf(const double *origin, const double *vec, double eps, DataArrayIdType *&cellIds) const
{
  checkFullyDefined();
  if(getMeshDimension()!=2 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSlice3DSurf works on umeshes with meshdim equal to 2 and spaceDim equal to 3 !");
  MCAuto<DataArrayIdType> candidates(getCellIdsCrossingPlane(origin,vec,eps));
  if(candidates->empty())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSlice3DSurf : No 3D surf cells in this intercepts the specified plane considering bounding boxes !");
  std::vector<mcIdType> nodes;
  DataArrayIdType *cellIds1D(0);
  MCAuto<MEDCouplingUMesh> subMesh(buildPartOfMySelf(candidates->begin(),candidates->end(),false));
  subMesh->findNodesOnPlane(origin,vec,eps,nodes);
  MCAuto<DataArrayIdType> desc1(DataArrayIdType::New()),descIndx1(DataArrayIdType::New()),revDesc1(DataArrayIdType::New()),revDescIndx1(DataArrayIdType::New());
  MCAuto<MEDCouplingUMesh> mDesc1(subMesh->buildDescendingConnectivity(desc1,descIndx1,revDesc1,revDescIndx1));//meshDim==1 spaceDim==3
  mDesc1->fillCellIdsToKeepFromNodeIds(&nodes[0],&nodes[0]+nodes.size(),true,cellIds1D);
  MCAuto<DataArrayIdType> cellIds1DTmp(cellIds1D);
  //
  std::vector<mcIdType> cut3DCurve(mDesc1->getNumberOfCells(),-2);
  for(const mcIdType *it=cellIds1D->begin();it!=cellIds1D->end();it++)
    cut3DCurve[*it]=-1;
  mDesc1->split3DCurveWithPlane(origin,vec,eps,cut3DCurve);
  mcIdType ncellsSub=subMesh->getNumberOfCells();
  std::vector< std::pair<mcIdType,mcIdType> > cut3DSurf(ncellsSub);
  AssemblyForSplitFrom3DCurve(cut3DCurve,nodes,subMesh->getNodalConnectivity()->getConstPointer(),subMesh->getNodalConnectivityIndex()->getConstPointer(),
                              mDesc1->getNodalConnectivity()->getConstPointer(),mDesc1->getNodalConnectivityIndex()->getConstPointer(),
                              desc1->getConstPointer(),descIndx1->getConstPointer(),cut3DSurf);
  MCAuto<DataArrayIdType> conn(DataArrayIdType::New()),connI(DataArrayIdType::New()),cellIds2(DataArrayIdType::New()); connI->pushBackSilent(0);
  conn->alloc(0,1);
  const mcIdType *nodal=subMesh->getNodalConnectivity()->getConstPointer();
  const mcIdType *nodalI=subMesh->getNodalConnectivityIndex()->getConstPointer();
  for(mcIdType i=0;i<ncellsSub;i++)
    {
      if(cut3DSurf[i].first!=-1 && cut3DSurf[i].second!=-1)
        {
          if(cut3DSurf[i].first!=-2)
            {
              conn->pushBackSilent(ToIdType(INTERP_KERNEL::NORM_SEG2)); conn->pushBackSilent(cut3DSurf[i].first); conn->pushBackSilent(cut3DSurf[i].second);
              connI->pushBackSilent(conn->getNumberOfTuples());
              cellIds2->pushBackSilent(i);
            }
          else
            {
              mcIdType cellId3DSurf=cut3DSurf[i].second;
              mcIdType offset=nodalI[cellId3DSurf]+1;
              mcIdType nbOfEdges=nodalI[cellId3DSurf+1]-offset;
              for(mcIdType j=0;j<nbOfEdges;j++)
                {
                  conn->pushBackSilent(ToIdType(INTERP_KERNEL::NORM_SEG2)); conn->pushBackSilent(nodal[offset+j]); conn->pushBackSilent(nodal[offset+(j+1)%nbOfEdges]);
                  connI->pushBackSilent(conn->getNumberOfTuples());
                  cellIds2->pushBackSilent(cellId3DSurf);
                }
            }
        }
    }
  if(cellIds2->empty())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSlice3DSurf : No 3DSurf cells in this intercepts the specified plane !");
  MCAuto<MEDCouplingUMesh> ret=MEDCouplingUMesh::New("Slice3DSurf",1);
  ret->setCoords(mDesc1->getCoords());
  ret->setConnectivity(conn,connI,true);
  cellIds=candidates->selectByTupleId(cellIds2->begin(),cellIds2->end());
  return ret.retn();
}

MCAuto<MEDCouplingUMesh> MEDCouplingUMesh::clipSingle3DCellByPlane(const double origin[3], const double vec[3], double eps) const
{
  checkFullyDefined();
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::clipSingle3DCellByPlane works on umeshes with meshdim equal to 3 and spaceDim equal to 3 too!");
  if(getNumberOfCells()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::clipSingle3DCellByPlane works only on mesh containing exactly one cell !");
  //
  std::vector<mcIdType> nodes;
  findNodesOnPlane(origin,vec,eps,nodes);
  MCAuto<DataArrayIdType> desc1(DataArrayIdType::New()),desc2(DataArrayIdType::New()),descIndx1(DataArrayIdType::New()),descIndx2(DataArrayIdType::New()),revDesc1(DataArrayIdType::New()),revDesc2(DataArrayIdType::New()),revDescIndx1(DataArrayIdType::New()),revDescIndx2(DataArrayIdType::New());
  MCAuto<MEDCouplingUMesh> mDesc2(buildDescendingConnectivity(desc2,descIndx2,revDesc2,revDescIndx2));//meshDim==2 spaceDim==3
  revDesc2=0; revDescIndx2=0;
  MCAuto<MEDCouplingUMesh> mDesc1(mDesc2->buildDescendingConnectivity(desc1,descIndx1,revDesc1,revDescIndx1));//meshDim==1 spaceDim==3
  revDesc1=0; revDescIndx1=0;
  DataArrayIdType *cellIds1D(0);
  mDesc1->fillCellIdsToKeepFromNodeIds(&nodes[0],&nodes[0]+nodes.size(),true,cellIds1D);
  MCAuto<DataArrayIdType> cellIds1DTmp(cellIds1D);
  std::vector<mcIdType> cut3DCurve(mDesc1->getNumberOfCells(),-2);
  for(const mcIdType *it=cellIds1D->begin();it!=cellIds1D->end();it++)
    cut3DCurve[*it]=-1;
  bool sameNbNodes;
  {
    mcIdType oldNbNodes(mDesc1->getNumberOfNodes());
    mDesc1->split3DCurveWithPlane(origin,vec,eps,cut3DCurve);
    sameNbNodes=(mDesc1->getNumberOfNodes()==oldNbNodes);
  }
  std::vector< std::pair<mcIdType,mcIdType> > cut3DSurf(mDesc2->getNumberOfCells());
  AssemblyForSplitFrom3DCurve(cut3DCurve,nodes,mDesc2->getNodalConnectivity()->begin(),mDesc2->getNodalConnectivityIndex()->begin(),
                              mDesc1->getNodalConnectivity()->begin(),mDesc1->getNodalConnectivityIndex()->begin(),
                              desc1->begin(),descIndx1->begin(),cut3DSurf);
  MCAuto<DataArrayIdType> conn(DataArrayIdType::New()),connI(DataArrayIdType::New());
  connI->pushBackSilent(0); conn->alloc(0,1);
  {
    MCAuto<DataArrayIdType> cellIds2(DataArrayIdType::New()); cellIds2->alloc(0,1);
    assemblyForSplitFrom3DSurf(cut3DSurf,desc2->begin(),descIndx2->begin(),conn,connI,cellIds2);
    if(cellIds2->empty())
      throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSlice3D : No 3D cells in this intercepts the specified plane !");
  }
  std::vector<std::vector<mcIdType> > res;
  buildSubCellsFromCut(cut3DSurf,desc2->begin(),descIndx2->begin(),mDesc1->getCoords()->begin(),eps,res);
  std::size_t sz(res.size());
  if(ToIdType(res.size())==mDesc1->getNumberOfCells() && sameNbNodes)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::clipSingle3DCellByPlane : cell is not clipped !");
  for(std::size_t i=0;i<sz;i++)
    {
      conn->pushBackSilent(ToIdType(INTERP_KERNEL::NORM_POLYGON));
      conn->insertAtTheEnd(res[i].begin(),res[i].end());
      connI->pushBackSilent(conn->getNumberOfTuples());
    }
  MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::New("",2));
  ret->setCoords(mDesc1->getCoords());
  ret->setConnectivity(conn,connI,true);
  mcIdType nbCellsRet(ret->getNumberOfCells());
  //
  MCAuto<DataArrayDouble> vec2(DataArrayDouble::New()); vec2->alloc(1,3); std::copy(vec,vec+3,vec2->getPointer());
  MCAuto<MEDCouplingFieldDouble> ortho(ret->buildOrthogonalField());
  MCAuto<DataArrayDouble> ortho2(ortho->getArray()->selectByTupleIdSafeSlice(0,1,1));
  MCAuto<DataArrayDouble> dott(DataArrayDouble::Dot(ortho2,vec2));
  MCAuto<DataArrayDouble> ccm(ret->computeCellCenterOfMass());
  MCAuto<DataArrayDouble> occm;
  {
    MCAuto<DataArrayDouble> pt(DataArrayDouble::New()); pt->alloc(1,3); std::copy(origin,origin+3,pt->getPointer());
    occm=DataArrayDouble::Substract(ccm,pt);
  }
  vec2=DataArrayDouble::New(); vec2->alloc(nbCellsRet,3);
  vec2->setPartOfValuesSimple1(vec[0],0,nbCellsRet,1,0,1,1); vec2->setPartOfValuesSimple1(vec[1],0,nbCellsRet,1,1,2,1); vec2->setPartOfValuesSimple1(vec[2],0,nbCellsRet,1,2,3,1);
  MCAuto<DataArrayDouble> dott2(DataArrayDouble::Dot(occm,vec2));
  //
  const mcIdType *cPtr(ret->getNodalConnectivity()->begin()),*ciPtr(ret->getNodalConnectivityIndex()->begin());
  MCAuto<MEDCouplingUMesh> ret2(MEDCouplingUMesh::New("Clip3D",3));
  ret2->setCoords(mDesc1->getCoords());
  MCAuto<DataArrayIdType> conn2(DataArrayIdType::New()),conn2I(DataArrayIdType::New());
  conn2I->pushBackSilent(0); conn2->alloc(0,1);
  std::vector<mcIdType> cell0(1,ToIdType(INTERP_KERNEL::NORM_POLYHED));
  std::vector<mcIdType> cell1(1,ToIdType(INTERP_KERNEL::NORM_POLYHED));
  if(dott->getIJ(0,0)>0)
    {
      cell0.insert(cell0.end(),cPtr+1,cPtr+ciPtr[1]);
      std::reverse_copy(cPtr+1,cPtr+ciPtr[1],std::inserter(cell1,cell1.end()));
    }
  else
    {
      cell1.insert(cell1.end(),cPtr+1,cPtr+ciPtr[1]);
      std::reverse_copy(cPtr+1,cPtr+ciPtr[1],std::inserter(cell0,cell0.end()));
    }
  for(mcIdType i=1;i<nbCellsRet;i++)
    {
      if(dott2->getIJ(i,0)<0)
        {
          if(ciPtr[i+1]-ciPtr[i]>=4)
            {
              cell0.push_back(-1);
              cell0.insert(cell0.end(),cPtr+ciPtr[i]+1,cPtr+ciPtr[i+1]);
            }
        }
      else
        {
          if(ciPtr[i+1]-ciPtr[i]>=4)
            {
              cell1.push_back(-1);
              cell1.insert(cell1.end(),cPtr+ciPtr[i]+1,cPtr+ciPtr[i+1]);
            }
        }
    }
  conn2->insertAtTheEnd(cell0.begin(),cell0.end());
  conn2I->pushBackSilent(conn2->getNumberOfTuples());
  conn2->insertAtTheEnd(cell1.begin(),cell1.end());
  conn2I->pushBackSilent(conn2->getNumberOfTuples());
  ret2->setConnectivity(conn2,conn2I,true);
  ret2->checkConsistencyLight();
  ret2->orientCorrectlyPolyhedrons();
  return ret2;
}

/*!
 * Finds cells whose bounding boxes intersect a given plane.
 *  \param [in] origin - 3 components of a point defining location of the plane.
 *  \param [in] vec - 3 components of a vector normal to the plane. Vector magnitude
 *         must be greater than 1e-6.
 *  \param [in] eps - half-thickness of the plane.
 *  \return DataArrayIdType * - a new instance of DataArrayIdType holding ids of the found
 *         cells. The caller is to delete this array using decrRef() as it is no more
 *         needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If \a this->getSpaceDimension() != 3.
 *  \throw If magnitude of \a vec is less than 1e-6.
 *  \sa buildSlice3D()
 */
DataArrayIdType *MEDCouplingUMesh::getCellIdsCrossingPlane(const double *origin, const double *vec, double eps) const
{
  checkFullyDefined();
  if(getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSlice3D works on umeshes with spaceDim equal to 3 !");
  double normm=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  if(normm<1e-6)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getCellIdsCrossingPlane : parameter 'vec' should have a norm2 greater than 1e-6 !");
  double vec2[3];
  vec2[0]=vec[1]; vec2[1]=-vec[0]; vec2[2]=0.;//vec2 is the result of cross product of vec with (0,0,1)
  double angle=acos(vec[2]/normm);
  MCAuto<DataArrayIdType> cellIds;
  double bbox[6];
  if(angle>eps)
    {
      MCAuto<DataArrayDouble> coo=_coords->deepCopy();
      double normm2(sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2]));
      if(normm2/normm>1e-6)
        DataArrayDouble::Rotate3DAlg(origin,vec2,angle,coo->getNumberOfTuples(),coo->getPointer(),coo->getPointer());
      MCAuto<MEDCouplingUMesh> mw=clone(false);//false -> shallow copy
      mw->setCoords(coo);
      mw->getBoundingBox(bbox);
      bbox[4]=origin[2]-eps; bbox[5]=origin[2]+eps;
      cellIds=mw->getCellsInBoundingBox(bbox,eps);
    }
  else
    {
      getBoundingBox(bbox);
      bbox[4]=origin[2]-eps; bbox[5]=origin[2]+eps;
      cellIds=getCellsInBoundingBox(bbox,eps);
    }
  return cellIds.retn();
}

/*!
 * This method checks that \a this is a contiguous mesh. The user is expected to call this method on a mesh with meshdim==1.
 * If not an exception will thrown. If this is an empty mesh with no cell an exception will be thrown too.
 * No consideration of coordinate is done by this method.
 * A 1D mesh is said contiguous if : a cell i with nodal connectivity (k,p) the cell i+1 the nodal connectivity should be (p,m)
 * If not false is returned. In case that false is returned a call to MEDCoupling::MEDCouplingUMesh::mergeNodes could be useful.
 */
bool MEDCouplingUMesh::isContiguous1D() const
{
  if(getMeshDimension()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::isContiguous1D : this method has a sense only for 1D mesh !");
  mcIdType nbCells=getNumberOfCells();
  if(nbCells<1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::isContiguous1D : this method has a sense for non empty mesh !");
  const mcIdType *connI(_nodal_connec_index->begin()),*conn(_nodal_connec->begin());
  mcIdType ref=conn[connI[0]+2];
  for(mcIdType i=1;i<nbCells;i++)
    {
      if(conn[connI[i]+1]!=ref)
        return false;
      ref=conn[connI[i]+2];
    }
  return true;
}

/*!
 * This method is only callable on mesh with meshdim == 1 containing only SEG2 and spaceDim==3.
 * This method projects this on the 3D line defined by (pt,v). This methods first checks that all SEG2 are along v vector.
 * \param pt reference point of the line
 * \param v normalized director vector of the line
 * \param eps max precision before throwing an exception
 * \param res output of size this->getNumberOfCells
 */
void MEDCouplingUMesh::project1D(const double *pt, const double *v, double eps, double *res) const
{
  if(getMeshDimension()!=1)
    throw INTERP_KERNEL::Exception("Expected a umesh with meshDim == 1 for project1D !");
  if(_types.size()!=1 || *(_types.begin())!=INTERP_KERNEL::NORM_SEG2)
    throw INTERP_KERNEL::Exception("Expected a umesh with only NORM_SEG2 type of elements for project1D !");
  if(getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("Expected a umesh with spaceDim==3 for project1D !");
  MCAuto<MEDCouplingFieldDouble> f=buildDirectionVectorField();
  const double *fPtr=f->getArray()->getConstPointer();
  double tmp[3];
  for(mcIdType i=0;i<getNumberOfCells();i++)
    {
      const double *tmp1=fPtr+3*i;
      tmp[0]=tmp1[1]*v[2]-tmp1[2]*v[1];
      tmp[1]=tmp1[2]*v[0]-tmp1[0]*v[2];
      tmp[2]=tmp1[0]*v[1]-tmp1[1]*v[0];
      double n1=INTERP_KERNEL::norm<3>(tmp);
      n1/=INTERP_KERNEL::norm<3>(tmp1);
      if(n1>eps)
        throw INTERP_KERNEL::Exception("UMesh::Projection 1D failed !");
    }
  const double *coo=getCoords()->getConstPointer();
  for(mcIdType i=0;i<getNumberOfNodes();i++)
    {
      std::transform(coo+i*3,coo+i*3+3,pt,tmp,std::minus<double>());
      std::transform(tmp,tmp+3,v,tmp,std::multiplies<double>());
      res[i]=std::accumulate(tmp,tmp+3,0.);
    }
}

/*!
 * This method computes the distance from a point \a pt to \a this and the first \a cellId in \a this corresponding to the returned distance.
 * \a this is expected to be a mesh so that its space dimension is equal to its
 * mesh dimension + 1. Furthermore only mesh dimension 1 and 2 are supported for the moment.
 * Distance from \a ptBg to \a ptEnd is expected to be equal to the space dimension. \a this is also expected to be fully defined (connectivity and coordinates).
 *
 * WARNING, if there is some orphan nodes in \a this (nodes not fetched by any cells in \a this ( see MEDCouplingUMesh::zipCoords ) ) these nodes will ** not ** been taken
 * into account in this method. Only cells and nodes lying on them are considered in the algorithm (even if one of these orphan nodes is closer than returned distance).
 * A user that needs to consider orphan nodes should invoke DataArrayDouble::minimalDistanceTo method on the coordinates array of \a this.
 *
 * So this method is more accurate (so, more costly) than simply searching for the closest point in \a this.
 * If only this information is enough for you simply call \c getCoords()->distanceToTuple on \a this.
 *
 * \param [in] ptBg the start pointer (included) of the coordinates of the point
 * \param [in] ptEnd the end pointer (not included) of the coordinates of the point
 * \param [out] cellId that corresponds to minimal distance. If the closer node is not linked to any cell in \a this -1 is returned.
 * \return the positive value of the distance.
 * \throw if distance from \a ptBg to \a ptEnd is not equal to the space dimension. An exception is also thrown if mesh dimension of \a this is not equal to space
 * dimension - 1.
 * \sa DataArrayDouble::distanceToTuple, MEDCouplingUMesh::distanceToPoints
 */
double MEDCouplingUMesh::distanceToPoint(const double *ptBg, const double *ptEnd, mcIdType& cellId) const
{
  int meshDim=getMeshDimension(),spaceDim=getSpaceDimension();
  if(meshDim!=spaceDim-1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::distanceToPoint works only for spaceDim=meshDim+1 !");
  if(meshDim!=2 && meshDim!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::distanceToPoint : only mesh dimension 2 and 1 are implemented !");
  checkFullyDefined();
  if(ToIdType(std::distance(ptBg,ptEnd))!=spaceDim)
    { std::ostringstream oss; oss << "MEDCouplingUMesh::distanceToPoint : input point has to have dimension equal to the space dimension of this (" << spaceDim << ") !"; throw INTERP_KERNEL::Exception(oss.str()); }
  DataArrayIdType *ret1=0;
  MCAuto<DataArrayDouble> pts=DataArrayDouble::New(); pts->useArray(ptBg,false,DeallocType::C_DEALLOC,1,spaceDim);
  MCAuto<DataArrayDouble> ret0=distanceToPoints(pts,ret1);
  MCAuto<DataArrayIdType> ret1Safe(ret1);
  cellId=*ret1Safe->begin();
  return *ret0->begin();
}

/*!
 * This method computes the distance from each point of points serie \a pts (stored in a DataArrayDouble in which each tuple represents a point)
 *  to \a this  and the first \a cellId in \a this corresponding to the returned distance.
 * WARNING, if there is some orphan nodes in \a this (nodes not fetched by any cells in \a this ( see MEDCouplingUMesh::zipCoords ) ) these nodes will ** not ** been taken
 * into account in this method. Only cells and nodes lying on them are considered in the algorithm (even if one of these orphan nodes is closer than returned distance).
 * A user that needs to consider orphan nodes should invoke DataArrayDouble::minimalDistanceTo method on the coordinates array of \a this.
 *
 * \a this is expected to be a mesh so that its space dimension is equal to its
 * mesh dimension + 1. Furthermore only mesh dimension 1 and 2 are supported for the moment.
 * Number of components of \a pts is expected to be equal to the space dimension. \a this is also expected to be fully defined (connectivity and coordinates).
 *
 * So this method is more accurate (so, more costly) than simply searching for each point in \a pts the closest point in \a this.
 * If only this information is enough for you simply call \c getCoords()->distanceToTuple on \a this.
 *
 * \param [in] pts the list of points in which each tuple represents a point
 * \param [out] cellIds a newly allocated object that tells for each point in \a pts the first cell id in \a this that minimizes the distance.
 * \return a newly allocated object to be dealed by the caller that tells for each point in \a pts the distance to \a this.
 * \throw if number of components of \a pts is not equal to the space dimension.
 * \throw if mesh dimension of \a this is not equal to space dimension - 1.
 * \sa DataArrayDouble::distanceToTuple, MEDCouplingUMesh::distanceToPoint
 */
DataArrayDouble *MEDCouplingUMesh::distanceToPoints(const DataArrayDouble *pts, DataArrayIdType *& cellIds) const
{
  if(!pts)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::distanceToPoints : input points pointer is NULL !");
  pts->checkAllocated();
  int meshDim=getMeshDimension(),spaceDim=getSpaceDimension();
  if(meshDim!=spaceDim-1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::distanceToPoints works only for spaceDim=meshDim+1 !");
  if(meshDim!=2 && meshDim!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::distanceToPoints : only mesh dimension 2 and 1 are implemented !");
  if(ToIdType(pts->getNumberOfComponents())!=spaceDim)
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::distanceToPoints : input pts DataArrayDouble has " << pts->getNumberOfComponents() << " components whereas it should be equal to " << spaceDim << " (mesh spaceDimension) !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  checkFullyDefined();
  mcIdType nbCells=getNumberOfCells();
  if(nbCells==0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::distanceToPoints : no cells in this !");
  mcIdType nbOfPts=pts->getNumberOfTuples();
  MCAuto<DataArrayDouble> ret0=DataArrayDouble::New(); ret0->alloc(nbOfPts,1);
  MCAuto<DataArrayIdType> ret1=DataArrayIdType::New(); ret1->alloc(nbOfPts,1);
  const mcIdType *nc=_nodal_connec->begin(),*ncI=_nodal_connec_index->begin(); const double *coords=_coords->begin();
  double *ret0Ptr=ret0->getPointer(); mcIdType *ret1Ptr=ret1->getPointer(); const double *ptsPtr=pts->begin();
  MCAuto<DataArrayDouble> bboxArr(getBoundingBoxForBBTree());
  const double *bbox(bboxArr->begin());
  switch(spaceDim)
  {
    case 3:
      {
        BBTreeDst<3> myTree(bbox,0,0,nbCells);
        for(mcIdType i=0;i<nbOfPts;i++,ret0Ptr++,ret1Ptr++,ptsPtr+=3)
          {
            double x=std::numeric_limits<double>::max();
            std::vector<mcIdType> elems;
            myTree.getMinDistanceOfMax(ptsPtr,x);
            myTree.getElemsWhoseMinDistanceToPtSmallerThan(ptsPtr,x,elems);
            DistanceToPoint3DSurfAlg(ptsPtr,&elems[0],&elems[0]+elems.size(),coords,nc,ncI,*ret0Ptr,*ret1Ptr);
          }
        break;
      }
    case 2:
      {
        BBTreeDst<2> myTree(bbox,0,0,nbCells);
        for(mcIdType i=0;i<nbOfPts;i++,ret0Ptr++,ret1Ptr++,ptsPtr+=2)
          {
            double x=std::numeric_limits<double>::max();
            std::vector<mcIdType> elems;
            myTree.getMinDistanceOfMax(ptsPtr,x);
            myTree.getElemsWhoseMinDistanceToPtSmallerThan(ptsPtr,x,elems);
            DistanceToPoint2DCurveAlg(ptsPtr,&elems[0],&elems[0]+elems.size(),coords,nc,ncI,*ret0Ptr,*ret1Ptr);
          }
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingUMesh::distanceToPoints : only spacedim 2 and 3 supported !");
  }
  cellIds=ret1.retn();
  return ret0.retn();
}

/// @cond INTERNAL

/// @endcond

/*!
 * Finds cells in contact with a ball (i.e. a point with precision).
 * For speed reasons, the INTERP_KERNEL::NORM_QUAD4, INTERP_KERNEL::NORM_TRI6 and INTERP_KERNEL::NORM_QUAD8 cells are considered as convex cells to detect if a point is IN or OUT.
 * If it is not the case, please change their types to INTERP_KERNEL::NORM_POLYGON or INTERP_KERNEL::NORM_QPOLYG before invoking this method.
 *
 * \warning This method is suitable if the caller intends to evaluate only one
 *          point, for more points getCellsContainingPoints() is recommended as it is
 *          faster.
 *  \param [in] pos - array of coordinates of the ball central point.
 *  \param [in] eps - ball radius.
 *  \return mcIdType - a smallest id of cells being in contact with the ball, -1 in case
 *         if there are no such cells.
 *  \throw If the coordinates array is not set.
 *  \throw If \a this->getMeshDimension() != \a this->getSpaceDimension().
 */
mcIdType MEDCouplingUMesh::getCellContainingPoint(const double *pos, double eps) const
{
  std::vector<mcIdType> elts;
  getCellsContainingPoint(pos,eps,elts);
  if(elts.empty())
    return -1;
  return elts.front();
}

/*!
 * Finds cells in contact with a ball (i.e. a point with precision).
 * For speed reasons, the INTERP_KERNEL::NORM_QUAD4, INTERP_KERNEL::NORM_TRI6 and INTERP_KERNEL::NORM_QUAD8 cells are considered as convex cells to detect if a point is IN or OUT.
 * If it is not the case, please change their types to INTERP_KERNEL::NORM_POLYGON or INTERP_KERNEL::NORM_QPOLYG before invoking this method.
 * \warning This method is suitable if the caller intends to evaluate only one
 *          point, for more points getCellsContainingPoints() is recommended as it is
 *          faster.
 *  \param [in] pos - array of coordinates of the ball central point.
 *  \param [in] eps - ball radius (i.e. the precision) relative to max dimension along X, Y or Z of the the considered cell bounding box.
 *  \param [out] elts - vector returning ids of the found cells. It is cleared
 *         before inserting ids.
 *  \throw If the coordinates array is not set.
 *  \throw If \a this->getMeshDimension() != \a this->getSpaceDimension().
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_getCellsContainingPoint "Here is a C++ example".<br>
 *  \ref  py_mcumesh_getCellsContainingPoint "Here is a Python example".
 *  \endif
 */
void MEDCouplingUMesh::getCellsContainingPoint(const double *pos, double eps, std::vector<mcIdType>& elts) const
{
  MCAuto<DataArrayIdType> eltsUg,eltsIndexUg;
  getCellsContainingPoints(pos,1,eps,eltsUg,eltsIndexUg);
  elts.clear(); elts.insert(elts.end(),eltsUg->begin(),eltsUg->end());
}

void MEDCouplingUMesh::getCellsContainingPointsZeAlg(const double *pos, mcIdType nbOfPoints, double eps,
                                                     MCAuto<DataArrayIdType>& elts, MCAuto<DataArrayIdType>& eltsIndex,
                                                     std::function<bool(INTERP_KERNEL::NormalizedCellType,mcIdType)> sensibilityTo2DQuadraticLinearCellsFunc) const
{
  int spaceDim(getSpaceDimension()),mDim(getMeshDimension());
  if(spaceDim==3)
    {
      if(mDim==3)
        {
          const double *coords=_coords->getConstPointer();
          getCellsContainingPointsAlg<3>(coords,pos,nbOfPoints,eps,elts,eltsIndex,sensibilityTo2DQuadraticLinearCellsFunc);
        }
      else
        throw INTERP_KERNEL::Exception("For spaceDim==3 only meshDim==3 implemented for getelementscontainingpoints !");
    }
  else if(spaceDim==2)
    {
      if(mDim==2)
        {
          const double *coords=_coords->getConstPointer();
          getCellsContainingPointsAlg<2>(coords,pos,nbOfPoints,eps,elts,eltsIndex,sensibilityTo2DQuadraticLinearCellsFunc);
        }
      else
        throw INTERP_KERNEL::Exception("For spaceDim==2 only meshDim==2 implemented for getelementscontainingpoints !");
    }
  else if(spaceDim==1)
    {
      if(mDim==1)
        {
          const double *coords=_coords->getConstPointer();
          getCellsContainingPointsAlg<1>(coords,pos,nbOfPoints,eps,elts,eltsIndex,sensibilityTo2DQuadraticLinearCellsFunc);
        }
      else
        throw INTERP_KERNEL::Exception("For spaceDim==1 only meshDim==1 implemented for getelementscontainingpoints !");
    }
  else
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getCellsContainingPoints : not managed for mdim not in [1,2,3] !");
}

/*!
 * Finds cells in contact with several balls (i.e. points with precision).
 * This method is an extension of getCellContainingPoint() and
 * getCellsContainingPoint() for the case of multiple points.
 * For speed reasons, the INTERP_KERNEL::NORM_QUAD4, INTERP_KERNEL::NORM_TRI6 and INTERP_KERNEL::NORM_QUAD8 cells are considered as convex cells to detect if a point is IN or OUT.
 * If it is not the case, please change their types to INTERP_KERNEL::NORM_POLYGON or INTERP_KERNEL::NORM_QPOLYG before invoking this method.
 *  \param [in] pos - an array of coordinates of points in full interlace mode :
 *         X0,Y0,Z0,X1,Y1,Z1,... Size of the array must be \a
 *         this->getSpaceDimension() * \a nbOfPoints
 *  \param [in] nbOfPoints - number of points to locate within \a this mesh.
 *  \param [in] eps - ball radius (i.e. the precision) relative to max dimension along X, Y or Z of the the considered cell bounding box.
 *  \param [out] elts - vector returning ids of found cells.
 *  \param [out] eltsIndex - an array, of length \a nbOfPoints + 1,
 *         dividing cell ids in \a elts into groups each referring to one
 *         point. Its every element (except the last one) is an index pointing to the
 *         first id of a group of cells. For example cells in contact with the *i*-th
 *         point are described by following range of indices:
 *         [ \a eltsIndex[ *i* ], \a eltsIndex[ *i*+1 ] ) and the cell ids are
 *         \a elts[ \a eltsIndex[ *i* ]], \a elts[ \a eltsIndex[ *i* ] + 1 ], ...
 *         Number of cells in contact with the *i*-th point is
 *         \a eltsIndex[ *i*+1 ] - \a eltsIndex[ *i* ].
 *  \throw If the coordinates array is not set.
 *  \throw If \a this->getMeshDimension() != \a this->getSpaceDimension().
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_getCellsContainingPoints "Here is a C++ example".<br>
 *  \ref  py_mcumesh_getCellsContainingPoints "Here is a Python example".
 *  \endif
 */
void MEDCouplingUMesh::getCellsContainingPoints(const double *pos, mcIdType nbOfPoints, double eps,
                                                MCAuto<DataArrayIdType>& elts, MCAuto<DataArrayIdType>& eltsIndex) const
{
  auto yesImSensibleTo2DQuadraticLinearCellsFunc([](INTERP_KERNEL::NormalizedCellType ct, int mdim) { return INTERP_KERNEL::CellModel::GetCellModel(ct).isQuadratic() && mdim == 2; } );
  this->getCellsContainingPointsZeAlg(pos,nbOfPoints,eps,elts,eltsIndex,yesImSensibleTo2DQuadraticLinearCellsFunc);
}

/*!
 * Behaves like MEDCouplingMesh::getCellsContainingPoints for cells in \a this that are linear.
 * For quadratic cells in \a this, this method behaves by just considering linear part of cells.
 * This method is here only for backward compatibility (interpolation GaussPoints to GaussPoints).
 * 
 * \sa MEDCouplingUMesh::getCellsContainingPoints, MEDCouplingRemapper::prepareNotInterpKernelOnlyGaussGauss
 */
void MEDCouplingUMesh::getCellsContainingPointsLinearPartOnlyOnNonDynType(const double *pos, mcIdType nbOfPoints, double eps, MCAuto<DataArrayIdType>& elts, MCAuto<DataArrayIdType>& eltsIndex) const
{
  auto noImNotSensibleTo2DQuadraticLinearCellsFunc([](INTERP_KERNEL::NormalizedCellType,mcIdType) { return false; } );
  this->getCellsContainingPointsZeAlg(pos,nbOfPoints,eps,elts,eltsIndex,noImNotSensibleTo2DQuadraticLinearCellsFunc);
}

/*!
 * Finds butterfly cells in \a this mesh. A 2D cell is considered to be butterfly if at
 * least two its edges intersect each other anywhere except their extremities. An
 * INTERP_KERNEL::NORM_NORI3 cell can \b not be butterfly.
 *  \param [in,out] cells - a vector returning ids of the found cells. It is not
 *         cleared before filling in.
 *  \param [in] eps - precision.
 *  \throw If \a this->getMeshDimension() != 2.
 *  \throw If \a this->getSpaceDimension() != 2 && \a this->getSpaceDimension() != 3.
 */
void MEDCouplingUMesh::checkButterflyCells(std::vector<mcIdType>& cells, double eps) const
{
  const char msg[]="Butterfly detection work only for 2D cells with spaceDim==2 or 3!";
  if(getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception(msg);
  int spaceDim=getSpaceDimension();
  if(spaceDim!=2 && spaceDim!=3)
    throw INTERP_KERNEL::Exception(msg);
  const mcIdType *conn=_nodal_connec->getConstPointer();
  const mcIdType *connI=_nodal_connec_index->getConstPointer();
  mcIdType nbOfCells=getNumberOfCells();
  std::vector<double> cell2DinS2;
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      mcIdType offset=connI[i];
      mcIdType nbOfNodesForCell=connI[i+1]-offset-1;
      if(nbOfNodesForCell<=3)
        continue;
      bool isQuad=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[offset]).isQuadratic();
      project2DCellOnXY(conn+offset+1,conn+connI[i+1],cell2DinS2);
      if(isButterfly2DCell(cell2DinS2,isQuad,eps))
        cells.push_back(i);
      cell2DinS2.clear();
    }
}

/*!
 * This method is typically requested to unbutterfly 2D linear cells in \b this.
 *
 * This method expects that space dimension is equal to 2 and mesh dimension is equal to 2 too. If it is not the case an INTERP_KERNEL::Exception will be thrown.
 * This method works only for linear 2D cells. If there is any of non linear cells (INTERP_KERNEL::NORM_QUAD8 for example) an INTERP_KERNEL::Exception will be thrown too.
 *
 * For each 2D linear cell in \b this, this method builds the convex envelop (or the convex hull) of the current cell.
 * This convex envelop is computed using Jarvis march algorithm.
 * The coordinates and the number of cells of \b this remain unchanged on invocation of this method.
 * Only connectivity of some cells could be modified if those cells were not representing a convex envelop. If a cell already equals its convex envelop (regardless orientation)
 * its connectivity will remain unchanged. If the computation leads to a modification of nodal connectivity of a cell its geometric type will be modified to INTERP_KERNEL::NORM_POLYGON.
 *
 * \return a newly allocated array containing cellIds that have been modified if any. If no cells have been impacted by this method NULL is returned.
 * \sa MEDCouplingUMesh::colinearize2D
 */
DataArrayIdType *MEDCouplingUMesh::convexEnvelop2D()
{
  if(getMeshDimension()!=2 || getSpaceDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::convexEnvelop2D  works only for meshDim=2 and spaceDim=2 !");
  checkFullyDefined();
  const double *coords=getCoords()->getConstPointer();
  mcIdType nbOfCells=getNumberOfCells();
  MCAuto<DataArrayIdType> nodalConnecIndexOut=DataArrayIdType::New();
  nodalConnecIndexOut->alloc(nbOfCells+1,1);
  MCAuto<DataArrayIdType> nodalConnecOut(DataArrayIdType::New());
  mcIdType *workIndexOut=nodalConnecIndexOut->getPointer();
  *workIndexOut=0;
  const mcIdType *nodalConnecIn=_nodal_connec->getConstPointer();
  const mcIdType *nodalConnecIndexIn=_nodal_connec_index->getConstPointer();
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  MCAuto<DataArrayIdType> isChanged(DataArrayIdType::New());
  isChanged->alloc(0,1);
  for(mcIdType i=0;i<nbOfCells;i++,workIndexOut++)
    {
      mcIdType pos=nodalConnecOut->getNumberOfTuples();
      if(BuildConvexEnvelopOf2DCellJarvis(coords,nodalConnecIn+nodalConnecIndexIn[i],nodalConnecIn+nodalConnecIndexIn[i+1],nodalConnecOut))
        isChanged->pushBackSilent(i);
      types.insert((INTERP_KERNEL::NormalizedCellType)nodalConnecOut->getIJ(pos,0));
      workIndexOut[1]=nodalConnecOut->getNumberOfTuples();
    }
  if(isChanged->empty())
    return 0;
  setConnectivity(nodalConnecOut,nodalConnecIndexOut,false);
  _types=types;
  return isChanged.retn();
}

/*!
 * This method is \b NOT const because it can modify \a this.
 * \a this is expected to be an unstructured mesh with meshDim==2 and spaceDim==3. If not an exception will be thrown.
 * \param mesh1D is an unstructured mesh with MeshDim==1 and spaceDim==3. If not an exception will be thrown.
 * \param policy specifies the type of extrusion chosen:
 *   - \b 0 for translation only (most simple): the cells of the 1D mesh represent the vectors along which the 2D mesh
 *   will be repeated to build each level
 *   - \b 1 for translation and rotation: the translation is done as above. For each level, an arc of circle is fitted on
 *   the 3 preceding points of the 1D mesh. The center of the arc is the center of rotation for each level, the rotation is done
 *   along an axis normal to the plane containing the arc, and finally the angle of rotation is defined by the first two points on the
 *   arc.
 * \return an unstructured mesh with meshDim==3 and spaceDim==3. The returned mesh has the same coords than \a this.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildExtrudedMesh(const MEDCouplingUMesh *mesh1D, int policy)
{
  checkFullyDefined();
  mesh1D->checkFullyDefined();
  if(!mesh1D->isContiguous1D())
    throw INTERP_KERNEL::Exception("buildExtrudedMesh : 1D mesh passed in parameter is not contiguous !");
  if(getSpaceDimension()!=mesh1D->getSpaceDimension())
    throw INTERP_KERNEL::Exception("Invalid call to buildExtrudedMesh this and mesh1D must have same space dimension !");
  if((getMeshDimension()!=2 || getSpaceDimension()!=3) && (getMeshDimension()!=1 || getSpaceDimension()!=2))
    throw INTERP_KERNEL::Exception("Invalid 'this' for buildExtrudedMesh method : must be (meshDim==2 and spaceDim==3) or (meshDim==1 and spaceDim==2) !");
  if(mesh1D->getMeshDimension()!=1)
    throw INTERP_KERNEL::Exception("Invalid 'mesh1D' for buildExtrudedMesh method : must be meshDim==1 !");
  bool isQuad=false;
  if(isPresenceOfQuadratic())
    {
      if(mesh1D->isFullyQuadratic())
        isQuad=true;
      else
        throw INTERP_KERNEL::Exception("Invalid 2D mesh and 1D mesh because 2D mesh has quadratic cells and 1D is not fully quadratic !");
    }
  mcIdType oldNbOfNodes(getNumberOfNodes());
  MCAuto<DataArrayDouble> newCoords;
  switch(policy)
  {
    case 0:
      {
        newCoords=fillExtCoordsUsingTranslation(mesh1D,isQuad);
        break;
      }
    case 1:
      {
        newCoords=fillExtCoordsUsingTranslAndAutoRotation(mesh1D,isQuad);
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("Not implemented extrusion policy : must be in (0) !");
  }
  setCoords(newCoords);
  MCAuto<MEDCouplingUMesh> ret(buildExtrudedMeshFromThisLowLev(oldNbOfNodes,isQuad));
  updateTime();
  return ret.retn();
}


/*!
 * Checks if \a this mesh is constituted by only quadratic cells.
 *  \return bool - \c true if there are only quadratic cells in \a this mesh.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 */
bool MEDCouplingUMesh::isFullyQuadratic() const
{
  checkFullyDefined();
  bool ret=true;
  mcIdType nbOfCells=getNumberOfCells();
  for(mcIdType i=0;i<nbOfCells && ret;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      ret=cm.isQuadratic();
    }
  return ret;
}

/*!
 * Checks if \a this mesh includes any quadratic cell.
 *  \return bool - \c true if there is at least one quadratic cells in \a this mesh.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 */
bool MEDCouplingUMesh::isPresenceOfQuadratic() const
{
  checkFullyDefined();
  bool ret=false;
  mcIdType nbOfCells=getNumberOfCells();
  for(mcIdType i=0;i<nbOfCells && !ret;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      ret=cm.isQuadratic();
    }
  return ret;
}

/*!
 * Converts all quadratic cells to linear ones. If there are no quadratic cells in \a
 * this mesh, it remains unchanged.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 */
void MEDCouplingUMesh::convertQuadraticCellsToLinear()
{
  checkFullyDefined();
  mcIdType nbOfCells=getNumberOfCells();
  mcIdType delta=0;
  const mcIdType *iciptr=_nodal_connec_index->begin();
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      if(cm.isQuadratic())
        {
          INTERP_KERNEL::NormalizedCellType typel=cm.getLinearType();
          const INTERP_KERNEL::CellModel& cml=INTERP_KERNEL::CellModel::GetCellModel(typel);
          if(!cml.isDynamic())
            delta+=cm.getNumberOfNodes()-cml.getNumberOfNodes();
          else
            delta+=(iciptr[i+1]-iciptr[i]-1)/2;
        }
    }
  if(delta==0)
    return ;
  MCAuto<DataArrayIdType> newConn(DataArrayIdType::New()),newConnI(DataArrayIdType::New());
  const mcIdType *icptr(_nodal_connec->begin());
  newConn->alloc(getNodalConnectivityArrayLen()-delta,1);
  newConnI->alloc(nbOfCells+1,1);
  mcIdType *ocptr(newConn->getPointer()),*ociptr(newConnI->getPointer());
  *ociptr=0;
  _types.clear();
  for(mcIdType i=0;i<nbOfCells;i++,ociptr++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)icptr[iciptr[i]];
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      if(!cm.isQuadratic())
        {
          _types.insert(type);
          ocptr=std::copy(icptr+iciptr[i],icptr+iciptr[i+1],ocptr);
          ociptr[1]=ociptr[0]+iciptr[i+1]-iciptr[i];
        }
      else
        {
          INTERP_KERNEL::NormalizedCellType typel=cm.getLinearType();
          _types.insert(typel);
          const INTERP_KERNEL::CellModel& cml=INTERP_KERNEL::CellModel::GetCellModel(typel);
          mcIdType newNbOfNodes=cml.getNumberOfNodes();
          if(cml.isDynamic())
            newNbOfNodes=(iciptr[i+1]-iciptr[i]-1)/2;
          *ocptr++=ToIdType(typel);
          ocptr=std::copy(icptr+iciptr[i]+1,icptr+iciptr[i]+newNbOfNodes+1,ocptr);
          ociptr[1]=ociptr[0]+newNbOfNodes+1;
        }
    }
  setConnectivity(newConn,newConnI,false);
}

/*!
 * This method converts all linear cell in \a this to quadratic one.
 * Contrary to MEDCouplingUMesh::convertQuadraticCellsToLinear method, here it is needed to specify the target
 * type of cells expected. For example INTERP_KERNEL::NORM_TRI3 can be converted to INTERP_KERNEL::NORM_TRI6 if \a conversionType is equal to 0 (the default)
 * or to INTERP_KERNEL::NORM_TRI7 if \a conversionType is equal to 1. All non linear cells and polyhedron in \a this are let untouched.
 * Contrary to MEDCouplingUMesh::convertQuadraticCellsToLinear method, the coordinates in \a this can be become bigger. All created nodes will be put at the
 * end of the existing coordinates.
 *
 * \param [in] conversionType specifies the type of conversion expected. Only 0 (default) and 1 are supported presently. 0 those that creates the 'most' simple
 *             corresponding quadratic cells. 1 is those creating the 'most' complex.
 * \return a newly created DataArrayIdType instance that the caller should deal with containing cell ids of converted cells.
 *
 * \throw if \a this is not fully defined. It throws too if \a conversionType is not in [0,1].
 *
 * \sa MEDCouplingUMesh::convertQuadraticCellsToLinear
 */
DataArrayIdType *MEDCouplingUMesh::convertLinearCellsToQuadratic(int conversionType)
{
  DataArrayIdType *conn=0,*connI=0;
  DataArrayDouble *coords=0;
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  checkFullyDefined();
  MCAuto<DataArrayIdType> ret,connSafe,connISafe;
  MCAuto<DataArrayDouble> coordsSafe;
  int meshDim=getMeshDimension();
  switch(conversionType)
  {
    case 0:
      switch(meshDim)
      {
        case 1:
          ret=convertLinearCellsToQuadratic1D0(conn,connI,coords,types);
          connSafe=conn; connISafe=connI; coordsSafe=coords;
          break;
        case 2:
          ret=convertLinearCellsToQuadratic2D0(conn,connI,coords,types);
          connSafe=conn; connISafe=connI; coordsSafe=coords;
          break;
        case 3:
          ret=convertLinearCellsToQuadratic3D0(conn,connI,coords,types);
          connSafe=conn; connISafe=connI; coordsSafe=coords;
          break;
        default:
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::convertLinearCellsToQuadratic : conversion of type 0 mesh dimensions available are [1,2,3] !");
      }
      break;
        case 1:
          {
            switch(meshDim)
            {
              case 1:
                ret=convertLinearCellsToQuadratic1D0(conn,connI,coords,types);//it is not a bug. In 1D policy 0 and 1 are equals
                connSafe=conn; connISafe=connI; coordsSafe=coords;
                break;
              case 2:
                ret=convertLinearCellsToQuadratic2D1(conn,connI,coords,types);
                connSafe=conn; connISafe=connI; coordsSafe=coords;
                break;
              case 3:
                ret=convertLinearCellsToQuadratic3D1(conn,connI,coords,types);
                connSafe=conn; connISafe=connI; coordsSafe=coords;
                break;
              default:
                throw INTERP_KERNEL::Exception("MEDCouplingUMesh::convertLinearCellsToQuadratic : conversion of type 1 mesh dimensions available are [1,2,3] !");
            }
            break;
          }
        default:
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::convertLinearCellsToQuadratic : conversion type available are 0 (default, the simplest) and 1 (the most complex) !");
  }
  setConnectivity(connSafe,connISafe,false);
  _types=types;
  setCoords(coordsSafe);
  return ret.retn();
}

/*!
 * Tessellates \a this 2D mesh by dividing not straight edges of quadratic faces,
 * so that the number of cells remains the same. Quadratic faces are converted to
 * polygons. This method works only for 2D meshes in
 * 2D space. If no cells are quadratic (INTERP_KERNEL::NORM_QUAD8,
 * INTERP_KERNEL::NORM_TRI6, INTERP_KERNEL::NORM_QPOLYG ), \a this mesh remains unchanged.
 * \warning This method can lead to a huge amount of nodes if \a eps is very low.
 *  \param [in] eps - specifies the maximal angle (in radians) between 2 sub-edges of
 *         a polylinized edge constituting the input polygon.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If \a this->getMeshDimension() != 2.
 *  \throw If \a this->getSpaceDimension() != 2.
 */
void MEDCouplingUMesh::tessellate2D(double eps)
{
  int meshDim(getMeshDimension()),spaceDim(getSpaceDimension());
  if(spaceDim!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::tessellate2D : works only with space dimension equal to 2 !");
  switch(meshDim)
    {
    case 1:
      return tessellate2DCurveInternal(eps);
    case 2:
      return tessellate2DInternal(eps);
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingUMesh::tessellate2D : mesh dimension must be in [1,2] !");
    }
}

#if 0
/*!
 * This method only works if \a this has spaceDimension equal to 2 and meshDimension also equal to 2.
 * This method allows to modify connectivity of cells in \a this that shares some edges in \a edgeIdsToBeSplit.
 * The nodes to be added in those 2D cells are defined by the pair of \a  nodeIdsToAdd and \a nodeIdsIndexToAdd.
 * Length of \a nodeIdsIndexToAdd is expected to equal to length of \a edgeIdsToBeSplit + 1.
 * The node ids in \a nodeIdsToAdd should be valid. Those nodes have to be sorted exactly following exactly the direction of the edge.
 * This method can be seen as the opposite method of colinearize2D.
 * This method can be lead to create some new nodes if quadratic polygon cells have to be split. In this case the added nodes will be put at the end
 * to avoid to modify the numbering of existing nodes.
 *
 * \param [in] nodeIdsToAdd - the list of node ids to be added (\a nodeIdsIndexToAdd array allows to walk on this array)
 * \param [in] nodeIdsIndexToAdd - the entry point of \a nodeIdsToAdd to point to the corresponding nodes to be added.
 * \param [in] mesh1Desc - 1st output of buildDescendingConnectivity2 on \a this.
 * \param [in] desc - 2nd output of buildDescendingConnectivity2 on \a this.
 * \param [in] descI - 3rd output of buildDescendingConnectivity2 on \a this.
 * \param [in] revDesc - 4th output of buildDescendingConnectivity2 on \a this.
 * \param [in] revDescI - 5th output of buildDescendingConnectivity2 on \a this.
 *
 * \sa buildDescendingConnectivity2
 */
void MEDCouplingUMesh::splitSomeEdgesOf2DMesh(const DataArrayIdType *nodeIdsToAdd, const DataArrayIdType *nodeIdsIndexToAdd, const DataArrayIdType *edgeIdsToBeSplit,
                                              const MEDCouplingUMesh *mesh1Desc, const DataArrayIdType *desc, const DataArrayIdType *descI, const DataArrayIdType *revDesc, const DataArrayIdType *revDescI)
{
  if(!nodeIdsToAdd || !nodeIdsIndexToAdd || !edgeIdsToBeSplit || !mesh1Desc || !desc || !descI || !revDesc || !revDescI)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::splitSomeEdgesOf2DMesh : input pointers must be not NULL !");
  nodeIdsToAdd->checkAllocated(); nodeIdsIndexToAdd->checkAllocated(); edgeIdsToBeSplit->checkAllocated(); desc->checkAllocated(); descI->checkAllocated(); revDesc->checkAllocated(); revDescI->checkAllocated();
  if(getSpaceDimension()!=2 || getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::splitSomeEdgesOf2DMesh : this must have spacedim=meshdim=2 !");
  if(mesh1Desc->getSpaceDimension()!=2 || mesh1Desc->getMeshDimension()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::splitSomeEdgesOf2DMesh : mesh1Desc must be the explosion of this with spaceDim=2 and meshDim = 1 !");
  //DataArrayIdType *out0(0),*outi0(0);
  //MEDCouplingUMesh::ExtractFromIndexedArrays(idsInDesc2DToBeRefined->begin(),idsInDesc2DToBeRefined->end(),dd3,dd4,out0,outi0);
  //MCAuto<DataArrayIdType> out0s(out0),outi0s(outi0);
  //out0s=out0s->buildUnique(); out0s->sort(true);
}
#endif


/*!
 * Divides every cell of \a this mesh into simplices (triangles in 2D and tetrahedra in 3D).
 * In addition, returns an array mapping new cells to old ones. <br>
 * This method typically increases the number of cells in \a this mesh
 * but the number of nodes remains \b unchanged.
 * That's why the 3D splitting policies
 * INTERP_KERNEL::GENERAL_24 and INTERP_KERNEL::GENERAL_48 are not available here.
 *  \param [in] policy - specifies a pattern used for splitting.
 * The semantic of \a policy is:
 * - 0 - to split QUAD4 by cutting it along 0-2 diagonal (for 2D mesh only).
 * - 1 - to split QUAD4 by cutting it along 1-3 diagonal (for 2D mesh only).
 * - INTERP_KERNEL::PLANAR_FACE_5 - to split HEXA8  into 5 TETRA4 (for 3D mesh only - see INTERP_KERNEL::SplittingPolicy for an image).
 * - INTERP_KERNEL::PLANAR_FACE_6 - to split HEXA8  into 6 TETRA4 (for 3D mesh only - see INTERP_KERNEL::SplittingPolicy for an image).
 *
 *
 *  \return DataArrayIdType * - a new instance of DataArrayIdType holding, for each new cell,
 *          an id of old cell producing it. The caller is to delete this array using
 *         decrRef() as it is no more needed.
 *
 *  \throw If \a policy is 0 or 1 and \a this->getMeshDimension() != 2.
 *  \throw If \a policy is INTERP_KERNEL::PLANAR_FACE_5 or INTERP_KERNEL::PLANAR_FACE_6
 *          and \a this->getMeshDimension() != 3.
 *  \throw If \a policy is not one of the four discussed above.
 *  \throw If the nodal connectivity of cells is not defined.
 * \sa MEDCouplingUMesh::tetrahedrize, MEDCoupling1SGTUMesh::sortHexa8EachOther
 */
DataArrayIdType *MEDCouplingUMesh::simplexize(int policy)
{
  switch(policy)
  {
    case 0:
      return simplexizePol0();
    case 1:
      return simplexizePol1();
    case INTERP_KERNEL::PLANAR_FACE_5:
        return simplexizePlanarFace5();
    case INTERP_KERNEL::PLANAR_FACE_6:
        return simplexizePlanarFace6();
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingUMesh::simplexize : unrecognized policy ! Must be :\n  - 0 or 1 (only available for meshdim=2) \n  - PLANAR_FACE_5, PLANAR_FACE_6  (only for meshdim=3)");
  }
}

/*!
 * Checks if \a this mesh is constituted by simplex cells only. Simplex cells are:
 * - 1D: INTERP_KERNEL::NORM_SEG2
 * - 2D: INTERP_KERNEL::NORM_TRI3
 * - 3D: INTERP_KERNEL::NORM_TETRA4.
 *
 * This method is useful for users that need to use P1 field services as
 * MEDCouplingFieldDouble::getValueOn(), MEDCouplingField::buildMeasureField() etc.
 * All these methods need mesh support containing only simplex cells.
 *  \return bool - \c true if there are only simplex cells in \a this mesh.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If \a this->getMeshDimension() < 1.
 */
bool MEDCouplingUMesh::areOnlySimplexCells() const
{
  checkFullyDefined();
  int mdim=getMeshDimension();
  if(mdim<1 || mdim>3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::areOnlySimplexCells : only available with meshes having a meshdim 1, 2 or 3 !");
  mcIdType nbCells=getNumberOfCells();
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  for(mcIdType i=0;i<nbCells;i++)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[connI[i]]);
      if(!cm.isSimplex())
        return false;
    }
  return true;
}



/*!
 * Converts degenerated 2D or 3D linear cells of \a this mesh into cells of simpler
 * type. For example an INTERP_KERNEL::NORM_QUAD4 cell having only three unique nodes in
 * its connectivity is transformed into an INTERP_KERNEL::NORM_TRI3 cell.
 * Quadratic cells in 2D are also handled. In those cells edges where start=end=midpoint are removed.
 * This method does \b not perform geometrical checks and checks only nodal connectivity of cells,
 * so it can be useful to call mergeNodes() before calling this method.
 *  \throw If \a this->getMeshDimension() <= 1.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 */
void MEDCouplingUMesh::convertDegeneratedCells()
{
  checkFullyDefined();
  if(getMeshDimension()<=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::convertDegeneratedCells works on umeshes with meshdim equals to 2 or 3 !");
  mcIdType nbOfCells=getNumberOfCells();
  if(nbOfCells<1)
    return ;
  mcIdType initMeshLgth=getNodalConnectivityArrayLen();
  mcIdType *conn=_nodal_connec->getPointer();
  mcIdType *index=_nodal_connec_index->getPointer();
  mcIdType posOfCurCell=0;
  mcIdType newPos=0;
  mcIdType lgthOfCurCell;
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      lgthOfCurCell=index[i+1]-posOfCurCell;
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)conn[posOfCurCell];
      mcIdType newLgth;
      INTERP_KERNEL::NormalizedCellType newType=INTERP_KERNEL::CellSimplify::simplifyDegeneratedCell(type,conn+posOfCurCell+1,lgthOfCurCell-1,
                                                                                                     conn+newPos+1,newLgth);
      conn[newPos]=newType;
      newPos+=newLgth+1;
      posOfCurCell=index[i+1];
      index[i+1]=newPos;
    }
  if(newPos!=initMeshLgth)
    _nodal_connec->reAlloc(newPos);
  computeTypes();
}

/*!
 * Same as MEDCouplingUMesh::convertDegeneratedCells() plus deletion of the flat cells.
 * A cell is flat in the following cases:
 *   - for a linear cell, all points in the connectivity are equal
 *   - for a quadratic cell, either the above, or a quadratic polygon with two (linear) points and two
 *   identical quadratic points
 * \return a new instance of DataArrayIdType holding ids of removed cells. The caller is to delete
 *      this array using decrRef() as it is no more needed.
 */
DataArrayIdType *MEDCouplingUMesh::convertDegeneratedCellsAndRemoveFlatOnes()
{
  checkFullyDefined();
  if(getMeshDimension()<=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::convertDegeneratedCells works on umeshes with meshdim equals to 2 or 3 !");
  mcIdType nbOfCells=getNumberOfCells();
  MCAuto<DataArrayIdType> ret(DataArrayIdType::New()); ret->alloc(0,1);
  if(nbOfCells<1)
    return ret.retn();
  mcIdType initMeshLgth=getNodalConnectivityArrayLen();
  mcIdType *conn=_nodal_connec->getPointer();
  mcIdType *index=_nodal_connec_index->getPointer();
  mcIdType posOfCurCell=0;
  mcIdType newPos=0;
  mcIdType lgthOfCurCell, nbDelCells(0);
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      lgthOfCurCell=index[i+1]-posOfCurCell;
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)conn[posOfCurCell];
      mcIdType newLgth;
      INTERP_KERNEL::NormalizedCellType newType=INTERP_KERNEL::CellSimplify::simplifyDegeneratedCell(type,conn+posOfCurCell+1,lgthOfCurCell-1,
                                                                                                     conn+newPos+1,newLgth);
      // Shall we delete the cell if it is completely degenerated:
      bool delCell=INTERP_KERNEL::CellSimplify::isFlatCell(conn, newPos, newLgth, newType);
      if (delCell)
        {
          nbDelCells++;
          ret->pushBackSilent(i);
        }
      else   //if the cell is to be deleted, simply stay at the same place
        {
          conn[newPos]=newType;
          newPos+=newLgth+1;
        }
      posOfCurCell=index[i+1];
      index[i+1-nbDelCells]=newPos;
    }
  if(newPos!=initMeshLgth)
    _nodal_connec->reAlloc(newPos);
  const mcIdType nCellDel=ret->getNumberOfTuples();
  if (nCellDel)
    _nodal_connec_index->reAlloc(nbOfCells-nCellDel+1);
  computeTypes();
  return ret.retn();
}

/*!
 * This method remove null 1D cells from \a this. A 1D cell is considered null if start node is equal to end node.
 * Only connectivity is considered here.
 */
bool MEDCouplingUMesh::removeDegenerated1DCells()
{
  checkConnectivityFullyDefined();
  if(getMeshDimension()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::removeDegenerated1DCells works on umeshes with meshdim equals to 1 !");
  std::size_t nbCells(getNumberOfCells()),newSize(0),newSize2(0);
  const mcIdType *conn(getNodalConnectivity()->begin()),*conni(getNodalConnectivityIndex()->begin());
  {
    for(std::size_t i=0;i<nbCells;i++)
      {
        INTERP_KERNEL::NormalizedCellType ct((INTERP_KERNEL::NormalizedCellType)conn[conni[i]]);
        if(ct==INTERP_KERNEL::NORM_SEG2 || ct==INTERP_KERNEL::NORM_SEG3)
          {
            if(conn[conni[i]+1]!=conn[conni[i]+2])
              {
                newSize++;
                newSize2+=conni[i+1]-conni[i];
              }
          }
        else
          {
            std::ostringstream oss; oss << "MEDCouplingUMesh::removeDegenerated1DCells : cell #" << i << " in this is not of type SEG2/SEG3 !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
  }
  if(newSize==nbCells)//no cells has been removed -> do nothing
    return false;
  MCAuto<DataArrayIdType> newConn(DataArrayIdType::New()),newConnI(DataArrayIdType::New()); newConnI->alloc(newSize+1,1); newConn->alloc(newSize2,1);
  mcIdType *newConnPtr(newConn->getPointer()),*newConnIPtr(newConnI->getPointer()); newConnIPtr[0]=0;
  for(std::size_t i=0;i<nbCells;i++)
    {
      if(conn[conni[i]+1]!=conn[conni[i]+2])
        {
          newConnIPtr[1]=newConnIPtr[0]+conni[i+1]-conni[i];
          newConnPtr=std::copy(conn+conni[i],conn+conni[i+1],newConnPtr);
          newConnIPtr++;
        }
    }
  setConnectivity(newConn,newConnI,true);
  return true;
}

/*!
 * Finds incorrectly oriented cells of this 2D mesh in 3D space.
 * A cell is considered to be oriented correctly if an angle between its
 * normal vector and a given vector is less than \c PI / \c 2.
 *  \param [in] vec - 3 components of the vector specifying the correct orientation of
 *         cells.
 *  \param [in] polyOnly - if \c true, only polygons are checked, else, all cells are
 *         checked.
 *  \param [in,out] cells - a vector returning ids of incorrectly oriented cells. It
 *         is not cleared before filling in.
 *  \throw If \a this->getMeshDimension() != 2.
 *  \throw If \a this->getSpaceDimension() != 3.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_are2DCellsNotCorrectlyOriented "Here is a C++ example".<br>
 *  \ref  py_mcumesh_are2DCellsNotCorrectlyOriented "Here is a Python example".
 *  \endif
 */
void MEDCouplingUMesh::are2DCellsNotCorrectlyOriented(const double *vec, bool polyOnly, std::vector<mcIdType>& cells) const
{
  if(getMeshDimension()!=2 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply are2DCellsNotCorrectlyOriented on it : must be meshDim==2 and spaceDim==3 !");
  mcIdType nbOfCells=getNumberOfCells();
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  const double *coordsPtr=_coords->begin();
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)conn[connI[i]];
      if(!polyOnly || (type==INTERP_KERNEL::NORM_POLYGON || type==INTERP_KERNEL::NORM_QPOLYG))
        {
          bool isQuadratic=INTERP_KERNEL::CellModel::GetCellModel(type).isQuadratic();
          if(!IsPolygonWellOriented(isQuadratic,vec,conn+connI[i]+1,conn+connI[i+1],coordsPtr))
            cells.push_back(i);
        }
    }
}

/*!
 * Reverse connectivity of 2D cells whose orientation is not correct. A cell is
 * considered to be oriented correctly if an angle between its normal vector and a
 * given vector is less than \c PI / \c 2.
 *  \param [in] vec - 3 components of the vector specifying the correct orientation of
 *         cells.
 *  \param [in] polyOnly - if \c true, only polygons are checked, else, all cells are
 *         checked.
 *  \throw If \a this->getMeshDimension() != 2.
 *  \throw If \a this->getSpaceDimension() != 3.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_are2DCellsNotCorrectlyOriented "Here is a C++ example".<br>
 *  \ref  py_mcumesh_are2DCellsNotCorrectlyOriented "Here is a Python example".
 *  \endif
 *
 *  \sa changeOrientationOfCells
 */
void MEDCouplingUMesh::orientCorrectly2DCells(const double *vec, bool polyOnly)
{
  if(getMeshDimension()!=2 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply orientCorrectly2DCells on it : must be meshDim==2 and spaceDim==3 !");
  mcIdType nbOfCells=getNumberOfCells();
  mcIdType *conn(_nodal_connec->getPointer());
  const mcIdType *connI(_nodal_connec_index->begin());
  const double *coordsPtr(_coords->begin());
  bool isModified(false);
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type((INTERP_KERNEL::NormalizedCellType)conn[connI[i]]);
      if(!polyOnly || (type==INTERP_KERNEL::NORM_POLYGON || type==INTERP_KERNEL::NORM_QPOLYG))
        {
          const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(type));
          bool isQuadratic(cm.isQuadratic());
          if(!IsPolygonWellOriented(isQuadratic,vec,conn+connI[i]+1,conn+connI[i+1],coordsPtr))
            {
              isModified=true;
              cm.changeOrientationOf2D(conn+connI[i]+1,(unsigned int)(connI[i+1]-connI[i]-1));
            }
        }
    }
  if(isModified)
    _nodal_connec->declareAsNew();
  updateTime();
}

/*!
 * This method change the orientation of cells in \a this without any consideration of coordinates. Only connectivity is impacted.
 *
 * \sa orientCorrectly2DCells
 */
void MEDCouplingUMesh::changeOrientationOfCells()
{
  int mdim(getMeshDimension());
  if(mdim!=2 && mdim!=1)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply changeOrientationOfCells on it : must be meshDim==2 or meshDim==1 !");
  mcIdType nbOfCells=getNumberOfCells();
  mcIdType *conn(_nodal_connec->getPointer());
  const mcIdType *connI(_nodal_connec_index->begin());
  if(mdim==2)
    {//2D
      for(mcIdType i=0;i<nbOfCells;i++)
        {
          INTERP_KERNEL::NormalizedCellType type((INTERP_KERNEL::NormalizedCellType)conn[connI[i]]);
          const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(type));
          cm.changeOrientationOf2D(conn+connI[i]+1,(unsigned int)(connI[i+1]-connI[i]-1));
        }
    }
  else
    {//1D
      for(mcIdType i=0;i<nbOfCells;i++)
        {
          INTERP_KERNEL::NormalizedCellType type((INTERP_KERNEL::NormalizedCellType)conn[connI[i]]);
          const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(type));
          cm.changeOrientationOf1D(conn+connI[i]+1,(unsigned int)(connI[i+1]-connI[i]-1));
        }
    }
}

/*!
 * Finds incorrectly oriented polyhedral cells, i.e. polyhedrons having correctly
 * oriented facets. The normal vector of the facet should point out of the cell.
 *  \param [in,out] cells - a vector returning ids of incorrectly oriented cells. It
 *         is not cleared before filling in.
 *  \throw If \a this->getMeshDimension() != 3.
 *  \throw If \a this->getSpaceDimension() != 3.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_arePolyhedronsNotCorrectlyOriented "Here is a C++ example".<br>
 *  \ref  py_mcumesh_arePolyhedronsNotCorrectlyOriented "Here is a Python example".
 *  \endif
 */
void MEDCouplingUMesh::arePolyhedronsNotCorrectlyOriented(std::vector<mcIdType>& cells) const
{
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply arePolyhedronsNotCorrectlyOriented on it : must be meshDim==3 and spaceDim==3 !");
  mcIdType nbOfCells=getNumberOfCells();
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  const double *coordsPtr=_coords->begin();
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)conn[connI[i]];
      if(type==INTERP_KERNEL::NORM_POLYHED)
        {
          if(!IsPolyhedronWellOriented(conn+connI[i]+1,conn+connI[i+1],coordsPtr))
            cells.push_back(i);
        }
    }
}

/*!
 * Tries to fix connectivity of polyhedra, so that normal vector of all facets to point
 * out of the cell.
 *  \throw If \a this->getMeshDimension() != 3.
 *  \throw If \a this->getSpaceDimension() != 3.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If the reparation fails.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_arePolyhedronsNotCorrectlyOriented "Here is a C++ example".<br>
 *  \ref  py_mcumesh_arePolyhedronsNotCorrectlyOriented "Here is a Python example".
 *  \endif
 * \sa MEDCouplingUMesh::findAndCorrectBadOriented3DCells
 */
void MEDCouplingUMesh::orientCorrectlyPolyhedrons()
{
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply orientCorrectlyPolyhedrons on it : must be meshDim==3 and spaceDim==3 !");
  mcIdType nbOfCells=getNumberOfCells();
  mcIdType *conn=_nodal_connec->getPointer();
  const mcIdType *connI=_nodal_connec_index->begin();
  const double *coordsPtr=_coords->begin();
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)conn[connI[i]];
      if(type==INTERP_KERNEL::NORM_POLYHED)
        {
          try
          {
              if(!IsPolyhedronWellOriented(conn+connI[i]+1,conn+connI[i+1],coordsPtr))
                TryToCorrectPolyhedronOrientation(conn+connI[i]+1,conn+connI[i+1],coordsPtr);
          }
          catch(INTERP_KERNEL::Exception& e)
          {
              std::ostringstream oss; oss << "Something wrong in polyhedron #" << i << " : " << e.what();
              throw INTERP_KERNEL::Exception(oss.str());
          }
        }
    }
  updateTime();
}

/*!
 * This method invert orientation of all cells in \a this.
 * After calling this method the absolute value of measure of cells in \a this are the same than before calling.
 * This method only operates on the connectivity so coordinates are not touched at all.
 */
void MEDCouplingUMesh::invertOrientationOfAllCells()
{
  checkConnectivityFullyDefined();
  std::set<INTERP_KERNEL::NormalizedCellType> gts(getAllGeoTypes());
  mcIdType *conn(_nodal_connec->getPointer());
  const mcIdType *conni(_nodal_connec_index->begin());
  for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator gt=gts.begin();gt!=gts.end();gt++)
    {
      INTERP_KERNEL::AutoCppPtr<INTERP_KERNEL::OrientationInverter> oi(INTERP_KERNEL::OrientationInverter::BuildInstanceFrom(*gt));
      MCAuto<DataArrayIdType> cwt(giveCellsWithType(*gt));
      for(const mcIdType *it=cwt->begin();it!=cwt->end();it++)
        oi->operate(conn+conni[*it]+1,conn+conni[*it+1]);
    }
  updateTime();
}

/*!
 * Finds and fixes incorrectly oriented linear extruded volumes (INTERP_KERNEL::NORM_HEXA8,
 * INTERP_KERNEL::NORM_PENTA6, INTERP_KERNEL::NORM_HEXGP12 etc) to respect the MED convention
 * according to which the first facet of the cell should be oriented to have the normal vector
 * pointing out of cell.
 *  \return DataArrayIdType * - a new instance of DataArrayIdType holding ids of fixed
 *         cells. The caller is to delete this array using decrRef() as it is no more
 *         needed.
 *  \throw If \a this->getMeshDimension() != 3.
 *  \throw If \a this->getSpaceDimension() != 3.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_findAndCorrectBadOriented3DExtrudedCells "Here is a C++ example".<br>
 *  \ref  py_mcumesh_findAndCorrectBadOriented3DExtrudedCells "Here is a Python example".
 *  \endif
 * \sa MEDCouplingUMesh::findAndCorrectBadOriented3DCells
 */
DataArrayIdType *MEDCouplingUMesh::findAndCorrectBadOriented3DExtrudedCells()
{
  const char msg[]="check3DCellsWellOriented detection works only for 3D cells !";
  if(getMeshDimension()!=3)
    throw INTERP_KERNEL::Exception(msg);
  int spaceDim=getSpaceDimension();
  if(spaceDim!=3)
    throw INTERP_KERNEL::Exception(msg);
  //
  mcIdType nbOfCells=getNumberOfCells();
  mcIdType *conn=_nodal_connec->getPointer();
  const mcIdType *connI=_nodal_connec_index->begin();
  const double *coo=getCoords()->begin();
  MCAuto<DataArrayIdType> cells(DataArrayIdType::New()); cells->alloc(0,1);
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[connI[i]]);
      if(cm.isExtruded() && !cm.isDynamic() && !cm.isQuadratic())
        {
          if(!Is3DExtrudedStaticCellWellOriented(conn+connI[i]+1,conn+connI[i+1],coo))
            {
              CorrectExtrudedStaticCell(conn+connI[i]+1,conn+connI[i+1]);
              cells->pushBackSilent(i);
            }
        }
    }
  return cells.retn();
}

/*!
 * This method is a faster method to correct orientation of all 3D cells in \a this.
 * This method works only if \a this is a 3D mesh, that is to say a mesh with mesh dimension 3 and a space dimension 3.
 * This method makes the hypothesis that \a this a coherent that is to say MEDCouplingUMesh::checkConsistency should throw no exception.
 *
 * \return a newly allocated mcIdType array with one components containing cell ids renumbered to fit the convention of MED (MED file and MEDCoupling)
 * \sa MEDCouplingUMesh::orientCorrectlyPolyhedrons,
 */
DataArrayIdType *MEDCouplingUMesh::findAndCorrectBadOriented3DCells()
{
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply findAndCorrectBadOriented3DCells on it : must be meshDim==3 and spaceDim==3 !");
  mcIdType nbOfCells=getNumberOfCells();
  mcIdType *conn=_nodal_connec->getPointer();
  const mcIdType *connI=_nodal_connec_index->begin();
  const double *coordsPtr=_coords->begin();
  MCAuto<DataArrayIdType> ret=DataArrayIdType::New(); ret->alloc(0,1);
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)conn[connI[i]];
      switch(type)
      {
        case INTERP_KERNEL::NORM_TETRA4:
          {
            if(!IsTetra4WellOriented(conn+connI[i]+1,conn+connI[i+1],coordsPtr))
              {
                std::swap(*(conn+connI[i]+2),*(conn+connI[i]+3));
                ret->pushBackSilent(i);
              }
            break;
          }
        case INTERP_KERNEL::NORM_PYRA5:
          {
            if(!IsPyra5WellOriented(conn+connI[i]+1,conn+connI[i+1],coordsPtr))
              {
                std::swap(*(conn+connI[i]+2),*(conn+connI[i]+4));
                ret->pushBackSilent(i);
              }
            break;
          }
        case INTERP_KERNEL::NORM_PENTA6:
        case INTERP_KERNEL::NORM_HEXA8:
        case INTERP_KERNEL::NORM_HEXGP12:
          {
            if(!Is3DExtrudedStaticCellWellOriented(conn+connI[i]+1,conn+connI[i+1],coordsPtr))
              {
                CorrectExtrudedStaticCell(conn+connI[i]+1,conn+connI[i+1]);
                ret->pushBackSilent(i);
              }
            break;
          }
        case INTERP_KERNEL::NORM_POLYHED:
          {
            if(!IsPolyhedronWellOriented(conn+connI[i]+1,conn+connI[i+1],coordsPtr))
              {
                TryToCorrectPolyhedronOrientation(conn+connI[i]+1,conn+connI[i+1],coordsPtr);
                ret->pushBackSilent(i);
              }
            break;
          }
        default:
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::orientCorrectly3DCells : Your mesh contains type of cell not supported yet ! send mail to anthony.geay@cea.fr to add it !");
      }
    }
  updateTime();
  return ret.retn();
}

/*!
 * This method has a sense for meshes with spaceDim==3 and meshDim==2.
 * If it is not the case an exception will be thrown.
 * This method is fast because the first cell of \a this is used to compute the plane.
 * \param vec output of size at least 3 used to store the normal vector (with norm equal to Area ) of searched plane.
 * \param pos output of size at least 3 used to store a point owned of searched plane.
 */
void MEDCouplingUMesh::getFastAveragePlaneOfThis(double *vec, double *pos) const
{
  if(getMeshDimension()!=2 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply getFastAveragePlaneOfThis on it : must be meshDim==2 and spaceDim==3 !");
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  const double *coordsPtr=_coords->begin();
  INTERP_KERNEL::areaVectorOfPolygon<mcIdType,INTERP_KERNEL::ALL_C_MODE>(conn+1,connI[1]-connI[0]-1,coordsPtr,vec);
  std::copy(coordsPtr+3*conn[1],coordsPtr+3*conn[1]+3,pos);
}

/*!
 * Creates a new MEDCouplingFieldDouble holding Edge Ratio values of all
 * cells. Currently cells of the following types are treated:
 * INTERP_KERNEL::NORM_TRI3, INTERP_KERNEL::NORM_QUAD4 and INTERP_KERNEL::NORM_TETRA4.
 * For a cell of other type an exception is thrown.
 * Space dimension of a 2D mesh can be either 2 or 3.
 * The Edge Ratio of a cell \f$t\f$ is:
 *  \f$\frac{|t|_\infty}{|t|_0}\f$,
 *  where \f$|t|_\infty\f$ and \f$|t|_0\f$ respectively denote the greatest and
 *  the smallest edge lengths of \f$t\f$.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble on
 *          cells and one time, lying on \a this mesh. The caller is to delete this
 *          field using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If \a this mesh contains elements of dimension different from the mesh dimension.
 *  \throw If the connectivity data array has more than one component.
 *  \throw If the connectivity data array has a named component.
 *  \throw If the connectivity index data array has more than one component.
 *  \throw If the connectivity index data array has a named component.
 *  \throw If \a this->getMeshDimension() is neither 2 nor 3.
 *  \throw If \a this->getSpaceDimension() is neither 2 nor 3.
 *  \throw If \a this mesh includes cells of type different from the ones enumerated above.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::getEdgeRatioField() const
{
  checkConsistencyLight();
  int spaceDim=getSpaceDimension();
  int meshDim=getMeshDimension();
  if(spaceDim!=2 && spaceDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getEdgeRatioField : SpaceDimension must be equal to 2 or 3 !");
  if(meshDim!=2 && meshDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getEdgeRatioField : MeshDimension must be equal to 2 or 3 !");
  MCAuto<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  ret->setMesh(this);
  mcIdType nbOfCells=getNumberOfCells();
  MCAuto<DataArrayDouble> arr=DataArrayDouble::New();
  arr->alloc(nbOfCells,1);
  double *pt=arr->getPointer();
  ret->setArray(arr);//In case of throw to avoid mem leaks arr will be used after decrRef.
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  const double *coo=_coords->begin();
  double tmp[12];
  for(mcIdType i=0;i<nbOfCells;i++,pt++)
    {
      INTERP_KERNEL::NormalizedCellType t=(INTERP_KERNEL::NormalizedCellType)*conn;
      switch(t)
      {
        case INTERP_KERNEL::NORM_TRI3:
          {
            FillInCompact3DMode(spaceDim,3,conn+1,coo,tmp);
            *pt=INTERP_KERNEL::triEdgeRatio(tmp);
            break;
          }
        case INTERP_KERNEL::NORM_QUAD4:
          {
            FillInCompact3DMode(spaceDim,4,conn+1,coo,tmp);
            *pt=INTERP_KERNEL::quadEdgeRatio(tmp);
            break;
          }
        case INTERP_KERNEL::NORM_TETRA4:
          {
            FillInCompact3DMode(spaceDim,4,conn+1,coo,tmp);
            *pt=INTERP_KERNEL::tetraEdgeRatio(tmp);
            break;
          }
        default:
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getEdgeRatioField : A cell with not manged type (NORM_TRI3, NORM_QUAD4 and NORM_TETRA4) has been detected !");
      }
      conn+=connI[i+1]-connI[i];
    }
  ret->setName("EdgeRatio");
  ret->synchronizeTimeWithSupport();
  return ret.retn();
}

/*!
 * Creates a new MEDCouplingFieldDouble holding Aspect Ratio values of all
 * cells. Currently cells of the following types are treated:
 * INTERP_KERNEL::NORM_TRI3, INTERP_KERNEL::NORM_QUAD4 and INTERP_KERNEL::NORM_TETRA4.
 * For a cell of other type an exception is thrown.
 * Space dimension of a 2D mesh can be either 2 or 3.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble on
 *          cells and one time, lying on \a this mesh. The caller is to delete this
 *          field using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If \a this mesh contains elements of dimension different from the mesh dimension.
 *  \throw If the connectivity data array has more than one component.
 *  \throw If the connectivity data array has a named component.
 *  \throw If the connectivity index data array has more than one component.
 *  \throw If the connectivity index data array has a named component.
 *  \throw If \a this->getMeshDimension() is neither 2 nor 3.
 *  \throw If \a this->getSpaceDimension() is neither 2 nor 3.
 *  \throw If \a this mesh includes cells of type different from the ones enumerated above.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::getAspectRatioField() const
{
  checkConsistencyLight();
  int spaceDim=getSpaceDimension();
  int meshDim=getMeshDimension();
  if(spaceDim!=2 && spaceDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getAspectRatioField : SpaceDimension must be equal to 2 or 3 !");
  if(meshDim!=2 && meshDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getAspectRatioField : MeshDimension must be equal to 2 or 3 !");
  MCAuto<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  ret->setMesh(this);
  mcIdType nbOfCells=getNumberOfCells();
  MCAuto<DataArrayDouble> arr=DataArrayDouble::New();
  arr->alloc(nbOfCells,1);
  double *pt=arr->getPointer();
  ret->setArray(arr);//In case of throw to avoid mem leaks arr will be used after decrRef.
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  const double *coo=_coords->begin();
  double tmp[12];
  for(mcIdType i=0;i<nbOfCells;i++,pt++)
    {
      INTERP_KERNEL::NormalizedCellType t=(INTERP_KERNEL::NormalizedCellType)*conn;
      switch(t)
      {
        case INTERP_KERNEL::NORM_TRI3:
          {
            FillInCompact3DMode(spaceDim,3,conn+1,coo,tmp);
            *pt=INTERP_KERNEL::triAspectRatio(tmp);
            break;
          }
        case INTERP_KERNEL::NORM_QUAD4:
          {
            FillInCompact3DMode(spaceDim,4,conn+1,coo,tmp);
            *pt=INTERP_KERNEL::quadAspectRatio(tmp);
            break;
          }
        case INTERP_KERNEL::NORM_TETRA4:
          {
            FillInCompact3DMode(spaceDim,4,conn+1,coo,tmp);
            *pt=INTERP_KERNEL::tetraAspectRatio(tmp);
            break;
          }
        default:
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getAspectRatioField : A cell with not manged type (NORM_TRI3, NORM_QUAD4 and NORM_TETRA4) has been detected !");
      }
      conn+=connI[i+1]-connI[i];
    }
  ret->setName("AspectRatio");
  ret->synchronizeTimeWithSupport();
  return ret.retn();
}

/*!
 * Creates a new MEDCouplingFieldDouble holding Warping factor values of all
 * cells of \a this 2D mesh in 3D space. It is a measure of the "planarity" of 2D cell
 * in 3D space. Currently only cells of the following types are
 * treated: INTERP_KERNEL::NORM_QUAD4.
 * For a cell of other type an exception is thrown.
 * The warp field is computed as follows: let (a,b,c,d) be the points of the quad.
 * Defining
 * \f$t=\vec{da}\times\vec{ab}\f$,
 * \f$u=\vec{ab}\times\vec{bc}\f$
 * \f$v=\vec{bc}\times\vec{cd}\f$
 * \f$w=\vec{cd}\times\vec{da}\f$, the warp is defined as \f$W^3\f$ with
 *  \f[
 *     W=min(\frac{t}{|t|}\cdot\frac{v}{|v|}, \frac{u}{|u|}\cdot\frac{w}{|w|})
 *  \f]
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble on
 *          cells and one time, lying on \a this mesh. The caller is to delete this
 *          field using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If \a this mesh contains elements of dimension different from the mesh dimension.
 *  \throw If the connectivity data array has more than one component.
 *  \throw If the connectivity data array has a named component.
 *  \throw If the connectivity index data array has more than one component.
 *  \throw If the connectivity index data array has a named component.
 *  \throw If \a this->getMeshDimension() != 2.
 *  \throw If \a this->getSpaceDimension() != 3.
 *  \throw If \a this mesh includes cells of type different from the ones enumerated above.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::getWarpField() const
{
  checkConsistencyLight();
  int spaceDim=getSpaceDimension();
  int meshDim=getMeshDimension();
  if(spaceDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getWarpField : SpaceDimension must be equal to 3 !");
  if(meshDim!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getWarpField : MeshDimension must be equal to 2 !");
  MCAuto<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  ret->setMesh(this);
  mcIdType nbOfCells=getNumberOfCells();
  MCAuto<DataArrayDouble> arr=DataArrayDouble::New();
  arr->alloc(nbOfCells,1);
  double *pt=arr->getPointer();
  ret->setArray(arr);//In case of throw to avoid mem leaks arr will be used after decrRef.
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  const double *coo=_coords->begin();
  double tmp[12];
  for(mcIdType i=0;i<nbOfCells;i++,pt++)
    {
      INTERP_KERNEL::NormalizedCellType t=(INTERP_KERNEL::NormalizedCellType)*conn;
      switch(t)
      {
        case INTERP_KERNEL::NORM_QUAD4:
          {
            FillInCompact3DMode(3,4,conn+1,coo,tmp);
            *pt=INTERP_KERNEL::quadWarp(tmp);
            break;
          }
        default:
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getWarpField : A cell with not manged type (NORM_QUAD4) has been detected !");
      }
      conn+=connI[i+1]-connI[i];
    }
  ret->setName("Warp");
  ret->synchronizeTimeWithSupport();
  return ret.retn();
}


/*!
 * Creates a new MEDCouplingFieldDouble holding Skew factor values of all
 * cells of \a this 2D mesh in 3D space. Currently cells of the following types are
 * treated: INTERP_KERNEL::NORM_QUAD4.
 * The skew is computed as follow for a quad with points (a,b,c,d): let
 * \f$u=\vec{ab}+\vec{dc}\f$ and \f$v=\vec{ac}+\vec{bd}\f$
 * then the skew is computed as:
 *  \f[
 *    s=\frac{u}{|u|}\cdot\frac{v}{|v|}
 *  \f]
 *
 * For a cell of other type an exception is thrown.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble on
 *          cells and one time, lying on \a this mesh. The caller is to delete this
 *          field using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If \a this mesh contains elements of dimension different from the mesh dimension.
 *  \throw If the connectivity data array has more than one component.
 *  \throw If the connectivity data array has a named component.
 *  \throw If the connectivity index data array has more than one component.
 *  \throw If the connectivity index data array has a named component.
 *  \throw If \a this->getMeshDimension() != 2.
 *  \throw If \a this->getSpaceDimension() != 3.
 *  \throw If \a this mesh includes cells of type different from the ones enumerated above.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::getSkewField() const
{
  checkConsistencyLight();
  int spaceDim=getSpaceDimension();
  int meshDim=getMeshDimension();
  if(spaceDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getSkewField : SpaceDimension must be equal to 3 !");
  if(meshDim!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getSkewField : MeshDimension must be equal to 2 !");
  MCAuto<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  ret->setMesh(this);
  mcIdType nbOfCells=getNumberOfCells();
  MCAuto<DataArrayDouble> arr=DataArrayDouble::New();
  arr->alloc(nbOfCells,1);
  double *pt=arr->getPointer();
  ret->setArray(arr);//In case of throw to avoid mem leaks arr will be used after decrRef.
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  const double *coo=_coords->begin();
  double tmp[12];
  for(mcIdType i=0;i<nbOfCells;i++,pt++)
    {
      INTERP_KERNEL::NormalizedCellType t=(INTERP_KERNEL::NormalizedCellType)*conn;
      switch(t)
      {
        case INTERP_KERNEL::NORM_QUAD4:
          {
            FillInCompact3DMode(3,4,conn+1,coo,tmp);
            *pt=INTERP_KERNEL::quadSkew(tmp);
            break;
          }
        default:
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getSkewField : A cell with not manged type (NORM_QUAD4) has been detected !");
      }
      conn+=connI[i+1]-connI[i];
    }
  ret->setName("Skew");
  ret->synchronizeTimeWithSupport();
  return ret.retn();
}

/*!
 * Returns the cell field giving for each cell in \a this its diameter. Diameter means the max length of all possible SEG2 in the cell.
 *
 * \return a new instance of field containing the result. The returned instance has to be deallocated by the caller.
 *
 * \sa getSkewField, getWarpField, getAspectRatioField, getEdgeRatioField
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::computeDiameterField() const
{
  checkConsistencyLight();
  MCAuto<MEDCouplingFieldDouble> ret(MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME));
  ret->setMesh(this);
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  ComputeAllTypesInternal(types,_nodal_connec,_nodal_connec_index);
  int spaceDim(getSpaceDimension());
  mcIdType nbCells(getNumberOfCells());
  MCAuto<DataArrayDouble> arr(DataArrayDouble::New());
  arr->alloc(nbCells,1);
  for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator it=types.begin();it!=types.end();it++)
    {
      INTERP_KERNEL::AutoCppPtr<INTERP_KERNEL::DiameterCalculator> dc(INTERP_KERNEL::CellModel::GetCellModel(*it).buildInstanceOfDiameterCalulator(spaceDim));
      MCAuto<DataArrayIdType> cellIds(giveCellsWithType(*it));
      dc->computeForListOfCellIdsUMeshFrmt(cellIds->begin(),cellIds->end(),_nodal_connec_index->begin(),_nodal_connec->begin(),getCoords()->begin(),arr->getPointer());
    }
  ret->setArray(arr);
  ret->setName("Diameter");
  return ret.retn();
}

/*!
 * This method aggregate the bbox of each cell and put it into bbox parameter (xmin,xmax,ymin,ymax,zmin,zmax).
 *
 * \param [in] arcDetEps - a parameter specifying in case of 2D quadratic polygon cell the detection limit between linear and arc circle. (By default 1e-12)
 *                         For all other cases this input parameter is ignored.
 * \return DataArrayDouble * - newly created object (to be managed by the caller) \a this number of cells tuples and 2*spacedim components.
 *
 * \throw If \a this is not fully set (coordinates and connectivity).
 * \throw If a cell in \a this has no valid nodeId.
 * \sa MEDCouplingUMesh::getBoundingBoxForBBTreeFast, MEDCouplingUMesh::getBoundingBoxForBBTree2DQuadratic
 */
DataArrayDouble *MEDCouplingUMesh::getBoundingBoxForBBTree(double arcDetEps) const
{
  int mDim(getMeshDimension()),sDim(getSpaceDimension());
  if((mDim==3 && sDim==3) || (mDim==2 && sDim==3) || (mDim==1 && sDim==1) || ( mDim==1 && sDim==3))  // Compute refined boundary box for quadratic elements only in 2D.
    return getBoundingBoxForBBTreeFast();
  if((mDim==2 && sDim==2) || (mDim==1 && sDim==2))
    {
      bool presenceOfQuadratic(false);
      for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator it=_types.begin();it!=_types.end();it++)
        {
          const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(*it));
          if(cm.isQuadratic())
            presenceOfQuadratic=true;
        }
      if(!presenceOfQuadratic)
        return getBoundingBoxForBBTreeFast();
      if(mDim==2 && sDim==2)
        return getBoundingBoxForBBTree2DQuadratic(arcDetEps);
      else
        return getBoundingBoxForBBTree1DQuadratic(arcDetEps);
    }
  throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getBoundingBoxForBBTree : Managed dimensions are (mDim=1,sDim=1), (mDim=1,sDim=2), (mDim=1,sDim=3), (mDim=2,sDim=2), (mDim=2,sDim=3) and (mDim=3,sDim=3) !");
}

/*!
 * This method aggregate the bbox of each cell only considering the nodes constituting each cell and put it into bbox parameter.
 * So meshes having quadratic cells the computed bounding boxes can be invalid !
 *
 * \return DataArrayDouble * - newly created object (to be managed by the caller) \a this number of cells tuples and 2*spacedim components.
 *
 * \throw If \a this is not fully set (coordinates and connectivity).
 * \throw If a cell in \a this has no valid nodeId.
 */
DataArrayDouble *MEDCouplingUMesh::getBoundingBoxForBBTreeFast() const
{
  checkFullyDefined();
  int spaceDim(getSpaceDimension());
  mcIdType nbOfCells(getNumberOfCells()), nbOfNodes(getNumberOfNodes());
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New()); ret->alloc(nbOfCells,2*spaceDim);
  double *bbox(ret->getPointer());
  for(mcIdType i=0;i<nbOfCells*spaceDim;i++)
    {
      bbox[2*i]=std::numeric_limits<double>::max();
      bbox[2*i+1]=-std::numeric_limits<double>::max();
    }
  const double *coordsPtr(_coords->begin());
  const mcIdType *conn(_nodal_connec->begin()),*connI(_nodal_connec_index->begin());
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      mcIdType offset=connI[i]+1;
      mcIdType nbOfNodesForCell(connI[i+1]-offset),kk(0);
      for(mcIdType j=0;j<nbOfNodesForCell;j++)
        {
          mcIdType nodeId=conn[offset+j];
          if(nodeId>=0 && nodeId<nbOfNodes)
            {
              for(int k=0;k<spaceDim;k++)
                {
                  bbox[2*spaceDim*i+2*k]=std::min(bbox[2*spaceDim*i+2*k],coordsPtr[spaceDim*nodeId+k]);
                  bbox[2*spaceDim*i+2*k+1]=std::max(bbox[2*spaceDim*i+2*k+1],coordsPtr[spaceDim*nodeId+k]);
                }
              kk++;
            }
        }
      if(kk==0)
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::getBoundingBoxForBBTree : cell #" << i << " contains no valid nodeId !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  return ret.retn();
}

/*!
 * This method aggregates the bbox of each 2D cell in \a this considering the whole shape. This method is particularly
 * useful for 2D meshes having quadratic cells
 * because for this type of cells getBoundingBoxForBBTreeFast method may return invalid bounding boxes (since it just considers
 * the two extremities of the arc of circle).
 *
 * \param [in] arcDetEps - a parameter specifying in case of 2D quadratic polygon cell the detection limit between linear and arc circle. (By default 1e-12)
 * \return DataArrayDouble * - newly created object (to be managed by the caller) \a this number of cells tuples and 2*spacedim components.
 * \throw If \a this is not fully defined.
 * \throw If \a this is not a mesh with meshDimension equal to 2.
 * \throw If \a this is not a mesh with spaceDimension equal to 2.
 * \sa MEDCouplingUMesh::getBoundingBoxForBBTree1DQuadratic
 */
DataArrayDouble *MEDCouplingUMesh::getBoundingBoxForBBTree2DQuadratic(double arcDetEps) const
{
  checkFullyDefined();
  INTERP_KERNEL::QuadraticPlanarPrecision arcPrec(arcDetEps);

  int spaceDim(getSpaceDimension()),mDim(getMeshDimension());
  mcIdType nbOfCells=getNumberOfCells();
  if(spaceDim!=2 || mDim!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getBoundingBoxForBBTree2DQuadratic : This method should be applied on mesh with mesh dimension equal to 2 and space dimension also equal to 2!");
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New()); ret->alloc(nbOfCells,2*spaceDim);
  double *bbox(ret->getPointer());
  const double *coords(_coords->begin());
  const mcIdType *conn(_nodal_connec->begin()),*connI(_nodal_connec_index->begin());
  for(mcIdType i=0;i<nbOfCells;i++,bbox+=4,connI++)
    {
      const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[*connI]));
      mcIdType sz(connI[1]-connI[0]-1);
      std::vector<INTERP_KERNEL::Node *> nodes(sz);
      INTERP_KERNEL::QuadraticPolygon *pol(0);
      for(mcIdType j=0;j<sz;j++)
        {
          mcIdType nodeId(conn[*connI+1+j]);
          nodes[j]=new INTERP_KERNEL::Node(coords[nodeId*2],coords[nodeId*2+1]);
        }
      if(!cm.isQuadratic())
        pol=INTERP_KERNEL::QuadraticPolygon::BuildLinearPolygon(nodes);
      else
        pol=INTERP_KERNEL::QuadraticPolygon::BuildArcCirclePolygon(nodes);
      INTERP_KERNEL::Bounds b; b.prepareForAggregation(); pol->fillBounds(b); delete pol;
      bbox[0]=b.getXMin(); bbox[1]=b.getXMax(); bbox[2]=b.getYMin(); bbox[3]=b.getYMax();
    }
  return ret.retn();
}

/*!
 * This method aggregates the bbox of each 1D cell in \a this considering the whole shape. This method is particularly
 * useful for 2D meshes having quadratic cells
 * because for this type of cells getBoundingBoxForBBTreeFast method may return invalid bounding boxes (since it just considers
 * the two extremities of the arc of circle).
 *
 * \param [in] arcDetEps - a parameter specifying in case of 2D quadratic polygon cell the detection limit between linear and arc circle. (By default 1e-12)
 * \return DataArrayDouble * - newly created object (to be managed by the caller) \a this number of cells tuples and 2*spacedim components.
 * \throw If \a this is not fully defined.
 * \throw If \a this is not a mesh with meshDimension equal to 1.
 * \throw If \a this is not a mesh with spaceDimension equal to 2.
 * \sa MEDCouplingUMesh::getBoundingBoxForBBTree2DQuadratic
 */
DataArrayDouble *MEDCouplingUMesh::getBoundingBoxForBBTree1DQuadratic(double arcDetEps) const
{
  checkFullyDefined();
  int spaceDim(getSpaceDimension()),mDim(getMeshDimension());
  mcIdType nbOfCells=getNumberOfCells();
  if(spaceDim!=2 || mDim!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getBoundingBoxForBBTree1DQuadratic : This method should be applied on mesh with mesh dimension equal to 1 and space dimension also equal to 2!");
  INTERP_KERNEL::QuadraticPlanarPrecision arcPrec(arcDetEps);
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New()); ret->alloc(nbOfCells,2*spaceDim);
  double *bbox(ret->getPointer());
  const double *coords(_coords->begin());
  const mcIdType *conn(_nodal_connec->begin()),*connI(_nodal_connec_index->begin());
  for(mcIdType i=0;i<nbOfCells;i++,bbox+=4,connI++)
    {
      const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[*connI]));
      mcIdType sz(connI[1]-connI[0]-1);
      std::vector<INTERP_KERNEL::Node *> nodes(sz);
      INTERP_KERNEL::Edge *edge(0);
      for(mcIdType j=0;j<sz;j++)
        {
          mcIdType nodeId(conn[*connI+1+j]);
          nodes[j]=new INTERP_KERNEL::Node(coords[nodeId*2],coords[nodeId*2+1]);
        }
      if(!cm.isQuadratic())
        edge=INTERP_KERNEL::QuadraticPolygon::BuildLinearEdge(nodes);
      else
        edge=INTERP_KERNEL::QuadraticPolygon::BuildArcCircleEdge(nodes);
      const INTERP_KERNEL::Bounds& b(edge->getBounds());
      bbox[0]=b.getXMin(); bbox[1]=b.getXMax(); bbox[2]=b.getYMin(); bbox[3]=b.getYMax(); edge->decrRef();
    }
  return ret.retn();
}

/// @cond INTERNAL

namespace MEDCouplingImpl
{
  class ConnReader
  {
  public:
    ConnReader(const mcIdType *c, mcIdType val):_conn(c),_val(val) { }
    bool operator() (const mcIdType& pos) { return _conn[pos]!=_val; }
  private:
    const mcIdType *_conn;
    mcIdType _val;
  };

  class ConnReader2
  {
  public:
    ConnReader2(const mcIdType *c, mcIdType val):_conn(c),_val(val) { }
    bool operator() (const mcIdType& pos) { return _conn[pos]==_val; }
  private:
    const mcIdType *_conn;
    mcIdType _val;
  };
}

/// @endcond

/*!
 * This method expects that \a this is sorted by types. If not an exception will be thrown.
 * This method returns in the same format as code (see MEDCouplingUMesh::checkTypeConsistencyAndContig or MEDCouplingUMesh::splitProfilePerType) how
 * \a this is composed in cell types.
 * The returned array is of size 3*n where n is the number of different types present in \a this.
 * For every k in [0,n] ret[3*k+2]==-1 because it has no sense here.
 * This parameter is kept only for compatibility with other method listed above.
 */
std::vector<mcIdType> MEDCouplingUMesh::getDistributionOfTypes() const
{
  checkConnectivityFullyDefined();
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  const mcIdType *work=connI;
  mcIdType nbOfCells=getNumberOfCells();
  std::size_t n=getAllGeoTypes().size();
  std::vector<mcIdType> ret(3*n,-1); //ret[3*k+2]==-1 because it has no sense here
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  for(std::size_t i=0;work!=connI+nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType)conn[*work];
      if(types.find(typ)!=types.end())
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::getDistributionOfTypes : Type " << INTERP_KERNEL::CellModel::GetCellModel(typ).getRepr();
          oss << " is not contiguous !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      types.insert(typ);
      ret[3*i]=typ;
      const mcIdType *work2=std::find_if(work+1,connI+nbOfCells,MEDCouplingImpl::ConnReader(conn,typ));
      ret[3*i+1]=ToIdType(std::distance(work,work2));
      work=work2;
    }
  return ret;
}

/*!
 * This method is used to check that this has contiguous cell type in same order than described in \a code.
 * only for types cell, type node is not managed.
 * Format of \a code is the following. \a code should be of size 3*n and non empty. If not an exception is thrown.
 * foreach k in [0,n) on 3*k pos represent the geometric type and 3*k+1 number of elements of type 3*k.
 * 3*k+2 refers if different from -1 the pos in 'idsPerType' to get the corresponding array.
 * If 2 or more same geometric type is in \a code and exception is thrown too.
 *
 * This method firstly checks
 * If it exists k so that 3*k geometric type is not in geometric types of this an exception will be thrown.
 * If it exists k so that 3*k geometric type exists but the number of consecutive cell types does not match,
 * an exception is thrown too.
 *
 * If all geometric types in \a code are exactly those in \a this null pointer is returned.
 * If it exists a geometric type in \a this \b not in \a code \b no exception is thrown
 * and a DataArrayIdType instance is returned that the user has the responsibility to deallocate.
 */
DataArrayIdType *MEDCouplingUMesh::checkTypeConsistencyAndContig(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const
{
  if(code.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::checkTypeConsistencyAndContig : code is empty, should not !");
  std::size_t sz=code.size();
  std::size_t n=sz/3;
  if(sz%3!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::checkTypeConsistencyAndContig : code size is NOT %3 !");
  std::vector<INTERP_KERNEL::NormalizedCellType> types;
  mcIdType nb=0;
  bool isNoPflUsed=true;
  for(std::size_t i=0;i<n;i++)
    if(std::find(types.begin(),types.end(),(INTERP_KERNEL::NormalizedCellType)code[3*i])==types.end())
      {
        types.push_back((INTERP_KERNEL::NormalizedCellType)code[3*i]);
        nb+=code[3*i+1];
        if(_types.find((INTERP_KERNEL::NormalizedCellType)code[3*i])==_types.end())
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::checkTypeConsistencyAndContig : expected geo types not in this !");
        isNoPflUsed=isNoPflUsed && (code[3*i+2]==-1);
      }
  if(types.size()!=n)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::checkTypeConsistencyAndContig : code contains duplication of types in unstructured mesh !");
  if(isNoPflUsed)
    {
      if(!checkConsecutiveCellTypesAndOrder(&types[0],&types[0]+types.size()))
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::checkTypeConsistencyAndContig : non contiguous type !");
      if(types.size()==_types.size())
        return 0;
    }
  MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
  ret->alloc(nb,1);
  mcIdType *retPtr=ret->getPointer();
  const mcIdType *connI=_nodal_connec_index->begin();
  const mcIdType *conn=_nodal_connec->begin();
  mcIdType nbOfCells=getNumberOfCells();
  const mcIdType *i=connI;
  int kk=0;
  for(std::vector<INTERP_KERNEL::NormalizedCellType>::const_iterator it=types.begin();it!=types.end();it++,kk++)
    {
      i=std::find_if(i,connI+nbOfCells,MEDCouplingImpl::ConnReader2(conn,ToIdType((*it))));
      mcIdType offset=ToIdType(std::distance(connI,i));
      const mcIdType *j=std::find_if(i+1,connI+nbOfCells,MEDCouplingImpl::ConnReader(conn,ToIdType((*it))));
      mcIdType nbOfCellsOfCurType=ToIdType(std::distance(i,j));
      if(code[3*kk+2]==-1)
        for(mcIdType k=0;k<nbOfCellsOfCurType;k++)
          *retPtr++=k+offset;
      else
        {
          mcIdType idInIdsPerType=code[3*kk+2];
          if(idInIdsPerType>=0 && idInIdsPerType<ToIdType(idsPerType.size()))
            {
              const DataArrayIdType *zePfl=idsPerType[idInIdsPerType];
              if(zePfl)
                {
                  zePfl->checkAllocated();
                  if(zePfl->getNumberOfComponents()==1)
                    {
                      for(const mcIdType *k=zePfl->begin();k!=zePfl->end();k++,retPtr++)
                        {
                          if(*k>=0 && *k<nbOfCellsOfCurType)
                            *retPtr=(*k)+offset;
                          else
                            {
                              std::ostringstream oss; oss << "MEDCouplingUMesh::checkTypeConsistencyAndContig : the section " << kk << " points to the profile #" << idInIdsPerType;
                              oss << ", and this profile contains a value " << *k << " should be in [0," << nbOfCellsOfCurType << ") !";
                              throw INTERP_KERNEL::Exception(oss.str());
                            }
                        }
                    }
                  else
                    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::checkTypeConsistencyAndContig : presence of a profile with nb of compo != 1 !");
                }
              else
                throw INTERP_KERNEL::Exception("MEDCouplingUMesh::checkTypeConsistencyAndContig : presence of null profile !");
            }
          else
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::checkTypeConsistencyAndContig : at section " << kk << " of code it points to the array #" << idInIdsPerType;
              oss << " should be in [0," << idsPerType.size() << ") !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
        }
      i=j;
    }
  return ret.retn();
}

/*!
 * This method makes the hypothesis that \a this is sorted by type. If not an exception will be thrown.
 * This method is the opposite of MEDCouplingUMesh::checkTypeConsistencyAndContig method. Given a list of cells in \a profile it returns a list of sub-profiles sorted by geo type.
 * The result is put in the array \a idsPerType. In the returned parameter \a code, foreach i \a code[3*i+2] refers (if different from -1) to a location into the \a idsPerType.
 * This method has 1 input \a profile and 3 outputs \a code \a idsInPflPerType and \a idsPerType.
 *
 * \param [in] profile list of IDs constituing the profile
 * \param [out] code is a vector of size 3*n where n is the number of different geometric type in \a this \b reduced to the profile \a profile. \a code has exactly the same semantic than in MEDCouplingUMesh::checkTypeConsistencyAndContig method.
 * \param [out] idsInPflPerType is a vector of size of different geometric type in the subpart defined by \a profile of \a this ( equal to \a code.size()/3). For each i,
 *              \a idsInPflPerType[i] stores the tuple ids in \a profile that correspond to the geometric type code[3*i+0]
 * \param [out] idsPerType is a vector of size of different sub profiles needed to be defined to represent the profile \a profile for a given geometric type.
 *              This vector can be empty in case of all geometric type cells are fully covered in ascending in the given input \a profile.
 * \throw if \a profile has not exactly one component. It throws too, if \a profile contains some values not in [0,getNumberOfCells()) or if \a this is not fully defined
 */
void MEDCouplingUMesh::splitProfilePerType(const DataArrayIdType *profile, std::vector<mcIdType>& code, std::vector<DataArrayIdType *>& idsInPflPerType, std::vector<DataArrayIdType *>& idsPerType, bool smartPflKiller) const
{
  if(!profile)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::splitProfilePerType : input profile is NULL !");
  if(profile->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::splitProfilePerType : input profile should have exactly one component !");
  checkConnectivityFullyDefined();
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  mcIdType nbOfCells=getNumberOfCells();
  std::vector<INTERP_KERNEL::NormalizedCellType> types;
  std::vector<mcIdType> typeRangeVals(1);
  for(const mcIdType *i=connI;i!=connI+nbOfCells;)
    {
      INTERP_KERNEL::NormalizedCellType curType=(INTERP_KERNEL::NormalizedCellType)conn[*i];
      if(std::find(types.begin(),types.end(),curType)!=types.end())
        {
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::splitProfilePerType : current mesh is not sorted by type !");
        }
      types.push_back(curType);
      i=std::find_if(i+1,connI+nbOfCells,MEDCouplingImpl::ConnReader(conn,ToIdType(curType)));
      typeRangeVals.push_back(ToIdType(std::distance(connI,i)));
    }
  //
  DataArrayIdType *castArr=0,*rankInsideCast=0,*castsPresent=0;
  profile->splitByValueRange(&typeRangeVals[0],&typeRangeVals[0]+typeRangeVals.size(),castArr,rankInsideCast,castsPresent);
  MCAuto<DataArrayIdType> tmp0=castArr;
  MCAuto<DataArrayIdType> tmp1=rankInsideCast;
  MCAuto<DataArrayIdType> tmp2=castsPresent;
  //
  mcIdType nbOfCastsFinal=castsPresent->getNumberOfTuples();
  code.resize(3*nbOfCastsFinal);
  std::vector< MCAuto<DataArrayIdType> > idsInPflPerType2;
  std::vector< MCAuto<DataArrayIdType> > idsPerType2;
  for(mcIdType i=0;i<nbOfCastsFinal;i++)
    {
      mcIdType castId=castsPresent->getIJ(i,0);
      MCAuto<DataArrayIdType> tmp3=castArr->findIdsEqual(castId);
      idsInPflPerType2.push_back(tmp3);
      code[3*i]=ToIdType(types[castId]);
      code[3*i+1]=tmp3->getNumberOfTuples();
      MCAuto<DataArrayIdType> tmp4=rankInsideCast->selectByTupleId(tmp3->begin(),tmp3->begin()+tmp3->getNumberOfTuples());
      if(!smartPflKiller || !tmp4->isIota(typeRangeVals[castId+1]-typeRangeVals[castId]))
        {
          tmp4->copyStringInfoFrom(*profile);
          idsPerType2.push_back(tmp4);
          code[3*i+2]=ToIdType(idsPerType2.size())-1;
        }
      else
        {
          code[3*i+2]=-1;
        }
    }
  std::size_t sz2=idsInPflPerType2.size();
  idsInPflPerType.resize(sz2);
  for(std::size_t i=0;i<sz2;i++)
    {
      DataArrayIdType *locDa=idsInPflPerType2[i];
      locDa->incrRef();
      idsInPflPerType[i]=locDa;
    }
  std::size_t sz=idsPerType2.size();
  idsPerType.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      DataArrayIdType *locDa=idsPerType2[i];
      locDa->incrRef();
      idsPerType[i]=locDa;
    }
}

/*!
 * This method is here too emulate the MEDMEM behaviour on BDC (buildDescendingConnectivity). Hoping this method becomes deprecated very soon.
 * This method make the assumption that \a this and 'nM1LevMesh' mesh lyies on same coords (same pointer) as MED and MEDMEM does.
 * The following equality should be verified 'nM1LevMesh->getMeshDimension()==this->getMeshDimension()-1'
 * This method returns 5+2 elements. 'desc', 'descIndx', 'revDesc', 'revDescIndx' and 'meshnM1' behaves exactly as MEDCoupling::MEDCouplingUMesh::buildDescendingConnectivity except the content as described after. The returned array specifies the n-1 mesh reordered by type as MEDMEM does. 'nM1LevMeshIds' contains the ids in returned 'meshnM1'. Finally 'meshnM1Old2New' contains numbering old2new that is to say the cell #k in coarse 'nM1LevMesh' will have the number ret[k] in returned mesh 'nM1LevMesh' MEDMEM reordered.
 */
MEDCouplingUMesh *MEDCouplingUMesh::emulateMEDMEMBDC(const MEDCouplingUMesh *nM1LevMesh, DataArrayIdType *desc, DataArrayIdType *descIndx, DataArrayIdType *&revDesc, DataArrayIdType *&revDescIndx, DataArrayIdType *& nM1LevMeshIds, DataArrayIdType *&meshnM1Old2New) const
{
  checkFullyDefined();
  nM1LevMesh->checkFullyDefined();
  if(getMeshDimension()-1!=nM1LevMesh->getMeshDimension())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::emulateMEDMEMBDC : The mesh passed as first argument should have a meshDim equal to this->getMeshDimension()-1 !" );
  if(_coords!=nM1LevMesh->getCoords())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::emulateMEDMEMBDC : 'this' and mesh in first argument should share the same coords : Use tryToShareSameCoords method !");
  MCAuto<DataArrayIdType> tmp0=DataArrayIdType::New();
  MCAuto<DataArrayIdType> tmp1=DataArrayIdType::New();
  MCAuto<MEDCouplingUMesh> ret1=buildDescendingConnectivity(desc,descIndx,tmp0,tmp1);
  MCAuto<DataArrayIdType> ret0=ret1->sortCellsInMEDFileFrmt();
  desc->transformWithIndArr(ret0->begin(),ret0->begin()+ret0->getNbOfElems());
  MCAuto<MEDCouplingUMesh> tmp=MEDCouplingUMesh::New();
  tmp->setConnectivity(tmp0,tmp1);
  tmp->renumberCells(ret0->begin(),false);
  revDesc=tmp->getNodalConnectivity();
  revDescIndx=tmp->getNodalConnectivityIndex();
  DataArrayIdType *ret=0;
  if(!ret1->areCellsIncludedIn(nM1LevMesh,2,ret))
    {
      mcIdType tmp2;
      ret->getMaxValue(tmp2);
      ret->decrRef();
      std::ostringstream oss; oss << "MEDCouplingUMesh::emulateMEDMEMBDC : input N-1 mesh present a cell not in descending mesh ... Id of cell is " << tmp2 << " !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  nM1LevMeshIds=ret;
  //
  revDesc->incrRef();
  revDescIndx->incrRef();
  ret1->incrRef();
  ret0->incrRef();
  meshnM1Old2New=ret0;
  return ret1;
}

/*!
 * Permutes the nodal connectivity arrays so that the cells are sorted by type, which is
 * necessary for writing the mesh to MED file. Additionally returns a permutation array
 * in "Old to New" mode.
 *  \return DataArrayIdType * - a new instance of DataArrayIdType. The caller is to delete
 *          this array using decrRef() as it is no more needed.
 *  \throw If the nodal connectivity of cells is not defined.
 */
DataArrayIdType *MEDCouplingUMesh::sortCellsInMEDFileFrmt()
{
  checkConnectivityFullyDefined();
  MCAuto<DataArrayIdType> ret=getRenumArrForMEDFileFrmt();
  renumberCells(ret->begin(),false);
  return ret.retn();
}

/*!
 * This methods checks that cells are sorted by their types.
 * This method makes asumption (no check) that connectivity is correctly set before calling.
 */
bool MEDCouplingUMesh::checkConsecutiveCellTypes() const
{
  checkFullyDefined();
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  mcIdType nbOfCells=getNumberOfCells();
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  for(const mcIdType *i=connI;i!=connI+nbOfCells;)
    {
      INTERP_KERNEL::NormalizedCellType curType=(INTERP_KERNEL::NormalizedCellType)conn[*i];
      if(types.find(curType)!=types.end())
        return false;
      types.insert(curType);
      i=std::find_if(i+1,connI+nbOfCells,MEDCouplingImpl::ConnReader(conn,ToIdType(curType)));
    }
  return true;
}

/*!
 * This method is a specialization of MEDCouplingUMesh::checkConsecutiveCellTypesAndOrder method that is called here.
 * The geometric type order is specified by MED file.
 *
 * \sa  MEDCouplingUMesh::checkConsecutiveCellTypesAndOrder
 */
bool MEDCouplingUMesh::checkConsecutiveCellTypesForMEDFileFrmt() const
{
  return checkConsecutiveCellTypesAndOrder(MEDMEM_ORDER,MEDMEM_ORDER+N_MEDMEM_ORDER);
}

/*!
 * This method performs the same job as checkConsecutiveCellTypes except that the order of types sequence is analyzed to check
 * that the order is specified in array defined by [ \a orderBg , \a orderEnd ).
 * If there is some geo types in \a this \b NOT in [ \a orderBg, \a orderEnd ) it is OK (return true) if contiguous.
 * If there is some geo types in [ \a orderBg, \a orderEnd ) \b NOT in \a this it is OK too (return true) if contiguous.
 */
bool MEDCouplingUMesh::checkConsecutiveCellTypesAndOrder(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd) const
{
  checkFullyDefined();
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  mcIdType nbOfCells=getNumberOfCells();
  if(nbOfCells==0)
    return true;
  mcIdType lastPos=-1;
  std::set<INTERP_KERNEL::NormalizedCellType> sg;
  for(const mcIdType *i=connI;i!=connI+nbOfCells;)
    {
      INTERP_KERNEL::NormalizedCellType curType=(INTERP_KERNEL::NormalizedCellType)conn[*i];
      const INTERP_KERNEL::NormalizedCellType *isTypeExists=std::find(orderBg,orderEnd,curType);
      if(isTypeExists!=orderEnd)
        {
          mcIdType pos=ToIdType(std::distance(orderBg,isTypeExists));
          if(pos<=lastPos)
            return false;
          lastPos=pos;
          i=std::find_if(i+1,connI+nbOfCells,MEDCouplingImpl::ConnReader(conn,ToIdType(curType)));
        }
      else
        {
          if(sg.find(curType)==sg.end())
            {
              i=std::find_if(i+1,connI+nbOfCells,MEDCouplingImpl::ConnReader(conn,ToIdType(curType)));
              sg.insert(curType);
            }
          else
            return false;
        }
    }
  return true;
}

/*!
 * This method returns 2 newly allocated DataArrayIdType instances. The first is an array of size 'this->getNumberOfCells()' with one component,
 * that tells for each cell the pos of its type in the array on type given in input parameter. The 2nd output parameter is an array with the same
 * number of tuples than input type array and with one component. This 2nd output array gives type by type the number of occurrence of type in 'this'.
 */
DataArrayIdType *MEDCouplingUMesh::getLevArrPerCellTypes(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd, DataArrayIdType *&nbPerType) const
{
  checkConnectivityFullyDefined();
  mcIdType nbOfCells=getNumberOfCells();
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  MCAuto<DataArrayIdType> tmpa=DataArrayIdType::New();
  MCAuto<DataArrayIdType> tmpb=DataArrayIdType::New();
  tmpa->alloc(nbOfCells,1);
  tmpb->alloc(std::distance(orderBg,orderEnd),1);
  tmpb->fillWithZero();
  mcIdType *tmp=tmpa->getPointer();
  mcIdType *tmp2=tmpb->getPointer();
  for(const mcIdType *i=connI;i!=connI+nbOfCells;i++)
    {
      const INTERP_KERNEL::NormalizedCellType *where=std::find(orderBg,orderEnd,(INTERP_KERNEL::NormalizedCellType)conn[*i]);
      if(where!=orderEnd)
        {
          mcIdType pos=ToIdType(std::distance(orderBg,where));
          tmp2[pos]++;
          tmp[std::distance(connI,i)]=pos;
        }
      else
        {
          const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[*i]);
          std::ostringstream oss; oss << "MEDCouplingUMesh::getLevArrPerCellTypes : Cell #" << std::distance(connI,i);
          oss << " has a type " << cm.getRepr() << " not in input array of type !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  nbPerType=tmpb.retn();
  return tmpa.retn();
}

/*!
 * This method behaves exactly as MEDCouplingUMesh::getRenumArrForConsecutiveCellTypesSpec but the order is those defined in MED file spec.
 *
 * \return a new object containing the old to new correspondence.
 *
 * \sa MEDCouplingUMesh::getRenumArrForConsecutiveCellTypesSpec, MEDCouplingUMesh::sortCellsInMEDFileFrmt.
 */
DataArrayIdType *MEDCouplingUMesh::getRenumArrForMEDFileFrmt() const
{
  return getRenumArrForConsecutiveCellTypesSpec(MEDMEM_ORDER,MEDMEM_ORDER+N_MEDMEM_ORDER);
}

/*!
 * This method is similar to method MEDCouplingUMesh::rearrange2ConsecutiveCellTypes except that the type order is specified by [ \a orderBg , \a orderEnd ) (as MEDCouplingUMesh::checkConsecutiveCellTypesAndOrder method) and that this method is \b const and performs \b NO permutation in \a this.
 * This method returns an array of size getNumberOfCells() that gives a renumber array old2New that can be used as input of MEDCouplingMesh::renumberCells.
 * The mesh after this call to MEDCouplingMesh::renumberCells will pass the test of MEDCouplingUMesh::checkConsecutiveCellTypesAndOrder with the same inputs.
 * The returned array minimizes the permutations that is to say the order of cells inside same geometric type remains the same.
 */
DataArrayIdType *MEDCouplingUMesh::getRenumArrForConsecutiveCellTypesSpec(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd) const
{
  DataArrayIdType *nbPerType=0;
  MCAuto<DataArrayIdType> tmpa=getLevArrPerCellTypes(orderBg,orderEnd,nbPerType);
  nbPerType->decrRef();
  return tmpa->buildPermArrPerLevel();
}

/*!
 * This method reorganize the cells of \a this so that the cells with same geometric types are put together.
 * The number of cells remains unchanged after the call of this method.
 * This method tries to minimizes the number of needed permutations. So, this method behaves not exactly as
 * MEDCouplingUMesh::sortCellsInMEDFileFrmt.
 *
 * \return the array giving the correspondence old to new.
 */
DataArrayIdType *MEDCouplingUMesh::rearrange2ConsecutiveCellTypes()
{
  checkFullyDefined();
  computeTypes();
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  mcIdType nbOfCells=getNumberOfCells();
  std::vector<INTERP_KERNEL::NormalizedCellType> types;
  for(const mcIdType *i=connI;i!=connI+nbOfCells && (types.size()!=_types.size());)
    if(std::find(types.begin(),types.end(),(INTERP_KERNEL::NormalizedCellType)conn[*i])==types.end())
      {
        INTERP_KERNEL::NormalizedCellType curType=(INTERP_KERNEL::NormalizedCellType)conn[*i];
        types.push_back(curType);
        for(i++;i!=connI+nbOfCells && (INTERP_KERNEL::NormalizedCellType)conn[*i]==curType;i++);
      }
  DataArrayIdType *ret=DataArrayIdType::New();
  ret->alloc(nbOfCells,1);
  mcIdType *retPtr=ret->getPointer();
  std::fill(retPtr,retPtr+nbOfCells,-1);
  mcIdType newCellId=0;
  for(std::vector<INTERP_KERNEL::NormalizedCellType>::const_iterator iter=types.begin();iter!=types.end();iter++)
    {
      for(const mcIdType *i=connI;i!=connI+nbOfCells;i++)
        if((INTERP_KERNEL::NormalizedCellType)conn[*i]==(*iter))
          retPtr[std::distance(connI,i)]=newCellId++;
    }
  renumberCells(retPtr,false);
  return ret;
}

/*!
 * This method splits \a this into as mush as untructured meshes that consecutive set of same type cells.
 * So this method has typically a sense if MEDCouplingUMesh::checkConsecutiveCellTypes has a sense.
 * This method makes asumption that connectivity is correctly set before calling.
 */
std::vector<MEDCouplingUMesh *> MEDCouplingUMesh::splitByType() const
{
  checkConnectivityFullyDefined();
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  mcIdType nbOfCells=getNumberOfCells();
  std::vector<MEDCouplingUMesh *> ret;
  for(const mcIdType *i=connI;i!=connI+nbOfCells;)
    {
      INTERP_KERNEL::NormalizedCellType curType=(INTERP_KERNEL::NormalizedCellType)conn[*i];
      mcIdType beginCellId=ToIdType(std::distance(connI,i));
      i=std::find_if(i+1,connI+nbOfCells,MEDCouplingImpl::ConnReader(conn,ToIdType(curType)));
      mcIdType endCellId=ToIdType(std::distance(connI,i));
      mcIdType sz=endCellId-beginCellId;
      mcIdType *cells=new mcIdType[sz];
      for(mcIdType j=0;j<sz;j++)
        cells[j]=beginCellId+j;
      MEDCouplingUMesh *m=(MEDCouplingUMesh *)buildPartOfMySelf(cells,cells+sz,true);
      delete [] cells;
      ret.push_back(m);
    }
  return ret;
}

/*!
 * This method performs the opposite operation than those in MEDCoupling1SGTUMesh::buildUnstructured.
 * If \a this is a single geometric type unstructured mesh, it will be converted into a more compact data structure,
 * MEDCoupling1GTUMesh instance. The returned instance will aggregate the same DataArrayDouble instance of coordinates than \a this.
 *
 * \return a newly allocated instance, that the caller must manage.
 * \throw If \a this contains more than one geometric type.
 * \throw If the nodal connectivity of \a this is not fully defined.
 * \throw If the internal data is not coherent.
 */
MEDCoupling1GTUMesh *MEDCouplingUMesh::convertIntoSingleGeoTypeMesh() const
{
  checkConnectivityFullyDefined();
  if(_types.size()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::convertIntoSingleGeoTypeMesh : current mesh does not contain exactly one geometric type !");
  INTERP_KERNEL::NormalizedCellType typ=*_types.begin();
  MCAuto<MEDCoupling1GTUMesh> ret=MEDCoupling1GTUMesh::New(getName(),typ);
  ret->setCoords(getCoords());
  MEDCoupling1SGTUMesh *retC=dynamic_cast<MEDCoupling1SGTUMesh *>((MEDCoupling1GTUMesh*)ret);
  if(retC)
    {
      MCAuto<DataArrayIdType> c=convertNodalConnectivityToStaticGeoTypeMesh();
      retC->setNodalConnectivity(c);
    }
  else
    {
      MEDCoupling1DGTUMesh *retD=dynamic_cast<MEDCoupling1DGTUMesh *>((MEDCoupling1GTUMesh*)ret);
      if(!retD)
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::convertIntoSingleGeoTypeMesh : Internal error !");
      DataArrayIdType *c=0,*ci=0;
      convertNodalConnectivityToDynamicGeoTypeMesh(c,ci);
      MCAuto<DataArrayIdType> cs(c),cis(ci);
      retD->setNodalConnectivity(cs,cis);
    }
  return ret.retn();
}

DataArrayIdType *MEDCouplingUMesh::convertNodalConnectivityToStaticGeoTypeMesh() const
{
  checkConnectivityFullyDefined();
  if(_types.size()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::convertNodalConnectivityToStaticGeoTypeMesh : current mesh does not contain exactly one geometric type !");
  INTERP_KERNEL::NormalizedCellType typ=*_types.begin();
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(typ);
  if(cm.isDynamic())
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::convertNodalConnectivityToStaticGeoTypeMesh : this contains a single geo type (" << cm.getRepr() << ") but ";
      oss << "this type is dynamic ! Only static geometric type is possible for that type ! call convertNodalConnectivityToDynamicGeoTypeMesh instead !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  mcIdType nbCells=getNumberOfCells();
  mcIdType typi=ToIdType(typ);
  mcIdType nbNodesPerCell=ToIdType(cm.getNumberOfNodes());
  MCAuto<DataArrayIdType> connOut=DataArrayIdType::New(); connOut->alloc(nbCells*nbNodesPerCell,1);
  mcIdType *outPtr=connOut->getPointer();
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  nbNodesPerCell++;
  for(mcIdType i=0;i<nbCells;i++,connI++)
    {
      if(conn[connI[0]]==typi && connI[1]-connI[0]==nbNodesPerCell)
        outPtr=std::copy(conn+connI[0]+1,conn+connI[1],outPtr);
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::convertNodalConnectivityToStaticGeoTypeMesh : there something wrong in cell #" << i << " ! The type of cell is not those expected, or the length of nodal connectivity is not those expected (" << nbNodesPerCell-1 << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  return connOut.retn();
}

/*!
 * Convert the nodal connectivity of the mesh so that all the cells are of dynamic types (polygon or quadratic
 * polygon). This returns the corresponding new nodal connectivity in \ref numbering-indirect format.
 * \param nodalConn nodal connectivity
 * \param nodalConnIndex nodal connectivity indices
 */
void MEDCouplingUMesh::convertNodalConnectivityToDynamicGeoTypeMesh(DataArrayIdType *&nodalConn, DataArrayIdType *&nodalConnIndex) const
{
  static const char msg0[]="MEDCouplingUMesh::convertNodalConnectivityToDynamicGeoTypeMesh : nodal connectivity in this are invalid ! Call checkConsistency !";
  checkConnectivityFullyDefined();
  if(_types.size()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::convertNodalConnectivityToDynamicGeoTypeMesh : current mesh does not contain exactly one geometric type !");
  mcIdType nbCells=getNumberOfCells(),
           lgth=_nodal_connec->getNumberOfTuples();
  if(lgth<nbCells)
    throw INTERP_KERNEL::Exception(msg0);
  MCAuto<DataArrayIdType> c(DataArrayIdType::New()),ci(DataArrayIdType::New());
  c->alloc(lgth-nbCells,1); ci->alloc(nbCells+1,1);
  mcIdType *cp(c->getPointer()),*cip(ci->getPointer());
  const mcIdType *incp(_nodal_connec->begin()),*incip(_nodal_connec_index->begin());
  cip[0]=0;
  for(mcIdType i=0;i<nbCells;i++,cip++,incip++)
    {
      mcIdType strt(incip[0]+1),stop(incip[1]);//+1 to skip geo type
      mcIdType delta(stop-strt);
      if(delta>=1)
        {
          if((strt>=0 && strt<lgth) && (stop>=0 && stop<=lgth))
            cp=std::copy(incp+strt,incp+stop,cp);
          else
            throw INTERP_KERNEL::Exception(msg0);
        }
      else
        throw INTERP_KERNEL::Exception(msg0);
      cip[1]=cip[0]+delta;
    }
  nodalConn=c.retn(); nodalConnIndex=ci.retn();
}

/*!
 * This method takes in input a vector of MEDCouplingUMesh instances lying on the same coordinates with same mesh dimensions.
 * Each mesh in \b ms must be sorted by type with the same order (typically using MEDCouplingUMesh::sortCellsInMEDFileFrmt).
 * This method is particularly useful for MED file interaction. It allows to aggregate several meshes and keeping the type sorting
 * and the track of the permutation by chunk of same geotype cells to retrieve it. The traditional formats old2new and new2old
 * are not used here to avoid the build of big permutation array.
 *
 * \param [in] ms meshes with same mesh dimension lying on the same coords and sorted by type following de the same geometric type order than
 *                those specified in MEDCouplingUMesh::sortCellsInMEDFileFrmt method.
 * \param [out] szOfCellGrpOfSameType is a newly allocated DataArrayIdType instance whose number of tuples is equal to the number of chunks of same geotype
 *              in all meshes in \b ms. The accumulation of all values of this array is equal to the number of cells of returned mesh.
 * \param [out] idInMsOfCellGrpOfSameType is a newly allocated DataArrayIdType instance having the same size than \b szOfCellGrpOfSameType. This
 *              output array gives for each chunck of same type the corresponding mesh id in \b ms.
 * \return A newly allocated unstructured mesh that is the result of the aggregation on same coords of all meshes in \b ms. This returned mesh
 *         is sorted by type following the geo cell types order of MEDCouplingUMesh::sortCellsInMEDFileFrmt method.
 */
MEDCouplingUMesh *MEDCouplingUMesh::AggregateSortedByTypeMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& ms,
                                                                            DataArrayIdType *&szOfCellGrpOfSameType,
                                                                            DataArrayIdType *&idInMsOfCellGrpOfSameType)
{
  std::vector<const MEDCouplingUMesh *> ms2;
  for(std::vector<const MEDCouplingUMesh *>::const_iterator it=ms.begin();it!=ms.end();it++)
    if(*it)
      {
        (*it)->checkConnectivityFullyDefined();
        ms2.push_back(*it);
      }
  if(ms2.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::AggregateSortedByTypeMeshesOnSameCoords : input vector is empty !");
  const DataArrayDouble *refCoo=ms2[0]->getCoords();
  int meshDim=ms2[0]->getMeshDimension();
  std::vector<const MEDCouplingUMesh *> m1ssm;
  std::vector< MCAuto<MEDCouplingUMesh> > m1ssmAuto;
  //
  std::vector<const MEDCouplingUMesh *> m1ssmSingle;
  std::vector< MCAuto<MEDCouplingUMesh> > m1ssmSingleAuto;
  mcIdType fake=0,rk=0;
  MCAuto<DataArrayIdType> ret1(DataArrayIdType::New()),ret2(DataArrayIdType::New());
  ret1->alloc(0,1); ret2->alloc(0,1);
  for(std::vector<const MEDCouplingUMesh *>::const_iterator it=ms2.begin();it!=ms2.end();it++,rk++)
    {
      if(meshDim!=(*it)->getMeshDimension())
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::AggregateSortedByTypeMeshesOnSameCoords : meshdims mismatch !");
      if(refCoo!=(*it)->getCoords())
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::AggregateSortedByTypeMeshesOnSameCoords : meshes are not shared by a single coordinates coords !");
      std::vector<MEDCouplingUMesh *> sp=(*it)->splitByType();
      std::copy(sp.begin(),sp.end(),std::back_insert_iterator< std::vector<const MEDCouplingUMesh *> >(m1ssm));
      std::copy(sp.begin(),sp.end(),std::back_insert_iterator< std::vector<MCAuto<MEDCouplingUMesh> > >(m1ssmAuto));
      for(std::vector<MEDCouplingUMesh *>::const_iterator it2=sp.begin();it2!=sp.end();it2++)
        {
          MEDCouplingUMesh *singleCell=static_cast<MEDCouplingUMesh *>((*it2)->buildPartOfMySelf(&fake,&fake+1,true));
          m1ssmSingleAuto.push_back(singleCell);
          m1ssmSingle.push_back(singleCell);
          ret1->pushBackSilent((*it2)->getNumberOfCells()); ret2->pushBackSilent(rk);
        }
    }
  MCAuto<MEDCouplingUMesh> m1ssmSingle2=MEDCouplingUMesh::MergeUMeshesOnSameCoords(m1ssmSingle);
  MCAuto<DataArrayIdType> renum=m1ssmSingle2->sortCellsInMEDFileFrmt();
  std::vector<const MEDCouplingUMesh *> m1ssmfinal(m1ssm.size());
  for(mcIdType i=0;i<ToIdType(m1ssm.size());i++)
    m1ssmfinal[renum->getIJ(i,0)]=m1ssm[i];
  MCAuto<MEDCouplingUMesh> ret0=MEDCouplingUMesh::MergeUMeshesOnSameCoords(m1ssmfinal);
  szOfCellGrpOfSameType=ret1->renumber(renum->begin());
  idInMsOfCellGrpOfSameType=ret2->renumber(renum->begin());
  return ret0.retn();
}

/*!
 * This method returns a newly created DataArrayIdType instance.
 * This method retrieves cell ids in [ \a begin, \a end ) that have the type \a type.
 */
DataArrayIdType *MEDCouplingUMesh::keepCellIdsByType(INTERP_KERNEL::NormalizedCellType type, const mcIdType *begin, const mcIdType *end) const
{
  checkFullyDefined();
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connIndex=_nodal_connec_index->begin();
  MCAuto<DataArrayIdType> ret(DataArrayIdType::New()); ret->alloc(0,1);
  for(const mcIdType *w=begin;w!=end;w++)
    if((INTERP_KERNEL::NormalizedCellType)conn[connIndex[*w]]==type)
      ret->pushBackSilent(*w);
  return ret.retn();
}

/*!
 * This method makes the assumption that da->getNumberOfTuples()<this->getNumberOfCells(). This method makes the assumption that ids contained in 'da'
 * are in [0:getNumberOfCells())
 */
DataArrayIdType *MEDCouplingUMesh::convertCellArrayPerGeoType(const DataArrayIdType *da) const
{
  checkFullyDefined();
  const mcIdType *conn=_nodal_connec->begin();
  const mcIdType *connI=_nodal_connec_index->begin();
  mcIdType nbOfCells=getNumberOfCells();
  std::set<INTERP_KERNEL::NormalizedCellType> types(getAllGeoTypes());
  mcIdType *tmp=new mcIdType[nbOfCells];
  for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iter=types.begin();iter!=types.end();iter++)
    {
      mcIdType j=0;
      for(const mcIdType *i=connI;i!=connI+nbOfCells;i++)
        if((INTERP_KERNEL::NormalizedCellType)conn[*i]==(*iter))
          tmp[std::distance(connI,i)]=j++;
    }
  DataArrayIdType *ret=DataArrayIdType::New();
  ret->alloc(da->getNumberOfTuples(),da->getNumberOfComponents());
  ret->copyStringInfoFrom(*da);
  mcIdType *retPtr=ret->getPointer();
  const mcIdType *daPtr=da->begin();
  mcIdType nbOfElems=da->getNbOfElems();
  for(mcIdType k=0;k<nbOfElems;k++)
    retPtr[k]=tmp[daPtr[k]];
  delete [] tmp;
  return ret;
}

/*!
 * This method reduced number of cells of this by keeping cells whose type is different from 'type' and if type=='type'
 * This method \b works \b for mesh sorted by type.
 * cells whose ids is in 'idsPerGeoType' array.
 * This method conserves coords and name of mesh.
 */
MEDCouplingUMesh *MEDCouplingUMesh::keepSpecifiedCells(INTERP_KERNEL::NormalizedCellType type, const mcIdType *idsPerGeoTypeBg, const mcIdType *idsPerGeoTypeEnd) const
{
  std::vector<mcIdType> code=getDistributionOfTypes();
  std::size_t nOfTypesInThis=code.size()/3;
  mcIdType sz=0,szOfType=0;
  for(std::size_t i=0;i<nOfTypesInThis;i++)
    {
      if(code[3*i]!=type)
        sz+=code[3*i+1];
      else
        szOfType=code[3*i+1];
    }
  for(const mcIdType *work=idsPerGeoTypeBg;work!=idsPerGeoTypeEnd;work++)
    if(*work<0 || *work>=szOfType)
      {
        std::ostringstream oss; oss << "MEDCouplingUMesh::keepSpecifiedCells : Request on type " << type << " at place #" << std::distance(idsPerGeoTypeBg,work) << " value " << *work;
        oss << ". It should be in [0," << szOfType << ") !";
        throw INTERP_KERNEL::Exception(oss.str());
      }
  MCAuto<DataArrayIdType> idsTokeep=DataArrayIdType::New(); idsTokeep->alloc(sz+std::distance(idsPerGeoTypeBg,idsPerGeoTypeEnd),1);
  mcIdType *idsPtr=idsTokeep->getPointer();
  mcIdType offset=0;
  for(std::size_t i=0;i<nOfTypesInThis;i++)
    {
      if(code[3*i]!=type)
        for(mcIdType j=0;j<code[3*i+1];j++)
          *idsPtr++=offset+j;
      else
        idsPtr=std::transform(idsPerGeoTypeBg,idsPerGeoTypeEnd,idsPtr,std::bind(std::plus<mcIdType>(),std::placeholders::_1,offset));
      offset+=code[3*i+1];
    }
  MCAuto<MEDCouplingUMesh> ret=static_cast<MEDCouplingUMesh *>(buildPartOfMySelf(idsTokeep->begin(),idsTokeep->end(),true));
  ret->copyTinyInfoFrom(this);
  return ret.retn();
}

/*!
 * This method returns a vector of size 'this->getNumberOfCells()'.
 * This method retrieves for each cell in \a this if it is linear (false) or quadratic(true).
 */
std::vector<bool> MEDCouplingUMesh::getQuadraticStatus() const
{
  mcIdType ncell=getNumberOfCells();
  std::vector<bool> ret(ncell);
  const mcIdType *cI=getNodalConnectivityIndex()->begin();
  const mcIdType *c=getNodalConnectivity()->begin();
  for(mcIdType i=0;i<ncell;i++)
    {
      INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType)c[cI[i]];
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(typ);
      ret[i]=cm.isQuadratic();
    }
  return ret;
}

/*!
 * Returns a newly created mesh (with ref count ==1) that contains merge of \a this and \a other.
 */
MEDCouplingMesh *MEDCouplingUMesh::mergeMyselfWith(const MEDCouplingMesh *other) const
{
  if(other->getType()!=UNSTRUCTURED)
    throw INTERP_KERNEL::Exception("Merge of umesh only available with umesh each other !");
  const MEDCouplingUMesh *otherC=static_cast<const MEDCouplingUMesh *>(other);
  return MergeUMeshes(this,otherC);
}

/*!
 * Returns a new DataArrayDouble holding barycenters of all cells. The barycenter is
 * computed by averaging coordinates of cell nodes, so this method is not a right
 * choice for degenerated meshes (not well oriented, cells with measure close to zero).
 * Beware also that for quadratic meshes, degenerated arc of circles are turned into linear edges for the computation.
 * This happens with a default detection precision of eps=1.0e-14. If you need control over this use computeCellCenterOfMassWithPrecision().
 *  \return DataArrayDouble * - a new instance of DataArrayDouble, of size \a
 *          this->getNumberOfCells() tuples per \a this->getSpaceDimension()
 *          components. The caller is to delete this array using decrRef() as it is
 *          no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \sa MEDCouplingUMesh::computeIsoBarycenterOfNodesPerCell
 *  \sa MEDCouplingUMesh::computeCellCenterOfMassWithPrecision
 */
DataArrayDouble *MEDCouplingUMesh::computeCellCenterOfMass() const
{
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  int spaceDim=getSpaceDimension();
  mcIdType nbOfCells=getNumberOfCells();
  ret->alloc(nbOfCells,spaceDim);
  ret->copyStringInfoFrom(*getCoords());
  double *ptToFill=ret->getPointer();
  const mcIdType *nodal=_nodal_connec->begin();
  const mcIdType *nodalI=_nodal_connec_index->begin();
  const double *coor=_coords->begin();
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)nodal[nodalI[i]];
      INTERP_KERNEL::computeBarycenter2<mcIdType,INTERP_KERNEL::ALL_C_MODE>(type,nodal+nodalI[i]+1,nodalI[i+1]-nodalI[i]-1,coor,spaceDim,ptToFill);
      ptToFill+=spaceDim;
    }
  return ret.retn();
}


/*!
 * See computeCellCenterOfMass().
 *  \param eps a precision for the detection of degenerated arc of circles.
 *  \return DataArrayDouble * - a new instance of DataArrayDouble, of size \a
 *          this->getNumberOfCells() tuples per \a this->getSpaceDimension()
 *          components. The caller is to delete this array using decrRef() as it is
 *          no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \sa MEDCouplingUMesh::computeIsoBarycenterOfNodesPerCell
 *  \sa MEDCouplingUMesh::computeCellCenterOfMassWithPrecision
 */
DataArrayDouble *MEDCouplingUMesh::computeCellCenterOfMassWithPrecision(double eps) const
{
  INTERP_KERNEL::QuadraticPlanarPrecision prec(eps);
  MCAuto<DataArrayDouble> ret = computeCellCenterOfMass();
  return ret.retn();
}


/*!
 * This method computes for each cell in \a this, the location of the iso barycenter of nodes constituting
 * the cell. Contrary to badly named MEDCouplingUMesh::computeCellCenterOfMass method that returns the center of inertia of the
 *
 * \return a newly allocated DataArrayDouble instance that the caller has to deal with. The returned
 *          DataArrayDouble instance will have \c this->getNumberOfCells() tuples and \c this->getSpaceDimension() components.
 *
 * \sa MEDCouplingUMesh::computeCellCenterOfMass
 * \throw If \a this is not fully defined (coordinates and connectivity)
 * \throw If there is presence in nodal connectivity in \a this of node ids not in [0, \c this->getNumberOfNodes() )
 */
DataArrayDouble *MEDCouplingUMesh::computeIsoBarycenterOfNodesPerCell() const
{
  checkFullyDefined();
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  int spaceDim=getSpaceDimension();
  mcIdType nbOfCells=getNumberOfCells();
  mcIdType nbOfNodes=getNumberOfNodes();
  ret->alloc(nbOfCells,spaceDim);
  double *ptToFill=ret->getPointer();
  const mcIdType *nodal=_nodal_connec->begin();
  const mcIdType *nodalI=_nodal_connec_index->begin();
  const double *coor=_coords->begin();
  for(mcIdType i=0;i<nbOfCells;i++,ptToFill+=spaceDim)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)nodal[nodalI[i]];
      std::fill(ptToFill,ptToFill+spaceDim,0.);
      if(type!=INTERP_KERNEL::NORM_POLYHED)
        {
          for(const mcIdType *conn=nodal+nodalI[i]+1;conn!=nodal+nodalI[i+1];conn++)
            {
              if(*conn>=0 && *conn<nbOfNodes)
                std::transform(coor+spaceDim*conn[0],coor+spaceDim*(conn[0]+1),ptToFill,ptToFill,std::plus<double>());
              else
                {
                  std::ostringstream oss; oss << "MEDCouplingUMesh::computeIsoBarycenterOfNodesPerCell : on cell #" << i << " presence of nodeId #" << *conn << " should be in [0," <<   nbOfNodes << ") !";
                  throw INTERP_KERNEL::Exception(oss.str());
                }
            }
          mcIdType nbOfNodesInCell=nodalI[i+1]-nodalI[i]-1;
          if(nbOfNodesInCell>0)
            std::transform(ptToFill,ptToFill+spaceDim,ptToFill,std::bind(std::multiplies<double>(),std::placeholders::_1,1./(double)nbOfNodesInCell));
          else
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::computeIsoBarycenterOfNodesPerCell : on cell #" << i << " presence of cell with no nodes !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
        }
      else
        {
          std::set<mcIdType> s(nodal+nodalI[i]+1,nodal+nodalI[i+1]);
          s.erase(-1);
          for(std::set<mcIdType>::const_iterator it=s.begin();it!=s.end();it++)
            {
              if(*it>=0 && *it<nbOfNodes)
                std::transform(coor+spaceDim*(*it),coor+spaceDim*((*it)+1),ptToFill,ptToFill,std::plus<double>());
              else
                {
                  std::ostringstream oss; oss << "MEDCouplingUMesh::computeIsoBarycenterOfNodesPerCell : on cell polyhedron cell #" << i << " presence of nodeId #" << *it << " should be in [0," <<   nbOfNodes << ") !";
                  throw INTERP_KERNEL::Exception(oss.str());
                }
            }
          if(!s.empty())
            std::transform(ptToFill,ptToFill+spaceDim,ptToFill,std::bind(std::multiplies<double>(),std::placeholders::_1,1./(double)s.size()));
          else
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::computeIsoBarycenterOfNodesPerCell : on polyhedron cell #" << i << " there are no nodes !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
        }
    }
  return ret.retn();
}

/*!
 * Returns a new DataArrayDouble holding barycenters of specified cells. The
 * barycenter is computed by averaging coordinates of cell nodes. The cells to treat
 * are specified via an array of cell ids.
 *  \warning Validity of the specified cell ids is not checked!
 *           Valid range is [ 0, \a this->getNumberOfCells() ).
 *  \param [in] begin - an array of cell ids of interest.
 *  \param [in] end - the end of \a begin, i.e. a pointer to its (last+1)-th element.
 *  \return DataArrayDouble * - a new instance of DataArrayDouble, of size ( \a
 *          end - \a begin ) tuples per \a this->getSpaceDimension() components. The
 *          caller is to delete this array using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_getPartBarycenterAndOwner "Here is a C++ example".<br>
 *  \ref  py_mcumesh_getPartBarycenterAndOwner "Here is a Python example".
 *  \endif
 */
DataArrayDouble *MEDCouplingUMesh::getPartBarycenterAndOwner(const mcIdType *begin, const mcIdType *end) const
{
  DataArrayDouble *ret=DataArrayDouble::New();
  int spaceDim=getSpaceDimension();
  std::size_t nbOfTuple=std::distance(begin,end);
  ret->alloc(nbOfTuple,spaceDim);
  double *ptToFill=ret->getPointer();
  double *tmp=new double[spaceDim];
  const mcIdType *nodal=_nodal_connec->begin();
  const mcIdType *nodalI=_nodal_connec_index->begin();
  const double *coor=_coords->begin();
  for(const mcIdType *w=begin;w!=end;w++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)nodal[nodalI[*w]];
      INTERP_KERNEL::computeBarycenter2<mcIdType,INTERP_KERNEL::ALL_C_MODE>(type,nodal+nodalI[*w]+1,nodalI[*w+1]-nodalI[*w]-1,coor,spaceDim,ptToFill);
      ptToFill+=spaceDim;
    }
  delete [] tmp;
  return ret;
}

/*!
 * Returns a DataArrayDouble instance giving for each cell in \a this the equation of plane given by "a*X+b*Y+c*Z+d=0".
 * So the returned instance will have 4 components and \c this->getNumberOfCells() tuples.
 * So this method expects that \a this has a spaceDimension equal to 3 and meshDimension equal to 2.
 * The computation of the plane equation is done using each time the 3 first nodes of 2D cells.
 * This method is useful to detect 2D cells in 3D space that are not coplanar.
 *
 * \return DataArrayDouble * - a new instance of DataArrayDouble having 4 components and a number of tuples equal to number of cells in \a this.
 * \throw If spaceDim!=3 or meshDim!=2.
 * \throw If connectivity of \a this is invalid.
 * \throw If connectivity of a cell in \a this points to an invalid node.
 */
DataArrayDouble *MEDCouplingUMesh::computePlaneEquationOf3DFaces() const
{
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  mcIdType nbOfCells=getNumberOfCells();
  mcIdType nbOfNodes(getNumberOfNodes());
  if(getSpaceDimension()!=3 || getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::computePlaneEquationOf3DFaces : This method must be applied on a mesh having meshDimension equal 2 and a spaceDimension equal to 3 !");
  ret->alloc(nbOfCells,4);
  double *retPtr(ret->getPointer());
  const mcIdType *nodal(_nodal_connec->begin()),*nodalI(_nodal_connec_index->begin());
  const double *coor(_coords->begin());
  for(mcIdType i=0;i<nbOfCells;i++,nodalI++,retPtr+=4)
    {
      double matrix[16]={0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,0},matrix2[16];
      if(nodalI[1]-nodalI[0]>=4)
        {
          double aa[3]={coor[nodal[nodalI[0]+1+1]*3+0]-coor[nodal[nodalI[0]+1+0]*3+0],
                        coor[nodal[nodalI[0]+1+1]*3+1]-coor[nodal[nodalI[0]+1+0]*3+1],
                        coor[nodal[nodalI[0]+1+1]*3+2]-coor[nodal[nodalI[0]+1+0]*3+2]}
          ,bb[3]={coor[nodal[nodalI[0]+1+2]*3+0]-coor[nodal[nodalI[0]+1+0]*3+0],
                        coor[nodal[nodalI[0]+1+2]*3+1]-coor[nodal[nodalI[0]+1+0]*3+1],
                        coor[nodal[nodalI[0]+1+2]*3+2]-coor[nodal[nodalI[0]+1+0]*3+2]};
          double cc[3]={aa[1]*bb[2]-aa[2]*bb[1],aa[2]*bb[0]-aa[0]*bb[2],aa[0]*bb[1]-aa[1]*bb[0]};
          double aa_norm(sqrt(aa[0]*aa[0]+aa[1]*aa[1]+aa[2]*aa[2])),bb_norm(sqrt(bb[0]*bb[0]+bb[1]*bb[1]+bb[2]*bb[2]));
          for(int j=0;j<3;j++)
            {
              mcIdType nodeId(nodal[nodalI[0]+1+j]);
              if(nodeId>=0 && nodeId<nbOfNodes)
                std::copy(coor+nodeId*3,coor+(nodeId+1)*3,matrix+4*j);
              else
                {
                  std::ostringstream oss; oss << "MEDCouplingUMesh::computePlaneEquationOf3DFaces : invalid 2D cell #" << i << " ! This cell points to an invalid nodeId : " << nodeId << " !";
                  throw INTERP_KERNEL::Exception(oss.str());
                }
            }
          if(sqrt(cc[0]*cc[0]+cc[1]*cc[1]+cc[2]*cc[2])>(1e-3*aa_norm*bb_norm))
            {
              INTERP_KERNEL::inverseMatrix(matrix,4,matrix2);
              retPtr[0]=matrix2[3]; retPtr[1]=matrix2[7]; retPtr[2]=matrix2[11]; retPtr[3]=matrix2[15];
            }
          else
            {
              if(nodalI[1]-nodalI[0]==4)
                {
                  std::ostringstream oss; oss << "MEDCouplingUMesh::computePlaneEquationOf3DFaces : cell" << i << " : Presence of The 3 colinear points !";
                  throw INTERP_KERNEL::Exception(oss.str());
                }
              //
              double dd[3]={0.,0.,0.};
              for(mcIdType offset=nodalI[0]+1;offset<nodalI[1];offset++)
                std::transform(coor+3*nodal[offset],coor+3*(nodal[offset]+1),dd,dd,std::plus<double>());
              mcIdType nbOfNodesInCell(nodalI[1]-nodalI[0]-1);
              std::transform(dd,dd+3,dd,std::bind(std::multiplies<double>(),std::placeholders::_1,1./(double)nbOfNodesInCell));
              std::copy(dd,dd+3,matrix+4*2);
              INTERP_KERNEL::inverseMatrix(matrix,4,matrix2);
              retPtr[0]=matrix2[3]; retPtr[1]=matrix2[7]; retPtr[2]=matrix2[11]; retPtr[3]=matrix2[15];
            }
        }
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::computePlaneEquationOf3DFaces : invalid 2D cell #" << i << " ! Must be constitued by more than 3 nodes !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  return ret.retn();
}

/*!
 * This method expects as input a DataArrayDouble non nul instance 'da' that should be allocated. If not an exception is thrown.
 *
 */
MEDCouplingUMesh *MEDCouplingUMesh::Build0DMeshFromCoords(DataArrayDouble *da)
{
  if(!da)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::Build0DMeshFromCoords : instance of DataArrayDouble must be not null !");
  da->checkAllocated();
  std::string name(da->getName());
  MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::New(name,0));
  if(name.empty())
    ret->setName("Mesh");
  ret->setCoords(da);
  mcIdType nbOfTuples(da->getNumberOfTuples());
  MCAuto<DataArrayIdType> c(DataArrayIdType::New()),cI(DataArrayIdType::New());
  c->alloc(2*nbOfTuples,1);
  cI->alloc(nbOfTuples+1,1);
  mcIdType *cp(c->getPointer()),*cip(cI->getPointer());
  *cip++=0;
  for(mcIdType i=0;i<nbOfTuples;i++)
    {
      *cp++=INTERP_KERNEL::NORM_POINT1;
      *cp++=i;
      *cip++=2*(i+1);
    }
  ret->setConnectivity(c,cI,true);
  return ret.retn();
}

MCAuto<MEDCouplingUMesh> MEDCouplingUMesh::Build1DMeshFromCoords(DataArrayDouble *da)
{
  if(!da)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::Build01MeshFromCoords : instance of DataArrayDouble must be not null !");
  da->checkAllocated();
  std::string name(da->getName());
  MCAuto<MEDCouplingUMesh> ret;
  {
    MCAuto<MEDCouplingCMesh> tmp(MEDCouplingCMesh::New());
    MCAuto<DataArrayDouble> arr(DataArrayDouble::New());
    arr->alloc(da->getNumberOfTuples());
    tmp->setCoordsAt(0,arr);
    ret=tmp->buildUnstructured();
  }
  ret->setCoords(da);
  if(name.empty())
    ret->setName("Mesh");
  else
    ret->setName(name);
  return ret;
}

/*!
 * Creates a new MEDCouplingUMesh by concatenating two given meshes of the same dimension.
 * Cells and nodes of
 * the first mesh precede cells and nodes of the second mesh within the result mesh.
 *  \param [in] mesh1 - the first mesh.
 *  \param [in] mesh2 - the second mesh.
 *  \return MEDCouplingUMesh * - the result mesh. It is a new instance of
 *          MEDCouplingUMesh. The caller is to delete this mesh using decrRef() as it
 *          is no more needed.
 *  \throw If \a mesh1 == NULL or \a mesh2 == NULL.
 *  \throw If the coordinates array is not set in none of the meshes.
 *  \throw If \a mesh1->getMeshDimension() < 0 or \a mesh2->getMeshDimension() < 0.
 *  \throw If \a mesh1->getMeshDimension() != \a mesh2->getMeshDimension().
 */
MEDCouplingUMesh *MEDCouplingUMesh::MergeUMeshes(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2)
{
  std::vector<const MEDCouplingUMesh *> tmp(2);
  tmp[0]=const_cast<MEDCouplingUMesh *>(mesh1); tmp[1]=const_cast<MEDCouplingUMesh *>(mesh2);
  return MergeUMeshes(tmp);
}

/*!
 * Creates a new MEDCouplingUMesh by concatenating all given meshes of the same dimension.
 * Cells and nodes of
 * the *i*-th mesh precede cells and nodes of the (*i*+1)-th mesh within the result mesh.
 *  \param [in] a - a vector of meshes (MEDCouplingUMesh) to concatenate.
 *  \return MEDCouplingUMesh * - the result mesh. It is a new instance of
 *          MEDCouplingUMesh. The caller is to delete this mesh using decrRef() as it
 *          is no more needed.
 *  \throw If \a a.size() == 0.
 *  \throw If \a a[ *i* ] == NULL.
 *  \throw If the coordinates array is not set in none of the meshes.
 *  \throw If \a a[ *i* ]->getMeshDimension() < 0.
 *  \throw If the meshes in \a a are of different dimension (getMeshDimension()).
 */
MEDCouplingUMesh *MEDCouplingUMesh::MergeUMeshes(const std::vector<const MEDCouplingUMesh *>& a)
{
  std::size_t sz=a.size();
  if(sz==0)
    return MergeUMeshesLL(a);
  for(std::size_t ii=0;ii<sz;ii++)
    if(!a[ii])
      {
        std::ostringstream oss; oss << "MEDCouplingUMesh::MergeUMeshes : item #" << ii << " in input array of size "<< sz << " is empty !";
        throw INTERP_KERNEL::Exception(oss.str());
      }
  std::vector< MCAuto<MEDCouplingUMesh> > bb(sz);
  std::vector< const MEDCouplingUMesh * > aa(sz);
  int spaceDim=-3;
  for(std::size_t i=0;i<sz && spaceDim==-3;i++)
    {
      const MEDCouplingUMesh *cur=a[i];
      const DataArrayDouble *coo=cur->getCoords();
      if(coo)
        spaceDim=int(coo->getNumberOfComponents());
    }
  if(spaceDim==-3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::MergeUMeshes : no spaceDim specified ! unable to perform merge !");
  for(std::size_t i=0;i<sz;i++)
    {
      bb[i]=a[i]->buildSetInstanceFromThis(spaceDim);
      aa[i]=bb[i];
    }
  return MergeUMeshesLL(aa);
}

/*!
 * Creates a new MEDCouplingUMesh by concatenating cells of two given meshes of same
 * dimension and sharing the node coordinates array.
 * All cells of the first mesh precede all cells of the second mesh
 * within the result mesh.
 *  \param [in] mesh1 - the first mesh.
 *  \param [in] mesh2 - the second mesh.
 *  \return MEDCouplingUMesh * - the result mesh. It is a new instance of
 *          MEDCouplingUMesh. The caller is to delete this mesh using decrRef() as it
 *          is no more needed.
 *  \throw If \a mesh1 == NULL or \a mesh2 == NULL.
 *  \throw If the meshes do not share the node coordinates array.
 *  \throw If \a mesh1->getMeshDimension() < 0 or \a mesh2->getMeshDimension() < 0.
 *  \throw If \a mesh1->getMeshDimension() != \a mesh2->getMeshDimension().
 */
MEDCouplingUMesh *MEDCouplingUMesh::MergeUMeshesOnSameCoords(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2)
{
  std::vector<const MEDCouplingUMesh *> tmp(2);
  tmp[0]=mesh1; tmp[1]=mesh2;
  return MergeUMeshesOnSameCoords(tmp);
}

/*!
 * Creates a new MEDCouplingUMesh by concatenating cells of all given meshes of same
 * dimension and sharing the node coordinates array.
 * All cells of the *i*-th mesh precede all cells of the
 * (*i*+1)-th mesh within the result mesh.
 *  \param [in] meshes - a vector of meshes (MEDCouplingUMesh) to concatenate.
 *  \return MEDCouplingUMesh * - the result mesh. It is a new instance of
 *          MEDCouplingUMesh. The caller is to delete this mesh using decrRef() as it
 *          is no more needed.
 *  \throw If \a a.size() == 0.
 *  \throw If \a a[ *i* ] == NULL.
 *  \throw If the meshes do not share the node coordinates array.
 *  \throw If \a a[ *i* ]->getMeshDimension() < 0.
 *  \throw If the meshes in \a a are of different dimension (getMeshDimension()).
 */
MEDCouplingUMesh *MEDCouplingUMesh::MergeUMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& meshes)
{
  if(meshes.empty())
    throw INTERP_KERNEL::Exception("meshes input parameter is expected to be non empty.");
  for(std::size_t ii=0;ii<meshes.size();ii++)
    if(!meshes[ii])
      {
        std::ostringstream oss; oss << "MEDCouplingUMesh::MergeUMeshesOnSameCoords : item #" << ii << " in input array of size "<< meshes.size() << " is empty !";
        throw INTERP_KERNEL::Exception(oss.str());
      }
  const DataArrayDouble *coords=meshes.front()->getCoords();
  int meshDim=meshes.front()->getMeshDimension();
  std::vector<const MEDCouplingUMesh *>::const_iterator iter=meshes.begin();
  mcIdType meshLgth=0;
  mcIdType meshIndexLgth=0;
  for(;iter!=meshes.end();iter++)
    {
      if(coords!=(*iter)->getCoords())
        throw INTERP_KERNEL::Exception("meshes does not share the same coords ! Try using tryToShareSameCoords method !");
      if(meshDim!=(*iter)->getMeshDimension())
        throw INTERP_KERNEL::Exception("Mesh dimensions mismatches, FuseUMeshesOnSameCoords impossible !");
      meshLgth+=(*iter)->getNodalConnectivityArrayLen();
      meshIndexLgth+=(*iter)->getNumberOfCells();
    }
  MCAuto<DataArrayIdType> nodal=DataArrayIdType::New();
  nodal->alloc(meshLgth,1);
  mcIdType *nodalPtr=nodal->getPointer();
  MCAuto<DataArrayIdType> nodalIndex=DataArrayIdType::New();
  nodalIndex->alloc(meshIndexLgth+1,1);
  mcIdType *nodalIndexPtr=nodalIndex->getPointer();
  mcIdType offset=0;
  for(iter=meshes.begin();iter!=meshes.end();iter++)
    {
      const mcIdType *nod=(*iter)->getNodalConnectivity()->begin();
      const mcIdType *index=(*iter)->getNodalConnectivityIndex()->begin();
      mcIdType nbOfCells=(*iter)->getNumberOfCells();
      mcIdType meshLgth2=(*iter)->getNodalConnectivityArrayLen();
      nodalPtr=std::copy(nod,nod+meshLgth2,nodalPtr);
      if(iter!=meshes.begin())
        nodalIndexPtr=std::transform(index+1,index+nbOfCells+1,nodalIndexPtr,std::bind(std::plus<mcIdType>(),std::placeholders::_1,offset));
      else
        nodalIndexPtr=std::copy(index,index+nbOfCells+1,nodalIndexPtr);
      offset+=meshLgth2;
    }
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New();
  ret->setName("merge");
  ret->setMeshDimension(meshDim);
  ret->setConnectivity(nodal,nodalIndex,true);
  ret->setCoords(coords);
  return ret;
}

/*!
 * Creates a new MEDCouplingUMesh by concatenating cells of all given meshes of same
 * dimension and sharing the node coordinates array. Cells of the *i*-th mesh precede
 * cells of the (*i*+1)-th mesh within the result mesh. Duplicates of cells are
 * removed from \a this mesh and arrays mapping between new and old cell ids in "Old to
 * New" mode are returned for each input mesh.
 *  \param [in] meshes - a vector of meshes (MEDCouplingUMesh) to concatenate.
 *  \param [in] compType - specifies a cell comparison technique. For meaning of its
 *          valid values [0,1,2], see zipConnectivityTraducer().
 *  \param [in,out] corr - an array of DataArrayIdType, of the same size as \a
 *          meshes. The *i*-th array describes cell ids mapping for \a meshes[ *i* ]
 *          mesh. The caller is to delete each of the arrays using decrRef() as it is
 *          no more needed.
 *  \return MEDCouplingUMesh * - the result mesh. It is a new instance of
 *          MEDCouplingUMesh. The caller is to delete this mesh using decrRef() as it
 *          is no more needed.
 *  \throw If \a meshes.size() == 0.
 *  \throw If \a meshes[ *i* ] == NULL.
 *  \throw If the meshes do not share the node coordinates array.
 *  \throw If \a meshes[ *i* ]->getMeshDimension() < 0.
 *  \throw If the \a meshes are of different dimension (getMeshDimension()).
 *  \throw If the nodal connectivity of cells of any of \a meshes is not defined.
 *  \throw If the nodal connectivity any of \a meshes includes an invalid id.
 */
MEDCouplingUMesh *MEDCouplingUMesh::FuseUMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& meshes, int compType, std::vector<DataArrayIdType *>& corr)
{
  //All checks are delegated to MergeUMeshesOnSameCoords
  MCAuto<MEDCouplingUMesh> ret=MergeUMeshesOnSameCoords(meshes);
  MCAuto<DataArrayIdType> o2n=ret->zipConnectivityTraducer(compType);
  corr.resize(meshes.size());
  std::size_t nbOfMeshes=meshes.size();
  mcIdType offset=0;
  const mcIdType *o2nPtr=o2n->begin();
  for(std::size_t i=0;i<nbOfMeshes;i++)
    {
      DataArrayIdType *tmp=DataArrayIdType::New();
      mcIdType curNbOfCells=meshes[i]->getNumberOfCells();
      tmp->alloc(curNbOfCells,1);
      std::copy(o2nPtr+offset,o2nPtr+offset+curNbOfCells,tmp->getPointer());
      offset+=curNbOfCells;
      tmp->setName(meshes[i]->getName());
      corr[i]=tmp;
    }
  return ret.retn();
}

/*!
 * Makes all given meshes share the nodal connectivity array. The common connectivity
 * array is created by concatenating the connectivity arrays of all given meshes. All
 * the given meshes must be of the same space dimension but dimension of cells **can
 * differ**. This method is particularly useful in MEDLoader context to build a \ref
 * MEDCoupling::MEDFileUMesh "MEDFileUMesh" instance that expects that underlying
 * MEDCouplingUMesh'es of different dimension share the same nodal connectivity array.
 *  \param [in,out] meshes - a vector of meshes to update.
 *  \throw If any of \a meshes is NULL.
 *  \throw If the coordinates array is not set in any of \a meshes.
 *  \throw If the nodal connectivity of cells is not defined in any of \a meshes.
 *  \throw If \a meshes are of different space dimension.
 */
void MEDCouplingUMesh::PutUMeshesOnSameAggregatedCoords(const std::vector<MEDCouplingUMesh *>& meshes)
{
  std::size_t sz=meshes.size();
  if(sz==0 || sz==1)
    return;
  std::vector< const DataArrayDouble * > coords(meshes.size());
  std::vector< const DataArrayDouble * >::iterator it2=coords.begin();
  for(std::vector<MEDCouplingUMesh *>::const_iterator it=meshes.begin();it!=meshes.end();it++,it2++)
    {
      if((*it))
        {
          (*it)->checkConnectivityFullyDefined();
          const DataArrayDouble *coo=(*it)->getCoords();
          if(coo)
            *it2=coo;
          else
            {
              std::ostringstream oss; oss << " MEDCouplingUMesh::PutUMeshesOnSameAggregatedCoords : Item #" << std::distance(meshes.begin(),it) << " inside the vector of length " << meshes.size();
              oss << " has no coordinate array defined !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
        }
      else
        {
          std::ostringstream oss; oss << " MEDCouplingUMesh::PutUMeshesOnSameAggregatedCoords : Item #" << std::distance(meshes.begin(),it) << " inside the vector of length " << meshes.size();
          oss << " is null !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  MCAuto<DataArrayDouble> res=DataArrayDouble::Aggregate(coords);
  std::vector<MEDCouplingUMesh *>::const_iterator it=meshes.begin();
  mcIdType offset=(*it)->getNumberOfNodes();
  (*it++)->setCoords(res);
  for(;it!=meshes.end();it++)
    {
      mcIdType oldNumberOfNodes=(*it)->getNumberOfNodes();
      (*it)->setCoords(res);
      (*it)->shiftNodeNumbersInConn(offset);
      offset+=oldNumberOfNodes;
    }
}

/*!
 * Merges nodes coincident with a given precision within all given meshes that share
 * the nodal connectivity array. The given meshes **can be of different** mesh
 * dimension. This method is particularly useful in MEDLoader context to build a \ref
 * MEDCoupling::MEDFileUMesh "MEDFileUMesh" instance that expects that underlying
 * MEDCouplingUMesh'es of different dimension share the same nodal connectivity array.
 *  \param [in,out] meshes - a vector of meshes to update.
 *  \param [in] eps - the precision used to detect coincident nodes (infinite norm).
 *  \throw If any of \a meshes is NULL.
 *  \throw If the \a meshes do not share the same node coordinates array.
 *  \throw If the nodal connectivity of cells is not defined in any of \a meshes.
 */
void MEDCouplingUMesh::MergeNodesOnUMeshesSharingSameCoords(const std::vector<MEDCouplingUMesh *>& meshes, double eps)
{
  if(meshes.empty())
    return ;
  std::set<const DataArrayDouble *> s;
  for(std::vector<MEDCouplingUMesh *>::const_iterator it=meshes.begin();it!=meshes.end();it++)
    {
      if(*it)
        s.insert((*it)->getCoords());
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::MergeNodesOnUMeshesSharingSameCoords : In input vector of unstructured meshes of size " << meshes.size() << " the element #" << std::distance(meshes.begin(),it) << " is null !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  if(s.size()!=1)
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::MergeNodesOnUMeshesSharingSameCoords : In input vector of unstructured meshes of size " << meshes.size() << ", it appears that they do not share the same instance of DataArrayDouble for coordiantes ! tryToShareSameCoordsPermute method can help to reach that !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  const DataArrayDouble *coo=*(s.begin());
  if(!coo)
    return;
  //
  DataArrayIdType *comm,*commI;
  coo->findCommonTuples(eps,-1,comm,commI);
  MCAuto<DataArrayIdType> tmp1(comm),tmp2(commI);
  mcIdType oldNbOfNodes=coo->getNumberOfTuples();
  mcIdType newNbOfNodes;
  MCAuto<DataArrayIdType> o2n=DataArrayIdType::ConvertIndexArrayToO2N(oldNbOfNodes,comm->begin(),commI->begin(),commI->end(),newNbOfNodes);
  if(oldNbOfNodes==newNbOfNodes)
    return ;
  MCAuto<DataArrayDouble> newCoords=coo->renumberAndReduce(o2n->begin(),newNbOfNodes);
  for(std::vector<MEDCouplingUMesh *>::const_iterator it=meshes.begin();it!=meshes.end();it++)
    {
      (*it)->renumberNodesInConn(o2n->begin());
      (*it)->setCoords(newCoords);
    }
}


/*!
 * This static operates only for coords in 3D. The polygon is specified by its connectivity nodes in [ \a begin , \a end ).
 */
bool MEDCouplingUMesh::IsPolygonWellOriented(bool isQuadratic, const double *vec, const mcIdType *begin, const mcIdType *end, const double *coords)
{
  std::size_t i, ip1;
  double v[3]={0.,0.,0.};
  std::size_t sz=std::distance(begin,end);
  if(!isQuadratic)
    for(i=0;i<sz;i++)
      {
        // Algorithm: sum in v the cross products of (e1, e2) where e_i it the vector between (0,0,0) and point i
        // and e2 is linear point directly following e1 in the connectivity. All points are used.
        v[0]+=coords[3*begin[i]+1]*coords[3*begin[(i+1)%sz]+2]-coords[3*begin[i]+2]*coords[3*begin[(i+1)%sz]+1];
        v[1]+=coords[3*begin[i]+2]*coords[3*begin[(i+1)%sz]]-coords[3*begin[i]]*coords[3*begin[(i+1)%sz]+2];
        v[2]+=coords[3*begin[i]]*coords[3*begin[(i+1)%sz]+1]-coords[3*begin[i]+1]*coords[3*begin[(i+1)%sz]];
      }
  else
    {
      // Same algorithm as above but also using intermediate quadratic points.
      // (taking only linear points might lead to issues if the linearized version of the
      // polygon is not convex or self-intersecting ... see testCellOrientation4)
      std::size_t hsz = sz/2;
      for(std::size_t j=0;j<sz;j++)
        {
          if (j%2)  // current point i is quadratic, next point i+1 is standard
            {
              i = hsz+(j-1)/2;
              ip1 = ((j-1)/2 + 1)%hsz; // ip1 means "i+1", i.e. next point
            }
          else      // current point i is standard, next point i+1 is quadratic
            {
              i = j/2;
              ip1 = j/2+hsz;
            }
          v[0]+=coords[3*begin[i]+1]*coords[3*begin[ip1]+2]-coords[3*begin[i]+2]*coords[3*begin[ip1]+1];
          v[1]+=coords[3*begin[i]+2]*coords[3*begin[ip1]]-coords[3*begin[i]]*coords[3*begin[ip1]+2];
          v[2]+=coords[3*begin[i]]*coords[3*begin[ip1]+1]-coords[3*begin[i]+1]*coords[3*begin[ip1]];
        }
    }
  double ret = vec[0]*v[0]+vec[1]*v[1]+vec[2]*v[2];
  return (ret>0.);
}

/*!
 * The polyhedron is specified by its connectivity nodes in [ \a begin , \a end ).
 */
bool MEDCouplingUMesh::IsPolyhedronWellOriented(const mcIdType *begin, const mcIdType *end, const double *coords)
{
  std::vector<std::pair<mcIdType,mcIdType> > edges;
  std::size_t nbOfFaces=std::count(begin,end,-1)+1;
  const mcIdType *bgFace=begin;
  for(std::size_t i=0;i<nbOfFaces;i++)
    {
      const mcIdType *endFace=std::find(bgFace+1,end,-1);
      std::size_t nbOfEdgesInFace=std::distance(bgFace,endFace);
      for(std::size_t j=0;j<nbOfEdgesInFace;j++)
        {
          std::pair<mcIdType,mcIdType> p1(bgFace[j],bgFace[(j+1)%nbOfEdgesInFace]);
          if(std::find(edges.begin(),edges.end(),p1)!=edges.end())
            return false;
          edges.push_back(p1);
        }
      bgFace=endFace+1;
    }
  return INTERP_KERNEL::calculateVolumeForPolyh2<mcIdType,INTERP_KERNEL::ALL_C_MODE>(begin,ToIdType(std::distance(begin,end)),coords)>-EPS_FOR_POLYH_ORIENTATION;
}

/*!
 * The 3D extruded static cell (PENTA6,HEXA8,HEXAGP12...) its connectivity nodes in [ \a begin , \a end ).
 */
bool MEDCouplingUMesh::Is3DExtrudedStaticCellWellOriented(const mcIdType *begin, const mcIdType *end, const double *coords)
{
  double vec0[3],vec1[3];
  std::size_t sz=std::distance(begin,end);
  if(sz%2!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::Is3DExtrudedStaticCellWellOriented : the length of nodal connectivity of extruded cell is not even !");
  mcIdType nbOfNodes=ToIdType(sz/2);
  INTERP_KERNEL::areaVectorOfPolygon<mcIdType,INTERP_KERNEL::ALL_C_MODE>(begin,nbOfNodes,coords,vec0);
  const double *pt0=coords+3*begin[0];
  const double *pt1=coords+3*begin[nbOfNodes];
  vec1[0]=pt1[0]-pt0[0]; vec1[1]=pt1[1]-pt0[1]; vec1[2]=pt1[2]-pt0[2];
  return (vec0[0]*vec1[0]+vec0[1]*vec1[1]+vec0[2]*vec1[2])<0.;
}

void MEDCouplingUMesh::CorrectExtrudedStaticCell(mcIdType *begin, mcIdType *end)
{
  std::size_t sz=std::distance(begin,end);
  INTERP_KERNEL::AutoPtr<mcIdType> tmp=new mcIdType[sz];
  std::size_t nbOfNodes(sz/2);
  std::copy(begin,end,(mcIdType *)tmp);
  for(std::size_t j=1;j<nbOfNodes;j++)
    {
      begin[j]=tmp[nbOfNodes-j];
      begin[j+nbOfNodes]=tmp[nbOfNodes+nbOfNodes-j];
    }
}

bool MEDCouplingUMesh::IsTetra4WellOriented(const mcIdType *begin, const mcIdType *end, const double *coords)
{
  std::size_t sz=std::distance(begin,end);
  if(sz!=4)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::IsTetra4WellOriented : Tetra4 cell with not 4 nodes ! Call checkConsistency !");
  double vec0[3],vec1[3];
  const double *pt0=coords+3*begin[0],*pt1=coords+3*begin[1],*pt2=coords+3*begin[2],*pt3=coords+3*begin[3];
  vec0[0]=pt1[0]-pt0[0]; vec0[1]=pt1[1]-pt0[1]; vec0[2]=pt1[2]-pt0[2]; vec1[0]=pt2[0]-pt0[0]; vec1[1]=pt2[1]-pt0[1]; vec1[2]=pt2[2]-pt0[2];
  return ((vec0[1]*vec1[2]-vec0[2]*vec1[1])*(pt3[0]-pt0[0])+(vec0[2]*vec1[0]-vec0[0]*vec1[2])*(pt3[1]-pt0[1])+(vec0[0]*vec1[1]-vec0[1]*vec1[0])*(pt3[2]-pt0[2]))<0;
}

bool MEDCouplingUMesh::IsPyra5WellOriented(const mcIdType *begin, const mcIdType *end, const double *coords)
{
  std::size_t sz=std::distance(begin,end);
  if(sz!=5)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::IsPyra5WellOriented : Pyra5 cell with not 5 nodes ! Call checkConsistency !");
  double vec0[3];
  INTERP_KERNEL::areaVectorOfPolygon<mcIdType,INTERP_KERNEL::ALL_C_MODE>(begin,4,coords,vec0);
  const double *pt0=coords+3*begin[0],*pt1=coords+3*begin[4];
  return (vec0[0]*(pt1[0]-pt0[0])+vec0[1]*(pt1[1]-pt0[1])+vec0[2]*(pt1[2]-pt0[2]))<0.;
}

/*!
 * This method performs a simplyfication of a single polyedron cell. To do that each face of cell whose connectivity is defined by [ \b begin , \b end )
 * is compared with the others in order to find faces in the same plane (with approx of eps). If any, the cells are grouped together and projected to
 * a 2D space.
 *
 * \param [in] eps is a relative precision that allows to establish if some 3D plane are coplanar or not.
 * \param [in] coords the coordinates with nb of components exactly equal to 3
 * \param [in] index begin of the nodal connectivity (geometric type included) of a single polyhedron cell
 * \param [out] res the result is put at the end of the vector without any alteration of the data.
 */
void MEDCouplingUMesh::SimplifyPolyhedronCell(double eps, const DataArrayDouble *coords, mcIdType index, DataArrayIdType *res, MEDCouplingUMesh *faces,
                                              DataArrayIdType *E_Fi, DataArrayIdType *E_F, DataArrayIdType *F_Ei, DataArrayIdType *F_E)
{
  mcIdType nbFaces = E_Fi->getIJ(index + 1, 0) - E_Fi->getIJ(index, 0);
  MCAuto<DataArrayDouble> v=DataArrayDouble::New(); v->alloc(nbFaces,3);
  double *vPtr=v->getPointer();
  MCAuto<DataArrayDouble> p=DataArrayDouble::New(); p->alloc(nbFaces,2);
  double *pPtr=p->getPointer();
  mcIdType *e_fi = E_Fi->getPointer(), *e_f = E_F->getPointer(), *f_ei = F_Ei->getPointer(), *f_e = F_E->getPointer();
  const mcIdType *f_idx = faces->getNodalConnectivityIndex()->getPointer(), *f_cnn = faces->getNodalConnectivity()->getPointer();
  for(mcIdType i=0;i<nbFaces;i++,vPtr+=3,pPtr++)
    {
      mcIdType face = e_f[e_fi[index] + i];
      ComputeVecAndPtOfFace(eps, coords->begin(), f_cnn + f_idx[face] + 1, f_cnn + f_idx[face + 1], vPtr, pPtr);
      // to differentiate faces going to different cells:
      pPtr++, *pPtr = 0;
      for (mcIdType j = f_ei[face]; j < f_ei[face + 1]; j++)
        *pPtr += FromIdType<double>(f_e[j]);
    }
  pPtr=p->getPointer(); vPtr=v->getPointer();
  DataArrayIdType *comm1=0,*commI1=0;
  v->findCommonTuples(eps,-1,comm1,commI1);
  for (mcIdType i = 0; i < nbFaces; i++)
    if (comm1->findIdFirstEqual(i) < 0)
      {
        comm1->pushBackSilent(i);
        commI1->pushBackSilent(comm1->getNumberOfTuples());
      }
  MCAuto<DataArrayIdType> comm1Auto(comm1),commI1Auto(commI1);
  const mcIdType *comm1Ptr=comm1->begin();
  const mcIdType *commI1Ptr=commI1->begin();
  mcIdType nbOfGrps1=commI1Auto->getNumberOfTuples()-1;
  res->pushBackSilent(ToIdType(INTERP_KERNEL::NORM_POLYHED));
  //
  for(mcIdType i=0;i<nbOfGrps1;i++)
    {
      mcIdType vecId=comm1Ptr[commI1Ptr[i]];
      MCAuto<DataArrayDouble> tmpgrp2=p->selectByTupleId(comm1Ptr+commI1Ptr[i],comm1Ptr+commI1Ptr[i+1]);
      DataArrayIdType *comm2=0,*commI2=0;
      tmpgrp2->findCommonTuples(eps,-1,comm2,commI2);
      for (mcIdType j = 0; j < commI1Ptr[i+1] - commI1Ptr[i]; j++)
        if (comm2->findIdFirstEqual(j) < 0)
          {
            comm2->pushBackSilent(j);
            commI2->pushBackSilent(comm2->getNumberOfTuples());
          }
      MCAuto<DataArrayIdType> comm2Auto(comm2),commI2Auto(commI2);
      const mcIdType *comm2Ptr=comm2->begin();
      const mcIdType *commI2Ptr=commI2->begin();
      mcIdType nbOfGrps2=commI2Auto->getNumberOfTuples()-1;
      for(mcIdType j=0;j<nbOfGrps2;j++)
        {
          if(commI2Ptr[j+1] == commI2Ptr[j] + 1)
            {
              mcIdType face = e_f[e_fi[index] + comm1Ptr[commI1Ptr[i] + comm2Ptr[commI2Ptr[j]]]]; //hmmm
              res->insertAtTheEnd(f_cnn + f_idx[face] + 1, f_cnn + f_idx[face + 1]);
              res->pushBackSilent(-1);
            }
          else
            {
              mcIdType pointId=comm1Ptr[commI1Ptr[i]+comm2Ptr[commI2Ptr[j]]];
              MCAuto<DataArrayIdType> ids2=comm2->selectByTupleIdSafeSlice(commI2Ptr[j],commI2Ptr[j+1],1);
              ids2->transformWithIndArr(comm1Ptr+commI1Ptr[i],comm1Ptr+commI1Ptr[i+1]);
              ids2->transformWithIndArr(e_f + e_fi[index], e_f + e_fi[index + 1]);
              MCAuto<MEDCouplingUMesh> mm3=static_cast<MEDCouplingUMesh *>(faces->buildPartOfMySelf(ids2->begin(),ids2->end(),true));
              MCAuto<DataArrayIdType> idsNodeTmp=mm3->zipCoordsTraducer();
              MCAuto<DataArrayIdType> idsNode=idsNodeTmp->invertArrayO2N2N2O(mm3->getNumberOfNodes());
              const mcIdType *idsNodePtr=idsNode->begin();
              double center[3]; center[0]=pPtr[2*pointId]*vPtr[3*vecId]; center[1]=pPtr[2*pointId]*vPtr[3*vecId+1]; center[2]=pPtr[2*pointId]*vPtr[3*vecId+2];
              double vec[3]; vec[0]=vPtr[3*vecId+1]; vec[1]=-vPtr[3*vecId]; vec[2]=0.;
              double norm=vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
              if(std::abs(norm)>eps)
                {
                  double angle=INTERP_KERNEL::EdgeArcCircle::SafeAsin(norm);
                  mm3->rotate(center,vec,angle);
                }
              mm3->changeSpaceDimension(2);
              MCAuto<MEDCouplingUMesh> mm4=mm3->buildSpreadZonesWithPoly();
              const mcIdType *conn4=mm4->getNodalConnectivity()->begin();
              const mcIdType *connI4=mm4->getNodalConnectivityIndex()->begin();
              mcIdType nbOfCells=mm4->getNumberOfCells();
              for(mcIdType k=0;k<nbOfCells;k++)
                {
                  int l=0;
                  for(const mcIdType *work=conn4+connI4[k]+1;work!=conn4+connI4[k+1];work++,l++)
                    res->pushBackSilent(idsNodePtr[*work]);
                  res->pushBackSilent(-1);
                }
            }
        }
    }
  res->popBackSilent();
}

/*!
 * This method computes the normalized vector of the plane and the pos of the point belonging to the plane and the line defined by the vector going
 * through origin. The plane is defined by its nodal connectivity [ \b begin, \b end ).
 *
 * \param [in] eps below that value the dot product of 2 vectors is considered as colinears
 * \param [in] coords coordinates expected to have 3 components.
 * \param [in] begin start of the nodal connectivity of the face.
 * \param [in] end end of the nodal connectivity (excluded) of the face.
 * \param [out] v the normalized vector of size 3
 * \param [out] p the pos of plane
 */
void MEDCouplingUMesh::ComputeVecAndPtOfFace(double eps, const double *coords, const mcIdType *begin, const mcIdType *end, double *v, double *p)
{
  std::size_t nbPoints=std::distance(begin,end);
  if(nbPoints<3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ComputeVecAndPtOfFace : < of 3 points in face ! not able to find a plane on that face !");
  double vec[3]={0.,0.,0.};
  std::size_t j=0;
  bool refFound=false;
  for(;j<nbPoints-1 && !refFound;j++)
    {
      vec[0]=coords[3*begin[j+1]]-coords[3*begin[j]];
      vec[1]=coords[3*begin[j+1]+1]-coords[3*begin[j]+1];
      vec[2]=coords[3*begin[j+1]+2]-coords[3*begin[j]+2];
      double norm=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
      if(norm>eps)
        {
          refFound=true;
          vec[0]/=norm; vec[1]/=norm; vec[2]/=norm;
        }
    }
  for(std::size_t i=j;i<nbPoints-1;i++)
    {
      double curVec[3];
      curVec[0]=coords[3*begin[i+1]]-coords[3*begin[i]];
      curVec[1]=coords[3*begin[i+1]+1]-coords[3*begin[i]+1];
      curVec[2]=coords[3*begin[i+1]+2]-coords[3*begin[i]+2];
      double norm=sqrt(curVec[0]*curVec[0]+curVec[1]*curVec[1]+curVec[2]*curVec[2]);
      if(norm<eps)
        continue;
      curVec[0]/=norm; curVec[1]/=norm; curVec[2]/=norm;
      v[0]=vec[1]*curVec[2]-vec[2]*curVec[1]; v[1]=vec[2]*curVec[0]-vec[0]*curVec[2]; v[2]=vec[0]*curVec[1]-vec[1]*curVec[0];
      norm=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
      if(norm>eps)
        {
          v[0]/=norm; v[1]/=norm; v[2]/=norm;
          *p=v[0]*coords[3*begin[i]]+v[1]*coords[3*begin[i]+1]+v[2]*coords[3*begin[i]+2];
          return ;
        }
    }
  throw INTERP_KERNEL::Exception("Not able to find a normal vector of that 3D face !");
}

/*!
 * This method tries to obtain a well oriented polyhedron.
 * If the algorithm fails, an exception will be thrown.
 */
void MEDCouplingUMesh::TryToCorrectPolyhedronOrientation(mcIdType *begin, mcIdType *end, const double *coords)
{
  std::list< std::pair<mcIdType,mcIdType> > edgesOK,edgesFinished;
  std::size_t nbOfFaces=std::count(begin,end,-1)+1;
  std::vector<bool> isPerm(nbOfFaces,false);//field on faces False: I don't know, True : oriented
  isPerm[0]=true;
  mcIdType *bgFace=begin,*endFace=std::find(begin+1,end,-1);
  std::size_t nbOfEdgesInFace=std::distance(bgFace,endFace);
  for(std::size_t l=0;l<nbOfEdgesInFace;l++) { std::pair<mcIdType,mcIdType> p1(bgFace[l],bgFace[(l+1)%nbOfEdgesInFace]); edgesOK.push_back(p1); }
  //
  while(std::find(isPerm.begin(),isPerm.end(),false)!=isPerm.end())
    {
      bgFace=begin;
      std::size_t smthChanged=0;
      for(std::size_t i=0;i<nbOfFaces;i++)
        {
          endFace=std::find(bgFace+1,end,-1);
          nbOfEdgesInFace=std::distance(bgFace,endFace);
          if(!isPerm[i])
            {
              bool b=false;
              for(std::size_t j=0;j<nbOfEdgesInFace;j++)
                {
                  std::pair<mcIdType,mcIdType> p1(bgFace[j],bgFace[(j+1)%nbOfEdgesInFace]);
                  std::pair<mcIdType,mcIdType> p2(p1.second,p1.first);
                  bool b1=std::find(edgesOK.begin(),edgesOK.end(),p1)!=edgesOK.end();
                  bool b2=std::find(edgesOK.begin(),edgesOK.end(),p2)!=edgesOK.end();
                  if(b1 || b2) { b=b2; isPerm[i]=true; smthChanged++; break; }
                }
              if(isPerm[i])
                {
                  if(!b)
                    std::reverse(bgFace+1,endFace);
                  for(std::size_t j=0;j<nbOfEdgesInFace;j++)
                    {
                      std::pair<mcIdType,mcIdType> p1(bgFace[j],bgFace[(j+1)%nbOfEdgesInFace]);
                      std::pair<mcIdType,mcIdType> p2(p1.second,p1.first);
                      if(std::find(edgesOK.begin(),edgesOK.end(),p1)!=edgesOK.end())
                        { std::ostringstream oss; oss << "Face #" << j << " of polyhedron looks bad !"; throw INTERP_KERNEL::Exception(oss.str()); }
                      if(std::find(edgesFinished.begin(),edgesFinished.end(),p1)!=edgesFinished.end() || std::find(edgesFinished.begin(),edgesFinished.end(),p2)!=edgesFinished.end())
                        { std::ostringstream oss; oss << "Face #" << j << " of polyhedron looks bad !"; throw INTERP_KERNEL::Exception(oss.str()); }
                      std::list< std::pair<mcIdType,mcIdType> >::iterator it=std::find(edgesOK.begin(),edgesOK.end(),p2);
                      if(it!=edgesOK.end())
                        {
                          edgesOK.erase(it);
                          edgesFinished.push_back(p1);
                        }
                      else
                        edgesOK.push_back(p1);
                    }
                }
            }
          bgFace=endFace+1;
        }
      if(smthChanged==0)
        { throw INTERP_KERNEL::Exception("The polyhedron looks too bad to be repaired !"); }
    }
  if(!edgesOK.empty())
    { throw INTERP_KERNEL::Exception("The polyhedron looks too bad to be repaired : Some edges are shared only once !"); }
  if(INTERP_KERNEL::calculateVolumeForPolyh2<mcIdType,INTERP_KERNEL::ALL_C_MODE>(begin,ToIdType(std::distance(begin,end)),coords)<-EPS_FOR_POLYH_ORIENTATION)
    {//not lucky ! The first face was not correctly oriented : reorient all faces...
      bgFace=begin;
      for(std::size_t i=0;i<nbOfFaces;i++)
        {
          endFace=std::find(bgFace+1,end,-1);
          std::reverse(bgFace+1,endFace);
          bgFace=endFace+1;
        }
    }
}


/*!
 * This method makes the assumption spacedimension == meshdimension == 2.
 * This method works only for linear cells.
 *
 * \return a newly allocated array containing the connectivity of a polygon type enum included (NORM_POLYGON in pos#0)
 */
DataArrayIdType *MEDCouplingUMesh::buildUnionOf2DMesh() const
{
  if(getMeshDimension()!=2 || getSpaceDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildUnionOf2DMesh : meshdimension, spacedimension must be equal to 2 !");
  MCAuto<MEDCouplingUMesh> skin(computeSkin());
  mcIdType oldNbOfNodes(skin->getNumberOfNodes());
  MCAuto<DataArrayIdType> o2n(skin->zipCoordsTraducer());
  mcIdType nbOfNodesExpected(skin->getNumberOfNodes());
  MCAuto<DataArrayIdType> n2o(o2n->invertArrayO2N2N2O(oldNbOfNodes));
  mcIdType nbCells=skin->getNumberOfCells();
  if(nbCells==nbOfNodesExpected)
    return buildUnionOf2DMeshLinear(skin,n2o);
  else if(2*nbCells==nbOfNodesExpected)
    return buildUnionOf2DMeshQuadratic(skin,n2o);
  else
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildUnionOf2DMesh : the mesh 2D in input appears to be not in a single part of a 2D mesh !");
}

/*!
 * This method makes the assumption spacedimension == meshdimension == 3.
 * This method works only for linear cells.
 *
 * \return a newly allocated array containing the connectivity of a polygon type enum included (NORM_POLYHED in pos#0)
 */
DataArrayIdType *MEDCouplingUMesh::buildUnionOf3DMesh() const
{
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildUnionOf3DMesh : meshdimension, spacedimension must be equal to 2 !");
  MCAuto<MEDCouplingUMesh> m=computeSkin();
  const mcIdType *conn=m->getNodalConnectivity()->begin();
  const mcIdType *connI=m->getNodalConnectivityIndex()->begin();
  mcIdType nbOfCells=m->getNumberOfCells();
  MCAuto<DataArrayIdType> ret=DataArrayIdType::New(); ret->alloc(m->getNodalConnectivity()->getNumberOfTuples(),1);
  mcIdType *work=ret->getPointer();  *work++=INTERP_KERNEL::NORM_POLYHED;
  if(nbOfCells<1)
    return ret.retn();
  work=std::copy(conn+connI[0]+1,conn+connI[1],work);
  for(mcIdType i=1;i<nbOfCells;i++)
    {
      *work++=-1;
      work=std::copy(conn+connI[i]+1,conn+connI[i+1],work);
    }
  return ret.retn();
}

/*!
 * \brief Creates a graph of cell neighbors
 *  \return MEDCouplingSkyLineArray * - an sky line array the user should delete.
 *  In the sky line array, graph arcs are stored in terms of (index,value) notation.
 *  For example
 *  - index:  0 3 5 6 6
 *  - value:  1 2 3 2 3 3
 *  means 6 arcs (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
 *  Arcs are not doubled but reflexive (1,1) arcs are present for each cell
 */
MEDCouplingSkyLineArray* MEDCouplingUMesh::generateGraph() const
{
  checkConnectivityFullyDefined();

  int meshDim = this->getMeshDimension();
  MEDCoupling::DataArrayIdType* indexr=MEDCoupling::DataArrayIdType::New();
  MEDCoupling::DataArrayIdType* revConn=MEDCoupling::DataArrayIdType::New();
  this->getReverseNodalConnectivity(revConn,indexr);
  const mcIdType* indexr_ptr=indexr->begin();
  const mcIdType* revConn_ptr=revConn->begin();

  const MEDCoupling::DataArrayIdType* index;
  const MEDCoupling::DataArrayIdType* conn;
  conn=this->getNodalConnectivity(); // it includes a type as the 1st element!!!
  index=this->getNodalConnectivityIndex();
  mcIdType nbCells=this->getNumberOfCells();
  const mcIdType* index_ptr=index->begin();
  const mcIdType* conn_ptr=conn->begin();

  //creating graph arcs (cell to cell relations)
  //arcs are stored in terms of (index,value) notation
  // 0 3 5 6 6
  // 1 2 3 2 3 3
  // means 6 arcs (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
  // in present version arcs are not doubled but reflexive (1,1) arcs are present for each cell

  //warning here one node have less than or equal effective number of cell with it
  //but cell could have more than effective nodes
  //because other equals nodes in other domain (with other global inode)
  std::vector <mcIdType> cell2cell_index(nbCells+1,0);
  std::vector <mcIdType> cell2cell;
  cell2cell.reserve(3*nbCells);

  for (mcIdType icell=0; icell<nbCells;icell++)
    {
      std::map<mcIdType,mcIdType > counter;
      for (mcIdType iconn=index_ptr[icell]+1; iconn<index_ptr[icell+1];iconn++)
        {
          mcIdType inode=conn_ptr[iconn];
          for (mcIdType iconnr=indexr_ptr[inode]; iconnr<indexr_ptr[inode+1];iconnr++)
            {
              mcIdType icell2=revConn_ptr[iconnr];
              std::map<mcIdType,mcIdType>::iterator iter=counter.find(icell2);
              if (iter!=counter.end()) (iter->second)++;
              else counter.insert(std::make_pair(icell2,1));
            }
        }
      for (std::map<mcIdType,mcIdType>::const_iterator iter=counter.begin();
           iter!=counter.end(); iter++)
        if (iter->second >= meshDim)
          {
            cell2cell_index[icell+1]++;
            cell2cell.push_back(iter->first);
          }
    }
  indexr->decrRef();
  revConn->decrRef();
  cell2cell_index[0]=0;
  for (mcIdType icell=0; icell<nbCells;icell++)
    cell2cell_index[icell+1]=cell2cell_index[icell]+cell2cell_index[icell+1];

  //filling up index and value to create skylinearray structure
  MEDCouplingSkyLineArray * array(MEDCouplingSkyLineArray::New(cell2cell_index,cell2cell));
  return array;
}


void MEDCouplingUMesh::writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const
{
  mcIdType nbOfCells=getNumberOfCells();
  if(nbOfCells<=0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::writeVTK : the unstructured mesh has no cells !");
  ofs << "  <" << getVTKDataSetType() << ">\n";
  ofs << "    <Piece NumberOfPoints=\"" << getNumberOfNodes() << "\" NumberOfCells=\"" << nbOfCells << "\">\n";
  ofs << "      <PointData>\n" << pointData << std::endl;
  ofs << "      </PointData>\n";
  ofs << "      <CellData>\n" << cellData << std::endl;
  ofs << "      </CellData>\n";
  ofs << "      <Points>\n";
  if(getSpaceDimension()==3)
    _coords->writeVTK(ofs,8,"Points",byteData);
  else
    {
      MCAuto<DataArrayDouble> coo=_coords->changeNbOfComponents(3,0.);
      coo->writeVTK(ofs,8,"Points",byteData);
    }
  ofs << "      </Points>\n";
  ofs << "      <Cells>\n";
  const mcIdType *cPtr=_nodal_connec->begin();
  const mcIdType *cIPtr=_nodal_connec_index->begin();
  MCAuto<DataArrayIdType> faceoffsets=DataArrayIdType::New(); faceoffsets->alloc(nbOfCells,1);
  MCAuto<DataArrayIdType> types=DataArrayIdType::New(); types->alloc(nbOfCells,1);
  MCAuto<DataArrayIdType> offsets=DataArrayIdType::New(); offsets->alloc(nbOfCells,1);
  MCAuto<DataArrayIdType> connectivity=DataArrayIdType::New(); connectivity->alloc(_nodal_connec->getNumberOfTuples()-nbOfCells,1);
  mcIdType *w1=faceoffsets->getPointer(),*w2=types->getPointer(),*w3=offsets->getPointer(),*w4=connectivity->getPointer();
  mcIdType szFaceOffsets=0,szConn=0;
  for(mcIdType i=0;i<nbOfCells;i++,w1++,w2++,w3++)
    {
      *w2=cPtr[cIPtr[i]];
      if((INTERP_KERNEL::NormalizedCellType)cPtr[cIPtr[i]]!=INTERP_KERNEL::NORM_POLYHED)
        {
          *w1=-1;
          *w3=szConn+cIPtr[i+1]-cIPtr[i]-1; szConn+=cIPtr[i+1]-cIPtr[i]-1;
          w4=std::copy(cPtr+cIPtr[i]+1,cPtr+cIPtr[i+1],w4);
        }
      else
        {
          mcIdType deltaFaceOffset=cIPtr[i+1]-cIPtr[i]+1;
          *w1=szFaceOffsets+deltaFaceOffset; szFaceOffsets+=deltaFaceOffset;
          std::set<mcIdType> c(cPtr+cIPtr[i]+1,cPtr+cIPtr[i+1]); c.erase(-1);
          *w3=szConn+ToIdType(c.size()); szConn+=ToIdType(c.size());
          w4=std::copy(c.begin(),c.end(),w4);
        }
    }
  std::unique_ptr<mcIdType[]> medcoupling2vtkTypeTraducer_mcIdType(new mcIdType[MEDCOUPLING2VTKTYPETRADUCER_LGTH]);
  for(auto ii = 0; ii<MEDCOUPLING2VTKTYPETRADUCER_LGTH ; ++ii)
    medcoupling2vtkTypeTraducer_mcIdType[ii] = MEDCOUPLING2VTKTYPETRADUCER[ii]!=MEDCOUPLING2VTKTYPETRADUCER_NONE?MEDCOUPLING2VTKTYPETRADUCER[ii] : -1;
  types->transformWithIndArr(medcoupling2vtkTypeTraducer_mcIdType.get(),medcoupling2vtkTypeTraducer_mcIdType.get()+MEDCOUPLING2VTKTYPETRADUCER_LGTH);
  types->writeVTK(ofs,8,"UInt8","types",byteData);
  std::string vtkTypeName = Traits<mcIdType>::VTKReprStr;
  offsets->writeVTK(ofs,8,vtkTypeName,"offsets",byteData);
  if(szFaceOffsets!=0)
    {//presence of Polyhedra
      connectivity->reAlloc(szConn);
      faceoffsets->writeVTK(ofs,8,vtkTypeName,"faceoffsets",byteData);
      MCAuto<DataArrayIdType> faces=DataArrayIdType::New(); faces->alloc(szFaceOffsets,1);
      w1=faces->getPointer();
      for(mcIdType i=0;i<nbOfCells;i++)
        if((INTERP_KERNEL::NormalizedCellType)cPtr[cIPtr[i]]==INTERP_KERNEL::NORM_POLYHED)
          {
            mcIdType nbFaces=ToIdType(std::count(cPtr+cIPtr[i]+1,cPtr+cIPtr[i+1],-1))+1;
            *w1++=nbFaces;
            const mcIdType *w6=cPtr+cIPtr[i]+1,*w5=0;
            for(mcIdType j=0;j<nbFaces;j++)
              {
                w5=std::find(w6,cPtr+cIPtr[i+1],-1);
                *w1++=ToIdType(std::distance(w6,w5));
                w1=std::copy(w6,w5,w1);
                w6=w5+1;
              }
          }
      faces->writeVTK(ofs,8,vtkTypeName,"faces",byteData);
    }
  connectivity->writeVTK(ofs,8,vtkTypeName,"connectivity",byteData);
  ofs << "      </Cells>\n";
  ofs << "    </Piece>\n";
  ofs << "  </" << getVTKDataSetType() << ">\n";
}

void MEDCouplingUMesh::reprQuickOverview(std::ostream& stream) const
{
  stream << "MEDCouplingUMesh C++ instance at " << this << ". Name : \"" << getName() << "\".";
  if(_mesh_dim==-2)
    { stream << " Not set !"; return ; }
  stream << " Mesh dimension : " << _mesh_dim << ".";
  if(_mesh_dim==-1)
    return ;
  if(!_coords)
    { stream << " No coordinates set !"; return ; }
  if(!_coords->isAllocated())
    { stream << " Coordinates set but not allocated !"; return ; }
  stream << " Space dimension : " << _coords->getNumberOfComponents() << "." << std::endl;
  stream << "Number of nodes : " << _coords->getNumberOfTuples() << ".";
  if(!_nodal_connec_index)
    { stream << std::endl << "Nodal connectivity NOT set !"; return ; }
  if(!_nodal_connec_index->isAllocated())
    { stream << std::endl << "Nodal connectivity set but not allocated !"; return ; }
  mcIdType lgth=_nodal_connec_index->getNumberOfTuples();
  std::size_t cpt=_nodal_connec_index->getNumberOfComponents();
  if(cpt!=1 || lgth<1)
    return ;
  stream << std::endl << "Number of cells : " << lgth-1 << ".";
}

std::string MEDCouplingUMesh::getVTKDataSetType() const
{
  return std::string("UnstructuredGrid");
}

std::string MEDCouplingUMesh::getVTKFileExtension() const
{
  return std::string("vtu");
}



/**
 * Provides a renumbering of the cells of this (which has to be a piecewise connected 1D line), so that
 * the segments of the line are indexed in consecutive order (i.e. cells \a i and \a i+1 are neighbors).
 * This doesn't modify the mesh. This method only works using nodal connectivity consideration. Coordinates of nodes are ignored here.
 * The caller is to deal with the resulting DataArrayIdType.
 *  \throw If the coordinate array is not set.
 *  \throw If the nodal connectivity of the cells is not defined.
 *  \throw If m1 is not a mesh of dimension 2, or m1 is not a mesh of dimension 1
 *  \throw If m2 is not a (piecewise) line (i.e. if a point has more than 2 adjacent segments)
 *
 * \sa DataArrayIdType::sortEachPairToMakeALinkedList
 */
DataArrayIdType *MEDCouplingUMesh::orderConsecutiveCells1D() const
{
  checkFullyDefined();
  if(getMeshDimension()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::orderConsecutiveCells1D works on unstructured mesh with meshdim = 1 !");

  // Check that this is a line (and not a more complex 1D mesh) - each point is used at most by 2 segments:
  MCAuto<DataArrayIdType> _d(DataArrayIdType::New()),_dI(DataArrayIdType::New());
  MCAuto<DataArrayIdType> _rD(DataArrayIdType::New()),_rDI(DataArrayIdType::New());
  MCAuto<MEDCouplingUMesh> m_points(buildDescendingConnectivity(_d, _dI, _rD, _rDI));
  const mcIdType *d(_d->begin()), *dI(_dI->begin());
  const mcIdType *rD(_rD->begin()), *rDI(_rDI->begin());
  MCAuto<DataArrayIdType> _dsi(_rDI->deltaShiftIndex());
  const mcIdType * dsi(_dsi->begin());
  MCAuto<DataArrayIdType> dsii = _dsi->findIdsNotInRange(0,3);
  m_points=0;
  if (dsii->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::orderConsecutiveCells1D only work with a mesh being a (piecewise) connected line!");

  mcIdType nc=getNumberOfCells();
  MCAuto<DataArrayIdType> result(DataArrayIdType::New());
  result->alloc(nc,1);

  // set of edges not used so far
  std::set<mcIdType> edgeSet;
  for (mcIdType i=0; i<nc; edgeSet.insert(i), i++);

  mcIdType startSeg=0;
  mcIdType newIdx=0;
  // while we have points with only one neighbor segments
  do
    {
      std::list<mcIdType> linePiece;
      // fills a list of consecutive segment linked to startSeg. This can go forward or backward.
      for (int direction=0;direction<2;direction++) // direction=0 --> forward, direction=1 --> backward
        {
          // Fill the list forward (resp. backward) from the start segment:
          mcIdType activeSeg = startSeg;
          mcIdType prevPointId = -20;
          mcIdType ptId;
          while (!edgeSet.empty())
            {
              if (!(direction == 1 && prevPointId==-20)) // prevent adding twice startSeg
                {
                  if (direction==0)
                    linePiece.push_back(activeSeg);
                  else
                    linePiece.push_front(activeSeg);
                  edgeSet.erase(activeSeg);
                }

              mcIdType ptId1 = d[dI[activeSeg]], ptId2 = d[dI[activeSeg]+1];
              ptId = direction ? (ptId1 == prevPointId ? ptId2 : ptId1) : (ptId2 == prevPointId ? ptId1 : ptId2);
              if (dsi[ptId] == 1) // hitting the end of the line
                break;

              prevPointId = ptId;
              mcIdType seg1 = rD[rDI[ptId]], seg2 = rD[rDI[ptId]+1];
              activeSeg = (seg1 == activeSeg) ? seg2 : seg1;

              //for piecewise meshes made up of closed parts
              bool segmentAlreadyTreated = (std::find(linePiece.begin(), linePiece.end(), activeSeg) != linePiece.end());
              if(segmentAlreadyTreated)
                break;
            }
        }
      // Done, save final piece into DA:
      std::copy(linePiece.begin(), linePiece.end(), result->getPointer()+newIdx);
      newIdx += ToIdType(linePiece.size());

      // identify next valid start segment (one which is not consumed)
      if(!edgeSet.empty())
        startSeg = *(edgeSet.begin());

    }
  while (!edgeSet.empty());
  return result.retn();
}

/**
 * This method split some of edges of 2D cells in \a this. The edges to be split are specified in \a subNodesInSeg
 * and in \a subNodesInSegI using \ref numbering-indirect storage mode.
 * To do the work this method can optionally needs information about middle of subedges for quadratic cases if
 * a minimal creation of new nodes is wanted.
 * So this method try to reduce at most the number of new nodes. The only case that can lead this method to add
 * nodes if a SEG3 is split without information of middle.
 * \b WARNING : is returned value is different from 0 a call to MEDCouplingUMesh::mergeNodes is necessary to
 * avoid to have a non conform mesh.
 *
 * \return mcIdType - the number of new nodes created (in most of cases 0).
 *
 * \throw If \a this is not coherent.
 * \throw If \a this has not spaceDim equal to 2.
 * \throw If \a this has not meshDim equal to 2.
 * \throw If some subcells needed to be split are orphan.
 * \sa MEDCouplingUMesh::conformize2D
 */
mcIdType MEDCouplingUMesh::split2DCells(const DataArrayIdType *desc, const DataArrayIdType *descI, const DataArrayIdType *subNodesInSeg, const DataArrayIdType *subNodesInSegI, const DataArrayIdType *midOpt, const DataArrayIdType *midOptI)
{
  if(!desc || !descI || !subNodesInSeg || !subNodesInSegI)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::split2DCells : the 4 first arrays must be not null !");
  desc->checkAllocated(); descI->checkAllocated(); subNodesInSeg->checkAllocated(); subNodesInSegI->checkAllocated();
  if(getSpaceDimension()!=2 || getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::split2DCells : This method only works for meshes with spaceDim=2 and meshDim=2 !");
  if(midOpt==0 && midOptI==0)
    {
      split2DCellsLinear(desc,descI,subNodesInSeg,subNodesInSegI);
      return 0;
    }
  else if(midOpt!=0 && midOptI!=0)
    return split2DCellsQuadratic(desc,descI,subNodesInSeg,subNodesInSegI,midOpt,midOptI);
  else
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::split2DCells : middle parameters must be set to null for all or not null for all.");
}

/*!
 * This method compute the convex hull of a single 2D cell. This method tries to conserve at maximum the given input connectivity. In particular, if the orientation of cell is not clockwise
 * as in MED format norm. If definitely the result of Jarvis algorithm is not matchable with the input connectivity, the result will be copied into \b nodalConnecOut parameter and
 * the geometric cell type set to INTERP_KERNEL::NORM_POLYGON.
 * This method excepts that \b coords parameter is expected to be in dimension 2. [ \b nodalConnBg , \b nodalConnEnd ) is the nodal connectivity of the input
 * cell (geometric cell type included at the position 0). If the meshdimension of the input cell is not equal to 2 an INTERP_KERNEL::Exception will be thrown.
 *
 * \return false if the input connectivity represents already the convex hull, true if the input cell needs to be reordered.
 */
bool MEDCouplingUMesh::BuildConvexEnvelopOf2DCellJarvis(const double *coords, const mcIdType *nodalConnBg, const mcIdType *nodalConnEnd, DataArrayIdType *nodalConnecOut)
{
  std::size_t sz=std::distance(nodalConnBg,nodalConnEnd);
  if(sz>=4)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)*nodalConnBg);
      if(cm.getDimension()==2)
        {
          const mcIdType *node=nodalConnBg+1;
          mcIdType startNode=*node++;
          double refX=coords[2*startNode];
          for(;node!=nodalConnEnd;node++)
            {
              if(coords[2*(*node)]<refX)
                {
                  startNode=*node;
                  refX=coords[2*startNode];
                }
            }
          std::vector<mcIdType> tmpOut; tmpOut.reserve(sz); tmpOut.push_back(startNode);
          refX=1e300;
          double tmp1;
          double tmp2[2];
          double angle0=-M_PI/2;
          //
          mcIdType nextNode=-1;
          mcIdType prevNode=-1;
          double resRef;
          double angleNext=0.;
          while(nextNode!=startNode)
            {
              nextNode=-1;
              resRef=1e300;
              for(node=nodalConnBg+1;node!=nodalConnEnd;node++)
                {
                  if(*node!=tmpOut.back() && *node!=prevNode)
                    {
                      tmp2[0]=coords[2*(*node)]-coords[2*tmpOut.back()]; tmp2[1]=coords[2*(*node)+1]-coords[2*tmpOut.back()+1];
                      double angleM=INTERP_KERNEL::EdgeArcCircle::GetAbsoluteAngle(tmp2,tmp1);
                      double res;
                      if(angleM<=angle0)
                        res=angle0-angleM;
                      else
                        res=angle0-angleM+2.*M_PI;
                      if(res<resRef)
                        {
                          nextNode=*node;
                          resRef=res;
                          angleNext=angleM;
                        }
                    }
                }
              if(nextNode!=startNode)
                {
                  angle0=angleNext-M_PI;
                  if(angle0<-M_PI)
                    angle0+=2*M_PI;
                  prevNode=tmpOut.back();
                  tmpOut.push_back(nextNode);
                }
            }
          std::vector<mcIdType> tmp3(2*(sz-1));
          std::vector<mcIdType>::iterator it=std::copy(nodalConnBg+1,nodalConnEnd,tmp3.begin());
          std::copy(nodalConnBg+1,nodalConnEnd,it);
          if(std::search(tmp3.begin(),tmp3.end(),tmpOut.begin(),tmpOut.end())!=tmp3.end())
            {
              nodalConnecOut->insertAtTheEnd(nodalConnBg,nodalConnEnd);
              return false;
            }
          if(std::search(tmp3.rbegin(),tmp3.rend(),tmpOut.begin(),tmpOut.end())!=tmp3.rend())
            {
              nodalConnecOut->insertAtTheEnd(nodalConnBg,nodalConnEnd);
              return false;
            }
          else
            {
              nodalConnecOut->pushBackSilent(ToIdType(INTERP_KERNEL::NORM_POLYGON));
              nodalConnecOut->insertAtTheEnd(tmpOut.begin(),tmpOut.end());
              return true;
            }
        }
      else
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::BuildConvexEnvelopOf2DCellJarvis : invalid 2D cell connectivity !");
    }
  else
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::BuildConvexEnvelopOf2DCellJarvis : invalid 2D cell connectivity !");
}

/*!
 * This method works on a pair input (\b arrIn, \b arrIndxIn) where \b arr indexes is in \b arrIndxIn.
 * This method expects that these two input arrays come from the output of MEDCouplingUMesh::computeNeighborsOfCells method.
 * This method start from id 0 that will be contained in output DataArrayIdType. It searches then all neighbors of id0 looking at arrIn[arrIndxIn[0]:arrIndxIn[0+1]].
 * Then it is repeated recursively until either all ids are fetched or no more ids are reachable step by step.
 * A negative value in \b arrIn means that it is ignored.
 * This method is useful to see if a mesh is contiguous regarding its connectivity. If it is not the case the size of returned array is different from arrIndxIn->getNumberOfTuples()-1.
 *
 * \param [in] arrIn arr origin array from which the extraction will be done.
 * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
 * \return a newly allocated DataArray that stores all ids fetched by the gradually spread process.
 * \sa MEDCouplingUMesh::ComputeSpreadZoneGraduallyFromSeed, MEDCouplingUMesh::partitionBySpreadZone
 */
DataArrayIdType *MEDCouplingUMesh::ComputeSpreadZoneGradually(const DataArrayIdType *arrIn, const DataArrayIdType *arrIndxIn)
{
  mcIdType seed=0,nbOfDepthPeelingPerformed=0;
  return ComputeSpreadZoneGraduallyFromSeed(&seed,&seed+1,arrIn,arrIndxIn,-1,nbOfDepthPeelingPerformed);
}

/*!
 * This method works on a pair input (\b arrIn, \b arrIndxIn) where \b arr indexes is in \b arrIndxIn.
 * This method expects that these two input arrays come from the output of MEDCouplingUMesh::computeNeighborsOfCells method.
 * This method start from id 0 that will be contained in output DataArrayIdType. It searches then all neighbors of id0 regarding arrIn[arrIndxIn[0]:arrIndxIn[0+1]].
 * Then it is repeated recursively until either all ids are fetched or no more ids are reachable step by step.
 * A negative value in \b arrIn means that it is ignored.
 * This method is useful to see if a mesh is contiguous regarding its connectivity. If it is not the case the size of returned array is different from arrIndxIn->getNumberOfTuples()-1.
 * \param [in] seedBg the begin pointer (included) of an array containing the seed of the search zone
 * \param [in] seedEnd the end pointer (not included) of an array containing the seed of the search zone
 * \param [in] arrIn arr origin array from which the extraction will be done.
 * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
 * \param [in] nbOfDepthPeeling the max number of peels requested in search. By default -1, that is to say, no limit.
 * \param [out] nbOfDepthPeelingPerformed the number of peels effectively performed. May be different from \a nbOfDepthPeeling
 * \return a newly allocated DataArray that stores all ids fetched by the gradually spread process.
 * \sa MEDCouplingUMesh::partitionBySpreadZone
 */
DataArrayIdType *MEDCouplingUMesh::ComputeSpreadZoneGraduallyFromSeed(const mcIdType *seedBg, const mcIdType *seedEnd, const DataArrayIdType *arrIn, const DataArrayIdType *arrIndxIn, mcIdType nbOfDepthPeeling, mcIdType& nbOfDepthPeelingPerformed)
{
  nbOfDepthPeelingPerformed=0;
  if(!arrIndxIn)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ComputeSpreadZoneGraduallyFromSeed : arrIndxIn input pointer is NULL !");
  mcIdType nbOfTuples=arrIndxIn->getNumberOfTuples()-1;
  if(nbOfTuples<=0)
    {
      DataArrayIdType *ret=DataArrayIdType::New(); ret->alloc(0,1);
      return ret;
    }
  //
  std::vector<bool> fetched(nbOfTuples,false);
  return ComputeSpreadZoneGraduallyFromSeedAlg(fetched,seedBg,seedEnd,arrIn,arrIndxIn,nbOfDepthPeeling,nbOfDepthPeelingPerformed);
}



/*!
 * \b this is expected to be a mesh fully defined whose spaceDim==meshDim.
 * It returns a new allocated mesh having the same mesh dimension and lying on same coordinates.
 * The returned mesh contains as poly cells as number of contiguous zone (regarding connectivity).
 * A spread contiguous zone is built using poly cells (polyhedra in 3D, polygons in 2D and polyline in 1D).
 * The sum of measure field of returned mesh is equal to the sum of measure field of this.
 *
 * \return a newly allocated mesh lying on the same coords than \b this with same meshdimension than \b this.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildSpreadZonesWithPoly() const
{
  checkFullyDefined();
  int mdim=getMeshDimension();
  int spaceDim=getSpaceDimension();
  if(mdim!=spaceDim)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSpreadZonesWithPoly : meshdimension and spacedimension do not match !");
  std::vector<DataArrayIdType *> partition=partitionBySpreadZone();
  std::vector< MCAuto<DataArrayIdType> > partitionAuto; partitionAuto.reserve(partition.size());
  std::copy(partition.begin(),partition.end(),std::back_insert_iterator<std::vector< MCAuto<DataArrayIdType> > >(partitionAuto));
  MCAuto<MEDCouplingUMesh> ret=MEDCouplingUMesh::New(getName(),mdim);
  ret->setCoords(getCoords());
  ret->allocateCells(ToIdType(partition.size()));
  //
  for(std::vector<DataArrayIdType *>::const_iterator it=partition.begin();it!=partition.end();it++)
    {
      MCAuto<MEDCouplingUMesh> tmp=static_cast<MEDCouplingUMesh *>(buildPartOfMySelf((*it)->begin(),(*it)->end(),true));
      MCAuto<DataArrayIdType> cell;
      switch(mdim)
      {
        case 2:
          cell=tmp->buildUnionOf2DMesh();
          break;
        case 3:
          cell=tmp->buildUnionOf3DMesh();
          break;
        default:
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSpreadZonesWithPoly : meshdimension supported are [2,3] ! Not implemented yet for others !");
      }

      ret->insertNextCell((INTERP_KERNEL::NormalizedCellType)cell->getIJSafe(0,0),cell->getNumberOfTuples()-1,cell->begin()+1);
    }
  //
  ret->finishInsertingCells();
  return ret.retn();
}

/*!
 * This method partitions \b this into contiguous zone.
 * This method only needs a well defined connectivity. Coordinates are not considered here.
 * This method returns a vector of \b newly allocated arrays that the caller has to deal with.
 */
std::vector<DataArrayIdType *> MEDCouplingUMesh::partitionBySpreadZone() const
{
  DataArrayIdType *neigh=0,*neighI=0;
  computeNeighborsOfCells(neigh,neighI);
  MCAuto<DataArrayIdType> neighAuto(neigh),neighIAuto(neighI);
  return PartitionBySpreadZone(neighAuto,neighIAuto);
}

std::vector<DataArrayIdType *> MEDCouplingUMesh::PartitionBySpreadZone(const DataArrayIdType *arrIn, const DataArrayIdType *arrIndxIn)
{
  if(!arrIn || !arrIndxIn)
    throw INTERP_KERNEL::Exception("PartitionBySpreadZone : null input pointers !");
  arrIn->checkAllocated(); arrIndxIn->checkAllocated();
  mcIdType nbOfTuples(arrIndxIn->getNumberOfTuples());
  if(arrIn->getNumberOfComponents()!=1 || arrIndxIn->getNumberOfComponents()!=1 || nbOfTuples<1)
    throw INTERP_KERNEL::Exception("PartitionBySpreadZone : invalid arrays in input !");
  mcIdType nbOfCellsCur(nbOfTuples-1);
  std::vector<DataArrayIdType *> ret;
  if(nbOfCellsCur<=0)
    return ret;
  std::vector<bool> fetchedCells(nbOfCellsCur,false);
  std::vector< MCAuto<DataArrayIdType> > ret2;
  mcIdType seed=0;
  while(seed<nbOfCellsCur)
    {
      mcIdType nbOfPeelPerformed=0;
      ret2.push_back(ComputeSpreadZoneGraduallyFromSeedAlg(fetchedCells,&seed,&seed+1,arrIn,arrIndxIn,-1,nbOfPeelPerformed));
      seed=ToIdType(std::distance(fetchedCells.begin(),std::find(fetchedCells.begin()+seed,fetchedCells.end(),false)));
    }
  for(std::vector< MCAuto<DataArrayIdType> >::iterator it=ret2.begin();it!=ret2.end();it++)
    ret.push_back((*it).retn());
  return ret;
}

/*!
 * This method returns given a distribution of cell type (returned for example by MEDCouplingUMesh::getDistributionOfTypes method and customized after) a
 * newly allocated DataArrayIdType instance with 2 components ready to be interpreted as input of DataArrayIdType::findRangeIdForEachTuple method.
 *
 * \param [in] code a code with the same format than those returned by MEDCouplingUMesh::getDistributionOfTypes except for the code[3*k+2] that should contain start id of chunck.
 * \return a newly allocated DataArrayIdType to be managed by the caller.
 * \throw In case of \a code has not the right format (typically of size 3*n)
 */
DataArrayIdType *MEDCouplingUMesh::ComputeRangesFromTypeDistribution(const std::vector<mcIdType>& code)
{
  MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
  std::size_t nb=code.size()/3;
  if(code.size()%3!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ComputeRangesFromTypeDistribution : invalid input code !");
  ret->alloc(nb,2);
  mcIdType *retPtr=ret->getPointer();
  for(std::size_t i=0;i<nb;i++,retPtr+=2)
    {
      retPtr[0]=code[3*i+2];
      retPtr[1]=code[3*i+2]+code[3*i+1];
    }
  return ret.retn();
}

/*!
 * This method expects that \a this a 3D mesh (spaceDim=3 and meshDim=3) with all coordinates and connectivities set.
 * All cells in \a this are expected to be linear 3D cells.
 * This method will split **all** 3D cells in \a this into INTERP_KERNEL::NORM_TETRA4 cells and put them in the returned mesh.
 * It leads to an increase to number of cells.
 * This method contrary to MEDCouplingUMesh::simplexize can append coordinates in \a this to perform its work.
 * The \a nbOfAdditionalPoints returned value informs about it. If > 0, the coordinates array in returned mesh will have \a nbOfAdditionalPoints
 * more tuples (nodes) than in \a this. Anyway, all the nodes in \a this (with the same order) will be in the returned mesh.
 *
 * \param [in] policy - the policy of splitting that must be in (PLANAR_FACE_5, PLANAR_FACE_6, GENERAL_24, GENERAL_48). The policy will be used only for INTERP_KERNEL::NORM_HEXA8 cells.
 *                      For all other cells, the splitting policy will be ignored. See INTERP_KERNEL::SplittingPolicy for the images.
 * \param [out] nbOfAdditionalPoints - number of nodes added to \c this->_coords. If > 0 a new coordinates object will be constructed result of the aggregation of the old one and the new points added.
 * \param [out] n2oCells - A new instance of DataArrayIdType holding, for each new cell,
 *          an id of old cell producing it. The caller is to delete this array using
 *         decrRef() as it is no more needed.
 * \return MEDCoupling1SGTUMesh * - the mesh containing only INTERP_KERNEL::NORM_TETRA4 cells.
 *
 * \warning This method operates on each cells in this independently ! So it can leads to non conform mesh in returned value ! If you expect to have a conform mesh in output
 * the policy PLANAR_FACE_6 should be used on a mesh sorted with MEDCoupling1SGTUMesh::sortHexa8EachOther.
 *
 * \throw If \a this is not a 3D mesh (spaceDim==3 and meshDim==3).
 * \throw If \a this is not fully constituted with linear 3D cells.
 * \sa MEDCouplingUMesh::simplexize, MEDCoupling1SGTUMesh::sortHexa8EachOther
 */
MEDCoupling1SGTUMesh *MEDCouplingUMesh::tetrahedrize(int policy, DataArrayIdType *& n2oCells, mcIdType& nbOfAdditionalPoints) const
{
  INTERP_KERNEL::SplittingPolicy pol((INTERP_KERNEL::SplittingPolicy)policy);
  checkConnectivityFullyDefined();
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::tetrahedrize : only available for mesh with meshdim == 3 and spacedim == 3 !");
  mcIdType nbOfCells=getNumberOfCells();
  mcIdType nbNodes(getNumberOfNodes());
  MCAuto<MEDCoupling1SGTUMesh> ret0(MEDCoupling1SGTUMesh::New(getName(),INTERP_KERNEL::NORM_TETRA4));
  MCAuto<DataArrayIdType> ret(DataArrayIdType::New()); ret->alloc(nbOfCells,1);
  mcIdType *retPt(ret->getPointer());
  MCAuto<DataArrayIdType> newConn(DataArrayIdType::New()); newConn->alloc(0,1);
  MCAuto<DataArrayDouble> addPts(DataArrayDouble::New()); addPts->alloc(0,1);
  const mcIdType *oldc(_nodal_connec->begin());
  const mcIdType *oldci(_nodal_connec_index->begin());
  const double *coords(_coords->begin());
  for(mcIdType i=0;i<nbOfCells;i++,oldci++,retPt++)
    {
      std::vector<mcIdType> a; std::vector<double> b;
      INTERP_KERNEL::SplitIntoTetras(pol,(INTERP_KERNEL::NormalizedCellType)oldc[oldci[0]],oldc+oldci[0]+1,oldc+oldci[1],coords,a,b);
      std::size_t nbOfTet(a.size()/4); *retPt=ToIdType(nbOfTet);
      const mcIdType *aa(&a[0]);
      if(!b.empty())
        {
          for(std::vector<mcIdType>::iterator it=a.begin();it!=a.end();it++)
            if(*it<0)
              *it=(-(*(it))-1+nbNodes);
          addPts->insertAtTheEnd(b.begin(),b.end());
          nbNodes+=ToIdType(b.size()/3);
        }
      for(std::size_t j=0;j<nbOfTet;j++,aa+=4)
        newConn->insertAtTheEnd(aa,aa+4);
    }
  if(!addPts->empty())
    {
      addPts->rearrange(3);
      nbOfAdditionalPoints=addPts->getNumberOfTuples();
      addPts=DataArrayDouble::Aggregate(getCoords(),addPts);
      ret0->setCoords(addPts);
    }
  else
    {
      nbOfAdditionalPoints=0;
      ret0->setCoords(getCoords());
    }
  ret0->setNodalConnectivity(newConn);
  //
  ret->computeOffsetsFull();
  n2oCells=ret->buildExplicitArrOfSliceOnScaledArr(0,nbOfCells,1);
  return ret0.retn();
}

MEDCouplingUMeshCellIterator::MEDCouplingUMeshCellIterator(MEDCouplingUMesh *mesh):_mesh(mesh),_cell(new MEDCouplingUMeshCell(mesh)),
    _own_cell(true),_cell_id(-1),_nb_cell(0)
{
  if(mesh)
    {
      mesh->incrRef();
      _nb_cell=mesh->getNumberOfCells();
    }
}

MEDCouplingUMeshCellIterator::~MEDCouplingUMeshCellIterator()
{
  if(_mesh)
    _mesh->decrRef();
  if(_own_cell)
    delete _cell;
}

MEDCouplingUMeshCellIterator::MEDCouplingUMeshCellIterator(MEDCouplingUMesh *mesh, MEDCouplingUMeshCell *itc, mcIdType bg, mcIdType end):_mesh(mesh),_cell(itc),
    _own_cell(false),_cell_id(bg-1),
    _nb_cell(end)
{
  if(mesh)
    mesh->incrRef();
}

MEDCouplingUMeshCell *MEDCouplingUMeshCellIterator::nextt()
{
  _cell_id++;
  if(_cell_id<_nb_cell)
    {
      _cell->next();
      return _cell;
    }
  else
    return 0;
}

MEDCouplingUMeshCellByTypeEntry::MEDCouplingUMeshCellByTypeEntry(MEDCouplingUMesh *mesh):_mesh(mesh)
{
  if(_mesh)
    _mesh->incrRef();
}

MEDCouplingUMeshCellByTypeIterator *MEDCouplingUMeshCellByTypeEntry::iterator()
{
  return new MEDCouplingUMeshCellByTypeIterator(_mesh);
}

MEDCouplingUMeshCellByTypeEntry::~MEDCouplingUMeshCellByTypeEntry()
{
  if(_mesh)
    _mesh->decrRef();
}

MEDCouplingUMeshCellEntry::MEDCouplingUMeshCellEntry(MEDCouplingUMesh *mesh,  INTERP_KERNEL::NormalizedCellType type, MEDCouplingUMeshCell *itc, mcIdType bg, mcIdType end):_mesh(mesh),_type(type),
    _itc(itc),
    _bg(bg),_end(end)
{
  if(_mesh)
    _mesh->incrRef();
}

MEDCouplingUMeshCellEntry::~MEDCouplingUMeshCellEntry()
{
  if(_mesh)
    _mesh->decrRef();
}

INTERP_KERNEL::NormalizedCellType MEDCouplingUMeshCellEntry::getType() const
{
  return _type;
}

mcIdType MEDCouplingUMeshCellEntry::getNumberOfElems() const
{
  return _end-_bg;
}

MEDCouplingUMeshCellIterator *MEDCouplingUMeshCellEntry::iterator()
{
  return new MEDCouplingUMeshCellIterator(_mesh,_itc,_bg,_end);
}

MEDCouplingUMeshCellByTypeIterator::MEDCouplingUMeshCellByTypeIterator(MEDCouplingUMesh *mesh):_mesh(mesh),_cell(new MEDCouplingUMeshCell(mesh)),_cell_id(0),_nb_cell(0)
{
  if(mesh)
    {
      mesh->incrRef();
      _nb_cell=mesh->getNumberOfCells();
    }
}

MEDCouplingUMeshCellByTypeIterator::~MEDCouplingUMeshCellByTypeIterator()
{
  if(_mesh)
    _mesh->decrRef();
  delete _cell;
}

MEDCouplingUMeshCellEntry *MEDCouplingUMeshCellByTypeIterator::nextt()
{
  const mcIdType *c=_mesh->getNodalConnectivity()->begin();
  const mcIdType *ci=_mesh->getNodalConnectivityIndex()->begin();
  if(_cell_id<_nb_cell)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)c[ci[_cell_id]];
      mcIdType nbOfElems=ToIdType(std::distance(ci+_cell_id,std::find_if(ci+_cell_id,ci+_nb_cell,MEDCouplingImpl::ConnReader(c,type))));
      mcIdType startId=_cell_id;
      _cell_id+=nbOfElems;
      return new MEDCouplingUMeshCellEntry(_mesh,type,_cell,startId,_cell_id);
    }
  else
    return 0;
}

MEDCouplingUMeshCell::MEDCouplingUMeshCell(MEDCouplingUMesh *mesh):_conn(0),_conn_indx(0),_conn_lgth(NOTICABLE_FIRST_VAL)
{
  if(mesh)
    {
      _conn=mesh->getNodalConnectivity()->getPointer();
      _conn_indx=mesh->getNodalConnectivityIndex()->getPointer();
    }
}

void MEDCouplingUMeshCell::next()
{
  if(_conn_lgth!=NOTICABLE_FIRST_VAL)
    {
      _conn+=_conn_lgth;
      _conn_indx++;
    }
  _conn_lgth=_conn_indx[1]-_conn_indx[0];
}

std::string MEDCouplingUMeshCell::repr() const
{
  if(_conn_lgth!=NOTICABLE_FIRST_VAL)
    {
      std::ostringstream oss; oss << "Cell Type " << INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)_conn[0]).getRepr();
      oss << " : ";
      std::copy(_conn+1,_conn+_conn_lgth,std::ostream_iterator<mcIdType>(oss," "));
      return oss.str();
    }
  else
    return std::string("MEDCouplingUMeshCell::repr : Invalid pos");
}

INTERP_KERNEL::NormalizedCellType MEDCouplingUMeshCell::getType() const
{
  if(_conn_lgth!=NOTICABLE_FIRST_VAL)
    return (INTERP_KERNEL::NormalizedCellType)_conn[0];
  else
    return INTERP_KERNEL::NORM_ERROR;
}

const mcIdType *MEDCouplingUMeshCell::getAllConn(mcIdType& lgth) const
{
  lgth=_conn_lgth;
  if(_conn_lgth!=NOTICABLE_FIRST_VAL)
    return _conn;
  else
    return 0;
}

/// @cond INTERNAL

namespace MEDCouplingImpl
{
  const mcIdType theUndefID = std::numeric_limits< mcIdType >::max(); //!< undefined cell id

  //================================================================================
  /*!
   * \brief Encode a cell id and a mesh index into a code
   *  \param [in] id - cell id
   *  \param [in] iMesh - mesh index [0,1]
   *  \return mcIdType - code
   */
  //================================================================================

  mcIdType encodeID( mcIdType id, int iMesh )
  {
    return ( id + 1 ) * ( iMesh ? -1 : 1 );
  }
  //================================================================================
  /*!
   * \brief Return cell id and mesh index by a given id
   *  \param [in] id - code of a cell in a mesh
   *  \param [out] iMesh - returned mesh index
   *  \return mcIdType - cell id
   */
  //================================================================================

  mcIdType decodeID( mcIdType id, int& iMesh )
  {
    iMesh = ( id < 0 );
    return std::abs( id ) - 1;
  }

  //================================================================================
  /*!
   * \brief return another face sharing two given nodes of a face edge
   *  \param [in] n0 - 1st node of the edge
   *  \param [in] n1 - 2nd node of the edge
   *  \param [in] inputFaceID - face including \a n0 andf \a n2
   *  \param [in] mesh - object and reference meshes
   *  \param [in] revNodal - reverse nodal connectivity of the two meshes
   *  \param [in] revNodalIndx - index of reverse nodal connectivity of the two meshes
   *  \param [out] facesByEdge - return another face including \a n0 andf \a n2
   *  \param [out] equalFaces - return faces equal to facesByEdge
   */
  //================================================================================

  void getFacesOfEdge( mcIdType n0,
                       mcIdType n1,
                       mcIdType inputFaceID,
                       MEDCouplingUMesh* mesh[],
                       MCAuto<DataArrayIdType> revNodal[],
                       MCAuto<DataArrayIdType> revNodalIndx[],
                       std::vector< mcIdType >& facesByEdge,
                       std::vector< mcIdType >& equalFaces)
  {
    // find faces sharing the both nodes of edge

    facesByEdge.clear();
    size_t prevNbF; // nb faces found in 0-th mesh
    for ( int iM = 0; iM < 2; ++iM )
      {
        const mcIdType * revInd = revNodalIndx[ iM ]->begin();
        const mcIdType * rev    = revNodal    [ iM ]->begin();

        mcIdType nbRevFaces0 = revInd[ n0 + 1 ] - revInd[ n0 ];
        mcIdType nbRevFaces1 = revInd[ n1 + 1 ] - revInd[ n1 ];

        prevNbF = facesByEdge.size();
        facesByEdge.resize( prevNbF + std::max( nbRevFaces0, nbRevFaces1 ));

        auto it = std::set_intersection( rev + revInd[ n0 ],
                                         rev + revInd[ n0 ] + nbRevFaces0,
                                         rev + revInd[ n1 ],
                                         rev + revInd[ n1 ] + nbRevFaces1,
                                         facesByEdge.begin() + prevNbF );
        facesByEdge.resize( it - facesByEdge.begin() );
      }

    // facesByEdge now contains at least the 'inputFaceID'
    // check if there are other faces

    size_t nbF = facesByEdge.size();
    if ( nbF > 1 )
      {
        if ( prevNbF > 0 && prevNbF < nbF ) // faces found in both meshes
          {
            // remove from facesByEdge equal faces in different meshes
            const mcIdType *conn [2] = { mesh[0]->getNodalConnectivity()->getConstPointer(),
                                         mesh[1]->getNodalConnectivity()->getConstPointer() };
            const mcIdType *connI[2] = { mesh[0]->getNodalConnectivityIndex()->getConstPointer(),
                                         mesh[1]->getNodalConnectivityIndex()->getConstPointer() };
            for ( size_t i0 = 0; i0 < prevNbF; ++i0 )
              {
                if ( facesByEdge[ i0 ] == theUndefID )
                  continue;
                mcIdType objFaceID = MEDCouplingImpl::encodeID( facesByEdge[ i0 ], 0 );
                bool   isInputFace = ( objFaceID == inputFaceID );

                for ( size_t i1 = prevNbF; i1 < facesByEdge.size(); ++i1 )
                  {
                    if ( facesByEdge[ i1 ] == theUndefID )
                      continue;

                    mcIdType f0 = facesByEdge[ i0 ];
                    mcIdType f1 = facesByEdge[ i1 ];
                    size_t nbNodes0 = connI[0][ f0 + 1 ] - connI[0][ f0 ] - 1;
                    size_t nbNodes1 = connI[1][ f1 + 1 ] - connI[1][ f1 ] - 1;
                    if ( nbNodes0 != nbNodes1 )
                      continue;

                    const mcIdType * fConn0 = conn[0] + connI[0][ f0 ] + 1;
                    const mcIdType * fConn1 = conn[1] + connI[1][ f1 ] + 1;
                    if ( std::equal( fConn0, fConn0 + nbNodes0, fConn1 ))
                      {
                        // equal faces; remove an object one
                        mcIdType refFaceID = MEDCouplingImpl::encodeID( facesByEdge[ i1 ], 1 );
                        if ( refFaceID == inputFaceID )
                          isInputFace = true;

                        if ( std::find( equalFaces.begin(),
                                        equalFaces.end(), objFaceID ) == equalFaces.end() )
                          equalFaces.push_back( objFaceID );

                        facesByEdge[ i0 ] = theUndefID;
                        if ( isInputFace )
                          facesByEdge[ i1 ] = theUndefID;
                        break;
                      }
                  }
                if ( isInputFace )
                  facesByEdge[ i0 ] = theUndefID;
              }
          }
      }

    nbF = facesByEdge.size();
    for ( size_t i = 0; i < facesByEdge.size(); ++i )
      {
        if ( facesByEdge[ i ] != theUndefID )
          {
            facesByEdge[ i ] = MEDCouplingImpl::encodeID( facesByEdge[ i ], i >= prevNbF );
            if ( facesByEdge[ i ] == inputFaceID )
              facesByEdge[ i ] = theUndefID;
          }
        nbF -= ( facesByEdge[ i ] == theUndefID );
      }

    if ( nbF > 1 )
      return; // non-manifold

    if ( nbF < 1 )
      {
        facesByEdge.clear();
      }
    else // nbF == 1, set a found face first
      {
        if ( facesByEdge[ 0 ] == theUndefID )
          {
            for ( size_t i = 1; i < facesByEdge.size(); ++i )
              if ( facesByEdge[ i ] != theUndefID )
                {
                  facesByEdge[ 0 ] = facesByEdge[ i ];
                  break;
                }
          }
        facesByEdge.resize( 1 );
      }
    return;
  }

  //================================================================================
  /*!
   * \brief Remove a face from nodal reversed connectivity
   *  \param [in] node - a node of the face
   *  \param [in] face - the face
   *  \param [in.out] revNodal - reverse nodal connectivity
   *  \param [in,out] revNodalIndx - reverse nodal connectivity index
   */
  //================================================================================

  void removeFromRevNodal( mcIdType node,
                           mcIdType face,
                           MCAuto<DataArrayIdType>& revNodal,
                           MCAuto<DataArrayIdType>& revNodalIndx)
  {
    mcIdType* fBeg = revNodal->getPointer() + revNodalIndx->getIJ( node, 0 );
    mcIdType* fEnd = revNodal->getPointer() + revNodalIndx->getIJ( node + 1, 0);
    auto it = std::find( fBeg, fEnd, face );
    if ( it != fEnd )
      {
        for ( auto it2 = it + 1; it2 < fEnd; ++it2 ) // keep faces sorted
          *( it2 - 1 ) = *it2;

        *( fEnd - 1 ) = theUndefID;
      }
  }

  //================================================================================
  /*!
   * \brief Check order of two nodes in a given face
   *  \param [inout] n0 - node 1
   *  \param [inout] n1 - node 2
   *  \param [inout] iFEnc - face
   *  \param [inout] mesh - mesh
   *  \return bool - true if the nodes are in [ .., n1, n0, ..] order in face
   */
  //================================================================================

  bool isReverseOrder( mcIdType n0,
                       mcIdType n1,
                       mcIdType iFEnc,
                       MEDCouplingUMesh* mesh[] )
  {
    int iMesh;
    mcIdType iF = decodeID( iFEnc, iMesh );

    const mcIdType *conn  = mesh[ iMesh ]->getNodalConnectivity()->getConstPointer();
    const mcIdType *connI = mesh[ iMesh ]->getNodalConnectivityIndex()->getConstPointer();

    auto it0 = std::find( conn + connI[ iF ] + 1,
                          conn + connI[ iF + 1 ],
                          n0 );
    auto it1 = std::find( conn + connI[ iF ] + 1,
                          conn + connI[ iF + 1 ],
                          n1 );
    long i0 = it0 - conn;
    long i1 = it1 - conn;

    bool isRev = ( std::abs( i1 - i0 ) == 1 ) ?  i1 < i0 :  i0 < i1;
    return isRev;
  }

  //================================================================================
  /*!
   * \brief Change orientation of a face in one of given meshes
   *  \param [in] iFEnc - face ID also encoding a mesh index
   *  \param [in,out] mesh - object and reference meshes
   */
  //================================================================================

  void reverseFace( mcIdType iFEnc, MEDCouplingUMesh* mesh[] )
  {
    int iMesh;
    mcIdType face = decodeID( iFEnc, iMesh );

    mcIdType *conn  = mesh[ iMesh ]->getNodalConnectivity()->getPointer();
    mcIdType *connI = mesh[ iMesh ]->getNodalConnectivityIndex()->getPointer();

    const INTERP_KERNEL::CellModel& cm =
      INTERP_KERNEL::CellModel::GetCellModel( mesh[iMesh]->getTypeOfCell( face ));

    cm.changeOrientationOf2D( conn + connI[ face ] + 1,
                              (unsigned int)( connI[ face + 1 ] - connI[ face ] - 1 ));
    return;
  }
}

/// @endcond

//================================================================================
/*!
 * \brief Orient cells of \a this 2D mesh equally to \a refFaces
 *  \param [in] refFaces - 2D mesh containing correctly oriented faces. It is optional.
 *              If there are no cells in \a refFaces or it is nullptr, then any face
 *              in \a this mesh is used as a reference
 *  \throw If \a this mesh is not well defined.
 *  \throw If \a this mesh or \refFaces are not 2D.
 *  \throw If \a this mesh and \refFaces do not share nodes.
 *  \throw If \a refFaces are not equally oriented.
 *  \throw If \a this mesh plus \a refFaces together form a non-manifold mesh.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_are2DCellsNotCorrectlyOriented "Here is a C++ example".<br>
 *  \ref  py_mcumesh_are2DCellsNotCorrectlyOriented "Here is a Python example".
 *  \endif
 */
//================================================================================

void MEDCouplingUMesh::orientCorrectly2DCells(const MEDCouplingUMesh* refFaces)
{
  checkConsistencyLight();
  if ( getMeshDimension() != 2 )
    throw INTERP_KERNEL::Exception("The mesh dimension must be 2");
  if ( refFaces )
    {
      refFaces->checkConsistencyLight();
      if ( refFaces->getMeshDimension() != 2 )
        throw INTERP_KERNEL::Exception("The reference mesh dimension must be 2");
      if ( getCoords() != refFaces->getCoords() )
        throw INTERP_KERNEL::Exception("Object and reference meshes must share nodes ");
      if ( refFaces->getNumberOfCells() == 0 )
        refFaces = nullptr;
    }
  if ( getNumberOfCells() == 0 )
    return;

  enum { _OBJ, _REF };
  MEDCouplingUMesh* mesh[2] = { this, const_cast< MEDCouplingUMesh* >( refFaces ) };
  MCAuto<MEDCouplingUMesh> meshPtr;
  if ( !mesh[_REF] )
    {
      meshPtr = mesh[_REF] = MEDCouplingUMesh::New();
      mesh[_REF]->setCoords( mesh[_OBJ]->getCoords() );
      mesh[_REF]->allocateCells(0);
      mesh[_REF]->finishInsertingCells();
    }
  mcIdType nbFacesToCheck[2] = { mesh[_OBJ]->getNumberOfCells(),
                                 mesh[_REF]->getNumberOfCells() };
  std::vector< bool > isFaceQueued[ 2 ]; // enqueued faces of 2 meshes
  isFaceQueued[_OBJ].resize( nbFacesToCheck[_OBJ] );
  isFaceQueued[_REF].resize( nbFacesToCheck[_REF] );

  MCAuto<DataArrayIdType> revNodal    [2] = { DataArrayIdType::New(), DataArrayIdType::New() };
  MCAuto<DataArrayIdType> revNodalIndx[2] = { DataArrayIdType::New(), DataArrayIdType::New() };
  mesh[_OBJ]->getReverseNodalConnectivity( revNodal[_OBJ], revNodalIndx[_OBJ] );
  mesh[_REF]->getReverseNodalConnectivity( revNodal[_REF], revNodalIndx[_REF] );

  std::vector< mcIdType > faceNodes(4);
  std::vector< mcIdType > facesByEdge(4), equalFaces;
  std::vector< mcIdType > faceQueue; // starting faces with IDs counted from 1; negative ID mean a face in ref mesh

  while ( nbFacesToCheck[_OBJ] + nbFacesToCheck[_REF] > 0 ) // until all faces checked
    {
      if ( faceQueue.empty() ) // all neighbors checked, find more faces to check
        {
          for ( int iMesh = 1; iMesh >= 0; --iMesh ) // on [ _REF, _OBJ ]
            if ( nbFacesToCheck[iMesh] > 0 )
              for ( mcIdType f = 0, nbF = mesh[iMesh]->getNumberOfCells(); f < nbF; ++f )
                if ( !isFaceQueued[iMesh][f] )
                  {
                    faceQueue.push_back( MEDCouplingImpl::encodeID( f, iMesh ));
                    isFaceQueued[ iMesh ][ f ] = true;
                    iMesh = 0;
                    break;
                  }
          if ( faceQueue.empty() )
            break;
        }

      mcIdType fID = faceQueue.back();
      faceQueue.pop_back();

      int iMesh, iMesh2;
      mcIdType refFace = MEDCouplingImpl::decodeID( fID, iMesh );

      nbFacesToCheck[iMesh]--;

      equalFaces.clear();
      faceNodes.clear();
      mesh[iMesh]->getNodeIdsOfCell( refFace, faceNodes );
      const INTERP_KERNEL::CellModel& cm = INTERP_KERNEL::CellModel::GetCellModel( mesh[iMesh]->getTypeOfCell( refFace ));
      const int nbEdges = cm.getNumberOfSons();

      // loop on edges of the refFace
      mcIdType n0 = faceNodes[ nbEdges - 1 ]; // 1st node of edge
      for ( int edge = 0; edge < nbEdges; ++edge )
        {
          mcIdType n1 = faceNodes[ edge ]; // 2nd node of edge

          // get faces sharing the edge
          MEDCouplingImpl::getFacesOfEdge( n0, n1, fID, mesh, revNodal, revNodalIndx,
                                           facesByEdge, equalFaces );

          if ( facesByEdge.size() > 1 )
            THROW_IK_EXCEPTION("Non-manifold mesh at edge " << n0+1 << " - " << n1+1);

          if ( facesByEdge.size() == 1 )
            {
              // compare orientation of two faces
              //
              if ( !MEDCouplingImpl::isReverseOrder( n0, n1, facesByEdge[0], mesh ))
                {
                  if ( facesByEdge[0] < 0 ) // in the ref mesh
                    throw INTERP_KERNEL::Exception("Different orientation of reference faces");

                  MEDCouplingImpl::reverseFace( facesByEdge[0], mesh );
                }
              mcIdType face2 = MEDCouplingImpl::decodeID( facesByEdge[0], iMesh2 );
              if ( !isFaceQueued[iMesh2][face2] )
                {
                  isFaceQueued[iMesh2][face2] = true;
                  faceQueue.push_back( facesByEdge[0] );
                }
            }
          n0 = n1;
        }

      // remove face and equalFaces from revNodal in order not to treat them again
      equalFaces.push_back( fID );
      for ( mcIdType face : equalFaces )
        {
          mcIdType f            = MEDCouplingImpl::decodeID( face, iMesh2 );
          const mcIdType *conn  = mesh[iMesh2]->getNodalConnectivity()->getConstPointer();
          const mcIdType *connI = mesh[iMesh2]->getNodalConnectivityIndex()->getConstPointer();
          mcIdType nbNodes      = connI[ f + 1 ] - connI[ f ] - 1;
          for ( const mcIdType* n = conn + connI[ f ] + 1, *nEnd = n + nbNodes; n < nEnd; ++n )

            MEDCouplingImpl::removeFromRevNodal( *n, f,  // not to treat f again
                                                 revNodal[ iMesh2 ], revNodalIndx[ iMesh2 ] );
        }

    } // while() until all faces checked

  return;
}
