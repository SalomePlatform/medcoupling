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

#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingMemArray.txx"
#include "MEDCouplingFieldDouble.hxx"
#include "CellModel.hxx"
#include "VolSurfUser.txx"
#include "InterpolationUtils.hxx"
#include "PointLocatorAlgos.txx"
#include "BBTree.txx"
#include "DirectedBoundingBox.hxx"
#include "InterpKernelMeshQuality.hxx"
#include "InterpKernelCellSimplify.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "InterpKernelAutoPtr.hxx"
#include "InterpKernelGeo2DNode.hxx"
#include "InterpKernelGeo2DEdgeLin.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelGeo2DQuadraticPolygon.hxx"

#include <sstream>
#include <fstream>
#include <numeric>
#include <cstring>
#include <limits>
#include <list>

using namespace ParaMEDMEM;

const char MEDCouplingUMesh::PART_OF_NAME[]="PartOf_";

double MEDCouplingUMesh::EPS_FOR_POLYH_ORIENTATION=1.e-14;

const INTERP_KERNEL::NormalizedCellType MEDCouplingUMesh::MEDMEM_ORDER[N_MEDMEM_ORDER] = { INTERP_KERNEL::NORM_POINT1, INTERP_KERNEL::NORM_SEG2, INTERP_KERNEL::NORM_SEG3, INTERP_KERNEL::NORM_TRI3, INTERP_KERNEL::NORM_QUAD4, INTERP_KERNEL::NORM_TRI6, INTERP_KERNEL::NORM_QUAD8, INTERP_KERNEL::NORM_TETRA4, INTERP_KERNEL::NORM_PYRA5, INTERP_KERNEL::NORM_PENTA6, INTERP_KERNEL::NORM_HEXA8, INTERP_KERNEL::NORM_HEXGP12, INTERP_KERNEL::NORM_TETRA10, INTERP_KERNEL::NORM_PYRA13, INTERP_KERNEL::NORM_PENTA15, INTERP_KERNEL::NORM_HEXA20, INTERP_KERNEL::NORM_POLYL, INTERP_KERNEL::NORM_POLYGON, INTERP_KERNEL::NORM_QPOLYG, INTERP_KERNEL::NORM_POLYHED };

MEDCouplingUMesh *MEDCouplingUMesh::New()
{
  return new MEDCouplingUMesh;
}

MEDCouplingUMesh *MEDCouplingUMesh::New(const char *meshName, int meshDim)
{
  MEDCouplingUMesh *ret=new MEDCouplingUMesh;
  ret->setName(meshName);
  ret->setMeshDimension(meshDim);
  return ret;
}

MEDCouplingMesh *MEDCouplingUMesh::deepCpy() const
{
  return clone(true);
}

MEDCouplingUMesh *MEDCouplingUMesh::clone(bool recDeepCpy) const
{
  return new MEDCouplingUMesh(*this,recDeepCpy);
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

MEDCouplingUMesh::MEDCouplingUMesh():_iterator(-1),_mesh_dim(-2),
                                     _nodal_connec(0),_nodal_connec_index(0)
{
}

/*!
 * This method checks that this is correctly designed. For example le coordinates are set, nodal connectivity.
 * When this method returns without throwing any exception, 'this' is expected to be writable, exchangeable and to be 
 * available for most of algorithm. When a mesh has been constructed from scratch it is a good habits to call this method to check
 * that all is in order in 'this'.
 */
void MEDCouplingUMesh::checkCoherency() const throw(INTERP_KERNEL::Exception)
{
  if(_mesh_dim<-1)
    throw INTERP_KERNEL::Exception("No mesh dimension specified !");
  for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iter=_types.begin();iter!=_types.end();iter++)
    {
      if((int)INTERP_KERNEL::CellModel::GetCellModel(*iter).getDimension()!=_mesh_dim)
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
  if(_nodal_connec_index)
    {
      if(_nodal_connec_index->getNumberOfComponents()!=1)
        throw INTERP_KERNEL::Exception("Nodal connectivity index array is expected to be with number of components set to one !");
      if(_nodal_connec_index->getInfoOnComponent(0)!="")
        throw INTERP_KERNEL::Exception("Nodal connectivity index array is expected to have no info on its single component !");
    }
  if(_iterator!=-1)
    {
      throw INTERP_KERNEL::Exception("It appears that finishInsertingCells method has not been invoked after a insertNextCell session !");
    }
}

/*!
 * This method performs deeper checking in 'this' than MEDCouplingUMesh::checkCoherency does.
 * So this method is more time-consuming. This method checks that nodal connectivity points to valid node ids.
 * No geometrical aspects are checked here. These aspects are done in MEDCouplingUMesh::checkCoherency2.
 */
void MEDCouplingUMesh::checkCoherency1(double eps) const throw(INTERP_KERNEL::Exception)
{
  checkCoherency();
  if(_mesh_dim==-1)
    return ;
  int meshDim=getMeshDimension();
  int nbOfNodes=getNumberOfNodes();
  int nbOfCells=getNumberOfCells();
  const int *ptr=_nodal_connec->getConstPointer();
  const int *ptrI=_nodal_connec_index->getConstPointer();
  for(int i=0;i<nbOfCells;i++)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)ptr[ptrI[i]]);
      if((int)cm.getDimension()!=meshDim)
        {
          std::ostringstream oss;
          oss << "MEDCouplingUMesh::checkCoherency1 : cell << #" << i<< " with type Type " << cm.getRepr() << " in 'this' whereas meshdim == " << meshDim << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      int nbOfNodesInCell=ptrI[i+1]-ptrI[i]-1;
      if(!cm.isDynamic())
        if(nbOfNodesInCell!=(int)cm.getNumberOfNodes())
          {
            std::ostringstream oss;
            oss << "MEDCouplingUMesh::checkCoherency1 : cell #" << i << " with static Type '" << cm.getRepr() << "' has " <<  cm.getNumberOfNodes();
            oss << " nodes whereas in connectivity there is " << nbOfNodesInCell << " nodes ! Looks very bad !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      for(const int *w=ptr+ptrI[i]+1;w!=ptr+ptrI[i+1];w++)
        {
          int nodeId=*w;
          if(nodeId>=0)
            {
              if(nodeId>=nbOfNodes)
                {
                  std::ostringstream oss; oss << "Cell #" << i << " is consituted of node #" << nodeId << " whereas there are only " << nbOfNodes << " nodes !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
            }
          else if(nodeId<-1)
            {
              std::ostringstream oss; oss << "Cell #" << i << " is consituted of node #" << nodeId << " in connectivity ! sounds bad !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          else
            {
              if((INTERP_KERNEL::NormalizedCellType)(ptr[ptrI[i]])!=INTERP_KERNEL::NORM_POLYHED)
                {
                  std::ostringstream oss; oss << "Cell #" << i << " is consituted of node #-1 in connectivity ! sounds bad !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
            }
        }
    }
}

void MEDCouplingUMesh::checkCoherency2(double eps) const throw(INTERP_KERNEL::Exception)
{
  checkCoherency1(eps);
}

void MEDCouplingUMesh::setMeshDimension(int meshDim)
{
  if(meshDim<-1)
    throw INTERP_KERNEL::Exception("Invalid meshDim specified ! Must be greater or equal to -1 !");
  _mesh_dim=meshDim;
  declareAsNew();
}

void MEDCouplingUMesh::allocateCells(int nbOfCells)
{
  if(_nodal_connec_index)
    {
      _nodal_connec_index->decrRef();
    }
  if(_nodal_connec)
    {
      _nodal_connec->decrRef();
    }

  _nodal_connec_index=DataArrayInt::New();
  _nodal_connec_index->alloc(nbOfCells+1,1);
  int *pt=_nodal_connec_index->getPointer();
  pt[0]=0;
  _nodal_connec=DataArrayInt::New();
  _nodal_connec->alloc(2*nbOfCells,1);
  _iterator=0;
  _types.clear();
  declareAsNew();
}

/*!
 * Appends a cell in connectivity array.
 * @param type type of cell to add.
 * @param size number of nodes constituting this cell.
 * @param nodalConnOfCell the connectivity of the cell to add.
 */
void MEDCouplingUMesh::insertNextCell(INTERP_KERNEL::NormalizedCellType type, int size, const int *nodalConnOfCell) throw(INTERP_KERNEL::Exception)
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  if(_nodal_connec_index==0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::insertNextCell : nodal connectivity not set ! invoke allocateCells before calling insertNextCell !");
  if((int)cm.getDimension()==_mesh_dim)
    {
      int nbOfElems=_nodal_connec_index->getNbOfElems()-1;
      if(_iterator>=nbOfElems)
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::insertNextCell : allocation of cells was wide enough ! Call insertNextCell with higher value or call finishInsertingCells !");
      int *pt=_nodal_connec_index->getPointer();
      int idx=pt[_iterator];
      
      _nodal_connec->writeOnPlace(idx,type,nodalConnOfCell,size);
      _types.insert(type);
      pt[++_iterator]=idx+size+1;
    }
  else
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::insertNextCell : cell type " << cm.getRepr() << " has a dimension " << cm.getDimension();
      oss << " whereas Mesh Dimension of current UMesh instance is set to " << _mesh_dim << " ! Please invoke \"setMeshDimension\" method before or invoke ";
      oss << "\"MEDCouplingUMesh::New\" static method with 2 parameters name and meshDimension !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * Method to be called to cloture the insertion of cells using this->insertNextCell.
 */
void MEDCouplingUMesh::finishInsertingCells()
{
  const int *pt=_nodal_connec_index->getConstPointer();
  int idx=pt[_iterator];

  _nodal_connec->reAlloc(idx);
  _nodal_connec_index->reAlloc(_iterator+1);
  _iterator=-1;
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
 * If 'this' is not so that that cells are grouped by geo types this method will throw an exception.
 * In this case MEDCouplingUMesh::sortCellsInMEDFileFrmt or MEDCouplingUMesh::rearrange2ConsecutiveCellTypes methods for example can be called before invoking this method.
 * Useful for python users.
 */
MEDCouplingUMeshCellByTypeEntry *MEDCouplingUMesh::cellsByType() throw(INTERP_KERNEL::Exception)
{
  if(!checkConsecutiveCellTypes())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::cellsByType : this mesh is not sorted by type !");
  return new MEDCouplingUMeshCellByTypeEntry(this);
}

std::set<INTERP_KERNEL::NormalizedCellType> MEDCouplingUMesh::getAllGeoTypes() const
{
  return _types;
}

/*!
 * This method is a method that compares 'this' and 'other'.
 * This method compares \b all attributes, even names and component names.
 */
bool MEDCouplingUMesh::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const throw(INTERP_KERNEL::Exception)
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
 * This method looks if 'this' and 'other' are geometrically equivalent that is to say if each cell in 'other' correspond to one cell and only one
 * in 'this' is found regarding 'prec' parameter and 'cellCompPol' parameter.
 * 
 * In case of success cellCor and nodeCor are informed both. 
 * @param cellCompPol values are described in MEDCouplingUMesh::zipConnectivityTraducer method.
 * @param cellCor output array giving the correspondance of cells from 'other' to 'this'.
 * @param nodeCor output array giving the correspondance of nodes from 'other' to 'this'.
 */
void MEDCouplingUMesh::checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                            DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingUMesh *otherC=dynamic_cast<const MEDCouplingUMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("checkDeepEquivalWith : Two meshes are not not unstructured !");
  MEDCouplingMesh::checkFastEquivalWith(other,prec);
  if(_types!=otherC->_types)
    throw INTERP_KERNEL::Exception("checkDeepEquivalWith : Types are not equal !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m=MergeUMeshes(this,otherC);
  bool areNodesMerged;
  int newNbOfNodes;
  int oldNbOfNodes=getNumberOfNodes();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da=m->buildPermArrayForMergeNode(prec,oldNbOfNodes,areNodesMerged,newNbOfNodes);
  //mergeNodes
  if(!areNodesMerged)
    throw INTERP_KERNEL::Exception("checkDeepEquivalWith : Nodes are incompatible ! ");
  const int *pt=std::find_if(da->getConstPointer()+oldNbOfNodes,da->getConstPointer()+da->getNbOfElems(),std::bind2nd(std::greater<int>(),oldNbOfNodes-1));
  if(pt!=da->getConstPointer()+da->getNbOfElems())
    throw INTERP_KERNEL::Exception("checkDeepEquivalWith : some nodes in other are not in this !");
  m->renumberNodes(da->getConstPointer(),newNbOfNodes);
  //
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> nodeCor2=da->substr(oldNbOfNodes);
  da=m->mergeNodes(prec,areNodesMerged,newNbOfNodes);
  
  //
  da=m->zipConnectivityTraducer(cellCompPol);
  int maxId=*std::max_element(da->getConstPointer(),da->getConstPointer()+getNumberOfCells());
  pt=std::find_if(da->getConstPointer()+getNumberOfCells(),da->getConstPointer()+da->getNbOfElems(),std::bind2nd(std::greater<int>(),maxId));
  if(pt!=da->getConstPointer()+da->getNbOfElems())
    throw INTERP_KERNEL::Exception("checkDeepEquivalWith : some cells in other are not in this !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cellCor2=DataArrayInt::New();
  cellCor2->alloc(otherC->getNumberOfCells(),1);
  std::copy(da->getConstPointer()+getNumberOfCells(),da->getConstPointer()+da->getNbOfElems(),cellCor2->getPointer());
  bool nident=nodeCor2->isIdentity();
  bool cident=cellCor2->isIdentity();
  if(!nident) { nodeCor=nodeCor2; nodeCor2->incrRef(); } else nodeCor=0;
  if(!cident) { cellCor=cellCor2; cellCor2->incrRef(); } else cellCor=0;
}

/*!
 * This method looks if 'this' and 'other' are geometrically equivalent that is to say if each cell in 'other' correspond to one cell and only one
 * in 'this' is found regarding 'prec' parameter and 'cellCompPol' parameter. The difference with MEDCouplingUMesh::checkDeepEquivalWith method is that
 * coordinates of 'this' and 'other' are expected to be the same. If not an exception will be thrown.
 * This method is close to MEDCouplingUMesh::areCellsIncludedIn except that this method throws exception !
 * 
 * In case of success cellCor are informed both. 
 * @param cellCompPol values are described in MEDCouplingUMesh::zipConnectivityTraducer method.
 * @param cellCor output array giving the correspondance of cells from 'other' to 'this'.
 */
void MEDCouplingUMesh::checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                       DataArrayInt *&cellCor) const throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingUMesh *otherC=dynamic_cast<const MEDCouplingUMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("checkDeepEquivalOnSameNodesWith : Two meshes are not not unstructured !");
  MEDCouplingMesh::checkFastEquivalWith(other,prec);
  if(_types!=otherC->_types)
    throw INTERP_KERNEL::Exception("checkDeepEquivalOnSameNodesWith : Types are not equal !");
  if(_coords!=otherC->_coords)
    throw INTERP_KERNEL::Exception("checkDeepEquivalOnSameNodesWith : meshes do not share the same coordinates ! Use tryToShareSameCoordinates or call checkDeepEquivalWith !");
  std::vector<const MEDCouplingUMesh *> ms(2);
  ms[0]=this;
  ms[1]=otherC;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m=MergeUMeshesOnSameCoords(ms);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da=m->zipConnectivityTraducer(cellCompPol);
  int maxId=*std::max_element(da->getConstPointer(),da->getConstPointer()+getNumberOfCells());
  const int *pt=std::find_if(da->getConstPointer()+getNumberOfCells(),da->getConstPointer()+da->getNbOfElems(),std::bind2nd(std::greater<int>(),maxId));
  if(pt!=da->getConstPointer()+da->getNbOfElems())
    {
      throw INTERP_KERNEL::Exception("checkDeepEquivalOnSameNodesWith : some cells in other are not in this !");
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cellCor2=DataArrayInt::New();
  cellCor2->alloc(otherC->getNumberOfCells(),1);
  std::copy(da->getConstPointer()+getNumberOfCells(),da->getConstPointer()+da->getNbOfElems(),cellCor2->getPointer());
  if(!cellCor2->isIdentity()) { cellCor=cellCor2; cellCor2->incrRef(); } else cellCor=0;
}

/*!
 * This method checks fastly that 'this' and 'other' are equal. 
 */
void MEDCouplingUMesh::checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingUMesh *otherC=dynamic_cast<const MEDCouplingUMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("checkFastEquivalWith : Two meshes are not not unstructured !");
  MEDCouplingPointSet::checkFastEquivalWith(other,prec);
  int nbOfCells=getNumberOfCells();
  if(nbOfCells<1)
    return ;
  bool status=true;
  status&=areCellsFrom2MeshEqual(otherC,0,prec);
  status&=areCellsFrom2MeshEqual(otherC,nbOfCells/2,prec);
  status&=areCellsFrom2MeshEqual(otherC,nbOfCells-1,prec);
  if(!status)
    throw INTERP_KERNEL::Exception("checkFastEquivalWith : Two meshes are not equal because on 3 test cells some difference have been detected !");
}

/*!
 * \b WARNING this method do the assumption that connectivity lies on the coordinates set.
 * For speed reasons no check of this will be done.
 */
void MEDCouplingUMesh::getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  int nbOfNodes=getNumberOfNodes();
  int *revNodalIndxPtr=new int[nbOfNodes+1];
  revNodalIndx->useArray(revNodalIndxPtr,true,CPP_DEALLOC,nbOfNodes+1,1);
  std::fill(revNodalIndxPtr,revNodalIndxPtr+nbOfNodes+1,0);
  const int *conn=_nodal_connec->getConstPointer();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  int nbOfCells=getNumberOfCells();
  int nbOfEltsInRevNodal=0;
  for(int eltId=0;eltId<nbOfCells;eltId++)
    {
      const int *strtNdlConnOfCurCell=conn+connIndex[eltId]+1;
      const int *endNdlConnOfCurCell=conn+connIndex[eltId+1];
      for(const int *iter=strtNdlConnOfCurCell;iter!=endNdlConnOfCurCell;iter++)
        if(*iter>=0)//for polyhedrons
          {
            nbOfEltsInRevNodal++;
            revNodalIndxPtr[(*iter)+1]++;
          }
    }
  std::transform(revNodalIndxPtr+1,revNodalIndxPtr+nbOfNodes+1,revNodalIndxPtr,revNodalIndxPtr+1,std::plus<int>());
  int *revNodalPtr=new int[nbOfEltsInRevNodal];
  revNodal->useArray(revNodalPtr,true,CPP_DEALLOC,nbOfEltsInRevNodal,1);
  std::fill(revNodalPtr,revNodalPtr+nbOfEltsInRevNodal,-1);
  for(int eltId=0;eltId<nbOfCells;eltId++)
    {
      const int *strtNdlConnOfCurCell=conn+connIndex[eltId]+1;
      const int *endNdlConnOfCurCell=conn+connIndex[eltId+1];
      for(const int *iter=strtNdlConnOfCurCell;iter!=endNdlConnOfCurCell;iter++)
        if(*iter>=0)//for polyhedrons
          *std::find_if(revNodalPtr+revNodalIndxPtr[*iter],revNodalPtr+revNodalIndxPtr[*iter+1],std::bind2nd(std::equal_to<int>(),-1))=eltId;
    }
}

/// @cond INTERNAL

int MEDCouplingFastNbrer(int id, unsigned nb, const INTERP_KERNEL::CellModel& cm, bool compute, const int *conn1, const int *conn2)
{
  return id;
}

int MEDCouplingOrientationSensitiveNbrer(int id, unsigned nb, const INTERP_KERNEL::CellModel& cm, bool compute, const int *conn1, const int *conn2)
{
  if(!compute)
    return id+1;
  else
    {
      if(cm.getOrientationStatus(nb,conn1,conn2))
        return id+1;
      else
        return -(id+1);
    }
}

/// @endcond

/*!
 * \b WARNING this method do the assumption that connectivity lies on the coordinates set.
 * For speed reasons no check of this will be done.
 * Given 'this' with spacedim equal to s and meshdim equal to p, this method returns a new allocated mesh
 * lying on the same coordinates than 'this' and having a meshdim equal to p-1.
 * The algorithm to compute this p-1 mesh is the following :
 * For each cell in 'this' it splits into p-1 elements.
 *   If this p-1 element does not already exists it is appended to the returned mesh
 *   If this p-1 element already exists, it is not appended.
 * This method returns or 4 arrays plus the returned mesh.
 * 'desc' and 'descIndx' are the descending connectivity. These 2 arrays tell for each cell in 'this', to wich p-1 dimension cell in returned mesh it refers.
 * For a cell with a cellid c in 'this' it is constituted of cells in [desc+descIndx[c],desc+descIndex[c+1])
 *
 * Reversely 'revDesc' and 'revDescIndx' are the reverse descending connectivity. These 2 arrays tell for each cell in returned mesh, to wich cell in 'this' it refers.
 * For a cell with a cellid d in returned p-1 mesh it is shared by the following cells in 'this' [revDesc+revDescIndx[d],revDesc+revDescIndx[d+1])
 *
 * \warning This method returns a mesh whose geometric type order in are \b not sorted.
 * In view of the MED file writing, a renumbering of cells in returned mesh (using MEDCouplingUMesh::sortCellsInMEDFileFrmt) should be necessary.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildDescendingConnectivity(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const throw(INTERP_KERNEL::Exception)
{
  return buildDescendingConnectivityGen(desc,descIndx,revDesc,revDescIndx,MEDCouplingFastNbrer);
}

/*!
 * WARNING this method do the assumption that connectivity lies on the coordinates set.
 * For speed reasons no check of this will be done.
 * This method differs from MEDCouplingUMesh::buildDescendingConnectivity method in that 'desc' is in different format.
 * This method is more precise because it returns in descending connectivity giving the direction. If value is positive the n-1 dim element is taken in the same direction,
 * if it is in the opposite direction it is retrieved negative. So the problem is for elemt #0 in C convention. That's why this method is the only one that retrieves 
 * an array in relative "FORTRAN" mode.
 *
 * \warning This method returns a mesh whose geometric type order in are \b not sorted.
 * In view of the MED file writing, a renumbering of cells in returned mesh (using MEDCouplingUMesh::sortCellsInMEDFileFrmt) should be necessary.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildDescendingConnectivity2(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const throw(INTERP_KERNEL::Exception)
{
  return buildDescendingConnectivityGen(desc,descIndx,revDesc,revDescIndx,MEDCouplingOrientationSensitiveNbrer);
}

/*!
 * \b WARNING this method do the assumption that connectivity lies on the coordinates set.
 * For speed reasons no check of this will be done. This method calls MEDCouplingUMesh::buildDescendingConnectivity to compute the result.
 * This method lists cell by cell in \b this which are its neighbors. To compute the result only connectivities are considered.
 * The a cell with id 'cellId' its neighbors are neighbors[neighborsIndx[cellId]:neighborsIndx[cellId+1]].
 *
 * \param [out] neighbors is an array storing all the neighbors of all cells in \b this. This array is newly allocated and should be dealt by the caller. \b neighborsIndx 2nd output
 *                        parameter allows to select the right part in this array. The number of tuples is equal to the last values in \b neighborsIndx.
 * \param [out] neighborsIndx is an array of size this->getNumberOfCells()+1 newly allocated and should be dealt by the caller. This arrays allow to use the first output parameter \b neighbors.
 */
void MEDCouplingUMesh::computeNeighborsOfCells(DataArrayInt *&neighbors, DataArrayInt *&neighborsIndx) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> desc=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> descIndx=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> revDesc=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> revDescIndx=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> meshDM1=buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  meshDM1=0;
  ComputeNeighborsOfCellsAdv(desc,descIndx,revDesc,revDescIndx,neighbors,neighborsIndx);
}

/*!
 * This method is called by MEDCouplingUMesh::computeNeighborsOfCells. This methods performs the algorithm of MEDCouplingUMesh::computeNeighborsOfCells.
 * This method is useful for users that want to reduce along a criterion the set of neighbours cell. This is typically the case to extract a set a neighbours,
 * excluding a set of meshdim-1 cells in input descending connectivity.
 * Typically \b desc, \b descIndx, \b revDesc and \b revDescIndx input params are the result of MEDCouplingUMesh::buildDescendingConnectivity.
 * This method lists cell by cell in \b this which are its neighbors. To compute the result only connectivities are considered.
 * The a cell with id 'cellId' its neighbors are neighbors[neighborsIndx[cellId]:neighborsIndx[cellId+1]].
 *
 * \param [in] desc descending connectivity array.
 * \param [in] descIndx descending connectivity index array used to walk through \b desc.
 * \param [in] revDesc reverse descending connectivity array.
 * \param [in] revDescIndx reverse descending connectivity index array used to walk through \b revDesc.
 * \param [out] neighbors is an array storing all the neighbors of all cells in \b this. This array is newly allocated and should be dealt by the caller. \b neighborsIndx 2nd output
 *                        parameter allows to select the right part in this array. The number of tuples is equal to the last values in \b neighborsIndx.
 * \param [out] neighborsIndx is an array of size this->getNumberOfCells()+1 newly allocated and should be dealt by the caller. This arrays allow to use the first output parameter \b neighbors.
 */
void MEDCouplingUMesh::ComputeNeighborsOfCellsAdv(const DataArrayInt *desc, const DataArrayInt *descIndx, const DataArrayInt *revDesc, const DataArrayInt *revDescIndx,
                                                  DataArrayInt *&neighbors, DataArrayInt *&neighborsIndx) throw(INTERP_KERNEL::Exception)
{
  if(!desc || !descIndx || !revDesc || !revDescIndx)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ComputeNeighborsOfCellsAdv some input array is empty !");
  const int *descPtr=desc->getConstPointer();
  const int *descIPtr=descIndx->getConstPointer();
  const int *revDescPtr=revDesc->getConstPointer();
  const int *revDescIPtr=revDescIndx->getConstPointer();
  //
  int nbCells=descIndx->getNumberOfTuples()-1;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> out0=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> out1=DataArrayInt::New(); out1->alloc(nbCells+1,1);
  int *out1Ptr=out1->getPointer();
  *out1Ptr++=0;
  std::vector<int> out0v;
  out0v.reserve(desc->getNumberOfTuples());
  for(int i=0;i<nbCells;i++,descIPtr++,out1Ptr++)
    {
      for(const int *w1=descPtr+descIPtr[0];w1!=descPtr+descIPtr[1];w1++)
        {
          std::set<int> s(revDescPtr+revDescIPtr[*w1],revDescPtr+revDescIPtr[(*w1)+1]);
          s.erase(i);
          out0v.insert(out0v.end(),s.begin(),s.end());
        }
      *out1Ptr=out0v.size();
    }
  out0->alloc((int)out0v.size(),1);
  std::copy(out0v.begin(),out0v.end(),out0->getPointer());
  neighbors=out0; out0->incrRef();
  neighborsIndx=out1; out1->incrRef();
}

/// @cond INTERNAL

/*!
 * \b WARNING this method do the assumption that connectivity lies on the coordinates set.
 * For speed reasons no check of this will be done.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildDescendingConnectivityGen(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx, DimM1DescNbrer nbrer) const throw(INTERP_KERNEL::Exception)
{
  checkConnectivityFullyDefined();
  int nbOfCells=getNumberOfCells();
  int nbOfNodes=getNumberOfNodes();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  std::vector< std::vector<int> > descMeshConnB(nbOfCells);
  std::vector< std::vector<int> > revDescMeshConnB;
  std::vector< std::vector<int> > revNodalB(nbOfNodes);
  std::vector<int> meshDM1Conn;
  std::vector<int> meshDM1ConnIndex(1); meshDM1ConnIndex[0]=0;
  std::vector<int> meshDM1Type;
  for(int eltId=0;eltId<nbOfCells;eltId++)
    {
      int pos=connIndex[eltId];
      int posP1=connIndex[eltId+1];
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[pos]);
      unsigned nbOfSons=cm.getNumberOfSons2(conn+pos+1,posP1-pos-1);
      int *tmp=new int[posP1-pos];
      for(unsigned i=0;i<nbOfSons;i++)
        {
          INTERP_KERNEL::NormalizedCellType cmsId;
          unsigned nbOfNodesSon=cm.fillSonCellNodalConnectivity2(i,conn+pos+1,posP1-pos-1,tmp,cmsId);
          const INTERP_KERNEL::CellModel& cms=INTERP_KERNEL::CellModel::GetCellModel(cmsId);
          std::set<int> shareableCells(revNodalB[tmp[0]].begin(),revNodalB[tmp[0]].end());
          for(unsigned j=1;j<nbOfNodesSon && !shareableCells.empty();j++)
            {
              std::set<int> tmp2(revNodalB[tmp[j]].begin(),revNodalB[tmp[j]].end());
              std::set<int> tmp3;
              std::set_intersection(tmp2.begin(),tmp2.end(),shareableCells.begin(),shareableCells.end(),inserter(tmp3,tmp3.begin()));
              shareableCells=tmp3;
            }
          std::list<int> shareableCellsL(shareableCells.begin(),shareableCells.end());
          std::set<int> ref(tmp,tmp+nbOfNodesSon);
          for(std::list<int>::iterator iter=shareableCellsL.begin();iter!=shareableCellsL.end();)
            {
              if(cms.isCompatibleWith((INTERP_KERNEL::NormalizedCellType)meshDM1Type[*iter]))
                {
                  std::set<int> ref2(meshDM1Conn.begin()+meshDM1ConnIndex[*iter],meshDM1Conn.begin()+meshDM1ConnIndex[(*iter)+1]);
                  if(ref==ref2)
                    break;
                  else
                    iter=shareableCellsL.erase(iter);
                }
              else
                iter=shareableCellsL.erase(iter);
            }
          if(shareableCellsL.empty())
            {
              meshDM1Conn.insert(meshDM1Conn.end(),tmp,tmp+nbOfNodesSon);
              meshDM1ConnIndex.push_back(meshDM1ConnIndex.back()+nbOfNodesSon);
              int cellDM1Id=(int)meshDM1Type.size();
              meshDM1Type.push_back((int)cmsId);
              for(unsigned k=0;k<nbOfNodesSon;k++)
                revNodalB[tmp[k]].push_back(cellDM1Id);
              revDescMeshConnB.resize(cellDM1Id+1);
              revDescMeshConnB.back().push_back(eltId);
              descMeshConnB[eltId].push_back(nbrer(cellDM1Id,0,cms,false,0,0));
            }
          else
            {
              int DM1cellId=shareableCellsL.front();
              revDescMeshConnB[DM1cellId].push_back(eltId);
              descMeshConnB[eltId].push_back(nbrer(DM1cellId,nbOfNodesSon,cms,true,tmp,&meshDM1Conn[meshDM1ConnIndex[DM1cellId]]));
            }
        }
      delete [] tmp;
    }
  revNodalB.clear();
  //
  std::string name="Mesh constituent of "; name+=getName();
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New(name.c_str(),getMeshDimension()-1);
  ret->setCoords(getCoords());
  int nbOfCellsInConstituent=(int)meshDM1Type.size();
  ret->allocateCells(nbOfCellsInConstituent);
  revDescIndx->alloc(nbOfCellsInConstituent+1,1);
  int *tmp3=revDescIndx->getPointer(); tmp3[0]=0;
  for(int ii=0;ii<nbOfCellsInConstituent;ii++)
    {
      ret->insertNextCell((INTERP_KERNEL::NormalizedCellType)meshDM1Type[ii],meshDM1ConnIndex[ii+1]-meshDM1ConnIndex[ii],&meshDM1Conn[meshDM1ConnIndex[ii]]);
      tmp3[ii+1]=tmp3[ii]+((int)revDescMeshConnB[ii].size());
    }
  ret->finishInsertingCells();
  revDesc->alloc(tmp3[nbOfCellsInConstituent],1);
  tmp3=revDesc->getPointer();
  for(std::vector< std::vector<int> >::const_iterator iter2=revDescMeshConnB.begin();iter2!=revDescMeshConnB.end();iter2++)
    tmp3=std::copy((*iter2).begin(),(*iter2).end(),tmp3);
  meshDM1Type.clear(); meshDM1ConnIndex.clear(); meshDM1Conn.clear();
  descIndx->alloc(nbOfCells+1,1);
  tmp3=descIndx->getPointer(); tmp3[0]=0;
  for(int jj=0;jj<nbOfCells;jj++)
    tmp3[jj+1]=tmp3[jj]+((int)descMeshConnB[jj].size());
  desc->alloc(tmp3[nbOfCells],1);
  tmp3=desc->getPointer();
  for(std::vector< std::vector<int> >::const_iterator iter3=descMeshConnB.begin();iter3!=descMeshConnB.end();iter3++)
    tmp3=std::copy((*iter3).begin(),(*iter3).end(),tmp3);
  //
  return ret;
}

struct MEDCouplingAccVisit
{
  MEDCouplingAccVisit():_new_nb_of_nodes(0) { }
  int operator()(int val) { if(val!=-1) return _new_nb_of_nodes++; else return -1; }
  int _new_nb_of_nodes;
};

/// @endcond


/*!
 * This method convert cell with ids in ['cellIdsToConvertBg','cellIdsToConvertEnd') into 'this' into dynamic types without changing geometry.
 * That is to say if 'this' is a 2D, mesh after the invocation of this method it will contain only polygons.
 * If 'this' is a 3D mesh after the invocation of this method it will contain only polyhedra.
 * If mesh dimension is not in [2,3] an exception is thrown.
 * Of course pay attention that the resulting mesh is slower than previous one.
 * If in ['cellIdsToConvertBg','cellIdsToConvertEnd') there is a cell id not in [0,'this->getNumberOfCells()') an exception will be thrown.
 * In this case if meshDim==2 the mesh is still valid and only cells treated before throw will be converted into polygon.
 * If mesh==3, after throw the mesh is \b unconsistent !
 * This method is above all designed to test more extensively algorithms able to deal with polygons/polyhedra.
 * 
 * \warning This method modifies can modify significantly the geometric type order in \a this.
 * In view of the MED file writing, a renumbering of cells in \a this (using MEDCouplingUMesh::sortCellsInMEDFileFrmt) should be necessary.
 */
void MEDCouplingUMesh::convertToPolyTypes(const int *cellIdsToConvertBg, const int *cellIdsToConvertEnd)
{
  checkFullyDefined();
  int dim=getMeshDimension();
  if(dim<2 || dim>3)
    throw INTERP_KERNEL::Exception("Invalid mesh dimension : must be 2 or 3 !");
  int nbOfCells=getNumberOfCells();
  if(dim==2)
    {
      const int *connIndex=_nodal_connec_index->getConstPointer();
      int *conn=_nodal_connec->getPointer();
      for(const int *iter=cellIdsToConvertBg;iter!=cellIdsToConvertEnd;iter++)
        {
          if(*iter>=0 && *iter<nbOfCells)
            {
              const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[connIndex[*iter]]);
              if(!cm.isDynamic())
                conn[connIndex[*iter]]=INTERP_KERNEL::NORM_POLYGON;
              else
                conn[connIndex[*iter]]=INTERP_KERNEL::NORM_QPOLYG;
            }
          else
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::convertToPolyTypes : On rank #" << std::distance(cellIdsToConvertBg,iter) << " value is " << *iter << " which is not";
              oss << " in range [0," << nbOfCells << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
    }
  else
    {
      int *connIndex=_nodal_connec_index->getPointer();
      int connIndexLgth=_nodal_connec_index->getNbOfElems();
      const int *connOld=_nodal_connec->getConstPointer();
      int connOldLgth=_nodal_connec->getNbOfElems();
      std::vector<int> connNew(connOld,connOld+connOldLgth);
      for(const int *iter=cellIdsToConvertBg;iter!=cellIdsToConvertEnd;iter++)
        {
          if(*iter>=0 && *iter<nbOfCells)
            {
              int pos=connIndex[*iter];
              int posP1=connIndex[(*iter)+1];
              int lgthOld=posP1-pos-1;
              const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)connNew[pos]);
              connNew[pos]=INTERP_KERNEL::NORM_POLYHED;
              unsigned nbOfFaces=cm.getNumberOfSons2(&connNew[pos+1],lgthOld);
              int *tmp=new int[nbOfFaces*lgthOld];
              int *work=tmp;
              for(int j=0;j<(int)nbOfFaces;j++)
                {
                  INTERP_KERNEL::NormalizedCellType type;
                  unsigned offset=cm.fillSonCellNodalConnectivity2(j,&connNew[pos+1],lgthOld,work,type);
                  work+=offset;
                  *work++=-1;
                }
              std::size_t newLgth=std::distance(tmp,work)-1;
              std::size_t delta=newLgth-lgthOld;
              std::transform(connIndex+(*iter)+1,connIndex+connIndexLgth,connIndex+(*iter)+1,std::bind2nd(std::plus<int>(),delta));
              connNew.insert(connNew.begin()+posP1,tmp+lgthOld,tmp+newLgth);
              std::copy(tmp,tmp+lgthOld,connNew.begin()+pos+1);
              delete [] tmp;
            }
          else
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::convertToPolyTypes : On rank #" << std::distance(cellIdsToConvertBg,iter) << " value is " << *iter << " which is not";
              oss << " in range [0," << nbOfCells << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      _nodal_connec->alloc((int)connNew.size(),1);
      int *newConnPtr=_nodal_connec->getPointer();
      std::copy(connNew.begin(),connNew.end(),newConnPtr);
    }
  computeTypes();
}

/*!
 * This method converts all cells into poly type if possible.
 * This method is purely for userfriendliness.
 * As this method can be costly in Memory, no optimization is done to avoid construction of useless vector.
 */
void MEDCouplingUMesh::convertAllToPoly()
{
  int nbOfCells=getNumberOfCells();
  std::vector<int> cellIds(nbOfCells);
  for(int i=0;i<nbOfCells;i++)
    cellIds[i]=i;
  convertToPolyTypes(&cellIds[0],&cellIds[0]+cellIds.size());
}

/*!
 * This method expects that 'this' has a spacedim equal to 3 and a mesh dimension equal to 3 too, if not an exception will be thrown.
 * This method work only on cells with type NORM_POLYHED, all other cells with different type, are remains unchanged.
 * For such polyhedra, they are expected to have only 1 face (containing 2 faces in opposition), having 2*n number of nodes (n nodes on
 * each 2 faces hidden in the single face of polyhedron).
 * The first face is expected to be right oriented because all faces of this polyhedron will be deduced.
 * When called 'this' is an invalid mesh on MED sense. This method will correct that for polyhedra.
 * In case of presence of polyhedron that has not the extruded aspect (2 faces with the same number of nodes) an exception is thrown and 'this'
 * remains unchanged.
 * This method is usefull only for users that wants to build extruded unstructured mesh.
 * This method is a convenient one that avoids boring polyhedra setting during insertNextCell process.
 * In case of success, 'this' has be corrected contains the same number of cells and is valid in MED sense.
 */
void MEDCouplingUMesh::convertExtrudedPolyhedra() throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::convertExtrudedPolyhedra works on umeshes with meshdim equal to 3 and spaceDim equal to 3 too!");
  int nbOfCells=getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newCi=DataArrayInt::New();
  newCi->alloc(nbOfCells+1,1);
  int *newci=newCi->getPointer();
  const int *ci=_nodal_connec_index->getConstPointer();
  const int *c=_nodal_connec->getConstPointer();
  newci[0]=0;
  for(int i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)c[ci[i]];
      if(type==INTERP_KERNEL::NORM_POLYHED)
        {
          if(std::count(c+ci[i]+1,c+ci[i+1],-1)!=0)
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::convertExtrudedPolyhedra : cell # " << i << " is a polhedron BUT it has NOT exactly 1 face !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          std::size_t n2=std::distance(c+ci[i]+1,c+ci[i+1]);
          if(n2%2!=0)
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::convertExtrudedPolyhedra : cell # " << i << " is a polhedron with 1 face but there is a mismatch of number of nodes in face should be even !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          int n1=(int)(n2/2);
          newci[i+1]=7*n1+2+newci[i];//6*n1 (nodal length) + n1+2 (number of faces) - 1 (number of '-1' separator is equal to number of faces -1) + 1 (for cell type)
        }
      else
        newci[i+1]=(ci[i+1]-ci[i])+newci[i];
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newC=DataArrayInt::New();
  newC->alloc(newci[nbOfCells],1);
  int *newc=newC->getPointer();
  for(int i=0;i<nbOfCells;i++)
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
              newc[n1+5*j+2]=c[ci[i]+1+(j+1)%n1];
              newc[n1+5*j+3]=c[ci[i]+1+(j+1)%n1+n1];
              newc[n1+5*j+4]=c[ci[i]+1+j+n1];
            }
          newc+=n1*6;
        }
      else
        newc=std::copy(c+ci[i],c+ci[i+1],newc);
    }
  _nodal_connec_index->decrRef(); _nodal_connec_index=newCi;
  _nodal_connec->decrRef(); _nodal_connec=newC;
  newC->incrRef(); newCi->incrRef();
}

/*!
 * This method is the opposite of ParaMEDMEM::MEDCouplingUMesh::convertToPolyTypes method.
 * The aim is to take all polygons or polyhedrons cell and to try to traduce them into classical cells.
 *
 *  \return If true at least one cell has been unpolyzed.
            \n If false has been returned the nodal connectivity of \a this has **not** been altered and \a this has remains unchanged.
 *
 * \warning This method modifies can modify significantly the geometric type order in \a this.
 * In view of the MED file writing, a renumbering of cells in \a this (using MEDCouplingUMesh::sortCellsInMEDFileFrmt) should be necessary.
 */
bool MEDCouplingUMesh::unPolyze()
{
  checkFullyDefined();
  int mdim=getMeshDimension();
  if(mdim<0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::unPolyze works on umeshes with meshdim equals to 0, 1 2 or 3 !");
  if(mdim<=1)
    return false;
  int nbOfCells=getNumberOfCells();
  if(nbOfCells<1)
    return false;
  int initMeshLgth=getMeshLength();
  int *conn=_nodal_connec->getPointer();
  int *index=_nodal_connec_index->getPointer();
  int posOfCurCell=0;
  int newPos=0;
  int lgthOfCurCell;
  bool ret=false;
  for(int i=0;i<nbOfCells;i++)
    {
      lgthOfCurCell=index[i+1]-posOfCurCell;
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)conn[posOfCurCell];
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      INTERP_KERNEL::NormalizedCellType newType=INTERP_KERNEL::NORM_ERROR;
      int newLgth;
      if(cm.isDynamic())
        {
          switch(cm.getDimension())
            {
            case 2:
              {
                INTERP_KERNEL::AutoPtr<int> tmp=new int[lgthOfCurCell-1];
                std::copy(conn+posOfCurCell+1,conn+posOfCurCell+lgthOfCurCell,(int *)tmp);
                newType=INTERP_KERNEL::CellSimplify::tryToUnPoly2D(cm.isQuadratic(),tmp,lgthOfCurCell-1,conn+newPos+1,newLgth);
                break;
              }
            case 3:
              {
                int nbOfFaces,lgthOfPolyhConn;
                INTERP_KERNEL::AutoPtr<int> zipFullReprOfPolyh=INTERP_KERNEL::CellSimplify::getFullPolyh3DCell(type,conn+posOfCurCell+1,lgthOfCurCell-1,nbOfFaces,lgthOfPolyhConn);
                newType=INTERP_KERNEL::CellSimplify::tryToUnPoly3D(zipFullReprOfPolyh,nbOfFaces,lgthOfPolyhConn,conn+newPos+1,newLgth);
                break;
              }
            case 1:
              {
                newType=(lgthOfCurCell==3)?INTERP_KERNEL::NORM_SEG2:INTERP_KERNEL::NORM_POLYL;
                break;
              }
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
 * \param [in] eps is a relative precision that allows to establish if some 3D plane are coplanar or not. This epsilon is used to recenter around origin to have maximal 
 *             precision.
 */
void MEDCouplingUMesh::simplifyPolyhedra(double eps) throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::simplifyPolyhedra : works on meshdimension 3 and spaceDimension 3 !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> coords=getCoords()->deepCpy();
  coords->recenterForMaxPrecision(eps);
  const double *coordsPtr=coords->getConstPointer();
  //
  int nbOfCells=getNumberOfCells();
  const int *conn=_nodal_connec->getConstPointer();
  const int *index=_nodal_connec_index->getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> connINew=DataArrayInt::New();
  connINew->alloc(nbOfCells+1,1);
  int *connINewPtr=connINew->getPointer(); *connINewPtr++=0;
  std::vector<int> connNew;
  bool changed=false;
  for(int i=0;i<nbOfCells;i++,connINewPtr++)
    {
      if(conn[index[i]]==(int)INTERP_KERNEL::NORM_POLYHED)
        {
          SimplifyPolyhedronCell(eps,coords,conn+index[i],conn+index[i+1],connNew);
          changed=true;
        }
      else
        connNew.insert(connNew.end(),conn+index[i],conn+index[i+1]);
      *connINewPtr=(int)connNew.size();
    }
  if(changed)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> connNew2=DataArrayInt::New();
      connNew2->alloc((int)connNew.size(),1);
      std::copy(connNew.begin(),connNew.end(),connNew2->getPointer());
      setConnectivity(connNew2,connINew,false);
    }
}

/*!
 * This method returns all node ids used in \b this. The data array returned has to be dealt by the caller.
 * The returned node ids are sortes ascendingly. This method is closed to MEDCouplingUMesh::getNodeIdsInUse except
 * the format of returned DataArrayInt instance.
 * 
 * @return a newly allocated DataArrayInt sorted ascendingly of fetched node ids.
 * \sa MEDCouplingUMesh::getNodeIdsInUse
 */
DataArrayInt *MEDCouplingUMesh::computeFetchedNodeIds() const throw(INTERP_KERNEL::Exception)
{
  checkConnectivityFullyDefined();
  std::set<int> retS;
  int nbOfCells=getNumberOfCells();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  const int *conn=_nodal_connec->getConstPointer();
  for(int i=0;i<nbOfCells;i++)
    for(int j=connIndex[i]+1;j<connIndex[i+1];j++)
      if(conn[j]>=0)
        retS.insert(conn[j]);
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc((int)retS.size(),1);
  std::copy(retS.begin(),retS.end(),ret->getPointer());
  return ret;
}

/*!
 * Array returned is the correspondance in \b old \b to \b new format (that's why 'nbrOfNodesInUse' is returned too).
 * The returned array is newly created and should be dealt by the caller.
 * To retrieve the new to old format the user can use DataArrayInt::invertArrayO2N2N2O method.
 * The size of returned array is the number of nodes of 'this'.
 * -1 values in returned array means that the corresponding node never appear in any nodal connectivity of cells constituting 'this'.
 * @param [out] nbrOfNodesInUse out parameter that specifies how many of nodes in 'this' is really used in nodal connectivity.
 * @return a newly allocated DataArrayInt that tells for each nodeid in \b this if it is unused (-1) or used (the corresponding new id)
 */
DataArrayInt *MEDCouplingUMesh::getNodeIdsInUse(int& nbrOfNodesInUse) const throw(INTERP_KERNEL::Exception)
{
  nbrOfNodesInUse=-1;
  int nbOfNodes=getNumberOfNodes();
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbOfNodes,1);
  int *traducer=ret->getPointer();
  std::fill(traducer,traducer+nbOfNodes,-1);
  int nbOfCells=getNumberOfCells();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  const int *conn=_nodal_connec->getConstPointer();
  for(int i=0;i<nbOfCells;i++)
    for(int j=connIndex[i]+1;j<connIndex[i+1];j++)
      if(conn[j]>=0)
        traducer[conn[j]]=1;
  nbrOfNodesInUse=(int)std::count(traducer,traducer+nbOfNodes,1);
  std::transform(traducer,traducer+nbOfNodes,traducer,MEDCouplingAccVisit());
  return ret;
}

/*!
 * This method returns a newly allocated array containing this->getNumberOfCells() tuples and 1 component.
 * For each cell in \b this the number of nodes constituting cell is computed.
 * Excepted for poyhedrons, the result can be deduced by performing a deltaShiftIndex on the nodal connectivity index in \b this minus 1.
 * For polyhedrons, the face separation (-1) are excluded from the couting.
 * 
 * \return a newly allocated array
 */
DataArrayInt *MEDCouplingUMesh::computeNbOfNodesPerCell() const throw(INTERP_KERNEL::Exception)
{
  checkConnectivityFullyDefined();
  int nbOfCells=getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbOfCells,1);
  int *retPtr=ret->getPointer();
  const int *conn=getNodalConnectivity()->getConstPointer();
  const int *connI=getNodalConnectivityIndex()->getConstPointer();
  for(int i=0;i<nbOfCells;i++,retPtr++)
    {
      if(conn[connI[i]]!=(int)INTERP_KERNEL::NORM_POLYHED)
        *retPtr=connI[i+1]-connI[i]-1;
      else
        *retPtr=connI[i+1]-connI[i]-1-std::count(conn+connI[i]+1,conn+connI[i+1],-1);
    }
  ret->incrRef(); return ret;
}

/*!
 * Array returned is the correspondance in \b old \b to \b new format. The returned array is newly created and should be dealt by the caller.
 * The maximum value stored in returned array is the number of nodes of 'this' minus 1 after call of this method.
 * The size of returned array is the number of nodes of the old (previous to the call of this method) number of nodes.
 * -1 values in returned array means that the corresponding old node is no more used.
 */
DataArrayInt *MEDCouplingUMesh::zipCoordsTraducer() throw(INTERP_KERNEL::Exception)
{
  int newNbOfNodes=-1;
  DataArrayInt *traducer=getNodeIdsInUse(newNbOfNodes);
  renumberNodes(traducer->getConstPointer(),newNbOfNodes);
  return traducer;
}

/*!
 * This method stands if 'cell1' and 'cell2' are equals regarding 'compType' policy.
 * The semantic of 'compType' is specified in MEDCouplingUMesh::zipConnectivityTraducer method.
 */
int MEDCouplingUMesh::areCellsEqual(int cell1, int cell2, int compType) const
{
  switch(compType)
    {
    case 0:
      return areCellsEqual0(cell1,cell2);
    case 1:
      return areCellsEqual1(cell1,cell2);
    case 2:
      return areCellsEqual2(cell1,cell2);
    case 7:
      return areCellsEqual7(cell1,cell2);
    }
  throw INTERP_KERNEL::Exception("Unknown comparison asked ! Must be in 0,1 or 2.");
}

/*!
 * This method is the last step of the MEDCouplingUMesh::zipConnectivityTraducer with policy 0.
 */
int MEDCouplingUMesh::areCellsEqual0(int cell1, int cell2) const
{
  const int *conn=getNodalConnectivity()->getConstPointer();
  const int *connI=getNodalConnectivityIndex()->getConstPointer();
  if(connI[cell1+1]-connI[cell1]==connI[cell2+1]-connI[cell2])
    return std::equal(conn+connI[cell1]+1,conn+connI[cell1+1],conn+connI[cell2]+1)?1:0;
  return 0;
}

/*!
 * This method is the last step of the MEDCouplingUMesh::zipConnectivityTraducer with policy 1.
 */
int MEDCouplingUMesh::areCellsEqual1(int cell1, int cell2) const
{
  const int *conn=getNodalConnectivity()->getConstPointer();
  const int *connI=getNodalConnectivityIndex()->getConstPointer();
  int sz=connI[cell1+1]-connI[cell1];
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
                  int sz1=2*(sz-1);
                  int *tmp=new int[sz1];
                  int *work=std::copy(conn+connI[cell1]+1,conn+connI[cell1+1],tmp);
                  std::copy(conn+connI[cell1]+1,conn+connI[cell1+1],work);
                  work=std::search(tmp,tmp+sz1,conn+connI[cell2]+1,conn+connI[cell2+1]);
                  delete [] tmp;
                  return work!=tmp+sz1?1:0;
                }
              else
                return std::equal(conn+connI[cell1]+1,conn+connI[cell1+1],conn+connI[cell2]+1)?1:0;//case of SEG2 and SEG3
            }
          else
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::areCellsEqual1 : not implemented yet for meshdim == 3 !");
        }
    }
  return 0;
}

/*!
 * This method is the last step of the MEDCouplingUMesh::zipConnectivityTraducer with policy 2.
 */
int MEDCouplingUMesh::areCellsEqual2(int cell1, int cell2) const
{
  const int *conn=getNodalConnectivity()->getConstPointer();
  const int *connI=getNodalConnectivityIndex()->getConstPointer();
  if(connI[cell1+1]-connI[cell1]==connI[cell2+1]-connI[cell2])
    {
      if(conn[connI[cell1]]==conn[connI[cell2]])
        {
          std::set<int> s1(conn+connI[cell1]+1,conn+connI[cell1+1]);
          std::set<int> s2(conn+connI[cell2]+1,conn+connI[cell2+1]);
          return s1==s2?1:0;
        }
    }
  return 0;
}

/*!
 * This method is the last step of the MEDCouplingUMesh::zipConnectivityTraducer with policy 7.
 */
int MEDCouplingUMesh::areCellsEqual7(int cell1, int cell2) const
{
  const int *conn=getNodalConnectivity()->getConstPointer();
  const int *connI=getNodalConnectivityIndex()->getConstPointer();
  int sz=connI[cell1+1]-connI[cell1];
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
                  int sz1=2*(sz-1);
                  int *tmp=new int[sz1];
                  int *work=std::copy(conn+connI[cell1]+1,conn+connI[cell1+1],tmp);
                  std::copy(conn+connI[cell1]+1,conn+connI[cell1+1],work);
                  work=std::search(tmp,tmp+sz1,conn+connI[cell2]+1,conn+connI[cell2+1]);
                  if(work!=tmp+sz1)
                    {
                      delete [] tmp;
                      return 1;
                    }
                  else
                    {
                      std::reverse_iterator<int *> it1(tmp+sz1);
                      std::reverse_iterator<int *> it2(tmp);
                      if(std::search(it1,it2,conn+connI[cell2]+1,conn+connI[cell2+1])!=it2)
                        {
                          delete [] tmp;
                          return 2;
                        }
                      else
                        {
                          delete [] tmp;
                          return 0;
                        }
                    }
                  
                  return work!=tmp+sz1?1:0;
                }
              else
                {//case of SEG2 and SEG3
                  if(std::equal(conn+connI[cell1]+1,conn+connI[cell1+1],conn+connI[cell2]+1))
                    return 1;
                  if(!cm.isQuadratic())
                    {
                      std::reverse_iterator<const int *> it1(conn+connI[cell1+1]);
                      std::reverse_iterator<const int *> it2(conn+connI[cell1]+1);
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
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::areCellsEqual7 : not implemented yet for meshdim == 3 !");
        }
    }
  return 0;
}


/*!
 * This method compares 2 cells coming from two unstructured meshes : 'this' and 'other'.
 * This method compares 2 cells having the same id 'cellId' in 'this' and 'other'.
 */
bool MEDCouplingUMesh::areCellsFrom2MeshEqual(const MEDCouplingUMesh *other, int cellId, double prec) const
{
  if(getTypeOfCell(cellId)!=other->getTypeOfCell(cellId))
    return false;
  std::vector<int> c1,c2;
  getNodeIdsOfCell(cellId,c1);
  other->getNodeIdsOfCell(cellId,c2);
  std::size_t sz=c1.size();
  if(sz!=c2.size())
    return false;
  for(std::size_t i=0;i<sz;i++)
    {
      std::vector<double> n1,n2;
      getCoordinatesOfNode(c1[0],n1);
      other->getCoordinatesOfNode(c2[0],n2);
      std::transform(n1.begin(),n1.end(),n2.begin(),n1.begin(),std::minus<double>());
      std::transform(n1.begin(),n1.end(),n1.begin(),std::ptr_fun<double,double>(fabs));
      if(*std::max_element(n1.begin(),n1.end())>prec)
        return false;
    }
  return true;
}

/*!
 * This method find in candidate pool defined by 'candidates' the cells equal following the polycy 'compType'.
 * If any true is returned and the results will be put at the end of 'result' output parameter. If not false is returned
 * and result remains unchanged.
 * The semantic of 'compType' is specified in MEDCouplingUMesh::zipConnectivityTraducer method.
 * If in 'candidates' pool -1 value is considered as an empty value.
 * WARNING this method returns only ONE set of result !
 */
bool MEDCouplingUMesh::areCellsEqualInPool(const std::vector<int>& candidates, int compType, std::vector<int>& result) const
{
  std::set<int> cand(candidates.begin(),candidates.end());
  cand.erase(-1);
  if(cand.size()<=1)
    return false;
  bool ret=false;
  std::set<int>::const_iterator iter=cand.begin();
  int start=(*iter++);
  for(;iter!=cand.end();iter++)
    {
      int status=areCellsEqual(start,*iter,compType);
      if(status!=0)
        {
          if(!ret)
            {
              result.push_back(start);
              ret=true;
            }
          if(status==1)
            result.push_back(*iter);
          else
            result.push_back(status==2?(*iter+1):-(*iter+1));
        }
    }
  return ret;
}

/*!
 * This method common cells base regarding 'compType' comparison policy described in ParaMEDMEM::MEDCouplingUMesh::zipConnectivityTraducer for details.
 * This method returns 2 values 'res' and 'resI'.
 * If 'res' and 'resI' are not empty before calling this method they will be cleared before set.
 * The format of 'res' and 'resI' is as explained here.
 * resI.size()-1 is the number of set of cells equal.
 * The nth set is [res.begin()+resI[n];res.begin()+resI[n+1]) with 0<=n<resI.size()-1 
 */
template<int SPACEDIM>
void MEDCouplingUMesh::findCommonCellsBase(int compType, std::vector<int>& res, std::vector<int>& resI) const
{
  res.clear(); resI.clear();
  resI.push_back(0);
  std::vector<double> bbox;
  int nbOfCells=getNumberOfCells();
  getBoundingBoxForBBTree(bbox);
  double bb[2*SPACEDIM];
  double eps=getCaracteristicDimension();
  eps*=1.e-12;
  BBTree<SPACEDIM,int> myTree(&bbox[0],0,0,nbOfCells,-eps);
  const int *conn=getNodalConnectivity()->getConstPointer();
  const int *connI=getNodalConnectivityIndex()->getConstPointer();
  const double *coords=getCoords()->getConstPointer();
  std::vector<bool> isFetched(nbOfCells);
  for(int k=0;k<nbOfCells;k++)
    {
      if(!isFetched[k])
        {
          for(int j=0;j<SPACEDIM;j++)
            { bb[2*j]=std::numeric_limits<double>::max(); bb[2*j+1]=-std::numeric_limits<double>::max(); }
          for(const int *pt=conn+connI[k]+1;pt!=conn+connI[k+1];pt++)
            if(*pt>-1)
              {
                for(int j=0;j<SPACEDIM;j++)
                  {
                    bb[2*j]=std::min(bb[2*j],coords[SPACEDIM*(*pt)+j]);
                    bb[2*j+1]=std::max(bb[2*j+1],coords[SPACEDIM*(*pt)+j]);
                  }
              }
          std::vector<int> candidates1;
          myTree.getIntersectingElems(bb,candidates1);
          std::vector<int> candidates;
          for(std::vector<int>::const_iterator iter=candidates1.begin();iter!=candidates1.end();iter++)
            if(!isFetched[*iter])
              candidates.push_back(*iter);
          if(areCellsEqualInPool(candidates,compType,res))
            {
              int pos=resI.back();
              resI.push_back((int)res.size());
              for(std::vector<int>::const_iterator it=res.begin()+pos;it!=res.end();it++)
                isFetched[*it]=true;
            }
          isFetched[k]=true;
        }
    }
}

/*!
 * This method could potentially modify 'this'. This method merges cells if there are cells equal in 'this'. The comparison is specified by 'compType'.
 * This method keeps the coordiantes of 'this'.
 *
 * @param compType input specifying the technique used to compare cells each other.
 *   - 0 : exactly. A cell is detected to be the same if and only if the connectivity is exactly the same without permutation and types same too. This is the strongest policy.
 *   - 1 : permutation same orientation. cell1 and cell2 are considered equal if the connectivity of cell2 can be deduced by those of cell1 by direct permutation (with exactly the same orientation)
 * and their type equal. For 1D mesh the policy 1 is equivalent to 0.
 *   - 2 : nodal. cell1 and cell2 are equal if and only if cell1 and cell2 have same type and have the same nodes constituting connectivity. This is the laziest policy. This policy
 * can be used for users not sensitive to orientation of cell
 * @return the correspondance array old to new.
 * 
 * \warning This method modifies can modify significantly the geometric type order in \a this.
 * In view of the MED file writing, a renumbering of cells in \a this (using MEDCouplingUMesh::sortCellsInMEDFileFrmt) should be necessary.
 */
DataArrayInt *MEDCouplingUMesh::zipConnectivityTraducer(int compType) throw(INTERP_KERNEL::Exception)
{
  int spaceDim=getSpaceDimension();
  int nbOfCells=getNumberOfCells();
  std::vector<int> commonCells;
  std::vector<int> commonCellsI;
  switch(spaceDim)
    {
    case 3:
      {
        findCommonCellsBase<3>(compType,commonCells,commonCellsI);
        break;
      }
    case 2:
      {
        findCommonCellsBase<2>(compType,commonCells,commonCellsI);
        break;
      }
    case 1:
      {
        findCommonCellsBase<1>(compType,commonCells,commonCellsI);
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("Invalid spaceDimension : must be 1, 2 or 3.");
    }
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbOfCells,1);
  int *retPtr=ret->getPointer();
  std::fill(retPtr,retPtr+nbOfCells,0);
  const std::size_t nbOfTupleSmCells=commonCellsI.size()-1;
  int id=-1;
  std::vector<int> cellsToKeep;
  for(std::size_t i=0;i<nbOfTupleSmCells;i++)
    {
      for(std::vector<int>::const_iterator it=commonCells.begin()+commonCellsI[i];it!=commonCells.begin()+commonCellsI[i+1];it++)
        retPtr[*it]=id;
      id--;
    }
  id=0;
  std::map<int,int> m;
  for(int i=0;i<nbOfCells;i++)
    {
      int val=retPtr[i];
      if(val==0)
        {
          retPtr[i]=id++;
          cellsToKeep.push_back(i);
        }
      else
        {
          std::map<int,int>::const_iterator iter=m.find(val);
          if(iter==m.end())
            {
              m[val]=id;
              retPtr[i]=id++;
              cellsToKeep.push_back(i);
            }
          else
            retPtr[i]=(*iter).second;
        }
    }
  MEDCouplingUMesh *self=(MEDCouplingUMesh *)buildPartOfMySelf(&cellsToKeep[0],&cellsToKeep[0]+cellsToKeep.size(),true);
  setConnectivity(self->getNodalConnectivity(),self->getNodalConnectivityIndex(),true);
  self->decrRef();
  return ret;
}

/*!
 * This method makes the assumption that 'this' and 'other' share the same coords. If not an exception will be thrown !
 * This method tries to determine if 'other' is fully included in 'this'. To compute that, this method works with connectivity as MEDCouplingUMesh::zipConnectivityTraducer method does. 
 * This method is close to MEDCouplingUMesh::checkDeepEquivalOnSameNodesWith or MEDCouplingMesh::checkGeoEquivalWith with policy 20,21,or 22.
 * The main difference is that this method is not expected to throw exception.
 * This method has two outputs :
 *
 * @param compType is the comparison type. The possible values of this parameter are described in ParaMEDMEM::MEDCouplingUMesh::zipConnectivityTraducer method
 * @param arr is an output parameter that returns a \b newly created instance. This array is of size 'other->getNumberOfCells()'.
 * @return If 'other' is fully included in 'this 'true is returned. If not false is returned.
 */
bool MEDCouplingUMesh::areCellsIncludedIn(const MEDCouplingUMesh *other, int compType, DataArrayInt *& arr) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mesh=MergeUMeshesOnSameCoords(this,other);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> o2n=mesh->zipConnectivityTraducer(compType);
  int nbOfCells=getNumberOfCells();
  arr=o2n->substr(nbOfCells);
  arr->setName(other->getName());
  int tmp;
  if(other->getNumberOfCells()==0)
    return true;
  return arr->getMaxValue(tmp)<nbOfCells;
}

/*!
 * This method makes the assumption that 'this' and 'other' share the same coords. If not an exception will be thrown !
 * This method tries to determine if \b other is fully included in \b this.
 * The main difference is that this method is not expected to throw exception.
 * This method has two outputs :
 *
 * @param arr is an output parameter that returns a \b newly created instance. This array is of size 'other->getNumberOfCells()'.
 * @return If 'other' is fully included in 'this 'true is returned. If not false is returned.
 */
bool MEDCouplingUMesh::areCellsIncludedIn2(const MEDCouplingUMesh *other, DataArrayInt *& arr) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mesh=MergeUMeshesOnSameCoords(this,other);
  int spaceDim=mesh->getSpaceDimension();
  std::vector<int> commonCells;
  std::vector<int> commonCellsI;
  switch(spaceDim)
    {
    case 3:
      {
        findCommonCellsBase<3>(7,commonCells,commonCellsI);
        break;
      }
    case 2:
      {
        findCommonCellsBase<2>(7,commonCells,commonCellsI);
        break;
      }
    case 1:
      {
        findCommonCellsBase<1>(7,commonCells,commonCellsI);
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("Invalid spaceDimension : must be 1, 2 or 3.");
    }
  int thisNbCells=getNumberOfCells();
  int otherNbCells=other->getNumberOfCells();
  int nbOfCells=mesh->getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr2=DataArrayInt::New();
  arr2->alloc(otherNbCells,1);
  arr2->fillWithZero();
  int *arr2Ptr=arr2->getPointer();
  int nbOfCommon=(int)commonCellsI.size()-1;
  for(int i=0;i<nbOfCommon;i++)
    {
      int start=commonCells[commonCellsI[i]];
      if(start<thisNbCells)
        {
          for(int j=commonCellsI[i]+1;j!=commonCellsI[i+1];j++)
            {
              int sig=commonCells[j]>0?1:-1;
              int val=std::abs(commonCells[j])-1;
              if(val>=thisNbCells)
                arr2Ptr[val-thisNbCells]=sig*(start+1);
            }
        }
    }
  arr2->setName(other->getName());
  if(arr2->presenceOfValue(0))
    return false;
  arr=arr2;
  arr2->incrRef();
  return true;
}

/*!
 * @param areNodesMerged if at least two nodes have been merged.
 * @return old to new node correspondance.
 */
DataArrayInt *MEDCouplingUMesh::mergeNodes(double precision, bool& areNodesMerged, int& newNbOfNodes)
{
  DataArrayInt *ret=buildPermArrayForMergeNode(precision,-1,areNodesMerged,newNbOfNodes);
  if(areNodesMerged)
    renumberNodes(ret->getConstPointer(),newNbOfNodes);
  return ret;
}

/*!
 * Idem ParaMEDMEM::MEDCouplingUMesh::mergeNodes method except that the merged nodes are meld into the barycenter of them.
 */
DataArrayInt *MEDCouplingUMesh::mergeNodes2(double precision, bool& areNodesMerged, int& newNbOfNodes)
{
  DataArrayInt *ret=buildPermArrayForMergeNode(precision,-1,areNodesMerged,newNbOfNodes);
  if(areNodesMerged)
    renumberNodes2(ret->getConstPointer(),newNbOfNodes);
  return ret;
}

/*!
 * This method tries to use 'other' coords and use it for 'this'. If no exception was thrown after the call of this method :
 * this->_coords==other->_coords. If an exception is thrown 'this' remains unchanged.
 * Contrary to MEDCouplingUMesh::tryToShareSameCoords method this method makes a deeper analyze of coordinates (and so more expensive) than simple equality.
 * Two nodes one in 'this' and other in 'other' are considered equal if the distance between the two is lower than epsilon.
 */
void MEDCouplingUMesh::tryToShareSameCoordsPermute(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception)
{
  const DataArrayDouble *coords=other.getCoords();
  if(!coords)
    throw INTERP_KERNEL::Exception("tryToShareSameCoordsPermute : No coords specified in other !");
  if(!_coords)
    throw INTERP_KERNEL::Exception("tryToShareSameCoordsPermute : No coords specified in this whereas there is any in other !");
  int otherNbOfNodes=other.getNumberOfNodes();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> newCoords=MergeNodesArray(&other,this);
  _coords->incrRef();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> oldCoords=_coords;
  setCoords(newCoords);
  bool areNodesMerged;
  int newNbOfNodes;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da=buildPermArrayForMergeNode(epsilon,otherNbOfNodes,areNodesMerged,newNbOfNodes);
  if(!areNodesMerged)
    {
      setCoords(oldCoords);
      throw INTERP_KERNEL::Exception("tryToShareSameCoordsPermute fails : no nodes are mergeable with specified given epsilon !");
    }
  int maxId=*std::max_element(da->getConstPointer(),da->getConstPointer()+otherNbOfNodes);
  const int *pt=std::find_if(da->getConstPointer()+otherNbOfNodes,da->getConstPointer()+da->getNbOfElems(),std::bind2nd(std::greater<int>(),maxId));
  if(pt!=da->getConstPointer()+da->getNbOfElems())
    {
      setCoords(oldCoords);
      throw INTERP_KERNEL::Exception("tryToShareSameCoordsPermute fails : some nodes in this are not in other !");
    }
  setCoords(oldCoords);
  renumberNodesInConn(da->getConstPointer()+otherNbOfNodes);
  setCoords(coords);
}

/*!
 * Build a sub part of \b this lying or not on the same coordinates than \b this (regarding value of \b keepCoords).
 * By default coordinates are kept. This method is close to MEDCouplingUMesh::buildPartOfMySelf except that here input
 * cellIds is not given explicitely but by a range python like.
 * 
 * \param keepCoords that specifies if you want or not to keep coords as this or zip it (see ParaMEDMEM::MEDCouplingUMesh::zipCoords). If true zipCoords is \b NOT called, if false, zipCoords is called.
 * \return a newly allocated
 * 
 * \warning This method modifies can generate an unstructured mesh whose cells are not sorted by geometric type order.
 * In view of the MED file writing, a renumbering of cells of returned unstructured mesh (using MEDCouplingUMesh::sortCellsInMEDFileFrmt) should be necessary.
 */
MEDCouplingPointSet *MEDCouplingUMesh::buildPartOfMySelf2(int start, int end, int step, bool keepCoords) const throw(INTERP_KERNEL::Exception)
{
  if(getMeshDimension()!=-1)
    {
      MEDCouplingUMesh *ret=buildPartOfMySelfKeepCoords2(start,end,step);
      if(!keepCoords)
        ret->zipCoords();
      return ret;
    }
  else
    {
      int newNbOfCells=DataArray::GetNumberOfItemGivenBESRelative(start,end,step,"MEDCouplingUMesh::buildPartOfMySelf2 for -1 dimension mesh ");
      if(newNbOfCells!=1)
        throw INTERP_KERNEL::Exception("-1D mesh has only one cell !");
      if(start!=0)
        throw INTERP_KERNEL::Exception("-1D mesh has only one cell : 0 !");
      incrRef();
      return const_cast<MEDCouplingUMesh *>(this);
    }
}

/*!
 * build a sub part of \b this. This sub part is defined by the cell ids contained in the array in [begin,end).
 * @param begin begin of array containing the cell ids to keep.
 * @param end end of array of cell ids to keep. \b WARNING end param is \b not included ! Idem STL standard definitions.
 * @param keepCoords that specifies if you want or not to keep coords as this or zip it (see ParaMEDMEM::MEDCouplingUMesh::zipCoords). If true zipCoords is \b NOT called, if false, zipCoords is called.
 * 
 * \warning This method modifies can generate an unstructured mesh whose cells are not sorted by geometric type order.
 * In view of the MED file writing, a renumbering of cells of returned unstructured mesh (using MEDCouplingUMesh::sortCellsInMEDFileFrmt) should be necessary.
 */
MEDCouplingPointSet *MEDCouplingUMesh::buildPartOfMySelf(const int *begin, const int *end, bool keepCoords) const
{
  if(getMeshDimension()!=-1)
    {
      MEDCouplingUMesh *ret=buildPartOfMySelfKeepCoords(begin,end);
      if(!keepCoords)
        ret->zipCoords();
      return ret;
    }
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
 * This method allows to partially modify some cells in \b this (whose list is specified by [\b cellIdsBg, \b cellIdsEnd) ) with cells coming in \b otherOnSameCoordsThanThis.
 * Size of [\b cellIdsBg, \b cellIdsEnd) ) must be equal to the number of cells of otherOnSameCoordsThanThis.
 * The number of cells of \b this will remain the same with this method.
 *
 * \param [in] begin begin of cell ids (included) of cells in this to assign
 * \param [in] end end of cell ids (excluded) of cells in this to assign
 * \param [in] otherOnSameCoordsThanThis an another mesh with same meshdimension than \b this with exactly the same number of cells than cell ids list in [\b cellIdsBg, \b cellIdsEnd).
 *             Coordinate pointer of \b this and those of \b otherOnSameCoordsThanThis must be the same
 */
void MEDCouplingUMesh::setPartOfMySelf(const int *cellIdsBg, const int *cellIdsEnd, const MEDCouplingUMesh& otherOnSameCoordsThanThis) throw(INTERP_KERNEL::Exception)
{
  checkConnectivityFullyDefined();
  otherOnSameCoordsThanThis.checkConnectivityFullyDefined();
  if(getCoords()!=otherOnSameCoordsThanThis.getCoords())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::setPartOfMySelf : coordinates pointer are not the same ! Invoke setCoords or call tryToShareSameCoords method !");
  if(getMeshDimension()!=otherOnSameCoordsThanThis.getMeshDimension())
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::setPartOfMySelf : Mismatch of meshdimensions ! this is equal to " << getMeshDimension();
      oss << ", whereas other mesh dimension is set equal to " << otherOnSameCoordsThanThis.getMeshDimension() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  int nbOfCellsToModify=(int)std::distance(cellIdsBg,cellIdsEnd);
  if(nbOfCellsToModify!=otherOnSameCoordsThanThis.getNumberOfCells())
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::setPartOfMySelf : cells ids length (" <<  nbOfCellsToModify << ") do not match the number of cells of other mesh (" << otherOnSameCoordsThanThis.getNumberOfCells() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  int nbOfCells=getNumberOfCells();
  bool easyAssign=true;
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const int *connOther=otherOnSameCoordsThanThis._nodal_connec->getConstPointer();
  const int *connIOther=otherOnSameCoordsThanThis._nodal_connec_index->getConstPointer();
  for(const int *it=cellIdsBg;it!=cellIdsEnd && easyAssign;it++,connIOther++)
    {
      if(*it>=0 && *it<nbOfCells)
        {
          easyAssign=(connIOther[1]-connIOther[0])==(connI[*it+1]-connI[*it]);
        }
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::setPartOfMySelf : On pos #" << std::distance(cellIdsBg,it) << " id is equal to " << *it << " which is not in [0," << nbOfCells << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  if(easyAssign)
    {
      MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx(cellIdsBg,cellIdsEnd,_nodal_connec,_nodal_connec_index,otherOnSameCoordsThanThis._nodal_connec,otherOnSameCoordsThanThis._nodal_connec_index);
      computeTypes();
    }
  else
    {
      DataArrayInt *arrOut=0,*arrIOut=0;
      MEDCouplingUMesh::SetPartOfIndexedArrays(cellIdsBg,cellIdsEnd,_nodal_connec,_nodal_connec_index,otherOnSameCoordsThanThis._nodal_connec,otherOnSameCoordsThanThis._nodal_connec_index,
                                               arrOut,arrIOut);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arrOutAuto(arrOut),arrIOutAuto(arrIOut);
      setConnectivity(arrOut,arrIOut,true);
    }
}

void MEDCouplingUMesh::setPartOfMySelf2(int start, int end, int step, const MEDCouplingUMesh& otherOnSameCoordsThanThis) throw(INTERP_KERNEL::Exception)
{
  checkConnectivityFullyDefined();
  otherOnSameCoordsThanThis.checkConnectivityFullyDefined();
  if(getCoords()!=otherOnSameCoordsThanThis.getCoords())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::setPartOfMySelf2 : coordinates pointer are not the same ! Invoke setCoords or call tryToShareSameCoords method !");
  if(getMeshDimension()!=otherOnSameCoordsThanThis.getMeshDimension())
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::setPartOfMySelf2 : Mismatch of meshdimensions ! this is equal to " << getMeshDimension();
      oss << ", whereas other mesh dimension is set equal to " << otherOnSameCoordsThanThis.getMeshDimension() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  int nbOfCellsToModify=DataArray::GetNumberOfItemGivenBESRelative(start,end,step,"MEDCouplingUMesh::setPartOfMySelf2 : ");
  if(nbOfCellsToModify!=otherOnSameCoordsThanThis.getNumberOfCells())
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::setPartOfMySelf2 : cells ids length (" <<  nbOfCellsToModify << ") do not match the number of cells of other mesh (" << otherOnSameCoordsThanThis.getNumberOfCells() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  int nbOfCells=getNumberOfCells();
  bool easyAssign=true;
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const int *connOther=otherOnSameCoordsThanThis._nodal_connec->getConstPointer();
  const int *connIOther=otherOnSameCoordsThanThis._nodal_connec_index->getConstPointer();
  int it=start;
  for(int i=0;i<nbOfCellsToModify && easyAssign;i++,it+=step,connIOther++)
    {
      if(it>=0 && it<nbOfCells)
        {
          easyAssign=(connIOther[1]-connIOther[0])==(connI[it+1]-connI[it]);
        }
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::setPartOfMySelf2 : On pos #" << i << " id is equal to " << it << " which is not in [0," << nbOfCells << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  if(easyAssign)
    {
      MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx2(start,end,step,_nodal_connec,_nodal_connec_index,otherOnSameCoordsThanThis._nodal_connec,otherOnSameCoordsThanThis._nodal_connec_index);
      computeTypes();
    }
  else
    {
      DataArrayInt *arrOut=0,*arrIOut=0;
      MEDCouplingUMesh::SetPartOfIndexedArrays2(start,end,step,_nodal_connec,_nodal_connec_index,otherOnSameCoordsThanThis._nodal_connec,otherOnSameCoordsThanThis._nodal_connec_index,
                                                arrOut,arrIOut);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arrOutAuto(arrOut),arrIOutAuto(arrIOut);
      setConnectivity(arrOut,arrIOut,true);
    }
}                      

DataArrayInt *MEDCouplingUMesh::getCellIdsFullyIncludedInNodeIds(const int *partBg, const int *partEnd) const
{
  std::vector<int> cellIdsKept;
  fillCellIdsToKeepFromNodeIds(partBg,partEnd,true,cellIdsKept);
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc((int)cellIdsKept.size(),1);
  std::copy(cellIdsKept.begin(),cellIdsKept.end(),ret->getPointer());
  return ret;
}

/*!
 * Keeps from 'this' only cells which constituing point id are in the ids specified by ['begin','end').
 * The resulting cell ids are stored at the end of the 'cellIdsKept' parameter.
 * Parameter 'fullyIn' specifies if a cell that has part of its nodes in ids array is kept or not.
 * If 'fullyIn' is true only cells whose ids are \b fully contained in ['begin','end') tab will be kept.
 *
 * @param begin input start of array of node ids.
 * @param end input end of array of node ids.
 * @param fullyIn input that specifies if all node ids must be in ['begin','end') array to consider cell to be in.
 * @param cellIdsKept in/out array where all candidate cell ids are put at the end.
 */
void MEDCouplingUMesh::fillCellIdsToKeepFromNodeIds(const int *begin, const int *end, bool fullyIn, std::vector<int>& cellIdsKept) const
{
  std::set<int> fastFinder(begin,end);
  int nbOfCells=getNumberOfCells();
  const int *conn=getNodalConnectivity()->getConstPointer();
  const int *connIndex=getNodalConnectivityIndex()->getConstPointer();
  for(int i=0;i<nbOfCells;i++)
    {
      std::set<int> connOfCell(conn+connIndex[i]+1,conn+connIndex[i+1]);
      connOfCell.erase(-1);//polyhedron separator
      int refLgth=(int)connOfCell.size();
      std::set<int> locMerge;
      std::insert_iterator< std::set<int> > it(locMerge,locMerge.begin());
      std::set_intersection(connOfCell.begin(),connOfCell.end(),fastFinder.begin(),fastFinder.end(),it);
      if(((int)locMerge.size()==refLgth && fullyIn) || (locMerge.size()!=0 && !fullyIn))
        cellIdsKept.push_back(i);
    }
}

/*!
 * This method is very close too MEDCouplingUMesh::buildPartOfMySelfNode. The difference is that it returns directly ids.
 */
DataArrayInt *MEDCouplingUMesh::getCellIdsLyingOnNodes(const int *begin, const int *end, bool fullyIn) const
{
  std::vector<int> cellIdsKept;
  fillCellIdsToKeepFromNodeIds(begin,end,fullyIn,cellIdsKept);
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc((int)cellIdsKept.size(),1);
  std::copy(cellIdsKept.begin(),cellIdsKept.end(),ret->getPointer());
  ret->setName(getName());
  return ret;
}

/*!
 * Keeps from 'this' only cells which constituing point id are in the ids specified by ['begin','end').
 * The return newly allocated mesh will share the same coordinates as 'this'.
 * Parameter 'fullyIn' specifies if a cell that has part of its nodes in ids array is kept or not.
 * If 'fullyIn' is true only cells whose ids are \b fully contained in ['begin','end') tab will be kept.
 */
MEDCouplingPointSet *MEDCouplingUMesh::buildPartOfMySelfNode(const int *begin, const int *end, bool fullyIn) const
{
  std::vector<int> cellIdsKept;
  fillCellIdsToKeepFromNodeIds(begin,end,fullyIn,cellIdsKept);
  return buildPartOfMySelf(&cellIdsKept[0],&cellIdsKept[0]+cellIdsKept.size(),true);
}

/*!
 * Contrary to MEDCouplingUMesh::buildPartOfMySelfNode method this method builds a mesh with a meshDimension equal to
 * this->getMeshDimension()-1. The return newly allocated mesh will share the same coordinates as 'this'.
 * Parameter 'fullyIn' specifies if a face that has part of its nodes in ids array is kept or not.
 * If 'fullyIn' is true only faces whose ids are \b fully contained in ['begin','end') tab will be kept.
 */
MEDCouplingPointSet *MEDCouplingUMesh::buildFacePartOfMySelfNode(const int *begin, const int *end, bool fullyIn) const
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> desc,descIndx,revDesc,revDescIndx;
  desc=DataArrayInt::New(); descIndx=DataArrayInt::New(); revDesc=DataArrayInt::New(); revDescIndx=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> subMesh=buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  desc=0; descIndx=0; revDesc=0; revDescIndx=0;
  return subMesh->buildPartOfMySelfNode(begin,end,fullyIn);
}

/*!
 * This method returns a mesh with meshDim=this->getMeshDimension()-1.
 * This returned mesh contains cells that are linked with one and only one cell of this.
 * @param keepCoords specifies if ParaMEDMEM::MEDCouplingUMesh::zipCoords is called on returned mesh before being returned. If true zipCoords is \b NOT called, if false, zipCoords is called.
 * @return mesh with ref counter equal to 1.
 */
MEDCouplingPointSet *MEDCouplingUMesh::buildBoundaryMesh(bool keepCoords) const
{
  DataArrayInt *desc=DataArrayInt::New();
  DataArrayInt *descIndx=DataArrayInt::New();
  DataArrayInt *revDesc=DataArrayInt::New();
  DataArrayInt *revDescIndx=DataArrayInt::New();
  //
  MEDCouplingUMesh *meshDM1=buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  revDesc->decrRef();
  desc->decrRef();
  descIndx->decrRef();
  int nbOfCells=meshDM1->getNumberOfCells();
  const int *revDescIndxC=revDescIndx->getConstPointer();
  std::vector<int> boundaryCells;
  for(int i=0;i<nbOfCells;i++)
    if(revDescIndxC[i+1]-revDescIndxC[i]==1)
      boundaryCells.push_back(i);
  revDescIndx->decrRef();
  MEDCouplingPointSet *ret=meshDM1->buildPartOfMySelf(&boundaryCells[0],&boundaryCells[0]+boundaryCells.size(),keepCoords);
  meshDM1->decrRef();
  return ret;
}

/*!
 * This method returns a newly created DataArrayInt instance containing ids of cells located in boundary.
 * A cell is detected to be on boundary if it contains one or more than one face having only one father.
 * This method makes the assumption that 'this' is fully defined (coords,connectivity). If not an exception will be thrown. 
 */
DataArrayInt *MEDCouplingUMesh::findCellIdsOnBoundary() const throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  DataArrayInt *desc=DataArrayInt::New();
  DataArrayInt *descIndx=DataArrayInt::New();
  DataArrayInt *revDesc=DataArrayInt::New();
  DataArrayInt *revDescIndx=DataArrayInt::New();
  //
  MEDCouplingUMesh *meshDM1=buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  meshDM1->decrRef();
  desc->decrRef();
  descIndx->decrRef();
  //
  DataArrayInt *tmp=revDescIndx->deltaShiftIndex();
  DataArrayInt *faceIds=tmp->getIdsEqual(1);
  tmp->decrRef();
  int nbOfFaces=faceIds->getNumberOfTuples();
  const int *faces=faceIds->getConstPointer();
  std::set<int> ret;
  for(const int *w=faces;w!=faces+nbOfFaces;w++)
    ret.insert(revDesc->getIJ(revDescIndx->getIJ(*w,0),0));
  faceIds->decrRef();
  //
  revDescIndx->decrRef();
  revDesc->decrRef();
  //
  DataArrayInt *ret2=DataArrayInt::New();
  ret2->alloc((int)ret.size(),1);
  std::copy(ret.begin(),ret.end(),ret2->getPointer());
  ret2->setName("BoundaryCells");
  return ret2;
}

/*!
 * This method find in \b this cells ids that lie on mesh \b otherDimM1OnSameCoords.
 * \b this and \b otherDimM1OnSameCoords have to lie on the same coordinate array pointer. The coherency of that coords array with connectivity
 * of \b this and \b otherDimM1OnSameCoords is not important here because this method works only on connectivity.
 * this->getMeshDimension() - 1 must be equal to otherDimM1OnSameCoords.getMeshDimension()
 *
 * s0 is the cells ids set in \b this lying on at least one node in fetched nodes in \b otherDimM1OnSameCoords.
 * This method method returns cells ids set s = s1 + s2 where :
 * 
 *  - s1 are cells ids in \b this whose dim-1 constituent equals a cell in \b otherDimM1OnSameCoords.
 *  - s2 are cells ids in \b s0 - \b s1 whose at least two neighbors are in s1.
 *
 * \throw if \b otherDimM1OnSameCoords is not part of constituent of \b this, or if coordinate pointer of \b this and \b otherDimM1OnSameCoords
 *        are not same, or if this->getMeshDimension()-1!=otherDimM1OnSameCoords.getMeshDimension()
 *
 * \param [out] cellIdsRk0 a newly allocated array containing cells ids in \b this containg s0 in above algorithm.
 * \param [out] cellIdsRk1 a newly allocated array containing cells ids of s1+s2 \b into \b cellIdsRk0 subset. To get absolute ids of s1+s2 simply invoke
 *              cellIdsRk1->transformWithIndArr(cellIdsRk0->begin(),cellIdsRk0->end());
 */
void MEDCouplingUMesh::findCellIdsLyingOn(const MEDCouplingUMesh& otherDimM1OnSameCoords, DataArrayInt *&cellIdsRk0, DataArrayInt *&cellIdsRk1) const throw(INTERP_KERNEL::Exception)
{
  if(getCoords()!=otherDimM1OnSameCoords.getCoords())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findCellIdsLyingOn : coordinates pointer are not the same ! Use tryToShareSameCoords method !");
  checkConnectivityFullyDefined();
  otherDimM1OnSameCoords.checkConnectivityFullyDefined();
  if(getMeshDimension()-1!=otherDimM1OnSameCoords.getMeshDimension())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findCellIdsLyingOn : invalid mesh dimension of input mesh regarding meshdimesion of this !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> fetchedNodeIds1=otherDimM1OnSameCoords.computeFetchedNodeIds();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> s0arr=getCellIdsLyingOnNodes(fetchedNodeIds1->begin(),fetchedNodeIds1->end(),false);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> thisPart=static_cast<MEDCouplingUMesh *>(buildPartOfMySelf(s0arr->begin(),s0arr->end(),true));
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> descThisPart=DataArrayInt::New(),descIThisPart=DataArrayInt::New(),revDescThisPart=DataArrayInt::New(),revDescIThisPart=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> thisPartConsti=thisPart->buildDescendingConnectivity(descThisPart,descIThisPart,revDescThisPart,revDescIThisPart);
  const int *revDescThisPartPtr=revDescThisPart->getConstPointer(),*revDescIThisPartPtr=revDescIThisPart->getConstPointer();
  DataArrayInt *idsOtherInConsti=0;
  bool b=thisPartConsti->areCellsIncludedIn(&otherDimM1OnSameCoords,2,idsOtherInConsti);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> idsOtherInConstiAuto(idsOtherInConsti);
  if(!b)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findCellIdsLyingOn : the given mdim-1 mesh in other is not a constituent of this !");
  std::set<int> s1;
  for(const int *idOther=idsOtherInConsti->begin();idOther!=idsOtherInConsti->end();idOther++)
    s1.insert(revDescThisPartPtr+revDescIThisPartPtr[*idOther],revDescThisPartPtr+revDescIThisPartPtr[*idOther+1]);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> s1arr_renum1=DataArrayInt::New(); s1arr_renum1->alloc((int)s1.size(),1); std::copy(s1.begin(),s1.end(),s1arr_renum1->getPointer());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> s1Comparr_renum1=s1arr_renum1->buildComplement(s0arr->getNumberOfTuples());
  DataArrayInt *neighThisPart=0,*neighIThisPart=0;
  ComputeNeighborsOfCellsAdv(descThisPart,descIThisPart,revDescThisPart,revDescIThisPart,neighThisPart,neighIThisPart);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> neighThisPartAuto(neighThisPart),neighIThisPartAuto(neighIThisPart);
  ExtractFromIndexedArrays(s1Comparr_renum1->begin(),s1Comparr_renum1->end(),neighThisPart,neighIThisPart,neighThisPart,neighIThisPart);// reuse of neighThisPart and neighIThisPart
  neighThisPartAuto=neighThisPart; neighIThisPartAuto=neighIThisPart;
  RemoveIdsFromIndexedArrays(s1Comparr_renum1->begin(),s1Comparr_renum1->end(),neighThisPart,neighIThisPart);
  neighThisPartAuto=0;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> s2_tmp=neighIThisPart->deltaShiftIndex();
  const int li[2]={0,1};
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> s2_renum2=s2_tmp->getIdsNotEqualList(li,li+2);
  s2_renum2->transformWithIndArr(s1Comparr_renum1->begin(),s1Comparr_renum1->end());//s2_renum2==s2_renum1
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> s_renum1=DataArrayInt::Aggregate(s2_renum2,s1arr_renum1,0);
  s_renum1->sort();
  //
  s0arr->incrRef(); cellIdsRk0=s0arr;
  s_renum1->incrRef(); cellIdsRk1=s_renum1;
}

/*!
 * This method computes the skin of \b this. That is to say the consituting meshdim-1 mesh is built and only the boundary subpart is
 * returned. This subpart of meshdim-1 mesh is built using meshdim-1 cells in it shared only one cell in \b this.
 * 
 * \return a newly allocated mesh lying on the same coordinates than \b this. The caller has to deal with returned mesh.
 */
MEDCouplingUMesh *MEDCouplingUMesh::computeSkin() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> desc=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> descIndx=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> revDesc=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> revDescIndx=DataArrayInt::New();
  //
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> meshDM1=buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  revDesc=0; desc=0; descIndx=0;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> revDescIndx2=revDescIndx->deltaShiftIndex();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> part=revDescIndx2->getIdsEqual(1);
  return static_cast<MEDCouplingUMesh *>(meshDM1->buildPartOfMySelf(part->begin(),part->end(),true));
}

/*!
 * This methods returns set of nodes in a newly allocated array that the caller has to deal with.
 * The returned nodes ids are those lying on the boundary of \b this.
 */
DataArrayInt *MEDCouplingUMesh::findBoundaryNodes() const
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> skin=computeSkin();
  return skin->computeFetchedNodeIds();
}

MEDCouplingUMesh *MEDCouplingUMesh::buildUnstructured() const throw(INTERP_KERNEL::Exception)
{
  incrRef();
  return const_cast<MEDCouplingUMesh *>(this);
}

/*
 * This method renumber 'this' using 'newNodeNumbers' array of size this->getNumberOfNodes.
 * newNbOfNodes specifies the *std::max_element(newNodeNumbers,newNodeNumbers+this->getNumberOfNodes())
 * This value is asked because often known by the caller of this method.
 * This method, contrary to MEDCouplingMesh::renumberCells does NOT conserve the number of nodes before and after.
 *
 * @param newNodeNumbers array specifying the new numbering in old2New convention.
 * @param newNbOfNodes the new number of nodes.
 */
void MEDCouplingUMesh::renumberNodes(const int *newNodeNumbers, int newNbOfNodes)
{
  MEDCouplingPointSet::renumberNodes(newNodeNumbers,newNbOfNodes);
  renumberNodesInConn(newNodeNumbers);
}

/*
 * This method renumber 'this' using 'newNodeNumbers' array of size this->getNumberOfNodes.
 * newNbOfNodes specifies the *std::max_element(newNodeNumbers,newNodeNumbers+this->getNumberOfNodes())
 * This value is asked because often known by the caller of this method.
 * This method, contrary to MEDCouplingMesh::renumberCells does NOT conserve the number of nodes before and after.
 * The difference with ParaMEDMEM::MEDCouplingUMesh::renumberNodes method is in the fact that the barycenter of merged nodes is computed here.
 *
 * @param newNodeNumbers array specifying the new numbering.
 * @param newNbOfNodes the new number of nodes.
 *
 */
void MEDCouplingUMesh::renumberNodes2(const int *newNodeNumbers, int newNbOfNodes)
{
  MEDCouplingPointSet::renumberNodes2(newNodeNumbers,newNbOfNodes);
  renumberNodesInConn(newNodeNumbers);
}

/*!
 * This method expects that \b this and \b otherDimM1OnSameCoords share the same coordinates array.
 * otherDimM1OnSameCoords->getMeshDimension() is expected to be equal to this->getMeshDimension()-1.
 * This method searches for nodes needed to be duplicated. These nodes are nodes fetched by \b otherDimM1OnSameCoords which are not part of the boundary of \b otherDimM1OnSameCoords.
 * If a node is in the boundary of \b this \b and in the boundary of \b otherDimM1OnSameCoords this node is considerd as needed to be duplicated.
 * When the set of node ids \b nodeIdsToDuplicate is computed, cell ids in \b this is searched so that their connectivity includes at least 1 node in \b nodeIdsToDuplicate.
 *
 * \param [in] otherDimM1OnSameCoords a mesh lying on the same coords than \b this and with a mesh dimension equal to those of \b this minus 1. WARNING this input
 *             parameter is altered during the call.
 * \param [out] nodeIdsToDuplicate node ids needed to be duplicated following the algorithm explain above.
 * \param [out] cellIdsNeededToBeRenum cell ids in \b this in which the renumber of nodes should be performed.
 * \param [out] cellIdsNotModified cell ids int \b this that lies on \b otherDimM1OnSameCoords mesh whose connectivity do \b not need to be modified as it is the case for \b cellIdsNeededToBeRenum.
 *
 * \warning This method modifies param \b otherDimM1OnSameCoords (for speed reasons).
 */
void MEDCouplingUMesh::findNodesToDuplicate(const MEDCouplingUMesh& otherDimM1OnSameCoords, DataArrayInt *& nodeIdsToDuplicate,
                                            DataArrayInt *& cellIdsNeededToBeRenum, DataArrayInt *& cellIdsNotModified) const throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  otherDimM1OnSameCoords.checkFullyDefined();
  if(getCoords()!=otherDimM1OnSameCoords.getCoords())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findNodesToDuplicate : meshes do not share the same coords array !");
  if(otherDimM1OnSameCoords.getMeshDimension()!=getMeshDimension()-1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findNodesToDuplicate : the mesh given in other parameter must have this->getMeshDimension()-1 !");
  DataArrayInt *cellIdsRk0=0,*cellIdsRk1=0;
  findCellIdsLyingOn(otherDimM1OnSameCoords,cellIdsRk0,cellIdsRk1);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cellIdsRk0Auto(cellIdsRk0),cellIdsRk1Auto(cellIdsRk1);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> s0=cellIdsRk1->buildComplement(cellIdsRk0->getNumberOfTuples());
  s0->transformWithIndArr(cellIdsRk0Auto->begin(),cellIdsRk0Auto->end());
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m0Part=static_cast<MEDCouplingUMesh *>(buildPartOfMySelf(s0->begin(),s0->end(),true));
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> s1=m0Part->computeFetchedNodeIds();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> s2=otherDimM1OnSameCoords.computeFetchedNodeIds();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> s3=s2->buildSubstraction(s1);
  cellIdsRk1->transformWithIndArr(cellIdsRk0Auto->begin(),cellIdsRk0Auto->end());
  //
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m0Part2=static_cast<MEDCouplingUMesh *>(buildPartOfMySelf(cellIdsRk1->begin(),cellIdsRk1->end(),true));
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> desc00=DataArrayInt::New(),descI00=DataArrayInt::New(),revDesc00=DataArrayInt::New(),revDescI00=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m01=m0Part2->buildDescendingConnectivity(desc00,descI00,revDesc00,revDescI00);
  DataArrayInt *idsTmp=0;
  bool b=m01->areCellsIncludedIn(&otherDimM1OnSameCoords,2,idsTmp);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ids(idsTmp);
  if(!b)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findNodesToDuplicate : the given mdim-1 mesh in other is not a constituent of this !");
  MEDCouplingUMesh::RemoveIdsFromIndexedArrays(ids->begin(),ids->end(),desc00,descI00);
  DataArrayInt *tmp0=0,*tmp1=0;
  ComputeNeighborsOfCellsAdv(desc00,descI00,revDesc00,revDescI00,tmp0,tmp1);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> neigh00(tmp0);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> neighI00(tmp1);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cellsToModifyConn0_torenum=MEDCouplingUMesh::ComputeSpreadZoneGradually(neigh00,neighI00);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cellsToModifyConn1_torenum=cellsToModifyConn0_torenum->buildComplement(neighI00->getNumberOfTuples()-1);
  cellsToModifyConn0_torenum->transformWithIndArr(cellIdsRk1->begin(),cellIdsRk1->end());
  cellsToModifyConn1_torenum->transformWithIndArr(cellIdsRk1->begin(),cellIdsRk1->end());
  //
  cellIdsNeededToBeRenum=cellsToModifyConn0_torenum; cellsToModifyConn0_torenum->incrRef();
  cellIdsNotModified=cellsToModifyConn1_torenum; cellsToModifyConn1_torenum->incrRef();
  nodeIdsToDuplicate=s3; s3->incrRef();
}

/*!
 * This method operates a modification of the connectivity and coords in \b this.
 * Every time that a node id in [\b nodeIdsToDuplicateBg, \b nodeIdsToDuplicateEnd) will append in nodal connectivity of \b this 
 * its ids will be modified to id this->getNumberOfNodes()+std::distance(nodeIdsToDuplicateBg,std::find(nodeIdsToDuplicateBg,nodeIdsToDuplicateEnd,id)).
 * More explicitely the renumber array in nodes is not explicitely given in old2new to avoid to build a big array of renumbering whereas typically few node ids needs to be
 * renumbered. The node id nodeIdsToDuplicateBg[0] will have id this->getNumberOfNodes()+0, node id nodeIdsToDuplicateBg[1] will have id this->getNumberOfNodes()+1,
 * node id nodeIdsToDuplicateBg[2] will have id this->getNumberOfNodes()+2...
 * 
 * As a consequence nodal connectivity array length will remain unchanged by this method, and nodal connectivity index array will remain unchanged by this method.
 * 
 * \param [in] nodeIdsToDuplicateBg begin of node ids (included) to be duplicated in connectivity only
 * \param [in] nodeIdsToDuplicateEnd end of node ids (excluded) to be duplicated in connectivity only
 */
void MEDCouplingUMesh::duplicateNodes(const int *nodeIdsToDuplicateBg, const int *nodeIdsToDuplicateEnd) throw(INTERP_KERNEL::Exception)
{
  int nbOfNodes=getNumberOfNodes();
  duplicateNodesInCoords(nodeIdsToDuplicateBg,nodeIdsToDuplicateEnd);
  duplicateNodesInConn(nodeIdsToDuplicateBg,nodeIdsToDuplicateEnd,nbOfNodes);
}

/*!
 * This method renumbers nodes \b in \b connectivity \b only \b without \b any \b reference \b to \b coords.
 * This method performs no check on the fact that new coordinate ids are valid. \b Use \b it \b with \b care !
 * This method is an generalization of \ref ParaMEDMEM::MEDCouplingUMesh::shiftNodeNumbersInConn "shiftNodeNumbersInConn method".
 * @param [in] newNodeNumbers in old2New convention
 */
void MEDCouplingUMesh::renumberNodesInConn(const int *newNodeNumbersO2N)
{
  checkConnectivityFullyDefined();
  int *conn=getNodalConnectivity()->getPointer();
  const int *connIndex=getNodalConnectivityIndex()->getConstPointer();
  int nbOfCells=getNumberOfCells();
  for(int i=0;i<nbOfCells;i++)
    for(int iconn=connIndex[i]+1;iconn!=connIndex[i+1];iconn++)
      {
        int& node=conn[iconn];
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
 * This method is an specialization of \ref ParaMEDMEM::MEDCouplingUMesh::renumberNodesInConn "renumberNodesInConn method".
 * 
 * @param [in] delta specifies the shift size applied to nodeId in nodal connectivity in \b this.
 */
void MEDCouplingUMesh::shiftNodeNumbersInConn(int delta) throw(INTERP_KERNEL::Exception)
{
  checkConnectivityFullyDefined();
  int *conn=getNodalConnectivity()->getPointer();
  const int *connIndex=getNodalConnectivityIndex()->getConstPointer();
  int nbOfCells=getNumberOfCells();
  for(int i=0;i<nbOfCells;i++)
    for(int iconn=connIndex[i]+1;iconn!=connIndex[i+1];iconn++)
      {
        int& node=conn[iconn];
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
 * Every time that a node id in [\b nodeIdsToDuplicateBg, \b nodeIdsToDuplicateEnd) will append in nodal connectivity of \b this 
 * its ids will be modified to id offset+std::distance(nodeIdsToDuplicateBg,std::find(nodeIdsToDuplicateBg,nodeIdsToDuplicateEnd,id)).
 * More explicitely the renumber array in nodes is not explicitely given in old2new to avoid to build a big array of renumbering whereas typically few node ids needs to be
 * renumbered. The node id nodeIdsToDuplicateBg[0] will have id offset+0, node id nodeIdsToDuplicateBg[1] will have id offset+1,
 * node id nodeIdsToDuplicateBg[2] will have id offset+2...
 * 
 * As a consequence nodal connectivity array length will remain unchanged by this method, and nodal connectivity index array will remain unchanged by this method.
 * As an another consequense after the call of this method \b this can be transiently non cohrent.
 * 
 * \param [in] nodeIdsToDuplicateBg begin of node ids (included) to be duplicated in connectivity only
 * \param [in] nodeIdsToDuplicateEnd end of node ids (excluded) to be duplicated in connectivity only
 * \param [in] offset the offset applied to all node ids in connectivity that are in [nodeIdsToDuplicateBg,nodeIdsToDuplicateEnd). 
 */
void MEDCouplingUMesh::duplicateNodesInConn(const int *nodeIdsToDuplicateBg, const int *nodeIdsToDuplicateEnd, int offset) throw(INTERP_KERNEL::Exception)
{
  checkConnectivityFullyDefined();
  std::map<int,int> m;
  int val=offset;
  for(const int *work=nodeIdsToDuplicateBg;work!=nodeIdsToDuplicateEnd;work++,val++)
    m[*work]=val;
  int *conn=getNodalConnectivity()->getPointer();
  const int *connIndex=getNodalConnectivityIndex()->getConstPointer();
  int nbOfCells=getNumberOfCells();
  for(int i=0;i<nbOfCells;i++)
    for(int iconn=connIndex[i]+1;iconn!=connIndex[i+1];iconn++)
      {
        int& node=conn[iconn];
        if(node>=0)//avoid polyhedron separator
          {
            std::map<int,int>::iterator it=m.find(node);
            if(it!=m.end())
              node=(*it).second;
          }
      }
  updateTime();
}

/*!
 * This method renumbers cells of 'this' using the array specified by [old2NewBg;old2NewBg+getNumberOfCells())
 *
 * Contrary to MEDCouplingPointSet::renumberNodes, this method makes a permutation without any fuse of cell.
 * After the call of this method the number of cells remains the same as before.
 *
 * If 'check' equals true the method will check that any elements in [old2NewBg;old2NewEnd) is unique ; if not
 * an INTERP_KERNEL::Exception will be thrown. When 'check' equals true [old2NewBg;old2NewEnd) is not expected to
 * be strictly in [0;this->getNumberOfCells()).
 *
 * If 'check' equals false the method will not check the content of [old2NewBg;old2NewEnd).
 * To avoid any throw of SIGSEGV when 'check' equals false, the elements in [old2NewBg;old2NewEnd) should be unique and
 * should be contained in[0;this->getNumberOfCells()).
 * 
 * \param [in] old2NewBg is expected to be a dynamically allocated pointer of size at least equal to this->getNumberOfCells()
 */
void MEDCouplingUMesh::renumberCells(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
  checkConnectivityFullyDefined();
  int nbCells=getNumberOfCells();
  const int *array=old2NewBg;
  if(check)
    array=DataArrayInt::CheckAndPreparePermutation(old2NewBg,old2NewBg+nbCells);
  //
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newConn=DataArrayInt::New();
  newConn->alloc(_nodal_connec->getNumberOfTuples(),_nodal_connec->getNumberOfComponents());
  newConn->copyStringInfoFrom(*_nodal_connec);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newConnI=DataArrayInt::New();
  newConnI->alloc(_nodal_connec_index->getNumberOfTuples(),_nodal_connec_index->getNumberOfComponents());
  newConnI->copyStringInfoFrom(*_nodal_connec_index);
  //
  int *newC=newConn->getPointer();
  int *newCI=newConnI->getPointer();
  int loc=0;
  newCI[0]=loc;
  for(int i=0;i<nbCells;i++)
    {
      std::size_t pos=std::distance(array,std::find(array,array+nbCells,i));
      int nbOfElts=connI[pos+1]-connI[pos];
      newC=std::copy(conn+connI[pos],conn+connI[pos+1],newC);
      loc+=nbOfElts;
      newCI[i+1]=loc;
    }
  //
  setConnectivity(newConn,newConnI);
  if(check)
    delete [] const_cast<int *>(array);
}

/*!
 * Given a boundary box 'bbox' returns elements 'elems' contained in this 'bbox'.
 * Warning 'elems' is incremented during the call so if elems is not empty before call returned elements will be
 * added in 'elems' parameter.
 */
void MEDCouplingUMesh::getCellsInBoundingBox(const double *bbox, double eps, std::vector<int>& elems) const
{
  if(getMeshDimension()==-1)
    {
      elems.push_back(0);
      return;
    }
  int dim=getSpaceDimension();
  double* elem_bb=new double[2*dim];
  const int* conn      = getNodalConnectivity()->getConstPointer();
  const int* conn_index= getNodalConnectivityIndex()->getConstPointer();
  const double* coords = getCoords()->getConstPointer();
  int nbOfCells=getNumberOfCells();
  for ( int ielem=0; ielem<nbOfCells;ielem++ )
    {
      for (int i=0; i<dim; i++)
        {
          elem_bb[i*2]=std::numeric_limits<double>::max();
          elem_bb[i*2+1]=-std::numeric_limits<double>::max();
        }

      for (int inode=conn_index[ielem]+1; inode<conn_index[ielem+1]; inode++)//+1 due to offset of cell type.
        {
          int node= conn[inode];
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
        {
          elems.push_back(ielem);
        }
    }
  delete [] elem_bb;
}

/*!
 * Given a boundary box 'bbox' returns elements 'elems' contained in this 'bbox' or touching 'bbox' (within 'eps' distance).
 * Warning 'elems' is incremented during the call so if elems is not empty before call returned elements will be
 * added in 'elems' parameter.
 */
void MEDCouplingUMesh::getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps, std::vector<int>& elems)
{
  if(getMeshDimension()==-1)
    {
      elems.push_back(0);
      return;
    }
  int dim=getSpaceDimension();
  double* elem_bb=new double[2*dim];
  const int* conn      = getNodalConnectivity()->getConstPointer();
  const int* conn_index= getNodalConnectivityIndex()->getConstPointer();
  const double* coords = getCoords()->getConstPointer();
  int nbOfCells=getNumberOfCells();
  for ( int ielem=0; ielem<nbOfCells;ielem++ )
    {
      for (int i=0; i<dim; i++)
        {
          elem_bb[i*2]=std::numeric_limits<double>::max();
          elem_bb[i*2+1]=-std::numeric_limits<double>::max();
        }

      for (int inode=conn_index[ielem]+1; inode<conn_index[ielem+1]; inode++)//+1 due to offset of cell type.
        {
          int node= conn[inode];
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
      if (intersectsBoundingBox(bbox, elem_bb, dim, eps))
        {
          elems.push_back(ielem);
        }
    }
  delete [] elem_bb;
}

/*!
 * Returns the cell type of cell with id 'cellId'.
 */
INTERP_KERNEL::NormalizedCellType MEDCouplingUMesh::getTypeOfCell(int cellId) const
{
  const int *ptI=_nodal_connec_index->getConstPointer();
  const int *pt=_nodal_connec->getConstPointer();
  return (INTERP_KERNEL::NormalizedCellType) pt[ptI[cellId]];
}

/*!
 * This method returns a newly allocated array containing cell ids (ascendingly sorted) whose geometric type are equal to type.
 * This method throws an INTERP_KERNEL::Exception if meshdimension of \b this is not equal to those of \b type.
 * The coordinates array is not considered here.
 *
 * \param [in] type the geometric type
 * \return the 
 */
DataArrayInt *MEDCouplingUMesh::giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception)
{
  
  std::vector<int> v;
  checkConnectivityFullyDefined();
  int nbCells=getNumberOfCells();
  int mdim=getMeshDimension();
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  if(mdim!=(int)cm.getDimension())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::giveCellsWithType : Mismatch between mesh dimension and dimension of the cell !");
  const int *ptI=_nodal_connec_index->getConstPointer();
  const int *pt=_nodal_connec->getConstPointer();
  for(int i=0;i<nbCells;i++)
    {
      if((INTERP_KERNEL::NormalizedCellType)pt[ptI[i]]==type)
        v.push_back(i);
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New(); ret->alloc((int)v.size(),1);
  std::copy(v.begin(),v.end(),ret->getPointer());
  ret->incrRef();
  return ret;
}

/*!
 * Returns nb of cells having the geometric type 'type'.
 */
int MEDCouplingUMesh::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const
{
  const int *ptI=_nodal_connec_index->getConstPointer();
  const int *pt=_nodal_connec->getConstPointer();
  int nbOfCells=getNumberOfCells();
  int ret=0;
  for(int i=0;i<nbOfCells;i++)
    if((INTERP_KERNEL::NormalizedCellType) pt[ptI[i]]==type)
      ret++;
  return ret;
}

/*!
 * Appends the nodal connectivity in 'conn' of cell with id 'cellId'.
 * All elements added in conn can be used by MEDCouplingUMesh::getCoordinatesOfNode method.
 * That is to say -1 separator is omitted in returned conn.
 */
void MEDCouplingUMesh::getNodeIdsOfCell(int cellId, std::vector<int>& conn) const
{
  const int *ptI=_nodal_connec_index->getConstPointer();
  const int *pt=_nodal_connec->getConstPointer();
  for(const int *w=pt+ptI[cellId]+1;w!=pt+ptI[cellId+1];w++)
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
  ret << "Mesh dimension : " << _mesh_dim << "\nSpace dimension : ";
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

std::string MEDCouplingUMesh::reprConnectivityOfThis() const
{
  std::ostringstream ret;
  reprConnectivityOfThisLL(ret);
  return ret.str();
}

/*!
 * This method builds a newly allocated instance (with the same name than 'this') that the caller has the responsability to deal with.
 * This method returns an instance with all arrays allocated (connectivity, connectivity index, coordinates)
 * but with length of these arrays set to 0. It allows to define an "empty" mesh (with nor cells nor nodes but compliant with
 * some algos).
 * 
 * This method expects that 'this' has a mesh dimension set and higher or equal to 0. If not an exception will be thrown.
 * This method analyzes the 3 arrays of 'this'. For each the following behaviour is done : if the array is null a newly one is created
 * with number of tuples set to 0, if not the array is taken as this in the returned instance.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildSetInstanceFromThis(int spaceDim) const throw(INTERP_KERNEL::Exception)
{
  int mdim=getMeshDimension();
  if(mdim<0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSetInstanceFromThis : invalid mesh dimension ! Should be >= 0 !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> ret=MEDCouplingUMesh::New(getName(),mdim);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp1,tmp2;
  bool needToCpyCT=true;
  if(!_nodal_connec)
    {
      tmp1=DataArrayInt::New(); tmp1->alloc(0,1);
      needToCpyCT=false;
    }
  else
    {
      tmp1=_nodal_connec;
      tmp1->incrRef();
    }
  if(!_nodal_connec_index)
    {
      tmp2=DataArrayInt::New(); tmp2->alloc(1,1); tmp2->setIJ(0,0,0);
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
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> coords=DataArrayDouble::New(); coords->alloc(0,spaceDim);
      ret->setCoords(coords);
    }
  else
    ret->setCoords(_coords);
  ret->incrRef();
  return ret;
}

void MEDCouplingUMesh::reprConnectivityOfThisLL(std::ostringstream& stream) const
{
  if(_nodal_connec!=0 && _nodal_connec_index!=0)
    {
      int nbOfCells=getNumberOfCells();
      const int *c=_nodal_connec->getConstPointer();
      const int *ci=_nodal_connec_index->getConstPointer();
      for(int i=0;i<nbOfCells;i++)
        {
          const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)c[ci[i]]);
          stream << "Cell #" << i << " " << cm.getRepr() << " : ";
          std::copy(c+ci[i]+1,c+ci[i+1],std::ostream_iterator<int>(stream," "));
          stream << "\n";
        }
    }
  else
    stream << "Connectivity not defined !\n";
}

int MEDCouplingUMesh::getNumberOfNodesInCell(int cellId) const
{
  const int *ptI=_nodal_connec_index->getConstPointer();
  const int *pt=_nodal_connec->getConstPointer();
  if(pt[ptI[cellId]]!=INTERP_KERNEL::NORM_POLYHED)
    return ptI[cellId+1]-ptI[cellId]-1;
  else
    return (int)std::count_if(pt+ptI[cellId]+1,pt+ptI[cellId+1],std::bind2nd(std::not_equal_to<int>(),-1));
}

/*!
 * This method is equivalent to MEDCouplingUMesh::getAllTypes excecpt that it returns only types of submesh which cell ids are in [begin,end).
 * This method avoids to compute explicitely submesh to get its types.
 */
std::set<INTERP_KERNEL::NormalizedCellType> MEDCouplingUMesh::getTypesOfPart(const int *begin, const int *end) const throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  std::set<INTERP_KERNEL::NormalizedCellType> ret;
  const int *conn=_nodal_connec->getConstPointer();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  for(const int *w=begin;w!=end;w++)
    ret.insert((INTERP_KERNEL::NormalizedCellType)conn[connIndex[*w]]);
  return ret;
}

/*!
 * Method reserved for advanced users having prepared their connectivity before.
 * Arrays 'conn' and 'connIndex' will be aggregated without any copy and their counter will be incremented.
 */
void MEDCouplingUMesh::setConnectivity(DataArrayInt *conn, DataArrayInt *connIndex, bool isComputingTypes)
{
  DataArrayInt::SetArrayIn(conn,_nodal_connec);
  DataArrayInt::SetArrayIn(connIndex,_nodal_connec_index);
  if(isComputingTypes)
    computeTypes();
  declareAsNew();
}

/*!
 * Copy constructor. If 'deepCpy' is false 'this' is a shallow copy of other.
 * If 'deeCpy' is true all arrays (coordinates and connectivities) are deeply copied.
 */
MEDCouplingUMesh::MEDCouplingUMesh(const MEDCouplingUMesh& other, bool deepCopy):MEDCouplingPointSet(other,deepCopy),_iterator(-1),_mesh_dim(other._mesh_dim),
                                                                                 _nodal_connec(0),_nodal_connec_index(0),
                                                                                _types(other._types)
{
  if(other._nodal_connec)
    _nodal_connec=other._nodal_connec->performCpy(deepCopy);
  if(other._nodal_connec_index)
    _nodal_connec_index=other._nodal_connec_index->performCpy(deepCopy);
}

MEDCouplingUMesh::~MEDCouplingUMesh()
{
  if(_nodal_connec)
    _nodal_connec->decrRef();
  if(_nodal_connec_index)
    _nodal_connec_index->decrRef();
}

/*!
 * This method recomputes all cell types of 'this'.
 */
void MEDCouplingUMesh::computeTypes()
{
  if(_nodal_connec && _nodal_connec_index)
    {
      _types.clear();
      const int *conn=_nodal_connec->getConstPointer();
      const int *connIndex=_nodal_connec_index->getConstPointer();
      int nbOfElem=_nodal_connec_index->getNbOfElems()-1;
      for(const int *pt=connIndex;pt!=connIndex+nbOfElem;pt++)
        _types.insert((INTERP_KERNEL::NormalizedCellType)conn[*pt]);
    }
}

/*!
 * This method checks that all arrays are set. If yes nothing done if no an exception is thrown.
 */
void MEDCouplingUMesh::checkFullyDefined() const throw(INTERP_KERNEL::Exception)
{
  if(!_nodal_connec_index || !_nodal_connec || !_coords)
    throw INTERP_KERNEL::Exception("Reverse nodal connectivity computation requires full connectivity and coordinates set in unstructured mesh.");
}

/*!
 * This method checks that all connectivity arrays are set. If yes nothing done if no an exception is thrown.
 */
void MEDCouplingUMesh::checkConnectivityFullyDefined() const throw(INTERP_KERNEL::Exception)
{
  if(!_nodal_connec_index || !_nodal_connec)
    throw INTERP_KERNEL::Exception("Reverse nodal connectivity computation requires full connectivity set in unstructured mesh.");
}

int MEDCouplingUMesh::getNumberOfCells() const
{ 
  if(_nodal_connec_index)
    if(_iterator==-1)
      return _nodal_connec_index->getNumberOfTuples()-1;
    else
      return _iterator;
  else
    if(_mesh_dim==-1)
      return 1;
    else
      throw INTERP_KERNEL::Exception("Unable to get number of cells because no connectivity specified !");
}

int MEDCouplingUMesh::getMeshDimension() const
{
  if(_mesh_dim<-1)
    throw INTERP_KERNEL::Exception("No mesh dimension specified !");
  return _mesh_dim;
}

/*!
 * This method is for test reason. Normally the integer returned is not useable by user.
 */
int MEDCouplingUMesh::getMeshLength() const
{
  return _nodal_connec->getNbOfElems();
}

/*!
 * First step of serialization process. Used by ParaMEDMEM and MEDCouplingCorba to transfert data between process.
 */
void MEDCouplingUMesh::getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const
{
  MEDCouplingPointSet::getTinySerializationInformation(tinyInfoD,tinyInfo,littleStrings);
  tinyInfo.push_back(getMeshDimension());
  tinyInfo.push_back(getNumberOfCells());
  if(_nodal_connec)
    tinyInfo.push_back(getMeshLength());
  else
    tinyInfo.push_back(-1);
}

/*!
 * First step of unserialization process.
 */
bool MEDCouplingUMesh::isEmptyMesh(const std::vector<int>& tinyInfo) const
{
  return tinyInfo[6]<=0;
}

/*!
 * Second step of serialization process.
 * @param tinyInfo must be equal to the result given by getTinySerializationInformation method.
 */
void MEDCouplingUMesh::resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const
{
  MEDCouplingPointSet::resizeForUnserialization(tinyInfo,a1,a2,littleStrings);
  if(tinyInfo[5]!=-1)
    a1->alloc(tinyInfo[7]+tinyInfo[6]+1,1);
}

/*!
 * Third and final step of serialization process.
 */
void MEDCouplingUMesh::serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const
{
  MEDCouplingPointSet::serialize(a1,a2);
  if(getMeshDimension()>-1)
    {
      a1=DataArrayInt::New();
      a1->alloc(getMeshLength()+getNumberOfCells()+1,1);
      int *ptA1=a1->getPointer();
      const int *conn=getNodalConnectivity()->getConstPointer();
      const int *index=getNodalConnectivityIndex()->getConstPointer();
      ptA1=std::copy(index,index+getNumberOfCells()+1,ptA1);
      std::copy(conn,conn+getMeshLength(),ptA1);
    }
  else
    a1=0;
}

/*!
 * Second and final unserialization process.
 * @param tinyInfo must be equal to the result given by getTinySerializationInformation method.
 */
void MEDCouplingUMesh::unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2, const std::vector<std::string>& littleStrings)
{
  MEDCouplingPointSet::unserialization(tinyInfoD,tinyInfo,a1,a2,littleStrings);
  setMeshDimension(tinyInfo[5]);
  if(tinyInfo[7]!=-1)
    {
      // Connectivity
      const int *recvBuffer=a1->getConstPointer();
      DataArrayInt* myConnecIndex=DataArrayInt::New();
      myConnecIndex->alloc(tinyInfo[6]+1,1);
      std::copy(recvBuffer,recvBuffer+tinyInfo[6]+1,myConnecIndex->getPointer());
      DataArrayInt* myConnec=DataArrayInt::New();
      myConnec->alloc(tinyInfo[7],1);
      std::copy(recvBuffer+tinyInfo[6]+1,recvBuffer+tinyInfo[6]+1+tinyInfo[7],myConnec->getPointer());
      setConnectivity(myConnec, myConnecIndex) ;
      myConnec->decrRef();
      myConnecIndex->decrRef();
    }
}

/*!
 * This is the low algorithm of MEDCouplingUMesh::buildPartOfMySelf2.
 * CellIds are given using range specified by a start an end and step.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildPartOfMySelfKeepCoords2(int start, int end, int step) const
{
  checkFullyDefined();
  int ncell=getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> ret=MEDCouplingUMesh::New();
  ret->_mesh_dim=_mesh_dim;
  ret->setCoords(_coords);
  int newNbOfCells=DataArray::GetNumberOfItemGivenBESRelative(start,end,step,"MEDCouplingUMesh::buildPartOfMySelfKeepCoords2 : ");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newConnI=DataArrayInt::New(); newConnI->alloc(newNbOfCells+1,1);
  int *newConnIPtr=newConnI->getPointer(); *newConnIPtr=0;
  int work=start;
  const int *conn=_nodal_connec->getConstPointer();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  for(int i=0;i<newNbOfCells;i++,newConnIPtr++,work+=step)
    {
      if(work>=0 && work<ncell)
        {
          newConnIPtr[1]=newConnIPtr[0]+connIndex[work+1]-connIndex[work];
        }
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::buildPartOfMySelfKeepCoords2 : On pos #" << i << " input cell id =" << work << " should be in [0," << ncell << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newConn=DataArrayInt::New(); newConn->alloc(newConnIPtr[0],1);
  int *newConnPtr=newConn->getPointer();
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  work=start;
  for(int i=0;i<newNbOfCells;i++,newConnIPtr++,work+=step)
    {
      types.insert((INTERP_KERNEL::NormalizedCellType)conn[connIndex[work]]);
      newConnPtr=std::copy(conn+connIndex[work],conn+connIndex[work+1],newConnPtr);
    }
  ret->setConnectivity(newConn,newConnI,false);
  ret->_types=types;
  ret->copyTinyInfoFrom(this);
  std::string name(getName());
  std::size_t sz=strlen(PART_OF_NAME);
  if(name.length()>=sz)
    name=name.substr(0,sz);
  if(name!=PART_OF_NAME)
    {
      std::ostringstream stream; stream << PART_OF_NAME << getName();
      ret->setName(stream.str().c_str());
    }
  else
    ret->setName(getName());
  ret->incrRef();
  return ret;
}

/*!
 * This is the low algorithm of MEDCouplingUMesh::buildPartOfMySelf.
 * Keeps from 'this' only cells which constituing point id are in the ids specified by ['begin','end').
 * The return newly allocated mesh will share the same coordinates as 'this'.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildPartOfMySelfKeepCoords(const int *begin, const int *end) const
{
  checkFullyDefined();
  int ncell=getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> ret=MEDCouplingUMesh::New();
  ret->_mesh_dim=_mesh_dim;
  ret->setCoords(_coords);
  std::size_t nbOfElemsRet=std::distance(begin,end);
  int *connIndexRet=new int[nbOfElemsRet+1];
  connIndexRet[0]=0;
  const int *conn=_nodal_connec->getConstPointer();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  int newNbring=0;
  for(const int *work=begin;work!=end;work++,newNbring++)
    {
      if(*work>=0 && *work<ncell)
        connIndexRet[newNbring+1]=connIndexRet[newNbring]+connIndex[*work+1]-connIndex[*work];
      else
        {
          delete [] connIndexRet;
          std::ostringstream oss; oss << "MEDCouplingUMesh::buildPartOfMySelfKeepCoords : On pos #" << std::distance(begin,work) << " input cell id =" << *work << " should be in [0," << ncell << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  int *connRet=new int[connIndexRet[nbOfElemsRet]];
  int *connRetWork=connRet;
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  for(const int *work=begin;work!=end;work++)
    {
      types.insert((INTERP_KERNEL::NormalizedCellType)conn[connIndex[*work]]);
      connRetWork=std::copy(conn+connIndex[*work],conn+connIndex[*work+1],connRetWork);
    }
  DataArrayInt *connRetArr=DataArrayInt::New();
  connRetArr->useArray(connRet,true,CPP_DEALLOC,connIndexRet[nbOfElemsRet],1);
  DataArrayInt *connIndexRetArr=DataArrayInt::New();
  connIndexRetArr->useArray(connIndexRet,true,CPP_DEALLOC,(int)nbOfElemsRet+1,1);
  ret->setConnectivity(connRetArr,connIndexRetArr,false);
  ret->_types=types;
  connRetArr->decrRef();
  connIndexRetArr->decrRef();
  ret->copyTinyInfoFrom(this);
  std::string name(getName());
  std::size_t sz=strlen(PART_OF_NAME);
  if(name.length()>=sz)
    name=name.substr(0,sz);
  if(name!=PART_OF_NAME)
    {
      std::ostringstream stream; stream << PART_OF_NAME << getName();
      ret->setName(stream.str().c_str());
    }
  else
    ret->setName(getName());
  ret->incrRef();
  return ret;
}

/*!
 * brief returns the volumes of the cells underlying the field \a field
 *
 * For 2D geometries, the returned field contains the areas.
 * For 3D geometries, the returned field contains the volumes.
 *
 * param field field on which cells the volumes are required
 * return field containing the volumes, area or length depending the meshdimension.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::getMeasureField(bool isAbs) const
{
  std::string name="MeasureOfMesh_";
  name+=getName();
  int nbelem=getNumberOfCells();
  MEDCouplingFieldDouble *field=MEDCouplingFieldDouble::New(ON_CELLS);
  field->setName(name.c_str());
  DataArrayDouble* array=DataArrayDouble::New();
  array->alloc(nbelem,1);
  double *area_vol=array->getPointer();
  field->setArray(array) ;
  array->decrRef();
  field->setMesh(const_cast<MEDCouplingUMesh *>(this));
  if(getMeshDimension()!=-1)
    {
      int ipt;
      INTERP_KERNEL::NormalizedCellType type;
      int dim_space=getSpaceDimension();
      const double *coords=getCoords()->getConstPointer();
      const int *connec=getNodalConnectivity()->getConstPointer();
      const int *connec_index=getNodalConnectivityIndex()->getConstPointer();
      for(int iel=0;iel<nbelem;iel++)
        {
          ipt=connec_index[iel];
          type=(INTERP_KERNEL::NormalizedCellType)connec[ipt];
          area_vol[iel]=INTERP_KERNEL::computeVolSurfOfCell2<int,INTERP_KERNEL::ALL_C_MODE>(type,connec+ipt+1,connec_index[iel+1]-ipt-1,coords,dim_space);
        }
      if(isAbs)
        std::transform(area_vol,area_vol+nbelem,area_vol,std::ptr_fun<double,double>(fabs));
    }
  else
    {
      area_vol[0]=std::numeric_limits<double>::max();
    }
  return field;
}

/*!
 * This method is equivalent to MEDCouplingUMesh::getMeasureField except that only part defined by [begin,end) is returned !
 * This method avoids to build explicitely part of this to perform the work.
 */
DataArrayDouble *MEDCouplingUMesh::getPartMeasureField(bool isAbs, const int *begin, const int *end) const
{
  std::string name="PartMeasureOfMesh_";
  name+=getName();
  int nbelem=(int)std::distance(begin,end);
  DataArrayDouble* array=DataArrayDouble::New();
  array->setName(name.c_str());
  array->alloc(nbelem,1);
  double *area_vol=array->getPointer();
  if(getMeshDimension()!=-1)
    {
      int ipt;
      INTERP_KERNEL::NormalizedCellType type;
      int dim_space=getSpaceDimension();
      const double *coords=getCoords()->getConstPointer();
      const int *connec=getNodalConnectivity()->getConstPointer();
      const int *connec_index=getNodalConnectivityIndex()->getConstPointer();
      for(const int *iel=begin;iel!=end;iel++)
        {
          ipt=connec_index[*iel];
          type=(INTERP_KERNEL::NormalizedCellType)connec[ipt];
          *area_vol++=INTERP_KERNEL::computeVolSurfOfCell2<int,INTERP_KERNEL::ALL_C_MODE>(type,connec+ipt+1,connec_index[*iel+1]-ipt-1,coords,dim_space);
        }
      if(isAbs)
        std::transform(array->getPointer(),area_vol,array->getPointer(),std::ptr_fun<double,double>(fabs));
    }
  else
    {
      area_vol[0]=std::numeric_limits<double>::max();
    }
  return array;
}

/*!
 * This methods returns a field on nodes and no time. This method is usefull to check "P1*" conservative interpolators.
 * This field returns the getMeasureField of the dualMesh in P1 sens of 'this'.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::getMeasureFieldOnNode(bool isAbs) const
{
  MEDCouplingFieldDouble *tmp=getMeasureField(isAbs);
  std::string name="MeasureOnNodeOfMesh_";
  name+=getName();
  int nbNodes=getNumberOfNodes();
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(ON_NODES);
  double cst=1./((double)getMeshDimension()+1.);
  DataArrayDouble* array=DataArrayDouble::New();
  array->alloc(nbNodes,1);
  double *valsToFill=array->getPointer();
  std::fill(valsToFill,valsToFill+nbNodes,0.);
  const double *values=tmp->getArray()->getConstPointer();
  DataArrayInt *da=DataArrayInt::New();
  DataArrayInt *daInd=DataArrayInt::New();
  getReverseNodalConnectivity(da,daInd);
  const int *daPtr=da->getConstPointer();
  const int *daIPtr=daInd->getConstPointer();
  for(int i=0;i<nbNodes;i++)
    for(const int *cell=daPtr+daIPtr[i];cell!=daPtr+daIPtr[i+1];cell++)
      valsToFill[i]+=cst*values[*cell];
  ret->setMesh(this);
  da->decrRef();
  daInd->decrRef();
  ret->setArray(array);
  array->decrRef();
  tmp->decrRef();
  return ret;
}

/*!
 * This methods returns a vector field on cells that represents the orthogonal vector normalized of each 2D cell of this.
 * This method is only callable on mesh with meshdim == 2 and spacedim==2 or 3.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::buildOrthogonalField() const
{
  if((getMeshDimension()!=2) && (getMeshDimension()!=1 || getSpaceDimension()!=2))
    throw INTERP_KERNEL::Exception("Expected a umesh with ( meshDim == 2 spaceDim == 2 or 3 ) or ( meshDim == 1 spaceDim == 2 ) !");
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  DataArrayDouble *array=DataArrayDouble::New();
  int nbOfCells=getNumberOfCells();
  int nbComp=getMeshDimension()+1;
  array->alloc(nbOfCells,nbComp);
  double *vals=array->getPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const int *conn=_nodal_connec->getConstPointer();
  const double *coords=_coords->getConstPointer();
  if(getMeshDimension()==2)
    {
      if(getSpaceDimension()==3)
        {
          DataArrayDouble *loc=getBarycenterAndOwner();
          const double *locPtr=loc->getConstPointer();
          for(int i=0;i<nbOfCells;i++,vals+=3)
            {
              int offset=connI[i];
              INTERP_KERNEL::crossprod<3>(locPtr+3*i,coords+3*conn[offset+1],coords+3*conn[offset+2],vals);
              double n=INTERP_KERNEL::norm<3>(vals);
              std::transform(vals,vals+3,vals,std::bind2nd(std::multiplies<double>(),1./n));
            }
          loc->decrRef();
        }
      else
        {
          for(int i=0;i<nbOfCells;i++)
            { vals[3*i]=0.; vals[3*i+1]=0.; vals[3*i+2]=1.; }
        }
    }
  else//meshdimension==1
    {
      double tmp[2];
      for(int i=0;i<nbOfCells;i++)
        {
          int offset=connI[i];
          std::transform(coords+2*conn[offset+2],coords+2*conn[offset+2]+2,coords+2*conn[offset+1],tmp,std::minus<double>());
          double n=INTERP_KERNEL::norm<2>(tmp);
          std::transform(tmp,tmp+2,tmp,std::bind2nd(std::multiplies<double>(),1./n));
          *vals++=-tmp[1];
          *vals++=tmp[0];
        }
    }
  ret->setArray(array);
  array->decrRef();
  ret->setMesh(this);
  return ret;
}

/*!
 * This method is equivalent to MEDCouplingUMesh::buildOrthogonalField except that only part defined by [begin,end) is returned !
 * This method avoids to build explicitely part of this to perform the work.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::buildPartOrthogonalField(const int *begin, const int *end) const
{
  if((getMeshDimension()!=2) && (getMeshDimension()!=1 || getSpaceDimension()!=2))
    throw INTERP_KERNEL::Exception("Expected a umesh with ( meshDim == 2 spaceDim == 2 or 3 ) or ( meshDim == 1 spaceDim == 2 ) !");
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  DataArrayDouble *array=DataArrayDouble::New();
  std::size_t nbelems=std::distance(begin,end);
  int nbComp=getMeshDimension()+1;
  array->alloc((int)nbelems,nbComp);
  double *vals=array->getPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const int *conn=_nodal_connec->getConstPointer();
  const double *coords=_coords->getConstPointer();
  if(getMeshDimension()==2)
    {
      if(getSpaceDimension()==3)
        {
          DataArrayDouble *loc=getPartBarycenterAndOwner(begin,end);
          const double *locPtr=loc->getConstPointer();
          for(const int *i=begin;i!=end;i++,vals+=3,locPtr+=3)
            {
              int offset=connI[*i];
              INTERP_KERNEL::crossprod<3>(locPtr,coords+3*conn[offset+1],coords+3*conn[offset+2],vals);
              double n=INTERP_KERNEL::norm<3>(vals);
              std::transform(vals,vals+3,vals,std::bind2nd(std::multiplies<double>(),1./n));
            }
          loc->decrRef();
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
      for(const int *i=begin;i!=end;i++)
        {
          int offset=connI[*i];
          std::transform(coords+2*conn[offset+2],coords+2*conn[offset+2]+2,coords+2*conn[offset+1],tmp,std::minus<double>());
          double n=INTERP_KERNEL::norm<2>(tmp);
          std::transform(tmp,tmp+2,tmp,std::bind2nd(std::multiplies<double>(),1./n));
          *vals++=-tmp[1];
          *vals++=tmp[0];
        }
    }
  ret->setArray(array);
  array->decrRef();
  ret->setMesh(this);
  return ret;
}

/*!
 * This methods returns a vector newly created field on cells that represents the direction vector of each 1D cell of this.
 * This method is only callable on mesh with meshdim == 1 containing only SEG2.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::buildDirectionVectorField() const
{
   if(getMeshDimension()!=1)
    throw INTERP_KERNEL::Exception("Expected a umesh with meshDim == 1 for buildDirectionVectorField !");
   if(_types.size()!=1 || *(_types.begin())!=INTERP_KERNEL::NORM_SEG2)
     throw INTERP_KERNEL::Exception("Expected a umesh with only NORM_SEG2 type of elements for buildDirectionVectorField !");
   MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
   DataArrayDouble *array=DataArrayDouble::New();
   int nbOfCells=getNumberOfCells();
   int spaceDim=getSpaceDimension();
   array->alloc(nbOfCells,spaceDim);
   double *pt=array->getPointer();
   const double *coo=getCoords()->getConstPointer();
   std::vector<int> conn;
   conn.reserve(2);
   for(int i=0;i<nbOfCells;i++)
     {
       conn.resize(0);
       getNodeIdsOfCell(i,conn);
       pt=std::transform(coo+conn[1]*spaceDim,coo+(conn[1]+1)*spaceDim,coo+conn[0]*spaceDim,pt,std::minus<double>());
     }
   ret->setArray(array);
   array->decrRef();
   ret->setMesh(this);
   return ret;   
}

/*!
 * This method expects that 'this' is fully defined and has a spaceDim==3 and a meshDim==3. If it is not the case an exception will be thrown.
 * This method returns 2 objects : 
 * - a newly created mesh instance containing the result of the slice lying on different coords than 'this' and with a meshdim == 2
 * - a newly created dataarray having number of tuples equal to the number of cells in returned mesh that tells for each 2D cell in returned
 *   mesh the 3D cell id is 'this' it comes from.
 * This method works only for linear meshes (non quadratic).
 * If plane crosses within 'eps' a face in 'this' shared by more than 1 cell, 2 output faces will be generated. The 2 faces having the same geometry than intersecting
 * face. Only 'cellIds' parameter can distinguish the 2.
 * @param origin is the origin of the plane. It should be an array of length 3.
 * @param vec is the direction vector of the plane. It should be an array of length 3. Norm of 'vec' should be > 1e-6.
 * @param eps is the precision. It is used by called method MEDCouplingUMesh::getCellIdsCrossingPlane for the first 3D cell selection (in absolute). 'eps' is
 * also used to state if new points should be created or already existing points are reused. 'eps' is also used to tells if plane overlaps a face, edge or nodes (in absolute).
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildSlice3D(const double *origin, const double *vec, double eps, DataArrayInt *&cellIds) const throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSlice3D works on umeshes with meshdim equal to 3 and spaceDim equal to 3 too!");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> candidates=getCellIdsCrossingPlane(origin,vec,eps);
  if(candidates->empty())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSlice3D : No 3D cells in this intercepts the specified plane considering bounding boxes !");
  std::vector<int> nodes;
  std::vector<int> cellIds2D,cellIds1D;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> subMesh=static_cast<MEDCouplingUMesh*>(buildPartOfMySelf(candidates->begin(),candidates->end(),false));
  subMesh->findNodesOnPlane(origin,vec,eps,nodes);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> desc1=DataArrayInt::New(),desc2=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> descIndx1=DataArrayInt::New(),descIndx2=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> revDesc1=DataArrayInt::New(),revDesc2=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> revDescIndx1=DataArrayInt::New(),revDescIndx2=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mDesc2=subMesh->buildDescendingConnectivity(desc2,descIndx2,revDesc2,revDescIndx2);//meshDim==2 spaceDim==3
  revDesc2=0; revDescIndx2=0;
  mDesc2->fillCellIdsToKeepFromNodeIds(&nodes[0],&nodes[0]+nodes.size(),true,cellIds2D);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mDesc1=mDesc2->buildDescendingConnectivity(desc1,descIndx1,revDesc1,revDescIndx1);//meshDim==1 spaceDim==3
  revDesc1=0; revDescIndx1=0;
  mDesc1->fillCellIdsToKeepFromNodeIds(&nodes[0],&nodes[0]+nodes.size(),true,cellIds1D);
  //
  std::vector<int> cut3DCurve(mDesc1->getNumberOfCells(),-2);
  for(std::vector<int>::const_iterator it=cellIds1D.begin();it!=cellIds1D.end();it++)
    cut3DCurve[*it]=-1;
  mDesc1->split3DCurveWithPlane(origin,vec,eps,cut3DCurve);
  std::vector< std::pair<int,int> > cut3DSurf(mDesc2->getNumberOfCells());
  AssemblyForSplitFrom3DCurve(cut3DCurve,nodes,mDesc2->getNodalConnectivity()->getConstPointer(),mDesc2->getNodalConnectivityIndex()->getConstPointer(),
                              mDesc1->getNodalConnectivity()->getConstPointer(),mDesc1->getNodalConnectivityIndex()->getConstPointer(),
                              desc1->getConstPointer(),descIndx1->getConstPointer(),cut3DSurf);
  std::vector<int> conn,connI,cellIds2;
  connI.push_back(0);
  subMesh->assemblyForSplitFrom3DSurf(cut3DSurf,desc2->getConstPointer(),descIndx2->getConstPointer(),conn,connI,cellIds2);
  if(cellIds2.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSlice3D : No 3D cells in this intercepts the specified plane !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> ret=MEDCouplingUMesh::New("Slice3D",2);
  ret->setCoords(mDesc1->getCoords());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> c=DataArrayInt::New();
  c->alloc((int)conn.size(),1); std::copy(conn.begin(),conn.end(),c->getPointer());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cI=DataArrayInt::New();
  cI->alloc((int)connI.size(),1); std::copy(connI.begin(),connI.end(),cI->getPointer());
  ret->setConnectivity(c,cI,true);
  cellIds=candidates->selectByTupleId(&cellIds2[0],&cellIds2[0]+cellIds2.size());
  ret->incrRef();
  return ret;
}

/*!
 * This method expects that 'this' is fully defined and has a spaceDim==3 and a meshDim==2. If it is not the case an exception will be thrown.
 * This method returns 2 objects : 
 * - a newly created mesh instance containing the result of the slice lying on different coords than 'this' and with a meshdim == 1
 * - a newly created dataarray having number of tuples equal to the number of cells in returned mesh that tells for each 2D cell in returned
 *   mesh the 3DSurf cell id is 'this' it comes from.
 * This method works only for linear meshes (non quadratic).
 * If plane crosses within 'eps' a face in 'this' shared by more than 1 cell, 2 output faces will be generated. The 2 faces having the same geometry than intersecting
 * face. Only 'cellIds' parameter can distinguish the 2.
 * @param origin is the origin of the plane. It should be an array of length 3.
 * @param vec is the direction vector of the plane. It should be an array of length 3. Norm of 'vec' should be > 1e-6.
 * @param eps is the precision. It is used by called method MEDCouplingUMesh::getCellIdsCrossingPlane for the first 3DSurf cell selection (in absolute). 'eps' is
 * also used to state if new points should be created or already existing points are reused. 'eps' is also used to tells if plane overlaps a face, edge or nodes (in absolute).
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildSlice3DSurf(const double *origin, const double *vec, double eps, DataArrayInt *&cellIds) const throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  if(getMeshDimension()!=2 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSlice3D works on umeshes with meshdim equal to 2 and spaceDim equal to 3 !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> candidates=getCellIdsCrossingPlane(origin,vec,eps);
  if(candidates->empty())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSlice3D : No 3D cells in this intercepts the specified plane considering bounding boxes !");
  std::vector<int> nodes;
  std::vector<int> cellIds1D;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> subMesh=static_cast<MEDCouplingUMesh*>(buildPartOfMySelf(candidates->begin(),candidates->end(),false));
  subMesh->findNodesOnPlane(origin,vec,eps,nodes);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> desc1=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> descIndx1=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> revDesc1=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> revDescIndx1=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mDesc1=subMesh->buildDescendingConnectivity(desc1,descIndx1,revDesc1,revDescIndx1);//meshDim==1 spaceDim==3
  mDesc1->fillCellIdsToKeepFromNodeIds(&nodes[0],&nodes[0]+nodes.size(),true,cellIds1D);
  //
  std::vector<int> cut3DCurve(mDesc1->getNumberOfCells(),-2);
  for(std::vector<int>::const_iterator it=cellIds1D.begin();it!=cellIds1D.end();it++)
    cut3DCurve[*it]=-1;
  mDesc1->split3DCurveWithPlane(origin,vec,eps,cut3DCurve);
  int ncellsSub=subMesh->getNumberOfCells();
  std::vector< std::pair<int,int> > cut3DSurf(ncellsSub);
  AssemblyForSplitFrom3DCurve(cut3DCurve,nodes,subMesh->getNodalConnectivity()->getConstPointer(),subMesh->getNodalConnectivityIndex()->getConstPointer(),
                              mDesc1->getNodalConnectivity()->getConstPointer(),mDesc1->getNodalConnectivityIndex()->getConstPointer(),
                              desc1->getConstPointer(),descIndx1->getConstPointer(),cut3DSurf);
  std::vector<int> conn,connI,cellIds2; connI.push_back(0);
  const int *nodal=subMesh->getNodalConnectivity()->getConstPointer();
  const int *nodalI=subMesh->getNodalConnectivityIndex()->getConstPointer();
  for(int i=0;i<ncellsSub;i++)
    {
      if(cut3DSurf[i].first!=-1 && cut3DSurf[i].second!=-1)
        {
          if(cut3DSurf[i].first!=-2)
            {
              conn.push_back((int)INTERP_KERNEL::NORM_SEG2); conn.push_back(cut3DSurf[i].first); conn.push_back(cut3DSurf[i].second);
              connI.push_back((int)conn.size());
              cellIds2.push_back(i);
            }
          else
            {
              int cellId3DSurf=cut3DSurf[i].second;
              int offset=nodalI[cellId3DSurf]+1;
              int nbOfEdges=nodalI[cellId3DSurf+1]-offset;
              for(int j=0;j<nbOfEdges;j++)
                {
                  conn.push_back((int)INTERP_KERNEL::NORM_SEG2); conn.push_back(nodal[offset+j]); conn.push_back(nodal[offset+(j+1)%nbOfEdges]);
                  connI.push_back((int)conn.size());
                  cellIds2.push_back(cellId3DSurf);
                }
            }
        }
    }
  if(cellIds2.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSlice3DSurf : No 3DSurf cells in this intercepts the specified plane !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> ret=MEDCouplingUMesh::New("Slice3DSurf",1);
  ret->setCoords(mDesc1->getCoords());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> c=DataArrayInt::New();
  c->alloc((int)conn.size(),1); std::copy(conn.begin(),conn.end(),c->getPointer());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cI=DataArrayInt::New();
  cI->alloc((int)connI.size(),1); std::copy(connI.begin(),connI.end(),cI->getPointer());
  ret->setConnectivity(c,cI,true);
  cellIds=candidates->selectByTupleId(&cellIds2[0],&cellIds2[0]+cellIds2.size());
  ret->incrRef();
  return ret;
}

/*!
 * This method expects that 'this' is fully defined and has a spaceDim==3. If it is not the case an exception will be thrown.
 * This method returns a newly created dataarray containing cellsids in 'this' that potentially crosses the plane specified by 'origin' and 'vec'.
 * @param origin is the origin of the plane. It should be an array of length 3.
 * @param vec is the direction vector of the plane. It should be an array of length 3. Norm of 'vec' should be > 1e-6.
 */
DataArrayInt *MEDCouplingUMesh::getCellIdsCrossingPlane(const double *origin, const double *vec, double eps) const throw(INTERP_KERNEL::Exception)
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
  std::vector<int> cellIds;
  double bbox[6];
  if(angle>eps)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> coo=_coords->deepCpy();
      MEDCouplingPointSet::Rotate3DAlg(origin,vec2,angle,coo->getNumberOfTuples(),coo->getPointer());
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mw=clone(false);//false -> shallow copy
      mw->setCoords(coo);
      mw->getBoundingBox(bbox);
      bbox[4]=origin[2]-eps; bbox[5]=origin[2]+eps;
      mw->getCellsInBoundingBox(bbox,eps,cellIds);
    }
  else
    {
      getBoundingBox(bbox);
      bbox[4]=origin[2]-eps; bbox[5]=origin[2]+eps;
      getCellsInBoundingBox(bbox,eps,cellIds);
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc((int)cellIds.size(),1);
  std::copy(cellIds.begin(),cellIds.end(),ret->getPointer());
  ret->incrRef();
  return ret;
}

/*!
 * This method checks that 'this' is a contiguous mesh. The user is expected to call this method on a mesh with meshdim==1.
 * If not an exception will thrown. If this is an empty mesh with no cell an exception will be thrown too.
 * No consideration of coordinate is done by this method.
 * A 1D mesh is said contiguous if : a cell i with nodal connectivity (k,p) the cell i+1 the nodal connectivity should be (p,m)
 * If not false is returned. In case that false is returned a call to ParaMEDMEM::MEDCouplingUMesh::mergeNodes could be usefull.
 */
bool MEDCouplingUMesh::isContiguous1D() const throw(INTERP_KERNEL::Exception)
{
  if(getMeshDimension()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::isContiguous1D : this method has a sense only for 1D mesh !");
  int nbCells=getNumberOfCells();
  if(nbCells<1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::isContiguous1D : this method has a sense for non empty mesh !");
  const int *connI=_nodal_connec_index->getConstPointer();
  const int *conn=_nodal_connec->getConstPointer();
  int ref=conn[connI[0]+2];
  for(int i=1;i<nbCells;i++)
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
 * @param pt reference point of the line
 * @param v normalized director vector of the line
 * @param eps max precision before throwing an exception
 * @param res output of size this->getNumberOfCells
 */
void MEDCouplingUMesh::project1D(const double *pt, const double *v, double eps, double *res) const
{
  if(getMeshDimension()!=1)
    throw INTERP_KERNEL::Exception("Expected a umesh with meshDim == 1 for project1D !");
   if(_types.size()!=1 || *(_types.begin())!=INTERP_KERNEL::NORM_SEG2)
     throw INTERP_KERNEL::Exception("Expected a umesh with only NORM_SEG2 type of elements for project1D !");
   if(getSpaceDimension()!=3)
     throw INTERP_KERNEL::Exception("Expected a umesh with spaceDim==3 for project1D !");
   MEDCouplingFieldDouble *f=buildDirectionVectorField();
   const double *fPtr=f->getArray()->getConstPointer();
   double tmp[3];
   for(int i=0;i<getNumberOfCells();i++)
     {
       const double *tmp1=fPtr+3*i;
       tmp[0]=tmp1[1]*v[2]-tmp1[2]*v[1];
       tmp[1]=tmp1[2]*v[0]-tmp1[0]*v[2];
       tmp[2]=tmp1[0]*v[1]-tmp1[1]*v[0];
       double n1=INTERP_KERNEL::norm<3>(tmp);
       n1/=INTERP_KERNEL::norm<3>(tmp1);
       if(n1>eps)
         {
           f->decrRef();
           throw INTERP_KERNEL::Exception("UMesh::Projection 1D failed !");
         }
     }
   const double *coo=getCoords()->getConstPointer();
   for(int i=0;i<getNumberOfNodes();i++)
     {
       std::transform(coo+i*3,coo+i*3+3,pt,tmp,std::minus<double>());
       std::transform(tmp,tmp+3,v,tmp,std::multiplies<double>());
       res[i]=std::accumulate(tmp,tmp+3,0.);
     }
   f->decrRef();
}

/*!
 * Returns a cell if any that contains the point located on 'pos' with precison eps.
 * If 'pos' is outside 'this' -1 is returned. If several cells contain this point the cell with the smallest id is returned.
 * \b Warning this method is good if the caller intends to evaluate only one point. But if more than one point is requested on 'this'
 * it is better to use MEDCouplingUMesh::getCellsContainingPoints method because in this case, the acceleration structure will be computed only once.
 */
int MEDCouplingUMesh::getCellContainingPoint(const double *pos, double eps) const
{
  std::vector<int> elts;
  getCellsContainingPoint(pos,eps,elts);
  if(elts.empty())
    return -1;
  return elts.front();
}

/*!
 * Returns all cellIds in 'elts' of point 'pos' with eps accuracy.
 * \b Warning this method is good if the caller intends to evaluate only one point. But if more than one point is requested on 'this'
 * it is better to use MEDCouplingUMesh::getCellsContainingPoints method because in this case, the acceleration structure will be computed only once.
 */
void MEDCouplingUMesh::getCellsContainingPoint(const double *pos, double eps, std::vector<int>& elts) const
{
  std::vector<int> eltsIndex;
  getCellsContainingPoints(pos,1,eps,elts,eltsIndex);
}

/// @cond INTERNAL

namespace ParaMEDMEM
{
  template<const int SPACEDIMM>
  class DummyClsMCUG
  {
  public:
    static const int MY_SPACEDIM=SPACEDIMM;
    static const int MY_MESHDIM=8;
    typedef int MyConnType;
    static const INTERP_KERNEL::NumberingPolicy My_numPol=INTERP_KERNEL::ALL_C_MODE;
    // begin
    // useless, but for windows compilation ...
    const double* getCoordinatesPtr() const { return 0; }
    const int* getConnectivityPtr() const { return 0; }
    const int* getConnectivityIndexPtr() const { return 0; }
    INTERP_KERNEL::NormalizedCellType getTypeOfElement(int) const { return (INTERP_KERNEL::NormalizedCellType)0; }
    // end
  };

  INTERP_KERNEL::Edge *MEDCouplingUMeshBuildQPFromEdge(INTERP_KERNEL::NormalizedCellType typ, std::map<int, std::pair<INTERP_KERNEL::Node *,bool> >& mapp2, const int *bg)
  {
    INTERP_KERNEL::Edge *ret=0;
    switch(typ)
      {
      case INTERP_KERNEL::NORM_SEG2:
        {
          ret=new INTERP_KERNEL::EdgeLin(mapp2[bg[0]].first,mapp2[bg[1]].first);
          break;
        }
      case INTERP_KERNEL::NORM_SEG3:
        {
          INTERP_KERNEL::EdgeLin *e1=new INTERP_KERNEL::EdgeLin(mapp2[bg[0]].first,mapp2[bg[2]].first);
          INTERP_KERNEL::EdgeLin *e2=new INTERP_KERNEL::EdgeLin(mapp2[bg[2]].first,mapp2[bg[1]].first);
          INTERP_KERNEL::SegSegIntersector inters(*e1,*e2);
          bool colinearity=inters.areColinears();
          delete e1; delete e2;
          if(colinearity)
            ret=new INTERP_KERNEL::EdgeLin(mapp2[bg[0]].first,mapp2[bg[1]].first);
          else
            ret=new INTERP_KERNEL::EdgeArcCircle(mapp2[bg[0]].first,mapp2[bg[2]].first,mapp2[bg[1]].first);
          mapp2[bg[2]].second=false;
          break;
        }
      default:
        throw INTERP_KERNEL::Exception("MEDCouplingUMeshBuildQPFromEdge : Expecting a mesh with spaceDim==2 and meshDim==1 !");
      }
    return ret;
  }

  /*!
   * This method creates a sub mesh in Geometric2D DS. The sub mesh is composed be the sub set of cells in 'candidates' and the global mesh 'mDesc'.
   * The input meth 'mDesc' must be so that mDim==1 et spaceDim==3.
   * 'mapp' contains a mapping between local numbering in submesh and the global node numbering in 'mDesc'.
   */
  INTERP_KERNEL::QuadraticPolygon *MEDCouplingUMeshBuildQPFromMesh(const MEDCouplingUMesh *mDesc, const std::vector<int>& candidates, std::map<INTERP_KERNEL::Node *,int>& mapp) throw(INTERP_KERNEL::Exception)
  {
    mapp.clear();
    std::map<int, std::pair<INTERP_KERNEL::Node *,bool> > mapp2;//bool is for a flag specifying if node is boundary (true) or only a middle for SEG3.
    const double *coo=mDesc->getCoords()->getConstPointer();
    const int *c=mDesc->getNodalConnectivity()->getConstPointer();
    const int *cI=mDesc->getNodalConnectivityIndex()->getConstPointer();
    std::set<int> s;
    for(std::vector<int>::const_iterator it=candidates.begin();it!=candidates.end();it++)
      s.insert(c+cI[*it]+1,c+cI[(*it)+1]);
    for(std::set<int>::const_iterator it2=s.begin();it2!=s.end();it2++)
      {
        INTERP_KERNEL::Node *n=new INTERP_KERNEL::Node(coo[2*(*it2)],coo[2*(*it2)+1]);
        mapp2[*it2]=std::pair<INTERP_KERNEL::Node *,bool>(n,true);
      }
    INTERP_KERNEL::QuadraticPolygon *ret=new INTERP_KERNEL::QuadraticPolygon;
    for(std::vector<int>::const_iterator it=candidates.begin();it!=candidates.end();it++)
      {
        INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType)c[cI[*it]];
        ret->pushBack(MEDCouplingUMeshBuildQPFromEdge(typ,mapp2,c+cI[*it]+1));
      }
    for(std::map<int, std::pair<INTERP_KERNEL::Node *,bool> >::const_iterator it2=mapp2.begin();it2!=mapp2.end();it2++)
      {
        if((*it2).second.second)
          mapp[(*it2).second.first]=(*it2).first;
        ((*it2).second.first)->decrRef();
      }
    return ret;
  }

  INTERP_KERNEL::Node *MEDCouplingUMeshBuildQPNode(int nodeId, const double *coo1, int offset1, const double *coo2, int offset2, const std::vector<double>& addCoo)
  {
    if(nodeId>=offset2)
      {
        int locId=nodeId-offset2;
        return new INTERP_KERNEL::Node(addCoo[2*locId],addCoo[2*locId+1]);
      }
    if(nodeId>=offset1)
      {
        int locId=nodeId-offset1;
        return new INTERP_KERNEL::Node(coo2[2*locId],coo2[2*locId+1]);
      }
    return new INTERP_KERNEL::Node(coo1[2*nodeId],coo1[2*nodeId+1]);
  }

  void MEDCouplingUMeshBuildQPFromMesh3(const double *coo1, int offset1, const double *coo2, int offset2, const std::vector<double>& addCoo,
                                        const int *desc1Bg, const int *desc1End, const std::vector<std::vector<int> >& intesctEdges1,
                                        /*output*/std::map<INTERP_KERNEL::Node *,int>& mapp, std::map<int,INTERP_KERNEL::Node *>& mappRev)
  {
    for(const int *desc1=desc1Bg;desc1!=desc1End;desc1++)
      {
        int eltId1=abs(*desc1)-1;
        for(std::vector<int>::const_iterator it1=intesctEdges1[eltId1].begin();it1!=intesctEdges1[eltId1].end();it1++)
          {
            std::map<int,INTERP_KERNEL::Node *>::const_iterator it=mappRev.find(*it1);
            if(it==mappRev.end())
              {
                INTERP_KERNEL::Node *node=MEDCouplingUMeshBuildQPNode(*it1,coo1,offset1,coo2,offset2,addCoo);
                mapp[node]=*it1;
                mappRev[*it1]=node;
              }
          }
      }
  }
}

/// @endcond

template<int SPACEDIM>
void MEDCouplingUMesh::getCellsContainingPointsAlg(const double *coords, const double *pos, int nbOfPoints,
                                                   double eps, std::vector<int>& elts, std::vector<int>& eltsIndex) const
{
  std::vector<double> bbox;
  eltsIndex.resize(nbOfPoints+1);
  eltsIndex[0]=0;
  elts.clear();
  getBoundingBoxForBBTree(bbox);
  int nbOfCells=getNumberOfCells();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  double bb[2*SPACEDIM];
  BBTree<SPACEDIM,int> myTree(&bbox[0],0,0,nbOfCells,-eps);
  for(int i=0;i<nbOfPoints;i++)
    {
      eltsIndex[i+1]=eltsIndex[i];
      for(int j=0;j<SPACEDIM;j++)
        {
          bb[2*j]=pos[SPACEDIM*i+j];
          bb[2*j+1]=pos[SPACEDIM*i+j];
        }
      std::vector<int> candidates;
      myTree.getIntersectingElems(bb,candidates);
      for(std::vector<int>::const_iterator iter=candidates.begin();iter!=candidates.end();iter++)
        {
          int sz=connI[(*iter)+1]-connI[*iter]-1;
          if(INTERP_KERNEL::PointLocatorAlgos<DummyClsMCUG<SPACEDIM> >::isElementContainsPoint(pos+i*SPACEDIM,
                                                                                               (INTERP_KERNEL::NormalizedCellType)conn[connI[*iter]],
                                                                                               coords,conn+connI[*iter]+1,sz,eps))
            {
              eltsIndex[i+1]++;
              elts.push_back(*iter);
            }
        }
    }
}

/*!
 * This method is an extension of MEDCouplingUMesh::getCellContainingPoint and MEDCouplingUMesh::getCellsContainingPoint.
 * This method performs 'nbOfPoints' time the getCellsContainingPoint request. This method is recommended rather than the 2 others
 * in case of multi points searching.
 * This method returns 2 arrays 'elts' and 'eltsIndex'. 'eltsIndex' is of size 'nbOfPoints+1' and 'elts' is of size 'eltsIndex[nbOfPoints-1]'.
 * For point j in [0,nbOfPoints), (eltsIndex[j+1]-eltsIndex[j]) cells contain this point. These cells are : [elts.begin()+eltsIndex[j],elts.begin():eltsIndex[j+1]).
 * 
 * \param pos input parameter that points to an array of size 'getSpaceDim()*nbOfPoints' points stored in full interlace mode : X0,Y0,Z0,X1,Y1,Z1...
 */
void MEDCouplingUMesh::getCellsContainingPoints(const double *pos, int nbOfPoints, double eps,
                                                std::vector<int>& elts, std::vector<int>& eltsIndex) const
{
  int spaceDim=getSpaceDimension();
  int mDim=getMeshDimension();
  if(spaceDim==3)
    {
      if(mDim==3)
        {
          const double *coords=_coords->getConstPointer();
          getCellsContainingPointsAlg<3>(coords,pos,nbOfPoints,eps,elts,eltsIndex);
        }
      /*else if(mDim==2)
        {
          
        }*/
      else
        throw INTERP_KERNEL::Exception("For spaceDim==3 only meshDim==3 implemented for getelementscontainingpoints !");
    }
  else if(spaceDim==2)
    {
      if(mDim==2)
        {
          const double *coords=_coords->getConstPointer();
          getCellsContainingPointsAlg<2>(coords,pos,nbOfPoints,eps,elts,eltsIndex);
        }
      else
        throw INTERP_KERNEL::Exception("For spaceDim==2 only meshDim==2 implemented for getelementscontainingpoints !");
    }
  else if(spaceDim==1)
    {
      if(mDim==1)
        {
          const double *coords=_coords->getConstPointer();
          getCellsContainingPointsAlg<1>(coords,pos,nbOfPoints,eps,elts,eltsIndex);
        }
      else
        throw INTERP_KERNEL::Exception("For spaceDim==1 only meshDim==1 implemented for getelementscontainingpoints !");
    }
  else
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getCellsContainingPoints : not managed for mdim not in [1,2,3] !");
}

/*!
 * This method is only available for a mesh with meshDim==2 and spaceDim==2||spaceDim==3.
 * This method returns a vector 'cells' where all detected butterfly cells have been added to cells.
 * A 2D cell is considered to be butterfly if it exists at least one pair of distinct edges of it that intersect each other
 * anywhere excepted their extremities. An INTERP_KERNEL::NORM_NORI3 could \b not be butterfly.
 */
void MEDCouplingUMesh::checkButterflyCells(std::vector<int>& cells, double eps) const
{
  const char msg[]="Butterfly detection work only for 2D cells with spaceDim==2 or 3!";
  if(getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception(msg);
  int spaceDim=getSpaceDimension();
  if(spaceDim!=2 && spaceDim!=3)
    throw INTERP_KERNEL::Exception(msg);
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  int nbOfCells=getNumberOfCells();
  std::vector<double> cell2DinS2;
  for(int i=0;i<nbOfCells;i++)
    {
      int offset=connI[i];
      int nbOfNodesForCell=connI[i+1]-offset-1;
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
 * @return a newly allocated array containing cellIds that have been modified if any. If no cells have been impacted by this method NULL is returned.
 */
DataArrayInt *MEDCouplingUMesh::convexEnvelop2D() throw(INTERP_KERNEL::Exception)
{
  if(getMeshDimension()!=2 || getSpaceDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::convexEnvelop2D  works only for meshDim=2 and spaceDim=2 !");
  checkFullyDefined();
  const double *coords=getCoords()->getConstPointer();
  int nbOfCells=getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> nodalConnecIndexOut=DataArrayInt::New();
  nodalConnecIndexOut->alloc(nbOfCells+1,1);
  std::vector<int> nodalConnecOut;
  int *workIndexOut=nodalConnecIndexOut->getPointer();
  *workIndexOut=0;
  const int *nodalConnecIn=_nodal_connec->getConstPointer();
  const int *nodalConnecIndexIn=_nodal_connec_index->getConstPointer();
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  std::vector<int> isChanged;
  for(int i=0;i<nbOfCells;i++,workIndexOut++)
    {
      std::size_t pos=nodalConnecOut.size();
      if(BuildConvexEnvelopOf2DCellJarvis(coords,nodalConnecIn+nodalConnecIndexIn[i],nodalConnecIn+nodalConnecIndexIn[i+1],nodalConnecOut))
        isChanged.push_back(i);
      types.insert((INTERP_KERNEL::NormalizedCellType)nodalConnecOut[pos]);
      workIndexOut[1]=(int)nodalConnecOut.size();
    }
  if(isChanged.empty())
    return 0;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> nodalConnecOut2=DataArrayInt::New();
  nodalConnecOut2->alloc((int)nodalConnecOut.size(),1);
  std::copy(nodalConnecOut.begin(),nodalConnecOut.end(),nodalConnecOut2->getPointer());
  setConnectivity(nodalConnecOut2,nodalConnecIndexOut,false);
  _types=types;
  DataArrayInt *ret=DataArrayInt::New(); ret->alloc((int)isChanged.size(),1);
  std::copy(isChanged.begin(),isChanged.end(),ret->getPointer());
  return ret;
}

/*!
 * This method is expected to be applied on a mesh with spaceDim==3 and meshDim==3. If not an exception will be thrown.
 * This method analyzes only linear extruded 3D cells (NORM_HEXA8,NORM_PENTA6,NORM_HEXGP12...)
 * If some extruded cells does not fulfill the MED norm for extruded cells (first face of 3D cell should be oriented to the exterior of the 3D cell).
 * Some viewers are very careful of that (SMESH), but ParaVis ignore that.
 */
void MEDCouplingUMesh::findAndCorrectBadOriented3DExtrudedCells(std::vector<int>& cells) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="check3DCellsWellOriented detection works only for 3D cells !";
  if(getMeshDimension()!=3)
    throw INTERP_KERNEL::Exception(msg);
  int spaceDim=getSpaceDimension();
  if(spaceDim!=3)
    throw INTERP_KERNEL::Exception(msg);
  //
  int nbOfCells=getNumberOfCells();
  int *conn=_nodal_connec->getPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const double *coo=getCoords()->getConstPointer();
  double vec0[3],vec1[3];
  for(int i=0;i<nbOfCells;i++)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[connI[i]]);
      if(cm.isExtruded() && !cm.isDynamic() && !cm.isQuadratic())
        {
          INTERP_KERNEL::AutoPtr<int> tmp=new int[connI[i+1]-connI[i]-1];
          int nbOfNodes=cm.fillSonCellNodalConnectivity(0,conn+connI[i]+1,tmp);
          INTERP_KERNEL::areaVectorOfPolygon<int,INTERP_KERNEL::ALL_C_MODE>(tmp,nbOfNodes,coo,vec0);
          const double *pt0=coo+3*conn[connI[i]+1];
          const double *pt1=coo+3*conn[connI[i]+nbOfNodes+1];
          vec1[0]=pt0[0]-pt1[0]; vec1[1]=pt0[1]-pt1[1]; vec1[2]=pt0[2]-pt1[2];
          double dot=vec0[0]*vec1[0]+vec0[1]*vec1[1]+vec0[2]*vec1[2];
          if(dot<0)
            {
              cells.push_back(i);
              std::copy(conn+connI[i]+1,conn+connI[i+1],(int *)tmp);
              for(int j=1;j<nbOfNodes;j++)
                {
                  conn[connI[i]+1+j]=tmp[nbOfNodes-j];
                  conn[connI[i]+1+j+nbOfNodes]=tmp[nbOfNodes+nbOfNodes-j];
                }
            }
        }
    }
}

/*!
 * This method is \b NOT const because it can modify 'this'.
 * 'this' is expected to be an unstructured mesh with meshDim==2 and spaceDim==3. If not an exception will be thrown.
 * @param mesh1D is an unstructured mesh with MeshDim==1 and spaceDim==3. If not an exception will be thrown.
 * @param policy specifies the type of extrusion chosen. \b 0 for translation (most simple),
 * \b 1 for translation and rotation around point of 'mesh1D'.
 * @return an unstructured mesh with meshDim==3 and spaceDim==3. The returned mesh has the same coords than 'this'.  
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildExtrudedMesh(const MEDCouplingUMesh *mesh1D, int policy)
{
  checkFullyDefined();
  mesh1D->checkFullyDefined();
  if(!mesh1D->isContiguous1D())
    throw INTERP_KERNEL::Exception("buildExtrudedMesh : 1D mesh passed in parameter is not contiguous !");
  if(getSpaceDimension()!=mesh1D->getSpaceDimension())
    throw INTERP_KERNEL::Exception("Invalid call to buildExtrudedMesh this and mesh1D must have same dimension !");
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
  zipCoords();
  int oldNbOfNodes=getNumberOfNodes();
  DataArrayDouble *newCoords=0;
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
  newCoords->decrRef();
  MEDCouplingUMesh *ret=buildExtrudedMeshFromThisLowLev(oldNbOfNodes,isQuad);
  updateTime();
  return ret;
}

/*!
 * This method works on a 3D curve linear mesh that is to say (meshDim==1 and spaceDim==3).
 * If it is not the case an exception will be thrown.
 * This method is non const because the coordinate of 'this' can be appended with some new points issued from
 * intersection of plane defined by ('origin','vec').
 * This method has one in/out parameter : 'cut3DCurve'.
 * Param 'cut3DCurve' is expected to be of size 'this->getNumberOfCells()'. For each i in [0,'this->getNumberOfCells()')
 * if cut3DCurve[i]==-2, it means that for cell #i in 'this' nothing has been detected previously.
 * if cut3DCurve[i]==-1, it means that cell#i has been already detected to be fully part of plane defined by ('origin','vec').
 * This method will throw an exception if 'this' contains a non linear segment.
 */
void MEDCouplingUMesh::split3DCurveWithPlane(const double *origin, const double *vec, double eps, std::vector<int>& cut3DCurve) throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  if(getMeshDimension()!=1 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::split3DCurveWithPlane works on umeshes with meshdim equal to 1 and spaceDim equal to 3 !");
  int ncells=getNumberOfCells();
  int nnodes=getNumberOfNodes();
  double vec2[3],vec3[3],vec4[3];
  double normm=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  if(normm<1e-6)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::split3DCurveWithPlane : parameter 'vec' should have a norm2 greater than 1e-6 !");
  vec2[0]=vec[0]/normm; vec2[1]=vec[1]/normm; vec2[2]=vec[2]/normm;
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const double *coo=_coords->getConstPointer();
  std::vector<double> addCoo;
  for(int i=0;i<ncells;i++)
    {
      if(conn[connI[i]]==(int)INTERP_KERNEL::NORM_SEG2)
        {
          if(cut3DCurve[i]==-2)
            {
              int st=conn[connI[i]+1],endd=conn[connI[i]+2];
              vec3[0]=coo[3*endd]-coo[3*st]; vec3[1]=coo[3*endd+1]-coo[3*st+1]; vec3[2]=coo[3*endd+2]-coo[3*st+2];
              double normm2=sqrt(vec3[0]*vec3[0]+vec3[1]*vec3[1]+vec3[2]*vec3[2]);
              double colin=std::abs((vec3[0]*vec2[0]+vec3[1]*vec2[1]+vec3[2]*vec2[2])/normm2);
              if(colin>eps)//if colin<=eps -> current SEG2 is colinear to the input plane
                {
                  const double *st2=coo+3*st;
                  vec4[0]=st2[0]-origin[0]; vec4[1]=st2[1]-origin[1]; vec4[2]=st2[2]-origin[2];
                  double pos=-(vec4[0]*vec2[0]+vec4[1]*vec2[1]+vec4[2]*vec2[2])/((vec3[0]*vec2[0]+vec3[1]*vec2[1]+vec3[2]*vec2[2]));
                  if(pos>eps && pos<1-eps)
                    {
                      int nNode=((int)addCoo.size())/3;
                      vec4[0]=st2[0]+pos*vec3[0]; vec4[1]=st2[1]+pos*vec3[1]; vec4[2]=st2[2]+pos*vec3[2];
                      addCoo.insert(addCoo.end(),vec4,vec4+3);
                      cut3DCurve[i]=nnodes+nNode;
                    }
                }
            }
        }
      else
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::split3DCurveWithPlane : this method is only available for linear cell (NORM_SEG2) !");
    }
  if(!addCoo.empty())
    {
      int newNbOfNodes=nnodes+((int)addCoo.size())/3;
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> coo2=DataArrayDouble::New();
      coo2->alloc(newNbOfNodes,3);
      double *tmp=coo2->getPointer();
      tmp=std::copy(_coords->begin(),_coords->end(),tmp);
      std::copy(addCoo.begin(),addCoo.end(),tmp);
      DataArrayDouble::SetArrayIn(coo2,_coords);
    }
}

/*!
 * This method incarnates the policy 0 for MEDCouplingUMesh::buildExtrudedMesh method.
 * @param mesh1D is the input 1D mesh used for translation computation.
 * @return newCoords new coords filled by this method. 
 */
DataArrayDouble *MEDCouplingUMesh::fillExtCoordsUsingTranslation(const MEDCouplingUMesh *mesh1D, bool isQuad) const
{
  int oldNbOfNodes=getNumberOfNodes();
  int nbOf1DCells=mesh1D->getNumberOfCells();
  int spaceDim=getSpaceDimension();
  DataArrayDouble *ret=DataArrayDouble::New();
  std::vector<bool> isQuads;
  int nbOfLevsInVec=isQuad?2*nbOf1DCells+1:nbOf1DCells+1;
  ret->alloc(oldNbOfNodes*nbOfLevsInVec,spaceDim);
  double *retPtr=ret->getPointer();
  const double *coords=getCoords()->getConstPointer();
  double *work=std::copy(coords,coords+spaceDim*oldNbOfNodes,retPtr);
  std::vector<int> v;
  std::vector<double> c;
  double vec[3];
  v.reserve(3);
  c.reserve(6);
  for(int i=0;i<nbOf1DCells;i++)
    {
      v.resize(0);
      mesh1D->getNodeIdsOfCell(i,v);
      c.resize(0);
      mesh1D->getCoordinatesOfNode(v[isQuad?2:1],c);
      mesh1D->getCoordinatesOfNode(v[0],c);
      std::transform(c.begin(),c.begin()+spaceDim,c.begin()+spaceDim,vec,std::minus<double>());
      for(int j=0;j<oldNbOfNodes;j++)
        work=std::transform(vec,vec+spaceDim,retPtr+spaceDim*(i*oldNbOfNodes+j),work,std::plus<double>());
      if(isQuad)
        {
          c.resize(0);
          mesh1D->getCoordinatesOfNode(v[1],c);
          mesh1D->getCoordinatesOfNode(v[0],c);
          std::transform(c.begin(),c.begin()+spaceDim,c.begin()+spaceDim,vec,std::minus<double>());
          for(int j=0;j<oldNbOfNodes;j++)
            work=std::transform(vec,vec+spaceDim,retPtr+spaceDim*(i*oldNbOfNodes+j),work,std::plus<double>());
        }
    }
  ret->copyStringInfoFrom(*getCoords());
  return ret;
}

/*!
 * This method incarnates the policy 1 for MEDCouplingUMesh::buildExtrudedMesh method.
 * @param mesh1D is the input 1D mesh used for translation and automatic rotation computation.
 * @return newCoords new coords filled by this method. 
 */
DataArrayDouble *MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation(const MEDCouplingUMesh *mesh1D, bool isQuad) const throw(INTERP_KERNEL::Exception)
{
  if(mesh1D->getSpaceDimension()==2)
    return fillExtCoordsUsingTranslAndAutoRotation2D(mesh1D,isQuad);
  if(mesh1D->getSpaceDimension()==3)
    return fillExtCoordsUsingTranslAndAutoRotation3D(mesh1D,isQuad);
  throw INTERP_KERNEL::Exception("Not implemented rotation and translation alg. for spacedim other than 2 and 3 !");
}

/*!
 * This method incarnates the policy 1 for MEDCouplingUMesh::buildExtrudedMesh method.
 * @param mesh1D is the input 1D mesh used for translation and automatic rotation computation.
 * @return newCoords new coords filled by this method. 
 */
DataArrayDouble *MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation2D(const MEDCouplingUMesh *mesh1D, bool isQuad) const throw(INTERP_KERNEL::Exception)
{
  if(isQuad)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation2D : not implemented for quadratic cells !");
  int oldNbOfNodes=getNumberOfNodes();
  int nbOf1DCells=mesh1D->getNumberOfCells();
  if(nbOf1DCells<2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation2D : impossible to detect any angle of rotation ! Change extrusion policy 1->0 !");
  DataArrayDouble *ret=DataArrayDouble::New();
  int nbOfLevsInVec=nbOf1DCells+1;
  ret->alloc(oldNbOfNodes*nbOfLevsInVec,2);
  double *retPtr=ret->getPointer();
  retPtr=std::copy(getCoords()->getConstPointer(),getCoords()->getConstPointer()+getCoords()->getNbOfElems(),retPtr);
  MEDCouplingUMesh *tmp=MEDCouplingUMesh::New();
  DataArrayDouble *tmp2=getCoords()->deepCpy();
  tmp->setCoords(tmp2);
  tmp2->decrRef();
  const double *coo1D=mesh1D->getCoords()->getConstPointer();
  const int *conn1D=mesh1D->getNodalConnectivity()->getConstPointer();
  const int *connI1D=mesh1D->getNodalConnectivityIndex()->getConstPointer();
  for(int i=1;i<nbOfLevsInVec;i++)
    {
      const double *begin=coo1D+2*conn1D[connI1D[i-1]+1];
      const double *end=coo1D+2*conn1D[connI1D[i-1]+2];
      const double *third=i+1<nbOfLevsInVec?coo1D+2*conn1D[connI1D[i]+2]:coo1D+2*conn1D[connI1D[i-2]+1];
      const double vec[2]={end[0]-begin[0],end[1]-begin[1]};
      tmp->translate(vec);
      double tmp3[2],radius,alpha,alpha0;
      const double *p0=i+1<nbOfLevsInVec?begin:third;
      const double *p1=i+1<nbOfLevsInVec?end:begin;
      const double *p2=i+1<nbOfLevsInVec?third:end;
      INTERP_KERNEL::EdgeArcCircle::GetArcOfCirclePassingThru(p0,p1,p2,tmp3,radius,alpha,alpha0);
      double cosangle=i+1<nbOfLevsInVec?(p0[0]-tmp3[0])*(p1[0]-tmp3[0])+(p0[1]-tmp3[1])*(p1[1]-tmp3[1]):(p2[0]-tmp3[0])*(p1[0]-tmp3[0])+(p2[1]-tmp3[1])*(p1[1]-tmp3[1]);
      double angle=acos(cosangle/(radius*radius));
      tmp->rotate(end,0,angle);
      retPtr=std::copy(tmp2->getConstPointer(),tmp2->getConstPointer()+tmp2->getNbOfElems(),retPtr);
    }
  tmp->decrRef();
  return ret;
}

/*!
 * This method incarnates the policy 1 for MEDCouplingUMesh::buildExtrudedMesh method.
 * @param mesh1D is the input 1D mesh used for translation and automatic rotation computation.
 * @return newCoords new coords filled by this method. 
 */
DataArrayDouble *MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation3D(const MEDCouplingUMesh *mesh1D, bool isQuad) const throw(INTERP_KERNEL::Exception)
{
  if(isQuad)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation3D : not implemented for quadratic cells !");
  int oldNbOfNodes=getNumberOfNodes();
  int nbOf1DCells=mesh1D->getNumberOfCells();
  if(nbOf1DCells<2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation3D : impossible to detect any angle of rotation ! Change extrusion policy 1->0 !");
  DataArrayDouble *ret=DataArrayDouble::New();
  int nbOfLevsInVec=nbOf1DCells+1;
  ret->alloc(oldNbOfNodes*nbOfLevsInVec,3);
  double *retPtr=ret->getPointer();
  retPtr=std::copy(getCoords()->getConstPointer(),getCoords()->getConstPointer()+getCoords()->getNbOfElems(),retPtr);
  MEDCouplingUMesh *tmp=MEDCouplingUMesh::New();
  DataArrayDouble *tmp2=getCoords()->deepCpy();
  tmp->setCoords(tmp2);
  tmp2->decrRef();
  const double *coo1D=mesh1D->getCoords()->getConstPointer();
  const int *conn1D=mesh1D->getNodalConnectivity()->getConstPointer();
  const int *connI1D=mesh1D->getNodalConnectivityIndex()->getConstPointer();
  for(int i=1;i<nbOfLevsInVec;i++)
    {
      const double *begin=coo1D+3*conn1D[connI1D[i-1]+1];
      const double *end=coo1D+3*conn1D[connI1D[i-1]+2];
      const double *third=i+1<nbOfLevsInVec?coo1D+3*conn1D[connI1D[i]+2]:coo1D+3*conn1D[connI1D[i-2]+1];
      const double vec[3]={end[0]-begin[0],end[1]-begin[1],end[2]-begin[2]};
      tmp->translate(vec);
      double tmp3[2],radius,alpha,alpha0;
      const double *p0=i+1<nbOfLevsInVec?begin:third;
      const double *p1=i+1<nbOfLevsInVec?end:begin;
      const double *p2=i+1<nbOfLevsInVec?third:end;
      double vecPlane[3]={
        (p1[1]-p0[1])*(p2[2]-p1[2])-(p1[2]-p0[2])*(p2[1]-p1[1]),
        (p1[2]-p0[2])*(p2[0]-p1[0])-(p1[0]-p0[0])*(p2[2]-p1[2]),
        (p1[0]-p0[0])*(p2[1]-p1[1])-(p1[1]-p0[1])*(p2[0]-p1[0]),
      };
      double norm=sqrt(vecPlane[0]*vecPlane[0]+vecPlane[1]*vecPlane[1]+vecPlane[2]*vecPlane[2]);
      if(norm>1.e-7)
        {
          vecPlane[0]/=norm; vecPlane[1]/=norm; vecPlane[2]/=norm;
          double norm2=sqrt(vecPlane[0]*vecPlane[0]+vecPlane[1]*vecPlane[1]);
          double vec2[2]={vecPlane[1]/norm2,-vecPlane[0]/norm2};
          double s2=norm2;
          double c2=cos(asin(s2));
          double m[3][3]={
            {vec2[0]*vec2[0]*(1-c2)+c2, vec2[0]*vec2[1]*(1-c2), vec2[1]*s2},
            {vec2[0]*vec2[1]*(1-c2), vec2[1]*vec2[1]*(1-c2)+c2, -vec2[0]*s2},
            {-vec2[1]*s2, vec2[0]*s2, c2}
          };
          double p0r[3]={m[0][0]*p0[0]+m[0][1]*p0[1]+m[0][2]*p0[2], m[1][0]*p0[0]+m[1][1]*p0[1]+m[1][2]*p0[2], m[2][0]*p0[0]+m[2][1]*p0[1]+m[2][2]*p0[2]};
          double p1r[3]={m[0][0]*p1[0]+m[0][1]*p1[1]+m[0][2]*p1[2], m[1][0]*p1[0]+m[1][1]*p1[1]+m[1][2]*p1[2], m[2][0]*p1[0]+m[2][1]*p1[1]+m[2][2]*p1[2]};
          double p2r[3]={m[0][0]*p2[0]+m[0][1]*p2[1]+m[0][2]*p2[2], m[1][0]*p2[0]+m[1][1]*p2[1]+m[1][2]*p2[2], m[2][0]*p2[0]+m[2][1]*p2[1]+m[2][2]*p2[2]};
          INTERP_KERNEL::EdgeArcCircle::GetArcOfCirclePassingThru(p0r,p1r,p2r,tmp3,radius,alpha,alpha0);
          double cosangle=i+1<nbOfLevsInVec?(p0r[0]-tmp3[0])*(p1r[0]-tmp3[0])+(p0r[1]-tmp3[1])*(p1r[1]-tmp3[1]):(p2r[0]-tmp3[0])*(p1r[0]-tmp3[0])+(p2r[1]-tmp3[1])*(p1r[1]-tmp3[1]);
          double angle=acos(cosangle/(radius*radius));
          tmp->rotate(end,vecPlane,angle);
          
        }
      retPtr=std::copy(tmp2->getConstPointer(),tmp2->getConstPointer()+tmp2->getNbOfElems(),retPtr);
    }
  tmp->decrRef();
  return ret;
}

/*!
 * This method is private because not easy to use for end user. This method is const contrary to
 * MEDCouplingUMesh::buildExtrudedMesh method because this->_coords are expected to contain
 * the coords sorted slice by slice.
 * @param isQuad specifies presence of quadratic cells.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildExtrudedMeshFromThisLowLev(int nbOfNodesOf1Lev, bool isQuad) const
{
  int nbOf1DCells=getNumberOfNodes()/nbOfNodesOf1Lev-1;
  int nbOf2DCells=getNumberOfCells();
  int nbOf3DCells=nbOf2DCells*nbOf1DCells;
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New("Extruded",getMeshDimension()+1);
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  DataArrayInt *newConn=DataArrayInt::New();
  DataArrayInt *newConnI=DataArrayInt::New();
  newConnI->alloc(nbOf3DCells+1,1);
  int *newConnIPtr=newConnI->getPointer();
  *newConnIPtr++=0;
  std::vector<int> newc;
  for(int j=0;j<nbOf2DCells;j++)
    {
      AppendExtrudedCell(conn+connI[j],conn+connI[j+1],nbOfNodesOf1Lev,isQuad,newc);
      *newConnIPtr++=(int)newc.size();
    }
  newConn->alloc((int)(newc.size())*nbOf1DCells,1);
  int *newConnPtr=newConn->getPointer();
  int deltaPerLev=isQuad?2*nbOfNodesOf1Lev:nbOfNodesOf1Lev;
  newConnIPtr=newConnI->getPointer();
  for(int iz=0;iz<nbOf1DCells;iz++)
    {
      if(iz!=0)
        std::transform(newConnIPtr+1,newConnIPtr+1+nbOf2DCells,newConnIPtr+1+iz*nbOf2DCells,std::bind2nd(std::plus<int>(),newConnIPtr[iz*nbOf2DCells]));
      for(std::vector<int>::const_iterator iter=newc.begin();iter!=newc.end();iter++,newConnPtr++)
        {
          int icell=(int)(iter-newc.begin());
          if(std::find(newConnIPtr,newConnIPtr+nbOf2DCells,icell)==newConnIPtr+nbOf2DCells)
            {
              if(*iter!=-1)
                *newConnPtr=(*iter)+iz*deltaPerLev;
              else
                *newConnPtr=-1;
            }
          else
            *newConnPtr=(*iter);
        }
    }
  ret->setConnectivity(newConn,newConnI,true);
  newConn->decrRef();
  newConnI->decrRef();
  ret->setCoords(getCoords());
  return ret;
}

/*!
 * This method returns if 'this' is constituted by only quadratic cells.
 */
bool MEDCouplingUMesh::isFullyQuadratic() const
{
  checkFullyDefined();
  bool ret=true;
  int nbOfCells=getNumberOfCells();
  for(int i=0;i<nbOfCells && ret;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      ret=cm.isQuadratic();
    }
  return ret;
}

/*!
 * This method returns if there is at least one quadratic cell.
 */
bool MEDCouplingUMesh::isPresenceOfQuadratic() const
{
  checkFullyDefined();
  bool ret=false;
  int nbOfCells=getNumberOfCells();
  for(int i=0;i<nbOfCells && !ret;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      ret=cm.isQuadratic();
    }
  return ret;
}

/*!
 * This method convert quadratic cells to linear cells if any was found.
 * If no such cells exists 'this' remains unchanged.
 */
void MEDCouplingUMesh::convertQuadraticCellsToLinear() throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  int nbOfCells=getNumberOfCells();
  int delta=0;
  for(int i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      if(cm.isQuadratic())
        {
          INTERP_KERNEL::NormalizedCellType typel=cm.getLinearType();
          const INTERP_KERNEL::CellModel& cml=INTERP_KERNEL::CellModel::GetCellModel(typel);
          delta+=cm.getNumberOfNodes()-cml.getNumberOfNodes();
        }
    }
  if(delta==0)
    return ;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newConn=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newConnI=DataArrayInt::New();
  newConn->alloc(getMeshLength()-delta,1);
  newConnI->alloc(nbOfCells+1,1);
  const int *icptr=_nodal_connec->getConstPointer();
  const int *iciptr=_nodal_connec_index->getConstPointer();
  int *ocptr=newConn->getPointer();
  int *ociptr=newConnI->getPointer();
  *ociptr=0;
  _types.clear();
  for(int i=0;i<nbOfCells;i++,ociptr++)
    {
      INTERP_KERNEL::NormalizedCellType type=getTypeOfCell(i);
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
          int newNbOfNodes=cml.getNumberOfNodes();
          *ocptr++=(int)typel;
          ocptr=std::copy(icptr+iciptr[i]+1,icptr+iciptr[i]+newNbOfNodes+1,ocptr);
          ociptr[1]=ociptr[0]+newNbOfNodes+1;
        }
    }
  setConnectivity(newConn,newConnI,false);
}

/*!
 * This method tessallates 'this' so that the number of cells remains the same.
 * This method works only for meshes with spaceDim equal to 2 and meshDim equal to 2.
 * If no cells are quadratic in 'this' (INTERP_KERNEL::NORM_QUAD8, INTERP_KERNEL::NORM_TRI6, INTERP_KERNEL::NORM_QPOLYG ) this method will remain unchanged.
 * 
 * \b WARNING this method can lead to a uge amount of nodes if eps is very low.
 * @param eps specifies the maximal angle (in radian) between 2 subedges of polylinized edge constituting the input polygon.
 */
void MEDCouplingUMesh::tessellate2D(double eps) throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  if(getMeshDimension()!=2 || getSpaceDimension()!=2)  
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::tessellate2D works on umeshes with meshdim equal to 2 and spaceDim equal to 2 too!");
  double epsa=fabs(eps);
  if(epsa<std::numeric_limits<double>::min())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::tessellate2DCurve : epsilon is null ! Please specify a higher epsilon. If too tiny it can lead to a huge amount of nodes and memory !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> desc1=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> descIndx1=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> revDesc1=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> revDescIndx1=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mDesc=buildDescendingConnectivity2(desc1,descIndx1,revDesc1,revDescIndx1);
  revDesc1=0; revDescIndx1=0;
  mDesc->tessellate2DCurve(eps);
  subDivide2DMesh(mDesc->_nodal_connec->getConstPointer(),mDesc->_nodal_connec_index->getConstPointer(),desc1->getConstPointer(),descIndx1->getConstPointer());
  setCoords(mDesc->getCoords());
}

/*!
 * This method tessallates 'this' so that the number of cells remains the same.
 * This method works only for meshes with spaceDim equal to 2 and meshDim equal to 1.
 * If no cells are quadratic in 'this' (INTERP_KERNEL::NORM_QUAD8, INTERP_KERNEL::NORM_TRI6, INTERP_KERNEL::NORM_QPOLYG ) this method will remain unchanged.
 * 
 * \b WARNING this method can lead to a uge amount of nodes if eps is very low.
 * @param eps specifies the maximal angle (in radian) between 2 subedges of polylinized edge constituting the input polygon.
 */
void MEDCouplingUMesh::tessellate2DCurve(double eps) throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  if(getMeshDimension()!=1 || getSpaceDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::tessellate2DCurve works on umeshes with meshdim equal to 1 and spaceDim equal to 2 too!");
  double epsa=fabs(eps);
  if(epsa<std::numeric_limits<double>::min())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::tessellate2DCurve : epsilon is null ! Please specify a higher epsilon. If too tiny it can lead to a huge amount of nodes and memory !");
  INTERP_KERNEL::QUADRATIC_PLANAR::_arc_detection_precision=1.e-10;
  int nbCells=getNumberOfCells();
  int nbNodes=getNumberOfNodes();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const double *coords=_coords->getConstPointer();
  std::vector<double> addCoo;
  std::vector<int> newConn;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newConnI=DataArrayInt::New();
  newConnI->alloc(nbCells+1,1);
  int *newConnIPtr=newConnI->getPointer();
  *newConnIPtr=0;
  int tmp1[3];
  INTERP_KERNEL::Node *tmp2[3];
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  for(int i=0;i<nbCells;i++,newConnIPtr++)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[connI[i]]);
      if(cm.isQuadratic())
        {//assert(connI[i+1]-connI[i]-1==3)
          tmp1[0]=conn[connI[i]+1+0]; tmp1[1]=conn[connI[i]+1+1]; tmp1[2]=conn[connI[i]+1+2];
          tmp2[0]=new INTERP_KERNEL::Node(coords[2*tmp1[0]],coords[2*tmp1[0]+1]);
          tmp2[1]=new INTERP_KERNEL::Node(coords[2*tmp1[1]],coords[2*tmp1[1]+1]);
          tmp2[2]=new INTERP_KERNEL::Node(coords[2*tmp1[2]],coords[2*tmp1[2]+1]);
          INTERP_KERNEL::EdgeArcCircle *eac=INTERP_KERNEL::EdgeArcCircle::BuildFromNodes(tmp2[0],tmp2[2],tmp2[1]);
          if(eac)
            {
              eac->tesselate(tmp1,nbNodes,epsa,newConn,addCoo);
              types.insert((INTERP_KERNEL::NormalizedCellType)newConn[newConnIPtr[0]]);
              delete eac;
              newConnIPtr[1]=(int)newConn.size();
            }
          else
            {
              types.insert(INTERP_KERNEL::NORM_SEG2);
              newConn.push_back(INTERP_KERNEL::NORM_SEG2);
              newConn.insert(newConn.end(),conn+connI[i]+1,conn+connI[i]+3);
              newConnIPtr[1]=newConnIPtr[0]+3;
            }
        }
      else
        {
          types.insert((INTERP_KERNEL::NormalizedCellType)conn[connI[i]]);
          newConn.insert(newConn.end(),conn+connI[i],conn+connI[i+1]);
          newConnIPtr[1]=newConnIPtr[0]+3;
        }
    }
  if(addCoo.empty() && ((int)newConn.size())==_nodal_connec->getNumberOfTuples())//nothing happens during tasselation : no update needed
    return ;
  _types=types;
  DataArrayInt::SetArrayIn(newConnI,_nodal_connec_index);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newConnArr=DataArrayInt::New();
  newConnArr->alloc((int)newConn.size(),1);
  std::copy(newConn.begin(),newConn.end(),newConnArr->getPointer());
  DataArrayInt::SetArrayIn(newConnArr,_nodal_connec);
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> newCoords=DataArrayDouble::New();
  newCoords->alloc(nbNodes+((int)addCoo.size())/2,2);
  double *work=std::copy(_coords->begin(),_coords->end(),newCoords->getPointer());
  std::copy(addCoo.begin(),addCoo.end(),work);
  DataArrayDouble::SetArrayIn(newCoords,_coords);
  updateTime();
}

/*!
 * This methods modify this by converting each cells into simplex cell, that is too say triangle for meshdim==2 or tetra for meshdim==3.
 * This cut into simplex is performed following the parameter 'policy'. This method so typically increases the number of cells of this.
 * This method returns new2old array that specifies a each cell of 'this' after the call what was its id it comes.
 * 
 * The semantic of 'policy' parameter :
 * - 1 only QUAD4. For QUAD4 the cut is done along 0-2 diagonal for QUAD4
 * - 2 only QUAD4. For QUAD4 the cut is done along 1-3 diagonal for QUAD4
 */
DataArrayInt *MEDCouplingUMesh::simplexize(int policy) throw(INTERP_KERNEL::Exception)
{
  switch(policy)
    {
    case 0:
      return simplexizePol0();
    case 1:
      return simplexizePol1();
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingUMesh::simplexize : unrecognized policy ! Must be 0 or 1 !");
    }
}

bool MEDCouplingUMesh::areOnlySimplexCells() const throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  if(getMeshDimension()<1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::areOnlySimplexCells : only available with meshes having a meshdim >= 1 !");
  int nbCells=getNumberOfCells();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  for(int i=0;i<nbCells;i++)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[connI[i]]);
      if(!cm.isSimplex())
        return false;
    }
  return true;
}

/*!
 * This method implements policy 0 of virtual method ParaMEDMEM::MEDCouplingUMesh::simplexize.
 */
DataArrayInt *MEDCouplingUMesh::simplexizePol0() throw(INTERP_KERNEL::Exception)
{
  if(getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::simplexizePol0 : this policy is only available for mesh with meshdim == 2 !");
  int nbOfCells=getNumberOfCells();
  DataArrayInt *ret=DataArrayInt::New();
  int nbOfCutCells=getNumberOfCellsWithType(INTERP_KERNEL::NORM_QUAD4);
  ret->alloc(nbOfCells+nbOfCutCells,1);
  if(nbOfCutCells==0)
    {
      ret->iota(0);
      return ret;
    }
  int *retPt=ret->getPointer();
  DataArrayInt *newConn=DataArrayInt::New();
  DataArrayInt *newConnI=DataArrayInt::New();
  newConnI->alloc(nbOfCells+nbOfCutCells+1,1);
  newConn->alloc(getMeshLength()+3*nbOfCutCells,1);
  int *pt=newConn->getPointer();
  int *ptI=newConnI->getPointer();
  ptI[0]=0;
  const int *oldc=_nodal_connec->getConstPointer();
  const int *ci=_nodal_connec_index->getConstPointer();
  for(int i=0;i<nbOfCells;i++,ci++)
    {
      if((INTERP_KERNEL::NormalizedCellType)oldc[ci[0]]==INTERP_KERNEL::NORM_QUAD4)
        {
          const int tmp[8]={(int)INTERP_KERNEL::NORM_TRI3,oldc[ci[0]+1],oldc[ci[0]+2],oldc[ci[0]+3],
                            (int)INTERP_KERNEL::NORM_TRI3,oldc[ci[0]+1],oldc[ci[0]+3],oldc[ci[0]+4]};
          pt=std::copy(tmp,tmp+8,pt);
          ptI[1]=ptI[0]+4;
          ptI[2]=ptI[0]+8;
          *retPt++=i;
          *retPt++=i;
          ptI+=2;
        }
      else
        {
          pt=std::copy(oldc+ci[0],oldc+ci[1],pt);
          ptI[1]=ptI[0]+ci[1]-ci[0];
          ptI++;
          *retPt++=i;
        }
    }
  _nodal_connec->decrRef();
  _nodal_connec=newConn;
  _nodal_connec_index->decrRef();
  _nodal_connec_index=newConnI;
  computeTypes();
  updateTime();
  return ret;
}

/*!
 * This method implements policy 1 of virtual method ParaMEDMEM::MEDCouplingUMesh::simplexize.
 */
DataArrayInt *MEDCouplingUMesh::simplexizePol1() throw(INTERP_KERNEL::Exception)
{
  if(getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::simplexizePol0 : this policy is only available for mesh with meshdim == 2 !");
  int nbOfCells=getNumberOfCells();
  DataArrayInt *ret=DataArrayInt::New();
  int nbOfCutCells=getNumberOfCellsWithType(INTERP_KERNEL::NORM_QUAD4);
  ret->alloc(nbOfCells+nbOfCutCells,1);
  if(nbOfCutCells==0)
    {
      ret->iota(0);
      return ret;
    }
  int *retPt=ret->getPointer();
  DataArrayInt *newConn=DataArrayInt::New();
  DataArrayInt *newConnI=DataArrayInt::New();
  newConnI->alloc(nbOfCells+nbOfCutCells+1,1);
  newConn->alloc(getMeshLength()+3*nbOfCutCells,1);
  int *pt=newConn->getPointer();
  int *ptI=newConnI->getPointer();
  ptI[0]=0;
  const int *oldc=_nodal_connec->getConstPointer();
  const int *ci=_nodal_connec_index->getConstPointer();
  for(int i=0;i<nbOfCells;i++,ci++)
    {
      if((INTERP_KERNEL::NormalizedCellType)oldc[ci[0]]==INTERP_KERNEL::NORM_QUAD4)
        {
          const int tmp[8]={(int)INTERP_KERNEL::NORM_TRI3,oldc[ci[0]+1],oldc[ci[0]+2],oldc[ci[0]+4],
                            (int)INTERP_KERNEL::NORM_TRI3,oldc[ci[0]+2],oldc[ci[0]+3],oldc[ci[0]+4]};
          pt=std::copy(tmp,tmp+8,pt);
          ptI[1]=ptI[0]+4;
          ptI[2]=ptI[0]+8;
          *retPt++=i;
          *retPt++=i;
          ptI+=2;
        }
      else
        {
          pt=std::copy(oldc+ci[0],oldc+ci[1],pt);
          ptI[1]=ptI[0]+ci[1]-ci[0];
          ptI++;
          *retPt++=i;
        }
    }
  _nodal_connec->decrRef();
  _nodal_connec=newConn;
  _nodal_connec_index->decrRef();
  _nodal_connec_index=newConnI;
  computeTypes();
  updateTime();
  return ret;
}

/*!
 * This private method is used to subdivide edges of a mesh with meshdim==2. If 'this' has no a meshdim equal to 2 an exception will be thrown.
 * This method completly ignore coordinates.
 * @param nodeSubdived is the nodal connectivity of subdivision of edges
 * @param nodeIndxSubdived is the nodal connectivity index of subdivision of edges
 * @param desc is descending connectivity in format specified in MEDCouplingUMesh::buildDescendingConnectivity2
 * @param descIndex is descending connectivity index in format specified in MEDCouplingUMesh::buildDescendingConnectivity2
 */
void MEDCouplingUMesh::subDivide2DMesh(const int *nodeSubdived, const int *nodeIndxSubdived, const int *desc, const int *descIndex) throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  if(getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::subDivide2DMesh : works only on umesh with meshdim==2 !");
  int nbOfCells=getNumberOfCells();
  int *connI=_nodal_connec_index->getPointer();
  int newConnLgth=0;
  for(int i=0;i<nbOfCells;i++,connI++)
    {
      int offset=descIndex[i];
      int nbOfEdges=descIndex[i+1]-offset;
      //
      bool ddirect=desc[offset+nbOfEdges-1]>0;
      int eedgeId=std::abs(desc[offset+nbOfEdges-1])-1;
      int ref=ddirect?nodeSubdived[nodeIndxSubdived[eedgeId+1]-1]:nodeSubdived[nodeIndxSubdived[eedgeId]+1];
      for(int j=0;j<nbOfEdges;j++)
        {
          bool direct=desc[offset+j]>0;
          int edgeId=std::abs(desc[offset+j])-1;
          if(!INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)nodeSubdived[nodeIndxSubdived[edgeId]]).isQuadratic())
            {
              int id1=nodeSubdived[nodeIndxSubdived[edgeId]+1];
              int id2=nodeSubdived[nodeIndxSubdived[edgeId+1]-1];
              int ref2=direct?id1:id2;
              if(ref==ref2)
                {
                  int nbOfSubNodes=nodeIndxSubdived[edgeId+1]-nodeIndxSubdived[edgeId]-1;
                  newConnLgth+=nbOfSubNodes-1;
                  ref=direct?id2:id1;
                }
              else
                {
                  std::ostringstream oss; oss << "MEDCouplingUMesh::subDivide2DMesh : On polygon #" << i << " edgeid #" << j << " subedges mismatch : end subedge k!=start subedge k+1 !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
            }
          else
            {
              throw INTERP_KERNEL::Exception("MEDCouplingUMesh::subDivide2DMesh : this method only subdivides into linear edges !");
            }
        }
      newConnLgth++;//+1 is for cell type
      connI[1]=newConnLgth;
    }
  //
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newConn=DataArrayInt::New();
  newConn->alloc(newConnLgth,1);
  int *work=newConn->getPointer();
  for(int i=0;i<nbOfCells;i++)
    {
      *work++=INTERP_KERNEL::NORM_POLYGON;
      int offset=descIndex[i];
      int nbOfEdges=descIndex[i+1]-offset;
      for(int j=0;j<nbOfEdges;j++)
        {
          bool direct=desc[offset+j]>0;
          int edgeId=std::abs(desc[offset+j])-1;
          if(direct)
            work=std::copy(nodeSubdived+nodeIndxSubdived[edgeId]+1,nodeSubdived+nodeIndxSubdived[edgeId+1]-1,work);
          else
            {
              int nbOfSubNodes=nodeIndxSubdived[edgeId+1]-nodeIndxSubdived[edgeId]-1;
              std::reverse_iterator<const int *> it(nodeSubdived+nodeIndxSubdived[edgeId+1]);
              work=std::copy(it,it+nbOfSubNodes-1,work);
            }
        }
    }
  DataArrayInt::SetArrayIn(newConn,_nodal_connec);
  _types.clear();
  if(nbOfCells>0)
    _types.insert(INTERP_KERNEL::NORM_POLYGON);
}

/*!
 * This method converts all degenerated cells to simpler cells. For example a NORM_QUAD4 cell consituted from 2 same node id in its
 * nodal connectivity will be transform to a NORM_TRI3 cell.
 * This method works \b only \b on \b linear cells.
 * This method works on nodes ids, that is to say a call to ParaMEDMEM::MEDCouplingUMesh::mergeNodes
 * method could be usefull before calling this method in case of presence of several pair of nodes located on same position.
 * This method throws an exception if 'this' is not fully defined (connectivity).
 * This method throws an exception too if a "too" degenerated cell is detected. For example a NORM_TRI3 with 3 times the same node id.
 */
void MEDCouplingUMesh::convertDegeneratedCells() throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  if(getMeshDimension()<=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::convertDegeneratedCells works on umeshes with meshdim equals to 2 or 3 !");
  int nbOfCells=getNumberOfCells();
  if(nbOfCells<1)
    return ;
  int initMeshLgth=getMeshLength();
  int *conn=_nodal_connec->getPointer();
  int *index=_nodal_connec_index->getPointer();
  int posOfCurCell=0;
  int newPos=0;
  int lgthOfCurCell;
  for(int i=0;i<nbOfCells;i++)
    {
      lgthOfCurCell=index[i+1]-posOfCurCell;
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)conn[posOfCurCell];
      int newLgth;
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
 * This method checks that all or only polygons (depending 'polyOnly' parameter) 2D cells are correctly oriented relative to 'vec' vector.
 * The 'vec' vector has to have a non nul norm.
 * If not 'cells' parameter will be appended with cellIds of incorrect cells.
 * @throw when 'this' is not a mesh with meshdim==2 and spacedim==3
 */
void MEDCouplingUMesh::are2DCellsNotCorrectlyOriented(const double *vec, bool polyOnly, std::vector<int>& cells) const throw(INTERP_KERNEL::Exception)
{
  if(getMeshDimension()!=2 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply are2DCellsNotCorrectlyOriented on it : must be meshDim==2 and spaceDim==3 !");
  int nbOfCells=getNumberOfCells();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const double *coordsPtr=_coords->getConstPointer();
  for(int i=0;i<nbOfCells;i++)
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
 * This method orient correctly (if needed) all or only polygons (depending 'polyOnly' parameter)  2D cells are correctly oriented relative to 'vec' vector.
 * The 'vec' vector has to have a non nul norm.
 * @throw when 'this' is not a mesh with meshdim==2 and spacedim==3
 */
void MEDCouplingUMesh::orientCorrectly2DCells(const double *vec, bool polyOnly) throw(INTERP_KERNEL::Exception)
{
  if(getMeshDimension()!=2 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply orientCorrectly2DCells on it : must be meshDim==2 and spaceDim==3 !");
  int nbOfCells=getNumberOfCells();
  int *conn=_nodal_connec->getPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const double *coordsPtr=_coords->getConstPointer();
  bool isModified=false;
  for(int i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)conn[connI[i]];
      if(!polyOnly || (type==INTERP_KERNEL::NORM_POLYGON || type==INTERP_KERNEL::NORM_QPOLYG))
        {
          bool isQuadratic=INTERP_KERNEL::CellModel::GetCellModel(type).isQuadratic();
          if(!IsPolygonWellOriented(isQuadratic,vec,conn+connI[i]+1,conn+connI[i+1],coordsPtr))
            {
              isModified=true;
              std::vector<int> tmp(connI[i+1]-connI[i]-2);
              std::copy(conn+connI[i]+2,conn+connI[i+1],tmp.rbegin());
              std::copy(tmp.begin(),tmp.end(),conn+connI[i]+2);
            }
        }
    }
  if(isModified)
    _nodal_connec->declareAsNew();
  updateTime();
}

/*!
 * This method checks that all polyhedrons cells have correctly oriented faces.
 * If not, 'cells' parameter will be appended with cellIds of incorrect cells.
 * @throw when 'this' is not a mesh with meshdim==3 and spacedim==3
 */
void MEDCouplingUMesh::arePolyhedronsNotCorrectlyOriented(std::vector<int>& cells) const throw(INTERP_KERNEL::Exception)
{
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply arePolyhedronsNotCorrectlyOriented on it : must be meshDim==3 and spaceDim==3 !");
  int nbOfCells=getNumberOfCells();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const double *coordsPtr=_coords->getConstPointer();
  for(int i=0;i<nbOfCells;i++)
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
 * This method tries to orient correctly polhedrons cells.
 * @throw when 'this' is not a mesh with meshdim==3 and spacedim==3. An exception is also thrown when the attempt of reparation fails.
 */
void MEDCouplingUMesh::orientCorrectlyPolyhedrons() throw(INTERP_KERNEL::Exception)
{
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply orientCorrectlyPolyhedrons on it : must be meshDim==3 and spaceDim==3 !");
  int nbOfCells=getNumberOfCells();
  int *conn=_nodal_connec->getPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const double *coordsPtr=_coords->getConstPointer();
  bool isModified=false;
  for(int i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)conn[connI[i]];
      if(type==INTERP_KERNEL::NORM_POLYHED)
        if(!IsPolyhedronWellOriented(conn+connI[i]+1,conn+connI[i+1],coordsPtr))
          {
            TryToCorrectPolyhedronOrientation(conn+connI[i]+1,conn+connI[i+1],coordsPtr);
            isModified=true;
          }
    }
  if(isModified)
    _nodal_connec->declareAsNew();
  updateTime();
}

/*!
 * This method has a sense for meshes with spaceDim==3 and meshDim==2.
 * If it is not the case an exception will be thrown.
 * This method is fast because the first cell of 'this' is used to compute the plane.
 * @param vec output of size at least 3 used to store the normal vector (with norm equal to Area ) of searched plane.
 * @param pos output of size at least 3 used to store a point owned of searched plane.
 */
void MEDCouplingUMesh::getFastAveragePlaneOfThis(double *vec, double *pos) const throw(INTERP_KERNEL::Exception)
{
  if(getMeshDimension()!=2 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply getFastAveragePlaneOfThis on it : must be meshDim==2 and spaceDim==3 !");
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const double *coordsPtr=_coords->getConstPointer();
  INTERP_KERNEL::areaVectorOfPolygon<int,INTERP_KERNEL::ALL_C_MODE>(conn+1,connI[1]-connI[0]-1,coordsPtr,vec);
  std::copy(coordsPtr+3*conn[1],coordsPtr+3*conn[1]+3,pos);
}

/*!
 * The returned newly created field has to be managed by the caller.
 * This method returns a field on cell with no time lying on 'this'. The meshdimension and spacedimension of this are expected to be both in [2,3]. If not an exception will be thrown.
 * This method for the moment only deals with NORM_TRI3, NORM_QUAD4 and NORM_TETRA4 geometric types.
 * If a cell has an another type an exception will be thrown.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::getEdgeRatioField() const throw(INTERP_KERNEL::Exception)
{
  checkCoherency();
  int spaceDim=getSpaceDimension();
  int meshDim=getMeshDimension();
  if(spaceDim!=2 && spaceDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getEdgeRatioField : SpaceDimension must be equal to 2 or 3 !");
  if(meshDim!=2 && meshDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getEdgeRatioField : MeshDimension must be equal to 2 or 3 !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  ret->setMesh(this);
  int nbOfCells=getNumberOfCells();
  DataArrayDouble *arr=DataArrayDouble::New();
  arr->alloc(nbOfCells,1);
  double *pt=arr->getPointer();
  ret->setArray(arr);//In case of throw to avoid mem leaks arr will be used after decrRef.
  arr->decrRef();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const double *coo=_coords->getConstPointer();
  double tmp[12];
  for(int i=0;i<nbOfCells;i++,pt++)
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
  ret->incrRef();
  return ret;
}

/*!
 * The returned newly created field has to be managed by the caller.
 * This method returns a field on cell with no time lying on 'this'. The meshdimension and spacedimension of this are expected to be both in [2,3]. If not an exception will be thrown.
 * This method for the moment only deals with NORM_TRI3, NORM_QUAD4 and NORM_TETRA4 geometric types.
 * If a cell has an another type an exception will be thrown.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::getAspectRatioField() const throw(INTERP_KERNEL::Exception)
{
  checkCoherency();
  int spaceDim=getSpaceDimension();
  int meshDim=getMeshDimension();
  if(spaceDim!=2 && spaceDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getAspectRatioField : SpaceDimension must be equal to 2 or 3 !");
  if(meshDim!=2 && meshDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getAspectRatioField : MeshDimension must be equal to 2 or 3 !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  ret->setMesh(this);
  int nbOfCells=getNumberOfCells();
  DataArrayDouble *arr=DataArrayDouble::New();
  arr->alloc(nbOfCells,1);
  double *pt=arr->getPointer();
  ret->setArray(arr);//In case of throw to avoid mem leaks arr will be used after decrRef.
  arr->decrRef();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const double *coo=_coords->getConstPointer();
  double tmp[12];
  for(int i=0;i<nbOfCells;i++,pt++)
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
  ret->incrRef();
  return ret;
}

/*!
 * The returned newly created field has to be managed by the caller.
 * This method returns a field on cell with no time lying on 'this'. The meshdimension must be equal to 2 and the spacedimension must be equal to 3. If not an exception will be thrown.
 * This method for the moment only deals with NORM_QUAD4 geometric type.
 * If a cell has an another type an exception will be thrown.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::getWarpField() const throw(INTERP_KERNEL::Exception)
{
  checkCoherency();
  int spaceDim=getSpaceDimension();
  int meshDim=getMeshDimension();
  if(spaceDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getWarpField : SpaceDimension must be equal to 3 !");
  if(meshDim!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getWarpField : MeshDimension must be equal to 2 !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  ret->setMesh(this);
  int nbOfCells=getNumberOfCells();
  DataArrayDouble *arr=DataArrayDouble::New();
  arr->alloc(nbOfCells,1);
  double *pt=arr->getPointer();
  ret->setArray(arr);//In case of throw to avoid mem leaks arr will be used after decrRef.
  arr->decrRef();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const double *coo=_coords->getConstPointer();
  double tmp[12];
  for(int i=0;i<nbOfCells;i++,pt++)
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
  ret->incrRef();
  return ret;
}

/*!
 * The returned newly created field has to be managed by the caller.
 * This method returns a field on cell with no time lying on 'this'. The meshdimension must be equal to 2 and the spacedimension must be equal to 3. If not an exception will be thrown.
 * This method for the moment only deals with NORM_QUAD4 geometric type.
 * If a cell has an another type an exception will be thrown.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::getSkewField() const throw(INTERP_KERNEL::Exception)
{
  checkCoherency();
  int spaceDim=getSpaceDimension();
  int meshDim=getMeshDimension();
  if(spaceDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getSkewField : SpaceDimension must be equal to 3 !");
  if(meshDim!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getSkewField : MeshDimension must be equal to 2 !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  ret->setMesh(this);
  int nbOfCells=getNumberOfCells();
  DataArrayDouble *arr=DataArrayDouble::New();
  arr->alloc(nbOfCells,1);
  double *pt=arr->getPointer();
  ret->setArray(arr);//In case of throw to avoid mem leaks arr will be used after decrRef.
  arr->decrRef();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const double *coo=_coords->getConstPointer();
  double tmp[12];
  for(int i=0;i<nbOfCells;i++,pt++)
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
  ret->incrRef();
  return ret;
}

/*!
 * This method aggregate the bbox of each cell and put it into bbox parameter.
 * @param bbox out parameter of size 2*spacedim*nbOfcells.
 */
void MEDCouplingUMesh::getBoundingBoxForBBTree(std::vector<double>& bbox) const
{
  int spaceDim=getSpaceDimension();
  int nbOfCells=getNumberOfCells();
  bbox.resize(2*nbOfCells*spaceDim);
  for(int i=0;i<nbOfCells*spaceDim;i++)
    {
      bbox[2*i]=std::numeric_limits<double>::max();
      bbox[2*i+1]=-std::numeric_limits<double>::max();
    }
  const double *coordsPtr=_coords->getConstPointer();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  for(int i=0;i<nbOfCells;i++)
    {
      int offset=connI[i]+1;
      int nbOfNodesForCell=connI[i+1]-offset;
      for(int j=0;j<nbOfNodesForCell;j++)
        {
          int nodeId=conn[offset+j];
          if(nodeId>=0)
            for(int k=0;k<spaceDim;k++)
              {
                bbox[2*spaceDim*i+2*k]=std::min(bbox[2*spaceDim*i+2*k],coordsPtr[spaceDim*nodeId+k]);
                bbox[2*spaceDim*i+2*k+1]=std::max(bbox[2*spaceDim*i+2*k+1],coordsPtr[spaceDim*nodeId+k]);
              }
        }
    }
}

/// @cond INTERNAL

namespace ParaMEDMEMImpl
{
  class ConnReader
  {
  public:
    ConnReader(const int *c, int val):_conn(c),_val(val) { }
    bool operator() (const int& pos) { return _conn[pos]!=_val; }
  private:
    const int *_conn;
    int _val;
  };

  class ConnReader2
  {
  public:
    ConnReader2(const int *c, int val):_conn(c),_val(val) { }
    bool operator() (const int& pos) { return _conn[pos]==_val; }
  private:
    const int *_conn;
    int _val;
  };
}

/// @endcond

/*!
 * This method expects that 'this' is sorted by types. If not an exception will be thrown.
 * This method returns in the same format as code (see MEDCouplingUMesh::checkTypeConsistencyAndContig or MEDCouplingUMesh::splitProfilePerType) how
 * 'this' is composed in cell types.
 * The returned array is of size 3*n where n is the number of different types present in 'this'. 
 * For every k in [0,n] ret[3*k+2]==0 because it has no sense here. 
 * This parameter is kept only for compatibility with other methode listed above.
 */
std::vector<int> MEDCouplingUMesh::getDistributionOfTypes() const throw(INTERP_KERNEL::Exception)
{
  checkConnectivityFullyDefined();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const int *work=connI;
  int nbOfCells=getNumberOfCells();
  std::size_t n=getAllTypes().size();
  std::vector<int> ret(3*n,0); //ret[3*k+2]==0 because it has no sense here
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  for(std::size_t i=0;work!=connI+nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType)conn[*work];
      if(types.find(typ)!=types.end())
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::getDistributionOfTypes : Type " << INTERP_KERNEL::CellModel::GetCellModel(typ).getRepr();
          oss << " is not contiguous !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      types.insert(typ);
      ret[3*i]=typ;
      const int *work2=std::find_if(work+1,connI+nbOfCells,ParaMEDMEMImpl::ConnReader(conn,typ));
      ret[3*i+1]=(int)std::distance(work,work2);
      work=work2;
    }
  return ret;
}

/*!
 * This method is used to check that this has contiguous cell type in same order than described in 'code'.
 * only for types cell, type node is not managed.
 * Format of 'code' is the following. 'code' should be of size 3*n and non empty. If not an exception is thrown.
 * foreach k in [0,n) on 3*k pos represent the geometric type and 3*k+1 number of elements of type 3*k.
 * 3*k+2 refers if different from -1 the pos in 'idsPerType' to get the corresponding array.
 * If 2 or more same geometric type is in 'code' and exception is thrown too.
 *
 * This method firstly checks
 * If it exists k so that 3*k geometric type is not in geometric types of this an exception will be thrown.
 * If it exists k so that 3*k geometric type exists but the number of consecutive cell types does not match,
 * an exception is thrown too.
 * 
 * If all geometric types in 'code' are exactly those in 'this' null pointer is returned.
 * If it exists a geometric type in 'this' \b not in 'code' \b no exception is thrown 
 * and a DataArrayInt instance is returned that the user has the responsability to deallocate.
 */
DataArrayInt *MEDCouplingUMesh::checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception)
{
  if(code.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::checkTypeConsistencyAndContig : code is empty, should not !");
  std::size_t sz=code.size();
  std::size_t n=sz/3;
  if(sz%3!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::checkTypeConsistencyAndContig : code size is NOT %3 !");
  std::vector<INTERP_KERNEL::NormalizedCellType> types;
  int nb=0;
  for(std::size_t i=0;i<n;i++)
    if(std::find(types.begin(),types.end(),(INTERP_KERNEL::NormalizedCellType)code[3*i])==types.end())
      {
        types.push_back((INTERP_KERNEL::NormalizedCellType)code[3*i]);
        nb+=code[3*i+1];
        if(_types.find((INTERP_KERNEL::NormalizedCellType)code[3*i])==_types.end())
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::checkTypeConsistencyAndContig : expected geo types not in this !");
      }
  if(types.size()!=n)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::checkTypeConsistencyAndContig : code contains duplication of types in unstructured mesh !");
  if(idsPerType.empty())
    {
      if(!checkConsecutiveCellTypesAndOrder(&types[0],&types[0]+types.size()))
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::checkTypeConsistencyAndContig : non contiguous type !");
      if(types.size()==_types.size())
        return 0;
    }
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nb,1);
  int *retPtr=ret->getPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const int *conn=_nodal_connec->getConstPointer();
  int nbOfCells=getNumberOfCells();
  const int *i=connI;
  int kk=0;
  for(std::vector<INTERP_KERNEL::NormalizedCellType>::const_iterator it=types.begin();it!=types.end();it++,kk++)
    {
      i=std::find_if(i,connI+nbOfCells,ParaMEDMEMImpl::ConnReader2(conn,(int)(*it)));
      int offset=(int)std::distance(connI,i);
      if(code[3*kk+2]==-1)
        {
          const int *j=std::find_if(i+1,connI+nbOfCells,ParaMEDMEMImpl::ConnReader(conn,(int)(*it)));
          std::size_t pos2=std::distance(i,j);
          for(std::size_t k=0;k<pos2;k++)
            *retPtr++=(int)k+offset;
          i=j;
        }
      else
        {
          retPtr=std::transform(idsPerType[code[3*kk+2]]->getConstPointer(),idsPerType[code[3*kk+2]]->getConstPointer()+idsPerType[code[3*kk+2]]->getNbOfElems(),
                                retPtr,std::bind2nd(std::plus<int>(),offset));
        }
    }
  return ret;
}

/*!
 * This method makes the hypothesis that 'this' is sorted by type. If not an exception will be thrown.
 * This method is the opposite of MEDCouplingUMesh::checkTypeConsistencyAndContig method. Given a list of cells in 'profile' it returns a list of profiles sorted by geo type.
 * This method has 1 input 'profile' and 2 outputs 'code' and 'idsPerType'.
 * @throw if 'profile' has not exactly one component. It throws too, if 'profile' contains some values not in [0,getNumberOfCells()) or if 'this' is not fully defined
 */
void MEDCouplingUMesh::splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception)
{
  if(profile->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::splitProfilePerType : input profile should have exactly one component !");
  checkConnectivityFullyDefined();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  int nbOfCells=getNumberOfCells();
  std::vector<INTERP_KERNEL::NormalizedCellType> types;
  std::vector<int> typeRangeVals(1);
  for(const int *i=connI;i!=connI+nbOfCells;)
    {
      INTERP_KERNEL::NormalizedCellType curType=(INTERP_KERNEL::NormalizedCellType)conn[*i];
      if(std::find(types.begin(),types.end(),curType)!=types.end())
        {
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::splitProfilePerType : current mesh is not sorted by type !");
        }
      types.push_back(curType);
      i=std::find_if(i+1,connI+nbOfCells,ParaMEDMEMImpl::ConnReader(conn,(int)curType));
      typeRangeVals.push_back((int)std::distance(connI,i));
    }
  //
  DataArrayInt *castArr=0,*rankInsideCast=0,*castsPresent=0;
  profile->splitByValueRange(&typeRangeVals[0],&typeRangeVals[0]+typeRangeVals.size(),castArr,rankInsideCast,castsPresent);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp0=castArr;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp1=rankInsideCast;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp2=castsPresent;
  //
  int nbOfCastsFinal=castsPresent->getNumberOfTuples();
  code.resize(3*nbOfCastsFinal);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > idsInPflPerType2;
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > idsPerType2;
  for(int i=0;i<nbOfCastsFinal;i++)
    {
      int castId=castsPresent->getIJ(i,0);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp3=castArr->getIdsEqual(castId);
      idsInPflPerType2.push_back(tmp3);
      code[3*i]=(int)types[castId];
      code[3*i+1]=tmp3->getNumberOfTuples();
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp4=rankInsideCast->selectByTupleId(tmp3->getConstPointer(),tmp3->getConstPointer()+tmp3->getNumberOfTuples());
      if(tmp4->getNumberOfTuples()!=typeRangeVals[castId+1]-typeRangeVals[castId] || !tmp4->isIdentity())
        {
          tmp4->copyStringInfoFrom(*profile);
          idsPerType2.push_back(tmp4);
          code[3*i+2]=(int)idsPerType2.size()-1;
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
      DataArrayInt *locDa=idsInPflPerType2[i];
      locDa->incrRef();
      idsInPflPerType[i]=locDa;
    }
  std::size_t sz=idsPerType2.size();
  idsPerType.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      DataArrayInt *locDa=idsPerType2[i];
      locDa->incrRef();
      idsPerType[i]=locDa;
    }
}

/*!
 * This method is here too emulate the MEDMEM behaviour on BDC (buildDescendingConnectivity). Hoping this method becomes deprecated very soon.
 * This method make the assumption that 'this' and 'nM1LevMesh' mesh lyies on same coords (same pointer) as MED and MEDMEM does.
 * The following equality should be verified 'nM1LevMesh->getMeshDimension()==this->getMeshDimension()-1'
 * This method returns 5+2 elements. 'desc', 'descIndx', 'revDesc', 'revDescIndx' and 'meshnM1' behaves exactly as ParaMEDMEM::MEDCouplingUMesh::buildDescendingConnectivity except the content as described after. The returned array specifies the n-1 mesh reordered by type as MEDMEM does. 'nM1LevMeshIds' contains the ids in returned 'meshnM1'. Finally 'meshnM1Old2New' contains numbering old2new that is to say the cell #k in coarse 'nM1LevMesh' will have the number ret[k] in returned mesh 'nM1LevMesh' MEDMEM reordered.
 */
MEDCouplingUMesh *MEDCouplingUMesh::emulateMEDMEMBDC(const MEDCouplingUMesh *nM1LevMesh, DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *&revDesc, DataArrayInt *&revDescIndx, DataArrayInt *& nM1LevMeshIds, DataArrayInt *&meshnM1Old2New) const throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  nM1LevMesh->checkFullyDefined();
  if(getMeshDimension()-1!=nM1LevMesh->getMeshDimension())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::emulateMEDMEMBDC : The mesh passed as first argument should have a meshDim equal to this->getMeshDimension()-1 !" );
  if(_coords!=nM1LevMesh->getCoords())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::emulateMEDMEMBDC : 'this' and mesh in first argument should share the same coords : Use tryToShareSameCoords method !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp0=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp1=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> ret1=buildDescendingConnectivity(desc,descIndx,tmp0,tmp1);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret0=ret1->sortCellsInMEDFileFrmt();
  desc->transformWithIndArr(ret0->getConstPointer(),ret0->getConstPointer()+ret0->getNbOfElems());
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> tmp=MEDCouplingUMesh::New();
  tmp->setConnectivity(tmp0,tmp1);
  tmp->renumberCells(ret0->getConstPointer(),false);
  revDesc=tmp->getNodalConnectivity();
  revDescIndx=tmp->getNodalConnectivityIndex();
  DataArrayInt *ret=0;
  if(!ret1->areCellsIncludedIn(nM1LevMesh,2,ret))
    {
      int tmp2;
      ret->getMaxValue(tmp2);
      ret->decrRef();
      std::ostringstream oss; oss << "MEDCouplingUMesh::emulateMEDMEMBDC : input N-1 mesh present a cell not in descending mesh ... Id of cell is " << tmp2 << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
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
 * This method sorts cell in this so that cells are sorted by cell type specified by MEDMEM and so for MED file.
 * It avoids to deal with renum in MEDLoader so it is usefull for MED file R/W with multi types.
 * This method returns a newly allocated array old2New.
 * This method expects that connectivity of this is set. If not an exception will be thrown. Coordinates are not taken into account.
 */
DataArrayInt *MEDCouplingUMesh::sortCellsInMEDFileFrmt() throw(INTERP_KERNEL::Exception)
{
  checkConnectivityFullyDefined();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=getRenumArrForConsecutiveCellTypesSpec(MEDMEM_ORDER,MEDMEM_ORDER+N_MEDMEM_ORDER);
  renumberCells(ret->getConstPointer(),false);
  ret->incrRef();
  return ret;
}

/*!
 * This methods checks that cells are sorted by their types.
 * This method makes asumption (no check) that connectivity is correctly set before calling.
 */
bool MEDCouplingUMesh::checkConsecutiveCellTypes() const
{
  checkFullyDefined();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  int nbOfCells=getNumberOfCells();
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  for(const int *i=connI;i!=connI+nbOfCells;)
    {
      INTERP_KERNEL::NormalizedCellType curType=(INTERP_KERNEL::NormalizedCellType)conn[*i];
      if(types.find(curType)!=types.end())
        return false;
      types.insert(curType);
      i=std::find_if(i+1,connI+nbOfCells,ParaMEDMEMImpl::ConnReader(conn,(int)curType));
    }
  return true;
}

/*!
 * This method performs the same job as checkConsecutiveCellTypes except that the order of types sequence is analyzed to check
 * that the order is specified in array defined by [orderBg,orderEnd). 
 */
bool MEDCouplingUMesh::checkConsecutiveCellTypesAndOrder(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd) const
{
  checkFullyDefined();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  int nbOfCells=getNumberOfCells();
  int lastPos=-1;
  for(const int *i=connI;i!=connI+nbOfCells;)
    {
      INTERP_KERNEL::NormalizedCellType curType=(INTERP_KERNEL::NormalizedCellType)conn[*i];
      int pos=(int)std::distance(orderBg,std::find(orderBg,orderEnd,curType));
      if(pos<=lastPos)
        return false;
      lastPos=pos;
      i=std::find_if(i+1,connI+nbOfCells,ParaMEDMEMImpl::ConnReader(conn,(int)curType));
    }
  return true;
}

/*!
 * This method returns 2 newly allocated DataArrayInt instances. The first is an array of size 'this->getNumberOfCells()' with one component,
 * that tells for each cell the pos of its type in the array on type given in input parameter. The 2nd output parameter is an array with the same
 * number of tuples than input type array and with one component. This 2nd output array gives type by type the number of occurence of type in 'this'.
 */
DataArrayInt *MEDCouplingUMesh::getLevArrPerCellTypes(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd, DataArrayInt *&nbPerType) const throw(INTERP_KERNEL::Exception)
{
  checkConnectivityFullyDefined();
  int nbOfCells=getNumberOfCells();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmpa=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmpb=DataArrayInt::New();
  tmpa->alloc(nbOfCells,1);
  tmpb->alloc((int)std::distance(orderBg,orderEnd),1);
  tmpb->fillWithZero();
  int *tmp=tmpa->getPointer();
  int *tmp2=tmpb->getPointer();
  for(const int *i=connI;i!=connI+nbOfCells;i++)
    {
      const INTERP_KERNEL::NormalizedCellType *where=std::find(orderBg,orderEnd,(INTERP_KERNEL::NormalizedCellType)conn[*i]);
      if(where!=orderEnd)
        {
          int pos=(int)std::distance(orderBg,where);
          tmp2[pos]++;
          tmp[std::distance(connI,i)]=pos;
        }
      else
        {
          const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[*i]);
          std::ostringstream oss; oss << "MEDCouplingUMesh::getLevArrPerCellTypes : Cell #" << std::distance(connI,i);
          oss << " has a type " << cm.getRepr() << " not in input array of type !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  nbPerType=tmpb;
  tmpa->incrRef();
  tmpb->incrRef();
  return tmpa;
}

/*!
 * This method is similar to method MEDCouplingUMesh::rearrange2ConsecutiveCellTypes except that the type order is specfied by [orderBg,orderEnd) (as MEDCouplingUMesh::checkConsecutiveCellTypesAndOrder method) and that this method is \b const and performs \b NO permutation in 'this'.
 * This method returns an array of size getNumberOfCells() that gives a renumber array old2New that can be used as input of MEDCouplingMesh::renumberCells.
 * The mesh after this call to MEDCouplingMesh::renumberCells will pass the test of MEDCouplingUMesh::checkConsecutiveCellTypesAndOrder with the same inputs.
 * The returned array minimizes the permutations that is to say the order of cells inside same geometric type remains the same.
 */
DataArrayInt *MEDCouplingUMesh::getRenumArrForConsecutiveCellTypesSpec(const INTERP_KERNEL::NormalizedCellType *orderBg, const INTERP_KERNEL::NormalizedCellType *orderEnd) const throw(INTERP_KERNEL::Exception)
{
  DataArrayInt *nbPerType=0;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmpa=getLevArrPerCellTypes(orderBg,orderEnd,nbPerType);
  nbPerType->decrRef();
  return tmpa->buildPermArrPerLevel();
}

/*!
 * This method reorganize the cells of 'this' so that the cells with same geometric types are put together.
 * The number of cells remains unchanged after the call of this method.
 * This method tries to minimizes the number of needed permutations. So, this method behaves not exactly as
 * MEDCouplingUMesh::sortCellsInMEDFileFrmt.
 *
 * @return the array giving the correspondance old to new.
 */
DataArrayInt *MEDCouplingUMesh::rearrange2ConsecutiveCellTypes()
{
  checkFullyDefined();
  computeTypes();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  int nbOfCells=getNumberOfCells();
  std::vector<INTERP_KERNEL::NormalizedCellType> types;
  for(const int *i=connI;i!=connI+nbOfCells && (types.size()!=_types.size());)
    if(std::find(types.begin(),types.end(),(INTERP_KERNEL::NormalizedCellType)conn[*i])==types.end())
      {
        INTERP_KERNEL::NormalizedCellType curType=(INTERP_KERNEL::NormalizedCellType)conn[*i];
        types.push_back(curType);
        for(i++;i!=connI+nbOfCells && (INTERP_KERNEL::NormalizedCellType)conn[*i]==curType;i++);
      }
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbOfCells,1);
  int *retPtr=ret->getPointer();
  std::fill(retPtr,retPtr+nbOfCells,-1);
  int newCellId=0;
  for(std::vector<INTERP_KERNEL::NormalizedCellType>::const_iterator iter=types.begin();iter!=types.end();iter++)
    {
      for(const int *i=connI;i!=connI+nbOfCells;i++)
        if((INTERP_KERNEL::NormalizedCellType)conn[*i]==(*iter))
          retPtr[std::distance(connI,i)]=newCellId++;
    }
  renumberCells(retPtr,false);
  return ret;
}

/*!
 * This method splits 'this' into as mush as untructured meshes that consecutive set of same type cells.
 * So this method has typically a sense if MEDCouplingUMesh::checkConsecutiveCellTypes has a sense.
 * This method makes asumption that connectivity is correctly set before calling.
 */
std::vector<MEDCouplingUMesh *> MEDCouplingUMesh::splitByType() const
{
  checkFullyDefined();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  int nbOfCells=getNumberOfCells();
  std::vector<MEDCouplingUMesh *> ret;
  for(const int *i=connI;i!=connI+nbOfCells;)
    {
      INTERP_KERNEL::NormalizedCellType curType=(INTERP_KERNEL::NormalizedCellType)conn[*i];
      int beginCellId=(int)std::distance(connI,i);
      i=std::find_if(i+1,connI+nbOfCells,ParaMEDMEMImpl::ConnReader(conn,(int)curType));
      int endCellId=(int)std::distance(connI,i);
      int sz=endCellId-beginCellId;
      int *cells=new int[sz];
      for(int j=0;j<sz;j++)
        cells[j]=beginCellId+j;
      MEDCouplingUMesh *m=(MEDCouplingUMesh *)buildPartOfMySelf(cells,cells+sz,true);
      delete [] cells;
      ret.push_back(m);
    }
  return ret;
}

/*!
 * This method takes in input a vector of MEDCouplingUMesh instances lying on the same coordinates with same mesh dimensions.
 * Each mesh in \b ms must be sorted by type with the same order (typically using MEDCouplingUMesh::sortCellsInMEDFileFrmt).
 * This method is particulary useful for MED file interaction. It allows to aggregate several meshes and keeping the type sorting
 * and the track of the permutation by chunk of same geotype cells to retrieve it. The traditional formats old2new and new2old
 * are not used here to avoid the build of big permutation array.
 *
 * \param [in] ms meshes with same mesh dimension lying on the same coords and sorted by type following de the same geometric type order than
 *                those specified in MEDCouplingUMesh::sortCellsInMEDFileFrmt method.
 * \param [out] szOfCellGrpOfSameType is a newly allocated DataArrayInt instance whose number of tuples is equal to the number of chunks of same geotype
 *              in all meshes in \b ms. The accumulation of all values of this array is equal to the number of cells of returned mesh.
 * \param [out] idInMsOfCellGrpOfSameType is a newly allocated DataArrayInt instance having the same size than \b szOfCellGrpOfSameType. This
 *              output array gives for each chunck of same type the corresponding mesh id in \b ms.
 * \return A newly allocated unstructured mesh that is the result of the aggregation on same coords of all meshes in \b ms. This returned mesh
 *         is sorted by type following the geo cell types order of MEDCouplingUMesh::sortCellsInMEDFileFrmt method.
 */
MEDCouplingUMesh *MEDCouplingUMesh::AggregateSortedByTypeMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& ms,
                                                                            DataArrayInt *&szOfCellGrpOfSameType,
                                                                            DataArrayInt *&idInMsOfCellGrpOfSameType) throw(INTERP_KERNEL::Exception)
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
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> > m1ssmAuto;
  //
  std::vector<const MEDCouplingUMesh *> m1ssmSingle;
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> > m1ssmSingleAuto;
  int fake=0,rk=0;
  std::vector<int> ret1Data;
  std::vector<int> ret2Data;
  for(std::vector<const MEDCouplingUMesh *>::const_iterator it=ms2.begin();it!=ms2.end();it++,rk++)
    {
      if(meshDim!=(*it)->getMeshDimension())
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::AggregateSortedByTypeMeshesOnSameCoords : meshdims mismatch !");
      if(refCoo!=(*it)->getCoords())
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::AggregateSortedByTypeMeshesOnSameCoords : meshes are not shared by a single coordinates coords !");
      std::vector<MEDCouplingUMesh *> sp=(*it)->splitByType();
      std::copy(sp.begin(),sp.end(),std::back_insert_iterator< std::vector<const MEDCouplingUMesh *> >(m1ssm));
      std::copy(sp.begin(),sp.end(),std::back_insert_iterator< std::vector<MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> > >(m1ssmAuto));
      for(std::vector<MEDCouplingUMesh *>::const_iterator it2=sp.begin();it2!=sp.end();it2++)
        {
          MEDCouplingUMesh *singleCell=static_cast<MEDCouplingUMesh *>((*it2)->buildPartOfMySelf(&fake,&fake+1,true));
          m1ssmSingleAuto.push_back(singleCell);
          m1ssmSingle.push_back(singleCell);
          ret1Data.push_back((*it2)->getNumberOfCells()); ret2Data.push_back(rk);
        }
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret1=DataArrayInt::New(); ret1->alloc((int)m1ssmSingle.size(),1); std::copy(ret1Data.begin(),ret1Data.end(),ret1->getPointer());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret2=DataArrayInt::New(); ret2->alloc((int)m1ssmSingle.size(),1); std::copy(ret2Data.begin(),ret2Data.end(),ret2->getPointer());
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m1ssmSingle2=MEDCouplingUMesh::MergeUMeshesOnSameCoords(m1ssmSingle);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> renum=m1ssmSingle2->sortCellsInMEDFileFrmt();
  std::vector<const MEDCouplingUMesh *> m1ssmfinal(m1ssm.size());
  for(std::size_t i=0;i<m1ssm.size();i++)
    m1ssmfinal[renum->getIJ(i,0)]=m1ssm[i];
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> ret0=MEDCouplingUMesh::MergeUMeshesOnSameCoords(m1ssmfinal);
  szOfCellGrpOfSameType=ret1->renumber(renum->getConstPointer());
  idInMsOfCellGrpOfSameType=ret2->renumber(renum->getConstPointer());
  ret0->incrRef();
  return ret0;
}

/*!
 * This method returns a newly created DataArrayInt instance.
 * This method retrieves cell ids in [begin,end) that have the type 'type'.
 */
DataArrayInt *MEDCouplingUMesh::keepCellIdsByType(INTERP_KERNEL::NormalizedCellType type, const int *begin, const int *end) const throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  std::vector<int> r;
  const int *conn=_nodal_connec->getConstPointer();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  for(const int *w=begin;w!=end;w++)
    if((INTERP_KERNEL::NormalizedCellType)conn[connIndex[*w]]==type)
      r.push_back(*w);
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc((int)r.size(),1);
  std::copy(r.begin(),r.end(),ret->getPointer());
  return ret;
}

/*!
 * This method makes the assumption that da->getNumberOfTuples()<this->getNumberOfCells(). This method makes the assumption that ids contained in 'da'
 * are in [0:getNumberOfCells())
 */
DataArrayInt *MEDCouplingUMesh::convertCellArrayPerGeoType(const DataArrayInt *da) const throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  int nbOfCells=getNumberOfCells();
  std::set<INTERP_KERNEL::NormalizedCellType> types=getAllTypes();
  int *tmp=new int[nbOfCells];
  for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iter=types.begin();iter!=types.end();iter++)
    {
      int j=0;
      for(const int *i=connI;i!=connI+nbOfCells;i++)
        if((INTERP_KERNEL::NormalizedCellType)conn[*i]==(*iter))
          tmp[std::distance(connI,i)]=j++;
    }
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(da->getNumberOfTuples(),da->getNumberOfComponents());
  ret->copyStringInfoFrom(*da);
  int *retPtr=ret->getPointer();
  const int *daPtr=da->getConstPointer();
  int nbOfElems=da->getNbOfElems();
  for(int k=0;k<nbOfElems;k++)
    retPtr[k]=tmp[daPtr[k]];
  delete [] tmp;
  return ret;
}

/*!
 * This method reduced number of cells of this by keeping cells whose type is different from 'type' and if type=='type'
 * cells whose ids is in 'idsPerGeoType' array.
 * This method conserves coords and name of mesh.
 */
MEDCouplingUMesh *MEDCouplingUMesh::keepSpecifiedCells(INTERP_KERNEL::NormalizedCellType type, const int *idsPerGeoTypeBg, const int *idsPerGeoTypeEnd) const
{
  std::vector<int> idsTokeep;
  int nbOfCells=getNumberOfCells();
  int j=0;
  for(int i=0;i<nbOfCells;i++)
    if(getTypeOfCell(i)!=type)
      idsTokeep.push_back(i);
    else
      {
        if(std::find(idsPerGeoTypeBg,idsPerGeoTypeEnd,j)!=idsPerGeoTypeEnd)
          idsTokeep.push_back(i);
        j++;
      }
  MEDCouplingPointSet *ret=buildPartOfMySelf(&idsTokeep[0],&idsTokeep[0]+idsTokeep.size(),true);
  MEDCouplingUMesh *ret2=dynamic_cast<MEDCouplingUMesh *>(ret);
  if(!ret2)
    {
      ret->decrRef();
      return 0;
    }
  ret2->copyTinyInfoFrom(this);
  return ret2;
}

/*!
 * This method returns a vector of size 'this->getNumberOfCells()'.
 * This method retrieves for each cell in 'this' if it is linear (false) or quadratic(true).
 */
std::vector<bool> MEDCouplingUMesh::getQuadraticStatus() const throw(INTERP_KERNEL::Exception)
{
  int ncell=getNumberOfCells();
  std::vector<bool> ret(ncell);
  const int *cI=getNodalConnectivityIndex()->getConstPointer();
  const int *c=getNodalConnectivity()->getConstPointer();
  for(int i=0;i<ncell;i++)
    {
      INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType)c[cI[i]];
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(typ);
      ret[i]=cm.isQuadratic();
    }
  return ret;
}

/*!
 * Returns a newly created mesh (with ref count ==1) that contains merge of 'this' and 'other'.
 */
MEDCouplingMesh *MEDCouplingUMesh::mergeMyselfWith(const MEDCouplingMesh *other) const
{
  if(other->getType()!=UNSTRUCTURED)
    throw INTERP_KERNEL::Exception("Merge of umesh only available with umesh each other !");
  const MEDCouplingUMesh *otherC=static_cast<const MEDCouplingUMesh *>(other);
  return MergeUMeshes(this,otherC);
}

/*!
 * Returns an array with this->getNumberOfCells() tuples and this->getSpaceDimension() dimension.
 * The false barycenter is computed that is to say barycenter of a cell is computed using average on each
 * components of coordinates of the cell.
 */
DataArrayDouble *MEDCouplingUMesh::getBarycenterAndOwner() const
{
  DataArrayDouble *ret=DataArrayDouble::New();
  int spaceDim=getSpaceDimension();
  int nbOfCells=getNumberOfCells();
  ret->alloc(nbOfCells,spaceDim);
  ret->copyStringInfoFrom(*getCoords());
  double *ptToFill=ret->getPointer();
  double *tmp=new double[spaceDim];
  const int *nodal=_nodal_connec->getConstPointer();
  const int *nodalI=_nodal_connec_index->getConstPointer();
  const double *coor=_coords->getConstPointer();
  for(int i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)nodal[nodalI[i]];
      INTERP_KERNEL::computeBarycenter2<int,INTERP_KERNEL::ALL_C_MODE>(type,nodal+nodalI[i]+1,nodalI[i+1]-nodalI[i]-1,coor,spaceDim,ptToFill);
      ptToFill+=spaceDim;
    }
  delete [] tmp;
  return ret;
}

/*!
 * This method is similar to MEDCouplingUMesh::getBarycenterAndOwner except that it works on subPart of 'this' without
 * building explicitely it. The input part is defined by an array [begin,end). All ids contained in this array should be less than this->getNumberOfCells().
 * No check of that will be done !
 */
DataArrayDouble *MEDCouplingUMesh::getPartBarycenterAndOwner(const int *begin, const int *end) const
{
  DataArrayDouble *ret=DataArrayDouble::New();
  int spaceDim=getSpaceDimension();
  int nbOfTuple=(int)std::distance(begin,end);
  ret->alloc(nbOfTuple,spaceDim);
  double *ptToFill=ret->getPointer();
  double *tmp=new double[spaceDim];
  const int *nodal=_nodal_connec->getConstPointer();
  const int *nodalI=_nodal_connec_index->getConstPointer();
  const double *coor=_coords->getConstPointer();
  for(const int *w=begin;w!=end;w++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)nodal[nodalI[*w]];
      INTERP_KERNEL::computeBarycenter2<int,INTERP_KERNEL::ALL_C_MODE>(type,nodal+nodalI[*w]+1,nodalI[*w+1]-nodalI[*w]-1,coor,spaceDim,ptToFill);
      ptToFill+=spaceDim;
    }
  delete [] tmp;
  return ret;
}

/*!
 * This method expects as input a DataArrayDouble non nul instance 'da' that should be allocated. If not an exception is thrown.
 * 
 */
MEDCouplingUMesh *MEDCouplingUMesh::Build0DMeshFromCoords(DataArrayDouble *da) throw(INTERP_KERNEL::Exception)
{
  if(!da)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::Build0DMeshFromCoords : instance of DataArrayDouble must be not null !");
  da->checkAllocated();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> ret=MEDCouplingUMesh::New(da->getName().c_str(),0);
  ret->setCoords(da);
  int nbOfTuples=da->getNumberOfTuples();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> c=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cI=DataArrayInt::New();
  c->alloc(2*nbOfTuples,1);
  cI->alloc(nbOfTuples+1,1);
  int *cp=c->getPointer();
  int *cip=cI->getPointer();
  *cip++=0;
  for(int i=0;i<nbOfTuples;i++)
    {
      *cp++=INTERP_KERNEL::NORM_POINT1;
      *cp++=i;
      *cip++=2*(i+1);
    }
  ret->setConnectivity(c,cI,true);
  ret->incrRef();
  return ret;
}

/*!
 * Returns a newly created mesh (with ref count ==1) that contains merge of 'mesh1' and 'other'.
 * The coords of 'mesh2' are added at the end of coords of 'mesh1'.
 */
MEDCouplingUMesh *MEDCouplingUMesh::MergeUMeshes(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2) throw(INTERP_KERNEL::Exception)
{
  std::vector<const MEDCouplingUMesh *> tmp(2);
  tmp[0]=const_cast<MEDCouplingUMesh *>(mesh1); tmp[1]=const_cast<MEDCouplingUMesh *>(mesh2);
  return MergeUMeshes(tmp);
}

/*!
 * This method returns in case of success a mesh constitued from union of all meshes in 'a'.
 * There should be \b no presence of null pointer into 'a'. If any an INTERP_KERNEL::Exception will be thrown.
 * The returned mesh will contain aggregation of nodes in 'a' (in the same order) and aggregation of
 * cells in meshes in 'a' (in the same order too).
 */
MEDCouplingUMesh *MEDCouplingUMesh::MergeUMeshes(std::vector<const MEDCouplingUMesh *>& a) throw(INTERP_KERNEL::Exception)
{
  std::size_t sz=a.size();
  if(sz==0)
    return MergeUMeshesLL(a);
  for(std::size_t ii=0;ii<sz;ii++)
    if(!a[ii])
      {
        std::ostringstream oss; oss << "MEDCouplingUMesh::MergeUMeshes : item #" << ii << " in input array of size "<< sz << " is empty !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> > bb(sz);
  std::vector< const MEDCouplingUMesh * > aa(sz);
  int spaceDim=-3;
  for(std::size_t i=0;i<sz && spaceDim==-3;i++)
    {
      const MEDCouplingUMesh *cur=a[i];
      const DataArrayDouble *coo=cur->getCoords();
      if(coo)
        spaceDim=coo->getNumberOfComponents();
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

/// @cond INTERNAL

MEDCouplingUMesh *MEDCouplingUMesh::MergeUMeshesLL(std::vector<const MEDCouplingUMesh *>& a) throw(INTERP_KERNEL::Exception)
{
  if(a.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::MergeUMeshes : input array must be NON EMPTY !");
  std::vector<const MEDCouplingUMesh *>::const_iterator it=a.begin();
  int meshDim=(*it)->getMeshDimension();
  int nbOfCells=(*it)->getNumberOfCells();
  int meshLgth=(*it++)->getMeshLength();
  for(;it!=a.end();it++)
    {
      if(meshDim!=(*it)->getMeshDimension())
        throw INTERP_KERNEL::Exception("Mesh dimensions mismatches, MergeUMeshes impossible !");
      nbOfCells+=(*it)->getNumberOfCells();
      meshLgth+=(*it)->getMeshLength();
    }
  std::vector<const MEDCouplingPointSet *> aps(a.size());
  std::copy(a.begin(),a.end(),aps.begin());
  DataArrayDouble *pts=MergeNodesArray(aps);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> ret=MEDCouplingUMesh::New("merge",meshDim);
  ret->setCoords(pts);
  pts->decrRef();
  DataArrayInt *c=DataArrayInt::New();
  c->alloc(meshLgth,1);
  int *cPtr=c->getPointer();
  DataArrayInt *cI=DataArrayInt::New();
  cI->alloc(nbOfCells+1,1);
  int *cIPtr=cI->getPointer();
  *cIPtr++=0;
  int offset=0;
  int offset2=0;
  for(it=a.begin();it!=a.end();it++)
    {
      int curNbOfCell=(*it)->getNumberOfCells();
      const int *curCI=(*it)->_nodal_connec_index->getConstPointer();
      const int *curC=(*it)->_nodal_connec->getConstPointer();
      cIPtr=std::transform(curCI+1,curCI+curNbOfCell+1,cIPtr,std::bind2nd(std::plus<int>(),offset));
      for(int j=0;j<curNbOfCell;j++)
        {
          const int *src=curC+curCI[j];
          *cPtr++=*src++;
          for(;src!=curC+curCI[j+1];src++,cPtr++)
            {
              if(*src!=-1)
                *cPtr=*src+offset2;
              else
                *cPtr=-1;
            }
        }
      offset+=curCI[curNbOfCell];
      offset2+=(*it)->getNumberOfNodes();
    }
  //
  ret->setConnectivity(c,cI,true);
  c->decrRef();
  cI->decrRef();
  ret->incrRef();
  return ret;
}

/// @endcond

/*!
 * Idem MergeUMeshes except that 'meshes' are expected to lyie on the same coords and 'meshes' have the same meshdim.
 * 'meshes' must be a non empty vector.
 */
MEDCouplingUMesh *MEDCouplingUMesh::MergeUMeshesOnSameCoords(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2) throw(INTERP_KERNEL::Exception)
{
  std::vector<const MEDCouplingUMesh *> tmp(2);
  tmp[0]=mesh1; tmp[1]=mesh2;
  return MergeUMeshesOnSameCoords(tmp);
}

/*!
 * Idem MergeUMeshes except that 'meshes' are expected to lyie on the same coords and 'meshes' have the same meshdim.
 * 'meshes' must be a non empty vector.
 */
MEDCouplingUMesh *MEDCouplingUMesh::MergeUMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& meshes)
{
  if(meshes.empty())
    throw INTERP_KERNEL::Exception("meshes input parameter is expected to be non empty.");
  for(std::size_t ii=0;ii<meshes.size();ii++)
    if(!meshes[ii])
      {
        std::ostringstream oss; oss << "MEDCouplingUMesh::MergeUMeshesOnSameCoords : item #" << ii << " in input array of size "<< meshes.size() << " is empty !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  const DataArrayDouble *coords=meshes.front()->getCoords();
  int meshDim=meshes.front()->getMeshDimension();
  std::vector<const MEDCouplingUMesh *>::const_iterator iter=meshes.begin();
  int meshLgth=0;
  int meshIndexLgth=0;
  for(;iter!=meshes.end();iter++)
    {
      if(coords!=(*iter)->getCoords())
        throw INTERP_KERNEL::Exception("meshes does not share the same coords ! Try using tryToShareSameCoords method !");
      if(meshDim!=(*iter)->getMeshDimension())
        throw INTERP_KERNEL::Exception("Mesh dimensions mismatches, FuseUMeshesOnSameCoords impossible !");
      meshLgth+=(*iter)->getMeshLength();
      meshIndexLgth+=(*iter)->getNumberOfCells();
    }
  DataArrayInt *nodal=DataArrayInt::New();
  nodal->alloc(meshLgth,1);
  int *nodalPtr=nodal->getPointer();
  DataArrayInt *nodalIndex=DataArrayInt::New();
  nodalIndex->alloc(meshIndexLgth+1,1);
  int *nodalIndexPtr=nodalIndex->getPointer();
  int offset=0;
  for(iter=meshes.begin();iter!=meshes.end();iter++)
    {
      const int *nod=(*iter)->getNodalConnectivity()->getConstPointer();
      const int *index=(*iter)->getNodalConnectivityIndex()->getConstPointer();
      int nbOfCells=(*iter)->getNumberOfCells();
      int meshLgth2=(*iter)->getMeshLength();
      nodalPtr=std::copy(nod,nod+meshLgth2,nodalPtr);
      if(iter!=meshes.begin())
        nodalIndexPtr=std::transform(index+1,index+nbOfCells+1,nodalIndexPtr,std::bind2nd(std::plus<int>(),offset));
      else
        nodalIndexPtr=std::copy(index,index+nbOfCells+1,nodalIndexPtr);
      offset+=meshLgth2;
    }
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New();
  ret->setName("merge");
  ret->setMeshDimension(meshDim);
  ret->setConnectivity(nodal,nodalIndex,true);
  ret->setCoords(coords);
  nodalIndex->decrRef();
  nodal->decrRef();
  return ret;
}

/*!
 * This method fuses meshes 'meshes' and returns the fused mesh and the correspondances arrays for each mesh in 'meshes' in returned mesh.
 * If a same cell is detected in several meshes in 'meshes', this cell will appear only once in returned mesh (see ParaMEDMEM::MEDCouplingUMesh::zipConnectivityTraducer for more details)
 *
 * @param meshes input non empty vector containing meshes having same coordiantes array and same mesh dimension.
 * @param compType see MEDCouplingUMesh::zipConnectivityTraducer
 * @param corr output vector with same size as 'meshes' parameter. corr[i] is the correspondance array of mesh meshes[i] in returned mesh.
 *             The arrays contained in 'corr' parameter are returned with refcounter set to one.
 *             To avoid memory leaks the caller have to deal with each instances of DataArrayInt contained in 'corr' parameter.
 * @return The mesh lying on the same coordinates than those in meshes. All cells in 'meshes' are in returned mesh with 
 * @exception if meshes is a empty vector or meshes are not lying on same coordinates or meshes not have the same dimension.
 */
MEDCouplingUMesh *MEDCouplingUMesh::FuseUMeshesOnSameCoords(const std::vector<const MEDCouplingUMesh *>& meshes, int compType, std::vector<DataArrayInt *>& corr)
{
  //All checks are delegated to MergeUMeshesOnSameCoords
  MEDCouplingUMesh *ret=MergeUMeshesOnSameCoords(meshes);
  DataArrayInt *o2n=ret->zipConnectivityTraducer(compType);
  corr.resize(meshes.size());
  std::size_t nbOfMeshes=meshes.size();
  int offset=0;
  const int *o2nPtr=o2n->getConstPointer();
  for(std::size_t i=0;i<nbOfMeshes;i++)
    {
      DataArrayInt *tmp=DataArrayInt::New();
      int curNbOfCells=meshes[i]->getNumberOfCells();
      tmp->alloc(curNbOfCells,1);
      std::copy(o2nPtr+offset,o2nPtr+offset+curNbOfCells,tmp->getPointer());
      offset+=curNbOfCells;
      tmp->setName(meshes[i]->getName());
      corr[i]=tmp;
    }
  o2n->decrRef();
  return ret;
}

/*!
 * This method takes in input meshes \b meshes containing no null reference. If any an INTERP_KERNEL::Exception will be thrown.
 * \b meshes should have a good coherency (connectivity and coordinates well defined).
 * All mesh in \b meshes must have the same space dimension. If not an INTERP_KERNEL:Exception will be thrown.
 * But mesh in \b meshes \b can \b have \b different \b mesh \b dimension \b each \b other.
 *
 * This method performs nothing if size of \b meshes is in [0,1].
 * This method is particulary usefull in MEDLoader context to build a \ref ParaMEDMEM::MEDFileUMesh "MEDFileUMesh" instance that expects that underlying
 * coordinates DataArrayDouble instance.
 *
 * \param [in,out] meshes : vector containing no null instance of MEDCouplingUMesh that in case of success of this method will be modified.
 */
void MEDCouplingUMesh::PutUMeshesOnSameAggregatedCoords(const std::vector<MEDCouplingUMesh *>& meshes) throw(INTERP_KERNEL::Exception)
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
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << " MEDCouplingUMesh::PutUMeshesOnSameAggregatedCoords : Item #" << std::distance(meshes.begin(),it) << " inside the vector of length " << meshes.size();
          oss << " is null !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> res=DataArrayDouble::Aggregate(coords);
  std::vector<MEDCouplingUMesh *>::const_iterator it=meshes.begin();
  int offset=(*it)->getNumberOfNodes();
  (*it++)->setCoords(res);
  for(;it!=meshes.end();it++)
    {
      int oldNumberOfNodes=(*it)->getNumberOfNodes();
      (*it)->setCoords(res);
      (*it)->shiftNodeNumbersInConn(offset);
      offset+=oldNumberOfNodes;
    }
}

/*!
 * This method takes in input meshes \b meshes containing no null reference. If any an INTERP_KERNEL::Exception will be thrown.
 * \b meshes should have a good coherency (connectivity and coordinates well defined).
 * All mesh in \b meshes must have the same space dimension. If not an INTERP_KERNEL:Exception will be thrown.
 * But mesh in \b meshes \b can \b have \b different \b mesh \b dimension \b each \b other.
 * If \b meshes share the same instance of DataArrayDouble as coordinates and that this instance is null, this method do nothing and no exception will be thrown.
 *
 * This method performs nothing if size of \b meshes is empty.
 * This method is particulary usefull in MEDLoader context to perform a treatment of a MEDFileUMesh instance on different levels.
 * coordinates DataArrayDouble instance.
 *
 * \param [in,out] meshes :vector containing no null instance of MEDCouplingUMesh sharing the same DataArrayDouble instance of coordinates, that in case of success of this method will be modified.
 * \param [in] eps is the distance in absolute (that should be positive !), so that 2 or more points within a distance of eps will be merged into a single point.
 */
void MEDCouplingUMesh::MergeNodesOnUMeshesSharingSameCoords(const std::vector<MEDCouplingUMesh *>& meshes, double eps) throw(INTERP_KERNEL::Exception)
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
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  if(s.size()!=1)
    {
      std::ostringstream oss; oss << "MEDCouplingUMesh::MergeNodesOnUMeshesSharingSameCoords : In input vector of unstructured meshes of size " << meshes.size() << ", it appears that they do not share the same instance of DataArrayDouble for coordiantes ! tryToShareSameCoordsPermute method can help to reach that !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const DataArrayDouble *coo=*(s.begin());
  if(!coo)
    return;
  //
  DataArrayInt *comm,*commI;
  coo->findCommonTuples(eps,-1,comm,commI);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp1(comm),tmp2(commI);
  int oldNbOfNodes=coo->getNumberOfTuples();
  int newNbOfNodes;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> o2n=DataArrayInt::BuildOld2NewArrayFromSurjectiveFormat2(oldNbOfNodes,comm,commI,newNbOfNodes);
  if(oldNbOfNodes==newNbOfNodes)
    return ;
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> newCoords=coo->renumberAndReduce(o2n->getConstPointer(),newNbOfNodes);
  for(std::vector<MEDCouplingUMesh *>::const_iterator it=meshes.begin();it!=meshes.end();it++)
    {
      (*it)->renumberNodesInConn(o2n->getConstPointer());
      (*it)->setCoords(newCoords);
    } 
}

/*!
 * This method takes in input a cell defined by its MEDcouplingUMesh connectivity [connBg,connEnd) and returns its extruded cell by inserting the result at the end of ret.
 * @param nbOfNodesPerLev in parameter that specifies the number of nodes of one slice of global dataset
 * @param isQuad specifies the policy of connectivity.
 * @ret in/out parameter in which the result will be append
 */
void MEDCouplingUMesh::AppendExtrudedCell(const int *connBg, const int *connEnd, int nbOfNodesPerLev, bool isQuad, std::vector<int>& ret)
{
  INTERP_KERNEL::NormalizedCellType flatType=(INTERP_KERNEL::NormalizedCellType)connBg[0];
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(flatType);
  ret.push_back(cm.getExtrudedType());
  int deltaz=isQuad?2*nbOfNodesPerLev:nbOfNodesPerLev;
  switch(flatType)
    {
    case INTERP_KERNEL::NORM_POINT1:
      {
        ret.push_back(connBg[1]);
        ret.push_back(connBg[1]+nbOfNodesPerLev);
        break;
      }
    case INTERP_KERNEL::NORM_SEG2:
      {
        int conn[4]={connBg[1],connBg[2],connBg[2]+deltaz,connBg[1]+deltaz};
        ret.insert(ret.end(),conn,conn+4);
        break;
      }
    case INTERP_KERNEL::NORM_SEG3:
      {
        int conn[8]={connBg[1],connBg[3],connBg[3]+deltaz,connBg[1]+deltaz,connBg[2],connBg[3]+nbOfNodesPerLev,connBg[2]+deltaz,connBg[1]+nbOfNodesPerLev};
        ret.insert(ret.end(),conn,conn+8);
        break;
      }
    case INTERP_KERNEL::NORM_QUAD4:
      {
        int conn[8]={connBg[1],connBg[2],connBg[3],connBg[4],connBg[1]+deltaz,connBg[2]+deltaz,connBg[3]+deltaz,connBg[4]+deltaz};
        ret.insert(ret.end(),conn,conn+8);
        break;
      }
    case INTERP_KERNEL::NORM_TRI3:
      {
        int conn[6]={connBg[1],connBg[2],connBg[3],connBg[1]+deltaz,connBg[2]+deltaz,connBg[3]+deltaz};
        ret.insert(ret.end(),conn,conn+6);
        break;
      }
    case INTERP_KERNEL::NORM_TRI6:
      {
        int conn[15]={connBg[1],connBg[2],connBg[3],connBg[1]+deltaz,connBg[2]+deltaz,connBg[3]+deltaz,connBg[4],connBg[5],connBg[6],connBg[4]+deltaz,connBg[5]+deltaz,connBg[6]+deltaz,
                      connBg[1]+nbOfNodesPerLev,connBg[2]+nbOfNodesPerLev,connBg[3]+nbOfNodesPerLev};
        ret.insert(ret.end(),conn,conn+15);
        break;
      }
    case INTERP_KERNEL::NORM_QUAD8:
      {
        int conn[20]={
          connBg[1],connBg[2],connBg[3],connBg[4],connBg[1]+deltaz,connBg[2]+deltaz,connBg[3]+deltaz,connBg[4]+deltaz,
          connBg[5],connBg[6],connBg[7],connBg[8],connBg[5]+deltaz,connBg[6]+deltaz,connBg[7]+deltaz,connBg[8]+deltaz,
          connBg[1]+nbOfNodesPerLev,connBg[2]+nbOfNodesPerLev,connBg[3]+nbOfNodesPerLev,connBg[4]+nbOfNodesPerLev
        };
        ret.insert(ret.end(),conn,conn+20);
        break;
      }
    case INTERP_KERNEL::NORM_POLYGON:
      {
        std::back_insert_iterator< std::vector<int> > ii(ret);
        std::copy(connBg+1,connEnd,ii);
        *ii++=-1;
        std::reverse_iterator<const int *> rConnBg(connEnd);
        std::reverse_iterator<const int *> rConnEnd(connBg+1);
        std::transform(rConnBg,rConnEnd,ii,std::bind2nd(std::plus<int>(),deltaz));
        std::size_t nbOfRadFaces=std::distance(connBg+1,connEnd);
        for(std::size_t i=0;i<nbOfRadFaces;i++)
          {
            *ii++=-1;
            int conn[4]={connBg[(i+1)%nbOfRadFaces+1],connBg[i+1],connBg[i+1]+deltaz,connBg[(i+1)%nbOfRadFaces+1]+deltaz};
            std::copy(conn,conn+4,ii);
          }
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("A flat type has been detected that has not its extruded representation !");
    }
}

/*!
 * This static operates only for coords in 3D. The polygon is specfied by its connectivity nodes in [begin,end).
 */
bool MEDCouplingUMesh::IsPolygonWellOriented(bool isQuadratic, const double *vec, const int *begin, const int *end, const double *coords)
{
  double v[3]={0.,0.,0.};
  std::size_t sz=std::distance(begin,end);
  if(isQuadratic)
    sz/=2;
  for(std::size_t i=0;i<sz;i++)
    {
      v[0]+=coords[3*begin[i]+1]*coords[3*begin[(i+1)%sz]+2]-coords[3*begin[i]+2]*coords[3*begin[(i+1)%sz]+1];
      v[1]+=coords[3*begin[i]+2]*coords[3*begin[(i+1)%sz]]-coords[3*begin[i]]*coords[3*begin[(i+1)%sz]+2];
      v[2]+=coords[3*begin[i]]*coords[3*begin[(i+1)%sz]+1]-coords[3*begin[i]+1]*coords[3*begin[(i+1)%sz]];
    }
  return vec[0]*v[0]+vec[1]*v[1]+vec[2]*v[2]>0.;
}

/*!
 * The polyhedron is specfied by its connectivity nodes in [begin,end).
 */
bool MEDCouplingUMesh::IsPolyhedronWellOriented(const int *begin, const int *end, const double *coords)
{
  std::vector<std::pair<int,int> > edges;
  std::size_t nbOfFaces=std::count(begin,end,-1)+1;
  const int *bgFace=begin;
  for(std::size_t i=0;i<nbOfFaces;i++)
    {
      const int *endFace=std::find(bgFace+1,end,-1);
      std::size_t nbOfEdgesInFace=std::distance(bgFace,endFace);
      for(std::size_t j=0;j<nbOfEdgesInFace;j++)
        {
          std::pair<int,int> p1(bgFace[j],bgFace[(j+1)%nbOfEdgesInFace]);
          if(std::find(edges.begin(),edges.end(),p1)!=edges.end())
            return false;
          edges.push_back(p1);
        }
      bgFace=endFace+1;
    }
  return INTERP_KERNEL::calculateVolumeForPolyh2<int,INTERP_KERNEL::ALL_C_MODE>(begin,(int)std::distance(begin,end),coords)>-EPS_FOR_POLYH_ORIENTATION;
}

/*!
 * This method performs a simplyfication of a single polyedron cell. To do that each face of cell whose connectivity is defined by [\b begin, \b end) 
 * is compared with the others in order to find faces in the same plane (with approx of eps). If any, the cells are grouped together and projected to
 * a 2D space.
 *
 * \param [in] eps is a relative precision that allows to establish if some 3D plane are coplanar or not.
 * \param [in] coords the coordinates with nb of components exactly equal to 3
 * \param [in] begin begin of the nodal connectivity (geometric type included) of a single polyhedron cell
 * \param [in] end end of nodal connectivity of a single polyhedron cell (excluded)
 * \param [out] the result is put at the end of the vector without any alteration of the data.
 */
void MEDCouplingUMesh::SimplifyPolyhedronCell(double eps, const DataArrayDouble *coords, const int *begin, const int *end, std::vector<int>& res) throw(INTERP_KERNEL::Exception)
{
  int nbFaces=std::count(begin+1,end,-1)+1;
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> v=DataArrayDouble::New(); v->alloc(nbFaces,3);
  double *vPtr=v->getPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> p=DataArrayDouble::New(); p->alloc(nbFaces,1);
  double *pPtr=p->getPointer();
  const int *stFaceConn=begin+1;
  for(int i=0;i<nbFaces;i++,vPtr+=3,pPtr++)
    {
      const int *endFaceConn=std::find(stFaceConn,end,-1);
      ComputeVecAndPtOfFace(eps,coords->getConstPointer(),stFaceConn,endFaceConn,vPtr,pPtr);
      stFaceConn=endFaceConn+1;
    }
  pPtr=p->getPointer(); vPtr=v->getPointer();
  DataArrayInt *comm1=0,*commI1=0;
  v->findCommonTuples(eps,-1,comm1,commI1);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> comm1Auto(comm1),commI1Auto(commI1);
  const int *comm1Ptr=comm1->getConstPointer();
  const int *commI1Ptr=commI1->getConstPointer();
  int nbOfGrps1=commI1Auto->getNumberOfTuples()-1;
  res.push_back((int)INTERP_KERNEL::NORM_POLYHED);
  //
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mm=MEDCouplingUMesh::New("",3);
  mm->setCoords(const_cast<DataArrayDouble *>(coords)); mm->allocateCells(1); mm->insertNextCell(INTERP_KERNEL::NORM_POLYHED,(int)std::distance(begin+1,end),begin+1);
  mm->finishInsertingCells();
  //
  for(int i=0;i<nbOfGrps1;i++)
    {
      int vecId=comm1Ptr[commI1Ptr[i]];
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> tmpgrp2=p->selectByTupleId(comm1Ptr+commI1Ptr[i],comm1Ptr+commI1Ptr[i+1]);
      DataArrayInt *comm2=0,*commI2=0;
      tmpgrp2->findCommonTuples(eps,-1,comm2,commI2);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> comm2Auto(comm2),commI2Auto(commI2);
      const int *comm2Ptr=comm2->getConstPointer();
      const int *commI2Ptr=commI2->getConstPointer();
      int nbOfGrps2=commI2Auto->getNumberOfTuples()-1;
      for(int j=0;j<nbOfGrps2;j++)
        {
          if(commI2Ptr[j+1]-commI2Ptr[j]<=1)
            {
              res.insert(res.end(),begin,end);
              res.push_back(-1);
            }
          else
            {
              int pointId=comm1Ptr[commI1Ptr[i]+comm2Ptr[commI2Ptr[j]]];
              MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ids2=comm2->selectByTupleId2(commI2Ptr[j],commI2Ptr[j+1],1);
              ids2->transformWithIndArr(comm1Ptr+commI1Ptr[i],comm1Ptr+commI1Ptr[i+1]);
              DataArrayInt *tmp0=DataArrayInt::New(),*tmp1=DataArrayInt::New(),*tmp2=DataArrayInt::New(),*tmp3=DataArrayInt::New();
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mm2=mm->buildDescendingConnectivity(tmp0,tmp1,tmp2,tmp3); tmp0->decrRef(); tmp1->decrRef(); tmp2->decrRef(); tmp3->decrRef();
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mm3=static_cast<MEDCouplingUMesh *>(mm2->buildPartOfMySelf(ids2->begin(),ids2->end(),true));
              MEDCouplingAutoRefCountObjectPtr<DataArrayInt> idsNodeTmp=mm3->zipCoordsTraducer();
              MEDCouplingAutoRefCountObjectPtr<DataArrayInt> idsNode=idsNodeTmp->invertArrayO2N2N2O(mm3->getNumberOfNodes());
              const int *idsNodePtr=idsNode->getConstPointer();
              double center[3]; center[0]=pPtr[pointId]*vPtr[3*vecId]; center[1]=pPtr[pointId]*vPtr[3*vecId+1]; center[2]=pPtr[pointId]*vPtr[3*vecId+2];
              double vec[3]; vec[0]=vPtr[3*vecId+1]; vec[1]=-vPtr[3*vecId]; vec[2]=0.;
              double norm=vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
              if(std::abs(norm)>eps)
                {
                  double angle=INTERP_KERNEL::EdgeArcCircle::SafeAsin(norm);
                  mm3->rotate(center,vec,angle);
                }
              mm3->changeSpaceDimension(2);
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mm4=mm3->buildSpreadZonesWithPoly();
              const int *conn4=mm4->getNodalConnectivity()->getConstPointer();
              const int *connI4=mm4->getNodalConnectivityIndex()->getConstPointer();
              int nbOfCells=mm4->getNumberOfCells();
              for(int k=0;k<nbOfCells;k++)
                {
                  int l=0;
                  for(const int *work=conn4+connI4[k]+1;work!=conn4+connI4[k+1];work++,l++)
                    res.push_back(idsNodePtr[*work]);
                  res.push_back(-1);
                }
            }
        }
    }
  res.pop_back();
}

/*!
 * This method computes the normalized vector of the plane and the pos of the point belonging to the plane and the line defined by the vector going
 * through origin. The plane is defined by its nodal connectivity [\b begin, \b end).
 * 
 * \param [in] eps below that value the dot product of 2 vectors is considered as colinears
 * \param [in] coords coordinates expected to have 3 components.
 * \param [in] begin start of the nodal connectivity of the face.
 * \param [in] end end of the nodal connectivity (excluded) of the face.
 * \param [out] v the normalized vector of size 3
 * \param [out] p the pos of plane
 */
void MEDCouplingUMesh::ComputeVecAndPtOfFace(double eps, const double *coords, const int *begin, const int *end, double *v, double *p) throw(INTERP_KERNEL::Exception)
{
  std::size_t nbPoints=std::distance(begin,end);
  if(nbPoints<3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ComputeVecAndPtOfFace : < of 3 points in face ! not able to find a plane on that face !");
  double vec[3];
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
void MEDCouplingUMesh::TryToCorrectPolyhedronOrientation(int *begin, int *end, const double *coords) throw(INTERP_KERNEL::Exception)
{
  std::vector<std::pair<int,int> > edges;
  std::size_t nbOfFaces=std::count(begin,end,-1)+1;
  int *bgFace=begin;
  std::vector<bool> isPerm(nbOfFaces);
  for(std::size_t i=0;i<nbOfFaces;i++)
    {
      int *endFace=std::find(bgFace+1,end,-1);
      std::size_t nbOfEdgesInFace=std::distance(bgFace,endFace);
      for(std::size_t l=0;l<nbOfEdgesInFace;l++)
        {
          std::pair<int,int> p1(bgFace[l],bgFace[(l+1)%nbOfEdgesInFace]);
          edges.push_back(p1);
        }
      int *bgFace2=endFace+1;
      for(std::size_t k=i+1;k<nbOfFaces;k++)
        {
          int *endFace2=std::find(bgFace2+1,end,-1);
          std::size_t nbOfEdgesInFace2=std::distance(bgFace2,endFace2);
          for(std::size_t j=0;j<nbOfEdgesInFace2;j++)
            {
              std::pair<int,int> p2(bgFace2[j],bgFace2[(j+1)%nbOfEdgesInFace2]);
              if(std::find(edges.begin(),edges.end(),p2)!=edges.end())
                {
                  if(isPerm[k])
                    throw INTERP_KERNEL::Exception("Fail to repare polyhedron ! Polyedron looks bad !");
                  std::vector<int> tmp(nbOfEdgesInFace2-1);
                  std::copy(bgFace2+1,endFace2,tmp.rbegin());
                  std::copy(tmp.begin(),tmp.end(),bgFace2+1);
                  isPerm[k]=true;
                  continue;
                }
            }
          bgFace2=endFace2+1;
        }
      bgFace=endFace+1;
    }
  if(INTERP_KERNEL::calculateVolumeForPolyh2<int,INTERP_KERNEL::ALL_C_MODE>(begin,(int)std::distance(begin,end),coords)<-EPS_FOR_POLYH_ORIENTATION)
    {//not lucky ! The first face was not correctly oriented : reorient all faces...
      bgFace=begin;
      for(std::size_t i=0;i<nbOfFaces;i++)
        {
          int *endFace=std::find(bgFace+1,end,-1);
          std::size_t nbOfEdgesInFace=std::distance(bgFace,endFace);
          std::vector<int> tmp(nbOfEdgesInFace-1);
          std::copy(bgFace+1,endFace,tmp.rbegin());
          std::copy(tmp.begin(),tmp.end(),bgFace+1);
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
DataArrayInt *MEDCouplingUMesh::buildUnionOf2DMesh() const throw(INTERP_KERNEL::Exception)
{
  if(getMeshDimension()!=2 || getSpaceDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildUnionOf2DMesh : meshdimension, spacedimension must be equal to 2 !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m=computeSkin();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> o2n=m->zipCoordsTraducer();
  int nbOfNodesExpected=m->getNumberOfNodes();
  if(m->getNumberOfCells()!=nbOfNodesExpected)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildUnionOf2DMesh : the mesh 2D in input appears to be not in a single part or a quadratic 2D mesh !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> n2o=o2n->invertArrayO2N2N2O(m->getNumberOfNodes());
  const int *n2oPtr=n2o->getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> revNodal(DataArrayInt::New()),revNodalI(DataArrayInt::New());
  m->getReverseNodalConnectivity(revNodal,revNodalI);
  const int *revNodalPtr=revNodal->getConstPointer(),*revNodalIPtr=revNodalI->getConstPointer();
  const int *nodalPtr=m->getNodalConnectivity()->getConstPointer();
  const int *nodalIPtr=m->getNodalConnectivityIndex()->getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(nbOfNodesExpected+1,1);
  int *work=ret->getPointer();  *work++=INTERP_KERNEL::NORM_POLYGON;
  if(nbOfNodesExpected<1)
    { ret->incrRef(); return ret; }
  int prevCell=0;
  int prevNode=nodalPtr[nodalIPtr[0]+1];
  *work++=n2oPtr[prevNode];
  for(int i=1;i<nbOfNodesExpected;i++)
    {
      if(nodalIPtr[prevCell+1]-nodalIPtr[prevCell]==3)
        {
          std::set<int> conn(nodalPtr+nodalIPtr[prevCell]+1,nodalPtr+nodalIPtr[prevCell]+3);
          conn.erase(prevNode);
          if(conn.size()==1)
            {
              int curNode=*(conn.begin());
              *work++=n2oPtr[curNode];
              std::set<int> shar(revNodalPtr+revNodalIPtr[curNode],revNodalPtr+revNodalIPtr[curNode+1]);
              shar.erase(prevCell);
              if(shar.size()==1)
                {
                  prevCell=*(shar.begin());
                  prevNode=curNode;
                }
              else
                throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildUnionOf2DMesh : presence of unexpected 2 !");
            }
          else
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildUnionOf2DMesh : presence of unexpected 1 !");
        }
      else
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildUnionOf2DMesh : presence of unexpected cell !");
    }
  ret->incrRef(); return ret;
}

/*!
 * This method makes the assumption spacedimension == meshdimension == 3.
 * This method works only for linear cells.
 * 
 * \return a newly allocated array containing the connectivity of a polygon type enum included (NORM_POLYHED in pos#0)
 */
DataArrayInt *MEDCouplingUMesh::buildUnionOf3DMesh() const throw(INTERP_KERNEL::Exception)
{
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildUnionOf3DMesh : meshdimension, spacedimension must be equal to 2 !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m=computeSkin();
  const int *conn=m->getNodalConnectivity()->getConstPointer();
  const int *connI=m->getNodalConnectivityIndex()->getConstPointer();
  int nbOfCells=m->getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(m->getNodalConnectivity()->getNumberOfTuples(),1);
  int *work=ret->getPointer();  *work++=INTERP_KERNEL::NORM_POLYHED;
  if(nbOfCells<1)
    { ret->incrRef(); return ret; }
  work=std::copy(conn+connI[0]+1,conn+connI[1],work);
  for(int i=1;i<nbOfCells;i++)
    {
      *work++=-1;
      work=std::copy(conn+connI[i]+1,conn+connI[i+1],work);
    }
  ret->incrRef();
  return ret;
}

/*!
 * This method put in zip format into parameter 'zipFrmt' in full interlace mode.
 * This format is often asked by INTERP_KERNEL algorithms to avoid many indirections into coordinates array.
 */
void MEDCouplingUMesh::FillInCompact3DMode(int spaceDim, int nbOfNodesInCell, const int *conn, const double *coo, double *zipFrmt) throw(INTERP_KERNEL::Exception)
{
  double *w=zipFrmt;
  if(spaceDim==3)
    for(int i=0;i<nbOfNodesInCell;i++)
      w=std::copy(coo+3*conn[i],coo+3*conn[i]+3,w);
  else if(spaceDim==2)
    {
      for(int i=0;i<nbOfNodesInCell;i++)
        {
          w=std::copy(coo+2*conn[i],coo+2*conn[i]+2,w);
          *w++=0.;
        }
    }
  else
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::FillInCompact3DMode : Invalid spaceDim specified : must be 2 or 3 !");
}

void MEDCouplingUMesh::writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfCells()<=0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::writeVTK : the unstructured mesh has no cells !");
  static const int PARAMEDMEM2VTKTYPETRADUCER[INTERP_KERNEL::NORM_MAXTYPE+1]={1,3,21,5,9,7,22,-1,23,-1,-1,-1,-1,-1,10,14,13,-1,12,-1,24,-1,16,27,-1,26,-1,-1,-1,-1,25,42,-1,4};
  ofs << "  <" << getVTKDataSetType() << ">\n";
  ofs << "    <Piece NumberOfPoints=\"" << getNumberOfNodes() << "\" NumberOfCells=\"" << getNumberOfCells() << "\">\n";
  ofs << "      <PointData>\n" << pointData << std::endl;
  ofs << "      </PointData>\n";
  ofs << "      <CellData>\n" << cellData << std::endl;
  ofs << "      </CellData>\n";
  ofs << "      <Points>\n";
  if(getSpaceDimension()==3)
    _coords->writeVTK(ofs,8,"Points");
  else
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> coo=_coords->changeNbOfComponents(3,0.);
      coo->writeVTK(ofs,8,"Points");
    }
  ofs << "      </Points>\n";
  ofs << "      <Cells>\n";
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> c0=_nodal_connec_index->buildComplement(_nodal_connec->getNumberOfTuples()+1);
  c0=_nodal_connec->selectByTupleId(c0->begin(),c0->end());
  c0->writeVTK(ofs,8,"Int64","connectivity");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> c1=_nodal_connec_index->deltaShiftIndex();
  c1->applyLin(1,-1);
  c1->computeOffsets2();
  c1=c1->selectByTupleId2(1,c1->getNumberOfTuples(),1);
  c1->writeVTK(ofs,8,"Int64","offsets");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> c2=_nodal_connec->selectByTupleId(_nodal_connec_index->getConstPointer(),_nodal_connec_index->getConstPointer()+getNumberOfCells());
  c2->transformWithIndArr(PARAMEDMEM2VTKTYPETRADUCER,PARAMEDMEM2VTKTYPETRADUCER+INTERP_KERNEL::NORM_MAXTYPE);
  c2->writeVTK(ofs,8,"UInt8","types");
  ofs << "      </Cells>\n";
  ofs << "    </Piece>\n";
  ofs << "  </" << getVTKDataSetType() << ">\n";
}

std::string MEDCouplingUMesh::getVTKDataSetType() const throw(INTERP_KERNEL::Exception)
{
  return std::string("UnstructuredGrid");
}

/// @cond INTERNAL

MEDCouplingUMesh *MEDCouplingUMesh::Intersect2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps, DataArrayInt *&cellNb1, DataArrayInt *&cellNb2) throw(INTERP_KERNEL::Exception)
{
  m1->checkFullyDefined();
  m2->checkFullyDefined();
  if(m1->getMeshDimension()!=2 || m1->getSpaceDimension()!=2 || m2->getMeshDimension()!=2 || m2->getSpaceDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::Intersect2DMeshes works on umeshes m1 AND m2  with meshdim equal to 2 and spaceDim equal to 2 too!");
  std::vector< std::vector<int> > intersectEdge1, colinear2, subDiv2;
  MEDCouplingUMesh *m1Desc=0,*m2Desc=0;
  DataArrayInt *desc1=0,*descIndx1=0,*revDesc1=0,*revDescIndx1=0,*desc2=0,*descIndx2=0,*revDesc2=0,*revDescIndx2=0;
  std::vector<double> addCoo,addCoordsQuadratic;
  INTERP_KERNEL::QUADRATIC_PLANAR::_precision=eps;
  INTERP_KERNEL::QUADRATIC_PLANAR::_arc_detection_precision=eps;
  IntersectDescending2DMeshes(m1,m2,eps,intersectEdge1,colinear2, subDiv2,m1Desc,desc1,descIndx1,revDesc1,revDescIndx1,
                              m2Desc,desc2,descIndx2,revDesc2,revDescIndx2,addCoo);
  revDesc1->decrRef(); revDescIndx1->decrRef(); revDesc2->decrRef(); revDescIndx2->decrRef();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> dd1(desc1),dd2(descIndx1),dd3(desc2),dd4(descIndx2);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> dd5(m1Desc),dd6(m2Desc);
  std::vector< std::vector<int> > intersectEdge2;
  BuildIntersectEdges(m1Desc,m2Desc,addCoo,subDiv2,intersectEdge2);
  subDiv2.clear(); dd5=0; dd6=0;
  std::vector<int> cr,crI;
  std::vector<int> cNb1,cNb2;
  BuildIntersecting2DCellsFromEdges(eps,m1,desc1->getConstPointer(),descIndx1->getConstPointer(),intersectEdge1,colinear2,m2,desc2->getConstPointer(),descIndx2->getConstPointer(),intersectEdge2,addCoo,
                                    /* outputs -> */addCoordsQuadratic,cr,crI,cNb1,cNb2);
  //
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> addCooDa=DataArrayDouble::New();
  addCooDa->alloc((int)(addCoo.size())/2,2);
  std::copy(addCoo.begin(),addCoo.end(),addCooDa->getPointer());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> addCoordsQuadraticDa=DataArrayDouble::New();
  addCoordsQuadraticDa->alloc((int)(addCoordsQuadratic.size())/2,2);
  std::copy(addCoordsQuadratic.begin(),addCoordsQuadratic.end(),addCoordsQuadraticDa->getPointer());
  std::vector<const DataArrayDouble *> coordss(4);
  coordss[0]=m1->getCoords(); coordss[1]=m2->getCoords(); coordss[2]=addCooDa; coordss[3]=addCoordsQuadraticDa;
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> coo=DataArrayDouble::Aggregate(coordss);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> ret=MEDCouplingUMesh::New("Intersect2D",2);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> conn=DataArrayInt::New(); conn->alloc((int)cr.size(),1); std::copy(cr.begin(),cr.end(),conn->getPointer());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> connI=DataArrayInt::New(); connI->alloc((int)crI.size(),1); std::copy(crI.begin(),crI.end(),connI->getPointer());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> c1=DataArrayInt::New(); c1->alloc((int)cNb1.size(),1); std::copy(cNb1.begin(),cNb1.end(),c1->getPointer()); cellNb1=c1;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> c2=DataArrayInt::New(); c2->alloc((int)cNb2.size(),1); std::copy(cNb2.begin(),cNb2.end(),c2->getPointer()); cellNb2=c2;
  ret->setConnectivity(conn,connI,true);
  ret->setCoords(coo);
  ret->incrRef(); c1->incrRef(); c2->incrRef();
  return ret;
}

/// @endcond

void MEDCouplingUMesh::BuildIntersecting2DCellsFromEdges(double eps, const MEDCouplingUMesh *m1, const int *desc1, const int *descIndx1,
                                                         const std::vector<std::vector<int> >& intesctEdges1, const std::vector< std::vector<int> >& colinear2,
                                                         const MEDCouplingUMesh *m2, const int *desc2, const int *descIndx2, const std::vector<std::vector<int> >& intesctEdges2,
                                                         const std::vector<double>& addCoords,
                                                         std::vector<double>& addCoordsQuadratic, std::vector<int>& cr, std::vector<int>& crI, std::vector<int>& cNb1, std::vector<int>& cNb2)
{
  static const int SPACEDIM=2;
  std::vector<double> bbox1,bbox2;
  const double *coo1=m1->getCoords()->getConstPointer();
  const int *conn1=m1->getNodalConnectivity()->getConstPointer();
  const int *connI1=m1->getNodalConnectivityIndex()->getConstPointer();
  int offset1=m1->getNumberOfNodes();
  const double *coo2=m2->getCoords()->getConstPointer();
  const int *conn2=m2->getNodalConnectivity()->getConstPointer();
  const int *connI2=m2->getNodalConnectivityIndex()->getConstPointer();
  int offset2=offset1+m2->getNumberOfNodes();
  int offset3=offset2+((int)addCoords.size())/2;
  m1->getBoundingBoxForBBTree(bbox1);
  m2->getBoundingBoxForBBTree(bbox2);
  BBTree<SPACEDIM,int> myTree(&bbox2[0],0,0,m2->getNumberOfCells(),eps);
  int ncell1=m1->getNumberOfCells();
  crI.push_back(0);
  for(int i=0;i<ncell1;i++)
    {
      std::vector<int> candidates2;
      myTree.getIntersectingElems(&bbox1[i*2*SPACEDIM],candidates2);
      std::map<INTERP_KERNEL::Node *,int> mapp;
      std::map<int,INTERP_KERNEL::Node *> mappRev;
      INTERP_KERNEL::QuadraticPolygon pol1;
      INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType)conn1[connI1[i]];
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(typ);
      MEDCouplingUMeshBuildQPFromMesh3(coo1,offset1,coo2,offset2,addCoords,desc1+descIndx1[i],desc1+descIndx1[i+1],intesctEdges1,/* output */mapp,mappRev);
      pol1.buildFromCrudeDataArray(mappRev,cm.isQuadratic(),conn1+connI1[i]+1,coo1,
                                   desc1+descIndx1[i],desc1+descIndx1[i+1],intesctEdges1);
      std::vector<int> crTmp,crITmp;
      crITmp.push_back(crI.back());
      for(std::vector<int>::const_iterator it2=candidates2.begin();it2!=candidates2.end();it2++)
        {
          INTERP_KERNEL::QuadraticPolygon pol2;
          pol1.initLocations();
          MEDCouplingUMeshBuildQPFromMesh3(coo1,offset1,coo2,offset2,addCoords,desc2+descIndx2[*it2],desc2+descIndx2[*it2+1],intesctEdges2,/* output */mapp,mappRev);
          INTERP_KERNEL::NormalizedCellType typ2=(INTERP_KERNEL::NormalizedCellType)conn2[connI2[*it2]];
          const INTERP_KERNEL::CellModel& cm2=INTERP_KERNEL::CellModel::GetCellModel(typ2);
          pol2.buildFromCrudeDataArray2(mappRev,cm2.isQuadratic(),conn2+connI2[*it2]+1,coo2,desc2+descIndx2[*it2],desc2+descIndx2[*it2+1],intesctEdges2,
                                        pol1,desc1+descIndx1[i],desc1+descIndx1[i+1],intesctEdges1,colinear2);
          //MEDCouplingUMeshAssignOnLoc(pol1,pol2,desc1+descIndx1[i],desc1+descIndx1[i+1],intesctEdges1,desc2+descIndx2[*it2],desc2+descIndx2[*it2+1],intesctEdges2,colinear2);
          pol1.buildPartitionsAbs(pol2,mapp,i,*it2,offset3,addCoordsQuadratic,cr,crI,cNb1,cNb2);
        }
      if(!crTmp.empty())
        {
          cr.insert(cr.end(),crTmp.begin(),crTmp.end());
          crI.insert(crI.end(),crITmp.begin()+1,crITmp.end());
        }
      for(std::map<int,INTERP_KERNEL::Node *>::const_iterator it=mappRev.begin();it!=mappRev.end();it++)
        (*it).second->decrRef();
    }
}

/*!
 * This method is private and is the first step of Partition of 2D mesh (spaceDim==2 and meshDim==2).
 * 
 */
void MEDCouplingUMesh::IntersectDescending2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps,
                                                   std::vector< std::vector<int> >& intersectEdge1, std::vector< std::vector<int> >& colinear2, std::vector< std::vector<int> >& subDiv2,
                                                   MEDCouplingUMesh *& m1Desc, DataArrayInt *&desc1, DataArrayInt *&descIndx1, DataArrayInt *&revDesc1, DataArrayInt *&revDescIndx1,
                                                   MEDCouplingUMesh *& m2Desc, DataArrayInt *&desc2, DataArrayInt *&descIndx2, DataArrayInt *&revDesc2, DataArrayInt *&revDescIndx2,
                                                   std::vector<double>& addCoo) throw(INTERP_KERNEL::Exception)
{
  static const int SPACEDIM=2;
  desc1=DataArrayInt::New(); descIndx1=DataArrayInt::New(); revDesc1=DataArrayInt::New(); revDescIndx1=DataArrayInt::New();
  desc2=DataArrayInt::New();
  descIndx2=DataArrayInt::New();
  revDesc2=DataArrayInt::New();
  revDescIndx2=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> dd1(desc1),dd2(descIndx1),dd3(revDesc1),dd4(revDescIndx1);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> dd5(desc2),dd6(descIndx2),dd7(revDesc2),dd8(revDescIndx2);
  m1Desc=m1->buildDescendingConnectivity2(desc1,descIndx1,revDesc1,revDescIndx1);
  m2Desc=m2->buildDescendingConnectivity2(desc2,descIndx2,revDesc2,revDescIndx2);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> dd9(m1Desc),dd10(m2Desc);
  const int *c1=m1Desc->getNodalConnectivity()->getConstPointer();
  const int *ci1=m1Desc->getNodalConnectivityIndex()->getConstPointer();
  std::vector<double> bbox1,bbox2;
  m1Desc->getBoundingBoxForBBTree(bbox1);
  m2Desc->getBoundingBoxForBBTree(bbox2);
  int ncell1=m1Desc->getNumberOfCells();
  int ncell2=m2Desc->getNumberOfCells();
  intersectEdge1.resize(ncell1);
  colinear2.resize(ncell2);
  subDiv2.resize(ncell2);
  BBTree<SPACEDIM,int> myTree(&bbox2[0],0,0,m2Desc->getNumberOfCells(),-eps);
  std::vector<int> candidates1(1);
  int offset1=m1->getNumberOfNodes();
  int offset2=offset1+m2->getNumberOfNodes();
  for(int i=0;i<ncell1;i++)
    {
      std::vector<int> candidates2;
      myTree.getIntersectingElems(&bbox1[i*2*SPACEDIM],candidates2);
      if(!candidates2.empty())
        {
          std::map<INTERP_KERNEL::Node *,int> map1,map2;
          INTERP_KERNEL::QuadraticPolygon *pol2=MEDCouplingUMeshBuildQPFromMesh(m2Desc,candidates2,map2);
          candidates1[0]=i;
          INTERP_KERNEL::QuadraticPolygon *pol1=MEDCouplingUMeshBuildQPFromMesh(m1Desc,candidates1,map1);
          pol1->splitAbs(*pol2,map1,map2,offset1,offset2,candidates2,intersectEdge1[i],i,colinear2,subDiv2,addCoo);
          delete pol2;
          delete pol1;
        }
      else
        intersectEdge1[i].insert(intersectEdge1[i].end(),c1+ci1[i]+1,c1+ci1[i+1]);
    }
  m1Desc->incrRef(); desc1->incrRef(); descIndx1->incrRef(); revDesc1->incrRef(); revDescIndx1->incrRef();
  m2Desc->incrRef(); desc2->incrRef(); descIndx2->incrRef(); revDesc2->incrRef(); revDescIndx2->incrRef();
}

/*!
 * This method performs the 2nd step of Partition of 2D mesh.
 * This method has 4 inputs :
 *  - a mesh 'm1' with meshDim==1 and a SpaceDim==2
 *  - a mesh 'm2' with meshDim==1 and a SpaceDim==2
 *  - subDiv of size 'm2->getNumberOfCells()' that lists for each seg cell in 'm' the splitting node ids in randomly sorted.
 * The aim of this method is to sort the splitting nodes, if any, and to put in 'intersectEdge' output paramter based on edges of mesh 'm2'
 * @param m1 is expected to be a mesh of meshDimension equal to 1 and spaceDim equal to 2. No check of that is performed by this method. Only present for its coords in case of 'subDiv' shares some nodes of 'm1'
 * @param m2 is expected to be a mesh of meshDimension equal to 1 and spaceDim equal to 2. No check of that is performed by this method.
 * @param addCoo input parameter with additionnal nodes linked to intersection of the 2 meshes.
 */
void MEDCouplingUMesh::BuildIntersectEdges(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, const std::vector<double>& addCoo, const std::vector< std::vector<int> >& subDiv, std::vector< std::vector<int> >& intersectEdge) throw(INTERP_KERNEL::Exception)
{
  int offset1=m1->getNumberOfNodes();
  int ncell=m2->getNumberOfCells();
  const int *c=m2->getNodalConnectivity()->getConstPointer();
  const int *cI=m2->getNodalConnectivityIndex()->getConstPointer();
  const double *coo=m2->getCoords()->getConstPointer();
  const double *cooBis=m1->getCoords()->getConstPointer();
  int offset2=offset1+m2->getNumberOfNodes();
  intersectEdge.resize(ncell);
  for(int i=0;i<ncell;i++,cI++)
    {
      const std::vector<int>& divs=subDiv[i];
      int nnode=cI[1]-cI[0]-1;
      std::map<int, std::pair<INTERP_KERNEL::Node *,bool> > mapp2;
      std::map<INTERP_KERNEL::Node *, int> mapp22;
      for(int j=0;j<nnode;j++)
        {
          INTERP_KERNEL::Node *nn=new INTERP_KERNEL::Node(coo[2*c[(*cI)+j+1]],coo[2*c[(*cI)+j+1]+1]);
          int nnid=c[(*cI)+j+1];
          mapp2[nnid]=std::pair<INTERP_KERNEL::Node *,bool>(nn,true);
          mapp22[nn]=nnid+offset1;
        }
      INTERP_KERNEL::Edge *e=MEDCouplingUMeshBuildQPFromEdge((INTERP_KERNEL::NormalizedCellType)c[*cI],mapp2,c+(*cI)+1);
      for(std::map<int, std::pair<INTERP_KERNEL::Node *,bool> >::const_iterator it=mapp2.begin();it!=mapp2.end();it++)
        ((*it).second.first)->decrRef();
      std::vector<INTERP_KERNEL::Node *> addNodes(divs.size());
      std::map<INTERP_KERNEL::Node *,int> mapp3;
      for(std::size_t j=0;j<divs.size();j++)
        {
          int id=divs[j];
          INTERP_KERNEL::Node *tmp=0;
          if(id<offset1)
            tmp=new INTERP_KERNEL::Node(cooBis[2*id],cooBis[2*id+1]);
          else if(id<offset2)
            tmp=new INTERP_KERNEL::Node(coo[2*(id-offset1)],coo[2*(id-offset1)+1]);//if it happens, bad news mesh 'm2' is non conform.
          else
            tmp=new INTERP_KERNEL::Node(addCoo[2*(id-offset2)],addCoo[2*(id-offset2)+1]);
          addNodes[j]=tmp;
          mapp3[tmp]=id;
        }
      e->sortIdsAbs(addNodes,mapp22,mapp3,intersectEdge[i]);
      for(std::vector<INTERP_KERNEL::Node *>::const_iterator it=addNodes.begin();it!=addNodes.end();it++)
        (*it)->decrRef();
      e->decrRef();
    }
}

/*!
 * This method is part of the Slice3D algorithm. It is the first step of assembly process, ones coordinates have been computed (by MEDCouplingUMesh::split3DCurveWithPlane method).
 * This method allows to compute given the status of 3D curve cells and the descending connectivity 3DSurf->3DCurve to deduce the intersection of each 3D surf cells
 * with a plane. The result will be put in 'cut3DSuf' out parameter.
 * @param cut3DCurve  input paramter that gives for each 3DCurve cell if it owns fully to the plane or partially.
 * @param nodesOnPlane, returns all the nodes that are on the plane.
 * @param nodal3DSurf is the nodal connectivity of 3D surf mesh.
 * @param nodalIndx3DSurf is the nodal connectivity index of 3D surf mesh.
 * @param nodal3DCurve is the nodal connectivity of 3D curve mesh.
 * @param nodal3DIndxCurve is the nodal connectivity index of 3D curve mesh.
 * @param desc is the descending connectivity 3DSurf->3DCurve
 * @param descIndx is the descending connectivity index 3DSurf->3DCurve
 * @param cut3DSuf input/output param.
 */
void MEDCouplingUMesh::AssemblyForSplitFrom3DCurve(const std::vector<int>& cut3DCurve, std::vector<int>& nodesOnPlane, const int *nodal3DSurf, const int *nodalIndx3DSurf,
                                                   const int *nodal3DCurve, const int *nodalIndx3DCurve,
                                                   const int *desc, const int *descIndx, 
                                                   std::vector< std::pair<int,int> >& cut3DSurf) throw(INTERP_KERNEL::Exception)
{
  std::set<int> nodesOnP(nodesOnPlane.begin(),nodesOnPlane.end());
  int nbOf3DSurfCell=(int)cut3DSurf.size();
  for(int i=0;i<nbOf3DSurfCell;i++)
    {
      std::vector<int> res;
      int offset=descIndx[i];
      int nbOfSeg=descIndx[i+1]-offset;
      for(int j=0;j<nbOfSeg;j++)
        {
          int edgeId=desc[offset+j];
          int status=cut3DCurve[edgeId];
          if(status!=-2)
            {
              if(status>-1)
                res.push_back(status);
              else
                {
                  res.push_back(nodal3DCurve[nodalIndx3DCurve[edgeId]+1]);
                  res.push_back(nodal3DCurve[nodalIndx3DCurve[edgeId]+2]);
                }
            }
        }
      switch(res.size())
        {
        case 2:
          {
            cut3DSurf[i].first=res[0]; cut3DSurf[i].second=res[1];
            break;
          }
        case 1:
        case 0:
          {
            std::set<int> s1(nodal3DSurf+nodalIndx3DSurf[i]+1,nodal3DSurf+nodalIndx3DSurf[i+1]);
            std::set_intersection(nodesOnP.begin(),nodesOnP.end(),s1.begin(),s1.end(),std::back_insert_iterator< std::vector<int> >(res));
            if(res.size()==2)
              {
                cut3DSurf[i].first=res[0]; cut3DSurf[i].second=res[1];
              }
            else
              {
                cut3DSurf[i].first=-1; cut3DSurf[i].second=-1;
              }
            break;
          }
        default:
          {// case when plane is on a multi colinear edge of a polyhedron
            if((int)res.size()==2*nbOfSeg)
              {
                cut3DSurf[i].first=-2; cut3DSurf[i].second=i;
              }
            else
              throw INTERP_KERNEL::Exception("MEDCouplingUMesh::AssemblyPointsFrom3DCurve : unexpected situation !");
          }
        }
    }
}

/*!
 * 'this' is expected to be a mesh with spaceDim==3 and meshDim==3. If not an exception will be thrown.
 * This method is part of the Slice3D algorithm. It is the second step of assembly process, ones coordinates have been computed (by MEDCouplingUMesh::split3DCurveWithPlane method).
 * This method allows to compute given the result of 3D surf cells with plane and the descending connectivity 3D->3DSurf to deduce the intersection of each 3D cells
 * with a plane. The result will be put in 'nodalRes' 'nodalResIndx' and 'cellIds' out parameters.
 * @param cut3DSurf  input paramter that gives for each 3DSurf its intersection with plane (result of MEDCouplingUMesh::AssemblyForSplitFrom3DCurve).
 * @param desc is the descending connectivity 3D->3DSurf
 * @param descIndx is the descending connectivity index 3D->3DSurf
 */
void MEDCouplingUMesh::assemblyForSplitFrom3DSurf(const std::vector< std::pair<int,int> >& cut3DSurf,
                                                  const int *desc, const int *descIndx,
                                                  std::vector<int>& nodalRes, std::vector<int>& nodalResIndx, std::vector<int>& cellIds) const throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::assemblyForSplitFrom3DSurf works on umeshes with meshdim equal to 3 and spaceDim equal to 3 too!");
  const int *nodal3D=_nodal_connec->getConstPointer();
  const int *nodalIndx3D=_nodal_connec_index->getConstPointer();
  int nbOfCells=getNumberOfCells();
  for(int i=0;i<nbOfCells;i++)
    {
      std::map<int, std::set<int> > m;
      int offset=descIndx[i];
      int nbOfFaces=descIndx[i+1]-offset;
      int start=-1;
      int end=-1;
      for(int j=0;j<nbOfFaces;j++)
        {
          const std::pair<int,int>& p=cut3DSurf[desc[offset+j]];
          if(p.first!=-1 && p.second!=-1)
            {
              if(p.first!=-2)
                {
                  start=p.first; end=p.second;
                  m[p.first].insert(p.second);
                  m[p.second].insert(p.first);
                }
              else
                {
                  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)nodal3D[nodalIndx3D[i]]);
                  int sz=nodalIndx3D[i+1]-nodalIndx3D[i]-1;
                  INTERP_KERNEL::AutoPtr<int> tmp=new int[sz];
                  INTERP_KERNEL::NormalizedCellType cmsId;
                  unsigned nbOfNodesSon=cm.fillSonCellNodalConnectivity2(j,nodal3D+nodalIndx3D[i]+1,sz,tmp,cmsId);
                  start=tmp[0]; end=tmp[nbOfNodesSon-1];
                  for(unsigned k=0;k<nbOfNodesSon;k++)
                    {
                      m[tmp[k]].insert(tmp[(k+1)%nbOfNodesSon]);
                      m[tmp[(k+1)%nbOfNodesSon]].insert(tmp[k]);
                    }
                }
            }
        }
      if(m.empty())
        continue;
      std::vector<int> conn(1,(int)INTERP_KERNEL::NORM_POLYGON);
      int prev=end;
      while(end!=start)
        {
          std::map<int, std::set<int> >::const_iterator it=m.find(start);
          const std::set<int>& s=(*it).second;
          std::set<int> s2; s2.insert(prev);
          std::set<int> s3;
          std::set_difference(s.begin(),s.end(),s2.begin(),s2.end(),inserter(s3,s3.begin()));
          if(s3.size()==1)
            {
              int val=*s3.begin();
              conn.push_back(start);
              prev=start;
              start=val;
            }
          else
            start=end;
        }
      conn.push_back(end);
      if(conn.size()>3)
        {
          nodalRes.insert(nodalRes.end(),conn.begin(),conn.end());
          nodalResIndx.push_back((int)nodalRes.size());
          cellIds.push_back(i);
        }
    }
}

/*!
 * This method compute the convex hull of a single 2D cell. This method tries to conserve at maximum the given input connectivity. In particular, if the orientation of cell is not clockwise
 * as in MED format norm. If definitely the result of Jarvis algorithm is not matchable with the input connectivity, the result will be copied into \b nodalConnecOut parameter and
 * the geometric cell type set to INTERP_KERNEL::NORM_POLYGON.
 * This method excepts that \b coords parameter is expected to be in dimension 2. [\b nodalConnBg, \b nodalConnEnd) is the nodal connectivity of the input
 * cell (geometric cell type included at the position 0). If the meshdimension of the input cell is not equal to 2 an INTERP_KERNEL::Exception will be thrown.
 * 
 * @return false if the input connectivity represents already the convex hull, true if the input cell needs to be reordered.
 */
bool MEDCouplingUMesh::BuildConvexEnvelopOf2DCellJarvis(const double *coords, const int *nodalConnBg, const int *nodalConnEnd, std::vector<int>& nodalConnecOut) throw(INTERP_KERNEL::Exception)
{
  std::size_t sz=std::distance(nodalConnBg,nodalConnEnd);
  if(sz>=4)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)*nodalConnBg);
      if(cm.getDimension()==2)
        {
          const int *node=nodalConnBg+1;
          int startNode=*node++;
          double refX=coords[2*startNode];
          for(;node!=nodalConnEnd;node++)
            {
              if(coords[2*(*node)]<refX)
                {
                  startNode=*node;
                  refX=coords[2*startNode];
                }
            }
          std::vector<int> tmpOut; tmpOut.reserve(sz); tmpOut.push_back(startNode);
          refX=1e300;
          double tmp1;
          double tmp2[2];
          double angle0=-M_PI/2;
          //
          int nextNode=-1;
          int prevNode=-1;
          double resRef;
          double angleNext;
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
          std::vector<int> tmp3(2*(sz-1));
          std::vector<int>::iterator it=std::copy(nodalConnBg+1,nodalConnEnd,tmp3.begin());
          std::copy(nodalConnBg+1,nodalConnEnd,it);
          if(std::search(tmp3.begin(),tmp3.end(),tmpOut.begin(),tmpOut.end())!=tmp3.end())
            {
              nodalConnecOut.insert(nodalConnecOut.end(),nodalConnBg,nodalConnEnd);
              return false;
            }
          if(std::search(tmp3.rbegin(),tmp3.rend(),tmpOut.begin(),tmpOut.end())!=tmp3.rend())
            {
              nodalConnecOut.insert(nodalConnecOut.end(),nodalConnBg,nodalConnEnd);
              return false;
            }
          else
            {
              nodalConnecOut.push_back((int)INTERP_KERNEL::NORM_POLYGON);
              nodalConnecOut.insert(nodalConnecOut.end(),tmpOut.begin(),tmpOut.end());
              return true;
            }
        }
    }
  else
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::BuildConvexEnvelopOf2DCellJarvis : invalid 2D cell connectivity !");
}

/*!
 * This method works on an input pair (\b arr, \b arrIndx) where \b arr indexes is in \b arrIndx.
 * This method will not impact the size of inout parameter \b arrIndx but the size of \b arr will be modified in case of suppression.
 * 
 * \param [in] idsToRemoveBg begin of set of ids to remove in \b arr (included)
 * \param [in] idsToRemoveEnd end of set of ids to remove in \b arr (excluded)
 * \param [in,out] arr array in which the remove operation will be done.
 * \param [in,out] arrIndx array in the remove operation will modify
 * \param [in] offsetForRemoval (by default 0) offset so that for each i in [0,arrIndx->getNumberOfTuples()-1) removal process will be performed in the following range [arr+arrIndx[i]+offsetForRemoval,arr+arr[i+1])
 * \return true if \b arr and \b arrIndx have been modified, false if not.
 */
bool MEDCouplingUMesh::RemoveIdsFromIndexedArrays(const int *idsToRemoveBg, const int *idsToRemoveEnd, DataArrayInt *arr, DataArrayInt *arrIndx, int offsetForRemoval) throw(INTERP_KERNEL::Exception)
{
  if(!arrIndx || !arr)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::RemoveIdsFromIndexedArrays : some input arrays are empty !");
  if(offsetForRemoval<0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::RemoveIdsFromIndexedArrays : offsetForRemoval should be >=0 !");
  std::set<int> s(idsToRemoveBg,idsToRemoveEnd);
  int nbOfGrps=arrIndx->getNumberOfTuples()-1;
  int *arrIPtr=arrIndx->getPointer();
  *arrIPtr++=0;
  int previousArrI=0;
  const int *arrPtr=arr->getConstPointer();
  std::vector<int> arrOut;
  for(int i=0;i<nbOfGrps;i++,arrIPtr++)
    {
      if(*arrIPtr-previousArrI>offsetForRemoval)
        {
          for(const int *work=arrPtr+previousArrI+offsetForRemoval;work!=arrPtr+*arrIPtr;work++)
            {
              if(s.find(*work)==s.end())
                arrOut.push_back(*work);
            }
        }
      previousArrI=*arrIPtr;
      *arrIPtr=(int)arrOut.size();
    }
  if(arr->getNumberOfTuples()==(int)arrOut.size())
    return false;
  arr->alloc((int)arrOut.size(),1);
  std::copy(arrOut.begin(),arrOut.end(),arr->getPointer());
  return true;
}

/*!
 * This method works on a pair input (\b arrIn, \b arrIndxIn) where \b arrIn indexes is in \b arrIndxIn.
 * This method returns the result of the extraction ( specified by a set of ids in [\b idsOfSelectBg , \b idsOfSelectEnd ) ).
 * The selection of extraction is done standardly in new2old format.
 * This method returns indexed arrays using 2 arrays (arrOut,arrIndexOut).
 *
 * \param [in] idsOfSelectBg begin of set of ids of the input extraction (included)
 * \param [in] idsOfSelectEnd end of set of ids of the input extraction (excluded)
 * \param [in] arrIn arr origin array from which the extraction will be done.
 * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
 * \param [out] arrOut the resulting array
 * \param [out] arrIndexOut the index array of the resulting array \b arrOut
 */
void MEDCouplingUMesh::ExtractFromIndexedArrays(const int *idsOfSelectBg, const int *idsOfSelectEnd, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn,
                                                DataArrayInt* &arrOut, DataArrayInt* &arrIndexOut) throw(INTERP_KERNEL::Exception)
{
  if(!arrIn || !arrIndxIn)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ExtractFromIndexedArrays : input pointer is NULL !");
  std::size_t sz=std::distance(idsOfSelectBg,idsOfSelectEnd);
  const int *arrInPtr=arrIn->getConstPointer();
  const int *arrIndxPtr=arrIndxIn->getConstPointer();
  int nbOfGrps=arrIndxIn->getNumberOfTuples()-1;
  int maxSizeOfArr=arrIn->getNumberOfTuples();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arro=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arrIo=DataArrayInt::New();
  arrIo->alloc((int)(sz+1),1);
  const int *idsIt=idsOfSelectBg;
  int *work=arrIo->getPointer();
  *work++=0;
  int lgth=0;
  for(std::size_t i=0;i<sz;i++,work++,idsIt++)
    {
      if(*idsIt>=0 && *idsIt<nbOfGrps)
        lgth+=arrIndxPtr[*idsIt+1]-arrIndxPtr[*idsIt];
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::ExtractFromIndexedArrays : id located on pos #" << i << " value is " << *idsIt << " ! Must be in [0," << nbOfGrps << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      if(lgth>=work[-1])
        *work=lgth;
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::ExtractFromIndexedArrays : id located on pos #" << i << " value is " << *idsIt << " and at this pos arrIndxIn[" << *idsIt;
          oss << "+1]-arrIndxIn[" << *idsIt << "] < 0 ! The input index array is bugged !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  arro->alloc(lgth,1);
  work=arro->getPointer();
  idsIt=idsOfSelectBg;
  for(std::size_t i=0;i<sz;i++,idsIt++)
    {
      if(arrIndxPtr[*idsIt]>=0 && arrIndxPtr[*idsIt+1]<=maxSizeOfArr)
        work=std::copy(arrInPtr+arrIndxPtr[*idsIt],arrInPtr+arrIndxPtr[*idsIt+1],work);
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::ExtractFromIndexedArrays : id located on pos #" << i << " value is " << *idsIt << " arrIndx[" << *idsIt << "] must be >= 0 and arrIndx[";
          oss << *idsIt << "+1] <= " << maxSizeOfArr << " (the size of arrIn)!";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  arrOut=arro;
  arrIndexOut=arrIo;
  arro->incrRef();
  arrIo->incrRef();
}

/*!
 * This method works on an input pair (\b arrIn, \b arrIndxIn) where \b arrIn indexes is in \b arrIndxIn.
 * This method builds an output pair (\b arrOut,\b arrIndexOut) that is a copy from \b arrIn for all cell ids \b not \b in [\b idsOfSelectBg, \b idsOfSelectEnd) and for
 * cellIds \b in [\b idsOfSelectBg, \b idsOfSelectEnd) a copy coming from the corresponding values in input pair (\b srcArr, \b srcArrIndex).
 * This method is an generalization of MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx that performs the same thing but by without building explicitely a result output arrays.
 *
 * \param [in] idsOfSelectBg begin of set of ids of the input extraction (included)
 * \param [in] idsOfSelectEnd end of set of ids of the input extraction (excluded)
 * \param [in] arrIn arr origin array from which the extraction will be done.
 * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
 * \param [in] srcArr input array that will be used as source of copy for ids in [\b idsOfSelectBg, \b idsOfSelectEnd)
 * \param [in] srcArrIndex index array of \b srcArr
 * \param [out] arrOut the resulting array
 * \param [out] arrIndexOut the index array of the resulting array \b arrOut
 * 
 * \sa MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx
 */
void MEDCouplingUMesh::SetPartOfIndexedArrays(const int *idsOfSelectBg, const int *idsOfSelectEnd, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn,
                                              const DataArrayInt *srcArr, const DataArrayInt *srcArrIndex,
                                              DataArrayInt* &arrOut, DataArrayInt* &arrIndexOut) throw(INTERP_KERNEL::Exception)
{
  if(arrIn==0 || arrIndxIn==0 || srcArr==0 || srcArrIndex==0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::SetPartOfIndexedArrays : presence of null pointer in input parameter !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arro=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arrIo=DataArrayInt::New();
  int nbOfTuples=arrIndxIn->getNumberOfTuples()-1;
  std::vector<bool> v(nbOfTuples,true);
  int offset=0;
  const int *arrIndxInPtr=arrIndxIn->getConstPointer();
  const int *srcArrIndexPtr=srcArrIndex->getConstPointer();
  for(const int *it=idsOfSelectBg;it!=idsOfSelectEnd;it++,srcArrIndexPtr++)
    {
      if(*it>=0 && *it<nbOfTuples)
        {
          v[*it]=false;
          offset+=(srcArrIndexPtr[1]-srcArrIndexPtr[0])-(arrIndxInPtr[*it+1]-arrIndxInPtr[*it]);
        }
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::SetPartOfIndexedArrays : On pos #" << std::distance(idsOfSelectBg,it) << " value is " << *it << " not in [0," << nbOfTuples << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  srcArrIndexPtr=srcArrIndex->getConstPointer();
  arrIo->alloc(nbOfTuples+1,1);
  arro->alloc(arrIn->getNumberOfTuples()+offset,1);
  const int *arrInPtr=arrIn->getConstPointer();
  const int *srcArrPtr=srcArr->getConstPointer();
  int *arrIoPtr=arrIo->getPointer(); *arrIoPtr++=0;
  int *arroPtr=arro->getPointer();
  for(int ii=0;ii<nbOfTuples;ii++,arrIoPtr++)
    {
      if(v[ii])
        {
          arroPtr=std::copy(arrInPtr+arrIndxInPtr[ii],arrInPtr+arrIndxInPtr[ii+1],arroPtr);
          *arrIoPtr=arrIoPtr[-1]+(arrIndxInPtr[ii+1]-arrIndxInPtr[ii]);
        }
      else
        {
          std::size_t pos=std::distance(idsOfSelectBg,std::find(idsOfSelectBg,idsOfSelectEnd,ii));
          arroPtr=std::copy(srcArrPtr+srcArrIndexPtr[pos],srcArrPtr+srcArrIndexPtr[pos+1],arroPtr);
          *arrIoPtr=arrIoPtr[-1]+(srcArrIndexPtr[pos+1]-srcArrIndexPtr[pos]);
        }
    }
  arrOut=arro; arro->incrRef();
  arrIndexOut=arrIo; arrIo->incrRef();
}

/*!
 * This method works on an input pair (\b arrIn, \b arrIndxIn) where \b arrIn indexes is in \b arrIndxIn.
 * This method is an specialization of MEDCouplingUMesh::SetPartOfIndexedArrays in the case of assignement do not modify the index in \b arrIndxIn.
 *
 * \param [in] idsOfSelectBg begin of set of ids of the input extraction (included)
 * \param [in] idsOfSelectEnd end of set of ids of the input extraction (excluded)
 * \param [in,out] arrInOut arr origin array from which the extraction will be done.
 * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
 * \param [in] srcArr input array that will be used as source of copy for ids in [\b idsOfSelectBg, \b idsOfSelectEnd)
 * \param [in] srcArrIndex index array of \b srcArr
 * 
 * \sa MEDCouplingUMesh::SetPartOfIndexedArrays
 */
void MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx(const int *idsOfSelectBg, const int *idsOfSelectEnd, DataArrayInt *arrInOut, const DataArrayInt *arrIndxIn,
                                                     const DataArrayInt *srcArr, const DataArrayInt *srcArrIndex) throw(INTERP_KERNEL::Exception)
{
  if(arrInOut==0 || arrIndxIn==0 || srcArr==0 || srcArrIndex==0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx : presence of null pointer in input parameter !");
  int nbOfTuples=arrIndxIn->getNumberOfTuples()-1;
  const int *arrIndxInPtr=arrIndxIn->getConstPointer();
  const int *srcArrIndexPtr=srcArrIndex->getConstPointer();
  int *arrInOutPtr=arrInOut->getPointer();
  const int *srcArrPtr=srcArr->getConstPointer();
  for(const int *it=idsOfSelectBg;it!=idsOfSelectEnd;it++,srcArrIndexPtr++)
    {
      if(*it>=0 && *it<nbOfTuples)
        {
          if(srcArrIndexPtr[1]-srcArrIndexPtr[0]==arrIndxInPtr[*it+1]-arrIndxInPtr[*it])
            std::copy(srcArrPtr+srcArrIndexPtr[0],srcArrPtr+srcArrIndexPtr[1],arrInOutPtr+arrIndxInPtr[*it]);
          else
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx : On pos #" << std::distance(idsOfSelectBg,it) << " id (idsOfSelectBg[" << std::distance(idsOfSelectBg,it)<< "]) is " << *it << " arrIndxIn[id+1]-arrIndxIn[id]!=srcArrIndex[pos+1]-srcArrIndex[pos] !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx : On pos #" << std::distance(idsOfSelectBg,it) << " value is " << *it << " not in [0," << nbOfTuples << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}

/*!
 * This method works on a pair input (\b arrIn, \b arrIndxIn) where \b arr indexes is in \b arrIndxIn.
 * This method expects that these two input arrays come from the output of MEDCouplingUMesh::computeNeighborsOfCells method.
 * This method start from id 0 that will be contained in output DataArrayInt. It searches then all neighbors of id0 regarding arrIn[arrIndxIn[0]:arrIndxIn[0+1]].
 * Then it is repeated recursively until either all ids are fetched or no more ids are reachable step by step.
 * A negative value in \b arrIn means that it is ignored.
 * This method is usefull to see if a mesh is contiguous regarding its connectivity. If it is not the case the size of returned array is different from arrIndxIn->getNumberOfTuples()-1.
 * 
 * \param [in] arrIn arr origin array from which the extraction will be done.
 * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
 * \return a newly allocated DataArray that stores all ids fetched by the gradually spread process.
 */
DataArrayInt *MEDCouplingUMesh::ComputeSpreadZoneGradually(const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn) throw(INTERP_KERNEL::Exception)
{
  if(!arrIn || !arrIndxIn)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ExtractFromIndexedArrays : input pointer is NULL !");
  int nbOfTuples=arrIndxIn->getNumberOfTuples()-1;
  if(nbOfTuples<=0)
    {
      DataArrayInt *ret=DataArrayInt::New(); ret->alloc(0,1);
      return ret;
    }
  const int *arrInPtr=arrIn->getConstPointer();
  const int *arrIndxPtr=arrIndxIn->getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arro=DataArrayInt::New();
  arro->alloc(nbOfTuples,1);
  arro->fillWithValue(-1);
  int *arroPtr=arro->getPointer();
  std::set<int> s; s.insert(0);
  while(!s.empty())
    {
      std::set<int> s2;
      for(std::set<int>::const_iterator it=s.begin();it!=s.end();it++)
        {
          for(const int *work=arrInPtr+arrIndxPtr[*it];work!=arrInPtr+arrIndxPtr[*it+1];work++)
            {
              if(*work>=0 && arroPtr[*work]<0)
                {
                  arroPtr[*work]=1;
                  s2.insert(*work);
                }
            }
        }
      s=s2;
    }
  return arro->getIdsEqual(1);
}

/*!
 * This method works on an input pair (\b arrIn, \b arrIndxIn) where \b arrIn indexes is in \b arrIndxIn.
 * This method builds an output pair (\b arrOut,\b arrIndexOut) that is a copy from \b arrIn for all cell ids \b not \b in [\b idsOfSelectBg, \b idsOfSelectEnd) and for
 * cellIds \b in [\b idsOfSelectBg, \b idsOfSelectEnd) a copy coming from the corresponding values in input pair (\b srcArr, \b srcArrIndex).
 * This method is an generalization of MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx that performs the same thing but by without building explicitely a result output arrays.
 *
 * \param [in] start begin of set of ids of the input extraction (included)
 * \param [in] end end of set of ids of the input extraction (excluded)
 * \param [in] step step of the set of ids in range mode.
 * \param [in] arrIn arr origin array from which the extraction will be done.
 * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
 * \param [in] srcArr input array that will be used as source of copy for ids in [\b idsOfSelectBg, \b idsOfSelectEnd)
 * \param [in] srcArrIndex index array of \b srcArr
 * \param [out] arrOut the resulting array
 * \param [out] arrIndexOut the index array of the resulting array \b arrOut
 * 
 * \sa MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx MEDCouplingUMesh::SetPartOfIndexedArrays
 */
void MEDCouplingUMesh::SetPartOfIndexedArrays2(int start, int end, int step, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn,
                                               const DataArrayInt *srcArr, const DataArrayInt *srcArrIndex,
                                               DataArrayInt* &arrOut, DataArrayInt* &arrIndexOut) throw(INTERP_KERNEL::Exception)
{
  if(arrIn==0 || arrIndxIn==0 || srcArr==0 || srcArrIndex==0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::SetPartOfIndexedArrays2 : presence of null pointer in input parameter !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arro=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arrIo=DataArrayInt::New();
  int nbOfTuples=arrIndxIn->getNumberOfTuples()-1;
  int offset=0;
  const int *arrIndxInPtr=arrIndxIn->getConstPointer();
  const int *srcArrIndexPtr=srcArrIndex->getConstPointer();
  int nbOfElemsToSet=DataArray::GetNumberOfItemGivenBESRelative(start,end,step,"MEDCouplingUMesh::SetPartOfIndexedArrays2 : ");
  int it=start;
  for(int i=0;i<nbOfElemsToSet;i++,srcArrIndexPtr++,it+=step)
    {
      if(it>=0 && it<nbOfTuples)
        offset+=(srcArrIndexPtr[1]-srcArrIndexPtr[0])-(arrIndxInPtr[it+1]-arrIndxInPtr[it]);
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::SetPartOfIndexedArrays2 : On pos #" << i << " value is " << it << " not in [0," << nbOfTuples << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  srcArrIndexPtr=srcArrIndex->getConstPointer();
  arrIo->alloc(nbOfTuples+1,1);
  arro->alloc(arrIn->getNumberOfTuples()+offset,1);
  const int *arrInPtr=arrIn->getConstPointer();
  const int *srcArrPtr=srcArr->getConstPointer();
  int *arrIoPtr=arrIo->getPointer(); *arrIoPtr++=0;
  int *arroPtr=arro->getPointer();
  for(int ii=0;ii<nbOfTuples;ii++,arrIoPtr++)
    {
      int pos=DataArray::GetPosOfItemGivenBESRelativeNoThrow(ii,start,end,step);
      if(pos<0)
        {
          arroPtr=std::copy(arrInPtr+arrIndxInPtr[ii],arrInPtr+arrIndxInPtr[ii+1],arroPtr);
          *arrIoPtr=arrIoPtr[-1]+(arrIndxInPtr[ii+1]-arrIndxInPtr[ii]);
        }
      else
        {
          arroPtr=std::copy(srcArrPtr+srcArrIndexPtr[pos],srcArrPtr+srcArrIndexPtr[pos+1],arroPtr);
          *arrIoPtr=arrIoPtr[-1]+(srcArrIndexPtr[pos+1]-srcArrIndexPtr[pos]);
        }
    }
  arrOut=arro; arro->incrRef();
  arrIndexOut=arrIo; arrIo->incrRef();
}

/*!
 * This method works on an input pair (\b arrIn, \b arrIndxIn) where \b arrIn indexes is in \b arrIndxIn.
 * This method is an specialization of MEDCouplingUMesh::SetPartOfIndexedArrays in the case of assignement do not modify the index in \b arrIndxIn.
 *
 * \param [in] start begin of set of ids of the input extraction (included)
 * \param [in] end end of set of ids of the input extraction (excluded)
 * \param [in] step step of the set of ids in range mode.
 * \param [in,out] arrInOut arr origin array from which the extraction will be done.
 * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
 * \param [in] srcArr input array that will be used as source of copy for ids in [\b idsOfSelectBg, \b idsOfSelectEnd)
 * \param [in] srcArrIndex index array of \b srcArr
 * 
 * \sa MEDCouplingUMesh::SetPartOfIndexedArrays2 MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx
 */
void MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx2(int start, int end, int step, DataArrayInt *arrInOut, const DataArrayInt *arrIndxIn,
                                                      const DataArrayInt *srcArr, const DataArrayInt *srcArrIndex) throw(INTERP_KERNEL::Exception)
{
  if(arrInOut==0 || arrIndxIn==0 || srcArr==0 || srcArrIndex==0)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx2 : presence of null pointer in input parameter !");
  int nbOfTuples=arrIndxIn->getNumberOfTuples()-1;
  const int *arrIndxInPtr=arrIndxIn->getConstPointer();
  const int *srcArrIndexPtr=srcArrIndex->getConstPointer();
  int *arrInOutPtr=arrInOut->getPointer();
  const int *srcArrPtr=srcArr->getConstPointer();
  int nbOfElemsToSet=DataArray::GetNumberOfItemGivenBESRelative(start,end,step,"MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx2 : ");
  int it=start;
  for(int i=0;i<nbOfElemsToSet;i++,srcArrIndexPtr++,it+=step)
    {
      if(it>=0 && it<nbOfTuples)
        {
          if(srcArrIndexPtr[1]-srcArrIndexPtr[0]==arrIndxInPtr[it+1]-arrIndxInPtr[it])
            std::copy(srcArrPtr+srcArrIndexPtr[0],srcArrPtr+srcArrIndexPtr[1],arrInOutPtr+arrIndxInPtr[it]);
          else
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx2 : On pos #" << i << " id (idsOfSelectBg[" << i << "]) is " << it << " arrIndxIn[id+1]-arrIndxIn[id]!=srcArrIndex[pos+1]-srcArrIndex[pos] !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx2 : On pos #" << i << " value is " << it << " not in [0," << nbOfTuples << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
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
MEDCouplingUMesh *MEDCouplingUMesh::buildSpreadZonesWithPoly() const throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
  int mdim=getMeshDimension();
  int spaceDim=getSpaceDimension();
  if(mdim!=spaceDim)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSpreadZonesWithPoly : meshdimension and spacedimension do not match !");
  int nbCells=getNumberOfCells();
  std::vector<DataArrayInt *> partition=partitionBySpreadZone();
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > partitionAuto; partitionAuto.reserve(partition.size());
  std::copy(partition.begin(),partition.end(),std::back_insert_iterator<std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > >(partitionAuto));
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> ret=MEDCouplingUMesh::New(getName(),mdim);
  ret->setCoords(getCoords());
  ret->allocateCells((int)partition.size());
  //
  for(std::vector<DataArrayInt *>::const_iterator it=partition.begin();it!=partition.end();it++)
    {
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> tmp=static_cast<MEDCouplingUMesh *>(buildPartOfMySelf((*it)->begin(),(*it)->end(),true));
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cell;
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
      
      ret->insertNextCell((INTERP_KERNEL::NormalizedCellType)cell->getIJSafe(0,0),cell->getNumberOfTuples()-1,cell->getConstPointer()+1);
    }
  //
  ret->finishInsertingCells();
  ret->incrRef(); return ret;
}

/*!
 * This method partitions \b this into contiguous zone.
 * This method only needs a well defined connectivity. Coordinates are not considered here.
 * This method returns a vector of \b newly allocated arrays that the caller has to deal with.
 */
std::vector<DataArrayInt *> MEDCouplingUMesh::partitionBySpreadZone() const throw(INTERP_KERNEL::Exception)
{
  int nbOfCellsCur=getNumberOfCells();
  DataArrayInt *neigh=0,*neighI=0;
  computeNeighborsOfCells(neigh,neighI);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> neighAuto(neigh),neighIAuto(neighI);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ids=DataArrayInt::New(); ids->alloc(nbOfCellsCur,1); ids->iota();
  std::vector<DataArrayInt *> ret;
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > ret2;
  while(nbOfCellsCur>0)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp=MEDCouplingUMesh::ComputeSpreadZoneGradually(neighAuto,neighIAuto);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp3=tmp->buildComplement(nbOfCellsCur);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp2=ids->selectByTupleId(tmp->begin(),tmp->end());
      ret2.push_back(tmp2);  ret.push_back(tmp2);
      nbOfCellsCur=tmp3->getNumberOfTuples();
      if(nbOfCellsCur>0)
        {
          ids=ids->selectByTupleId(tmp3->begin(),tmp3->end());
          MEDCouplingUMesh::ExtractFromIndexedArrays(tmp3->begin(),tmp3->end(),neighAuto,neighIAuto,neigh,neighI);
          neighAuto=neigh;
          neighIAuto=neighI;
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> renum=tmp3->invertArrayN2O2O2N(nbOfCellsCur+tmp->getNumberOfTuples());
          neighAuto->transformWithIndArr(renum->begin(),renum->end());
        }
    }
  for(std::vector<DataArrayInt *>::const_iterator it=ret.begin();it!=ret.end();it++)
    (*it)->incrRef();
  return ret;
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

MEDCouplingUMeshCellIterator::MEDCouplingUMeshCellIterator(MEDCouplingUMesh *mesh, MEDCouplingUMeshCell *itc, int bg, int end):_mesh(mesh),_cell(itc),
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

MEDCouplingUMeshCellEntry::MEDCouplingUMeshCellEntry(MEDCouplingUMesh *mesh,  INTERP_KERNEL::NormalizedCellType type, MEDCouplingUMeshCell *itc, int bg, int end):_mesh(mesh),_type(type),
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

int MEDCouplingUMeshCellEntry::getNumberOfElems() const
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
  const int *c=_mesh->getNodalConnectivity()->getConstPointer();
  const int *ci=_mesh->getNodalConnectivityIndex()->getConstPointer();
  if(_cell_id<_nb_cell)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)c[ci[_cell_id]];
      int nbOfElems=(int)std::distance(ci+_cell_id,std::find_if(ci+_cell_id,ci+_nb_cell,ParaMEDMEMImpl::ConnReader(c,type)));
      int startId=_cell_id;
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
      std::copy(_conn+1,_conn+_conn_lgth,std::ostream_iterator<int>(oss," "));
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

const int *MEDCouplingUMeshCell::getAllConn(int& lgth) const
{
  lgth=_conn_lgth;
  if(_conn_lgth!=NOTICABLE_FIRST_VAL)
    return _conn;
  else
    return 0;
}
