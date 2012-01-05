// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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
bool MEDCouplingUMesh::isEqual(const MEDCouplingMesh *other, double prec) const
{
  const MEDCouplingUMesh *otherC=dynamic_cast<const MEDCouplingUMesh *>(other);
  if(!otherC)
    return false;
  if(!MEDCouplingPointSet::isEqual(other,prec))
    return false;
  if(_mesh_dim!=otherC->_mesh_dim)
    return false;
  if(_types!=otherC->_types)
    return false;
  if(_nodal_connec!=0 || otherC->_nodal_connec!=0)
    if(_nodal_connec==0 || otherC->_nodal_connec==0)
      return false;
  if(_nodal_connec!=otherC->_nodal_connec)
    if(!_nodal_connec->isEqual(*otherC->_nodal_connec))
      return false;
  if(_nodal_connec_index!=0 || otherC->_nodal_connec_index!=0)
    if(_nodal_connec_index==0 || otherC->_nodal_connec_index==0)
      return false;
  if(_nodal_connec_index!=otherC->_nodal_connec_index)
    if(!_nodal_connec_index->isEqual(*otherC->_nodal_connec_index))
      return false;
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
  nodeCor=da->substr(oldNbOfNodes);
  da=m->mergeNodes(prec,areNodesMerged,newNbOfNodes);
  if(nodeCor->isIdentity())
    {
      nodeCor->decrRef();
      nodeCor=0;
    }
  //
  da=m->zipConnectivityTraducer(cellCompPol);
  int maxId=*std::max_element(da->getConstPointer(),da->getConstPointer()+getNumberOfCells());
  pt=std::find_if(da->getConstPointer()+getNumberOfCells(),da->getConstPointer()+da->getNbOfElems(),std::bind2nd(std::greater<int>(),maxId));
  if(pt!=da->getConstPointer()+da->getNbOfElems())
    {
      nodeCor->decrRef(); nodeCor=0;
      throw INTERP_KERNEL::Exception("checkDeepEquivalWith : some cells in other are not in this !");
    }
  cellCor=DataArrayInt::New();
  cellCor->alloc(otherC->getNumberOfCells(),1);
  std::copy(da->getConstPointer()+getNumberOfCells(),da->getConstPointer()+da->getNbOfElems(),cellCor->getPointer());
  if(cellCor->isIdentity())
    {
      cellCor->decrRef();
      cellCor=0;
    }
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
  cellCor=DataArrayInt::New();
  cellCor->alloc(otherC->getNumberOfCells(),1);
  std::copy(da->getConstPointer()+getNumberOfCells(),da->getConstPointer()+da->getNbOfElems(),cellCor->getPointer());
  if(cellCor->isIdentity())
    {
      cellCor->decrRef();
      cellCor=0;
    }
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
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildDescendingConnectivity2(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const throw(INTERP_KERNEL::Exception)
{
  return buildDescendingConnectivityGen(desc,descIndx,revDesc,revDescIndx,MEDCouplingOrientationSensitiveNbrer);
}

/// @cond INTERNAL

/*!
 * \b WARNING this method do the assumption that connectivity lies on the coordinates set.
 * For speed reasons no check of this will be done.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildDescendingConnectivityGen(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx, DimM1DescNbrer nbrer) const throw(INTERP_KERNEL::Exception)
{
  checkFullyDefined();
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
              int cellDM1Id=meshDM1Type.size();
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
  int nbOfCellsInConstituent=meshDM1Type.size();
  ret->allocateCells(nbOfCellsInConstituent);
  revDescIndx->alloc(nbOfCellsInConstituent+1,1);
  int *tmp3=revDescIndx->getPointer(); tmp3[0]=0;
  for(int ii=0;ii<nbOfCellsInConstituent;ii++)
    {
      ret->insertNextCell((INTERP_KERNEL::NormalizedCellType)meshDM1Type[ii],meshDM1ConnIndex[ii+1]-meshDM1ConnIndex[ii],&meshDM1Conn[meshDM1ConnIndex[ii]]);
      tmp3[ii+1]=tmp3[ii]+revDescMeshConnB[ii].size();
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
    tmp3[jj+1]=tmp3[jj]+descMeshConnB[jj].size();
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
 * This method convert this into dynamic types without changing geometry.
 * That is to say if 'this' is a 2D, mesh after the invocation of this method it will contain only polygons.
 * If 'this' is a 3D mesh after the invocation of this method it will contain only polyhedra.
 * If mesh dimension is not in [2,3] an exception is thrown.
 * Of course pay attention that the resulting mesh is slower than previous one.
 * This method is above all designed to test more extensively algorithms able to deal with polygons/polyhedra.
 */
void MEDCouplingUMesh::convertToPolyTypes(const std::vector<int>& cellIdsToConvert)
{
  checkFullyDefined();
  int dim=getMeshDimension();
  if(dim<2 || dim>3)
    throw INTERP_KERNEL::Exception("Invalid mesh dimension : must be 2 or 3 !");
  if(dim==2)
    {
      const int *connIndex=_nodal_connec_index->getConstPointer();
      int *conn=_nodal_connec->getPointer();
      for(std::vector<int>::const_iterator iter=cellIdsToConvert.begin();iter!=cellIdsToConvert.end();iter++)
        {
          const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[connIndex[*iter]]);
          if(!cm.isDynamic())
            conn[connIndex[*iter]]=INTERP_KERNEL::NORM_POLYGON;
          else
            conn[connIndex[*iter]]=INTERP_KERNEL::NORM_QPOLYG;
        }
    }
  else
    {
      int *connIndex=_nodal_connec_index->getPointer();
      int connIndexLgth=_nodal_connec_index->getNbOfElems();
      const int *connOld=_nodal_connec->getConstPointer();
      int connOldLgth=_nodal_connec->getNbOfElems();
      std::vector<int> connNew(connOld,connOld+connOldLgth);
      for(std::vector<int>::const_iterator iter=cellIdsToConvert.begin();iter!=cellIdsToConvert.end();iter++)
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
          unsigned newLgth=work-tmp-1;
          unsigned delta=newLgth-lgthOld;
          std::transform(connIndex+(*iter)+1,connIndex+connIndexLgth,connIndex+(*iter)+1,std::bind2nd(std::plus<int>(),delta));
          connNew.insert(connNew.begin()+posP1,tmp+lgthOld,tmp+newLgth);
          std::copy(tmp,tmp+lgthOld,connNew.begin()+pos+1);
          delete [] tmp;
        }
      _nodal_connec->alloc(connNew.size(),1);
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
  convertToPolyTypes(cellIds);
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
          int n2=std::distance(c+ci[i]+1,c+ci[i+1]);
          if(n2%2!=0)
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::convertExtrudedPolyhedra : cell # " << i << " is a polhedron with 1 face but there is a mismatch of number of nodes in face should be even !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          int n1=n2/2;
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
          int n1=std::distance(c+ci[i]+1,c+ci[i+1])/2;
          newc=std::copy(c+ci[i],c+ci[i]+n1+1,newc);
          *newc++=-1;
          for(int j=0;j<n1;j++)
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
 */
void MEDCouplingUMesh::unPolyze()
{
  checkFullyDefined();
  if(getMeshDimension()<=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::unPolyze works on umeshes with meshdim equals to 2 or 3 !");
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
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      INTERP_KERNEL::NormalizedCellType newType=INTERP_KERNEL::NORM_ERROR;
      int newLgth;
      if(cm.isDynamic())
        {
          if(cm.getDimension()==2)
            {
              INTERP_KERNEL::AutoPtr<int> tmp=new int[lgthOfCurCell-1];
              std::copy(conn+posOfCurCell+1,conn+posOfCurCell+lgthOfCurCell,(int *)tmp);
              newType=INTERP_KERNEL::CellSimplify::tryToUnPoly2D(cm.isQuadratic(),tmp,lgthOfCurCell-1,conn+newPos+1,newLgth);
            }
          if(cm.getDimension()==3)
            {
              int nbOfFaces,lgthOfPolyhConn;
              INTERP_KERNEL::AutoPtr<int> zipFullReprOfPolyh=INTERP_KERNEL::CellSimplify::getFullPolyh3DCell(type,conn+posOfCurCell+1,lgthOfCurCell-1,nbOfFaces,lgthOfPolyhConn);
              newType=INTERP_KERNEL::CellSimplify::tryToUnPoly3D(zipFullReprOfPolyh,nbOfFaces,lgthOfPolyhConn,conn+newPos+1,newLgth);
            }
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
  computeTypes();
}

/*!
 * Array returned is the correspondance old to new.
 * The maximum value stored in returned array is the number of nodes of 'this' minus 1 after call of this method.
 * The size of returned array is the number of nodes of the old (previous to the call of this method) number of nodes.
 * -1 values in returned array means that the corresponding old node is no more used.
 */
DataArrayInt *MEDCouplingUMesh::zipCoordsTraducer() throw(INTERP_KERNEL::Exception)
{
  int nbOfNodes=getNumberOfNodes();
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbOfNodes,1);
  int *traducer=ret->getPointer();
  std::fill(traducer,traducer+nbOfNodes,-1);
  int nbOfCells=getNumberOfCells();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  int *conn=_nodal_connec->getPointer();
  for(int i=0;i<nbOfCells;i++)
    for(int j=connIndex[i]+1;j<connIndex[i+1];j++)
      if(conn[j]>=0)
        traducer[conn[j]]=1;
  int newNbOfNodes=std::count(traducer,traducer+nbOfNodes,1);
  std::transform(traducer,traducer+nbOfNodes,traducer,MEDCouplingAccVisit());
  for(int i=0;i<nbOfCells;i++)
    for(int j=connIndex[i]+1;j<connIndex[i+1];j++)
      if(conn[j]>=0)
        conn[j]=traducer[conn[j]];
  DataArrayDouble *newCoords=_coords->renumberAndReduce(traducer,newNbOfNodes);
  setCoords(newCoords);
  newCoords->decrRef();
  return ret;
}

/*!
 * This method stands if 'cell1' and 'cell2' are equals regarding 'compType' policy.
 * The semantic of 'compType' is specified in MEDCouplingUMesh::zipConnectivityTraducer method.
 */
bool MEDCouplingUMesh::areCellsEqual(int cell1, int cell2, int compType) const
{
  switch(compType)
    {
    case 0:
      return areCellsEqual0(cell1,cell2);
    case 1:
      return areCellsEqual1(cell1,cell2);
    case 2:
      return areCellsEqual2(cell1,cell2);
    }
  throw INTERP_KERNEL::Exception("Unknown comparison asked ! Must be in 0,1 or 2.");
}

/*!
 * This method is the last step of the MEDCouplingUMesh::zipConnectivityTraducer with policy 0.
 */
bool MEDCouplingUMesh::areCellsEqual0(int cell1, int cell2) const
{
  const int *conn=getNodalConnectivity()->getConstPointer();
  const int *connI=getNodalConnectivityIndex()->getConstPointer();
  return std::equal(conn+connI[cell1],conn+connI[cell1+1],conn+connI[cell2]);
}

/*!
 * This method is the last step of the MEDCouplingUMesh::zipConnectivityTraducer with policy 1.
 */
bool MEDCouplingUMesh::areCellsEqual1(int cell1, int cell2) const
{
  throw INTERP_KERNEL::Exception("Policy comparison, not implemented yet !");
}

/*!
 * This method is the last step of the MEDCouplingUMesh::zipConnectivityTraducer with policy 2.
 */
bool MEDCouplingUMesh::areCellsEqual2(int cell1, int cell2) const
{
  const int *conn=getNodalConnectivity()->getConstPointer();
  const int *connI=getNodalConnectivityIndex()->getConstPointer();
  std::set<int> s1(conn+connI[cell1],conn+connI[cell1+1]);
  std::set<int> s2(conn+connI[cell2],conn+connI[cell2+1]);
  return s1==s2;
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
  int sz=c1.size();
  if(sz!=(int)c2.size())
    return false;
  for(int i=0;i<sz;i++)
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
  std::set<int>::const_iterator end=cand.end(); end--;
  bool ret=false;
  for(std::set<int>::const_iterator iter=cand.begin();iter!=end && !ret;iter++)
    {
      std::set<int>::const_iterator begin2=iter; begin2++;
      for(std::set<int>::const_iterator iter2=begin2;iter2!=cand.end();iter2++)
        {
          if(areCellsEqual(*iter,*iter2,compType))
            {
              if(!ret)
                {
                  result.push_back(*iter);
                  ret=true;
                }
              result.push_back(*iter2);
            }
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
              resI.push_back(res.size());
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
 *                 0 : exactly. A cell is detected to be the same if and only if the connectivity is exactly the same without permutation and types same too. This is the strongest policy.
 *                 1 : permutation. cell1 and cell2 are equal if and the connectivity of cell2 can be deduced by those of cell1 by direct permutation and their type equal.
 *                 2 : nodal. cell1 and cell2 are equal if and only if cell1 and cell2 have same type and have the same nodes constituting connectivity. This is the laziest policy.
 * @return the correspondance array old to new.
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
  const int nbOfTupleSmCells=commonCellsI.size()-1;
  int id=-1;
  std::vector<int> cellsToKeep;
  for(int i=0;i<nbOfTupleSmCells;i++)
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
 * build a sub part of 'this'. This sub part is defined by the cell ids contained in the array in [begin,end).
 * @param begin begin of array containing the cell ids to keep.
 * @param end end of array of cell ids to keep. \b WARNING end param is \b not included ! Idem STL standard definitions.
 * @param keepCoords that specifies if you want or not to keep coords as this or zip it (see ParaMEDMEM::MEDCouplingUMesh::zipCoords). If true zipCoords is \b NOT called, if false, zipCoords is called.
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
      return (MEDCouplingUMesh *)this;
    }
}

DataArrayInt *MEDCouplingUMesh::getCellIdsFullyIncludedInNodeIds(const int *partBg, const int *partEnd) const
{
  std::vector<int> cellIdsKept;
  fillCellIdsToKeepFromNodeIds(partBg,partEnd,true,cellIdsKept);
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(cellIdsKept.size(),1);
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
      int refLgth=std::min(connOfCell.size(),fastFinder.size());
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
  ret->alloc(cellIdsKept.size(),1);
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
  DataArrayInt *desc,*descIndx,*revDesc,*revDescIndx;
  desc=DataArrayInt::New(); descIndx=DataArrayInt::New(); revDesc=DataArrayInt::New(); revDescIndx=DataArrayInt::New();
  MEDCouplingUMesh *subMesh=buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  desc->decrRef(); descIndx->decrRef(); revDesc->decrRef(); revDescIndx->decrRef();
  MEDCouplingUMesh *ret=(MEDCouplingUMesh *)subMesh->buildPartOfMySelfNode(begin,end,fullyIn);
  subMesh->decrRef();
  return ret;
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
DataArrayInt *MEDCouplingUMesh::findCellsIdsOnBoundary() const throw(INTERP_KERNEL::Exception)
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
  ret2->alloc(ret.size(),1);
  std::copy(ret.begin(),ret.end(),ret2->getPointer());
  ret2->setName("BoundaryCells");
  return ret2;
}

/*!
 * This methods returns set of nodes lying on the boundary of this.
 */
void MEDCouplingUMesh::findBoundaryNodes(std::vector<int>& nodes) const
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
  std::set<int> ret;
  int nbOfCells=meshDM1->getNumberOfCells();
  const int *revDescIndxC=revDescIndx->getConstPointer();
  std::vector<int> boundaryCells;
  for(int i=0;i<nbOfCells;i++)
    if(revDescIndxC[i+1]-revDescIndxC[i]==1)
      boundaryCells.push_back(i);
  revDescIndx->decrRef();
  const int *conn=meshDM1->getNodalConnectivity()->getConstPointer();
  const int *connIndx=meshDM1->getNodalConnectivityIndex()->getConstPointer();
  for(std::vector<int>::const_iterator iter=boundaryCells.begin();iter!=boundaryCells.end();iter++)
    for(int k=connIndx[*iter]+1;k<connIndx[(*iter)+1];k++)
      ret.insert(conn[k]);
  nodes.resize(ret.size());
  std::copy(ret.begin(),ret.end(),nodes.begin());
  //
  meshDM1->decrRef();
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
 */
void MEDCouplingUMesh::renumberNodes2(const int *newNodeNumbers, int newNbOfNodes)
{
  MEDCouplingPointSet::renumberNodes2(newNodeNumbers,newNbOfNodes);
  renumberNodesInConn(newNodeNumbers);
}

/*!
 * This method renumbers nodes in connectivity only without any reference with coords.
 * Use it with care !
 * @param 'newNodeNumbers' in old2New convention
 */
void MEDCouplingUMesh::renumberNodesInConn(const int *newNodeNumbers)
{
  int *conn=getNodalConnectivity()->getPointer();
  const int *connIndex=getNodalConnectivityIndex()->getConstPointer();
  int nbOfCells=getNumberOfCells();
  for(int i=0;i<nbOfCells;i++)
    for(int iconn=connIndex[i]+1;iconn!=connIndex[i+1];iconn++)
      {
        int& node=conn[iconn];
        if(node>=0)//avoid polyhedron separator
          {
            node=newNodeNumbers[node];
          }
      }
  _nodal_connec->declareAsNew();
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
  DataArrayInt *newConn=DataArrayInt::New();
  newConn->alloc(_nodal_connec->getNumberOfTuples(),_nodal_connec->getNumberOfComponents());
  newConn->copyStringInfoFrom(*_nodal_connec);
  DataArrayInt *newConnI=DataArrayInt::New();
  newConnI->alloc(_nodal_connec_index->getNumberOfTuples(),_nodal_connec_index->getNumberOfComponents());
  newConnI->copyStringInfoFrom(*_nodal_connec_index);
  //
  int *newC=newConn->getPointer();
  int *newCI=newConnI->getPointer();
  int loc=0;
  newCI[0]=loc;
  for(int i=0;i<nbCells;i++)
    {
      int pos=std::distance(array,std::find(array,array+nbCells,i));
      int nbOfElts=connI[pos+1]-connI[pos];
      newC=std::copy(conn+connI[pos],conn+connI[pos+1],newC);
      loc+=nbOfElts;
      newCI[i+1]=loc;
    }
  //
  setConnectivity(newConn,newConnI);
  //
  newConn->decrRef();
  newConnI->decrRef();
  if(check)
    delete [] (int *)array;
}

/*!
 * Given a boundary box 'bbox' returns elements 'elems' contained in this 'bbox'.
 * Warning 'elems' is incremented during the call so if elems is not empty before call returned elements will be
 * added in 'elems' parameter.
 */
void MEDCouplingUMesh::getCellsInBoundingBox(const double *bbox, double eps, std::vector<int>& elems)
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

/*!
 * Returns coordinates of node with id 'nodeId' and append it in 'coo'.
 */
void MEDCouplingUMesh::getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const
{
  const double *cooPtr=_coords->getConstPointer();
  int spaceDim=getSpaceDimension();
  coo.insert(coo.end(),cooPtr+spaceDim*nodeId,cooPtr+spaceDim*(nodeId+1));
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
    return std::count_if(pt+ptI[cellId]+1,pt+ptI[cellId+1],std::bind2nd(std::not_equal_to<int>(),-1));
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
MEDCouplingUMesh::MEDCouplingUMesh(const MEDCouplingUMesh& other, bool deepCpy):MEDCouplingPointSet(other,deepCpy),_iterator(-1),_mesh_dim(other._mesh_dim),
                                                                                _nodal_connec(0),_nodal_connec_index(0),
                                                                                _types(other._types)
{
  if(other._nodal_connec)
    _nodal_connec=other._nodal_connec->performCpy(deepCpy);
  if(other._nodal_connec_index)
    _nodal_connec_index=other._nodal_connec_index->performCpy(deepCpy);
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
 * This is the low algorithm of buildPartOfMySelf. 
 * Keeps from 'this' only cells which constituing point id are in the ids specified by ['begin','end').
 * The return newly allocated mesh will share the same coordinates as 'this'.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildPartOfMySelfKeepCoords(const int *begin, const int *end) const
{
  checkFullyDefined();
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New();
  ret->_mesh_dim=_mesh_dim;
  ret->setCoords(_coords);
  int nbOfElemsRet=end-begin;
  int *connIndexRet=new int[nbOfElemsRet+1];
  connIndexRet[0]=0;
  const int *conn=_nodal_connec->getConstPointer();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  int newNbring=0;
  for(const int *work=begin;work!=end;work++,newNbring++)
    connIndexRet[newNbring+1]=connIndexRet[newNbring]+connIndex[*work+1]-connIndex[*work];
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
  connIndexRetArr->useArray(connIndexRet,true,CPP_DEALLOC,nbOfElemsRet+1,1);
  ret->setConnectivity(connRetArr,connIndexRetArr,false);
  ret->_types=types;
  connRetArr->decrRef();
  connIndexRetArr->decrRef();
  ret->copyTinyInfoFrom(this);
  std::string name(getName());
  int sz=strlen(PART_OF_NAME);
  if((int)name.length()>=sz)
    name=name.substr(0,sz);
  if(name!=PART_OF_NAME)
    {
      std::ostringstream stream; stream << PART_OF_NAME << getName();
      ret->setName(stream.str().c_str());
    }
  else
    ret->setName(getName());
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
  int nbelem=std::distance(begin,end);
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
  int nbelems=std::distance(begin,end);
  int nbComp=getMeshDimension()+1;
  array->alloc(nbelems,nbComp);
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
          for(int i=0;i<nbelems;i++)
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

  INTERP_KERNEL::Edge *MEDCouplingUMeshBuildQPFromEdge(INTERP_KERNEL::NormalizedCellType typ, std::map<int, INTERP_KERNEL::Node *>& mapp2, const int *bg)
  {
    INTERP_KERNEL::Edge *ret=0;
    switch(typ)
      {
      case INTERP_KERNEL::NORM_SEG2:
        {
          ret=new INTERP_KERNEL::EdgeLin(mapp2[bg[0]],mapp2[bg[1]]);
          break;
        }
      case INTERP_KERNEL::NORM_SEG3:
        {
          ret=new INTERP_KERNEL::EdgeArcCircle(mapp2[bg[0]],mapp2[bg[2]],mapp2[bg[1]]);
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
    std::map<int, INTERP_KERNEL::Node *> mapp2;
    const double *coo=mDesc->getCoords()->getConstPointer();
    const int *c=mDesc->getNodalConnectivity()->getConstPointer();
    const int *cI=mDesc->getNodalConnectivityIndex()->getConstPointer();
    std::set<int> s;
    for(std::vector<int>::const_iterator it=candidates.begin();it!=candidates.end();it++)
      s.insert(c+cI[*it]+1,c+cI[(*it)+1]);
    for(std::set<int>::const_iterator it2=s.begin();it2!=s.end();it2++)
      {
        INTERP_KERNEL::Node *n=new INTERP_KERNEL::Node(coo[2*(*it2)],coo[2*(*it2)+1]);
        mapp[n]=(*it2); mapp2[*it2]=n;
      }
    INTERP_KERNEL::QuadraticPolygon *ret=new INTERP_KERNEL::QuadraticPolygon;
    for(std::vector<int>::const_iterator it=candidates.begin();it!=candidates.end();it++)
      {
        INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType)c[cI[*it]];
        ret->pushBack(MEDCouplingUMeshBuildQPFromEdge(typ,mapp2,c+cI[*it]+1));
      }
    for(std::set<int>::const_iterator it2=s.begin();it2!=s.end();it2++)
      mapp2[*it2]->decrRef();
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

  void MEDCouplingUMeshBuildQPFromMesh2(const double *coo1, int offset1, const double *coo2, int offset2, const std::vector<double>& addCoo,
                                        bool isQuad1, const int *desc1Bg, const int *desc1End, const std::vector<std::vector<int> >& intesctEdges1,
                                        bool isQuad2, const int *desc2Bg, const int *desc2End, const std::vector<std::vector<int> >& intesctEdges2,
                                        /*output*/INTERP_KERNEL::QuadraticPolygon& pol1, INTERP_KERNEL::QuadraticPolygon& pol2, std::map<INTERP_KERNEL::Node *,int>& mapp)
  {
    std::map<int,INTERP_KERNEL::Node *> mappRev;
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
    for(const int *desc2=desc2Bg;desc2!=desc2End;desc2++)
      {
        int eltId2=abs(*desc2)-1;
        for(std::vector<int>::const_iterator it2=intesctEdges2[eltId2].begin();it2!=intesctEdges2[eltId2].end();it2++)
          {
            int curNodeId2=*it2;
            std::map<int,INTERP_KERNEL::Node *>::const_iterator it=mappRev.find(curNodeId2);
            if(it==mappRev.end())
              {
                INTERP_KERNEL::Node *node=MEDCouplingUMeshBuildQPNode(curNodeId2,coo1,offset1,coo2,offset2,addCoo);
                mapp[node]=curNodeId2;
                mappRev[curNodeId2]=node;
              }
          }
      }
    //
    pol1.buildFromCrudeDataArray(mappRev,isQuad1,desc1Bg,desc1End,intesctEdges1);
    pol2.buildFromCrudeDataArray(mappRev,isQuad2,desc2Bg,desc2End,intesctEdges2);
    for(std::map<int,INTERP_KERNEL::Node *>::const_iterator it=mappRev.begin();it!=mappRev.end();it++)
      (*it).second->decrRef();
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
}

/*!
 * This method is only available for a mesh with meshDim==2 and spaceDim==2||spaceDim==3.
 * This method returns a vector 'cells' where all detected butterfly cells have been added to cells.
 * A 2D cell is considered to be butterfly if it exists at least one pair of distinct edges of it that intersect each other
 * anywhere excepted their extremities. An INTERP_KERNEL::NORM_NORI3 could \b not be butterfly.
 */
void MEDCouplingUMesh::checkButterflyCells(std::vector<int>& cells) const
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
      if(isButterfly2DCell(cell2DinS2,isQuad))
        cells.push_back(i);
      cell2DinS2.clear();
    }
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
      INTERP_KERNEL::EdgeArcCircle::getArcOfCirclePassingThru(p0,p1,p2,tmp3,radius,alpha,alpha0);
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
          INTERP_KERNEL::EdgeArcCircle::getArcOfCirclePassingThru(p0r,p1r,p2r,tmp3,radius,alpha,alpha0);
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
      *newConnIPtr++=newc.size();
    }
  newConn->alloc(newc.size()*nbOf1DCells,1);
  int *newConnPtr=newConn->getPointer();
  int deltaPerLev=isQuad?2*nbOfNodesOf1Lev:nbOfNodesOf1Lev;
  newConnIPtr=newConnI->getPointer();
  for(int iz=0;iz<nbOf1DCells;iz++)
    {
      if(iz!=0)
        std::transform(newConnIPtr+1,newConnIPtr+1+nbOf2DCells,newConnIPtr+1+iz*nbOf2DCells,std::bind2nd(std::plus<int>(),newConnIPtr[iz*nbOf2DCells]));
      for(std::vector<int>::const_iterator iter=newc.begin();iter!=newc.end();iter++,newConnPtr++)
        {
          int icell=iter-newc.begin();
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
 * The returned array is of size 3*n where n is the number of different types present in 'this'. For every k in [0,n] ret[3*k+2]==0 because it has no
 * sense here. This parameter is kept only for compatibility with other methode listed above.
 */
std::vector<int> MEDCouplingUMesh::getDistributionOfTypes() const throw(INTERP_KERNEL::Exception)
{
  checkConnectivityFullyDefined();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const int *work=connI;
  int nbOfCells=getNumberOfCells();
  std::size_t n=getAllTypes().size();
  std::vector<int> ret(3*n);
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
      ret[3*i+1]=std::distance(work,work2);
      work=work2;
    }
  return ret;
}

/*!
 * This method is used to check that this has contiguous cell type in same order than described in 'code'.
 * Format of 'code' is the following. 'code' should be of size 3*n and non empty. If not an exception is thrown.
 * foreach k in [0,n) on 3*k pos represent the geometric type and 3*k+1 number of elements of type 3*k.
 * 3*k+2 refers if different from -1 the pos in 'idsPerType' to get the corresponding array.
 * If 2 or more same geometric type is in 'code' and exception is thrown too.
 *
 * This method fistly checks
 * If it exists k so that 3*k geometric type is not in geometric types of this an exception will be thrown.
 * If it exists k so that 3*k geometric type exists but the number of consecutive cell types does not match,
 * an exception is thrown too.
 * 
 * If all geometric types in 'code' are exactly those in 'this' null pointer is returned.
 * If it exists a geometric type in 'this' \b not in 'code' \b no exception is thrown and a DataArrayInt instance is returned that the user has the responsability
 * to deallocate.
 */
DataArrayInt *MEDCouplingUMesh::checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception)
{
  if(code.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::checkTypeConsistencyAndContig : code is empty, should not !");
  int sz=code.size();
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
      int offset=std::distance(connI,i);
      if(code[3*kk+2]==-1)
        {
          const int *j=std::find_if(i+1,connI+nbOfCells,ParaMEDMEMImpl::ConnReader(conn,(int)(*it)));
          int pos2=std::distance(i,j);
          for(int k=0;k<pos2;k++)
            *retPtr++=k+offset;
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
      typeRangeVals.push_back(std::distance(connI,i));
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
  int sz2=idsInPflPerType2.size();
  idsInPflPerType.resize(sz2);
  for(int i=0;i<sz2;i++)
    {
      DataArrayInt *locDa=idsInPflPerType2[i];
      locDa->incrRef();
      idsInPflPerType[i]=locDa;
    }
  int sz=idsPerType2.size();
  idsPerType.resize(sz);
  for(int i=0;i<sz;i++)
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
      int tmp;
      ret->getMaxValue(tmp);
      ret->decrRef();
      std::ostringstream oss; oss << "MEDCouplingUMesh::emulateMEDMEMBDC : input N-1 mesh present a cell not in descending mesh ... Id of cell is " << tmp << " !";
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
  static const int N=19;
  static const INTERP_KERNEL::NormalizedCellType MEDMEM_ORDER[N] = { INTERP_KERNEL::NORM_POINT1, INTERP_KERNEL::NORM_SEG2, INTERP_KERNEL::NORM_SEG3, INTERP_KERNEL::NORM_TRI3, INTERP_KERNEL::NORM_QUAD4, INTERP_KERNEL::NORM_TRI6, INTERP_KERNEL::NORM_QUAD8, INTERP_KERNEL::NORM_TETRA4, INTERP_KERNEL::NORM_PYRA5, INTERP_KERNEL::NORM_PENTA6, INTERP_KERNEL::NORM_HEXA8, INTERP_KERNEL::NORM_HEXGP12, INTERP_KERNEL::NORM_TETRA10, INTERP_KERNEL::NORM_PYRA13, INTERP_KERNEL::NORM_PENTA15, INTERP_KERNEL::NORM_HEXA20, INTERP_KERNEL::NORM_POLYGON, INTERP_KERNEL::NORM_QPOLYG, INTERP_KERNEL::NORM_POLYHED };
  checkConnectivityFullyDefined();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=getRenumArrForConsecutiveCellTypesSpec(MEDMEM_ORDER,MEDMEM_ORDER+N);
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
      int pos=std::distance(orderBg,std::find(orderBg,orderEnd,curType));
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
  tmpb->alloc(std::distance(orderBg,orderEnd),1);
  tmpb->fillWithZero();
  int *tmp=tmpa->getPointer();
  int *tmp2=tmpb->getPointer();
  for(const int *i=connI;i!=connI+nbOfCells;i++)
    {
      const INTERP_KERNEL::NormalizedCellType *where=std::find(orderBg,orderEnd,(INTERP_KERNEL::NormalizedCellType)conn[*i]);
      if(where!=orderEnd)
        {
          int pos=std::distance(orderBg,where);
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
      int beginCellId=std::distance(connI,i);
      i=std::find_if(i+1,connI+nbOfCells,ParaMEDMEMImpl::ConnReader(conn,(int)curType));
      int endCellId=std::distance(connI,i);
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
  ret->alloc(r.size(),1);
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
  int nbOfTuple=std::distance(begin,end);
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
 * There should be \b no presence of null pointer into 'a'.
 * The returned mesh will contain aggregation of nodes in 'a' (in the same order) and aggregation of
 * cells in meshes in 'a' (in the same order too).
 */
MEDCouplingUMesh *MEDCouplingUMesh::MergeUMeshes(std::vector<const MEDCouplingUMesh *>& a) throw(INTERP_KERNEL::Exception)
{
  std::size_t sz=a.size();
  if(sz==0)
    return MergeUMeshesLL(a);
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> > bb(sz);
  std::vector< const MEDCouplingUMesh * > aa(sz);
  int spaceDim=-3;
  for(std::size_t i=0;i<sz && spaceDim==-3;i++)
    {
      const MEDCouplingUMesh *cur=a[i];
      if(!cur)
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::MergeUMeshes : item #" << i << " in input array is empty !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
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
      int meshLgth=(*iter)->getMeshLength();
      nodalPtr=std::copy(nod,nod+meshLgth,nodalPtr);
      if(iter!=meshes.begin())
        nodalIndexPtr=std::transform(index+1,index+nbOfCells+1,nodalIndexPtr,std::bind2nd(std::plus<int>(),offset));
      else
        nodalIndexPtr=std::copy(index,index+nbOfCells+1,nodalIndexPtr);
      offset+=meshLgth;
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
  int nbOfMeshes=meshes.size();
  int offset=0;
  const int *o2nPtr=o2n->getConstPointer();
  for(int i=0;i<nbOfMeshes;i++)
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
        int nbOfRadFaces=std::distance(connBg+1,connEnd);
        for(int i=0;i<nbOfRadFaces;i++)
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
  int sz=std::distance(begin,end);
  if(isQuadratic)
    sz/=2;
  for(int i=0;i<sz;i++)
    {
      v[0]+=coords[3*begin[i]+1]*coords[3*begin[(i+1)%sz]+2]-coords[3*begin[i]+2]*coords[3*begin[(i+1)%sz]+1];
      v[1]+=coords[3*begin[i]+2]*coords[3*begin[(i+1)%sz]]-coords[3*begin[i]]*coords[3*begin[(i+1)%sz]+2];
      v[2]+=coords[3*begin[i]]*coords[3*begin[(i+1)%sz]+1]-coords[3*begin[i]+1]*coords[3*begin[(i+1)%sz]];
    }
  return vec[0]*v[0]+vec[1]*v[1]+vec[2]*v[2]<0.;
}

/*!
 * The polyhedron is specfied by its connectivity nodes in [begin,end).
 */
bool MEDCouplingUMesh::IsPolyhedronWellOriented(const int *begin, const int *end, const double *coords)
{
  std::vector<std::pair<int,int> > edges;
  int nbOfFaces=std::count(begin,end,-1)+1;
  const int *bgFace=begin;
  for(int i=0;i<nbOfFaces;i++)
    {
      const int *endFace=std::find(bgFace+1,end,-1);
      int nbOfEdgesInFace=std::distance(bgFace,endFace);
      for(int j=0;j<nbOfEdgesInFace;j++)
        {
          std::pair<int,int> p1(bgFace[j],bgFace[(j+1)%nbOfEdgesInFace]);
          if(std::find(edges.begin(),edges.end(),p1)!=edges.end())
            return false;
          edges.push_back(p1);
        }
      bgFace=endFace+1;
    }
  return INTERP_KERNEL::calculateVolumeForPolyh2<int,INTERP_KERNEL::ALL_C_MODE>(begin,std::distance(begin,end),coords)>-EPS_FOR_POLYH_ORIENTATION;
}

/*!
 * This method tries to obtain a well oriented polyhedron.
 * If the algorithm fails, an exception will be thrown.
 */
void MEDCouplingUMesh::TryToCorrectPolyhedronOrientation(int *begin, int *end, const double *coords) throw(INTERP_KERNEL::Exception)
{
  std::vector<std::pair<int,int> > edges;
  int nbOfFaces=std::count(begin,end,-1)+1;
  int *bgFace=begin;
  std::vector<bool> isPerm(nbOfFaces);
  for(int i=0;i<nbOfFaces;i++)
    {
      int *endFace=std::find(bgFace+1,end,-1);
      int nbOfEdgesInFace=std::distance(bgFace,endFace);
      for(int l=0;l<nbOfEdgesInFace;l++)
        {
          std::pair<int,int> p1(bgFace[l],bgFace[(l+1)%nbOfEdgesInFace]);
          edges.push_back(p1);
        }
      int *bgFace2=endFace+1;
      for(int k=i+1;k<nbOfFaces;k++)
        {
          int *endFace2=std::find(bgFace2+1,end,-1);
          int nbOfEdgesInFace2=std::distance(bgFace2,endFace2);
          for(int j=0;j<nbOfEdgesInFace2;j++)
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
  if(INTERP_KERNEL::calculateVolumeForPolyh2<int,INTERP_KERNEL::ALL_C_MODE>(begin,std::distance(begin,end),coords)<-EPS_FOR_POLYH_ORIENTATION)
    {//not lucky ! The first face was not correctly oriented : reorient all faces...
      bgFace=begin;
      for(int i=0;i<nbOfFaces;i++)
        {
          int *endFace=std::find(bgFace+1,end,-1);
          int nbOfEdgesInFace=std::distance(bgFace,endFace);
          std::vector<int> tmp(nbOfEdgesInFace-1);
          std::copy(bgFace+1,endFace,tmp.rbegin());
          std::copy(tmp.begin(),tmp.end(),bgFace+1);
          bgFace=endFace+1;
        }
    }
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
  static const int PARAMEDMEM2VTKTYPETRADUCER[INTERP_KERNEL::NORM_MAXTYPE+1]={1,3,21,5,9,7,22,-1,23,-1,-1,-1,-1,-1,10,14,13,-1,12,-1,24,-1,16,27,-1,26,-1,-1,-1,-1,25,42,-1};
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
  std::vector< std::vector<int> > intersectEdge1, subDiv2;
  MEDCouplingUMesh *m1Desc=0,*m2Desc=0;
  DataArrayInt *desc1=0,*descIndx1=0,*revDesc1=0,*revDescIndx1=0,*desc2=0,*descIndx2=0,*revDesc2=0,*revDescIndx2=0;
  std::vector<double> addCoo;
  INTERP_KERNEL::QUADRATIC_PLANAR::_precision=eps;
  IntersectDescending2DMeshes(m1,m2,eps,intersectEdge1,subDiv2,m1Desc,desc1,descIndx1,revDesc1,revDescIndx1,
                              m2Desc,desc2,descIndx2,revDesc2,revDescIndx2,addCoo);
  revDesc1->decrRef(); revDescIndx1->decrRef(); revDesc2->decrRef(); revDescIndx2->decrRef();
  std::vector< std::vector<int> > intersectEdge2;
  BuildIntersectEdges(m1Desc,m2Desc,addCoo,subDiv2,intersectEdge2);
  std::vector<bool> b1=m1Desc->getQuadraticStatus();
  std::vector<bool> b2=m1Desc->getQuadraticStatus();
  subDiv2.clear(); m1Desc->decrRef(); m2Desc->decrRef();
  std::vector<int> cr,crI;
  std::vector<int> cNb1,cNb2;
  BuildIntersecting2DCellsFromEdges(eps,m1,b1,desc1->getConstPointer(),descIndx1->getConstPointer(),intersectEdge1,m2,b2,desc2->getConstPointer(),descIndx2->getConstPointer(),intersectEdge2,addCoo,
                                    /* outputs -> */cr,crI,cNb1,cNb2);
  //
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> addCooDa=DataArrayDouble::New();
  addCooDa->alloc(addCoo.size()/2,2);
  std::copy(addCoo.begin(),addCoo.end(),addCooDa->getPointer());
  std::vector<const DataArrayDouble *> coordss(3);
  coordss[0]=m1->getCoords(); coordss[1]=m2->getCoords(); coordss[2]=addCooDa;
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> coo=DataArrayDouble::Aggregate(coordss);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> ret=MEDCouplingUMesh::New("Intersect2D",2);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> conn=DataArrayInt::New(); conn->alloc(cr.size(),1); std::copy(cr.begin(),cr.end(),conn->getPointer());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> connI=DataArrayInt::New(); connI->alloc(crI.size(),1); std::copy(crI.begin(),crI.end(),connI->getPointer());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> c1=DataArrayInt::New(); c1->alloc(cNb1.size(),1); std::copy(cNb1.begin(),cNb1.end(),c1->getPointer()); cellNb1=c1;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> c2=DataArrayInt::New(); c2->alloc(cNb2.size(),1); std::copy(cNb2.begin(),cNb2.end(),c2->getPointer()); cellNb2=c2;
  ret->setConnectivity(conn,connI,true);
  ret->setCoords(coo);
  ret->incrRef(); c1->incrRef(); c2->incrRef(); desc1->decrRef(); descIndx1->decrRef(); desc2->decrRef(); descIndx2->decrRef();
  return ret;
}

/// @endcond

void MEDCouplingUMesh::BuildIntersecting2DCellsFromEdges(double eps, const MEDCouplingUMesh *m1, const std::vector<bool>& b1, const int *desc1, const int *descIndx1, const std::vector<std::vector<int> >& intesctEdges1,
                                                         const MEDCouplingUMesh *m2, const std::vector<bool>& b2, const int *desc2, const int *descIndx2, const std::vector<std::vector<int> >& intesctEdges2,
                                                         const std::vector<double>& addCoords,
                                                         std::vector<int>& cr, std::vector<int>& crI, std::vector<int>& cNb1, std::vector<int>& cNb2)
{
  static const int SPACEDIM=2;
  std::vector<double> bbox1,bbox2;
  const double *coo1=m1->getCoords()->getConstPointer();
  int offset1=m1->getNumberOfNodes();
  const double *coo2=m2->getCoords()->getConstPointer();
  int offset2=offset1+m2->getNumberOfNodes();
  m1->getBoundingBoxForBBTree(bbox1);
  m2->getBoundingBoxForBBTree(bbox2);
  BBTree<SPACEDIM,int> myTree(&bbox2[0],0,0,m2->getNumberOfCells(),-eps);
  int ncell1=m1->getNumberOfCells();
  crI.push_back(0);
  for(int i=0;i<ncell1;i++)
    {
      std::vector<int> candidates2;
      myTree.getIntersectingElems(&bbox1[i*2*SPACEDIM],candidates2);
      for(std::vector<int>::const_iterator it2=candidates2.begin();it2!=candidates2.end();it2++)
        {
          INTERP_KERNEL::QuadraticPolygon pol1,pol2;
          std::map<INTERP_KERNEL::Node *,int> mapp;
          MEDCouplingUMeshBuildQPFromMesh2(coo1,offset1,coo2,offset2,addCoords,
                                           b1[i],desc1+descIndx1[i],desc1+descIndx1[i+1],intesctEdges1,
                                           b2[*it2],desc2+descIndx2[*it2],desc2+descIndx2[*it2+1],intesctEdges2,
                                           /* output */pol1,pol2,mapp);
          pol1.buildPartitionsAbs(pol2,mapp,i,*it2,cr,crI,cNb1,cNb2);
        }
    }
}

/*!
 * This method is private and is the first step of Partition of 2D mesh (spaceDim==2 and meshDim==2).
 * 
 */
void MEDCouplingUMesh::IntersectDescending2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps,
                                                   std::vector< std::vector<int> >& intersectEdge1, std::vector< std::vector<int> >& subDiv2,
                                                   MEDCouplingUMesh *& m1Desc, DataArrayInt *&desc1, DataArrayInt *&descIndx1, DataArrayInt *&revDesc1, DataArrayInt *&revDescIndx1,
                                                   MEDCouplingUMesh *& m2Desc, DataArrayInt *&desc2, DataArrayInt *&descIndx2, DataArrayInt *&revDesc2, DataArrayInt *&revDescIndx2,
                                                   std::vector<double>& addCoo) throw(INTERP_KERNEL::Exception)
{
  static const int SPACEDIM=2;
  desc1=DataArrayInt::New(); descIndx1=DataArrayInt::New(); revDesc1=DataArrayInt::New(); revDescIndx1=DataArrayInt::New();
  desc2=DataArrayInt::New(); descIndx2=DataArrayInt::New(); revDesc2=DataArrayInt::New(); revDescIndx2=DataArrayInt::New();
  m1Desc=m1->buildDescendingConnectivity2(desc1,descIndx1,revDesc1,revDescIndx1);
  m2Desc=m2->buildDescendingConnectivity2(desc2,descIndx2,revDesc2,revDescIndx2);
  const int *c1=m1Desc->getNodalConnectivity()->getConstPointer();
  const int *ci1=m1Desc->getNodalConnectivityIndex()->getConstPointer();
  std::vector<double> bbox1,bbox2;
  m1Desc->getBoundingBoxForBBTree(bbox1);
  m2Desc->getBoundingBoxForBBTree(bbox2);
  int ncell1=m1Desc->getNumberOfCells();
  int ncell2=m2Desc->getNumberOfCells();
  intersectEdge1.resize(ncell1);
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
          pol1->splitAbs(*pol2,map1,map2,offset1,offset2,candidates2,intersectEdge1[i],subDiv2,addCoo);
          delete pol2;
          delete pol1;
        }
      else
        intersectEdge1[i].insert(intersectEdge1[i].end(),c1+ci1[i]+1,c1+ci1[i+1]);
    }
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
      std::map<int, INTERP_KERNEL::Node *> mapp2;
      std::map<INTERP_KERNEL::Node *, int> mapp22;
      for(int j=0;j<nnode;j++)
        {
          INTERP_KERNEL::Node *nn=new INTERP_KERNEL::Node(coo[2*c[(*cI)+j+1]],coo[2*c[(*cI)+j+1]+1]);
          int nnid=c[(*cI)+j+1];
          mapp2[nnid]=nn;
          mapp22[nn]=nnid+offset1;
        }
      INTERP_KERNEL::Edge *e=MEDCouplingUMeshBuildQPFromEdge((INTERP_KERNEL::NormalizedCellType)c[*cI],mapp2,c+(*cI)+1);
      for(std::map<int, INTERP_KERNEL::Node *>::const_iterator it=mapp2.begin();it!=mapp2.end();it++)
        (*it).second->decrRef();
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
      int nbOfElems=std::distance(ci+_cell_id,std::find_if(ci+_cell_id,ci+_nb_cell,ParaMEDMEMImpl::ConnReader(c,type)));
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
