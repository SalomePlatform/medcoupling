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

#include "MEDCouplingExtrudedMesh.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "CellModel.hxx"

#include "InterpolationUtils.hxx"

#include <limits>
#include <algorithm>
#include <functional>
#include <iterator>
#include <sstream>
#include <cmath>
#include <list>
#include <set>

using namespace ParaMEDMEM;

/*!
 * Build an extruded mesh instance from 3D and 2D unstructured mesh lying on the \b same \b coords.
 * @param mesh3D 3D unstructured mesh.
 * @param mesh2D 2D unstructured mesh lying on the same coordinates than mesh3D. \b Warning mesh2D is \b not \b const
 * because the mesh is aggregated and potentially modified by rotate or translate method.
 * @param cell2DId Id of cell in mesh2D mesh where the computation of 1D mesh will be done.
 */
MEDCouplingExtrudedMesh *MEDCouplingExtrudedMesh::New(const MEDCouplingUMesh *mesh3D, const MEDCouplingUMesh *mesh2D, int cell2DId)
{
  return new MEDCouplingExtrudedMesh(mesh3D,mesh2D,cell2DId);
}

/*!
 * This constructor is here only for unserialisation process.
 * This constructor is normally completely useless for end user.
 */
MEDCouplingExtrudedMesh *MEDCouplingExtrudedMesh::New()
{
  return new MEDCouplingExtrudedMesh;
}

MEDCouplingMeshType MEDCouplingExtrudedMesh::getType() const
{
  return EXTRUDED;
}

std::size_t MEDCouplingExtrudedMesh::getHeapMemorySizeWithoutChildren() const
{
  return MEDCouplingMesh::getHeapMemorySizeWithoutChildren();
}

std::vector<const BigMemoryObject *> MEDCouplingExtrudedMesh::getDirectChildren() const
{
  std::vector<const BigMemoryObject *> ret;
  if(_mesh2D)
    ret.push_back(_mesh2D);
  if(_mesh1D)
    ret.push_back(_mesh1D);
  if(_mesh3D_ids)
    ret.push_back(_mesh3D_ids);
  return ret;
}

/*!
 * This method copyies all tiny strings from other (name and components name).
 * @throw if other and this have not same mesh type.
 */
void MEDCouplingExtrudedMesh::copyTinyStringsFrom(const MEDCouplingMesh *other)
{
  const MEDCouplingExtrudedMesh *otherC=dynamic_cast<const MEDCouplingExtrudedMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCouplingExtrudedMesh::copyTinyStringsFrom : meshes have not same type !");
  MEDCouplingMesh::copyTinyStringsFrom(other);
  _mesh2D->copyTinyStringsFrom(otherC->_mesh2D);
  _mesh1D->copyTinyStringsFrom(otherC->_mesh1D);
}

MEDCouplingExtrudedMesh::MEDCouplingExtrudedMesh(const MEDCouplingUMesh *mesh3D, const MEDCouplingUMesh *mesh2D, int cell2DId)
try:_mesh2D(const_cast<MEDCouplingUMesh *>(mesh2D)),_mesh1D(MEDCouplingUMesh::New()),_mesh3D_ids(0),_cell_2D_id(cell2DId)
{
  if(_mesh2D!=0)
    _mesh2D->incrRef();
  computeExtrusion(mesh3D);
  setName(mesh3D->getName().c_str());
}
catch(INTERP_KERNEL::Exception& e)
  {
    if(_mesh2D)
      _mesh2D->decrRef();
    if(_mesh1D)
      _mesh1D->decrRef();
    if(_mesh3D_ids)
      _mesh3D_ids->decrRef();
    throw e;
  }

MEDCouplingExtrudedMesh::MEDCouplingExtrudedMesh():_mesh2D(0),_mesh1D(0),_mesh3D_ids(0),_cell_2D_id(-1)
{
}

MEDCouplingExtrudedMesh::MEDCouplingExtrudedMesh(const MEDCouplingExtrudedMesh& other, bool deepCopy):MEDCouplingMesh(other),_cell_2D_id(other._cell_2D_id)
{
  if(deepCopy)
    {
      _mesh2D=other._mesh2D->clone(true);
      _mesh1D=other._mesh1D->clone(true);
      _mesh3D_ids=other._mesh3D_ids->deepCpy();
    }
  else
    {
      _mesh2D=other._mesh2D;
      if(_mesh2D)
        _mesh2D->incrRef();
      _mesh1D=other._mesh1D;
      if(_mesh1D)
        _mesh1D->incrRef();
      _mesh3D_ids=other._mesh3D_ids;
      if(_mesh3D_ids)
        _mesh3D_ids->incrRef();
    }
}

int MEDCouplingExtrudedMesh::getNumberOfCells() const
{
  return _mesh2D->getNumberOfCells()*_mesh1D->getNumberOfCells();
}

int MEDCouplingExtrudedMesh::getNumberOfNodes() const
{
  return _mesh2D->getNumberOfNodes();
}

int MEDCouplingExtrudedMesh::getSpaceDimension() const
{
  return 3;
}

int MEDCouplingExtrudedMesh::getMeshDimension() const
{
  return 3;
}

MEDCouplingMesh *MEDCouplingExtrudedMesh::deepCpy() const
{
  return clone(true);
}

MEDCouplingExtrudedMesh *MEDCouplingExtrudedMesh::clone(bool recDeepCpy) const
{
  return new MEDCouplingExtrudedMesh(*this,recDeepCpy);
}

bool MEDCouplingExtrudedMesh::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingExtrudedMesh::isEqualIfNotWhy : input other pointer is null !");
  const MEDCouplingExtrudedMesh *otherC=dynamic_cast<const MEDCouplingExtrudedMesh *>(other);
  std::ostringstream oss;
  if(!otherC)
    {
      reason="mesh given in input is not castable in MEDCouplingExtrudedMesh !";
      return false;
    }
  if(!MEDCouplingMesh::isEqualIfNotWhy(other,prec,reason))
    return false;
  if(!_mesh2D->isEqualIfNotWhy(otherC->_mesh2D,prec,reason))
    {
      reason.insert(0,"Mesh2D unstructured meshes differ : ");
      return false;
    }
  if(!_mesh1D->isEqualIfNotWhy(otherC->_mesh1D,prec,reason))
    {
      reason.insert(0,"Mesh1D unstructured meshes differ : ");
      return false;
    }
  if(!_mesh3D_ids->isEqualIfNotWhy(*otherC->_mesh3D_ids,reason))
    {
      reason.insert(0,"Mesh3D ids DataArrayInt instances differ : ");
      return false;
    }
  if(_cell_2D_id!=otherC->_cell_2D_id)
    {
      oss << "Cell 2D id of the two extruded mesh differ : this = " << _cell_2D_id << " other = " <<  otherC->_cell_2D_id;
      reason=oss.str();
      return false;
    }
  return true;
}

bool MEDCouplingExtrudedMesh::isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const
{
  const MEDCouplingExtrudedMesh *otherC=dynamic_cast<const MEDCouplingExtrudedMesh *>(other);
  if(!otherC)
    return false;
  if(!_mesh2D->isEqualWithoutConsideringStr(otherC->_mesh2D,prec))
    return false;
  if(!_mesh1D->isEqualWithoutConsideringStr(otherC->_mesh1D,prec))
    return false;
  if(!_mesh3D_ids->isEqualWithoutConsideringStr(*otherC->_mesh3D_ids))
    return false;
  if(_cell_2D_id!=otherC->_cell_2D_id)
    return false;
  return true;
}

void MEDCouplingExtrudedMesh::checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                   DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("MEDCouplingExtrudedMesh::checkDeepEquivalWith : not implemented yet !");
}

void MEDCouplingExtrudedMesh::checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                              DataArrayInt *&cellCor) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("MEDCouplingExtrudedMesh::checkDeepEquivalOnSameNodesWith : not implemented yet !");
}

INTERP_KERNEL::NormalizedCellType MEDCouplingExtrudedMesh::getTypeOfCell(int cellId) const
{
  const int *ids=_mesh3D_ids->getConstPointer();
  int nbOf3DCells=_mesh3D_ids->getNumberOfTuples();
  const int *where=std::find(ids,ids+nbOf3DCells,cellId);
  if(where==ids+nbOf3DCells)
    throw INTERP_KERNEL::Exception("Invalid cellId specified >= getNumberOfCells() !");
  int nbOfCells2D=_mesh2D->getNumberOfCells();
  int locId=((int)std::distance(ids,where))%nbOfCells2D;
  INTERP_KERNEL::NormalizedCellType tmp=_mesh2D->getTypeOfCell(locId);
  return INTERP_KERNEL::CellModel::GetCellModel(tmp).getExtrudedType();
}

std::set<INTERP_KERNEL::NormalizedCellType> MEDCouplingExtrudedMesh::getAllGeoTypes() const
{
  std::set<INTERP_KERNEL::NormalizedCellType> ret2D(_mesh2D->getAllGeoTypes());
  std::set<INTERP_KERNEL::NormalizedCellType> ret;
  for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator it=ret2D.begin();it!=ret2D.end();it++)
    ret.insert(INTERP_KERNEL::CellModel::GetCellModel(*it).getExtrudedType());
  return ret;
}

DataArrayInt *MEDCouplingExtrudedMesh::giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  INTERP_KERNEL::NormalizedCellType revExtTyp=cm.getReverseExtrudedType();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  if(revExtTyp==INTERP_KERNEL::NORM_ERROR)
    {
      ret->alloc(0,1);
      return ret.retn();
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp=_mesh2D->giveCellsWithType(revExtTyp);
  int nbOfLevs=_mesh1D->getNumberOfCells();
  int nbOfCells2D=_mesh2D->getNumberOfCells();
  int nbOfTuples=tmp->getNumberOfTuples();
  ret->alloc(nbOfLevs*nbOfTuples,1);
  int *pt=ret->getPointer();
  for(int i=0;i<nbOfLevs;i++,pt+=nbOfTuples)
    std::transform(tmp->begin(),tmp->end(),pt,std::bind2nd(std::plus<int>(),i*nbOfCells2D));
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret2=ret->renumberR(_mesh3D_ids->begin());
  ret2->sort();
  return ret2.retn();
}

DataArrayInt *MEDCouplingExtrudedMesh::computeNbOfNodesPerCell() const
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret2D=_mesh2D->computeNbOfNodesPerCell();
  int nbOfLevs=_mesh1D->getNumberOfCells();
  int nbOfCells2D=_mesh2D->getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret3D=DataArrayInt::New(); ret3D->alloc(nbOfLevs*nbOfCells2D,1);
  int *pt=ret3D->getPointer();
  for(int i=0;i<nbOfLevs;i++,pt+=nbOfCells2D)
     std::copy(ret2D->begin(),ret2D->end(),pt);
  ret3D->applyLin(2,0,0);
  return ret3D->renumberR(_mesh3D_ids->begin());
}

DataArrayInt *MEDCouplingExtrudedMesh::computeNbOfFacesPerCell() const
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret2D=_mesh2D->computeNbOfNodesPerCell();
  int nbOfLevs=_mesh1D->getNumberOfCells();
  int nbOfCells2D=_mesh2D->getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret3D=DataArrayInt::New(); ret3D->alloc(nbOfLevs*nbOfCells2D,1);
  int *pt=ret3D->getPointer();
  for(int i=0;i<nbOfLevs;i++,pt+=nbOfCells2D)
     std::copy(ret2D->begin(),ret2D->end(),pt);
  ret3D->applyLin(2,2,0);
  return ret3D->renumberR(_mesh3D_ids->begin());
}

DataArrayInt *MEDCouplingExtrudedMesh::computeEffectiveNbOfNodesPerCell() const
{
  return computeNbOfNodesPerCell();
}

int MEDCouplingExtrudedMesh::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const
{
  int ret=0;
  int nbOfCells2D=_mesh2D->getNumberOfCells();
  for(int i=0;i<nbOfCells2D;i++)
    {
      INTERP_KERNEL::NormalizedCellType t=_mesh2D->getTypeOfCell(i);
      if(INTERP_KERNEL::CellModel::GetCellModel(t).getExtrudedType()==type)
        ret++;
    }
  return ret*_mesh1D->getNumberOfCells();
}

void MEDCouplingExtrudedMesh::getNodeIdsOfCell(int cellId, std::vector<int>& conn) const
{
  int nbOfCells2D=_mesh2D->getNumberOfCells();
  int nbOfNodes2D=_mesh2D->getNumberOfNodes();
  int locId=cellId%nbOfCells2D;
  int lev=cellId/nbOfCells2D;
  std::vector<int> tmp,tmp2;
  _mesh2D->getNodeIdsOfCell(locId,tmp);
  tmp2=tmp;
  std::transform(tmp.begin(),tmp.end(),tmp.begin(),std::bind2nd(std::plus<int>(),nbOfNodes2D*lev));
  std::transform(tmp2.begin(),tmp2.end(),tmp2.begin(),std::bind2nd(std::plus<int>(),nbOfNodes2D*(lev+1)));
  conn.insert(conn.end(),tmp.begin(),tmp.end());
  conn.insert(conn.end(),tmp2.begin(),tmp2.end());
}

void MEDCouplingExtrudedMesh::getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const
{
  int nbOfNodes2D=_mesh2D->getNumberOfNodes();
  int locId=nodeId%nbOfNodes2D;
  int lev=nodeId/nbOfNodes2D;
  std::vector<double> tmp,tmp2;
  _mesh2D->getCoordinatesOfNode(locId,tmp);
  tmp2=tmp;
  int spaceDim=_mesh1D->getSpaceDimension();
  const double *z=_mesh1D->getCoords()->getConstPointer();
  std::transform(tmp.begin(),tmp.end(),z+lev*spaceDim,tmp.begin(),std::plus<double>());
  std::transform(tmp2.begin(),tmp2.end(),z+(lev+1)*spaceDim,tmp2.begin(),std::plus<double>());
  coo.insert(coo.end(),tmp.begin(),tmp.end());
  coo.insert(coo.end(),tmp2.begin(),tmp2.end());
}

std::string MEDCouplingExtrudedMesh::simpleRepr() const
{
  std::ostringstream ret;
  ret << "3D Extruded mesh from a 2D Surf Mesh with name : \"" << getName() << "\"\n";
  ret << "Description of mesh : \"" << getDescription() << "\"\n";
  int tmpp1,tmpp2;
  double tt=getTime(tmpp1,tmpp2);
  ret << "Time attached to the mesh [unit] : " << tt << " [" << getTimeUnit() << "]\n";
  ret << "Iteration : " << tmpp1  << " Order : " << tmpp2 << "\n";
  ret << "Cell id where 1D mesh has been deduced : " << _cell_2D_id << "\n";
  ret << "Number of cells : " << getNumberOfCells() << "(" << _mesh2D->getNumberOfCells() << "x" << _mesh1D->getNumberOfCells() << ")\n";
  ret << "1D Mesh info : _____________________\n\n\n";
  ret << _mesh1D->simpleRepr();
  ret << "\n\n\n2D Mesh info : _____________________\n\n\n" << _mesh2D->simpleRepr() << "\n\n\n";
  return ret.str();
}

std::string MEDCouplingExtrudedMesh::advancedRepr() const
{
  std::ostringstream ret;
  ret << "3D Extruded mesh from a 2D Surf Mesh with name : \"" << getName() << "\"\n";
  ret << "Description of mesh : \"" << getDescription() << "\"\n";
  int tmpp1,tmpp2;
  double tt=getTime(tmpp1,tmpp2);
  ret << "Time attached to the mesh (unit) : " << tt << " (" << getTimeUnit() << ")\n";
  ret << "Iteration : " << tmpp1  << " Order : " << tmpp2 << "\n";
  ret << "Cell id where 1D mesh has been deduced : " << _cell_2D_id << "\n";
  ret << "Number of cells : " << getNumberOfCells() << "(" << _mesh2D->getNumberOfCells() << "x" << _mesh1D->getNumberOfCells() << ")\n";
  ret << "1D Mesh info : _____________________\n\n\n";
  ret << _mesh1D->advancedRepr();
  ret << "\n\n\n2D Mesh info : _____________________\n\n\n" << _mesh2D->advancedRepr() << "\n\n\n";
  ret << "3D cell ids per level :\n";
  return ret.str();
}

void MEDCouplingExtrudedMesh::checkCoherency() const throw (INTERP_KERNEL::Exception)
{
}

void MEDCouplingExtrudedMesh::checkCoherency1(double eps) const
{
  checkCoherency();
}

void MEDCouplingExtrudedMesh::checkCoherency2(double eps) const
{
  checkCoherency1(eps);
}

void MEDCouplingExtrudedMesh::getBoundingBox(double *bbox) const
{
  double bbox2D[6];
  _mesh2D->getBoundingBox(bbox2D);
  const double *nodes1D=_mesh1D->getCoords()->getConstPointer();
  int nbOfNodes1D=_mesh1D->getNumberOfNodes();
  double bbox1DMin[3],bbox1DMax[3],tmp[3];
  std::fill(bbox1DMin,bbox1DMin+3,std::numeric_limits<double>::max());
  std::fill(bbox1DMax,bbox1DMax+3,-(std::numeric_limits<double>::max()));
  for(int i=0;i<nbOfNodes1D;i++)
    {
      std::transform(nodes1D+3*i,nodes1D+3*(i+1),bbox1DMin,bbox1DMin,static_cast<const double& (*)(const double&, const double&)>(std::min<double>));
      std::transform(nodes1D+3*i,nodes1D+3*(i+1),bbox1DMax,bbox1DMax,static_cast<const double& (*)(const double&, const double&)>(std::max<double>));
    }
  std::transform(bbox1DMax,bbox1DMax+3,bbox1DMin,tmp,std::minus<double>());
  int id=(int)std::distance(tmp,std::max_element(tmp,tmp+3));
  bbox[0]=bbox1DMin[0]; bbox[1]=bbox1DMax[0];
  bbox[2]=bbox1DMin[1]; bbox[3]=bbox1DMax[1];
  bbox[4]=bbox1DMin[2]; bbox[5]=bbox1DMax[2];
  bbox[2*id+1]+=tmp[id];
}

void MEDCouplingExtrudedMesh::updateTime() const
{
  if(_mesh2D)
    {
      updateTimeWith(*_mesh2D);
    }
  if(_mesh1D)
    {
      updateTimeWith(*_mesh1D);
    }
}

void MEDCouplingExtrudedMesh::renumberCells(const int *old2NewBg, bool check)
{
  throw INTERP_KERNEL::Exception("Functionnality of renumbering cells unavailable for ExtrudedMesh");
}

MEDCouplingUMesh *MEDCouplingExtrudedMesh::build3DUnstructuredMesh() const
{
  MEDCouplingUMesh *ret=_mesh2D->buildExtrudedMesh(_mesh1D,0);
  const int *renum=_mesh3D_ids->getConstPointer();
  ret->renumberCells(renum,false);
  ret->setName(getName().c_str());
  return ret;
}

MEDCouplingUMesh *MEDCouplingExtrudedMesh::buildUnstructured() const
{
  return build3DUnstructuredMesh();
}

MEDCouplingFieldDouble *MEDCouplingExtrudedMesh::getMeasureField(bool) const
{
  std::string name="MeasureOfMesh_";
  name+=getName();
  MEDCouplingFieldDouble *ret2D=_mesh2D->getMeasureField(true);
  MEDCouplingFieldDouble *ret1D=_mesh1D->getMeasureField(true);
  const double *ret2DPtr=ret2D->getArray()->getConstPointer();
  const double *ret1DPtr=ret1D->getArray()->getConstPointer();
  int nbOf2DCells=_mesh2D->getNumberOfCells();
  int nbOf1DCells=_mesh1D->getNumberOfCells();
  int nbOf3DCells=nbOf2DCells*nbOf1DCells;
  const int *renum=_mesh3D_ids->getConstPointer();
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  ret->setMesh(this);
  ret->synchronizeTimeWithMesh();
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(nbOf3DCells,1);
  double *retPtr=da->getPointer();
  for(int i=0;i<nbOf1DCells;i++)
    for(int j=0;j<nbOf2DCells;j++)
      retPtr[renum[i*nbOf2DCells+j]]=ret2DPtr[j]*ret1DPtr[i];
  ret->setArray(da);
  da->decrRef();
  ret->setName(name.c_str());
  ret2D->decrRef();
  ret1D->decrRef();
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingExtrudedMesh::getMeasureFieldOnNode(bool isAbs) const
{
  //not implemented yet
  return 0;
}

MEDCouplingFieldDouble *MEDCouplingExtrudedMesh::buildOrthogonalField() const
{
  throw INTERP_KERNEL::Exception("MEDCouplingExtrudedMesh::buildOrthogonalField : This method has no sense for MEDCouplingExtrudedMesh that is 3D !");
}

int MEDCouplingExtrudedMesh::getCellContainingPoint(const double *pos, double eps) const
{
  throw INTERP_KERNEL::Exception("MEDCouplingExtrudedMesh::getCellContainingPoint : not implemented yet !");
}

MEDCouplingExtrudedMesh::~MEDCouplingExtrudedMesh()
{
  if(_mesh2D)
    _mesh2D->decrRef();
  if(_mesh1D)
    _mesh1D->decrRef();
  if(_mesh3D_ids)
    _mesh3D_ids->decrRef();
}

void MEDCouplingExtrudedMesh::computeExtrusion(const MEDCouplingUMesh *mesh3D)
{
  const char errMsg1[]="2D mesh is empty unable to compute extrusion !";
  const char errMsg2[]="Coords between 2D and 3D meshes are not the same ! Try MEDCouplingPointSet::tryToShareSameCoords method";
  const char errMsg3[]="No chance to find extrusion pattern in mesh3D,mesh2D couple because nbCells3D%nbCells2D!=0 !";
  if(_mesh2D==0 || mesh3D==0)
    throw INTERP_KERNEL::Exception(errMsg1);
  if(_mesh2D->getCoords()!=mesh3D->getCoords())
    throw INTERP_KERNEL::Exception(errMsg2);
  if(mesh3D->getNumberOfCells()%_mesh2D->getNumberOfCells()!=0)
    throw INTERP_KERNEL::Exception(errMsg3);
  if(!_mesh3D_ids)
    _mesh3D_ids=DataArrayInt::New();
  if(!_mesh1D)
    _mesh1D=MEDCouplingUMesh::New();
  computeExtrusionAlg(mesh3D);
}

void MEDCouplingExtrudedMesh::build1DExtrusion(int idIn3DDesc, int newId, int nbOf1DLev, MEDCouplingUMesh *subMesh,
                                               const int *desc3D, const int *descIndx3D,
                                               const int *revDesc3D, const int *revDescIndx3D,
                                               bool computeMesh1D) throw(INTERP_KERNEL::Exception)
{
  int nbOf2DCells=_mesh2D->getNumberOfCells();
  int start=revDescIndx3D[idIn3DDesc];
  int end=revDescIndx3D[idIn3DDesc+1];
  if(end-start!=1)
    {
      std::ostringstream ost; ost << "Invalid bases 2D mesh specified : 2D cell # " <<  idIn3DDesc;
      ost << " shared by more than 1 3D cell !!!";
      throw INTERP_KERNEL::Exception(ost.str().c_str());
    }
  int current3DCell=revDesc3D[start];
  int current2DCell=idIn3DDesc;
  int *mesh3DIDs=_mesh3D_ids->getPointer();
  mesh3DIDs[newId]=current3DCell;
  const int *conn2D=subMesh->getNodalConnectivity()->getConstPointer();
  const int *conn2DIndx=subMesh->getNodalConnectivityIndex()->getConstPointer();
  for(int i=1;i<nbOf1DLev;i++)
    {
      std::vector<int> conn(conn2D+conn2DIndx[current2DCell]+1,conn2D+conn2DIndx[current2DCell+1]);
      std::sort(conn.begin(),conn.end());
      if(computeMesh1D)
        computeBaryCenterOfFace(conn,i-1);
      current2DCell=findOppositeFaceOf(current2DCell,current3DCell,conn,
                                       desc3D,descIndx3D,conn2D,conn2DIndx);
      start=revDescIndx3D[current2DCell];
      end=revDescIndx3D[current2DCell+1];
      if(end-start!=2)
        {
          std::ostringstream ost; ost << "Expecting to have 2 3D cells attached to 2D cell " << current2DCell << "!";
          ost << " : Impossible or call tryToShareSameCoords method !";
          throw INTERP_KERNEL::Exception(ost.str().c_str());
        }
      if(revDesc3D[start]!=current3DCell)
        current3DCell=revDesc3D[start];
      else
        current3DCell=revDesc3D[start+1];
      mesh3DIDs[i*nbOf2DCells+newId]=current3DCell;
    }
  if(computeMesh1D)
    {
      std::vector<int> conn(conn2D+conn2DIndx[current2DCell]+1,conn2D+conn2DIndx[current2DCell+1]);
      std::sort(conn.begin(),conn.end());
      computeBaryCenterOfFace(conn,nbOf1DLev-1);
      current2DCell=findOppositeFaceOf(current2DCell,current3DCell,conn,
                                       desc3D,descIndx3D,conn2D,conn2DIndx);
      conn.clear();
      conn.insert(conn.end(),conn2D+conn2DIndx[current2DCell]+1,conn2D+conn2DIndx[current2DCell+1]);
      std::sort(conn.begin(),conn.end());
      computeBaryCenterOfFace(conn,nbOf1DLev);
    }
}

int MEDCouplingExtrudedMesh::findOppositeFaceOf(int current2DCell, int current3DCell, const std::vector<int>& connSorted,
                                                const int *desc3D, const int *descIndx3D,
                                                const int *conn2D, const int *conn2DIndx) throw(INTERP_KERNEL::Exception)
{
  int start=descIndx3D[current3DCell];
  int end=descIndx3D[current3DCell+1];
  bool found=false;
  for(const int *candidate2D=desc3D+start;candidate2D!=desc3D+end && !found;candidate2D++)
    {
      if(*candidate2D!=current2DCell)
        {
          std::vector<int> conn2(conn2D+conn2DIndx[*candidate2D]+1,conn2D+conn2DIndx[*candidate2D+1]);
          std::sort(conn2.begin(),conn2.end());
          std::list<int> intersect;
          std::set_intersection(connSorted.begin(),connSorted.end(),conn2.begin(),conn2.end(),
                                std::insert_iterator< std::list<int> >(intersect,intersect.begin()));
          if(intersect.empty())
            return *candidate2D;
        }
    }
  std::ostringstream ost; ost << "Impossible to find an opposite 2D face of face # " <<  current2DCell;
  ost << " in 3D cell # " << current3DCell << " : Impossible or call tryToShareSameCoords method !";
  throw INTERP_KERNEL::Exception(ost.str().c_str());
}

void MEDCouplingExtrudedMesh::computeBaryCenterOfFace(const std::vector<int>& nodalConnec, int lev1DId)
{
  double *zoneToUpdate=_mesh1D->getCoords()->getPointer()+lev1DId*3;
  std::fill(zoneToUpdate,zoneToUpdate+3,0.);
  const double *coords=_mesh2D->getCoords()->getConstPointer();
  for(std::vector<int>::const_iterator iter=nodalConnec.begin();iter!=nodalConnec.end();iter++)
    std::transform(zoneToUpdate,zoneToUpdate+3,coords+3*(*iter),zoneToUpdate,std::plus<double>());
  std::transform(zoneToUpdate,zoneToUpdate+3,zoneToUpdate,std::bind2nd(std::multiplies<double>(),(double)(1./(int)nodalConnec.size())));
}

int MEDCouplingExtrudedMesh::FindCorrespCellByNodalConn(const std::vector<int>& nodalConnec, const int *revNodalPtr, const int *revNodalIndxPtr) throw(INTERP_KERNEL::Exception)
{
  std::vector<int>::const_iterator iter=nodalConnec.begin();
  std::set<int> s1(revNodalPtr+revNodalIndxPtr[*iter],revNodalPtr+revNodalIndxPtr[*iter+1]);
  iter++;
  for(;iter!=nodalConnec.end();iter++)
    {
      std::set<int> s2(revNodalPtr+revNodalIndxPtr[*iter],revNodalPtr+revNodalIndxPtr[*iter+1]);
      std::set<int> s3;
      std::set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::insert_iterator< std::set<int> >(s3,s3.end()));
      s1=s3;
    }
  if(s1.size()==1)
    return *(s1.begin());
  std::ostringstream ostr;
  ostr << "Cell with nodal connec : ";
  std::copy(nodalConnec.begin(),nodalConnec.end(),std::ostream_iterator<int>(ostr," "));
  ostr << " is not part of mesh";
  throw INTERP_KERNEL::Exception(ostr.str().c_str());
}

/*!
 * This method is callable on 1Dmeshes (meshDim==1 && spaceDim==3) returned by MEDCouplingExtrudedMesh::getMesh1D typically.
 * These 1Dmeshes (meshDim==1 && spaceDim==3) have a special semantic because these meshes do not specify a static location but a translation along a path.
 * This method checks that 'm1' and 'm2' are compatible, if not an exception is thrown. In case these meshes ('m1' and 'm2') are compatible 2 corresponding meshes
 * are created ('m1r' and 'm2r') that can be used for interpolation.
 * @param m1 input mesh with meshDim==1 and spaceDim==3
 * @param m2 input mesh with meshDim==1 and spaceDim==3
 * @param eps tolerance acceptable to determine compatibility
 * @param m1r output mesh with ref count equal to 1 with meshDim==1 and spaceDim==1
 * @param m2r output mesh with ref count equal to 1 with meshDim==1 and spaceDim==1
 * @param v is the output normalized vector of the common direction of 'm1' and 'm2'  
 * @throw in case that m1 and m2 are not compatible each other.
 */
void MEDCouplingExtrudedMesh::Project1DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps,
                                              MEDCouplingUMesh *&m1r, MEDCouplingUMesh *&m2r, double *v) throw(INTERP_KERNEL::Exception)
{
  if(m1->getSpaceDimension()!=3 || m1->getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("Input meshes are expected to have a spaceDim==3 for Projec1D !");
  m1r=m1->clone(true);
  m2r=m2->clone(true);
  m1r->changeSpaceDimension(1);
  m2r->changeSpaceDimension(1);
  std::vector<int> c;
  std::vector<double> ref,ref2;
  m1->getNodeIdsOfCell(0,c);
  m1->getCoordinatesOfNode(c[0],ref);
  m1->getCoordinatesOfNode(c[1],ref2);
  std::transform(ref2.begin(),ref2.end(),ref.begin(),v,std::minus<double>());
  double n=INTERP_KERNEL::norm<3>(v);
  std::transform(v,v+3,v,std::bind2nd(std::multiplies<double>(),1/n));
  m1->project1D(&ref[0],v,eps,m1r->getCoords()->getPointer());
  m2->project1D(&ref[0],v,eps,m2r->getCoords()->getPointer());
  
}

void MEDCouplingExtrudedMesh::rotate(const double *center, const double *vector, double angle)
{
  _mesh2D->rotate(center,vector,angle);
  _mesh1D->rotate(center,vector,angle);
}

void MEDCouplingExtrudedMesh::translate(const double *vector)
{
  _mesh2D->translate(vector);
  _mesh1D->translate(vector);
}

void MEDCouplingExtrudedMesh::scale(const double *point, double factor)
{
  _mesh2D->scale(point,factor);
  _mesh1D->scale(point,factor);
}

std::vector<int> MEDCouplingExtrudedMesh::getDistributionOfTypes() const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

DataArrayInt *MEDCouplingExtrudedMesh::checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

void MEDCouplingExtrudedMesh::splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

MEDCouplingMesh *MEDCouplingExtrudedMesh::buildPart(const int *start, const int *end) const
{
  // not implemented yet !
  return 0;
}

MEDCouplingMesh *MEDCouplingExtrudedMesh::buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const
{
  // not implemented yet !
  return 0;
}

DataArrayInt *MEDCouplingExtrudedMesh::simplexize(int policy)
{
  throw INTERP_KERNEL::Exception("MEDCouplingExtrudedMesh::simplexize : unavailable for such a type of mesh : Extruded !");
}

MEDCouplingMesh *MEDCouplingExtrudedMesh::mergeMyselfWith(const MEDCouplingMesh *other) const
{
  // not implemented yet !
  return 0;
}

DataArrayDouble *MEDCouplingExtrudedMesh::getCoordinatesAndOwner() const
{
  DataArrayDouble *arr2D=_mesh2D->getCoords();
  DataArrayDouble *arr1D=_mesh1D->getCoords();
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(getNumberOfNodes(),3);
  int nbOf1DLev=_mesh1D->getNumberOfNodes();
  int nbOf2DNodes=_mesh2D->getNumberOfNodes();
  const double *ptSrc=arr2D->getConstPointer();
  double *pt=ret->getPointer();
  std::copy(ptSrc,ptSrc+3*nbOf2DNodes,pt);
  for(int i=1;i<nbOf1DLev;i++)
    {
      std::copy(ptSrc,ptSrc+3*nbOf2DNodes,pt+3*i*nbOf2DNodes);
      double vec[3];
      std::copy(arr1D->getConstPointer()+3*i,arr1D->getConstPointer()+3*(i+1),vec);
      std::transform(arr1D->getConstPointer()+3*(i-1),arr1D->getConstPointer()+3*i,vec,vec,std::minus<double>());
      for(int j=0;j<nbOf2DNodes;j++)
        std::transform(vec,vec+3,pt+3*(i*nbOf2DNodes+j),pt+3*(i*nbOf2DNodes+j),std::plus<double>());
    }
  return ret;
}

DataArrayDouble *MEDCouplingExtrudedMesh::getBarycenterAndOwner() const
{
  throw INTERP_KERNEL::Exception("MEDCouplingExtrudedMesh::getBarycenterAndOwner : not yet implemented !");
}

DataArrayDouble *MEDCouplingExtrudedMesh::computeIsoBarycenterOfNodesPerCell() const
{
  throw INTERP_KERNEL::Exception("MEDCouplingExtrudedMesh::computeIsoBarycenterOfNodesPerCell: not yet implemented !");
}

void MEDCouplingExtrudedMesh::computeExtrusionAlg(const MEDCouplingUMesh *mesh3D)
{
  _mesh3D_ids->alloc(mesh3D->getNumberOfCells(),1);
  int nbOf1DLev=mesh3D->getNumberOfCells()/_mesh2D->getNumberOfCells();
  _mesh1D->setMeshDimension(1);
  _mesh1D->allocateCells(nbOf1DLev);
  int tmpConn[2];
  for(int i=0;i<nbOf1DLev;i++)
    {
      tmpConn[0]=i;
      tmpConn[1]=i+1;
      _mesh1D->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,tmpConn);
    }
  _mesh1D->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(nbOf1DLev+1,3);
  _mesh1D->setCoords(myCoords);
  myCoords->decrRef();
  DataArrayInt *desc,*descIndx,*revDesc,*revDescIndx;
  desc=DataArrayInt::New(); descIndx=DataArrayInt::New(); revDesc=DataArrayInt::New(); revDescIndx=DataArrayInt::New();
  MEDCouplingUMesh *subMesh=mesh3D->buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  DataArrayInt *revNodal2D,*revNodalIndx2D;
  revNodal2D=DataArrayInt::New(); revNodalIndx2D=DataArrayInt::New();
  subMesh->getReverseNodalConnectivity(revNodal2D,revNodalIndx2D);
  const int *nodal2D=_mesh2D->getNodalConnectivity()->getConstPointer();
  const int *nodal2DIndx=_mesh2D->getNodalConnectivityIndex()->getConstPointer();
  const int *revNodal2DPtr=revNodal2D->getConstPointer();
  const int *revNodalIndx2DPtr=revNodalIndx2D->getConstPointer();
  const int *descP=desc->getConstPointer();
  const int *descIndxP=descIndx->getConstPointer();
  const int *revDescP=revDesc->getConstPointer();
  const int *revDescIndxP=revDescIndx->getConstPointer();
  //
  int nbOf2DCells=_mesh2D->getNumberOfCells();
  for(int i=0;i<nbOf2DCells;i++)
    {
      int idInSubMesh;
       std::vector<int> nodalConnec(nodal2D+nodal2DIndx[i]+1,nodal2D+nodal2DIndx[i+1]);
       try
        {
          idInSubMesh=FindCorrespCellByNodalConn(nodalConnec,revNodal2DPtr,revNodalIndx2DPtr);
        }
       catch(INTERP_KERNEL::Exception& e)
         {
           std::ostringstream ostr; ostr << "mesh2D cell # " << i << " is not part of any cell of 3D mesh !\n";
           ostr << e.what();
           throw INTERP_KERNEL::Exception(ostr.str().c_str());
         }
       build1DExtrusion(idInSubMesh,i,nbOf1DLev,subMesh,descP,descIndxP,revDescP,revDescIndxP,i==_cell_2D_id);
    }
  //
  revNodal2D->decrRef();
  revNodalIndx2D->decrRef();
  subMesh->decrRef();
  desc->decrRef();
  descIndx->decrRef();
  revDesc->decrRef();
  revDescIndx->decrRef();
}

void MEDCouplingExtrudedMesh::getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const
{
  std::vector<int> tinyInfo1;
  std::vector<std::string> ls1;
  std::vector<double> ls3;
  _mesh2D->getTinySerializationInformation(ls3,tinyInfo1,ls1);
  std::vector<int> tinyInfo2;
  std::vector<std::string> ls2;
  std::vector<double> ls4;
  _mesh1D->getTinySerializationInformation(ls4,tinyInfo2,ls2);
  tinyInfo.clear(); littleStrings.clear();
  tinyInfo.insert(tinyInfo.end(),tinyInfo1.begin(),tinyInfo1.end());
  littleStrings.insert(littleStrings.end(),ls1.begin(),ls1.end());
  tinyInfo.insert(tinyInfo.end(),tinyInfo2.begin(),tinyInfo2.end());
  littleStrings.insert(littleStrings.end(),ls2.begin(),ls2.end());
  tinyInfo.push_back(_cell_2D_id);
  tinyInfo.push_back((int)tinyInfo1.size());
  tinyInfo.push_back(_mesh3D_ids->getNbOfElems());
  littleStrings.push_back(getName());
  littleStrings.push_back(getDescription());
}

void MEDCouplingExtrudedMesh::resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const
{
  std::size_t sz=tinyInfo.size();
  int sz1=tinyInfo[sz-2];
  std::vector<int> ti1(tinyInfo.begin(),tinyInfo.begin()+sz1);
  std::vector<int> ti2(tinyInfo.begin()+sz1,tinyInfo.end()-3);
  MEDCouplingUMesh *um=MEDCouplingUMesh::New();
  DataArrayInt *a1tmp=DataArrayInt::New();
  DataArrayDouble *a2tmp=DataArrayDouble::New();
  int la1=0,la2=0;
  std::vector<std::string> ls1,ls2;
  um->resizeForUnserialization(ti1,a1tmp,a2tmp,ls1);
  la1+=a1tmp->getNbOfElems(); la2+=a2tmp->getNbOfElems();
  a1tmp->decrRef(); a2tmp->decrRef();
  a1tmp=DataArrayInt::New(); a2tmp=DataArrayDouble::New();
  um->resizeForUnserialization(ti2,a1tmp,a2tmp,ls2);
  la1+=a1tmp->getNbOfElems(); la2+=a2tmp->getNbOfElems();
  a1tmp->decrRef(); a2tmp->decrRef();
  um->decrRef();
  //
  a1->alloc(la1+tinyInfo[sz-1],1);
  a2->alloc(la2,1);
  littleStrings.resize(ls1.size()+ls2.size()+2);
}

void MEDCouplingExtrudedMesh::serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const
{
  a1=DataArrayInt::New(); a2=DataArrayDouble::New();
  DataArrayInt *a1_1=0,*a1_2=0;
  DataArrayDouble *a2_1=0,*a2_2=0;
  _mesh2D->serialize(a1_1,a2_1);
  _mesh1D->serialize(a1_2,a2_2);
  a1->alloc(a1_1->getNbOfElems()+a1_2->getNbOfElems()+_mesh3D_ids->getNbOfElems(),1);
  int *ptri=a1->getPointer();
  ptri=std::copy(a1_1->getConstPointer(),a1_1->getConstPointer()+a1_1->getNbOfElems(),ptri);
  a1_1->decrRef();
  ptri=std::copy(a1_2->getConstPointer(),a1_2->getConstPointer()+a1_2->getNbOfElems(),ptri);
  a1_2->decrRef();
  std::copy(_mesh3D_ids->getConstPointer(),_mesh3D_ids->getConstPointer()+_mesh3D_ids->getNbOfElems(),ptri);
  a2->alloc(a2_1->getNbOfElems()+a2_2->getNbOfElems(),1);
  double *ptrd=a2->getPointer();
  ptrd=std::copy(a2_1->getConstPointer(),a2_1->getConstPointer()+a2_1->getNbOfElems(),ptrd);
  a2_1->decrRef();
  std::copy(a2_2->getConstPointer(),a2_2->getConstPointer()+a2_2->getNbOfElems(),ptrd);
  a2_2->decrRef();
}

void MEDCouplingExtrudedMesh::unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2, const std::vector<std::string>& littleStrings)
{
  setName(littleStrings[littleStrings.size()-2].c_str());
  setDescription(littleStrings.back().c_str());
  std::size_t sz=tinyInfo.size();
  int sz1=tinyInfo[sz-2];
  _cell_2D_id=tinyInfo[sz-3];
  std::vector<int> ti1(tinyInfo.begin(),tinyInfo.begin()+sz1);
  std::vector<int> ti2(tinyInfo.begin()+sz1,tinyInfo.end()-3);
  DataArrayInt *a1tmp=DataArrayInt::New();
  DataArrayDouble *a2tmp=DataArrayDouble::New();
  const int *a1Ptr=a1->getConstPointer();
  const double *a2Ptr=a2->getConstPointer();
  _mesh2D=MEDCouplingUMesh::New();
  std::vector<std::string> ls1,ls2;
  _mesh2D->resizeForUnserialization(ti1,a1tmp,a2tmp,ls1);
  std::copy(a2Ptr,a2Ptr+a2tmp->getNbOfElems(),a2tmp->getPointer());
  std::copy(a1Ptr,a1Ptr+a1tmp->getNbOfElems(),a1tmp->getPointer());
  a2Ptr+=a2tmp->getNbOfElems();
  a1Ptr+=a1tmp->getNbOfElems();
  ls2.insert(ls2.end(),littleStrings.begin(),littleStrings.begin()+ls1.size());
  std::vector<double> d1(1);
  _mesh2D->unserialization(d1,ti1,a1tmp,a2tmp,ls2);
  a1tmp->decrRef(); a2tmp->decrRef();
  //
  ls2.clear();
  ls2.insert(ls2.end(),littleStrings.begin()+ls1.size(),littleStrings.end()-2);
  _mesh1D=MEDCouplingUMesh::New();
  a1tmp=DataArrayInt::New(); a2tmp=DataArrayDouble::New();
  _mesh1D->resizeForUnserialization(ti2,a1tmp,a2tmp,ls1);
  std::copy(a2Ptr,a2Ptr+a2tmp->getNbOfElems(),a2tmp->getPointer());
  std::copy(a1Ptr,a1Ptr+a1tmp->getNbOfElems(),a1tmp->getPointer());
  a1Ptr+=a1tmp->getNbOfElems();
  _mesh1D->unserialization(d1,ti2,a1tmp,a2tmp,ls2);
  a1tmp->decrRef(); a2tmp->decrRef();
  //
  _mesh3D_ids=DataArrayInt::New();
  int szIds=(int)std::distance(a1Ptr,a1->getConstPointer()+a1->getNbOfElems());
  _mesh3D_ids->alloc(szIds,1);
  std::copy(a1Ptr,a1Ptr+szIds,_mesh3D_ids->getPointer());
}

void MEDCouplingExtrudedMesh::writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m=buildUnstructured();
  m->writeVTKLL(ofs,cellData,pointData,byteData);
}

void MEDCouplingExtrudedMesh::reprQuickOverview(std::ostream& stream) const
{
  stream << "MEDCouplingExtrudedMesh C++ instance at " << this << ". Name : \"" << getName() << "\".";
}

std::string MEDCouplingExtrudedMesh::getVTKDataSetType() const
{
  return _mesh2D->getVTKDataSetType();
}
