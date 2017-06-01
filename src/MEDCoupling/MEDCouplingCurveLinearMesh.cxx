// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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

#include "MEDCouplingCurveLinearMesh.hxx"
#include "MEDCouplingPointSet.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include "VolSurfUser.txx"
#include "PointLocatorAlgos.txx"

#include <functional>
#include <algorithm>
#include <sstream>
#include <numeric>

using namespace MEDCoupling;

MEDCouplingCurveLinearMesh::MEDCouplingCurveLinearMesh():_coords(0),_structure(0)
{
}

MEDCouplingCurveLinearMesh::MEDCouplingCurveLinearMesh(const MEDCouplingCurveLinearMesh& other, bool deepCpy):MEDCouplingStructuredMesh(other,deepCpy),_structure(other._structure)
{
  if(deepCpy)
    {
      if((const DataArrayDouble *)other._coords)
        _coords=other._coords->deepCopy();
    }
  else
    _coords=other._coords;
}

MEDCouplingCurveLinearMesh::~MEDCouplingCurveLinearMesh()
{
}

MEDCouplingCurveLinearMesh *MEDCouplingCurveLinearMesh::New()
{
  return new MEDCouplingCurveLinearMesh;
}

MEDCouplingCurveLinearMesh *MEDCouplingCurveLinearMesh::New(const std::string& meshName)
{
  MEDCouplingCurveLinearMesh *ret=new MEDCouplingCurveLinearMesh;
  ret->setName(meshName);
  return ret;
}

MEDCouplingCurveLinearMesh *MEDCouplingCurveLinearMesh::deepCopy() const
{
  return clone(true);
}

MEDCouplingCurveLinearMesh *MEDCouplingCurveLinearMesh::clone(bool recDeepCpy) const
{
  return new MEDCouplingCurveLinearMesh(*this,recDeepCpy);
}

void MEDCouplingCurveLinearMesh::updateTime() const
{
  if((const DataArrayDouble *)_coords)
    updateTimeWith(*_coords);
}

std::size_t MEDCouplingCurveLinearMesh::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(MEDCouplingStructuredMesh::getHeapMemorySizeWithoutChildren());
  ret+=_structure.capacity()*sizeof(int);
  return ret;
}

std::vector<const BigMemoryObject *> MEDCouplingCurveLinearMesh::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  ret.push_back((const DataArrayDouble *)_coords);
  return ret;
}

/*!
 * This method copyies all tiny strings from other (name and components name).
 * @throw if other and this have not same mesh type.
 */
void MEDCouplingCurveLinearMesh::copyTinyStringsFrom(const MEDCouplingMesh *other)
{ 
  const MEDCouplingCurveLinearMesh *otherC=dynamic_cast<const MEDCouplingCurveLinearMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::copyTinyStringsFrom : meshes have not same type !");
  MEDCouplingStructuredMesh::copyTinyStringsFrom(other);
  if((DataArrayDouble *)_coords && (const DataArrayDouble *)otherC->_coords)
    _coords->copyStringInfoFrom(*otherC->_coords);
}

bool MEDCouplingCurveLinearMesh::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::isEqualIfNotWhy : input other pointer is null !");
  const MEDCouplingCurveLinearMesh *otherC=dynamic_cast<const MEDCouplingCurveLinearMesh *>(other);
  if(!otherC)
    {
      reason="mesh given in input is not castable in MEDCouplingCurveLinearMesh !";
      return false;
    }
  if(!MEDCouplingStructuredMesh::isEqualIfNotWhy(other,prec,reason))
    return false;
  std::ostringstream oss; oss.precision(15);
  if(((const DataArrayDouble *)_coords && ((const DataArrayDouble *)otherC->_coords)==0) || (((const DataArrayDouble *)_coords)==0 && (const DataArrayDouble *)otherC->_coords))
    {
      oss << "Only one CurveLinearMesh between the two this and other has its coordinates defined !";
      reason=oss.str();
      return false;
    }
  if((const DataArrayDouble *)_coords)
    {
      if(!_coords->isEqualIfNotWhy(*(otherC->_coords),prec,reason))
        {
          oss << "Coordinates DataArrayDouble of differ :";
          reason.insert(0,oss.str());
          return false;
        }
      if(_structure!=otherC->_structure)
        { reason="CurveLinearMesh structures differ !"; return false; }
    }
  return true;
}

bool MEDCouplingCurveLinearMesh::isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const
{
  const MEDCouplingCurveLinearMesh *otherC=dynamic_cast<const MEDCouplingCurveLinearMesh *>(other);
  if(!otherC)
    return false;
  if(((const DataArrayDouble *)_coords && ((const DataArrayDouble *)otherC->_coords)==0) || (((const DataArrayDouble *)_coords)==0 && (const DataArrayDouble *)otherC->_coords))
    return false;
  if((const DataArrayDouble *)_coords)
    {
      if(!_coords->isEqualWithoutConsideringStr(*(otherC->_coords),prec))
        return false;
      if(_structure!=otherC->_structure)
        return false;
    }
  return true;
}

void MEDCouplingCurveLinearMesh::checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                      DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const
{
  if(!isEqualWithoutConsideringStr(other,prec))
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::checkDeepEquivalWith : Meshes are not the same !");
}

/*!
 * Nothing is done here (except to check that the other is a MEDCoupling::MEDCouplingCurveLinearMesh instance too).
 * The user intend that the nodes are the same, so by construction of MEDCoupling::MEDCouplingCurveLinearMesh, \a this and \a other are the same !
 */
void MEDCouplingCurveLinearMesh::checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                                 DataArrayInt *&cellCor) const
{
  if(!isEqualWithoutConsideringStr(other,prec))
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::checkDeepEquivalOnSameNodesWith : Meshes are not the same !");
}

void MEDCouplingCurveLinearMesh::checkConsistencyLight() const
{
  std::size_t sz=_structure.size(),i=0,nbOfNodes=1;
  if(sz<1)
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::checkConsistencyLight : structure should have a lgth of size 1 at least !");
  for(std::vector<int>::const_iterator it=_structure.begin();it!=_structure.end();it++,i++)
    {
      if((*it)<1)
        { std::ostringstream oss; oss << "MEDCouplingCurveLinearMesh::checkConsistencyLight : At pos #" << i << " of structure value is " << *it << "should be >= 1 !"; throw INTERP_KERNEL::Exception(oss.str().c_str()); }
      nbOfNodes*=*it;
    }
  if(!((const DataArrayDouble *)_coords))
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::checkConsistencyLight : the array is not set !");
  if(!_coords->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::checkConsistencyLight : the array is not allocated !");
  if(_coords->getNumberOfComponents()<1)
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::checkConsistencyLight : the array should have >= 1 components !");
  if(_coords->getNumberOfTuples()!=(int)nbOfNodes)
    {
      std::ostringstream oss; oss << "MEDCouplingCurveLinearMesh::checkConsistencyLight : structure said that number of nodes should be equal to " << nbOfNodes << " but number of tuples in array is equal to " << _coords->getNumberOfTuples() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void MEDCouplingCurveLinearMesh::checkConsistency(double eps) const
{
  checkConsistencyLight();
}

std::size_t MEDCouplingCurveLinearMesh::getNumberOfCells() const
{
  checkConsistencyLight();
  return MEDCouplingStructuredMesh::getNumberOfCells();
}

int MEDCouplingCurveLinearMesh::getNumberOfNodes() const
{
  checkConsistencyLight();
  return MEDCouplingStructuredMesh::getNumberOfNodes();
}

void MEDCouplingCurveLinearMesh::getNodeGridStructure(int *res) const
{
  std::copy(_structure.begin(),_structure.end(),res);
}

/*!
 * MEDCouplingCurveLinearMesh has the property to define 2 space dimensions. One coming from its coordinates. The other coming from the node structure.
 * Normally they should be equal ! This method returns the space dimension from coordinates. If the other one is requested call getSpaceDimensionOnNodeStruct.
 *
 * \sa MEDCouplingStructuredMesh::getSpaceDimensionOnNodeStruct
 */
int MEDCouplingCurveLinearMesh::getSpaceDimension() const
{
  if(!((const DataArrayDouble *)_coords))
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::getSpaceDimension : no array set ! impossible to deduce a space dimension !");
  return _coords->getNumberOfComponents();
}

void MEDCouplingCurveLinearMesh::getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const
{
  if(!((const DataArrayDouble *)_coords))
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::getCoordinatesOfNode : Coordinates not set !");
  int nbOfCompo=_coords->getNumberOfComponents();
  if(nodeId>=0 && nodeId<_coords->getNumberOfTuples())
    coo.insert(coo.end(),_coords->begin()+nodeId*nbOfCompo,_coords->begin()+(nodeId+1)*nbOfCompo);
  else
    { std::ostringstream oss; oss << "MEDCouplingCurveLinearMesh::getCoordinatesOfNode : nodeId has to be in [0," << _coords->getNumberOfTuples() << ") !"; throw INTERP_KERNEL::Exception(oss.str().c_str()); }
}

std::string MEDCouplingCurveLinearMesh::simpleRepr() const
{
  std::ostringstream ret;
  ret << "Curve linear mesh with name : \"" << getName() << "\"\n";
  ret << "Description of mesh : \"" << getDescription() << "\"\n";
  int tmpp1,tmpp2;
  double tt=getTime(tmpp1,tmpp2);
  ret << "Time attached to the mesh [unit] : " << tt << " [" << getTimeUnit() << "]\n";
  ret << "Iteration : " << tmpp1  << " Order : " << tmpp2 << "\n";
  ret << "The nodal stucture of curve linear mesh is : [";
  std::copy(_structure.begin(),_structure.end(),std::ostream_iterator<int>(ret,",")); ret << "]\n";
  ret << "The coords array is this : ";
  if((const DataArrayDouble *)_coords)
    _coords->reprZipWithoutNameStream(ret);
  else
    ret << "no array specified !";
  return ret.str();
}

std::string MEDCouplingCurveLinearMesh::advancedRepr() const
{
  return simpleRepr();
}

const DataArrayDouble *MEDCouplingCurveLinearMesh::getDirectAccessOfCoordsArrIfInStructure() const
{
  return _coords;
}

DataArrayDouble *MEDCouplingCurveLinearMesh::getCoords()
{
  return _coords;
}

const DataArrayDouble *MEDCouplingCurveLinearMesh::getCoords() const
{
  return _coords;
}

void MEDCouplingCurveLinearMesh::setCoords(const DataArrayDouble *coords)
{
  if(coords!=(const DataArrayDouble *)_coords)
    {
      _coords=const_cast<DataArrayDouble *>(coords);
      if(coords)
        coords->incrRef();
      declareAsNew();
    }
}

void MEDCouplingCurveLinearMesh::setNodeGridStructure(const int *gridStructBg, const int *gridStructEnd)
{
  std::size_t sz=std::distance(gridStructBg,gridStructEnd);
  if(sz>=1 && sz<=3)
    {
      _structure.resize(0);
      _structure.insert(_structure.end(),gridStructBg,gridStructEnd);
    }
  else
    {
      std::ostringstream oss; oss << "MEDCouplingCurveLinearMesh::setNodeGridStructure : size of input nodal grid structure (" << sz << ") should be in 1, 2 or 3 !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

std::vector<int> MEDCouplingCurveLinearMesh::getNodeGridStructure() const
{
  return _structure;
}

MEDCouplingStructuredMesh *MEDCouplingCurveLinearMesh::buildStructuredSubPart(const std::vector< std::pair<int,int> >& cellPart) const
{
  checkConsistencyLight();
  int dim(getSpaceDimension());
  std::vector<int> dims(getMeshDimension());
  if(dim!=(int)cellPart.size())
    {
      std::ostringstream oss; oss << "MEDCouplingCurveLinearMesh::buildStructuredSubPart : the space dimension is " << dim << " and cell part size is " << cellPart.size() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::vector< std::pair<int,int> > nodePartFormat(cellPart);
  for(std::vector< std::pair<int,int> >::iterator it=nodePartFormat.begin();it!=nodePartFormat.end();it++)
    (*it).second++;
  MCAuto<DataArrayInt> tmp1(BuildExplicitIdsFrom(getNodeGridStructure(),nodePartFormat));
  MCAuto<MEDCouplingCurveLinearMesh> ret(dynamic_cast<MEDCouplingCurveLinearMesh *>(deepCopy()));
  const DataArrayDouble *coo(ret->getCoords());
  if(coo)
    {
      MCAuto<DataArrayDouble> coo2(coo->selectByTupleIdSafe(tmp1->begin(),tmp1->end()));
      ret->setCoords(coo2);
    }
  for(int i=0;i<dim;i++)
    {
      dims[i]=cellPart[i].second-cellPart[i].first+1;
      if(dims[i]<1)
        throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::buildStructuredSubPart : invalid input cellPart !");
    }
  ret->setNodeGridStructure(&dims[0],&dims[0]+dims.size());
  return ret.retn();
}

void MEDCouplingCurveLinearMesh::getBoundingBox(double *bbox) const
{
  if(!((const DataArrayDouble *)_coords))
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::getBoundingBox : Coordinates not set !");
  _coords->getMinMaxPerComponent(bbox);
}

MEDCouplingFieldDouble *MEDCouplingCurveLinearMesh::getMeasureField(bool isAbs) const
{
  checkConsistencyLight();
  int meshDim=getMeshDimension();
  std::string name="MeasureOfMesh_"; name+=getName();
  MCAuto<MEDCouplingFieldDouble> field=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  field->setName(name); field->setMesh(const_cast<MEDCouplingCurveLinearMesh *>(this)); field->synchronizeTimeWithMesh();
  switch(meshDim)
  {
    case 3:
      { getMeasureFieldMeshDim3(isAbs,field); return field.retn(); }
    case 2:
      { getMeasureFieldMeshDim2(isAbs,field); return field.retn(); }
    case 1:
      { getMeasureFieldMeshDim1(isAbs,field); return field.retn(); }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::getMeasureField : mesh dimension must be in [1,2,3] !");
  }
}

/*!
 * \param [in,out] f field feeded with good values.
 * \sa MEDCouplingCurveLinearMesh::getMeasureField
 */
void MEDCouplingCurveLinearMesh::getMeasureFieldMeshDim1(bool isAbs, MEDCouplingFieldDouble *field) const
{
  int nbnodes=getNumberOfNodes();
  int spaceDim=getSpaceDimension();
  MCAuto<DataArrayDouble> arr=DataArrayDouble::New(); field->setArray(arr);
  if(nbnodes==0)
    { arr->alloc(0,1); return; }
  if(spaceDim==1)
    {
      arr->alloc(nbnodes-1,1);
      std::transform(_coords->begin()+1,_coords->end(),_coords->begin(),arr->getPointer(),std::minus<double>());
      if(isAbs)
        arr->abs();
    }
  else
    {
      MCAuto<DataArrayDouble> tmp=DataArrayDouble::New(); tmp->alloc(nbnodes-1,spaceDim);
      std::transform(_coords->begin()+spaceDim,_coords->end(),_coords->begin(),tmp->getPointer(),std::minus<double>());
      MCAuto<DataArrayDouble> tmp2=tmp->magnitude(); field->setArray(tmp2);
    }
}

/*!
 * \param [in,out] f field feeded with good values.
 * \sa MEDCouplingCurveLinearMesh::getMeasureField
 */
void MEDCouplingCurveLinearMesh::getMeasureFieldMeshDim2(bool isAbs, MEDCouplingFieldDouble *field) const
{
  int nbcells=getNumberOfCells();
  int spaceDim=getSpaceDimension();
  if(spaceDim!=2 && spaceDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::getMeasureFieldMeshDim2 : with meshDim 2 only space dimension 2 and 3 are possible !");
  MCAuto<DataArrayDouble> arr=DataArrayDouble::New(); field->setArray(arr);
  arr->alloc(nbcells,1);
  double *pt=arr->getPointer();
  const double *coords=_coords->begin();
  int nX=_structure[0]-1;
  int conn[4];
  for(int i=0;i<nbcells;i++,pt++)
    {
      int cy=i/nX,cx=i-cy*nX;
      conn[0]=cy*(nX+1)+cx; conn[1]=(cy+1)*(nX+1)+cx; conn[2]=(cy+1)*(nX+1)+1+cx; conn[3]=cy*(nX+1)+cx+1;
      *pt=INTERP_KERNEL::computeVolSurfOfCell2<int,INTERP_KERNEL::ALL_C_MODE>(INTERP_KERNEL::NORM_QUAD4,conn,4,coords,spaceDim);
    }
  if(isAbs)
    arr->abs();
}

/*!
 * \param [in,out] f field feeded with good values.
 * \sa MEDCouplingCurveLinearMesh::getMeasureField
 */
void MEDCouplingCurveLinearMesh::getMeasureFieldMeshDim3(bool isAbs, MEDCouplingFieldDouble *field) const
{
  int nbcells=getNumberOfCells();
  int spaceDim=getSpaceDimension();
  if(spaceDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::getMeasureFieldMeshDim3 : with meshDim 3 only space dimension 3 is possible !");
  MCAuto<DataArrayDouble> arr=DataArrayDouble::New(); field->setArray(arr);
  arr->alloc(nbcells,1);
  double *pt=arr->getPointer();
  const double *coords=_coords->begin();
  int nX=_structure[0]-1,nY=(_structure[0]-1)*(_structure[1]-1);
  int nY1=_structure[0]*_structure[1];
  int conn[8];
  for(int i=0;i<nbcells;i++,pt++)
    {
      int cz=i/nY;
      int cy=(i-cz*nY)/nX;
      int cx=(i-cz*nY)-nX*cy;
      conn[0]=cz*nY1+cy*(nX+1)+cx; conn[1]=cz*nY1+(cy+1)*(nX+1)+cx; conn[2]=cz*nY1+(cy+1)*(nX+1)+1+cx; conn[3]=cz*nY1+cy*(nX+1)+cx+1;
      conn[4]=(cz+1)*nY1+cy*(nX+1)+cx; conn[5]=(cz+1)*nY1+(cy+1)*(nX+1)+cx; conn[6]=(cz+1)*nY1+(cy+1)*(nX+1)+1+cx; conn[7]=(cz+1)*nY1+cy*(nX+1)+cx+1;
      *pt=INTERP_KERNEL::computeVolSurfOfCell2<int,INTERP_KERNEL::ALL_C_MODE>(INTERP_KERNEL::NORM_HEXA8,conn,8,coords,3);
    }
  if(isAbs)
    arr->abs();
}

/*!
 * not implemented yet !
 */
MEDCouplingFieldDouble *MEDCouplingCurveLinearMesh::getMeasureFieldOnNode(bool isAbs) const
{
  throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::getMeasureFieldOnNode : not implemented yet !");
}

MEDCouplingFieldDouble *MEDCouplingCurveLinearMesh::buildOrthogonalField() const
{
  if(getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("Expected a cmesh with meshDim == 2 !");
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  DataArrayDouble *array=DataArrayDouble::New();
  int nbOfCells=getNumberOfCells();
  array->alloc(nbOfCells,3);
  double *vals=array->getPointer();
  for(int i=0;i<nbOfCells;i++)
    { vals[3*i]=0.; vals[3*i+1]=0.; vals[3*i+2]=1.; }
  ret->setArray(array);
  array->decrRef();
  ret->setMesh(this);
  return ret;
}

/// @cond INTERNAL

namespace MEDCoupling
{
  template<const int SPACEDIMM>
  class DummyClsMCL
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
}

/// @endcond

int MEDCouplingCurveLinearMesh::getCellContainingPoint(const double *pos, double eps) const
{
  checkConsistencyLight();
  int spaceDim=getSpaceDimension();
  const double *coords=_coords->getConstPointer();
  int nodeId=-1;
  _coords->distanceToTuple(pos,pos+spaceDim,nodeId);
  if(nodeId<0)
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::getCellContainingPoint : internal problem 1 !");
  int conn[8];
  int nbOfNodes=getNumberOfNodes();
  if(nbOfNodes==1)
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::getCellContainingPoint : No cells in this !");
  switch(getMeshDimension())
  {
    case 1:
      if(spaceDim==1)
        {
          if(nodeId>0)
            {
              conn[0]=nodeId-1; conn[1]=nodeId;
              if(INTERP_KERNEL::PointLocatorAlgos<DummyClsMCL<1> >::isElementContainsPoint(pos,INTERP_KERNEL::NORM_SEG2,coords,conn,2,eps))
                return nodeId-1;
            }
          if(nodeId<nbOfNodes-1)
            {
              conn[0]=nodeId; conn[1]=nodeId+1;
              if(INTERP_KERNEL::PointLocatorAlgos<DummyClsMCL<1> >::isElementContainsPoint(pos,INTERP_KERNEL::NORM_SEG2,coords,conn,2,eps))
                return nodeId;
            }
        }
    case 2:
      if(spaceDim==2)
        {
          int ny=nodeId/_structure[0],nx=nodeId-ny*_structure[0];
          if(nx>0 && ny>0)
            {
              conn[0]=nx-1+_structure[0]*(ny-1); conn[1]=nx-1+_structure[0]*ny; conn[2]=nx+_structure[0]*ny; conn[3]=nx+_structure[0]*(ny-1);
              if(INTERP_KERNEL::PointLocatorAlgos<DummyClsMCL<2> >::isElementContainsPoint(pos,INTERP_KERNEL::NORM_QUAD4,coords,conn,4,eps))
                return nx-1+(ny-1)*_structure[0];
            }
          if(nx<_structure[0]-1 && ny>0)
            {
              conn[0]=nx+_structure[0]*(ny-1); conn[1]=nx+_structure[0]*ny; conn[2]=nx+1+_structure[0]*ny; conn[3]=nx+1+_structure[0]*(ny-1);
              if(INTERP_KERNEL::PointLocatorAlgos<DummyClsMCL<2> >::isElementContainsPoint(pos,INTERP_KERNEL::NORM_QUAD4,coords,conn,4,eps))
                return nx+(ny-1)*_structure[0];
            }
          if(nx>0 && ny<_structure[1]-1)
            {
              conn[0]=nx-1+_structure[0]*ny; conn[1]=nx-1+_structure[0]*(ny+1); conn[2]=nx+_structure[0]*(ny+1); conn[3]=nx+_structure[0]*ny;
              if(INTERP_KERNEL::PointLocatorAlgos<DummyClsMCL<2> >::isElementContainsPoint(pos,INTERP_KERNEL::NORM_QUAD4,coords,conn,4,eps))
                return nx-1+ny*_structure[0];
            }
          if(nx<_structure[0]-1 && ny<_structure[1]-1)
            {
              conn[0]=nx+_structure[0]*ny; conn[1]=nx+_structure[0]*(ny+1); conn[2]=nx+1+_structure[0]*(ny+1); conn[3]=nx+1+_structure[0]*ny;
              if(INTERP_KERNEL::PointLocatorAlgos<DummyClsMCL<2> >::isElementContainsPoint(pos,INTERP_KERNEL::NORM_QUAD4,coords,conn,4,eps))
                return nx+ny*_structure[0];
            }
        }
    case 3:
      {
        if(spaceDim==3)
          {
            int nY=_structure[0]*_structure[1];
            int nz=nodeId/_structure[1]; int ny=(nodeId-nz*nY)/_structure[0]; int nx=(nodeId-nz*nY)-_structure[0]*ny;
            if(nx>0 && ny>0 && nz>0)
              {
                conn[0]=nx-1+_structure[0]*(ny-1)+nY*(nz-1); conn[1]=nx-1+_structure[2]*ny+nY*(nz-1); conn[2]=nx+_structure[2]*ny+nY*(nz-1); conn[3]=nx+_structure[0]*(ny-1)+nY*(nz-1);
                conn[4]=nx-1+_structure[0]*(ny-1)+nY*nz; conn[5]=nx-1+_structure[0]*ny+nY*nz; conn[6]=nx+_structure[0]*ny+nY*nz; conn[7]=nx+_structure[0]*(ny-1)+nY*nz;
                if(INTERP_KERNEL::PointLocatorAlgos<DummyClsMCL<3> >::isElementContainsPoint(pos,INTERP_KERNEL::NORM_HEXA8,coords,conn,8,eps))
                  return nx-1+(ny-1)*_structure[0]+(nz-1)*nY;
              }
            if(nx<_structure[0]-1 && ny>0 && nz>0)
              {
                conn[0]=nx+_structure[0]*(ny-1)+nY*(nz-1); conn[1]=nx+_structure[0]*ny+nY*(nz-1); conn[2]=nx+1+_structure[0]*ny+nY*(nz-1); conn[3]=nx+1+_structure[0]*(ny-1)+nY*(nz-1);
                conn[4]=nx+_structure[0]*(ny-1)+nY*nz; conn[5]=nx+_structure[0]*ny+nY*nz; conn[6]=nx+1+_structure[0]*ny+nY*nz; conn[7]=nx+1+_structure[0]*(ny-1)+nY*nz;
                if(INTERP_KERNEL::PointLocatorAlgos<DummyClsMCL<3> >::isElementContainsPoint(pos,INTERP_KERNEL::NORM_HEXA8,coords,conn,8,eps))
                  return nx+(ny-1)*_structure[0]+(nz-1)*nY;
              }
            if(nx>0 && ny<_structure[1]-1 && nz>0)
              {
                conn[0]=nx-1+_structure[0]*ny+nY*(nz-1); conn[1]=nx-1+_structure[0]*(ny+1)+nY*(nz-1); conn[2]=nx+_structure[0]*(ny+1)+nY*(nz-1); conn[3]=nx+_structure[0]*ny+nY*(nz-1);
                conn[4]=nx-1+_structure[0]*ny+nY*nz; conn[5]=nx-1+_structure[0]*(ny+1)+nY*nz; conn[6]=nx+_structure[0]*(ny+1)+nY*nz; conn[7]=nx+_structure[0]*ny+nY*nz;
                if(INTERP_KERNEL::PointLocatorAlgos<DummyClsMCL<3> >::isElementContainsPoint(pos,INTERP_KERNEL::NORM_HEXA8,coords,conn,8,eps))
                  return nx-1+ny*_structure[0]+(nz-1)*nY;
              }
            if(nx<_structure[0]-1 && ny<_structure[1]-1 && nz>0)
              {
                conn[0]=nx+_structure[0]*ny+nY*(nz-1); conn[1]=nx+_structure[0]*(ny+1)+nY*(nz-1); conn[2]=nx+1+_structure[0]*(ny+1)+nY*(nz-1); conn[3]=nx+1+_structure[0]*ny+nY*(nz-1);
                conn[4]=nx+_structure[0]*ny+nY*nz; conn[5]=nx+_structure[0]*(ny+1)+nY*nz; conn[6]=nx+1+_structure[0]*(ny+1)+nY*nz; conn[7]=nx+1+_structure[0]*ny+nY*nz;
                if(INTERP_KERNEL::PointLocatorAlgos<DummyClsMCL<3> >::isElementContainsPoint(pos,INTERP_KERNEL::NORM_HEXA8,coords,conn,8,eps))
                  return nx+ny*_structure[0]+(nz-1)*nY;
              }
            if(nx>0 && ny>0 && nz<_structure[2]-1)
              {
                conn[0]=nx-1+_structure[0]*(ny-1)+nY*(nz-1); conn[1]=nx-1+_structure[2]*ny+nY*(nz-1); conn[2]=nx+_structure[2]*ny+nY*(nz-1); conn[3]=nx+_structure[0]*(ny-1)+nY*(nz-1);
                conn[4]=nx-1+_structure[0]*(ny-1)+nY*nz; conn[5]=nx-1+_structure[0]*ny+nY*nz; conn[6]=nx+_structure[0]*ny+nY*nz; conn[7]=nx+_structure[0]*(ny-1)+nY*nz;
                if(INTERP_KERNEL::PointLocatorAlgos<DummyClsMCL<3> >::isElementContainsPoint(pos,INTERP_KERNEL::NORM_HEXA8,coords,conn,8,eps))
                  return nx-1+(ny-1)*_structure[0]+nz*nY;
              }
            if(nx<_structure[0]-1 && ny>0 && nz<_structure[2]-1)
              {
                conn[0]=nx+_structure[0]*(ny-1)+nY*(nz-1); conn[1]=nx+_structure[0]*ny+nY*(nz-1); conn[2]=nx+1+_structure[0]*ny+nY*(nz-1); conn[3]=nx+1+_structure[0]*(ny-1)+nY*(nz-1);
                conn[4]=nx+_structure[0]*(ny-1)+nY*nz; conn[5]=nx+_structure[0]*ny+nY*nz; conn[6]=nx+1+_structure[0]*ny+nY*nz; conn[7]=nx+1+_structure[0]*(ny-1)+nY*nz;
                if(INTERP_KERNEL::PointLocatorAlgos<DummyClsMCL<3> >::isElementContainsPoint(pos,INTERP_KERNEL::NORM_HEXA8,coords,conn,8,eps))
                  return nx+(ny-1)*_structure[0]+nz*nY;
              }
            if(nx>0 && ny<_structure[1]-1 && nz<_structure[2]-1)
              {
                conn[0]=nx-1+_structure[0]*ny+nY*(nz-1); conn[1]=nx-1+_structure[0]*(ny+1)+nY*(nz-1); conn[2]=nx+_structure[0]*(ny+1)+nY*(nz-1); conn[3]=nx+_structure[0]*ny+nY*(nz-1);
                conn[4]=nx-1+_structure[0]*ny+nY*nz; conn[5]=nx-1+_structure[0]*(ny+1)+nY*nz; conn[6]=nx+_structure[0]*(ny+1)+nY*nz; conn[7]=nx+_structure[0]*ny+nY*nz;
                if(INTERP_KERNEL::PointLocatorAlgos<DummyClsMCL<3> >::isElementContainsPoint(pos,INTERP_KERNEL::NORM_HEXA8,coords,conn,8,eps))
                  return nx-1+ny*_structure[0]+nz*nY;
              }
            if(nx<_structure[0]-1 && ny<_structure[1]-1 && nz<_structure[2]-1)
              {
                conn[0]=nx+_structure[0]*ny+nY*(nz-1); conn[1]=nx+_structure[0]*(ny+1)+nY*(nz-1); conn[2]=nx+1+_structure[0]*(ny+1)+nY*(nz-1); conn[3]=nx+1+_structure[0]*ny+nY*(nz-1);
                conn[4]=nx+_structure[0]*ny+nY*nz; conn[5]=nx+_structure[0]*(ny+1)+nY*nz; conn[6]=nx+1+_structure[0]*(ny+1)+nY*nz; conn[7]=nx+1+_structure[0]*ny+nY*nz;
                if(INTERP_KERNEL::PointLocatorAlgos<DummyClsMCL<3> >::isElementContainsPoint(pos,INTERP_KERNEL::NORM_HEXA8,coords,conn,8,eps))
                  return nx+ny*_structure[0]+nz*nY;
              }
          }
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::getCellContainingPoint : mesh dimension managed are 1, 2 or 3 !");
  }
}

void MEDCouplingCurveLinearMesh::getCellsContainingPoint(const double *pos, double eps, std::vector<int>& elts) const
{
  int ret(getCellContainingPoint(pos,eps));
  elts.push_back(ret);
}

void MEDCouplingCurveLinearMesh::rotate(const double *center, const double *vector, double angle)
{
  if(!((DataArrayDouble *)_coords))
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::rotate : no coordinates set !");
  int spaceDim=getSpaceDimension();
  int nbNodes=_coords->getNumberOfTuples();
  double *coords=_coords->getPointer();
  if(spaceDim==3)
    DataArrayDouble::Rotate3DAlg(center,vector,angle,nbNodes,coords,coords);
  else if(spaceDim==2)
    DataArrayDouble::Rotate2DAlg(center,angle,nbNodes,coords,coords);
  else
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::rotate : invalid space dim for rotation must be 2 or 3");
  _coords->declareAsNew();
  updateTime();
}

void MEDCouplingCurveLinearMesh::translate(const double *vector)
{
  if(!vector)
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::translate : NULL input point !");
  if(!((DataArrayDouble *)_coords))
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::translate : no coordinates set !");
  double *coords=_coords->getPointer();
  int nbNodes=getNumberOfNodes();
  int dim=getSpaceDimension();
  for(int i=0; i<nbNodes; i++)
    for(int idim=0; idim<dim;idim++)
      coords[i*dim+idim]+=vector[idim];
  _coords->declareAsNew();
  updateTime();
}

void MEDCouplingCurveLinearMesh::scale(const double *point, double factor)
{
  if(!point)
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::scale : NULL input point !");
  if(!((DataArrayDouble *)_coords))
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::scale : no coordinates set !");
  double *coords=_coords->getPointer();
  int nbNodes=_coords->getNumberOfTuples();
  int dim=_coords->getNumberOfComponents();
  for(int i=0;i<nbNodes;i++)
    {
      std::transform(coords+i*dim,coords+(i+1)*dim,point,coords+i*dim,std::minus<double>());
      std::transform(coords+i*dim,coords+(i+1)*dim,coords+i*dim,std::bind2nd(std::multiplies<double>(),factor));
      std::transform(coords+i*dim,coords+(i+1)*dim,point,coords+i*dim,std::plus<double>());
    }
  _coords->declareAsNew();
  updateTime();
}

MEDCouplingMesh *MEDCouplingCurveLinearMesh::mergeMyselfWith(const MEDCouplingMesh *other) const
{
  throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::mergeMyselfWith : not available for CurveLinear Mesh !");
}

DataArrayDouble *MEDCouplingCurveLinearMesh::getCoordinatesAndOwner() const
{
  DataArrayDouble *ret=const_cast<DataArrayDouble *>((const DataArrayDouble *)_coords);
  if(ret)
    ret->incrRef();
  return ret;
}

DataArrayDouble *MEDCouplingCurveLinearMesh::computeCellCenterOfMass() const
{
  checkConsistencyLight();
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  int spaceDim=getSpaceDimension();
  int meshDim=getMeshDimension();
  int nbOfCells=getNumberOfCells();
  ret->alloc(nbOfCells,spaceDim);
  ret->copyStringInfoFrom(*getCoords());
  switch(meshDim)
  {
    case 3:
      { getBarycenterAndOwnerMeshDim3(ret); return ret.retn(); }
    case 2:
      { getBarycenterAndOwnerMeshDim2(ret); return ret.retn(); }
    case 1:
      { getBarycenterAndOwnerMeshDim1(ret); return ret.retn(); }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::computeCellCenterOfMass : mesh dimension must be in [1,2,3] !");
  }
}

DataArrayDouble *MEDCouplingCurveLinearMesh::computeIsoBarycenterOfNodesPerCell() const
{
  return MEDCouplingCurveLinearMesh::computeCellCenterOfMass();
}

/*!
 * \param [in,out] bary Barycenter array feeded with good values.
 * \sa MEDCouplingCurveLinearMesh::computeCellCenterOfMass
 */
void MEDCouplingCurveLinearMesh::getBarycenterAndOwnerMeshDim3(DataArrayDouble *bary) const
{
  int nbOfCells=getNumberOfCells();
  double *ptToFill=bary->getPointer();
  const double *coor=_coords->getConstPointer();
  if(getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::getBarycenterAndOwnerMeshDim3 : with meshDim 3 only space dimension 3 is possible !");
  int nX=_structure[0]-1,nY=(_structure[0]-1)*(_structure[1]-1);
  int nY1=_structure[0]*_structure[1];
  int conn[8];
  for(int i=0;i<nbOfCells;i++)
    {
      int cz=i/nY;
      int cy=(i-cz*nY)/nX;
      int cx=(i-cz*nY)-nX*cy;
      conn[0]=cz*nY1+cy*(nX+1)+cx+1; conn[1]=cz*nY1+cy*(nX+1)+cx; conn[2]=cz*nY1+(cy+1)*(nX+1)+cx; conn[3]=cz*nY1+(cy+1)*(nX+1)+1+cx;
      conn[4]=(cz+1)*nY1+cy*(nX+1)+cx+1; conn[5]=(cz+1)*nY1+cy*(nX+1)+cx; conn[6]=(cz+1)*nY1+(cy+1)*(nX+1)+cx; conn[7]=(cz+1)*nY1+(cy+1)*(nX+1)+1+cx;
      INTERP_KERNEL::computeBarycenter2<int,INTERP_KERNEL::ALL_C_MODE>(INTERP_KERNEL::NORM_HEXA8,conn,8,coor,3,ptToFill);
      ptToFill+=3;
    }
}

/*!
 * \param [in,out] bary Barycenter array feeded with good values.
 * \sa MEDCouplingCurveLinearMesh::computeCellCenterOfMass
 */
void MEDCouplingCurveLinearMesh::getBarycenterAndOwnerMeshDim2(DataArrayDouble *bary) const
{
  int nbcells=getNumberOfCells();
  int spaceDim=getSpaceDimension();
  double *ptToFill=bary->getPointer();
  const double *coor=_coords->getConstPointer();
  if(spaceDim!=2 && spaceDim!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::getBarycenterAndOwnerMeshDim2 : with meshDim 2 only space dimension 2 and 3 are possible !");
  int nX=_structure[0]-1;
  int conn[4];
  for(int i=0;i<nbcells;i++)
    {
      int cy=i/nX,cx=i-cy*nX;
      conn[0]=cy*(nX+1)+cx; conn[1]=(cy+1)*(nX+1)+cx; conn[2]=(cy+1)*(nX+1)+1+cx; conn[3]=cy*(nX+1)+cx+1;
      INTERP_KERNEL::computeBarycenter2<int,INTERP_KERNEL::ALL_C_MODE>(INTERP_KERNEL::NORM_QUAD4,conn,4,coor,spaceDim,ptToFill);
      ptToFill+=spaceDim;
    }
}

/*!
 * \param [in,out] bary Barycenter array feeded with good values.
 * \sa MEDCouplingCurveLinearMesh::computeCellCenterOfMass
 */
void MEDCouplingCurveLinearMesh::getBarycenterAndOwnerMeshDim1(DataArrayDouble *bary) const
{
  int spaceDim=getSpaceDimension();
  std::transform(_coords->begin()+spaceDim,_coords->end(),_coords->begin(),bary->getPointer(),std::plus<double>());
  std::transform(bary->begin(),bary->end(),bary->getPointer(),std::bind2nd(std::multiplies<double>(),0.5));
}

void MEDCouplingCurveLinearMesh::renumberCells(const int *old2NewBg, bool check)
{
  throw INTERP_KERNEL::Exception("Functionnality of renumbering cell not available for CurveLinear Mesh !");
}

void MEDCouplingCurveLinearMesh::getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const
{
  int it,order;
  double time=getTime(it,order);
  tinyInfo.clear();
  tinyInfoD.clear();
  littleStrings.clear();
  littleStrings.push_back(getName());
  littleStrings.push_back(getDescription());
  littleStrings.push_back(getTimeUnit());
  //
  std::vector<std::string> littleStrings2;
  if((const DataArrayDouble *)_coords)
    _coords->getTinySerializationStrInformation(littleStrings2);
  littleStrings.insert(littleStrings.end(),littleStrings2.begin(),littleStrings2.end());
  //
  tinyInfo.push_back(it);
  tinyInfo.push_back(order);
  tinyInfo.push_back((int)_structure.size());
  for(std::vector<int>::const_iterator itt=_structure.begin();itt!=_structure.end();itt++)
    tinyInfo.push_back(*itt);
  std::vector<int> tinyInfo2;
  if((const DataArrayDouble *)_coords)
    _coords->getTinySerializationIntInformation(tinyInfo2);
  tinyInfo.insert(tinyInfo.end(),tinyInfo2.begin(),tinyInfo2.end());
  //
  tinyInfoD.push_back(time);
}

void MEDCouplingCurveLinearMesh::resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const
{
  a1->alloc(tinyInfo[2],1);
  std::vector<int> tinyInfo2(tinyInfo.begin()+3+tinyInfo[2],tinyInfo.end());
  a2->resizeForUnserialization(tinyInfo2);
}

void MEDCouplingCurveLinearMesh::serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const
{
  a1=DataArrayInt::New();
  a1->alloc((int)_structure.size(),1);
  int *ptr=a1->getPointer();
  for(std::vector<int>::const_iterator it=_structure.begin();it!=_structure.end();it++,ptr++)
    *ptr=(*it);
  int sz=0;
  if((const DataArrayDouble *)_coords)
    if(_coords->isAllocated())
      sz=_coords->getNbOfElems();
  a2=DataArrayDouble::New();
  a2->alloc(sz,1);
  if(sz!=0 && (const DataArrayDouble *)_coords)
    std::copy(_coords->begin(),_coords->end(),a2->getPointer());
}

void MEDCouplingCurveLinearMesh::unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                                                 const std::vector<std::string>& littleStrings)
{
  setName(littleStrings[0]);
  setDescription(littleStrings[1]);
  setTimeUnit(littleStrings[2]);
  setTime(tinyInfoD[0],tinyInfo[0],tinyInfo[1]);
  int sz=tinyInfo[2];
  _structure.resize(sz);
  for(int i=0;i<sz;i++)
    _structure[i]=tinyInfo[3+i];
  if((int)tinyInfo.size()>sz+3)
    {
      _coords=DataArrayDouble::New();
      std::vector<int> tinyInfo2(tinyInfo.begin()+3+sz,tinyInfo.end());
      _coords->resizeForUnserialization(tinyInfo2);
      std::copy(a2->begin(),a2->end(),_coords->getPointer());
      std::vector<std::string> littleStrings2(littleStrings.begin()+3,littleStrings.end());
      _coords->finishUnserialization(tinyInfo2,littleStrings2);
    }
}

void MEDCouplingCurveLinearMesh::writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const
{
  std::ostringstream extent;
  int meshDim=(int)_structure.size();
  if(meshDim<=0 || meshDim>3)
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::writeVTKLL : meshDim invalid ! must be in [1,2,3] !");
  for(int i=0;i<3;i++)
    { int val=i<meshDim?_structure[i]-1:0; extent << "0 " <<  val << " "; }
  ofs << "  <" << getVTKDataSetType() << " WholeExtent=\"" << extent.str() << "\">\n";
  ofs << "    <Piece Extent=\"" << extent.str() << "\">\n";
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
  ofs << "    </Piece>\n";
  ofs << "  </" << getVTKDataSetType() << ">\n";
}

void MEDCouplingCurveLinearMesh::reprQuickOverview(std::ostream& stream) const
{
  stream << "MEDCouplingCurveLinearMesh C++ instance at " << this << ". Name : \"" << getName() << "\".";
  stream << " Nodal structure : [";
  for(std::size_t i=0;i<_structure.size();i++)
    {
      char tmp='X'+i;
      stream << " " << tmp << "=" << _structure[i];
      if(i!=_structure.size()-1)
        stream << ", ";
    }
  stream << " ].";
  const DataArrayDouble *coo(_coords);
  if(!coo)
    { stream << std::endl << "No coordinates set !"; return ; }
  if(!coo->isAllocated())
    { stream << std::endl << "Coordinates set but not allocated !"; return ; }
  int nbOfCompo(coo->getNumberOfComponents());
  int nbOfCompoExp(-1);
  try
    {
      nbOfCompoExp=getSpaceDimension();
    }
  catch(INTERP_KERNEL::Exception& e)
    {
    }
  if(nbOfCompo!=nbOfCompoExp)
    { stream << std::endl << "Coordinates set and allocated but mismatch number of components !"; return ; }
  stream << std::endl << "Coordinates ( number of tuples = " << coo->getNumberOfTuples() << " ) : ";
  coo->reprQuickOverviewData(stream,200);
}

std::string MEDCouplingCurveLinearMesh::getVTKFileExtension() const
{
  return std::string("vts");
}

std::string MEDCouplingCurveLinearMesh::getVTKDataSetType() const
{
  return std::string("StructuredGrid");
}
