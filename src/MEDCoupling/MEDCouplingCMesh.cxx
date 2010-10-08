//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include <functional>
#include <algorithm>
#include <numeric>

using namespace ParaMEDMEM;

MEDCouplingCMesh::MEDCouplingCMesh():_x_array(0),_y_array(0),_z_array(0)
{
}

MEDCouplingCMesh::MEDCouplingCMesh(const MEDCouplingCMesh& other, bool deepCpy):MEDCouplingMesh(other)
{
  if(deepCpy)
    {
      if(other._x_array)
        _x_array=_x_array->deepCopy();
      if(other._y_array)
        _y_array=_y_array->deepCopy();
      if(other._z_array)
        _z_array=_z_array->deepCopy();
    }
  else
    {
      _x_array=other._x_array;
      if(_x_array)
        _x_array->incrRef();
      _y_array=other._y_array;
      if(_y_array)
        _y_array->incrRef();
      _z_array=other._z_array;
      if(_z_array)
        _z_array->incrRef();
    }
}

MEDCouplingCMesh::~MEDCouplingCMesh()
{
  if(_x_array)
    _x_array->decrRef();
  if(_y_array)
    _y_array->decrRef();
  if(_z_array)
    _z_array->decrRef();
}

MEDCouplingCMesh *MEDCouplingCMesh::New()
{
  return new MEDCouplingCMesh;
}

MEDCouplingMesh *MEDCouplingCMesh::deepCpy() const
{
  return clone(true);
}

MEDCouplingCMesh *MEDCouplingCMesh::clone(bool recDeepCpy) const
{
  return new MEDCouplingCMesh(*this,recDeepCpy);
}

void MEDCouplingCMesh::updateTime()
{
  if(_x_array)
    updateTimeWith(*_x_array);
  if(_y_array)
    updateTimeWith(*_y_array);
  if(_z_array)
    updateTimeWith(*_z_array);
}

/*!
 * This method copyies all tiny strings from other (name and components name).
 * @throw if other and this have not same mesh type.
 */
void MEDCouplingCMesh::copyTinyStringsFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception)
{ 
  const MEDCouplingCMesh *otherC=dynamic_cast<const MEDCouplingCMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCouplingCMesh::copyTinyStringsFrom : meshes have not same type !");
  MEDCouplingMesh::copyTinyStringsFrom(other);
  if(_x_array && otherC->_x_array)
    _x_array->copyStringInfoFrom(*otherC->_x_array);
  if(_y_array && otherC->_y_array)
    _y_array->copyStringInfoFrom(*otherC->_y_array);
  if(_z_array && otherC->_z_array)
    _z_array->copyStringInfoFrom(*otherC->_z_array);
}

bool MEDCouplingCMesh::isEqual(const MEDCouplingMesh *other, double prec) const
{
  const MEDCouplingCMesh *otherC=dynamic_cast<const MEDCouplingCMesh *>(other);
  if(!otherC)
    return false;
  if(!MEDCouplingMesh::isEqual(other,prec))
    return false;
  const DataArrayDouble *thisArr[3]={_x_array,_y_array,_z_array};
  const DataArrayDouble *otherArr[3]={otherC->_x_array,otherC->_y_array,otherC->_z_array};
  for(int i=0;i<3;i++)
    {
      if((thisArr[i]!=0 && otherArr[i]==0) || (thisArr[i]==0 && otherArr[i]!=0))
        return false;
      if(thisArr[i])
        if(!thisArr[i]->isEqual(*otherArr[i],prec))
          return false;
    }
  return true;
}

bool MEDCouplingCMesh::isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const
{
  const MEDCouplingCMesh *otherC=dynamic_cast<const MEDCouplingCMesh *>(other);
  if(!otherC)
    return false;
  const DataArrayDouble *thisArr[3]={_x_array,_y_array,_z_array};
  const DataArrayDouble *otherArr[3]={otherC->_x_array,otherC->_y_array,otherC->_z_array};
  for(int i=0;i<3;i++)
    {
      if((thisArr[i]!=0 && otherArr[i]==0) || (thisArr[i]==0 && otherArr[i]!=0))
        return false;
      if(thisArr[i])
        if(!thisArr[i]->isEqualWithoutConsideringStr(*otherArr[i],prec))
          return false;
    }
  return true;
}

void MEDCouplingCMesh::checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                            DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

void MEDCouplingCMesh::checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                       DataArrayInt *&cellCor) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

void MEDCouplingCMesh::checkCoherency() const throw(INTERP_KERNEL::Exception)
{
  const char msg0[]="Invalid ";
  const char msg1[]=" array ! Must contain more than 1 element.";
  const char msg2[]=" array ! Must be with only one component.";
  if(_x_array)
    {
      if(_x_array->getNbOfElems()<2)
        {
          std::ostringstream os; os << msg0 << 'X' << msg1;
          throw INTERP_KERNEL::Exception(os.str().c_str());
        }
      if(_x_array->getNumberOfComponents()!=1)
        {
          std::ostringstream os; os << msg0 << 'X' << msg2;
          throw INTERP_KERNEL::Exception(os.str().c_str());
        }
    }
  if(_y_array)
    {
      if(_y_array->getNbOfElems()<2)
        {
          std::ostringstream os; os << msg0 << 'Y' << msg1;
          throw INTERP_KERNEL::Exception(os.str().c_str());
        }
      if(_y_array->getNumberOfComponents()!=1)
        {
          std::ostringstream os; os << msg0 << 'Y' << msg2;
          throw INTERP_KERNEL::Exception(os.str().c_str());
        }
      
    }
  if(_z_array)
    {
      if(_z_array->getNbOfElems()<2)
        {
          std::ostringstream os; os << msg0 << 'Z' << msg1;
          throw INTERP_KERNEL::Exception(os.str().c_str());
        }
      if(_z_array->getNumberOfComponents()!=1)
        {
          std::ostringstream os; os << msg0 << 'Z' << msg2;
          throw INTERP_KERNEL::Exception(os.str().c_str());
        }
    }
}

int MEDCouplingCMesh::getNumberOfCells() const
{
  int ret=1;
  if(_x_array)
    ret*=_x_array->getNbOfElems()-1;
  if(_y_array)
    ret*=_y_array->getNbOfElems()-1;
  if(_z_array)
    ret*=_z_array->getNbOfElems()-1;
  return ret;
}

int MEDCouplingCMesh::getNumberOfNodes() const
{
  int ret=1;
  if(_x_array)
    ret*=_x_array->getNbOfElems();
  if(_y_array)
    ret*=_y_array->getNbOfElems();
  if(_z_array)
    ret*=_z_array->getNbOfElems();
  return ret;
}

void MEDCouplingCMesh::getSplitCellValues(int *res) const
{
  int spaceDim=getSpaceDimension();
  for(int l=0;l<spaceDim;l++)
    {
      int val=1;
      for(int p=l;p<spaceDim-1;p++)
        val*=getCoordsAt(p)->getNbOfElems()-1;
      res[spaceDim-l-1]=val;
    }
}

void MEDCouplingCMesh::getSplitNodeValues(int *res) const
{
  int spaceDim=getSpaceDimension();
  for(int l=0;l<spaceDim;l++)
    {
      int val=1;
      for(int p=l;p<spaceDim-1;p++)
        val*=getCoordsAt(p)->getNbOfElems();
      res[spaceDim-l-1]=val;
    }
}

int MEDCouplingCMesh::getCellIdFromPos(int i, int j, int k) const
{
  int tmp[3]={i,j,k};
  int tmp2[3];
  int spaceDim=getSpaceDimension();
  getSplitCellValues(tmp2);
  std::transform(tmp,tmp+spaceDim,tmp2,tmp,std::multiplies<int>());
  return std::accumulate(tmp,tmp+spaceDim,0);
}

int MEDCouplingCMesh::getNodeIdFromPos(int i, int j, int k) const
{
  int tmp[3]={i,j,k};
  int tmp2[3];
  int spaceDim=getSpaceDimension();
  getSplitNodeValues(tmp2);
  std::transform(tmp,tmp+spaceDim,tmp2,tmp,std::multiplies<int>());
  return std::accumulate(tmp,tmp+spaceDim,0);
}

void MEDCouplingCMesh::getPosFromId(int nodeId, int spaceDim, const int *split, int *res)
{
  int work=nodeId;
  for(int i=spaceDim-1;i>=0;i--)
    {
      int pos=work/split[i];
      work=work%split[i];
      res[i]=pos;
    }
}

int MEDCouplingCMesh::getSpaceDimension() const
{
  int ret=0;
  if(_x_array)
    ret++;
  if(_y_array)
    ret++;
  if(_z_array)
    ret++;
  return ret;
}

int MEDCouplingCMesh::getMeshDimension() const
{
  return getSpaceDimension();
}

INTERP_KERNEL::NormalizedCellType MEDCouplingCMesh::getTypeOfCell(int cellId) const
{
  switch(getMeshDimension())
    {
    case 3:
      return INTERP_KERNEL::NORM_HEXA8;
    case 2:
      return INTERP_KERNEL::NORM_QUAD4;
    case 1:
      return INTERP_KERNEL::NORM_SEG2;
    default:
      throw INTERP_KERNEL::Exception("Unexpected dimension for MEDCouplingCMesh::getTypeOfCell !");
    }
}

int MEDCouplingCMesh::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const
{
  int ret=getNumberOfCells();
  int dim=getMeshDimension();
  switch(type)
    {
    case INTERP_KERNEL::NORM_HEXA8:
      if(dim==3)
        return ret;
    case INTERP_KERNEL::NORM_QUAD4:
      if(dim==2)
      return ret;
    case INTERP_KERNEL::NORM_SEG2:
      if(dim==1)
        return ret;
    default:
      throw INTERP_KERNEL::Exception("Unexpected dimension for MEDCouplingCMesh::getTypeOfCell !");
    }
  return 0;
}

void MEDCouplingCMesh::getNodeIdsOfCell(int cellId, std::vector<int>& conn) const
{
  //not implemented yet
}

void MEDCouplingCMesh::getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const
{
  //not implemented yet
}

std::string MEDCouplingCMesh::simpleRepr() const
{
  std::ostringstream ret;
  ret << "Cartesian mesh with name : \"" << getName() << "\"\n";
  ret << "Mesh and SpaceDimension dimension : " << getSpaceDimension() << "\n\nArrays :\n________\n\n";
  if(_x_array)
    {
      ret << "X Array :\n";
      _x_array->reprZipWithoutNameStream(ret);
    }
  if(_y_array)
    {
      ret << "Y Array :\n";
      _y_array->reprZipWithoutNameStream(ret);
    }
  if(_z_array)
    {
      ret << "Z Array :\n";
      _z_array->reprZipWithoutNameStream(ret);
    }
  return ret.str();
}

std::string MEDCouplingCMesh::advancedRepr() const
{
  return simpleRepr();
}

DataArrayDouble *MEDCouplingCMesh::getCoordsAt(int i) const throw(INTERP_KERNEL::Exception)
{
  switch(i)
    {
    case 0:
      return _x_array;
    case 1:
      return _y_array;
    case 2:
      return _z_array;
    default:
      throw INTERP_KERNEL::Exception("Invalid rank specified must be 0 or 1 or 2.");
    }
}

void MEDCouplingCMesh::setCoordsAt(int i, DataArrayDouble *arr) throw(INTERP_KERNEL::Exception)
{
  DataArrayDouble **thisArr[3]={&_x_array,&_y_array,&_z_array};
  if(i<0 || i>2)
    throw INTERP_KERNEL::Exception("Invalid rank specified must be 0 or 1 or 2.");
  if(arr!=*(thisArr[i]))
    {
      if(*(thisArr[i]))
        (*(thisArr[i]))->decrRef();
      (*(thisArr[i]))=arr;
      if(*(thisArr[i]))
        (*(thisArr[i]))->incrRef();
      declareAsNew();
    }
}

void MEDCouplingCMesh::setCoords(DataArrayDouble *coordsX, DataArrayDouble *coordsY, DataArrayDouble *coordsZ)
{
  if(_x_array)
    _x_array->decrRef();
  _x_array=coordsX;
  if(_x_array)
    _x_array->incrRef();
  if(_y_array)
    _y_array->decrRef();
  _y_array=coordsY;
  if(_y_array)
    _y_array->incrRef();
  if(_z_array)
    _z_array->decrRef();
  _z_array=coordsZ;
  if(_z_array)
    _z_array->incrRef();
  declareAsNew();
}

MEDCouplingMesh *MEDCouplingCMesh::buildPart(const int *start, const int *end) const
{
  //not implemented yet !
  return 0;
}

MEDCouplingMesh *MEDCouplingCMesh::buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const
{
  //not implemented yet !
  return 0;
}

void MEDCouplingCMesh::getBoundingBox(double *bbox) const
{
  //not implemented yet !
}

MEDCouplingFieldDouble *MEDCouplingCMesh::getMeasureField(bool isAbs) const
{
  //not implemented yet !
  return 0;
}

MEDCouplingFieldDouble *MEDCouplingCMesh::getMeasureFieldOnNode(bool isAbs) const
{
  //not implemented yet !
  return 0;
}

MEDCouplingFieldDouble *MEDCouplingCMesh::buildOrthogonalField() const
{
  if(getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("Expected a cmesh with meshDim == 2 !");
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  DataArrayDouble *array=DataArrayDouble::New();
  int nbOfCells=getNumberOfCells();
  array->alloc(nbOfCells,3);
  double *vals=array->getPointer();
  for(int i=0;i<nbOfCells;i++)
    { vals[3*i]=1.; vals[3*i+1]=1.; vals[3*i+2]=1.; }
  ret->setArray(array);
  array->decrRef();
  ret->setMesh(this);
  return ret;
}

int MEDCouplingCMesh::getCellContainingPoint(const double *pos, double eps) const
{
  //not implemented yet !
  return -1;
}

void MEDCouplingCMesh::rotate(const double *center, const double *vector, double angle)
{
  throw INTERP_KERNEL::Exception("No rotation available on CMesh : Traduce it to StructuredMesh to apply it !");
}

void MEDCouplingCMesh::translate(const double *vector)
{
  if(_x_array)
    std::transform(_x_array->getConstPointer(),_x_array->getConstPointer()+_x_array->getNbOfElems(),
                   _x_array->getPointer(),std::bind2nd(std::plus<double>(),vector[0]));
  if(_y_array)
    std::transform(_y_array->getConstPointer(),_y_array->getConstPointer()+_y_array->getNbOfElems(),
                   _y_array->getPointer(),std::bind2nd(std::plus<double>(),vector[1]));
  if(_z_array)
    std::transform(_z_array->getConstPointer(),_z_array->getConstPointer()+_z_array->getNbOfElems(),
                   _z_array->getPointer(),std::bind2nd(std::plus<double>(),vector[2]));
}

MEDCouplingMesh *MEDCouplingCMesh::mergeMyselfWith(const MEDCouplingMesh *other) const
{
  //not implemented yet !
  return 0;
}

DataArrayDouble *MEDCouplingCMesh::getCoordinatesAndOwner() const
{
  DataArrayDouble *ret=DataArrayDouble::New();
  int spaceDim=getSpaceDimension();
  int nbNodes=getNumberOfNodes();
  ret->alloc(nbNodes,spaceDim);
  double *pt=ret->getPointer();
  int tmp[3];
  getSplitNodeValues(tmp);
  DataArrayDouble *tabs[3]={getCoordsAt(0),getCoordsAt(1),getCoordsAt(2)};
  const double *tabsPtr[3];
  for(int j=0;j<spaceDim;j++)
    tabsPtr[j]=tabs[j]->getConstPointer();
  int tmp2[3];
  for(int i=0;i<nbNodes;i++)
    {
      getPosFromId(i,spaceDim,tmp,tmp2);
      for(int j=0;j<spaceDim;j++)
        pt[i*spaceDim+j]=tabsPtr[j][tmp2[j]];
    }
  return ret;
}

DataArrayDouble *MEDCouplingCMesh::getBarycenterAndOwner() const
{
  DataArrayDouble *ret=DataArrayDouble::New();
  int spaceDim=getSpaceDimension();
  int nbCells=getNumberOfCells();
  ret->alloc(nbCells,spaceDim);
  double *pt=ret->getPointer();
  int tmp[3];
  getSplitCellValues(tmp);
  DataArrayDouble *tabs[3]={getCoordsAt(0),getCoordsAt(1),getCoordsAt(2)};
  std::vector<double> tabsPtr[3];
  for(int j=0;j<spaceDim;j++)
    {
      int sz=tabs[j]->getNbOfElems()-1;
      const double *srcPtr=tabs[j]->getConstPointer();
      tabsPtr[j].insert(tabsPtr[j].end(),srcPtr,srcPtr+sz);
      std::transform(tabsPtr[j].begin(),tabsPtr[j].end(),srcPtr+1,tabsPtr[j].begin(),std::plus<double>());
      std::transform(tabsPtr[j].begin(),tabsPtr[j].end(),tabsPtr[j].begin(),std::bind2nd(std::multiplies<double>(),0.5));
    }
  int tmp2[3];
  for(int i=0;i<nbCells;i++)
    {
      getPosFromId(i,spaceDim,tmp,tmp2);
      for(int j=0;j<spaceDim;j++)
        pt[i*spaceDim+j]=tabsPtr[j][tmp2[j]];
    }
  return ret;
}

void MEDCouplingCMesh::renumberCells(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Functionnality of renumbering cell not available for CMesh !");
}

void MEDCouplingCMesh::getTinySerializationInformation(std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const
{
  tinyInfo.clear();
  littleStrings.clear();
  littleStrings.push_back(getName());
  const DataArrayDouble *thisArr[3]={_x_array,_y_array,_z_array};
  for(int i=0;i<3;i++)
    {
      int val=-1;
      std::string st;
      if(thisArr[i])
        {
          val=thisArr[i]->getNumberOfTuples();
          st=thisArr[i]->getInfoOnComponent(0);
        }
      tinyInfo.push_back(val);
      littleStrings.push_back(st);
    }
}

void MEDCouplingCMesh::resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const
{
  a1->alloc(0,1);
  int sum=0;
  for(int i=0;i<3;i++)
    if(tinyInfo[i]!=-1)
      sum+=tinyInfo[i];
  a2->alloc(sum,1);
}

void MEDCouplingCMesh::serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const
{
  a1=DataArrayInt::New();
  a1->alloc(0,1);
  const DataArrayDouble *thisArr[3]={_x_array,_y_array,_z_array};
  int sz=0;
  for(int i=0;i<3;i++)
    {
      if(thisArr[i])
        sz+=thisArr[i]->getNumberOfTuples();
    }
  a2=DataArrayDouble::New();
  a2->alloc(sz,1);
  double *a2Ptr=a2->getPointer();
  for(int i=0;i<3;i++)
    if(thisArr[i])
      a2Ptr=std::copy(thisArr[i]->getConstPointer(),thisArr[i]->getConstPointer()+thisArr[i]->getNumberOfTuples(),a2Ptr);
}

void MEDCouplingCMesh::unserialization(const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                                       const std::vector<std::string>& littleStrings)
{
  setName(littleStrings[0].c_str());
  DataArrayDouble **thisArr[3]={&_x_array,&_y_array,&_z_array};
  const double *data=a2->getConstPointer();
  for(int i=0;i<3;i++)
    {
      if(tinyInfo[i]!=-1)
        {
          (*(thisArr[i]))=DataArrayDouble::New();
          (*(thisArr[i]))->alloc(tinyInfo[i],1);
          (*(thisArr[i]))->setInfoOnComponent(0,littleStrings[i+1].c_str());
          std::copy(data,data+tinyInfo[i],(*(thisArr[i]))->getPointer());
          data+=tinyInfo[i];
        }
    }
}

