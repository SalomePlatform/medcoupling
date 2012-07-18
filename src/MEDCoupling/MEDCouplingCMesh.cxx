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

#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include <functional>
#include <algorithm>
#include <sstream>
#include <numeric>

using namespace ParaMEDMEM;

MEDCouplingCMesh::MEDCouplingCMesh():_x_array(0),_y_array(0),_z_array(0)
{
}

MEDCouplingCMesh::MEDCouplingCMesh(const MEDCouplingCMesh& other, bool deepCopy):MEDCouplingMesh(other)
{
  if(deepCopy)
    {
      if(other._x_array)
        _x_array=other._x_array->deepCpy();
      else
        _x_array=0;
      if(other._y_array)
        _y_array=other._y_array->deepCpy();
      else
        _y_array=0;
      if(other._z_array)
        _z_array=other._z_array->deepCpy();
      else
        _z_array=0;
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

MEDCouplingCMesh *MEDCouplingCMesh::New(const char *meshName)
{
  MEDCouplingCMesh *ret=new MEDCouplingCMesh;
  ret->setName(meshName);
  return ret;
}

MEDCouplingMesh *MEDCouplingCMesh::deepCpy() const
{
  return clone(true);
}

MEDCouplingCMesh *MEDCouplingCMesh::clone(bool recDeepCpy) const
{
  return new MEDCouplingCMesh(*this,recDeepCpy);
}

void MEDCouplingCMesh::updateTime() const
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

bool MEDCouplingCMesh::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingCMesh::isEqualIfNotWhy : input other pointer is null !");
  const MEDCouplingCMesh *otherC=dynamic_cast<const MEDCouplingCMesh *>(other);
  if(!otherC)
    {
      reason="mesh given in input is not castable in MEDCouplingCMesh !";
      return false;
    }
  if(!MEDCouplingMesh::isEqualIfNotWhy(other,prec,reason))
    return false;
  const DataArrayDouble *thisArr[3]={_x_array,_y_array,_z_array};
  const DataArrayDouble *otherArr[3]={otherC->_x_array,otherC->_y_array,otherC->_z_array};
  std::ostringstream oss; oss.precision(15);
  for(int i=0;i<3;i++)
    {
      if((thisArr[i]!=0 && otherArr[i]==0) || (thisArr[i]==0 && otherArr[i]!=0))
        {
          oss << "Only one CMesh between the two this and other has its coordinates of rank" << i << " defined !";
          reason=oss.str();
          return false;
        }
      if(thisArr[i])
        if(!thisArr[i]->isEqualIfNotWhy(*otherArr[i],prec,reason))
          {
            oss << "Coordinates DataArrayDouble of rank #" << i << " differ :";
            reason.insert(0,oss.str());
            return false;
          }
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
  if(!isEqualWithoutConsideringStr(other,prec))
    throw INTERP_KERNEL::Exception("MEDCouplingCMesh::checkDeepEquivalWith : Meshes are not the same !");
}

/*!
 * Nothing is done here (except to check that the other is a ParaMEDMEM::MEDCouplingCMesh instance too).
 * The user intend that the nodes are the same, so by construction of ParaMEDMEM::MEDCouplingCMesh, 'this' and 'other' are the same !
 */
void MEDCouplingCMesh::checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                       DataArrayInt *&cellCor) const throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingCMesh *otherC=dynamic_cast<const MEDCouplingCMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCouplingCMesh::checkDeepEquivalOnSameNodesWith : other is NOT a cartesian mesh ! Impossible to check equivalence !");
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

void MEDCouplingCMesh::checkCoherency1(double eps) const throw(INTERP_KERNEL::Exception)
{
  checkCoherency();
  if(_x_array)
    _x_array->checkMonotonic(true, eps);
  if(_y_array)
    _y_array->checkMonotonic(true, eps);
  if(_z_array)
    _z_array->checkMonotonic(true, eps);
}

void MEDCouplingCMesh::checkCoherency2(double eps) const throw(INTERP_KERNEL::Exception)
{
  checkCoherency1(eps);
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
      for(int p=0;p<spaceDim-l-1;p++)
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
      for(int p=0;p<spaceDim-l-1;p++)
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

void MEDCouplingCMesh::GetPosFromId(int nodeId, int spaceDim, const int *split, int *res)
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

std::set<INTERP_KERNEL::NormalizedCellType> MEDCouplingCMesh::getAllGeoTypes() const
{
  INTERP_KERNEL::NormalizedCellType ret;
  switch(getMeshDimension())
    {
    case 3:
      ret=INTERP_KERNEL::NORM_HEXA8;
      break;
    case 2:
      ret=INTERP_KERNEL::NORM_QUAD4;
      break;
    case 1:
      ret=INTERP_KERNEL::NORM_SEG2;
      break;
    default:
      throw INTERP_KERNEL::Exception("Unexpected dimension for MEDCouplingCMesh::getAllGeoTypes !");
    }
  std::set<INTERP_KERNEL::NormalizedCellType> ret2;
  ret2.insert(ret);
  return ret2;
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
  int spaceDim=getSpaceDimension();
  int tmpCell[3],tmpNode[3];
  getSplitCellValues(tmpCell);
  getSplitNodeValues(tmpNode);
  int tmp2[3];
  GetPosFromId(cellId,spaceDim,tmpCell,tmp2);
  switch(spaceDim)
    {
    case 1:
      conn.push_back(tmp2[0]); conn.push_back(tmp2[0]+1);
      break;
    case 2:
      conn.push_back(tmp2[1]*tmpCell[1]+tmp2[0]); conn.push_back(tmp2[1]*tmpCell[1]+tmp2[0]+1);
      conn.push_back((tmp2[1]+1)*(tmpCell[1]+1)+tmp2[0]+1); conn.push_back((tmp2[1]+1)*(tmpCell[1]+1)+tmp2[0]);
      break;
    case 3:
      conn.push_back(tmp2[1]*tmpCell[1]+tmp2[0]+tmp2[2]*tmpNode[2]); conn.push_back(tmp2[1]*tmpCell[1]+tmp2[0]+1+tmp2[2]*tmpNode[2]);
      conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+1+tmp2[2]*tmpNode[2]); conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+tmp2[2]*tmpNode[2]);
      conn.push_back(tmp2[1]*tmpCell[1]+tmp2[0]+(tmp2[2]+1)*tmpNode[2]); conn.push_back(tmp2[1]*tmpCell[1]+tmp2[0]+1+(tmp2[2]+1)*tmpNode[2]);
      conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+1+(tmp2[2]+1)*tmpNode[2]); conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+(tmp2[2]+1)*tmpNode[2]);
      break;
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingCMesh::getNodeIdsOfCell : big problem spacedim must be in 1,2 or 3 !");
    };
}

void MEDCouplingCMesh::getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const throw(INTERP_KERNEL::Exception)
{
  int tmp[3];
  int spaceDim=getSpaceDimension();
  getSplitNodeValues(tmp);
  const DataArrayDouble *tabs[3]={getCoordsAt(0),getCoordsAt(1),getCoordsAt(2)};
  int tmp2[3];
  GetPosFromId(nodeId,spaceDim,tmp,tmp2);
  for(int j=0;j<spaceDim;j++)
    if(tabs[j])
      coo.push_back(tabs[j]->getConstPointer()[tmp2[j]]);
}

std::string MEDCouplingCMesh::simpleRepr() const
{
  std::ostringstream ret;
  ret << "Cartesian mesh with name : \"" << getName() << "\"\n";
  ret << "Description of mesh : \"" << getDescription() << "\"\n";
  int tmpp1,tmpp2;
  double tt=getTime(tmpp1,tmpp2);
  ret << "Time attached to the mesh [unit] : " << tt << " [" << getTimeUnit() << "]\n";
  ret << "Iteration : " << tmpp1  << " Order : " << tmpp2 << "\n";
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

const DataArrayDouble *MEDCouplingCMesh::getCoordsAt(int i) const throw(INTERP_KERNEL::Exception)
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

DataArrayDouble *MEDCouplingCMesh::getCoordsAt(int i) throw(INTERP_KERNEL::Exception)
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

void MEDCouplingCMesh::setCoordsAt(int i, const DataArrayDouble *arr) throw(INTERP_KERNEL::Exception)
{
  DataArrayDouble **thisArr[3]={&_x_array,&_y_array,&_z_array};
  if(i<0 || i>2)
    throw INTERP_KERNEL::Exception("Invalid rank specified must be 0 or 1 or 2.");
  if(arr!=*(thisArr[i]))
    {
      if(*(thisArr[i]))
        (*(thisArr[i]))->decrRef();
      (*(thisArr[i]))=const_cast<DataArrayDouble *>(arr);
      if(*(thisArr[i]))
        (*(thisArr[i]))->incrRef();
      declareAsNew();
    }
}

void MEDCouplingCMesh::setCoords(const DataArrayDouble *coordsX, const DataArrayDouble *coordsY, const DataArrayDouble *coordsZ)
{
  if(_x_array)
    _x_array->decrRef();
  _x_array=const_cast<DataArrayDouble *>(coordsX);
  if(_x_array)
    _x_array->incrRef();
  if(_y_array)
    _y_array->decrRef();
  _y_array=const_cast<DataArrayDouble *>(coordsY);
  if(_y_array)
    _y_array->incrRef();
  if(_z_array)
    _z_array->decrRef();
  _z_array=const_cast<DataArrayDouble *>(coordsZ);
  if(_z_array)
    _z_array->incrRef();
  declareAsNew();
}

/*!
 * See MEDCouplingUMesh::getDistributionOfTypes for more information
 */
std::vector<int> MEDCouplingCMesh::getDistributionOfTypes() const throw(INTERP_KERNEL::Exception)
{
  //only one type of cell
  std::vector<int> ret(3);
  ret[0]=getTypeOfCell(0);
  ret[1]=getNumberOfCells();
  ret[2]=0; //ret[3*k+2]==0 because it has no sense here
  return ret;
}

/*!
 * See MEDCouplingUMesh::checkTypeConsistencyAndContig for more information
 */
DataArrayInt *MEDCouplingCMesh::checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception)
{
  if(code.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingCMesh::checkTypeConsistencyAndContig : code is empty, should not !");
  std::size_t sz=code.size();
  if(sz!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingCMesh::checkTypeConsistencyAndContig : code should be of size 3 exactly !");

  int nbCells=getNumberOfCellsWithType((INTERP_KERNEL::NormalizedCellType)code[0]);
  if(code[2]==-1)
    {
      if(code[1]==nbCells)
        return 0;
      else
        throw INTERP_KERNEL::Exception("MEDCouplingCMesh::checkTypeConsistencyAndContig : number of cells mismatch !");
    }
  else
    {
      if(code[2]<-1) 
        throw INTERP_KERNEL::Exception("MEDCouplingCMesh::checkTypeConsistencyAndContig : code[2]<-1 mismatch !");
      if(code[2]>=(int)idsPerType.size()) 
        throw INTERP_KERNEL::Exception("MEDCouplingCMesh::checkTypeConsistencyAndContig : code[2]>size idsPerType !");
      return idsPerType[code[2]]->deepCpy();
    }
}

/*!
 * See MEDCouplingUMesh::splitProfilePerType for more information
 */
void MEDCouplingCMesh::splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception)
{
  int nbCells=getNumberOfCells();
  code.resize(3);
  code[0]=(int)getTypeOfCell(0);
  code[1]=nbCells;
  code[2]=0;
  idsInPflPerType.push_back(profile->deepCpy());
  idsPerType.push_back(profile->deepCpy());
}

MEDCouplingUMesh *MEDCouplingCMesh::buildUnstructured() const throw(INTERP_KERNEL::Exception)
{
  int spaceDim=getSpaceDimension();
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New(getName(),spaceDim);
  DataArrayDouble *coords=getCoordinatesAndOwner();
  ret->setCoords(coords);
  coords->decrRef();
  switch(spaceDim)
    {
    case 1:
      fill1DUnstructuredMesh(ret);
      break;
    case 2:
      fill2DUnstructuredMesh(ret);
      break;
    case 3:
      fill3DUnstructuredMesh(ret);
      break;
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingCMesh::buildUnstructured : big problem spacedim must be in 1,2 or 3 !");
    };
  return ret;
}

MEDCouplingMesh *MEDCouplingCMesh::buildPart(const int *start, const int *end) const
{
  MEDCouplingUMesh *um=buildUnstructured();
  MEDCouplingMesh *ret=um->buildPart(start,end);
  um->decrRef();
  return ret;
}

MEDCouplingMesh *MEDCouplingCMesh::buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const
{
  MEDCouplingUMesh *um=buildUnstructured();
  MEDCouplingMesh *ret=um->buildPartAndReduceNodes(start,end,arr);
  um->decrRef();
  return ret;
}

DataArrayInt *MEDCouplingCMesh::simplexize(int policy) throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("MEDCouplingCMesh::simplexize : not available for Cartesian mesh !");
}

void MEDCouplingCMesh::getBoundingBox(double *bbox) const
{
  int dim=getSpaceDimension();
  int j=0;
  for (int idim=0; idim<dim; idim++)
    {
      const DataArrayDouble *c=getCoordsAt(idim);
      if(c)
        {
          const double *coords=c->getConstPointer();
          int nb=c->getNbOfElems();
          bbox[2*j]=coords[0];
          bbox[2*j+1]=coords[nb-1];
          j++;
        }
    }
}

MEDCouplingFieldDouble *MEDCouplingCMesh::getMeasureField(bool isAbs) const
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
  field->setMesh(const_cast<MEDCouplingCMesh *>(this));
  int tmp[3];
  getSplitCellValues(tmp);
  int dim=getSpaceDimension();
  const double **thisArr=new const double *[dim];
  const DataArrayDouble *thisArr2[3]={_x_array,_y_array,_z_array};
  for(int i=0;i<dim;i++)
    thisArr[i]=thisArr2[i]->getConstPointer();
  for(int icell=0;icell<nbelem;icell++)
    {
      int tmp2[3];
      GetPosFromId(icell,dim,tmp,tmp2);
      area_vol[icell]=1.;
      for(int i=0;i<dim;i++)
        area_vol[icell]*=thisArr[i][tmp2[i]+1]-thisArr[i][tmp2[i]];
    }
  delete [] thisArr;
  return field;
}

/*!
 * not implemented yet !
 */
MEDCouplingFieldDouble *MEDCouplingCMesh::getMeasureFieldOnNode(bool isAbs) const
{
  throw INTERP_KERNEL::Exception("MEDCouplingCMesh::getMeasureFieldOnNode : not implemented yet !");
  //return 0;
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
    { vals[3*i]=0.; vals[3*i+1]=0.; vals[3*i+2]=1.; }
  ret->setArray(array);
  array->decrRef();
  ret->setMesh(this);
  return ret;
}

int MEDCouplingCMesh::getCellContainingPoint(const double *pos, double eps) const
{
  int dim=getSpaceDimension();
  int ret=0;
  int coeff=1;
  for(int i=0;i<dim;i++)
    {
      const double *d=getCoordsAt(i)->getConstPointer();
      int nbOfNodes=getCoordsAt(i)->getNbOfElems();
      double ref=pos[i];
      const double *w=std::find_if(d,d+nbOfNodes,std::bind2nd(std::greater_equal<double>(),ref));
      int w2=(int)std::distance(d,w);
      if(w2<nbOfNodes)
        {
          if(w2==0)
            {
              if(ref>d[0]-eps)
                w2=1;
              else
                return -1;
            }
          ret+=coeff*(w2-1);
          coeff*=nbOfNodes-1;
        }
      else
        return -1;
    }
  return ret;
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

void MEDCouplingCMesh::scale(const double *point, double factor)
{
  for(int i=0;i<3;i++)
    {
      DataArrayDouble *c=getCoordsAt(i);
      if(c)
        {
          double *coords=c->getPointer();
          int lgth=c->getNbOfElems();
          std::transform(coords,coords+lgth,coords,std::bind2nd(std::minus<double>(),point[i]));
          std::transform(coords,coords+lgth,coords,std::bind2nd(std::multiplies<double>(),factor));
          std::transform(coords,coords+lgth,coords,std::bind2nd(std::plus<double>(),point[i]));
          c->declareAsNew();
        }
    }
  updateTime();
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
  const DataArrayDouble *tabs[3]={getCoordsAt(0),getCoordsAt(1),getCoordsAt(2)};
  const double *tabsPtr[3];
  for(int j=0;j<spaceDim;j++)
    {
      tabsPtr[j]=tabs[j]->getConstPointer();
      ret->setInfoOnComponent(j,tabs[j]->getInfoOnComponent(0).c_str());
    }
  int tmp2[3];
  for(int i=0;i<nbNodes;i++)
    {
      GetPosFromId(i,spaceDim,tmp,tmp2);
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
  const DataArrayDouble *tabs[3]={getCoordsAt(0),getCoordsAt(1),getCoordsAt(2)};
  std::vector<double> tabsPtr[3];
  for(int j=0;j<spaceDim;j++)
    {
      int sz=tabs[j]->getNbOfElems()-1;
      ret->setInfoOnComponent(j,tabs[j]->getInfoOnComponent(0).c_str());
      const double *srcPtr=tabs[j]->getConstPointer();
      tabsPtr[j].insert(tabsPtr[j].end(),srcPtr,srcPtr+sz);
      std::transform(tabsPtr[j].begin(),tabsPtr[j].end(),srcPtr+1,tabsPtr[j].begin(),std::plus<double>());
      std::transform(tabsPtr[j].begin(),tabsPtr[j].end(),tabsPtr[j].begin(),std::bind2nd(std::multiplies<double>(),0.5));
    }
  int tmp2[3];
  for(int i=0;i<nbCells;i++)
    {
      GetPosFromId(i,spaceDim,tmp,tmp2);
      for(int j=0;j<spaceDim;j++)
        pt[i*spaceDim+j]=tabsPtr[j][tmp2[j]];
    }
  return ret;
}

void MEDCouplingCMesh::renumberCells(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Functionnality of renumbering cell not available for CMesh !");
}

void MEDCouplingCMesh::fill1DUnstructuredMesh(MEDCouplingUMesh *m) const
{
  const DataArrayDouble *c=getCoordsAt(0);
  int nbOfCells=c->getNbOfElems()-1;
  DataArrayInt *connI=DataArrayInt::New();
  connI->alloc(nbOfCells+1,1);
  int *ci=connI->getPointer();
  DataArrayInt *conn=DataArrayInt::New();
  conn->alloc(3*nbOfCells,1);
  ci[0]=0;
  int *cp=conn->getPointer();
  for(int i=0;i<nbOfCells;i++)
    {
      cp[3*i]=(int)INTERP_KERNEL::NORM_SEG2;
      cp[3*i+1]=i;
      cp[3*i+2]=i+1;
      ci[i+1]=3*(i+1);
    }
  m->setConnectivity(conn,connI,true);
  conn->decrRef();
  connI->decrRef();
}

void MEDCouplingCMesh::fill2DUnstructuredMesh(MEDCouplingUMesh *m) const
{
  const DataArrayDouble *c1=getCoordsAt(0);
  const DataArrayDouble *c2=getCoordsAt(1);
  int n1=c1->getNbOfElems()-1;
  int n2=c2->getNbOfElems()-1;
  DataArrayInt *connI=DataArrayInt::New();
  connI->alloc(n1*n2+1,1);
  int *ci=connI->getPointer();
  DataArrayInt *conn=DataArrayInt::New();
  conn->alloc(5*n1*n2,1);
  ci[0]=0;
  int *cp=conn->getPointer();
  int pos=0;
  for(int j=0;j<n2;j++)
    for(int i=0;i<n1;i++,pos++)
      {
        cp[5*pos]=(int)INTERP_KERNEL::NORM_QUAD4;
        cp[5*pos+1]=i+1+j*(n1+1);
        cp[5*pos+2]=i+j*(n1+1);
        cp[5*pos+3]=i+(j+1)*(n1+1);
        cp[5*pos+4]=i+1+(j+1)*(n1+1);
        ci[pos+1]=5*(pos+1);
    }
  m->setConnectivity(conn,connI,true);
  conn->decrRef();
  connI->decrRef();
}

void MEDCouplingCMesh::fill3DUnstructuredMesh(MEDCouplingUMesh *m) const
{
  const DataArrayDouble *c1=getCoordsAt(0);
  const DataArrayDouble *c2=getCoordsAt(1);
  const DataArrayDouble *c3=getCoordsAt(2);
  int n1=c1->getNbOfElems()-1;
  int n2=c2->getNbOfElems()-1;
  int n3=c3->getNbOfElems()-1;
  DataArrayInt *connI=DataArrayInt::New();
  connI->alloc(n1*n2*n3+1,1);
  int *ci=connI->getPointer();
  DataArrayInt *conn=DataArrayInt::New();
  conn->alloc(9*n1*n2*n3,1);
  ci[0]=0;
  int *cp=conn->getPointer();
  int pos=0;
  for(int k=0;k<n3;k++)
    for(int j=0;j<n2;j++)
      for(int i=0;i<n1;i++,pos++)
        {
          cp[9*pos]=(int)INTERP_KERNEL::NORM_HEXA8;
          int tmp=(n1+1)*(n2+1);
          cp[9*pos+1]=i+1+j*(n1+1)+k*tmp;
          cp[9*pos+2]=i+j*(n1+1)+k*tmp;
          cp[9*pos+3]=i+(j+1)*(n1+1)+k*tmp;
          cp[9*pos+4]=i+1+(j+1)*(n1+1)+k*tmp;
          cp[9*pos+5]=i+1+j*(n1+1)+(k+1)*tmp;
          cp[9*pos+6]=i+j*(n1+1)+(k+1)*tmp;
          cp[9*pos+7]=i+(j+1)*(n1+1)+(k+1)*tmp;
          cp[9*pos+8]=i+1+(j+1)*(n1+1)+(k+1)*tmp;
          ci[pos+1]=9*(pos+1);
        }
  m->setConnectivity(conn,connI,true);
  conn->decrRef();
  connI->decrRef();
}

void MEDCouplingCMesh::getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const
{
  int it,order;
  double time=getTime(it,order);
  tinyInfo.clear();
  tinyInfoD.clear();
  littleStrings.clear();
  littleStrings.push_back(getName());
  littleStrings.push_back(getDescription());
  littleStrings.push_back(getTimeUnit());
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
  tinyInfo.push_back(it);
  tinyInfo.push_back(order);
  tinyInfoD.push_back(time);
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

void MEDCouplingCMesh::unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                                       const std::vector<std::string>& littleStrings)
{
  setName(littleStrings[0].c_str());
  setDescription(littleStrings[1].c_str());
  setTimeUnit(littleStrings[2].c_str());
  DataArrayDouble **thisArr[3]={&_x_array,&_y_array,&_z_array};
  const double *data=a2->getConstPointer();
  for(int i=0;i<3;i++)
    {
      if(tinyInfo[i]!=-1)
        {
          (*(thisArr[i]))=DataArrayDouble::New();
          (*(thisArr[i]))->alloc(tinyInfo[i],1);
          (*(thisArr[i]))->setInfoOnComponent(0,littleStrings[i+3].c_str());
          std::copy(data,data+tinyInfo[i],(*(thisArr[i]))->getPointer());
          data+=tinyInfo[i];
        }
    }
  setTime(tinyInfoD[0],tinyInfo[3],tinyInfo[4]);
}

void MEDCouplingCMesh::writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("MEDCouplingCMesh::writeVTKLL : not implemented yet !");
}

std::string MEDCouplingCMesh::getVTKDataSetType() const throw(INTERP_KERNEL::Exception)
{
  return std::string("RectilinearGrid");
}
