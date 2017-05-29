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
// Author : Anthony Geay (CEA/DEN)

#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingCurveLinearMesh.hxx"

#include "InterpKernelAutoPtr.hxx"

#include <functional>
#include <algorithm>
#include <sstream>
#include <numeric>

using namespace MEDCoupling;

MEDCouplingCMesh::MEDCouplingCMesh():_x_array(0),_y_array(0),_z_array(0)
{
}

MEDCouplingCMesh::MEDCouplingCMesh(const MEDCouplingCMesh& other, bool deepCpy):MEDCouplingStructuredMesh(other,deepCpy)
{
  if(deepCpy)
    {
      if(other._x_array)
        _x_array=other._x_array->deepCopy();
      else
        _x_array=0;
      if(other._y_array)
        _y_array=other._y_array->deepCopy();
      else
        _y_array=0;
      if(other._z_array)
        _z_array=other._z_array->deepCopy();
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

MEDCouplingCMesh *MEDCouplingCMesh::New(const std::string& meshName)
{
  MEDCouplingCMesh *ret(new MEDCouplingCMesh);
  ret->setName(meshName);
  return ret;
}

MEDCouplingCMesh *MEDCouplingCMesh::deepCopy() const
{
  return clone(true);
}

MEDCouplingCMesh *MEDCouplingCMesh::clone(bool recDeepCpy) const
{
  return new MEDCouplingCMesh(*this,recDeepCpy);
}

const DataArrayDouble *MEDCouplingCMesh::getDirectAccessOfCoordsArrIfInStructure() const
{
  throw INTERP_KERNEL::Exception("MEDCouplingCMesh::getDirectAccessOfCoordsArrIfInStructure : MEDCouplingCMesh does not aggregate array of coordinates !");
}

MEDCouplingCurveLinearMesh *MEDCouplingCMesh::buildCurveLinear() const
{
  checkConsistencyLight();
  int dim(getSpaceDimension());
  MCAuto<MEDCouplingCurveLinearMesh> ret(MEDCouplingCurveLinearMesh::New());
  ret->MEDCouplingStructuredMesh::operator=(*this);
  INTERP_KERNEL::AutoPtr<int> ngs(new int[dim]);
  getNodeGridStructure(ngs);
  ret->setNodeGridStructure(ngs,ngs+dim);
  MCAuto<DataArrayDouble> coo(getCoordinatesAndOwner());
  ret->setCoords(coo);
  return ret.retn();
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

std::size_t MEDCouplingCMesh::getHeapMemorySizeWithoutChildren() const
{
  return MEDCouplingStructuredMesh::getHeapMemorySizeWithoutChildren();
}

std::vector<const BigMemoryObject *> MEDCouplingCMesh::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  ret.push_back(_x_array);
  ret.push_back(_y_array);
  ret.push_back(_z_array);
  return ret;
}

/*!
 * This method copyies all tiny strings from other (name and components name).
 * @throw if other and this have not same mesh type.
 */
void MEDCouplingCMesh::copyTinyStringsFrom(const MEDCouplingMesh *other)
{
  MEDCouplingStructuredMesh::copyTinyStringsFrom(other);
  const MEDCouplingCMesh *otherC(dynamic_cast<const MEDCouplingCMesh *>(other));
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCouplingCMesh::copyTinyStringsFrom : meshes have not same type !");
  if(_x_array && otherC->_x_array)
    _x_array->copyStringInfoFrom(*otherC->_x_array);
  if(_y_array && otherC->_y_array)
    _y_array->copyStringInfoFrom(*otherC->_y_array);
  if(_z_array && otherC->_z_array)
    _z_array->copyStringInfoFrom(*otherC->_z_array);
}

bool MEDCouplingCMesh::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingCMesh::isEqualIfNotWhy : input other pointer is null !");
  const MEDCouplingCMesh *otherC=dynamic_cast<const MEDCouplingCMesh *>(other);
  if(!otherC)
    {
      reason="mesh given in input is not castable in MEDCouplingCMesh !";
      return false;
    }
  if(!MEDCouplingStructuredMesh::isEqualIfNotWhy(other,prec,reason))
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
                                            DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const
{
  if(!isEqualWithoutConsideringStr(other,prec))
    throw INTERP_KERNEL::Exception("MEDCouplingCMesh::checkDeepEquivalWith : Meshes are not the same !");
}

/*!
 * Nothing is done here (except to check that the other is a MEDCoupling::MEDCouplingCMesh instance too).
 * The user intend that the nodes are the same, so by construction of MEDCoupling::MEDCouplingCMesh, \a this and \a other are the same !
 */
void MEDCouplingCMesh::checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                       DataArrayInt *&cellCor) const
{
  if(!isEqualWithoutConsideringStr(other,prec))
    throw INTERP_KERNEL::Exception("MEDCouplingCMesh::checkDeepEquivalOnSameNodesWith : Meshes are not the same !");
}

void MEDCouplingCMesh::checkConsistencyLight() const
{
  const char msg0[]="Invalid ";
  const char msg1[]=" array ! Must contain more than 1 element.";
  const char msg2[]=" array ! Must be with only one component.";
  getSpaceDimension();// here to check that no holes in arrays !
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

void MEDCouplingCMesh::checkConsistency(double eps) const
{
  checkConsistencyLight();
  if(_x_array)
    _x_array->checkMonotonic(true, eps);
  if(_y_array)
    _y_array->checkMonotonic(true, eps);
  if(_z_array)
    _z_array->checkMonotonic(true, eps);
}

void MEDCouplingCMesh::getNodeGridStructure(int *res) const
{
  std::vector<int> ret(getNodeGridStructure());
  std::copy(ret.begin(),ret.end(),res);
}

std::vector<int> MEDCouplingCMesh::getNodeGridStructure() const
{
  static const char MSG[]="MEDCouplingCMesh::getNodeGridStructure : mesh is invalid ! null vectors (X, Y or Z) must be put contiguously at the end !";
  std::vector<int> ret;
  bool isOK(true);
  if(_x_array)
    {
      if(!_x_array->isAllocated() || _x_array->getNumberOfComponents()!=1)
        throw INTERP_KERNEL::Exception("MEDCouplingCMesh::getNodeGridStructure : X array exits but it is not allocated or with nb of components equal to one !");
      ret.push_back(_x_array->getNumberOfTuples());
    }
  else
    isOK=false;
  if(_y_array)
    {
      if(!_y_array->isAllocated() || _y_array->getNumberOfComponents()!=1)
        throw INTERP_KERNEL::Exception("MEDCouplingCMesh::getNodeGridStructure : Y array exits but it is not allocated or with nb of components equal to one !");
      if(!isOK)
        throw INTERP_KERNEL::Exception(MSG);
      ret.push_back(_y_array->getNumberOfTuples());
    }
  else
    isOK=false;
  if(_z_array)
    {
      if(!_z_array->isAllocated() || _z_array->getNumberOfComponents()!=1)
        throw INTERP_KERNEL::Exception("MEDCouplingCMesh::getNodeGridStructure : Z array exits but it is not allocated or with nb of components equal to one !");
      if(!isOK)
        throw INTERP_KERNEL::Exception(MSG);
      ret.push_back(_z_array->getNumberOfTuples());
    }
  return ret;
}

MEDCouplingStructuredMesh *MEDCouplingCMesh::buildStructuredSubPart(const std::vector< std::pair<int,int> >& cellPart) const
{
  checkConsistencyLight();
  int dim(getSpaceDimension());
  if(dim!=(int)cellPart.size())
    {
      std::ostringstream oss; oss << "MEDCouplingCMesh::buildStructuredSubPart : the space dimension is " << dim << " and cell part size is " << cellPart.size() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MCAuto<MEDCouplingCMesh> ret(dynamic_cast<MEDCouplingCMesh *>(deepCopy()));
  for(int i=0;i<dim;i++)
    {
      MCAuto<DataArrayDouble> tmp(ret->getCoordsAt(i)->selectByTupleIdSafeSlice(cellPart[i].first,cellPart[i].second+1,1));
      ret->setCoordsAt(i,tmp);
    }
  return ret.retn();
}

/*!
 * Return the space dimension of \a this. It only considers the arrays along X, Y and Z to deduce that.
 * This method throws exceptions if the not null arrays defining this are not contiguously at the end. For example X!=0,Y==0,Z!=0 will throw.
 */
int MEDCouplingCMesh::getSpaceDimension() const
{
  return (int)getNodeGridStructure().size();
}

void MEDCouplingCMesh::getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const
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
  ret << "Space dimension : " << getSpaceDimension() << "\n\nArrays :\n________\n\n";
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

/*!
 * Returns a DataArrayDouble holding positions of nodes along a given axis.
 * For more info on Cartesian meshes, see \ref MEDCouplingCMeshPage.
 *  \param [in] i - an index of axis, a value from [0,1,2].
 *  \return const DataArrayDouble * - a pointer to the data array of node coordinates
 *         referred by \a this mesh.
 *  \throw If \a i is not one of [0,1,2].
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mccmesh_getCoordsAt "Here is a C++ example".<br>
 *  \ref  py_mccmesh_getCoordsAt "Here is a Python example".
 *  \endif
 */
const DataArrayDouble *MEDCouplingCMesh::getCoordsAt(int i) const
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

/*!
 * Returns a DataArrayDouble holding positions of nodes along a given axis.
 * For more info on Cartesian meshes, see \ref MEDCouplingCMeshPage.
 *  \param [in] i - an index of axis, a value from [0,1,2].
 *  \return const DataArrayDouble * - a pointer to the data array of node coordinates
 *         referred by \a this mesh.
 *  \throw If \a i is not one of [0,1,2].
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mccmesh_getCoordsAt "Here is a C++ example".<br>
 *  \ref  py_mccmesh_getCoordsAt "Here is a Python example".
 *  \endif
 */
DataArrayDouble *MEDCouplingCMesh::getCoordsAt(int i)
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

/*!
 * Sets node coordinates along a given axis. For more info on Cartesian meshes, see 
 * \ref MEDCouplingCMeshPage.
 *  \param [in] i - an index of axis, a value in range [0,1,2].
 *  \param [in] arr - DataArrayDouble holding positions of nodes along the i-th
 *         axis. It must be an array of one component.
 *  \throw If \a arr->getNumberOfComponents() != 1.
 *  \throw If \a i is not one of [0,1,2].
 *
 *  \if ENABLE_EXAMPLES
 *  \ref medcouplingcppexamplesCmeshStdBuild1 "Here is a C++ example".<br>
 *  \ref  medcouplingpyexamplesCmeshStdBuild1 "Here is a Python example".
 *  \endif
 */
void MEDCouplingCMesh::setCoordsAt(int i, const DataArrayDouble *arr)
{
  if(arr)
    arr->checkNbOfComps(1,"MEDCouplingCMesh::setCoordsAt");
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

/*!
 * Sets node coordinates along some of the tree axes. This method updates all the
 * three node coordinates arrays at once. For more info on Cartesian meshes, see 
 * \ref MEDCouplingCMeshPage.
 *  \param [in] coordsX - DataArrayDouble holding positions of nodes along the X
 *         axis. It must be an array of one component or \c NULL.
 *  \param [in] coordsY - DataArrayDouble holding positions of nodes along the Y
 *         axis. It must be an array of one component or \c NULL.
 *  \param [in] coordsZ - DataArrayDouble holding positions of nodes along the Z
 *         axis. It must be an array of one component or \c NULL.
 *  \throw If \a coords*->getNumberOfComponents() != 1.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref medcouplingcppexamplesCmeshStdBuild1 "Here is a C++ example".<br>
 *  \ref  medcouplingpyexamplesCmeshStdBuild1 "Here is a Python example".
 *  \endif
 */
void MEDCouplingCMesh::setCoords(const DataArrayDouble *coordsX, const DataArrayDouble *coordsY, const DataArrayDouble *coordsZ)
{
  if(coordsX)
    coordsX->checkNbOfComps(1,"MEDCouplingCMesh::setCoords : coordsX");
  if(coordsY)
    coordsY->checkNbOfComps(1,"MEDCouplingCMesh::setCoords : coordsY");
  if(coordsZ)
    coordsZ->checkNbOfComps(1,"MEDCouplingCMesh::setCoords : coordsZ");
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

/*!
 * Returns a new MEDCouplingFieldDouble containing volumes of cells constituting \a this
 * mesh.<br>
 * For 1D cells, the returned field contains lengths.<br>
 * For 2D cells, the returned field contains areas.<br>
 * For 3D cells, the returned field contains volumes.
 *  \param [in] isAbs - a not used parameter.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble on cells
 *         and one time . The caller is to delete this field using decrRef() as it is no
 *         more needed.
 */
MEDCouplingFieldDouble *MEDCouplingCMesh::getMeasureField(bool isAbs) const
{
  std::string name="MeasureOfMesh_";
  name+=getName();
  int nbelem=getNumberOfCells();
  MEDCouplingFieldDouble *field=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  field->setName(name);
  DataArrayDouble* array=DataArrayDouble::New();
  array->alloc(nbelem,1);
  double *area_vol=array->getPointer();
  field->setArray(array) ;
  array->decrRef();
  field->setMesh(const_cast<MEDCouplingCMesh *>(this));
  field->synchronizeTimeWithMesh();
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

void MEDCouplingCMesh::getCellsContainingPoint(const double *pos, double eps, std::vector<int>& elts) const
{
  int ret(getCellContainingPoint(pos,eps));
  elts.push_back(ret);
}

void MEDCouplingCMesh::rotate(const double *center, const double *vector, double angle)
{
  throw INTERP_KERNEL::Exception("No rotation available on CMesh : Traduce it to untructured mesh to apply it !");
}

/*!
 * Translates all nodes of \a this mesh by a given vector. Actually, it adds each
 * component of the \a vector to all node coordinates of a corresponding axis.
 *  \param [in] vector - the translation vector whose size must be not less than \a
 *         this->getSpaceDimension().
 */
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

/*!
 * Applies scaling transformation to all nodes of \a this mesh.
 *  \param [in] point - coordinates of a scaling center. This array is to be of
 *         size \a this->getSpaceDimension() at least.
 *  \param [in] factor - a scale factor.
 */
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

/*!
 * Returns a new DataArrayDouble holding coordinates of all nodes of \a this mesh.
 *  \return DataArrayDouble * - a new instance of DataArrayDouble, of size \a
 *          this->getNumberOfNodes() tuples per \a this->getSpaceDimension()
 *          components. The caller is to delete this array using decrRef() as it is
 *          no more needed.
 */
DataArrayDouble *MEDCouplingCMesh::getCoordinatesAndOwner() const
{
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  int spaceDim(getSpaceDimension()),nbNodes(getNumberOfNodes());
  ret->alloc(nbNodes,spaceDim);
  double *pt(ret->getPointer());
  int tmp[3];
  getSplitNodeValues(tmp);
  const DataArrayDouble *tabs[3]={getCoordsAt(0),getCoordsAt(1),getCoordsAt(2)};
  const double *tabsPtr[3];
  for(int j=0;j<spaceDim;j++)
    {
      tabsPtr[j]=tabs[j]->getConstPointer();
      ret->setInfoOnComponent(j,tabs[j]->getInfoOnComponent(0));
    }
  int tmp2[3];
  for(int i=0;i<nbNodes;i++)
    {
      GetPosFromId(i,spaceDim,tmp,tmp2);
      for(int j=0;j<spaceDim;j++)
        pt[i*spaceDim+j]=tabsPtr[j][tmp2[j]];
    }
  return ret.retn();
}

/*!
 * Returns a new DataArrayDouble holding barycenters of all cells. The barycenter is
 * computed by averaging coordinates of cell nodes.
 *  \return DataArrayDouble * - a new instance of DataArrayDouble, of size \a
 *          this->getNumberOfCells() tuples per \a this->getSpaceDimension()
 *          components. The caller is to delete this array using decrRef() as it is
 *          no more needed.
 */
DataArrayDouble *MEDCouplingCMesh::computeCellCenterOfMass() const
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
      ret->setInfoOnComponent(j,tabs[j]->getInfoOnComponent(0));
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

DataArrayDouble *MEDCouplingCMesh::computeIsoBarycenterOfNodesPerCell() const
{
  return MEDCouplingCMesh::computeCellCenterOfMass();
}

void MEDCouplingCMesh::renumberCells(const int *old2NewBg, bool check)
{
  throw INTERP_KERNEL::Exception("Functionnality of renumbering cell not available for CMesh !");
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
  setName(littleStrings[0]);
  setDescription(littleStrings[1]);
  setTimeUnit(littleStrings[2]);
  DataArrayDouble **thisArr[3]={&_x_array,&_y_array,&_z_array};
  const double *data=a2->getConstPointer();
  for(int i=0;i<3;i++)
    {
      if(tinyInfo[i]!=-1)
        {
          (*(thisArr[i]))=DataArrayDouble::New();
          (*(thisArr[i]))->alloc(tinyInfo[i],1);
          (*(thisArr[i]))->setInfoOnComponent(0,littleStrings[i+3]);
          std::copy(data,data+tinyInfo[i],(*(thisArr[i]))->getPointer());
          data+=tinyInfo[i];
        }
    }
  setTime(tinyInfoD[0],tinyInfo[3],tinyInfo[4]);
}

void MEDCouplingCMesh::writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const
{
  std::ostringstream extent;
  DataArrayDouble *thisArr[3]={_x_array,_y_array,_z_array};
  for(int i=0;i<3;i++)
    {
      if(thisArr[i])
        { extent << "0 " <<  thisArr[i]->getNumberOfTuples()-1 << " "; }
      else
        { extent << "0 0 "; }
    }
  ofs << "  <" << getVTKDataSetType() << " WholeExtent=\"" << extent.str() << "\">\n";
  ofs << "    <Piece Extent=\"" << extent.str() << "\">\n";
  ofs << "      <PointData>\n" << pointData << std::endl;
  ofs << "      </PointData>\n";
  ofs << "      <CellData>\n" << cellData << std::endl;
  ofs << "      </CellData>\n";
  ofs << "      <Coordinates>\n";
  for(int i=0;i<3;i++)
    {
      if(thisArr[i])
        thisArr[i]->writeVTK(ofs,8,"Array",byteData);
      else
        {
          MCAuto<DataArrayDouble> coo=DataArrayDouble::New(); coo->alloc(1,1);
          coo->setIJ(0,0,0.);
          coo->writeVTK(ofs,8,"Array",byteData);
        }
    }
  ofs << "      </Coordinates>\n";
  ofs << "    </Piece>\n";
  ofs << "  </" << getVTKDataSetType() << ">\n";
}

void MEDCouplingCMesh::reprQuickOverview(std::ostream& stream) const
{
  stream << "MEDCouplingCMesh C++ instance at " << this << ". Name : \"" << getName() << "\".";
  const DataArrayDouble *thisArr[3]={_x_array,_y_array,_z_array};
  std::ostringstream stream2[3];
  bool isDef[3];
  int nbOfCells=1,nbOfNodes=1;
  for(int i=0;i<3;i++)
    {
      isDef[i]=thisArr[i]!=0;
      if(isDef[i])
        {    
          char tmp='X'+i;
          stream2[i] << tmp << " positions array ";
          if(!thisArr[i]->isAllocated())
            stream2[i] << "set but not allocated.";
          else
            {
              int nbCompo=thisArr[i]->getNumberOfComponents();
              if(nbCompo==1)
                {
                  int nbTuples=thisArr[i]->getNumberOfTuples();
                  if(nbTuples<1)
                    { stream2[i] << "set and allocated - WARNING number of elements < 1 !"; nbOfCells=-1; nbOfNodes=-1; }
                  else
                    {
                      stream2[i] << "(length=" << nbTuples << ")" << ": ";
                      thisArr[i]->reprQuickOverviewData(stream2[i],200);
                      if(nbOfCells!=-1)
                        { nbOfNodes*=nbTuples; nbOfCells*=nbTuples-1; }
                    }
                }
              else
                { stream2[i] << "set and allocated - WARNING number of components != 1 !"; nbOfCells=-1; nbOfNodes=-1; }
            }
        }
    }
  if(!isDef[0] && !isDef[1] && !isDef[2])
    { stream << " No arrays set !"; return; }
  if(nbOfCells>=0)
    { stream << std::endl << "Number of cells : " << nbOfCells << ". Number of nodes : " << nbOfNodes << "."; }
  for(int i=0;i<3;i++)
    {
      if(isDef[i])
        stream << std::endl << stream2[i].str();
    }
}

std::string MEDCouplingCMesh::getVTKFileExtension() const
{
  return std::string("vtr");
}

std::string MEDCouplingCMesh::getVTKDataSetType() const
{
  return std::string("RectilinearGrid");
}
