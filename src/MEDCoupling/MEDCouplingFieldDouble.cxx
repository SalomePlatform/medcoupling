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

#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingPointSet.hxx"
#include "MEDCouplingTimeDiscretization.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include <sstream>
#include <functional>

using namespace ParaMEDMEM;

MEDCouplingFieldDouble *MEDCouplingFieldDouble::New(TypeOfField type, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldDouble(type,td);
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::clone(bool recDeepCpy) const
{
  return new MEDCouplingFieldDouble(*this,recDeepCpy);
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::buildNewTimeReprFromThis(TypeOfTimeDiscretization td, bool deepCpy) const
{
  MEDCouplingTimeDiscretization *tdo=_time_discr->buildNewTimeReprFromThis(_time_discr,td,deepCpy);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(getNature(),tdo,getTypeOfField());
  ret->setMesh(getMesh());
  ret->setName(getName());
  ret->setDescription(getDescription());
  return ret;
}

bool MEDCouplingFieldDouble::isEqual(const MEDCouplingField *other, double meshPrec, double valsPrec) const
{
  const MEDCouplingFieldDouble *otherC=dynamic_cast<const MEDCouplingFieldDouble *>(other);
  if(!otherC)
    return false;
  if(_nature!=otherC->_nature)
    return false;
  if(!MEDCouplingField::isEqual(other,meshPrec,valsPrec))
    return false;
  if(!_time_discr->isEqual(otherC->_time_discr,valsPrec))
    return false;
  return true;
}

/*!
 * This method states if 'this' and 'other' are compatibles each other before performing any treatment.
 * This method is good for methods like : mergeFields.
 * This method is not very demanding compared to areStrictlyCompatible that is better for operation on fields.
 */
bool MEDCouplingFieldDouble::areCompatibleForMerge(const MEDCouplingField *other) const
{
  if(!MEDCouplingField::areCompatibleForMerge(other))
    return false;
  const MEDCouplingFieldDouble *otherC=dynamic_cast<const MEDCouplingFieldDouble *>(other);
  if(!otherC)
    return false;
  if(_nature!=otherC->_nature)
    return false;
  if(!_time_discr->areCompatible(otherC->_time_discr))
    return false;
  return true;
}

/*!
 * This method is more strict than MEDCouplingField::areCompatible method.
 * This method is used for operation on fields to operate a first check before attempting operation.
 */
bool MEDCouplingFieldDouble::areStrictlyCompatible(const MEDCouplingField *other) const
{
  if(!MEDCouplingField::areStrictlyCompatible(other))
    return false;
  const MEDCouplingFieldDouble *otherC=dynamic_cast<const MEDCouplingFieldDouble *>(other);
  if(!otherC)
    return false;
  if(_nature!=otherC->_nature)
    return false;
  if(!_time_discr->areStrictlyCompatible(otherC->_time_discr))
    return false;
  return true;
}

bool MEDCouplingFieldDouble::areCompatibleForMul(const MEDCouplingField *other) const
{
  if(!MEDCouplingField::areStrictlyCompatible(other))
    return false;
  const MEDCouplingFieldDouble *otherC=dynamic_cast<const MEDCouplingFieldDouble *>(other);
  if(!otherC)
    return false;
  if(_nature!=otherC->_nature)
    return false;
  if(!_time_discr->areStrictlyCompatibleForMul(otherC->_time_discr))
    return false;
  return true;
}

/*!
 * This method performs a clone of mesh and a renumbering of underlying cells of it. The number of cells remains the same.
 * The values of field are impacted in consequence to have the same geometrical field.
 */
void MEDCouplingFieldDouble::renumberCells(const int *old2NewBg, const int *old2NewEnd, bool check) throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Expecting a defined mesh to be able to operate a renumbering !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=_mesh->deepCpy();
  m->renumberCells(old2NewBg,old2NewEnd,check);
  //
  _type->renumberCells(old2NewBg,old2NewEnd,check);
  std::vector<DataArrayDouble *> arrays;
  _time_discr->getArrays(arrays);
  _type->renumberArraysForCell(_mesh,arrays,old2NewBg,old2NewEnd,check);
  //
  setMesh(m);
  updateTime();
}

/*!
 * This method performs a clone of mesh and a renumbering of underlying nodes of it. The number of nodes remains not compulsory the same as renumberCells method.
 * The values of field are impacted in consequence to have the same geometrical field.
 */
void MEDCouplingFieldDouble::renumberNodes(const int *old2NewBg, const int *old2NewEnd) throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingPointSet *meshC=dynamic_cast<const MEDCouplingPointSet *>(_mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply renumberNodes on it !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingPointSet> meshC2((MEDCouplingPointSet *)meshC->deepCpy());
  std::vector<DataArrayDouble *> arrays;
  _time_discr->getArrays(arrays);
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    _type->renumberValuesOnNodes(old2NewBg,*iter);
  meshC2->renumberNodes(old2NewBg,*std::max_element(old2NewBg,old2NewEnd)+1);
  setMesh(meshC2);
}

TypeOfTimeDiscretization MEDCouplingFieldDouble::getTimeDiscretization() const
{
  return _time_discr->getEnum();
}

MEDCouplingFieldDouble::MEDCouplingFieldDouble(TypeOfField type, TypeOfTimeDiscretization td):MEDCouplingField(type),_nature(NoNature),
                                                                                              _time_discr(MEDCouplingTimeDiscretization::New(td))
{
}

MEDCouplingFieldDouble::MEDCouplingFieldDouble(const MEDCouplingFieldDouble& other, bool deepCpy):MEDCouplingField(other),_nature(other._nature),
                                                                                                  _time_discr(other._time_discr->performCpy(deepCpy))
{
}

MEDCouplingFieldDouble::MEDCouplingFieldDouble(NatureOfField n, MEDCouplingTimeDiscretization *td, TypeOfField type):MEDCouplingField(type),
                                                                                                                     _nature(n),_time_discr(td)
{
}

MEDCouplingFieldDouble::~MEDCouplingFieldDouble()
{
  delete _time_discr;
}

void MEDCouplingFieldDouble::checkCoherency() const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Field invalid because no mesh specified !");
  _time_discr->checkCoherency();
  _type->checkCoherencyBetween(_mesh,getArray());
}

/*!
 * Returns the accumulation (the sum) of comId_th component of each tuples of default array.
 */
double MEDCouplingFieldDouble::accumulate(int compId) const
{
  const double *ptr=getArray()->getConstPointer();
  int nbTuple=getArray()->getNumberOfTuples();
  int nbComps=getArray()->getNumberOfComponents();
  if(compId>=nbComps)
    throw INTERP_KERNEL::Exception("Invalid compId specified : No such nb of components !");
  double ret=0.;
  for(int i=0;i<nbTuple;i++)
    ret+=ptr[i*nbComps+compId];
  return ret;
}

/*!
 * Returns the accumulation (the sum) of all tuples of default array.
 * The res is expected to be of size getNumberOfComponents().
 */
void MEDCouplingFieldDouble::accumulate(double *res) const
{
  const double *ptr=getArray()->getConstPointer();
  int nbTuple=getArray()->getNumberOfTuples();
  int nbComps=getArray()->getNumberOfComponents();
  std::fill(res,res+nbComps,0.);
  for(int i=0;i<nbTuple;i++)
    std::transform(ptr+i*nbComps,ptr+(i+1)*nbComps,res,res,std::plus<double>());
}

/*!
 * Returns the normL1 of current field on compId component :
 * \f[
 * \frac{\sum_{0 \leq i < nbOfEntity}|val[i]*Vol[i]|}{\sum_{0 \leq i < nbOfEntity}|Vol[i]|}
 * \f]
 * If compId>=nbOfComponent an exception is thrown.
 */
double MEDCouplingFieldDouble::normL1(int compId, bool isWAbs) const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform normL1");
  int nbComps=getArray()->getNumberOfComponents();
  if(compId>=nbComps)
    throw INTERP_KERNEL::Exception("Invalid compId specified : No such nb of components !");
  double *res=new double[nbComps];
  try
    {
      _type->normL1(_mesh,getArray(),isWAbs,res);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      delete [] res;
      throw e;
    }
  double ret=res[compId];
  delete [] res;
  return ret;
}

/*!
 * Returns the normL1 of current field on each components :
 * \f[
 * \frac{\sum_{0 \leq i < nbOfEntity}|val[i]*Vol[i]|}{\sum_{0 \leq i < nbOfEntity}|Vol[i]|}
 * \f]
 * The res is expected to be of size getNumberOfComponents().
 */
void MEDCouplingFieldDouble::normL1(bool isWAbs, double *res) const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform normL1");
  _type->normL1(_mesh,getArray(),isWAbs,res);
}

/*!
 * Returns the normL2 of current field on compId component :
 * \f[
 * \sqrt{\frac{\sum_{0 \leq i < nbOfEntity}|val[i]^{2}*Vol[i]|}{\sum_{0 \leq i < nbOfEntity}|Vol[i]|}}
 * \f]
 * If compId>=nbOfComponent an exception is thrown.
 */
double MEDCouplingFieldDouble::normL2(int compId, bool isWAbs) const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform normL1");
  int nbComps=getArray()->getNumberOfComponents();
  if(compId>=nbComps)
    throw INTERP_KERNEL::Exception("Invalid compId specified : No such nb of components !");
  double *res=new double[nbComps];
  try
    {
      _type->normL2(_mesh,getArray(),isWAbs,res);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      delete [] res;
      throw e;
    }
  double ret=res[compId];
  delete [] res;
  return ret;
}

/*!
 * Returns the normL2 of current field on each components :
 * \f[
 * \sqrt{\frac{\sum_{0 \leq i < nbOfEntity}|val[i]^{2}*Vol[i]|}{\sum_{0 \leq i < nbOfEntity}|Vol[i]|}}
 * \f]
 * The res is expected to be of size getNumberOfComponents().
 */
void MEDCouplingFieldDouble::normL2(bool isWAbs, double *res) const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform normL1");
  _type->normL2(_mesh,getArray(),isWAbs,res);
}

/*!
 * Returns the accumulation (the sum) of comId_th component of each tuples weigthed by the field
 * returns by getWeightingField relative of the _type of field of default array.
 * This method is usefull to check the conservativity of interpolation method.
 */
double MEDCouplingFieldDouble::integral(int compId, bool isWAbs) const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform integral");
  int nbComps=getArray()->getNumberOfComponents();
  if(compId>=nbComps)
    throw INTERP_KERNEL::Exception("Invalid compId specified : No such nb of components !");
  double *res=new double[nbComps];
  try
    {
      _type->integral(_mesh,getArray(),isWAbs,res);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      delete [] res;
      throw e;
    }
  double ret=res[compId];
  delete [] res;
  return ret;
}

/*!
 * Returns the accumulation (the sum) of each tuples weigthed by the field
 * returns by getWeightingField relative of the _type of field of default array.
 * This method is usefull to check the conservativity of interpolation method.
 */
void MEDCouplingFieldDouble::integral(bool isWAbs, double *res) const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform integral2");
  _type->integral(_mesh,getArray(),isWAbs,res);
}

/*!
 * This method is reserved for field lying on structured mesh spatial support. It returns the value of cell localized by (i,j,k)
 * If spatial support is not structured mesh an exception will be thrown.
 * @param res out array expected to be equal to size getNumberOfComponents()
 */
void MEDCouplingFieldDouble::getValueOnPos(int i, int j, int k, double *res) const throw(INTERP_KERNEL::Exception)
{
  const DataArrayDouble *arr=_time_discr->getArray();
  _type->getValueOnPos(arr,_mesh,i,j,k,res);
}

/*!
 * Returns value of 'this' on default time of point 'spaceLoc' using spatial discretization.
 * If 'point' is outside the spatial discretization of this an exception will be thrown.
 */
void MEDCouplingFieldDouble::getValueOn(const double *spaceLoc, double *res) const throw(INTERP_KERNEL::Exception)
{
  const DataArrayDouble *arr=_time_discr->getArray();
  _type->getValueOn(arr,_mesh,spaceLoc,res);
}

/*!
 * Returns value of 'this' on time 'time' of point 'spaceLoc' using spatial discretization.
 * If 'time' is not covered by this->_time_discr an exception will be thrown.
 * If 'point' is outside the spatial discretization of this an exception will be thrown.
 */
void MEDCouplingFieldDouble::getValueOn(const double *spaceLoc, double time, double *res) const throw(INTERP_KERNEL::Exception)
{
  std::vector< const DataArrayDouble *> arrs=_time_discr->getArraysForTime(time);
  std::vector<double> res2;
  for(std::vector< const DataArrayDouble *>::const_iterator iter=arrs.begin();iter!=arrs.end();iter++)
    {
      int sz=res2.size();
      res2.resize(sz+(*iter)->getNumberOfComponents());
      _type->getValueOn(*iter,_mesh,spaceLoc,&res2[sz]);
    }
  _time_discr->getValueForTime(time,res2,res);
}

/*!
 * Applies a*x+b on 'compoId'th component of each cell.
 */
void MEDCouplingFieldDouble::applyLin(double a, double b, int compoId)
{
  _time_discr->applyLin(a,b,compoId);
}

/*!
 * Applyies the function specified by pointer 'func' on each tuples on all arrays contained in _time_discr.
 * If '*func' returns false during one evaluation an exception will be thrown.
 */
void MEDCouplingFieldDouble::applyFunc(int nbOfComp, FunctionToEvaluate func)
{
  _time_discr->applyFunc(nbOfComp,func);
}

/*!
 * Applyies the function specified by the string repr 'func' on each tuples on all arrays contained in _time_discr.
 * If '*func' fails in evaluation during one evaluation an exception will be thrown.
 * The field will contain 'nbOfComp' components after the call.
 */
void MEDCouplingFieldDouble::applyFunc(int nbOfComp, const char *func)
{
  _time_discr->applyFunc(nbOfComp,func);
}

/*!
 * Applyies the function specified by the string repr 'func' on each tuples on all arrays contained in _time_discr.
 * If '*func' fails in evaluation during one evaluation an exception will be thrown.
 * The field will contain exactly the same number of components after the call.
 */
void MEDCouplingFieldDouble::applyFunc(const char *func)
{
  _time_discr->applyFunc(func);
}

int MEDCouplingFieldDouble::getNumberOfComponents() const
{
  return getArray()->getNumberOfComponents();
}

int MEDCouplingFieldDouble::getNumberOfTuples() const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Impossible to retrieve number of tuples because no mesh specified !");
  return _type->getNumberOfTuples(_mesh);
}

void MEDCouplingFieldDouble::updateTime()
{
  MEDCouplingField::updateTime();
  updateTimeWith(*_time_discr);
}

void MEDCouplingFieldDouble::setNature(NatureOfField nat) throw(INTERP_KERNEL::Exception)
{
  _type->checkCompatibilityWithNature(nat);
  _nature=nat;
}

double MEDCouplingFieldDouble::getIJK(int cellId, int nodeIdInCell, int compoId) const
{
  return _type->getIJK(_mesh,getArray(),cellId,nodeIdInCell,compoId);
}

void MEDCouplingFieldDouble::setArray(DataArrayDouble *array)
{
  _time_discr->setArray(array,this);
}

void MEDCouplingFieldDouble::setEndArray(DataArrayDouble *array)
{
  _time_discr->setEndArray(array,this);
}

void MEDCouplingFieldDouble::getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const
{
  tinyInfo.clear();
  _time_discr->getTinySerializationStrInformation(tinyInfo);
  tinyInfo.push_back(_name);
  tinyInfo.push_back(_desc);
}

/*!
 * This method retrieves some critical values to resize and prepare remote instance.
 * The first two elements returned in tinyInfo correspond to the parameters to give in constructor.
 * @param tinyInfo out parameter resized correctly after the call. The length of this vector is tiny.
 */
void MEDCouplingFieldDouble::getTinySerializationIntInformation(std::vector<int>& tinyInfo) const
{
  tinyInfo.clear();
  tinyInfo.push_back((int)_type->getEnum());
  tinyInfo.push_back((int)_time_discr->getEnum());
  tinyInfo.push_back((int)_nature);
  _time_discr->getTinySerializationIntInformation(tinyInfo);
  std::vector<int> tinyInfo2;
  _type->getTinySerializationIntInformation(tinyInfo2);
  tinyInfo.insert(tinyInfo.end(),tinyInfo2.begin(),tinyInfo2.end());
  tinyInfo.push_back(tinyInfo2.size());
}

/*!
 * This method retrieves some critical values to resize and prepare remote instance.
 * @param tinyInfo out parameter resized correctly after the call. The length of this vector is tiny.
 */
void MEDCouplingFieldDouble::getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const
{
  tinyInfo.clear();
  _time_discr->getTinySerializationDbleInformation(tinyInfo);
  std::vector<double> tinyInfo2;
  _type->getTinySerializationDbleInformation(tinyInfo2);
  tinyInfo.insert(tinyInfo.end(),tinyInfo2.begin(),tinyInfo2.end());
  tinyInfo.push_back(tinyInfo2.size());
}

/*!
 * This method has to be called to the new instance filled by CORBA, MPI, File...
 * @param tinyInfoI is the value retrieves from distant result of getTinySerializationIntInformation on source instance to be copied.
 * @param dataInt out parameter. If not null the pointer is already owned by 'this' after the call of this method. In this case no decrRef must be applied.
 * @param arrays out parameter is a vector resized to the right size. The pointers in the vector is already owned by 'this' after the call of this method.
 *               No decrRef must be applied to every instances in returned vector.
 */
void MEDCouplingFieldDouble::resizeForUnserialization(const std::vector<int>& tinyInfoI, DataArrayInt *&dataInt, std::vector<DataArrayDouble *>& arrays)
{
  dataInt=0;
  std::vector<int> tinyInfoITmp(tinyInfoI);
  int sz=tinyInfoITmp.back();
  tinyInfoITmp.pop_back();
  std::vector<int> tinyInfoITmp2(tinyInfoITmp.begin(),tinyInfoITmp.end()-sz);
  std::vector<int> tinyInfoI2(tinyInfoITmp2.begin()+3,tinyInfoITmp2.end());
  _time_discr->resizeForUnserialization(tinyInfoI2,arrays);
  std::vector<int> tinyInfoITmp3(tinyInfoITmp.end()-sz,tinyInfoITmp.end());
  _type->resizeForUnserialization(tinyInfoITmp3,dataInt);
}

void MEDCouplingFieldDouble::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS)
{
  std::vector<int> tinyInfoI2(tinyInfoI.begin()+3,tinyInfoI.end());
  //
  std::vector<double> tmp(tinyInfoD);
  int sz=tinyInfoD.back();
  tmp.pop_back();
  std::vector<double> tmp1(tmp.begin(),tmp.end()-sz);
  std::vector<double> tmp2(tmp.end()-sz,tmp.end());
  //
  _time_discr->finishUnserialization(tinyInfoI2,tmp1,tinyInfoS);
  _nature=(NatureOfField)tinyInfoI[2];
  _type->finishUnserialization(tmp2);
  int nbOfElemS=tinyInfoS.size();
  _name=tinyInfoS[nbOfElemS-2];
  _desc=tinyInfoS[nbOfElemS-1];
}

/*!
 * Contrary to MEDCouplingPointSet class the returned arrays are \b not the responsabilities of the caller.
 * The values returned must be consulted only in readonly mode.
 */
void MEDCouplingFieldDouble::serialize(DataArrayInt *&dataInt, std::vector<DataArrayDouble *>& arrays) const
{
  _time_discr->getArrays(arrays);
  _type->getSerializationIntArray(dataInt);
}

/*!
 * Merge nodes of underlying mesh. In case of some node will be merged the underlying mesh instance will change.
 */
bool MEDCouplingFieldDouble::mergeNodes(double eps) throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingPointSet *meshC=dynamic_cast<const MEDCouplingPointSet *>(_mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply mergeNodes on it !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingPointSet> meshC2((MEDCouplingPointSet *)meshC->deepCpy());
  bool ret;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=meshC2->mergeNodes(eps,ret);
  if(!ret)//no nodes have been merged.
    return ret;
  std::vector<DataArrayDouble *> arrays;
  _time_discr->getArrays(arrays);
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    _type->renumberValuesOnNodes(arr->getConstPointer(),*iter);
  setMesh(meshC2);
  return true;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::mergeFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1->areCompatibleForMerge(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply mergeFields on them !");
  const MEDCouplingMesh *m1=f1->getMesh();
  const MEDCouplingMesh *m2=f2->getMesh();
  MEDCouplingMesh *m=m1->mergeMyselfWith(m2);
  MEDCouplingTimeDiscretization *td=f1->_time_discr->aggregate(f2->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->getTypeOfField());
  ret->setMesh(m);
  m->decrRef();
  ret->setName(f1->getName());
  ret->setDescription(f1->getDescription());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::addFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1->areStrictlyCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply addFields on them !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->add(f2->_time_discr);
  td->copyTinyAttrFrom(*f1->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->getTypeOfField());
  ret->setMesh(f1->getMesh());
  return ret;
}

const MEDCouplingFieldDouble &MEDCouplingFieldDouble::operator+=(const MEDCouplingFieldDouble& other)
{
  if(!areStrictlyCompatible(&other))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply += on them !");
  _time_discr->addEqual(other._time_discr);
  return *this;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::substractFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1->areStrictlyCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply substractFields on them !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->substract(f2->_time_discr);
  td->copyTinyAttrFrom(*f1->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->getTypeOfField());
  ret->setMesh(f1->getMesh());
  return ret;
}

const MEDCouplingFieldDouble &MEDCouplingFieldDouble::operator-=(const MEDCouplingFieldDouble& other)
{
  if(!areStrictlyCompatible(&other))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply -= on them !");
  _time_discr->substractEqual(other._time_discr);
  return *this;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::multiplyFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1->areCompatibleForMul(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply multiplyFields on them !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->multiply(f2->_time_discr);
  td->copyTinyAttrFrom(*f1->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->getTypeOfField());
  ret->setMesh(f1->getMesh());
  return ret;
}

const MEDCouplingFieldDouble &MEDCouplingFieldDouble::operator*=(const MEDCouplingFieldDouble& other)
{
  if(!areCompatibleForMul(&other))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply *= on them !");
  _time_discr->multiplyEqual(other._time_discr);
  return *this;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::divideFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1->areStrictlyCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply divideFields on them !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->divide(f2->_time_discr);
  td->copyTinyAttrFrom(*f1->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->getTypeOfField());
  ret->setMesh(f1->getMesh());
  return ret;
}

const MEDCouplingFieldDouble &MEDCouplingFieldDouble::operator/=(const MEDCouplingFieldDouble& other)
{
  if(!areStrictlyCompatible(&other))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply /= on them !");
  _time_discr->divideEqual(other._time_discr);
  return *this;
}
