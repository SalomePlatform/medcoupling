//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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

#include <sstream>

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

bool MEDCouplingFieldDouble::areCompatible(const MEDCouplingField *other) const
{
  if(!MEDCouplingField::areCompatible(other))
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
  if(!getArray())
    throw INTERP_KERNEL::Exception("Field invalid because no values set !");
  _type->checkCoherencyBetween(_mesh,getArray());
}

double MEDCouplingFieldDouble::accumulate(int compId) const
{
  const double *ptr=getArray()->getConstPointer();
  int nbTuple=getArray()->getNumberOfTuples();
  int nbComps=getArray()->getNumberOfComponents();
  double ret=0.;
  for(int i=0;i<nbTuple;i++)
    ret+=ptr[i*nbComps+compId];
  return ret;
}

void MEDCouplingFieldDouble::accumulate(double *res) const
{
  const double *ptr=getArray()->getConstPointer();
  int nbTuple=getArray()->getNumberOfTuples();
  int nbComps=getArray()->getNumberOfComponents();
  std::fill(res,res+nbComps,0.);
  for(int i=0;i<nbTuple;i++)
    std::transform(ptr+i*nbComps,ptr+(i+1)*nbComps,res,res,std::plus<double>());
}

double MEDCouplingFieldDouble::measureAccumulate(int compId, bool isWAbs) const
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform measureAccumulate");
  MEDCouplingFieldDouble *weight=_type->getWeightingField(_mesh,isWAbs);
  const double *ptr=weight->getArray()->getConstPointer();
  int nbOfValues=weight->getArray()->getNbOfElems();
  double ret=0.;
  for (int i=0; i<nbOfValues; i++)
    ret+=getIJ(i,compId)*ptr[i];
  weight->decrRef();
  return ret;
}

void MEDCouplingFieldDouble::measureAccumulate(bool isWAbs, double *res) const
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform measureAccumulate2");
  MEDCouplingFieldDouble *weight=_type->getWeightingField(_mesh,isWAbs);
  const double *ptr=weight->getArray()->getConstPointer();
  int nbOfValues=weight->getArray()->getNbOfElems();
  int nbComps=getArray()->getNumberOfComponents();
  const double *vals=getArray()->getConstPointer();
  std::fill(res,res+nbComps,0.);
  double *tmp=new double[nbComps];
  for (int i=0; i<nbOfValues; i++)
    {
      std::transform(vals+i*nbComps,vals+(i+1)*nbComps,tmp,std::bind2nd(std::multiplies<double>(),ptr[i]));
      std::transform(tmp,tmp+nbComps,res,res,std::plus<double>());
    }
  weight->decrRef();
  delete [] tmp;
}

void MEDCouplingFieldDouble::getValueOn(const double *spaceLoc, double *res) const throw(INTERP_KERNEL::Exception)
{
  _time_discr->checkNoTimePresence();
}

void MEDCouplingFieldDouble::getValueOn(const double *spaceLoc, double time, double *res) const throw(INTERP_KERNEL::Exception)
{
  _time_discr->checkTimePresence(time);
}

void MEDCouplingFieldDouble::applyLin(double a, double b, int compoId)
{
  _time_discr->applyLin(a,b,compoId);
}

void MEDCouplingFieldDouble::applyFunc(int nbOfComp, FunctionToEvaluate func)
{
  _time_discr->applyFunc(nbOfComp,func);
}

void MEDCouplingFieldDouble::applyFunc(int nbOfComp, const char *func)
{
  _time_discr->applyFunc(nbOfComp,func);
}

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
  if(getArray())
    updateTimeWith(*getArray());
}

void MEDCouplingFieldDouble::setNature(NatureOfField nat) throw(INTERP_KERNEL::Exception)
{
  _type->checkCompatibilityWithNature(nat);
  _nature=nat;
}

void MEDCouplingFieldDouble::setArray(DataArrayDouble *array)
{
  _time_discr->setArray(array,this);
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
}

/*!
 * This method retrieves some critical values to resize and prepare remote instance.
 * @param tinyInfo out parameter resized correctly after the call. The length of this vector is tiny.
 */
void MEDCouplingFieldDouble::getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const
{
  tinyInfo.clear();
  _time_discr->getTinySerializationDbleInformation(tinyInfo);
}

/*!
 * This method has to be called to the new instance filled by CORBA, MPI, File...
 * @param tinyInfoI is the value retrieves from distant result of getTinySerializationIntInformation on source instance to be copied.
 * @param arrays out parameter is a vector resized to the right size. The pointers in the vector is already owned by 'this' after the call of this method.
 *               No decrRef must be applied to every instances in returned vector.
 */
void MEDCouplingFieldDouble::resizeForUnserialization(const std::vector<int>& tinyInfoI, std::vector<DataArrayDouble *>& arrays)
{
  std::vector<int> tinyInfoI2(tinyInfoI.begin()+3,tinyInfoI.end());
  _time_discr->resizeForUnserialization(tinyInfoI2,arrays);
}

void MEDCouplingFieldDouble::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS)
{
  std::vector<int> tinyInfoI2(tinyInfoI.begin()+3,tinyInfoI.end());
  _time_discr->finishUnserialization(tinyInfoI2,tinyInfoD,tinyInfoS);
  _nature=(NatureOfField)tinyInfoI[2];
  int nbOfElemS=tinyInfoS.size();
  _name=tinyInfoS[nbOfElemS-2];
  _desc=tinyInfoS[nbOfElemS-1];
}

/*!
 * Contrary to MEDCouplingPointSet class the returned arrays are \b not the responsabilities of the caller.
 * The values returned must be consulted only in readonly mode.
 */
void MEDCouplingFieldDouble::serialize(std::vector<DataArrayDouble *>& arrays) const
{
  _time_discr->getArrays(arrays);
}

/*!
 * \b Warning ! This method potentially modifies the underlying mesh ! If the mesh is shared by other fields, these fields could be unavailable.
 */
bool MEDCouplingFieldDouble::mergeNodes(double eps)
{
  MEDCouplingPointSet *meshC=dynamic_cast<MEDCouplingPointSet *>((MEDCouplingMesh *)(_mesh));
  if(!meshC)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply mergeNodes on it !");
  bool ret;
  DataArrayInt *arr=meshC->mergeNodes(eps,ret);
  if(!ret)//no nodes have been merged.
    return ret;
  std::vector<DataArrayDouble *> arrays;
  _time_discr->getArrays(arrays);
  try
    {
      for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
        _type->renumberValuesOnNodes(arr,*iter);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      arr->decrRef();
      throw e;
    }
  arr->decrRef();
  return true;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::mergeFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1->areCompatible(f2))
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
  if(!f1->areCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply addFields on them !");
  if(f1->getMesh()!=f2->getMesh())
    throw INTERP_KERNEL::Exception("Fields are not lying on same mesh ; addFields impossible !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->add(f2->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->getTypeOfField());
  ret->setMesh(f1->getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::substractFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1->areCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply substractFields on them !");
  if(f1->getMesh()!=f2->getMesh())
    throw INTERP_KERNEL::Exception("Fields are not lying on same mesh ; substractFields impossible !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->substract(f2->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->getTypeOfField());
  ret->setMesh(f1->getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::multiplyFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1->areCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to applymultiplyFields  on them !");
  if(f1->getMesh()!=f2->getMesh())
    throw INTERP_KERNEL::Exception("Fields are not lying on same mesh ; multiplyFields impossible !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->multiply(f2->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->getTypeOfField());
  ret->setMesh(f1->getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::divideFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1->areCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply divideFields on them !");
  if(f1->getMesh()!=f2->getMesh())
    throw INTERP_KERNEL::Exception("Fields are not lying on same mesh ; divideFields impossible !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->divide(f2->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->getTypeOfField());
  ret->setMesh(f1->getMesh());
  return ret;
}
