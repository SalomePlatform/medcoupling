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

#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingTimeDiscretization.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "MEDCouplingNatureOfField.hxx"

#include <sstream>
#include <limits>
#include <functional>

using namespace ParaMEDMEM;

MEDCouplingFieldDouble *MEDCouplingFieldDouble::New(TypeOfField type, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldDouble(type,td);
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::New(const MEDCouplingFieldTemplate *ft, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldDouble(ft,td);
}

void MEDCouplingFieldDouble::setTimeUnit(const char *unit)
{
  _time_discr->setTimeUnit(unit);
}

const char *MEDCouplingFieldDouble::getTimeUnit() const
{
  return _time_discr->getTimeUnit();
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::clone(bool recDeepCpy) const
{
  return new MEDCouplingFieldDouble(*this,recDeepCpy);
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::cloneWithMesh(bool recDeepCpy) const
{
  MEDCouplingFieldDouble *ret=clone(recDeepCpy);
  if(_mesh)
    {
      MEDCouplingMesh *mCpy=_mesh->deepCpy();
      ret->setMesh(mCpy);
      mCpy->decrRef();
    }
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::deepCpy() const
{
  return cloneWithMesh(true);
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::buildNewTimeReprFromThis(TypeOfTimeDiscretization td, bool deepCpy) const
{
  MEDCouplingTimeDiscretization *tdo=_time_discr->buildNewTimeReprFromThis(td,deepCpy);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(getNature(),tdo,_type->clone());
  ret->setMesh(getMesh());
  ret->setName(getName());
  ret->setDescription(getDescription());
  return ret;
}

/*!
 * Copy tiny info (component names, name, description) but warning the underlying mesh is not renamed (for safety reason).
 */
void MEDCouplingFieldDouble::copyTinyStringsFrom(const MEDCouplingFieldDouble *other) throw(INTERP_KERNEL::Exception)
{
  if(other)
    {
      setName(other->_name.c_str());
      setDescription(other->_desc.c_str());
      _time_discr->copyTinyStringsFrom(*other->_time_discr);
    }
}

/*!
 * Copy only times, order, iteration from other. The underlying mesh is not impacted by this method.
 */
void MEDCouplingFieldDouble::copyTinyAttrFrom(const MEDCouplingFieldDouble *other) throw(INTERP_KERNEL::Exception)
{
  if(other)
    {
      _time_discr->copyTinyAttrFrom(*other->_time_discr);
    }
}

std::string MEDCouplingFieldDouble::simpleRepr() const
{
  std::ostringstream ret;
  ret << "FieldDouble with name : \"" << getName() << "\"\n";
  ret << "Description of field is : \"" << getDescription() << "\"\n";
  ret << "FieldDouble space discretization is : " << _type->getStringRepr() << "\n";
  ret << "FieldDouble time discretization is : " << _time_discr->getStringRepr() << "\n";
  ret << "FieldDouble nature of field is : " << MEDCouplingNatureOfField::getRepr(_nature) << "\n";
  if(getArray())
    {
      int nbOfCompo=getArray()->getNumberOfComponents();
      ret << "FieldDouble default array has " << nbOfCompo << " components and " << getArray()->getNumberOfTuples() << " tuples.\n";
      ret << "FieldDouble default array has following info on components : ";
      for(int i=0;i<nbOfCompo;i++)
        ret << "\"" << getArray()->getInfoOnComponent(i) << "\" ";
      ret << "\n";
    }
  if(_mesh)
    ret << "Mesh support information :\n__________________________\n" << _mesh->simpleRepr();
  else
    ret << "Mesh support information : No mesh set !\n";
  return ret.str();
}

std::string MEDCouplingFieldDouble::advancedRepr() const
{
  std::ostringstream ret;
  ret << "FieldDouble with name : \"" << getName() << "\"\n";
  ret << "Description of field is : \"" << getDescription() << "\"\n";
  ret << "FieldDouble space discretization is : " << _type->getStringRepr() << "\n";
  ret << "FieldDouble time discretization is : " << _time_discr->getStringRepr() << "\n";
  if(getArray())
    ret << "FieldDouble default array has " << getArray()->getNumberOfComponents() << " components and " << getArray()->getNumberOfTuples() << " tuples.\n";
  if(_mesh)
    ret << "Mesh support information :\n__________________________\n" << _mesh->simpleRepr();
  else
    ret << "Mesh support information : No mesh set !\n";
  std::vector<DataArrayDouble *> arrays;
  _time_discr->getArrays(arrays);
  int arrayId=0;
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++,arrayId++)
    {
      ret << "Array #" << arrayId << " :\n__________\n";
      if(*iter)
        (*iter)->reprWithoutNameStream(ret);
      else
        ret << "Array empty !";
      ret << "\n";
    }
  return ret.str();
}

bool MEDCouplingFieldDouble::isEqual(const MEDCouplingField *other, double meshPrec, double valsPrec) const
{
  const MEDCouplingFieldDouble *otherC=dynamic_cast<const MEDCouplingFieldDouble *>(other);
  if(!otherC)
    return false;
  if(!MEDCouplingField::isEqual(other,meshPrec,valsPrec))
    return false;
  if(!_time_discr->isEqual(otherC->_time_discr,valsPrec))
    return false;
  return true;
}

bool MEDCouplingFieldDouble::isEqualWithoutConsideringStr(const MEDCouplingField *other, double meshPrec, double valsPrec) const
{
  const MEDCouplingFieldDouble *otherC=dynamic_cast<const MEDCouplingFieldDouble *>(other);
  if(!otherC)
    return false;
  if(!MEDCouplingField::isEqualWithoutConsideringStr(other,meshPrec,valsPrec))
    return false;
  if(!_time_discr->isEqualWithoutConsideringStr(otherC->_time_discr,valsPrec))
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
  if(!_time_discr->areCompatible(otherC->_time_discr))
    return false;
  return true;
}

/*!
 * This method is more strict than MEDCouplingField::areCompatibleForMerge method.
 * This method is used for operation on fields to operate a first check before attempting operation.
 */
bool MEDCouplingFieldDouble::areStrictlyCompatible(const MEDCouplingField *other) const
{
  if(!MEDCouplingField::areStrictlyCompatible(other))
    return false;
  const MEDCouplingFieldDouble *otherC=dynamic_cast<const MEDCouplingFieldDouble *>(other);
  if(!otherC)
    return false;
  if(!_time_discr->areStrictlyCompatible(otherC->_time_discr))
    return false;
  return true;
}

/*!
 * Method with same principle than MEDCouplingFieldDouble::areStrictlyCompatible method except that
 * number of components between 'this' and 'other' can be different here (for operator*).
 */
bool MEDCouplingFieldDouble::areCompatibleForMul(const MEDCouplingField *other) const
{
  if(!MEDCouplingField::areStrictlyCompatible(other))
    return false;
  const MEDCouplingFieldDouble *otherC=dynamic_cast<const MEDCouplingFieldDouble *>(other);
  if(!otherC)
    return false;
  if(!_time_discr->areStrictlyCompatibleForMul(otherC->_time_discr))
    return false;
  return true;
}

/*!
 * Method with same principle than MEDCouplingFieldDouble::areStrictlyCompatible method except that
 * number of components between 'this' and 'other' can be different here (for operator/).
 */
bool MEDCouplingFieldDouble::areCompatibleForDiv(const MEDCouplingField *other) const
{
  if(!MEDCouplingField::areStrictlyCompatible(other))
    return false;
  const MEDCouplingFieldDouble *otherC=dynamic_cast<const MEDCouplingFieldDouble *>(other);
  if(!otherC)
    return false;
  if(!_time_discr->areStrictlyCompatibleForDiv(otherC->_time_discr))
    return false;
  return true;
}

/*!
 * This method is invocated before any attempt of melding. This method is very close to areStrictlyCompatible,
 * except that 'this' and other can have different number of components.
 */
bool MEDCouplingFieldDouble::areCompatibleForMeld(const MEDCouplingFieldDouble *other) const
{
  if(!MEDCouplingField::areStrictlyCompatible(other))
    return false;
  if(!_time_discr->areCompatibleForMeld(other->_time_discr))
    return false;
  return true;
}

/*!
 * This method performs a clone of mesh and a renumbering of underlying cells of it. The number of cells remains the same.
 * The values of field are impacted in consequence to have the same geometrical field.
 */
void MEDCouplingFieldDouble::renumberCells(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
  renumberCellsWithoutMesh(old2NewBg,check);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=_mesh->deepCpy();
  m->renumberCells(old2NewBg,check);
  setMesh(m);
  updateTime();
}

/*!
 * \b WARNING : use this method with lot of care !
 * This method performs half job of MEDCouplingFieldDouble::renumberCells. That is to say no permutation of cells is done on underlying mesh.
 * That is to say, the field content is changed by this method. The reason of this method is only for multi-field instances lying on the same mesh to
 * avoid a systematic duplication and renumbering of _mesh attribute.
 */
void MEDCouplingFieldDouble::renumberCellsWithoutMesh(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
   if(!_mesh)
    throw INTERP_KERNEL::Exception("Expecting a defined mesh to be able to operate a renumbering !");
  //
  _type->renumberCells(old2NewBg,check);
  std::vector<DataArrayDouble *> arrays;
  _time_discr->getArrays(arrays);
  _type->renumberArraysForCell(_mesh,arrays,old2NewBg,check);
  //
  updateTime();
}

/*!
 * This method performs a clone of mesh and a renumbering of underlying nodes of it. The number of nodes remains not compulsory the same as renumberCells method.
 * The values of field are impacted in consequence to have the same geometrical field.
 */
void MEDCouplingFieldDouble::renumberNodes(const int *old2NewBg) throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingPointSet *meshC=dynamic_cast<const MEDCouplingPointSet *>(_mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply renumberNodes on it !");
  int nbOfNodes=meshC->getNumberOfNodes();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingPointSet> meshC2((MEDCouplingPointSet *)meshC->deepCpy());
  renumberNodesWithoutMesh(old2NewBg);
  meshC2->renumberNodes(old2NewBg,*std::max_element(old2NewBg,old2NewBg+nbOfNodes)+1);
  setMesh(meshC2);
}

/*!
 * \b WARNING : use this method with lot of care !
 * This method performs half job of MEDCouplingFieldDouble::renumberNodes. That is to say no permutation of cells is done on underlying mesh.
 * That is to say, the field content is changed by this method.
 */
void MEDCouplingFieldDouble::renumberNodesWithoutMesh(const int *old2NewBg, double eps) throw(INTERP_KERNEL::Exception)
{
  std::vector<DataArrayDouble *> arrays;
  _time_discr->getArrays(arrays);
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    if(*iter)
      _type->renumberValuesOnNodes(eps,old2NewBg,*iter);
}

/*!
 * This method makes the assumption that the default array is set. If not an exception will be thrown.
 * This method is usable only if the default array has exactly one component. If not an exception will be thrown too.
 * This method returns all tuples ids that fit the range [vmin,vmax].
 * The caller has the responsability of the returned DataArrayInt.
 */
DataArrayInt *MEDCouplingFieldDouble::getIdsInRange(double vmin, double vmax) const throw(INTERP_KERNEL::Exception)
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::getIdsInRange : no default array set !");
  return getArray()->getIdsInRange(vmin,vmax);
}

/*!
 * Builds a newly created field, that the caller will have the responsability.
 * This method makes the assumption that the field is correctly defined when this method is called, no check of this will be done.
 * This method returns a restriction of 'this' so that only tuples id specified in 'part' will be contained in returned field. 
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::buildSubPart(const DataArrayInt *part) const throw(INTERP_KERNEL::Exception)
{
  if(part==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::buildSubPart : not empty array must be passed to this method !");
  const int *start=part->getConstPointer();
  const int *end=start+part->getNbOfElems();
  return buildSubPart(start,end);
}

/*!
 * Builds a newly created field, that the caller will have the responsability.
 * This method makes the assumption that the field is correctly defined when this method is called, no check of this will be done.
 * This method returns a restriction of 'this' so that only tuples id specified in ['partBg';'partEnd') will be contained in returned field. 
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::buildSubPart(const int *partBg, const int *partEnd) const throw(INTERP_KERNEL::Exception)
{
  DataArrayInt *cellRest;
  _type->computeMeshRestrictionFromTupleIds(_mesh,partBg,partEnd,cellRest);
  DataArrayInt *arrSelect;
  MEDCouplingMesh *m=_type->buildSubMeshData(_mesh,cellRest->getConstPointer(),cellRest->getConstPointer()+cellRest->getNbOfElems(),arrSelect);
  if(cellRest)
    cellRest->decrRef();
  MEDCouplingFieldDouble *ret=clone(false);//quick shallow copy.
  ret->setMesh(m);
  m->decrRef();
  std::vector<DataArrayDouble *> arrays;
  _time_discr->getArrays(arrays);
  std::vector<DataArrayDouble *> arrs;
  const int *arrSelBg=arrSelect->getConstPointer();
  const int *arrSelEnd=arrSelBg+arrSelect->getNbOfElems();
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    {
      DataArrayDouble *arr=0;
      if(*iter)
        arr=(*iter)->selectByTupleId(arrSelBg,arrSelEnd);
      arrs.push_back(arr);
    }
  ret->_time_discr->setArrays(arrs,0);
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrs.begin();iter!=arrs.end();iter++)
    if(*iter)
      (*iter)->decrRef();
  arrSelect->decrRef();
  return ret;
}

TypeOfTimeDiscretization MEDCouplingFieldDouble::getTimeDiscretization() const
{
  return _time_discr->getEnum();
}

MEDCouplingFieldDouble::MEDCouplingFieldDouble(TypeOfField type, TypeOfTimeDiscretization td):MEDCouplingField(type),
                                                                                              _time_discr(MEDCouplingTimeDiscretization::New(td))
{
}

MEDCouplingFieldDouble::MEDCouplingFieldDouble(const MEDCouplingFieldTemplate *ft, TypeOfTimeDiscretization td):MEDCouplingField(*ft),
                                                                                                                _time_discr(MEDCouplingTimeDiscretization::New(td))
{
}

MEDCouplingFieldDouble::MEDCouplingFieldDouble(const MEDCouplingFieldDouble& other, bool deepCpy):MEDCouplingField(other),
                                                                                                  _time_discr(other._time_discr->performCpy(deepCpy))
{
}

MEDCouplingFieldDouble::MEDCouplingFieldDouble(NatureOfField n, MEDCouplingTimeDiscretization *td, MEDCouplingFieldDiscretization *type):MEDCouplingField(type,n),_time_discr(td)
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
 * Returns the accumulation (the sum) of comId_th component of each tuples of \b default and \b only \b default array.
 */
double MEDCouplingFieldDouble::accumulate(int compId) const
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::accumulate : no default array defined !");
  return getArray()->accumulate(compId);
}

/*!
 * Returns the accumulation (the sum) of all tuples of \b default and \b only default array.
 * The res is expected to be of size getNumberOfComponents().
 */
void MEDCouplingFieldDouble::accumulate(double *res) const
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::accumulate : no default array defined !");
  getArray()->accumulate(res);
}

/*!
 * This method returns the max value in 'this'. 'This' is expected to be a field with exactly \b one component. If not an exception will be thrown.
 * To getMaxValue on vector field applyFunc is needed before. This method looks only on all arrays stored in 'this->_time_discr'.
 * If no arrays exists, an exception will be thrown.
 */
double MEDCouplingFieldDouble::getMaxValue() const throw(INTERP_KERNEL::Exception)
{
  std::vector<DataArrayDouble *> arrays;
  _time_discr->getArrays(arrays);
  double ret=-std::numeric_limits<double>::max();
  bool isExistingArr=false;
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    {
      if(*iter)
        {
          isExistingArr=true;
          int loc;
          ret=std::max(ret,(*iter)->getMaxValue(loc));
        }
    }
  if(!isExistingArr)
    throw INTERP_KERNEL::Exception("getMaxValue : No arrays defined !");
  return ret;
}

/*!
 * This method is an extension of ParaMEDMEM::MEDCouplingFieldDouble::getMaxValue method because the returned 
 * value is the same but this method also returns to you a tupleIds object which the caller have the responsibility
 * to deal with. The main difference is that the returned tupleIds is those corresponding the first set array.
 * If you have more than one array set (in LINEAR_TIME instance for example) only the first not null array will be used
 * to compute tupleIds.
 */
double MEDCouplingFieldDouble::getMaxValue2(DataArrayInt*& tupleIds) const throw(INTERP_KERNEL::Exception)
{
  std::vector<DataArrayDouble *> arrays;
  _time_discr->getArrays(arrays);
  double ret=-std::numeric_limits<double>::max();
  bool isExistingArr=false;
  tupleIds=0;
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    {
      if(*iter)
        {
          isExistingArr=true;
          DataArrayInt *tmp;
          ret=std::max(ret,(*iter)->getMaxValue2(tmp));
          if(!tupleIds)
            tupleIds=tmp;
          else
            tmp->decrRef();
        }
    }
  if(!isExistingArr)
    throw INTERP_KERNEL::Exception("getMaxValue2 : No arrays defined !");
  return ret;
}

/*!
 * This method returns the min value in 'this'. 'This' is expected to be a field with exactly \b one component. If not an exception will be thrown.
 * To getMinValue on vector field applyFunc is needed before. This method looks only on all arrays stored in 'this->_time_discr'.
 * If no arrays exists, an exception will be thrown.
 */
double MEDCouplingFieldDouble::getMinValue() const throw(INTERP_KERNEL::Exception)
{
  std::vector<DataArrayDouble *> arrays;
  _time_discr->getArrays(arrays);
  double ret=std::numeric_limits<double>::max();
  bool isExistingArr=false;
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    {
      if(*iter)
        {
          isExistingArr=true;
          int loc;
          ret=std::min(ret,(*iter)->getMinValue(loc));
        }
    }
  if(!isExistingArr)
    throw INTERP_KERNEL::Exception("getMinValue : No arrays defined !");
  return ret;
}

/*!
 * This method is an extension of ParaMEDMEM::MEDCouplingFieldDouble::getMinValue method because the returned 
 * value is the same but this method also returns to you a tupleIds object which the caller have the responsibility
 * to deal with. The main difference is that the returned tupleIds is those corresponding the first set array.
 * If you have more than one array set (in LINEAR_TIME instance for example) only the first not null array will be used
 * to compute tupleIds.
 */
double MEDCouplingFieldDouble::getMinValue2(DataArrayInt*& tupleIds) const throw(INTERP_KERNEL::Exception)
{
  std::vector<DataArrayDouble *> arrays;
  _time_discr->getArrays(arrays);
  double ret=-std::numeric_limits<double>::max();
  bool isExistingArr=false;
  tupleIds=0;
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    {
      if(*iter)
        {
          isExistingArr=true;
          DataArrayInt *tmp;
          ret=std::max(ret,(*iter)->getMinValue2(tmp));
          if(!tupleIds)
            tupleIds=tmp;
          else
            tmp->decrRef();
        }
    }
  if(!isExistingArr)
    throw INTERP_KERNEL::Exception("getMinValue2 : No arrays defined !");
  return ret;
}

/*!
 * This method returns the average value in 'this'. 'This' is expected to be a field with exactly \b one component. If not an exception will be thrown.
 * To getAverageValue on vector field applyFunc is needed before. This method looks only \b default array \b and \b only \b default.
 * If default array does not exist, an exception will be thrown.
 */
double MEDCouplingFieldDouble::getAverageValue() const throw(INTERP_KERNEL::Exception)
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::getAverageValue : no default array defined !");
  return getArray()->getAverageValue();
}

/*!
 * This method returns the euclidean norm of 'this'.
 * \f[
 * \sqrt{\sum_{0 \leq i < nbOfEntity}val[i]*val[i]}
 * \f]
 * If default array does not exist, an exception will be thrown.
 */
double MEDCouplingFieldDouble::norm2() const throw(INTERP_KERNEL::Exception)
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::norm2 : no default array defined !");
  return getArray()->norm2();
}

/*!
 * This method returns the max norm of 'this'.
 * \f[
 * \max_{0 \leq i < nbOfEntity}{abs(val[i])}
 * \f]
 * If default array does not exist, an exception will be thrown.
 */
double MEDCouplingFieldDouble::normMax() const throw(INTERP_KERNEL::Exception)
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::normMax : no default array defined !");
  return getArray()->normMax();
}

/*!
 * This method returns the average value in 'this' weighted by ParaMEDMEM::MEDCouplingField::buildMeasureField.
 * 'This' is expected to be a field with exactly \b one component. If not an exception will be thrown.
 * To getAverageValue on vector field applyFunc is needed before. This method looks only \b default array \b and \b only \b default.
 * If default array does not exist, an exception will be thrown.
 */
double MEDCouplingFieldDouble::getWeightedAverageValue() const throw(INTERP_KERNEL::Exception)
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::getWeightedAverageValue : no default array defined !");
  MEDCouplingFieldDouble *w=buildMeasureField(true);
  double deno=w->getArray()->accumulate(0);
  w->getArray()->multiplyEqual(getArray());
  double res=w->getArray()->accumulate(0);
  w->decrRef();
  return res/deno;
}

/*!
 * Returns the normL1 of current field on compId component :
 * \f[
 * \frac{\sum_{0 \leq i < nbOfEntity}|val[i]*Vol[i]|}{\sum_{0 \leq i < nbOfEntity}|Vol[i]|}
 * \f]
 * If compId>=nbOfComponent an exception is thrown.
 */
double MEDCouplingFieldDouble::normL1(int compId) const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform normL1");
  int nbComps=getArray()->getNumberOfComponents();
  if(compId>=nbComps)
    throw INTERP_KERNEL::Exception("Invalid compId specified : No such nb of components !");
  double *res=new double[nbComps];
  try
    {
      _type->normL1(_mesh,getArray(),res);
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
void MEDCouplingFieldDouble::normL1(double *res) const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform normL1");
  _type->normL1(_mesh,getArray(),res);
}

/*!
 * Returns the normL2 of current field on compId component :
 * \f[
 * \sqrt{\frac{\sum_{0 \leq i < nbOfEntity}|val[i]^{2}*Vol[i]|}{\sum_{0 \leq i < nbOfEntity}|Vol[i]|}}
 * \f]
 * If compId>=nbOfComponent an exception is thrown.
 */
double MEDCouplingFieldDouble::normL2(int compId) const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform normL2");
  int nbComps=getArray()->getNumberOfComponents();
  if(compId>=nbComps)
    throw INTERP_KERNEL::Exception("Invalid compId specified : No such nb of components !");
  double *res=new double[nbComps];
  try
    {
      _type->normL2(_mesh,getArray(),res);
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
void MEDCouplingFieldDouble::normL2(double *res) const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform normL2");
  _type->normL2(_mesh,getArray(),res);
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
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform getValueOnPos");
  _type->getValueOnPos(arr,_mesh,i,j,k,res);
}

/*!
 * Returns value of 'this' on default time of point 'spaceLoc' using spatial discretization.
 * If 'point' is outside the spatial discretization of this an exception will be thrown.
 */
void MEDCouplingFieldDouble::getValueOn(const double *spaceLoc, double *res) const throw(INTERP_KERNEL::Exception)
{
  const DataArrayDouble *arr=_time_discr->getArray();
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform getValueOn");
  _type->getValueOn(arr,_mesh,spaceLoc,res);
}

/*!
 * Returns a newly allocated array with 'nbOfPoints' tuples and nb of components equal to 'this->getNumberOfComponents()'.
 */
DataArrayDouble *MEDCouplingFieldDouble::getValueOnMulti(const double *spaceLoc, int nbOfPoints) const throw(INTERP_KERNEL::Exception)
{
  const DataArrayDouble *arr=_time_discr->getArray();
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform getValueOnMulti");
  return _type->getValueOnMulti(arr,_mesh,spaceLoc,nbOfPoints);
}

/*!
 * Returns value of 'this' on time 'time' of point 'spaceLoc' using spatial discretization.
 * If 'time' is not covered by this->_time_discr an exception will be thrown.
 * If 'point' is outside the spatial discretization of this an exception will be thrown.
 */
void MEDCouplingFieldDouble::getValueOn(const double *spaceLoc, double time, double *res) const throw(INTERP_KERNEL::Exception)
{
  std::vector< const DataArrayDouble *> arrs=_time_discr->getArraysForTime(time);
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform getValueOn");
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
 * This method sets 'this' to a uniform scalar field with one component.
 * All tuples will have the same value 'value'.
 * An exception is thrown if no underlying mesh is defined.
 */
MEDCouplingFieldDouble &MEDCouplingFieldDouble::operator=(double value) throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::operator= : no mesh defined !");
  int nbOfTuple=_type->getNumberOfTuples(_mesh);
  _time_discr->setUniformValue(nbOfTuple,1,value);
  return *this;
}

/*!
 * This method is very similar to this one MEDCouplingMesh::fillFromAnalytic.
 * See MEDCouplingMesh::fillFromAnalytic method doc to have more details.
 * The main difference is that the field as been started to be constructed here.
 * An exception is thrown if no underlying mesh is set before the call of this method.
 */
void MEDCouplingFieldDouble::fillFromAnalytic(int nbOfComp, FunctionToEvaluate func) throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::fillFromAnalytic : no mesh defined !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> loc=_type->getLocalizationOfDiscValues(_mesh);
  _time_discr->fillFromAnalytic(loc,nbOfComp,func);
}

/*!
 * This method is very similar to this one MEDCouplingMesh::fillFromAnalytic.
 * See MEDCouplingMesh::fillFromAnalytic method doc to have more details.
 * The main difference is that the field as been started to be constructed here.
 * An exception is thrown if no underlying mesh is set before the call of this method.
 */
void MEDCouplingFieldDouble::fillFromAnalytic(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::fillFromAnalytic : no mesh defined !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> loc=_type->getLocalizationOfDiscValues(_mesh);
  _time_discr->fillFromAnalytic(loc,nbOfComp,func);
}

/*!
 * This method is very similar to this one MEDCouplingMesh::fillFromAnalytic2.
 * See MEDCouplingMesh::fillFromAnalytic method doc to have more details.
 * The main difference is that the field as been started to be constructed here.
 * An exception is throw if no underlying mesh is set before the call of this method.
 */
void MEDCouplingFieldDouble::fillFromAnalytic2(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::fillFromAnalytic2 : no mesh defined !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> loc=_type->getLocalizationOfDiscValues(_mesh);
  _time_discr->fillFromAnalytic2(loc,nbOfComp,func);
}

/*!
 * This method is very similar to this one MEDCouplingMesh::fillFromAnalytic3.
 * See MEDCouplingMesh::fillFromAnalytic method doc to have more details.
 * The main difference is that the field as been started to be constructed here.
 * An exception is thrown if no underlying mesh is set before the call of this method.
 */
void MEDCouplingFieldDouble::fillFromAnalytic3(int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::fillFromAnalytic2 : no mesh defined !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> loc=_type->getLocalizationOfDiscValues(_mesh);
  _time_discr->fillFromAnalytic3(loc,nbOfComp,varsOrder,func);
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
 * This method is a specialization of other overloaded methods. When 'nbOfComp' equals 1 this method is equivalent to
 * ParaMEDMEM::MEDCouplingFieldDouble::operator=.
 */
void MEDCouplingFieldDouble::applyFunc(int nbOfComp, double val)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::applyFunc : no mesh defined !");
  int nbOfTuple=_type->getNumberOfTuples(_mesh);
  _time_discr->setUniformValue(nbOfTuple,nbOfComp,val);
}

/*!
 * Applyies the function specified by the string repr 'func' on each tuples on all arrays contained in _time_discr.
 * If '*func' fails in evaluation during one evaluation an exception will be thrown.
 * The field will contain 'nbOfComp' components after the call.
 */
void MEDCouplingFieldDouble::applyFunc(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception)
{
  _time_discr->applyFunc(nbOfComp,func);
}

/*!
 * This method is equivalent to MEDCouplingFieldDouble::applyFunc, except that here components info are used to determine variables position in 'func'.
 * If there is vars detected in 'func' that is not in an info on components an exception will be thrown.
 */
void MEDCouplingFieldDouble::applyFunc2(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception)
{
  _time_discr->applyFunc2(nbOfComp,func);
}

/*!
 * This method is equivalent to MEDCouplingFieldDouble::applyFunc, except that here 'varsOrder' is used to determine variables position in 'func'.
 * If there is vars detected in 'func' that is not in 'varsOrder' an exception will be thrown.
 */
void MEDCouplingFieldDouble::applyFunc3(int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) throw(INTERP_KERNEL::Exception)
{
  _time_discr->applyFunc3(nbOfComp,varsOrder,func);
}

/*!
 * Applyies the function specified by the string repr 'func' on each tuples on all arrays contained in _time_discr.
 * If '*func' fails in evaluation during one evaluation an exception will be thrown.
 * The field will contain exactly the same number of components after the call.
 */
void MEDCouplingFieldDouble::applyFunc(const char *func) throw(INTERP_KERNEL::Exception)
{
  _time_discr->applyFunc(func);
}

/*!
 * Applyies the function specified by the string repr 'func' on each tuples on all arrays contained in _time_discr.
 * The field will contain exactly the same number of components after the call.
 * Use is not warranted for the moment !
 */
void MEDCouplingFieldDouble::applyFuncFast32(const char *func) throw(INTERP_KERNEL::Exception)
{
  _time_discr->applyFuncFast32(func);
}

/*!
 * Applyies the function specified by the string repr 'func' on each tuples on all arrays contained in _time_discr.
 * The field will contain exactly the same number of components after the call.
 * Use is not warranted for the moment !
 */
void MEDCouplingFieldDouble::applyFuncFast64(const char *func) throw(INTERP_KERNEL::Exception)
{
  _time_discr->applyFuncFast64(func);
}

/*!
 * This method makes the assumption that the default array has been set before.
 * If not an exception will be sent.
 * If default array set, the number of components will be sent.
 */
int MEDCouplingFieldDouble::getNumberOfComponents() const throw(INTERP_KERNEL::Exception)
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::getNumberOfComponents : No array specified !");
  return getArray()->getNumberOfComponents();
}

/*!
 * This method makes the assumption that _mesh has be set before the call of this method and description of gauss
 * localizations in case of Gauss field. If not an exception will sent.
 * \b Contrary to MEDCouplingFieldDouble::getNumberOfComponents and MEDCouplingFieldDouble::getNumberOfValues is
 * \b not aware of the presence of the default array.
 * \b WARNING \b no coherency check is done here. MEDCouplingFieldDouble::checkCoherency method should be called to check that !
 */
int MEDCouplingFieldDouble::getNumberOfTuples() const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Impossible to retrieve number of tuples because no mesh specified !");
  return _type->getNumberOfTuples(_mesh);
}

/*!
 * This method makes the assumption that the default array has been set before.
 * If not an exception will be sent.
 * If default array set, the number of values present in the default array will be sent.
 */
int MEDCouplingFieldDouble::getNumberOfValues() const throw(INTERP_KERNEL::Exception)
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::getNumberOfValues : No array specified !");
  return getArray()->getNbOfElems();
}

void MEDCouplingFieldDouble::updateTime() const
{
  MEDCouplingField::updateTime();
  updateTimeWith(*_time_discr);
}

void MEDCouplingFieldDouble::setNature(NatureOfField nat) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingField::setNature(nat);
  _type->checkCompatibilityWithNature(nat);
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

void MEDCouplingFieldDouble::setArrays(const std::vector<DataArrayDouble *>& arrs) throw(INTERP_KERNEL::Exception)
{
  _time_discr->setArrays(arrs,this);
}

void MEDCouplingFieldDouble::getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const
{
  tinyInfo.clear();
  _time_discr->getTinySerializationStrInformation(tinyInfo);
  tinyInfo.push_back(_name);
  tinyInfo.push_back(_desc);
  tinyInfo.push_back(getTimeUnit());
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
  _name=tinyInfoS[nbOfElemS-3];
  _desc=tinyInfoS[nbOfElemS-2];
  setTimeUnit(tinyInfoS[nbOfElemS-1].c_str());
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
 * This method tries to to change the mesh support of 'this' following the parameter 'levOfCheck' and 'prec'.
 * Semantic of 'levOfCheck' is explained in MEDCouplingMesh::checkGeoEquivalWith method. This method is used to perform the job.
 * If this->_mesh is not defined or other an exeption will be throw.
 */
void MEDCouplingFieldDouble::changeUnderlyingMesh(const MEDCouplingMesh *other, int levOfCheck, double prec) throw(INTERP_KERNEL::Exception)
{
  if(_mesh==0 || other==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::changeUnderlyingMesh : is expected to operate on not null meshes !");
  DataArrayInt *cellCor,*nodeCor;
  other->checkGeoEquivalWith(_mesh,levOfCheck,prec,cellCor,nodeCor);
  if(cellCor)
    {
      renumberCellsWithoutMesh(cellCor->getConstPointer(),false);
      cellCor->decrRef();
    }
  if(nodeCor)
    {
      renumberNodesWithoutMesh(nodeCor->getConstPointer());
      nodeCor->decrRef();
    }
  setMesh((MEDCouplingMesh *)other);
}

/*!
 * This method is an extension of MEDCouplingFieldDouble::operator-=. It allows a user to operate a difference of 2 fields ('this' and 'f') even if they do not share same meshes.
 * No interpolation will be done here only an analyze of two underlying mesh will be done to see if the meshes are geometrically equivalent. If yes, the eventual renumbering will be done and operator-= applyed after.
 * This method requires that 'f' and 'this' are coherent (check coherency) and that 'f' and 'this' would be coherent for a merge.
 * Semantic of 'levOfCheck' is explained in MEDCouplingMesh::checkGeoEquivalWith method.
 */
void MEDCouplingFieldDouble::substractInPlaceDM(const MEDCouplingFieldDouble *f, int levOfCheck, double prec) throw(INTERP_KERNEL::Exception)
{
  checkCoherency();
  f->checkCoherency();
  if(!areCompatibleForMerge(f))
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::diffWith : Fields are not compatible ; unable to apply mergeFields on them !");
  changeUnderlyingMesh(f->getMesh(),levOfCheck,prec);
  operator-=(*f);
}

/*!
 * Merge nodes of underlying mesh. In case of some node will be merged the underlying mesh instance will change.
 * The first 'eps' stands for geometric approximation. The second 'epsOnVals' is for epsilon on values in case of node merging.
 * If 2 nodes distant from less than 'eps' and with value different with more than 'epsOnVals' an exception will be thrown.
 */
bool MEDCouplingFieldDouble::mergeNodes(double eps, double epsOnVals) throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingPointSet *meshC=dynamic_cast<const MEDCouplingPointSet *>(_mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("Invalid support mesh to apply mergeNodes on it : must be a MEDCouplingPointSet one !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingPointSet> meshC2((MEDCouplingPointSet *)meshC->deepCpy());
  bool ret;
  int ret2;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=meshC2->mergeNodes(eps,ret,ret2);
  if(!ret)//no nodes have been merged.
    return ret;
  std::vector<DataArrayDouble *> arrays;
  _time_discr->getArrays(arrays);
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    if(*iter)
      _type->renumberValuesOnNodes(epsOnVals,arr->getConstPointer(),*iter);
  setMesh(meshC2);
  return true;
}

/*!
 * Merge nodes with (barycenter computation) of underlying mesh. In case of some node will be merged the underlying mesh instance will change.
 * The first 'eps' stands for geometric approximation. The second 'epsOnVals' is for epsilon on values in case of node merging.
 * If 2 nodes distant from less than 'eps' and with value different with more than 'epsOnVals' an exception will be thrown.
 */
bool MEDCouplingFieldDouble::mergeNodes2(double eps, double epsOnVals) throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingPointSet *meshC=dynamic_cast<const MEDCouplingPointSet *>(_mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("Invalid support mesh to apply mergeNodes on it : must be a MEDCouplingPointSet one !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingPointSet> meshC2((MEDCouplingPointSet *)meshC->deepCpy());
  bool ret;
  int ret2;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=meshC2->mergeNodes2(eps,ret,ret2);
  if(!ret)//no nodes have been merged.
    return ret;
  std::vector<DataArrayDouble *> arrays;
  _time_discr->getArrays(arrays);
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    if(*iter)
      _type->renumberValuesOnNodes(epsOnVals,arr->getConstPointer(),*iter);
  setMesh(meshC2);
  return true;
}

/*!
 * This method applyies ParaMEDMEM::MEDCouplingPointSet::zipCoords method on 'this->_mesh' that should be set and of type ParaMEDMEM::MEDCouplingPointSet.
 * If some nodes have disappeared true is returned.
 * 'epsOnVals' stands for epsilon in case of merge of cells. This value is used as tolerance in case the corresponding values differ.
 */
bool MEDCouplingFieldDouble::zipCoords(double epsOnVals) throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingPointSet *meshC=dynamic_cast<const MEDCouplingPointSet *>(_mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::zipCoords : Invalid support mesh to apply zipCoords on it : must be a MEDCouplingPointSet one !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingPointSet> meshC2((MEDCouplingPointSet *)meshC->deepCpy());
  int oldNbOfNodes=meshC2->getNumberOfNodes();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=meshC2->zipCoordsTraducer();
  if(meshC2->getNumberOfNodes()!=oldNbOfNodes)
    {
      std::vector<DataArrayDouble *> arrays;
      _time_discr->getArrays(arrays);
      for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
        if(*iter)
          _type->renumberValuesOnNodes(epsOnVals,arr->getConstPointer(),*iter);
      setMesh(meshC2);
      return true;
    }
  return false;
}

/*!
 * This method applyies ParaMEDMEM::MEDCouplingUMesh::zipConnectivityTraducer on 'this->_mesh' that should be set and of type ParaMEDMEM::MEDCouplingUMesh.
 * The semantic of 'compType' is given in ParaMEDMEM::MEDCouplingUMesh::zipConnectivityTraducer method.
 * 'epsOnVals' stands for epsilon in case of merge of cells. This value is used as tolerance in case the corresponding values differ.
 */
bool MEDCouplingFieldDouble::zipConnectivity(int compType, double epsOnVals) throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingUMesh *meshC=dynamic_cast<const MEDCouplingUMesh *>(_mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::zipCoords : Invalid support mesh to apply zipCoords on it : must be a MEDCouplingPointSet one !");
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> meshC2((MEDCouplingUMesh *)meshC->deepCpy());
  int oldNbOfCells=meshC2->getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=meshC2->zipConnectivityTraducer(compType);
  if(meshC2->getNumberOfCells()!=oldNbOfCells)
    {
      std::vector<DataArrayDouble *> arrays;
      _time_discr->getArrays(arrays);
      for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
        if(*iter)
          _type->renumberValuesOnCells(epsOnVals,meshC,arr->getConstPointer(),*iter);
      setMesh(meshC2);
      return true;
    }
  return false;
}

/*!
 * This method applyies ParaMEDMEM::MEDCouplingUMesh::simplexize on 'this->_mesh'.
 * The semantic of 'policy' is given in ParaMEDMEM::MEDCouplingUMesh::simplexize method.
 */
bool MEDCouplingFieldDouble::simplexize(int policy) throw(INTERP_KERNEL::Exception)
{
  int oldNbOfCells=_mesh->getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> meshC2(_mesh->deepCpy());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=meshC2->simplexize(policy);
  int newNbOfCells=meshC2->getNumberOfCells();
  if(oldNbOfCells==newNbOfCells)
    return false;
  std::vector<DataArrayDouble *> arrays;
  _time_discr->getArrays(arrays);
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    if(*iter)
      _type->renumberValuesOnCellsR(_mesh,arr->getConstPointer(),arr->getNbOfElems(),*iter);
  setMesh(meshC2);
  return true;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::doublyContractedProduct() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingTimeDiscretization *td=_time_discr->doublyContractedProduct();
  td->copyTinyAttrFrom(*_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(getNature(),td,_type->clone());
  ret->setName("DoublyContractedProduct");
  ret->setMesh(getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::determinant() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingTimeDiscretization *td=_time_discr->determinant();
  td->copyTinyAttrFrom(*_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(getNature(),td,_type->clone());
  ret->setName("Determinant");
  ret->setMesh(getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::eigenValues() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingTimeDiscretization *td=_time_discr->eigenValues();
  td->copyTinyAttrFrom(*_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(getNature(),td,_type->clone());
  ret->setName("EigenValues");
  ret->setMesh(getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::eigenVectors() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingTimeDiscretization *td=_time_discr->eigenVectors();
  td->copyTinyAttrFrom(*_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(getNature(),td,_type->clone());
  ret->setName("EigenVectors");
  ret->setMesh(getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::inverse() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingTimeDiscretization *td=_time_discr->inverse();
  td->copyTinyAttrFrom(*_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(getNature(),td,_type->clone());
  ret->setName("Inversion");
  ret->setMesh(getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::trace() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingTimeDiscretization *td=_time_discr->trace();
  td->copyTinyAttrFrom(*_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(getNature(),td,_type->clone());
  ret->setName("Trace");
  ret->setMesh(getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::deviator() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingTimeDiscretization *td=_time_discr->deviator();
  td->copyTinyAttrFrom(*_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(getNature(),td,_type->clone());
  ret->setName("Trace");
  ret->setMesh(getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::magnitude() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingTimeDiscretization *td=_time_discr->magnitude();
  td->copyTinyAttrFrom(*_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(getNature(),td,_type->clone());
  ret->setName("Magnitude");
  ret->setMesh(getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::maxPerTuple() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingTimeDiscretization *td=_time_discr->maxPerTuple();
  td->copyTinyAttrFrom(*_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(getNature(),td,_type->clone());
  std::ostringstream oss;
  oss << "Max_" << getName();
  ret->setName(oss.str().c_str());
  ret->setMesh(getMesh());
  return ret;
}

void MEDCouplingFieldDouble::changeNbOfComponents(int newNbOfComp, double dftValue) throw(INTERP_KERNEL::Exception)
{
  _time_discr->changeNbOfComponents(newNbOfComp,dftValue);
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::keepSelectedComponents(const std::vector<int>& compoIds) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingTimeDiscretization *td=_time_discr->keepSelectedComponents(compoIds);
  td->copyTinyAttrFrom(*_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(getNature(),td,_type->clone());
  ret->setName(getName());
  ret->setMesh(getMesh());
  return ret;
}

void MEDCouplingFieldDouble::setSelectedComponents(const MEDCouplingFieldDouble *f, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception)
{
  _time_discr->setSelectedComponents(f->_time_discr,compoIds);
}

void MEDCouplingFieldDouble::sortPerTuple(bool asc) throw(INTERP_KERNEL::Exception)
{
  _time_discr->sortPerTuple(asc);
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::MergeFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception)
{
  if(!f1->areCompatibleForMerge(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply MergeFields on them !");
  const MEDCouplingMesh *m1=f1->getMesh();
  const MEDCouplingMesh *m2=f2->getMesh();
  MEDCouplingMesh *m=m1->mergeMyselfWith(m2);
  MEDCouplingTimeDiscretization *td=f1->_time_discr->aggregate(f2->_time_discr);
  td->copyTinyAttrFrom(*f1->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone());
  ret->setMesh(m);
  m->decrRef();
  ret->setName(f1->getName());
  ret->setDescription(f1->getDescription());
  return ret;
}

/*!
 * This method returns a newly created field that is the union of all fields in input array 'a'.
 * This method expects that 'a' is non empty. If not an exception will be thrown.
 * If there is only one field in 'a' a deepCopy (except time information of mesh and field) of the unique field instance in 'a' will be returned.
 * Generally speaking the first instance field in 'a' will be used to assign tiny attributes of returned field.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::MergeFields(const std::vector<const MEDCouplingFieldDouble *>& a) throw(INTERP_KERNEL::Exception)
{
  if(a.size()<1)
    throw INTERP_KERNEL::Exception("FieldDouble::MergeFields : size of array must be >= 1 !");
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> > ms(a.size());
  std::vector< const MEDCouplingUMesh *> ms2(a.size());
  std::vector< const MEDCouplingTimeDiscretization *> tds(a.size());
  std::vector<const MEDCouplingFieldDouble *>::const_iterator it=a.begin();
  const MEDCouplingFieldDouble *ref=(*it++);
  for(;it!=a.end();it++)
    if(!ref->areCompatibleForMerge(*it))
      throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply MergeFields on them !");
  for(int i=0;i<(int)a.size();i++)
    {
      if(!a[i]->getMesh())
        throw INTERP_KERNEL::Exception("MergeFields : A field as no underlying mesh !");
      ms[i]=a[i]->getMesh()->buildUnstructured();
      ms2[i]=ms[i];
      tds[i]=a[i]->_time_discr;
    }
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m=MEDCouplingUMesh::MergeUMeshes(ms2);
  m->setName(ms2[0]->getName()); m->setDescription(ms2[0]->getDescription());
  MEDCouplingTimeDiscretization *td=tds[0]->aggregate(tds);
  td->copyTinyAttrFrom(*(a[0]->_time_discr));
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(a[0]->getNature(),td,a[0]->_type->clone());
  ret->setMesh(m);
  ret->setName(a[0]->getName());
  ret->setDescription(a[0]->getDescription());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::MeldFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception)
{
  if(!f1->areCompatibleForMeld(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply MeldFields on them !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->meld(f2->_time_discr);
  td->copyTinyAttrFrom(*f1->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone());
  ret->setMesh(f1->getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::DotFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception)
{
  if(!f1->areStrictlyCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply DotFields on them !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->dot(f2->_time_discr);
  td->copyTinyAttrFrom(*f1->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone());
  ret->setMesh(f1->getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::CrossProductFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception)
{
  if(!f1->areStrictlyCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply CrossProductFields on them !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->crossProduct(f2->_time_discr);
  td->copyTinyAttrFrom(*f1->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone());
  ret->setMesh(f1->getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::MaxFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception)
{
  if(!f1->areStrictlyCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply MaxFields on them !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->max(f2->_time_discr);
  td->copyTinyAttrFrom(*f1->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone());
  ret->setMesh(f1->getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::MinFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception)
{
  if(!f1->areStrictlyCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply MinFields on them !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->min(f2->_time_discr);
  td->copyTinyAttrFrom(*f1->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone());
  ret->setMesh(f1->getMesh());
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::AddFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception)
{
  if(!f1->areStrictlyCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply AddFields on them !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->add(f2->_time_discr);
  td->copyTinyAttrFrom(*f1->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone());
  ret->setMesh(f1->getMesh());
  return ret;
}

const MEDCouplingFieldDouble &MEDCouplingFieldDouble::operator+=(const MEDCouplingFieldDouble& other) throw(INTERP_KERNEL::Exception)
{
  if(!areStrictlyCompatible(&other))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply += on them !");
  _time_discr->addEqual(other._time_discr);
  return *this;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::SubstractFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception)
{
  if(!f1->areStrictlyCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply SubstractFields on them !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->substract(f2->_time_discr);
  td->copyTinyAttrFrom(*f1->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone());
  ret->setMesh(f1->getMesh());
  return ret;
}

const MEDCouplingFieldDouble &MEDCouplingFieldDouble::operator-=(const MEDCouplingFieldDouble& other) throw(INTERP_KERNEL::Exception)
{
  if(!areStrictlyCompatible(&other))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply -= on them !");
  _time_discr->substractEqual(other._time_discr);
  return *this;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::MultiplyFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception)
{
  if(!f1->areCompatibleForMul(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply MultiplyFields on them !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->multiply(f2->_time_discr);
  td->copyTinyAttrFrom(*f1->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone());
  ret->setMesh(f1->getMesh());
  return ret;
}

const MEDCouplingFieldDouble &MEDCouplingFieldDouble::operator*=(const MEDCouplingFieldDouble& other) throw(INTERP_KERNEL::Exception)
{
  if(!areCompatibleForMul(&other))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply *= on them !");
  _time_discr->multiplyEqual(other._time_discr);
  return *this;
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::DivideFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception)
{
  if(!f1->areCompatibleForDiv(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply DivideFields on them !");
  MEDCouplingTimeDiscretization *td=f1->_time_discr->divide(f2->_time_discr);
  td->copyTinyAttrFrom(*f1->_time_discr);
  MEDCouplingFieldDouble *ret=new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone());
  ret->setMesh(f1->getMesh());
  return ret;
}

const MEDCouplingFieldDouble &MEDCouplingFieldDouble::operator/=(const MEDCouplingFieldDouble& other) throw(INTERP_KERNEL::Exception)
{
  if(!areCompatibleForDiv(&other))
    throw INTERP_KERNEL::Exception("Fields are not compatible ; unable to apply /= on them !");
  _time_discr->divideEqual(other._time_discr);
  return *this;
}

/*!
 * This method writes the field series 'fs' in the VTK file 'fileName'.
 * If 'fs' is empty no file is written. If fields lies on more than one mesh an exception will be thrown and no file will be written too.
 * If the single mesh is empty an exception will be thrown.
 * Finally there is a field in 'fs' with no name an exception will be thrown too.
 */
void MEDCouplingFieldDouble::WriteVTK(const char *fileName, const std::vector<const MEDCouplingFieldDouble *>& fs) throw(INTERP_KERNEL::Exception)
{
  if(fs.empty())
    return;
  std::size_t nfs=fs.size();
  const MEDCouplingMesh *m=fs[0]->getMesh();
  for(std::size_t i=1;i<nfs;i++)
    if(fs[i]->getMesh()!=m)
      throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::WriteVTK : Fields are not lying on a same mesh ! Expected by VTK ! MEDCouplingFieldDouble::setMesh or MEDCouplingFieldDouble::changeUnderlyingMesh can help to that.");
  if(!m)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::WriteVTK : Fields are lying on a same mesh but it is empty !");
  std::ostringstream coss,noss;
  for(std::size_t i=0;i<nfs;i++)
    {
      const MEDCouplingFieldDouble *cur=fs[i];
      std::string name(cur->getName());
      if(name.empty())
        {
          std::ostringstream oss; oss << "MEDCouplingFieldDouble::WriteVTK : Field in pos #" << i << " has no name !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      TypeOfField typ=cur->getTypeOfField();
      if(typ==ON_CELLS)
        cur->getArray()->writeVTK(coss,8,cur->getName());
      else if(typ==ON_NODES)
        cur->getArray()->writeVTK(noss,8,cur->getName());
    }
  m->writeVTKAdvanced(fileName,coss.str(),noss.str());
}
