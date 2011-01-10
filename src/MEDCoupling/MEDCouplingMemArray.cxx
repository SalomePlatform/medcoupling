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

#include "MEDCouplingMemArray.txx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include "GenMathFormulae.hxx"
#include "InterpKernelExprParser.hxx"

#include <set>
#include <cmath>
#include <limits>
#include <numeric>
#include <functional>

typedef double (*MYFUNCPTR)(double);

using namespace ParaMEDMEM;

void DataArray::setName(const char *name)
{
  _name=name;
}

void DataArray::copyStringInfoFrom(const DataArray& other) throw(INTERP_KERNEL::Exception)
{
  if(_info_on_compo.size()!=other._info_on_compo.size())
    throw INTERP_KERNEL::Exception("Size of arrays mismatches on copyStringInfoFrom !");
  _name=other._name;
  _info_on_compo=other._info_on_compo;
}

void DataArray::copyPartOfStringInfoFrom(const DataArray& other, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception)
{
  int nbOfCompoOth=other.getNumberOfComponents();
  int newNbOfCompo=compoIds.size();
  for(int i=0;i<newNbOfCompo;i++)
    if(compoIds[i]>=nbOfCompoOth)
      {
        std::ostringstream oss; oss << "Specified component ids is to high (" << compoIds[i] << ") compared with nb of actual components (" << nbOfCompoOth << ")";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  for(int i=0;i<newNbOfCompo;i++)
    setInfoOnComponent(i,other.getInfoOnComponent(compoIds[i]).c_str());
}

void DataArray::copyPartOfStringInfoFrom2(const std::vector<int>& compoIds, const DataArray& other) throw(INTERP_KERNEL::Exception)
{
  int nbOfCompo=getNumberOfComponents();
  int partOfCompoToSet=compoIds.size();
  if(partOfCompoToSet!=other.getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Given compoIds has not the same size as number of components of given array !");
  for(int i=0;i<partOfCompoToSet;i++)
    if(compoIds[i]>=nbOfCompo)
      {
        std::ostringstream oss; oss << "Specified component ids is to high (" << compoIds[i] << ") compared with nb of actual components (" << nbOfCompo << ")";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  for(int i=0;i<partOfCompoToSet;i++)
    setInfoOnComponent(compoIds[i],other.getInfoOnComponent(i).c_str());
}

bool DataArray::areInfoEquals(const DataArray& other) const
{
  if(_nb_of_tuples!=other._nb_of_tuples)
    return false;
  if(_name!=other._name)
    return false;
  return _info_on_compo==other._info_on_compo;
}

void DataArray::reprWithoutNameStream(std::ostream& stream) const
{
  stream << "Number of components : "<< getNumberOfComponents() << "\n";
  stream << "Info of these components : ";
  for(std::vector<std::string>::const_iterator iter=_info_on_compo.begin();iter!=_info_on_compo.end();iter++)
    stream << "\"" << *iter << "\"   ";
  stream << "\n";
}

std::string DataArray::getInfoOnComponent(int i) const throw(INTERP_KERNEL::Exception)
{
  if(i<(int)_info_on_compo.size())
    return _info_on_compo[i];
  else
    {
      std::ostringstream oss; oss << "getInfoOnComponent : Invalid component id transmitted (" << i << ") >= " << (int) _info_on_compo.size();
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void DataArray::setInfoOnComponent(int i, const char *info) throw(INTERP_KERNEL::Exception)
{
  if(i<(int)_info_on_compo.size())
    _info_on_compo[i]=info;
  else
    {
      std::ostringstream oss; oss << "setInfoOnComponent : Invalid component id transmitted (" << i << ") >= " << (int) _info_on_compo.size();
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

DataArrayDouble *DataArrayDouble::New()
{
  return new DataArrayDouble;
}

bool DataArrayDouble::isAllocated() const
{
  return getConstPointer()!=0;
}

void DataArrayDouble::checkAllocated() const throw(INTERP_KERNEL::Exception)
{
  if(!isAllocated())
    throw INTERP_KERNEL::Exception("DataArrayDouble::checkAllocated : Array is defined but not allocated ! Call alloc or setValues method first !");
}

DataArrayDouble *DataArrayDouble::deepCpy() const
{
  return new DataArrayDouble(*this);
}

DataArrayDouble *DataArrayDouble::performCpy(bool dCpy) const
{
  if(dCpy)
    return deepCpy();
  else
    {
      incrRef();
      return const_cast<DataArrayDouble *>(this);
    }
}

void DataArrayDouble::alloc(int nbOfTuple, int nbOfCompo)
{
  _nb_of_tuples=nbOfTuple;
  _info_on_compo.resize(nbOfCompo);
  _mem.alloc(nbOfCompo*_nb_of_tuples);
  declareAsNew();
}

void DataArrayDouble::fillWithZero() throw(INTERP_KERNEL::Exception)
{
  fillWithValue(0.);
}

void DataArrayDouble::fillWithValue(double val) throw(INTERP_KERNEL::Exception)
{
   if(!getPointer())
    throw INTERP_KERNEL::Exception("DataArrayDouble::fillWithValue : allocate first !");
   _mem.fillWithValue(val);
  declareAsNew();
}

void DataArrayDouble::iota(double init) throw(INTERP_KERNEL::Exception)
{
  double *ptr=getPointer();
  if(!ptr)
    throw INTERP_KERNEL::Exception("DataArrayDouble::iota : allocate first !");
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::iota : works only for arrays with only one component!");
  int ntuples=getNumberOfTuples();
  for(int i=0;i<ntuples;i++)
    ptr[i]=init+double(i);
  declareAsNew();
}

bool DataArrayDouble::isUniform(double val, double eps) const
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::isUniform : must be applied on DataArrayDouble with only one component !");
  int nbOfTuples=getNumberOfTuples();
  const double *w=getConstPointer();
  const double *end=w+nbOfTuples;
  const double vmin=val-eps;
  const double vmax=val+eps;
  for(;w!=end;w++)
    if(*w<vmin || *w>vmax)
      return false;
  return true;
}

void DataArrayDouble::sort() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::sort : only supported with 'this' array with ONE component !");
  _mem.sort();
}

std::string DataArrayDouble::repr() const
{
  std::ostringstream ret;
  reprStream(ret);
  return ret.str();
}

std::string DataArrayDouble::reprZip() const
{
  std::ostringstream ret;
  reprZipStream(ret);
  return ret.str();
}

void DataArrayDouble::reprStream(std::ostream& stream) const
{
  stream << "Name of double array : \"" << _name << "\"\n";
  reprWithoutNameStream(stream);
}

void DataArrayDouble::reprZipStream(std::ostream& stream) const
{
  stream << "Name of double array : \"" << _name << "\"\n";
  reprZipWithoutNameStream(stream);
}

void DataArrayDouble::reprWithoutNameStream(std::ostream& stream) const
{
  DataArray::reprWithoutNameStream(stream);
  stream.precision(15);
  _mem.repr(getNumberOfComponents(),stream);
}

void DataArrayDouble::reprZipWithoutNameStream(std::ostream& stream) const
{
  DataArray::reprWithoutNameStream(stream);
  stream.precision(15);
  _mem.reprZip(getNumberOfComponents(),stream);
}

bool DataArrayDouble::isEqual(const DataArrayDouble& other, double prec) const
{
  if(!areInfoEquals(other))
    return false;
  return _mem.isEqual(other._mem,prec);
}

bool DataArrayDouble::isEqualWithoutConsideringStr(const DataArrayDouble& other, double prec) const
{
  return _mem.isEqual(other._mem,prec);
}

void DataArrayDouble::reAlloc(int nbOfTuples) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  _mem.reAlloc(_info_on_compo.size()*nbOfTuples);
  _nb_of_tuples=nbOfTuples;
  declareAsNew();
}

DataArrayInt *DataArrayDouble::convertToIntArr() const
{
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(getNumberOfTuples(),getNumberOfComponents());
  int nbOfVals=getNbOfElems();
  const double *src=getConstPointer();
  int *dest=ret->getPointer();
  std::copy(src,src+nbOfVals,dest);
  ret->copyStringInfoFrom(*this);
  return ret;
}

DataArrayDouble *DataArrayDouble::fromNoInterlace() const throw(INTERP_KERNEL::Exception)
{
  if(_mem.isNull())
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromNoInterlace : Not defined array !");
  double *tab=_mem.fromNoInterlace(getNumberOfComponents());
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->useArray(tab,true,CPP_DEALLOC,getNumberOfTuples(),getNumberOfComponents());
  return ret;
}

DataArrayDouble *DataArrayDouble::toNoInterlace() const throw(INTERP_KERNEL::Exception)
{
  if(_mem.isNull())
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromNoInterlace : Not defined array !");
  double *tab=_mem.toNoInterlace(getNumberOfComponents());
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->useArray(tab,true,CPP_DEALLOC,getNumberOfTuples(),getNumberOfComponents());
  return ret;
}

/*!
 * This method does \b not change the number of tuples after this call.
 * Only a permutation is done. If a permutation reduction is needed substr, or selectByTupleId should be used.
 */
void DataArrayDouble::renumberInPlace(const int *old2New)
{
  int nbTuples=getNumberOfTuples();
  int nbOfCompo=getNumberOfComponents();
  double *tmp=new double[nbTuples*nbOfCompo];
  const double *iptr=getConstPointer();
  for(int i=0;i<nbTuples;i++)
    std::copy(iptr+nbOfCompo*i,iptr+nbOfCompo*(i+1),tmp+nbOfCompo*old2New[i]);
  std::copy(tmp,tmp+nbTuples*nbOfCompo,getPointer());
  delete [] tmp;
  declareAsNew();
}

/*!
 * This method does \b not change the number of tuples after this call.
 * Only a permutation is done.
 */
void DataArrayDouble::renumberInPlaceR(const int *new2Old)
{
  int nbTuples=getNumberOfTuples();
  int nbOfCompo=getNumberOfComponents();
  double *tmp=new double[nbTuples*nbOfCompo];
  const double *iptr=getConstPointer();
  for(int i=0;i<nbTuples;i++)
    std::copy(iptr+nbOfCompo*new2Old[i],iptr+nbOfCompo*(new2Old[i]+1),tmp+nbOfCompo*i);
  std::copy(tmp,tmp+nbTuples*nbOfCompo,getPointer());
  delete [] tmp;
  declareAsNew();
}

/*!
 * This method does \b not change the number of tuples after this call.
 * Only a permutation is done. If a permutation reduction is needed substr, or selectByTupleId should be used.
 */
DataArrayDouble *DataArrayDouble::renumber(const int *old2New) const
{
  int nbTuples=getNumberOfTuples();
  int nbOfCompo=getNumberOfComponents();
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbTuples,nbOfCompo);
  ret->copyStringInfoFrom(*this);
  const double *iptr=getConstPointer();
  double *optr=ret->getPointer();
  for(int i=0;i<nbTuples;i++)
    std::copy(iptr+nbOfCompo*i,iptr+nbOfCompo*(i+1),optr+nbOfCompo*old2New[i]);
  ret->copyStringInfoFrom(*this);
  return ret;
}

/*!
 * This method does \b not change the number of tuples after this call.
 * Only a permutation is done.
 */
DataArrayDouble *DataArrayDouble::renumberR(const int *new2Old) const
{
  int nbTuples=getNumberOfTuples();
  int nbOfCompo=getNumberOfComponents();
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbTuples,nbOfCompo);
  ret->copyStringInfoFrom(*this);
  const double *iptr=getConstPointer();
  double *optr=ret->getPointer();
  for(int i=0;i<nbTuples;i++)
    std::copy(iptr+nbOfCompo*new2Old[i],iptr+nbOfCompo*(new2Old[i]+1),optr+i*nbOfCompo);
  ret->copyStringInfoFrom(*this);
  return ret;
}

/*!
 * Idem DataArrayDouble::renumber method except that the number of tuples is reduced.
 * That is to say that it is expected that newNbOfTuple<this->getNumberOfTuples().
 * ['old2New','old2New'+getNumberOfTuples()) defines a range containing old to new array. For every negative value in ['old2NewBg','old2New'+getNumberOfTuples()) the corresponding tuple is
 * omitted.
 */
DataArrayDouble *DataArrayDouble::renumberAndReduce(const int *old2New, int newNbOfTuple) const
{
  int nbTuples=getNumberOfTuples();
  int nbOfCompo=getNumberOfComponents();
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(newNbOfTuple,nbOfCompo);
  const double *iptr=getConstPointer();
  double *optr=ret->getPointer();
  for(int i=0;i<nbTuples;i++)
    {
      int w=old2New[i];
      if(w>=0)
        std::copy(iptr+i*nbOfCompo,iptr+(i+1)*nbOfCompo,optr+w*nbOfCompo);
    }
  ret->copyStringInfoFrom(*this);
  return ret;
}

/*!
 * This method is a generalization of DataArrayDouble::substr method because a not contigous range can be specified here.
 * This method is equavalent to DataArrayDouble::renumberAndReduce except that convention in input is new2old and \b not old2new.
 */
DataArrayDouble *DataArrayDouble::selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const
{
  DataArrayDouble *ret=DataArrayDouble::New();
  int nbComp=getNumberOfComponents();
  ret->alloc(std::distance(new2OldBg,new2OldEnd),nbComp);
  ret->copyStringInfoFrom(*this);
  double *pt=ret->getPointer();
  const double *srcPt=getConstPointer();
  int i=0;
  for(const int *w=new2OldBg;w!=new2OldEnd;w++,i++)
    std::copy(srcPt+(*w)*nbComp,srcPt+((*w)+1)*nbComp,pt+i*nbComp);
  ret->copyStringInfoFrom(*this);
  return ret;
}

/*!
 * This method is equivalent to DataArrayInt::selectByTupleId except that an analyze to the content of input range to check that it will not lead to memory corruption !
 */
DataArrayDouble *DataArrayDouble::selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
  int nbComp=getNumberOfComponents();
  int oldNbOfTuples=getNumberOfTuples();
  ret->alloc(std::distance(new2OldBg,new2OldEnd),nbComp);
  ret->copyStringInfoFrom(*this);
  double *pt=ret->getPointer();
  const double *srcPt=getConstPointer();
  int i=0;
  for(const int *w=new2OldBg;w!=new2OldEnd;w++,i++)
    if(*w>=0 && *w<oldNbOfTuples)
      std::copy(srcPt+(*w)*nbComp,srcPt+((*w)+1)*nbComp,pt+i*nbComp);
    else
      throw INTERP_KERNEL::Exception("DataArrayInt::selectByTupleIdSafe : some ids has been detected to be out of [0,this->getNumberOfTuples) !");
  ret->copyStringInfoFrom(*this);
  ret->incrRef();
  return ret;
}

/*!
 * This methods has a similar behaviour than std::string::substr. This method returns a newly created DataArrayInt that is part of this with same number of components.
 * The intervall is specified by [tupleIdBg,tupleIdEnd) except if tupleIdEnd ==-1 in this case the [tupleIdBg,this->end()) will be kept.
 * This method check that interval is valid regarding this, if not an exception will be thrown.
 */
DataArrayDouble *DataArrayDouble::substr(int tupleIdBg, int tupleIdEnd) const throw(INTERP_KERNEL::Exception)
{
  int nbt=getNumberOfTuples();
  if(tupleIdBg<0)
    throw INTERP_KERNEL::Exception("DataArrayInt::substr : The tupleIdBg parameter must be greater than 0 !");
  if(tupleIdBg>nbt)
    throw INTERP_KERNEL::Exception("DataArrayInt::substr : The tupleIdBg parameter is greater than number of tuples !");
  int trueEnd=tupleIdEnd;
  if(tupleIdEnd!=-1)
    {
      if(tupleIdEnd>nbt)
        throw INTERP_KERNEL::Exception("DataArrayInt::substr : The tupleIdBg parameter is greater or equal than number of tuples !");
    }
  else
    trueEnd=nbt;
  int nbComp=getNumberOfComponents();
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(trueEnd-tupleIdBg,nbComp);
  ret->copyStringInfoFrom(*this);
  std::copy(getConstPointer()+tupleIdBg*nbComp,getConstPointer()+trueEnd*nbComp,ret->getPointer());
  return ret;
}

/*!
 * This method builds a new instance of DataArrayDouble (to deal with) that is reduction or an extension of 'this'.
 * if 'newNbOfComp' < this->getNumberOfComponents() a reduction is done and for each tuple 'newNbOfComp' first components are kept.
 * If 'newNbOfComp' > this->getNumberOfComponents() an extension is done, and for each components i such that i > getNumberOfComponents() 'dftValue' parameter is taken.
 */
DataArrayDouble *DataArrayDouble::changeNbOfComponents(int newNbOfComp, double dftValue) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(getNumberOfTuples(),newNbOfComp);
  const double *oldc=getConstPointer();
  double *nc=ret->getPointer();
  int nbOfTuples=getNumberOfTuples();
  int oldNbOfComp=getNumberOfComponents();
  int dim=std::min(oldNbOfComp,newNbOfComp);
  for(int i=0;i<nbOfTuples;i++)
    {
      int j=0;
      for(;j<dim;j++)
        nc[newNbOfComp*i+j]=oldc[i*oldNbOfComp+j];
      for(;j<newNbOfComp;j++)
        nc[newNbOfComp*i+j]=dftValue;
    }
  ret->setName(getName().c_str());
  for(int i=0;i<dim;i++)
    ret->setInfoOnComponent(i,getInfoOnComponent(i).c_str());
  ret->setName(getName().c_str());
  return ret;
}

/*!
 * Contrary to DataArrayDouble::changeNbOfComponents method this method is \b not const. The content 
 * This method \b do \b not change the content of data but changes the splitting of this data seen by the caller.
 * This method makes the assumption that 'this' is already allocated. If not an exception will be thrown.
 * This method checks that getNbOfElems()%newNbOfCompo==0. If not an exception will be throw !
 * This method erases all components info set before call !
 */
void DataArrayDouble::rearrange(int newNbOfCompo) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfElems=getNbOfElems();
  if(nbOfElems%newNbOfCompo!=0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::rearrange : nbOfElems%newNbOfCompo!=0 !");
  _nb_of_tuples=nbOfElems/newNbOfCompo;
  _info_on_compo.clear();
  _info_on_compo.resize(newNbOfCompo);
  declareAsNew();
}

DataArrayDouble *DataArrayDouble::keepSelectedComponents(const std::vector<int>& compoIds) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret(DataArrayDouble::New());
  int newNbOfCompo=compoIds.size();
  int oldNbOfCompo=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  ret->alloc(nbOfTuples,newNbOfCompo);
  ret->copyPartOfStringInfoFrom(*this,compoIds);
  const double *oldc=getConstPointer();
  double *nc=ret->getPointer();
  for(int i=0;i<nbOfTuples;i++)
    for(int j=0;j<newNbOfCompo;j++,nc++)
      *nc=oldc[i*oldNbOfCompo+compoIds[j]];
  ret->incrRef();
  return ret;
}

/*!
 * This method melds the components of 'this' with components of 'other'.
 * After this call in case of success, 'this' will contain a number of components equal to the sum of 'this'
 * before the call and the number of components of 'other'.
 * This method expects that 'this' and 'other' have exactly the same number of tuples. If not an exception is thrown.
 */
void DataArrayDouble::meldWith(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  other->checkAllocated();
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples!=other->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("DataArrayDouble::meldWith : mismatch of number of tuples !");
  int nbOfComp1=getNumberOfComponents();
  int nbOfComp2=other->getNumberOfComponents();
  double *newArr=new double[nbOfTuples*(nbOfComp1+nbOfComp2)];
  double *w=newArr;
  const double *inp1=getConstPointer();
  const double *inp2=other->getConstPointer();
  for(int i=0;i<nbOfTuples;i++,inp1+=nbOfComp1,inp2+=nbOfComp2)
    {
      w=std::copy(inp1,inp1+nbOfComp1,w);
      w=std::copy(inp2,inp2+nbOfComp2,w);
    }
  useArray(newArr,true,CPP_DEALLOC,nbOfTuples,nbOfComp1+nbOfComp2);
  std::vector<int> compIds(nbOfComp2);
  for(int i=0;i<nbOfComp2;i++)
    compIds[i]=nbOfComp1+i;
  copyPartOfStringInfoFrom2(compIds,*other);
}

void DataArrayDouble::setSelectedComponents(const DataArrayDouble *a, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception)
{
  copyPartOfStringInfoFrom2(compoIds,*a);
  int partOfCompoSz=compoIds.size();
  int nbOfCompo=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  const double *ac=a->getConstPointer();
  double *nc=getPointer();
  for(int i=0;i<nbOfTuples;i++)
    for(int j=0;j<partOfCompoSz;j++,ac++)
      nc[nbOfCompo*i+compoIds[j]]=*ac;
}

void DataArrayDouble::SetArrayIn(DataArrayDouble *newArray, DataArrayDouble* &arrayToSet)
{
  if(newArray!=arrayToSet)
    {
      if(arrayToSet)
        arrayToSet->decrRef();
      arrayToSet=newArray;
      if(arrayToSet)
        arrayToSet->incrRef();
    }
}

void DataArrayDouble::useArray(const double *array, bool ownership,  DeallocType type, int nbOfTuple, int nbOfCompo)
{
  _nb_of_tuples=nbOfTuple;
  _info_on_compo.resize(nbOfCompo);
  _mem.useArray(array,ownership,type,nbOfTuple*nbOfCompo);
  declareAsNew();
}

void DataArrayDouble::checkNoNullValues() const throw(INTERP_KERNEL::Exception)
{
  const double *tmp=getConstPointer();
  int nbOfElems=getNbOfElems();
  const double *where=std::find(tmp,tmp+nbOfElems,0.);
  if(where!=tmp+nbOfElems)
    throw INTERP_KERNEL::Exception("A value 0.0 have been detected !");
}

double DataArrayDouble::getMaxValue(int& tupleId) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getMaxValue : must be applied on DataArrayDouble with only one component !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<=0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getMaxValue : array exists but number of tuples must be > 0 !");
  const double *vals=getConstPointer();
  const double *loc=std::max_element(vals,vals+nbOfTuples);
  tupleId=std::distance(vals,loc);
  return *loc;
}

double DataArrayDouble::getMaxValue2(DataArrayInt*& tupleIds) const throw(INTERP_KERNEL::Exception)
{
  int tmp;
  tupleIds=0;
  double ret=getMaxValue(tmp);
  tupleIds=getIdsInRange(ret,ret);
  return ret;
}

double DataArrayDouble::getMinValue(int& tupleId) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getMinValue : must be applied on DataArrayDouble with only one component !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<=0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getMinValue : array exists but number of tuples must be > 0 !");
  const double *vals=getConstPointer();
  const double *loc=std::min_element(vals,vals+nbOfTuples);
  tupleId=std::distance(vals,loc);
  return *loc;
}

double DataArrayDouble::getMinValue2(DataArrayInt*& tupleIds) const throw(INTERP_KERNEL::Exception)
{
  int tmp;
  tupleIds=0;
  double ret=getMinValue(tmp);
  tupleIds=getIdsInRange(ret,ret);
  return ret;
}

double DataArrayDouble::getAverageValue() const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getAverageValue : must be applied on DataArrayDouble with only one component !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<=0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getAverageValue : array exists but number of tuples must be > 0 !");
  const double *vals=getConstPointer();
  double ret=std::accumulate(vals,vals+nbOfTuples,0.);
  return ret/nbOfTuples;
}

void DataArrayDouble::accumulate(double *res) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  const double *ptr=getConstPointer();
  int nbTuple=getNumberOfTuples();
  int nbComps=getNumberOfComponents();
  std::fill(res,res+nbComps,0.);
  for(int i=0;i<nbTuple;i++)
    std::transform(ptr+i*nbComps,ptr+(i+1)*nbComps,res,res,std::plus<double>());
}

double DataArrayDouble::accumulate(int compId) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  const double *ptr=getConstPointer();
  int nbTuple=getNumberOfTuples();
  int nbComps=getNumberOfComponents();
  if(compId>=nbComps)
    throw INTERP_KERNEL::Exception("DataArrayDouble::accumulate : Invalid compId specified : No such nb of components !");
  double ret=0.;
  for(int i=0;i<nbTuple;i++)
    ret+=ptr[i*nbComps+compId];
  return ret;
}

DataArrayDouble *DataArrayDouble::fromPolarToCart() const throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromPolarToCart : must be an array with exactly 2 components !");
  int nbOfTuple=getNumberOfTuples();
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbOfTuple,2);
  double *w=ret->getPointer();
  const double *wIn=getConstPointer();
  for(int i=0;i<nbOfTuple;i++,w+=2,wIn+=2)
    {
      w[0]=wIn[0]*cos(wIn[1]);
      w[1]=wIn[0]*sin(wIn[1]);
    }
  return ret;
}

DataArrayDouble *DataArrayDouble::fromCylToCart() const throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=3)
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromCylToCart : must be an array with exactly 3 components !");
  int nbOfTuple=getNumberOfTuples();
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(getNumberOfTuples(),3);
  double *w=ret->getPointer();
  const double *wIn=getConstPointer();
  for(int i=0;i<nbOfTuple;i++,w+=3,wIn+=3)
    {
      w[0]=wIn[0]*cos(wIn[1]);
      w[1]=wIn[0]*sin(wIn[1]);
      w[2]=wIn[2];
    }
  ret->setInfoOnComponent(2,getInfoOnComponent(2).c_str());
  return ret;
}

DataArrayDouble *DataArrayDouble::fromSpherToCart() const throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=3)
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromSpherToCart : must be an array with exactly 3 components !");
  int nbOfTuple=getNumberOfTuples();
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(getNumberOfTuples(),3);
  double *w=ret->getPointer();
  const double *wIn=getConstPointer();
  for(int i=0;i<nbOfTuple;i++,w+=3,wIn+=3)
    {
      w[0]=wIn[0]*cos(wIn[2])*sin(wIn[1]);
      w[1]=wIn[0]*sin(wIn[2])*sin(wIn[1]);
      w[2]=wIn[0]*cos(wIn[1]);
    }
  return ret;
}

DataArrayDouble *DataArrayDouble::doublyContractedProduct() const throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=6)
    throw INTERP_KERNEL::Exception("DataArrayDouble::doublyContractedProduct : must be an array with exactly 6 components !");
  DataArrayDouble *ret=DataArrayDouble::New();
  int nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,1);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  for(int i=0;i<nbOfTuple;i++,dest++,src+=6)
    *dest=src[0]*src[0]+src[1]*src[1]+src[2]*src[2]+2.*src[3]*src[3]+2.*src[4]*src[4]+2.*src[5]*src[5];
  return ret;
}

DataArrayDouble *DataArrayDouble::determinant() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  DataArrayDouble *ret=DataArrayDouble::New();
  int nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,1);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  switch(getNumberOfComponents())
    {
    case 6:
      for(int i=0;i<nbOfTuple;i++,dest++,src+=6)
        *dest=src[0]*src[1]*src[2]+2.*src[4]*src[5]*src[3]-src[0]*src[4]*src[4]-src[2]*src[3]*src[3]-src[1]*src[5]*src[5];
        return ret;
    case 4:
      for(int i=0;i<nbOfTuple;i++,dest++,src+=4)
        *dest=src[0]*src[3]-src[1]*src[2];
      return ret;
    case 9:
      for(int i=0;i<nbOfTuple;i++,dest++,src+=9)
        *dest=src[0]*src[4]*src[8]+src[1]*src[5]*src[6]+src[2]*src[3]*src[7]-src[0]*src[5]*src[7]-src[1]*src[3]*src[8]-src[2]*src[4]*src[6];
      return ret;
    default:
      ret->decrRef();
      throw INTERP_KERNEL::Exception("DataArrayDouble::determinant : Invalid number of components ! must be in 4,6,9 !");
    }
}

DataArrayDouble *DataArrayDouble::eigenValues() const throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=6)
    throw INTERP_KERNEL::Exception("DataArrayDouble::eigenValues : must be an array with exactly 6 components !");
  DataArrayDouble *ret=DataArrayDouble::New();
  int nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,3);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  for(int i=0;i<nbOfTuple;i++,dest+=3,src+=6)
    INTERP_KERNEL::computeEigenValues6(src,dest);
  return ret;
}

DataArrayDouble *DataArrayDouble::eigenVectors() const throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=6)
    throw INTERP_KERNEL::Exception("DataArrayDouble::eigenVectors : must be an array with exactly 6 components !");
  DataArrayDouble *ret=DataArrayDouble::New();
  int nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,9);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  for(int i=0;i<nbOfTuple;i++,src+=6)
    {
      double tmp[3];
      INTERP_KERNEL::computeEigenValues6(src,tmp);
      for(int j=0;j<3;j++,dest+=3)
        INTERP_KERNEL::computeEigenVectorForEigenValue6(src,tmp[j],1e-12,dest);
    }
  return ret;
}

DataArrayDouble *DataArrayDouble::inverse() const throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=6 && nbOfComp!=9 && nbOfComp!=4)
    throw INTERP_KERNEL::Exception("DataArrayDouble::inversion : must be an array with 4,6 or 9 components !");
  DataArrayDouble *ret=DataArrayDouble::New();
  int nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,nbOfComp);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
if(nbOfComp==6)
    for(int i=0;i<nbOfTuple;i++,dest+=6,src+=6)
      {
        double det=src[0]*src[1]*src[2]+2.*src[4]*src[5]*src[3]-src[0]*src[4]*src[4]-src[2]*src[3]*src[3]-src[1]*src[5]*src[5];
        dest[0]=(src[1]*src[2]-src[4]*src[4])/det;
        dest[1]=(src[0]*src[2]-src[5]*src[5])/det;
        dest[2]=(src[0]*src[1]-src[3]*src[3])/det;
        dest[3]=(src[5]*src[4]-src[3]*src[2])/det;
        dest[4]=(src[5]*src[3]-src[0]*src[4])/det;
        dest[5]=(src[3]*src[4]-src[1]*src[5])/det;
      }
  else if(nbOfComp==4)
    for(int i=0;i<nbOfTuple;i++,dest+=4,src+=4)
      {
        double det=src[0]*src[3]-src[1]*src[2];
        dest[0]=src[3]/det;
        dest[1]=-src[1]/det;
        dest[2]=-src[2]/det;
        dest[3]=src[0]/det;
      }
  else
    for(int i=0;i<nbOfTuple;i++,dest+=9,src+=9)
      {
        double det=src[0]*src[4]*src[8]+src[1]*src[5]*src[6]+src[2]*src[3]*src[7]-src[0]*src[5]*src[7]-src[1]*src[3]*src[8]-src[2]*src[4]*src[6];
        dest[0]=(src[4]*src[8]-src[7]*src[5])/det;
        dest[1]=(src[7]*src[2]-src[1]*src[8])/det;
        dest[2]=(src[1]*src[5]-src[4]*src[2])/det;
        dest[3]=(src[6]*src[5]-src[3]*src[8])/det;
        dest[4]=(src[0]*src[8]-src[6]*src[2])/det;
        dest[5]=(src[2]*src[3]-src[0]*src[5])/det;
        dest[6]=(src[3]*src[7]-src[6]*src[4])/det;
        dest[7]=(src[6]*src[1]-src[0]*src[7])/det;
        dest[8]=(src[0]*src[4]-src[1]*src[3])/det;
      }
  return ret;
}

DataArrayDouble *DataArrayDouble::trace() const throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=6 && nbOfComp!=9 && nbOfComp!=4)
    throw INTERP_KERNEL::Exception("DataArrayDouble::trace : must be an array with 4,6 or 9 components !");
  DataArrayDouble *ret=DataArrayDouble::New();
  int nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,1);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  if(nbOfComp==6)
    for(int i=0;i<nbOfTuple;i++,dest++,src+=6)
      *dest=src[0]+src[1]+src[2];
  else if(nbOfComp==4)
    for(int i=0;i<nbOfTuple;i++,dest++,src+=4)
      *dest=src[0]+src[3];
  else
    for(int i=0;i<nbOfTuple;i++,dest++,src+=9)
      *dest=src[0]+src[4]+src[8];
  return ret;
}

DataArrayDouble *DataArrayDouble::deviator() const throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=6)
    throw INTERP_KERNEL::Exception("DataArrayDouble::deviator : must be an array with exactly 6 components !");
  DataArrayDouble *ret=DataArrayDouble::New();
  int nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,6);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  for(int i=0;i<nbOfTuple;i++,dest+=6,src+=6)
    {
      double tr=(src[0]+src[1]+src[2])/3.;
      dest[0]=src[0]-tr;
      dest[1]=src[1]-tr;
      dest[2]=src[2]-tr;
      dest[3]=src[3];
      dest[4]=src[4];
      dest[5]=src[5];
    }
  return ret;
}

DataArrayDouble *DataArrayDouble::magnitude() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfComp=getNumberOfComponents();
  DataArrayDouble *ret=DataArrayDouble::New();
  int nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,1);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  for(int i=0;i<nbOfTuple;i++,dest++)
    {
      double sum=0.;
      for(int j=0;j<nbOfComp;j++,src++)
        sum+=(*src)*(*src);
      *dest=sqrt(sum);
    }
  return ret;
}

DataArrayDouble *DataArrayDouble::maxPerTuple() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfComp=getNumberOfComponents();
  DataArrayDouble *ret=DataArrayDouble::New();
  int nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,1);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  for(int i=0;i<nbOfTuple;i++,dest++,src+=nbOfComp)
    *dest=*std::max_element(src,src+nbOfComp);
  return ret;
}

void DataArrayDouble::sortPerTuple(bool asc) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  double *pt=getPointer();
  int nbOfTuple=getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  if(asc)
    for(int i=0;i<nbOfTuple;i++,pt+=nbOfComp)
      std::sort(pt,pt+nbOfComp);
  else
    for(int i=0;i<nbOfTuple;i++,pt+=nbOfComp)
      std::sort(pt,pt+nbOfComp,std::greater<double>());
  declareAsNew();
}

void DataArrayDouble::applyLin(double a, double b, int compoId) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  double *ptr=getPointer()+compoId;
  int nbOfComp=getNumberOfComponents();
  int nbOfTuple=getNumberOfTuples();
  for(int i=0;i<nbOfTuple;i++,ptr+=nbOfComp)
    *ptr=a*(*ptr)+b;
  declareAsNew();
}

DataArrayDouble *DataArrayDouble::applyFunc(int nbOfComp, FunctionToEvaluate func) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  DataArrayDouble *newArr=DataArrayDouble::New();
  int nbOfTuples=getNumberOfTuples();
  int oldNbOfComp=getNumberOfComponents();
  newArr->alloc(nbOfTuples,nbOfComp);
  const double *ptr=getConstPointer();
  double *ptrToFill=newArr->getPointer();
  for(int i=0;i<nbOfTuples;i++)
    {
      if(!func(ptr+i*oldNbOfComp,ptrToFill+i*nbOfComp))
        {
          std::ostringstream oss; oss << "For tuple # " << i << " with value (";
          std::copy(ptr+oldNbOfComp*i,ptr+oldNbOfComp*(i+1),std::ostream_iterator<double>(oss,", "));
          oss << ") : Evaluation of function failed !";
          newArr->decrRef();
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return newArr;
}

DataArrayDouble *DataArrayDouble::applyFunc(int nbOfComp, const char *func) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  std::set<std::string> vars;
  expr.getTrueSetOfVars(vars);
  int oldNbOfComp=getNumberOfComponents();
  if((int)vars.size()>oldNbOfComp)
    {
      std::ostringstream oss; oss << "The field has a " << oldNbOfComp << " components and there are ";
      oss << vars.size() << " variables : ";
      std::copy(vars.begin(),vars.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::vector<std::string> varsV(vars.begin(),vars.end());
  expr.prepareExprEvaluation(varsV);
  //
  DataArrayDouble *newArr=DataArrayDouble::New();
  int nbOfTuples=getNumberOfTuples();
  newArr->alloc(nbOfTuples,nbOfComp);
  const double *ptr=getConstPointer();
  double *ptrToFill=newArr->getPointer();
  for(int i=0;i<nbOfTuples;i++)
    {
      try
        {
          expr.evaluateExpr(nbOfComp,ptr+i*oldNbOfComp,ptrToFill+i*nbOfComp);
        }
      catch(INTERP_KERNEL::Exception& e)
        {
          std::ostringstream oss; oss << "For tuple # " << i << " with value (";
          std::copy(ptr+oldNbOfComp*i,ptr+oldNbOfComp*(i+1),std::ostream_iterator<double>(oss,", "));
          oss << ") : Evaluation of function failed !" << e.what();
          newArr->decrRef();
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return newArr;
}

DataArrayDouble *DataArrayDouble::applyFunc(const char *func) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  expr.prepareExprEvaluationVec();
  //
  DataArrayDouble *newArr=DataArrayDouble::New();
  int nbOfTuples=getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  newArr->alloc(nbOfTuples,nbOfComp);
  const double *ptr=getConstPointer();
  double *ptrToFill=newArr->getPointer();
  for(int i=0;i<nbOfTuples;i++)
    {
      try
        {
          expr.evaluateExpr(nbOfComp,ptr+i*nbOfComp,ptrToFill+i*nbOfComp);
        }
      catch(INTERP_KERNEL::Exception& e)
        {
          std::ostringstream oss; oss << "For tuple # " << i << " with value (";
          std::copy(ptr+nbOfComp*i,ptr+nbOfComp*(i+1),std::ostream_iterator<double>(oss,", "));
          oss << ") : Evaluation of function failed ! " << e.what();
          newArr->decrRef();
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return newArr;
}

void DataArrayDouble::applyFuncFast32(const char *func) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  char *funcStr=expr.compileX86();
  MYFUNCPTR funcPtr=(MYFUNCPTR)funcStr;//he he...
  //
  double *ptr=getPointer();
  int nbOfComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  int nbOfElems=nbOfTuples*nbOfComp;
  for(int i=0;i<nbOfElems;i++,ptr++)
    *ptr=funcPtr(*ptr);
  declareAsNew();
}

void DataArrayDouble::applyFuncFast64(const char *func) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  char *funcStr=expr.compileX86_64();
  MYFUNCPTR funcPtr=(MYFUNCPTR)funcStr;//he he...
  //
  double *ptr=getPointer();
  int nbOfComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  int nbOfElems=nbOfTuples*nbOfComp;
  for(int i=0;i<nbOfElems;i++,ptr++)
    *ptr=funcPtr(*ptr);
  declareAsNew();
}

DataArrayInt *DataArrayDouble::getIdsInRange(double vmin, double vmax) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getIdsInRange : this must have exactly one component !");
  const double *cptr=getConstPointer();
  std::vector<int> res;
  int nbOfTuples=getNumberOfTuples();
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(*cptr>=vmin && *cptr<=vmax)
      res.push_back(i);
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(res.size(),1);
  std::copy(res.begin(),res.end(),ret->getPointer());
  return ret;
}

DataArrayDouble *DataArrayDouble::Aggregate(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayDouble *> tmp(2);
  tmp[0]=a1; tmp[1]=a2;
  return Aggregate(tmp);
}

DataArrayDouble *DataArrayDouble::Aggregate(const std::vector<const DataArrayDouble *>& a) throw(INTERP_KERNEL::Exception)
{
  if(a.empty())
    throw INTERP_KERNEL::Exception("DataArrayDouble::Aggregate : input list must be NON EMPTY !");
  std::vector<const DataArrayDouble *>::const_iterator it=a.begin();
  int nbOfComp=(*it)->getNumberOfComponents();
  int nbt=(*it++)->getNumberOfTuples();
  for(int i=1;it!=a.end();it++,i++)
    {
      if((*it)->getNumberOfComponents()!=nbOfComp)
        throw INTERP_KERNEL::Exception("DataArrayDouble::Aggregate : Nb of components mismatch for array aggregation !");
      nbt+=(*it)->getNumberOfTuples();
    }
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbt,nbOfComp);
  double *pt=ret->getPointer();
  for(it=a.begin();it!=a.end();it++)
    pt=std::copy((*it)->getConstPointer(),(*it)->getConstPointer()+(*it)->getNbOfElems(),pt);
  ret->copyStringInfoFrom(*(a[0]));
  return ret;
}

DataArrayDouble *DataArrayDouble::Meld(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayDouble *> arr(2);
  arr[0]=a1; arr[1]=a2;
  return Meld(arr);
}

DataArrayDouble *DataArrayDouble::Meld(const std::vector<const DataArrayDouble *>& a) throw(INTERP_KERNEL::Exception)
{
  if(a.empty())
    throw INTERP_KERNEL::Exception("DataArrayDouble::Meld : array must be NON empty !");
  std::vector<const DataArrayDouble *>::const_iterator it;
  for(it=a.begin();it!=a.end();it++)
    (*it)->checkAllocated();
  it=a.begin();
  int nbOfTuples=(*it)->getNumberOfTuples();
  std::vector<int> nbc(a.size());
  std::vector<const double *> pts(a.size());
  nbc[0]=(*it)->getNumberOfComponents();
  pts[0]=(*it++)->getConstPointer();
  for(int i=1;it!=a.end();it++,i++)
    {
      if(nbOfTuples!=(*it)->getNumberOfTuples())
        throw INTERP_KERNEL::Exception("DataArrayDouble::Meld : mismatch of number of tuples !");
      nbc[i]=(*it)->getNumberOfComponents();
      pts[i]=(*it)->getConstPointer();
    }
  int totalNbOfComp=std::accumulate(nbc.begin(),nbc.end(),0);
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbOfTuples,totalNbOfComp);
  double *retPtr=ret->getPointer();
  for(int i=0;i<nbOfTuples;i++)
    for(int j=0;j<(int)a.size();j++)
      {
        retPtr=std::copy(pts[j],pts[j]+nbc[j],retPtr);
        pts[j]+=nbc[j];
      }
  int k=0;
  for(int i=0;i<(int)a.size();i++)
    for(int j=0;j<nbc[i];j++,k++)
      ret->setInfoOnComponent(k,a[i]->getInfoOnComponent(j).c_str());
  return ret;
}

DataArrayDouble *DataArrayDouble::Dot(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  a1->checkAllocated();
  a2->checkAllocated();
  int nbOfComp=a1->getNumberOfComponents();
  if(nbOfComp!=a2->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array Dot !");
  int nbOfTuple=a1->getNumberOfTuples();
  if(nbOfTuple!=a2->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array Dot !");
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbOfTuple,1);
  double *retPtr=ret->getPointer();
  const double *a1Ptr=a1->getConstPointer();
  const double *a2Ptr=a2->getConstPointer();
  for(int i=0;i<nbOfTuple;i++)
    {
      double sum=0.;
      for(int j=0;j<nbOfComp;j++)
        sum+=a1Ptr[i*nbOfComp+j]*a2Ptr[i*nbOfComp+j];
      retPtr[i]=sum;
    }
  ret->setInfoOnComponent(0,a1->getInfoOnComponent(0).c_str());
  ret->setName(a1->getName().c_str());
  return ret;
}

DataArrayDouble *DataArrayDouble::CrossProduct(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=a1->getNumberOfComponents();
  if(nbOfComp!=a2->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array crossProduct !");
  if(nbOfComp!=3)
    throw INTERP_KERNEL::Exception("Nb of components must be equal to 3 for array crossProduct !");
  int nbOfTuple=a1->getNumberOfTuples();
  if(nbOfTuple!=a2->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array crossProduct !");
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbOfTuple,3);
  double *retPtr=ret->getPointer();
  const double *a1Ptr=a1->getConstPointer();
  const double *a2Ptr=a2->getConstPointer();
  for(int i=0;i<nbOfTuple;i++)
    {
      retPtr[3*i]=a1Ptr[3*i+1]*a2Ptr[3*i+2]-a1Ptr[3*i+2]*a2Ptr[3*i+1];
      retPtr[3*i+1]=a1Ptr[3*i+2]*a2Ptr[3*i]-a1Ptr[3*i]*a2Ptr[3*i+2];
      retPtr[3*i+2]=a1Ptr[3*i]*a2Ptr[3*i+1]-a1Ptr[3*i+1]*a2Ptr[3*i];
    }
  ret->copyStringInfoFrom(*a1);
  return ret;
}

DataArrayDouble *DataArrayDouble::Max(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=a1->getNumberOfComponents();
  if(nbOfComp!=a2->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array Max !");
  int nbOfTuple=a1->getNumberOfTuples();
  if(nbOfTuple!=a2->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array Max !");
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbOfTuple,nbOfComp);
  double *retPtr=ret->getPointer();
  const double *a1Ptr=a1->getConstPointer();
  const double *a2Ptr=a2->getConstPointer();
  int nbElem=nbOfTuple*nbOfComp;
  for(int i=0;i<nbElem;i++)
    retPtr[i]=std::max(a1Ptr[i],a2Ptr[i]);
  ret->copyStringInfoFrom(*a1);
  return ret;
}

DataArrayDouble *DataArrayDouble::Min(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=a1->getNumberOfComponents();
  if(nbOfComp!=a2->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array min !");
  int nbOfTuple=a1->getNumberOfTuples();
  if(nbOfTuple!=a2->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array min !");
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbOfTuple,nbOfComp);
  double *retPtr=ret->getPointer();
  const double *a1Ptr=a1->getConstPointer();
  const double *a2Ptr=a2->getConstPointer();
  int nbElem=nbOfTuple*nbOfComp;
  for(int i=0;i<nbElem;i++)
    retPtr[i]=std::min(a1Ptr[i],a2Ptr[i]);
  ret->copyStringInfoFrom(*a1);
  return ret;
}

DataArrayDouble *DataArrayDouble::Add(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=a1->getNumberOfComponents();
  if(nbOfComp!=a2->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array add !");
  int nbOfTuple=a1->getNumberOfTuples();
  if(nbOfTuple!=a2->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array add !");
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbOfTuple,nbOfComp);
  std::transform(a1->getConstPointer(),a1->getConstPointer()+nbOfTuple*nbOfComp,a2->getConstPointer(),ret->getPointer(),std::plus<double>());
  ret->copyStringInfoFrom(*a1);
  return ret;
}

void DataArrayDouble::addEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=other->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array add !");
  int nbOfTuple=getNumberOfTuples();
  if(nbOfTuple!=other->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array add !");
  std::transform(getConstPointer(),getConstPointer()+nbOfTuple*nbOfComp,other->getConstPointer(),getPointer(),std::plus<double>());
  declareAsNew();
}

DataArrayDouble *DataArrayDouble::Substract(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=a1->getNumberOfComponents();
  if(nbOfComp!=a2->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array Substract !");
  int nbOfTuple=a1->getNumberOfTuples();
  if(nbOfTuple!=a2->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array Substract !");
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbOfTuple,nbOfComp);
  std::transform(a1->getConstPointer(),a1->getConstPointer()+nbOfTuple*nbOfComp,a2->getConstPointer(),ret->getPointer(),std::minus<double>());
  ret->copyStringInfoFrom(*a1);
  return ret;
}

void DataArrayDouble::substractEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception)
{
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=other->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array substract !");
  int nbOfTuple=getNumberOfTuples();
  if(nbOfTuple!=other->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array substract !");
  std::transform(getConstPointer(),getConstPointer()+nbOfTuple*nbOfComp,other->getConstPointer(),getPointer(),std::minus<double>());
  declareAsNew();
}

DataArrayDouble *DataArrayDouble::Multiply(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  int nbOfTuple=a1->getNumberOfTuples();
  int nbOfTuple2=a2->getNumberOfTuples();
  int nbOfComp=a1->getNumberOfComponents();
  int nbOfComp2=a2->getNumberOfComponents();
  if(nbOfTuple!=nbOfTuple2)
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array multiply !");
  DataArrayDouble *ret=0;
  if(nbOfComp==nbOfComp2)
    {
      ret=DataArrayDouble::New();
      ret->alloc(nbOfTuple,nbOfComp);
      std::transform(a1->getConstPointer(),a1->getConstPointer()+nbOfTuple*nbOfComp,a2->getConstPointer(),ret->getPointer(),std::multiplies<double>());
      ret->copyStringInfoFrom(*a1);
    }
  else
    {
      int nbOfCompMin,nbOfCompMax;
      const DataArrayDouble *aMin, *aMax;
      if(nbOfComp>nbOfComp2)
        {
          nbOfCompMin=nbOfComp2; nbOfCompMax=nbOfComp;
          aMin=a2; aMax=a1;
        }
      else
        {
          nbOfCompMin=nbOfComp; nbOfCompMax=nbOfComp2;
          aMin=a1; aMax=a2;
        }
      if(nbOfCompMin==1)
        {
          ret=DataArrayDouble::New();
          ret->alloc(nbOfTuple,nbOfCompMax);
          const double *aMinPtr=aMin->getConstPointer();
          const double *aMaxPtr=aMax->getConstPointer();
          double *res=ret->getPointer();
          for(int i=0;i<nbOfTuple;i++)
            res=std::transform(aMaxPtr+i*nbOfCompMax,aMaxPtr+(i+1)*nbOfCompMax,res,std::bind2nd(std::multiplies<double>(),aMinPtr[i]));
          ret->copyStringInfoFrom(*aMax);
        }
      else
        throw INTERP_KERNEL::Exception("Nb of components mismatch for array multiply !");
    }
  return ret;
}

void DataArrayDouble::multiplyEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception)
{
  int nbOfTuple=getNumberOfTuples();
  int nbOfTuple2=other->getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  int nbOfComp2=other->getNumberOfComponents();
  if(nbOfTuple!=nbOfTuple2)
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array multiplyEqual !");
  if(nbOfComp==nbOfComp2)
    {
      std::transform(getConstPointer(),getConstPointer()+nbOfTuple*nbOfComp,other->getConstPointer(),getPointer(),std::multiplies<double>());
    }
  else
    {
      if(nbOfComp2==1)
        {
          const double *ptr=other->getConstPointer();
          double *myPtr=getPointer();
          for(int i=0;i<nbOfTuple;i++)
            myPtr=std::transform(myPtr,myPtr+nbOfComp,myPtr,std::bind2nd(std::multiplies<double>(),ptr[i]));
        }
      else
        throw INTERP_KERNEL::Exception("Nb of components mismatch for array multiplyEqual !");
    }
  declareAsNew();
}

DataArrayDouble *DataArrayDouble::Divide(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  int nbOfTuple=a1->getNumberOfTuples();
  int nbOfTuple2=a2->getNumberOfTuples();
  int nbOfComp=a1->getNumberOfComponents();
  int nbOfComp2=a2->getNumberOfComponents();
  if(nbOfTuple!=nbOfTuple2)
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array divide !");
  DataArrayDouble *ret=0;
  if(nbOfComp==nbOfComp2)
    {
      ret=DataArrayDouble::New();
      ret->alloc(nbOfTuple,nbOfComp);
      std::transform(a1->getConstPointer(),a1->getConstPointer()+nbOfTuple*nbOfComp,a2->getConstPointer(),ret->getPointer(),std::divides<double>());
      ret->copyStringInfoFrom(*a1);
    }
  else
    {
      if(nbOfComp2==1)
        {
          ret=DataArrayDouble::New();
          ret->alloc(nbOfTuple,nbOfComp);
          const double *a2Ptr=a2->getConstPointer();
          const double *a1Ptr=a1->getConstPointer();
          double *res=ret->getPointer();
          for(int i=0;i<nbOfTuple;i++)
            res=std::transform(a1Ptr+i*nbOfComp,a1Ptr+(i+1)*nbOfComp,res,std::bind2nd(std::divides<double>(),a2Ptr[i]));
          ret->copyStringInfoFrom(*a1);
        }
      else
        throw INTERP_KERNEL::Exception("Nb of components mismatch for array divide !");
    }
  return ret;
}

void DataArrayDouble::divideEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception)
{
  int nbOfTuple=getNumberOfTuples();
  int nbOfTuple2=other->getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  int nbOfComp2=other->getNumberOfComponents();
  if(nbOfTuple!=nbOfTuple2)
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array divideEqual !");
  if(nbOfComp==nbOfComp2)
    {
      std::transform(getConstPointer(),getConstPointer()+nbOfTuple*nbOfComp,other->getConstPointer(),getPointer(),std::divides<double>());
    }
  else
    {
      if(nbOfComp2==1)
        {
          const double *ptr=other->getConstPointer();
          double *myPtr=getPointer();
          for(int i=0;i<nbOfTuple;i++)
            myPtr=std::transform(myPtr,myPtr+nbOfComp,myPtr,std::bind2nd(std::divides<double>(),ptr[i]));
        }
      else
        throw INTERP_KERNEL::Exception("Nb of components mismatch for array divideEqual !");
    }
  declareAsNew();
}

DataArrayInt *DataArrayInt::New()
{
  return new DataArrayInt;
}

bool DataArrayInt::isAllocated() const
{
  return getConstPointer()!=0;
}

void DataArrayInt::checkAllocated() const throw(INTERP_KERNEL::Exception)
{
  if(!isAllocated())
    throw INTERP_KERNEL::Exception("DataArrayInt::checkAllocated : Array is defined but not allocated ! Call alloc or setValues method first !");
}

DataArrayInt *DataArrayInt::deepCpy() const
{
  return new DataArrayInt(*this);
}

DataArrayInt *DataArrayInt::performCpy(bool dCpy) const
{
  if(dCpy)
    return deepCpy();
  else
    {
      incrRef();
      return const_cast<DataArrayInt *>(this);
    }
}

void DataArrayInt::alloc(int nbOfTuple, int nbOfCompo)
{
  _nb_of_tuples=nbOfTuple;
  _info_on_compo.resize(nbOfCompo);
  _mem.alloc(nbOfCompo*_nb_of_tuples);
  declareAsNew();
}

void DataArrayInt::fillWithZero()
{
  _mem.fillWithValue(0);
  declareAsNew();
}

void DataArrayInt::fillWithValue(int val)
{
  _mem.fillWithValue(val);
  declareAsNew();
}

void DataArrayInt::iota(int init) throw(INTERP_KERNEL::Exception)
{
  int *ptr=getPointer();
  if(!ptr)
    throw INTERP_KERNEL::Exception("DataArrayInt::iota : allocate first !");
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::iota : works only for arrays with only one component!");
  int ntuples=getNumberOfTuples();
  for(int i=0;i<ntuples;i++)
    ptr[i]=init+i;
  declareAsNew();
}

std::string DataArrayInt::repr() const
{
  std::ostringstream ret;
  reprStream(ret);
  return ret.str();
}

std::string DataArrayInt::reprZip() const
{
  std::ostringstream ret;
  reprZipStream(ret);
  return ret.str();
}

void DataArrayInt::reprStream(std::ostream& stream) const
{
  stream << "Name of int array : \"" << _name << "\"\n";
  reprWithoutNameStream(stream);
}

void DataArrayInt::reprZipStream(std::ostream& stream) const
{
  stream << "Name of int array : \"" << _name << "\"\n";
  reprZipWithoutNameStream(stream);
}

void DataArrayInt::reprWithoutNameStream(std::ostream& stream) const
{
  DataArray::reprWithoutNameStream(stream);
  _mem.repr(getNumberOfComponents(),stream);
}

void DataArrayInt::reprZipWithoutNameStream(std::ostream& stream) const
{
  DataArray::reprWithoutNameStream(stream);
  _mem.reprZip(getNumberOfComponents(),stream);
}

void DataArrayInt::transformWithIndArr(const int *indArr)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("Call transformWithIndArr method on DataArrayInt with only one component !");
  int nbOfTuples=getNumberOfTuples();
  int *pt=getPointer();
  for(int i=0;i<nbOfTuples;i++)
    pt[i]=indArr[pt[i]];
}

/*!
 * This method invert array 'di' that is a conversion map from Old to New numbering to New to Old numbering.
 */
DataArrayInt *DataArrayInt::invertArrayO2N2N2O(int newNbOfElem) const
{
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(newNbOfElem,1);
  int nbOfOldNodes=getNumberOfTuples();
  const int *old2New=getConstPointer();
  int *pt=ret->getPointer();
  for(int i=0;i!=nbOfOldNodes;i++)
    if(old2New[i]!=-1)
      pt[old2New[i]]=i;
  return ret;
}

/*!
 * This method invert array 'di' that is a conversion map from New to old numbering to Old to New numbering.
 */
DataArrayInt *DataArrayInt::invertArrayN2O2O2N(int oldNbOfElem) const
{
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(oldNbOfElem,1);
  const int *new2Old=getConstPointer();
  int *pt=ret->getPointer();
  std::fill(pt,pt+oldNbOfElem,-1);
  int nbOfNewElems=getNumberOfTuples();
  for(int i=0;i<nbOfNewElems;i++)
    pt[new2Old[i]]=i;
  return ret;
}

bool DataArrayInt::isEqual(const DataArrayInt& other) const
{
  if(!areInfoEquals(other))
    return false;
  return _mem.isEqual(other._mem,0);
}

bool DataArrayInt::isEqualWithoutConsideringStr(const DataArrayInt& other) const
{
  return _mem.isEqual(other._mem,0);
}

bool DataArrayInt::isEqualWithoutConsideringStrAndOrder(const DataArrayInt& other) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> a=deepCpy();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> b=other.deepCpy();
  a->sort();
  b->sort();
  return a->isEqualWithoutConsideringStr(*b);
}

void DataArrayInt::sort() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::sort : only supported with 'this' array with ONE component !");
  _mem.sort();
}

/*!
 * This method expects that 'this' and 'other' have the same number of tuples and exactly one component both. If not an exception will be thrown.
 * This method retrieves a newly created array with same number of tuples than 'this' and 'other' with one component.
 * The returned array 'ret' contains the correspondance from 'this' to 'other' that is to say for every i so that 0<=i<getNumberOfTuples()
 * other.getIJ(i,0)==this->getIJ(ret->getIJ(i),0)
 * If such permutation is not possible because it exists some elements in 'other' not in 'this', an exception will be thrown.
 */
DataArrayInt *DataArrayInt::buildPermutationArr(const DataArrayInt& other) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1 || other.getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::buildPermutationArr : 'this' and 'other' have to have exactly ONE component !");
  int nbTuple=getNumberOfTuples();
  if(nbTuple!=other.getNumberOfTuples())
    throw INTERP_KERNEL::Exception("DataArrayInt::buildPermutationArr : 'this' and 'other' must have the same number of tuple !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbTuple,1);
  ret->fillWithValue(-1);
  const int *pt=getConstPointer();
  std::map<int,int> mm;
  for(int i=0;i<nbTuple;i++)
    mm[pt[i]]=i;
  pt=other.getConstPointer();
  int *retToFill=ret->getPointer();
  for(int i=0;i<nbTuple;i++)
    {
      std::map<int,int>::const_iterator it=mm.find(pt[i]);
      if(it==mm.end())
        {
          std::ostringstream oss; oss << "DataArrayInt::buildPermutationArr : Arrays mismatch : element (" << pt[i] << ") in 'other' not findable in 'this' !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      retToFill[i]=(*it).second;
    }
  ret->incrRef();
  return ret;
}

void DataArrayInt::useArray(const int *array, bool ownership,  DeallocType type, int nbOfTuple, int nbOfCompo)
{
  _nb_of_tuples=nbOfTuple;
  _info_on_compo.resize(nbOfCompo);
  _mem.useArray(array,ownership,type,nbOfTuple*nbOfCompo);
  declareAsNew();
}

DataArrayInt *DataArrayInt::fromNoInterlace() const throw(INTERP_KERNEL::Exception)
{
  if(_mem.isNull())
    throw INTERP_KERNEL::Exception("DataArrayInt::fromNoInterlace : Not defined array !");
  int *tab=_mem.fromNoInterlace(getNumberOfComponents());
  DataArrayInt *ret=DataArrayInt::New();
  ret->useArray(tab,true,CPP_DEALLOC,getNumberOfTuples(),getNumberOfComponents());
  return ret;
}

DataArrayInt *DataArrayInt::toNoInterlace() const throw(INTERP_KERNEL::Exception)
{
  if(_mem.isNull())
    throw INTERP_KERNEL::Exception("DataArrayInt::toNoInterlace : Not defined array !");
  int *tab=_mem.toNoInterlace(getNumberOfComponents());
  DataArrayInt *ret=DataArrayInt::New();
  ret->useArray(tab,true,CPP_DEALLOC,getNumberOfTuples(),getNumberOfComponents());
  return ret;
}

void DataArrayInt::renumberInPlace(const int *old2New)
{
  int nbTuples=getNumberOfTuples();
  int nbOfCompo=getNumberOfComponents();
  int *tmp=new int[nbTuples*nbOfCompo];
  const int *iptr=getConstPointer();
  for(int i=0;i<nbTuples;i++)
    std::copy(iptr+nbOfCompo*i,iptr+nbOfCompo*(i+1),tmp+nbOfCompo*old2New[i]);
  std::copy(tmp,tmp+nbTuples*nbOfCompo,getPointer());
  delete [] tmp;
  declareAsNew();
}

void DataArrayInt::renumberInPlaceR(const int *new2Old)
{
  int nbTuples=getNumberOfTuples();
  int nbOfCompo=getNumberOfComponents();
  int *tmp=new int[nbTuples*nbOfCompo];
  const int *iptr=getConstPointer();
  for(int i=0;i<nbTuples;i++)
    std::copy(iptr+nbOfCompo*new2Old[i],iptr+nbOfCompo*(new2Old[i]+1),tmp+nbOfCompo*i);
  std::copy(tmp,tmp+nbTuples*nbOfCompo,getPointer());
  delete [] tmp;
  declareAsNew();
}

DataArrayInt *DataArrayInt::renumber(const int *old2New) const
{
  int nbTuples=getNumberOfTuples();
  int nbOfCompo=getNumberOfComponents();
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbTuples,nbOfCompo);
  ret->copyStringInfoFrom(*this);
  const int *iptr=getConstPointer();
  int *optr=ret->getPointer();
  for(int i=0;i<nbTuples;i++)
    std::copy(iptr+nbOfCompo*i,iptr+nbOfCompo*(i+1),optr+nbOfCompo*old2New[i]);
  ret->copyStringInfoFrom(*this);
  return ret;
}

DataArrayInt *DataArrayInt::renumberR(const int *new2Old) const
{
  int nbTuples=getNumberOfTuples();
  int nbOfCompo=getNumberOfComponents();
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbTuples,nbOfCompo);
  ret->copyStringInfoFrom(*this);
  const int *iptr=getConstPointer();
  int *optr=ret->getPointer();
  for(int i=0;i<nbTuples;i++)
    std::copy(iptr+nbOfCompo*new2Old[i],iptr+nbOfCompo*(new2Old[i]+1),optr+nbOfCompo*i);
  ret->copyStringInfoFrom(*this);
  return ret;
}

/*!
 * Idem DataArrayDouble::renumber method except that the number of tuples is reduced.
 * That is to say that it is expected that newNbOfTuple<this->getNumberOfTuples().
 * ['old2New','old2New'+getNumberOfTuples()) defines a range containing old to new array. For every negative value in ['old2NewBg','old2New'getNumberOfTuples()) the corresponding tuple is
 * omitted.
 */
DataArrayInt *DataArrayInt::renumberAndReduce(const int *old2New, int newNbOfTuple) const
{
  int nbTuples=getNumberOfTuples();
  int nbOfCompo=getNumberOfComponents();
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(newNbOfTuple,nbOfCompo);
  const int *iptr=getConstPointer();
  int *optr=ret->getPointer();
  for(int i=0;i<nbTuples;i++)
    {
      int w=old2New[i];
      if(w>=0)
        std::copy(iptr+i*nbOfCompo,iptr+(i+1)*nbOfCompo,optr+w*nbOfCompo);
    }
  ret->copyStringInfoFrom(*this);
  return ret;
}

/*!
 * This method is a generalization of DataArrayDouble::substr method because a not contigous range can be specified here.
 * This method is equavalent to DataArrayInt::renumberAndReduce except that convention in input is new2old and \b not old2new.
 */
DataArrayInt *DataArrayInt::selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const
{
  DataArrayInt *ret=DataArrayInt::New();
  int nbComp=getNumberOfComponents();
  ret->alloc(std::distance(new2OldBg,new2OldEnd),nbComp);
  ret->copyStringInfoFrom(*this);
  int *pt=ret->getPointer();
  const int *srcPt=getConstPointer();
  int i=0;
  for(const int *w=new2OldBg;w!=new2OldEnd;w++,i++)
    std::copy(srcPt+(*w)*nbComp,srcPt+((*w)+1)*nbComp,pt+i*nbComp);
  ret->copyStringInfoFrom(*this);
  return ret;
}

/*!
 * This method is equivalent to DataArrayInt::selectByTupleId except that an analyze to the content of input range to check that it will not lead to memory corruption !
 */
DataArrayInt *DataArrayInt::selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  int nbComp=getNumberOfComponents();
  int oldNbOfTuples=getNumberOfTuples();
  ret->alloc(std::distance(new2OldBg,new2OldEnd),nbComp);
  ret->copyStringInfoFrom(*this);
  int *pt=ret->getPointer();
  const int *srcPt=getConstPointer();
  int i=0;
  for(const int *w=new2OldBg;w!=new2OldEnd;w++,i++)
    if(*w>=0 && *w<oldNbOfTuples)
      std::copy(srcPt+(*w)*nbComp,srcPt+((*w)+1)*nbComp,pt+i*nbComp);
    else
      throw INTERP_KERNEL::Exception("DataArrayInt::selectByTupleIdSafe : some ids has been detected to be out of [0,this->getNumberOfTuples) !");
  ret->copyStringInfoFrom(*this);
  ret->incrRef();
  return ret;
}

/*!
 * This method checks that 'this' is with numberofcomponents == 1 and that it is equal to
 * stdext::iota() of size getNumberOfTuples. This method is particalary usefull for DataArrayInt instances
 * that represents a renumbering array to check the real need in renumbering. 
 */
bool DataArrayInt::isIdentity() const
{
  if(getNumberOfComponents()!=1)
    return false;
  int nbOfTuples=getNumberOfTuples();
  const int *pt=getConstPointer();
  for(int i=0;i<nbOfTuples;i++,pt++)
    if(*pt!=i)
      return false;
  return true;
}

bool DataArrayInt::isUniform(int val) const
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::isUniform : must be applied on DataArrayInt with only one component !");
  int nbOfTuples=getNumberOfTuples();
  const int *w=getConstPointer();
  const int *end=w+nbOfTuples;
  for(;w!=end;w++)
    if(*w!=val)
      return false;
  return true;
}

DataArrayDouble *DataArrayInt::convertToDblArr() const
{
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(getNumberOfTuples(),getNumberOfComponents());
  int nbOfVals=getNbOfElems();
  const int *src=getConstPointer();
  double *dest=ret->getPointer();
  std::copy(src,src+nbOfVals,dest);
  ret->copyStringInfoFrom(*this);
  return ret;
}

/*!
 * This methods has a similar behaviour than std::string::substr. This method returns a newly created DataArrayInt that is part of this with same number of components.
 * The intervall is specified by [tupleIdBg,tupleIdEnd) except if tupleIdEnd ==-1 in this case the [tupleIdBg,this->end()) will be kept.
 * This method check that interval is valid regarding this, if not an exception will be thrown.
 */
DataArrayInt *DataArrayInt::substr(int tupleIdBg, int tupleIdEnd) const throw(INTERP_KERNEL::Exception)
{
  int nbt=getNumberOfTuples();
  if(tupleIdBg<0)
    throw INTERP_KERNEL::Exception("DataArrayInt::substr : The tupleIdBg parameter must be greater than 0 !");
  if(tupleIdBg>nbt)
    throw INTERP_KERNEL::Exception("DataArrayInt::substr : The tupleIdBg parameter is greater than number of tuples !");
  int trueEnd=tupleIdEnd;
  if(tupleIdEnd!=-1)
    {
      if(tupleIdEnd>nbt)
        throw INTERP_KERNEL::Exception("DataArrayInt::substr : The tupleIdBg parameter is greater or equal than number of tuples !");
    }
  else
    trueEnd=nbt;
  int nbComp=getNumberOfComponents();
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(trueEnd-tupleIdBg,nbComp);
  ret->copyStringInfoFrom(*this);
  std::copy(getConstPointer()+tupleIdBg*nbComp,getConstPointer()+trueEnd*nbComp,ret->getPointer());
  return ret;
}

/*!
 * Contrary to DataArrayInt::changeNbOfComponents method this method is \b not const. The content 
 * This method \b do \b not change the content of data but changes the splitting of this data seen by the caller.
 * This method makes the assumption that 'this' is already allocated. If not an exception will be thrown.
 * This method checks that getNbOfElems()%newNbOfCompo==0. If not an exception will be throw !
 * This method erases all components info set before call !
 */
void DataArrayInt::rearrange(int newNbOfCompo) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfElems=getNbOfElems();
  if(nbOfElems%newNbOfCompo!=0)
    throw INTERP_KERNEL::Exception("DataArrayInt::rearrange : nbOfElems%newNbOfCompo!=0 !");
  _nb_of_tuples=nbOfElems/newNbOfCompo;
  _info_on_compo.clear();
  _info_on_compo.resize(newNbOfCompo);
  declareAsNew();
}

/*!
 * This method builds a new instance of DataArrayInt (to deal with) that is reduction or an extension of 'this'.
 * if 'newNbOfComp' < this->getNumberOfComponents() a reduction is done and for each tuple 'newNbOfComp' first components are kept.
 * If 'newNbOfComp' > this->getNumberOfComponents() an extension is done, and for each components i such that i > getNumberOfComponents() 'dftValue' parameter is taken.
 */
DataArrayInt *DataArrayInt::changeNbOfComponents(int newNbOfComp, int dftValue) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(getNumberOfTuples(),newNbOfComp);
  const int *oldc=getConstPointer();
  int *nc=ret->getPointer();
  int nbOfTuples=getNumberOfTuples();
  int oldNbOfComp=getNumberOfComponents();
  int dim=std::min(oldNbOfComp,newNbOfComp);
  for(int i=0;i<nbOfTuples;i++)
    {
      int j=0;
      for(;j<dim;j++)
        nc[newNbOfComp*i+j]=oldc[i*oldNbOfComp+j];
      for(;j<newNbOfComp;j++)
        nc[newNbOfComp*i+j]=dftValue;
    }
  ret->setName(getName().c_str());
  for(int i=0;i<dim;i++)
    ret->setInfoOnComponent(i,getInfoOnComponent(i).c_str());
  ret->setName(getName().c_str());
  return ret;
}

void DataArrayInt::reAlloc(int nbOfTuples) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  _mem.reAlloc(_info_on_compo.size()*nbOfTuples);
  _nb_of_tuples=nbOfTuples;
  declareAsNew();
}

DataArrayInt *DataArrayInt::keepSelectedComponents(const std::vector<int>& compoIds) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret(DataArrayInt::New());
  int newNbOfCompo=compoIds.size();
  int oldNbOfCompo=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  ret->alloc(nbOfTuples,newNbOfCompo);
  ret->copyPartOfStringInfoFrom(*this,compoIds);
  const int *oldc=getConstPointer();
  int *nc=ret->getPointer();
  for(int i=0;i<nbOfTuples;i++)
    for(int j=0;j<newNbOfCompo;j++,nc++)
      *nc=oldc[i*oldNbOfCompo+compoIds[j]];
  ret->incrRef();
  return ret;
}

/*!
 * This method melds the components of 'this' with components of 'other'.
 * After this call in case of success, 'this' will contain a number of components equal to the sum of 'this'
 * before the call and the number of components of 'other'.
 * This method expects that 'this' and 'other' have exactly the same number of tuples. If not an exception is thrown.
 */
void DataArrayInt::meldWith(const DataArrayInt *other) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  other->checkAllocated();
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples!=other->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("DataArrayInt::meldWith : mismatch of number of tuples !");
  int nbOfComp1=getNumberOfComponents();
  int nbOfComp2=other->getNumberOfComponents();
  int *newArr=new int[nbOfTuples*(nbOfComp1+nbOfComp2)];
  int *w=newArr;
  const int *inp1=getConstPointer();
  const int *inp2=other->getConstPointer();
  for(int i=0;i<nbOfTuples;i++,inp1+=nbOfComp1,inp2+=nbOfComp2)
    {
      w=std::copy(inp1,inp1+nbOfComp1,w);
      w=std::copy(inp2,inp2+nbOfComp2,w);
    }
  useArray(newArr,true,CPP_DEALLOC,nbOfTuples,nbOfComp1+nbOfComp2);
  std::vector<int> compIds(nbOfComp2);
  for(int i=0;i<nbOfComp2;i++)
    compIds[i]=nbOfComp1+i;
  copyPartOfStringInfoFrom2(compIds,*other);
}

void DataArrayInt::setSelectedComponents(const DataArrayInt *a, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception)
{
  copyPartOfStringInfoFrom2(compoIds,*a);
  int partOfCompoSz=compoIds.size();
  int nbOfCompo=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  const int *ac=a->getConstPointer();
  int *nc=getPointer();
  for(int i=0;i<nbOfTuples;i++)
    for(int j=0;j<partOfCompoSz;j++,ac++)
      nc[nbOfCompo*i+compoIds[j]]=*ac;
}

void DataArrayInt::SetArrayIn(DataArrayInt *newArray, DataArrayInt* &arrayToSet)
{
  if(newArray!=arrayToSet)
    {
      if(arrayToSet)
        arrayToSet->decrRef();
      arrayToSet=newArray;
      if(arrayToSet)
        arrayToSet->incrRef();
    }
}

DataArrayInt *DataArrayInt::getIdsEqual(int val) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::getIdsEqual : the array must have only one component !");
  const int *cptr=getConstPointer();
  std::vector<int> res;
  int nbOfTuples=getNumberOfTuples();
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(*cptr==val)
      res.push_back(i);
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(res.size(),1);
  std::copy(res.begin(),res.end(),ret->getPointer());
  return ret;
}

DataArrayInt *DataArrayInt::getIdsEqualList(const std::vector<int>& vals) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::getIdsEqualList : the array must have only one component !");
  std::set<int> vals2(vals.begin(),vals.end());
  const int *cptr=getConstPointer();
  std::vector<int> res;
  int nbOfTuples=getNumberOfTuples();
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(vals2.find(*cptr)!=vals2.end())
      res.push_back(i);
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(res.size(),1);
  std::copy(res.begin(),res.end(),ret->getPointer());
  return ret;
}

DataArrayInt *DataArrayInt::Aggregate(const DataArrayInt *a1, const DataArrayInt *a2, int offsetA2)
{
  int nbOfComp=a1->getNumberOfComponents();
  if(nbOfComp!=a2->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array Aggregation !");
  int nbOfTuple1=a1->getNumberOfTuples();
  int nbOfTuple2=a2->getNumberOfTuples();
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbOfTuple1+nbOfTuple2-offsetA2,nbOfComp);
  int *pt=std::copy(a1->getConstPointer(),a1->getConstPointer()+nbOfTuple1*nbOfComp,ret->getPointer());
  std::copy(a2->getConstPointer()+offsetA2*nbOfComp,a2->getConstPointer()+nbOfTuple2*nbOfComp,pt);
  ret->copyStringInfoFrom(*a1);
  return ret;
}

DataArrayInt *DataArrayInt::Aggregate(const std::vector<const DataArrayInt *>& a) throw(INTERP_KERNEL::Exception)
{
  if(a.empty())
    throw INTERP_KERNEL::Exception("DataArrayInt::Aggregate : input list must be NON EMPTY !");
  std::vector<const DataArrayInt *>::const_iterator it=a.begin();
  int nbOfComp=(*it)->getNumberOfComponents();
  int nbt=(*it++)->getNumberOfTuples();
  for(int i=1;it!=a.end();it++,i++)
    {
      if((*it)->getNumberOfComponents()!=nbOfComp)
        throw INTERP_KERNEL::Exception("DataArrayInt::Aggregate : Nb of components mismatch for array aggregation !");
      nbt+=(*it)->getNumberOfTuples();
    }
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbt,nbOfComp);
  int *pt=ret->getPointer();
  for(it=a.begin();it!=a.end();it++)
    pt=std::copy((*it)->getConstPointer(),(*it)->getConstPointer()+(*it)->getNbOfElems(),pt);
  ret->copyStringInfoFrom(*(a[0]));
  return ret;
}

int DataArrayInt::getMaxValue(int& tupleId) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getMaxValue : must be applied on DataArrayDouble with only one component !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<=0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getMaxValue : array exists but number of tuples must be > 0 !");
  const int *vals=getConstPointer();
  const int *loc=std::max_element(vals,vals+nbOfTuples);
  tupleId=std::distance(vals,loc);
  return *loc;
}

int DataArrayInt::getMinValue(int& tupleId) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getMaxValue : must be applied on DataArrayDouble with only one component !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<=0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getMaxValue : array exists but number of tuples must be > 0 !");
  const int *vals=getConstPointer();
  const int *loc=std::min_element(vals,vals+nbOfTuples);
  tupleId=std::distance(vals,loc);
  return *loc;
}

void DataArrayInt::applyLin(int a, int b, int compoId) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int *ptr=getPointer()+compoId;
  int nbOfComp=getNumberOfComponents();
  int nbOfTuple=getNumberOfTuples();
  for(int i=0;i<nbOfTuple;i++,ptr+=nbOfComp)
    *ptr=a*(*ptr)+b;
  declareAsNew();
}

DataArrayInt *DataArrayInt::Meld(const DataArrayInt *a1, const DataArrayInt *a2) throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayInt *> arr(2);
  arr[0]=a1; arr[1]=a2;
  return Meld(arr);
}

DataArrayInt *DataArrayInt::Meld(const std::vector<const DataArrayInt *>& a) throw(INTERP_KERNEL::Exception)
{
  if(a.empty())
    throw INTERP_KERNEL::Exception("DataArrayInt::Meld : array must be NON empty !");
  std::vector<const DataArrayInt *>::const_iterator it;
  for(it=a.begin();it!=a.end();it++)
    (*it)->checkAllocated();
  it=a.begin();
  int nbOfTuples=(*it)->getNumberOfTuples();
  std::vector<int> nbc(a.size());
  std::vector<const int *> pts(a.size());
  nbc[0]=(*it)->getNumberOfComponents();
  pts[0]=(*it++)->getConstPointer();
  for(int i=1;it!=a.end();it++,i++)
    {
      if(nbOfTuples!=(*it)->getNumberOfTuples())
        throw INTERP_KERNEL::Exception("DataArrayInt::meld : mismatch of number of tuples !");
      nbc[i]=(*it)->getNumberOfComponents();
      pts[i]=(*it)->getConstPointer();
    }
  int totalNbOfComp=std::accumulate(nbc.begin(),nbc.end(),0);
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbOfTuples,totalNbOfComp);
  int *retPtr=ret->getPointer();
  for(int i=0;i<nbOfTuples;i++)
    for(int j=0;j<(int)a.size();j++)
      {
        retPtr=std::copy(pts[j],pts[j]+nbc[j],retPtr);
        pts[j]+=nbc[j];
      }
  int k=0;
  for(int i=0;i<(int)a.size();i++)
    for(int j=0;j<nbc[i];j++,k++)
      ret->setInfoOnComponent(k,a[i]->getInfoOnComponent(j).c_str());
  return ret;
}

/*!
 * This method create a minimal partition of groups 'groups' the std::iota array of size 'newNb'.
 * This method returns an array of size 'newNb' that specifies for each item at which familyId it owns to, and this method returns
 * for each group the familyId it contains. If an id so that id<newNb and that appears in no groups will appears with 0 in return array.
 *
 * @param groups in arrays specifying ids of each groups.
 * @param newNb specifies size of whole set. Must be at least equal to max eltid in 'groups'.
 * @return an array of size newNb specifying fid of each item.
 */
DataArrayInt *DataArrayInt::MakePartition(const std::vector<const DataArrayInt *>& groups, int newNb, std::vector< std::vector<int> >& fidsOfGroups)
{
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(newNb,1);
  int *retPtr=ret->getPointer();
  std::fill(retPtr,retPtr+newNb,0);
  int fid=1;
  for(std::vector<const DataArrayInt *>::const_iterator iter=groups.begin();iter!=groups.end();iter++)
    {
      const int *ptr=(*iter)->getConstPointer();
      int nbOfElem=(*iter)->getNbOfElems();
      int sfid=fid;
      for(int j=0;j<sfid;j++)
        {
          bool found=false;
          for(int i=0;i<nbOfElem;i++)
            {
              if(retPtr[ptr[i]]==j)
                {
                  retPtr[ptr[i]]=fid;
                  found=true;
                }
            }
          if(found)
            fid++;
        }
    }
  fidsOfGroups.clear();
  fidsOfGroups.resize(groups.size());
  int grId=0;
  for(std::vector<const DataArrayInt *>::const_iterator iter=groups.begin();iter!=groups.end();iter++,grId++)
    {
      std::set<int> tmp;
      const int *ptr=(*iter)->getConstPointer();
      int nbOfElem=(*iter)->getNbOfElems();
      for(const int *p=ptr;p!=ptr+nbOfElem;p++)
        tmp.insert(retPtr[*p]);
      fidsOfGroups[grId].insert(fidsOfGroups[grId].end(),tmp.begin(),tmp.end());
    }
  return ret;
}

DataArrayInt *DataArrayInt::BuildUnion(const std::vector<const DataArrayInt *>& a) throw(INTERP_KERNEL::Exception)
{
  int valm=std::numeric_limits<int>::max();
  for(std::vector<const DataArrayInt *>::const_iterator it=a.begin();it!=a.end();it++)
    {
      (*it)->checkAllocated();
      if((*it)->getNumberOfComponents()!=1)
        throw INTERP_KERNEL::Exception("DataArrayInt::BuildUnion : only single component allowed !");
      int tmp1;
      valm=std::min((*it)->getMinValue(tmp1),valm);
    }
  if(valm<0)
    throw INTERP_KERNEL::Exception("DataArrayInt::BuildUnion : a negative value has been detected !");
  //
  std::set<int> r;
  for(std::vector<const DataArrayInt *>::const_iterator it=a.begin();it!=a.end();it++)
    {
      const int *pt=(*it)->getConstPointer();
      int nbOfTuples=(*it)->getNumberOfTuples();
      r.insert(pt,pt+nbOfTuples);
    }
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(r.size(),1);
  std::copy(r.begin(),r.end(),ret->getPointer());
  return ret;
}

DataArrayInt *DataArrayInt::BuildIntersection(const std::vector<const DataArrayInt *>& a) throw(INTERP_KERNEL::Exception)
{
  int valm=std::numeric_limits<int>::max();
  for(std::vector<const DataArrayInt *>::const_iterator it=a.begin();it!=a.end();it++)
    {
      (*it)->checkAllocated();
      if((*it)->getNumberOfComponents()!=1)
        throw INTERP_KERNEL::Exception("DataArrayInt::BuildUnion : only single component allowed !");
      int tmp1;
      valm=std::min((*it)->getMinValue(tmp1),valm);
    }
  if(valm<0)
    throw INTERP_KERNEL::Exception("DataArrayInt::BuildUnion : a negative value has been detected !");
  //
  std::set<int> r;
  for(std::vector<const DataArrayInt *>::const_iterator it=a.begin();it!=a.end();it++)
    {
      const int *pt=(*it)->getConstPointer();
      int nbOfTuples=(*it)->getNumberOfTuples();
      std::set<int> s1(pt,pt+nbOfTuples);
      if(it!=a.begin())
        {
          std::set<int> r2;
          std::set_intersection(r.begin(),r.end(),s1.begin(),s1.end(),inserter(r2,r2.end()));
          r=r2;
        }
      else
        r=s1;
    }
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(r.size(),1);
  std::copy(r.begin(),r.end(),ret->getPointer());
  return ret;
}

DataArrayInt *DataArrayInt::buildComplement(int nbOfElement) const throw(INTERP_KERNEL::Exception)
{
   checkAllocated();
   if(getNumberOfComponents()!=1)
     throw INTERP_KERNEL::Exception("DataArrayInt::buildComplement : only single component allowed !");
   std::vector<bool> tmp(nbOfElement);
   const int *pt=getConstPointer();
   int nbOfTuples=getNumberOfTuples();
   for(const int *w=pt;w!=pt+nbOfTuples;w++)
     if(*w>=0 && *w<nbOfElement)
       tmp[*w]=true;
     else
       throw INTERP_KERNEL::Exception("DataArrayInt::buildComplement : an element is not in valid range : [0,nbOfElement) !");
   int nbOfRetVal=std::count(tmp.begin(),tmp.end(),false);
   DataArrayInt *ret=DataArrayInt::New();
   ret->alloc(nbOfRetVal,1);
   int j=0;
   int *retPtr=ret->getPointer();
   for(int i=0;i<nbOfElement;i++)
     if(!tmp[i])
       retPtr[j++]=i;
   return ret;
}

DataArrayInt *DataArrayInt::buildSubstraction(const DataArrayInt *other) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  other->checkAllocated();
  if(getNumberOfComponents()!=1)
     throw INTERP_KERNEL::Exception("DataArrayInt::buildSubstraction : only single component allowed !");
  if(other->getNumberOfComponents()!=1)
     throw INTERP_KERNEL::Exception("DataArrayInt::buildSubstraction : only single component allowed for other type !");
  const int *pt=getConstPointer();
  int nbOfTuples=getNumberOfTuples();
  std::set<int> s1(pt,pt+nbOfTuples);
  pt=other->getConstPointer();
  nbOfTuples=other->getNumberOfTuples();
  std::set<int> s2(pt,pt+nbOfTuples);
  std::vector<int> r;
  std::set_difference(s1.begin(),s1.end(),s2.begin(),s2.end(),std::back_insert_iterator< std::vector<int> >(r));
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(r.size(),1);
  std::copy(r.begin(),r.end(),ret->getPointer());
  return ret;
}

DataArrayInt *DataArrayInt::buildUnion(const DataArrayInt *other) const throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayInt *>arrs(2);
  arrs[0]=this; arrs[1]=other;
  return BuildUnion(arrs);
}

DataArrayInt *DataArrayInt::buildIntersection(const DataArrayInt *other) const throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayInt *>arrs(2);
  arrs[0]=this; arrs[1]=other;
  return BuildIntersection(arrs);
}

/*!
 * This method could be usefull for returned DataArrayInt marked as index. Some methods that generate such DataArrayInt instances:
 * - ParaMEDMEM::MEDCouplingUMesh::buildDescendingConnectivity
 * - ParaMEDMEM::MEDCouplingUMesh::getNodalConnectivityIndex
 * This method makes the assumption that 'this' is allocated and has exactly one component and 2 or more tuples. If not an exception is thrown.
 * This method retrives a newly created DataArrayInt instance with 1 component and this->getNumberOfTuples()-1 tuples.
 * If this contains [1,3,6,7,7,9,15] -> returned array will contain [2,3,1,0,2,6].
 */
DataArrayInt *DataArrayInt::deltaShiftIndex() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
     throw INTERP_KERNEL::Exception("DataArrayInt::deltaShiftIndex : only single component allowed !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<2)
    throw INTERP_KERNEL::Exception("DataArrayInt::deltaShiftIndex : 1 tuple at least must be present in 'this' !");
  const int *ptr=getPointer();
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbOfTuples-1,1);
  int *out=ret->getPointer();
  std::transform(ptr+1,ptr+nbOfTuples,ptr,out,std::minus<int>());
  return ret;
}

/*!
 * This method returns all different values found in 'this'.
 */
std::set<int> DataArrayInt::getDifferentValues() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  std::set<int> ret;
  ret.insert(getConstPointer(),getConstPointer()+getNbOfElems());
  return ret;
}

int *DataArrayInt::CheckAndPreparePermutation(const int *start, const int *end)
{
  int sz=std::distance(start,end);
  int *ret=new int[sz];
  int *work=new int[sz];
  std::copy(start,end,work);
  std::sort(work,work+sz);
  if(std::unique(work,work+sz)!=work+sz)
    {
      delete [] work;
      delete [] ret;
      throw INTERP_KERNEL::Exception("Some elements are equals in the specified array !");
    }
  int *iter2=ret;
  for(const int *iter=start;iter!=end;iter++,iter2++)
    *iter2=std::distance(work,std::find(work,work+sz,*iter));
  delete [] work;
  return ret;
}
