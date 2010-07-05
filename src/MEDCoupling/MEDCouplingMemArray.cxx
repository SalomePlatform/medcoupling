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

#include <set>
#include <functional>

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

bool DataArray::areInfoEquals(const DataArray& other) const
{
  if(_nb_of_tuples!=other._nb_of_tuples)
    return false;
  if(_name!=other._name)
    return false;
  return _info_on_compo==other._info_on_compo;
}

DataArrayDouble *DataArrayDouble::New()
{
  return new DataArrayDouble;
}

DataArrayDouble *DataArrayDouble::deepCopy() const
{
  return new DataArrayDouble(*this);
}

DataArrayDouble *DataArrayDouble::performCpy(bool deepCpy) const
{
  if(deepCpy)
    return deepCopy();
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

bool DataArrayDouble::isEqual(const DataArrayDouble& other, double prec) const
{
  if(!areInfoEquals(other))
    return false;
  return _mem.isEqual(other._mem,prec);
}

void DataArrayDouble::reAlloc(int nbOfTuples)
{
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

void DataArrayDouble::setArrayIn(DataArrayDouble *newArray, DataArrayDouble* &arrayToSet)
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

DataArrayDouble *DataArrayDouble::aggregate(const DataArrayDouble *a1, const DataArrayDouble *a2)
{
  int nbOfComp=a1->getNumberOfComponents();
  if(nbOfComp!=a2->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array aggregation !");
  int nbOfTuple1=a1->getNumberOfTuples();
  int nbOfTuple2=a2->getNumberOfTuples();
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbOfTuple1+nbOfTuple2,nbOfComp);
  double *pt=std::copy(a1->getConstPointer(),a1->getConstPointer()+nbOfTuple1*nbOfComp,ret->getPointer());
  std::copy(a2->getConstPointer(),a2->getConstPointer()+nbOfTuple2*nbOfComp,pt);
  ret->copyStringInfoFrom(*a1);
  return ret;
}

DataArrayDouble *DataArrayDouble::add(const DataArrayDouble *a1, const DataArrayDouble *a2)
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

void DataArrayDouble::addEqual(const DataArrayDouble *other)
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

DataArrayDouble *DataArrayDouble::substract(const DataArrayDouble *a1, const DataArrayDouble *a2)
{
  int nbOfComp=a1->getNumberOfComponents();
  if(nbOfComp!=a2->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array substract !");
  int nbOfTuple=a1->getNumberOfTuples();
  if(nbOfTuple!=a2->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array substract !");
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbOfTuple,nbOfComp);
  std::transform(a1->getConstPointer(),a1->getConstPointer()+nbOfTuple*nbOfComp,a2->getConstPointer(),ret->getPointer(),std::minus<double>());
  ret->copyStringInfoFrom(*a1);
  return ret;
}

void DataArrayDouble::substractEqual(const DataArrayDouble *other)
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

DataArrayDouble *DataArrayDouble::multiply(const DataArrayDouble *a1, const DataArrayDouble *a2)
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

void DataArrayDouble::multiplyEqual(const DataArrayDouble *other)
{
  int nbOfTuple=getNumberOfTuples();
  int nbOfTuple2=other->getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  int nbOfComp2=other->getNumberOfComponents();
  if(nbOfTuple!=nbOfTuple2)
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array multiplyEqual !");
  DataArrayDouble *ret=0;
  if(nbOfComp==nbOfComp2)
    {
      ret=DataArrayDouble::New();
      ret->alloc(nbOfTuple,nbOfComp);
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

DataArrayDouble *DataArrayDouble::divide(const DataArrayDouble *a1, const DataArrayDouble *a2)
{
  int nbOfComp=a1->getNumberOfComponents();
  if(nbOfComp!=a2->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array divide !");
  int nbOfTuple=a1->getNumberOfTuples();
  if(nbOfTuple!=a2->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array divide !");
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbOfTuple,nbOfComp);
  std::transform(a1->getConstPointer(),a1->getConstPointer()+nbOfTuple*nbOfComp,a2->getConstPointer(),ret->getPointer(),std::divides<double>());
  ret->copyStringInfoFrom(*a1);
  return ret;
}

void DataArrayDouble::divideEqual(const DataArrayDouble *other)
{
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=other->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array divideEqual !");
  int nbOfTuple=getNumberOfTuples();
  if(nbOfTuple!=other->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array divideEqual !");
  std::transform(getConstPointer(),getConstPointer()+nbOfTuple*nbOfComp,other->getConstPointer(),getPointer(),std::divides<double>());
  declareAsNew();
}

DataArrayInt *DataArrayInt::New()
{
  return new DataArrayInt;
}

DataArrayInt *DataArrayInt::deepCopy() const
{
  return new DataArrayInt(*this);
}

DataArrayInt *DataArrayInt::performCpy(bool deepCpy) const
{
  if(deepCpy)
    return deepCopy();
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

bool DataArrayInt::isEqual(const DataArrayInt& other) const
{
  if(!areInfoEquals(other))
    return false;
  return _mem.isEqual(other._mem,0);
}

void DataArrayInt::useArray(const int *array, bool ownership,  DeallocType type, int nbOfTuple, int nbOfCompo)
{
  _nb_of_tuples=nbOfTuple;
  _info_on_compo.resize(nbOfCompo);
  _mem.useArray(array,ownership,type,nbOfTuple*nbOfCompo);
  declareAsNew();
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

void DataArrayInt::reAlloc(int nbOfTuples)
{
  _mem.reAlloc(_info_on_compo.size()*nbOfTuples);
  _nb_of_tuples=nbOfTuples;
  declareAsNew();
}

void DataArrayInt::setArrayIn(DataArrayInt *newArray, DataArrayInt* &arrayToSet)
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

DataArrayInt *DataArrayInt::aggregate(const DataArrayInt *a1, const DataArrayInt *a2, int offsetA2)
{
  int nbOfComp=a1->getNumberOfComponents();
  if(nbOfComp!=a2->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array aggregation !");
  int nbOfTuple1=a1->getNumberOfTuples();
  int nbOfTuple2=a2->getNumberOfTuples();
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbOfTuple1+nbOfTuple2-offsetA2,nbOfComp);
  int *pt=std::copy(a1->getConstPointer(),a1->getConstPointer()+nbOfTuple1*nbOfComp,ret->getPointer());
  std::copy(a2->getConstPointer()+offsetA2*nbOfComp,a2->getConstPointer()+nbOfTuple2*nbOfComp,pt);
  ret->copyStringInfoFrom(*a1);
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
DataArrayInt *DataArrayInt::makePartition(const std::vector<DataArrayInt *>& groups, int newNb, std::vector< std::vector<int> >& fidsOfGroups)
{
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(newNb,1);
  int *retPtr=ret->getPointer();
  std::fill(retPtr,retPtr+newNb,0);
  int fid=1;
  for(std::vector<DataArrayInt *>::const_iterator iter=groups.begin();iter!=groups.end();iter++)
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
  for(std::vector<DataArrayInt *>::const_iterator iter=groups.begin();iter!=groups.end();iter++,grId++)
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

int *DataArrayInt::checkAndPreparePermutation(const int *start, const int *end)
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
