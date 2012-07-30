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

template<int SPACEDIM>
void DataArrayDouble::findCommonTuplesAlg(const double *bbox, int nbNodes, int limitNodeId, double prec, std::vector<int>& c, std::vector<int>& cI) const
{
  const double *coordsPtr=getConstPointer();
  BBTree<SPACEDIM,int> myTree(bbox,0,0,nbNodes,prec/10);
  std::vector<bool> isDone(nbNodes);
  for(int i=0;i<nbNodes;i++)
    {
      if(!isDone[i])
        {
          std::vector<int> intersectingElems;
          myTree.getElementsAroundPoint(coordsPtr+i*SPACEDIM,intersectingElems);
          if(intersectingElems.size()>1)
            {
              std::vector<int> commonNodes;
              for(std::vector<int>::const_iterator it=intersectingElems.begin();it!=intersectingElems.end();it++)
                if(*it!=i)
                  if(*it>=limitNodeId)
                    {
                      commonNodes.push_back(*it);
                      isDone[*it]=true;
                    }
              if(!commonNodes.empty())
                {
                  cI.push_back(cI.back()+(int)commonNodes.size()+1);
                  c.push_back(i);
                  c.insert(c.end(),commonNodes.begin(),commonNodes.end());
                }
            }
        }
    }
}

template<int SPACEDIM>
void DataArrayDouble::findTupleIdsNearTuplesAlg(const BBTree<SPACEDIM,int>& myTree, const double *pos, int nbOfTuples, double eps,
                                                std::vector<int>& c, std::vector<int>& cI) const
{
  const double *coordsPtr=getConstPointer();
  for(int i=0;i<nbOfTuples;i++)
    {
      std::vector<int> intersectingElems;
      myTree.getElementsAroundPoint(pos+i*SPACEDIM,intersectingElems);
      std::vector<int> commonNodes;
      for(std::vector<int>::const_iterator it=intersectingElems.begin();it!=intersectingElems.end();it++)
        commonNodes.push_back(*it);
      cI.push_back(cI.back()+(int)commonNodes.size());
      c.insert(c.end(),commonNodes.begin(),commonNodes.end());
    }
}

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
  std::size_t newNbOfCompo=compoIds.size();
  for(std::size_t i=0;i<newNbOfCompo;i++)
    if(compoIds[i]>=nbOfCompoOth || compoIds[i]<0)
      {
        std::ostringstream oss; oss << "Specified component id is out of range (" << compoIds[i] << ") compared with nb of actual components (" << nbOfCompoOth << ")";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  for(std::size_t i=0;i<newNbOfCompo;i++)
    setInfoOnComponent((int)i,other.getInfoOnComponent(compoIds[i]).c_str());
}

void DataArray::copyPartOfStringInfoFrom2(const std::vector<int>& compoIds, const DataArray& other) throw(INTERP_KERNEL::Exception)
{
  int nbOfCompo=getNumberOfComponents();
  std::size_t partOfCompoToSet=compoIds.size();
  if((int)partOfCompoToSet!=other.getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Given compoIds has not the same size as number of components of given array !");
  for(std::size_t i=0;i<partOfCompoToSet;i++)
    if(compoIds[i]>=nbOfCompo || compoIds[i]<0)
      {
        std::ostringstream oss; oss << "Specified component id is out of range (" << compoIds[i] << ") compared with nb of actual components (" << nbOfCompo << ")";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  for(std::size_t i=0;i<partOfCompoToSet;i++)
    setInfoOnComponent(compoIds[i],other.getInfoOnComponent((int)i).c_str());
}

bool DataArray::areInfoEqualsIfNotWhy(const DataArray& other, std::string& reason) const
{
  std::ostringstream oss;
  if(_nb_of_tuples!=other._nb_of_tuples)
    {
      oss << "Number of tuples of DataArray mismatch : this number of tuples=" << _nb_of_tuples << " other number of tuples=" << other._nb_of_tuples;
      reason=oss.str();
      return false;
    }
  if(_name!=other._name)
    {
      oss << "Names DataArray mismatch : this name=\"" << _name << " other name=\"" << other._name << "\" !";
      reason=oss.str();
      return false;
    }
  if(_info_on_compo!=other._info_on_compo)
    {
      oss << "Components DataArray mismatch : \nThis components=";
      for(std::vector<std::string>::const_iterator it=_info_on_compo.begin();it!=_info_on_compo.end();it++)
        oss << "\"" << *it << "\",";
      oss << "\nOther components=";
      for(std::vector<std::string>::const_iterator it=other._info_on_compo.begin();it!=other._info_on_compo.end();it++)
        oss << "\"" << *it << "\",";
      reason=oss.str();
      return false;
    }
  return true;
}

bool DataArray::areInfoEquals(const DataArray& other) const
{
  std::string tmp;
  return areInfoEqualsIfNotWhy(other,tmp);
}

void DataArray::reprWithoutNameStream(std::ostream& stream) const
{
  stream << "Number of components : "<< getNumberOfComponents() << "\n";
  stream << "Info of these components : ";
  for(std::vector<std::string>::const_iterator iter=_info_on_compo.begin();iter!=_info_on_compo.end();iter++)
    stream << "\"" << *iter << "\"   ";
  stream << "\n";
}

void DataArray::setInfoOnComponents(const std::vector<std::string>& info) throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=(int)info.size())
    {
      std::ostringstream oss; oss << "DataArray::setInfoOnComponents : input is of size " << info.size() << " whereas number of components is equal to " << getNumberOfComponents() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  _info_on_compo=info;
}

std::vector<std::string> DataArray::getVarsOnComponent() const
{
  int nbOfCompo=(int)_info_on_compo.size();
  std::vector<std::string> ret(nbOfCompo);
  for(int i=0;i<nbOfCompo;i++)
    ret[i]=getVarOnComponent(i);
  return ret;
}

std::vector<std::string> DataArray::getUnitsOnComponent() const
{
  int nbOfCompo=(int)_info_on_compo.size();
  std::vector<std::string> ret(nbOfCompo);
  for(int i=0;i<nbOfCompo;i++)
    ret[i]=getUnitOnComponent(i);
  return ret;
}

std::string DataArray::getInfoOnComponent(int i) const throw(INTERP_KERNEL::Exception)
{
  if(i<(int)_info_on_compo.size() && i>=0)
    return _info_on_compo[i];
  else
    {
      std::ostringstream oss; oss << "DataArray::getInfoOnComponent : Specified component id is out of range (" << i << ") compared with nb of actual components (" << (int) _info_on_compo.size();
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * In the info part of i_th component this method returns the var part.
 * For example, if getInfoOnComponent(0) return "SIGXY (N/m^2)", getVarOnComponent(0) will return "SIGXY"
 */
std::string DataArray::getVarOnComponent(int i) const throw(INTERP_KERNEL::Exception)
{
  if(i<(int)_info_on_compo.size() && i>=0)
    {
      return GetVarNameFromInfo(_info_on_compo[i]);
    }
  else
    {
      std::ostringstream oss; oss << "DataArray::getVarOnComponent : Specified component id is out of range  (" << i << ") compared with nb of actual components (" << (int) _info_on_compo.size();
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * In the info part of i_th component this method returns the var part.
 * For example, if getInfoOnComponent(0) return "SIGXY (N/m^2)", getUnitOnComponent(0) will return "N/m^2"
 */
std::string DataArray::getUnitOnComponent(int i) const throw(INTERP_KERNEL::Exception)
{
  if(i<(int)_info_on_compo.size() && i>=0)
    {
      return GetUnitFromInfo(_info_on_compo[i]);
    }
  else
    {
      std::ostringstream oss; oss << "DataArray::getUnitOnComponent : Specified component id is out of range  (" << i << ") compared with nb of actual components (" << (int) _info_on_compo.size();
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

std::string DataArray::GetVarNameFromInfo(const std::string& info) throw(INTERP_KERNEL::Exception)
{
  std::size_t p1=info.find_last_of('[');
  std::size_t p2=info.find_last_of(']');
  if(p1==std::string::npos || p2==std::string::npos)
    return info;
  if(p1>p2)
    return info;
  if(p1==0)
    return std::string();
  std::size_t p3=info.find_last_not_of(' ',p1-1);
  return info.substr(0,p3+1);
}

std::string DataArray::GetUnitFromInfo(const std::string& info) throw(INTERP_KERNEL::Exception)
{
  std::size_t p1=info.find_last_of('[');
  std::size_t p2=info.find_last_of(']');
  if(p1==std::string::npos || p2==std::string::npos)
    return std::string();
  if(p1>p2)
    return std::string();
  return info.substr(p1+1,p2-p1-1);
}

void DataArray::setInfoOnComponent(int i, const char *info) throw(INTERP_KERNEL::Exception)
{
  if(i<(int)_info_on_compo.size() && i>=0)
    _info_on_compo[i]=info;
  else
    {
      std::ostringstream oss; oss << "DataArray::setInfoOnComponent : Specified component id is out of range  (" << i << ") compared with nb of actual components (" << (int) _info_on_compo.size();
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void DataArray::checkNbOfTuples(int nbOfTuples, const char *msg) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfTuples()!=nbOfTuples)
    {
      std::ostringstream oss; oss << msg << " : mismatch number of tuples : expected " <<  nbOfTuples << " having " << getNumberOfTuples() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void DataArray::checkNbOfComps(int nbOfCompo, const char *msg) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=nbOfCompo)
    {
      std::ostringstream oss; oss << msg << " : mismatch number of components : expected " << nbOfCompo << " having " << getNumberOfComponents() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void DataArray::checkNbOfElems(int nbOfElems, const char *msg) const throw(INTERP_KERNEL::Exception)
{
  if(getNbOfElems()!=nbOfElems)
    {
      std::ostringstream oss; oss << msg << " : mismatch number of elems : Expected " << nbOfElems << " having " << getNbOfElems() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void DataArray::checkNbOfTuplesAndComp(const DataArray& other, const char *msg) const throw(INTERP_KERNEL::Exception)
{
   if(getNumberOfTuples()!=other.getNumberOfTuples())
    {
      std::ostringstream oss; oss << msg << " : mismatch number of tuples : expected " <<  other.getNumberOfTuples() << " having " << getNumberOfTuples() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(getNumberOfComponents()!=other.getNumberOfComponents())
    {
      std::ostringstream oss; oss << msg << " : mismatch number of components : expected " << other.getNumberOfComponents() << " having " << getNumberOfComponents() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void DataArray::checkNbOfTuplesAndComp(int nbOfTuples, int nbOfCompo, const char *msg) const throw(INTERP_KERNEL::Exception)
{
  checkNbOfTuples(nbOfTuples,msg);
  checkNbOfComps(nbOfCompo,msg);
}

/*!
 * Simply this method checks that \b value is in [0,\b ref).
 */
void DataArray::CheckValueInRange(int ref, int value, const char *msg) throw(INTERP_KERNEL::Exception)
{
  if(value<0 || value>=ref)
    {
      std::ostringstream oss; oss << "DataArray::CheckValueInRange : " << msg  << " ! Expected in range [0," << ref << "[ having " << value << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * This method checks that [\b start, \b end) is compliant with ref length \b value.
 * typicaly start in [0,\b value) and end in [0,\b value). If value==start and start==end, it is supported.
 */
void DataArray::CheckValueInRangeEx(int value, int start, int end, const char *msg) throw(INTERP_KERNEL::Exception)
{
  if(start<0 || start>=value)
    {
      if(value!=start || end!=start)
        {
          std::ostringstream oss; oss << "DataArray::CheckValueInRange : " << msg  << " ! Expected start " << start << " of range in [0," << value << "[ !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  if(end<0 || end>value)
    {
      std::ostringstream oss; oss << "DataArray::CheckClosingParInRange : " << msg  << " ! Expected start " << end << " of range in [0," << value << "[ !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void DataArray::CheckClosingParInRange(int ref, int value, const char *msg) throw(INTERP_KERNEL::Exception)
{
  if(value<0 || value>ref)
    {
      std::ostringstream oss; oss << "DataArray::CheckClosingParInRange : " << msg  << " ! Expected a range in [0," << ref << "] having closing open parenthesis " << value << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

int DataArray::GetNumberOfItemGivenBES(int begin, int end, int step, const char *msg) throw(INTERP_KERNEL::Exception)
{
  if(end<begin)
    {
      std::ostringstream oss; oss << msg << " : end before begin !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(step<=0)
    {
      std::ostringstream oss; oss << msg << " : invalid step should be > 0 !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return (end-1-begin)/step+1;
}

int DataArray::GetNumberOfItemGivenBESRelative(int begin, int end, int step, const char *msg) throw(INTERP_KERNEL::Exception)
{
  if(step==0)
    throw INTERP_KERNEL::Exception("DataArray::GetNumberOfItemGivenBES : step=0 is not allowed !");
  if(end<begin && step>0)
    {
      std::ostringstream oss; oss << msg << " : end before begin whereas step is positive !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(begin<end && step<0)
    {
      std::ostringstream oss; oss << msg << " : invalid step should be > 0 !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(begin!=end)
    return (std::max(begin,end)-1-std::min(begin,end))/std::abs(step)+1;
  else
    return 0;
}

int DataArray::GetPosOfItemGivenBESRelativeNoThrow(int value, int begin, int end, int step) throw(INTERP_KERNEL::Exception)
{
  if(step!=0)
    {
      if(step>0)
        {
          if(begin<=value && value<end)
            {
              if((value-begin)%step==0)
                return (value-begin)/step;
              else
                return -1;
            }
          else
            return -1;
        }
      else
        {
          if(begin>=value && value>end)
            {
              if((begin-value)%(-step)==0)
                return (begin-value)/(-step);
              else
                return -1;
            }
          else
            return -1;
        }
    }
  else
    return -1;
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

/*!
 * This method differs from DataArray::setInfoOnComponents in the sense that if 'this->getNumberOfComponents()!=info.size()'
 * and if 'this' is not allocated it will change the number of components of 'this'.
 * If 'this->getNumberOfComponents()==info.size()' the behaviour is the same than DataArray::setInfoOnComponents method.
 * If 'this->getNumberOfComponents()!=info.size()' and the 'this' is already allocated an exception will be thrown.
 */
void DataArrayDouble::setInfoAndChangeNbOfCompo(const std::vector<std::string>& info) throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=(int)info.size())
    {
      if(!isAllocated())
        _info_on_compo=info;
      else
        {
          std::ostringstream oss; oss << "DataArrayDouble::setInfoAndChangeNbOfCompo : input is of size " << info.size() << " whereas number of components is equal to " << getNumberOfComponents() << "  and this is already allocated !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  else
    _info_on_compo=info;
}

/*!
 * This method returns the only one value in 'this', if and only if number of elements (nb of tuples * nb of components) is equal to 1, and that 'this' is allocated.
 * If one or more conditions is not fulfilled an exception will be thrown.
 */
double DataArrayDouble::doubleValue() const throw(INTERP_KERNEL::Exception)
{
  if(isAllocated())
    {
      if(getNbOfElems()==1)
        {
          return *getConstPointer();
        }
      else
        throw INTERP_KERNEL::Exception("DataArrayDouble::doubleValue : DataArrayDouble instance is allocated but number of elements is not equal to 1 !");
    }
  else
    throw INTERP_KERNEL::Exception("DataArrayDouble::doubleValue : DataArrayDouble instance is not allocated !");
}

/*!
 * This method should be called on an allocated DataArrayDouble instance. If not an exception will be throw !
 * This method checks the number of tupes. If it is equal to 0, it returns true, if not false is returned.
 */
bool DataArrayDouble::empty() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  return getNumberOfTuples()==0;
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

void DataArrayDouble::cpyFrom(const DataArrayDouble& other) throw(INTERP_KERNEL::Exception)
{
  other.checkAllocated();
  int nbOfTuples=other.getNumberOfTuples();
  int nbOfComp=other.getNumberOfComponents();
  allocIfNecessary(nbOfTuples,nbOfComp);
  int nbOfElems=nbOfTuples*nbOfComp;
  double *pt=getPointer();
  const double *ptI=other.getConstPointer();
  for(int i=0;i<nbOfElems;i++)
    pt[i]=ptI[i];
  copyStringInfoFrom(other);
}

void DataArrayDouble::allocIfNecessary(int nbOfTuple, int nbOfCompo)
{
  if(isAllocated())
    {
      if(nbOfTuple!=getNumberOfTuples() || nbOfCompo!=getNumberOfComponents())
        alloc(nbOfTuple,nbOfCompo);
    }
  else
    alloc(nbOfTuple,nbOfCompo);
}

void DataArrayDouble::alloc(int nbOfTuple, int nbOfCompo) throw(INTERP_KERNEL::Exception)
{
  if(nbOfTuple<0 || nbOfCompo<0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::alloc : request for negative length of data !");
  _nb_of_tuples=nbOfTuple;
  _info_on_compo.resize(nbOfCompo);
  _mem.alloc(nbOfCompo*_nb_of_tuples);
  declareAsNew();
}

void DataArrayDouble::fillWithZero() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  _mem.fillWithValue(0.);
  declareAsNew();
}

void DataArrayDouble::fillWithValue(double val) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  _mem.fillWithValue(val);
  declareAsNew();
}

void DataArrayDouble::iota(double init) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::iota : works only for arrays with only one component, you can call 'rearrange' method before !");
  double *ptr=getPointer();
  int ntuples=getNumberOfTuples();
  for(int i=0;i<ntuples;i++)
    ptr[i]=init+double(i);
  declareAsNew();
}

bool DataArrayDouble::isUniform(double val, double eps) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::isUniform : must be applied on DataArrayDouble with only one component, you can call 'rearrange' method before !");
  int nbOfTuples=getNumberOfTuples();
  const double *w=getConstPointer();
  const double *end2=w+nbOfTuples;
  const double vmin=val-eps;
  const double vmax=val+eps;
  for(;w!=end2;w++)
    if(*w<vmin || *w>vmax)
      return false;
  return true;
}

void DataArrayDouble::sort(bool asc) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::sort : only supported with 'this' array with ONE component !");
  _mem.sort(asc);
}

void DataArrayDouble::reverse() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::reverse : only supported with 'this' array with ONE component !");
  _mem.reverse();
}

/*!
 * This method check that (Maths) array consistently INCREASING or DECREASING in value,
 * with at least absolute difference value of |eps| at each step.
 * if not an exception will be thrown.
 */
 void DataArrayDouble::checkMonotonic(bool increasing, double eps) const throw(INTERP_KERNEL::Exception)
{
  if(!isMonotonic(increasing, eps))
    {
      if (increasing)
        {
          throw INTERP_KERNEL::Exception("DataArrayDouble::checkMonotonic : 'this' is not INCREASING monotonic !");
        }
      else
        {
          throw INTERP_KERNEL::Exception("DataArrayDouble::checkMonotonic : 'this' is not DECREASING monotonic !");
        }
    }
}

/*!
 * This method check that (Maths) array consistently INCREASING or DECREASING in value,
 * with at least absolute difference value of |eps| at each step.
 */
bool DataArrayDouble::isMonotonic(bool increasing, double eps) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::isMonotonic : only supported with 'this' array with ONE component !");
  int nbOfElements=getNumberOfTuples();
  const double *ptr=getConstPointer();
  if(nbOfElements==0)
    return true;
  double ref=ptr[0];
  double absEps=fabs(eps);
  if (increasing)
    {
      for(int i=1;i<nbOfElements;i++)
        {
          if(ptr[i]<(ref+absEps))
            return false;
          ref=ptr[i];
        }
      return true;
    }
  else
    {
      for(int i=1;i<nbOfElements;i++)
        {
          if(ptr[i]>(ref-absEps))
            return false;
          ref=ptr[i];
        }
      return true;
    }
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

void DataArrayDouble::writeVTK(std::ostream& ofs, int indent, const char *nameInFile) const throw(INTERP_KERNEL::Exception)
{
  std::string idt(indent,' ');
  ofs.precision(15);
  ofs << idt << "<DataArray type=\"Float32\" Name=\"" << nameInFile << "\" NumberOfComponents=\"" << getNumberOfComponents() << "\"";
  ofs << " format=\"ascii\" RangeMin=\"" << getMinValueInArray() << "\" RangeMax=\"" << getMaxValueInArray() << "\">\n" << idt;
  std::copy(begin(),end(),std::ostream_iterator<double>(ofs," "));
  ofs << std::endl << idt << "</DataArray>\n";
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

bool DataArrayDouble::isEqualIfNotWhy(const DataArrayDouble& other, double prec, std::string& reason) const
{
  if(!areInfoEqualsIfNotWhy(other,reason))
    return false;
  return _mem.isEqual(other._mem,prec,reason);
}

bool DataArrayDouble::isEqual(const DataArrayDouble& other, double prec) const
{
  std::string tmp;
  return isEqualIfNotWhy(other,prec,tmp);
}

bool DataArrayDouble::isEqualWithoutConsideringStr(const DataArrayDouble& other, double prec) const
{
  std::string tmp;
  return _mem.isEqual(other._mem,prec,tmp);
}

void DataArrayDouble::reAlloc(int nbOfTuples) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  _mem.reAlloc((int)(_info_on_compo.size())*nbOfTuples);
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
 * Only a permutation is done. If a permutation reduction is needed renumberAndReduce.
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
 * Only a permutation is done. If a permutation reduction is needed substr, or selectByTupleId should be used.
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
 * This method is equivalent to DataArrayDouble::renumberAndReduce except that convention in input is new2old and \b not old2new.
 */
DataArrayDouble *DataArrayDouble::selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const
{
  DataArrayDouble *ret=DataArrayDouble::New();
  int nbComp=getNumberOfComponents();
  ret->alloc((int)std::distance(new2OldBg,new2OldEnd),nbComp);
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
 * This method is equivalent to DataArrayDouble::selectByTupleId except that an analyze to the content of input range to check that it will not lead to memory corruption !
 */
DataArrayDouble *DataArrayDouble::selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
  int nbComp=getNumberOfComponents();
  int oldNbOfTuples=getNumberOfTuples();
  ret->alloc((int)std::distance(new2OldBg,new2OldEnd),nbComp);
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
 * Idem than DataArrayInt::selectByTupleIdSafe except that the input array is not constructed explicitely.
 * The convention is as python one. ['bg','end2') with steps of 'step'.
 * Returns a newly created array.
 * This method is a generalization of DataArrayDouble::substr.
 *
 * \sa DataArrayDouble::substr
 */
DataArrayDouble *DataArrayDouble::selectByTupleId2(int bg, int end2, int step) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
  int nbComp=getNumberOfComponents();
  int newNbOfTuples=GetNumberOfItemGivenBES(bg,end2,step,"DataArrayDouble::selectByTupleId2 : ");
  ret->alloc(newNbOfTuples,nbComp);
  double *pt=ret->getPointer();
  const double *srcPt=getConstPointer()+bg*nbComp;
  for(int i=0;i<newNbOfTuples;i++,srcPt+=step*nbComp)
    std::copy(srcPt,srcPt+nbComp,pt+i*nbComp);
  ret->copyStringInfoFrom(*this);
  ret->incrRef();
  return ret;
}

/*!
 * This method returns a newly allocated array that is the concatenation of all tuples ranges in param 'ranges'.
 * Each pair in input 'ranges' is in [begin,end) format. If there is a range in 'ranges' so that end is before begin an exception
 * will be thrown. If there is a range in 'ranges' so that end is greater than number of tuples of 'this', an exception will be thrown too.
 */
DataArrayDouble *DataArrayDouble::selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfComp=getNumberOfComponents();
  int nbOfTuplesThis=getNumberOfTuples();
  if(ranges.empty())
    {
      DataArrayDouble *ret=DataArrayDouble::New();
      ret->alloc(0,nbOfComp);
      ret->copyStringInfoFrom(*this);
      return ret;
    }
  int ref=ranges.front().first;
  int nbOfTuples=0;
  bool isIncreasing=true;
  for(std::vector<std::pair<int,int> >::const_iterator it=ranges.begin();it!=ranges.end();it++)
    {
      if((*it).first<=(*it).second)
        {
          if((*it).first>=0 && (*it).second<=nbOfTuplesThis)
            {
              nbOfTuples+=(*it).second-(*it).first;
              if(isIncreasing)
                isIncreasing=ref<=(*it).first;
              ref=(*it).second;
            }
          else
            {
              std::ostringstream oss; oss << "DataArrayDouble::selectByTupleRanges : on range #" << std::distance(ranges.begin(),it);
              oss << " (" << (*it).first << "," << (*it).second << ") is greater than number of tuples of this :" << nbOfTuples << " !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayDouble::selectByTupleRanges : on range #" << std::distance(ranges.begin(),it);
          oss << " (" << (*it).first << "," << (*it).second << ") end is before begin !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  if(isIncreasing && nbOfTuplesThis==nbOfTuples)
    return deepCpy();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
  ret->alloc(nbOfTuples,nbOfComp);
  ret->copyStringInfoFrom(*this);
  const double *src=getConstPointer();
  double *work=ret->getPointer();
  for(std::vector<std::pair<int,int> >::const_iterator it=ranges.begin();it!=ranges.end();it++)
    work=std::copy(src+(*it).first*nbOfComp,src+(*it).second*nbOfComp,work);
  ret->incrRef();
  return ret;
}

/*!
 * This methods has a similar behaviour than std::string::substr. This method returns a newly created DataArrayInt that is part of this with same number of components.
 * The intervall is specified by [tupleIdBg,tupleIdEnd) except if tupleIdEnd ==-1 in this case the [tupleIdBg,this->end()) will be kept.
 * This method check that interval is valid regarding this, if not an exception will be thrown.
 * This method is a specialization of method DataArrayDouble::selectByTupleId2.
 *
 * \sa DataArrayDouble::selectByTupleId2
 */
DataArrayDouble *DataArrayDouble::substr(int tupleIdBg, int tupleIdEnd) const throw(INTERP_KERNEL::Exception)
{
  int nbt=getNumberOfTuples();
  if(tupleIdBg<0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::substr : The tupleIdBg parameter must be greater than 0 !");
  if(tupleIdBg>nbt)
    throw INTERP_KERNEL::Exception("DataArrayDouble::substr : The tupleIdBg parameter is greater than number of tuples !");
  int trueEnd=tupleIdEnd;
  if(tupleIdEnd!=-1)
    {
      if(tupleIdEnd>nbt)
        throw INTERP_KERNEL::Exception("DataArrayDouble::substr : The tupleIdBg parameter is greater or equal than number of tuples !");
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

/*!
 * This method makes the assumption that \b this is allocated. If not an INTERP_KERNEL::Exception will be raised.
 * This method does not echange the values stored in \b this. Simply, the number of components before the call becomes the number of
 * tuples and inversely the number of tuples becomes the number of components. \b WARNING the info on components can be alterated by this method.
 */
void DataArrayDouble::transpose() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfTuples=getNumberOfTuples();
  rearrange(nbOfTuples);
}

DataArrayDouble *DataArrayDouble::keepSelectedComponents(const std::vector<int>& compoIds) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret(DataArrayDouble::New());
  std::size_t newNbOfCompo=compoIds.size();
  int oldNbOfCompo=getNumberOfComponents();
  for(std::vector<int>::const_iterator it=compoIds.begin();it!=compoIds.end();it++)
    if((*it)<0 || (*it)>=oldNbOfCompo)
      {
        std::ostringstream oss; oss << "DataArrayDouble::keepSelectedComponents : invalid requested component : " << *it << " whereas it should be in [0," << oldNbOfCompo << ") !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  int nbOfTuples=getNumberOfTuples();
  ret->alloc(nbOfTuples,(int)newNbOfCompo);
  ret->copyPartOfStringInfoFrom(*this,compoIds);
  const double *oldc=getConstPointer();
  double *nc=ret->getPointer();
  for(int i=0;i<nbOfTuples;i++)
    for(std::size_t j=0;j<newNbOfCompo;j++,nc++)
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

/*!
 * This methods searches for each tuple if there are any tuples in 'this' that are less far than 'prec' from n1. if any, these tuples are stored in out params
 * comm and commIndex. The distance is computed using norm2. This method expects that 'this' is allocated and that the number of components is in [1,2,3].
 * If not an exception will be thrown.
 * This method is typically used by MEDCouplingPointSet::findCommonNodes and MEDCouplingUMesh::mergeNodes.
 * In case of success, commIndex->getNumberOfTuples()-1 gives the number of tuples groupes that are within distance 'prec'.
 * comm->getNumberOfTuples()==commIndex->back()
 * The returned pair of DataArrayInt instances ('comm','commIndex') is called Surjectived Format 2 \sa DataArrayInt::BuildNew2OldArrayFromSurjectiveFormat2.
 * This format is more compact in surjective format because only all tuple ids not in 'comm' are remain unchanged.
 * 
 * @param prec is an absolute precision.
 * @param limitTupleId is the limit tuple id. All tuples which id is strictly lower than 'limiTupleId' will not be merged each other.
 * @param comm out parameter (not inout). Number of components is equal to 1.
 * @param commIndex out parameter (not inout). Number of components is equal to 1.
 */
void DataArrayDouble::findCommonTuples(double prec, int limitTupleId, DataArrayInt *&comm, DataArrayInt *&commIndex) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfCompo=getNumberOfComponents();
  if ((nbOfCompo<1) || (nbOfCompo>3)) //test before work
    throw INTERP_KERNEL::Exception("DataArrayDouble::findCommonTuples : Unexpected spacedim of coords. Must be 1, 2 or 3.");
  
  int nbOfTuples=getNumberOfTuples();
  comm=DataArrayInt::New();
  commIndex=DataArrayInt::New();
  //
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> bbox=computeBBoxPerTuple(prec);
  //
  std::vector<int> c,cI(1);
  switch(nbOfCompo)
    {
    case 3:
      findCommonTuplesAlg<3>(bbox->getConstPointer(),nbOfTuples,limitTupleId,prec,c,cI);
      break;
    case 2:
      findCommonTuplesAlg<2>(bbox->getConstPointer(),nbOfTuples,limitTupleId,prec,c,cI);
      break;
    case 1:
      findCommonTuplesAlg<1>(bbox->getConstPointer(),nbOfTuples,limitTupleId,prec,c,cI);
      break;
    default:
      throw INTERP_KERNEL::Exception("DataArrayDouble::findCommonTuples : nb of components managed are 1,2 and 3 ! not implemented for other number of components !");
    }
  commIndex->alloc((int)cI.size(),1);
  std::copy(cI.begin(),cI.end(),commIndex->getPointer());
  comm->alloc(cI.back(),1);
  std::copy(c.begin(),c.end(),comm->getPointer());
}

/*!
 * This method returns a newly allocated object the user should deal with.
 * This method works for arrays which have number of components into [1,2,3]. If not an exception will be thrown.
 * This method returns the different values in 'this' using 'prec'. The different values are kept in the same
 * order than 'this'. That is to say that returned DataArrayDouble instance is not systematically sorted.
 *
 * @param prec is an absolute precision.
 * @param limitTupleId is the limit tuple id. All tuples which id is strictly lower than 'limiTupleId' will not be merged each other.
 */
DataArrayDouble *DataArrayDouble::getDifferentValues(double prec, int limitTupleId) const throw(INTERP_KERNEL::Exception)
{
  DataArrayInt *c0=0,*cI0=0;
  findCommonTuples(prec,limitTupleId,c0,cI0);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> c(c0),cI(cI0);
  int newNbOfTuples=-1;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> o2n=DataArrayInt::BuildOld2NewArrayFromSurjectiveFormat2(getNumberOfTuples(),c0,cI0,newNbOfTuples);
  return renumberAndReduce(o2n->getConstPointer(),newNbOfTuples);
}

void DataArrayDouble::setSelectedComponents(const DataArrayDouble *a, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayDouble::setSelectedComponents : input DataArrayDouble is NULL !");
  copyPartOfStringInfoFrom2(compoIds,*a);
  std::size_t partOfCompoSz=compoIds.size();
  int nbOfCompo=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  const double *ac=a->getConstPointer();
  double *nc=getPointer();
  for(int i=0;i<nbOfTuples;i++)
    for(std::size_t j=0;j<partOfCompoSz;j++,ac++)
      nc[nbOfCompo*i+compoIds[j]]=*ac;
}

/*!
 * This method performs a partial assignment of 'this' using 'a' as input. Other input parameters specifies the subpart being considered by the assignment.
 * 'strictCompoCompare' specifies if DataArray 'a' should have exactly same number of components and tuples than 'this' (true) or not (false). By default set to true with maximal test.
 */
void DataArrayDouble::setPartOfValues1(const DataArrayDouble *a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare) throw(INTERP_KERNEL::Exception)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayDouble::setPartOfValues1 : input DataArrayDouble is NULL !");
  const char msg[]="DataArrayDouble::setPartOfValues1";
  checkAllocated();
  a->checkAllocated();
  int newNbOfTuples=DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg);
  int newNbOfComp=DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg);
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
  DataArray::CheckValueInRangeEx(nbComp,bgComp,endComp,"invalid component value");
  a->checkNbOfElems(newNbOfTuples*newNbOfComp,msg);
  if(strictCompoCompare)
    a->checkNbOfTuplesAndComp(newNbOfTuples,newNbOfComp,msg);
  double *pt=getPointer()+bgTuples*nbComp+bgComp;
  const double *srcPt=a->getConstPointer();
  for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
    for(int j=0;j<newNbOfComp;j++,srcPt++)
      pt[j*stepComp]=*srcPt;
}

/*!
 * This method performs a partial assignment of 'this' using 'a' as input. Other input parameters specifies the subpart being considered by the assignment.
 */
void DataArrayDouble::setPartOfValuesSimple1(double a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="DataArrayDouble::setPartOfValuesSimple1";
  checkAllocated();
  int newNbOfTuples=DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg);
  int newNbOfComp=DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg);
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
  DataArray::CheckValueInRangeEx(nbComp,bgComp,endComp,"invalid component value");
  double *pt=getPointer()+bgTuples*nbComp+bgComp;
  for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
    for(int j=0;j<newNbOfComp;j++)
      pt[j*stepComp]=a;
}

/*!
 * This method performs a partial assignment of 'this' using 'a' as input. Other input parameters specifies the subpart being considered by the assignment.
 * 'strictCompoCompare' specifies if DataArray 'a' should have exactly same number of components and tuples than 'this' (true) or not (false). By default set to true with maximal test.
 */
void DataArrayDouble::setPartOfValues2(const DataArrayDouble *a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp, bool strictCompoCompare) throw(INTERP_KERNEL::Exception)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayDouble::setPartOfValues2 : input DataArrayDouble is NULL !");
  const char msg[]="DataArrayDouble::setPartOfValues2";
  checkAllocated();
  a->checkAllocated();
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  for(const int *z=bgComp;z!=endComp;z++)
    DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
  int newNbOfTuples=(int)std::distance(bgTuples,endTuples);
  int newNbOfComp=(int)std::distance(bgComp,endComp);
  bool assignTech=true;
  if(a->getNbOfElems()==newNbOfTuples*newNbOfComp)
    {
      if(strictCompoCompare)
        a->checkNbOfTuplesAndComp(newNbOfTuples,newNbOfComp,msg);
    }
  else
    {
      a->checkNbOfTuplesAndComp(1,newNbOfComp,msg);
      assignTech=false;
    }
  double *pt=getPointer();
  const double *srcPt=a->getConstPointer();
  if(assignTech)
    {    
      for(const int *w=bgTuples;w!=endTuples;w++)
        {
          DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
          for(const int *z=bgComp;z!=endComp;z++,srcPt++)
            {    
              pt[(*w)*nbComp+(*z)]=*srcPt;
            }
        }
    }
  else
    {
      for(const int *w=bgTuples;w!=endTuples;w++)
        {
          const double *srcPt2=srcPt;
          DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
          for(const int *z=bgComp;z!=endComp;z++,srcPt2++)
            {    
              pt[(*w)*nbComp+(*z)]=*srcPt2;
            }
        }
    }
}

/*!
 * This method performs a partial assignment of 'this' using 'a' as input. Other input parameters specifies the subpart being considered by the assignment.
 */
void DataArrayDouble::setPartOfValuesSimple2(double a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  for(const int *z=bgComp;z!=endComp;z++)
    DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
  double *pt=getPointer();
  for(const int *w=bgTuples;w!=endTuples;w++)
    for(const int *z=bgComp;z!=endComp;z++)
      {
        DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
        pt[(*w)*nbComp+(*z)]=a;
      }
}

/*!
 * This method performs a partial assignment of 'this' using 'a' as input. Other input parameters specifies the subpart being considered by the assignment.
 * 'strictCompoCompare' specifies if DataArray 'a' should have exactly same number of components and tuples than 'this' (true) or not (false). By default set to true with maximal test.
 */
void DataArrayDouble::setPartOfValues3(const DataArrayDouble *a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare) throw(INTERP_KERNEL::Exception)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayDouble::setPartOfValues3 : input DataArrayDouble is NULL !");
  const char msg[]="DataArrayDouble::setPartOfValues3";
  checkAllocated();
  a->checkAllocated();
  int newNbOfComp=DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg);
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbComp,bgComp,endComp,"invalid component value");
  int newNbOfTuples=(int)std::distance(bgTuples,endTuples);
  bool assignTech=true;
  if(a->getNbOfElems()==newNbOfTuples*newNbOfComp)
    {
      if(strictCompoCompare)
        a->checkNbOfTuplesAndComp(newNbOfTuples,newNbOfComp,msg);
    }
  else
    {
      a->checkNbOfTuplesAndComp(1,newNbOfComp,msg);
      assignTech=false;
    }
  double *pt=getPointer()+bgComp;
  const double *srcPt=a->getConstPointer();
  if(assignTech)
    {
      for(const int *w=bgTuples;w!=endTuples;w++)
        for(int j=0;j<newNbOfComp;j++,srcPt++)
          {
            DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
            pt[(*w)*nbComp+j*stepComp]=*srcPt;
          }
    }
  else
    {
      for(const int *w=bgTuples;w!=endTuples;w++)
        {
          const double *srcPt2=srcPt;
          for(int j=0;j<newNbOfComp;j++,srcPt2++)
            {
              DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
              pt[(*w)*nbComp+j*stepComp]=*srcPt2;
            }
        }
    }
}

/*!
 * This method performs a partial assignment of 'this' using 'a' as input. Other input parameters specifies the subpart being considered by the assignment.
 */
void DataArrayDouble::setPartOfValuesSimple3(double a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="DataArrayDouble::setPartOfValuesSimple3";
  checkAllocated();
  int newNbOfComp=DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg);
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbComp,bgComp,endComp,"invalid component value");
  double *pt=getPointer()+bgComp;
  for(const int *w=bgTuples;w!=endTuples;w++)
    for(int j=0;j<newNbOfComp;j++)
      {
        DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
        pt[(*w)*nbComp+j*stepComp]=a;
      }
}

/*!
 * 'this', 'a' and 'tuplesSelec' are expected to be defined. If not an exception will be thrown.
 * @param a is an array having exactly the same number of components than 'this'
 * @param tuplesSelec is an array having exactly 2 components. The first one refers to the tuple ids of 'this' that will be set. The second one refers to the tuple ids of 'a' that will be used for setting.
 */
void DataArrayDouble::setPartOfValuesAdv(const DataArrayDouble *a, const DataArrayInt *tuplesSelec) throw(INTERP_KERNEL::Exception)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayDouble::setPartOfValuesAdv : input DataArrayDouble is NULL !");
  checkAllocated();
  a->checkAllocated();
  tuplesSelec->checkAllocated();
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=a->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("DataArrayDouble::setPartOfValuesAdv : This and a do not have the same number of components !");
  if(tuplesSelec->getNumberOfComponents()!=2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::setPartOfValuesAdv : Expecting to have a tuple selector DataArrayInt instance with exactly 2 components !");
  int thisNt=getNumberOfTuples();
  int aNt=a->getNumberOfTuples();
  double *valsToSet=getPointer();
  const double *valsSrc=a->getConstPointer();
  for(const int *tuple=tuplesSelec->begin();tuple!=tuplesSelec->end();tuple+=2)
    {
      if(tuple[1]>=0 && tuple[1]<aNt)
        {
          if(tuple[0]>=0 && tuple[0]<thisNt)
            std::copy(valsSrc+nbOfComp*tuple[1],valsSrc+nbOfComp*(tuple[1]+1),valsToSet+nbOfComp*tuple[0]);
          else
            {
              std::ostringstream oss; oss << "DataArrayDouble::setPartOfValuesAdv : Tuple #" << std::distance(tuplesSelec->begin(),tuple)/2;
              oss << " of 'tuplesSelec' request of tuple id #" << tuple[0] << " in 'this' ! It should be in [0," << thisNt << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayDouble::setPartOfValuesAdv : Tuple #" << std::distance(tuplesSelec->begin(),tuple)/2;
          oss << " of 'tuplesSelec' request of tuple id #" << tuple[1] << " in 'a' ! It should be in [0," << aNt << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}

/*!
 * 'this', 'a' and 'tuplesSelec' are expected to be defined. If not an exception will be thrown.
 * This is a method that is a specialization to DataArrayDouble::setPartOfValuesAdv method, except that here the tuple selection of 'a' is given by a range ('bg','end2' and 'step')
 * rather than an explicite array of tuple ids (given by the 2nd component) and the feeding is done in 'this' contiguously starting from 'tupleIdStart'.
 * @param a is an array having exactly the same number of components than 'this'
 */
void DataArrayDouble::setContigPartOfSelectedValues(int tupleIdStart, const DataArrayDouble *a, const DataArrayInt *tuplesSelec) throw(INTERP_KERNEL::Exception)
{
  if(!a || !tuplesSelec)
    throw INTERP_KERNEL::Exception("DataArrayDouble::setContigPartOfSelectedValues : input DataArray is NULL !");
  checkAllocated();
  a->checkAllocated();
  tuplesSelec->checkAllocated();
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=a->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("DataArrayDouble::setContigPartOfSelectedValues : This and a do not have the same number of components !");
  if(tuplesSelec->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::setContigPartOfSelectedValues : Expecting to have a tuple selector DataArrayInt instance with exactly 1 component !");
  int thisNt=getNumberOfTuples();
  int aNt=a->getNumberOfTuples();
  int nbOfTupleToWrite=tuplesSelec->getNumberOfTuples();
  double *valsToSet=getPointer()+tupleIdStart*nbOfComp;
  if(tupleIdStart+nbOfTupleToWrite>thisNt)
    throw INTERP_KERNEL::Exception("DataArrayDouble::setContigPartOfSelectedValues : invalid number range of values to write !");
  const double *valsSrc=a->getConstPointer();
  for(const int *tuple=tuplesSelec->begin();tuple!=tuplesSelec->end();tuple++,valsToSet+=nbOfComp)
    {
      if(*tuple>=0 && *tuple<aNt)
        {
          std::copy(valsSrc+nbOfComp*(*tuple),valsSrc+nbOfComp*(*tuple+1),valsToSet);
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayDouble::setContigPartOfSelectedValues : Tuple #" << std::distance(tuplesSelec->begin(),tuple);
          oss << " of 'tuplesSelec' request of tuple id #" << *tuple << " in 'a' ! It should be in [0," << aNt << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}

/*!
 * 'this' and 'a' are expected to be defined. If not an exception will be thrown.
 * This is a method that is a specialization to DataArrayDouble::setContigPartOfSelectedValues method, except that here the tuple selection is givenin a is done by a range ('bg','end2' and 'step')
 * rather than an explicite array of tuple ids.
 * @param a is an array having exactly the same number of components than 'this'
 */
void DataArrayDouble::setContigPartOfSelectedValues2(int tupleIdStart, const DataArrayDouble *a, int bg, int end2, int step) throw(INTERP_KERNEL::Exception)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayDouble::setContigPartOfSelectedValues2 : input DataArrayDouble is NULL !");
  checkAllocated();
  a->checkAllocated();
  int nbOfComp=getNumberOfComponents();
  const char msg[]="DataArrayDouble::setContigPartOfSelectedValues2";
  int nbOfTupleToWrite=DataArray::GetNumberOfItemGivenBES(bg,end2,step,msg);
  if(nbOfComp!=a->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("DataArrayDouble::setContigPartOfSelectedValues2 : This and a do not have the same number of components !");
  int thisNt=getNumberOfTuples();
  int aNt=a->getNumberOfTuples();
  double *valsToSet=getPointer()+tupleIdStart*nbOfComp;
  if(tupleIdStart+nbOfTupleToWrite>thisNt)
    throw INTERP_KERNEL::Exception("DataArrayDouble::setContigPartOfSelectedValues2 : invalid number range of values to write !");
  if(end2>aNt)
    throw INTERP_KERNEL::Exception("DataArrayDouble::setContigPartOfSelectedValues2 : invalid range of values to read !");
  const double *valsSrc=a->getConstPointer()+bg*nbOfComp;
  for(int i=0;i<nbOfTupleToWrite;i++,valsToSet+=nbOfComp,valsSrc+=step*nbOfComp)
    {
      std::copy(valsSrc,valsSrc+nbOfComp,valsToSet);
    }
}

/*!
 * This method is equivalent to DataArrayDouble::getIJ except that here \b tupleId is checked to be in [0,this->getNumberOfTuples()) and compoId to be in [0,this->getNumberOfComponents()).
 * If one of these check fails an INTERP_KERNEL::Exception will be thrown.
 * So this method is safe but expensive if used to go through all data of \b this.
 */
double DataArrayDouble::getIJSafe(int tupleId, int compoId) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(tupleId<0 || tupleId>=getNumberOfTuples())
    {
      std::ostringstream oss; oss << "DataArrayDouble::getIJSafe : request for tupleId " << tupleId << " should be in [0," << getNumberOfTuples() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(compoId<0 || compoId>=getNumberOfComponents())
    {
      std::ostringstream oss; oss << "DataArrayDouble::getIJSafe : request for compoId " << compoId << " should be in [0," << getNumberOfComponents() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return _mem[tupleId*((int)_info_on_compo.size())+compoId];
}

/*!
 * This method returns the last element in 'this'. So this method makes the hypothesis that 'this' is allocated.
 * This method works only for arrays that have exactly number of components equal to 1. If not an exception is thrown.
 * And to finish this method works for arrays that have number of tuples >= 1.
 */
double DataArrayDouble::back() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::back : number of components not equal to one !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::back : number of tuples must be >= 1 !");
  return *(getConstPointer()+nbOfTuples-1);
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

/*!
 * This method assume that \b this is allocated. If not an INTERP_KERNEL::Exception will be thrown.
 * This method fills \b bounds params like that : \b bounds[0]=XMin, \b bounds[1]=XMax, \b bounds[2]=YMin, \b bounds[3]=YMax...
 * Where X refers to component #0, and Y to component #1...
 * This method set 2*this->getNumberOfComponents() elements in \b bounds, so it is up to the caller to allocated enough space before calling this method.
 *
 * @param [out] bounds array of size 2*this->getNumberOfComponents().
 */
void DataArrayDouble::getMinMaxPerComponent(double *bounds) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int dim=getNumberOfComponents();
  for (int idim=0; idim<dim; idim++)
    {
      bounds[idim*2]=std::numeric_limits<double>::max();
      bounds[idim*2+1]=-std::numeric_limits<double>::max();
    } 
  const double *ptr=getConstPointer();
  int nbOfTuples=getNumberOfTuples();
  for(int i=0;i<nbOfTuples;i++)
    {
      for(int idim=0;idim<dim;idim++)
        {
          if(bounds[idim*2]>ptr[i*dim+idim])
            {
              bounds[idim*2]=ptr[i*dim+idim];
            }
          if(bounds[idim*2+1]<ptr[i*dim+idim])
            {
              bounds[idim*2+1]=ptr[i*dim+idim];
            }
        }
    }
}

/*!
 * This method retrieves a newly allocated DataArrayDouble instance having same number of tuples than \a this and twice number of components than \a this
 * to store both the min and max per component of each tuples. 
 * \param [in] epsilon the width of the bbox (identical in each direction) - 0.0 by default
 *
 * \return a newly created DataArrayDouble instance having \c this->getNumberOfTuples() tuples and 2 * \c this->getNumberOfComponent() components
 *
 * \throw If \a this is not allocated yet.
 */
DataArrayDouble *DataArrayDouble::computeBBoxPerTuple(double epsilon)const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  const double *dataPtr=getConstPointer();
  int nbOfCompo=getNumberOfComponents();
  int nbTuples=getNumberOfTuples();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> bbox=DataArrayDouble::New();
  bbox->alloc(nbTuples,2*nbOfCompo);
  double *bboxPtr=bbox->getPointer();
  for(int i=0;i<nbTuples;i++)
    {
      for(int j=0;j<nbOfCompo;j++)
        {
          bboxPtr[2*nbOfCompo*i+2*j]=dataPtr[nbOfCompo*i+j]-epsilon;
          bboxPtr[2*nbOfCompo*i+2*j+1]=dataPtr[nbOfCompo*i+j]+epsilon;
        }
    }
  bbox->incrRef();
  return bbox;
}

/*!
 * For each tuples **t** in \a other, this method retrieves tuples in \a this that are equal to **t**.
 * Two tuples are considered equal if the euclidian distance between the two tuples is lower than \a eps.
 * 
 * \param [in] other a DataArrayDouble having same number of components than \a this.
 * \param [in] eps absolute precision representing euclidian distance between 2 tuples behind which 2 tuples are considered equal.
 * \param [out] c will contain the set of tuple ids in \a this that are equal to to the tuple ids in \a other contiguously.
 *             \a cI allows to extract information in \a c.
 * \param [out] cI is an indirection array that allows to extract the data contained in \a c.
 *
 * \throw In case of:
 *  - \a this is not allocated
 *  - \a other is not allocated or null
 *  - \a this and \a other do not have the same number of components
 *  - if number of components of \a this is not in [1,2,3]
 *
 * \sa MEDCouplingPointSet::getNodeIdsNearPoints, DataArrayDouble::getDifferentValues
 */
void DataArrayDouble::computeTupleIdsNearTuples(const DataArrayDouble *other, double eps, std::vector<int>& c, std::vector<int>& cI) const throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::computeTupleIdsNearTuples : input pointer other is null !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> bbox=computeBBoxPerTuple(eps);
  other->checkAllocated();
  int nbOfCompo=getNumberOfComponents();
  int otherNbOfCompo=other->getNumberOfComponents();
  if(nbOfCompo!=otherNbOfCompo)
    throw INTERP_KERNEL::Exception("DataArrayDouble::computeTupleIdsNearTuples : number of components should be equal between this and other !");
  int nbOfTuplesOther=other->getNumberOfTuples();
  std::vector<int> ret;
  c.clear();
  cI.resize(1); cI[0]=0;
  switch(nbOfCompo)
    {
    case 3:
      {
        BBTree<3,int> myTree(bbox->getConstPointer(),0,0,getNumberOfTuples(),eps/10);
        findTupleIdsNearTuplesAlg<3>(myTree,other->getConstPointer(),nbOfTuplesOther,eps,c,cI);
        break;
      }
    case 2:
      {
        BBTree<2,int> myTree(bbox->getConstPointer(),0,0,getNumberOfTuples(),eps/10);
        findTupleIdsNearTuplesAlg<2>(myTree,other->getConstPointer(),nbOfTuplesOther,eps,c,cI);
        break;
      }
    case 1:
      {
        BBTree<1,int> myTree(bbox->getConstPointer(),0,0,getNumberOfTuples(),eps/10);
        findTupleIdsNearTuplesAlg<1>(myTree,other->getConstPointer(),nbOfTuplesOther,eps,c,cI);
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("Unexpected spacedim of coords for computeTupleIdsNearTuples. Must be 1, 2 or 3.");
    }
}

/*!
 * This method recenter tuples in \b this in order to be centered at the origin to benefit about the advantages of maximal precision to be around the box
 * around origin of 'radius' 1.
 * 
 * \param [in] eps absolute epsilon. under that value of delta between max and min no scale is performed.
 */
void DataArrayDouble::recenterForMaxPrecision(double eps) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int dim=getNumberOfComponents();
  std::vector<double> bounds(2*dim);
  getMinMaxPerComponent(&bounds[0]);
  for(int i=0;i<dim;i++)
    {
      double delta=bounds[2*i+1]-bounds[2*i];
      double offset=(bounds[2*i]+bounds[2*i+1])/2.;
      if(delta>eps)
        applyLin(1./delta,-offset/delta,i);
      else
        applyLin(1.,-offset,i);
    }
}

double DataArrayDouble::getMaxValue(int& tupleId) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getMaxValue : must be applied on DataArrayDouble with only one component, you can call 'rearrange' method before or call 'getMaxValueInArray' method !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<=0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getMaxValue : array exists but number of tuples must be > 0 !");
  const double *vals=getConstPointer();
  const double *loc=std::max_element(vals,vals+nbOfTuples);
  tupleId=(int)std::distance(vals,loc);
  return *loc;
}

/*!
 * Idem to DataArrayDouble::getMaxValue expect that here number of components can be >=1.
 */
double DataArrayDouble::getMaxValueInArray() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  const double *loc=std::max_element(begin(),end());
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
    throw INTERP_KERNEL::Exception("DataArrayDouble::getMinValue : must be applied on DataArrayDouble with only one component, you can call 'rearrange' method before call 'getMinValueInArray' method !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<=0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getMinValue : array exists but number of tuples must be > 0 !");
  const double *vals=getConstPointer();
  const double *loc=std::min_element(vals,vals+nbOfTuples);
  tupleId=(int)std::distance(vals,loc);
  return *loc;
}

/*!
 * Idem to DataArrayDouble::getMinValue expect that here number of components can be >=1.
 */
double DataArrayDouble::getMinValueInArray() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  const double *loc=std::min_element(begin(),end());
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
    throw INTERP_KERNEL::Exception("DataArrayDouble::getAverageValue : must be applied on DataArrayDouble with only one component, you can call 'rearrange' method before !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<=0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getAverageValue : array exists but number of tuples must be > 0 !");
  const double *vals=getConstPointer();
  double ret=std::accumulate(vals,vals+nbOfTuples,0.);
  return ret/nbOfTuples;
}

double DataArrayDouble::norm2() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  double ret=0.;
  int nbOfElems=getNbOfElems();
  const double *pt=getConstPointer();
  for(int i=0;i<nbOfElems;i++,pt++)
    ret+=(*pt)*(*pt);
  return sqrt(ret);
}

double DataArrayDouble::normMax() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  double ret=-1.;
  int nbOfElems=getNbOfElems();
  const double *pt=getConstPointer();
  for(int i=0;i<nbOfElems;i++,pt++)
    {
      double val=std::abs(*pt);
      if(val>ret)
        ret=val;
    }
  return ret;
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

/*!
 * This method returns a newly allocated DataArrayDouble instance having one component and \c this->getNumberOfTuples() * \c this->getNumberOfTuples() tuples.
 * \n This returned array contains the euclidian distance for each tuple in \a this. 
 * \n So the returned array can be seen as a dense symmetrical matrix whose diagonal elements are equal to 0.
 * \n The returned array has only one component (and **not** \c this->getNumberOfTuples() components to avoid the useless memory consumption due to components info in returned DataArrayDouble)
 *
 * \warning use this method with care because it can leads to big amount of consumed memory !
 * 
 * \return A newly allocated (huge) ParaMEDMEM::DataArrayDouble instance that the caller should deal with.
 *
 * \throw If \a this is not allocated.
 *
 * \sa DataArrayDouble::buildEuclidianDistanceDenseMatrixWith
 */
DataArrayDouble *DataArrayDouble::buildEuclidianDistanceDenseMatrix() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  const double *inData=getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
  ret->alloc(nbOfTuples*nbOfTuples,1);
  double *outData=ret->getPointer();
  for(int i=0;i<nbOfTuples;i++)
    {
      outData[i*nbOfTuples+i]=0.;
      for(int j=i+1;j<nbOfTuples;j++)
        {
          double dist=0.;
          for(int k=0;k<nbOfComp;k++)
            { double delta=inData[i*nbOfComp+k]-inData[j*nbOfComp+k]; dist+=delta*delta; }
          dist=sqrt(dist);
          outData[i*nbOfTuples+j]=dist;
          outData[j*nbOfTuples+i]=dist;
        }
    }
  ret->incrRef();
  return ret;
}

/*!
 * This method returns a newly allocated DataArrayDouble instance having one component and \c this->getNumberOfTuples() * \c other->getNumberOfTuples() tuples.
 * \n This returned array contains the euclidian distance for each tuple in \a other with each tuple in \a this. 
 * \n So the returned array can be seen as a dense rectangular matrix with \c other->getNumberOfTuples() rows and \c this->getNumberOfTuples() columns.
 * \n Output rectangular matrix is sorted along rows.
 * \n The returned array has only one component (and **not** \c this->getNumberOfTuples() components to avoid the useless memory consumption due to components info in returned DataArrayDouble)
 *
 * \warning use this method with care because it can leads to big amount of consumed memory !
 * 
 * \param [in] other DataArrayDouble instance having same number of components than \a this.
 * \return A newly allocated (huge) ParaMEDMEM::DataArrayDouble instance that the caller should deal with.
 *
 * \throw If \a this is not allocated, or if \a other is null or if \a other is not allocated, or if number of components of \a other and \a this differs.
 *
 * \sa DataArrayDouble::buildEuclidianDistanceDenseMatrix
 */
DataArrayDouble *DataArrayDouble::buildEuclidianDistanceDenseMatrixWith(const DataArrayDouble *other) const throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::buildEuclidianDistanceDenseMatrixWith : input parameter is null !");
  checkAllocated();
  other->checkAllocated();
  int nbOfComp=getNumberOfComponents();
  int otherNbOfComp=other->getNumberOfComponents();
  if(nbOfComp!=otherNbOfComp)
    {
      std::ostringstream oss; oss << "DataArrayDouble::buildEuclidianDistanceDenseMatrixWith : this nb of compo=" << nbOfComp << " and other nb of compo=" << otherNbOfComp << ". It should match !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  int nbOfTuples=getNumberOfTuples();
  int otherNbOfTuples=other->getNumberOfTuples();
  const double *inData=getConstPointer();
  const double *inDataOther=other->getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
  ret->alloc(otherNbOfTuples*nbOfTuples,1);
  double *outData=ret->getPointer();
  for(int i=0;i<otherNbOfTuples;i++,inDataOther+=nbOfComp)
    {
      for(int j=0;j<nbOfTuples;j++)
        {
          double dist=0.;
          for(int k=0;k<nbOfComp;k++)
            { double delta=inDataOther[k]-inData[j*nbOfComp+k]; dist+=delta*delta; }
          dist=sqrt(dist);
          outData[i*nbOfTuples+j]=dist;
        }
    }
  ret->incrRef();
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

void DataArrayDouble::abs() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  double *ptr=getPointer();
  int nbOfElems=getNbOfElems();
  std::transform(ptr,ptr+nbOfElems,ptr,std::ptr_fun<double,double>(fabs));
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

void DataArrayDouble::applyLin(double a, double b) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  double *ptr=getPointer();
  int nbOfElems=getNbOfElems();
  for(int i=0;i<nbOfElems;i++,ptr++)
    *ptr=a*(*ptr)+b;
  declareAsNew();
}

/*!
 * This method applies the operation 'numerator/x' for each element 'x' in 'this'.
 * If there is a value in 'this' exactly equal to 0. an exception is thrown.
 * Warning if presence of null this is modified for each values previous than place where exception was thrown !
 */
void DataArrayDouble::applyInv(double numerator) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  double *ptr=getPointer();
  int nbOfElems=getNbOfElems();
  for(int i=0;i<nbOfElems;i++,ptr++)
    {
      if(std::abs(*ptr)>std::numeric_limits<double>::min())
        {
          *ptr=numerator/(*ptr);
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayDouble::applyInv : presence of null value in tuple #" << i/getNumberOfComponents() << " component #" << i%getNumberOfComponents();
          oss << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  declareAsNew();
}

/*!
 * This method returns a newly allocated array containing the application of negate on \b this.
 * This method throws an INTERP_KERNEL::Exception if \b this is not allocated.
 */
DataArrayDouble *DataArrayDouble::negate() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  DataArrayDouble *newArr=DataArrayDouble::New();
  int nbOfTuples=getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  newArr->alloc(nbOfTuples,nbOfComp);
  const double *cptr=getConstPointer();
  std::transform(cptr,cptr+nbOfTuples*nbOfComp,newArr->getPointer(),std::negate<double>());
  newArr->copyStringInfoFrom(*this);
  return newArr;
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

/*!
 * This method returns a newly allocated array the caller should deal with.
 * The returned array will have 'nbOfComp' components (that can be different from this->getNumberOfComponents()) contrary to the other DataArrayDouble::applyFunc overload method.
 */
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
      std::ostringstream oss; oss << "The field has " << oldNbOfComp << " components and there are ";
      oss << vars.size() << " variables : ";
      std::copy(vars.begin(),vars.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::vector<std::string> varsV(vars.begin(),vars.end());
  expr.prepareExprEvaluation(varsV,oldNbOfComp,nbOfComp);
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

/*!
 * This method is equivalent than DataArrayDouble::applyFunc, except that here components names are used to determine vars orders.
 * If 'func' contains vars that are not in \c this->getInfoOnComponent() an exception will be thrown.
 */
DataArrayDouble *DataArrayDouble::applyFunc2(int nbOfComp, const char *func) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  std::set<std::string> vars;
  expr.getTrueSetOfVars(vars);
  int oldNbOfComp=getNumberOfComponents();
  if((int)vars.size()>oldNbOfComp)
    {
      std::ostringstream oss; oss << "The field has " << oldNbOfComp << " components and there are ";
      oss << vars.size() << " variables : ";
      std::copy(vars.begin(),vars.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  expr.prepareExprEvaluation(getVarsOnComponent(),oldNbOfComp,nbOfComp);
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

/*!
 * This method is equivalent than DataArrayDouble::applyFunc, except that here order of vars is passed explicitely in parameter.
 * In 'func' contains vars not in 'varsOrder' an exception will be thrown.
 */
DataArrayDouble *DataArrayDouble::applyFunc3(int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  std::set<std::string> vars;
  expr.getTrueSetOfVars(vars);
  int oldNbOfComp=getNumberOfComponents();
  if((int)vars.size()>oldNbOfComp)
    {
      std::ostringstream oss; oss << "The field has " << oldNbOfComp << " components and there are ";
      oss << vars.size() << " variables : ";
      std::copy(vars.begin(),vars.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  expr.prepareExprEvaluation(varsOrder,oldNbOfComp,nbOfComp);
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

void DataArrayDouble::applyFuncFast32(const char *func) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  char *funcStr=expr.compileX86();
  MYFUNCPTR funcPtr;
  *((void **)&funcPtr)=funcStr;//he he...
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
  MYFUNCPTR funcPtr;
  *((void **)&funcPtr)=funcStr;//he he...
  //
  double *ptr=getPointer();
  int nbOfComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  int nbOfElems=nbOfTuples*nbOfComp;
  for(int i=0;i<nbOfElems;i++,ptr++)
    *ptr=funcPtr(*ptr);
  declareAsNew();
}

DataArrayDoubleIterator *DataArrayDouble::iterator()
{
  return new DataArrayDoubleIterator(this);
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
  ret->alloc((int)res.size(),1);
  std::copy(res.begin(),res.end(),ret->getPointer());
  return ret;
}

DataArrayDouble *DataArrayDouble::Aggregate(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayDouble *> tmp(2);
  tmp[0]=a1; tmp[1]=a2;
  return Aggregate(tmp);
}

DataArrayDouble *DataArrayDouble::Aggregate(const std::vector<const DataArrayDouble *>& arr) throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayDouble *> a;
  for(std::vector<const DataArrayDouble *>::const_iterator it4=arr.begin();it4!=arr.end();it4++)
    if(*it4)
      a.push_back(*it4);
  if(a.empty())
    throw INTERP_KERNEL::Exception("DataArrayDouble::Aggregate : input list must contain at least one NON EMPTY DataArrayDouble !");
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

DataArrayDouble *DataArrayDouble::Meld(const std::vector<const DataArrayDouble *>& arr) throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayDouble *> a;
  for(std::vector<const DataArrayDouble *>::const_iterator it4=arr.begin();it4!=arr.end();it4++)
    if(*it4)
      a.push_back(*it4);
  if(a.empty())
    throw INTERP_KERNEL::Exception("DataArrayDouble::Meld : input list must contain at least one NON EMPTY DataArrayDouble !");
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
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::Dot : input DataArrayDouble instance is NULL !");
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
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::CrossProduct : input DataArrayDouble instance is NULL !");
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
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::Max : input DataArrayDouble instance is NULL !");
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
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::Min : input DataArrayDouble instance is NULL !");
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
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::Add : input DataArrayDouble instance is NULL !");
  int nbOfTuple=a1->getNumberOfTuples();
  int nbOfTuple2=a2->getNumberOfTuples();
  int nbOfComp=a1->getNumberOfComponents();
  int nbOfComp2=a2->getNumberOfComponents();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=0;
  if(nbOfTuple==nbOfTuple2)
    {
      if(nbOfComp==nbOfComp2)
        {
          ret=DataArrayDouble::New();
          ret->alloc(nbOfTuple,nbOfComp);
          std::transform(a1->begin(),a1->end(),a2->begin(),ret->getPointer(),std::plus<double>());
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
                res=std::transform(aMaxPtr+i*nbOfCompMax,aMaxPtr+(i+1)*nbOfCompMax,res,std::bind2nd(std::plus<double>(),aMinPtr[i]));
              ret->copyStringInfoFrom(*aMax);
            }
          else
            throw INTERP_KERNEL::Exception("Nb of components mismatch for array Add !");
        }
    }
  else if((nbOfTuple==1 && nbOfTuple2>1) || (nbOfTuple>1 && nbOfTuple2==1))
    {
      if(nbOfComp==nbOfComp2)
        {
          int nbOfTupleMax=std::max(nbOfTuple,nbOfTuple2);
          const DataArrayDouble *aMin=nbOfTuple>nbOfTuple2?a2:a1;
          const DataArrayDouble *aMax=nbOfTuple>nbOfTuple2?a1:a2;
          const double *aMinPtr=aMin->getConstPointer(),*aMaxPtr=aMax->getConstPointer();
          ret=DataArrayDouble::New();
          ret->alloc(nbOfTupleMax,nbOfComp);
          double *res=ret->getPointer();
          for(int i=0;i<nbOfTupleMax;i++)
            res=std::transform(aMaxPtr+i*nbOfComp,aMaxPtr+(i+1)*nbOfComp,aMinPtr,res,std::plus<double>());
          ret->copyStringInfoFrom(*aMax);
        }
      else
        throw INTERP_KERNEL::Exception("Nb of components mismatch for array Add !");
    }
  else
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array Add !");
  ret->incrRef();
  return ret;
}

void DataArrayDouble::addEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::addEqual : input DataArrayDouble instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayDouble::addEqual  !";
  int nbOfTuple=getNumberOfTuples();
  int nbOfTuple2=other->getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  int nbOfComp2=other->getNumberOfComponents();
  if(nbOfTuple==nbOfTuple2)
    {
      if(nbOfComp==nbOfComp2)
        {
          std::transform(begin(),end(),other->begin(),getPointer(),std::plus<double>());
        }
      else if(nbOfComp2==1)
        {
          double *ptr=getPointer();
          const double *ptrc=other->getConstPointer();
          for(int i=0;i<nbOfTuple;i++)
            std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptr+i*nbOfComp,std::bind2nd(std::plus<double>(),*ptrc++));
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else if(nbOfTuple2==1)
    {
      if(nbOfComp2==nbOfComp)
        {
          double *ptr=getPointer();
          const double *ptrc=other->getConstPointer();
          for(int i=0;i<nbOfTuple;i++)
            std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptrc,ptr+i*nbOfComp,std::plus<double>());
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else
    throw INTERP_KERNEL::Exception(msg);
  declareAsNew();
}

DataArrayDouble *DataArrayDouble::Substract(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::Substract : input DataArrayDouble instance is NULL !");
  int nbOfTuple1=a1->getNumberOfTuples();
  int nbOfTuple2=a2->getNumberOfTuples();
  int nbOfComp1=a1->getNumberOfComponents();
  int nbOfComp2=a2->getNumberOfComponents();
  if(nbOfTuple2==nbOfTuple1)
    {
      if(nbOfComp1==nbOfComp2)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
          ret->alloc(nbOfTuple2,nbOfComp1);
          std::transform(a1->begin(),a1->end(),a2->begin(),ret->getPointer(),std::minus<double>());
          ret->copyStringInfoFrom(*a1);
          ret->incrRef();
          return ret;
        }
      else if(nbOfComp2==1)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
          ret->alloc(nbOfTuple1,nbOfComp1);
          const double *a2Ptr=a2->getConstPointer();
          const double *a1Ptr=a1->getConstPointer();
          double *res=ret->getPointer();
          for(int i=0;i<nbOfTuple1;i++)
            res=std::transform(a1Ptr+i*nbOfComp1,a1Ptr+(i+1)*nbOfComp1,res,std::bind2nd(std::minus<double>(),a2Ptr[i]));
          ret->copyStringInfoFrom(*a1);
          ret->incrRef();
          return ret;
        }
      else
        {
          a1->checkNbOfComps(nbOfComp2,"Nb of components mismatch for array Substract !");
          return 0;
        }
    }
  else if(nbOfTuple2==1)
    {
      a1->checkNbOfComps(nbOfComp2,"Nb of components mismatch for array Substract !");
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
      ret->alloc(nbOfTuple1,nbOfComp1);
      const double *a1ptr=a1->getConstPointer(),*a2ptr=a2->getConstPointer();
      double *pt=ret->getPointer();
      for(int i=0;i<nbOfTuple1;i++)
        pt=std::transform(a1ptr+i*nbOfComp1,a1ptr+(i+1)*nbOfComp1,a2ptr,pt,std::minus<double>());
      ret->copyStringInfoFrom(*a1);
      ret->incrRef();
      return ret;
    }
  else
    {
      a1->checkNbOfTuples(nbOfTuple2,"Nb of tuples mismatch for array Substract !");//will always throw an exception
      return 0;
    }
}

void DataArrayDouble::substractEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::substractEqual : input DataArrayDouble instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayDouble::substractEqual  !";
  int nbOfTuple=getNumberOfTuples();
  int nbOfTuple2=other->getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  int nbOfComp2=other->getNumberOfComponents();
  if(nbOfTuple==nbOfTuple2)
    {
      if(nbOfComp==nbOfComp2)
        {
          std::transform(begin(),end(),other->begin(),getPointer(),std::minus<double>());
        }
      else if(nbOfComp2==1)
        {
          double *ptr=getPointer();
          const double *ptrc=other->getConstPointer();
          for(int i=0;i<nbOfTuple;i++)
            std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptr+i*nbOfComp,std::bind2nd(std::minus<double>(),*ptrc++)); 
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else if(nbOfTuple2==1)
    {
      if(nbOfComp2==nbOfComp)
        {
          double *ptr=getPointer();
          const double *ptrc=other->getConstPointer();
          for(int i=0;i<nbOfTuple;i++)
            std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptrc,ptr+i*nbOfComp,std::minus<double>());
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else
    throw INTERP_KERNEL::Exception(msg);
  declareAsNew();
}

DataArrayDouble *DataArrayDouble::Multiply(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::Multiply : input DataArrayDouble instance is NULL !");
  int nbOfTuple=a1->getNumberOfTuples();
  int nbOfTuple2=a2->getNumberOfTuples();
  int nbOfComp=a1->getNumberOfComponents();
  int nbOfComp2=a2->getNumberOfComponents();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=0;
  if(nbOfTuple==nbOfTuple2)
    {
      if(nbOfComp==nbOfComp2)
        {
          ret=DataArrayDouble::New();
          ret->alloc(nbOfTuple,nbOfComp);
          std::transform(a1->begin(),a1->end(),a2->begin(),ret->getPointer(),std::multiplies<double>());
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
            throw INTERP_KERNEL::Exception("Nb of components mismatch for array Multiply !");
        }
    }
  else if((nbOfTuple==1 && nbOfTuple2>1) || (nbOfTuple>1 && nbOfTuple2==1))
    {
      if(nbOfComp==nbOfComp2)
        {
          int nbOfTupleMax=std::max(nbOfTuple,nbOfTuple2);
          const DataArrayDouble *aMin=nbOfTuple>nbOfTuple2?a2:a1;
          const DataArrayDouble *aMax=nbOfTuple>nbOfTuple2?a1:a2;
          const double *aMinPtr=aMin->getConstPointer(),*aMaxPtr=aMax->getConstPointer();
          ret=DataArrayDouble::New();
          ret->alloc(nbOfTupleMax,nbOfComp);
          double *res=ret->getPointer();
          for(int i=0;i<nbOfTupleMax;i++)
            res=std::transform(aMaxPtr+i*nbOfComp,aMaxPtr+(i+1)*nbOfComp,aMinPtr,res,std::multiplies<double>());
          ret->copyStringInfoFrom(*aMax);
        }
      else
        throw INTERP_KERNEL::Exception("Nb of components mismatch for array Multiply !");
    }
  else
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array Multiply !");
  ret->incrRef();
  return ret;
}

void DataArrayDouble::multiplyEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::multiplyEqual : input DataArrayDouble instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayDouble::multiplyEqual !";
  int nbOfTuple=getNumberOfTuples();
  int nbOfTuple2=other->getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  int nbOfComp2=other->getNumberOfComponents();
  if(nbOfTuple==nbOfTuple2)
    {
      if(nbOfComp==nbOfComp2)
        {
          std::transform(begin(),end(),other->begin(),getPointer(),std::multiplies<double>());
        }
      else if(nbOfComp2==1)
        {
          double *ptr=getPointer();
          const double *ptrc=other->getConstPointer();
          for(int i=0;i<nbOfTuple;i++)
            std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptr+i*nbOfComp,std::bind2nd(std::multiplies<double>(),*ptrc++));
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else if(nbOfTuple2==1)
    {
      if(nbOfComp2==nbOfComp)
        {
          double *ptr=getPointer();
          const double *ptrc=other->getConstPointer();
          for(int i=0;i<nbOfTuple;i++)
            std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptrc,ptr+i*nbOfComp,std::multiplies<double>());
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else
    throw INTERP_KERNEL::Exception(msg);
  declareAsNew();
}

DataArrayDouble *DataArrayDouble::Divide(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::Divide : input DataArrayDouble instance is NULL !");
  int nbOfTuple1=a1->getNumberOfTuples();
  int nbOfTuple2=a2->getNumberOfTuples();
  int nbOfComp1=a1->getNumberOfComponents();
  int nbOfComp2=a2->getNumberOfComponents();
  if(nbOfTuple2==nbOfTuple1)
    {
      if(nbOfComp1==nbOfComp2)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
          ret->alloc(nbOfTuple2,nbOfComp1);
          std::transform(a1->begin(),a1->end(),a2->begin(),ret->getPointer(),std::divides<double>());
          ret->copyStringInfoFrom(*a1);
          ret->incrRef();
          return ret;
        }
      else if(nbOfComp2==1)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
          ret->alloc(nbOfTuple1,nbOfComp1);
          const double *a2Ptr=a2->getConstPointer();
          const double *a1Ptr=a1->getConstPointer();
          double *res=ret->getPointer();
          for(int i=0;i<nbOfTuple1;i++)
            res=std::transform(a1Ptr+i*nbOfComp1,a1Ptr+(i+1)*nbOfComp1,res,std::bind2nd(std::divides<double>(),a2Ptr[i]));
          ret->copyStringInfoFrom(*a1);
          ret->incrRef();
          return ret;
        }
      else
        {
          a1->checkNbOfComps(nbOfComp2,"Nb of components mismatch for array Divide !");
          return 0;
        }
    }
  else if(nbOfTuple2==1)
    {
      a1->checkNbOfComps(nbOfComp2,"Nb of components mismatch for array Divide !");
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
      ret->alloc(nbOfTuple1,nbOfComp1);
      const double *a1ptr=a1->getConstPointer(),*a2ptr=a2->getConstPointer();
      double *pt=ret->getPointer();
      for(int i=0;i<nbOfTuple1;i++)
        pt=std::transform(a1ptr+i*nbOfComp1,a1ptr+(i+1)*nbOfComp1,a2ptr,pt,std::divides<double>());
      ret->copyStringInfoFrom(*a1);
      ret->incrRef();
      return ret;
    }
  else
    {
      a1->checkNbOfTuples(nbOfTuple2,"Nb of tuples mismatch for array Divide !");//will always throw an exception
      return 0;
    }
}

void DataArrayDouble::divideEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::divideEqual : input DataArrayDouble instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayDouble::divideEqual !";
  int nbOfTuple=getNumberOfTuples();
  int nbOfTuple2=other->getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  int nbOfComp2=other->getNumberOfComponents();
  if(nbOfTuple==nbOfTuple2)
    {
      if(nbOfComp==nbOfComp2)
        {
          std::transform(begin(),end(),other->begin(),getPointer(),std::divides<double>());
        }
      else if(nbOfComp2==1)
        {
          double *ptr=getPointer();
          const double *ptrc=other->getConstPointer();
          for(int i=0;i<nbOfTuple;i++)
            std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptr+i*nbOfComp,std::bind2nd(std::divides<double>(),*ptrc++));
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else if(nbOfTuple2==1)
    {
      if(nbOfComp2==nbOfComp)
        {
          double *ptr=getPointer();
          const double *ptrc=other->getConstPointer();
          for(int i=0;i<nbOfTuple;i++)
            std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptrc,ptr+i*nbOfComp,std::divides<double>());
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else
    throw INTERP_KERNEL::Exception(msg);
  declareAsNew();
}

/*!
 * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
 * Server side.
 */
void DataArrayDouble::getTinySerializationIntInformation(std::vector<int>& tinyInfo) const
{
  tinyInfo.resize(2);
  if(isAllocated())
    {
      tinyInfo[0]=getNumberOfTuples();
      tinyInfo[1]=getNumberOfComponents();
    }
  else
    {
      tinyInfo[0]=-1;
      tinyInfo[1]=-1;
    }
}

/*!
 * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
 * Server side.
 */
void DataArrayDouble::getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const
{
  if(isAllocated())
    {
      int nbOfCompo=getNumberOfComponents();
      tinyInfo.resize(nbOfCompo+1);
      tinyInfo[0]=getName();
      for(int i=0;i<nbOfCompo;i++)
        tinyInfo[i+1]=getInfoOnComponent(i);
    }
  else
    {
      tinyInfo.resize(1);
      tinyInfo[0]=getName();
    }
}

/*!
 * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
 * This method returns if a feeding is needed.
 */
bool DataArrayDouble::resizeForUnserialization(const std::vector<int>& tinyInfoI)
{
  int nbOfTuple=tinyInfoI[0];
  int nbOfComp=tinyInfoI[1];
  if(nbOfTuple!=-1 || nbOfComp!=-1)
    {
      alloc(nbOfTuple,nbOfComp);
      return true;
    }
  return false;
}

/*!
 * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
 */
void DataArrayDouble::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<std::string>& tinyInfoS)
{
  setName(tinyInfoS[0].c_str());
  if(isAllocated())
    {
      int nbOfCompo=getNumberOfComponents();
      for(int i=0;i<nbOfCompo;i++)
        setInfoOnComponent(i,tinyInfoS[i+1].c_str());
    }
}

DataArrayDoubleIterator::DataArrayDoubleIterator(DataArrayDouble *da):_da(da),_tuple_id(0),_nb_comp(0),_nb_tuple(0)
{
  if(_da)
    {
      _da->incrRef();
      if(_da->isAllocated())
        {
          _nb_comp=da->getNumberOfComponents();
          _nb_tuple=da->getNumberOfTuples();
          _pt=da->getPointer();
        }
    }
}

DataArrayDoubleIterator::~DataArrayDoubleIterator()
{
  if(_da)
    _da->decrRef();
}

DataArrayDoubleTuple *DataArrayDoubleIterator::nextt()
{
  if(_tuple_id<_nb_tuple)
    {
      _tuple_id++;
      DataArrayDoubleTuple *ret=new DataArrayDoubleTuple(_pt,_nb_comp);
      _pt+=_nb_comp;
      return ret;
    }
  else
    return 0;
}

DataArrayDoubleTuple::DataArrayDoubleTuple(double *pt, int nbOfComp):_pt(pt),_nb_of_compo(nbOfComp)
{
}


std::string DataArrayDoubleTuple::repr() const
{
  std::ostringstream oss; oss.precision(15); oss << "(";
  for(int i=0;i<_nb_of_compo-1;i++)
    oss << _pt[i] << ", ";
  oss << _pt[_nb_of_compo-1] << ")";
  return oss.str();
}

double DataArrayDoubleTuple::doubleValue() const throw(INTERP_KERNEL::Exception)
{
  if(_nb_of_compo==1)
    return *_pt;
  throw INTERP_KERNEL::Exception("DataArrayDoubleTuple::doubleValue : DataArrayDoubleTuple instance has not exactly 1 component -> Not possible to convert it into a double precision float !");
}

/*!
 * This method returns a newly allocated instance the caller should dealed with by a ParaMEDMEM::DataArrayDouble::decrRef.
 * This method performs \b no copy of data. The content is only referenced using ParaMEDMEM::DataArrayDouble::useArray with ownership set to \b false.
 * This method throws an INTERP_KERNEL::Exception is it is impossible to match sizes of \b this that is too say \b nbOfCompo=this->_nb_of_elem and \bnbOfTuples==1 or
 * \b nbOfCompo=1 and \bnbOfTuples==this->_nb_of_elem.
 */
DataArrayDouble *DataArrayDoubleTuple::buildDADouble(int nbOfTuples, int nbOfCompo) const throw(INTERP_KERNEL::Exception)
{
  if((_nb_of_compo==nbOfCompo && nbOfTuples==1) || (_nb_of_compo==nbOfTuples && nbOfCompo==1))
    {
      DataArrayDouble *ret=DataArrayDouble::New();
      ret->useArray(_pt,false,CPP_DEALLOC,nbOfTuples,nbOfCompo);
      return ret;
    }
  else
    {
      std::ostringstream oss; oss << "DataArrayDoubleTuple::buildDADouble : unable to build a requested DataArrayDouble instance with nbofTuple=" << nbOfTuples << " and nbOfCompo=" << nbOfCompo;
      oss << ".\nBecause the number of elements in this is " << _nb_of_compo << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
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

/*!
 * This method differs from DataArray::setInfoOnComponents in the sense that if 'this->getNumberOfComponents()!=info.size()'
 * and if 'this' is not allocated it will change the number of components of 'this'.
 * If 'this->getNumberOfComponents()==info.size()' the behaviour is the same than DataArray::setInfoOnComponents method.
 * If 'this->getNumberOfComponents()!=info.size()' and the 'this' is already allocated an exception will be thrown.
 */
void DataArrayInt::setInfoAndChangeNbOfCompo(const std::vector<std::string>& info) throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=(int)info.size())
    {
      if(!isAllocated())
        _info_on_compo=info;
      else
        {
          std::ostringstream oss; oss << "DataArrayInt::setInfoAndChangeNbOfCompo : input is of size " << info.size() << " whereas number of components is equal to " << getNumberOfComponents() << "  and this is already allocated !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  else
    _info_on_compo=info;
}

/*!
 * This method returns the only one value in 'this', if and only if number of elements (nb of tuples * nb of components) is equal to 1, and that 'this' is allocated.
 * If one or more conditions is not fulfilled an exception will be thrown.
 */
int DataArrayInt::intValue() const throw(INTERP_KERNEL::Exception)
{
  if(isAllocated())
    {
      if(getNbOfElems()==1)
        {
          return *getConstPointer();
        }
      else
        throw INTERP_KERNEL::Exception("DataArrayInt::intValue : DataArrayInt instance is allocated but number of elements is not equal to 1 !");
    }
  else
    throw INTERP_KERNEL::Exception("DataArrayInt::intValue : DataArrayInt instance is not allocated !");
}

/*!
 * This method expects that \b this is well allocated. If not an INTERP_KERNEL::Exception will be thrown. This method is useful for a quick comparison of many instances of DataArrayInt.
 */
int DataArrayInt::getHashCode() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfElems=getNbOfElems();
  int ret=nbOfElems*65536;
  int delta=3;
  if(nbOfElems>48)
    delta=nbOfElems/8;
  int ret0=0;
  const int *pt=begin();
  for(int i=0;i<nbOfElems;i+=delta)
    ret0+=pt[i] & 0x1FFF;
  return ret+ret0;
}

/*!
 * This method should be called on an allocated DataArrayInt instance. If not an exception will be throw !
 * This method checks the number of tupes. If it is equal to 0, it returns true, if not false is returned.
 */
bool DataArrayInt::empty() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  return getNumberOfTuples()==0;
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

void DataArrayInt::cpyFrom(const DataArrayInt& other) throw(INTERP_KERNEL::Exception)
{
  other.checkAllocated();
  int nbOfTuples=other.getNumberOfTuples();
  int nbOfComp=other.getNumberOfComponents();
  allocIfNecessary(nbOfTuples,nbOfComp);
  int nbOfElems=nbOfTuples*nbOfComp;
  int *pt=getPointer();
  const int *ptI=other.getConstPointer();
  for(int i=0;i<nbOfElems;i++)
    pt[i]=ptI[i];
  copyStringInfoFrom(other);
}

void DataArrayInt::allocIfNecessary(int nbOfTuple, int nbOfCompo)
{
  if(isAllocated())
    {
      if(nbOfTuple!=getNumberOfTuples() || nbOfCompo!=getNumberOfComponents())
        alloc(nbOfTuple,nbOfCompo);
    }
  else
    alloc(nbOfTuple,nbOfCompo);
}

void DataArrayInt::alloc(int nbOfTuple, int nbOfCompo) throw(INTERP_KERNEL::Exception)
{
  if(nbOfTuple<0 || nbOfCompo<0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::alloc : request for negative length of data !");
  _nb_of_tuples=nbOfTuple;
  _info_on_compo.resize(nbOfCompo);
  _mem.alloc(nbOfCompo*_nb_of_tuples);
  declareAsNew();
}

void DataArrayInt::fillWithZero() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  _mem.fillWithValue(0);
  declareAsNew();
}

void DataArrayInt::fillWithValue(int val) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  _mem.fillWithValue(val);
  declareAsNew();
}

void DataArrayInt::iota(int init) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::iota : works only for arrays with only one component, you can call 'rearrange' method before !");
  int *ptr=getPointer();
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

void DataArrayInt::writeVTK(std::ostream& ofs, int indent, const char *type, const char *nameInFile) const throw(INTERP_KERNEL::Exception)
{
  std::string idt(indent,' ');
  ofs << idt << "<DataArray type=\"" << type << "\" Name=\"" << nameInFile << "\" NumberOfComponents=\"" << getNumberOfComponents() << "\"";
  ofs << " format=\"ascii\" RangeMin=\"" << getMinValueInArray() << "\" RangeMax=\"" << getMaxValueInArray() << "\">\n" << idt;
  std::copy(begin(),end(),std::ostream_iterator<int>(ofs," "));
  ofs << std::endl << idt << "</DataArray>\n";
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

/*!
 * This method expects a number of components equal to 1.
 * This method sweeps all the values (tuples) in 'this' (it should be allocated) and for each value v is replaced by
 * indArr[v] where 'indArr' is defined by ['indArrBg','indArrEnd').
 * This method is safe that is to say if there is a value in 'this' not in [0,std::distance('indArrBg','indArrEnd')) an exception
 * will be thrown.
 */
void DataArrayInt::transformWithIndArr(const int *indArrBg, const int *indArrEnd) throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("Call transformWithIndArr method on DataArrayInt with only one component, you can call 'rearrange' method before !");
  int nbElemsIn=(int)std::distance(indArrBg,indArrEnd);
  int nbOfTuples=getNumberOfTuples();
  int *pt=getPointer();
  for(int i=0;i<nbOfTuples;i++,pt++)
    {
      if(*pt>=0 && *pt<nbElemsIn)
        *pt=indArrBg[*pt];
      else
        {
          std::ostringstream oss; oss << "DataArrayInt::transformWithIndArr : error on tuple #" << i << " value is " << *pt << " and indirectionnal array as a size equal to " << nbElemsIn;
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  declareAsNew();
}

/*!
 * 'this' should be allocated and with numberOfComponents set to one. If not an exception will be thrown.
 * This method takes as input an array defined by ['arrBg','arrEnd'). The size of the array (std::distance(arrBg,arrEnd)) is equal to the number of cast + 1.
 * The values contained in ['arrBg','arrEnd') should be sorted ascendently. No check of this will be done. If not the result is not waranted.
 * For each cast j the value range that defines the cast is equal to [arrBg[j],arrBg[j+1]).
 * This method returns three arrays (to be managed by the caller).
 * This method is typically usefull for entity number spliting by types for example.
 * Example : If 'this' contains [6,5,0,3,2,7,8,1,4] and if ['arrBg','arrEnd') contains [0,4,9] then the output of this method will be :
 * - 'castArr'        : [1,1,0,0,0,1,1,0,1]
 * - 'rankInsideCast' : [2,1,0,3,2,3,4,1,0]
 * - 'return' : [0,1]
 *
 * @param castArr is a returned param has the same number of tuples than 'this' and number of components set to one. In case of sucess, this param contains for each tuple in 'this' in which cast it holds.
 * @param rankInsideCast is an another returned param has the same number of tuples than 'this' and number of components set to one too. In case of sucess, this param contains for each tuple its rank inside its cast.
 * @param castsPresent the casts that 'this' contains.
 * @throw if a value in 'this' is greater or equal to the last value of ['arrBg','arrEnd')
 */
void DataArrayInt::splitByValueRange(const int *arrBg, const int *arrEnd,
                                     DataArrayInt *& castArr, DataArrayInt *& rankInsideCast, DataArrayInt *& castsPresent) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("Call splitByValueRange  method on DataArrayInt with only one component, you can call 'rearrange' method before !");
  int nbOfTuples=getNumberOfTuples();
  std::size_t nbOfCast=std::distance(arrBg,arrEnd);
  if(nbOfCast<2)
    throw INTERP_KERNEL::Exception("DataArrayInt::splitByValueRange : The input array giving the cast range values should be of size >=2 !");
  nbOfCast--;
  const int *work=getConstPointer();
  typedef std::reverse_iterator<const int *> rintstart;
  rintstart bg(arrEnd);//OK no problem because size of 'arr' is greater of equal 2
  rintstart end2(arrBg);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret1=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret2=DataArrayInt::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret3=DataArrayInt::New();
  ret1->alloc(nbOfTuples,1);
  ret2->alloc(nbOfTuples,1);
  int *ret1Ptr=ret1->getPointer();
  int *ret2Ptr=ret2->getPointer();
  std::set<std::size_t> castsDetected;
  for(int i=0;i<nbOfTuples;i++)
    {
      rintstart res=std::find_if(bg,end2,std::bind2nd(std::less_equal<int>(), work[i]));
      std::size_t pos=std::distance(bg,res);
      std::size_t pos2=nbOfCast-pos;
      if(pos2<nbOfCast)
        {
          ret1Ptr[i]=(int)pos2;
          ret2Ptr[i]=work[i]-arrBg[pos2];
          castsDetected.insert(pos2);
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayInt::splitByValueRange : At rank #" << i << " the value is " << work[i] << " whereas the last value is " << *bg;
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  ret3->alloc((int)castsDetected.size(),1);
  std::copy(castsDetected.begin(),castsDetected.end(),ret3->getPointer());
  ret1->incrRef();
  castArr=ret1;
  ret2->incrRef();
  rankInsideCast=ret2;
  ret3->incrRef();
  castsPresent=ret3;
}

/*!
 * This method expects a number of components equal to 1.
 * This method sweeps all the values (tuples) in 'this' (it should be allocated) and for each value v on place i, place indArr[v] will have 
 * value i.
 * indArr[v] where 'indArr' is defined by ['indArrBg','indArrEnd').
 * This method is half/safe that is to say if there is location i so that indArr[v] is not in [0,this->getNumberOfTuples()) an exception
 * will be thrown.
 */
DataArrayInt *DataArrayInt::transformWithIndArrR(const int *indArrBg, const int *indArrEnd) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("Call transformWithIndArrR method on DataArrayInt with only one component, you can call 'rearrange' method before !");
  int nbElemsIn=(int)std::distance(indArrBg,indArrEnd);
  int nbOfTuples=getNumberOfTuples();
  const int *pt=getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbOfTuples,1);
  ret->fillWithValue(-1);
  int *tmp=ret->getPointer();
  for(int i=0;i<nbOfTuples;i++,pt++)
    {
      int pos=indArrBg[*pt];
      if(pos>=0 && pos<nbElemsIn)
        tmp[pos]=i;
      else
        {
          std::ostringstream oss; oss << "DataArrayInt::transformWithIndArrR : error on tuple #" << i << " value is " << *pt << " and indirectionnal array as a size equal to " << nbElemsIn;
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  ret->incrRef();
  return ret;
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

bool DataArrayInt::isEqualIfNotWhy(const DataArrayInt& other, std::string& reason) const
{
  if(!areInfoEqualsIfNotWhy(other,reason))
    return false;
  return _mem.isEqual(other._mem,0,reason);
}

bool DataArrayInt::isEqual(const DataArrayInt& other) const
{
  std::string tmp;
  return isEqualIfNotWhy(other,tmp);
}

bool DataArrayInt::isEqualWithoutConsideringStr(const DataArrayInt& other) const
{
  std::string tmp;
  return _mem.isEqual(other._mem,0,tmp);
}

bool DataArrayInt::isEqualWithoutConsideringStrAndOrder(const DataArrayInt& other) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> a=deepCpy();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> b=other.deepCpy();
  a->sort();
  b->sort();
  return a->isEqualWithoutConsideringStr(*b);
}

void DataArrayInt::sort(bool asc) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::sort : only supported with 'this' array with ONE component !");
  _mem.sort(asc);
}

void DataArrayInt::reverse() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::reverse : only supported with 'this' array with ONE component !");
  _mem.reverse();
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

/*!
 * This method expects that 'this' is allocated, if not an exception is thrown.
 * This method in case of success returns a newly created array the user should deal with.
 * In the case of having a renumber array in "old to new" format. More info on renumbering \ref MEDCouplingArrayRenumbering "here".
 */
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
  ret->alloc((int)std::distance(new2OldBg,new2OldEnd),nbComp);
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
  ret->alloc((int)std::distance(new2OldBg,new2OldEnd),nbComp);
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
 * Idem than DataArrayInt::selectByTupleIdSafe except that the input array is not constructed explicitely.
 * The convention is as python one. ['bg','end2') with steps of 'step'.
 * Returns a newly created array.
 * This method is an extension of DataArrayInt::substr method.
 * 
 * \sa DataArrayInt::substr
 */
DataArrayInt *DataArrayInt::selectByTupleId2(int bg, int end2, int step) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  int nbComp=getNumberOfComponents();
  int newNbOfTuples=GetNumberOfItemGivenBES(bg,end2,step,"DataArrayInt::selectByTupleId2 : ");
  ret->alloc(newNbOfTuples,nbComp);
  int *pt=ret->getPointer();
  const int *srcPt=getConstPointer()+bg*nbComp;
  for(int i=0;i<newNbOfTuples;i++,srcPt+=step*nbComp)
    std::copy(srcPt,srcPt+nbComp,pt+i*nbComp);
  ret->copyStringInfoFrom(*this);
  ret->incrRef();
  return ret;
}

/*!
 * This method returns a newly allocated array that is the concatenation of all tuples ranges in param 'ranges'.
 * Each pair in input 'ranges' is in [begin,end) format. If there is a range in 'ranges' so that end is before begin an exception
 * will be thrown. If there is a range in 'ranges' so that end is greater than number of tuples of 'this', an exception will be thrown too.
 */
DataArrayInt *DataArrayInt::selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfComp=getNumberOfComponents();
  int nbOfTuplesThis=getNumberOfTuples();
  if(ranges.empty())
    {
      DataArrayInt *ret=DataArrayInt::New();
      ret->alloc(0,nbOfComp);
      ret->copyStringInfoFrom(*this);
      return ret;
    }
  int ref=ranges.front().first;
  int nbOfTuples=0;
  bool isIncreasing=true;
  for(std::vector<std::pair<int,int> >::const_iterator it=ranges.begin();it!=ranges.end();it++)
    {
      if((*it).first<=(*it).second)
        {
          if((*it).first>=0 && (*it).second<=nbOfTuplesThis)
            {
              nbOfTuples+=(*it).second-(*it).first;
              if(isIncreasing)
                isIncreasing=ref<=(*it).first;
              ref=(*it).second;
            }
          else
            {
              std::ostringstream oss; oss << "DataArrayInt::selectByTupleRanges : on range #" << std::distance(ranges.begin(),it);
              oss << " (" << (*it).first << "," << (*it).second << ") is greater than number of tuples of this :" << nbOfTuples << " !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayInt::selectByTupleRanges : on range #" << std::distance(ranges.begin(),it);
          oss << " (" << (*it).first << "," << (*it).second << ") end is before begin !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  if(isIncreasing && nbOfTuplesThis==nbOfTuples)
    return deepCpy();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbOfTuples,nbOfComp);
  ret->copyStringInfoFrom(*this);
  const int *src=getConstPointer();
  int *work=ret->getPointer();
  for(std::vector<std::pair<int,int> >::const_iterator it=ranges.begin();it!=ranges.end();it++)
    work=std::copy(src+(*it).first*nbOfComp,src+(*it).second*nbOfComp,work);
  ret->incrRef();
  return ret;
}

/*!
 * This method works only for arrays having single component.
 * If this contains the array a1 containing [9,10,0,6,4,11,3,7] this method returns an array a2 [5,6,0,3,2,7,1,4].
 * By doing a1.renumber(a2) the user will obtain array a3 equal to a1 sorted.
 * This method is useful for renumbering (in MED file for example). This method is used by MEDCouplingFieldDouble::renumberCells when check is set to true.
 * This method throws an exception if more 2 or more elements in 'this' are same.
 */
DataArrayInt *DataArrayInt::checkAndPreparePermutation() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::checkAndPreparePermutation : number of components must == 1 !");
  int nbTuples=getNumberOfTuples();
  const int *pt=getConstPointer();
  int *pt2=CheckAndPreparePermutation(pt,pt+nbTuples);
  DataArrayInt *ret=DataArrayInt::New();
  ret->useArray(pt2,true,CPP_DEALLOC,nbTuples,1);
  return ret;
}

/*!
 * This method makes the assumption that 'this' is correctly set, and has exactly one component. If not an exception will be thrown.
 * Given a sujective application defined by 'this' from a set of size this->getNumberOfTuples() to a set of size targetNb.
 * 'targetNb'<this->getNumberOfTuples(). 'this' should be surjective that is to say for each id in [0,'targetNb') it exists at least one tupleId tid
 * so that this->getIJ(tid,0)==id.
 * If not an exception will be thrown.
 * This method returns 2 newly allocated arrays 'arr' and 'arrI', corresponding respectively to array and its corresponding index.
 * This method is usefull for methods that returns old2New numbering concecutive to a reduction ( MEDCouplingUMesh::zipConnectivityTraducer, MEDCouplingUMesh::zipConnectivityTraducer for example)
 * Example : if 'this' equals [0,3,2,3,2,2,1,2] this method will return arrI=[0,1,2,6,8] arr=[0,  6,  2,4,5,7,  1,3]
 * That is to say elt id 2 has arrI[2+1]-arrI[2]=4 places in 'this'. The corresponding tuple ids are [2,4,5,7].
 */
void DataArrayInt::changeSurjectiveFormat(int targetNb, DataArrayInt *&arr, DataArrayInt *&arrI) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::changeSurjectiveFormat : number of components must == 1 !");
  int nbOfTuples=getNumberOfTuples();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret(DataArrayInt::New());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> retI(DataArrayInt::New());
  retI->alloc(targetNb+1,1);
  const int *input=getConstPointer();
  std::vector< std::vector<int> > tmp(targetNb);
  for(int i=0;i<nbOfTuples;i++)
    {
      int tmp2=input[i];
      if(tmp2<targetNb)
        tmp[tmp2].push_back(i);
      else
        {
          std::ostringstream oss; oss << "DataArrayInt::changeSurjectiveFormat : At pos " << i << " presence of element " << tmp2 << " higher than " << targetNb;
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  int *retIPtr=retI->getPointer();
  *retIPtr=0;
  for(std::vector< std::vector<int> >::const_iterator it1=tmp.begin();it1!=tmp.end();it1++,retIPtr++)
    retIPtr[1]=retIPtr[0]+(int)((*it1).size());
  if(nbOfTuples!=retI->getIJ(targetNb,0))
    throw INTERP_KERNEL::Exception("DataArrayInt::changeSurjectiveFormat : big problem should never happen !");
  ret->alloc(nbOfTuples,1);
  int *retPtr=ret->getPointer();
  for(std::vector< std::vector<int> >::const_iterator it1=tmp.begin();it1!=tmp.end();it1++)
    retPtr=std::copy((*it1).begin(),(*it1).end(),retPtr);
  ret->incrRef();
  retI->incrRef();
  arr=ret;
  arrI=retI;
}

/*!
 * This static method computes a old 2 new format DataArrayInt instance from a zip representation of a surjective format (retrived by DataArrayDouble::findCommonTuples for example)
 * The retrieved array minimizes the permutation.
 * Let's take an example : 
 * If 'nbOfOldTuples'==10 and 'arr'==[0,3, 5,7,9] and 'arrI'==[0,2,5] it returns the following array [0,1,2,0,3,4,5,4,6,4] and newNbOfTuples==7.
 *
 * @param nbOfOldTuples is the number of tuples in initial array.
 * @param arr is the list of tuples ids grouped by 'arrI' array
 * @param arrI is the entry point of 'arr' array. arrI->getNumberOfTuples()-1 is the number of common groups > 1 tuple.
 * @param newNbOfTuples output parameter that retrieves the new number of tuples after surjection application
 */
DataArrayInt *DataArrayInt::BuildOld2NewArrayFromSurjectiveFormat2(int nbOfOldTuples, const DataArrayInt *arr, const DataArrayInt *arrI, int &newNbOfTuples) throw(INTERP_KERNEL::Exception)
{
  if(!arr || !arrI)
    throw INTERP_KERNEL::Exception("DataArrayInt::BuildOld2NewArrayFromSurjectiveFormat2 : presence of NULL ref of DataArrayInt in input !");
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbOfOldTuples,1);
  int *pt=ret->getPointer();
  std::fill(pt,pt+nbOfOldTuples,-1);
  int nbOfGrps=arrI->getNumberOfTuples()-1;
  const int *cIPtr=arrI->getConstPointer();
  const int *cPtr=arr->getConstPointer();
  for(int i=0;i<nbOfGrps;i++)
    pt[cPtr[cIPtr[i]]]=-(i+2);
  int newNb=0;
  for(int iNode=0;iNode<nbOfOldTuples;iNode++)
    {
      if(pt[iNode]<0)
        {
          if(pt[iNode]==-1)
            pt[iNode]=newNb++;
          else
            {
              int grpId=-(pt[iNode]+2);
              for(int j=cIPtr[grpId];j<cIPtr[grpId+1];j++)
                pt[cPtr[j]]=newNb;
              newNb++;
            }
        }
    }
  newNbOfTuples=newNb;
  return ret;
}

/*!
 * This method expects that 'this' is allocated and with only one component. If not an exception will be thrown.
 * This method returns a newly created array with 'this->getNumberOfTuples()' tuples and 1 component.
 * This methods returns an 'old2New' corresponding array that allows to follow the following rules :
 * - Lower a value in tuple in 'this' is, higher is its priority.
 * - If two tuples i and j have same value if i<j then ret[i]<ret[j]
 * - The first tuple with the lowest value will have the value 0, inversely the last tuple with highest value will have value 'this->getNumberOfTuples()-1'
 * 
 * Example if 'this' contains the following array : [2,0,1,1,0,1,2,0,1,1,0,0] this method returns [10,0,5,6,1,7,11,2,8,9,3,4]
 */
DataArrayInt *DataArrayInt::buildPermArrPerLevel() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::buildPermArrPerLevel : number of components must == 1 !");
  int nbOfTuples=getNumberOfTuples();
  const int *pt=getConstPointer();
  std::map<int,int> m;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbOfTuples,1);
  int *opt=ret->getPointer();
  for(int i=0;i<nbOfTuples;i++,pt++,opt++)
    {
      int val=*pt;
      std::map<int,int>::iterator it=m.find(val);
      if(it!=m.end())
        {
          *opt=(*it).second;
          (*it).second++;
        }
      else
        {
          *opt=0;
          m.insert(std::pair<int,int>(val,1));
        }
    }
  int sum=0;
  for(std::map<int,int>::iterator it=m.begin();it!=m.end();it++)
    {
      int vt=(*it).second;
      (*it).second=sum;
      sum+=vt;
    }
  pt=getConstPointer();
  opt=ret->getPointer();
  for(int i=0;i<nbOfTuples;i++,pt++,opt++)
    *opt+=m[*pt];
  //
  ret->incrRef();
  return ret;
}

/*!
 * This method checks that 'this' is with numberofcomponents == 1 and that it is equal to
 * stdext::iota() of size getNumberOfTuples. This method is particalary usefull for DataArrayInt instances
 * that represents a renumbering array to check the real need in renumbering. 
 */
bool DataArrayInt::isIdentity() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    return false;
  int nbOfTuples=getNumberOfTuples();
  const int *pt=getConstPointer();
  for(int i=0;i<nbOfTuples;i++,pt++)
    if(*pt!=i)
      return false;
  return true;
}

bool DataArrayInt::isUniform(int val) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::isUniform : must be applied on DataArrayInt with only one component, you can call 'rearrange' method before !");
  int nbOfTuples=getNumberOfTuples();
  const int *w=getConstPointer();
  const int *end2=w+nbOfTuples;
  for(;w!=end2;w++)
    if(*w!=val)
      return false;
  return true;
}

DataArrayDouble *DataArrayInt::convertToDblArr() const
{
  checkAllocated();
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
 * This method is a specialization of method DataArrayInt::selectByTupleId2.
 *
 * \sa DataArrayInt::selectByTupleId2
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
 * This method makes the assumption that \b this is allocated. If not an INTERP_KERNEL::Exception will be raised.
 * This method does not echange the values stored in \b this. Simply, the number of components before the call becomes the number of
 * tuples and inversely the number of tuples becomes the number of components. \b WARNING the info on components can be alterated by this method.
 */
void DataArrayInt::transpose() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbOfTuples=getNumberOfTuples();
  rearrange(nbOfTuples);
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
  _mem.reAlloc((int)_info_on_compo.size()*nbOfTuples);
  _nb_of_tuples=nbOfTuples;
  declareAsNew();
}

DataArrayInt *DataArrayInt::keepSelectedComponents(const std::vector<int>& compoIds) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret(DataArrayInt::New());
  int newNbOfCompo=(int)compoIds.size();
  int oldNbOfCompo=getNumberOfComponents();
  for(std::vector<int>::const_iterator it=compoIds.begin();it!=compoIds.end();it++)
    DataArray::CheckValueInRange(oldNbOfCompo,(*it),"keepSelectedComponents invalid requested component");
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
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayInt::meldWith : DataArrayInt pointer in input is NULL !");
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
  std::size_t partOfCompoSz=compoIds.size();
  int nbOfCompo=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  const int *ac=a->getConstPointer();
  int *nc=getPointer();
  for(int i=0;i<nbOfTuples;i++)
    for(std::size_t j=0;j<partOfCompoSz;j++,ac++)
      nc[nbOfCompo*i+compoIds[j]]=*ac;
}

/*!
 * This method performs a partial assignment of 'this' using 'a' as input. Other input parameters specifies the subpart being considered by the assignment.
 * 'strictCompoCompare' specifies if DataArray 'a' should have exactly same number of components and tuples than 'this' (true) or not (false). By default set to true with maximal test.
 */
void DataArrayInt::setPartOfValues1(const DataArrayInt *a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare) throw(INTERP_KERNEL::Exception)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayInt::setPartOfValues1 : DataArrayInt pointer in input is NULL !");
  const char msg[]="DataArrayInt::setPartOfValues1";
  checkAllocated();
  a->checkAllocated();
  int newNbOfTuples=DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg);
  int newNbOfComp=DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg);
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
  DataArray::CheckValueInRangeEx(nbComp,bgComp,endComp,"invalid component value");
  a->checkNbOfElems(newNbOfTuples*newNbOfComp,msg);
  if(strictCompoCompare)
    a->checkNbOfTuplesAndComp(newNbOfTuples,newNbOfComp,msg);
  int *pt=getPointer()+bgTuples*nbComp+bgComp;
  const int *srcPt=a->getConstPointer();
  for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
    for(int j=0;j<newNbOfComp;j++,srcPt++)
      pt[j*stepComp]=*srcPt;
}

/*!
 * This method performs a partial assignment of 'this' using 'a' as input. Other input parameters specifies the subpart being considered by the assignment.
 */
void DataArrayInt::setPartOfValuesSimple1(int a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="DataArrayInt::setPartOfValuesSimple1";
  checkAllocated();
  int newNbOfTuples=DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg);
  int newNbOfComp=DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg);
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
  DataArray::CheckValueInRangeEx(nbComp,bgComp,endComp,"invalid component value");
  int *pt=getPointer()+bgTuples*nbComp+bgComp;
  for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
    for(int j=0;j<newNbOfComp;j++)
      pt[j*stepComp]=a;
}

/*!
 * This method performs a partial assignment of 'this' using 'a' as input. Other input parameters specifies the subpart being considered by the assignment.
 * 'strictCompoCompare' specifies if DataArray 'a' should have exactly same number of components and tuples than 'this' (true) or not (false). By default set to true with maximal test.
 */
void DataArrayInt::setPartOfValues2(const DataArrayInt *a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp, bool strictCompoCompare) throw(INTERP_KERNEL::Exception)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayInt::setPartOfValues2 : DataArrayInt pointer in input is NULL !");
  const char msg[]="DataArrayInt::setPartOfValues2";
  checkAllocated();
  a->checkAllocated();
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  for(const int *z=bgComp;z!=endComp;z++)
    DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
  int newNbOfTuples=(int)std::distance(bgTuples,endTuples);
  int newNbOfComp=(int)std::distance(bgComp,endComp);
  bool assignTech=true;
  if(a->getNbOfElems()==newNbOfTuples*newNbOfComp)
    {
      if(strictCompoCompare)
        a->checkNbOfTuplesAndComp(newNbOfTuples,newNbOfComp,msg);
    }
  else
    {
      a->checkNbOfTuplesAndComp(1,newNbOfComp,msg);
      assignTech=false;
    }
  int *pt=getPointer();
  const int *srcPt=a->getConstPointer();
  if(assignTech)
    {    
      for(const int *w=bgTuples;w!=endTuples;w++)
        {
          DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
          for(const int *z=bgComp;z!=endComp;z++,srcPt++)
            {    
              pt[(*w)*nbComp+(*z)]=*srcPt;
            }
        }
    }
  else
    {
      for(const int *w=bgTuples;w!=endTuples;w++)
        {
          const int *srcPt2=srcPt;
          DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
          for(const int *z=bgComp;z!=endComp;z++,srcPt2++)
            {    
              pt[(*w)*nbComp+(*z)]=*srcPt2;
            }
        }
    }
}

/*!
 * This method performs a partial assignment of 'this' using 'a' as input. Other input parameters specifies the subpart being considered by the assignment.
 */
void DataArrayInt::setPartOfValuesSimple2(int a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  for(const int *z=bgComp;z!=endComp;z++)
    DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
  int *pt=getPointer();
  for(const int *w=bgTuples;w!=endTuples;w++)
    for(const int *z=bgComp;z!=endComp;z++)
      {
        DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
        pt[(*w)*nbComp+(*z)]=a;
      }
}

/*!
 * This method performs a partial assignment of 'this' using 'a' as input. Other input parameters specifies the subpart being considered by the assignment.
 * 'strictCompoCompare' specifies if DataArray 'a' should have exactly same number of components and tuples than 'this' (true) or not (false). By default set to true with maximal test.
 */
void DataArrayInt::setPartOfValues3(const DataArrayInt *a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare) throw(INTERP_KERNEL::Exception)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayInt::setPartOfValues3 : DataArrayInt pointer in input is NULL !");
  const char msg[]="DataArrayInt::setPartOfValues3";
  checkAllocated();
  a->checkAllocated();
  int newNbOfComp=DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg);
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbComp,bgComp,endComp,"invalid component value");
  int newNbOfTuples=(int)std::distance(bgTuples,endTuples);
  bool assignTech=true;
  if(a->getNbOfElems()==newNbOfTuples*newNbOfComp)
    {
      if(strictCompoCompare)
        a->checkNbOfTuplesAndComp(newNbOfTuples,newNbOfComp,msg);
    }
  else
    {
      a->checkNbOfTuplesAndComp(1,newNbOfComp,msg);
      assignTech=false;
    }
  int *pt=getPointer()+bgComp;
  const int *srcPt=a->getConstPointer();
  if(assignTech)
    {
      for(const int *w=bgTuples;w!=endTuples;w++)
        for(int j=0;j<newNbOfComp;j++,srcPt++)
          {
            DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
            pt[(*w)*nbComp+j*stepComp]=*srcPt;
          }
    }
  else
    {
      for(const int *w=bgTuples;w!=endTuples;w++)
        {
          const int *srcPt2=srcPt;
          for(int j=0;j<newNbOfComp;j++,srcPt2++)
            {
              DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
              pt[(*w)*nbComp+j*stepComp]=*srcPt2;
            }
        }
    }
}

/*!
 * This method performs a partial assignment of 'this' using 'a' as input. Other input parameters specifies the subpart being considered by the assignment.
 */
void DataArrayInt::setPartOfValuesSimple3(int a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="DataArrayInt::setPartOfValuesSimple3";
  checkAllocated();
  int newNbOfComp=DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg);
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbComp,bgComp,endComp,"invalid component value");
  int *pt=getPointer()+bgComp;
  for(const int *w=bgTuples;w!=endTuples;w++)
    for(int j=0;j<newNbOfComp;j++)
      {
        DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
        pt[(*w)*nbComp+j*stepComp]=a;
      }
}

/*!
 * 'this', 'a' and 'tuplesSelec' are expected to be defined. If not an exception will be thrown.
 * @param a is an array having exactly the same number of components than 'this'
 * @param tuplesSelec is an array having exactly 2 components. The first one refers to the tuple ids of 'this' that will be set. The second one refers to the tuple ids of 'a' that will be used for setting.
 */
void DataArrayInt::setPartOfValuesAdv(const DataArrayInt *a, const DataArrayInt *tuplesSelec) throw(INTERP_KERNEL::Exception)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayInt::setPartOfValuesAdv : DataArrayInt pointer in input is NULL !");
  checkAllocated();
  a->checkAllocated();
  tuplesSelec->checkAllocated();
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=a->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("DataArrayInt::setPartOfValuesAdv : This and a do not have the same number of components !");
  if(tuplesSelec->getNumberOfComponents()!=2)
    throw INTERP_KERNEL::Exception("DataArrayInt::setPartOfValuesAdv : Expecting to have a tuple selector DataArrayInt instance with exactly 2 components !");
  int thisNt=getNumberOfTuples();
  int aNt=a->getNumberOfTuples();
  int *valsToSet=getPointer();
  const int *valsSrc=a->getConstPointer();
  for(const int *tuple=tuplesSelec->begin();tuple!=tuplesSelec->end();tuple+=2)
    {
      if(tuple[1]>=0 && tuple[1]<aNt)
        {
          if(tuple[0]>=0 && tuple[0]<thisNt)
            std::copy(valsSrc+nbOfComp*tuple[1],valsSrc+nbOfComp*(tuple[1]+1),valsToSet+nbOfComp*tuple[0]);
          else
            {
              std::ostringstream oss; oss << "DataArrayInt::setPartOfValuesAdv : Tuple #" << std::distance(tuplesSelec->begin(),tuple)/2;
              oss << " of 'tuplesSelec' request of tuple id #" << tuple[0] << " in 'this' ! It should be in [0," << thisNt << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayInt::setPartOfValuesAdv : Tuple #" << std::distance(tuplesSelec->begin(),tuple)/2;
          oss << " of 'tuplesSelec' request of tuple id #" << tuple[1] << " in 'a' ! It should be in [0," << aNt << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}

/*!
 * 'this', 'a' and 'tuplesSelec' are expected to be defined. If not an exception will be thrown.
 * This is a method that is a specialization to DataArrayInt::setPartOfValuesAdv method, except that here the tuple selection of 'a' is given by a range ('bg','end2' and 'step')
 * rather than an explicite array of tuple ids (given by the 2nd component) and the feeding is done in 'this' contiguously starting from 'tupleIdStart'.
 * @param a is an array having exactly the same number of components than 'this'
 */
void DataArrayInt::setContigPartOfSelectedValues(int tupleIdStart, const DataArrayInt*a, const DataArrayInt *tuplesSelec) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  a->checkAllocated();
  tuplesSelec->checkAllocated();
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=a->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("DataArrayInt::setContigPartOfSelectedValues : This and a do not have the same number of components !");
  if(tuplesSelec->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::setContigPartOfSelectedValues : Expecting to have a tuple selector DataArrayInt instance with exactly 1 component !");
  int thisNt=getNumberOfTuples();
  int aNt=a->getNumberOfTuples();
  int nbOfTupleToWrite=tuplesSelec->getNumberOfTuples();
  int *valsToSet=getPointer()+tupleIdStart*nbOfComp;
  if(tupleIdStart+nbOfTupleToWrite>thisNt)
    throw INTERP_KERNEL::Exception("DataArrayInt::setContigPartOfSelectedValues : invalid number range of values to write !");
  const int *valsSrc=a->getConstPointer();
  for(const int *tuple=tuplesSelec->begin();tuple!=tuplesSelec->end();tuple++,valsToSet+=nbOfComp)
    {
      if(*tuple>=0 && *tuple<aNt)
        {
          std::copy(valsSrc+nbOfComp*(*tuple),valsSrc+nbOfComp*(*tuple+1),valsToSet);
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayInt::setContigPartOfSelectedValues : Tuple #" << std::distance(tuplesSelec->begin(),tuple);
          oss << " of 'tuplesSelec' request of tuple id #" << *tuple << " in 'a' ! It should be in [0," << aNt << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}

/*!
 * 'this' and 'a' are expected to be defined. If not an exception will be thrown.
 * This is a method that is a specialization to DataArrayInt::setContigPartOfSelectedValues method, except that here the tuple selection is givenin a is done by a range ('bg','end2' and 'step')
 * rather than an explicite array of tuple ids.
 * @param a is an array having exactly the same number of components than 'this'
 */
void DataArrayInt::setContigPartOfSelectedValues2(int tupleIdStart, const DataArrayInt *a, int bg, int end2, int step) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  a->checkAllocated();
  int nbOfComp=getNumberOfComponents();
  const char msg[]="DataArrayInt::setContigPartOfSelectedValues2";
  int nbOfTupleToWrite=DataArray::GetNumberOfItemGivenBES(bg,end2,step,msg);
  if(nbOfComp!=a->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("DataArrayInt::setContigPartOfSelectedValues2 : This and a do not have the same number of components !");
  int thisNt=getNumberOfTuples();
  int aNt=a->getNumberOfTuples();
  int *valsToSet=getPointer()+tupleIdStart*nbOfComp;
  if(tupleIdStart+nbOfTupleToWrite>thisNt)
    throw INTERP_KERNEL::Exception("DataArrayInt::setContigPartOfSelectedValues2 : invalid number range of values to write !");
  if(end2>aNt)
    throw INTERP_KERNEL::Exception("DataArrayInt::setContigPartOfSelectedValues2 : invalid range of values to read !");
  const int *valsSrc=a->getConstPointer()+bg*nbOfComp;
  for(int i=0;i<nbOfTupleToWrite;i++,valsToSet+=nbOfComp,valsSrc+=step*nbOfComp)
    {
      std::copy(valsSrc,valsSrc+nbOfComp,valsToSet);
    }
}

/*!
 * This method is equivalent to DataArrayInt::getIJ except that here \b tupleId is checked to be in [0,this->getNumberOfTuples()) and compoId to be in [0,this->getNumberOfComponents()).
 * If one of these check fails an INTERP_KERNEL::Exception will be thrown.
 * So this method is safe but expensive if used to go through all data of \b this.
 */
int DataArrayInt::getIJSafe(int tupleId, int compoId) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(tupleId<0 || tupleId>=getNumberOfTuples())
    {
      std::ostringstream oss; oss << "DataArrayInt::getIJSafe : request for tupleId " << tupleId << " should be in [0," << getNumberOfTuples() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(compoId<0 || compoId>=getNumberOfComponents())
    {
      std::ostringstream oss; oss << "DataArrayInt::getIJSafe : request for compoId " << compoId << " should be in [0," << getNumberOfComponents() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return _mem[tupleId*((int)_info_on_compo.size())+compoId];
}

/*!
 * This method returns the last element in 'this'. So this method makes the hypothesis that 'this' is allocated.
 * This method works only for arrays that have exactly number of components equal to 1. If not an exception is thrown.
 * And to finish this method works for arrays that have number of tuples >= 1.
 */
int DataArrayInt::back() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::back : number of components not equal to one !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<1)
    throw INTERP_KERNEL::Exception("DataArrayInt::back : number of tuples must be >= 1 !");
  return *(getConstPointer()+nbOfTuples-1);
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

DataArrayIntIterator *DataArrayInt::iterator()
{
  return new DataArrayIntIterator(this);
}

DataArrayInt *DataArrayInt::getIdsEqual(int val) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::getIdsEqual : the array must have only one component, you can call 'rearrange' method before !");
  const int *cptr=getConstPointer();
  std::vector<int> res;
  int nbOfTuples=getNumberOfTuples();
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(*cptr==val)
      res.push_back(i);
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc((int)res.size(),1);
  std::copy(res.begin(),res.end(),ret->getPointer());
  return ret;
}

DataArrayInt *DataArrayInt::getIdsNotEqual(int val) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::getIdsNotEqual : the array must have only one component, you can call 'rearrange' method before !");
  const int *cptr=getConstPointer();
  std::vector<int> res;
  int nbOfTuples=getNumberOfTuples();
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(*cptr!=val)
      res.push_back(i);
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc((int)res.size(),1);
  std::copy(res.begin(),res.end(),ret->getPointer());
  return ret;
}

/*!
 * This method expects that 'this' is allocated. If not an exception will be thrown.
 * This method expect that the number of components is exactly equal to 1. If not an exception will be thrown.
 * For each element in 'this' equal to 'oldValue' will take the value 'newValue'.
 * @return number of elements impacted by the modification.
 */
int DataArrayInt::changeValue(int oldValue, int newValue) throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::changeValue : the array must have only one component, you can call 'rearrange' method before !");
  checkAllocated();
  int *start=getPointer();
  int *end2=start+getNbOfElems();
  int ret=0;
  for(int *val=start;val!=end2;val++)
    {
      if(*val==oldValue)
        {
          *val=newValue;
          ret++;
        }
    }
  return ret;
}

DataArrayInt *DataArrayInt::getIdsEqualList(const int *valsBg, const int *valsEnd) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::getIdsEqualList : the array must have only one component, you can call 'rearrange' method before !");
  std::set<int> vals2(valsBg,valsEnd);
  const int *cptr=getConstPointer();
  std::vector<int> res;
  int nbOfTuples=getNumberOfTuples();
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(vals2.find(*cptr)!=vals2.end())
      res.push_back(i);
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc((int)res.size(),1);
  std::copy(res.begin(),res.end(),ret->getPointer());
  return ret;
}

DataArrayInt *DataArrayInt::getIdsNotEqualList(const int *valsBg, const int *valsEnd) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::getIdsNotEqualList : the array must have only one component, you can call 'rearrange' method before !");
  std::set<int> vals2(valsBg,valsEnd);
  const int *cptr=getConstPointer();
  std::vector<int> res;
  int nbOfTuples=getNumberOfTuples();
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(vals2.find(*cptr)==vals2.end())
      res.push_back(i);
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc((int)res.size(),1);
  std::copy(res.begin(),res.end(),ret->getPointer());
  return ret;
}

/*!
 * This method is an extension of DataArrayInt::locateValue method because this method works for DataArrayInt with
 * any number of components excepted 0 (an INTERP_KERNEL::Exception is thrown in this case).
 * This method searches in \b this is there is a tuple that matched the input parameter \b tupl.
 * If any the tuple id is returned. If not -1 is returned.
 * 
 * This method throws an INTERP_KERNEL::Exception if the number of components in \b this mismatches with the size of
 * the input vector. An INTERP_KERNEL::Exception is thrown too if \b this is not allocated.
 *
 * \return tuple id where \b tupl is. -1 if no such tuple exists in \b this.
 * \sa DataArrayInt::search, DataArrayInt::presenceOfTuple.
 */
int DataArrayInt::locateTuple(const std::vector<int>& tupl) const throw(INTERP_KERNEL::Exception)
{
  int nbOfCompo=getNumberOfComponents();
  if(nbOfCompo==0)
    throw INTERP_KERNEL::Exception("DataArrayInt::locateTuple : 0 components in 'this' !");
  if(nbOfCompo!=(int)tupl.size())
    {
      std::ostringstream oss; oss << "DataArrayInt::locateTuple : 'this' contains " << nbOfCompo << " components and searching for a tuple of length " << tupl.size() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const int *cptr=getConstPointer();
  int nbOfVals=getNbOfElems();
  for(const int *work=cptr;work!=cptr+nbOfVals;)
    {
      work=std::search(work,cptr+nbOfVals,tupl.begin(),tupl.end());
      if(work!=cptr+nbOfVals)
        {
          if(std::distance(cptr,work)%nbOfCompo!=0)
            work++;
          else
            return std::distance(cptr,work)/nbOfCompo;
        }
    }
  return -1;
}

/*!
 * This method searches the sequence specified in input parameter \b vals in \b this.
 * This works only for DataArrayInt having number of components equal to one (if not an INTERP_KERNEL::Exception will be thrown).
 * This method differs from DataArrayInt::locateTuple in that the position is internal raw data is not considered here contrary to DataArrayInt::locateTuple.
 * \sa DataArrayInt::locateTuple
 */
int DataArrayInt::search(const std::vector<int>& vals) const throw(INTERP_KERNEL::Exception)
{
  int nbOfCompo=getNumberOfComponents();
  if(nbOfCompo!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::search : works only for DataArrayInt instance with one component !");
  const int *cptr=getConstPointer();
  int nbOfVals=getNbOfElems();
  const int *loc=std::search(cptr,cptr+nbOfVals,vals.begin(),vals.end());
  if(loc!=cptr+nbOfVals)
    return std::distance(cptr,loc);
  return -1;
}

/*!
 * This method expects to be called when number of components of this is equal to one.
 * This method returns the tuple id, if it exists, of the first tuple equal to \b value.
 * If not any tuple contains \b value -1 is returned.
 * \sa DataArrayInt::presenceOfValue
 */
int DataArrayInt::locateValue(int value) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::presenceOfValue : the array must have only one component, you can call 'rearrange' method before !");
  const int *cptr=getConstPointer();
  int nbOfTuples=getNumberOfTuples();
  const int *ret=std::find(cptr,cptr+nbOfTuples,value);
  if(ret!=cptr+nbOfTuples)
    return std::distance(cptr,ret);
  return -1;
}

/*!
 * This method expects to be called when number of components of this is equal to one.
 * This method returns the tuple id, if it exists, of the first tuple so that the value is contained in \b vals.
 * If not any tuple contains one of the values contained in 'vals' false is returned.
 * \sa DataArrayInt::presenceOfValue
 */
int DataArrayInt::locateValue(const std::vector<int>& vals) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::presenceOfValue : the array must have only one component, you can call 'rearrange' method before !");
  std::set<int> vals2(vals.begin(),vals.end());
  const int *cptr=getConstPointer();
  int nbOfTuples=getNumberOfTuples();
  for(const int *w=cptr;w!=cptr+nbOfTuples;w++)
    if(vals2.find(*w)!=vals2.end())
      return std::distance(cptr,w);
  return -1;
}

/*!
 * This method is an extension of DataArrayInt::presenceOfValue method because this method works for DataArrayInt with
 * any number of components excepted 0 (an INTERP_KERNEL::Exception is thrown in this case).
 * This method searches in \b this is there is a tuple that matched the input parameter \b tupl.
 * This method throws an INTERP_KERNEL::Exception if the number of components in \b this mismatches with the size of
 * the input vector. An INTERP_KERNEL::Exception is thrown too if \b this is not allocated.
 * \sa DataArrayInt::locateTuple
 */
bool DataArrayInt::presenceOfTuple(const std::vector<int>& tupl) const throw(INTERP_KERNEL::Exception)
{
  return locateTuple(tupl)!=-1;
}

/*!
 * This method expects to be called when number of components of this is equal to one.
 * This method returns true if it exists a tuple equal to \b value.
 * If not any tuple contains \b value false is returned.
 * \sa DataArrayInt::locateValue
 */
bool DataArrayInt::presenceOfValue(int value) const throw(INTERP_KERNEL::Exception)
{
  return locateValue(value)!=-1;
}

/*!
 * This method expects to be called when number of components of this is equal to one.
 * This method returns true if it exists a tuple so that the value is contained in \b vals.
 * If not any tuple contains one of the values contained in 'vals' false is returned.
 * \sa DataArrayInt::locateValue
 */
bool DataArrayInt::presenceOfValue(const std::vector<int>& vals) const throw(INTERP_KERNEL::Exception)
{
  return locateValue(vals)!=-1;
}

DataArrayInt *DataArrayInt::Aggregate(const DataArrayInt *a1, const DataArrayInt *a2, int offsetA2)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayInt::Aggregate : input DataArrayInt instance is NULL !");
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

DataArrayInt *DataArrayInt::Aggregate(const std::vector<const DataArrayInt *>& arr) throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayInt *> a;
  for(std::vector<const DataArrayInt *>::const_iterator it4=arr.begin();it4!=arr.end();it4++)
    if(*it4)
      a.push_back(*it4);
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
    throw INTERP_KERNEL::Exception("DataArrayInt::getMaxValue : must be applied on DataArrayInt with only one component !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<=0)
    throw INTERP_KERNEL::Exception("DataArrayInt::getMaxValue : array exists but number of tuples must be > 0 !");
  const int *vals=getConstPointer();
  const int *loc=std::max_element(vals,vals+nbOfTuples);
  tupleId=(int)std::distance(vals,loc);
  return *loc;
}

/*!
 * Idem to DataArrayInt::getMaxValue expect that here number of components can be >=1.
 */
int DataArrayInt::getMaxValueInArray() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  const int *loc=std::max_element(begin(),end());
  return *loc;
}

int DataArrayInt::getMinValue(int& tupleId) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::getMaxValue : must be applied on DataArrayInt with only one component !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<=0)
    throw INTERP_KERNEL::Exception("DataArrayInt::getMaxValue : array exists but number of tuples must be > 0 !");
  const int *vals=getConstPointer();
  const int *loc=std::min_element(vals,vals+nbOfTuples);
  tupleId=(int)std::distance(vals,loc);
  return *loc;
}

/*!
 * Idem to DataArrayInt::getMinValue expect that here number of components can be >=1.
 */
int DataArrayInt::getMinValueInArray() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  const int *loc=std::min_element(begin(),end());
  return *loc;
}

void DataArrayInt::abs() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int *ptr=getPointer();
  int nbOfElems=getNbOfElems();
  std::transform(ptr,ptr+nbOfElems,ptr,std::ptr_fun<int,int>(std::abs));
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

void DataArrayInt::applyLin(int a, int b) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int *ptr=getPointer();
  int nbOfElems=getNbOfElems();
  for(int i=0;i<nbOfElems;i++,ptr++)
    *ptr=a*(*ptr)+b;
  declareAsNew();
}

/*!
 * This method returns a newly allocated array containing the application of negate on \b this.
 * This method throws an INTERP_KERNEL::Exception if \b this is not allocated.
 */
DataArrayInt *DataArrayInt::negate() const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  DataArrayInt *newArr=DataArrayInt::New();
  int nbOfTuples=getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  newArr->alloc(nbOfTuples,nbOfComp);
  const int *cptr=getConstPointer();
  std::transform(cptr,cptr+nbOfTuples*nbOfComp,newArr->getPointer(),std::negate<int>());
  newArr->copyStringInfoFrom(*this);
  return newArr;
}

/*!
 * This method applies the operation 'numerator/x' for each element 'x' in 'this'.
 * If there is a value in 'this' exactly equal to 0. an exception is thrown.
 * Warning if presence of null this is modified for each values previous than place where exception was thrown !
 */
void DataArrayInt::applyInv(int numerator) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int *ptr=getPointer();
  int nbOfElems=getNbOfElems();
  for(int i=0;i<nbOfElems;i++,ptr++)
    {
      if(*ptr!=0)
        {
          *ptr=numerator/(*ptr);
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayInt::applyInv : presence of null value in tuple #" << i/getNumberOfComponents() << " component #" << i%getNumberOfComponents();
          oss << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  declareAsNew();
}

void DataArrayInt::applyDivideBy(int val) throw(INTERP_KERNEL::Exception)
{
  if(val==0)
    throw INTERP_KERNEL::Exception("DataArrayInt::applyDivideBy : Trying to divide by 0 !");
  checkAllocated();
  int *ptr=getPointer();
  int nbOfElems=getNbOfElems();
  std::transform(ptr,ptr+nbOfElems,ptr,std::bind2nd(std::divides<int>(),val));
  declareAsNew();
}

void DataArrayInt::applyModulus(int val) throw(INTERP_KERNEL::Exception)
{
  if(val<=0)
    throw INTERP_KERNEL::Exception("DataArrayInt::applyDivideBy : Trying to operate modulus on value <= 0 !");
  checkAllocated();
  int *ptr=getPointer();
  int nbOfElems=getNbOfElems();
  std::transform(ptr,ptr+nbOfElems,ptr,std::bind2nd(std::modulus<int>(),val));
  declareAsNew();
}

/*!
 * This method works only on data array with one component.
 * This method returns a newly allocated array storing stored ascendantly tuple ids in \b this so that
 * this[*id] in [\b vmin,\b vmax)
 * 
 * \param [in] vmin begin of range. This value is included in range.
 * \param [out] vmax end of range. This value is \b not included in range.
 * \return a newly allocated data array that the caller should deal with.
 */
DataArrayInt *DataArrayInt::getIdsInRange(int vmin, int vmax) const throw(INTERP_KERNEL::Exception)
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::getIdsInRange : this must have exactly one component !");
  const int *cptr=getConstPointer();
  std::vector<int> res;
  int nbOfTuples=getNumberOfTuples();
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(*cptr>=vmin && *cptr<vmax)
      res.push_back(i);
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc((int)res.size(),1);
  std::copy(res.begin(),res.end(),ret->getPointer());
  return ret;
}

/*!
 * This method applies the operation 'numerator%x' for each element 'x' in 'this'.
 * If there is a value in 'this' exactly equals or lower than 0. an exception is thrown.
 * Warning if presence of null this is modified for each values previous than place where exception was thrown !
 */
void DataArrayInt::applyRModulus(int val) throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  int *ptr=getPointer();
  int nbOfElems=getNbOfElems();
  for(int i=0;i<nbOfElems;i++,ptr++)
    {
      if(*ptr>0)
        {
          *ptr=val%(*ptr);
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayInt::applyRModulus : presence of value <=0 in tuple #" << i/getNumberOfComponents() << " component #" << i%getNumberOfComponents();
          oss << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  declareAsNew();
}

DataArrayInt *DataArrayInt::Meld(const DataArrayInt *a1, const DataArrayInt *a2) throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayInt *> arr(2);
  arr[0]=a1; arr[1]=a2;
  return Meld(arr);
}

DataArrayInt *DataArrayInt::Meld(const std::vector<const DataArrayInt *>& arr) throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayInt *> a;
  for(std::vector<const DataArrayInt *>::const_iterator it4=arr.begin();it4!=arr.end();it4++)
    if(*it4)
      a.push_back(*it4);
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
DataArrayInt *DataArrayInt::MakePartition(const std::vector<const DataArrayInt *>& groups, int newNb, std::vector< std::vector<int> >& fidsOfGroups) throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayInt *> groups2;
  for(std::vector<const DataArrayInt *>::const_iterator it4=groups.begin();it4!=groups.end();it4++)
    if(*it4)
      groups2.push_back(*it4);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(newNb,1);
  int *retPtr=ret->getPointer();
  std::fill(retPtr,retPtr+newNb,0);
  int fid=1;
  for(std::vector<const DataArrayInt *>::const_iterator iter=groups2.begin();iter!=groups2.end();iter++)
    {
      const int *ptr=(*iter)->getConstPointer();
      int nbOfElem=(*iter)->getNbOfElems();
      int sfid=fid;
      for(int j=0;j<sfid;j++)
        {
          bool found=false;
          for(int i=0;i<nbOfElem;i++)
            {
              if(ptr[i]>=0 && ptr[i]<newNb)
                {
                  if(retPtr[ptr[i]]==j)
                    {
                      retPtr[ptr[i]]=fid;
                      found=true;
                    }
                }
              else
                {
                  std::ostringstream oss; oss << "DataArrayInt::MakePartition : In group \"" << (*iter)->getName() << "\" in tuple #" << i << " value = " << ptr[i] << " ! Should be in [0," << newNb;
                  oss << ") !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
            }
          if(found)
            fid++;
        }
    }
  fidsOfGroups.clear();
  fidsOfGroups.resize(groups2.size());
  int grId=0;
  for(std::vector<const DataArrayInt *>::const_iterator iter=groups2.begin();iter!=groups2.end();iter++,grId++)
    {
      std::set<int> tmp;
      const int *ptr=(*iter)->getConstPointer();
      int nbOfElem=(*iter)->getNbOfElems();
      for(const int *p=ptr;p!=ptr+nbOfElem;p++)
        tmp.insert(retPtr[*p]);
      fidsOfGroups[grId].insert(fidsOfGroups[grId].end(),tmp.begin(),tmp.end());
    }
  ret->incrRef();
  return ret;
}

DataArrayInt *DataArrayInt::BuildUnion(const std::vector<const DataArrayInt *>& arr) throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayInt *> a;
  for(std::vector<const DataArrayInt *>::const_iterator it4=arr.begin();it4!=arr.end();it4++)
    if(*it4)
      a.push_back(*it4);
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
  ret->alloc((int)r.size(),1);
  std::copy(r.begin(),r.end(),ret->getPointer());
  return ret;
}

DataArrayInt *DataArrayInt::BuildIntersection(const std::vector<const DataArrayInt *>& arr) throw(INTERP_KERNEL::Exception)
{
  std::vector<const DataArrayInt *> a;
  for(std::vector<const DataArrayInt *>::const_iterator it4=arr.begin();it4!=arr.end();it4++)
    if(*it4)
      a.push_back(*it4);
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
  ret->alloc((int)r.size(),1);
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
   int nbOfRetVal=(int)std::count(tmp.begin(),tmp.end(),false);
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
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayInt::buildSubstraction : DataArrayInt pointer in input is NULL !");
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
  ret->alloc((int)r.size(),1);
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
  const int *ptr=getConstPointer();
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbOfTuples-1,1);
  int *out=ret->getPointer();
  std::transform(ptr+1,ptr+nbOfTuples,ptr,out,std::minus<int>());
  return ret;
}

/*!
 * This method performs the work on itself. This method works on array with number of component equal to one and allocated. If not an exception is thrown.
 * This method conserves the number of tuples and number of components (1). No reallocation is done.
 * For an array [3,5,1,2,0,8] it becomes [0,3,8,9,11,11]. For each i>0 new[i]=new[i-1]+old[i-1] for i==0 new[i]=0.
 * This could be usefull for allToAllV in MPI with contiguous policy.
 */
void DataArrayInt::computeOffsets() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
     throw INTERP_KERNEL::Exception("DataArrayInt::computeOffsets : only single component allowed !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples==0)
    return ;
  int *work=getPointer();
  int tmp=work[0];
  work[0]=0;
  for(int i=1;i<nbOfTuples;i++)
    {
      int tmp2=work[i];
      work[i]=work[i-1]+tmp;
      tmp=tmp2;
    }
  declareAsNew();
}

/*!
 * Idem DataArrayInt::computeOffsets method execpt that 'this' changes its number of tuples.
 * After the call in case of success new number of tuples is equal to old number of tuples +1.
 * The content in 'this' for the first old number of tuples is exactly the same than those given by
 * DataArrayInt::computeOffsets method.
 * For an array [3,5,1,2,0,8] it becomes [0,3,8,9,11,11,19].
 */
void DataArrayInt::computeOffsets2() throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::computeOffsets2 : only single component allowed !");
  int nbOfTuples=getNumberOfTuples();
  int *ret=new int[nbOfTuples+1];
  if(nbOfTuples==0)
    return ;
  const int *work=getConstPointer();
  ret[0]=0;
  for(int i=0;i<nbOfTuples;i++)
    ret[i+1]=work[i]+ret[i];
  useArray(ret,true,CPP_DEALLOC,nbOfTuples+1,1);
  declareAsNew();
}

/*!
 * This method works on array with number of component equal to one and allocated. If not an exception is thrown.
 * 'offsets' should be monotic ascendently. If not, an exception will be thrown.
 * This method retrives a newly created DataArrayInt instance with 1 component and this->getNumberOfTuples()-1 tuples.
 * If 'this' contains [0,2,3] and 'offsets' [0,3,6,10,14,20] the returned array will contain [0,1,2,6,7,8,9,10,11,12,13]
 */
DataArrayInt *DataArrayInt::buildExplicitArrByRanges(const DataArrayInt *offsets) const throw(INTERP_KERNEL::Exception)
{
  if(!offsets)
    throw INTERP_KERNEL::Exception("DataArrayInt::buildExplicitArrByRanges : DataArrayInt pointer in input is NULL !");
  checkAllocated();
  if(getNumberOfComponents()!=1)
     throw INTERP_KERNEL::Exception("DataArrayInt::buildExplicitArrByRanges : only single component allowed !");
  offsets->checkAllocated();
  if(offsets->getNumberOfComponents()!=1)
     throw INTERP_KERNEL::Exception("DataArrayInt::buildExplicitArrByRanges : input array should have only single component !");
  int othNbTuples=offsets->getNumberOfTuples()-1;
  int nbOfTuples=getNumberOfTuples();
  int retNbOftuples=0;
  const int *work=getConstPointer();
  const int *offPtr=offsets->getConstPointer();
  for(int i=0;i<nbOfTuples;i++)
    {
      int val=work[i];
      if(val>=0 && val<othNbTuples)
        {
          int delta=offPtr[val+1]-offPtr[val];
          if(delta>=0)
            retNbOftuples+=delta;
          else
            {
              std::ostringstream oss; oss << "DataArrayInt::buildExplicitArrByRanges : Tuple #" << val << " of offset array has a delta < 0 !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayInt::buildExplicitArrByRanges : Tuple #" << i << " in this contains " << val;
          oss << " whereas offsets array is of size " << othNbTuples+1 << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(retNbOftuples,1);
  int *retPtr=ret->getPointer();
  for(int i=0;i<nbOfTuples;i++)
    {
      int val=work[i];
      int start=offPtr[val];
      int off=offPtr[val+1]-start;
      for(int j=0;j<off;j++,retPtr++)
        *retPtr=start+j;
    }
  ret->incrRef();
  return ret;
}

/*!
 * Given in input ranges \a ranges, it returns a newly allocated DataArrayInt instance having one component and the same number of tuples than \a this.
 * For each tuple at place **i** in \a this it tells which is the first range in \a ranges that contains value \c this->getIJ(i,0) and put the result
 * in tuple **i** of returned DataArrayInt.
 * If ranges overlapped (in theory it should not) this method do not detect it and always returns the first range.
 *
 * For example if \a this contains : [1,24,7,8,10,17] and \a ranges contains [(0,3),(3,8),(8,15),(15,22),(22,30)]
 * The return DataArrayInt will contain : **[0,4,1,2,2,3]**
 * 
 * \param [in] ranges typically come from output of MEDCouplingUMesh::ComputeRangesFromTypeDistribution. Each range is specified like this : 1st component is
 *             for lower value included and 2nd component is the upper value of corresponding range **excluded**.
 * \throw If offsets is a null pointer or does not have 2 components or if \a this is not allocated or \a this do not have exactly one component. To finish an exception
 *        is thrown if no ranges in \a ranges contains value in \a this.
 * 
 * \sa DataArrayInt::findIdInRangeForEachTuple
 */
DataArrayInt *DataArrayInt::findRangeIdForEachTuple(const DataArrayInt *ranges) const throw(INTERP_KERNEL::Exception)
{
  if(!ranges)
    throw INTERP_KERNEL::Exception("DataArrayInt::findRangeIdForEachTuple : null input pointer !");
  if(ranges->getNumberOfComponents()!=2)
    throw INTERP_KERNEL::Exception("DataArrayInt::findRangeIdForEachTuple : input DataArrayInt instance should have 2 components !");
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::findRangeIdForEachTuple : this should have only one component !");
  int nbTuples=getNumberOfTuples();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(nbTuples,1);
  int nbOfRanges=ranges->getNumberOfTuples();
  const int *rangesPtr=ranges->getConstPointer();
  int *retPtr=ret->getPointer();
  const int *inPtr=getConstPointer();
  for(int i=0;i<nbTuples;i++,retPtr++)
    {
      int val=inPtr[i];
      bool found=false;
      for(int j=0;j<nbOfRanges && !found;j++)
        if(val>=rangesPtr[2*j] && val<rangesPtr[2*j+1])
          { *retPtr=j; found=true; }
      if(found)
        continue;
      else
        {
          std::ostringstream oss; oss << "DataArrayInt::findRangeIdForEachTuple : tuple #" << i << " not found by any ranges !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  ret->incrRef();
  return ret;
}

/*!
 * Given in input ranges \a ranges, it returns a newly allocated DataArrayInt instance having one component and the same number of tuples than \a this.
 * For each tuple at place **i** in \a this it tells which is the sub position of the first range in \a ranges that contains value \c this->getIJ(i,0) and put the result
 * in tuple **i** of returned DataArrayInt.
 * If ranges overlapped (in theory it should not) this method do not detect it and always returns the sub position of the first range.
 *
 * For example if \a this contains : [1,24,7,8,10,17] and \a ranges contains [(0,3),(3,8),(8,15),(15,22),(22,30)]
 * The return DataArrayInt will contain : **[1,2,4,0,2,2]**
 * This method is often called in pair with DataArrayInt::findRangeIdForEachTuple method.
 * 
 * \param [in] ranges typically come from output of MEDCouplingUMesh::ComputeRangesFromTypeDistribution. Each range is specified like this : 1st component is
 *             for lower value included and 2nd component is the upper value of corresponding range **excluded**.
 * \throw If offsets is a null pointer or does not have 2 components or if \a this is not allocated or \a this do not have exactly one component. To finish an exception
 *        is thrown if no ranges in \a ranges contains value in \a this.
 * \sa DataArrayInt::findRangeIdForEachTuple
 */
DataArrayInt *DataArrayInt::findIdInRangeForEachTuple(const DataArrayInt *ranges) const throw(INTERP_KERNEL::Exception)
{
  if(!ranges)
    throw INTERP_KERNEL::Exception("DataArrayInt::findIdInRangeForEachTuple : null input pointer !");
  if(ranges->getNumberOfComponents()!=2)
    throw INTERP_KERNEL::Exception("DataArrayInt::findIdInRangeForEachTuple : input DataArrayInt instance should have 2 components !");
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::findIdInRangeForEachTuple : this should have only one component !");
  int nbTuples=getNumberOfTuples();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(nbTuples,1);
  int nbOfRanges=ranges->getNumberOfTuples();
  const int *rangesPtr=ranges->getConstPointer();
  int *retPtr=ret->getPointer();
  const int *inPtr=getConstPointer();
  for(int i=0;i<nbTuples;i++,retPtr++)
    {
      int val=inPtr[i];
      bool found=false;
      for(int j=0;j<nbOfRanges && !found;j++)
        if(val>=rangesPtr[2*j] && val<rangesPtr[2*j+1])
          { *retPtr=val-rangesPtr[2*j]; found=true; }
      if(found)
        continue;
      else
        {
          std::ostringstream oss; oss << "DataArrayInt::findIdInRangeForEachTuple : tuple #" << i << " not found by any ranges !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  ret->incrRef();
  return ret;
}

/*!
 * 
 * \param [in] nbTimes specifies the nb of times each tuples in \a this will be duplicated contiguouly in returned DataArrayInt instance.
 *             \a nbTimes  should be at least equal to 1.
 * \return a newly allocated DataArrayInt having one component and number of tuples equal to \a nbTimes * \c this->getNumberOfTuples.
 * \throw if \a this is not allocated or if \a this has not number of components set to one or if \a nbTimes is lower than 1.
 */
DataArrayInt *DataArrayInt::duplicateEachTupleNTimes(int nbTimes) const throw(INTERP_KERNEL::Exception)
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::duplicateEachTupleNTimes : this should have only one component !");
  if(nbTimes<1)
    throw INTERP_KERNEL::Exception("DataArrayInt::duplicateEachTupleNTimes : nb times should be >= 1 !");
  int nbTuples=getNumberOfTuples();
  const int *inPtr=getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(nbTimes*nbTuples,1);
  int *retPtr=ret->getPointer();
  for(int i=0;i<nbTuples;i++,inPtr++)
    {
      int val=*inPtr;
      for(int j=0;j<nbTimes;j++,retPtr++)
        *retPtr=val;
    }
  ret->copyStringInfoFrom(*this);
  ret->incrRef();
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

DataArrayInt *DataArrayInt::Add(const DataArrayInt *a1, const DataArrayInt *a2) throw(INTERP_KERNEL::Exception)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayInt::Add : input DataArrayInt instance is NULL !");
  int nbOfTuple=a1->getNumberOfTuples();
  int nbOfTuple2=a2->getNumberOfTuples();
  int nbOfComp=a1->getNumberOfComponents();
  int nbOfComp2=a2->getNumberOfComponents();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=0;
  if(nbOfTuple==nbOfTuple2)
    {
      if(nbOfComp==nbOfComp2)
        {
          ret=DataArrayInt::New();
          ret->alloc(nbOfTuple,nbOfComp);
          std::transform(a1->begin(),a1->end(),a2->begin(),ret->getPointer(),std::plus<int>());
          ret->copyStringInfoFrom(*a1);
        }
      else
        {
          int nbOfCompMin,nbOfCompMax;
          const DataArrayInt *aMin, *aMax;
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
              ret=DataArrayInt::New();
              ret->alloc(nbOfTuple,nbOfCompMax);
              const int *aMinPtr=aMin->getConstPointer();
              const int *aMaxPtr=aMax->getConstPointer();
              int *res=ret->getPointer();
              for(int i=0;i<nbOfTuple;i++)
                res=std::transform(aMaxPtr+i*nbOfCompMax,aMaxPtr+(i+1)*nbOfCompMax,res,std::bind2nd(std::plus<int>(),aMinPtr[i]));
              ret->copyStringInfoFrom(*aMax);
            }
          else
            throw INTERP_KERNEL::Exception("Nb of components mismatch for array Add !");
        }
    }
  else if((nbOfTuple==1 && nbOfTuple2>1) || (nbOfTuple>1 && nbOfTuple2==1))
    {
      if(nbOfComp==nbOfComp2)
        {
          int nbOfTupleMax=std::max(nbOfTuple,nbOfTuple2);
          const DataArrayInt *aMin=nbOfTuple>nbOfTuple2?a2:a1;
          const DataArrayInt *aMax=nbOfTuple>nbOfTuple2?a1:a2;
          const int *aMinPtr=aMin->getConstPointer(),*aMaxPtr=aMax->getConstPointer();
          ret=DataArrayInt::New();
          ret->alloc(nbOfTupleMax,nbOfComp);
          int *res=ret->getPointer();
          for(int i=0;i<nbOfTupleMax;i++)
            res=std::transform(aMaxPtr+i*nbOfComp,aMaxPtr+(i+1)*nbOfComp,aMinPtr,res,std::plus<int>());
          ret->copyStringInfoFrom(*aMax);
        }
      else
        throw INTERP_KERNEL::Exception("Nb of components mismatch for array Add !");
    }
  else
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array Add !");
  ret->incrRef();
  return ret;
}

void DataArrayInt::addEqual(const DataArrayInt *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayInt::addEqual : input DataArrayInt instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayInt::addEqual  !";
  int nbOfTuple=getNumberOfTuples();
  int nbOfTuple2=other->getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  int nbOfComp2=other->getNumberOfComponents();
  if(nbOfTuple==nbOfTuple2)
    {
      if(nbOfComp==nbOfComp2)
        {
          std::transform(begin(),end(),other->begin(),getPointer(),std::plus<int>());
        }
      else if(nbOfComp2==1)
        {
          int *ptr=getPointer();
          const int *ptrc=other->getConstPointer();
          for(int i=0;i<nbOfTuple;i++)
            std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptr+i*nbOfComp,std::bind2nd(std::plus<int>(),*ptrc++));
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else if(nbOfTuple2==1)
    {
      if(nbOfComp2==nbOfComp)
        {
          int *ptr=getPointer();
          const int *ptrc=other->getConstPointer();
          for(int i=0;i<nbOfTuple;i++)
            std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptrc,ptr+i*nbOfComp,std::plus<int>());
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else
    throw INTERP_KERNEL::Exception(msg);
  declareAsNew();
}

DataArrayInt *DataArrayInt::Substract(const DataArrayInt *a1, const DataArrayInt *a2) throw(INTERP_KERNEL::Exception)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayInt::Substract : input DataArrayInt instance is NULL !");
  int nbOfTuple1=a1->getNumberOfTuples();
  int nbOfTuple2=a2->getNumberOfTuples();
  int nbOfComp1=a1->getNumberOfComponents();
  int nbOfComp2=a2->getNumberOfComponents();
  if(nbOfTuple2==nbOfTuple1)
    {
      if(nbOfComp1==nbOfComp2)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
          ret->alloc(nbOfTuple2,nbOfComp1);
          std::transform(a1->begin(),a1->end(),a2->begin(),ret->getPointer(),std::minus<int>());
          ret->copyStringInfoFrom(*a1);
          ret->incrRef();
          return ret;
        }
      else if(nbOfComp2==1)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
          ret->alloc(nbOfTuple1,nbOfComp1);
          const int *a2Ptr=a2->getConstPointer();
          const int *a1Ptr=a1->getConstPointer();
          int *res=ret->getPointer();
          for(int i=0;i<nbOfTuple1;i++)
            res=std::transform(a1Ptr+i*nbOfComp1,a1Ptr+(i+1)*nbOfComp1,res,std::bind2nd(std::minus<int>(),a2Ptr[i]));
          ret->copyStringInfoFrom(*a1);
          ret->incrRef();
          return ret;
        }
      else
        {
          a1->checkNbOfComps(nbOfComp2,"Nb of components mismatch for array Substract !");
          return 0;
        }
    }
  else if(nbOfTuple2==1)
    {
      a1->checkNbOfComps(nbOfComp2,"Nb of components mismatch for array Substract !");
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
      ret->alloc(nbOfTuple1,nbOfComp1);
      const int *a1ptr=a1->getConstPointer(),*a2ptr=a2->getConstPointer();
      int *pt=ret->getPointer();
      for(int i=0;i<nbOfTuple1;i++)
        pt=std::transform(a1ptr+i*nbOfComp1,a1ptr+(i+1)*nbOfComp1,a2ptr,pt,std::minus<int>());
      ret->copyStringInfoFrom(*a1);
      ret->incrRef();
      return ret;
    }
  else
    {
      a1->checkNbOfTuples(nbOfTuple2,"Nb of tuples mismatch for array Substract !");//will always throw an exception
      return 0;
    }
}

void DataArrayInt::substractEqual(const DataArrayInt *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayInt::substractEqual : input DataArrayInt instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayInt::substractEqual  !";
  int nbOfTuple=getNumberOfTuples();
  int nbOfTuple2=other->getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  int nbOfComp2=other->getNumberOfComponents();
  if(nbOfTuple==nbOfTuple2)
    {
      if(nbOfComp==nbOfComp2)
        {
          std::transform(begin(),end(),other->begin(),getPointer(),std::minus<int>());
        }
      else if(nbOfComp2==1)
        {
          int *ptr=getPointer();
          const int *ptrc=other->getConstPointer();
          for(int i=0;i<nbOfTuple;i++)
            std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptr+i*nbOfComp,std::bind2nd(std::minus<int>(),*ptrc++));
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else if(nbOfTuple2==1)
    {
      int *ptr=getPointer();
      const int *ptrc=other->getConstPointer();
      for(int i=0;i<nbOfTuple;i++)
        std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptrc,ptr+i*nbOfComp,std::minus<int>());
    }
  else
    throw INTERP_KERNEL::Exception(msg);
  declareAsNew();
}

DataArrayInt *DataArrayInt::Multiply(const DataArrayInt *a1, const DataArrayInt *a2) throw(INTERP_KERNEL::Exception)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayInt::Multiply : input DataArrayInt instance is NULL !");
  int nbOfTuple=a1->getNumberOfTuples();
  int nbOfTuple2=a2->getNumberOfTuples();
  int nbOfComp=a1->getNumberOfComponents();
  int nbOfComp2=a2->getNumberOfComponents();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=0;
  if(nbOfTuple==nbOfTuple2)
    {
      if(nbOfComp==nbOfComp2)
        {
          ret=DataArrayInt::New();
          ret->alloc(nbOfTuple,nbOfComp);
          std::transform(a1->begin(),a1->end(),a2->begin(),ret->getPointer(),std::multiplies<int>());
          ret->copyStringInfoFrom(*a1);
        }
      else
        {
          int nbOfCompMin,nbOfCompMax;
          const DataArrayInt *aMin, *aMax;
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
              ret=DataArrayInt::New();
              ret->alloc(nbOfTuple,nbOfCompMax);
              const int *aMinPtr=aMin->getConstPointer();
              const int *aMaxPtr=aMax->getConstPointer();
              int *res=ret->getPointer();
              for(int i=0;i<nbOfTuple;i++)
                res=std::transform(aMaxPtr+i*nbOfCompMax,aMaxPtr+(i+1)*nbOfCompMax,res,std::bind2nd(std::multiplies<int>(),aMinPtr[i]));
              ret->copyStringInfoFrom(*aMax);
            }
          else
            throw INTERP_KERNEL::Exception("Nb of components mismatch for array Multiply !");
        }
    }
  else if((nbOfTuple==1 && nbOfTuple2>1) || (nbOfTuple>1 && nbOfTuple2==1))
    {
      if(nbOfComp==nbOfComp2)
        {
          int nbOfTupleMax=std::max(nbOfTuple,nbOfTuple2);
          const DataArrayInt *aMin=nbOfTuple>nbOfTuple2?a2:a1;
          const DataArrayInt *aMax=nbOfTuple>nbOfTuple2?a1:a2;
          const int *aMinPtr=aMin->getConstPointer(),*aMaxPtr=aMax->getConstPointer();
          ret=DataArrayInt::New();
          ret->alloc(nbOfTupleMax,nbOfComp);
          int *res=ret->getPointer();
          for(int i=0;i<nbOfTupleMax;i++)
            res=std::transform(aMaxPtr+i*nbOfComp,aMaxPtr+(i+1)*nbOfComp,aMinPtr,res,std::multiplies<int>());
          ret->copyStringInfoFrom(*aMax);
        }
      else
        throw INTERP_KERNEL::Exception("Nb of components mismatch for array Multiply !");
    }
  else
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array Multiply !");
  ret->incrRef();
  return ret;
}

void DataArrayInt::multiplyEqual(const DataArrayInt *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayInt::multiplyEqual : input DataArrayInt instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayInt::multiplyEqual !";
  int nbOfTuple=getNumberOfTuples();
  int nbOfTuple2=other->getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  int nbOfComp2=other->getNumberOfComponents();
  if(nbOfTuple==nbOfTuple2)
    {
      if(nbOfComp==nbOfComp2)
        {
          std::transform(begin(),end(),other->begin(),getPointer(),std::multiplies<int>());
        }
      else if(nbOfComp2==1)
        {
          int *ptr=getPointer();
          const int *ptrc=other->getConstPointer();
          for(int i=0;i<nbOfTuple;i++)
            std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptr+i*nbOfComp,std::bind2nd(std::multiplies<int>(),*ptrc++));    
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else if(nbOfTuple2==1)
    {
      if(nbOfComp2==nbOfComp)
        {
          int *ptr=getPointer();
          const int *ptrc=other->getConstPointer();
          for(int i=0;i<nbOfTuple;i++)
            std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptrc,ptr+i*nbOfComp,std::multiplies<int>());
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else
    throw INTERP_KERNEL::Exception(msg);
  declareAsNew();
}

DataArrayInt *DataArrayInt::Divide(const DataArrayInt *a1, const DataArrayInt *a2) throw(INTERP_KERNEL::Exception)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayInt::Divide : input DataArrayInt instance is NULL !");
  int nbOfTuple1=a1->getNumberOfTuples();
  int nbOfTuple2=a2->getNumberOfTuples();
  int nbOfComp1=a1->getNumberOfComponents();
  int nbOfComp2=a2->getNumberOfComponents();
  if(nbOfTuple2==nbOfTuple1)
    {
      if(nbOfComp1==nbOfComp2)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
          ret->alloc(nbOfTuple2,nbOfComp1);
          std::transform(a1->begin(),a1->end(),a2->begin(),ret->getPointer(),std::divides<int>());
          ret->copyStringInfoFrom(*a1);
          ret->incrRef();
          return ret;
        }
      else if(nbOfComp2==1)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
          ret->alloc(nbOfTuple1,nbOfComp1);
          const int *a2Ptr=a2->getConstPointer();
          const int *a1Ptr=a1->getConstPointer();
          int *res=ret->getPointer();
          for(int i=0;i<nbOfTuple1;i++)
            res=std::transform(a1Ptr+i*nbOfComp1,a1Ptr+(i+1)*nbOfComp1,res,std::bind2nd(std::divides<int>(),a2Ptr[i]));
          ret->copyStringInfoFrom(*a1);
          ret->incrRef();
          return ret;
        }
      else
        {
          a1->checkNbOfComps(nbOfComp2,"Nb of components mismatch for array Divide !");
          return 0;
        }
    }
  else if(nbOfTuple2==1)
    {
      a1->checkNbOfComps(nbOfComp2,"Nb of components mismatch for array Divide !");
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
      ret->alloc(nbOfTuple1,nbOfComp1);
      const int *a1ptr=a1->getConstPointer(),*a2ptr=a2->getConstPointer();
      int *pt=ret->getPointer();
      for(int i=0;i<nbOfTuple1;i++)
        pt=std::transform(a1ptr+i*nbOfComp1,a1ptr+(i+1)*nbOfComp1,a2ptr,pt,std::divides<int>());
      ret->copyStringInfoFrom(*a1);
      ret->incrRef();
      return ret;
    }
  else
    {
      a1->checkNbOfTuples(nbOfTuple2,"Nb of tuples mismatch for array Divide !");//will always throw an exception
      return 0;
    }
}

void DataArrayInt::divideEqual(const DataArrayInt *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayInt::divideEqual : input DataArrayInt instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayInt::divideEqual !";
  int nbOfTuple=getNumberOfTuples();
  int nbOfTuple2=other->getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  int nbOfComp2=other->getNumberOfComponents();
  if(nbOfTuple==nbOfTuple2)
    {
      if(nbOfComp==nbOfComp2)
        {
          std::transform(begin(),end(),other->begin(),getPointer(),std::divides<int>());
        }
      else if(nbOfComp2==1)
        {
          int *ptr=getPointer();
          const int *ptrc=other->getConstPointer();
          for(int i=0;i<nbOfTuple;i++)
            std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptr+i*nbOfComp,std::bind2nd(std::divides<int>(),*ptrc++));
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else if(nbOfTuple2==1)
    {
      if(nbOfComp2==nbOfComp)
        {
          int *ptr=getPointer();
          const int *ptrc=other->getConstPointer();
          for(int i=0;i<nbOfTuple;i++)
            std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptrc,ptr+i*nbOfComp,std::divides<int>());
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else
    throw INTERP_KERNEL::Exception(msg);
  declareAsNew();
}

DataArrayInt *DataArrayInt::Modulus(const DataArrayInt *a1, const DataArrayInt *a2) throw(INTERP_KERNEL::Exception)
{
    if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayInt::Modulus : input DataArrayInt instance is NULL !");
  int nbOfTuple1=a1->getNumberOfTuples();
  int nbOfTuple2=a2->getNumberOfTuples();
  int nbOfComp1=a1->getNumberOfComponents();
  int nbOfComp2=a2->getNumberOfComponents();
  if(nbOfTuple2==nbOfTuple1)
    {
      if(nbOfComp1==nbOfComp2)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
          ret->alloc(nbOfTuple2,nbOfComp1);
          std::transform(a1->begin(),a1->end(),a2->begin(),ret->getPointer(),std::modulus<int>());
          ret->copyStringInfoFrom(*a1);
          ret->incrRef();
          return ret;
        }
      else if(nbOfComp2==1)
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
          ret->alloc(nbOfTuple1,nbOfComp1);
          const int *a2Ptr=a2->getConstPointer();
          const int *a1Ptr=a1->getConstPointer();
          int *res=ret->getPointer();
          for(int i=0;i<nbOfTuple1;i++)
            res=std::transform(a1Ptr+i*nbOfComp1,a1Ptr+(i+1)*nbOfComp1,res,std::bind2nd(std::modulus<int>(),a2Ptr[i]));
          ret->copyStringInfoFrom(*a1);
          ret->incrRef();
          return ret;
        }
      else
        {
          a1->checkNbOfComps(nbOfComp2,"Nb of components mismatch for array Modulus !");
          return 0;
        }
    }
  else if(nbOfTuple2==1)
    {
      a1->checkNbOfComps(nbOfComp2,"Nb of components mismatch for array Modulus !");
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
      ret->alloc(nbOfTuple1,nbOfComp1);
      const int *a1ptr=a1->getConstPointer(),*a2ptr=a2->getConstPointer();
      int *pt=ret->getPointer();
      for(int i=0;i<nbOfTuple1;i++)
        pt=std::transform(a1ptr+i*nbOfComp1,a1ptr+(i+1)*nbOfComp1,a2ptr,pt,std::modulus<int>());
      ret->copyStringInfoFrom(*a1);
      ret->incrRef();
      return ret;
    }
  else
    {
      a1->checkNbOfTuples(nbOfTuple2,"Nb of tuples mismatch for array Modulus !");//will always throw an exception
      return 0;
    }
}

void DataArrayInt::modulusEqual(const DataArrayInt *other) throw(INTERP_KERNEL::Exception)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayInt::modulusEqual : input DataArrayInt instance is NULL !");
  const char *msg="Nb of tuples mismatch for DataArrayInt::modulusEqual !";
  int nbOfTuple=getNumberOfTuples();
  int nbOfTuple2=other->getNumberOfTuples();
  int nbOfComp=getNumberOfComponents();
  int nbOfComp2=other->getNumberOfComponents();
  if(nbOfTuple==nbOfTuple2)
    {
      if(nbOfComp==nbOfComp2)
        {
          std::transform(begin(),end(),other->begin(),getPointer(),std::modulus<int>());
        }
      else if(nbOfComp2==1)
        {
          if(nbOfComp2==nbOfComp)
            {
              int *ptr=getPointer();
              const int *ptrc=other->getConstPointer();
              for(int i=0;i<nbOfTuple;i++)
                std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptr+i*nbOfComp,std::bind2nd(std::modulus<int>(),*ptrc++));
            }
          else
            throw INTERP_KERNEL::Exception(msg);
        }
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  else if(nbOfTuple2==1)
    {
      int *ptr=getPointer();
      const int *ptrc=other->getConstPointer();
      for(int i=0;i<nbOfTuple;i++)
        std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptrc,ptr+i*nbOfComp,std::modulus<int>());
    }
  else
    throw INTERP_KERNEL::Exception(msg);
  declareAsNew();
}

int *DataArrayInt::CheckAndPreparePermutation(const int *start, const int *end)
{
  std::size_t sz=std::distance(start,end);
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
    *iter2=(int)std::distance(work,std::find(work,work+sz,*iter));
  delete [] work;
  return ret;
}

DataArrayInt *DataArrayInt::Range(int begin, int end, int step) throw(INTERP_KERNEL::Exception)
{
  int nbOfTuples=GetNumberOfItemGivenBESRelative(begin,end,step,"DataArrayInt::Range");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbOfTuples,1);
  int *ptr=ret->getPointer();
  if(step>0)
    {
      for(int i=begin;i<end;i+=step,ptr++)
        *ptr=i;
    }
  else
    {
      for(int i=begin;i>end;i+=step,ptr++)
        *ptr=i;
    }
  ret->incrRef();
  return ret;
}

/*!
 * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
 * Server side.
 */
void DataArrayInt::getTinySerializationIntInformation(std::vector<int>& tinyInfo) const
{
  tinyInfo.resize(2);
  if(isAllocated())
    {
      tinyInfo[0]=getNumberOfTuples();
      tinyInfo[1]=getNumberOfComponents();
    }
  else
    {
      tinyInfo[0]=-1;
      tinyInfo[1]=-1;
    }
}

/*!
 * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
 * Server side.
 */
void DataArrayInt::getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const
{
  if(isAllocated())
    {
      int nbOfCompo=getNumberOfComponents();
      tinyInfo.resize(nbOfCompo+1);
      tinyInfo[0]=getName();
      for(int i=0;i<nbOfCompo;i++)
        tinyInfo[i+1]=getInfoOnComponent(i);
    }
  else
    {
      tinyInfo.resize(1);
      tinyInfo[0]=getName();
    }
}

/*!
 * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
 * This method returns if a feeding is needed.
 */
bool DataArrayInt::resizeForUnserialization(const std::vector<int>& tinyInfoI)
{
  int nbOfTuple=tinyInfoI[0];
  int nbOfComp=tinyInfoI[1];
  if(nbOfTuple!=-1 || nbOfComp!=-1)
    {
      alloc(nbOfTuple,nbOfComp);
      return true;
    }
  return false;
}

/*!
 * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
 * This method returns if a feeding is needed.
 */
void DataArrayInt::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<std::string>& tinyInfoS)
{
  setName(tinyInfoS[0].c_str());
  if(isAllocated())
    {
      int nbOfCompo=getNumberOfComponents();
      for(int i=0;i<nbOfCompo;i++)
        setInfoOnComponent(i,tinyInfoS[i+1].c_str());
    }
}

DataArrayIntIterator::DataArrayIntIterator(DataArrayInt *da):_da(da),_pt(0),_tuple_id(0),_nb_comp(0),_nb_tuple(0)
{
  if(_da)
    {
      _da->incrRef();
      if(_da->isAllocated())
        {
          _nb_comp=da->getNumberOfComponents();
          _nb_tuple=da->getNumberOfTuples();
          _pt=da->getPointer();
        }
    }
}

DataArrayIntIterator::~DataArrayIntIterator()
{
  if(_da)
    _da->decrRef();
}

DataArrayIntTuple *DataArrayIntIterator::nextt()
{
  if(_tuple_id<_nb_tuple)
    {
      _tuple_id++;
      DataArrayIntTuple *ret=new DataArrayIntTuple(_pt,_nb_comp);
      _pt+=_nb_comp;
      return ret;
    }
  else
    return 0;
}

DataArrayIntTuple::DataArrayIntTuple(int *pt, int nbOfComp):_pt(pt),_nb_of_compo(nbOfComp)
{
}

std::string DataArrayIntTuple::repr() const
{
  std::ostringstream oss; oss << "(";
  for(int i=0;i<_nb_of_compo-1;i++)
    oss << _pt[i] << ", ";
  oss << _pt[_nb_of_compo-1] << ")";
  return oss.str();
}

int DataArrayIntTuple::intValue() const throw(INTERP_KERNEL::Exception)
{
  if(_nb_of_compo==1)
    return *_pt;
  throw INTERP_KERNEL::Exception("DataArrayIntTuple::intValue : DataArrayIntTuple instance has not exactly 1 component -> Not possible to convert it into an integer !");
}

/*!
 * This method returns a newly allocated instance the caller should dealed with by a ParaMEDMEM::DataArrayInt::decrRef.
 * This method performs \b no copy of data. The content is only referenced using ParaMEDMEM::DataArrayInt::useArray with ownership set to \b false.
 * This method throws an INTERP_KERNEL::Exception is it is impossible to match sizes of \b this that is too say \b nbOfCompo=this->_nb_of_elem and \bnbOfTuples==1 or
 * \b nbOfCompo=1 and \bnbOfTuples==this->_nb_of_elem.
 */
DataArrayInt *DataArrayIntTuple::buildDAInt(int nbOfTuples, int nbOfCompo) const throw(INTERP_KERNEL::Exception)
{
  if((_nb_of_compo==nbOfCompo && nbOfTuples==1) || (_nb_of_compo==nbOfTuples && nbOfCompo==1))
    {
      DataArrayInt *ret=DataArrayInt::New();
      ret->useArray(_pt,false,CPP_DEALLOC,nbOfTuples,nbOfCompo);
      return ret;
    }
  else
    {
      std::ostringstream oss; oss << "DataArrayIntTuple::buildDAInt : unable to build a requested DataArrayInt instance with nbofTuple=" << nbOfTuples << " and nbOfCompo=" << nbOfCompo;
      oss << ".\nBecause the number of elements in this is " << _nb_of_compo << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}
